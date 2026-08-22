// examples/example_reuse.cpp
// The three phases have three different lifetimes, and knowing which is which is most of what the
// split API is for.
//
//   analyze   depends on the PATTERN of A. Ordering, elimination forest, symbolic factorization,
//             all of them graph algorithms that never read a value
//   factor    depends on the VALUES. It reuses the analysis and produces the numeric factor
//   solve     depends on the RIGHT-HAND SIDE. It reuses the factor
//
// So a matrix whose structure is fixed and whose numbers change, which is what a Newton iteration
// and a time-stepping loop both produce, analyzes once and factors per step. And a step with
// several right-hand sides factors once and solves per side. This example is that loop: one
// pattern, three value sets, three right-hand sides each, and it counts the phases to show that
// nine solves cost one analysis.
//
// **The saving is not marginal.** benchmarks/pipeline measures analysis at 27 to 40 percent of one
// analyze-plus-factor on grid Laplacians, with the ordering about half of that, so a caller who
// re-analyzes per step is paying roughly a third more for every step after the first.
//
// The second half of the file is the other side of the same fact: which settings invalidate an
// analysis and which do not. A setting that changes the pattern's treatment throws the analysis
// away; one that changes only the arithmetic does not. The facade tracks this so a caller does not
// have to, and analyzed() reports it.
//
// Compile (macOS, from repo root):
//   g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude examples/example_reuse.cpp src/*.cpp -framework Accelerate -o example_reuse_cpp
// Linux: replace `-framework Accelerate` with `-llapack -lblas`.

#include "oblio/DirectSolver.h"

#include <cstdio>
#include <set>
#include <vector>

using namespace Oblio;

namespace {

using Val = double;

// A five-point grid Laplacian whose diagonal carries a parameter. Every value set below shares one
// pattern: the diagonal is stored whatever it holds, so changing it moves no structural entry, and
// that is exactly the condition under which an analysis can be kept.
SparseMatrix<Val> grid(std::size_t side, Val diagonal) {
    const std::size_t size = side * side;
    std::vector<std::set<std::int32_t>> adjacency(size);
    for (std::size_t r = 0; r < side; ++r) {
        for (std::size_t c = 0; c < side; ++c) {
            const std::int32_t v = static_cast<std::int32_t>(r * side + c);
            adjacency[v].insert(v);
            if (r > 0)        adjacency[v].insert(v - static_cast<std::int32_t>(side));
            if (r + 1 < side) adjacency[v].insert(v + static_cast<std::int32_t>(side));
            if (c > 0)        adjacency[v].insert(v - 1);
            if (c + 1 < side) adjacency[v].insert(v + 1);
        }
    }
    std::vector<std::size_t>  colPtr(size + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<Val>          val;
    for (std::size_t j = 0; j < size; ++j) {
        for (std::int32_t i : adjacency[j]) {
            rowIdx.push_back(i);
            val.push_back(i == static_cast<std::int32_t>(j) ? diagonal : Val(-1));
        }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<Val>(size, colPtr, rowIdx, val);
}

} // namespace

int main() {
    const std::size_t side = 20;
    const std::size_t size = side * side;

    // Three matrices sharing one pattern, as three steps of an iteration would produce.
    const Val diagonal[] = { Val(4.0), Val(4.5), Val(5.0) };

    // Three right-hand sides, as one step with several loads would.
    std::vector<Vector<Val>> b;
    for (int k = 0; k < 3; ++k) {
        Vector<Val> v(size);
        for (std::size_t i = 0; i < size; ++i)
            v[i] = Val(1) + Val(k) * Val(0.1) * static_cast<Val>(i % 7);
        b.push_back(v);
    }

    int analyses = 0, factorizations = 0, solves = 0;
    double worst = 0;

    DirectSolver<Val> solver(Ordering::MMD3, Factorization::Cholesky, Traversal::Multifrontal);

    printf("Grid %zux%zu, %zu unknowns. Three value sets, three right-hand sides each.\n\n",
           side, side, size);
    printf("  %-6s  %-8s  %-9s  %-6s  %s\n", "step", "diagonal", "analyzed?", "rhs", "residual");
    printf("  %-6s  %-8s  %-9s  %-6s  %s\n", "----", "--------", "---------", "---", "--------");

    for (int step = 0; step < 3; ++step) {
        const SparseMatrix<Val> A = grid(side, diagonal[step]);

        // Analyze only when there is no analysis to reuse. After the first step analyzed() is
        // still true, because nothing has invalidated it: factor and solve do not, and the
        // matrix, though it holds different numbers, has the same pattern.
        const bool hadAnalysis = solver.analyzed();
        if (!hadAnalysis) {
            if (!solver.analyze(A)) { printf("  analyze failed\n"); return 1; }
            ++analyses;
        }

        if (!solver.factor(A)) { printf("  factor failed\n"); return 1; }
        ++factorizations;

        for (int k = 0; k < 3; ++k) {
            Vector<Val> x(size);
            if (!solver.solve(b[k], x)) { printf("  solve failed\n"); return 1; }
            ++solves;

            const double residual = solver.relativeResidual(A, b[k], x);
            if (residual > worst) worst = residual;
            printf("  %-6d  %-8.1f  %-9s  %-6d  %.2e\n",
                   step, diagonal[step], hadAnalysis ? "reused" : "computed", k, residual);
        }
    }

    printf("\n  %d analysis, %d factorizations, %d solves. Worst residual %.2e.\n",
           analyses, factorizations, solves, worst);
    printf("  The analysis is computed once and reused twice; each factor serves three solves.\n");

    // ---------------------------------------------------------------------------------------
    // Which settings keep an analysis and which throw it away. The rule is one line: a setting
    // that changes what the structural phases compute invalidates it, and one that changes only
    // the arithmetic does not.
    // ---------------------------------------------------------------------------------------
    printf("\nWhat each setting does to a completed analysis:\n\n");
    printf("  %-42s  %s\n", "setting", "analysis");
    printf("  %-42s  %s\n", "-------", "--------");

    const SparseMatrix<Val> A = grid(side, diagonal[0]);
    const auto report = [&](const char* what) {
        printf("  %-42s  %s\n", what, solver.analyzed() ? "kept" : "discarded");
    };

    solver.setTraversal(Traversal::LeftLooking);   // out of multifrontal: the forest differs
    solver.analyze(A);

    // The factorization changes the arithmetic and nothing structural, so the analysis stands.
    // This is the one that matters in practice: switching from Cholesky to an indefinite method
    // after a step goes wrong costs a factorization, not an analysis.
    solver.setFactorization(Factorization::DynamicLDLT);
    report("setFactorization: Cholesky -> DynamicLDLT");

    // Left and right looking read the same forest, so the analysis stands for them too.
    solver.setTraversal(Traversal::RightLooking);
    report("setTraversal: LeftLooking -> RightLooking");

    // Multifrontal is the exception. It sorts each supernode's children and relabels them into a
    // postorder, so the forest itself differs and the analysis cannot be carried across.
    solver.setTraversal(Traversal::Multifrontal);
    report("setTraversal: RightLooking -> Multifrontal");

    solver.analyze(A);
    solver.setOrdering(Ordering::AMD3);
    report("setOrdering: MMD3 -> AMD3");

    solver.analyze(A);
    solver.setSupernodes(Supernodes::Nodal);
    report("setSupernodes: Fundamental -> Nodal");

    solver.analyze(A);
    solver.setAmalgamation(8);
    report("setAmalgamation: none -> 8");

    printf("\n  So the ordering and the two forest settings discard it, since the forest is what\n"
           "  they change and the analysis is the forest. The factorization never does. The\n"
           "  traversal does only across the multifrontal boundary, which is the one case where\n"
           "  a traversal changes the forest rather than merely reading it.\n");

    return 0;
}
