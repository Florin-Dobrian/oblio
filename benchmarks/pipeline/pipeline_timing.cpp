// pipeline_timing.cpp - where a solve's time actually goes, per ordering and per traversal.
//
// `benchmarks/ordering` measures one phase against itself: how long each ordering takes and how
// much it fills. This measures the phases against each other, which is a different question and
// the one that says whether the first question matters. An ordering that saves a millisecond of
// analysis and costs two per factorization is a loss for anyone who factors twice.
//
//   make run                        the default grid sides
//   ./pipeline_timing_cpp 200       any sides
//
// -O3 -DNDEBUG, best of three after a warm-up, as in the ordering folder and for the same reason.
//
// **Grid Laplacians only, and that is a real limitation rather than a starting point.** A grid is
// where the MMD fill gap is largest, because nearly every live vertex has the same degree and the
// tie-break decides almost every pick; `experiments/ordering` measured a 16 percent spread in fill
// across four filing orders of one algorithm. So a break-even computed here may be a property of
// grids rather than of the orderings, and widening the matrices is the first thing this folder
// needs. Cholesky, real, left-looking against right-looking against multifrontal.

#include "oblio/DirectSolver.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/OrderEngine.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"
#include "oblio/SymFactor.h"
#include "oblio/SymFactorEngine.h"
#include "oblio/Types.h"
#include "oblio/Vector.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace Oblio;

static SparseMatrix<double> grid2D(int m) {
    const int n = m * m;
    std::vector<std::vector<int>> adj(n);
    auto id = [&](int r, int c) { return r * m + c; };
    for (int r = 0; r < m; ++r)
        for (int c = 0; c < m; ++c) {
            const int u = id(r, c);
            if (r > 0)     adj[u].push_back(id(r - 1, c));
            if (c > 0)     adj[u].push_back(id(r, c - 1));
            if (r + 1 < m) adj[u].push_back(id(r + 1, c));
            if (c + 1 < m) adj[u].push_back(id(r, c + 1));
        }
    std::vector<std::size_t>  colPtr(n + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<double>       val;
    for (int j = 0; j < n; ++j) {
        std::vector<int> col = adj[j];
        col.push_back(j);
        std::sort(col.begin(), col.end());
        for (int i : col) {
            rowIdx.push_back(static_cast<std::int32_t>(i));
            val.push_back(i == j ? 4.0 : -1.0);
        }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<double>(static_cast<std::size_t>(n), colPtr, rowIdx, val);
}

// Best of three after a warm-up. The work is a callable so every phase is timed the same way.
template <class Work>
static double bestOfThree(Work work) {
    work();
    double best = 1e30;
    for (int trial = 0; trial < 3; ++trial) {
        const auto t0 = std::chrono::steady_clock::now();
        work();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

// nnz(L), from the symbolic factor rather than from any ordering's estimate: a supernode
// contributes its own triangle plus its update rows. Identical to the ordering folder's, and
// deliberately duplicated rather than shared, since a benchmark should stand alone.
static std::size_t fill(const SparseMatrix<double>& A, Ordering method) {
    const OrderEngine oe(method);
    Permutation       P;
    if (!oe.compute(A, P)) return 0;
    const ElmForestEngine fe;
    ElmForest             ef;
    if (!fe.compute(A, P, ef)) return 0;
    const SymFactorEngine se;
    SymFactor             sf;
    if (!se.compute(A, P, ef, sf)) return 0;

    std::size_t nnz = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.snodeSize()); ++kk) {
        const std::size_t f = sf.frontSize(kk), u = sf.updateSize(kk);
        nnz += f * (f + 1) / 2 + f * u;
    }
    return nnz;
}

struct Row {
    std::string name;
    double      order = 0;      // OrderEngine::compute alone
    double      analyze = 0;    // the whole facade call, left-looking
    double      analyzeMf = 0;  // the same, multifrontal: it sorts children and relabels
    double      factLl = 0;
    double      factRl = 0;
    double      factMf = 0;
    double      solve = 0;
    std::size_t nnzL = 0;
};

// One traversal's factor time, with the analysis already done. Reported separately because
// analyze depends on the traversal (multifrontal orders the children and relabels the supernodes)
// while factor depends on it far more.
static double factorTime(const SparseMatrix<double>& A, Ordering method, Traversal traversal,
                         bool& ok) {
    DirectSolver<double> solver(method, Factorization::Cholesky, traversal);
    ok = solver.analyze(A);
    if (!ok) return 0;
    return bestOfThree([&] { solver.factor(A); });
}

static Row measure(const SparseMatrix<double>& A, const std::string& name, Ordering method) {
    Row row;
    row.name = name;
    row.nnzL = fill(A, method);

    const OrderEngine oe(method);
    row.order = bestOfThree([&] { Permutation P; oe.compute(A, P); });

    row.analyze = bestOfThree([&] {
        DirectSolver<double> s(method, Factorization::Cholesky, Traversal::LeftLooking);
        s.analyze(A);
    });
    row.analyzeMf = bestOfThree([&] {
        DirectSolver<double> s(method, Factorization::Cholesky, Traversal::Multifrontal);
        s.analyze(A);
    });

    bool ok = false;
    row.factLl = factorTime(A, method, Traversal::LeftLooking, ok);
    row.factRl = factorTime(A, method, Traversal::RightLooking, ok);
    row.factMf = factorTime(A, method, Traversal::Multifrontal, ok);

    DirectSolver<double> solver(method, Factorization::Cholesky, Traversal::LeftLooking);
    if (solver.analyze(A) && solver.factor(A)) {
        Vector<double> b(A.size()), x(A.size());
        for (std::size_t i = 0; i < A.size(); ++i) b[i] = 1.0;
        row.solve = bestOfThree([&] { solver.solve(b, x); });
    }
    return row;
}

int main(int argc, char** argv) {
    std::vector<int> sides = {32, 64, 100, 140};
    if (argc > 1) {
        sides.clear();
        for (int k = 1; k < argc; ++k) sides.push_back(std::atoi(argv[k]));
    }

    const std::vector<std::pair<std::string, Ordering>> methods = {
        {"Natural", Ordering::Natural}, {"MMD", Ordering::MMD},
        {"MMD1", Ordering::MMD1},       {"MMD2", Ordering::MMD2},
        {"MMD3", Ordering::MMD3},
        {"AMD", Ordering::AMD},         {"AMD1", Ordering::AMD1},
        {"AMD1B", Ordering::AMD1B},     {"AMD2", Ordering::AMD2},
        {"AMD2B", Ordering::AMD2B},
    };

    for (int side : sides) {
        const SparseMatrix<double> A = grid2D(side);
        std::printf("\n=== grid %dx%d, n = %d, Cholesky, real ===\n\n", side, side, side * side);
        std::printf("%-8s %8s %9s %10s %10s %9s %9s %9s %8s\n", "ordering", "order", "analyze",
                    "analyzeMF", "nnz(L)", "factLL", "factRL", "factMF", "solve");

        std::vector<Row> rows;
        for (const auto& m : methods) {
            const Row r = measure(A, m.first, m.second);
            rows.push_back(r);
            std::printf("%-8s %8.2f %9.2f %10.2f %10zu %9.2f %9.2f %9.2f %8.2f\n", r.name.c_str(),
                        r.order, r.analyze, r.analyzeMf, r.nnzL, r.factLl, r.factRl, r.factMf,
                        r.solve);
        }

        // Break-even against the vendored AMD, which is the default and the one ours would have to
        // displace. An ordering that analyzes faster and factors slower pays only while the factor
        // count stays below the ratio, so this is the number a caller actually faces: analyze runs
        // once per pattern, factor runs per Newton step or per time step.
        const Row* ref = nullptr;
        for (const Row& r : rows)
            if (r.name == "AMD") ref = &r;
        if (ref == nullptr) continue;

        std::printf("\n%-8s %12s %12s %14s   %s\n", "against", "analyze", "factor LL",
                    "break-even", "reading");
        std::printf("%-8s %12s %12s %14s\n", "AMD", "saved ms", "cost ms", "factorizations");
        for (const Row& r : rows) {
            if (&r == ref) continue;
            const double saved = ref->analyze - r.analyze;   // positive: this ordering analyzes faster
            const double cost  = r.factLl - ref->factLl;     // positive: and factors slower
            std::printf("%-8s %12.2f %12.2f", r.name.c_str(), saved, cost);
            if (saved > 0 && cost > 0)
                std::printf(" %14.1f   pays below that many factorizations\n", saved / cost);
            else if (saved < 0 && cost < 0)
                std::printf(" %14.1f   pays above that many factorizations\n", saved / cost);
            else if (saved >= 0 && cost <= 0)
                std::printf(" %14s   wins on both, at every count\n", "always");
            else
                std::printf(" %14s   loses on both, at every count\n", "never");
        }
    }
    return 0;
}
