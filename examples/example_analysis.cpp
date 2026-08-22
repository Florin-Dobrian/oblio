// examples/example_analysis.cpp
// What the ordering buys, measured without factoring anything.
//
// The pipeline splits into analyze, factor and solve, and the first of those is where the fill is
// decided. This example runs analyze alone, over every ordering method, and reports what the
// symbolic factorization predicts: how large L will be, how many supernodes it falls into, and how
// tall the elimination forest is. Nothing numeric happens, and nothing needs to.
//
// **The analysis reads the pattern of A and never a value**, which is why the matrices below are
// filled with arbitrary numbers and why the same numbers serve every row. Ordering, the elimination
// forest and the symbolic factorization are graph algorithms; the values are the factor's business.
// That is also why one analysis serves a whole sequence of matrices sharing a pattern, which is the
// case a Newton iteration or a time-stepping loop presents.
//
// The columns:
//
//   nnz(L)      entries in the factor, from SymFactor::nnz(): each supernode's own lower
//               triangle plus its update rectangle
//   fill        nnz(L) minus the entries A already had below and on the diagonal: what elimination
//               adds, and the quantity the ordering exists to reduce
//   supernodes  how many blocks the columns group into. Fewer and larger is better for the numeric
//               phase, which hands each one to a dense BLAS kernel
//   height      the elimination forest's depth. A tall thin forest is a chain of dependencies; a
//               short bushy one has independent subtrees, which is where parallelism would come
//               from
//
// Read the arrow first: it is the whole argument in nine vertices.
//
// Compile (macOS, from repo root):
//   g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude examples/example_analysis.cpp src/*.cpp -framework Accelerate -o example_analysis_cpp
// Linux: replace `-framework Accelerate` with `-llapack -lblas`.

#include "oblio/DirectSolver.h"

#include <cstdio>
#include <set>
#include <vector>

using namespace Oblio;

namespace {

using Val = double;

const char* name(Ordering ordering) {
    switch (ordering) {
        case Ordering::Natural: return "Natural";
        case Ordering::MMD:     return "MMD";
        case Ordering::MMD3:    return "MMD3";
        case Ordering::AMD:     return "AMD";
        case Ordering::AMD3:    return "AMD3";
    }
    return "?";
}

// Build a full-symmetric CSC matrix from adjacency sets that already include the diagonal. The
// values are arbitrary and never read: 4 on the diagonal, -1 off it, so the matrix would also be
// factorable if anyone asked, which nothing here does.
SparseMatrix<Val> fromAdjacency(const std::vector<std::set<std::int32_t>>& adjacency) {
    const std::size_t size = adjacency.size();
    std::vector<std::size_t>  colPtr(size + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<Val>          val;
    for (std::size_t j = 0; j < size; ++j) {
        for (std::int32_t i : adjacency[j]) {
            rowIdx.push_back(i);
            val.push_back(i == static_cast<std::int32_t>(j) ? Val(4) : Val(-1));
        }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<Val>(size, colPtr, rowIdx, val);
}

// A star: one hub adjacent to every other vertex, the rest adjacent only to the hub. Whether the
// hub is numbered first or last is the entire content of this example's first two blocks.
SparseMatrix<Val> arrow(std::size_t size, bool hubFirst) {
    const std::int32_t hub = hubFirst ? 0 : static_cast<std::int32_t>(size) - 1;
    std::vector<std::set<std::int32_t>> adjacency(size);
    for (std::int32_t v = 0; v < static_cast<std::int32_t>(size); ++v) {
        adjacency[v].insert(v);
        if (v != hub) { adjacency[v].insert(hub); adjacency[hub].insert(v); }
    }
    return fromAdjacency(adjacency);
}

// A five-point grid Laplacian on a square mesh, numbered row by row.
SparseMatrix<Val> grid(std::size_t side) {
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
    return fromAdjacency(adjacency);
}

// One matrix, every ordering, analysis only.
void report(const char* what, const SparseMatrix<Val>& A) {
    // A is stored full, both triangles with the diagonal, so its lower triangle holds this many.
    const std::size_t nnzTrilA = (A.nnz() + A.size()) / 2;

    printf("\n%s: n = %zu, nnz(A) = %zu, of which %zu on and below the diagonal\n",
           what, A.size(), A.nnz(), nnzTrilA);
    printf("  %-8s  %8s  %8s  %10s  %6s\n",
           "ordering", "nnz(L)", "fill", "supernodes", "height");
    printf("  %-8s  %8s  %8s  %10s  %6s\n",
           "--------", "------", "----", "----------", "------");

    for (Ordering ordering : {Ordering::Natural, Ordering::MMD, Ordering::MMD3,
                              Ordering::AMD, Ordering::AMD3}) {
        // The factorization and traversal are left at their defaults and never used: analyze reads
        // neither. Only the traversal would matter, and only across the multifrontal boundary,
        // where the forest itself differs.
        DirectSolver<Val> solver(ordering);
        if (!solver.analyze(A)) { printf("  %-8s  analyze failed\n", name(ordering)); continue; }

        // The three intermediates the facade owns, which it exposes for exactly this.
        const ElmForest& ef = solver.elmForest();
        const SymFactor& sf = solver.symFactor();

        printf("  %-8s  %8zu  %8zu  %10zu  %6zu\n",
               name(ordering), sf.nnz(), sf.nnz() - nnzTrilA, sf.snodeSize(),
               ef.height());
    }
}

} // namespace

int main() {
    // The arrow, twice. Same graph, same nine vertices, one hub: only the numbering differs, and
    // that is the point. Eliminating the hub first makes every other vertex mutually adjacent, so
    // the factor is completely dense and lands in a single supernode. Eliminating it last costs
    // nothing at all, because each leaf's only neighbor is already gone by the time the hub is
    // reached. Every fill-reducing method finds the second arrangement from the first.
    report("arrow, hub numbered first", arrow(9, true));
    report("arrow, hub numbered last",  arrow(9, false));

    // A 3x3 grid, the smallest mesh worth ordering, and small enough to check by hand. Natural
    // costs about sixty percent more fill than any fill-reducing method, and all eight of those
    // tie exactly: at this size there is nothing left for them to disagree about.
    report("grid 3x3", grid(3));

    // An 8x8 grid, where they finally separate. The spread among the eight is a few percent,
    // against Natural's factor of two, which is the honest shape of the result: reordering matters
    // far more than choosing which reordering. Note the height column, where Natural's forest is a
    // chain of 56 and the others are around a dozen.
    report("grid 8x8", grid(8));

    return 0;
}
