// examples/example_matrix.cpp
// Getting a matrix into the solver, which is the first thing a caller has to do and the first
// place it can go wrong.
//
// SparseMatrix holds compressed sparse column: three arrays, one entry of colPtr per column plus a
// closing entry, and rowIdx and val running in parallel over the nonzeros. Column j occupies
// positions colPtr[j] through colPtr[j+1] - 1 in the other two. That is the standard format and
// most sparse codes speak it, so the usual job is a conversion rather than a construction.
//
// Below, the same 5x5 matrix built three ways, which must produce the same solution:
//
//   1. CSC written out directly, which is what the other examples do
//   2. from a dense array, the easy case, and the one to reach for while prototyping
//   3. from coordinate (triplet) form, which is what an assembly loop naturally produces and what
//      a file on disk usually holds
//
// **Three preconditions the solver requires and does not check**, which is what the conversion
// below has to satisfy and the reason it is worth reading. They are stated in SparseMatrix.h and
// tracked in docs/TODO.md under "Validate the input matrix":
//
//   FULL STORAGE. Both triangles, not one. A symmetric matrix stored as its lower triangle alone
//   is a different matrix as far as the structural phases are concerned, and they read the pattern
//   directly rather than symmetrizing it.
//
//   THE DIAGONAL IS STRUCTURALLY PRESENT, even where it is numerically zero. Symbolic
//   factorization builds a column's index set from that column's stored pattern, so a column
//   missing its own diagonal never enters its own index set. A zero diagonal is completely
//   ordinary in an indefinite problem, which is exactly what dynamic LDL exists for, so this is
//   not a corner case: a KKT matrix whose zero block is stored sparsely hits it immediately.
//
//   ROW INDICES SORTED ASCENDING within each column, and no duplicates. Assembly assigns rather
//   than accumulates, so a repeated entry would keep whichever value was written last.
//
// The triplet conversion below meets all three: it orders, it merges duplicates, and it puts the
// diagonal in whether or not anything landed there. Full storage is the caller's to respect, since
// only the caller knows which half it holds.
//
// Compile (macOS, from repo root):
//   g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude examples/example_matrix.cpp src/*.cpp -framework Accelerate -o example_matrix_cpp
// Linux: replace `-framework Accelerate` with `-llapack -lblas`.

#include "oblio/DirectSolver.h"

#include <algorithm>
#include <cstdio>
#include <map>
#include <vector>

using namespace Oblio;

namespace {

using Val = double;

// One nonzero, as an assembly loop or a file would give it: row, column, value.
struct Triplet {
    std::int32_t row;
    std::int32_t col;
    Val          val;
};

// Coordinate form to CSC. This is the conversion worth stealing, and the three preconditions above
// are exactly what it takes care of.
//
// A map keyed by (column, row) does all three at once: it orders by column and then by row, which
// is the sort, and one key per entry, which is the deduplication. Accumulating on collision is the
// choice an assembly needs, since a finite element method emits the same (row, column) from every
// element touching it and means them summed. A caller who instead wants "last one wins" assigns
// here rather than adding.
//
// A map is the readable version rather than the fast one. A production conversion counting-sorts
// into the arrays directly; at that point this loop is the specification of what it must produce.
SparseMatrix<Val> fromTriplets(std::size_t size, const std::vector<Triplet>& triplets) {
    std::map<std::pair<std::int32_t, std::int32_t>, Val> entries;

    for (const Triplet& t : triplets)
        entries[{t.col, t.row}] += t.val;

    // The diagonal, structurally present whether or not anything landed on it. Inserting a zero
    // where a key is absent is what keeps a numerically zero diagonal in the pattern; where one is
    // present this leaves it alone.
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(size); ++j)
        entries.emplace(std::make_pair(j, j), Val(0));

    // The map already walks in column-then-row order, so filling the arrays is one pass with a
    // cursor: take everything belonging to column j, then close j's slot in colPtr. A column with
    // no entries closes at the same position the previous one did, which is how CSC spells empty.
    std::vector<std::size_t>  colPtr(size + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<Val>          val;
    rowIdx.reserve(entries.size());
    val.reserve(entries.size());

    auto entry = entries.begin();
    for (std::size_t j = 0; j < size; ++j) {
        while (entry != entries.end() && entry->first.first == static_cast<std::int32_t>(j)) {
            rowIdx.push_back(entry->first.second);
            val.push_back(entry->second);
            ++entry;
        }
        colPtr[j + 1] = rowIdx.size();
    }

    return SparseMatrix<Val>(size, colPtr, rowIdx, val);
}

// A dense array to CSC, dropping the structural zeros but keeping the whole diagonal.
SparseMatrix<Val> fromDense(const std::vector<std::vector<Val>>& dense) {
    const std::size_t size = dense.size();
    std::vector<std::size_t>  colPtr(size + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<Val>          val;
    for (std::size_t j = 0; j < size; ++j) {
        for (std::size_t i = 0; i < size; ++i)
            if (dense[i][j] != Val(0) || i == j) {
                rowIdx.push_back(static_cast<std::int32_t>(i));
                val.push_back(dense[i][j]);
            }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<Val>(size, colPtr, rowIdx, val);
}

// Solve and report, so the three constructions can be compared on their answers.
bool solveAndPrint(const char* what, const SparseMatrix<Val>& A, Vector<Val>& x) {
    Vector<Val> b(A.size());
    for (std::size_t i = 0; i < A.size(); ++i) b[i] = Val(1);

    DirectSolver<Val> solver(Ordering::AmdCompacted, Factorization::DynamicLDLT);
    if (!solver.analyze(A)) { printf("  %-22s analyze refused\n", what); return false; }
    if (!solver.factor(A))  { printf("  %-22s factor refused\n",  what); return false; }
    if (!solver.solve(b, x)) { printf("  %-22s solve refused\n",  what); return false; }

    printf("  %-22s nnz %2zu   residual %.2e\n", what, A.nnz(), solver.relativeResidual(A, b, x));
    return true;
}

} // namespace

int main() {
    // The matrix, drawn once. Tridiagonal, and deliberately carrying a zero on the diagonal of
    // column 2, which is what makes the second demonstration below possible.
    //
    //   [  4  -1   0   0   0 ]
    //   [ -1   4  -1   0   0 ]
    //   [  0  -1   0  -1   0 ]      <- a structurally present, numerically zero diagonal
    //   [  0   0  -1   4  -1 ]
    //   [  0   0   0  -1   4 ]
    //
    const std::size_t size = 5;
    std::vector<std::vector<Val>> dense(size, std::vector<Val>(size, Val(0)));
    for (std::size_t i = 0; i < size; ++i) dense[i][i] = (i == 2) ? Val(0) : Val(4);
    for (std::size_t i = 0; i + 1 < size; ++i) dense[i][i + 1] = dense[i + 1][i] = Val(-1);

    printf("The same matrix, built three ways:\n\n");

    // 1. CSC by hand. Column j runs from colPtr[j] to colPtr[j+1] - 1 in the other two arrays, row
    // indices ascending within each column, and every diagonal present including column 2's zero.
    const std::vector<std::size_t>  colPtr = {0, 2, 5, 8, 11, 13};
    const std::vector<std::int32_t> rowIdx = {0, 1,  0, 1, 2,  1, 2, 3,  2, 3, 4,  3, 4};
    const std::vector<Val>          val    = {4, -1, -1, 4, -1, -1, 0, -1, -1, 4, -1, -1, 4};
    const SparseMatrix<Val> byHand(size, colPtr, rowIdx, val);

    // 2. From the dense array above.
    const SparseMatrix<Val> byDense = fromDense(dense);

    // 3. From triplets, given out of order and with one entry split in two, both of which the
    // conversion has to survive. The split pair sums to -1, which is the assembly case.
    const std::vector<Triplet> triplets = {
        {3, 4, Val(-1)}, {4, 4, Val(4)},  {1, 0, Val(-1)}, {0, 0, Val(4)},
        {2, 1, Val(-1)}, {1, 1, Val(4)},  {0, 1, Val(-1)}, {4, 3, Val(-1)},
        {3, 3, Val(4)},  {2, 3, Val(-1)}, {1, 2, Val(-0.5)}, {1, 2, Val(-0.5)},
        {3, 2, Val(-1)},
        // Column 2's diagonal is not listed at all, and the conversion inserts it as a zero.
    };
    const SparseMatrix<Val> byTriplets = fromTriplets(size, triplets);

    // Three separate calls rather than one chained expression, so that all three run and report
    // even if an earlier one fails. `&&` would short-circuit and hide the rest, and `&` on bools
    // avoids that at the price of a warning that Apple Clang raises and g++ does not.
    Vector<Val> x1(size), x2(size), x3(size);
    const bool ok1 = solveAndPrint("CSC by hand", byHand, x1);
    const bool ok2 = solveAndPrint("from a dense array", byDense, x2);
    const bool ok3 = solveAndPrint("from triplets", byTriplets, x3);

    if (ok1 && ok2 && ok3) {
        double worst = 0;
        for (std::size_t i = 0; i < size; ++i) {
            worst = std::max(worst, std::abs(x1[i] - x2[i]));
            worst = std::max(worst, std::abs(x1[i] - x3[i]));
        }
        printf("\n  Same nnz, same solution: the three agree to %.1e.\n", worst);
        printf("  x = [");
        for (std::size_t i = 0; i < size; ++i) printf(" %.4f", x1[i]);
        printf(" ]\n");
    }

    return 0;
}
