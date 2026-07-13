#include "oblio/SymFactorEngine.h"

#include <algorithm>

namespace Oblio {

// Adapter: the factor's pattern needs only A's pattern, so the matrix overload pulls it out
// and forwards. The implementation below is free of Val and compiled once.
template<class Val>
bool SymFactorEngine::compute(const SparseMatrix<Val>& A, const Permutation& p,
                            const ElmForest& f, SymFactor& s) const {
    return compute(A.colPtr(), A.rowIdx(), p, f, s);
}

bool SymFactorEngine::compute(const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            const Permutation& p, const ElmForest& f, SymFactor& s) const {
    if (colPtr.empty())
        return false;
    const std::size_t size = colPtr.size() - 1;
    if (p.size() != size || f.size() != size)
        return false;

    const std::size_t supSize = f.supSize();

    // Copy the forest attributes the factor needs, so it stands alone afterwards.
    s.mSize      = size;
    s.mSupSize   = supSize;
    s.mNumTrees  = f.numTrees();
    s.mHeight    = f.height();
    s.mFirstRoot = f.firstRoot();
    s.mLastRoot  = f.lastRoot();

    s.mIdxToSupIdx = f.idxToSupIdx();
    s.mParent      = f.parent();
    s.mFirstChild  = f.firstChild();
    s.mNextSibling = f.nextSibling();
    s.mFrontSize   = f.frontSize();
    s.mUpdateSize  = f.updateSize();

    // Offsets into the flat index array: each supernode holds its front indices
    // followed by its update indices, and both counts are known from the forest.
    s.mSupPtr.resize(supSize + 1);
    s.mSupPtr[0] = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(supSize); ++kk)
        s.mSupPtr[kk + 1] = s.mSupPtr[kk] + s.mFrontSize[kk] + s.mUpdateSize[kk];
    s.mNumRowIdx = s.mSupPtr[supSize];
    s.mRowIdx.resize(s.mNumRowIdx);

    std::vector<std::size_t>  frontSupPtr;
    std::vector<std::int32_t> frontRowIdx;
    gatherFrontalIndices(s, frontSupPtr, frontRowIdx);

    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // Naming, following 0.9. Supernodes are doubled letters: jj a child, kk its parent, and
    // jj < kk since a child's label is below its parent's. A single letter is a column, so
    // the doubling is the column-to-supernode map applied, and that reads literally here:
    // lk is a front column of kk, which is to say idxToSupIdx[lk] == kk. An l prefix is the
    // factor (permuted) ordering and an a prefix the ordering of A, so ak is lk's column in
    // A, and li and ai are the corresponding rows.
    //
    // Positions into a flat array carry a prefix naming the array they walk:
    //
    //   ap   into A's colPtr/rowIdx          (0.9 calls all of these p)
    //   fp   into frontSupPtr/frontRowIdx    the front indices alone
    //   sp   into supPtr/rowIdx              the whole index set, front and update
    //
    // sp1 reads the child's index set, sp2 writes the parent's, keeping 0.9's 1-is-child,
    // 2-is-parent pairing. The scratch pair mirrors SymFactor's own, with front marking the
    // part it holds: frontSupPtr : frontRowIdx as supPtr : rowIdx.
    //
    // Rows and columns share one index space, the structure being symmetric, so an index is
    // just an index and the l or a prefix says only which ordering it is in. Which is why lk
    // subscripts a column below (colPtr[ak]) and is compared against a row on the next line
    // (li < lk) with no conversion: they are the same numbers.
    //
    // Marker array, so each union is done in one pass: as an index is added to the
    // index set of supernode kk it is marked with kk, which makes "is it already
    // there?" a single lookup.
    std::vector<std::int32_t> mark(size, NIL);

    // For every supernode kk, in increasing order (a topological order, so the
    // children of kk are complete when kk is reached).
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(supSize); ++kk) {
        std::size_t sp2 = s.mSupPtr[kk];   // write cursor into the index set of kk

        // The contribution of kk itself. Without threshold-based compression it
        // would be enough to read the sparsity pattern of the first front column
        // of kk, but compression groups columns whose patterns are only nearly
        // identical, so every front column is read and the patterns unioned.
        for (std::size_t fp = frontSupPtr[kk]; fp < frontSupPtr[kk + 1]; ++fp) {
            const std::int32_t lk = frontRowIdx[fp];
            const std::int32_t ak = newToOld[lk];

            // For every factor row li with a structural nonzero in factor column lk.
            for (std::size_t ap = colPtr[ak]; ap < colPtr[ak + 1]; ++ap) {
                const std::int32_t ai = rowIdx[ap];
                const std::int32_t li = oldToNew[ai];

                if (li < lk)
                    continue;   // above the front column, so not in the index set of kk
                if (mark[li] == kk)
                    continue;   // already in the index set of kk

                mark[li] = kk;
                s.mRowIdx[sp2++] = li;
            }
        }

        // The contribution of the children of kk.
        for (std::int32_t jj = s.mFirstChild[kk]; jj != NIL; jj = s.mNextSibling[jj]) {

            for (std::size_t sp1 = s.mSupPtr[jj]; sp1 < s.mSupPtr[jj + 1]; ++sp1) {
                const std::int32_t li = s.mRowIdx[sp1];

                if (s.mIdxToSupIdx[li] == jj)
                    continue;   // a front index of jj, so it dies with jj: only the
                                // update indices of jj carry into the index set of kk
                if (mark[li] == kk)
                    continue;   // already in the index set of kk

                mark[li] = kk;
                s.mRowIdx[sp2++] = li;
            }
        }
    }

    sortIndices(s);

    return true;
}

void SymFactorEngine::gatherFrontalIndices(const SymFactor& s,
        std::vector<std::size_t>&  frontSupPtr,
        std::vector<std::int32_t>& frontRowIdx) const {
    // The offsets are the front sizes, accumulated.
    frontSupPtr.resize(s.mSupSize + 1);
    frontSupPtr[0] = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(s.mSupSize); ++kk)
        frontSupPtr[kk + 1] = frontSupPtr[kk] + s.mFrontSize[kk];

    // Scatter each column into its supernode's slot. Columns are walked in increasing
    // order, so each supernode's front indices come out sorted.
    frontRowIdx.resize(s.mSize);
    std::vector<std::size_t> cursor(frontSupPtr.begin(), frontSupPtr.end() - 1);
    for (std::int32_t lk = 0; lk < static_cast<std::int32_t>(s.mSize); ++lk) {
        const std::int32_t kk = s.mIdxToSupIdx[lk];   // the supernode of column lk
        frontRowIdx[cursor[kk]++] = lk;
    }
}

void SymFactorEngine::sortIndices(SymFactor& s) const {
    const std::size_t numIdx = s.mRowIdx.size();

    // Count the supernodes each index appears in, then accumulate, giving the offsets of
    // the transpose.
    std::vector<std::size_t> supPtr(s.mSize + 1, 0);
    for (std::size_t sp = 0; sp < numIdx; ++sp)
        ++supPtr[s.mRowIdx[sp] + 1];
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(s.mSize); ++i)
        supPtr[i + 1] += supPtr[i];

    // First pass: build the transpose, the supernodes each index appears in.
    std::vector<std::int32_t> sup(numIdx);
    std::vector<std::size_t> cursor(supPtr.begin(), supPtr.end() - 1);
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(s.mSupSize); ++kk)
        for (std::size_t sp = s.mSupPtr[kk]; sp < s.mSupPtr[kk + 1]; ++sp) {
            const std::int32_t i = s.mRowIdx[sp];
            sup[cursor[i]++] = kk;
        }

    // Second pass: read the transpose back, walking the indices in increasing order, which
    // writes each supernode's index set sorted.
    std::vector<std::size_t> writeCursor(s.mSupPtr.begin(), s.mSupPtr.end() - 1);
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(s.mSize); ++i)
        for (std::size_t sp = supPtr[i]; sp < supPtr[i + 1]; ++sp) {
            const std::int32_t kk = sup[sp];
            s.mRowIdx[writeCursor[kk]++] = i;
        }
}

template bool SymFactorEngine::compute(const SparseMatrix<double>&, const Permutation&, const ElmForest&, SymFactor&) const;
template bool SymFactorEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, const ElmForest&, SymFactor&) const;

} // namespace Oblio
