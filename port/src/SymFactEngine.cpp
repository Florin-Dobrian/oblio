#include "oblio/SymFactEngine.h"

#include <algorithm>

namespace Oblio {

template<class Val>
bool SymFactEngine::compute(const SparseMatrix<Val>& A, const Permutation& p,
                            const ElmForest& f, SymFact& s) const {
    const std::size_t size = A.size();
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
    s.mIdxPtr.resize(supSize + 1);
    s.mIdxPtr[0] = 0;
    for (std::size_t s1 = 0; s1 < supSize; ++s1)
        s.mIdxPtr[s1 + 1] = s.mIdxPtr[s1] + s.mFrontSize[s1] + s.mUpdateSize[s1];
    s.mNumIdx = s.mIdxPtr[supSize];
    s.mIdx.resize(s.mNumIdx);

    std::vector<std::size_t>  frontIdxPtr;
    std::vector<std::int32_t> frontIdx;
    gatherFrontIdx(size, supSize, s.mIdxToSupIdx, s.mFrontSize, frontIdxPtr, frontIdx);

    const std::vector<std::size_t>&  aColPtr = A.colPtr();
    const std::vector<std::int32_t>& aRowIdx = A.rowIdx();
    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // Marker array, so each union is done in one pass: as an index is added to the
    // index set of supernode s2 it is marked with s2, which makes "is it already
    // there?" a single lookup.
    std::vector<std::int32_t> mark(size, NIL);

    // For every supernode s2, in increasing order (a topological order, so the
    // children of s2 are complete when s2 is reached).
    for (std::size_t s2 = 0; s2 < supSize; ++s2) {
        std::size_t lk2 = s.mIdxPtr[s2];   // write cursor into the index set of s2

        // The contribution of s2 itself. Without threshold-based compression it
        // would be enough to read the sparsity pattern of the first front column
        // of s2, but compression groups columns whose patterns are only nearly
        // identical, so every front column is read and the patterns unioned.
        for (std::size_t lj = frontIdxPtr[s2]; lj < frontIdxPtr[s2 + 1]; ++lj) {
            const std::size_t lc = static_cast<std::size_t>(frontIdx[lj]);
            const std::size_t ac = static_cast<std::size_t>(newToOld[lc]);

            // For every factor row lr with a structural nonzero in factor column lc.
            for (std::size_t ai = aColPtr[ac]; ai < aColPtr[ac + 1]; ++ai) {
                const std::size_t ar = static_cast<std::size_t>(aRowIdx[ai]);
                const std::size_t lr = static_cast<std::size_t>(oldToNew[ar]);

                if (lr < lc)
                    continue;   // above the front column, so not in the index set of s2
                if (mark[lr] == static_cast<std::int32_t>(s2))
                    continue;   // already in the index set of s2

                mark[lr] = static_cast<std::int32_t>(s2);
                s.mIdx[lk2++] = static_cast<std::int32_t>(lr);
            }
        }

        // The contribution of the children of s2.
        for (std::int32_t s1 = s.mFirstChild[s2]; s1 != NIL; s1 = s.mNextSibling[s1]) {
            const std::size_t us1 = static_cast<std::size_t>(s1);

            for (std::size_t lk1 = s.mIdxPtr[us1]; lk1 < s.mIdxPtr[us1 + 1]; ++lk1) {
                const std::size_t lr = static_cast<std::size_t>(s.mIdx[lk1]);

                if (s.mIdxToSupIdx[lr] == s1)
                    continue;   // a front index of s1, so it dies with s1: only the
                                // update indices of s1 carry into the index set of s2
                if (mark[lr] == static_cast<std::int32_t>(s2))
                    continue;   // already in the index set of s2

                mark[lr] = static_cast<std::int32_t>(s2);
                s.mIdx[lk2++] = static_cast<std::int32_t>(lr);
            }
        }
    }

    sortIdx(size, supSize, s.mIdxPtr, s.mIdx);

    return true;
}

void SymFactEngine::gatherFrontIdx(std::size_t size, std::size_t supSize,
        const std::vector<std::int32_t>& idxToSupIdx,
        const std::vector<std::size_t>&  frontSize,
        std::vector<std::size_t>&  frontIdxPtr,
        std::vector<std::int32_t>& frontIdx) const {
    // The offsets are the front sizes, accumulated.
    frontIdxPtr.resize(supSize + 1);
    frontIdxPtr[0] = 0;
    for (std::size_t s = 0; s < supSize; ++s)
        frontIdxPtr[s + 1] = frontIdxPtr[s] + frontSize[s];

    // Scatter each column into its supernode's slot. Columns are walked in
    // increasing order, so each supernode's front indices come out sorted.
    frontIdx.resize(size);
    std::vector<std::size_t> cursor(frontIdxPtr.begin(), frontIdxPtr.end() - 1);
    for (std::size_t lc = 0; lc < size; ++lc) {
        const std::size_t s = static_cast<std::size_t>(idxToSupIdx[lc]);
        frontIdx[cursor[s]++] = static_cast<std::int32_t>(lc);
    }
}

void SymFactEngine::sortIdx(std::size_t size, std::size_t supSize,
        const std::vector<std::size_t>& idxPtr,
        std::vector<std::int32_t>& idx) const {
    const std::size_t numIdx = idx.size();

    // Count the supernodes each index appears in, then accumulate, giving the
    // offsets of the transpose.
    std::vector<std::size_t> supPtr(size + 1, 0);
    for (std::size_t k = 0; k < numIdx; ++k)
        ++supPtr[static_cast<std::size_t>(idx[k]) + 1];
    for (std::size_t i = 0; i < size; ++i)
        supPtr[i + 1] += supPtr[i];

    // First pass: build the transpose, the supernodes each index appears in.
    std::vector<std::int32_t> sup(numIdx);
    std::vector<std::size_t> cursor(supPtr.begin(), supPtr.end() - 1);
    for (std::size_t s = 0; s < supSize; ++s)
        for (std::size_t k = idxPtr[s]; k < idxPtr[s + 1]; ++k) {
            const std::size_t i = static_cast<std::size_t>(idx[k]);
            sup[cursor[i]++] = static_cast<std::int32_t>(s);
        }

    // Second pass: read the transpose back, walking the indices in increasing
    // order, which writes each supernode's index set sorted.
    std::vector<std::size_t> writeCursor(idxPtr.begin(), idxPtr.end() - 1);
    for (std::size_t i = 0; i < size; ++i)
        for (std::size_t k = supPtr[i]; k < supPtr[i + 1]; ++k) {
            const std::size_t s = static_cast<std::size_t>(sup[k]);
            idx[writeCursor[s]++] = static_cast<std::int32_t>(i);
        }
}

template bool SymFactEngine::compute(const SparseMatrix<double>&, const Permutation&, const ElmForest&, SymFact&) const;
template bool SymFactEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, const ElmForest&, SymFact&) const;

} // namespace Oblio
