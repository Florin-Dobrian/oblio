#include "oblio/SymFactEngine.h"

#include <algorithm>

namespace Oblio {

// Adapter: the factor's pattern needs only A's pattern, so the matrix overload pulls it out
// and forwards. The implementation below is free of Val and compiled once.
template<class Val>
bool SymFactEngine::compute(const SparseMatrix<Val>& A, const Permutation& p,
                            const ElmForest& f, SymFact& s) const {
    return compute(A.colPtr(), A.rowIdx(), p, f, s);
}

bool SymFactEngine::compute(const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            const Permutation& p, const ElmForest& f, SymFact& s) const {
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
    s.mIdxPtr.resize(supSize + 1);
    s.mIdxPtr[0] = 0;
    for (std::size_t s1 = 0; s1 < supSize; ++s1)
        s.mIdxPtr[s1 + 1] = s.mIdxPtr[s1] + s.mFrontSize[s1] + s.mUpdateSize[s1];
    s.mNumIdx = s.mIdxPtr[supSize];
    s.mIdx.resize(s.mNumIdx);

    std::vector<std::size_t>  frontIdxPtr;
    std::vector<std::int32_t> frontIdx;
    gatherFrontIdx(s, frontIdxPtr, frontIdx);

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
            for (std::size_t ai = colPtr[ac]; ai < colPtr[ac + 1]; ++ai) {
                const std::size_t ar = static_cast<std::size_t>(rowIdx[ai]);
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

    sortIdx(s);

    return true;
}

void SymFactEngine::gatherFrontIdx(const SymFact& s,
        std::vector<std::size_t>&  frontIdxPtr,
        std::vector<std::int32_t>& frontIdx) const {
    // The offsets are the front sizes, accumulated.
    frontIdxPtr.resize(s.mSupSize + 1);
    frontIdxPtr[0] = 0;
    for (std::size_t k = 0; k < s.mSupSize; ++k)
        frontIdxPtr[k + 1] = frontIdxPtr[k] + s.mFrontSize[k];

    // Scatter each column into its supernode's slot. Columns are walked in increasing
    // order, so each supernode's front indices come out sorted.
    frontIdx.resize(s.mSize);
    std::vector<std::size_t> cursor(frontIdxPtr.begin(), frontIdxPtr.end() - 1);
    for (std::size_t lc = 0; lc < s.mSize; ++lc) {
        const std::size_t k = static_cast<std::size_t>(s.mIdxToSupIdx[lc]);
        frontIdx[cursor[k]++] = static_cast<std::int32_t>(lc);
    }
}

void SymFactEngine::sortIdx(SymFact& s) const {
    const std::size_t numIdx = s.mIdx.size();

    // Count the supernodes each index appears in, then accumulate, giving the offsets of
    // the transpose.
    std::vector<std::size_t> supPtr(s.mSize + 1, 0);
    for (std::size_t t = 0; t < numIdx; ++t)
        ++supPtr[static_cast<std::size_t>(s.mIdx[t]) + 1];
    for (std::size_t i = 0; i < s.mSize; ++i)
        supPtr[i + 1] += supPtr[i];

    // First pass: build the transpose, the supernodes each index appears in.
    std::vector<std::int32_t> sup(numIdx);
    std::vector<std::size_t> cursor(supPtr.begin(), supPtr.end() - 1);
    for (std::size_t k = 0; k < s.mSupSize; ++k)
        for (std::size_t t = s.mIdxPtr[k]; t < s.mIdxPtr[k + 1]; ++t) {
            const std::size_t i = static_cast<std::size_t>(s.mIdx[t]);
            sup[cursor[i]++] = static_cast<std::int32_t>(k);
        }

    // Second pass: read the transpose back, walking the indices in increasing order, which
    // writes each supernode's index set sorted.
    std::vector<std::size_t> writeCursor(s.mIdxPtr.begin(), s.mIdxPtr.end() - 1);
    for (std::size_t i = 0; i < s.mSize; ++i)
        for (std::size_t t = supPtr[i]; t < supPtr[i + 1]; ++t) {
            const std::size_t k = static_cast<std::size_t>(sup[t]);
            s.mIdx[writeCursor[k]++] = static_cast<std::int32_t>(i);
        }
}

template bool SymFactEngine::compute(const SparseMatrix<double>&, const Permutation&, const ElmForest&, SymFact&) const;
template bool SymFactEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, const ElmForest&, SymFact&) const;

} // namespace Oblio
