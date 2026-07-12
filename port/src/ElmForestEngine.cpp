#include "oblio/ElmForestEngine.h"

#include <deque>

namespace Oblio {

template<class Val>
bool ElmForestEngine::compute(const SparseMatrix<Val>& A, const Permutation& p,
                              ElmForest& f) const {
    const std::size_t size = A.size();
    if (p.size() != size)
        return false;

    // A is full-symmetric: each column already holds its complete neighbour list,
    // so the etree reads it directly (diagonal entries map to lc1 == lc2 and are
    // skipped by the lc1 < lc2 test).
    f.mSize = size;
    computeParent(size, A.colPtr(), A.rowIdx(), p.oldToNew(), p.newToOld(), f.mParent);

    // Start nodal: one trivial supernode per column, so the supernode count equals the
    // column count and the map is the identity. compressFundamental below merges these
    // into fundamental supernodes; everything from here on is written to be indifferent
    // to which regime it is in.
    f.mSupSize = size;
    f.mIdxToSupIdx.resize(size);
    for (std::size_t j = 0; j < size; ++j)
        f.mIdxToSupIdx[j] = static_cast<std::int32_t>(j);

    // The child, sibling and root links, then the front and update sizes. Both must
    // precede compression, whose merge test reads a supernode's only-child status from
    // the links and its sparsity pattern from the sizes.
    finalizeLinks(f.mSupSize, f.mParent, f.mFirstChild, f.mLastChild,
                  f.mNextSibling, f.mPreviousSibling, f.mNumTrees,
                  f.mFirstRoot, f.mLastRoot);
    computeFrontAndUpdateSizes(size, f.mSupSize, A.colPtr(), A.rowIdx(),
                               p.oldToNew(), p.newToOld(), f.mParent,
                               f.mIdxToSupIdx, f.mFrontSize, f.mUpdateSize);

    // Merge the trivial supernodes into fundamental supernodes, unless asked to stay
    // nodal. This rewrites the map, the parent links, the sizes and the child/sibling
    // links; skipping it leaves the forest one supernode per column.
    if (mSupernodes == Supernodes::Fundamental)
        compressFundamental(size, f.mSupSize, f.mIdxToSupIdx, f.mParent,
                            f.mFirstChild, f.mLastChild, f.mNextSibling, f.mPreviousSibling,
                            f.mFrontSize, f.mUpdateSize,
                            f.mNumTrees, f.mFirstRoot, f.mLastRoot);

    // Height last, so it is computed once, on the final forest.
    f.mHeight = computeHeight(f.mSupSize, f.mLastRoot, f.mParent,
                              f.mLastChild, f.mPreviousSibling);

    return true;
}

void ElmForestEngine::computeParent(std::size_t size,
        const std::vector<std::size_t>&  colPtr,
        const std::vector<std::int32_t>& rowIdx,
        const std::vector<std::int32_t>& oldToNew,
        const std::vector<std::int32_t>& newToOld,
        std::vector<std::int32_t>& parent) const {
    parent.assign(size, NIL);
    std::vector<std::int32_t> ancestor(size, NIL);

    // For each factor column lc2 (increasing), look at its neighbours mapping to
    // earlier columns lc1 < lc2; path-compress to attach lc1's subtree under lc2.
    for (std::size_t lc2 = 0; lc2 < size; ++lc2) {
        const std::size_t ac2 = static_cast<std::size_t>(newToOld[lc2]);
        for (std::size_t sp = colPtr[ac2]; sp < colPtr[ac2 + 1]; ++sp) {
            const std::size_t lc1 = static_cast<std::size_t>(oldToNew[static_cast<std::size_t>(rowIdx[sp])]);
            if (lc1 >= lc2)
                continue;   // later column or the diagonal itself
            std::size_t lc3 = lc1;
            while (ancestor[lc3] != NIL && static_cast<std::size_t>(ancestor[lc3]) != lc2) {
                const std::size_t lc4 = static_cast<std::size_t>(ancestor[lc3]);
                ancestor[lc3] = static_cast<std::int32_t>(lc2);
                lc3 = lc4;
            }
            if (ancestor[lc3] == NIL) {
                ancestor[lc3] = static_cast<std::int32_t>(lc2);
                parent[lc3] = static_cast<std::int32_t>(lc2);
            }
        }
    }
}

void ElmForestEngine::finalizeLinks(std::size_t supSize,
        const std::vector<std::int32_t>& parent,
        std::vector<std::int32_t>& firstChild,
        std::vector<std::int32_t>& lastChild,
        std::vector<std::int32_t>& nextSibling,
        std::vector<std::int32_t>& previousSibling,
        std::size_t&  numTrees,
        std::int32_t& firstRoot,
        std::int32_t& lastRoot) const {
    firstChild.assign(supSize, NIL);
    lastChild.assign(supSize, NIL);
    nextSibling.assign(supSize, NIL);
    previousSibling.assign(supSize, NIL);
    numTrees  = 0;
    firstRoot = NIL;
    lastRoot  = NIL;

    // For each supernode in decreasing order, front-insert it into its parent's
    // child list, or into the root list if it has no parent. Front-insertion of
    // decreasing labels leaves both the child lists and the root list in
    // increasing label order.
    for (std::size_t t = supSize; t > 0; --t) {
        const std::int32_t s1 = static_cast<std::int32_t>(t - 1);
        const std::int32_t s2 = parent[t - 1];
        if (s2 == NIL) {
            // s1 is a tree root.
            nextSibling[t - 1] = firstRoot;
            if (firstRoot == NIL)
                lastRoot = s1;
            else
                previousSibling[static_cast<std::size_t>(firstRoot)] = s1;
            firstRoot = s1;
            ++numTrees;
        } else {
            // s1 becomes a child of s2.
            const std::size_t us2 = static_cast<std::size_t>(s2);
            nextSibling[t - 1] = firstChild[us2];
            if (firstChild[us2] == NIL)
                lastChild[us2] = s1;
            else
                previousSibling[static_cast<std::size_t>(firstChild[us2])] = s1;
            firstChild[us2] = s1;
        }
    }
}

std::size_t ElmForestEngine::computeHeight(std::size_t supSize, std::int32_t lastRoot,
        const std::vector<std::int32_t>& parent,
        const std::vector<std::int32_t>& lastChild,
        const std::vector<std::int32_t>& previousSibling) const {
    // Breadth-first from the roots: a supernode is processed after its parent, so
    // its depth is the parent's depth plus one, and the height is the largest
    // depth plus one. The root chain and each child list are walked backward via
    // previousSibling; order is irrelevant to a maximum.
    std::vector<std::size_t> depth(supSize, 0);
    std::deque<std::int32_t> queue;
    for (std::int32_t s = lastRoot; s != NIL; s = previousSibling[static_cast<std::size_t>(s)])
        queue.push_back(s);
    while (!queue.empty()) {
        const std::int32_t s2 = queue.front();
        queue.pop_front();
        const std::size_t us2 = static_cast<std::size_t>(s2);
        if (parent[us2] != NIL)
            depth[us2] = depth[static_cast<std::size_t>(parent[us2])] + 1;
        for (std::int32_t s1 = lastChild[us2]; s1 != NIL;
             s1 = previousSibling[static_cast<std::size_t>(s1)])
            queue.push_back(s1);
    }
    std::size_t height = 0;
    for (std::size_t s = 0; s < supSize; ++s)
        if (height < depth[s])
            height = depth[s];
    if (supSize > 0)
        ++height;
    return height;
}

void ElmForestEngine::computeFrontAndUpdateSizes(std::size_t size, std::size_t supSize,
        const std::vector<std::size_t>&  colPtr,
        const std::vector<std::int32_t>& rowIdx,
        const std::vector<std::int32_t>& oldToNew,
        const std::vector<std::int32_t>& newToOld,
        const std::vector<std::int32_t>& parent,
        const std::vector<std::int32_t>& idxToSupIdx,
        std::vector<std::size_t>& frontSize,
        std::vector<std::size_t>& updateSize) const {
    // Front size: the number of columns each supernode owns, counted from the map.
    // Trivial supernodes give all 1s; the map-count stays correct once
    // fundamental-supernode compression merges columns, unlike a fixed fill.
    frontSize.assign(supSize, 0);
    for (std::size_t j = 0; j < size; ++j)
        ++frontSize[static_cast<std::size_t>(idxToSupIdx[j])];

    // Update size: the subdiagonal nonzero count of each L column, by the pruned-
    // row-subtree walk. For factor column lc2, climb from each earlier neighbour
    // lc1 up the etree toward lc2, marking to avoid double counting; every column
    // passed below the root gains one update index. This is 0.9's columnSize with
    // the diagonal left out, so it yields updateSize directly (== columnSize - 1).
    updateSize.assign(supSize, 0);
    std::vector<std::int32_t> mark(size, NIL);
    for (std::size_t lc2 = 0; lc2 < size; ++lc2) {
        mark[lc2] = static_cast<std::int32_t>(lc2);
        const std::size_t ac2 = static_cast<std::size_t>(newToOld[lc2]);
        for (std::size_t sp = colPtr[ac2]; sp < colPtr[ac2 + 1]; ++sp) {
            const std::size_t lc1 =
                static_cast<std::size_t>(oldToNew[static_cast<std::size_t>(rowIdx[sp])]);
            if (lc1 >= lc2)
                continue;   // later column or the diagonal itself
            std::size_t lc3 = lc1;
            while (mark[lc3] != static_cast<std::int32_t>(lc2)) {
                ++updateSize[lc3];
                mark[lc3] = static_cast<std::int32_t>(lc2);
                lc3 = static_cast<std::size_t>(parent[lc3]);
            }
        }
    }
}

void ElmForestEngine::compressFundamental(std::size_t size, std::size_t& supSize,
        std::vector<std::int32_t>& idxToSupIdx,
        std::vector<std::int32_t>& parent,
        std::vector<std::int32_t>& firstChild,
        std::vector<std::int32_t>& lastChild,
        std::vector<std::int32_t>& nextSibling,
        std::vector<std::int32_t>& previousSibling,
        std::vector<std::size_t>&  frontSize,
        std::vector<std::size_t>&  updateSize,
        std::size_t&  numTrees,
        std::int32_t& firstRoot,
        std::int32_t& lastRoot) const {
    // Merge each supernode into its child where they belong to one fundamental
    // supernode, that is, where they form a path in the forest whose factor columns
    // share one sparsity pattern. Two conditions, for child s1 of parent s2:
    //
    //   (a) s1 is the only child of s2 (firstChild == lastChild), and
    //   (b) the columns of s1 and s2 have the same pattern.
    //
    // Condition (b) needs no pattern comparison, only arithmetic on sizes we already
    // hold. If s1 and s2 merge, then the rows below s1 are exactly the columns of s2
    // followed by the rows below s2, so
    //
    //   updateSize[s1] == frontSize[s2] + updateSize[s2]
    //
    // which is the test below, written as an addition so the unsigned arithmetic cannot
    // wrap. 0.9 writes it as updateSize[s2] == updateSize[s1] - 1, which is this same
    // identity with frontSize[s2] == 1, true only while supernodes are trivial. 10.12
    // does generalize it, but subtracts the child's front size where the derivation
    // calls for the parent's; that is latent (both are 1 on the nodal input it always
    // receives) but wrong, so we do not follow it. In the general form this function is
    // safe to run on an already-supernodal forest, which is what both references claim
    // for theirs.
    //
    // The scan runs in increasing order, so supOldToNew[s1] is always settled before s2
    // needs it: a child's label is smaller than its parent's.
    std::vector<std::int32_t> supOldToNew(supSize, NIL);
    std::size_t newSupSize = 0;
    for (std::size_t s2 = 0; s2 < supSize; ++s2) {
        const std::int32_t s1 = firstChild[s2];
        const bool onlyChild = (s1 != NIL) && (s1 == lastChild[s2]);
        const bool samePattern = onlyChild
            && (updateSize[s2] + frontSize[s2] == updateSize[static_cast<std::size_t>(s1)]);

        if (samePattern)
            supOldToNew[s2] = supOldToNew[static_cast<std::size_t>(s1)];   // s2 joins s1
        else
            supOldToNew[s2] = static_cast<std::int32_t>(newSupSize++);     // s2 starts anew
    }

    // Carry the column map through the merge. Note the columns of a supernode need not
    // be contiguous: the forest is in topological order, not necessarily postorder, so
    // merged supernodes can be non-consecutive. Consumers gather the front indices from
    // the map rather than assuming a contiguous range.
    for (std::size_t lc = 0; lc < size; ++lc)
        idxToSupIdx[lc] =
            supOldToNew[static_cast<std::size_t>(idxToSupIdx[lc])];

    // The parent links and the sizes of the merged supernodes. Scanning in decreasing
    // order, the first old supernode seen for a given new one is its highest-labeled
    // member, which is the top of the merged path: its parent link is the one that
    // leaves the merged supernode, and its update rows are the merged supernode's. Front
    // sizes instead accumulate over every member, so that runs before the mark test.
    std::vector<std::int32_t> newParent(newSupSize, NIL);
    std::vector<std::size_t>  newFrontSize(newSupSize, 0);
    std::vector<std::size_t>  newUpdateSize(newSupSize, 0);
    std::vector<bool>         done(newSupSize, false);

    for (std::size_t t = supSize; t > 0; --t) {
        const std::size_t s1 = t - 1;
        const std::size_t s1New = static_cast<std::size_t>(supOldToNew[s1]);

        newFrontSize[s1New] += frontSize[s1];

        if (done[s1New])
            continue;

        const std::int32_t s2 = parent[s1];
        if (s2 != NIL)
            newParent[s1New] = supOldToNew[static_cast<std::size_t>(s2)];
        newUpdateSize[s1New] = updateSize[s1];
        done[s1New] = true;
    }

    supSize = newSupSize;
    parent.swap(newParent);
    frontSize.swap(newFrontSize);
    updateSize.swap(newUpdateSize);

    // Rebuild the child, sibling and root links over the merged supernodes.
    finalizeLinks(supSize, parent, firstChild, lastChild, nextSibling, previousSibling,
                  numTrees, firstRoot, lastRoot);
}

template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
