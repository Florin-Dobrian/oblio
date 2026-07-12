#include "oblio/ElmForestEngine.h"

#include <deque>

namespace Oblio {

// Adapter: the forest needs only the sparsity pattern, so the matrix overload pulls it out
// and forwards. The implementation below is free of Val and compiled once.
template<class Val>
bool ElmForestEngine::compute(const SparseMatrix<Val>& A, const Permutation& p,
                              ElmForest& f) const {
    return compute(A.colPtr(), A.rowIdx(), p, f);
}

bool ElmForestEngine::compute(const std::vector<std::size_t>&  colPtr,
                              const std::vector<std::int32_t>& rowIdx,
                              const Permutation& p, ElmForest& f) const {
    if (colPtr.empty())
        return false;
    const std::size_t size = colPtr.size() - 1;
    if (p.size() != size)
        return false;

    computeParent(colPtr, rowIdx, p, f);

    // Start nodal: one trivial supernode per column, so the supernode count equals the
    // column count and the map is the identity. compressFundamental below merges these
    // into fundamental supernodes.
    f.mSupSize = size;
    f.mIdxToSupIdx.resize(size);
    for (std::size_t j = 0; j < size; ++j)
        f.mIdxToSupIdx[j] = static_cast<std::int32_t>(j);

    // The child, sibling and root links, then the front and update sizes. Both must
    // precede compression, whose merge test reads a column's only-child status from the
    // links and its sparsity pattern from the sizes.
    finalizeLinks(f);
    computeFrontAndUpdateSizes(colPtr, rowIdx, p, f);

    // Merge the columns into fundamental supernodes, unless asked to stay nodal. This
    // rewrites the map, the parent links and the sizes, leaving the child/sibling links
    // stale, so they are rebuilt straight after over the supernodes. finalizeLinks
    // therefore runs twice: once on the nodal forest, whose links the merge test reads,
    // and once on the compressed one.
    if (mSupernodes == Supernodes::Fundamental) {
        compressFundamental(f);
        finalizeLinks(f);
    }

    // Height last, so it is computed once, on the final forest. Compression shortens the
    // trees, every merged chain collapsing to a single level, so this cannot be carried
    // over from before the merge.
    f.mHeight = computeHeight(f);

    return true;
}

void ElmForestEngine::computeParent(
        const std::vector<std::size_t>&  colPtr,
        const std::vector<std::int32_t>& rowIdx,
        const Permutation& p,
        ElmForest& f) const {
    const std::size_t size = p.size();
    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    f.mSize = size;
    f.mParent.assign(size, NIL);
    std::vector<std::int32_t> ancestor(size, NIL);

    // For each factor column lc2 (increasing), look at its neighbours mapping to earlier
    // columns lc1 < lc2; path-compress to attach lc1's subtree under lc2.
    for (std::size_t lc2 = 0; lc2 < size; ++lc2) {
        const std::size_t ac2 = static_cast<std::size_t>(newToOld[lc2]);
        for (std::size_t sp = colPtr[ac2]; sp < colPtr[ac2 + 1]; ++sp) {
            const std::size_t lc1 =
                static_cast<std::size_t>(oldToNew[static_cast<std::size_t>(rowIdx[sp])]);
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
                f.mParent[lc3] = static_cast<std::int32_t>(lc2);
            }
        }
    }
}

void ElmForestEngine::finalizeLinks(ElmForest& f) const {
    f.mFirstChild.assign(f.mSupSize, NIL);
    f.mLastChild.assign(f.mSupSize, NIL);
    f.mNextSibling.assign(f.mSupSize, NIL);
    f.mPreviousSibling.assign(f.mSupSize, NIL);
    f.mNumTrees  = 0;
    f.mFirstRoot = NIL;
    f.mLastRoot  = NIL;

    // For each supernode in decreasing order, front-insert it into its f.mParent's
    // child list, or into the root list if it has no f.mParent. Front-insertion of
    // decreasing labels leaves both the child lists and the root list in
    // increasing label order.
    for (std::size_t t = f.mSupSize; t > 0; --t) {
        const std::int32_t s1 = static_cast<std::int32_t>(t - 1);
        const std::int32_t s2 = f.mParent[t - 1];
        if (s2 == NIL) {
            // s1 is a tree root.
            f.mNextSibling[t - 1] = f.mFirstRoot;
            if (f.mFirstRoot == NIL)
                f.mLastRoot = s1;
            else
                f.mPreviousSibling[static_cast<std::size_t>(f.mFirstRoot)] = s1;
            f.mFirstRoot = s1;
            ++f.mNumTrees;
        } else {
            // s1 becomes a child of s2.
            const std::size_t us2 = static_cast<std::size_t>(s2);
            f.mNextSibling[t - 1] = f.mFirstChild[us2];
            if (f.mFirstChild[us2] == NIL)
                f.mLastChild[us2] = s1;
            else
                f.mPreviousSibling[static_cast<std::size_t>(f.mFirstChild[us2])] = s1;
            f.mFirstChild[us2] = s1;
        }
    }
}

std::size_t ElmForestEngine::computeHeight(const ElmForest& f) const {
    // Breadth-first from the roots: a supernode is processed after its f.parent(), so
    // its depth is the f.parent()'s depth plus one, and the height is the largest
    // depth plus one. The root chain and each child list are walked backward via
    // f.previousSibling(); order is irrelevant to a maximum.
    std::vector<std::size_t> depth(f.supSize(), 0);
    std::deque<std::int32_t> queue;
    for (std::int32_t s = f.lastRoot(); s != NIL; s = f.previousSibling()[static_cast<std::size_t>(s)])
        queue.push_back(s);
    while (!queue.empty()) {
        const std::int32_t s2 = queue.front();
        queue.pop_front();
        const std::size_t us2 = static_cast<std::size_t>(s2);
        if (f.parent()[us2] != NIL)
            depth[us2] = depth[static_cast<std::size_t>(f.parent()[us2])] + 1;
        for (std::int32_t s1 = f.lastChild()[us2]; s1 != NIL;
             s1 = f.previousSibling()[static_cast<std::size_t>(s1)])
            queue.push_back(s1);
    }
    std::size_t height = 0;
    for (std::size_t s = 0; s < f.supSize(); ++s)
        if (height < depth[s])
            height = depth[s];
    if (f.supSize() > 0)
        ++height;
    return height;
}

void ElmForestEngine::computeFrontAndUpdateSizes(
        const std::vector<std::size_t>&  colPtr,
        const std::vector<std::int32_t>& rowIdx,
        const Permutation& p,
        ElmForest& f) const {
    const std::size_t size = f.mSize;
    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // Front size: the columns mapping to each supernode. The forest is nodal here, so the
    // map is the identity and this is all 1s; counting derives that rather than asserting
    // it, but buys no generality, since compression derives the merged front sizes itself.
    f.mFrontSize.assign(f.mSupSize, 0);
    for (std::size_t j = 0; j < size; ++j)
        ++f.mFrontSize[static_cast<std::size_t>(f.mIdxToSupIdx[j])];

    // Update size: the subdiagonal nonzero count of each L column, by the pruned-
    // row-subtree walk. For factor column lc2, climb from each earlier neighbour
    // lc1 up the etree toward lc2, marking to avoid double counting; every column
    // passed below the root gains one update index. This is 0.9's columnSize with
    // the diagonal left out, so it yields f.mUpdateSize directly (== columnSize - 1).
    f.mUpdateSize.assign(f.mSupSize, 0);
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
                ++f.mUpdateSize[lc3];
                mark[lc3] = static_cast<std::int32_t>(lc2);
                lc3 = static_cast<std::size_t>(f.mParent[lc3]);
            }
        }
    }
}

void ElmForestEngine::compressFundamental(ElmForest& f) const {
    // Precondition: f is nodal, so a supernode is a column, and the code names them as the
    // columns they are, j and k. The map is the identity and every front size is 1.
    const std::size_t size = f.mSize;

    // Assign every column to a supernode. Column k continues its child j's supernode when
    // the two form a fundamental supernode, otherwise it starts a new one. Increasing
    // order, so snode[j] is settled before k needs it: a child's column is numbered below
    // its parent's.
    std::vector<std::int32_t> snode(size, NIL);
    std::size_t numSnodes = 0;
    for (std::size_t k = 0; k < size; ++k) {
        const std::int32_t j = f.mFirstChild[k];

        // Does k merge into its child j? The two clauses after the guard are the
        // conditions of a fundamental supernode: j is the only child of k, and the two
        // share one sparsity pattern.
        //
        // The guard is not decoration. firstChild == lastChild also holds when both are
        // NIL, that is at every leaf, so without it a leaf would report exactly one child
        // and the pattern test would subscript updateSize with NIL. Short-circuit order
        // matters: the guard protects the clauses to its right, which is a property of
        // this expression, not something a later test could restore.
        //
        // The pattern test is |update(j)| = |front(k)| + |update(k)|, written as an
        // addition so the unsigned arithmetic cannot wrap. Every front size is 1 here, so
        // it is 0.9's |update(k)| == |update(j)| - 1; we keep the identity in full because
        // it is the form derived in the theory, and because a bare + 1 would go silently
        // wrong if the precondition were ever broken.
        const bool merge = (j != NIL)
            && (j == f.mLastChild[k])
            && (f.mUpdateSize[k] + f.mFrontSize[k]
                    == f.mUpdateSize[static_cast<std::size_t>(j)]);

        if (merge)
            snode[k] = snode[static_cast<std::size_t>(j)];        // k continues j's supernode
        else
            snode[k] = static_cast<std::int32_t>(numSnodes++);    // k starts a new one
    }

    // The parent links and the sizes of the supernodes. Scanning columns in decreasing
    // order, the first column seen for a given supernode is its topmost: its parent link is
    // the one that leaves the supernode, and its rows below are the supernode's update
    // rows. Front sizes instead accumulate over every column, so that runs before the test.
    std::vector<std::int32_t> parent(numSnodes, NIL);
    std::vector<std::size_t>  frontSize(numSnodes, 0);
    std::vector<std::size_t>  updateSize(numSnodes, 0);
    std::vector<bool>         done(numSnodes, false);

    for (std::size_t t = size; t > 0; --t) {
        const std::size_t k = t - 1;
        const std::size_t s = static_cast<std::size_t>(snode[k]);

        ++frontSize[s];   // every column of s adds one to its front size

        if (done[s])
            continue;

        const std::int32_t p = f.mParent[k];
        if (p != NIL)
            parent[s] = snode[static_cast<std::size_t>(p)];
        updateSize[s] = f.mUpdateSize[k];   // the rows below k are the rows below s
        done[s] = true;
    }

    f.mSupSize = numSnodes;
    f.mIdxToSupIdx.swap(snode);   // the map was the identity, and is now column -> supernode
    f.mParent.swap(parent);
    f.mFrontSize.swap(frontSize);
    f.mUpdateSize.swap(updateSize);

    // mFirstChild, mLastChild, mNextSibling, mPreviousSibling and the roots still describe
    // the nodal forest and are now stale. The caller rebuilds them with finalizeLinks.
}

template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
