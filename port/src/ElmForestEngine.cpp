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
    computeColumnSizes(colPtr, rowIdx, p, f);

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

    // Naming: an l prefix is the factor (permuted) ordering, an a prefix the original
    // ordering of A, and j and k carry the column roles, j the lower and k the higher. So lk
    // is the factor column being processed, lj an earlier neighbour of it (lj < lk, enforced
    // by the guard below), and ak is lk's column in A. The same prefix marks a position in
    // the storage: ap is an offset into A's flat arrays, which is 0.9's p (the name p being
    // taken here by the Permutation). Positions measure rather than name, so ap is a
    // std::size_t and cannot be NIL, unlike the columns. Beyond those, r is the node the
    // climb has reached and t the temporary that saves the next hop before compression
    // overwrites it.
    //
    // For each factor column lk (increasing), look at its neighbours mapping to earlier
    // columns lj < lk; path-compress to attach lj's subtree under lk.
    const std::int32_t n = static_cast<std::int32_t>(size);
    for (std::int32_t lk = 0; lk < n; ++lk) {
        const std::int32_t ak = newToOld[lk];
        for (std::size_t ap = colPtr[ak]; ap < colPtr[ak + 1]; ++ap) {
            const std::int32_t lj = oldToNew[rowIdx[ap]];
            if (lj >= lk)
                continue;   // later column or the diagonal itself
            std::int32_t r = lj;
            while (ancestor[r] != NIL && ancestor[r] != lk) {
                const std::int32_t t = ancestor[r];
                ancestor[r] = lk;
                r = t;
            }
            if (ancestor[r] == NIL) {
                ancestor[r] = lk;
                f.mParent[r] = lk;
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

    // For each supernode in decreasing order, front-insert it into its parent's child list,
    // or into the root list if it has no parent. Front-insertion of decreasing labels leaves
    // both the child lists and the root list in increasing label order.
    //
    // The loop counts down rather than indexing down. A std::size_t cannot express a negative
    // guard, so `for (jj = supSize - 1; jj >= 0; --jj)` never terminates: jj wraps to SIZE_MAX
    // instead of going below zero. Counting the supernodes remaining, from supSize down to 1,
    // keeps the guard in a range std::size_t can represent, and the index jj = t - 1 is then
    // always in bounds.
    for (std::size_t t = f.mSupSize; t > 0; --t) {
        const std::int32_t jj = static_cast<std::int32_t>(t - 1);   // supernode inserted
        const std::int32_t kk = f.mParent[jj];                      // its parent, or NIL

        if (kk == NIL) {
            // jj is a tree root.
            f.mNextSibling[jj] = f.mFirstRoot;
            if (f.mFirstRoot == NIL)
                f.mLastRoot = jj;
            else
                f.mPreviousSibling[f.mFirstRoot] = jj;
            f.mFirstRoot = jj;
            ++f.mNumTrees;
        } else {
            // jj becomes a child of kk.
            f.mNextSibling[jj] =
                f.mFirstChild[kk];
            if (f.mFirstChild[kk] == NIL)
                f.mLastChild[kk] = jj;
            else
                f.mPreviousSibling[f.mFirstChild[kk]] = jj;
            f.mFirstChild[kk] = jj;
        }
    }
}

std::size_t ElmForestEngine::computeHeight(const ElmForest& f) const {
    // Breadth-first from the roots: a supernode is processed after its parent, so its depth
    // is the parent's depth plus one, and the height is the largest depth plus one. The root
    // chain and each child list are walked backward via previousSibling; order is irrelevant
    // to a maximum.
    //
    // kk is the supernode being visited and jj each of its children, following the house
    // convention that a doubled letter names a supernode.
    std::vector<std::size_t> depth(f.mSupSize, 0);
    std::deque<std::int32_t> queue;
    for (std::int32_t kk = f.mLastRoot; kk != NIL; kk = f.mPreviousSibling[kk])
        queue.push_back(kk);
    while (!queue.empty()) {
        const std::int32_t kk = queue.front();
        queue.pop_front();
        if (f.mParent[kk] != NIL)
            depth[kk] = depth[f.mParent[kk]] + 1;
        for (std::int32_t jj = f.mLastChild[kk]; jj != NIL; jj = f.mPreviousSibling[jj])
            queue.push_back(jj);
    }
    std::size_t height = 0;
    for (std::size_t d : depth)
        if (height < d)
            height = d;
    if (f.mSupSize > 0)
        ++height;
    return height;
}

void ElmForestEngine::computeColumnSizes(
        const std::vector<std::size_t>&  colPtr,
        const std::vector<std::int32_t>& rowIdx,
        const Permutation& p,
        ElmForest& f) const {
    const std::size_t size = f.mSize;
    const std::int32_t n = static_cast<std::int32_t>(size);
    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // Front size. The forest is nodal here, so a supernode is a column and owns exactly
    // itself. Nothing to compute.
    f.mFrontSize.assign(size, 1);

    // Update size: the subdiagonal nonzero count of each column of L, by the pruned-row-
    // subtree walk. Fix a factor column lk. The columns holding a nonzero in row lk are
    // exactly those on the forest paths climbing from each earlier neighbour lj of lk up
    // to lk, so climb from each such lj and credit every column passed.
    //
    // The marker does two jobs. It prevents double counting, since two neighbours of lk may
    // climb into a shared upper path and the second must not credit it twice. And it halts
    // the climb: mark[lk] is set to lk before the neighbour loop, so a climb reaching lk
    // stops there without a separate r == lk test. That is the same trick symbolic
    // factorization uses to make a union idempotent, with the action reduced from storing
    // an index to counting it.
    //
    // 0.9 counts the diagonal too, giving |Struct(j)|, then subtracts one per column to get
    // the update size. Starting the count at zero yields the update size directly, which is
    // what the forest wants, and saves both the array and the pass.
    f.mUpdateSize.assign(size, 0);
    std::vector<std::int32_t> mark(size, NIL);
    for (std::int32_t lk = 0; lk < n; ++lk) {
        mark[lk] = lk;
        const std::int32_t ak = newToOld[lk];
        for (std::size_t ap = colPtr[ak]; ap < colPtr[ak + 1]; ++ap) {
            const std::int32_t lj = oldToNew[rowIdx[ap]];
            if (lj >= lk)
                continue;   // later column or the diagonal itself
            std::int32_t r = lj;
            while (mark[r] != lk) {
                ++f.mUpdateSize[r];   // column r gains row lk
                mark[r] = lk;
                r = f.mParent[r];
            }
        }
    }
}

void ElmForestEngine::compressFundamental(ElmForest& f) const {
    // Precondition: f is nodal, so a supernode is a column, and the code names them as the
    // columns they are, j and k. The map is the identity and every front size is 1.
    const std::size_t size = f.mSize;

    // Assign every column to a supernode. Column k continues its child j's supernode when the
    // two form a fundamental supernode, otherwise it starts a new one. Increasing order, so
    // idxToSupIdx[j] is settled before k needs it: a child's column is numbered below its
    // parent's.
    std::vector<std::int32_t> idxToSupIdx(size, NIL);
    std::int32_t supSize = 0;
    const std::int32_t n = static_cast<std::int32_t>(size);
    for (std::int32_t k = 0; k < n; ++k) {
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
        // The pattern test is |update(j)| = |front(k)| + |update(k)|, written as an addition
        // rather than a subtraction. A subtraction (0.9's |update(j)| - 1) is unsigned and
        // would wrap if its right side ever exceeded its left; the addition has nothing to
        // wrap. Nor can the addition overflow: the sum is the size of k's index set, which
        // is a subset of the factor's rows, so it is at most size, a number we already hold.
        //
        // Every front size is 1 here, so this is 0.9's |update(k)| == |update(j)| - 1. We
        // keep the identity in full because it is the form derived in the theory, and
        // because a bare + 1 would go silently wrong if the precondition were ever broken.
        const bool merge = (j != NIL)
            && (j == f.mLastChild[k])
            && (f.mFrontSize[k] + f.mUpdateSize[k]
                    == f.mUpdateSize[j]);

        if (merge)
            idxToSupIdx[k] = idxToSupIdx[j];   // k continues j's
        else
            idxToSupIdx[k] = supSize++;                                  // k starts a new one
    }

    // The parent links and the sizes of the supernodes. Scanning columns in decreasing
    // order, the first column seen for a given supernode is its topmost: its parent link is
    // the one that leaves the supernode, and its rows below are the supernode's update
    // rows. Front sizes instead accumulate over every column, so that runs before the test.
    const std::size_t numSup = static_cast<std::size_t>(supSize);
    std::vector<std::int32_t> parent(numSup, NIL);
    std::vector<std::size_t>  frontSize(numSup, 0);
    std::vector<std::size_t>  updateSize(numSup, 0);
    std::vector<bool>         seen(numSup, false);

    // Counting down, not indexing down: see finalizeLinks for why a std::size_t descending
    // loop must be written this way.
    //
    // Naming follows 0.9: a single letter is a column, a doubled one is that column's
    // supernode, so the doubling is the map applied. Column j is the one being scanned and
    // jj is its supernode; k is j's parent column and kk is the supernode that lands in.
    // The column convention j < k survives, since a parent's column is numbered above its
    // child's, and the assignment then reads as the fact it derives: the supernode of j's
    // parent column is the parent of j's supernode.
    for (std::size_t t = size; t > 0; --t) {
        const std::int32_t j  = static_cast<std::int32_t>(t - 1);   // column scanned
        const std::int32_t jj = idxToSupIdx[j];                     // its supernode

        ++frontSize[jj];   // every column of jj adds one to it

        if (seen[jj])
            continue;      // jj's topmost column came earlier; the rest is taken from it

        const std::int32_t k = f.mParent[j];
        if (k != NIL) {
            const std::int32_t kk = idxToSupIdx[k];
            parent[jj] = kk;
        }
        // the rows below j are the rows below jj
        updateSize[jj] = f.mUpdateSize[j];
        seen[jj] = true;
    }

    f.mSupSize = numSup;
    f.mIdxToSupIdx.swap(idxToSupIdx);   // was the identity, now column -> supernode
    f.mParent.swap(parent);
    f.mFrontSize.swap(frontSize);
    f.mUpdateSize.swap(updateSize);

    // mFirstChild, mLastChild, mNextSibling, mPreviousSibling and the roots still describe
    // the nodal forest and are now stale. The caller rebuilds them with finalizeLinks.
}

template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
