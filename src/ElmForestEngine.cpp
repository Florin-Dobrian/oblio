#include "oblio/ElmForestEngine.h"

#include <deque>

namespace Oblio {

// Adapter: the forest needs only the sparsity pattern, so the matrix overload pulls it out
// and forwards. The implementation below is free of Val and compiled once.
template<class Val>
bool ElmForestEngine::compute(const SparseMatrix<Val>& A, const Permutation& p,
                              ElmForest& ef) const {
    return compute(A.colPtr(), A.rowIdx(), p, ef);
}

bool ElmForestEngine::compute(const std::vector<std::size_t>&  colPtr,
                              const std::vector<std::int32_t>& rowIdx,
                              const Permutation& p, ElmForest& ef) const {
    if (colPtr.empty())
        return false;
    const std::size_t size = colPtr.size() - 1;
    if (p.size() != size)
        return false;

    computeParent(colPtr, rowIdx, p, ef);

    // Start nodal: one trivial supernode per column, so the supernode count equals the
    // column count and the map is the identity. compressFundamental below merges these
    // into fundamental supernodes.
    ef.mSnodeSize = size;
    ef.mNodeToSnode.resize(size);
    for (std::int32_t lj = 0; lj < static_cast<std::int32_t>(size); ++lj)
        ef.mNodeToSnode[lj] = lj;

    // The child, sibling and root links, then the front and update sizes. Both must
    // precede compression, whose merge test reads a column's only-child status from the
    // links and its sparsity pattern from the sizes.
    finalizeLinks(ef);
    computeColumnSizes(colPtr, rowIdx, p, ef);

    // Merge the columns into fundamental supernodes, unless asked to stay nodal, then
    // amalgamate if a threshold is set. Each rewrites the map, the parent links and the
    // sizes, leaving the child/sibling links stale, so they are rebuilt straight after.
    // Fundamental first: it does the free, canonical merging cheaply, leaving amalgamation
    // the tie-broken and the paid work on a smaller forest.
    if (mSupernodes == Supernodes::Fundamental) {
        compressFundamental(ef);
        finalizeLinks(ef);
    }
    if (mThreshold.has_value()) {
        compressThreshold(ef, *mThreshold);
        finalizeLinks(ef);
    }

    // Height last, so it is computed once, on the final forest. Compression shortens the
    // trees, every merged chain collapsing to a single level, so this cannot be carried
    // over from before the merge.
    ef.mHeight = computeHeight(ef);

    return true;
}

void ElmForestEngine::computeParent(
        const std::vector<std::size_t>&  colPtr,
        const std::vector<std::int32_t>& rowIdx,
        const Permutation& p,
        ElmForest& ef) const {
    const std::size_t size = p.size();
    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    ef.mSize = size;
    ef.mParent.assign(size, NIL);
    std::vector<std::int32_t> ancestor(size, NIL);

    // Naming: an l prefix is the factor (permuted) ordering, an a prefix the original
    // ordering of A, and j and k carry the column roles, j the lower and k the higher. So lk
    // is the factor column being processed, lj an earlier neighbor of it (lj < lk, enforced
    // by the guard below), and ak is lk's column in A. The same prefix marks a position in
    // the storage: cp is an offset into A's flat arrays, which is 0.9's p (the name p being
    // taken here by the Permutation). Positions measure rather than name, so cp is a
    // std::size_t and cannot be NIL, unlike the columns. Beyond those, r is the node the
    // climb has reached and t the temporary that saves the next hop before compression
    // overwrites it.
    //
    // For each factor column lk (increasing), look at its neighbors mapping to earlier
    // columns lj < lk; path-compress to attach lj's subtree under lk.
    for (std::int32_t lk = 0; lk < static_cast<std::int32_t>(size); ++lk) {
        const std::int32_t ak = newToOld[lk];
        for (std::size_t cp = colPtr[ak]; cp < colPtr[ak + 1]; ++cp) {
            const std::int32_t aj = rowIdx[cp];       // a row of A, hence a column of L
            const std::int32_t lj = oldToNew[aj];    // the same column, in L
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
                ef.mParent[r] = lk;
            }
        }
    }
}

void ElmForestEngine::finalizeLinks(ElmForest& ef) const {
    ef.mFirstChild.assign(ef.mSnodeSize, NIL);
    ef.mLastChild.assign(ef.mSnodeSize, NIL);
    ef.mNextSibling.assign(ef.mSnodeSize, NIL);
    ef.mPreviousSibling.assign(ef.mSnodeSize, NIL);
    ef.mNumTrees  = 0;
    ef.mFirstRoot = NIL;
    ef.mLastRoot  = NIL;

    // For each supernode in decreasing order, front-insert it into its parent's child list,
    // or into the root list if it has no parent. Front-insertion of decreasing labels leaves
    // both the child lists and the root list in increasing label order.
    for (std::int32_t jj = static_cast<std::int32_t>(ef.mSnodeSize) - 1; jj >= 0; --jj) {
        const std::int32_t kk = ef.mParent[jj];   // its parent, or NIL

        if (kk == NIL) {
            // jj is a tree root.
            ef.mNextSibling[jj] = ef.mFirstRoot;
            if (ef.mFirstRoot == NIL)
                ef.mLastRoot = jj;
            else
                ef.mPreviousSibling[ef.mFirstRoot] = jj;
            ef.mFirstRoot = jj;
            ++ef.mNumTrees;
        } else {
            // jj becomes a child of kk.
            ef.mNextSibling[jj] =
                ef.mFirstChild[kk];
            if (ef.mFirstChild[kk] == NIL)
                ef.mLastChild[kk] = jj;
            else
                ef.mPreviousSibling[ef.mFirstChild[kk]] = jj;
            ef.mFirstChild[kk] = jj;
        }
    }
}

std::size_t ElmForestEngine::computeHeight(const ElmForest& ef) const {
    // Breadth-first from the roots: a supernode is processed after its parent, so its depth
    // is the parent's depth plus one, and the height is the largest depth plus one. The root
    // chain and each child list are walked backward via previousSibling; order is irrelevant
    // to a maximum.
    //
    // kk is the supernode being visited and jj each of its children, following the house
    // convention that a doubled letter names a supernode.
    std::vector<std::size_t> depth(ef.mSnodeSize, 0);
    std::deque<std::int32_t> queue;
    for (std::int32_t kk = ef.mLastRoot; kk != NIL; kk = ef.mPreviousSibling[kk])
        queue.push_back(kk);
    while (!queue.empty()) {
        const std::int32_t kk = queue.front();
        queue.pop_front();
        if (ef.mParent[kk] != NIL)
            depth[kk] = depth[ef.mParent[kk]] + 1;
        for (std::int32_t jj = ef.mLastChild[kk]; jj != NIL; jj = ef.mPreviousSibling[jj])
            queue.push_back(jj);
    }
    std::size_t height = 0;
    for (std::size_t d : depth)
        if (height < d)
            height = d;
    if (ef.mSnodeSize > 0)
        ++height;
    return height;
}

void ElmForestEngine::computeColumnSizes(
        const std::vector<std::size_t>&  colPtr,
        const std::vector<std::int32_t>& rowIdx,
        const Permutation& p,
        ElmForest& ef) const {
    const std::size_t size = ef.mSize;
    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // Front size. The forest is nodal here, so a supernode is a column and owns exactly
    // itself. Nothing to compute.
    ef.mFrontSize.assign(size, 1);

    // Update size: the subdiagonal nonzero count of each column of L, by the pruned-row-
    // subtree walk. Fix a factor column lk. The columns holding a nonzero in row lk are
    // exactly those on the forest paths climbing from each earlier neighbor lj of lk up
    // to lk, so climb from each such lj and credit every column passed.
    //
    // The marker does two jobs. It prevents double counting, since two neighbors of lk may
    // climb into a shared upper path and the second must not credit it twice. And it halts
    // the climb: mark[lk] is set to lk before the neighbor loop, so a climb reaching lk
    // stops there without a separate r == lk test. That is the same trick symbolic
    // factorization uses to make a union idempotent, with the action reduced from storing
    // an index to counting it.
    //
    // 0.9 counts the diagonal too, giving |Struct(j)|, then subtracts one per column to get
    // the update size. Starting the count at zero yields the update size directly, which is
    // what the forest wants, and saves both the array and the pass.
    ef.mUpdateSize.assign(size, 0);
    std::vector<std::int32_t> mark(size, NIL);
    for (std::int32_t lk = 0; lk < static_cast<std::int32_t>(size); ++lk) {
        mark[lk] = lk;
        const std::int32_t ak = newToOld[lk];
        for (std::size_t cp = colPtr[ak]; cp < colPtr[ak + 1]; ++cp) {
            const std::int32_t aj = rowIdx[cp];       // a row of A, hence a column of L
            const std::int32_t lj = oldToNew[aj];    // the same column, in L
            if (lj >= lk)
                continue;   // later column or the diagonal itself
            std::int32_t r = lj;
            while (mark[r] != lk) {
                ++ef.mUpdateSize[r];   // column r gains row lk
                mark[r] = lk;
                r = ef.mParent[r];
            }
        }
    }
}

void ElmForestEngine::compressFundamental(ElmForest& ef) const {
    // Precondition: ef is nodal, so a supernode is a column, and the code names them as the
    // columns they are, lj and lk. The map is the identity and every front size is 1.
    const std::size_t size = ef.mSize;

    // Assign every column to a supernode. Column lk continues its child lj's supernode when the
    // two form a fundamental supernode, otherwise it starts a new one. Increasing order, so
    // newNodeToSnode[lj] is settled before lk needs it: a child's column is numbered below its
    // newParent's.
    std::vector<std::int32_t> newNodeToSnode(size, NIL);
    std::int32_t newSnodeSize = 0;
    for (std::int32_t lk = 0; lk < static_cast<std::int32_t>(size); ++lk) {
        const std::int32_t lj = ef.mFirstChild[lk];

        // Does lk merge into its child lj? The two clauses after the guard are the
        // conditions of a fundamental supernode: lj is the only child of lk, and the two
        // share one sparsity pattern.
        //
        // The guard is not decoration. firstChild == lastChild also holds when both are
        // NIL, that is at every leaf, so without it a leaf would report exactly one child
        // and the pattern test would subscript newUpdateSize with NIL. Short-circuit order
        // matters: the guard protects the clauses to its right, which is a property of
        // this expression, not something a later test could restore.
        //
        // The pattern test, in the notes' notation (Section 4.2), is
        //
        //     |update(J)| = |front(K)| + |update(K)|
        //
        // for a child supernode J and its newParent K. The forest is nodal here, so a supernode is
        // a column: J is lj and K is lk. It is written as an addition rather than a subtraction.
        // A subtraction (0.9's |update(J)| - 1) is unsigned and would wrap if its right side
        // ever exceeded its left; the addition has nothing to wrap. Nor can it overflow: the sum
        // is the size of lk's index set, a subset of the factor's rows, so at most size.
        //
        // Every front size is 1 here, so this reduces to 0.9's |update(K)| == |update(J)| - 1.
        // We keep the identity in full because it is the form derived in the theory, and because
        // a bare + 1 would go silently wrong if the precondition were ever broken.
        const bool merge = (lj != NIL)
            && (lj == ef.mLastChild[lk])
            && (ef.mFrontSize[lk] + ef.mUpdateSize[lk]
                    == ef.mUpdateSize[lj]);

        if (merge)
            newNodeToSnode[lk] = newNodeToSnode[lj];   // lk continues lj's
        else
            newNodeToSnode[lk] = newSnodeSize++;        // lk starts a new one
    }

    // The newParent links and the sizes of the supernodes. Scanning columns in decreasing
    // order, the first column seen for a given supernode is its topmost: its newParent link is
    // the one that leaves the supernode, and its rows below are the supernode's update
    // rows. Front sizes instead accumulate over every column, so that runs before the test.
    std::vector<std::int32_t> newParent(newSnodeSize, NIL);
    std::vector<std::size_t>  newFrontSize(newSnodeSize, 0);
    std::vector<std::size_t>  newUpdateSize(newSnodeSize, 0);
    std::vector<bool>         seen(newSnodeSize, false);

    // Counting down, not indexing down: see finalizeLinks for why a std::size_t descending
    // loop must be written this way.
    //
    // Naming follows 0.9: a single letter is a column, a doubled one is that column's
    // supernode, so the doubling is the map applied. Column lj is the one being scanned and
    // jj is its supernode; lk is lj's parent column and kk is the supernode that lands in.
    // The column convention lj < lk survives, since a parent's column is numbered above its
    // child's, and the assignment then reads as the fact it derives: the supernode of lj's
    // parent column is the parent of lj's supernode.
    for (std::int32_t lj = static_cast<std::int32_t>(size) - 1; lj >= 0; --lj) {
        const std::int32_t jj = newNodeToSnode[lj];                     // its supernode

        // Every column of jj adds one to its front. The forest here is nodal, so each column's
        // own front size is 1, left implicit in the `++` (threshold compression, running on a
        // supernodal forest, must instead add the child's actual front size). And unlike there,
        // these sizes are only tallied into the output, never read back to drive a decision, so
        // no mutable working copy of the front sizes is needed.
        ++newFrontSize[jj];

        if (seen[jj])
            continue;      // jj's topmost column came earlier; the rest is taken from it

        const std::int32_t lk = ef.mParent[lj];
        if (lk != NIL) {
            const std::int32_t kk = newNodeToSnode[lk];
            newParent[jj] = kk;
        }
        // the rows below lj are the rows below jj
        newUpdateSize[jj] = ef.mUpdateSize[lj];
        seen[jj] = true;
    }

    ef.mSnodeSize = newSnodeSize;   // the next free label is also the count
    ef.mNodeToSnode.swap(newNodeToSnode);   // was the identity, now column -> supernode
    ef.mParent.swap(newParent);
    ef.mFrontSize.swap(newFrontSize);
    ef.mUpdateSize.swap(newUpdateSize);

    // mFirstChild, mLastChild, mNextSibling, mPreviousSibling and the roots still describe
    // the nodal forest and are now stale. The caller rebuilds them with finalizeLinks.
}

void ElmForestEngine::compressThreshold(ElmForest& ef, std::size_t threshold) const {
    const std::size_t size    = ef.mSize;
    const std::size_t snodeSize = ef.mSnodeSize;

    // Where each old supernode ends up. Identity to begin with; snodeOldToNew[jj] = kk records
    // that jj was absorbed into kk.
    std::vector<std::int32_t> snodeOldToNew(snodeSize);
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj)
        snodeOldToNew[jj] = jj;

    // Children still worth testing. A child that fails the budget can never pass it later,
    // since its parent only ever grows, so striking it off is permanent and keeps the greedy
    // loop from rescanning it.
    std::vector<bool> candidate(snodeSize, true);

    // A working copy of the front sizes, because they grow. When supernode kk absorbs a
    // child, kk's front widens, and a later parent must price kk by its *current* width, not
    // its original one. Parents run in increasing order, so a child has already been a parent
    // itself by the time we look at it. 10.12 omits this update (the line is in its source,
    // commented out, because it had made the array const), which silently understates both
    // the fill and the resulting block.
    std::vector<std::size_t>        frontSize = ef.mFrontSize;
    const std::vector<std::size_t>& updateSize = ef.mUpdateSize;

    std::size_t newSnodeSize = snodeSize;

    // Whether any merge actually stored a zero. At threshold zero none can, so the merged
    // supernodes still have exactly matching patterns and the forest stays "exact"; above it,
    // the first zero stored breaks that.
    bool storedAZero = false;

    // For every supernode kk, absorb as many of its children as the budget allows.
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        std::size_t kkZerosFillInc = 0;   // zeros already bought for kk
        std::size_t kkFrontSizeInc = 0;   // columns already absorbed into kk

        // kk's index set before any of this parent's merges. It grows by kkFrontSizeInc as we go.
        const std::size_t kkNumNodeIdx = frontSize[kk] + updateSize[kk];

        for (;;) {
            // The best child so far and its zero-fill. "Best" is the full rule below, least fill
            // then widest front, so bestZerosFill is the winner's fill, not simply the least on
            // offer: at a tie the wider-front child wins with the same fill.
            std::int32_t bestChild = NIL;
            std::size_t  bestZerosFill = 0;

            // Scanned back to front, from the last child to the first, using the doubly-linked
            // sibling chain. This is a free choice: it does not change the fill or the block sizes,
            // only which of several equal-cost partitions is produced. Scanning backward strands the
            // *lower*-numbered child of a full tie as its own supernode and merges the higher one, so
            // the supernode labels come out more nearly monotone in column order, which is easier to
            // reason about downstream. (Forward scanning would strand the higher child instead.)
            for (std::int32_t jj = ef.mLastChild[kk]; jj != NIL; jj = ef.mPreviousSibling[jj]) {
                if (!candidate[jj])
                    continue;

                // The zeros each column of jj must store once it joins kk's front: it held
                // update(jj) rows below itself, and must now hold the whole of kk's index
                // set. By the containment theorem update(jj) is a subset of that index set,
                // so this is a set-difference size and cannot go negative.
                const std::size_t zerosFillPerCol = kkNumNodeIdx + kkFrontSizeInc - updateSize[jj];
                const std::size_t zerosFill = zerosFillPerCol * frontSize[jj];   // this child alone

                if (kkZerosFillInc + zerosFill > threshold) {
                    candidate[jj] = false;   // over budget, and kk only grows from here
                    continue;
                }

                // Least fill; ties to the widest front; ties again to the first child *reached*,
                // which with the backward scan is the last in the list. Arbitrariness made
                // deterministic: a canonical algorithm would need no such rule.
                if (bestChild == NIL || zerosFill < bestZerosFill
                        || (zerosFill == bestZerosFill && frontSize[jj] > frontSize[bestChild])) {
                    bestChild = jj;
                    bestZerosFill = zerosFill;
                }
            }

            if (bestChild == NIL)
                break;   // nothing more fits the budget

            snodeOldToNew[bestChild] = kk;
            if (bestZerosFill > 0)
                storedAZero = true;   // this merge pays in explicitly stored zeros
            kkZerosFillInc += bestZerosFill;
            kkFrontSizeInc += frontSize[bestChild];
            candidate[bestChild] = false;
            --newSnodeSize;
        }

        frontSize[kk] += kkFrontSizeInc;   // kk has grown; later parents must see it
    }

    // Resolve chains: a supernode may have been absorbed into one that was itself absorbed.
    // Decreasing order, so an absorber (whose label exceeds its children's) is resolved before
    // the children that point at it.
    for (std::int32_t jj = static_cast<std::int32_t>(snodeSize) - 1; jj >= 0; --jj) {
        const std::int32_t kk = snodeOldToNew[jj];
        if (kk != jj)
            snodeOldToNew[jj] = snodeOldToNew[kk];
    }

    // Compact the labels: the survivors keep topological order, but their old labels have
    // gaps where absorbed supernodes used to be.
    std::vector<std::int32_t> label(snodeSize, NIL);
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj)
        label[snodeOldToNew[jj]] = 0;               // mark the survivors
    std::int32_t next = 0;
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj)
        if (label[jj] != NIL)
            label[jj] = next++;
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj)
        snodeOldToNew[jj] = label[snodeOldToNew[jj]];

    // Carry the column map through.
    std::vector<std::int32_t> newNodeToSnode(size);
    for (std::int32_t lj = 0; lj < static_cast<std::int32_t>(size); ++lj)
        newNodeToSnode[lj] = snodeOldToNew[ef.mNodeToSnode[lj]];

    // Parent links and sizes of the merged supernodes. Scanning old supernodes in decreasing
    // order, the first one seen for a given survivor is the absorber itself (its label is the
    // largest in the merged set): its parent link leaves the merged supernode, and its update
    // rows are the merged supernode's. Front sizes instead accumulate over every member, from
    // the *input* sizes, not the working copy, which already holds the accumulated total.
    std::vector<std::int32_t> newParent(newSnodeSize, NIL);
    std::vector<std::size_t>  newFrontSize(newSnodeSize, 0);
    std::vector<std::size_t>  newUpdateSize(newSnodeSize, 0);
    std::vector<bool>         seen(newSnodeSize, false);

    for (std::int32_t jjOld = static_cast<std::int32_t>(snodeSize) - 1; jjOld >= 0; --jjOld) {
        const std::int32_t jjNew = snodeOldToNew[jjOld];

        newFrontSize[jjNew] += ef.mFrontSize[jjOld];

        if (seen[jjNew])
            continue;

        const std::int32_t kkOld = ef.mParent[jjOld];
        if (kkOld != NIL) {
            const std::int32_t kkNew = snodeOldToNew[kkOld];
            newParent[jjNew] = kkNew;
        }
        newUpdateSize[jjNew] = ef.mUpdateSize[jjOld];
        seen[jjNew] = true;
    }

    ef.mSnodeSize = newSnodeSize;
    ef.mNodeToSnode.swap(newNodeToSnode);
    ef.mParent.swap(newParent);
    ef.mFrontSize.swap(newFrontSize);
    ef.mUpdateSize.swap(newUpdateSize);

    if (storedAZero)
        ef.mExactPatterns = false;   // the merged columns now differ, by the zeros we bought

    // The child, sibling and root links still describe the old forest. The caller rebuilds
    // them with finalizeLinks.
}


template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
