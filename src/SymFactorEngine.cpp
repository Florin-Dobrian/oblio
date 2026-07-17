#include "oblio/SymFactorEngine.h"

#include <algorithm>

namespace Oblio {

// Adapter: the factor's pattern needs only A's pattern, so the matrix overload pulls it out
// and forwards. The implementation below is free of Val and compiled once.
template<class Val>
bool SymFactorEngine::compute(const SparseMatrix<Val>& A, const Permutation& p,
                            const ElmForest& ef, SymFactor& sf) const {
    return compute(A.colPtr(), A.rowIdx(), p, ef, sf);
}

bool SymFactorEngine::compute(const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            const Permutation& p, const ElmForest& ef, SymFactor& sf) const {
    if (colPtr.empty())
        return false;
    const std::size_t size = colPtr.size() - 1;
    if (p.size() != size || ef.size() != size)
        return false;

    const std::size_t snodeSize = ef.snodeSize();

    // Copy the forest attributes the factor needs, so it stands alone afterwards.
    sf.mSize      = size;
    sf.mSnodeSize   = snodeSize;
    sf.mNumTrees  = ef.numTrees();
    sf.mHeight    = ef.height();
    sf.mFirstRoot = ef.firstRoot();
    sf.mLastRoot  = ef.lastRoot();

    sf.mNodeToSnode = ef.nodeToSnode();
    sf.mParent      = ef.parent();
    sf.mFirstChild  = ef.firstChild();
    sf.mNextSibling = ef.nextSibling();
    sf.mFrontSize   = ef.frontSize();
    sf.mUpdateSize  = ef.updateSize();

    // Offsets into the flat index array: each supernode holds its front indices
    // followed by its update indices, and both counts are known from the forest.
    sf.mSnodePtr.resize(snodeSize + 1);
    sf.mSnodePtr[0] = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk)
        sf.mSnodePtr[kk + 1] = sf.mSnodePtr[kk] + sf.mFrontSize[kk] + sf.mUpdateSize[kk];
    sf.mNumNodeIdx = sf.mSnodePtr[snodeSize];
    sf.mNodeIdx.resize(sf.mNumNodeIdx);

    // The front columns each supernode must read from A. When the supernodes have exactly
    // matching patterns the lowest front column of each is enough (Section 4.6 of the notes);
    // otherwise every front column must be read. Only one of the two runs.
    const bool exact = ef.exactPatterns();

    std::vector<std::int32_t> firstFrontNodeIdx;   // exact: one column per supernode
    std::vector<std::size_t>  snodeFrontPtr;        // otherwise: all of them, with offsets
    std::vector<std::int32_t> frontNodeIdx;
    if (exact)
        getFirstFrontalIndex(sf, firstFrontNodeIdx);
    else
        getFrontalIndices(sf, snodeFrontPtr, frontNodeIdx);

    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // Naming, following 0.9. Supernodes are doubled letters: jj a child, kk its parent, and
    // jj < kk since a child's label is below its parent's. A single letter is a column, so
    // the doubling is the column-to-supernode map applied, and that reads literally here:
    // lk is a front column of kk, which is to say nodeToSnode[lk] == kk. An l prefix is the
    // factor (permuted) ordering and an a prefix the ordering of A, so ak is lk's column in
    // A, and li and ai are the corresponding rows.
    //
    // A position is named for the pointer array it walks: its initials, and nothing else.
    //
    //   cp    colPtr         one entry per column      (0.9 calls all of these p)
    //   sp    snodePtr       one entry per supernode
    //   sfp   snodeFrontPtr  one entry per supernode, the front indices alone
    //   np    nodePtr        one entry per node        (the transpose, in sortIndices)
    //
    // sp1 reads the child's index set, sp2 writes the parent's, keeping 0.9's 1-is-child,
    // 2-is-parent pairing. The scratch pair mirrors SymFactor's own, with front marking the
    // part it holds: snodeFrontPtr : frontNodeIdx as snodePtr : nodeIdx.
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

    // One of the two contributors to a supernode's index set: the indices that come from A.
    // Given a front column lk of kk, take the rows of A's column at or below lk, map them into
    // the factor's ordering, and add the ones not already there.
    //
    // Every index it adds is an A index, permuted. Nothing else here invents one. The other
    // contributor is the child loop below, and that is where fill enters: A's pattern is a
    // subset of L's, so the originals arrive here and the rest arrives from the children.
    //
    // This is also the only place A is read, which is what lets the two regimes differ in *how
    // many times it is called* rather than in what they do.
    auto addIndicesFromA = [&](std::int32_t lk, std::int32_t kk, std::size_t& sp2) {
        const std::int32_t ak = newToOld[lk];

        for (std::size_t cp = colPtr[ak]; cp < colPtr[ak + 1]; ++cp) {
            const std::int32_t ai = rowIdx[cp];
            const std::int32_t li = oldToNew[ai];

            if (li < lk)
                continue;   // above the front column, so not in the index set of kk
            if (mark[li] == kk)
                continue;   // already in the index set of kk

            mark[li] = kk;
            sf.mNodeIdx[sp2++] = li;
        }
    };

    // For every supernode kk, in increasing order (a topological order, so the children of kk
    // are complete when kk is reached).
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        std::size_t sp2 = sf.mSnodePtr[kk];   // write cursor into the index set of kk

        // The contribution of kk itself. One front column when the patterns match exactly, all
        // of them otherwise.
        if (exact) {
            addIndicesFromA(firstFrontNodeIdx[kk], kk, sp2);
        } else {
            for (std::size_t sfp = snodeFrontPtr[kk]; sfp < snodeFrontPtr[kk + 1]; ++sfp)
                addIndicesFromA(frontNodeIdx[sfp], kk, sp2);
        }

        // The contribution of the children of kk.
        for (std::int32_t jj = sf.mFirstChild[kk]; jj != NIL; jj = sf.mNextSibling[jj]) {
            for (std::size_t sp1 = sf.mSnodePtr[jj]; sp1 < sf.mSnodePtr[jj + 1]; ++sp1) {
                const std::int32_t li = sf.mNodeIdx[sp1];

                if (sf.mNodeToSnode[li] == jj)
                    continue;   // a front index of jj, so it dies with jj: only the update
                                // indices of jj carry into the index set of kk
                if (mark[li] == kk)
                    continue;   // already in the index set of kk

                mark[li] = kk;
                sf.mNodeIdx[sp2++] = li;
            }
        }
    }

    sortIndices(sf);

    return true;
}

void SymFactorEngine::getFrontalIndices(const SymFactor& sf,
        std::vector<std::size_t>&  snodeFrontPtr,
        std::vector<std::int32_t>& frontNodeIdx) const {
    // The offsets are the front sizes, accumulated.
    snodeFrontPtr.resize(sf.mSnodeSize + 1);
    snodeFrontPtr[0] = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.mSnodeSize); ++kk)
        snodeFrontPtr[kk + 1] = snodeFrontPtr[kk] + sf.mFrontSize[kk];

    // Scatter each column into its supernode's slot. Columns are walked in increasing order, so
    // each supernode's front indices come out sorted.
    //
    // sfp[kk] is the next free position in frontNodeIdx for supernode kk, so it is an sfp in the
    // usual sense, one per supernode rather than one in hand. It starts at the offsets and
    // advances as the columns land.
    frontNodeIdx.resize(sf.mSize);
    std::vector<std::size_t> sfp(snodeFrontPtr.begin(), snodeFrontPtr.end() - 1);
    for (std::int32_t lk = 0; lk < static_cast<std::int32_t>(sf.mSize); ++lk) {
        const std::int32_t kk = sf.mNodeToSnode[lk];   // the supernode of column lk
        frontNodeIdx[sfp[kk]++] = lk;
    }
}

void SymFactorEngine::getFirstFrontalIndex(const SymFactor& sf,
        std::vector<std::int32_t>& firstFrontNodeIdx) const {
    // The lowest column of each supernode. Columns are walked in increasing order, so the first
    // one seen for a supernode is its lowest; later ones are skipped. NIL doubles as "not seen
    // yet", and since every supernode owns at least one column, none survives the pass.
    //
    // No offsets: each supernode has exactly one entry, so the position is the supernode.
    firstFrontNodeIdx.assign(sf.mSnodeSize, NIL);
    for (std::int32_t lk = 0; lk < static_cast<std::int32_t>(sf.mSize); ++lk) {
        const std::int32_t kk = sf.mNodeToSnode[lk];   // the supernode of column lk
        if (firstFrontNodeIdx[kk] == NIL)
            firstFrontNodeIdx[kk] = lk;
    }
}

void SymFactorEngine::sortIndices(SymFactor& sf) const {
    // Sort each supernode's index set into increasing order, by transposing twice.
    //
    // The factor is a supernode-to-rows structure. Its transpose is a row-to-supernodes
    // structure, and the names mirror it exactly:
    //
    //   factor:     snodePtr -> nodeIdx     supernode kk, and the nodes it holds
    //   transpose:  nodePtr -> snodeIdx     node li, and the supernodes that hold it
    //
    // The positions mirror too, sp into the factor and np into the transpose, and so do the two
    // passes: each walks one structure with a scalar position and fills the other through a
    // cursor array, and the second pass is the first with the roles exchanged. Reading the
    // transpose back with rows in increasing order leaves every supernode's set sorted, with no
    // comparison anywhere. Two counting sorts, in opposite directions.

    // The transpose's offsets: how many supernodes each row appears in, accumulated.
    //
    // The + 1 is a trick, not an off-by-one, and it is worth naming properly.
    //
    // A prefix sum comes in two flavours. The *inclusive* one over [a, b, c] gives
    // [a, a+b, a+b+c]: each output counts itself. The *exclusive* one gives [0, a, a+b]: each
    // output counts only what comes before it, and the first is zero.
    //
    // Offsets are exclusive by nature. nodePtr[li] must be the total for all rows strictly below
    // li, which is where row li's run begins; that is what makes nodePtr[0] == 0 and makes
    // nodePtr[li + 1] - nodePtr[li] the length of row li's run.
    //
    // But the loop below, nodePtr[li + 1] += nodePtr[li] sweeping upward, is mechanically an
    // *inclusive* sum: each cell absorbs the running total beneath it, its own count included.
    // Storing each count one slot high is what reconciles the two. Row li's count lives at
    // nodePtr[li + 1], so it is never in nodePtr[li] to be absorbed, and an inclusive sweep lands
    // on exclusive results. nodePtr[0] is never written at all, and stays zero for free.
    //
    //   init         [0, 0, 0, 0]
    //   count        [0, 1, 2, 1]     counts, each at li + 1
    //   prefix sum   [0, 1, 3, 4]     offsets: row 0 -> [0,1), row 1 -> [1,3), row 2 -> [3,4)
    //
    // The shift converts an inclusive sweep into an exclusive result. Without it the counts would
    // sit at nodePtr[li], and getting exclusive offsets would need a genuinely exclusive sum: a
    // temporary to hold each count before it is overwritten, or a second backward pass.
    //
    // OrderEngine builds colPtrOff the same way. getFrontalIndices below does *not* need the
    // shift: its counts already exist, in mFrontSize, so it can write a plain exclusive sum
    // directly (snodeFrontPtr[kk + 1] = snodeFrontPtr[kk] + mFrontSize[kk]) with nothing to
    // overwrite. The shift is for when the counts must be tallied in place.
    std::vector<std::size_t> nodePtr(sf.mSize + 1, 0);
    for (std::size_t sp = 0; sp < sf.mNumNodeIdx; ++sp)
        ++nodePtr[sf.mNodeIdx[sp] + 1];
    for (std::int32_t li = 0; li < static_cast<std::int32_t>(sf.mSize); ++li)
        nodePtr[li + 1] += nodePtr[li];

    // Pass one: walk the factor, fill the transpose. np[li] is row li's next free position in it.
    std::vector<std::int32_t> snodeIdx(sf.mNumNodeIdx);
    {
        std::vector<std::size_t> np(nodePtr.begin(), nodePtr.end() - 1);
        for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.mSnodeSize); ++kk)
            for (std::size_t sp = sf.mSnodePtr[kk]; sp < sf.mSnodePtr[kk + 1]; ++sp) {
                const std::int32_t li = sf.mNodeIdx[sp];
                snodeIdx[np[li]++] = kk;
            }
    }

    // Pass two: walk the transpose, fill the factor. sp[kk] is supernode kk's next free position
    // in it. Rows are walked in increasing order, so each supernode's set comes out sorted.
    {
        std::vector<std::size_t> sp(sf.mSnodePtr.begin(), sf.mSnodePtr.end() - 1);
        for (std::int32_t li = 0; li < static_cast<std::int32_t>(sf.mSize); ++li)
            for (std::size_t np = nodePtr[li]; np < nodePtr[li + 1]; ++np) {
                const std::int32_t kk = snodeIdx[np];
                sf.mNodeIdx[sp[kk]++] = li;
            }
    }
}

template bool SymFactorEngine::compute(const SparseMatrix<double>&, const Permutation&, const ElmForest&, SymFactor&) const;
template bool SymFactorEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, const ElmForest&, SymFactor&) const;

} // namespace Oblio
