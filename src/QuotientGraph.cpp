#include "oblio/QuotientGraph.h"

#include <utility>

namespace Oblio {

QuotientGraph::QuotientGraph(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mSourcePtr(colPtr.empty() ? 1 : colPtr.size()),
      mAdjacencySize(colPtr.empty() ? 0 : colPtr.size() - 1),
      mIncidenceSize(mAdjacencySize.size(), 0),
      mCliquePtr(mAdjacencySize.size(), 0),
      mCliqueSize(mAdjacencySize.size(), 0),
      mEliminated(mAdjacencySize.size(), 0),
      mMark(mAdjacencySize.size(), NIL) {
    const std::int32_t size = static_cast<std::int32_t>(mAdjacencySize.size());

    // One pass to place each column's run, dropping the diagonal. The runs are laid out in column
    // order and never move afterwards, so the offsets are written once here and only the lengths
    // change. The run is u's whole allowance: A[u] fills it now and I[u] grows into the room A[u]
    // gives up, never past mSourcePtr[u + 1].
    mSource.reserve(colPtr.empty() ? 0 : colPtr.back());
    mSourcePtr[0] = 0;
    for (std::int32_t aj = 0; aj < size; ++aj) {
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSource.push_back(rowIdx[cp]);
        mSourcePtr[aj + 1]  = mSource.size();
        mAdjacencySize[aj]  = mSourcePtr[aj + 1] - mSourcePtr[aj];
    }

    // Every vertex begins as a supervariable of one: a chain holding itself, so next is NIL and
    // last is the vertex. Mass elimination splices these together and order() walks them.
    mSuperNext.assign(size, NIL);
    mSuperLast.resize(size);
    mWeight.assign(size, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

void QuotientGraph::reachableSet(std::int32_t u, std::vector<std::int32_t>& reached) {
    // reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}
    //
    // The mark array is the set: mMark[v] == mTag is the membership test, one comparison, and
    // mMark[v] = mTag is the insertion, one store. So the union costs one pass per source rather
    // than a hash per member, and the result comes out in the order the sources were walked.
    //
    // The buffer is the caller's so that a caller in a loop can keep one, which is the whole of
    // what the returning overload costs.
    // Eliminated vertices are skipped rather than purged. Mass elimination cannot leave one
    // behind, since a vertex it merges belongs to the pivot's clique alone and is removed from
    // it; a live merge can, since the vertex it folds away is left where it lies at weight zero
    // and every clique that named it still does. Amd.cpp reaches the same place from the other
    // side, testing `nvi = Nv[i]` before using an entry.
    const bool live = mLiveMerges;         // hoisted: see the member's note
    ++mTag;
    reached.clear();
    mMark[u] = mTag;                       // never its own neighbor
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per element. Measured at 300 ms of AMD1's
    // 6.31 s on alpamayo, in the one accessor that reads as though it were free.
    const std::int32_t* source        = mSource.data() + mSourcePtr[u];
    const std::size_t   adjacencySize = mAdjacencySize[u];
    const std::size_t   incidenceSize = mIncidenceSize[u];
    for (std::size_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v = source[k];
        if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) {
            mMark[v] = mTag;
            reached.push_back(v);
        }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::size_t i = 0; i < incidenceSize; ++i) {
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mCliqueArena.data() + mCliquePtr[c];
        const std::size_t   membersSize = mCliqueSize[c];
        for (std::size_t k = 0; k < membersSize; ++k) {
            const std::int32_t v = members[k];
            if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) {
                mMark[v] = mTag;
                reached.push_back(v);
            }
        }
    }
}

std::vector<std::int32_t> QuotientGraph::reachableSet(std::int32_t u) {
    std::vector<std::int32_t> reached;
    reachableSet(u, reached);
    return reached;
}

std::size_t QuotientGraph::reachableSize(std::int32_t u) {
    // The same two passes as reachableSet, counting rather than collecting.
    const bool live = mLiveMerges;
    ++mTag;
    std::size_t reached = 0;
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per element. Measured at 300 ms of AMD1's
    // 6.31 s on alpamayo, in the one accessor that reads as though it were free.
    const std::int32_t* source        = mSource.data() + mSourcePtr[u];
    const std::size_t   adjacencySize = mAdjacencySize[u];
    const std::size_t   incidenceSize = mIncidenceSize[u];
    for (std::size_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v = source[k];
        if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) { mMark[v] = mTag; ++reached; }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::size_t i = 0; i < incidenceSize; ++i) {
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mCliqueArena.data() + mCliquePtr[c];
        const std::size_t   membersSize = mCliqueSize[c];
        for (std::size_t k = 0; k < membersSize; ++k) {
            const std::int32_t v = members[k];
            if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) { mMark[v] = mTag; ++reached; }
        }
    }
    return reached;
}

std::size_t QuotientGraph::reachableWeight(std::int32_t u) {
    const bool live = mLiveMerges;
    ++mTag;
    std::size_t reached = 0;
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per element. Measured at 300 ms of AMD1's
    // 6.31 s on alpamayo, in the one accessor that reads as though it were free.
    const std::int32_t* source        = mSource.data() + mSourcePtr[u];
    const std::size_t   adjacencySize = mAdjacencySize[u];
    const std::size_t   incidenceSize = mIncidenceSize[u];
    for (std::size_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v = source[k];
        if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) { mMark[v] = mTag; reached += mWeight[v]; }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::size_t i = 0; i < incidenceSize; ++i) {
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mCliqueArena.data() + mCliquePtr[c];
        const std::size_t   membersSize = mCliqueSize[c];
        for (std::size_t k = 0; k < membersSize; ++k) {
            const std::int32_t v = members[k];
            if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) { mMark[v] = mTag; reached += mWeight[v]; }
        }
    }
    return reached;
}

void QuotientGraph::number(std::int32_t u) {
    mEliminated[u] = 1;
    mLiveMerges    = true;   // a numbered vertex lingers in lists, so the walks must filter
}

void QuotientGraph::beginElimination(std::int32_t pivot,
                                     std::int32_t& inClique, std::int32_t& absorbed) {
    // Fill the scratch, then copy it into the clique. Swapping instead would save the copy and
    // cost more than it saves: the scratch would come back empty every time and grow again from
    // nothing at the next pivot, several reallocations apiece, where keeping it lets its capacity
    // settle after a few pivots and never grow again. The clique's own allocation is then a single
    // right-sized one rather than a sequence of doublings. Measured both ways.
    reachableSet(pivot, mReached);
    const std::vector<std::int32_t>& reached = mReached;

    // The absorbed cliques are I[pivot], read where they lie. Nothing below writes the pivot's run
    // (the prune rewrites the runs of C[pivot]'s members, and the pivot is not one of them), so no
    // copy and no scratch is needed to keep them alive across the passes that follow.
    const std::int32_t* absorbedCliques = mSource.data() + mSourcePtr[pivot] + mAdjacencySize[pivot];
    const std::size_t   absorbedSize    = mIncidenceSize[pivot];
    for (std::size_t i = 0; i < absorbedSize; ++i)
        mCliqueSize[absorbedCliques[i]] = 0;    // dead, its block left behind

    // The pivot becomes a clique holding its own reach: one block at the end of the arena, its
    // size known now and never larger again.
    mCliquePtr[pivot]  = mCliqueArena.size();
    mCliqueSize[pivot] = reached.size();
    mCliqueArena.insert(mCliqueArena.end(), reached.begin(), reached.end());

    // Stamp the new clique once and the absorbed cliques once, each with its own tag. Both are
    // then queried for free by the two compaction passes below. Clique ids and vertex ids share
    // one space, so one mark array serves both, the tags keeping them apart.
    ++mTag;
    inClique = mTag;
    for (std::int32_t v : reached) mMark[v] = inClique;
    ++mTag;
    absorbed = mTag;
    for (std::size_t i = 0; i < absorbedSize; ++i) mMark[absorbedCliques[i]] = absorbed;
}

const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot) {
    std::int32_t inClique = NIL;
    std::int32_t absorbed = NIL;
    beginElimination(pivot, inClique, absorbed);
    const std::vector<std::int32_t>& reached = mReached;

    // Both lists are compacted in place rather than rebuilt into a scratch and swapped. Every
    // pass here only ever removes, so the survivors can be written over the entries already read,
    // the write cursor trailing the read one, and nothing is allocated at all. Rebuilding into a
    // shared scratch and swapping had each list inherit some other vertex's buffer, which then
    // had to grow again; that idiom is right where a pass can add, and this one cannot.
    // Both lists live in u's one run and are rewritten front to back, the adjacency first and the
    // incidence into whatever the adjacency has just given up. The pivot is appended rather than
    // pushed: the run never has to grow, because a source is destroyed here for each one created,
    // which is the conservation argument on the members. The two cursors are what makes the second
    // pass safe as well as the first, since the incidence write starts where the compacted
    // adjacency ends and its read starts where the original adjacency ended, which is never lower.
    for (std::int32_t u : reached) {
        std::int32_t*     source        = mSource.data() + mSourcePtr[u];
        const std::size_t adjacencySize = mAdjacencySize[u];
        std::size_t       kept          = 0;
        for (std::size_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (v == pivot) continue;                  // no longer a variable
            if (mMark[v] == inClique) continue;        // both ends inside the new clique
            if (mLiveMerges && mEliminated[v] != 0) continue;   // numbered by a prepass, gone for good
            source[kept++] = v;
        }
        mAdjacencySize[u] = kept;                      // A[u] - C[pivot] - {pivot}

        const std::size_t incidenceSize = mIncidenceSize[u];
        std::size_t       write         = kept;
        for (std::size_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = source[adjacencySize + i];
            if (mMark[c] != absorbed) source[write++] = c;
        }
        source[write++]   = pivot;                     // u joins the new clique, id = pivot
        mIncidenceSize[u] = write - kept;
    }

    return finishElimination(pivot);
}

// The same prune, with the driver's first scan folded into the two loops. The header carries the
// argument for why this is one call and why it cannot be folded further; the loops below are the
// plain ones with an accumulation added on each survivor, and nothing else differs.
const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot,
                                                         ApproximateScan& scan) {
    std::int32_t inClique = NIL;
    std::int32_t absorbed = NIL;
    beginElimination(pivot, inClique, absorbed);
    const std::vector<std::int32_t>& reached = mReached;

    for (std::int32_t u : reached) {
        const std::size_t weightU       = mWeight[u];
        std::int32_t*     source        = mSource.data() + mSourcePtr[u];
        const std::size_t adjacencySize = mAdjacencySize[u];
        std::size_t       kept          = 0;
        std::size_t       explicitPart  = 0;
        for (std::size_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (v == pivot) continue;
            if (mMark[v] == inClique) continue;
            if (mLiveMerges && mEliminated[v] != 0) continue;
            source[kept++] = v;
            explicitPart += mWeight[v];                // the bound's explicit term, in this visit
        }
        mAdjacencySize[u]       = kept;
        scan.explicitPart[u]    = explicitPart;

        const std::size_t incidenceSize = mIncidenceSize[u];
        std::size_t       write         = kept;
        for (std::size_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = source[adjacencySize + i];
            if (mMark[c] == absorbed) continue;
            source[write++] = c;
            if (c == pivot) continue;                  // the new clique subtracts from nothing
            if (scan.mark[c] != scan.tag) {            // first sighting: start from |C[c]|
                scan.mark[c] = scan.tag;
                scan.touchedCliques.push_back(c);
                scan.outside[c] = scan.cliqueDegree[c] - weightU;
            } else {
                scan.outside[c] -= weightU;
            }
        }
        source[write++]   = pivot;
        mIncidenceSize[u] = write - kept;
    }

    return finishElimination(pivot);
}

const std::vector<std::int32_t>& QuotientGraph::finishElimination(std::int32_t pivot) {
    const std::vector<std::int32_t>& reached = mReached;

    // Mass elimination. u is indistinguishable from the pivot when the two had the same closed
    // neighborhood before the step, equivalently when everything u can still reach now lies
    // inside the new clique, and eliminating it next then creates no fill at all. The test is a
    // cheap sufficient condition for that: nothing explicit left and no clique but the new one.
    // It is conservative, and deliberately so; the exact test costs a reachability query per
    // candidate. See the mass-elimination section of experiments/ordering/README.md.
    std::vector<std::int32_t>& merged = mMerged;   // scratch, kept for its capacity
    merged.clear();
    for (std::int32_t u : reached) {
        if (mAdjacencySize[u] == 0 && mIncidenceSize[u] == 1 &&
            mSource[mSourcePtr[u]] == pivot) {         // A[u] empty, so I[u] starts at the run
            mIncidenceSize[u] = 0;
            mEliminated[u]    = 1;
            merged.push_back(u);
        }
    }
    if (!merged.empty()) {                             // C[pivot] - merged, one compaction pass
        ++mTag;
        for (std::int32_t u : merged) mMark[u] = mTag;
        std::int32_t*     members     = mCliqueArena.data() + mCliquePtr[pivot];
        const std::size_t membersSize = mCliqueSize[pivot];
        std::size_t       kept        = 0;
        for (std::size_t k = 0; k < membersSize; ++k)
            if (mMark[members[k]] != mTag) members[kept++] = members[k];
        mCliqueSize[pivot] = kept;

        for (std::int32_t u : merged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }
    }

    mAdjacencySize[pivot] = 0;
    mIncidenceSize[pivot] = 0;
    mEliminated[pivot]    = 1;
    return merged;
}

void QuotientGraph::merge(std::int32_t u, std::int32_t v) {
    mLiveMerges = true;                            // the reachable set must filter from now on
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    mAdjacencySize[v] = 0;
    mIncidenceSize[v] = 0;
    mEliminated[v]    = 1;
}

void QuotientGraph::absorb(const std::vector<std::int32_t>& cliques,
                           const std::int32_t* vertices, std::size_t vertexCount) {
    if (cliques.empty()) return;

    ++mTag;
    for (std::int32_t c : cliques) { mCliqueSize[c] = 0; mMark[c] = mTag; }

    for (std::size_t k = 0; k < vertexCount; ++k) {   // I[u] - dead, compacted in place
        const std::int32_t u         = vertices[k];
        std::int32_t*      incidence = mSource.data() + mSourcePtr[u] + mAdjacencySize[u];
        const std::size_t  size      = mIncidenceSize[u];
        std::size_t        kept      = 0;
        for (std::size_t i = 0; i < size; ++i)
            if (mMark[incidence[i]] != mTag) incidence[kept++] = incidence[i];
        mIncidenceSize[u] = kept;
    }
}

std::vector<std::int32_t> QuotientGraph::order(const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order;
    order.reserve(size());
    for (std::int32_t pivot : pivots)
        for (std::int32_t u = pivot; u != NIL; u = mSuperNext[u]) order.push_back(u);
    return order;
}

} // namespace Oblio
