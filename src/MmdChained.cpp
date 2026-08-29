#include "oblio/MmdChained.h"

#include "oblio/QuotientGraphChained.h"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

// MmdChained.cpp - MmdFlat on a CHAINED CLIQUE STORE, everything else identical: every clique lives
// in the dead segment of the pivot that formed it, threaded through the runs it came from, so it
// costs no storage of its own.
//
// WHAT IT IS FOR, AND IT IS TWO THINGS.
//
// FIRST, IT IS THE ALIGNMENT VEHICLE FOR A DIFFERENTIAL. Comparing our ordering against an oracle
// is only clean when the two hold their cliques the same way, or every difference is confounded
// with layout. What then remains is either LAYOUT or an IMPROVEMENT to carry back into our own
// ladder.
//
// SECOND, IT IS A PREDICTABLE-SPACE VERSION OF MmdFlat. Given a machine you know whether A fits,
// but you cannot know whether L fits, nnz(L) depending on the ordering being computed. So a method
// that stays within `O(n + m)` carries a guarantee no amount of speed substitutes for: IF THE INPUT
// FITS, THE ANSWER IS REACHABLE. Chaining is that guarantee bought deliberately, not frugality for
// its own sake.
//
// WHAT IT CANNOT ANSWER. It prices our arena against a CHAINED store and says nothing about a
// COMPACTED one, which is a different thing again: one workspace reused, with a clique taking over
// the slots of the variable that formed it. That comparison is AmdCompacted's and MmdCompacted's.
//
// ITS OBLIGATION IS TWOFOLD. It must return MmdFlat's permutation exactly, and it must stay
// ENCODING-IDENTICAL to MmdFlat: a fold that lands in the shared quotient graph lands here too, or
// the comparison silently stops being about storage.
//
// `Buckets` is the shared type, nothing about degree buckets depending on how cliques are stored,
// and the quotient graph is `QuotientGraphChained`, header-only and so compiled into this
// translation unit, which is what keeps a comparison apples to apples.

namespace Oblio {

namespace {





}  // anonymous namespace


std::vector<std::int32_t> orderMmdChained(const std::vector<std::size_t>&  colPtr,
                                          const std::vector<std::int32_t>& rowIdx,
                                          std::int32_t delta) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraphChained qg(colPtr, rowIdx);
    // The deepest of the four walk-order conventions: it fixes the order of C[pivot], hence the
    // content order of every list built from it, hence which of several equal-degree vertices is
    // picked. Same sets and same cost either way, different winner among equals.
    qg.setReverseIncidence(true);
    std::vector<std::int32_t> pivots;
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;

    // NO `degrees` ARRAY AND NO `outmatched` ARRAY. Both live in `Buckets::mPrev`; the only reader
    // that needed a `degrees` array is `unfile`, which recovers the bucket from the link itself. We
    // file at the DEGREE at every site, and the running minimum is maintained here, at the moment
    // of filing, rather than recomputed after each round.
    Buckets buckets(size);
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        const std::uint32_t degree = qg.adjacencySize(u);
        buckets.file(degree, u);
        minDegree = std::min(minDegree, degree);
    }

    // NO `prepassVertices` LIST: the prepass reads each successor before unfiling and needs none.
    std::vector<std::int32_t> batchPivots, refreshedCliqueMembers, twoSourceQueue, manySourceQueue;

    // NO DRIVER MARK ARRAY. The two levels this refresh needs, one surviving a whole clique and one
    // fresh per vertex, are two tags rather than two arrays, drawn from the graph's own counter
    // through advanceTag/mark/setMark. One counter is what makes it safe, since two tags drawn from
    // it can never be equal.

    // ---- the prepass ------------------------------------------------------------ Buckets 0 and 1
    // hold the isolated and the degree-1 vertices, and both eliminate without fill. Number them and
    // leave both empty; nothing is eliminated in the quotient-graph sense.
    //
    // ONE PASS AND NO LIST. Unfiling a vertex while walking its bucket destroys the link the walk
    // is standing on, so the successor is read BEFORE the unfile.
    for (std::uint32_t b = 0; b < 2 && b < size; ++b)
    for (std::int32_t u = buckets.head(b); u != NIL; ) {
        const std::int32_t next = buckets.next(u);   // before the unfile invalidates it
        buckets.unfile(u);
        qg.number(u);
        pivots.push_back(u);
        ++numEliminated;
        u = next;
    }
    if (size > 2) minDegree = 2;                // buckets 0 and 1 are empty now

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;

        // ---- one batch, no degree refreshed inside it ---------------------------
        std::uint32_t batchLimit = minDegree;
        if (delta > 0)
            batchLimit = std::min(minDegree + static_cast<std::uint32_t>(delta),
                                  static_cast<std::uint32_t>(size) - 1);

        batchPivots.clear();
        while (true) {
            if (buckets.empty(minDegree)) {
                if (minDegree >= batchLimit) break;
                ++minDegree;
                continue;
            }
            const std::int32_t pivot = buckets.head(minDegree);
            buckets.unfile(pivot);

            const std::vector<std::int32_t>& merged = qg.eliminateMmd(pivot);
            batchPivots.push_back(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
            for (std::int32_t u : merged) buckets.unfile(u);

            qg.forEachMember(pivot, [&](std::int32_t u) {
                buckets.unfile(u);                  // evict, and put a withheld
                buckets.restore(u);                 //   vertex back: one operation, two calls
            });

            if (numEliminated >= size) break;       // nothing left to update
            if (delta < 0) break;
        }

        // ---- one refresh, walked clique by clique ----------------------------- The LAST pivot of
        // a batch is the FIRST clique refreshed, and that order decides which of several equal-
        // degree vertices ends at a bucket head.
        for (auto rcit = batchPivots.rbegin(); rcit != batchPivots.rend(); ++rcit) {
            const std::int32_t rc = *rcit;
            // The dead are dropped HERE, before anything stamps: `mMark` carries GONE as well as
            // the tag, so stamping a dead member would overwrite GONE and bring it back to life.
            // This store never compacts a clique, so the mass eliminated are in the list too.
            refreshedCliqueMembers.clear();
            qg.forEachMember(rc, [&](std::int32_t v) {
                if (!qg.eliminated(v)) refreshedCliqueMembers.push_back(v);
            });

            const std::int32_t cliqueTag = qg.advanceTag();   // marked once for the clique
            for (std::int32_t u : refreshedCliqueMembers) qg.setMark(u, cliqueTag);
            std::uint32_t refreshedCliqueWeight = 0;
            for (std::int32_t u : refreshedCliqueMembers) refreshedCliqueWeight += qg.weight(u);

            // reach(u) has |A[u]| + |I[u]| sources once the new clique is counted, so one other
            // source means everything u reaches is in this clique plus that one place.
            // refreshedCliqueWeight already counts the clique and the other source is walked
            // directly, so no union is formed at all.
            twoSourceQueue.clear();
            manySourceQueue.clear();
            for (std::int32_t u : refreshedCliqueMembers) {
                if (buckets.filed(u) || buckets.outmatched(u)) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? twoSourceQueue : manySourceQueue).push_back(u);
            }

            // The two-source queue, walked in reverse so the head-pushed order is preserved.
            for (auto uit = twoSourceQueue.rbegin(); uit != twoSourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;
                // merged or withheld by an earlier two-source vertex
                if (qg.eliminated(u) || buckets.outmatched(u)) continue;
                const std::int32_t vertexTag = qg.advanceTag();
                // refreshedCliqueWeight is kept WHOLE and u's own weight subtracted at the end,
                // which is not the same as subtracting it now: the walk below can MERGE a vertex
                // into u, growing that weight, and the post-merge value is the correct one to
                // remove. Subtracting first files a supervariable one bucket too high per merged
                // vertex, so it is not picked as early as its size has earned.
                std::uint32_t closedReachableSetWeight = refreshedCliqueWeight;

                // NO LOOPS. A two-source vertex has exactly TWO sources and one of them is the
                // refreshed clique, by the test that put it on this list, so the other is unique
                // and can be indexed. Our lists are split rather than interleaved, so the case
                // analysis reads off the two lengths: `incidenceSize` is at least 1, the prune
                // having appended the pivot, so the two sources are either one variable and the
                // clique, or two cliques of which one is the clique.
                if (qg.adjacencySize(u) == 1) {                    // a variable
                    const std::int32_t v = qg.adjacencyMmd(u)[0];
                    // ONE LOAD FOR BOTH QUESTIONS. `vertexTag` is the newest tag drawn, so anything
                    // at or above it is either this pass's own stamp or GONE, and both mean skip.
                    const std::int32_t vMark = qg.mark(v);
                    if (vMark < vertexTag && vMark != cliqueTag) {   // not seen, dead, or counted
                        qg.setMark(v, vertexTag);
                        closedReachableSetWeight += qg.weight(v);
                    }
                } else {                                           // two cliques; take the other
                    const std::int32_t* uIncidence = qg.incidenceMmd(u);
                    const std::int32_t  oc =
                        (uIncidence[0] == rc) ? uIncidence[1] : uIncidence[0];
                    qg.forEachMember(oc, [&](std::int32_t v) {
                        const std::int32_t vMark = qg.mark(v);
                        if (v == u || vMark >= vertexTag) return;  // seen this pass, or dead
                        if (vMark == cliqueTag) {
                            // v is in the new clique and in this same other source, so it sees
                            // at least what u sees.
                            if (buckets.filed(v) || buckets.outmatched(v)) return;
                            if (qg.adjacencySize(v) + qg.incidenceSize(v) - 1 == 1) {
                                qg.merge(u, v);      // identical reach: u absorbs it
                                ++numEliminated;
                            } else {
                                buckets.outmatch(v);    // reaches more, so never minimal first
                            }
                            return;
                        }
                        qg.setMark(v, vertexTag);
                        closedReachableSetWeight += qg.weight(v);
                    });
                }

                const std::uint32_t degree = closedReachableSetWeight - qg.weight(u);
                buckets.file(degree, u);
                minDegree = std::min(minDegree, degree);
            }

            // The many-source queue, the same stack.
            for (auto uit = manySourceQueue.rbegin(); uit != manySourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;                 // the full union
                if (qg.eliminated(u) || buckets.outmatched(u)) continue;
                const std::uint32_t degree = qg.reachableSetWeight(u); // reach excludes u already
                buckets.file(degree, u);
                minDegree = std::min(minDegree, degree);
            }
        }

    }

    gPeakCliqueMembers = qg.numPeakCliqueMembers();
    return qg.orderAscending(pivots);
}


}  // namespace Oblio
