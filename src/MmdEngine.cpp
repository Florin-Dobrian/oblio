#include "oblio/MmdEngine.h"

#include "oblio/Buckets.h"
#include "oblio/ElmOrder.h"
#include "oblio/QuotientGraphFlat.h"
#include "oblio/QuotientGraphCompacted.h"
#include "oblio/Types.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

// MmdEngine.cpp - the body of the mmd ordering, written once and instantiated for each clique
// store. Both instantiations and both graph classes are in this one unit, so each is compiled with
// its store's bodies visible and inlined into the pivot loop.
//
// THE STORES DIFFER IN WHAT THEY CAN BE ASKED. Only a store with a bounded pool can run out and be
// compacted, so only that one answers `numCompactions`. The overload pair below is where the body
// stops having to know which store it holds: the template answers zero for a store that does not
// publish the figure, and the compacted store's own overload answers for it.

namespace Oblio {
namespace {

template<class QuotientGraph> std::size_t numCompactionsOf(const QuotientGraph&) { return 0; }

std::size_t numCompactionsOf(const QuotientGraphCompacted& qg) { return qg.numCompactions(); }

} // namespace


template<class QuotientGraph>
void MmdEngine<QuotientGraph>::compute(const std::vector<std::size_t>&  colPtr,
                               const std::vector<std::int32_t>& rowIdx,
                               ElmOrder& eo) const {
    eo = ElmOrder();
    if (colPtr.empty()) return;
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return;

    QuotientGraph qg(colPtr, rowIdx);
    // The deepest of the four walk-order conventions: it fixes the order of C[pivot], hence the
    // content order of every list built from it, hence which of several equal-degree vertices is
    // picked. Same sets and same cost either way, different winner among equals.
    qg.setReverseIncidence(true);
    qg.enableMarks();   // this branch needs the tag array; the amd one does not
    std::vector<std::int32_t> pivots;
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;

    // WE FILE AT THE DEGREE, AT EVERY SITE. Filing a refreshed vertex one higher would penalise it
    // against one no pivot has reached yet, and the minimum selected would not always be the
    // minimum.
    //
    // NO `degrees` ARRAY: the only reader that needed one is `unfile`, which recovers the bucket
    // from the link. The running minimum is maintained here, at the moment of filing, rather than
    // recomputed after each round.
    Buckets buckets(size);
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        const std::uint32_t degree = qg.adjacencySize(u);
        buckets.file(degree, u);
        minDegree = std::min(minDegree, degree);
    }

    // NO `outmatched` ARRAY either: withholding is a value in the same link, `Buckets::outmatch`.

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
        if (mDelta > 0)
            batchLimit = std::min(minDegree + static_cast<std::uint32_t>(mDelta),
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

            const std::int32_t* newClique     = qg.clique(pivot);
            const std::uint32_t newCliqueSize = qg.cliqueSize(pivot);
            for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
                const std::int32_t u = newClique[uk];
                buckets.unfile(u);                  // evict, and put a withheld
                buckets.restore(u);                 //   vertex back: one operation, two calls
            }

            if (numEliminated >= size) break;       // nothing left to update
            if (mDelta < 0) break;
        }

        // ---- one refresh, walked clique by clique ----------------------------- The LAST pivot of
        // a batch is the FIRST clique refreshed, and that order decides which of several equal-
        // degree vertices ends at a bucket head.
        for (auto rcit = batchPivots.rbegin(); rcit != batchPivots.rend(); ++rcit) {
            const std::int32_t rc = *rcit;
            const std::int32_t* refreshedClique     = qg.clique(rc);
            const std::uint32_t refreshedCliqueSize = qg.cliqueSize(rc);

            // The dead are dropped HERE, before anything stamps: `mMark` carries GONE as well as
            // the tag, so stamping a dead member would overwrite GONE and bring it back to life.
            refreshedCliqueMembers.clear();
            for (std::uint32_t uk = 0; uk < refreshedCliqueSize; ++uk)
                if (!qg.eliminatedMmd(refreshedClique[uk]))
                    refreshedCliqueMembers.push_back(refreshedClique[uk]);

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
                if (qg.eliminatedMmd(u) || buckets.outmatched(u)) continue;
                const std::int32_t vertexTag = qg.advanceTag();
                // refreshedCliqueWeight is kept WHOLE and u's own weight subtracted at the end,
                // which is not the same as subtracting it now: the walk below can MERGE a vertex
                // into u, growing that weight, and the post-merge value is the correct one to
                // remove. Subtracting first files a supervariable one bucket too high per merged
                // vertex, so it is not picked as early as its size has earned.
                // Bounded by n: a sum of weights over DISJOINT sets, which is what the guards
                // below make them. It reads as a fixed value here and accumulates further down.
                std::uint32_t closedReachableSetWeight = refreshedCliqueWeight;

                // Not hoisted, deliberately. A two-source vertex has adjacencySize + incidenceSize
                // == 2 by the test that put it on this list, so these two loops run over at most
                // two cliques between them and a length loaded up front is overhead rather than a
                // saving.
                const std::int32_t* uAdjacency = qg.adjacencyMmd(u);
                for (std::uint32_t vk = 0; vk < qg.adjacencySize(u); ++vk) {
                    const std::int32_t v = uAdjacency[vk];
                    // ONE LOAD FOR BOTH QUESTIONS. `vertexTag` is the newest tag drawn, so anything
                    // at or above it is either this pass's own stamp or GONE, and both mean skip.
                    const std::int32_t vMark = qg.mark(v);
                    if (vMark >= vertexTag) continue;              // seen this pass, or dead
                    if (vMark == cliqueTag) continue;   // already in refreshedCliqueWeight
                    qg.setMark(v, vertexTag);
                    closedReachableSetWeight += qg.weight(v);
                }
                const std::int32_t* uIncidence = qg.incidenceMmd(u);
                for (std::uint32_t ock = 0; ock < qg.incidenceSize(u); ++ock) {
                    const std::int32_t oc = uIncidence[ock];
                    if (oc == rc) continue;
                    const std::int32_t* otherClique     = qg.clique(oc);
                    const std::uint32_t otherCliqueSize = qg.cliqueSize(oc);
                    for (std::uint32_t vk = 0; vk < otherCliqueSize; ++vk) {
                        const std::int32_t v = otherClique[vk];
                        const std::int32_t vMark = qg.mark(v);
                        if (v == u || vMark >= vertexTag) continue;    // seen this pass, or dead
                        if (vMark == cliqueTag) {
                            // v is in the new clique and in this same other source, so it sees
                            // at least what u sees.
                            if (buckets.filed(v) || buckets.outmatched(v)) continue;
                            if (qg.adjacencySize(v) + qg.incidenceSize(v) - 1 == 1) {
                                qg.merge(u, v);      // identical reach: u absorbs it
                                ++numEliminated;
                            } else {
                                buckets.outmatch(v);    // reaches more, so never minimal first
                            }
                            continue;
                        }
                        qg.setMark(v, vertexTag);
                        closedReachableSetWeight += qg.weight(v);
                    }
                }

                const std::uint32_t degree = closedReachableSetWeight - qg.weight(u);
                buckets.file(degree, u);
                minDegree = std::min(minDegree, degree);
            }

            // The many-source queue, the same stack.
            for (auto uit = manySourceQueue.rbegin(); uit != manySourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;                 // the full union
                if (qg.eliminatedMmd(u) || buckets.outmatched(u)) continue;
                const std::uint32_t degree = qg.reachableSetWeight(u); // reach excludes u already
                buckets.file(degree, u);
                minDegree = std::min(minDegree, degree);
            }
        }

    }

    // THE COUNTER CROSS-CHECKED AGAINST A RECOMPUTATION, which the driver can do exactly because it
    // holds the pivot list and a clique's owner is a pivot. NOT AN ASSERT THAT IT IS ZERO: at the
    // close of a run the last cliques can have had every member MASS ELIMINATED into the pivot
    // instead, leaving no one to absorb them, so a few entries legitimately survive.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");
    eo.mOrder                = qg.orderAscending(pivots);
    eo.mNumPeakCliqueMembers = qg.numPeakCliqueMembers();
    eo.mNumBornCliqueMembers = qg.numBornCliqueMembers();
    eo.mNumCompactions       = numCompactionsOf(qg);
}


template class MmdEngine<QuotientGraphFlat>;
template class MmdEngine<QuotientGraphCompacted>;

} // namespace Oblio
