#include "oblio/MmdCompacted.h"

#include "oblio/QuotientGraph.h"          // Buckets
#include "oblio/QuotientGraphCompacted.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

// MmdCompacted.cpp - MmdFlat on a COMPACTED CLIQUE STORE: one pooled workspace with a free cursor
// and a compaction, where the flat driver keeps a separate append-only arena in elimination order.
//
// IT IS A CELL OF THE LAYOUT MATRIX. Two algorithms and three clique layouts give six cells, and
// with only the diagonal filled every reading confounds the layout with the algorithm that happens
// to use it. This is the mmd counterpart of AmdCompacted and pairs with it down a column.
//
// ITS OBLIGATION IS TWOFOLD. It must return MmdFlat's permutation exactly, and it must stay
// ENCODING-IDENTICAL to MmdFlat: a fold that lands in the shared quotient graph lands here too, or
// its time column stops being about storage and starts being about whatever drifted.
//
// `Buckets` is the shared type, nothing about it depending on how cliques are stored, and the
// quotient graph is the shared `QuotientGraphCompacted` with a suffix on each method the two
// branches disagree about. This driver calls the `...Mmd` half of each pair.

namespace Oblio {

// HOW OFTEN THE POOL RAN OUT. It is the one number that says whether the elbow room is the right
// size for this branch's cliques: the two branches fill the pool at different rates, amd absorbing
// cliques aggressively where mmd does not, so the sizing is the part of this port that is not a
// copy. Non-zero and growing would mean the room is too small and the compactions are being paid
// for repeatedly.
//
// AND IT IS WHY THE COMPACTOR HAS TO BE TESTED BY SHRINKING THE POOL. At the shipped room it fires
// about once per ordering, so the mid-walk path is nearly untested by any ordinary run.
std::size_t gMmdCompactions = 0;

namespace {



struct TaggedScanCompacted {
    Buckets*                          buckets;
    std::vector<std::int32_t>&        w;             // per clique, the tagged workspace
    const std::vector<std::uint32_t>& degree;
    std::vector<std::int32_t>&        touchedCliques;// the cliques this step reached, once each
    std::int32_t                      wflg;          // the tag for this elimination
    std::int32_t                      modulus;       // the driver's bucket count, n + 1
};





} // namespace

std::vector<std::int32_t> orderMmdCompacted(const std::vector<std::size_t>&  colPtr,
                                            const std::vector<std::int32_t>& rowIdx,
                                            std::int32_t delta) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraphCompacted qg(colPtr, rowIdx);
    qg.setReverseIncidence(true);
    qg.enableMarks();   // this branch needs the tag array; the amd one does not
    std::vector<std::int32_t> pivots;
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;

    // WE FILE AT THE DEGREE, AT EVERY SITE. Filing a refreshed vertex one higher would penalise it
    // against one no pivot has reached yet, and the minimum selected would not always be the
    // minimum.
    Buckets buckets(size);
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        const std::uint32_t degree = qg.adjacencySize(u);
        buckets.file(degree, u);
        minDegree = std::min(minDegree, degree);
    }

    std::vector<std::int32_t> batchPivots, refreshedCliqueMembers, twoSourceQueue, manySourceQueue;

    // The prepass. Buckets 0 and 1 hold the isolated and the degree-1 vertices, both of which
    // eliminate without fill.
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

            const std::int32_t* newClique     = qg.clique(pivot);
            const std::uint32_t newCliqueSize = qg.cliqueSize(pivot);
            for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
                const std::int32_t u = newClique[uk];
                buckets.unfile(u);                  // evict, and put a withheld
                buckets.restore(u);                 //   vertex back: one operation, two calls
            }

            if (numEliminated >= size) break;       // nothing left to update
            if (delta < 0) break;
        }

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

            twoSourceQueue.clear();
            manySourceQueue.clear();
            for (std::int32_t u : refreshedCliqueMembers) {
                if (buckets.filed(u) || buckets.outmatched(u)) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? twoSourceQueue : manySourceQueue).push_back(u);
            }

            for (auto uit = twoSourceQueue.rbegin(); uit != twoSourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;
                // merged or withheld by an earlier two-source vertex
                if (qg.eliminatedMmd(u) || buckets.outmatched(u)) continue;
                const std::int32_t vertexTag = qg.advanceTag();
                // Bounded by n: a sum of weights over DISJOINT sets, which is what the guards
                // below make them. It reads as a fixed value here and accumulates further down.
                std::uint32_t closedReachableSetWeight = refreshedCliqueWeight;

                const std::int32_t* uAdjacency = qg.adjacencyMmd(u);
                for (std::uint32_t vk = 0; vk < qg.adjacencySize(u); ++vk) {
                    const std::int32_t v = uAdjacency[vk];
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

            for (auto uit = manySourceQueue.rbegin(); uit != manySourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;                 // the full union, as md5 computes it
                if (qg.eliminatedMmd(u) || buckets.outmatched(u)) continue;
                const std::uint32_t degree = qg.reachableSetWeight(u); // reach excludes u already
                buckets.file(degree, u);
                minDegree = std::min(minDegree, degree);
            }
        }

    }

    // NOT AN ASSERT THAT THE LIVE COUNT IS ZERO. A clique dies when one of its members becomes a
    // pivot, and at the close of a run the last cliques can have had every member mass eliminated
    // into the pivot instead, leaving no one to absorb them, so a few entries legitimately survive.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");
    gMmdCompactions  = qg.numCompactions();
    gPeakCliqueMembers = qg.numPeakCliqueMembers();
    return qg.orderAscending(pivots);
}

} // namespace Oblio
