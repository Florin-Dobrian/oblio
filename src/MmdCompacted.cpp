#include "oblio/MmdCompacted.h"

#include "oblio/QuotientGraph.h"          // Buckets
#include "oblio/QuotientGraphCompacted.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

// MmdCompacted.cpp - MmdFlat on AMD_2'S CLIQUE LAYOUT: one pooled workspace with a free cursor and
// a compaction, where production keeps a separate append-only arena in elimination order.
//
// IT IS A CELL OF THE LAYOUT MATRIX, and the matrix is why it exists. Two algorithms and three
// clique layouts give six cells; with only the diagonal filled, every reading confounds the layout
// with the algorithm that happens to use it. B is a driver on its OWN branch's vendored layout and
// C is a driver on the OTHER branch's, so this is the mmd counterpart of AmdCompacted and pairs
// with it down a column. See docs/DESIGN_DECISIONS.md (2026-08-16, the layout matrix, and
// 2026-08-17).
//
// WHAT IT IS FOR, precisely. With encoding held equal, genmmd's dead-segment scheme COSTS MmdFlat
// about 8 percent on squares and 29 on cubes, and AMD_2's pool EARNS AmdFlat about 9 percent,
// rising to 15 at large n. Each of those is a single reading on a single algorithm, so neither can
// be told apart from something about how that algorithm walks. This file is the second reading of
// the pool: if mmd on it also comes in below 1.0, the pool wins on both algorithms and our arena is
// the thing to change.
//
// IT DOES, BY 5 TO 8 PER CENT. `MmdCompacted / MmdFlat` on alpamayo reads 0.92 to 0.95 on square
// grids from 801 to 1601 a side and about 0.91 at 81 cubed, so the pooled workspace beats our
// append-only arena on this branch as it does on the amd one.
//
// THE FIGURE WAS 0.90 UNTIL 2026-08-19 AND THAT WAS NOT A FAIR COMPARISON. This file held its own
// quotient graph then, in an anonymous namespace, so the class was inlined into the pivot loop
// while `MmdFlat`'s was called across a translation-unit boundary. Build both the same way, either
// way round, and the answer is 0.92 to 0.95. See docs/DESIGN_DECISIONS.md (2026-08-19).
//
// A PREVIOUS MmdCompacted HELD THIS NAME AND WAS SOMETHING ELSE. It was mmd on the PRODUCTION
// layout, a transitional vehicle for working the amd array folds out on the mmd side without
// disturbing the class six drivers run. It did that; the folds are in QuotientGraph and the file
// was replaced. Nothing of it survives here but the name.
//
// ITS OBLIGATION. It must return MmdFlat's permutation exactly, which is genmmd's. `make digest` in
// benchmarks/ordering hashes every driver's permutation over 73 grids and names which one moved;
// `make mmdorder` in experiments/ordering is what says correct. And it must stay ENCODING-IDENTICAL
// to MmdFlat: a fold that lands in QuotientGraph lands here too, or its time column stops being
// about storage and starts being about whatever drifted. That is the same obligation MmdChained
// carries, and MmdChained was found to have broken it on 2026-08-17.
//
// NOTHING HERE IS A COPY ANY MORE. This file held its own `QuotientGraphCompacted` and its own
// `Buckets`, generated from the production ones and then hand edited in the storage layer, until
// 2026-08-19. The buckets turned out to be `Oblio::Buckets` verbatim, so they are used directly;
// the graph became the shared `QuotientGraphCompacted`, which `AmdCompacted` also uses, with a
// suffix on each method the two vendored routines disagree about. This driver calls the `...Mmd`
// half of each pair. Design notes for both types live with them and are authoritative there.
//
// THE STORAGE IS COMPLETE AS OF 2026-08-19, in four steps, each verified alone and none of which
// moved a permutation: the elbow room off `nzaat` rather than off `nnz` with the diagonal, the
// clique-count recomputation, the walk converted to positions and cursors, and the mid-walk
// compactor with the absorbed-clique capture and the contraction report that go with it. The line
// this replaces said the storage had not landed and the file was a verbatim copy of MmdFlat, which
// is how it was first committed and stopped being true well before then.

namespace Oblio {

// HOW OFTEN THE POOL RAN OUT, read by tmp/ probes as gAmdCompactions is on the amd side. It is
// the one number that says whether `AMD_2`'s elbow room, the pattern plus a fifth plus n, is the
// right size for MMD's cliques: the two branches fill the pool at different rates, amd absorbing
// cliques aggressively where mmd does not, so the sizing is the part of this port that is not a
// copy. Non-zero and growing would mean the room is too small and the compactions are being paid
// for repeatedly.
//
// IT READS 1 AT EVERY SQUARE SIZE FROM 8 TO 401 A SIDE, which is `AmdCompacted`'s figure exactly,
// so the vendored headroom suits mmd's cliques as well as amd's. That comparison only became
// possible on 2026-08-19: the room here was computed from `nnz` with the diagonal and then rounded
// up to whatever capacity the allocator had handed back, so it ran about a fifth large and was not
// a stated figure at all.
//
// AND IT IS WHY THE COMPACTOR HAS TO BE TESTED BY SHRINKING THE POOL. At the shipped room it fires
// once per ordering, so the mid-walk path is nearly untested by any normal run; sized to exactly
// `nzaat` it fires 8 to 13 times over the same sizes and every permutation is byte identical.
std::size_t gMmdCompactions = 0;

namespace {

// `Buckets` IS PRODUCTION'S, NOT A COPY. This file held a `Buckets` that was
// include/oblio/QuotientGraph.h's class verbatim apart from the name. Nothing about it
// depends on how cliques are stored, being three link arrays over vertices, so the storage
// question this file exists to ask never reaches it. The copy went on 2026-08-19; only the
// quotient graph is a separate class.


struct TaggedScanCompacted {
    Buckets*                          buckets;
    std::vector<std::int32_t>&        w;             // per clique, Amd.cpp's tagged W
    const std::vector<std::uint32_t>& degree;
    std::vector<std::int32_t>&        touchedCliques;// the cliques this step reached, once each
    std::int32_t                      wflg;          // the tag for this elimination
    std::int32_t                      modulus;       // the driver's bucket count, n + 1
};

// THE QUOTIENT GRAPH IS THE SHARED `QuotientGraphCompacted`, and this file no longer
// carries a copy of it. It held one until 2026-08-19, taken into an anonymous namespace so
// that `AmdCompacted` could hold another; the two then evolved apart for months, which is exactly
// what a private copy is bad at. The class is now one body with a suffix on each method the
// two vendored routines disagree about, so this driver calls the `...Mmd` half of each pair
// and nothing else in it changed. See include/oblio/QuotientGraphCompacted.h.




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

    // WE FILE AT THE DEGREE, AT EVERY SITE, where genmmd files at the degree initially and at the
    // degree plus one on refresh. See MmdFlat.cpp for what that costs it and
    // private/MmdCorrected.cpp for the repaired reference this matches.
    Buckets buckets(size);
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        const std::uint32_t degree = qg.adjacencySize(u);
        buckets.file(degree, u);
        minDegree = std::min(minDegree, degree);
    }

    std::vector<std::int32_t> batch, cliqueMembers, twoSourceQueue, manySourceQueue;

    // The prepass. Buckets 0 and 1 hold the isolated and the degree-1 vertices, both of which
    // eliminate without fill; genmmd has one bucket to drain here, its `if(dg==0)dg=1` putting the
    // two together, and filing at the true degree separates them.
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

        batch.clear();
        while (true) {
            if (buckets.empty(minDegree)) {
                if (minDegree >= batchLimit) break;
                ++minDegree;
                continue;
            }
            const std::int32_t pivot = buckets.head(minDegree);
            buckets.unfile(pivot);

            const std::vector<std::int32_t>& merged = qg.eliminateMmd(pivot);
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
            for (std::int32_t u : merged) buckets.unfile(u);

            const std::int32_t* newClique     = qg.clique(pivot);
            const std::uint32_t newCliqueSize = qg.cliqueSize(pivot);
            for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
                const std::int32_t u = newClique[uk];
                buckets.unfile(u);                  // evict; mmdelm's bwd[rn] = 0 is both this
                buckets.restore(u);                 //   and putting a withheld vertex back
            }

            if (numEliminated >= size) break;       // genmmd: nothing left to update
            if (delta < 0) break;
        }

        for (auto cliqueIt = batch.rbegin(); cliqueIt != batch.rend(); ++cliqueIt) {
            const std::int32_t clique = *cliqueIt;
            const std::int32_t* refreshedClique     = qg.clique(clique);
            const std::uint32_t refreshedCliqueSize = qg.cliqueSize(clique);

            cliqueMembers.clear();
            for (std::uint32_t uk = 0; uk < refreshedCliqueSize; ++uk)
                if (!qg.eliminatedMmd(refreshedClique[uk]))
                    cliqueMembers.push_back(refreshedClique[uk]);

            const std::int32_t cliqueTag = qg.advanceTag();   // marked once for the clique
            for (std::int32_t u : cliqueMembers) qg.setMark(u, cliqueTag);
            std::uint32_t refreshedCliqueWeight = 0;
            for (std::int32_t u : cliqueMembers) refreshedCliqueWeight += qg.weight(u);

            twoSourceQueue.clear();
            manySourceQueue.clear();
            for (std::int32_t u : cliqueMembers) {
                if (buckets.filed(u) || buckets.outmatched(u)) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? twoSourceQueue : manySourceQueue).push_back(u);
            }

            for (auto uit = twoSourceQueue.rbegin(); uit != twoSourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;
                // merged or withheld by an earlier two-source vertex
                if (qg.eliminatedMmd(u) || buckets.outmatched(u)) continue;
                const std::int32_t vertexTag = qg.advanceTag();
                std::uint32_t degree = refreshedCliqueWeight;

                const std::int32_t* uAdjacency = qg.adjacencyMmd(u);
                for (std::uint32_t vk = 0; vk < qg.adjacencySize(u); ++vk) {
                    const std::int32_t v = uAdjacency[vk];
                    const std::int32_t vMark = qg.mark(v);
                    if (vMark >= vertexTag) continue;              // seen this pass, or dead
                    if (vMark == cliqueTag) continue;   // already in refreshedCliqueWeight
                    qg.setMark(v, vertexTag);
                    degree += qg.weight(v);
                }
                const std::int32_t* uIncidence = qg.incidenceMmd(u);
                for (std::uint32_t ck = 0; ck < qg.incidenceSize(u); ++ck) {
                    const std::int32_t c = uIncidence[ck];
                    if (c == clique) continue;
                    const std::int32_t* otherClique     = qg.clique(c);
                    const std::uint32_t otherCliqueSize = qg.cliqueSize(c);
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
                        degree += qg.weight(v);
                    }
                }

                const std::uint32_t filed = degree - qg.weight(u);
                buckets.file(filed, u);
                minDegree = std::min(minDegree, filed);
            }

            for (auto uit = manySourceQueue.rbegin(); uit != manySourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;                 // the full union, as md5 computes it
                if (qg.eliminatedMmd(u) || buckets.outmatched(u)) continue;
                const std::uint32_t degree = qg.reachableSetWeight(u); // reach excludes u already
                const std::uint32_t filed = degree;
                buckets.file(filed, u);
                minDegree = std::min(minDegree, filed);
            }
        }

    }

    // NOT AN ASSERT THAT THE LIVE COUNT IS ZERO. A clique dies when one of its members becomes a
    // pivot, and at the close of a run the last cliques can have had every member mass eliminated
    // into the pivot instead, leaving no one to absorb them, so a few entries legitimately survive.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");
    gMmdCompactions  = qg.numCompactions();
    gPeakCliqueMembers = qg.numPeakCliqueMembers();   // see include/oblio/QuotientGraph.h
    return qg.orderAscending(pivots);   // genmmd's mmdnum. See the ledger, entry 6.
}

} // namespace Oblio
