#include "oblio/MmdFlat.h"

#include <cassert>
#include <algorithm>
#include <cstdint>

namespace Oblio {
namespace {

// THE BODY, WITH ONE OPTIONAL OUT-PARAMETER. The public forms are an overload pair rather than one
// function with another default argument; see MmdFlat.h for why the type must not change.
std::vector<std::int32_t> orderMmdFlatImpl(const std::vector<std::size_t>&  colPtr,
                                           const std::vector<std::int32_t>& rowIdx,
                                           std::int32_t delta,
                                           std::size_t* numBornCliqueMembers) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);
    // The fourth walk, and the deepest: it fixes the order of C[pivot], hence the content order of
    // every list built from it. The other three are below. All four mirror genmmd holding a linked
    // list pushed at the head and read from the head, `list[nb] = h; h = nb`, where we append to a
    // vector. Same sets, same cost, different winner among equals, and minimum degree is settled by
    // exactly that. See experiments/ordering/mmd3.py.
    qg.setReverseIncidence(true);
    qg.enableMarks();   // this branch needs the tag array; the amd one does not
    std::vector<std::int32_t> pivots;
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;

    // WE FILE AT THE DEGREE, AT EVERY SITE, which genmmd does not: `mmdint` files at the degree
    // and `mmdupd` at the degree PLUS ONE, so a refreshed vertex is penalised by exactly one
    // against one no pivot has reached yet and the minimum is not always the minimum. On a 4 by 4
    // grid that puts a degree 4 vertex in the same bucket as four degree 3 vertices and takes the
    // degree 4 one. `MmdCorrected` is genmmd with that repaired and is what this driver matches;
    // `MmdVendored` keeps the original and is reference only.
    //
    // NO `degrees` ARRAY. The only reader that needed one kept was `unfile`, which recovers the
    // bucket from the link. See Buckets.
    //
    // The running minimum moves with it: it was recomputed after each round from a `refreshed`
    // list that existed only to be walked once for this, where genmmd does `if(dg<*mdeg)*mdeg=dg`
    // at the moment it files, `private/Mmd.cpp` line 164.
    Buckets buckets(size);
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        const std::uint32_t degree = qg.adjacencySize(u);
        buckets.file(degree, u);
        minDegree = std::min(minDegree, degree);
    }

    // NO `outmatched` ARRAY either: withholding is a value in the same link, Buckets::outmatch,
    // which is genmmd's bwd[nd] = -maxint.

    // NO `prepassVertices` LIST. It held the degree-1 bucket so the prepass could walk it after
    // emptying it; the prepass now reads each successor before unfiling and needs no list at all.
    // A `touched` list with a `touchedRound` stamp beside it went the same way on 2026-08-15,
    // having been filled once per clique member per pivot and never read on this layer.
    std::vector<std::int32_t> batchPivots, refreshedCliqueMembers, twoSourceQueue, manySourceQueue;

    // NO DRIVER MARK ARRAY. The two levels this refresh needs, one surviving a whole clique and
    // one fresh per vertex, are two tags rather than two arrays, and they go into the graph's own
    // stamp array through advanceTag/mark/setMark. That is genmmd's `marker` exactly: `mmdelm`
    // stamps it at level `tag` and `mmdupd` at level `mt = tag + md0`, one array between them.
    // One counter is what makes it safe, since two tags drawn from it can never be equal.

    // ---- the prepass ------------------------------------------------------------
    // Buckets 0 and 1 hold the isolated and the degree-1 vertices, and both eliminate without
    // fill. Number them and leave both empty. Nothing is eliminated in the quotient-graph sense.
    // genmmd has ONE bucket to drain here, its `if(dg==0)dg=1` putting degree 0 and degree 1
    // together; filing at the true degree separates them and this loop takes both.
    //
    // ONE PASS AND NO LIST. This was a collect loop into `prepassVertices` followed by a walk of
    // it, and the list existed only because unfiling a vertex while walking the bucket destroys
    // the link the walk is standing on. Reading the successor BEFORE the unfile removes the need,
    // which is what genmmd does at `private/Mmd.cpp`:
    //
    //     while(nextmd>0){int mn=nextmd; nextmd=invp[mn]; marker[mn]=maxint; invp[mn]=-num; num++;}
    //
    // It matters most where it looks least worth doing. On a matrix with no off-diagonal entries
    // every vertex goes through here and nothing else runs, so this loop IS the ordering:
    // Boeing/bcsstm39, n = 46772 and nnz(A) = 46772, is one of the rows where MmdFlat reads worst
    // against genmmd. See benchmarks/matrices, `make ordering`.
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

        // ---- one refresh, walked clique by clique -----------------------------
        // THE ACCESSORS STAY. `qg.mark(v)` and `qg.weight(v)` are `mArray[v]`, a base load plus
        // the clique, and the obvious idea is to hoist both bases once per round and index them
        // directly, which is what genmmd gets for free by taking its arrays as parameters. Built
        // and measured 2026-08-15, BOTH WAYS. Cachegrind on a 100x100 grid: 371403 more
        // instructions and 161414 more data reads, the two extra live values costing more in the
        // register allocator than the reloads they remove. Then timed on alpamayo across seven
        // square grids and six cubic ones, with MMD1 and MMD2 as controls carrying the same class
        // and not the change: FLAT, three sizes down and three up, inside the noise everywhere.
        //
        // So it was reverted for costing something and buying nothing. The wider lesson is worth
        // more than the change: a two to three percent movement in instruction count is not
        // visible on this machine at all, in either direction, so it is not on its own a reason
        // to keep a change OR to drop one. Two other candidates were reverted on counters alone
        // that day, the stamping fold and the arena cursor, and neither was ever timed.
        // The driver's clique list, genmmd's `list[mn] = ehead; ehead = mn`, so the LAST pivot of
        // a batch is the FIRST clique refreshed.
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
            // refreshedCliqueWeight already counts the clique, and the other source is walked
            // directly, so no union is formed at all.
            twoSourceQueue.clear();
            manySourceQueue.clear();
            for (std::int32_t u : refreshedCliqueMembers) {
                if (buckets.filed(u) || buckets.outmatched(u)) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? twoSourceQueue : manySourceQueue).push_back(u);
            }

            // mmdupd's q2h list, `list[nb] = q2h; q2h = nb`.
            for (auto uit = twoSourceQueue.rbegin(); uit != twoSourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;
                // by an earlier two-source vertex
                if (qg.eliminatedMmd(u) || buckets.outmatched(u)) continue;
                const std::int32_t vertexTag = qg.advanceTag();
                // refreshedCliqueWeight is kept WHOLE and u's own weight subtracted at the end,
                // which is genmmd's `dg - qsize[en] + 1` and not the same as subtracting it now.
                // The walk below can MERGE a vertex into u, and genmmd's merge does `qsize[en] +=
                // qsize[nd]` in that same walk, so the weight it subtracts is the one AFTER the
                // merge. Subtracting first files a supervariable one bucket too high per merged
                // vertex, so it is not picked as early as its size has earned. See
                // experiments/ordering/mmd3.py, ledger entry 5.
                // Bounded by n: a sum of weights over DISJOINT sets, which is what the guards
                // below make them. It reads as a fixed value here and accumulates further down.
                std::uint32_t closedReachableSetWeight = refreshedCliqueWeight;

                // Not hoisted, deliberately. A two-source vertex has adjacencySize +
                // incidenceSize == 2 by the test that put it on this list, so these two loops run
                // over at most two
                // cliques between them and a length loaded up front is overhead rather than a
                // saving. Hoist where a loop is long; leave it where the loop is short or exits
                // early. Measured both ways.
                const std::int32_t* uAdjacency = qg.adjacencyMmd(u);
                for (std::uint32_t vk = 0; vk < qg.adjacencySize(u); ++vk) {
                    const std::int32_t v = uAdjacency[vk];
                    // ONE LOAD FOR BOTH QUESTIONS. `vertexTag` is the newest tag drawn, so
                    // anything at or above it is either this pass's own stamp or GONE, and both
                    // mean skip. This was `qg.eliminatedMmd(v) || mark[v] == vertexTag`, two
                    // arrays.
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

            // mmdupd's qxh list, the same stack.
            for (auto uit = manySourceQueue.rbegin(); uit != manySourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;                 // the full union, as md5 computes it
                if (qg.eliminatedMmd(u) || buckets.outmatched(u)) continue;
                const std::uint32_t degree = qg.reachableSetWeight(u); // reach excludes u already
                buckets.file(degree, u);
                minDegree = std::min(minDegree, degree);
            }
        }

    }

    // THE COUNTER CROSS-CHECKED AGAINST A RECOMPUTATION, which the driver can do exactly because
    // it holds the pivot list and a clique's owner is a pivot. Births and deaths are spread over
    // four call sites in the shared class and nothing else in the suite would notice if they
    // stopped balancing.
    //
    // NOT AN ASSERT THAT IT IS ZERO, which is what this said first and which is wrong even though
    // it passes here. A clique dies when a member of it becomes a pivot, and at the close of a run
    // the last cliques can have had every member MASS ELIMINATED into the pivot instead, leaving
    // no one to absorb them; the amd side leaves 1 to 3 entries on small grids. Holding on these
    // graphs would have made it a trap rather than a check.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");
    gPeakCliqueMembers = qg.numPeakCliqueMembers();   // see include/oblio/QuotientGraph.h
    if (numBornCliqueMembers != nullptr) *numBornCliqueMembers = qg.numBornCliqueMembers();
    return qg.orderAscending(pivots);   // genmmd's mmdnum. See the ledger, entry 6.
}

} // namespace

std::vector<std::int32_t> orderMmdFlat(const std::vector<std::size_t>&  colPtr,
                                       const std::vector<std::int32_t>& rowIdx,
                                       std::int32_t delta) {
    return orderMmdFlatImpl(colPtr, rowIdx, delta, nullptr);
}

// The same ordering, reporting every member ever put into a clique. `benchmarks/matrices` prints
// it beside nnz(L) as `cC`; see QuotientGraph::numBornCliqueMembers.
std::vector<std::int32_t> orderMmdFlat(const std::vector<std::size_t>&  colPtr,
                                       const std::vector<std::int32_t>& rowIdx,
                                       std::int32_t delta,
                                       std::size_t& numBornCliqueMembers) {
    return orderMmdFlatImpl(colPtr, rowIdx, delta, &numBornCliqueMembers);
}

} // namespace Oblio
