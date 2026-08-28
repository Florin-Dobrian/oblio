#include "oblio/MmdChained.h"

#include "oblio/QuotientGraphChained.h"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

// MmdChained.cpp - MmdFlat on GENMMD'S clique storage, everything else identical.
//
// WHAT IT IS FOR, AND IT IS TWO THINGS. Both are permanent; this file is not an experiment awaiting
// a verdict and its old stop condition is withdrawn.
//
// FIRST, IT IS THE ALIGNMENT VEHICLE FOR A DIFFERENTIAL. Comparing our ordering against genmmd is
// only clean when the two hold their cliques the same way; otherwise every difference is confounded
// with layout. This file removes that confound. What then remains between it and genmmd is either
// LAYOUT, whose price is now measured, or an IMPROVEMENT, which is carried back into our own
// ladder. That is not hypothetical: with storage held equal on the amd side, the same arrangement
// surfaced five array folds that had nothing to do with layout at all, and they are being ported to
// AmdFlat and, where applicable, to Amd1 and Amd2. The 2n mark shape came out of this file the same
// way. A verdict on storage therefore does NOT retire it, storage never having been the only thing
// it was for.
//
// SECOND, IT IS THE PREDICTABLE-SPACE VERSION OF MmdFlat, and that is worth having on its own. From
// a conversation with Alex Pothen: given a machine you know whether A fits, but you cannot know
// whether L fits, nnz(L) depending on the ordering being computed. So a method that stays within
// `O(n + m)` carries a guarantee no amount of speed substitutes for: IF THE INPUT FITS, THE ANSWER
// IS REACHABLE. genmmd's chaining and AMD_2's compaction are that guarantee bought
// deliberately, not frugality for its own sake. Our arena is the right default for a known shape on
// a known machine solved repeatedly; this is the right one when whether an answer exists is the
// open question. See docs/DESIGN_DECISIONS.md (2026-08-16).
//
// THE PRICE, measured: 4 to 10 percent slower than MmdFlat at every size on the square ladder, and
// 1.15 to 1.38x genmmd where MmdFlat reads 1.02 to 1.19x. On 16.61M instructions against 14.22M and
// 123510 D1 read misses against 119331. So the second arena buys speed and costs the guarantee.
//
// IT CARRIES A PRIVATE COPY OF QuotientGraph, named QuotientGraphChained, so the storage can be
// changed without touching the class six drivers share. THE COST IS EVERY SHARED FOLD LANDING
// TWICE, and it
// is accepted on purpose; `make digest` catches a copy that has stopped reproducing its original in
// half a second, which is what makes carrying it tolerable.
//
// WHAT THIS FILE CANNOT ANSWER, stated because it was briefly assumed to, 2026-08-16. It prices our
// arena against GENMMD's dead-segment scheme. It says nothing about `AMD_2`'s, which is a different
// thing again: one workspace that is compacted and reused, with a clique taking over the slots of
// the variable that formed it. The amd branch has no equivalent of this file, so the arena has
// never been compared against that storage at all. When the vendored AMD turned out to zigzag on a
// power-of-two scaling ladder while three other codes did not, storage was the first suspect and
// this file was cited in support; it does not support it. See docs/DESIGN_DECISIONS.md
// (2026-08-16, later).
//
// What it is NOT is a second implementation to be maintained for its own sake. Its obligation is
// to stay encoding-identical to MmdFlat: a fold that lands in QuotientGraph lands here too, or the
// comparison silently stops being about storage. See docs/DESIGN_DECISIONS.md (2026-08-15).
//
// THE COPY IS EXACT AS OF THIS COMMIT and the comments were not duplicated with it. Every design
// note for these types lives in include/oblio/QuotientGraph.h and is authoritative there; reading
// this file means reading it alongside this one. `Buckets` is that header's `Buckets` too, which is
// why promoting `QuotientGraphChained` needs no second bucket class. Only the DIFFERENCES are
// commented here, and at this commit there are none: MmdChained must return MmdFlat's permutation,
// which is genmmd's, so `make mmdorder` and `make test` are the acceptance check unchanged.
//
// The scheme itself: see the section "The vendored storage scheme, and what it is worth" in
// experiments/ordering/README.md. In one line, genmmd keeps every clique in the dead segment of
// the pivot that formed it, which places blocks in vertex-id order and costs no storage at all.
// That placement was once credited with 0.41 of the time ratio on 2D grids; THAT CLAIM IS
// WITHDRAWN, this file being what refuted it.

namespace Oblio {

namespace {

// `Buckets` IS PRODUCTION'S, NOT A COPY. This file held a `Buckets` that was
// include/oblio/QuotientGraph.h's class verbatim apart from the name, one method it lacked
// and one it added. Nothing about degree buckets depends on how cliques are stored, they are
// three link arrays over vertices, so the storage question this file exists to ask never
// reaches them. The copy went on 2026-08-19.
//
// THE ONE METHOD IT ADDED WAS `evict`, which was `unfile` followed by resetting an outmatched
// vertex to unfiled, so it is `unfile(u); restore(u);` in production's spelling and that is
// what src/MmdFlat.cpp already writes. Nothing had to be added to the shared class.


// What an approximate-degree driver accumulates while the eliminator is already walking the lists,
// handed to the second `eliminate` overload so that the walk serves both. The members are the
// driver's own arrays, held by reference: the graph fills them and owns none of them.
//
// `tag` is the only member that moves, and the driver sets it before each elimination, exactly as
// it would before its own scan. The rest are bound once and reused for the whole ordering.
// THE QUOTIENT GRAPH IS `QuotientGraphChained`, in include/oblio/QuotientGraphChained.h. It
// lived in this file, in an anonymous namespace, until 2026-08-19. Being header-only it is
// still compiled into this translation unit, which is the point: a driver has to be in the
// same unit as its graph or a comparison against genmmd is not apples to apples. What moved
// is where the source sits, not what the compiler sees.
//
// TWO DEAD DECLARATIONS WENT WITH IT. The class declared `eliminate` overloads taking an
// `ApproximateScanChained` and a `TaggedScanChained`, neither ever defined and neither ever
// called: they are the amd branch's fused scans, carried across by the copy and never
// reached by an mmd driver. An undefined declaration costs nothing until someone links a
// call to it, which is why nothing caught them.


}  // anonymous namespace


std::vector<std::int32_t> orderMmdChained(const std::vector<std::size_t>&  colPtr,
                                          const std::vector<std::int32_t>& rowIdx,
                                          std::int32_t delta) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraphChained qg(colPtr, rowIdx);
    // The fourth walk, and the deepest: it fixes the order of C[pivot], hence the content order of
    // every list built from it. The other three are below. All four mirror genmmd holding a linked
    // list pushed at the head and read from the head, `list[nb] = h; h = nb`, where we append to a
    // vector. Same sets, same cost, different winner among equals, and minimum degree is settled by
    // exactly that. See experiments/ordering/mmd3.py.
    qg.setReverseIncidence(true);
    std::vector<std::int32_t> pivots;
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;

    // NO `degrees` ARRAY AND NO `outmatched` ARRAY. Both live in Buckets::mPrev now, on
    // genmmd's `bwd` encoding; see the class. The only reader that needed a `degrees` array kept
    // was `unfile`, which recovers the bucket from the link itself. We file at the DEGREE at every
    // site, where genmmd files at the degree initially and at the degree plus one on refresh; see
    // MmdFlat.cpp and private/MmdCorrected.cpp.
    //
    // The running minimum moves with it. It was recomputed after each round from a `refreshed`
    // list, which existed only to be walked once for this; genmmd instead does `if(dg<*mdeg)
    // *mdeg=dg` at the moment it files, `private/Mmd.cpp` line 164, so the value is maintained
    // where it is produced and the list has nothing left to do.
    Buckets buckets(size);
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        const std::uint32_t degree = qg.adjacencySize(u);
        buckets.file(degree, u);
        minDegree = std::min(minDegree, degree);
    }

    // NO `prepassVertices` LIST. It held the degree-1 bucket so the prepass could walk it after
    // emptying it; the prepass now reads each successor before unfiling and needs no list at all.
    // A `touched` list with a `touchedRound` stamp beside it went the same way on 2026-08-15,
    // having been filled once per clique member per pivot and never read on this layer.
    std::vector<std::int32_t> batch, cliqueMembers, twoSourceQueue, manySourceQueue;

    // NO DRIVER MARK ARRAY. The two levels this refresh needs, one surviving a whole clique and
    // one fresh per vertex, are two tags rather than two arrays, and they go into the graph's own
    // stamp array through advanceTag/mark/setMark. That is genmmd's `marker` exactly: `mmdelm`
    // stamps it at level `tag` and `mmdupd` at level `mt = tag + md0`, one array between them.
    // One counter is what makes it safe, since two tags drawn from it can never be equal, so a
    // comparison against a captured tag means what it says and nothing else.

    // ---- the prepass ------------------------------------------------------------
    // Buckets 0 and 1 hold the isolated and the degree-1 vertices, and both eliminate without
    // fill. Number them and leave both empty. Nothing is eliminated in the quotient-graph sense.
    // genmmd has ONE bucket to drain here, its `if(dg==0)dg=1` putting degree 0 and degree 1
    // together; filing at the true degree separates them and this loop takes both.
    //
    // ONE PASS AND NO LIST. This was a compact loop into `prepassVertices` followed by a walk of
    // it, and the list existed only because unfiling a vertex while walking the bucket destroys
    // the link the walk is standing on. Reading the successor BEFORE the unfile removes the need,
    // which is what genmmd does at `private/Mmd.cpp`:
    //
    //     while(nextmd>0){int mn=nextmd; nextmd=invp[mn]; marker[mn]=maxint; invp[mn]=-num; num++;}
    //
    // It matters most where it looks least worth doing. On a matrix with no off-diagonal entries
    // every vertex goes through here and nothing else runs, so this loop IS the ordering. Measured
    // on MmdFlat: 8.6 percent of a pure-diagonal ordering at n = 46772, and 0.2 percent of a
    // 100x100 grid. See benchmarks/matrices, `make ordering`.
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

            qg.forEachMember(pivot, [&](std::int32_t u) {
                buckets.unfile(u);                  // evict; mmdelm's bwd[rn] = 0 is both this
                buckets.restore(u);                 //   and putting a withheld vertex back
            });

            if (numEliminated >= size) break;       // genmmd: nothing left to update
            if (delta < 0) break;
        }

        // ---- one refresh, walked clique by clique -----------------------------
        // The driver's clique list, genmmd's `list[mn] = ehead; ehead = mn`, so the LAST pivot of
        // a batch is the FIRST clique refreshed.
        for (auto cliqueIt = batch.rbegin(); cliqueIt != batch.rend(); ++cliqueIt) {
            const std::int32_t clique = *cliqueIt;
            cliqueMembers.clear();
            qg.forEachMember(clique, [&](std::int32_t v) {
                if (!qg.eliminated(v)) cliqueMembers.push_back(v);
            });

            const std::int32_t cliqueTag = qg.advanceTag();   // marked once for the clique
            for (std::int32_t u : cliqueMembers) qg.setMark(u, cliqueTag);
            std::uint32_t refreshedCliqueWeight = 0;
            for (std::int32_t u : cliqueMembers) refreshedCliqueWeight += qg.weight(u);

            // reach(u) has |A[u]| + |I[u]| sources once the new clique is counted, so one other
            // source means everything u reaches is in this clique plus that one place.
            // refreshedCliqueWeight already counts the clique, and the other source is walked
            // directly, so no union is formed at all.
            twoSourceQueue.clear();
            manySourceQueue.clear();
            for (std::int32_t u : cliqueMembers) {
                if (buckets.filed(u) || buckets.outmatched(u)) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? twoSourceQueue : manySourceQueue).push_back(u);
            }

            // mmdupd's q2h list, `list[nb] = q2h; q2h = nb`.
            for (auto uit = twoSourceQueue.rbegin(); uit != twoSourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;
                // merged or withheld by an earlier two-source vertex
                if (qg.eliminated(u) || buckets.outmatched(u)) continue;
                const std::int32_t vertexTag = qg.advanceTag();
                // refreshedCliqueWeight is kept WHOLE and u's own weight subtracted at the end,
                // which is genmmd's `dg - qsize[en] + 1` and not the same as subtracting it now.
                // The walk below can MERGE a vertex into u, and genmmd's merge does `qsize[en] +=
                // qsize[nd]` in that same walk, so the weight it subtracts is the one AFTER the
                // merge. Subtracting first files a supervariable one bucket too high per merged
                // vertex, so it is not picked as early as its size has earned. See
                // experiments/ordering/mmd3.py, ledger entry 5.
                std::uint32_t degree = refreshedCliqueWeight;

                // NO LOOPS. A two-source vertex has exactly TWO sources and one of them is the new
                // clique, by the test that put it on this list, so the other one is unique and
                // can be indexed. genmmd does the same and does it in three lines: it reads the
                // first entry of the segment, steps past it if that is the clique, and branches
                // once on whether what remains is a variable or a clique,
                //
                //     is=xadj[en]; nb=adjncy[is]; if(nb==el) nb=adjncy[is+1]; lk=nb;
                //     if(fwd[nb]>=0){ dg+=qsize[nb]; goto n2100; }
                //
                // What stood here was two loops, each re-reading its bound through an accessor on
                // every iteration and one of them testing `c == clique` per entry, to find a
                // single entry already known to be there. The comment defending that said the
                // loops were short enough not to be worth hoisting, which was true and beside the
                // point: the loops themselves are what genmmd does not have. This pass measured
                // 869993 instructions against its 341996.
                //
                // Our lists are split where genmmd's interleave, so the case analysis reads off
                // the two lengths rather than off a sign: `incidenceSize` is at least 1, the prune
                // having appended the pivot, so the two sources are either one variable and the
                // clique, or two cliques of which one is the clique.
                if (qg.adjacencySize(u) == 1) {                    // a variable, genmmd's fwd >= 0
                    const std::int32_t v = qg.adjacencyMmd(u)[0];
                    // ONE LOAD FOR BOTH QUESTIONS. `vertexTag` is the newest tag drawn, so
                    // anything at or above it is either this pass's own stamp or GONE, and both
                    // mean skip.
                    const std::int32_t vMark = qg.mark(v);
                    if (vMark < vertexTag && vMark != cliqueTag) {   // not seen, dead, or counted
                        qg.setMark(v, vertexTag);
                        degree += qg.weight(v);
                    }
                } else {                                           // two cliques; take the other
                    const std::int32_t* uIncidence = qg.incidenceMmd(u);
                    const std::int32_t  c =
                        (uIncidence[0] == clique) ? uIncidence[1] : uIncidence[0];
                    qg.forEachMember(c, [&](std::int32_t v) {
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
                        degree += qg.weight(v);
                    });
                }

                const std::uint32_t filed = degree - qg.weight(u);
                buckets.file(filed, u);
                minDegree = std::min(minDegree, filed);
            }

            // mmdupd's qxh list, the same stack.
            for (auto uit = manySourceQueue.rbegin(); uit != manySourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;                 // the full union, as md5 computes it
                if (qg.eliminated(u) || buckets.outmatched(u)) continue;
                const std::uint32_t degree = qg.reachableSetWeight(u); // reach excludes u already
                const std::uint32_t filed = degree;
                buckets.file(filed, u);
                minDegree = std::min(minDegree, filed);
            }
        }

    }

    gPeakCliqueMembers = qg.numPeakCliqueMembers();   // see include/oblio/QuotientGraph.h
    return qg.orderAscending(pivots);   // genmmd's mmdnum. See the ledger, entry 6.
}


}  // namespace Oblio
