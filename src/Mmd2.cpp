#include "oblio/Mmd2.h"

#include <algorithm>
#include <cstdint>

namespace Oblio {

std::vector<std::int32_t> orderMmd2(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    std::int32_t delta) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);
    std::vector<std::int32_t> pivots;
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;

    // NO `degrees` ARRAY. The filed value was never the degree (mmdint files a degree-0 vertex
    // under 1 and the refresh files under degree plus one), and the only reader that needed it
    // kept was `unfile`, which now recovers the bucket from the link. See Buckets.
    //
    // The running minimum moves with it: it was recomputed after each round from a `refreshed`
    // list that existed only to be walked once for this, where genmmd does `if(dg<*mdeg)*mdeg=dg`
    // at the moment it files, `private/Mmd.cpp` line 164.
    Buckets buckets(size);
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        const std::uint32_t degree = std::max<std::uint32_t>(qg.adjacencySize(u), 1);
        buckets.file(degree, u);
        minDegree = std::min(minDegree, degree);
    }

    // NO `outmatched` ARRAY either: withholding is a value in the same link, Buckets::outmatch,
    // which is genmmd's bwd[nd] = -maxint.

    // `prepassVertices` is the prepass's own list and nothing else. A `touched` list with a
    // `touchedRound` stamp beside it used to collect every member of every C[pivot] in a round.
    // Mmd1 READS that list, its refresh walking the touched vertices; this layer refreshes element
    // by element and walks `batch`, so it was filled once per clique member per pivot and never
    // read again. Inherited from the layer below and made redundant by the pass that replaced it,
    // which is the shape `refile`, the evicted list and the inert ternary all had.
    std::vector<std::int32_t> prepassVertices, batch, elementMembers, q2h, qxh;

    // NO DRIVER MARK ARRAY. The two levels this refresh needs, one surviving a whole element and
    // one fresh per vertex, are two tags rather than two arrays, and they go into the graph's own
    // stamp array through advanceTag/mark/setMark. That is genmmd's `marker` exactly: `mmdelm`
    // stamps it at level `tag` and `mmdupd` at level `mt = tag + md0`, one array between them.
    // One counter is what makes it safe, since two tags drawn from it can never be equal.

    // ---- the prepass ------------------------------------------------------------
    // Bucket 1 holds the isolated and the degree-1 vertices together, by the convention above.
    // Number them and leave the bucket empty. Nothing is eliminated in the quotient-graph sense.
    for (std::int32_t u = buckets.head(1); u != NIL; u = buckets.next(u))
        prepassVertices.push_back(u);
    for (std::int32_t u : prepassVertices) {
        buckets.unfile(u);
        qg.number(u);
        pivots.push_back(u);
        ++numEliminated;
    }
    if (size > 2) minDegree = 2;                // head[1] is empty now, and mdeg starts at 2

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

            const std::vector<std::int32_t>& merged = qg.eliminate(pivot);
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
            for (std::int32_t u : merged) buckets.unfile(u);

            const std::int32_t* clique     = qg.clique(pivot);
            const std::uint32_t cliqueSize = qg.cliqueSize(pivot);
            for (std::uint32_t k = 0; k < cliqueSize; ++k) {
                const std::int32_t u = clique[k];
                buckets.unfile(u);                  // evict; mmdelm's bwd[rn] = 0 is both this
                buckets.restore(u);                 //   and putting a withheld vertex back
            }

            if (numEliminated >= size) break;       // genmmd: nothing left to update
            if (delta < 0) break;
        }

        // ---- one refresh, walked element by element -----------------------------
        for (std::int32_t element : batch) {
            const std::int32_t* members     = qg.clique(element);
            const std::uint32_t membersSize = qg.cliqueSize(element);

            elementMembers.clear();
            for (std::uint32_t k = 0; k < membersSize; ++k)
                if (!qg.eliminated(members[k])) elementMembers.push_back(members[k]);

            const std::int32_t elementTag = qg.advanceTag();   // marked once for the element
            for (std::int32_t v : elementMembers) qg.setMark(v, elementTag);
            std::uint32_t dg0 = 0;
            for (std::int32_t v : elementMembers) dg0 += qg.weight(v);

            // reach(u) has |A[u]| + |I[u]| sources once the new element is counted, so one other
            // source means everything u reaches is in this element plus that one place. dg0
            // already counts the element, and the other source is walked directly, so no union is
            // formed at all.
            q2h.clear();
            qxh.clear();
            for (std::int32_t u : elementMembers) {
                if (buckets.filed(u) || buckets.outmatched(u)) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? q2h : qxh).push_back(u);
            }

            for (std::int32_t u : q2h) {
                if (qg.eliminated(u) || buckets.outmatched(u)) continue;   // by an earlier q2h vertex
                const std::int32_t vertexTag = qg.advanceTag();
                // dg0 is kept WHOLE and u's own weight subtracted at the end, which is
                // genmmd's `dg - qsize[en] + 1` and NOT the same as subtracting it now. The
                // walk below can MERGE a vertex into u, and genmmd's merge does
                // `qsize[en] += qsize[nd]` in that same walk, so the weight it subtracts is the
                // one AFTER the merge. Subtracting first files a supervariable one bucket too
                // high per merged vertex, so it is never picked as early as its size has earned.
                // This was a DEFECT here until 2026-08-07, found by aligning Mmd3 against
                // genmmd; see experiments/ordering/mmd3.py, ledger entry 5.
                std::uint32_t degree = dg0;

                // Not hoisted, deliberately. A q2h vertex has adjacencySize + incidenceSize == 2
                // by the test that put it on this list, so these two loops run over at most two
                // elements between them and a length loaded up front is overhead rather than a
                // saving. Hoist where a loop is long; leave it where the loop is short or exits
                // early. Measured both ways.
                const std::int32_t* adjacency = qg.adjacency(u);
                for (std::uint32_t a = 0; a < qg.adjacencySize(u); ++a) {
                    const std::int32_t v = adjacency[a];
                    // ONE LOAD FOR BOTH QUESTIONS. `vertexTag` is the newest tag drawn, so
                    // anything at or above it is either this pass's own stamp or GONE, and both
                    // mean skip. This was `qg.eliminated(v) || mark[v] == vertexTag`, two arrays.
                    const std::int32_t m = qg.mark(v);
                    if (m >= vertexTag) continue;                  // seen this pass, or dead
                    if (m == elementTag) continue;                 // already counted in dg0
                    qg.setMark(v, vertexTag);
                    degree += qg.weight(v);
                }
                const std::int32_t* incidence = qg.incidence(u);
                for (std::uint32_t i = 0; i < qg.incidenceSize(u); ++i) {
                    const std::int32_t c = incidence[i];
                    if (c == element) continue;
                    const std::int32_t* other     = qg.clique(c);
                    const std::uint32_t otherSize = qg.cliqueSize(c);
                    for (std::uint32_t k = 0; k < otherSize; ++k) {
                        const std::int32_t v = other[k];
                        const std::int32_t m = qg.mark(v);
                        if (v == u || m >= vertexTag) continue;    // seen this pass, or dead
                        if (m == elementTag) {
                            // v is in the new element and in this same other source, so it sees
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
                        degree += qg.weight(v);
                    }
                }

                const std::uint32_t filed = std::max<std::uint32_t>(degree - qg.weight(u) + 1, 1);
                buckets.file(filed, u);
                minDegree = std::min(minDegree, filed);
            }

            for (std::int32_t u : qxh) {                 // the full union, as md5 computes it
                if (qg.eliminated(u) || buckets.outmatched(u)) continue;
                const std::uint32_t degree = qg.reachableWeight(u); // reach excludes u already
                const std::uint32_t filed = std::max<std::uint32_t>(degree + 1, 1);
                buckets.file(filed, u);
                minDegree = std::min(minDegree, filed);
            }
        }

    }

    return qg.order(pivots);
}

} // namespace Oblio
