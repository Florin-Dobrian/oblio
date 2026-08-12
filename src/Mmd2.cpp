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

    // The filed value, not the degree. mmdint files a degree-0 vertex under 1, `if(dg==0)dg=1`,
    // and the refresh below files under degree + 1, so from here on this array holds what MMD
    // compares and files by rather than the degree itself.
    std::vector<std::uint32_t> degrees(size);
    Buckets buckets(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        degrees[u] = std::max<std::uint32_t>(qg.adjacencySize(u), 1);
        buckets.file(degrees[u], u);
    }
    std::uint32_t minDegree = *std::min_element(degrees.begin(), degrees.end());

    // Withheld from the buckets rather than refiled, which is genmmd's bwd[nd] = -maxint. Not the
    // same as merged: the vertex is still live and still reachable, it simply cannot be the
    // minimum before the vertex that outmatched it. An elimination that reaches it puts it back.
    // A byte per vertex, not std::vector<bool>; see the note on Buckets::mFiled. Constructed once
    // per ordering, which is where that specialization costs.
    std::vector<std::uint8_t> outmatched(size, 0);

    std::vector<std::int32_t> touchedRound(size, NIL);
    std::vector<std::int32_t> touched, batch, elementMembers, q2h, qxh, refreshed;
    std::int32_t numRounds = 0;

    // The driver's own membership scratch. Two levels, as mmdupd has: one tag survives a whole
    // element and says "already counted in dg0", the other is fresh per vertex, so one q2h vertex
    // cannot hide a neighbor from the next.
    std::vector<std::int32_t> mark(size, NIL);
    std::int32_t              tag = 0;

    // ---- the prepass ------------------------------------------------------------
    // Bucket 1 holds the isolated and the degree-1 vertices together, by the convention above.
    // Number them and leave the bucket empty. Nothing is eliminated in the quotient-graph sense.
    for (std::int32_t u = buckets.head(1); u != NIL; u = buckets.next(u)) touched.push_back(u);
    for (std::int32_t u : touched) {
        buckets.unfile(degrees[u], u);
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
        touched.clear();
        while (true) {
            if (buckets.empty(minDegree)) {
                if (minDegree >= batchLimit) break;
                ++minDegree;
                continue;
            }
            const std::int32_t pivot = buckets.head(minDegree);
            buckets.unfile(degrees[pivot], pivot);

            const std::vector<std::int32_t>& merged = qg.eliminate(pivot);
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
            for (std::int32_t u : merged) {
                buckets.unfile(degrees[u], u);
                degrees[u] = 0;
            }
            degrees[pivot] = 0;

            const std::int32_t* clique     = qg.clique(pivot);
            const std::uint32_t cliqueSize = qg.cliqueSize(pivot);
            for (std::uint32_t k = 0; k < cliqueSize; ++k) {
                const std::int32_t u = clique[k];
                outmatched[u] = 0;                  // mmdelm's bwd[rn] = 0, back in the running
                buckets.unfile(degrees[u], u);      // evict, with a stale degree
                if (touchedRound[u] != numRounds) {
                    touchedRound[u] = numRounds;
                    touched.push_back(u);
                }
            }

            if (numEliminated >= size) break;       // genmmd: nothing left to update
            if (delta < 0) break;
        }

        // ---- one refresh, walked element by element -----------------------------
        refreshed.clear();
        for (std::int32_t element : batch) {
            const std::int32_t* members     = qg.clique(element);
            const std::uint32_t membersSize = qg.cliqueSize(element);

            elementMembers.clear();
            for (std::uint32_t k = 0; k < membersSize; ++k)
                if (!qg.eliminated(members[k])) elementMembers.push_back(members[k]);

            ++tag;
            const std::int32_t elementTag = tag;    // dg0's members, marked once for the element
            for (std::int32_t v : elementMembers) mark[v] = elementTag;
            std::uint32_t dg0 = 0;
            for (std::int32_t v : elementMembers) dg0 += qg.weight(v);

            // reach(u) has |A[u]| + |I[u]| sources once the new element is counted, so one other
            // source means everything u reaches is in this element plus that one place. dg0
            // already counts the element, and the other source is walked directly, so no union is
            // formed at all.
            q2h.clear();
            qxh.clear();
            for (std::int32_t u : elementMembers) {
                if (buckets.filed(u) || outmatched[u] != 0) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? q2h : qxh).push_back(u);
            }

            for (std::int32_t u : q2h) {
                if (qg.eliminated(u) || outmatched[u] != 0) continue;   // by an earlier q2h vertex
                ++tag;
                const std::int32_t vertexTag = tag;
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
                    if (qg.eliminated(v) || mark[v] == vertexTag) continue;
                    if (mark[v] == elementTag) continue;           // already counted in dg0
                    mark[v] = vertexTag;
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
                        if (v == u || qg.eliminated(v) || mark[v] == vertexTag) continue;
                        if (mark[v] == elementTag) {
                            // v is in the new element and in this same other source, so it sees
                            // at least what u sees.
                            if (buckets.filed(v) || outmatched[v] != 0) continue;
                            if (qg.adjacencySize(v) + qg.incidenceSize(v) - 1 == 1) {
                                qg.merge(u, v);      // identical reach: u absorbs it
                                ++numEliminated;
                            } else {
                                outmatched[v] = 1;      // reaches more, so never minimal first
                            }
                            continue;
                        }
                        mark[v] = vertexTag;
                        degree += qg.weight(v);
                    }
                }

                degrees[u] = std::max<std::uint32_t>(degree - qg.weight(u) + 1, 1);
                buckets.file(degrees[u], u);
                refreshed.push_back(u);
            }

            for (std::int32_t u : qxh) {                 // the full union, as md5 computes it
                if (qg.eliminated(u) || outmatched[u] != 0) continue;
                const std::uint32_t degree = qg.reachableWeight(u); // reach excludes u already
                degrees[u] = std::max<std::uint32_t>(degree + 1, 1);
                buckets.file(degrees[u], u);
                refreshed.push_back(u);
            }
        }

        for (std::int32_t u : refreshed) minDegree = std::min(minDegree, degrees[u]);
        ++numRounds;
    }

    return qg.order(pivots);
}

} // namespace Oblio
