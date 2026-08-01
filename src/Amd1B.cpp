#include "oblio/Amd1B.h"

#include "oblio/QuotientGraph.h"

#include <algorithm>
#include <utility>

namespace Oblio {

std::vector<std::int32_t> orderAmd1B(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);
    std::vector<std::int32_t> pivots;               // the order over supervariables
    pivots.reserve(size);
    std::size_t numEliminated = 0;                  // a counter, not a scan of the flags

    // The cache. Exact at construction, since no clique exists yet and the whole neighborhood is
    // still explicit, and a bound from the first elimination onward. That the cached value is
    // itself a bound is what makes the second cap below hold inductively: an upper bound on an
    // earlier degree is still an upper bound now.
    std::vector<std::size_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    Buckets buckets(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        buckets.file(degrees[u], u);
    std::size_t minDegree = *std::min_element(degrees.begin(), degrees.end());

    // The driver's own membership scratch, separate from the quotient graph's: which vertices lie
    // in the new clique, and which cliques the step has already listed. Both are sets built by
    // stamping and queried by comparison, never allocated.
    std::vector<std::int32_t> mark(size, NIL);
    std::int32_t              tag = 0;

    // |C[c] - C[p]| per clique, indexed by clique id and hoisted out of the loop. The prototype
    // allocates and zeroes it per pivot, which reads better and is O(n) per step, O(n^2) over the
    // run in bookkeeping alone, independent of the graph. Only the entries this step wrote are
    // touched, and they are exactly the ones it will read.
    std::vector<std::size_t>  outside(size, 0);
    std::vector<std::int32_t> touchedCliques;

    // |C[c]| per live clique, weighted, which is what the scan below subtracts from. Exact rather
    // than an estimate, and the invariants that keep it so are worth stating because they are not
    // obvious. A live clique never holds an eliminated vertex, since eliminating v absorbs every
    // clique in I[v]. Mass elimination only removes a vertex whose I[u] is exactly {pivot}, so no
    // other clique is touched. And the value is written once, when the clique is formed, from the
    // already-trimmed member list.
    std::vector<std::size_t> cliqueDegree(size, 0);

    // The bound's explicit term, accumulated by the eliminator while it prunes A[u] rather than by
    // a second walk afterwards. Written for every reached vertex and read for every survivor.
    std::vector<std::size_t> explicitPart(size, 0);

    // Bound once and reused: only `tag` moves, and the driver sets it before each elimination
    // exactly as it would before its own scan.
    ApproximateScan scan{explicitPart, outside, cliqueDegree, touchedCliques, mark, 0};

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        // The scan runs inside the elimination, so its stamp and its list are prepared before the
        // call rather than after it. Everything else about the step is unchanged.
        touchedCliques.clear();
        ++tag;
        scan.tag = tag;

        const std::vector<std::int32_t>& merged = qg.eliminate(pivot, scan);
        pivots.push_back(pivot);
        numEliminated += 1 + merged.size();

        buckets.unfile(degrees[pivot], pivot);      // unfile before zeroing: the bucket index is
        degrees[pivot] = 0;                         //   read from the degree
        for (std::int32_t u : merged) {
            buckets.unfile(degrees[u], u);
            degrees[u] = 0;
        }

        // ---- the bound, in place of an exact refresh -----------------------------------
        // Everything the new clique reached needs a new degree, and nothing else can have
        // changed. C[p] is therefore both the refresh set and the domain the decomposition is a
        // statement about, which is the coincidence the whole placement rests on.
        // Read at the moment of use: the arena holding it grows as cliques are formed, so a
        // pointer taken before the next elimination is the only one that is safe, and this is
        // that window.
        const std::int32_t* pivotClique     = qg.clique(pivot);
        const std::size_t   pivotCliqueSize = qg.cliqueSize(pivot);

        // No membership stamp for C[p] is needed any more: the subtraction below asks which
        // cliques the new clique's members belong to, never which vertices a clique's members
        // are, so the one query that used it is gone. The amd2 prototype still carries the stamp,
        // inherited from amd1 and dead there for the same reason.
        std::size_t degme = 0;                      // |C[p]|, weighted
        for (std::size_t k = 0; k < pivotCliqueSize; ++k) degme += qg.weight(pivotClique[k]);
        cliqueDegree[pivot] = degme;                // what the scan below subtracts from

        // The clique-degree scan is not here any more: it ran inside the elimination above, in the
        // same visit that pruned each incidence list. What it computed is unchanged, and the
        // reasoning for why the fusion is sound is on the eliminate overload in QuotientGraph.h.
        //
        // What it computed, for the record, since the loop no longer says it:
        //
        //     |C[c] - C[p]| = |C[c]| - sum of weight(u) over u in C[c] & C[p]
        //
        // obtained by subtraction and never by looking at C[c] at all, which is the whole reason
        // the bound is cheap: the quantity depends on the clique and not on the vertex, so every
        // vertex whose incidence list holds c reads it rather than recomputing it.

        const std::size_t numLeft = size - numEliminated;
        for (std::size_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u = pivotClique[k];
            // bound(u) = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]| over I[u] - {p},
            // against the exact |( A[u] | C[c] for c in I[u] ) - {u}|. The first term needs no
            // subtraction: when a clique is formed every member has its explicit adjacency
            // pruned against it, and neither set grows afterwards, so A[u] and C[p] are already
            // disjoint. All the overcounting is therefore clique against clique, outside C[p],
            // which is the smallest place it could have been put.
            // Every length is hoisted out of its condition. These read as free accessors and
            // are not: the loops store, so the compiler cannot prove the store does not alias the
            // size, and the bound is re-loaded per element. `incidenceSize` alone measured 300 ms
            // of this driver's 6.31 s on alpamayo, walked once per member of C[p].
            const std::int32_t* incidence     = qg.incidence(u);
            const std::size_t   incidenceSize = qg.incidenceSize(u);

            // A[u] is not walked here at all: its weight sum was accumulated while the eliminator
            // pruned it, which is what takes that list from two visits to one.
            std::size_t bound = explicitPart[u] + degme - qg.weight(u);
            for (std::size_t i = 0; i < incidenceSize; ++i)
                if (incidence[i] != pivot) bound += outside[incidence[i]];

            // The two caps, both exact and both cheap, and load-bearing rather than defensive:
            // they are what stops the loose term accumulating over a run, which is also why this
            // branch cannot batch (the bound stays tight only while every reached vertex is
            // recomputed at each pivot).
            bound = std::min(bound, numLeft - qg.weight(u));
            bound = std::min(bound, degrees[u] + degme - qg.weight(u));

            buckets.refile(degrees, u, bound);
            minDegree = std::min(minDegree, bound);
        }

        // Nothing above reads an entry this step did not write, since the cliques read are the
        // cliques listed. Clearing anyway keeps that a property of the loop rather than of the
        // reader's memory, and it costs one pass over what was touched.
        for (std::int32_t c : touchedCliques) outside[c] = 0;
    }

    return qg.order(pivots);
}


} // namespace Oblio
