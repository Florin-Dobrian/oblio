#include "oblio/Amd1.h"

#include <algorithm>
#include <utility>

namespace Oblio {

std::vector<std::int32_t> orderAmd1(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);
    std::vector<std::int32_t> pivots;               // the order over supervariables
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;                // a counter, not a scan of the flags

    // The cache. Exact at construction, since no clique exists yet and the whole neighborhood is
    // still explicit, and a bound from the first elimination onward. That the cached value is
    // itself a bound is what makes the second cap below hold inductively: an upper bound on an
    // earlier degree is still an upper bound now.
    std::vector<std::uint32_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    Buckets buckets(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        buckets.file(degrees[u], u);
    std::uint32_t minDegree = *std::min_element(degrees.begin(), degrees.end());

    // The driver's own membership scratch, separate from the quotient graph's: which vertices lie
    // in the new clique, and which cliques the step has already listed. Both are sets built by
    // stamping and queried by comparison, never allocated.
    std::vector<std::int32_t> mark(size, NIL);
    std::int32_t              tag = 0;

    // |C[c] - C[p]| per clique, indexed by clique id and hoisted out of the loop. The prototype
    // allocates and zeroes it per pivot, which reads better and is O(n) per step, O(n^2) over the
    // run in bookkeeping alone, independent of the graph. Only the entries this step wrote are
    // touched, and they are exactly the ones it will read.
    std::vector<std::uint32_t> outside(size, 0);
    std::vector<std::int32_t> touchedCliques;

    // |C[c]| per live clique, weighted, which is what the scan below subtracts from. Exact rather
    // than an estimate, and the invariants that keep it so are worth stating because they are not
    // obvious. A live clique never holds an eliminated vertex, since eliminating v absorbs every
    // clique in I[v]. Mass elimination only removes a vertex whose I[u] is exactly {pivot}, so no
    // other clique is touched. And the value is written once, when the clique is formed, from the
    // already-trimmed member list.
    std::vector<std::uint32_t> cliqueDegree(size, 0);

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        const std::vector<std::int32_t>& merged = qg.eliminate(pivot);
        pivots.push_back(pivot);
        numEliminated += 1 + static_cast<std::uint32_t>(merged.size());

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
        const std::uint32_t pivotCliqueSize = qg.cliqueSize(pivot);

        // No membership stamp for C[p] is needed any more: the subtraction below asks which
        // cliques the new clique's members belong to, never which vertices a clique's members
        // are, so the one query that used it is gone. The amd2 prototype still carries the stamp,
        // inherited from amd1 and dead there for the same reason.
        // |C[p]| weighted, off the eliminator rather than from a pass of our own. AMD_2
        // accumulates `degme += nvi` while building the element and this is that; the pass this
        // replaces cost one scattered weight load per member per pivot, which is about 6 in 2D
        // and 13 on cubes.
        const std::uint32_t degme = qg.cliqueWeight();
        cliqueDegree[pivot] = degme;                // what the scan below subtracts from

        // |C[c] - C[p]| once per clique. This is the whole reason the bound is cheap: the
        // quantity depends on c alone, so every vertex whose incidence list holds c reads it
        // rather than recomputing it. A second tag makes the clique list a set too, so a clique
        // reached by several vertices is listed once.
        //
        // And it is obtained by SUBTRACTION, never by looking at C[c] at all:
        //
        //     |C[c] - C[p]| = |C[c]| - sum of weight(u) over u in C[c] & C[p]
        //
        // cliqueDegree[c] supplies the first term and the members of C[p] supply the second,
        // since c is in I[u] exactly when u is in C[c]. So the scan walks the incidence lists of
        // the new clique's members and pays sum |I[u]|, where walking the member lists of every
        // touched clique pays sum |C[c]|. Measured on a 100x100 grid, 74281 elements against
        // 272646, which is most of the reason this branch used to run three times slower than the
        // vendored routine. `Amd.cpp` does the same thing at `we = Degree[e] + wnvi`, then
        // `we -= nvi`, and it is the amd2 layer's pass 3.
        touchedCliques.clear();
        ++tag;
        const std::int32_t seenClique = tag;
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u       = pivotClique[k];
            const std::uint32_t weightU = qg.weight(u);
            const std::int32_t* incidence     = qg.incidence(u);
            const std::uint32_t incidenceSize = qg.incidenceSize(u);
            for (std::uint32_t i = 0; i < incidenceSize; ++i) {
                const std::int32_t c = incidence[i];
                if (c == pivot) continue;
                if (mark[c] != seenClique) {        // first sighting: start from |C[c]|
                    mark[c] = seenClique;
                    touchedCliques.push_back(c);
                    outside[c] = cliqueDegree[c] - weightU;
                } else {                            // every later member just subtracts
                    outside[c] -= weightU;
                }
            }
        }

        const std::uint32_t numLeft = static_cast<std::uint32_t>(size) - numEliminated;
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) {
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
            const std::int32_t* adjacency     = qg.adjacency(u);
            const std::uint32_t adjacencySize = qg.adjacencySize(u);
            const std::int32_t* incidence     = qg.incidence(u);
            const std::uint32_t incidenceSize = qg.incidenceSize(u);

            std::uint32_t explicitPart = 0;
            for (std::uint32_t a = 0; a < adjacencySize; ++a)
                explicitPart += qg.weight(adjacency[a]);

            // The seed is not the hazard: A[u] and C[p] are disjoint after the prune, so
            // `explicitPart + degme` is itself at most n. The cap below is the sum that would
            // reach 2n, and it is formed wide; see the note there.
            // WIDE, AND THE LOOP BELOW IS WHY. `bound` accumulates `outside[c]` over I[u], each
            // term up to n and O(n) of them, so the intermediate reaches O(n^2) exactly as Amd3's
            // `deg` does. The two caps afterwards are what bring it back to at most n, and the
            // narrowing therefore belongs after them and not here. The seed is not the hazard:
            // A[u] and C[p] are disjoint after the prune, so `explicitPart + degme` is at most n.
            std::size_t   bound = explicitPart + degme - qg.weight(u);
            for (std::uint32_t i = 0; i < incidenceSize; ++i)
                if (incidence[i] != pivot) bound += outside[incidence[i]];

            // The two caps, both exact and both cheap, and load-bearing rather than defensive:
            // they are what stops the loose term accumulating over a run, which is also why this
            // branch cannot batch (the bound stays tight only while every reached vertex is
            // recomputed at each pivot).
            bound = std::min<std::size_t>(bound, numLeft - qg.weight(u));
            // ONE OPERAND WIDENED, so the sum is formed in `std::size_t`. Widening cannot be
            // done afterwards the way narrowing can: `static_cast<std::size_t>(a + b)` adds
            // in 32 bits and casts the wreckage. One term is enough, the other promoting
            // to meet it. Without it `degrees[u] + degme` reaches 2n and fits only because
            // n is capped at 2^31 - 1, which is a dependency worth not having.
            bound = std::min<std::size_t>(bound,
                                          static_cast<std::size_t>(degrees[u]) + degme
                                              - qg.weight(u));

            // THE NARROWING POINT. Both caps are at most n, so the minimum is too, and the
            // cast is the losing direction, which is the one that gets written.
            const std::uint32_t filed = static_cast<std::uint32_t>(bound);

            buckets.refile(degrees, u, filed);
            minDegree = std::min(minDegree, filed);
        }

        // Nothing above reads an entry this step did not write, since the cliques read are the
        // cliques listed. Clearing anyway keeps that a property of the loop rather than of the
        // reader's memory, and it costs one pass over what was touched.
        for (std::int32_t c : touchedCliques) outside[c] = 0;
    }

    return qg.order(pivots);
}


} // namespace Oblio
