#include "oblio/Amd2.h"

#include <algorithm>
#include <utility>

namespace Oblio {

std::vector<std::int32_t> orderAmd2(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);
    std::vector<std::int32_t> pivots;               // the order over supervariables
    pivots.reserve(size);
    std::size_t numEliminated = 0;                  // a counter, not a scan of the flags

    // The live ORIGINAL vertices, which is what the first cap below is about. It stops agreeing
    // with size - numEliminated the moment a hash merge happens: the merged vertex is eliminated,
    // so it leaves that count, but its weight lives on in the vertex that absorbed it, so it has
    // not left the graph. Deriving the cap from the wrong one of these makes it too tight and the
    // ordering silently worse.
    std::size_t numLive = size;

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
    // Sized for twice the vertex space: the hash comparison below stamps a clique id at c + size,
    // so vertices and cliques can be tested against one stamp without two arrays.
    std::vector<std::int32_t> mark(2 * size, NIL);
    std::int32_t              tag = 0;

    // The hash groups, an array indexed by the hash value rather than a map: the key is already
    // an index into 0 .. size, so a map would cost a log per insertion and a node per group for
    // nothing. Allocated once and cleared only where it was used, which is Amd.cpp's Head[hval].
    // The hash buckets, as head and next arrays rather than a container per bucket. Same idiom as
    // Buckets, and the same reason: a bucket is a list to be pushed onto and walked once, which two
    // flat arrays do in O(1) with nothing allocated. A vector per bucket cost n + 1 headers
    // constructed and destroyed per ordering plus one allocation per bucket the step used, which
    // was 16855 of AMD2's 16911 allocations at 140x140.
    //
    // hashNext is indexed by VERTEX and hashHead by hash, which is why they have different lengths;
    // the same asymmetry Buckets carries between its links and its heads.
    std::vector<std::int32_t> hashHead(size + 1, NIL);
    std::vector<std::int32_t> hashNext(size, NIL);
    std::vector<std::size_t>  usedKeys;

    // |C[c] - C[p]| per clique, indexed by clique id and hoisted out of the loop. The prototype
    // allocates and zeroes it per pivot, which reads better and is O(n) per step, O(n^2) over the
    // run in bookkeeping alone, independent of the graph. Only the entries this step wrote are
    // touched, and they are exactly the ones it will read.
    std::vector<std::size_t>  outside(size, 0);
    std::vector<std::int32_t> touchedCliques;
    std::vector<std::int32_t> deadCliques;

    // |C[c]| per live clique, weighted, which is what the scan below subtracts from. Exact rather
    // than an estimate, and the invariants that keep it so are worth stating because they are not
    // obvious. A live clique never holds an eliminated vertex, since eliminating v absorbs every
    // clique in I[v]. Mass elimination only removes a vertex whose I[u] is exactly {pivot}, so no
    // other clique is touched. And the value is written once, when the clique is formed, from the
    // already-trimmed member list.
    std::vector<std::size_t> cliqueDegree(size, 0);

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        const std::vector<std::int32_t>& merged = qg.eliminate(pivot);
        pivots.push_back(pivot);
        numEliminated += 1 + merged.size();
        numLive -= qg.weight(pivot);                // every original the pivot stands for

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
        for (std::size_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u       = pivotClique[k];
            const std::size_t  weightU = qg.weight(u);
            const std::int32_t* incidence     = qg.incidence(u);
            const std::size_t   incidenceSize = qg.incidenceSize(u);
            for (std::size_t i = 0; i < incidenceSize; ++i) {
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

        // AGGRESSIVE ABSORPTION. outside[c] == 0 says C[c] lies wholly inside the new clique, so
        // it can never contribute anything again and its entries in the incidence lists are pure
        // cost. Amd.cpp's `if (aggressive && we == 0)`. It is worth doing here and nowhere else
        // because the quantity was computed for the bound anyway, so the test is free; and it
        // pays twice over, shortening the lists the bound walks and the lists a later scan walks.
        deadCliques.clear();
        for (std::int32_t c : touchedCliques)
            if (outside[c] == 0) deadCliques.push_back(c);
        qg.absorb(deadCliques, pivotClique, pivotCliqueSize);

        const std::size_t numLeft = numLive;
        for (std::size_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u = pivotClique[k];
            // bound(u) = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]| over I[u] - {p},
            // against the exact |( A[u] | C[c] for c in I[u] ) - {u}|. The first term needs no
            // subtraction: when a clique is formed every member has its explicit adjacency
            // pruned against it, and neither set grows afterwards, so A[u] and C[p] are already
            // disjoint. All the overcounting is therefore clique against clique, outside C[p],
            // which is the smallest place it could have been put.
            // Every length hoisted out of its condition; see the note in Amd1.
            const std::int32_t* adjacency     = qg.adjacency(u);
            const std::size_t   adjacencySize = qg.adjacencySize(u);
            const std::int32_t* incidence     = qg.incidence(u);
            const std::size_t   incidenceSize = qg.incidenceSize(u);

            std::size_t explicitPart = 0;
            for (std::size_t a = 0; a < adjacencySize; ++a)
                explicitPart += qg.weight(adjacency[a]);

            std::size_t bound = explicitPart + degme - qg.weight(u);
            for (std::size_t i = 0; i < incidenceSize; ++i)
                if (incidence[i] != pivot) bound += outside[incidence[i]];

            // The two caps, both exact and both cheap, and load-bearing rather than defensive:
            // they are what stops the loose term accumulating over a run, which is also why this
            // branch cannot batch (the bound stays tight only while every reached vertex is
            // recomputed at each pivot).
            bound = std::min(bound, numLeft - qg.weight(u));
            bound = std::min(bound, degrees[u] + degme - qg.weight(u));

            buckets.refile(degrees, u, bound);
        }

        // HASH SUPERVARIABLE DETECTION. Vertices indistinguishable from EACH OTHER, which the
        // pivot test cannot see: mass elimination only ever finds a vertex indistinguishable from
        // the pivot, and two vertices can become interchangeable with one another without either
        // being interchangeable with it. Hash first so the exact comparison runs only within a
        // group; the hash is a filter and never the decision, so a collision costs a comparison
        // rather than a wrong merge.
        // Filled in REVERSE, and that is load bearing rather than a preference. A chain pushed at
        // the head comes out reversed, and the order within a bucket decides which of two
        // indistinguishable vertices absorbs the other, so filling forward would have been a
        // tie-break change wearing a data-structure change's clothes. It moved the permutation on
        // four of the test graphs before this loop was turned around. Same hazard the degree
        // buckets carry, and the tie-break section of experiments/ordering/README.md describes it.
        usedKeys.clear();
        for (std::size_t k = pivotCliqueSize; k-- > 0;) {
            const std::int32_t u = pivotClique[k];
            if (qg.eliminated(u)) continue;

            // A SUM, because addition has no order and neither do the sets: sorting to build a
            // key would be a log factor for nothing. The two halves are separated by a stride so
            // that a vertex and a clique of the same index cannot cancel.
            std::size_t key = 0;
            const std::int32_t* adjacency = qg.adjacency(u);
            for (std::size_t a = 0; a < qg.adjacencySize(u); ++a)
                if (!qg.eliminated(adjacency[a]))
                    key += static_cast<std::size_t>(adjacency[a]) + 1;
            const std::int32_t* incidence = qg.incidence(u);
            for (std::size_t i = 0; i < qg.incidenceSize(u); ++i)
                key += (static_cast<std::size_t>(incidence[i]) + 1) * (size + 1);

            const std::size_t hash = key % (size + 1);
            if (hashHead[hash] == NIL) usedKeys.push_back(hash);
            hashNext[u]    = hashHead[hash];
            hashHead[hash] = u;
        }

        for (std::size_t hash : usedKeys) {
            for (std::int32_t u = hashHead[hash]; u != NIL; u = hashNext[u]) {
                if (qg.eliminated(u)) continue;
                for (std::int32_t v = hashNext[u]; v != NIL; v = hashNext[v]) {
                    if (qg.eliminated(v)) continue;

                    // The exact test the hash only filters for:
                    //     A[u] - {v} == A[v] - {u}   and   I[u] == I[v]
                    // Decided by stamping one side and counting matches on the other, one pass
                    // and no sort, as every other membership test here is.
                    ++tag;
                    const std::int32_t other = tag;
                    std::size_t sizeV = 0;
                    const std::int32_t* adjacencyV = qg.adjacency(v);
                    for (std::size_t a = 0; a < qg.adjacencySize(v); ++a) {
                        const std::int32_t w = adjacencyV[a];
                        if (w != u && !qg.eliminated(w)) { mark[w] = other; ++sizeV; }
                    }
                    const std::int32_t* incidenceV = qg.incidence(v);
                    for (std::size_t i = 0; i < qg.incidenceSize(v); ++i) {
                        mark[incidenceV[i] + static_cast<std::int32_t>(size)] = other;
                        ++sizeV;
                    }

                    std::size_t sizeU = 0;
                    bool        same  = true;
                    const std::int32_t* adjacencyU = qg.adjacency(u);
                    for (std::size_t a = 0; a < qg.adjacencySize(u) && same; ++a) {
                        const std::int32_t w = adjacencyU[a];
                        if (w == v || qg.eliminated(w)) continue;
                        ++sizeU;
                        if (mark[w] != other) same = false;
                    }
                    if (same) {
                        const std::int32_t* incidenceU = qg.incidence(u);
                        for (std::size_t i = 0; i < qg.incidenceSize(u) && same; ++i) {
                            ++sizeU;
                            if (mark[incidenceU[i] + static_cast<std::int32_t>(size)] != other)
                                same = false;
                        }
                    }
                    if (!same || sizeU != sizeV) continue;

                    buckets.unfile(degrees[v], v);
                    qg.merge(u, v);                 // v folded into u, left where it lies
                    ++numEliminated;                // out of the count, not out of the graph
                }
            }
        }
        for (std::size_t hash : usedKeys) hashHead[hash] = NIL;      // only what was used

        for (std::size_t k = 0; k < pivotCliqueSize; ++k)
            if (!qg.eliminated(pivotClique[k]))
                minDegree = std::min(minDegree, degrees[pivotClique[k]]);

        // Nothing above reads an entry this step did not write, since the cliques read are the
        // cliques listed. Clearing anyway keeps that a property of the loop rather than of the
        // reader's memory, and it costs one pass over what was touched.
        for (std::int32_t c : touchedCliques) outside[c] = 0;
    }

    return qg.order(pivots);
}


} // namespace Oblio
