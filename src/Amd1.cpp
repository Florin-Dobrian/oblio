#include "oblio/Amd1.h"

#include "oblio/QuotientGraph.h"

#include <algorithm>
#include <limits>
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
    // AMD_2's `Degree`, and it answers TWO questions from one array. For a live vertex it is the
    // cached degree; for a dead one it is the WEIGHTED SIZE of the clique that vertex's
    // elimination formed, which the scan below subtracts from. The two never overlap: a clique id
    // IS the id of the pivot that made it, and that vertex is dead from the moment the clique
    // exists. A separate `cliqueDegree` was one of the seven n-arrays AMD_2 allocates none of.
    std::vector<std::uint32_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    Buckets buckets(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        buckets.file(degrees[u], u);
    std::uint32_t minDegree = *std::min_element(degrees.begin(), degrees.end());

    // AMD.CPP'S W ARRAY, in place of `mark[c]` plus `outside[c]`, exactly as Amd1 now carries it.
    // The fused schedule was carried in a separate Amd1B layer until 2026-08-16, when it
    // measured permutation-identical and faster and was moved here. `mark` goes entirely: this
    // layer has no hash, so the stamp had no other reader. See src/Amd1.cpp for the encoding and
    // why the absorbed state is unused here.
    std::vector<std::int32_t> w(size, 1);       // every clique alive and unseen, Amd.cpp's W
    std::int32_t wflg  = 2;                     // the tag, Amd.cpp's wflg
    std::int32_t lemax = 0;                     // the largest clique so far, Amd.cpp's lemax
    const std::int32_t wbig = std::numeric_limits<std::int32_t>::max()
                              - static_cast<std::int32_t>(size);
    const auto clearFlag = [&]() {
        if (wflg < 2 || wflg >= wbig) {
            for (std::int32_t x = 0; x < static_cast<std::int32_t>(size); ++x) w[x] = 1;
            wflg = 2;
        }
    };
    std::vector<std::int32_t> touchedCliques;

    // |C[c]| per live clique, weighted, which is what the scan below subtracts from. Exact rather
    // than an estimate, and the invariants that keep it so are worth stating because they are not
    // obvious. A live clique never holds an eliminated vertex, since eliminating v absorbs every
    // clique in I[v]. Mass elimination only removes a vertex whose I[u] is exactly {pivot}, so no
    // other clique is touched. And the value is written once, when the clique is formed, from the
    // already-trimmed member list.

    // The bound's explicit term, accumulated by the eliminator while it prunes A[u] rather than by
    // a second walk afterwards. Written for every reached vertex and read for every survivor.
    // The adjacency half of each bound now rides in `w[u]`, free for a live vertex because `w` is
    // indexed by clique id and a clique id is a dead pivot's. The array that carried it is gone;
    // the obligation is the reset at the end of the bound pass. See src/Amd3.cpp.

    // THE TAGGED SCAN, not the ApproximateScan this layer used to hand the eliminator. Same
    // fusion, the driver's first scan folded into the eliminator's walk, on the tagged W encoding
    // instead of a stamp array plus a value array. `nullptr` for the buckets: that arrangement
    // takes every member of C[pivot] out of the degree lists and parks a hash key in the link it
    // frees, and this layer has no hash and refiles inside its own bound pass. See TaggedScan in
    // QuotientGraph.h.
    //
    // Bound once and reused: only `wflg` moves, and the driver sets it before each elimination
    // exactly as it would before its own scan.
    TaggedScan scan{nullptr, w, degrees, touchedCliques, wflg,
                    static_cast<std::int32_t>(size)};

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        // The scan runs inside the elimination, so its stamp and its list are prepared before the
        // call rather than after it. Everything else about the step is unchanged.
        clearFlag();                                   // Amd.cpp's clear_flag; almost never fires
        touchedCliques.clear();
        scan.wflg = wflg;

        const std::vector<std::int32_t>& merged = qg.eliminate(pivot, scan);
        pivots.push_back(pivot);
        numEliminated += 1 + static_cast<std::uint32_t>(merged.size());

        // The pivot leaves the lists. The zeroing that used to follow is gone: under the fold above
        // `degrees[pivot]` is the slot the new clique's weight is written into a few lines down,
        // so it was a store nobody read. The old comment warned to unfile before zeroing because
        // the bucket index came from the degree; Buckets reads it out of mPrev, so that ordering
        // was already vestigial.
        buckets.unfile(pivot);
        for (std::int32_t u : merged) {
            buckets.unfile(u);
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
        std::uint32_t degme = 0;                    // |C[p]|, weighted
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) degme += qg.weight(pivotClique[k]);
        degrees[pivot] = degme;                     // what the scan below subtracts from
        lemax = std::max(lemax, static_cast<std::int32_t>(degme));

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
            const std::int32_t* incidence     = qg.incidence(u);
            const std::uint32_t incidenceSize = qg.incidenceSize(u);

            // A[u] is not walked here at all: its weight sum was accumulated while the eliminator
            // pruned it, which is what takes that list from two visits to one.
            // The seed is at most n; the cap below is the sum that reaches 2n. See Amd1.
            // WIDE, AND THE LOOP BELOW IS WHY. `bound` accumulates `outside[c]` over I[u], each
            // term up to n and O(n) of them, so the intermediate reaches O(n^2) exactly as Amd3's
            // `deg` does. The two caps afterwards are what bring it back to at most n, and the
            // narrowing therefore belongs after them and not here. The seed is not the hazard:
            // A[u] and C[p] are disjoint after the prune, so `explicitPart + degme` is at most n.
            std::size_t   bound = static_cast<std::size_t>(w[u]) + degme - qg.weight(u);
            for (std::uint32_t i = 0; i < incidenceSize; ++i)
                if (incidence[i] != pivot)
                    bound += static_cast<std::size_t>(w[incidence[i]] - wflg);

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

            // THE SLOT GOES BACK TO ALIVE-AND-UNSEEN, the bound half having just been read out
            // of it. Without this a survivor later chosen as pivot would form a clique whose w
            // already held a bound. See src/Amd3.cpp.
            w[u] = 1;
            buckets.refile(degrees, u, filed);
            minDegree = std::min(minDegree, filed);
        }

        // Nothing above reads an entry this step did not write, since the cliques read are the
        // cliques listed. Clearing anyway keeps that a property of the loop rather than of the
        // reader's memory, and it costs one pass over what was touched.
        // THE TAG ADVANCES, which is what replaces the clearing pass this step used to end with.
        wflg += lemax;
    }

    return qg.order(pivots);
}


} // namespace Oblio
