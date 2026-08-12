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
    std::uint32_t numEliminated = 0;                // a counter, not a scan of the flags

    // The live ORIGINAL vertices, which is what the first cap below is about. It stops agreeing
    // with size - numEliminated the moment a hash merge happens: the merged vertex is eliminated,
    // so it leaves that count, but its weight lives on in the vertex that absorbed it, so it has
    // not left the graph. Deriving the cap from the wrong one of these makes it too tight and the
    // ordering silently worse.
    std::uint32_t numLive = static_cast<std::uint32_t>(size);

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
    // Sized for twice the vertex space: the hash comparison below stamps a clique id at
    // c + cliqueStamp, so vertices and cliques can be tested against one stamp without two arrays.
    std::vector<std::int32_t> mark(2 * size, NIL);
    std::int32_t              tag = 0;

    // The stride separating the two halves of `mark`, and the one place a COUNT becomes an offset
    // in the INDEX space. `size` is the matrix order, one dimensional and bounded by n, so it is a
    // count and is held as std::size_t like every other count here; the sum below is an index into
    // mark and so is std::int32_t. Naming the crossing once beats writing the cast at each of the
    // sites that make it, which is what this file used to do: the cast is not a hazard here, since
    // n is an int32_t by construction, but four unexplained casts of the same quantity read as
    // four separate events rather than one convention.
    //
    // It exists because there is no type for a count. docs/DESIGN_DECISIONS.md (2026-08-08) and
    // experiments/ordering/REPORT.md carry that: an index names an entity and may be NIL, a
    // position offsets into an n x n object and may exceed 2^31, and a count is bounded by a SIDE
    // rather than an AREA and has neither category. With one, `size` would already be the right
    // width and this line would not be needed.
    const std::int32_t cliqueStamp = static_cast<std::int32_t>(size);

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
    std::vector<std::uint32_t> usedKeys;

    // |C[c] - C[p]| per clique, indexed by clique id and hoisted out of the loop. The prototype
    // allocates and zeroes it per pivot, which reads better and is O(n) per step, O(n^2) over the
    // run in bookkeeping alone, independent of the graph. Only the entries this step wrote are
    // touched, and they are exactly the ones it will read.
    std::vector<std::uint32_t> outside(size, 0);
    std::vector<std::int32_t> touchedCliques;
    std::vector<std::int32_t> deadCliques;

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

        // AGGRESSIVE ABSORPTION. outside[c] == 0 says C[c] lies wholly inside the new clique, so
        // it can never contribute anything again and its entries in the incidence lists are pure
        // cost. Amd.cpp's `if (aggressive && we == 0)`. It is worth doing here and nowhere else
        // because the quantity was computed for the bound anyway, so the test is free; and it
        // pays twice over, shortening the lists the bound walks and the lists a later scan walks.
        deadCliques.clear();
        for (std::int32_t c : touchedCliques)
            if (outside[c] == 0) deadCliques.push_back(c);
        qg.absorb(deadCliques, pivotClique, pivotCliqueSize);

        const std::uint32_t numLeft = numLive;
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u = pivotClique[k];
            // bound(u) = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]| over I[u] - {p},
            // against the exact |( A[u] | C[c] for c in I[u] ) - {u}|. The first term needs no
            // subtraction: when a clique is formed every member has its explicit adjacency
            // pruned against it, and neither set grows afterwards, so A[u] and C[p] are already
            // disjoint. All the overcounting is therefore clique against clique, outside C[p],
            // which is the smallest place it could have been put.
            // Every length hoisted out of its condition; see the note in Amd1.
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
            // The minimum, taken HERE rather than in a pass of its own after the hash. `bound` is
            // in a register, where that pass paid a scattered read per survivor to recover it.
            // AMD_2 takes its minimum inside the same loop, `if (deg < mindeg)`.
            //
            // TWO SITES, and that is not tidiness: a hash merge below LOWERS a survivor's degree
            // after this loop has run, which is exactly why the pass being removed sat at the end.
            // Folding it here alone would miss those, so the merge takes the minimum too, off the
            // value it computes anyway. Amd3 needs one site only, its refile pass running after
            // the hash rather than before it.
            //
            // MEASURED AT ZERO, with the clique-weight fusion beside it: useful cycles unchanged
            // within half a percent in both families. A port and a simplification, not a speed
            // fix. See benchmarks/ordering/README.md (2026-08-09).
            minDegree = std::min(minDegree, filed);
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
        for (std::uint32_t k = pivotCliqueSize; k-- > 0;) {
            const std::int32_t u = pivotClique[k];
            if (qg.eliminated(u)) continue;

            // A SUM, because addition has no order and neither do the sets: sorting to build a
            // key would be a log factor for nothing. The two halves are separated by a stride so
            // that a vertex and a clique of the same index cannot cancel.
            //
            // Built in a pass of its own, and it was FUSED INTO THE BOUND LOOP ABOVE AND REVERTED
            // on 2026-08-08. Amd.cpp accumulates its key in the walks it is already making,
            // `hval += e` and `hval += j`, and REPORT.md had measured this separate traversal at
            // 72 percent of AMD2's overhead in 2D and 92 in 3D, so the fusion looked like the
            // whole answer. Measured on alpamayo it bought NOTHING at 140 a side and cost 2
            // percent at 400, which is the footprint trade REPORT attached as its own caution:
            // the key has to be carried in an array of size n, the same stream that made Amd1B
            // slower at large n after being faster at small. A tenth of the driver's element
            // visits went and the array ate it.
            std::size_t key = 0;
            const std::int32_t* adjacency = qg.adjacency(u);
            for (std::uint32_t a = 0; a < qg.adjacencySize(u); ++a)
                if (!qg.eliminated(adjacency[a]))
                    key += static_cast<std::size_t>(adjacency[a]) + 1;
            // ONE SUM, WITH NO STRIDE, and that is ledger entry 8. This added the incidence half
            // as `(c + 1) * (size + 1)` until 2026-08-09, so that a vertex and a clique of the same
            // index could not cancel. True of the KEY and false of the BUCKET: the modulus below
            // is the same number as the stride, so the incidence term is annihilated exactly and
            // the hash came out a function of the ADJACENCY ALONE. As the elimination proceeds
            // A[u] empties and everything a vertex reaches becomes cliques, so the surviving key
            // carried less and less, and cubic grids reach that state sooner than square ones.
            //
            // Measured, for the SAME MERGES: 19.0 pairs tested per pivot at 140 a side against the
            // vendored routine's 0.33, and 155.3 at 26 cubed against its 0.48. Amd.cpp accumulates
            // `hval += e` and `hval += j` into one running value and takes it mod n, letting a
            // vertex and a clique collide on purpose, because the hash is a FILTER and never the
            // decision: a collision costs one exact comparison and cannot produce a wrong merge.
            // The invariant the two lines have to hold TOGETHER is that the modulus must not
            // divide the stride, and having no stride is the cheapest way to hold it.
            const std::int32_t* incidence = qg.incidence(u);
            for (std::uint32_t i = 0; i < qg.incidenceSize(u); ++i)
                key += static_cast<std::size_t>(incidence[i]) + 1;

            // THE NARROWING POINT for the key. `key` is one of the five wide accumulators, this
            // one summing `c + 1` over A[u] and I[u]; the remainder is under `size + 1` and so
            // at most n, which is what makes the bucket index one dimensional.
            const std::uint32_t hash = static_cast<std::uint32_t>(key % (size + 1));
            if (hashHead[hash] == NIL) usedKeys.push_back(hash);
            hashNext[u]    = hashHead[hash];
            hashHead[hash] = u;
        }

        for (std::uint32_t hash : usedKeys) {
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
                    std::uint32_t sizeV = 0;   // list entries, so at most deg(v)
                    const std::int32_t* adjacencyV = qg.adjacency(v);
                    for (std::uint32_t a = 0; a < qg.adjacencySize(v); ++a) {
                        const std::int32_t w = adjacencyV[a];
                        if (w != u && !qg.eliminated(w)) { mark[w] = other; ++sizeV; }
                    }
                    const std::int32_t* incidenceV = qg.incidence(v);
                    for (std::uint32_t i = 0; i < qg.incidenceSize(v); ++i) {
                        mark[incidenceV[i] + cliqueStamp] = other;
                        ++sizeV;
                    }

                    std::uint32_t sizeU = 0;   // list entries, so at most deg(u)
                    bool        same  = true;
                    const std::int32_t* adjacencyU = qg.adjacency(u);
                    for (std::uint32_t a = 0; a < qg.adjacencySize(u) && same; ++a) {
                        const std::int32_t w = adjacencyU[a];
                        if (w == v || qg.eliminated(w)) continue;
                        ++sizeU;
                        if (mark[w] != other) same = false;
                    }
                    if (same) {
                        const std::int32_t* incidenceU = qg.incidence(u);
                        for (std::uint32_t i = 0; i < qg.incidenceSize(u) && same; ++i) {
                            ++sizeU;
                            if (mark[incidenceU[i] + cliqueStamp] != other)
                                same = false;
                        }
                    }
                    if (!same || sizeU != sizeV) continue;

                    // The TARGET is a live vertex and never the pivot, which is the whole
                    // difference from mass elimination. Both draw from C[p], the clique this
                    // step formed, and the pivot is not a member of its own clique. So mass
                    // elimination folds into a vertex that is leaving; this folds into one that
                    // stays, carrying the combined weight and still a candidate.
                    //
                    // u's REACHABLE SET is unchanged by the merge, and that half is worth
                    // keeping: the two were adjacent to each other, v leaves the graph, and an
                    // external degree excludes u's own supervariable, so what u can reach is
                    // exactly what it was. Every other member of C[p] is unaffected too, since v
                    // was in C[p] and its weight has moved to u, so |C[p]| weighted is unchanged
                    // and so is the middle term of their bounds.
                    //
                    // And the degree IS recomputed, by one subtraction. The bound a few lines
                    // above was written with u's weight as it stood BEFORE this merge, and the
                    // weight appears inside the bound, in `degme - weight(u)` and in the
                    // `numLeft - weight(u)` cap. So absorbing v leaves the filed value one
                    // bucket too high per original vertex taken, and u is never picked as early
                    // as its size has earned.
                    //
                    // This comment used to say the opposite, that nothing was stale because an
                    // external degree excludes u's own supervariable and only the WEIGHT
                    // changes. The first half is true and the second is the whole difficulty:
                    // the buckets are keyed on a degree that has the weight subtracted in it.
                    // AMD_2 has no such problem because it subtracts `nvi` in the pass that
                    // restores the degree lists, which runs AFTER supervariable detection and so
                    // reads the post-merge weight. Found by aligning the amd3 prototype against
                    // it, where the same timing is ledger entry 4.
                    //
                    // One subtraction rather than a recomputation, because all three terms of
                    // the bound's minimum shift by the same amount when the weight grows, so the
                    // minimum shifts with them. It remains an upper bound for the same reason:
                    // the true external degree drops by exactly weight(v) as v stops being
                    // outside u's supervariable and becomes part of it.
                    const std::uint32_t weightV = qg.weight(v);
                    buckets.unfile(degrees[v], v);
                    qg.merge(u, v);                 // v folded into u, left where it lies
                    const std::uint32_t merged = degrees[u] - weightV;
                    buckets.refile(degrees, u, merged);
                    minDegree = std::min(minDegree, merged);   // see the bound loop above
                    ++numEliminated;                // out of the count, not out of the graph
                }
            }
        }
        for (std::uint32_t hash : usedKeys) hashHead[hash] = NIL;     // only what was used


        // Nothing above reads an entry this step did not write, since the cliques read are the
        // cliques listed. Clearing anyway keeps that a property of the loop rather than of the
        // reader's memory, and it costs one pass over what was touched.
        for (std::int32_t c : touchedCliques) outside[c] = 0;
    }

    return qg.order(pivots);
}


} // namespace Oblio
