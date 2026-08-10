#include "oblio/Amd3.h"

#include <algorithm>
#include <limits>
#include <utility>

namespace Oblio {

std::vector<std::int32_t> orderAmd3(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);

    // The two shared-class conventions this layer differs by. Both are off for every other driver
    // and neither changes which sets are computed, only which permutation comes out. See the
    // setters, and experiments/ordering/AMD3.md for ledger entries 2, 3 and 5.
    qg.setVendoredListOrder(true);      // cliques before adjacency; the new clique at the front
    qg.setLateMassElimination(true);    // and mass elimination becomes this driver's, below

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
    std::vector<std::size_t>  usedKeys;

    // |C[c] - C[p]| per clique, indexed by clique id and hoisted out of the loop. The prototype
    // allocates and zeroes it per pivot, which reads better and is O(n) per step, O(n^2) over the
    // run in bookkeeping alone, independent of the graph. Only the entries this step wrote are
    // touched, and they are exactly the ones it will read.
    std::vector<std::int32_t> touchedCliques;
    std::vector<std::int32_t> deadCliques;

    // |C[c]| per live clique, weighted, which is what the scan below subtracts from. Exact rather
    // than an estimate, and the invariants that keep it so are worth stating because they are not
    // obvious. A live clique never holds an eliminated vertex, since eliminating v absorbs every
    // clique in I[v]. Mass elimination only removes a vertex whose I[u] is exactly {pivot}, so no
    // other clique is touched. And the value is written once, when the clique is formed, from the
    // already-trimmed member list.
    std::vector<std::size_t> cliqueDegree(size, 0);

    // Amd.cpp's W ARRAY, and this file used to keep the same three facts in three places.
    //
    // The old shape kept THREE facts about a clique in three places: `mark[c]` says whether this step has
    // seen it, `outside[c]` carries |C[c] - C[p]|, and a clique is dead when it has been stripped
    // from the incidence lists. Its scan 1 therefore loads and stores `mark[c]` AND loads and
    // stores `outside[c]` for every incidence element, two cache lines per element, and clears
    // `outside` over the touched list at the end of every step.
    //
    // Amd.cpp keeps all three in ONE array, with a tag:
    //
    //     w[c] == 0            the clique is absorbed and gone
    //     0 < w[c] < wflg      alive, not seen this step; the value is stale
    //     w[c] >= wflg         seen this step, and w[c] - wflg is |C[c] - C[p]|
    //
    // so its scan 1 is one load and one store into one array, and `cliqueDegree` is touched only
    // on first sighting. Its inner body is four lines and this is a transcription of them:
    //
    //     we = W [e] ;
    //     if      (we >= wflg) we -= nvi ;
    //     else if (we != 0)    we = Degree [e] + wnvi ;
    //     W [e] = we ;
    //
    // The tag advances by `lemax` at the end of each step, which is the largest clique this run
    // has built. That is what makes the stale range safe without a clearing pass: after scan 1 no
    // entry exceeds wflg + lemax, so advancing by lemax puts every one of them below the new wflg,
    // and the whole array is invalidated in a single addition.
    //
    // SIGNED, and int32 to match. `wnvi = wflg - weight(u)` is negative whenever the tag is still
    // small and the weight is not, which is exactly the case on the first eliminations, and the
    // arithmetic only comes right again at `w[c] - wflg`. An unsigned type wraps there and the
    // bound comes out enormous. Amd.cpp's `Int` is signed for the same reason.
    std::vector<std::int32_t> w(size, 1);       // every clique alive and unseen, Amd.cpp's W
    std::int32_t wflg  = 2;                     // the tag, Amd.cpp's wflg
    std::int32_t lemax = 0;                     // the largest clique so far, Amd.cpp's lemax
    // wflg + n must not overflow, which is the whole of Amd.cpp's wbig.
    const std::int32_t wbig = std::numeric_limits<std::int32_t>::max() - static_cast<std::int32_t>(size);
    std::size_t numFlagSweeps = 0;              // how often the guard below actually fires

    // Amd.cpp's clear_flag: reset the array and the tag when the tag can no longer be advanced
    // safely. Every live clique goes back to 1, which is the alive-and-unseen state, and the dead
    // ones stay 0. Called once per elimination, and almost never does anything.
    const auto clearFlag = [&]() {
        if (wflg < 2 || wflg >= wbig) {
            for (std::size_t x = 0; x < size; ++x)
                if (w[x] != 0) w[x] = 1;
            wflg = 2;
            ++numFlagSweeps;
        }
    };

    // The half of each bound that does not involve the vertex's own weight, carried from the pass
    // that forms it across supervariable detection to the pass that finishes it. Amd.cpp keeps the
    // same quantity in Degree[i] between its scan 2 and its degree-list pass. See ledger entry 4.
    std::vector<std::size_t> partial(size, 0);

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        // Under setLateMassElimination this returns an empty list and C[pivot] is reach(pivot)
        // exactly. The merge happens below, once the absorption has run.
        qg.eliminate(pivot);
        pivots.push_back(pivot);

        buckets.unfile(degrees[pivot], pivot);      // unfile before zeroing: the bucket index is
        degrees[pivot] = 0;                         //   read from the degree

        // ---- the bound, in place of an exact refresh -----------------------------------
        // Everything the new clique reached needs a new degree, and nothing else can have
        // changed. C[p] is therefore both the refresh set and the domain the decomposition is a
        // statement about, which is the coincidence the whole placement rests on.
        // Read at the moment of use: the arena holding it grows as cliques are formed, so a
        // pointer taken before the next elimination is the only one that is safe, and this is
        // that window.
        const std::int32_t* pivotClique     = qg.clique(pivot);
        std::size_t         pivotCliqueSize = qg.cliqueSize(pivot);

        // No membership stamp for C[p] is needed any more: the subtraction below asks which
        // cliques the new clique's members belong to, never which vertices a clique's members
        // are, so the one query that used it is gone. The amd2 prototype still carries the stamp,
        // inherited from amd1 and dead there for the same reason.
        // |C[p]| weighted is NOT taken here. It is the value after mass elimination, which now runs
        // below the absorption, so it is read there. The scan that follows is deliberately over
        // the UNTRIMMED clique, which is what Amd.cpp's scan 1 walks: it runs over the whole of
        // Lme before any of it has been mass eliminated. Nothing is lost by that, since a vertex
        // the merge will take belongs to no clique but the new one and so cannot appear in any
        // touched clique's member list either way.
        // |C[p]| weighted, off the eliminator rather than from a pass of our own. AMD_2
        // accumulates `degme += nvi` while building the element and this is that. The second
        // computation below is NOT removable: it runs after mass elimination has trimmed the
        // clique, which is ledger entry 7, and this one is deliberately over the untrimmed one.
        std::size_t degme = qg.cliqueWeight();
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
        clearFlag();                            // Amd.cpp calls this here too, before scan 1
        lemax = std::max(lemax, static_cast<std::int32_t>(degme));

        touchedCliques.clear();
        for (std::size_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u    = pivotClique[k];
            const std::int32_t nvi  = static_cast<std::int32_t>(qg.weight(u));
            const std::int32_t wnvi = wflg - nvi;
            const std::int32_t* incidence     = qg.incidence(u);
            const std::size_t   incidenceSize = qg.incidenceSize(u);
            for (std::size_t i = 0; i < incidenceSize; ++i) {
                const std::int32_t c = incidence[i];
                if (c == pivot) continue;
                std::int32_t we = w[c];
                if (we >= wflg) {                   // already seen this step: just subtract
                    we -= nvi;
                } else if (we != 0) {               // first sighting: start from |C[c]|, tagged
                    we = static_cast<std::int32_t>(cliqueDegree[c]) + wnvi;
                    touchedCliques.push_back(c);    // only for the absorption pass below
                }
                w[c] = we;
            }
        }

        // AGGRESSIVE ABSORPTION. w[c] - wflg == 0 says C[c] lies wholly inside the new clique, so
        // it can never contribute anything again and its entries in the incidence lists are pure
        // cost. Amd.cpp's `if (aggressive && we == 0)`. It is worth doing here and nowhere else
        // because the quantity was computed for the bound anyway, so the test is free; and it
        // pays twice over, shortening the lists the bound walks and the lists a later scan walks.
        deadCliques.clear();
        for (std::int32_t c : touchedCliques)
            if (w[c] == wflg) { deadCliques.push_back(c); w[c] = 0; }   // |C[c] - C[p]| == 0
        qg.absorb(deadCliques, pivotClique, pivotCliqueSize);

        // MASS ELIMINATION, and it runs HERE rather than inside the eliminator. Absorption is what
        // makes the cheap structural test agree with the true one: a clique whose members all lie
        // inside C[p] contributes nothing to what u can reach, yet its presence in I[u] makes the
        // test fail. Amd.cpp says so itself, making the same test in its scan 2 over an element
        // list absorption has already compacted: with aggressive absorption, `deg == 0` is
        // identical to `Elen[i] == 1 && p3 == pn`. Asking first, as every other driver here does,
        // declines merges the vendored routine makes. experiments/ordering/AMD3.md, ledger entry 3.
        const std::vector<std::int32_t>& merged = qg.massEliminate(pivot);
        numEliminated += 1 + merged.size();
        numLive -= qg.weight(pivot);                // every original the pivot stands for
        for (std::int32_t u : merged) {
            buckets.unfile(degrees[u], u);
            degrees[u] = 0;
        }

        // The clique and its weight are re-read, both having moved: massEliminate trims C[p] of
        // the merged members and folds their weight into the pivot. Amd.cpp reaches the same
        // value by decrementing degme inside scan 2, but it does not CONSUME it there: the term
        // enters a survivor's degree only in the pass that restores the degree lists,
        // `deg = Degree[i] + degme - nvi`, by which point degme is final. So every survivor sees
        // the same number, which is what re-taking it here gives.
        pivotClique     = qg.clique(pivot);
        pivotCliqueSize = qg.cliqueSize(pivot);
        degme = 0;
        for (std::size_t k = 0; k < pivotCliqueSize; ++k) degme += qg.weight(pivotClique[k]);

        // AND THE STORED CLIQUE DEGREE IS WRITTEN AGAIN, which is not a tidy-up. `Amd.cpp` writes
        // `Degree [me] = degme` TWICE, at its lines 1676 and 1940, and the second write is the
        // durable one: by then scan 2 has run `degme -= nvi` for every vertex mass elimination
        // took, so what a later step reads as |C[me]| is the post-merge size. The write above the
        // scan holds the pre-merge size, which is what the scan itself must subtract from, since
        // it runs over the untrimmed clique.
        //
        // We had only the first write, so any pivot that mass-eliminated left a clique degree
        // permanently too large by the merged weight, and every later `dext = cliqueDegree[c] -
        // ...` taken through that clique inherited it. That inflates a bound, and an inflated
        // bound moves the ordering only when it moves the head of the minimum bucket, which is
        // why it is invisible on 2D grids at every size to 140 a side and first surfaces on a 3D
        // grid at 16.
        //
        // It is this driver's alone. Amd1 and Amd2 mass-eliminate inside the eliminator, so their
        // clique is already trimmed when they take `degme` and their single write is correct. The
        // defect arrived with ledger entry 3, which moved mass elimination out and did not carry
        // the second write that placement is the whole reason for. Half a mechanism, as entry 6
        // was. See experiments/ordering/AMD3.md.
        cliqueDegree[pivot] = degme;

        const std::size_t numLeft = numLive;

        // The bound is formed in TWO halves here, where Amd2 forms it in one, and that split is
        // ledger entry 4. This pass computes the part that does not involve u's own weight and
        // stores it; the pass after the hash adds `degme` and subtracts the weight. Amd.cpp does
        // the same: scan 2 forms `deg = sum dext(e) + sum nvj` and keeps
        // `Degree[i] = MIN(Degree[i], deg)`, and only the later pass finishes it. The reason it
        // is load-bearing rather than a rearrangement is that supervariable detection runs BETWEEN
        // the two and a hash merge does `Nv[i] += Nv[j]`, so `nvi` there is the POST-merge weight.
        // Subtracting it in this pass would use the weight u had before absorbing v, which files a
        // supervariable one bucket too high per vertex taken. Amd2 and Amd2B carried exactly that
        // and it was costing 3 to 9 percent of fill on grids.
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

            std::size_t deg = explicitPart;
            for (std::size_t i = 0; i < incidenceSize; ++i)
                if (incidence[i] != pivot)
                    deg += static_cast<std::size_t>(w[incidence[i]] - wflg);

            // Amd.cpp's `Degree[i] = MIN (Degree[i], deg)`. The stored degree is a full one from
            // an earlier step and this is a partial, so the two are comparable only once the pass
            // below adds `degme - weight(u)` to whichever won. That is why the minimum is taken
            // here and the common term added there rather than the other way round. The second
            // cap, `numLeft - weight(u)`, is also that pass's, for the same weight reason.
            partial[u] = std::min(deg, degrees[u]);
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
        for (std::size_t k = 0; k < pivotCliqueSize; ++k) {
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
            for (std::size_t a = 0; a < qg.adjacencySize(u); ++a)
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
            for (std::size_t i = 0; i < qg.incidenceSize(u); ++i)
                key += static_cast<std::size_t>(incidence[i]) + 1;

            const std::size_t hash = key % (size + 1);
            if (hashHead[hash] == NIL) usedKeys.push_back(hash);
            hashNext[u]    = hashHead[hash];
            hashHead[hash] = u;
        }

        for (std::size_t hash : usedKeys) {
            for (std::int32_t u = hashHead[hash]; u != NIL; u = hashNext[u]) {
                if (qg.eliminated(u)) continue;
                // AMD_2 enters this loop only for a bucket member with a successor,
                // `while (i != EMPTY && Next [i] != EMPTY)`, and the guard is the other half of
                // the hoisted stamp below: without it a member with nothing after it pays a full
                // list of random writes for a pair that will never be tested. It cannot change
                // the answer, since the inner loop is empty in exactly the cases it skips. Inert
                // while the buckets were enormous, which is why it was never missed, and most of
                // the pass once entry 8 made them singletons. Amd2 and Amd2B need no counterpart:
                // they stamp INSIDE the pair loop, so a member with no successor already costs
                // them nothing.
                if (hashNext[u] == NIL) continue;

                // THE STAMP IS HOISTED, which is the second thing taken from Amd.cpp here. Its supervariable detection stamps the OUTER vertex once, before the
                // inner loop, and then tests every candidate against that one stamp:
                //
                //     for (p = Pe [i] + 1 ; ... ) W [Iw [p]] = wflg ;   /* i, once */
                //     while (j != EMPTY) { ok = ... ; for (p = Pe [j] + 1 ; ok && ... ) ... }
                //
                // This file used to stamp the INNER vertex, once per PAIR, over its whole list and with no
                // short-circuit, so every pair pays a full list of random writes even though
                // entry 6 got the comparison itself down to 1.08 iterations. Measured on a 140x140
                // grid: 639083 stamp writes before against 290473 after, with 263032 compare
                // iterations either way. We were stamping 2.4 times more elements than we compared.
                //
                // WHY IT CAN BE HOISTED, and why it took a day to see. The stamping carried
                // `w != u` and the walk carried `w == v`, exclusions that look pair-dependent and
                // so look to pin the stamp inside the loop. They are vestigial: u and v are both
                // members of C[pivot], and the prune drops every neighbour lying inside the new
                // clique, `if (mMark[v] == inClique) continue`, so A[u] cannot contain v and A[v]
                // cannot contain u. Nothing to exclude. That is exactly why Amd.cpp's stamp has no
                // such guard either. This was landed as a separate driver first, so that
                // `AMD3C == AMD3` could say the reasoning held rather than merely sounding right,
                // and folded in once it did.
                //
                // Roles swapped with the hoist: u is stamped and v is walked, where before it was
                // v stamped and u walked. The test is symmetric so the outcome does not move, and the SURVIVOR
                // does not either, u being the outer vertex in both and merge(u, v) folding v into
                // it. Amd.cpp merges j into i the same way round.
                ++tag;
                const std::int32_t other = tag;
                std::size_t sizeU = 0;
                const std::int32_t* adjacencyU = qg.adjacency(u);
                for (std::size_t a = 0; a < qg.adjacencySize(u); ++a) {
                    const std::int32_t w = adjacencyU[a];
                    if (!qg.eliminated(w)) { mark[w] = other; ++sizeU; }
                }
                // Index 1: the new clique is at the front of every I[u] and is shared by every
                // member of C[pivot], so it can never discriminate. Ledger entry 6.
                const std::int32_t* incidenceU = qg.incidence(u);
                for (std::size_t i = 1; i < qg.incidenceSize(u); ++i) {
                    mark[incidenceU[i] + cliqueStamp] = other;
                    ++sizeU;
                }

                for (std::int32_t v = hashNext[u]; v != NIL; v = hashNext[v]) {
                    if (qg.eliminated(v)) continue;

                    // The exact test the hash only filters for:
                    //     A[u] == A[v]   and   I[u] == I[v]
                    // against the stamp of u laid down once above. Both walks short-circuit on the
                    // first mismatch, which is what made the comparison cheap and the stamping the
                    // thing that had to move: see experiments/ordering/AMD3.md, iteration 15.
                    std::size_t sizeV = 0;
                    bool        same  = true;
                    const std::int32_t* adjacencyV = qg.adjacency(v);
                    for (std::size_t a = 0; a < qg.adjacencySize(v) && same; ++a) {
                        const std::int32_t w = adjacencyV[a];
                        if (qg.eliminated(w)) continue;
                        ++sizeV;
                        if (mark[w] != other) same = false;
                    }
                    if (same) {
                        const std::int32_t* incidenceV = qg.incidence(v);
                        for (std::size_t i = 1; i < qg.incidenceSize(v) && same; ++i) {
                            ++sizeV;
                            if (mark[incidenceV[i] + cliqueStamp] != other) same = false;
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
                    const std::size_t weightV = qg.weight(v);
                    buckets.unfile(degrees[v], v);
                    qg.merge(u, v);                 // v folded into u, left where it lies
                    buckets.refile(degrees, u, degrees[u] - weightV);
                    ++numEliminated;                // out of the count, not out of the graph
                }
            }
        }
        for (std::size_t hash : usedKeys) hashHead[hash] = NIL;      // only what was used

        // THE FOURTH PASS, which finishes the bounds and files them. Amd.cpp spells it
        // `deg = Degree[i] + degme - nvi` then `deg = MIN (deg, nleft - nvi)`, under its RESTORE
        // DEGREE LISTS heading, and it runs here for the reason in the note above: `nvi` is read
        // after supervariable detection, so a vertex that absorbed another subtracts the combined
        // weight. `degme` and `numLeft` were settled before the hash and are the same for
        // everyone, so this is one addition and two comparisons per survivor.
        for (std::size_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u = pivotClique[k];
            if (qg.eliminated(u)) continue;         // absorbed by the hash a moment ago
            const std::size_t  weightU = qg.weight(u);   // POST-merge, which is the whole point
            std::size_t bound = partial[u] + degme - weightU;
            bound = std::min(bound, numLeft - weightU);
            buckets.refile(degrees, u, bound);
            // The minimum, taken HERE rather than in a pass of its own. `bound` is in a register
            // and `degrees[u]` has just been written from it, so the pass this replaces was one
            // scattered read per survivor per pivot to recover a value it had already had.
            // AMD_2 does the same inside its restore-degree-lists loop, `if (deg < mindeg)`.
            // Amd1 has always done it this way; Amd2, Amd2B and Amd3 did not.
            //
            // MEASURED AT ZERO, with the clique-weight fusion beside it: useful cycles unchanged
            // within half a percent in both families. A port and a simplification, not a speed
            // fix. See benchmarks/ordering/README.md (2026-08-09).
            minDegree = std::min(minDegree, bound);
        }

        // The whole array is invalidated in ONE ADDITION, where Amd3 walks the touched list and
        // zeroes each entry. After scan 1 no entry exceeds wflg + lemax, so advancing the tag by
        // lemax puts every one of them into the stale range. Amd.cpp's `wflg += lemax`.
        wflg += lemax;
    }

    return qg.order(pivots);
}


} // namespace Oblio
