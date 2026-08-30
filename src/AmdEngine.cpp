#include "oblio/AmdEngine.h"

#include "oblio/Buckets.h"
#include "oblio/ElmOrder.h"
#include "oblio/QuotientGraphCompacted.h"
#include "oblio/QuotientGraphFlat.h"
#include "oblio/Types.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

// AmdEngine.cpp - the body of the amd ordering, written once and instantiated for each clique
// store. Both instantiations and both graph classes are in this one unit, so each is compiled with
// its store's bodies visible and inlined into the pivot loop.
//
// THE STORES DIFFER IN WHAT THEY CAN BE ASKED. Only a store with a bounded pool can run out and be
// compacted, so only that one answers `numCompactions`. The overload pair below is where the body
// stops having to know which store it holds: the template answers zero for a store that does not
// publish the figure, and the compacted store's own overload answers for it.

namespace Oblio {
namespace {

template<class QuotientGraph> std::size_t numCompactionsOf(const QuotientGraph&) { return 0; }

std::size_t numCompactionsOf(const QuotientGraphCompacted& qg) { return qg.numCompactions(); }

} // namespace


template<class QuotientGraph>
void AmdEngine<QuotientGraph>::compute(const std::vector<std::size_t>&  colPtr,
                                       const std::vector<std::int32_t>& rowIdx,
                                       ElmOrder& eo) const {
    eo = ElmOrder();
    if (colPtr.empty()) return;
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return;

    QuotientGraph qg(colPtr, rowIdx);   // detection stamps into `markAmd`, so no clique marks

    // Mass elimination runs in this driver rather than inside the eliminator; see below.
    qg.setLateMassElimination(true);    // see the note above

    std::vector<std::int32_t> pivots;               // the order over supervariables
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;                // a counter, not a scan of the flags

    // The live ORIGINAL vertices. NOT `size - numEliminated`: a hash merge eliminates a vertex
    // whose weight lives on in the one that absorbed it, so it leaves that count without leaving
    // the graph. The first cap below reads this.
    std::uint32_t numLive = static_cast<std::uint32_t>(size);

    // Two questions from one array. For a LIVE vertex this is the cached degree, exact at
    // construction and a bound afterwards; for a DEAD one it is the weighted size of the clique
    // that vertex's elimination formed. The two cannot overlap: a clique id is the id of the pivot
    // that made it, and that vertex is dead from the moment the clique exists.
    std::vector<std::uint32_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    // THE PREPASS, riding in the filing loop below. A degree-zero vertex has nothing to eliminate
    // and nothing to update, so it is numbered where it stands and not filed: no clique is formed
    // and no list is pruned. It does not call `number()`, which is for a vertex numbered while its
    // neighbors still name it; a degree-zero vertex is in nobody's adjacency. `numLive` loses it
    // too, or the `numLeft - weight(u)` cap runs one too large per empty row.
    //
    // AND THE DENSE-ROW RULE. A row above the threshold is SET ASIDE: not eliminated, not
    // available, kept out of every reachable set by a zero weight, and appended to the permutation
    // at the end. The threshold is fixed rather than exposed.
    const std::uint32_t dense = static_cast<std::uint32_t>(std::max<double>(
        16.0, 10.0 * std::sqrt(static_cast<double>(size))));
    std::vector<std::int32_t> denseRows;             // ascending by construction
    Buckets buckets(size);
    // THE MINIMUM IS TAKEN AS VERTICES ARE FILED, over the vertices actually filed. It seeds the
    // upward walk below and nothing else reads it, so it has only to be a LOWER BOUND on the first
    // live bucket. A pass over every degree would answer for buckets that do not exist, the empty
    // and the dense rows being numbered or set aside rather than filed.
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        if (degrees[u] == 0) {
            pivots.push_back(u);
            ++numEliminated;
            --numLive;
            continue;
        }
        if (degrees[u] > dense) {                    // a hub
            qg.setAside(u);
            denseRows.push_back(u);
            ++numEliminated;
            --numLive;
            continue;
        }
        buckets.file(degrees[u], u);
        minDegree = std::min(minDegree, degrees[u]);
    }


    // The hash groups, head and next arrays rather than a container per bucket. `hashHead` is
    // indexed by a hash and the chain rides in `Buckets`' own links, which are free for exactly the
    // span that needs them because `eliminateAmd` takes every member of C[pivot] out of the lists.
    // The key is reduced `% size`, so a hash lies in [0, n) and one array of that size serves.
    std::vector<std::int32_t> hashHead(size, NIL);

    // The cliques this step reached, and those aggressive absorption killed. Hoisted out of the
    // loop: only the entries a step writes are read, so there is no clearing pass.
    std::vector<std::int32_t> touchedCliques;
    std::vector<std::int32_t> deadCliques;


    // THE TAGGED CLIQUE ARRAY, three facts in one slot:
    //
    //     markAmd[c] == 0             absorbed and gone
    //     0 < markAmd[c] < tagAmd     alive, not seen this step; the value is stale
    //     markAmd[c] >= tagAmd        seen this step, and markAmd[c] - tagAmd is |C[c] - C[p]|
    //
    // The tag advances by `maxCliqueWeight` at the end of each step, and that is what makes the
    // stale range safe with no clearing pass: after the scan no entry exceeds
    // `tagAmd + maxCliqueWeight`, so one addition invalidates the whole array.
    //
    // SIGNED, and `std::int32_t` for it. `tagAmd - weight(u)` is negative whenever the tag is
    // still small and the weight is not, and the arithmetic only comes right again at `markAmd[c] -
    // tagAmd`; an unsigned type wraps there and the bound comes out enormous.
    std::vector<std::int32_t> markAmd(size, 1);   // every clique alive and unseen
    std::int32_t tagAmd = 2;                      // the tag
    // SUPERVARIABLE DETECTION STAMPS INTO `markAmd` TOO, so the two scales must interleave rather
    // than collide. `stampAmd` starts above the values this step's scan wrote, `tagAmd +
    // maxCliqueWeight`, and rises by one per candidate; `tagAmd` is then set past all of them at
    // the end of the step, so next step every stamped entry reads BELOW `tagAmd`, which is the
    // alive-and-unseen state. The zeros survive, only entries of a live vertex's list being
    // stamped.
    std::int32_t stampAmd = 2;                    // detection's marks, above tagAmd
    std::int32_t maxCliqueWeight = 0;
    // ONE n IS THE EXACT MARGIN, and the two guards below are what make it exact. Each raise adds
    // at most n to the value the guard before it left: `tagAmd + maxCliqueWeight + 1 <= tagAmd + n`
    // since a clique excludes its own pivot, and at most n stamps follow. So a value AT the ceiling
    // is admissible, `ceiling + n == INT32_MAX`, which is why both tests are `>` and not `>=`.
    const std::int32_t tagCeilingAmd =
        std::numeric_limits<std::int32_t>::max() - static_cast<std::int32_t>(size);
    std::size_t numTagResets = 0;               // how often the guard below actually fires

    // THE TAG GUARD, run first in every step. `AMD_2`'s `clear_flag` at Amd.cpp:1694. It resets the
    // array and the tag when the tag can no longer be advanced safely: every live clique goes back
    // to 1, the alive-and-unseen state, and the dead ones stay 0.
    //
    // `AMD_2` ALSO TESTS `wflg < 2` HERE AND WE DO NOT, because that test is its INITIALIZATION.
    // It calls `clear_flag(0, ...)` once at Amd.cpp:1350 to put `W` at 1 and `wflg` at 2, so the
    // low test is what performs the setup. We do that setup in the declarations above, and the tag
    // only ever climbs, so a value below 2 is unreachable and a test for one would be dead.
    const auto resetAtTag = [&]() {
        if (tagAmd > tagCeilingAmd) {
            for (std::int32_t k = 0; k < static_cast<std::int32_t>(size); ++k)
                if (markAmd[k] != 0) markAmd[k] = 1;
            tagAmd   = 2;
            stampAmd = 2;         // same array, same scale
            ++numTagResets;
        }
    };

    // THE STAMP GUARD, run in the middle of the step. `AMD_2`'s second `clear_flag`, at
    // Amd.cpp:1949. ONE CHECK CANNOT COVER BOTH RAISES: the first is bounded by `maxCliqueWeight`
    // and the second by the candidate count, each reaching n, and they belong to DIFFERENT cliques,
    // `maxCliqueWeight` being a maximum over all previous steps, so their sum is not bounded by n.
    // Measured worst climb is 1.21n, and with one guard the `stampAmd++` below overflows a tag that
    // entered the step at the ceiling. See docs/NEXT.md.
    //
    // It sweeps on the RAISED STAMP rather than on the tag, that being the value about to be
    // incremented.
    const auto resetAtStamp = [&]() {
        if (stampAmd > tagCeilingAmd) {
            for (std::int32_t k = 0; k < static_cast<std::int32_t>(size); ++k)
                if (markAmd[k] != 0) markAmd[k] = 1;
            tagAmd   = 2;
            stampAmd = 2;
            ++numTagResets;
        }
    };

    // The half of each bound that does not involve the vertex's own weight rides in `markAmd[u]`,
    // free for a live vertex for exactly the span required. The obligation is the reset at the end
    // of the bound pass below.

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        // Under `setLateMassElimination` this returns an empty list and C[pivot] is reach(pivot)
        // exactly; the merge happens below, once absorption has run. The first scan is folded into
        // the prune, so the eliminator accumulates |C[c] - C[p]| into the tagged `markAmd` on the
        // walk it already makes.
        //
        // THE TAG GUARD RUNS FIRST because that scan is inside the eliminator, so the tag has to be
        // valid before it. `maxCliqueWeight` still advances after, needing the clique weight.
        resetAtTag();
        touchedCliques.clear();
        TaggedScan scan{&buckets, markAmd, degrees, touchedCliques, tagAmd,
                        static_cast<std::int32_t>(size + 1)};
        qg.eliminateAmd(pivot, scan);
        pivots.push_back(pivot);

        // The pivot leaves the lists.
        buckets.unfile(pivot);

        // ---- the bound, in place of an exact refresh -----------------------------------
        // Everything the new clique reached needs a new degree and nothing else can have changed,
        // so C[p] is both the refresh set and the domain the decomposition below is a statement
        // about. Read at the moment of use: a clique's handle does not survive the next
        // elimination, whether because the store grew or because it was compacted underneath.
        const std::int32_t* newClique     = qg.clique(pivot);
        std::uint32_t       newCliqueSize = qg.cliqueSize(pivot);

        // |C[p]| weighted, off the eliminator. The scan below is deliberately over the UNTRIMMED
        // clique, and the second computation further down is not removable, running after mass
        // elimination has trimmed it. Weighted rather than `newCliqueSize`: a degree counts
        // original vertices and a clique holds supervariables.
        std::uint32_t newCliqueWeight = qg.cliqueWeight();
        degrees[pivot] = newCliqueWeight;       // what the scan below subtracts from

        // |C[c] - C[p]| is obtained by SUBTRACTION and never by looking at C[c]:
        //
        //     |C[c] - C[p]| = |C[c]| - sum of weight(u) over u in C[c] & C[p]
        //
        // `degrees[c]` supplies the first term and the members of C[p] the second, c being in I[u]
        // exactly when u is in C[c]. So the scan walks the incidence lists of the new clique's
        // members and pays sum |I[u]|, where walking member lists would pay sum |C[c]|.
        maxCliqueWeight = std::max(maxCliqueWeight, static_cast<std::int32_t>(newCliqueWeight));

        // AGGRESSIVE ABSORPTION. `markAmd[c] - tagAmd == 0` says C[c] lies wholly inside the new
        // clique, so it can never contribute again and its entries in the incidence lists are pure
        // cost. The quantity was computed for the bound anyway, so the test is free.
        deadCliques.clear();
        for (std::int32_t c : touchedCliques)
            // |C[c] - C[p]| == 0
            if (markAmd[c] == tagAmd) { deadCliques.push_back(c); markAmd[c] = 0; }
        qg.absorbAggressively(deadCliques, newClique, newCliqueSize);

        // MASS ELIMINATION, HERE rather than inside the eliminator, because absorption is what
        // makes the cheap structural test agree with the true one: a clique whose members all lie
        // inside C[p] contributes nothing to what u can reach, yet its presence in I[u] makes the
        // test fail. Asking before absorption declines merges that should be made.
        const std::vector<std::int32_t>& merged = qg.massEliminate(pivot);
        numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
        numLive -= qg.weight(pivot);                // every original the pivot stands for
        for (std::int32_t u : merged) {
            degrees[u] = 0;                         // already out of the lists
        }

        // The clique and its weight are re-read, both having moved: mass elimination trims C[p] of
        // the merged members and folds their weight into the pivot. Every survivor must see the
        // same final value.
        newClique       = qg.clique(pivot);
        newCliqueSize   = qg.cliqueSize(pivot);
        newCliqueWeight = 0;
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk)
            newCliqueWeight += qg.weight(newClique[uk]);

        // AND THE STORED CLIQUE DEGREE IS WRITTEN AGAIN, which is not a tidy-up. THIS write is the
        // durable one, holding the post-merge size that a later step reads as |C[pivot]|; the write
        // above the scan holds the pre-merge size, which is what the scan itself must subtract
        // from, running as it does over the untrimmed clique. Without this one a pivot that mass-
        // eliminated leaves a clique degree permanently too large by the merged weight, and every
        // later bound taken through that clique inherits it.
        degrees[pivot] = newCliqueWeight;

        const std::uint32_t numLeft = numLive;

        // The bound is formed in TWO halves. This pass computes the part that does not involve u's
        // own weight; the pass after the hash adds the clique weight and subtracts u's. The split
        // is load-bearing rather than a rearrangement: supervariable detection runs BETWEEN the two
        // and a merge grows u's weight, so subtracting it here would use the pre-merge value and
        // file a supervariable one bucket too high per vertex absorbed.
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
            const std::int32_t u = newClique[uk];
            // bound(u) = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]| over I[u] - {p}, against
            // the exact |( A[u] | C[c] for c in I[u] ) - {u}|. The first term needs no subtraction:
            // when a clique is formed every member has its adjacency pruned against it, and neither
            // set grows afterwards, so A[u] and C[p] are already disjoint. All the overcounting is
            // clique against clique, outside C[p].
            const std::int32_t* uIncidence     = qg.incidenceAmd(u);
            const std::uint32_t uIncidenceSize = qg.incidenceSize(u);

            // THE THIRD TERM OF THE BOUND, and the only inexact one. Each addend is the exact
            // weight |C[c] - C[p]|, but the C[c] can overlap EACH OTHER outside C[p], so a vertex
            // in two of them is counted twice. That over-count is the whole gap between the bound
            // and the true degree, which is why this is a `Bound` where the other two terms are
            // `Weight`s.
            //
            // WIDE, and it is the one accumulator in the ordering that no disjointness argument
            // bounds: each addend reaches n and there are O(n) of them, so the intermediate reaches
            // O(n^2). The minimum below is what brings it back into range.
            //
            // The first term is already in `markAmd[u]`, put there by the prune as
            // `uAdjacencyWeight` over exactly the sets it produced, and it is added once at the end
            // rather than used as a seed, so this accumulator means ONE thing for its whole life.
            std::size_t otherCliqueBound = 0;   // sum |C[c] - C[p]| over c in I[u] - {p}
            // The ADJACENCY HALF of the key, already reduced. The other half is accumulated below,
            // in the walk this pass makes anyway, and cannot move into the prune: absorption runs
            // between the two and compacts I[u], so the list to sum over does not exist there yet.
            std::uint32_t uHashKey = buckets.hashKey(u);
            // A SUM, because addition has no order and neither do the sets. NO STRIDE: multiplying
            // the incidence half by `size + 1` would be annihilated by the modulus at the filing
            // site, which is the same number, leaving the hash a function of the ADJACENCY ALONE.
            // The invariant the two lines have to hold together is that the modulus must not divide
            // the stride, and having no stride is the cheapest way to hold it.
            //
            // THE PIVOT IS EXCLUDED FROM BOTH, and the two `c != pivot` tests are one rule rather
            // than two. Excluding it is free either way: every member of C[p] carries the pivot, so
            // taking it would shift every key in this clique by the same constant and leave the
            // partition into buckets identical.
            for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
                const std::int32_t c = uIncidence[ck];
                if (c != pivot) otherCliqueBound += static_cast<std::size_t>(markAmd[c] - tagAmd);
                if (c != pivot) uHashKey += static_cast<std::uint32_t>(c);   // not the pivot
            }

            // TERMS ONE AND THREE MEET HERE, `markAmd[u]` carrying the adjacency weight in and
            // `degrees[u]` taking the sum of both back out. The adjacency half is SPENT at this
            // line, which is what leaves the tag array holding nothing but tags from here on.
            //
            // THE PARTIAL BOUND GOES INTO `degrees` AND NOT BACK INTO `markAmd`, which is `AMD_2`'s
            // arrangement: `Degree[i] = MIN(Degree[i], deg)` there, read back after detection as
            // `Degree[i] + degme - nvi`. It is not a preference. A sweep of the tag array has to be
            // possible between the two terms that raise the tag, and a live vertex's bound sitting
            // in that array is exactly what a sweep would destroy. See the entry in docs/NEXT.md.
            //
            // AND IT IS THE ONE PLACE THE WIDE ACCUMULATOR MEETS A NARROW DEGREE. The stored degree
            // is a full one from an earlier step and this is two terms of three, so the two are
            // comparable only once the pass below adds the clique weight and subtracts u's. The
            // minimum is taken WIDE and is what makes the result representable, being at most
            // `degrees[u]` and so at most n. The slot is read as the cap and written as the
            // destination in one statement, the read happening first.
            const std::size_t twoTerms = static_cast<std::size_t>(markAmd[u]) + otherCliqueBound;
            degrees[u] = static_cast<std::uint32_t>(std::min<std::size_t>(twoTerms, degrees[u]));

            // AND THE VERTEX IS FILED HERE, in the pass that completes its key. The skip for an
            // eliminated vertex moves onto the FILING alone: a bound computed for one is written
            // and never read. This loop walks C[p] FORWARD and pushes at the head, and that order
            // decides which of two indistinguishable vertices absorbs the other, so reversing it
            // moves the permutation.
            if (!qg.eliminatedAmd(u)) {
                // One reduction over a key that has wrapped in uint32 rather than been reduced per
                // term.
                const std::int32_t uHashBucket = static_cast<std::int32_t>(
                                                     uHashKey % static_cast<std::uint32_t>(size));
                buckets.setChain(u, hashHead[uHashBucket]);
                hashHead[uHashBucket] = u;
                // THE REDUCED FORM GOES BACK INTO THE SLOT, over the unreduced sum the prune left
                // there, which is what lets the detection pass below find a vertex's bucket from
                // the vertex and so walk C[pivot] instead of a list of the keys it used.
                buckets.setHashBucket(u, uHashBucket);
            }
        }

        // HASH SUPERVARIABLE DETECTION. Two members of C[pivot] indistinguishable from EACH OTHER,
        // which mass elimination cannot see, only ever comparing a member against the pivot. The
        // hash is a filter and never the decision, so a collision costs a comparison rather than a
        // wrong merge.
        //
        // DRIVEN BY C[pivot], and a bucket is EMPTIED the moment it is reached, so a later member
        // of the same bucket finds nothing and there is no clearing pass. Every bucket that was
        // filled is reached: only a principal member of C[pivot] is filed, and the survivor of a
        // merge is one too.
        //
        // THE STAMP BASE IS RAISED FIRST, and this is CORRECTNESS. `markAmd` holds the scan's
        // values, which occupy `[tagAmd, tagAmd + maxCliqueWeight]` INCLUSIVE, and detection's
        // stamps. A stamp at or below a scan value makes that clique read as marked, so two
        // vertices that are not duplicates compare equal. The base is therefore ONE PAST the top of
        // that block, and it is the FIRST STAMP rather than the value below it: `stampAmd` names
        // the next stamp to hand out, so every value from `tagAmd` to the next tag is either
        // written by the scan or handed out as a stamp, and none is skipped.
        stampAmd = std::max(stampAmd, tagAmd + maxCliqueWeight + 1);
        resetAtStamp();

        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
            const std::int32_t seed = newClique[uk];
            if (qg.eliminatedAmd(seed)) continue;
            const std::int32_t seedHashBucket = buckets.hashBucket(seed);
            const std::int32_t headOfBucket   = hashHead[seedHashBucket];
            if (headOfBucket == NIL) continue;      // an earlier member already emptied it
            hashHead[seedHashBucket] = NIL;

            // A vertex at the END of its chain has nothing after it to compare against, so a
            // SINGLETON BUCKET, which most are, costs no iteration at all: the condition skips
            // exactly the cases where the inner loop would be empty.
            for (std::int32_t u = headOfBucket; u != NIL && buckets.chain(u) != NIL;
                 u = buckets.chain(u)) {
                if (qg.eliminatedAmd(u)) continue;

                // THE STAMP IS HOISTED: u is stamped once and every candidate read against it. No
                // `v != u` guard is needed, both being members of C[pivot] and the prune having
                // dropped every neighbor that lies inside it. STORED LENGTHS, not live ones: a
                // vertex the hash absorbed earlier in this loop is still listed by its neighbors
                // and both sides count it, so the answer is consistent.
                //
                // ONE LOOP, forced by the layout: the run is contiguous, so I[u] and A[u] are one
                // span, and the prune's rotation puts the new clique at index 0. That entry is
                // shared by every member of C[pivot], so it can never discriminate and is skipped
                // by starting at 1.
                const std::int32_t  other          = stampAmd++;
                const std::uint32_t uAdjacencySize = qg.adjacencySize(u);
                const std::uint32_t uIncidenceSize = qg.incidenceSize(u);
                const std::int32_t* uSegment       = qg.incidenceAmd(u);   // the segment's start
                const std::uint32_t uSegmentSize   = uAdjacencySize + uIncidenceSize;
                for (std::uint32_t a = 1; a < uSegmentSize; ++a) markAmd[uSegment[a]] = other;

                for (std::int32_t v = buckets.chain(u); v != NIL; v = buckets.chain(v)) {
                    if (qg.eliminatedAmd(v)) continue;

                    // The lengths reject before either list is touched.
                    if (qg.adjacencySize(v) != uAdjacencySize) continue;
                    if (qg.incidenceSize(v) != uIncidenceSize) continue;

                    // The exact test the hash only filters for, A[u] == A[v] and I[u] == I[v],
                    // read against u's stampAmd and short-circuiting on the first mismatch.
                    bool                same     = true;
                    const std::int32_t* vSegment = qg.incidenceAmd(v);
                    for (std::uint32_t a = 1; a < uSegmentSize && same; ++a)
                        if (markAmd[vSegment[a]] != other) same = false;
                    if (!same) continue;

                    // The TARGET IS LIVE and never the pivot, which is the whole difference from
                    // mass elimination: that folds into a vertex which is leaving, this into one
                    // that stays. u's reach is unchanged, the two being adjacent and an external
                    // degree excluding u's own supervariable, and |C[p]| weighted is unchanged,
                    // v's weight having moved to u.
                    //
                    // WHAT MOVES IS u's WEIGHT, and NOTHING IS CORRECTED HERE. The fourth pass
                    // reads that weight post-merge, which is the whole reason the bound is finished
                    // after detection rather than before it.
                    qg.merge(u, v);                 // v folded into u, left where it lies
                    // NO BUCKET TRAFFIC. Both are already out of the lists and u's mPrev holds its
                    // hash bucket, so an unfile would read that as a link.
                    ++numEliminated;                // out of the count, not out of the graph
                }
            }
        }

        // THE FOURTH PASS, which finishes the bounds and files them. It runs after detection so
        // that a vertex which absorbed another subtracts the COMBINED weight. The clique weight and
        // `numLeft` were settled before the hash and are the same for everyone. AND THE CLIQUE IS
        // TRIMMED AS THIS PASS WALKS IT: survivors are written back over the front and what
        // detection absorbed falls off the end, one store on a walk this pass makes anyway.
        // Without it those vertices are visited by every later walk of this clique.

        std::int32_t* cliqueOut = qg.clique(pivot);
        std::uint32_t kept      = 0;
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
            const std::int32_t u = newClique[uk];
            if (qg.eliminatedAmd(u)) continue;         // absorbed by the hash a moment ago
            const std::uint32_t uWeight = qg.weight(u);          // POST-merge
            // ONE OPERAND WIDENED, so the sum is formed in `std::size_t`; widening cannot be done
            // after the addition the way narrowing is done after the subtraction. `degrees[u]` and
            // the clique weight each reach n, so the sum reaches 2n.
            std::size_t bound = static_cast<std::size_t>(degrees[u]) + newCliqueWeight - uWeight;
            // THE TAG SLOT GOES BACK TO ALIVE-AND-UNSEEN. It holds this step's detection stamp, and
            // before that the spent adjacency weight the prune left; either can exceed a small tag
            // and would then read as seen-this-step. A vertex the hash eliminated is skipped and
            // never reset, which is right: it is dead and no clique is ever named after it.
            markAmd[u] = 1;
            bound = std::min<std::size_t>(bound, numLeft - uWeight);
            // THE NARROWING POINT, after the cap, which is where the value is at most n.
            const std::uint32_t degreeBound = static_cast<std::uint32_t>(bound);
            // FILE, NOT REFILE. The vertex has been out of the lists since C[pivot] was formed and
            // its mPrev holds a hash key rather than a link, so the unfile inside refile would read
            // that key as one.
            degrees[u] = degreeBound;
            buckets.file(degreeBound, u);
            // The minimum, taken here rather than in a pass of its own.
            minDegree = std::min(minDegree, degreeBound);
            cliqueOut[kept++] = u;
        }
        qg.trimClique(pivot, kept);

        // The whole array is invalidated in ONE ASSIGNMENT rather than by walking the touched list
        // and zeroing each entry. `stampAmd` is already one past the last stamp handed out, and
        // every stamp was above every scan value, so taking it as the next tag puts the whole of
        // this step's range into the stale-and-therefore-alive band. NO `+ 1` IS NEEDED: the
        // post-increment left `stampAmd` on the first UNUSED value, which is exactly what a tag
        // must be.
        tagAmd = stampAmd;                     // one past the last stamp this step laid down
    }

    // THE COUNTER CROSS-CHECKED AGAINST A RECOMPUTATION, which the driver can do exactly because it
    // holds the pivot list and a clique's owner is a pivot. NOT AN ASSERT THAT IT IS ZERO: at the
    // close of a run the last cliques can have had every member mass eliminated into the pivot,
    // leaving no one to absorb them, so a few entries legitimately survive.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");
    // THE ROWS THE DENSE RULE SET ASIDE GO LAST, in index order. They were collected in an
    // ascending pass, and each stands only for itself, having been set aside before it could absorb
    // anything.
    pivots.insert(pivots.end(), denseRows.begin(), denseRows.end());
    eo.mOrder                = qg.orderAsMerged(pivots);
    eo.mNumPeakCliqueMembers = qg.numPeakCliqueMembers();
    eo.mNumBornCliqueMembers = qg.numBornCliqueMembers();
    eo.mNumCompactions       = numCompactionsOf(qg);
    eo.mNumTagResets         = numTagResets;
}


template class AmdEngine<QuotientGraphFlat>;
template class AmdEngine<QuotientGraphCompacted>;

} // namespace Oblio
