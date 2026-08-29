#include "oblio/AmdFlat.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

namespace Oblio {
namespace {

// The public forms are an overload pair rather than one function with a default argument: a default
// argument is not part of a function's type, and `orderAmdFlat` has to bind to a plain two-argument
// function pointer.
std::vector<std::int32_t> orderAmdFlatImpl(const std::vector<std::size_t>&  colPtr,
                                           const std::vector<std::int32_t>& rowIdx,
                                           std::size_t* numBornCliqueMembers) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);   // detection stamps into `work`, so no clique marks

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
    //     work[c] == 0            absorbed and gone
    //     0 < work[c] < workTag      alive, not seen this step; the value is stale
    //     work[c] >= workTag         seen this step, and work[c] - workTag is |C[c] - C[p]|
    //
    // The tag advances by `maxCliqueWeight` at the end of each step, and that is what makes the
    // stale range safe with no clearing pass: after the scan no entry exceeds
    // `workTag + maxCliqueWeight`, so one addition invalidates the whole array.
    //
    // SIGNED, and `std::int32_t` for it. `workTag - weight(u)` is negative whenever the tag is
    // still small and the weight is not, and the arithmetic only comes right again at `work[c] -
    // workTag`; an unsigned type wraps there and the bound comes out enormous.
    std::vector<std::int32_t> work(size, 1);       // every clique alive and unseen
    std::int32_t workTag  = 2;                     // the tag
    // SUPERVARIABLE DETECTION STAMPS INTO `work` TOO, so the two scales must interleave rather than
    // collide. `stamp` starts above the values this step's scan wrote, `workTag + maxCliqueWeight`,
    // and rises by one per candidate; `workTag` is then set past all of them at the end of the
    // step, so next step every stamped entry reads BELOW `workTag`, which is the alive-and-unseen
    // state. The zeros survive, only entries of a live vertex's list being stamped.
    std::int32_t stamp = 2;                     // detection's marks, above workTag
    std::int32_t maxCliqueWeight = 0;
    // `workTag + n` must not overflow.
    const std::int32_t wbig = std::numeric_limits<std::int32_t>::max() - static_cast<std::int32_t>(size);
    std::size_t numFlagSweeps = 0;              // how often the guard below actually fires

    // Reset the array and the tag when the tag can no longer be advanced safely. Every live clique
    // goes back to 1, the alive-and-unseen state, and the dead ones stay 0.
    const auto clearFlag = [&]() {
        if (workTag < 2 || workTag >= wbig) {
            for (std::int32_t x = 0; x < static_cast<std::int32_t>(size); ++x)
                if (work[x] != 0) work[x] = 1;
            workTag  = 2;
            stamp = 2;         // same array, same scale
            ++numFlagSweeps;
        }
    };

    // The half of each bound that does not involve the vertex's own weight rides in `work[u]`, free
    // for a live vertex for exactly the span required. The obligation is the reset at the end of
    // the bound pass below.

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        // Under `setLateMassElimination` this returns an empty list and C[pivot] is reach(pivot)
        // exactly; the merge happens below, once absorption has run. The first scan is folded into
        // the prune, so the eliminator accumulates |C[c] - C[p]| into the tagged `work` on the walk
        // it already makes.
        //
        // clearFlag RUNS FIRST because that scan is inside the eliminator, so the tag has to be
        // valid before it. `maxCliqueWeight` still advances after, needing the clique weight.
        clearFlag();
        touchedCliques.clear();
        TaggedScan scan{&buckets, work, degrees, touchedCliques, workTag,
                        static_cast<std::int32_t>(size + 1)};
        qg.eliminateAmd(pivot, scan);
        pivots.push_back(pivot);

        // The pivot leaves the lists.
        buckets.unfile(pivot);

        // ---- the bound, in place of an exact refresh -----------------------------------
        // Everything the new clique reached needs a new degree and nothing else can have changed,
        // so C[p] is both the refresh set and the domain the decomposition below is a statement
        // about. Read at the moment of use: the arena grows as cliques are formed, so a pointer
        // taken before the next elimination is the only safe one.
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

        // AGGRESSIVE ABSORPTION. `work[c] - workTag == 0` says C[c] lies wholly inside the new
        // clique, so it can never contribute again and its entries in the incidence lists are pure
        // cost. The quantity was computed for the bound anyway, so the test is free.
        deadCliques.clear();
        for (std::int32_t c : touchedCliques)
            if (work[c] == workTag) { deadCliques.push_back(c); work[c] = 0; }  // |C[c]-C[p]| == 0
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
            // The first term is already in `work[u]`, put there by the prune as `uAdjacencyWeight`
            // over exactly the sets it produced, and it is added once at the end rather than used
            // as a seed, so this accumulator means ONE thing for its whole life.
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
                if (c != pivot) otherCliqueBound += static_cast<std::size_t>(work[c] - workTag);
                if (c != pivot) uHashKey += static_cast<std::uint32_t>(c);   // not the pivot
            }

            // TERMS ONE AND THREE MEET HERE, `work[u]` carrying the adjacency weight in and the
            // sum of both back out. What the slot holds therefore changes at this line, from one
            // term to two.
            //
            // AND IT IS THE ONE PLACE THE WIDE ACCUMULATOR MEETS A NARROW DEGREE. The stored degree
            // is a full one from an earlier step and this is two terms of three, so the two are
            // comparable only once the pass below adds the clique weight and subtracts u's. The
            // minimum is taken WIDE and is what makes the result representable, being at most
            // `degrees[u]` and so at most n.
            const std::size_t twoTerms = static_cast<std::size_t>(work[u]) + otherCliqueBound;
            work[u] = static_cast<std::int32_t>(std::min<std::size_t>(twoTerms, degrees[u]));

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
        // THE STAMP BASE IS RAISED FIRST, and this is CORRECTNESS. `work` holds the scan's values,
        // up to `workTag + maxCliqueWeight`, and detection's stamps. A stamp at or below a scan
        // value makes that clique read as marked, so two vertices that are not duplicates compare
        // equal.
        stamp = std::max(stamp, workTag + maxCliqueWeight);

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
                const std::int32_t  other          = ++stamp;
                const std::uint32_t uAdjacencySize = qg.adjacencySize(u);
                const std::uint32_t uIncidenceSize = qg.incidenceSize(u);
                const std::int32_t* uSegment       = qg.incidenceAmd(u);   // the segment's start
                const std::uint32_t uSegmentSize   = uAdjacencySize + uIncidenceSize;
                for (std::uint32_t a = 1; a < uSegmentSize; ++a) work[uSegment[a]] = other;

                for (std::int32_t v = buckets.chain(u); v != NIL; v = buckets.chain(v)) {
                    if (qg.eliminatedAmd(v)) continue;

                    // The lengths reject before either list is touched.
                    if (qg.adjacencySize(v) != uAdjacencySize) continue;
                    if (qg.incidenceSize(v) != uIncidenceSize) continue;

                    // The exact test the hash only filters for, A[u] == A[v] and I[u] == I[v],
                    // read against u's stamp and short-circuiting on the first mismatch.
                    bool                same     = true;
                    const std::int32_t* vSegment = qg.incidenceAmd(v);
                    for (std::uint32_t a = 1; a < uSegmentSize && same; ++a)
                        if (work[vSegment[a]] != other) same = false;
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
        // detection absorbed falls off the end. Without it those vertices are visited by every
        // later walk of this clique.

        std::int32_t* cliqueOut = qg.clique(pivot);
        std::uint32_t kept      = 0;
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
            const std::int32_t u = newClique[uk];
            if (qg.eliminatedAmd(u)) continue;         // absorbed by the hash a moment ago
            const std::uint32_t uWeight = qg.weight(u);          // POST-merge
            // ONE OPERAND WIDENED, so the sum is formed in `std::size_t`; widening cannot be done
            // after the addition the way narrowing is done after the subtraction. `work[u]` and the
            // clique weight each reach n, so the sum reaches 2n.
            std::size_t bound = static_cast<std::size_t>(work[u]) + newCliqueWeight - uWeight;
            // THE SLOT GOES BACK TO ALIVE-AND-UNSEEN, the last read having just happened. Without
            // it a survivor later chosen as pivot would form a clique whose `work` already held a
            // bound. A vertex the hash eliminated is skipped and never reset, which is right: it is
            // dead and no clique is ever named after it.
            work[u] = 1;
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

        // The whole array is invalidated in ONE ADDITION rather than by walking the touched list
        // and zeroing each entry. After the scan no entry exceeds `workTag + maxCliqueWeight`, so
        // advancing past that puts every one of them into the stale range. PAST EVERY STAMP THIS
        // STEP LAID DOWN, not merely past the scan's values, since detection writes into this same
        // array above `workTag + maxCliqueWeight`.
        workTag = stamp + 1;                 // past every stamp this step laid down
    }

    // THE COUNTER CROSS-CHECKED AGAINST A RECOMPUTATION, which the driver can do exactly because it
    // holds the pivot list and a clique's owner is a pivot. NOT AN ASSERT THAT IT IS ZERO: at the
    // close of a run the last cliques can have had every member mass eliminated into the pivot,
    // leaving no one to absorb them, so a few entries legitimately survive.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");
    gPeakCliqueMembers = qg.numPeakCliqueMembers();
    if (numBornCliqueMembers != nullptr) *numBornCliqueMembers = qg.numBornCliqueMembers();
    // THE ROWS THE DENSE RULE SET ASIDE GO LAST, in index order. They were collected in an
    // ascending pass, and each stands only for itself, having been set aside before it could absorb
    // anything.
    pivots.insert(pivots.end(), denseRows.begin(), denseRows.end());
    return qg.orderAsMerged(pivots);
}


} // namespace

std::vector<std::int32_t> orderAmdFlat(const std::vector<std::size_t>&  colPtr,
                                       const std::vector<std::int32_t>& rowIdx) {
    return orderAmdFlatImpl(colPtr, rowIdx, nullptr);
}

// The same ordering, reporting every member ever put into a clique.
std::vector<std::int32_t> orderAmdFlat(const std::vector<std::size_t>&  colPtr,
                                       const std::vector<std::int32_t>& rowIdx,
                                       std::size_t& numBornCliqueMembers) {
    return orderAmdFlatImpl(colPtr, rowIdx, &numBornCliqueMembers);
}

} // namespace Oblio
