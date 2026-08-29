#include "oblio/AmdCompacted.h"

#include "oblio/QuotientGraph.h"          // Buckets and TaggedScan
#include "oblio/QuotientGraphCompacted.h"
#include "oblio/Types.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// AmdCompacted.cpp - AmdFlat on a COMPACTED CLIQUE STORE: one pool with a free cursor and a
// compaction, in place of an arena that only grows.
//
// WHAT IT IS FOR, AND IT IS TWO THINGS.
//
// FIRST, IT IS THE ALIGNMENT VEHICLE FOR A DIFFERENTIAL. Comparing our ordering against a vendored
// one is only clean when the two hold their cliques the same way, or every difference is confounded
// with layout. What then remains is either LAYOUT or an IMPROVEMENT to carry back into our own
// ladder.
//
// SECOND, IT IS THE PREDICTABLE-SPACE VERSION OF AmdFlat. Given a machine you know whether A fits,
// but you cannot know whether L fits, nnz(L) depending on the ordering being computed. So a method
// that stays within `O(n + m)` carries a guarantee no amount of speed substitutes for: IF THE INPUT
// FITS, THE ANSWER IS REACHABLE. The compaction below is that guarantee bought deliberately, not
// frugality. Our arena is the right default for a known shape on a known machine solved repeatedly;
// this is the right one when whether an answer exists is the open question.
//
// ITS OBLIGATION. It must return AmdFlat's permutation exactly, which is what makes it an
// instrument rather than a second ordering.
//
// THREE THINGS THE STORAGE REQUIRES, and they are why this file's passes read as they do:
//
//   1  THE RUN IS INCIDENCE FIRST, adjacency behind it, where the flat driver is the reverse. A
//      consumed prefix is then a PREFIX, so what remains is still contiguous, and the new clique
//      goes in by a three-move rotation rather than by holding a vertex back and swapping
//      afterwards. It is also what makes detection one loop rather than two.
//   2  THE WALK IS IN POSITIONS, off a base hoisted once. The pool never reallocates, so a
//      compaction invalidates a block's OFFSET and never the array's base; a position survives
//      one and a pointer does not.
//   3  THE COMPACTOR RUNS PER ENTRY and carries the half-built clique, so `beginElimination` does
//      not have to reserve room for a worst-case reach of n before the walk.
//
// The list order and the reverse incidence walk are flags this layout cannot serve two values of,
// and it adds one pass the flat driver does not make: a slide in `absorbAggressively`, the
// adjacency having to follow the incidence part down when that part shrinks.
//
// `Buckets` and `TaggedScan` are the shared types, nothing about either depending on how cliques
// are stored, and the quotient graph is the shared `QuotientGraphCompacted` with a suffix on each
// method the two branches disagree about. This driver calls the `...Amd` half of each pair.

namespace Oblio {

std::size_t gAmdCompactions = 0;   // read by the compaction probes


std::vector<std::int32_t> orderAmdCompacted(const std::vector<std::size_t>&  colPtr,
                                            const std::vector<std::int32_t>& rowIdx) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraphCompacted qg(colPtr, rowIdx);

    qg.setLateMassElimination(true);    // see the note below

    std::vector<std::int32_t> pivots;               // the order over supervariables
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;                // a counter, not a scan of the flags

    std::uint32_t numLive = static_cast<std::uint32_t>(size);

    std::vector<std::uint32_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    // THE PREPASS, riding in the filing loop below, and the DENSE-ROW RULE with it. A degree-zero
    // vertex is numbered where it stands and not filed; a row above the threshold is SET ASIDE,
    // kept out of every reachable set by a zero weight and appended to the permutation at the end.
    // The threshold is fixed rather than exposed. AmdFlat must produce the same permutation and
    // `test_order` asserts it.
    const std::uint32_t dense = static_cast<std::uint32_t>(std::max<double>(
        16.0, 10.0 * std::sqrt(static_cast<double>(size))));
    std::vector<std::int32_t> denseRows;             // ascending by construction
    Buckets buckets(size);
    // THE MINIMUM IS TAKEN AS VERTICES ARE FILED, over the vertices actually filed. It seeds the
    // upward walk below and nothing else reads it, so it has only to be a LOWER BOUND on the first
    // live bucket. The empty and the dense rows are numbered or set aside rather than filed, so a
    // pass over every degree would answer for buckets that do not exist.
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

    std::vector<std::int32_t> hashHead(size, NIL);

    std::vector<std::int32_t> touchedCliques;
    std::vector<std::int32_t> deadCliques;

    std::vector<std::int32_t> work(size, 1);       // every clique alive and unseen
    std::int32_t workTag  = 2;                     // the tag
    std::int32_t stamp = 2;                     // detection's marks, above workTag
    std::int32_t maxCliqueWeight = 0;
    const std::int32_t wbig = std::numeric_limits<std::int32_t>::max() - static_cast<std::int32_t>(size);
    std::size_t numFlagSweeps = 0;              // how often the guard below actually fires

    const auto clearFlag = [&]() {
        if (workTag < 2 || workTag >= wbig) {
            for (std::int32_t x = 0; x < static_cast<std::int32_t>(size); ++x)
                if (work[x] != 0) work[x] = 1;
            workTag  = 2;
            stamp = 2;         // same array, same scale
            ++numFlagSweeps;
        }
    };

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        clearFlag();
        touchedCliques.clear();
        TaggedScan scan{&buckets, work, degrees, touchedCliques, workTag,
                        static_cast<std::int32_t>(size + 1)};
        qg.eliminateAmd(pivot, scan);
        pivots.push_back(pivot);

        buckets.unfile(pivot);

        const std::int32_t* newClique     = qg.clique(pivot);
        std::uint32_t       newCliqueSize = qg.cliqueSize(pivot);

        // |C[p]| weighted, and weighted rather than `newCliqueSize`: a degree counts original
        // vertices and a clique holds supervariables. Taken twice per pivot, the second time after
        // mass elimination has trimmed the clique.
        std::uint32_t newCliqueWeight = qg.cliqueWeight();
        degrees[pivot] = newCliqueWeight;       // what the scan below subtracts from

        maxCliqueWeight = std::max(maxCliqueWeight, static_cast<std::int32_t>(newCliqueWeight));

        deadCliques.clear();
        for (std::int32_t c : touchedCliques)
            if (work[c] == workTag) { deadCliques.push_back(c); work[c] = 0; }  // |C[c]-C[p]| == 0
        qg.absorbAggressively(deadCliques, newClique, newCliqueSize);

        const std::vector<std::int32_t>& merged = qg.massEliminate(pivot);
        numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
        numLive -= qg.weight(pivot);                // every original the pivot stands for
        for (std::int32_t u : merged) {
            degrees[u] = 0;                         // already out of the lists
        }

        newClique       = qg.clique(pivot);
        newCliqueSize   = qg.cliqueSize(pivot);
        newCliqueWeight = 0;
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk)
            newCliqueWeight += qg.weight(newClique[uk]);

        degrees[pivot] = newCliqueWeight;

        const std::uint32_t numLeft = numLive;

        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
            const std::int32_t u = newClique[uk];
            const std::int32_t* uIncidence     = qg.incidenceAmd(u);
            const std::uint32_t uIncidenceSize = qg.incidenceSize(u);

            std::size_t partialBound = static_cast<std::size_t>(work[u]);   // the adjacency half
            std::uint32_t uHashKey = buckets.hashKey(u);
            for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
                const std::int32_t c = uIncidence[ck];
                if (c != pivot) uHashKey += static_cast<std::uint32_t>(c);   // not the pivot
                if (c != pivot) partialBound += static_cast<std::size_t>(work[c] - workTag);
            }

            work[u] = static_cast<std::int32_t>(std::min<std::size_t>(partialBound, degrees[u]));

            if (!qg.eliminatedAmd(u)) {
                const std::int32_t uHashBucket = static_cast<std::int32_t>(
                                                     uHashKey % static_cast<std::uint32_t>(size));
                buckets.setChain(u, hashHead[uHashBucket]);
                hashHead[uHashBucket] = u;
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

        // THE CLIQUE IS TRIMMED AS THIS PASS WALKS IT: survivors are written back over the front
        // and what detection absorbed falls off the end. One store on a walk this pass makes
        // anyway.
        std::int32_t* cliqueOut = qg.clique(pivot);
        std::uint32_t kept      = 0;
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
            const std::int32_t u = newClique[uk];
            if (qg.eliminatedAmd(u)) continue;         // absorbed by the hash a moment ago
            const std::uint32_t uWeight = qg.weight(u);          // POST-merge
            std::size_t bound = static_cast<std::size_t>(work[u]) + newCliqueWeight - uWeight;
            work[u] = 1;
            bound = std::min<std::size_t>(bound, numLeft - uWeight);
            const std::uint32_t degreeBound = static_cast<std::uint32_t>(bound);
            degrees[u] = degreeBound;
            buckets.file(degreeBound, u);
            minDegree = std::min(minDegree, degreeBound);
            cliqueOut[kept++] = u;
        }
        qg.trimClique(pivot, kept);

        // PAST EVERY STAMP THIS STEP LAID DOWN, not merely past the scan's values: a stamp must
        // read as alive-and-unseen next step, which means strictly below the new tag.
        workTag = stamp + 1;                 // past every stamp this step laid down
    }

    // NOT AN ASSERT THAT THE LIVE COUNT IS ZERO. A clique dies when one of its members becomes a
    // pivot, and at the close of a run the last cliques can have had every member mass eliminated
    // into the pivot instead, leaving no one to absorb them, so a handful of entries legitimately
    // survive.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");

    // How often the pool actually needed compacting.
    gAmdCompactions  = qg.numCompactions();
    gPeakCliqueMembers = qg.numPeakCliqueMembers();
    // THE ROWS THE DENSE RULE SET ASIDE GO LAST, in index order. They were collected in an
    // ascending pass, and each stands only for itself, having been set aside before it could absorb
    // anything.
    pivots.insert(pivots.end(), denseRows.begin(), denseRows.end());
    return qg.orderAsMerged(pivots);
}

} // namespace Oblio
