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

// AmdCompacted.cpp - AmdFlat on AMD_2'S CLIQUE STORAGE: one pool with a free cursor and a
// compaction.
//
// WHAT IT IS FOR, AND IT IS TWO THINGS. Both are permanent; this file is not an experiment awaiting
// a verdict and its old stop condition is withdrawn. It is the amd counterpart of MmdChained, which
// says the same of itself.
//
// FIRST, IT IS THE ALIGNMENT VEHICLE FOR A DIFFERENTIAL. Comparing our ordering against `AMD_2` is
// only clean when the two hold their cliques the same way; otherwise every difference is confounded
// with layout. This file removes that confound. What then remains is either LAYOUT, whose price is
// measured below, or an IMPROVEMENT, which is carried back into our own ladder.
//
// THAT IS WHAT ACTUALLY HAPPENED HERE, and it is worth stating because the storage answer alone
// would read as a null result. With storage held equal the differential surfaced FIVE ARRAY FOLDS
// that have nothing to do with layout: the sign of the weight as the membership mark, the restore
// riding in the bound pass, `eliminated()` off a zero weight, mass elimination merging before it
// compacts, and supervariable detection stamping into `w`.
//
// FOUR OF THE FIVE ARE PORTED, 2026-08-17, and the fold set is closed. The sign, the merge before
// the compaction and the restore all landed in `QuotientGraph`, so every driver has them; the
// restore rides in massEliminate's walk over C[pivot] rather than in a bound pass, which is a
// different existing pass serving the same purpose. The detection stamp landed in `AmdFlat` and in
// `Amd2`, both of which already carried the tagged `w`. `Amd1` needed nothing, having no
// supervariable detection at all. SO THIS FILE'S TIME COLUMN IS NOW THE STORAGE PRICE ALONE, as
// MmdChained's already was, and it reads about 6 to 8 percent below AmdFlat up to 256 a side and 15
// percent below it from 400 up, reproduced across two runs.
//
// THE FIFTH CANNOT PORT, AND THE REASON IS Mmd1 RATHER THAN THE PREPASS. `eliminated()` answers
// from `mWeight[u] == 0` here and from `mMark[u] == GONE` in the shared class, and the two differ
// on exactly one thing: an ELIMINATED PIVOT. Nothing zeroes a pivot's weight, `merge` and
// `massEliminate` zeroing the absorbed vertex instead, and `orderAscending` needs the pivot's
// weight to the very end. So the zero means ABSORBED, not eliminated, and this file gets away with
// it only because no list it walks can contain an eliminated pivot.
//
// A shared version would have to pick between the two answers, and the only flag available is
// `mHasNumbered`, which is false for `Mmd1`: that driver never calls `number()`, so it would take
// the weight branch, and its refresh filters a `touched` list that CAN hold a vertex a later pivot
// in the same batch eliminated. It would read live and be refiled. Gating on "is this driver amd"
// instead would be a worse thing to introduce than the load it saves, since the swap is
// `mMark[u]` for `mWeight[u]`, both scattered, and the array stays for mmd either way. Checked and
// declined 2026-08-17.
//
// SECOND, IT IS THE PREDICTABLE-SPACE VERSION OF AmdFlat. From a conversation with Alex Pothen:
// given a machine you know whether A fits, but you cannot know whether L fits, nnz(L) depending on the
// ordering being computed. So a method that stays within `O(n + m)` carries a guarantee no amount
// of speed substitutes for: IF THE INPUT FITS, THE ANSWER IS REACHABLE. The compaction
// below is that guarantee bought deliberately, not frugality. Our arena is the right default for a
// known shape on a known machine solved repeatedly; this is the right one when whether an answer
// exists is the open question. See docs/DESIGN_DECISIONS.md (2026-08-16).
//
// THE STORAGE PRICE, measured on its own before the folds went in: 2.6 percent fewer data reads,
// 6 to 8 percent more D1 read misses, both constant across the ladder. A wash. Two compactions at
// every size from 50 to 1600 a side, constant, which is the assumption AMD_2's own complexity bound
// rests on.
//
// IT CARRIES A PRIVATE COPY OF QuotientGraph, named QuotientGraphCompacted, so the storage can
// be changed without touching the class every driver shares. THE COST IS EVERY SHARED FOLD
// LANDING TWICE, and
// it is accepted on purpose; `make digest` catches a copy that has stopped reproducing its original
// in half a second.
//
// ITS OBLIGATION. It must return AmdFlat's permutation exactly, which is `AMD_2`'s raw order, so
// its nnz(L) column in benchmarks/ordering must equal AMDraw's on every row and its fill column carries
// nothing. That obligation is what makes it an instrument rather than a second ordering, and
// `make digest` checks it across every driver in half a second.
//
// THE STORAGE IS `AMD_2`'S, IN THREE STEPS TAKEN 2026-08-18, and the last of them closed the last
// divergence in it. None moved a permutation.
//
//   1  THE RUN IS ITS WAY ROUND: incidence first, adjacency behind it, where production is the
//      reverse. This is the one that unblocks the others: a consumed prefix is then a PREFIX, so
//      `Pe [me] = p ; Len [me] -= knt1` is expressible and what remains is still contiguous. It
//      also lets the new clique go in by `AMD_2`'s three-move rotation instead of by holding a
//      vertex back and swapping afterwards, which takes a test out of the adjacency loop.
//   2  THE WALK IS IN POSITIONS, off a base hoisted once. The pool never reallocates, so a
//      compaction invalidates a block's OFFSET and never the array's base; a position therefore
//      survives one and a pointer does not. Same cost, and it is how `AMD_2` walks.
//   3  THE COMPACTOR RUNS PER ENTRY and carries the half-built clique, so `beginElimination` no
//      longer reserves room for a worst-case reach of n before the walk. That reservation was the
//      last divergence and it was one-directional: we compacted where `AMD_2` would not.
//
// TWO THINGS WENT WITH THEM, both flags this layout cannot serve two values of and neither ever
// set here: the list order and the reverse incidence walk. See the notes where each setter was.
// One thing was added: a slide in `absorb`, the adjacency having to follow the incidence part down
// when that part shrinks. That is a pass `AMD_2` does not make at all, and the flip priced it
// rather than caused it. See the note there.
//
// VERIFIED BY VARYING THE HEADROOM, which is a better check than the digest alone. At the shipped
// reserve the compactor never runs on grids, so the mid-walk path would go untested; cut to
// `sum(Len)` with no elbow room at all it fires on nearly every large pivot, 3 to 14 compactions
// from 3 to 200 a side, and every permutation is byte identical to the shipped-headroom one, clean
// under ASan and UBSan. That also confirms `AMD_2`'s own claim that it runs with no elbow room,
// only slowly.
//
// WHAT IS STILL NOT ALIGNED, so that the list is one list. The elbow room: ours is
// `nnz + nnz/5 + n + 1` with nnz counting the diagonal, against `AMD_2`'s `nzaat + nzaat/5 + n`,
// about 1.2n larger, which is why we now compact once where it compacts once and would diverge
// under pressure. And the pass structure outside the storage, which a reading of `AMD_2`'s main
// loop on 2026-08-18 found several differences in; docs/NEXT.md carries them.
//
// NOTHING HERE IS A COPY ANY MORE. This file held its own `QuotientGraphCompacted` and its own
// `Buckets`, generated from the production ones and then hand edited in the storage layer, until
// 2026-08-19. The buckets turned out to be `Oblio::Buckets` verbatim, so they are used directly;
// the graph became the shared `QuotientGraphCompacted`, which `MmdCompacted` also uses, with a
// suffix on each method the two vendored routines disagree about. This driver calls the `...Amd` half of
// each pair. Design notes for both types live with them and are authoritative there.

namespace Oblio {

std::size_t gAmdCompactions = 0;   // read by tmp/ probes; see the note at the return below

// `Buckets` AND `TaggedScan` ARE PRODUCTION'S, NOT COPIES. This file held a `Buckets`
// and a `TaggedScan` that were include/oblio/QuotientGraph.h's two types verbatim,
// differing only in the name and in the one field that named the other. Nothing about either
// depends on how cliques are stored: the buckets are three link arrays over vertices and the
// scan is a bundle of references the driver assembles, so the storage question this file
// exists to ask never reaches them. The copies went on 2026-08-19 and the shared types are
// used directly; only the quotient graph is a separate class.


// THE QUOTIENT GRAPH IS THE SHARED `QuotientGraphCompacted`, and this file no longer carries a
// copy of it. It held one until 2026-08-19, in an anonymous namespace so that `MmdCompacted` could
// hold another; the two then evolved apart for months, which is what a private copy is bad
// at. The class is now one body with a suffix on each method the two vendored routines
// disagree about, so this driver calls the `...Amd` half of each pair.
// See include/oblio/QuotientGraphCompacted.h.


std::vector<std::int32_t> orderAmdCompacted(const std::vector<std::size_t>&  colPtr,
                                            const std::vector<std::int32_t>& rowIdx) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraphCompacted qg(colPtr, rowIdx);

    qg.setLateMassElimination(true);    // and mass elimination becomes this driver's, below

    std::vector<std::int32_t> pivots;               // the order over supervariables
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;                // a counter, not a scan of the flags

    std::uint32_t numLive = static_cast<std::uint32_t>(size);

    std::vector<std::uint32_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    // THE EMPTY-ROW PREPASS, riding in the filing loop as `AMD_2` does and as AmdFlat now does. See
    // src/AmdFlat.cpp for why it exists and what it fixes; the two must produce the same
    // permutation and `test_order` asserts it.
    // AND THE DENSE-ROW RULE, the other half of `AMD_2`'s initialization pass. A row whose degree
    // exceeds `max (16, 10 * sqrt (n))` is SET ASIDE: not eliminated, not available, kept out of
    // every reachable set by a zero weight, and appended to the permutation at the end. `AMD_2`:
    //
    //     ndense++ ; Nv [i] = 0 ; Elen [i] = EMPTY ; nel++ ; Pe [i] = EMPTY ;
    //
    // and at the output assembly, "This is a dense unordered variable, with no parent. Place it
    // last in the output order", `Next [i] = nel++` over i ascending.
    //
    // WHY IT MATTERS HERE AND NOT ON GRIDS. A grid has no vertex anywhere near the threshold, so
    // nothing in the digest or the scaling ladders can see this rule at all. On real matrices it
    // is the difference between our order and `AMD_2`'s on most social and power-law graphs, and
    // it is also where our worst timings on that set came from: a hub of degree in the thousands
    // that nobody set aside sits in every reachable set it touches. Measured before the rule went
    // in, benchmarks/matrices `make amdorder`: GHS_indef/bloweybq 0.36 ms for `amd_order` against
    // 20.4 for ours, bloweybl 0.90 against 41.0, QY/case9 1.05 against 12.7.
    //
    // THE THRESHOLD IS FIXED AT THE VENDORED DEFAULT rather than exposed. `AMD_2` reads
    // `Control [AMD_DENSE]`, defaulting to 10.0, and has a whole control structure to read it
    // from; this driver has none, and inventing one to hold a single constant would be the wrong
    // trade while the constant is the thing being matched.
    const std::uint32_t dense = static_cast<std::uint32_t>(std::max<double>(
        16.0, 10.0 * std::sqrt(static_cast<double>(size))));
    std::vector<std::int32_t> denseRows;             // ascending by construction; see the tail
    Buckets buckets(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        if (degrees[u] == 0) {
            pivots.push_back(u);
            ++numEliminated;
            --numLive;
            continue;
        }
        if (degrees[u] > dense) {                    // a hub; see the note above
            qg.setAside(u);
            denseRows.push_back(u);
            ++numEliminated;
            --numLive;
            continue;
        }
        buckets.file(degrees[u], u);
    }
    std::uint32_t minDegree = *std::min_element(degrees.begin(), degrees.end());


    std::vector<std::int32_t> hashHead(size, NIL);

    std::vector<std::int32_t> touchedCliques;
    std::vector<std::int32_t> deadCliques;

    std::vector<std::int32_t> w(size, 1);       // every clique alive and unseen, Amd.cpp's W
    std::int32_t wflg  = 2;                     // the tag, Amd.cpp's wflg
    std::int32_t stamp = 2;                     // detection's marks, above wflg; see below
    std::int32_t lemax = 0;                     // the largest clique so far, Amd.cpp's lemax
    const std::int32_t wbig = std::numeric_limits<std::int32_t>::max() - static_cast<std::int32_t>(size);
    std::size_t numFlagSweeps = 0;              // how often the guard below actually fires

    const auto clearFlag = [&]() {
        if (wflg < 2 || wflg >= wbig) {
            for (std::int32_t x = 0; x < static_cast<std::int32_t>(size); ++x)
                if (w[x] != 0) w[x] = 1;
            wflg  = 2;
            stamp = 2;         // the detection marks live in the same array and the same scale
            ++numFlagSweeps;
        }
    };

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        clearFlag();
        touchedCliques.clear();
        TaggedScan scan{&buckets, w, degrees, touchedCliques, wflg,
                        static_cast<std::int32_t>(size + 1)};
        qg.eliminate(pivot, scan);
        pivots.push_back(pivot);

        buckets.unfile(pivot);

        const std::int32_t* pivotClique     = qg.clique(pivot);
        std::uint32_t       pivotCliqueSize = qg.cliqueSize(pivot);

        std::uint32_t degme = qg.cliqueWeight();
        degrees[pivot] = degme;                     // what the scan below subtracts from

        lemax = std::max(lemax, static_cast<std::int32_t>(degme));

        deadCliques.clear();
        for (std::int32_t c : touchedCliques)
            if (w[c] == wflg) { deadCliques.push_back(c); w[c] = 0; }   // |C[c] - C[p]| == 0
        qg.absorb(deadCliques, pivotClique, pivotCliqueSize);

        const std::vector<std::int32_t>& merged = qg.massEliminate(pivot);
        numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
        numLive -= qg.weight(pivot);                // every original the pivot stands for
        for (std::int32_t u : merged) {
            degrees[u] = 0;                         // already out of the lists; see eliminate()
        }

        pivotClique     = qg.clique(pivot);
        pivotCliqueSize = qg.cliqueSize(pivot);
        degme = 0;
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) degme += qg.weight(pivotClique[k]);

        degrees[pivot] = degme;

        const std::uint32_t numLeft = numLive;

        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u = pivotClique[k];
            const std::int32_t* incidence     = qg.incidenceAmd(u);
            const std::uint32_t incidenceSize = qg.incidenceSize(u);

            std::size_t deg = static_cast<std::size_t>(w[u]);   // the adjacency half
            std::uint32_t key = static_cast<std::uint32_t>(buckets.key(u));
            for (std::uint32_t i = 0; i < incidenceSize; ++i) {
                const std::int32_t c = incidence[i];
                if (c != pivot) key += static_cast<std::uint32_t>(c);   // me is not in the key
                if (c != pivot) deg += static_cast<std::size_t>(w[c] - wflg);
            }

            w[u] = static_cast<std::int32_t>(std::min<std::size_t>(deg, degrees[u]));

            if (!qg.eliminatedAmd(u)) {
                const std::int32_t hash = static_cast<std::int32_t>(
                                              key % static_cast<std::uint32_t>(size));
                buckets.setChain(u, hashHead[hash]);
                hashHead[hash] = u;
                buckets.setKey(u, hash);
            }
        }

        // THE STAMP BASE IS RAISED BEFORE DETECTION, NOT AFTER IT, 2026-08-18, and this is a
        // CORRECTNESS requirement rather than tidiness. `w` holds two kinds of value: the scan's
        // `w[c] = degree[c] + wflg - nvi`, which reaches as high as `wflg + lemax`, and
        // detection's stamps. A stamp must be ABOVE every scan value of the same step, or a clique
        // whose scan value happens to land on the current stamp reads as marked and two vertices
        // that are not duplicates compare equal.
        //
        // Amd.cpp does exactly this, `wflg += lemax ; wflg = clear_flag (...)` between scan 2 and
        // SUPERVARIABLE DETECTION, and then stamps with `wflg` upward. Ours used to raise the base
        // at the END of the step, which left this step's stamps starting at `wflg` while this
        // step's scan values ran up to `wflg + lemax`: the two ranges overlapped exactly.
        //
        // WHAT IT COST. On Grund/meg4, n = 5860, vertices 5779 and 5780 were merged at pivot 5080
        // although their lists differ in six of sixteen entries, because one entry's scan value
        // equalled the stamp. That single false merge moved 109 positions of the permutation and
        // cost 297 entries of fill, 51809 against `AMD_2`'s 51512. No grid ever triggered it: the
        // overlap needs a clique degree that lands on the right value, and it fired on one matrix
        // in 246.
        stamp = std::max(stamp, wflg + lemax);

        for (std::uint32_t kk = 0; kk < pivotCliqueSize; ++kk) {
            const std::int32_t seed = pivotClique[kk];
            if (qg.eliminatedAmd(seed)) continue;
            const std::int32_t hash = buckets.key(seed);
            const std::int32_t headOfBucket = hashHead[hash];
            if (headOfBucket == NIL) continue;      // an earlier member already emptied it
            hashHead[hash] = NIL;

            for (std::int32_t u = headOfBucket; u != NIL && buckets.chain(u) != NIL;
                 u = buckets.chain(u)) {
                if (qg.eliminatedAmd(u)) continue;

                // THE STAMP GOES INTO `w`, WHICH RETIRES mMark. Amd.cpp does exactly this,
                // `W [Iw [p]] = wflg` over the whole of i's list, variables and cliques alike,
                // because both live in one id space and one array can hold a mark for either.
                //
                // ONE LOOP OVER THE RUN, FROM THE SECOND ENTRY, AND NO LIVENESS TEST, which is
                // Amd.cpp's `for (p = Pe[i]+1 ; p <= Pe[i]+ln-1 ; p++)` exactly. Two things make
                // it one loop rather than two: the run is contiguous, so the incidence and
                // adjacency parts are one span, and the pivot sits at the front from the prune's
                // rotation, so skipping it is skipping index 0. Amd.cpp needs no liveness test
                // because its lists never hold dead entries and neither do ours after the prune,
                // with ONE exception: a vertex the hash absorbed EARLIER IN THIS SAME LOOP is
                // still listed by its neighbors. Amd.cpp stamps it like any other member and
                // compares stored lengths, so both sides of a comparison count it and the answer
                // is consistent. Testing liveness instead is also consistent, but it is a
                // different quantity and it costs a test per entry.
                //
                // IT HAS TO INTERLEAVE WITH THE TAG PROTOCOL, not clobber it. `stamp` starts above
                // the values this step's scan wrote, `wflg + lemax`, and rises by one per
                // candidate; `wflg` is then set past all of them at the end of the step. So next
                // step every stamped entry reads BELOW wflg, which the prune's `we >= wflg` test
                // treats as alive-and-unseen and rebuilds from the clique degree. That is the same
                // reading Amd.cpp relies on, and it is why stamping is safe here at all.
                //
                // The zeros survive: only entries of a live vertex's list are stamped, and a dead
                // clique is not in one.
                const std::int32_t  other      = ++stamp;
                const std::int32_t* runU       = qg.incidenceAmd(u);        // the run's first entry
                const std::uint32_t incidenceU = qg.incidenceSize(u);
                const std::uint32_t adjacencyU = qg.adjacencySize(u);
                const std::uint32_t runSizeU   = incidenceU + adjacencyU;
                for (std::uint32_t a = 1; a < runSizeU; ++a) w[runU[a]] = other;

                for (std::int32_t v = buckets.chain(u); v != NIL; v = buckets.chain(v)) {
                    if (qg.eliminatedAmd(v)) continue;

                    // THE LENGTHS REJECT BEFORE THE LIST IS TOUCHED, which is Amd.cpp's
                    // `ok = (Len [j] == ln) && (Elen [j] == eln)`. Two compares throw out most
                    // candidates for nothing, where counting live entries as we walk cannot decide
                    // until the walk is over.
                    if (qg.incidenceSize(v) != incidenceU) continue;
                    if (qg.adjacencySize(v) != adjacencyU) continue;

                    bool                same = true;
                    const std::int32_t* runV = qg.incidenceAmd(v);
                    for (std::uint32_t a = 1; a < runSizeU && same; ++a)
                        if (w[runV[a]] != other) same = false;
                    if (!same) continue;

                    const std::uint32_t weightV = qg.weight(v);

                    qg.merge(u, v);                 // v folded into u, left where it lies
                    degrees[u] -= weightV;
                    ++numEliminated;                // out of the count, not out of the graph
                }
            }
        }

        // THE SIGNS COME BACK HERE, which is Amd.cpp's `Nv [i] = nvi` under RESTORE DEGREE LISTS
        // and is the same pass. Every member of C[pivot] has read negative since the clique was
        // built, which is what let the prune answer three questions from one load; from this point
        // on they are ordinary live vertices again.
        // THE SIGNS COME BACK INSIDE THIS PASS, not before it. Every member of C[pivot] has read
        // negative since the clique was built, which is what let the prune answer three questions
        // from one load; the store that ends that rides on the load this pass makes anyway.
        // Amd.cpp does the same, `nvi = -Nv [i]` then `Nv [i] = nvi`, under RESTORE DEGREE LISTS.
        //
        // A vertex the hash absorbed is skipped and never restored, which is right: `merge` zeroed
        // its weight, and zero is the absorbed state whatever sign it arrived with.
        qg.restorePivotWeight(pivot);

        // AND THE CLIQUE IS TRIMMED AS THIS PASS WALKS IT, which is Amd.cpp's `Iw [p++] = i` under
        // RESTORE DEGREE LISTS: the survivors are written back over the front of the clique and
        // what detection absorbed falls off the end. One store on a walk this pass makes anyway.
        std::int32_t* cliqueOut = qg.clique(pivot);
        std::uint32_t kept      = 0;
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u = pivotClique[k];
            if (qg.eliminatedAmd(u)) continue;         // absorbed by the hash a moment ago
            const std::uint32_t weightU = qg.restoreWeight(u);   // POST-merge, and un-negated here
            std::size_t bound = static_cast<std::size_t>(w[u]) + degme - weightU;
            w[u] = 1;
            bound = std::min<std::size_t>(bound, numLeft - weightU);
            const std::uint32_t filed = static_cast<std::uint32_t>(bound);
            degrees[u] = filed;
            buckets.file(filed, u);
            minDegree = std::min(minDegree, filed);
            cliqueOut[kept++] = u;
        }
        qg.trimClique(pivot, kept);

        // PAST EVERY STAMP THIS STEP LAID DOWN, not merely past the scan's values. Amd.cpp
        // advances wflg through detection for the same reason: a stamp must read as alive-unseen
        // next step, which means strictly below the new tag.
        wflg = stamp + 1;                 // past every stamp this step laid down
    }

    // NOT AN ASSERT THAT THE LIVE COUNT IS ZERO. A clique dies when one of its members becomes a
    // pivot, and at the close of a run the last cliques can have had every member mass eliminated
    // into the pivot instead, leaving no one to absorb them, so a handful of entries legitimately
    // survive. `AmdFlat` asserts the same thing on the same grounds.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");

    // How often the pool actually needed compacting. `AMD_2` reports the same figure as
    // Info[AMD_NCMPA] and its complexity bound assumes it stays constant; this is where that
    // assumption can be checked for OUR storage rather than for its.
    gAmdCompactions  = qg.numCompactions();
    gPeakCliqueMembers = qg.numPeakCliqueMembers();   // see include/oblio/QuotientGraph.h
    // THE ROWS THE DENSE RULE SET ASIDE GO LAST, in index order, which is where `AMD_2`'s output
    // assembly puts them. They were compacted in an ascending pass, so appending the vector is
    // that order; each stands only for itself, having been set aside before it could absorb
    // anything, so `order` expands a chain of one.
    pivots.insert(pivots.end(), denseRows.begin(), denseRows.end());
    return qg.order(pivots);
}

} // namespace Oblio
