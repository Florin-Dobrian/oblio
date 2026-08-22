#include "oblio/AmdFlat.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

namespace Oblio {
namespace {

// THE BODY, WITH ONE OPTIONAL OUT-PARAMETER. The public forms are an overload pair rather than one
// function with a default argument, because a default argument is not part of a function's type: a
// defaulted parameter here would stop `orderAmdFlat` binding to the plain two-argument function
// pointer that benchmarks/ordering and the digest harness take its address as. An overload leaves
// that type intact.
std::vector<std::int32_t> orderAmd3Impl(const std::vector<std::size_t>&  colPtr,
                                        const std::vector<std::int32_t>& rowIdx,
                                        std::size_t* arenaEntries) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);   // no clique marks: detection stamps into `w`; see below

    // The two shared-class conventions this layer differs by. Both are off for every other driver
    // and neither changes which sets are computed, only which permutation comes out. See the
    // setters, and experiments/ordering/AmdFlat.md for ledger entries 2, 3 and 5.
    qg.setLateMassElimination(true);    // and mass elimination becomes this driver's, below

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
    // AMD_2's `Degree`, and it answers TWO questions from one array. For a live vertex it is the
    // cached degree, exact at construction and a bound afterwards; for a dead one it is the
    // WEIGHTED SIZE of the clique that vertex's elimination formed, which the prune subtracts
    // from. The two never overlap: a clique id IS the id of the pivot that made it, and that
    // vertex is dead from the moment the clique exists, so the store that retires it is the store
    // that starts the clique's life. A separate `cliqueDegree` was one of the seven n-arrays
    // AMD_2 allocates none of.
    std::vector<std::uint32_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    // THE EMPTY-ROW PREPASS, AND IT RIDES IN THE FILING LOOP, 2026-08-18. `AMD_2` numbers every
    // degree-zero vertex where it stands during initialization, in one ascending pass that files
    // everything else:
    //
    //     for (i = 0 ; i < n ; i++) { deg = Degree [i] ;
    //         if (deg == 0) { Elen [i] = FLIP (1) ; nel++ ; Pe [i] = EMPTY ; W [i] = 0 ; }
    //         else { ...file i... } }
    //
    // A vertex with no off-diagonal entry has nothing to eliminate and nothing to update, so it is
    // simply given the next number. NUMBERED, NOT ELIMINATED: no clique is formed and no list is
    // pruned, which is what `number` means here and what the mmd prepass already does.
    //
    // WITHOUT IT THE ORDER IS REVERSED AMONG THOSE VERTICES and nothing else moves. Filed at
    // degree zero and popped from the head, they came out LIFO, so a pure diagonal gave `4 3 2 1
    // 0` where `AMD_2` gives `0 1 2 3 4`. Identical fill, different permutation, and twelve
    // matrices in benchmarks/matrices are entirely of this kind: every `m = 0` row there differed
    // for this reason alone.
    //
    // `numLive` LOSES THEM TOO, which is `nel++` above and `nleft = n - nel` at the degree bound.
    // Leaving them in would make the `numLeft - weight(u)` cap one too large per empty row.
    //
    // AND IT DOES NOT CALL `number`, which the mmd prepass does. That function exists for a vertex
    // that is numbered while still being NAMED by its neighbors, which is the degree-1 case; it
    // marks the vertex GONE and sets `mHasNumbered`, and the flag puts a test in every walk for
    // the rest of the run. A degree-ZERO vertex is in nobody's adjacency, so no walk can reach it
    // and there is nothing to mark. `AMD_2` writes `W [i] = 0` here for the same non-reason and
    // keeps `Nv [i] = 1`. Not filing it is the whole of what has to happen.
    //
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

    // NO SEPARATE MEMBERSHIP SCRATCH AT ALL, 2026-08-17. The driver's own `mark` and `tag` went
    // first, folded into the quotient graph's; then the graph's mark stopped answering membership,
    // the sign of the weight having taken it over; and now supervariable detection stamps into
    // `w` instead of into a clique half, so this driver no longer asks for clique marks and
    // `mMark` is n rather than 2n. `AMD_2` has no mark array whatever, which is the shape this
    // has been converging on.
    //
    // The `cliqueStamp` offset that used to sit here, `qg.cliqueBase()`, went with the clique
    // half: there is no second id space to bias into any more.

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
    // ONE ARRAY WHERE THERE WERE THREE. `hashNext` held the chain and the half-built key; both
    // now ride in the buckets' own links, free for exactly the span that needs them because
    // `eliminate` takes every member of C[pivot] out of the lists. `usedKeys` recorded which
    // buckets to clear afterwards; AMD_2 clears a bucket by EMPTYING IT AS IT FINDS IT, so there
    // is nothing to record and no clearing pass.
    //
    // `hashHead` STAYS, and that is where the alignment stops. AMD_2 overlays its hash heads on
    // Head[], the degree-list heads, with FLIP marking which kind a slot holds, and parks a second
    // head in Last[Head[hval]] when both kinds are live at one index. We cannot: Buckets carries
    // genmmd's `bwd` encoding, in which a head's mPrev holds -(degree + 1) rather than being free,
    // so the slot AMD_2 parks its second head in is already spoken for. Two encodings, each
    // coherent, that do not compose. DESIGN_DECISIONS.md calls FLIP an anti-model, and n int32 to
    // keep the two list kinds apart is the cheaper side of that trade.
    std::vector<std::int32_t> hashHead(size + 1, NIL);

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

    // Amd.cpp's W ARRAY, and this file used to keep the same three facts in three places.
    //
    // The old shape kept THREE facts about a clique in three places: `mark[c]` says whether this step has
    // seen it, `outside[c]` carries |C[c] - C[p]|, and a clique is dead when it has been stripped
    // from the incidence lists. Its scan 1 therefore loads and stores `mark[c]` AND loads and
    // stores `outside[c]` for every incidence clique, two cache lines per clique, and clears
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
    // SUPERVARIABLE DETECTION STAMPS INTO `w`, which is `AMD_2`'s `W [Iw [p]] = wflg` over the
    // whole of i's list, variables and cliques alike, both living in one id space. It replaces a
    // 2n `mMark`, so this driver no longer asks for clique marks at all.
    //
    // IT HAS TO INTERLEAVE WITH THE TAG PROTOCOL, not clobber it. `stamp` starts above the values
    // this step's scan wrote, `wflg + lemax`, and rises by one per candidate; `wflg` is then set
    // past all of them at the end of the step, so next step every stamped entry reads BELOW wflg,
    // which the prune's `we >= wflg` treats as alive-and-unseen and rebuilds from the clique
    // degree. That is the reading Amd.cpp relies on and it is why stamping here is safe at all.
    //
    // The zeros survive: only entries of a LIVE vertex's list are stamped, and a dead clique is
    // not in one.
    std::int32_t stamp = 2;                     // detection's marks, above wflg; see above
    std::int32_t lemax = 0;                     // the largest clique so far, Amd.cpp's lemax
    // wflg + n must not overflow, which is the whole of Amd.cpp's wbig.
    const std::int32_t wbig = std::numeric_limits<std::int32_t>::max() - static_cast<std::int32_t>(size);
    std::size_t numFlagSweeps = 0;              // how often the guard below actually fires

    // Amd.cpp's clear_flag: reset the array and the tag when the tag can no longer be advanced
    // safely. Every live clique goes back to 1, which is the alive-and-unseen state, and the dead
    // ones stay 0. Called once per elimination, and almost never does anything.
    const auto clearFlag = [&]() {
        if (wflg < 2 || wflg >= wbig) {
            for (std::int32_t x = 0; x < static_cast<std::int32_t>(size); ++x)
                if (w[x] != 0) w[x] = 1;
            wflg  = 2;
            stamp = 2;         // the detection marks live in the same array and the same scale
            ++numFlagSweeps;
        }
    };

    // The half of each bound that does not involve the vertex's own weight, carried from the pass
    // that forms it across supervariable detection to the pass that finishes it. Amd.cpp keeps the
    // same quantity in Degree[i] between its scan 2 and its degree-list pass. See ledger entry 4.
    // `partial` IS GONE, and with it the last driver array that had no counterpart in AMD_2. It
    // carried the half of each bound that does not involve the vertex's own weight, from the pass
    // that forms it across supervariable detection to the pass that finishes it. Amd.cpp keeps the
    // same quantity in `Degree[i]` between its scan 2 and its degree-list pass, and we cannot:
    // ledger entry 4 splits our bound so `nvi` is read AFTER the merge, and the min cap in the
    // middle still needs the old `degrees[u]`, so that slot is occupied.
    //
    // It rides in `w[u]` instead, free for a live vertex for exactly the span required; see the
    // store in QuotientGraph's tagged prune. The obligation is the reset at the end of the bound
    // pass below.

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket
        const std::int32_t pivot = buckets.head(minDegree);   // whatever was filed last

        // Under setLateMassElimination this returns an empty list and C[pivot] is reach(pivot)
        // exactly. The merge happens below, once the absorption has run.
        // THE FIRST SCAN IS FOLDED INTO THE PRUNE, landed 2026-08-10. This driver walked I[u]
        // three times per pivot, in the prune, in scan 1 and in the bound, and A[u] twice, in the
        // prune and in the bound; `AMD_2` walks I[u] twice and A[u] once, and the per-pass
        // inventory in benchmarks/ordering/README.md put that difference at 96 percent of what was
        // left of the cubic gap. The eliminator accumulates |C[c] - C[p]| into the tagged w array
        // on the walk it is already making, so scan 1 goes entirely and the bound's adjacency loop
        // with it, leaving the vendored routine's counts exactly.
        //
        // Worth 10 to 16 percent on cubic grids, measured with both codes down the same harness
        // path, and 0 to 8 percent in 2D. Over eight runs this driver reads 0.83 to 0.89 ms at 16
        // cubed where the vendored routine reads 0.74 to 0.86, so the two overlap there, and about
        // 1.2x at 32 a side. The same fusion exists in Amd1B and Amd2B and measured zero there,
        // five percent slower in Amd1B's case, which is why it took a re-run to find: that reading
        // was 2D, at one size, before ledger entry 8.
        //
        // clearFlag RUNS FIRST, where AmdFlat calls it after the elimination: the scan is inside
        // the eliminator now, so the tag has to be valid before it, and `Amd.cpp` calls clear_flag
        // before its own scan 1 for the same reason. `lemax` still advances after, needing degme,
        // and it is consumed only at the end of the step.
        clearFlag();
        touchedCliques.clear();
        // NO ARRAYS OF ITS OWN, which is the other half of the change and was worth more in 2D
        // than the fold itself. Two values have to cross from the prune to the bound, and the
        // first version carried them in two fresh vectors of size n: that measured 3 to 9 percent
        // faster on cubes and 12 percent SLOWER in 2D from 200 a side up, which is the footprint
        // trade REPORT.md names and the same one that sank the 2026-08-08 key fusion. Both fit in
        // arrays this driver already has and that are dead at this point: `partial[u]` is not
        // written until the end of the bound pass, and `hashNext[u]` holds nothing until the
        // vertex is filed, which happens in that same pass after the key has been read. With the
        // arrays gone the 2D penalty went with them.
        TaggedScan scan{&buckets, w, degrees, touchedCliques, wflg,
                        static_cast<std::int32_t>(size + 1)};
        qg.eliminate(pivot, scan);
        pivots.push_back(pivot);

        // The pivot leaves the lists. The zeroing that used to follow is gone: under the fold
        // above `degrees[pivot]` is the slot the new clique's weight is written into a few lines
        // down, so it was a store nobody read. The old comment warned to unfile before zeroing
        // because the bucket index came from the degree; Buckets reads it out of mPrev, so that
        // ordering was already vestigial.
        buckets.unfile(pivot);

        // ---- the bound, in place of an exact refresh -----------------------------------
        // Everything the new clique reached needs a new degree, and nothing else can have
        // changed. C[p] is therefore both the refresh set and the domain the decomposition is a
        // statement about, which is the coincidence the whole placement rests on.
        // Read at the moment of use: the arena holding it grows as cliques are formed, so a
        // pointer taken before the next elimination is the only one that is safe, and this is
        // that window.
        const std::int32_t* pivotClique     = qg.clique(pivot);
        std::uint32_t       pivotCliqueSize = qg.cliqueSize(pivot);

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
        // accumulates `degme += nvi` while building the clique and this is that. The second
        // computation below is NOT removable: it runs after mass elimination has trimmed the
        // clique, which is ledger entry 7, and this one is deliberately over the untrimmed one.
        std::uint32_t degme = qg.cliqueWeight();
        degrees[pivot] = degme;                     // what the scan below subtracts from

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
        // touched clique pays sum |C[c]|. Measured on a 100x100 grid, 74281 cliques against
        // 272646, which is most of the reason this branch used to run three times slower than the
        // vendored routine. `Amd.cpp` does the same thing at `we = Degree[e] + wnvi`, then
        // `we -= nvi`, and it is the amd2 layer's pass 3.
        lemax = std::max(lemax, static_cast<std::int32_t>(degme));

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
        // test fail. Amd.cpp says so itself, making the same test in its scan 2 over a clique
        // list absorption has already compacted: with aggressive absorption, `deg == 0` is
        // identical to `Elen[i] == 1 && p3 == pn`. Asking first, as every other driver here does,
        // declines merges the vendored routine makes. experiments/ordering/AmdFlat.md, ledger
        // entry 3.
        const std::vector<std::int32_t>& merged = qg.massEliminate(pivot);
        numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
        numLive -= qg.weight(pivot);                // every original the pivot stands for
        for (std::int32_t u : merged) {
            degrees[u] = 0;                         // already out of the lists; see eliminate()
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
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) degme += qg.weight(pivotClique[k]);

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
        // was. See experiments/ordering/AmdFlat.md.
        degrees[pivot] = degme;

        const std::uint32_t numLeft = numLive;

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
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u = pivotClique[k];
            // bound(u) = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]| over I[u] - {p},
            // against the exact |( A[u] | C[c] for c in I[u] ) - {u}|. The first term needs no
            // subtraction: when a clique is formed every member has its explicit adjacency
            // pruned against it, and neither set grows afterwards, so A[u] and C[p] are already
            // disjoint. All the overcounting is therefore clique against clique, outside C[p],
            // which is the smallest place it could have been put.
            // Every length hoisted out of its condition; see the note in Amd1.
            const std::int32_t* incidence     = qg.incidenceAmd(u);
            const std::uint32_t incidenceSize = qg.incidenceSize(u);

            // THE ADJACENCY TERM AND THE WHOLE KEY ARE ALREADY IN HAND, accumulated by the
            // prune over exactly the sets it produced. AmdFlat walks A[u] here for the bound's
            // explicit part and for the key's explicit half, and walks I[u] for the key's other
            // half; this reads both instead, which is one whole walk of A[u] removed and the key
            // out of the loop below entirely. The key's rule that an eliminated neighbor is
            // skipped is preserved by the prune's own filter: massEliminate and number both set
            // mLiveMerges, so from the first merge onward the prune drops eliminated vertices,
            // and before it there are none to drop.
            // WIDE: `deg` accumulates `w[c] - wflg` over I[u], each term up to n and O(n) of
            // them, so the intermediate reaches O(n^2). It is the one accumulator in the ordering
            // that no disjointness argument bounds, and the minimum below is what brings it back
            // into range. `partial[u]` seeds it and is itself at most n.
            std::size_t deg = static_cast<std::size_t>(w[u]);   // the adjacency half
            // The ADJACENCY HALF of the key, already reduced. The other half is accumulated below,
            // in the walk this pass makes anyway, and cannot move into the prune: absorption runs
            // between the two and compacts I[u], so the list to sum over does not exist there yet.
            std::uint32_t key = static_cast<std::uint32_t>(buckets.key(u));
            // THE HASH KEY IS ACCUMULATED HERE, in the walks the bound is already making, which
            // is what Amd.cpp does with `hval += e` and `hval += j` inside its scan 2. It had a
            // pass of its own until 2026-08-10, walking A[u] and I[u] a second time; see the
            // filing site below for what that cost and why an earlier attempt at this failed.
            //
            // A SUM, because addition has no order and neither do the sets: sorting to build a key
            // would be a log factor for nothing.
            // ONE SUM, WITH NO STRIDE, and that is ledger entry 8. The incidence half was added as
            // `(c + 1) * (size + 1)` until 2026-08-09, so that a vertex and a clique of the same
            // index could not cancel. True of the KEY and false of the BUCKET: the modulus at the
            // filing site is the same number as the stride, so the incidence term was annihilated
            // exactly and the hash came out a function of the ADJACENCY ALONE. As the elimination
            // proceeds A[u] empties and everything a vertex reaches becomes cliques, so the
            // surviving key carried less and less, and cubic grids reach that state sooner than
            // square ones.
            //
            // Measured, for the SAME MERGES: 19.0 pairs tested per pivot at 140 a side against the
            // vendored routine's 0.33, and 155.3 at 26 cubed against its 0.48. Amd.cpp accumulates
            // `hval += e` and `hval += j` into one running value and takes it mod n, letting a
            // vertex and a clique collide on purpose, because the hash is a FILTER and never the
            // decision: a collision costs one exact comparison and cannot produce a wrong merge.
            // The invariant the two lines have to hold TOGETHER is that the modulus must not
            // divide the stride, and having no stride is the cheapest way to hold it.
            //
            // EVERY ENTRY OF I[u] IS TAKEN, the pivot included, where the bound skips it. The two
            // rules differ and both are Amd.cpp's: its scan 2 accumulates `hval += e` over the
            // clique list it is compacting, which by then holds the new clique, while the degree
            // term for the new clique is added in a later pass. Fusing the walks does not fuse
            // the rules.
            for (std::uint32_t i = 0; i < incidenceSize; ++i) {
                const std::int32_t c = incidence[i];
                if (c != pivot) key += static_cast<std::uint32_t>(c);   // me is not in the key
                if (c != pivot) deg += static_cast<std::size_t>(w[c] - wflg);
            }

            // Amd.cpp's `Degree[i] = MIN (Degree[i], deg)`. The stored degree is a full one from
            // an earlier step and this is a partial, so the two are comparable only once the pass
            // below adds `degme - weight(u)` to whichever won. That is why the minimum is taken
            // here and the common term added there rather than the other way round. The second
            // cap, `numLeft - weight(u)`, is also that pass's, for the same weight reason.
            // THE ONE PLACE THE WIDE ACCUMULATOR MEETS A NARROW DEGREE. `deg` sums `w[c] - wflg`
            // over I[u], each term up to n and O(n) of them, so its intermediate reaches O(n^2)
            // and it stays `std::size_t`; `degrees[u]` is one dimensional and at most n. The
            // minimum is taken WIDE and is what makes the result representable, being at most
            // `degrees[u]` and so at most n.
            w[u] = static_cast<std::int32_t>(std::min<std::size_t>(deg, degrees[u]));

            // AND THE VERTEX IS FILED HERE, WHICH IS WHY THIS FUSION NEEDS NO ARRAY. It was tried
            // on 2026-08-08 and reverted: that version carried the key in a vector of size n and
            // measured nothing at 140 a side and minus two percent at 400, which REPORT.md had
            // already named as the footprint trade, the same stream that made Amd1B slower at
            // large n after being faster at small. Filing at the point the key completes stores
            // nothing extra, since hashNext is size n either way and hashHead is already
            // allocated. THE FAILURE WAS THE ARRAY, NOT THE FUSION.
            //
            // It was also measured while ledger entry 8 was live, when the exact comparison ran
            // 19.0 pairs per pivot against the vendored routine's 0.33, so the pass it shortens
            // was not the one the profile was standing on.
            //
            // Measured on alpamayo in a scratch Amd3B, 2026-08-10: about 4 to 7 percent faster at
            // six consecutive square grids from 64 to 400 a side, and 5 to 14 percent on cubic
            // grids from 12 to 32, with nnz(L) identical at every size. APPROXIMATE, both ranges:
            // the variant was timed as a free function and this driver through OrderEngine, which
            // also builds a Permutation, a bias of up to 2.4 percent. It removes a sweep over C[p]
            // and two walks, 26.70 of 149.96 element visits per pivot at 140 a side and 66.77 of
            // 352.57 at 26 cubed. See benchmarks/ordering/README.md and AmdFlat.md.
            //
            // THE GUARD AND THE DIRECTION ARE THE TWO THINGS THAT HAD TO COME ACROSS. The key pass
            // skipped an eliminated member and this bound pass does not, so the skip moves onto
            // the FILING alone: a bound computed for an eliminated vertex is written to partial[u]
            // and never read, exactly as before. And this loop walks C[p] FORWARD and pushes at
            // the head, which is the direction the key pass walked, so the chain comes out in the
            // same order. That order decides which of two indistinguishable vertices absorbs the
            // other, so a reversal would be a tie-break change wearing a schedule change's
            // clothes; it moved the permutation on four test graphs when this bucket was first
            // turned around.
            //
            // AND THIS DRIVER IS THE ONLY ONE THAT CAN DO IT, which is ledger entry 4's doing.
            // Amd2 and Amd2B form the bound in ONE pass and call `buckets.refile` inside it, so
            // the direction of their bound loop is ALREADY a tie-break input, deciding which
            // vertex sits at a degree bucket's head. Their key pass walks C[p] backward against
            // that forward bound, and head insertion into both structures wants opposite
            // directions, so one walk cannot serve both. Measured: fusing there changes the
            // permutation on all ten grids tried, so it is an ORDERING change there rather than a
            // schedule one. Tail insertion would preserve the order and is untried, needing a
            // hashTail array of size n, which is the footprint that made the 2026-08-08 version of
            // this fusion measure nothing. Entry 4 split this driver's bound in two and
            // moved the refile below the hash, for the post-merge weight, and that is what leaves
            // this loop free of tie-break duty. Second thing that split has bought by accident.
            if (!qg.eliminatedAmd(u)) {
                // Amd.cpp's `hval = hval % n`, one reduction over a key that has wrapped in uint32
                // rather than been reduced per term. See the prune for the other three halves of
                // this arithmetic.
                const std::int32_t hash = static_cast<std::int32_t>(
                                              key % static_cast<std::uint32_t>(size));
                buckets.setChain(u, hashHead[hash]);
                hashHead[hash] = u;
                // THE REDUCED KEY STAYS IN THE SLOT, which is `Last [i] = hval` in Amd.cpp. It is
                // what lets the detection pass below find a vertex's bucket from the vertex, and so
                // walk C[pivot] instead of a list of the keys it used.
                buckets.setKey(u, hash);
            }
        }

        // HASH SUPERVARIABLE DETECTION. Vertices indistinguishable from EACH OTHER, which the
        // pivot test cannot see: mass elimination only ever finds a vertex indistinguishable from
        // the pivot, and two vertices can become interchangeable with one another without either
        // being interchangeable with it. Hash first so the exact comparison runs only within a
        // group; the hash is a filter and never the decision, so a collision costs a comparison
        // rather than a wrong merge.
        //
        // THE BUCKETS ARE ALREADY FILLED, by the bound pass above, which accumulates each key in
        // the walks it is making anyway. Filling them was a sweep over C[p] and a second walk of
        // A[u] and I[u] until 2026-08-10. The chain is still built at the head from a forward
        // sweep, so it still comes out reversed against C[p], which is what the pair loop below
        // depends on.

        // DRIVEN BY C[pivot], NOT BY A LIST OF KEYS, and the bucket is EMPTIED the moment it is
        // reached. That is Amd.cpp's supervariable detection exactly: it walks Lme, reads each
        // member's key out of Last[i], takes the bucket and clears the head in one step, so a
        // later member of the same bucket finds nothing and the clearing pass does not exist.
        //
        // Every bucket that was filled is reached: filing only happens for a principal member of
        // C[pivot], and the survivor of every merge inside a bucket is principal and is itself a
        // member of C[pivot]. So no head is left dirty for the next step.
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

            // THE CONDITION IS AMD_2'S, `while (i != EMPTY && Next [i] != EMPTY)`. A vertex at the
            // END of its chain has nothing after it to compare against, so the body did nothing
            // and returned; testing it here means a SINGLETON BUCKET, which most are, costs no
            // iteration at all. Measured at 6.42 chain steps per pivot against AMD_2's 0.17 on a
            // 400 square before this, and 0.33 after.
            for (std::int32_t u = headOfBucket; u != NIL && buckets.chain(u) != NIL;
                 u = buckets.chain(u)) {
                if (qg.eliminatedAmd(u)) continue;
                // AMD_2 enters this loop only for a bucket member with a successor,
                // `while (i != EMPTY && Next [i] != EMPTY)`, and the guard is the other half of
                // the hoisted stamp below: without it a member with nothing after it pays a full
                // list of random writes for a pair that will never be tested. It cannot change
                // the answer, since the inner loop is empty in exactly the cases it skips. Inert
                // while the buckets were enormous, which is why it was never missed, and most of
                // the pass once entry 8 made them singletons. Amd2 and Amd2B need no counterpart:
                // they stamp INSIDE the pair loop, so a member with no successor already costs
                // them nothing.

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
                // iterations either way. We were stamping 2.4 times more cliques than we compared.
                //
                // WHY IT CAN BE HOISTED, and why it took a day to see. The stamping carried
                // `w != u` and the walk carried `w == v`, exclusions that look pair-dependent and
                // so look to pin the stamp inside the loop. They are vestigial: u and v are both
                // members of C[pivot], and the prune drops every neighbour lying inside the new
                // clique, `if (mMark[v] == inClique) continue`, so A[u] cannot contain v and A[v]
                // cannot contain u. Nothing to exclude. That is exactly why Amd.cpp's stamp has no
                // such guard either. This was landed as a separate driver first, so that
                // `AMD3C == AmdFlat` could say the reasoning held rather than merely sounding
                // right, and folded in once it did.
                //
                // Roles swapped with the hoist: u is stamped and v is walked, where before it was
                // v stamped and u walked. The test is symmetric so the outcome does not move, and the SURVIVOR
                // does not either, u being the outer vertex in both and merge(u, v) folding v into
                // it. Amd.cpp merges j into i the same way round.
                // THE STAMP GOES INTO `w`, which is what retires the clique half of mMark. Amd.cpp
                // does exactly this, `W [Iw [p]] = wflg` over the whole of i's list, variables and
                // cliques alike, both living in one id space so one array holds a mark for either.
                // See the declaration of `stamp` for the interleave with the tag protocol.
                // NO LIVENESS TEST, AND THE LENGTHS REJECT BEFORE THE LIST IS TOUCHED, which is
                // Amd.cpp's `for (p = Pe[i]+1 ; p <= Pe[i]+ln-1 ; p++) W [Iw[p]] = wflg` and its
                // `ok = (Len [j] == ln) && (Elen [j] == eln)`. Its lists never hold dead entries
                // and neither do ours after the prune, with ONE exception: a vertex the hash
                // absorbed EARLIER IN THIS SAME LOOP is still listed by its neighbors. Amd.cpp
                // stamps it like any other member and compares stored lengths, so both sides of a
                // comparison count it and the answer is consistent. Counting live entries instead
                // is also consistent, but it is a different quantity, it costs a test per entry,
                // and it cannot reject a candidate until its whole list has been walked.
                //
                // TWO LOOPS RATHER THAN Amd.cpp's ONE, and that is the layout rather than a
                // choice. Its run is cliques then variables with the new clique at the front, so
                // the whole list minus the first entry is one span. Ours is A[u] then I[u] with the
                // new clique at the front of I[u], so the entry to skip is in the middle. `Amd3B`,
                // which carries `AMD_2`'s order, does get the single loop.
                const std::int32_t  other      = ++stamp;
                const std::uint32_t adjacencyU = qg.adjacencySize(u);
                const std::uint32_t incidenceU = qg.incidenceSize(u);
                const std::int32_t* runU       = qg.adjacencyAmd(u);
                for (std::uint32_t a = 0; a < adjacencyU; ++a) w[runU[a]] = other;
                // Index 1: the new clique is at the front of every I[u] and is shared by every
                // member of C[pivot], so it can never discriminate. Ledger entry 6.
                const std::int32_t* incidenceRunU = qg.incidenceAmd(u);
                for (std::uint32_t i = 1; i < incidenceU; ++i) w[incidenceRunU[i]] = other;

                for (std::int32_t v = buckets.chain(u); v != NIL; v = buckets.chain(v)) {
                    if (qg.eliminatedAmd(v)) continue;

                    // The exact test the hash only filters for:
                    //     A[u] == A[v]   and   I[u] == I[v]
                    // against the stamp of u laid down once above. Both walks short-circuit on the
                    // first mismatch, which is what made the comparison cheap and the stamping the
                    // thing that had to move: see experiments/ordering/AmdFlat.md, iteration 15.
                    if (qg.adjacencySize(v) != adjacencyU) continue;
                    if (qg.incidenceSize(v) != incidenceU) continue;

                    bool                same       = true;
                    const std::int32_t* adjacencyV = qg.adjacencyAmd(v);
                    for (std::uint32_t a = 0; a < adjacencyU && same; ++a)
                        if (w[adjacencyV[a]] != other) same = false;
                    if (same) {
                        const std::int32_t* incidenceV = qg.incidenceAmd(v);
                        for (std::uint32_t i = 1; i < incidenceU && same; ++i)
                            if (w[incidenceV[i]] != other) same = false;
                    }
                    if (!same) continue;

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

                    qg.merge(u, v);                 // v folded into u, left where it lies
                    // NO BUCKET TRAFFIC HERE. Both are already out of the lists, taken out when
                    // C[pivot] was formed, and u's mPrev is holding its hash key, so an unfile
                    // would read that key as a link. AMD_2 does no list work here either: it
                    // writes Nv and moves on.
                    degrees[u] -= weightV;
                    ++numEliminated;                // out of the count, not out of the graph
                }
            }
        }

        // THE FOURTH PASS, which finishes the bounds and files them. Amd.cpp spells it
        // `deg = Degree[i] + degme - nvi` then `deg = MIN (deg, nleft - nvi)`, under its RESTORE
        // DEGREE LISTS heading, and it runs here for the reason in the note above: `nvi` is read
        // after supervariable detection, so a vertex that absorbed another subtracts the combined
        // weight. `degme` and `numLeft` were settled before the hash and are the same for
        // everyone, so this is one addition and two comparisons per survivor.
        // AND THE CLIQUE IS TRIMMED AS THIS PASS WALKS IT, which is Amd.cpp's `Iw [p++] = i` under
        // RESTORE DEGREE LISTS: the survivors are written back over the front of the clique and
        // what detection absorbed falls off the end. One store on a walk this pass makes anyway,
        // and without it those vertices are visited by every later walk of this clique.
        // THE SIGNS COME BACK INSIDE THIS PASS, not in mass elimination, which ran late and left
        // them negative. `AMD_2` does the same, `nvi = -Nv [i]` then `Nv [i] = nvi` under RESTORE
        // DEGREE LISTS, and the store rides on the load this pass makes anyway. A vertex the hash
        // absorbed is skipped and never restored, which is right: `merge` zeroed its weight, and
        // zero is the absorbed state whatever sign it arrived with.
        qg.restorePivotWeight(pivot);

        std::int32_t* cliqueOut = qg.clique(pivot);
        std::uint32_t kept      = 0;
        for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) {
            const std::int32_t u = pivotClique[k];
            if (qg.eliminatedAmd(u)) continue;         // absorbed by the hash a moment ago
            const std::uint32_t weightU = qg.restoreWeight(u);   // POST-merge, and un-negated here
            // ONE OPERAND WIDENED, so the sum is formed in `std::size_t`; see the note in Amd1
            // on why widening cannot be done after the addition the way narrowing is done after
            // the subtraction. `partial[u]` and `degme` each reach n, so the sum reaches 2n, and
            // in 32 bits that would fit only because n is capped at 2^31 - 1.
            std::size_t bound = static_cast<std::size_t>(w[u]) + degme - weightU;
            // THE SLOT GOES BACK TO ALIVE-AND-UNSEEN, the last read having just happened. Without
            // it a survivor later chosen as pivot would form a clique whose w already held a
            // bound, which the next step's prune would read as a running value above the tag or as
            // absorbed. A vertex the hash eliminated is skipped and never reset, which is right:
            // it is dead, no clique is ever named after it, and clear_flag resets it anyway.
            w[u] = 1;
            bound = std::min<std::size_t>(bound, numLeft - weightU);
            // THE NARROWING POINT, after the cap, which is where the value is at most n.
            const std::uint32_t filed = static_cast<std::uint32_t>(bound);
            // FILE, NOT REFILE. The vertex has been out of the lists since C[pivot] was formed
            // and its mPrev holds a hash key rather than a link, so the unfile inside refile would
            // read that key as one. It is also one pass less than refile did.
            degrees[u] = filed;
            buckets.file(filed, u);
            // The minimum, taken HERE rather than in a pass of its own. `bound` is in a register
            // and `degrees[u]` has just been written from it, so the pass this replaces was one
            // scattered read per survivor per pivot to recover a value it had already had.
            // AMD_2 does the same inside its restore-degree-lists loop, `if (deg < mindeg)`.
            // Amd1 has always done it this way; Amd2, Amd2B and AmdFlat did not.
            //
            // MEASURED AT ZERO, with the clique-weight fusion beside it: useful cycles unchanged
            // within half a percent in both families. A port and a simplification, not a speed
            // fix. See benchmarks/ordering/README.md (2026-08-09).
            minDegree = std::min(minDegree, filed);
            cliqueOut[kept++] = u;
        }
        qg.trimClique(pivot, kept);

        // The whole array is invalidated in ONE ADDITION rather than by walking the touched list and
        // zeroing each entry. After scan 1 no entry exceeds wflg + lemax, so advancing past that
        // puts every one of them into the stale range. Amd.cpp's `wflg += lemax`.
        //
        // PAST EVERY STAMP THIS STEP LAID DOWN, not merely past the scan's values, since detection
        // now writes into this same array above `wflg + lemax`. Amd.cpp advances wflg through
        // detection for the same reason: a stamp must read as alive-unseen next step, which means
        // strictly below the new tag.
        wflg = stamp + 1;                 // past every stamp this step laid down
    }

    // EVERY CLIQUE IS DEAD BY NOW, every vertex having been eliminated, so the live count must
    // have come back to zero. It is the whole check on the counter: births and deaths balance
    // or they do not, and nothing else in the suite would notice if they did not.
    // THE COUNTER CROSS-CHECKED AGAINST A RECOMPUTATION, which the driver can do exactly because
    // it holds the pivot list and a clique's owner is a pivot. Births and deaths are spread over
    // three call sites and nothing else in the suite would notice if they stopped balancing.
    //
    // IT IS NOT ZERO AT THE END, and that is correct rather than a leak. A clique dies when a
    // member of it becomes a pivot, and at the close of a run the last cliques have had every
    // member MASS ELIMINATED into the pivot instead, so no one is left to absorb them. A handful
    // of entries survive, 1 to 3 on grids from 2 to 5 a side.
    assert(qg.cliqueCountBalances() && "clique births and deaths do not balance");
    gPeakCliqueMembers = qg.numPeakCliqueMembers();   // see include/oblio/QuotientGraph.h
    if (arenaEntries != nullptr) *arenaEntries = qg.arenaEntries();
    // THE ROWS THE DENSE RULE SET ASIDE GO LAST, in index order, which is where `AMD_2`'s output
    // assembly puts them. They were collected in an ascending pass, so appending the vector is
    // that order; each stands only for itself, having been set aside before it could absorb
    // anything, so `order` expands a chain of one.
    pivots.insert(pivots.end(), denseRows.begin(), denseRows.end());
    return qg.order(pivots);
}


} // namespace

std::vector<std::int32_t> orderAmdFlat(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx) {
    return orderAmd3Impl(colPtr, rowIdx, nullptr);
}

// The same ordering, reporting how many entries the clique arena ended up holding. A SIZE, not a
// capacity, and for this scheme also the peak, the arena never shrinking. `benchmarks/matrices`
// prints it beside nnz(L); see QuotientGraph::arenaEntries.
std::vector<std::int32_t> orderAmdFlat(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    std::size_t& arenaEntries) {
    return orderAmd3Impl(colPtr, rowIdx, &arenaEntries);
}

} // namespace Oblio
