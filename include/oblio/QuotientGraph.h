#pragma once

// QuotientGraph.h - the representation Oblio's own minimum-degree orderings run on, and the
// degree buckets they pick from. Shared by Mmd1 and Amd1, which differ only in their drivers:
// one batches eliminations to make the degree refresh rare, the other bounds the degree to make
// each refresh cheap. Everything below that fork is here.
//
// The idea, in one line: an elimination does not create fill edges, it creates a CLIQUE, and a
// clique of d vertices is a d-element list rather than d(d-1)/2 edges. So the neighbor relation
// splits in two, and the true neighborhood is their union, formed on demand and never stored:
//
//     A[u]     the vertices u is still explicitly adjacent to
//     I[u]     the cliques that contain u
//     C[c]     the members of clique c
//
//     reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}
//
// This is George and Liu's reachable set, and it equals the neighborhood the filled graph would
// have. Section 5.3 of archive/sparse_factorization.md carries the derivation; the literature
// says "element" where we say "clique" and writes E_i for I[u].
//
// A[u] and I[u] are the two kinds of SOURCE the union is taken over, and their number falls
// monotonically: an elimination turns a source of the first kind into one of the second and never
// manufactures one. So the two lists share one block, sized once from u's column of A and never
// grown, which is what both vendored routines do and what section 5.3 of
// archive/sparse_factorization.md proves they may. C[c] has no such bound and gets an arena.
//
// A clique id is the pivot that created it, so cliques index into the vertex space and no
// separate id space is needed. A dead clique is an empty member list: absorption clears the
// members, and an empty clique contributes nothing to any reachable set or count, so no
// liveness flag is carried.
//
// No sets anywhere. Membership is a mark array stamped with a monotone tag, so a query is one
// comparison and nothing is allocated, which is what SymFactorEngine does with its index sets
// and what the vendored routines do with marker/tag. Every list edit is a compaction in place.
//
// The prototypes these were pulled from are experiments/ordering/mmd1 and amd1, where the same
// structures carry the tracing and counting that a production ordering has no use for.

#include "oblio/Types.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The degree buckets: one doubly linked list per degree, threaded through arrays of size n, so
// that filing, unfiling and taking the head are all O(1). MMD spells these fwd/bwd and AMD
// Next/Last. An ordered container cannot give O(1) removal from the middle, which is what a
// degree change needs and what happens far more often than a pick.
//
// The bodies are inline because they are single-statement pointer splices on the hot path, the
// same exception the tree makes for trivial accessors.
class Buckets {
public:
    // The heads are indexed by DEGREE and the links by VERTEX, which are not the same range. A
    // degree is at most size - 1 while a branch is filing degrees, but MMD's convention files a
    // vertex under its degree plus one, so the heads need one more slot than there are vertices.
    // One int32 either way, and getting it wrong is an out-of-bounds write on a one-vertex matrix.
    // NOT FILED, distinct from NIL, and it is what lets this class hold three arrays instead of
    // four. `unfile` used to leave mPrev[u] = NIL, which is also what a FILED vertex at the head
    // of its bucket has, so the state was ambiguous and a separate mFiled byte was needed to
    // disambiguate it. A sentinel of its own removes the ambiguity and the array with it.
    //
    // What that array was costing, measured on a 140x140 grid: amd3 wrote it 224054 times, once
    // per file and once per unfile, and read it ZERO times, `filed()` is MMD's accessor and no amd
    // driver calls it. mmd3 read it 115931 times. So one branch paid for a byte array it never
    // consulted, and the other kept information it could have derived from a link it already
    // touches. Amd.cpp has Head, Next and Last and no flag, for the same reason.
    // A DEGREE IS ONE DIMENSIONAL and so `std::uint32_t` in every signature here. The links stay
    // `std::int32_t`, carrying NIL and UNFILED; only the key does not. The heads are indexed by
    // degree and the links by vertex, which are not the same range: a true degree is at most
    // n - 1, but MMD files a vertex under its degree PLUS ONE, so the key reaches n and the heads
    // need one slot more than there are vertices.
    static constexpr std::int32_t UNFILED = -2;

    explicit Buckets(std::size_t size)
        : mHead(size + 1, NIL), mNext(size, NIL), mPrev(size, UNFILED) {}

    // buckets[degree].add(u), at the head. The head is the only O(1) end of a singly reachable
    // list, so the winner among equal degrees is whatever was filed last rather than the lowest
    // index. That is the vendored convention and it is why an ordering differs from an exact
    // scan's in its ties.
    void file(std::uint32_t degree, std::int32_t u) {
        mNext[u] = mHead[degree];
        mPrev[u] = NIL;                                // filed, and at the head
        if (mHead[degree] != NIL) mPrev[mHead[degree]] = u;
        mHead[degree] = u;
    }

    // buckets[degree].discard(u). Idempotent, which matters during a batch: a vertex evicted
    // early can be merged away by a later pivot in the same round, and unfiling it twice must
    // not splice a list it is no longer in.
    void unfile(std::uint32_t degree, std::int32_t u) {
        if (mPrev[u] == UNFILED) return;               // already gone; see the sentinel's note
        if (mPrev[u] != NIL) mNext[mPrev[u]] = mNext[u];
        else                 mHead[degree]   = mNext[u];
        if (mNext[u] != NIL) mPrev[mNext[u]] = mPrev[u];
        mNext[u] = NIL;
        mPrev[u] = UNFILED;
    }

    // Move u to the bucket for newDegree, carrying the cached degree with it. The three steps go
    // together, which is why this is one call: the bucket a vertex sits in is read from its
    // degree, so writing the degree first would erase it from the wrong list. A vertex whose
    // degree did not change is removed and reinserted into the same list, which is harmless.
    void refile(std::vector<std::uint32_t>& degrees, std::int32_t u, std::uint32_t newDegree) {
        unfile(degrees[u], u);
        degrees[u] = newDegree;
        file(newDegree, u);
    }

    // The next vertex in the same bucket, and whether u is filed at all. MMD reads both: it
    // walks a whole bucket in the prepass, and its refresh asks whether a vertex it reached has
    // already been dealt with this round.
    std::int32_t next(std::int32_t u) const      { return mNext[u]; }
    bool         filed(std::int32_t u) const     { return mPrev[u] != UNFILED; }

    std::int32_t head(std::uint32_t degree) const  { return mHead[degree]; }
    bool         empty(std::uint32_t degree) const { return mHead[degree] == NIL; }

private:
    std::vector<std::int32_t> mHead;   // mHead[d], the first live vertex of degree d
    std::vector<std::int32_t> mNext;   // mNext[u], toward the tail
    std::vector<std::int32_t> mPrev;   // mPrev[u], toward the head
    // A byte per vertex, not std::vector<bool>. That specialization packs one bit per element, and
    // what it costs here is CONSTRUCTION rather than access: a vector of n false is built through
    // the bit-reference machinery word by word, where a vector of n zero bytes is a memset. Both
    // this flag and QuotientGraph::mEliminated are built once per ordering, so a workload that
    // orders repeatedly pays it repeatedly.
    //
    // Measured on alpamayo, MMD2 at 140x140 over 3000 orderings: the graph constructor fell from
    // 710 ms to 186 ms, two traces agreeing, which is essentially the whole of a 12 percent
    // saving. Every branch gained about the same amount whatever its access pattern, which is the
    // signature of a per-construction cost rather than a per-read one.
    //
    // The access cost is real in principle, a load, a variable shift and a mask against one byte
    // load, and the write side is a read-modify-write of a byte eight vertices share. It did not
    // measure: the rows for these two flags moved by less than the run-to-run noise. So the
    // footprint, 39 KB at n = 19600 against a clique arena holding 115263 entries, buys the
    // construction and not the walk.
    //
    // The prototypes in experiments/ordering keep the bit-packed form, which is right: they
    // construct one graph per run.
};

// What an approximate-degree driver accumulates while the eliminator is already walking the lists,
// handed to the second `eliminate` overload so that the walk serves both. The members are the
// driver's own arrays, held by reference: the graph fills them and owns none of them.
//
// `tag` is the only member that moves, and the driver sets it before each elimination, exactly as
// it would before its own scan. The rest are bound once and reused for the whole ordering.
struct ApproximateScan {
    std::vector<std::uint32_t>&       explicitPart;  // per vertex, sum of weight over the pruned A[u]
    std::vector<std::uint32_t>&       outside;       // per clique, |C[c] - C[p]| weighted
    const std::vector<std::uint32_t>& cliqueDegree;  // per clique, |C[c]| weighted
    std::vector<std::int32_t>&        touchedCliques;// the cliques this step reached, once each
    std::vector<std::int32_t>&        mark;          // the driver's membership scratch
    std::int32_t                      tag;           // its stamp for this elimination
};

// The same thing for a driver that carries `Amd.cpp`'s TAGGED W ARRAY instead of a value array and
// a separate seen-this-step mark. One array holds three facts: `w[c] == 0` is absorbed, `0 < w[c] <
// wflg` is alive but stale, and `w[c] >= wflg` is seen this step with `w[c] - wflg` the value. So
// there is no mark to carry and no clearing pass, and first sighting is `w[c] != 0 && w[c] < wflg`
// rather than a tag comparison. `Amd3` uses that encoding; `Amd1` and `Amd2` use the other, which
// is why this is a second struct rather than a flag on the first.
//
// `key` is the odd one and is here because of where the walks are. `Amd3` builds its hash key from
// the PRUNED A[u] and the FINAL I[u], and a driver that folds its first scan in here keeps a walk
// of I[u] but loses its second walk of A[u], so the ADJACENCY HALF has nowhere else to go. This
// fills that half and the driver adds the other.
//
// **ONLY THAT HALF, and the reason is a phase boundary rather than a preference.** Aggressive
// absorption runs between this prune and the driver's bound pass and COMPACTS I[u] in place,
// dropping the cliques it killed. So the list the key must sum over does not exist yet here.
// Accumulating it anyway is correct on square grids, where absorption rarely fires, and wrong on
// cubic and random ones, which is how it was found.
//
// **REDUCED AS IT ACCUMULATES**, modulo the driver's bucket count, which is why it is an int32 and
// not a size_t. The key is a SUM and the driver only ever uses it modulo that number, so reducing
// early cannot change which bucket a vertex lands in. What it buys is the width: the running value
// stays below the modulus, so this fits in an array the driver already has rather than needing one
// of its own. Two extra arrays of size n across the phase boundary cost 12 percent in 2D at 400 a
// side when this fusion was first built, which is the footprint trade REPORT.md names and the
// same one that sank the 2026-08-08 key fusion.
struct TaggedScan {
    std::vector<std::uint32_t>&       explicitPart;  // per vertex, sum of weight over the pruned A[u]
    std::vector<std::int32_t>&        key;           // per vertex, the whole hash key, reduced
    std::vector<std::int32_t>&        w;             // per clique, Amd.cpp's tagged W
    const std::vector<std::uint32_t>& cliqueDegree;  // per clique, |C[c]| weighted
    std::vector<std::int32_t>&        touchedCliques;// the cliques this step reached, once each
    std::int32_t                      wflg;          // the tag for this elimination
    std::int32_t                      modulus;       // the driver's bucket count, n + 1
};

// The quotient graph itself: the three lists above, the liveness flags, and the supervariable
// members that mass elimination grows. A driver owns one of these, picks a pivot, calls
// eliminate, and refreshes whatever the elimination reached.
class QuotientGraph {
public:
    // Built straight from A's pattern. A is stored with both triangles, so a column's rows are
    // already its neighbors, and the only conversion is dropping the diagonal, which is a self
    // loop and says nothing about fill. Oblio's other input assumptions (sorted, unique, both
    // triangles, a structurally present diagonal) hold by construction, which is why nothing here
    // symmetrizes, deduplicates or sorts. See the pass-5 discussion in
    // experiments/ordering/README.md.
    QuotientGraph(const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const           { return mAdjacencySize.size(); }
    bool eliminated(std::int32_t u) const { return mEliminated[u] != 0; }

    // The members of clique c, which after eliminate(p) is the pattern of p's column of L: the
    // vertices the pivot reached, less those it absorbed. A pointer and a length, as the
    // adjacency is, for the same reason.
    //
    // **Valid until the next eliminate.** The members live in one arena that grows as cliques are
    // formed, so a reallocation moves them; the offsets are indices and survive it, but a pointer
    // taken beforehand does not. Read a clique at the moment of use, which is the same rule the
    // numeric factor's blocks live by.
    const std::int32_t* clique(std::int32_t c) const {
        return mCliqueArena.data() + mCliquePtr[c];
    }
    std::uint32_t cliqueSize(std::int32_t c) const { return mCliqueSize[c]; }

    // |C[pivot]| WEIGHTED, which every amd driver needs and all four used to compute for
    // themselves in a pass of their own, one scattered weight load per member per pivot. It is
    // accumulated in `beginElimination`'s stamping walk instead, which has the member in hand
    // already, so the pass is gone and nothing is added. AMD_2 does the same, `degme += nvi`
    // inside the loops that build the element, at its lines 1492 and 1636.
    //
    // MEASURED AT ZERO, and that is the point of saying so here. CPU Counters before and after,
    // same session: useful cycles unchanged within half a percent in both families, AMD3 still
    // 1.61x the vendored routine on cubes and 1.56x in 2D. The 25 visits per pivot this removes on
    // cubes were contiguous walks over hot arrays and cost near nothing. It is kept as a faithful
    // port that removes two passes and a stale-value hazard, NOT as a speed fix, and it should not
    // be cited as one. `benchmarks/ordering/README.md` (2026-08-09) has the numbers.
    //
    // VALID UNTIL THE NEXT ELIMINATION, and a scalar rather than an array for that reason: it is
    // read immediately after `eliminate` and never afterwards. Same contract as the pointer
    // `clique()` returns, and stated here because a scalar carrying per-pivot state is the kind
    // of thing that goes stale silently.
    //
    // Over the UNTRIMMED clique, which is what the drivers want: mass elimination runs after the
    // absorption in Amd3, and the value that pass needs is the one before any trimming. Amd3
    // recomputes it afterwards over the trimmed clique, which is ledger entry 7 and stays.
    std::uint32_t cliqueWeight() const { return mCliqueWeight; }

    // The two halves of the neighbor relation, read by an approximate degree, which decomposes
    // reach(u) rather than forming it: A[u] contributes its own term and each clique in I[u]
    // contributes one number computed for the clique rather than for u. An exact degree has no
    // use for either, since it unites them through reachableSet instead.
    //
    // Both are a pointer and a length rather than a container, because that is what their storage
    // is: the two lists share one run of one flat array. A[u] starts at the run and I[u] starts
    // immediately behind it, which is why the incidence lookup reads the adjacency's length. See
    // the note on the members below.
    const std::int32_t* adjacency(std::int32_t u) const {
        return mSource.data() + mSourcePtr[u];
    }
    std::uint32_t adjacencySize(std::int32_t u) const { return mAdjacencySize[u]; }

    const std::int32_t* incidence(std::int32_t u) const {
        return mSource.data() + mSourcePtr[u] + mAdjacencySize[u];
    }
    std::uint32_t incidenceSize(std::int32_t u) const { return mIncidenceSize[u]; }

    // How many original vertices u stands for. One until mass elimination merges into it, so a
    // degree that counts vertices has to count this rather than entries.
    //
    // Held in a flat array, and this is the condition under which the prototypes said it would
    // have to be: its members are a chain rather than a list, so the size is no longer free to
    // read. Caching it purely for locality was measured earlier and made no difference at all,
    // 2.77x against 2.76x, so the array is here for the chain and not for the cache.
    std::uint32_t weight(std::int32_t u) const { return mWeight[u]; }


    // reach(u), as above. Not const: the mark array and its tag are scratch, and threading them
    // through every call site is what the prototypes do only because their display functions
    // needed to borrow them.
    void reachableSet(std::int32_t u, std::vector<std::int32_t>& reached);
    std::vector<std::int32_t> reachableSet(std::int32_t u);

    // |reach(u)|, counted without materializing it. A maintained degree wants the size and never
    // the set, so the same two passes run with a counter where reachableSet has a push_back, and
    // an exact refresh stops allocating once per refreshed vertex. That is the dominant
    // allocation in a branch that keeps its degrees exact, and none at all in one that bounds
    // them, which is most of why the two differ in speed by more than they differ in work.
    std::uint32_t reachableSize(std::int32_t u);

    // The same set, weighted: the number of ORIGINAL vertices reach(u) stands for. Once a branch
    // merges into live vertices, a reached vertex can stand for several, and a degree that counts
    // vertices has to count those rather than entries.
    std::uint32_t reachableWeight(std::int32_t u);

    // Eliminate the pivot: turn it into a clique, absorb the cliques it belonged to, prune the
    // edges the new clique now implies, and merge in whatever it makes indistinguishable.
    // Returns the vertices merged into the pivot's supervariable, which the caller needs in
    // order to take them out of its own bookkeeping.
    //
    // **The returned reference is valid until the next eliminate**, being a scratch buffer whose
    // capacity survives from pivot to pivot. Returning by value cost one allocation per
    // elimination that merged anything, about 1700 per ordering at 140x140. A caller that needs
    // the list to outlive the next call copies it.
    //
    // This is now the ONLY scratch here. The reachable set used to have one too and no longer
    // does: the walk writes it straight into the clique arena, since C[pivot] is the reach and
    // there was never a reason for it to exist anywhere else first. See beginElimination.
    //
    // In set operations, and the order is the order below:
    //
    //     C[p] = reach(p)                    absorb into C[p]
    //     C    = C - I[p]                    reclaim I[p]
    //     for u in C[p]:
    //         A[u] = A[u] - C[p] - {p}       prune
    //         I[u] = ( I[u] - I[p] ) | {p}   absorb into C[p], reclaim I[p]
    //     merged = { u in C[p] : A[u] == {} and I[u] == {p} }
    //     C[p]   = C[p] - merged
    //
    // The new clique is C[p] and gets no name of its own, which is what makes the first line
    // read as what an elimination is: the pivot stops being a vertex with a reachable set and
    // becomes a clique holding that same set. Not one of the three differences builds a set;
    // each is a stamp of the subtrahend and one compaction pass over the minuend.
    const std::vector<std::int32_t>& eliminate(std::int32_t pivot);

    // The same elimination, with an approximate-degree driver's first scan folded into it.
    //
    // Why the two are one call rather than two loops. An approximate degree decomposes reach(u)
    // instead of forming it, so its driver walks exactly the lists the prune has just walked:
    // A[u] to sum the weights of what survived, I[u] to subtract this vertex's weight from each
    // clique's outside count. Run afterwards, that is a second and third visit to every element,
    // which is where AMD1 spends what the vendored routine does not: 483677 element visits against
    // 216662 on a 100x100 grid. Folded in, A[u] is visited once and I[u] twice, which is what
    // `amd_2`'s two scans cost.
    //
    // It cannot be folded further, and the reason is a property of the bound rather than of the
    // code: `outside[c]` is complete only once every member of C[p] has been seen, so the sum over
    // I[u] that the bound needs is a second pass by construction. `amd_2` has two scans for the
    // same reason.
    //
    // Two things make the fusion sound, and neither is obvious. The scan now runs over the
    // UNTRIMMED C[p], where the driver ran it over the trimmed one, and the difference is the
    // mass-eliminated vertices; they contribute nothing, since the merge test requires
    // I[u] == {pivot} and the scan skips the pivot. And the weights summed over A[u] are read
    // before mass elimination could change any of them, which is safe because the prune has
    // already removed every member of C[p] from A[u] and mass elimination merges only members of
    // C[p]. Neither would survive reordering the phases.
    const std::vector<std::int32_t>& eliminate(std::int32_t pivot, ApproximateScan& scan);

    // The same fusion for the tagged-W encoding, which is the one `Amd3` carries. Everything the
    // overload above says about why the fold is sound applies here unchanged; only the three-facts
    // array differs. It also fills the hash key's adjacency half, which that overload has no need
    // of because its callers keep a separate walk of A[u].
    const std::vector<std::int32_t>& eliminate(std::int32_t pivot, TaggedScan& scan);


    // Fold v into u, the two having been found indistinguishable from EACH OTHER rather than
    // from a pivot. Unlike mass elimination, which merges into a vertex that is being eliminated
    // in the same breath, this merges into one that stays live, so u carries v's weight onward.
    //
    // v is left exactly where it lies, at weight zero, rather than being purged from every clique
    // and adjacency that names it, which would cost a pass over the structure per merge. Nothing
    // is lost by leaving it: the caller's test required that every list holding v holds u as well,
    // so v is redundant wherever it appears and never the only way to reach anything, and a
    // weight of zero makes it invisible to every count.
    void merge(std::int32_t u, std::int32_t v);

    // Number u without eliminating it in the quotient-graph sense: no clique is formed, nothing is
    // pruned, and its neighbors keep degrees that still count it. That staleness is the point, and
    // it is what genmmd's prepass does with the degree-0 and degree-1 vertices. The lists are left
    // alone; every walk skips a numbered vertex from here on.
    void number(std::int32_t u);

    // Walk I[u] from the back in reachableSet(), matching genmmd's element stack. A tie-break
    // convention and nothing else: it changes which permutation comes out, never which sets are
    // computed. See the member's note.
    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    // Lay the lists out the way AMD_2 does, which is the opposite of genmmd's on both counts.
    // Two conventions under one switch because they are one fact, that AMD's lists run the other
    // way round, and are only ever wanted together. Like the flag above, this changes which
    // permutation comes out and never which sets are computed. Used by Amd3 alone.
    //
    //   reachableSet   walks the CLIQUES before the explicit adjacency, since AMD_2's
    //                  `for (knt1 = 1; knt1 <= elenme + 1; knt1++)` takes the elements of me and
    //                  the supervariables only on its last pass. genmmd is the other way round,
    //                  which is why the md ladder is laid out as it is.
    //   the prune      puts the new clique at the FRONT of I[u] rather than appending, with the
    //                  displaced entries ROTATED rather than shifted, which is AMD_2's
    //                  `Iw[pn] = Iw[p3]; Iw[p3] = Iw[p1]; Iw[p1] = me`. Our two lists share one
    //                  run exactly as its do, so the same three moves apply unchanged.
    //
    // See experiments/ordering/AMD3.md, ledger entries 2 and 5.
    void setVendoredListOrder(bool on) { mVendoredListOrder = on; }

    // Stop eliminate() at the prune, leaving mass elimination to the caller. AMD_2 makes the same
    // test in its scan 2, AFTER aggressive absorption has dropped every clique lying inside the
    // new one, and says why in its own comment: with aggressive absorption, `deg == 0` is
    // identical to the structural test. Asking first, which is what eliminate() does by default,
    // asks it of an I[u] that still holds cliques about to be removed, so the cheap test declines
    // vertices AMD merges.
    //
    // With this on, eliminate() returns an EMPTY merged list and C[pivot] is reach(pivot) exactly,
    // and the caller must call massEliminate() once it has absorbed. Used by Amd3 alone. See
    // experiments/ordering/AMD3.md, ledger entry 3.
    void setLateMassElimination(bool on) { mLateMassElimination = on; }

    // The half eliminate() no longer does under the flag above: fold into the pivot's supervariable
    // every member of C[pivot] that the new clique now accounts for entirely, and trim C[pivot] of
    // them. Returns the merged vertices, from a member scratch as eliminate() does.
    //
    //     merged   = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
    //     C[pivot] = C[pivot] - merged
    //
    // Calling it without the flag is a caller error and not guarded: eliminate() will already have
    // merged, so this would find nothing and cost a pass.
    const std::vector<std::int32_t>& massEliminate(std::int32_t pivot);

    // Kill these cliques and take them out of the incidence lists of these vertices. A clique
    // dies when it is found to lie wholly inside a newer one, which is aggressive absorption, and
    // the vertices to purge are the newer clique's members, since those are the only lists that
    // can still name it.
    // `vertexCount` is `|C[p]|`, so one dimensional; every caller passes `cliqueSize(pivot)`.
    void absorb(const std::vector<std::int32_t>& cliques,
                const std::int32_t* vertices, std::uint32_t vertexCount);

    // Expand a pivot sequence into an elimination order over the original vertices. A pivot
    // stands for its whole supervariable, whose members are eliminated consecutively, so this is
    // where a supervariable of size w becomes w columns.
    std::vector<std::int32_t> order(const std::vector<std::int32_t>& pivots) const;

    // The same permutation with each supervariable's members in ASCENDING VERTEX INDEX rather
    // than merge order. Indistinguishable members, so the fill and the forest are unchanged;
    // only the permutation is, and it is genmmd's. Used by Mmd3 alone. See the .cpp.
    std::vector<std::int32_t> orderAscending(const std::vector<std::int32_t>& pivots) const;

private:
    // The head and tail of an elimination, shared by the two overloads above and private because
    // between them the graph is half eliminated: the clique is written and stamped but the reached
    // vertices still name the pivot as a variable. Nothing outside may observe that state, which is
    // why the seam is two private calls rather than a public begin and end.
    void beginElimination(std::int32_t pivot, std::int32_t& inClique, std::int32_t& absorbed);

    const std::vector<std::int32_t>& finishElimination(std::int32_t pivot);

    // A[u] and I[u] share one run, and one array holds every run end to end. The two lists are
    // the SOURCES of reach(u), one per explicit neighbor and one per clique, and their number is
    // conserved: each elimination that reaches u replaces at least one source with the new
    // clique, the explicit neighbor `pivot` where u was reached through A[pivot], or an absorbed
    // clique where it was reached through one. So
    //
    //     |A[u]| + |I[u]| <= the number of off-diagonal entries in u's column of A
    //
    // holds for the whole run, and u's block can be sized once from the pattern and never grown.
    // Section 5.3 of archive/sparse_factorization.md carries the argument; 5.15 records that both
    // vendored routines rely on it and that neither states it.
    //
    // The order within the run is forced rather than chosen. The prune compacts A[u] and then
    // I[u], and A[u] shrinks by at least what I[u] gains, so writing the incidence behind the
    // adjacency always trails the read cursor. The other order would have the incidence overwrite
    // an adjacency entry it had not read yet.
    //
    // C[c] gets none of this and needs its own arena below: a clique's members are its pivot's
    // reach, which nothing about the pivot's own column bounds.
    // ONE DIMENSIONAL SIZES ARE `std::uint32_t`, POSITIONS ARE `std::size_t`. `mSourcePtr` offsets
    // into the pool and is bounded by nnz(A), so it is two dimensional and stays wide; the two
    // lengths are bounded by deg(u) and so by n, which the constructor caps at `MAX_IDX`. The two
    // kinds meet at exactly one place, the constructor's crossing, where a difference of positions
    // is written as a count, and that is where the cast belongs. `mIncidenceSize` has no such
    // crossing at all: it is only ever written from a cursor or from another length.
    //
    // The CONSERVATION LEMMA is why one uint32 covers both. Their sum is bounded by u's column of
    // A for the whole run, so neither can reach n on its own where the other is nonzero, and
    // `adjacencySize(u) + incidenceSize(u)` cannot overflow either.
    std::vector<std::int32_t>  mSource;         // every A[u] then I[u], run after run
    std::vector<std::size_t>   mSourcePtr;      // where u's run starts, fixed at construction
    std::vector<std::uint32_t> mAdjacencySize;  // A[u]'s length, from the run's start
    std::vector<std::uint32_t> mIncidenceSize;  // I[u]'s length, immediately behind A[u]

    // C[c] is flat too, and for a second property rather than the adjacency's. Its members are
    // not known in advance, so its block cannot be placed at construction; but they are known
    // exactly when the clique is formed, at the moment its pivot is eliminated, and the set only
    // ever shrinks afterwards. So a bump allocator suffices: take the next block, write the
    // members, and that block is the clique's for as long as it lives.
    //
    // Nothing is reclaimed. An absorbed clique leaves a hole, which is what the vendored pool
    // answers with a compaction pass, and we do not need the answer: the arena's total is the sum
    // of every clique ever formed, which is bounded by nnz(L) and is a few megabytes at the sizes
    // we run. The offsets being indices rather than pointers is what lets the arena simply grow
    // instead, since a reallocation leaves every offset valid.
    // The same split as the source pool above, and for the same reason: `mCliquePtr` offsets into
    // the arena, which grows toward nnz(L), so it is two dimensional; `mCliqueSize` is a member
    // count bounded by n. The one crossing is in `beginElimination`, where the arena's new length
    // less the block's start is written as a count.
    std::vector<std::int32_t>  mCliqueArena;  // every C[c] ever formed, end to end
    std::vector<std::size_t>   mCliquePtr;    // where c's block starts, fixed once written
    std::vector<std::uint32_t> mCliqueSize;   // how much of it is still live

    // The supervariable a vertex stands for, as a chain rather than a list per vertex. A list
    // meant one allocation per vertex before anything had happened, n of them for a structure
    // that is usually a singleton, plus a growth every time one absorbed another. The chain is
    // three arrays allocated once: the next member, the last one (so an absorption appends in
    // O(1) and the members keep their order, which the emitted permutation depends on), and the
    // count, which a chain no longer gives away for free.
    // `mWeight` is one dimensional and so `std::uint32_t`, and the bound is not the term count but
    // DISJOINTNESS: the weights partition the original vertices, so a sum of them over a set of
    // distinct vertices is at most n however many terms it has. That covers every accumulation of
    // weights in the ordering, in `reachableWeight`, in the clique weight below, and in the
    // drivers' bound and refresh passes. An accumulation over an unbounded number of one
    // dimensional terms would otherwise need a wider type, and this is the exception to that.
    std::vector<std::int32_t>  mSuperNext;
    std::vector<std::int32_t>  mSuperLast;
    std::vector<std::uint32_t> mWeight;
    std::uint32_t              mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex
    std::vector<std::uint8_t> mEliminated;   // a byte per vertex; see reachableSet on why the flag stays

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away
    // Whether any live merge has happened yet, which decides whether the reachable set has to
    // check for dead vertices at all. Mass elimination cannot leave one in a clique, so a branch
    // that only mass-eliminates never pays the test: the flag is false, the condition
    // short-circuits, and the load of the eliminated flag never happens. Measured at about a
    // quarter of the exact refresh, which is MMD1's hot loop, so it is worth the member.
    bool mLiveMerges = false;

    // Which end of I[u] reachableSet() walks from. genmmd threads its element list through an
    // integer array and pushes at the head, `list[nb] = el; el = nb`, then reads from the head, so
    // the element seen LAST is expanded FIRST; we hold a vector and append. Same set either way and
    // the same cost, but the order decides C[pivot]'s order, hence which of two equal-degree
    // candidates a later iteration finds first, and minimum degree is settled by exactly that.
    // Off by default, so every existing driver is unaffected; Mmd3 turns it on. See
    // experiments/ordering/mmd3.py, where the same four walks are reversed together.
    bool mReverseIncidence = false;

    // AMD_2's list conventions, off by default so the other five drivers are untouched. Read in
    // reachableSet() and in the prune, hoisted at both sites for the same reason mReverseIncidence
    // is: a member load the compiler cannot prove is unaliased by the stores in the loop. See the
    // setter for what each half does and why they are one flag.
    bool mVendoredListOrder = false;

    // Whether eliminate() stops at the prune and leaves mass elimination to the caller. Off by
    // default. See the setter, and massEliminate() for the half it hands over.
    bool mLateMassElimination = false;

    std::vector<std::int32_t> mMark;     // membership scratch, read against mTag
    std::int32_t              mTag = 0;
};

} // namespace Oblio
