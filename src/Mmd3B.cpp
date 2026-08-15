#include "oblio/Mmd3B.h"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

// Mmd3B.cpp - Mmd3 on a different clique STORAGE scheme, everything else identical.
//
// TEMPORARY, AND IT HAS A STOP CONDITION. This file carries a private copy of QuotientGraph,
// named QuotientGraphB, so the storage scheme can be changed and measured without touching the
// class that six drivers share. It exists until the scheme is measured against Mmd3; then either
// QuotientGraph adopts it and this file goes, or the scheme is rejected and this file goes. It is
// not a second implementation to be maintained.
//
// THE COPY IS EXACT AS OF THIS COMMIT and the comments were not duplicated with it. Every design
// note for these types lives in include/oblio/QuotientGraph.h and src/QuotientGraph.cpp and is
// authoritative there; reading this file means reading those alongside it. Only the DIFFERENCES
// are commented here, and at this commit there are none: Mmd3B must return Mmd3's permutation,
// which is genmmd's, so `make mmdorder` and `make test` are the acceptance check unchanged.
//
// What the scheme will become, and why: see the section "The vendored storage scheme, and what it
// is worth" in experiments/ordering/README.md. In one line, genmmd keeps every clique in the dead
// segment of the pivot that formed it, which places blocks in vertex-id order and costs no
// storage at all, and that placement is worth 0.41 of the time ratio on 2D grids.

namespace Oblio {
namespace {

class BucketsB {
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

    explicit BucketsB(std::size_t size)
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
    // this flag and QuotientGraphB::mEliminated are built once per ordering, so a workload that
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
struct ApproximateScanB {
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
struct TaggedScanB {
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
class QuotientGraphB {
public:
    // Built straight from A's pattern. A is stored with both triangles, so a column's rows are
    // already its neighbors, and the only conversion is dropping the diagonal, which is a self
    // loop and says nothing about fill. Oblio's other input assumptions (sorted, unique, both
    // triangles, a structurally present diagonal) hold by construction, which is why nothing here
    // symmetrizes, deduplicates or sorts. See the pass-5 discussion in
    // experiments/ordering/README.md.
    QuotientGraphB(const std::vector<std::size_t>&  colPtr,
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
        return mSource.data() + mSourcePtr[c];
    }
    std::uint32_t cliqueSize(std::int32_t c) const { return mCliqueSize[c]; }

    // Every read of a clique goes through here, because a clique is no longer one flat run. The
    // walk counts members down from mCliqueSize[c] and follows any negative entry it meets; a
    // link is not a member and does not count. The last entry of a segment is where a link is
    // written, so a walk that has members left always meets one before running off the end, and
    // one that does not stops on the count.
    //
    // Templated on the body so it inlines: these are the hottest loops in the ordering, and a
    // std::function or a virtual would cost more than the whole scheme saves.
    template <class F>
    void forEachMember(std::int32_t c, F f) const {
        std::size_t   p = mSourcePtr[c];
        std::uint32_t n = mCliqueSize[c];
        while (n != 0) {
            const std::int32_t v = mSource[p];
            if (v < 0) { p = mSourcePtr[-v - 1]; continue; }   // a link, see mCliqueSize's note
            f(v);
            --n;
            ++p;
        }
    }

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
    const std::vector<std::int32_t>& eliminate(std::int32_t pivot, ApproximateScanB& scan);

    // The same fusion for the tagged-W encoding, which is the one `Amd3` carries. Everything the
    // overload above says about why the fold is sound applies here unchanged; only the three-facts
    // array differs. It also fills the hash key's adjacency half, which that overload has no need
    // of because its callers keep a separate walk of A[u].
    const std::vector<std::int32_t>& eliminate(std::int32_t pivot, TaggedScanB& scan);


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
    // THE CHANGE UNDER TEST. There is no clique arena. C[c] lives in the segment of the vertex
    // that formed it, mSource[mSourcePtr[c] ...], which is dead the moment c is eliminated, so a
    // clique costs no storage and is addressed by its pivot. When the reach outgrows that segment
    // the list CONTINUES in the segment of a clique the pivot has just absorbed, and the last
    // entry of the segment being filled holds a LINK to it.
    //
    // A link is -(c + 1) and not genmmd's -c: ids here start at 0, so -c cannot distinguish
    // clique 0 from member 0. The encoding maps [0, 2^31-1] onto [-2^31, -1] exactly, no value
    // spare and none overlapping, which is what makes a sign test the whole discriminator.
    //
    // genmmd ends a list with a 0 entry, which works only because its ids start at 1. We keep the
    // length instead, in mCliqueSize, which every walk here counts down; a length is also what
    // several passes want without walking, and AMD_2 keeps one too, as `Len`.
    std::vector<std::uint32_t> mCliqueSize;   // members of C[c], links not counted
    std::vector<std::int32_t>  mAbsorbed;     // I[pivot], copied out before its segment is overwritten

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


QuotientGraphB::QuotientGraphB(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mSourcePtr(colPtr.empty() ? 1 : colPtr.size()),
      mAdjacencySize(colPtr.empty() ? 0 : colPtr.size() - 1),
      mIncidenceSize(mAdjacencySize.size(), 0),
      mCliqueSize(mAdjacencySize.size(), 0),
      mEliminated(mAdjacencySize.size(), 0),
      mMark(mAdjacencySize.size(), NIL) {
    const std::int32_t size = static_cast<std::int32_t>(mAdjacencySize.size());

    // One pass to place each column's run, dropping the diagonal. The runs are laid out in column
    // order and never move afterwards, so the offsets are written once here and only the lengths
    // change. The run is u's whole allowance: A[u] fills it now and I[u] grows into the room A[u]
    // gives up, never past mSourcePtr[u + 1].
    mSource.reserve(colPtr.empty() ? 0 : colPtr.back());
    mSourcePtr[0] = 0;
    for (std::int32_t aj = 0; aj < size; ++aj) {
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSource.push_back(rowIdx[cp]);
        mSourcePtr[aj + 1]  = mSource.size();
        // THE ONE CROSSING. A difference of two positions is a count, so this is the single place
        // in the class where a two-dimensional quantity is written into a one-dimensional one. It
        // is bounded by deg(aj) and so by n, which the SparseMatrix constructor has already capped
        // at MAX_IDX, but the cast is written rather than left implicit because that bound is an
        // argument and not something the types say.
        mAdjacencySize[aj]  = static_cast<std::uint32_t>(mSourcePtr[aj + 1] - mSourcePtr[aj]);
    }

    // There is no second arena to reserve: cliques live in mSource, in the segments their
    // pivots vacate. That is the change this file exists to measure.

    // Every vertex begins as a supervariable of one: a chain holding itself, so next is NIL and
    // last is the vertex. Mass elimination splices these together and order() walks them.
    mSuperNext.assign(size, NIL);
    mSuperLast.resize(size);
    mWeight.assign(size, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

std::uint32_t QuotientGraphB::reachableWeight(std::int32_t u) {
    const bool live = mLiveMerges;
    ++mTag;
    std::uint32_t reached = 0;   // a sum over DISTINCT vertices, so bounded by n; see the header
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per element. Measured at 300 ms of AMD1's
    // 6.31 s on alpamayo, in the one accessor that reads as though it were free.
    const std::int32_t* source        = mSource.data() + mSourcePtr[u];
    const std::uint32_t adjacencySize = mAdjacencySize[u];
    const std::uint32_t incidenceSize = mIncidenceSize[u];
    for (std::uint32_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v = source[k];
        if (!live || mEliminated[v] == 0) { mMark[v] = mTag; reached += mWeight[v]; }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::uint32_t i = 0; i < incidenceSize; ++i) {
        const std::int32_t c = incidence[i];
        forEachMember(c, [&](std::int32_t v) {
            if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) { mMark[v] = mTag; reached += mWeight[v]; }
        });
    }
    return reached;
}

void QuotientGraphB::number(std::int32_t u) {
    mEliminated[u] = 1;
    mLiveMerges    = true;   // a numbered vertex lingers in lists, so the walks must filter
}

void QuotientGraphB::beginElimination(std::int32_t pivot,
                                     std::int32_t& inClique, std::int32_t& absorbed) {
    // The reach is written STRAIGHT INTO THE ARENA, with no scratch and no copy. C[pivot] is the
    // reach, so the block the walk fills is already the clique's own block; there was never a
    // reason for the set to exist anywhere else first.
    //
    // The comment this replaces weighed the copy against SWAPPING a scratch into place, and picked
    // the copy because a swapped-out scratch comes back empty and grows again from nothing at the
    // next pivot. That was right about swapping and it missed the third option. What made the
    // third option safe is the reserve above: the arena no longer doubles, so appending to it does
    // not pay what a scratch would have.
    //
    // Measured before: 111 ms for the copy and 59 ms in the push_backs' capacity checks, of an
    // 8.38 s run.
    //
    // AND THE ARENA MUST NOT MOVE WHILE THAT WALK RUNS, which is what the reserve below is for and
    // is not an optimization. `reachableSet` reads each clique's members through a pointer into
    // the arena, `mCliqueArena.data() + mCliquePtr[c]`, while appending the reach to that same
    // arena. A push_back that outgrows the capacity reallocates, and every such pointer already
    // taken is then dangling for the rest of its clique. The constructor's reserve is nnz(A) and
    // the arena grows to the sum of |C[p]| over the run, 108705 against 97440 at 140 a side, so a
    // reallocation is ordinary rather than exceptional.
    //
    // It read as harmless for as long as it did because a vector growth COPIES and then frees, so
    // the stale pointer usually still finds the right values sitting in freed memory. That is a
    // property of the allocator and not of the program: on Apple Silicon a 6x6x6 grid came out
    // with a different ordering from the same source on the same input, which is what an ordering
    // with no floating point in it cannot legitimately do. Address sanitizer reports it on every
    // one of the six drivers, on 2D grids as well as 3D, so it is the shared class's and not any
    // driver's.
    //
    // The remedy is to make the arena unable to move rather than to re-fetch per element, which
    // would put a load in the innermost loop of the whole ordering for a hazard that occurs once
    // per elimination at most. A reach is at most `size()` entries, so room for one is room for
    // the whole walk, and the growth stays geometric so nothing is given back to the doubling this
    // reserve exists to avoid.
    //
    // The rule is already this tree's, stated for the dynamic factor in DESIGN_DECISIONS and
    // rehearsed in experiments/storage-options: structural growth invalidates every pointer taken
    // before it. What made it easy to miss here is that the pointer and the growth are in
    // DIFFERENT functions, and that the comment two lines below already names the hazard for the
    // one pointer it happens to be about.
    // C[pivot] IS BUILT IN THE PIVOT'S OWN SEGMENT, which is dead from this moment: A[pivot] and
    // I[pivot] are read here and never again. When the reach outgrows it the list continues in the
    // segment of a clique the pivot absorbs, linked by a negative entry. This is mmdelm, ported.
    //
    // THE ONE INVARIANT THAT MAKES IT SAFE, and it is worth stating because it is not obvious: the
    // link written at the current segment's last entry always names the clique ABOUT TO BE WALKED.
    // So the write cursor can only ever enter a segment we are reading, never one still to be
    // read, and within that segment writes trail reads because every entry written was read first.
    const bool live    = mLiveMerges;
    const bool reverse = mReverseIncidence;

    const std::size_t   base = mSourcePtr[pivot];
    const std::uint32_t adjN = mAdjacencySize[pivot];
    const std::uint32_t incN = mIncidenceSize[pivot];

    // I[pivot] must be copied out before the segment holding it is overwritten. genmmd threads the
    // same list through its `list[]` array; a small vector kept for its capacity costs no
    // allocation after the first few pivots and keeps the walk order explicit.
    mAbsorbed.clear();
    for (std::uint32_t ii = 0; ii < incN; ++ii)
        mAbsorbed.push_back(mSource[base + adjN + (reverse ? incN - 1 - ii : ii)]);

    ++mTag;
    mMark[pivot] = mTag;                       // never its own neighbor

    std::size_t   rl    = base;                        // write cursor
    std::size_t   rm    = mSourcePtr[pivot + 1] - 1;   // last entry of the segment being filled
    std::uint32_t count = 0;

    // One member written, following a link first if the segment is full. The link can only be
    // there to be followed: it is written before the walk that can reach it.
    const auto emit = [&](std::int32_t v) {
        while (rl >= rm) {
            const std::int32_t l = mSource[rm];         // -(c + 1); see mCliqueSize's note
            rl = mSourcePtr[-l - 1];
            rm = mSourcePtr[-l] - 1;
        }
        mSource[rl++] = v;
        ++count;
    };
    // THE LIVE NEIGHBORS FIRST, AND WITHOUT THE BOUND CHECK. They are a subset of A[pivot], which
    // was read from this same segment, so the cursor cannot pass the reader and cannot leave the
    // segment: at most adjN entries are written into a segment holding adjN + incN. genmmd's
    // first loop is unchecked for exactly this reason, and checking here is not merely wasteful,
    // it is wrong -- with incN == 0 the cursor legitimately lands on the last entry, and a check
    // would read it as a link that was never written. Found by ASan on a 2 by 2 grid.
    for (std::uint32_t k = 0; k < adjN; ++k) {
        const std::int32_t v = mSource[base + k];
        if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) {
            mMark[v] = mTag;
            mSource[rl++] = v;
            ++count;
        }
    }

    // Then the cliques, each writing its continuation before it can be reached. This loop is the
    // one that can outgrow the segment, and it is the only one that follows a link.
    for (const std::int32_t c : mAbsorbed) {
        mSource[rm] = -(c + 1);                     // the continuation, before it can be read
        forEachMember(c, [&](std::int32_t v) {
            if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) { mMark[v] = mTag; emit(v); }
        });
    }

    mCliqueSize[pivot] = count;

    // The absorbed cliques die only now: they were being READ until the loop above finished.
    for (const std::int32_t c : mAbsorbed) mCliqueSize[c] = 0;

    // Stamp the new clique once and the absorbed cliques once, each with its own tag. Both are
    // then queried for free by the two compaction passes below. Clique ids and vertex ids share
    // one space, so one mark array serves both, the tags keeping them apart.
    ++mTag;
    inClique = mTag;
    std::uint32_t cliqueWeight = 0;
    forEachMember(pivot, [&](std::int32_t v) {
        mMark[v] = inClique;
        cliqueWeight += mWeight[v];
    });
    mCliqueWeight = cliqueWeight;
    ++mTag;
    absorbed = mTag;
    for (const std::int32_t c : mAbsorbed) mMark[c] = absorbed;
}

const std::vector<std::int32_t>& QuotientGraphB::eliminate(std::int32_t pivot) {
    std::int32_t inClique = NIL;
    std::int32_t absorbed = NIL;
    beginElimination(pivot, inClique, absorbed);
    // C[pivot] IS the reach here: beginElimination wrote it there and nothing has trimmed it yet
    // (massEliminate does, and runs after). So the prune walks the arena block rather than a copy.
    const bool amdOrder = mVendoredListOrder;      // hoisted, as the other flags are

    // Both lists are compacted in place rather than rebuilt into a scratch and swapped. Every
    // pass here only ever removes, so the survivors can be written over the entries already read,
    // the write cursor trailing the read one, and nothing is allocated at all. Rebuilding into a
    // shared scratch and swapping had each list inherit some other vertex's buffer, which then
    // had to grow again; that idiom is right where a pass can add, and this one cannot.
    // Both lists live in u's one run and are rewritten front to back, the adjacency first and the
    // incidence into whatever the adjacency has just given up. The pivot is appended rather than
    // pushed: the run never has to grow, because a source is destroyed here for each one created,
    // which is the conservation argument on the members. The two cursors are what makes the second
    // pass safe as well as the first, since the incidence write starts where the compacted
    // adjacency ends and its read starts where the original adjacency ended, which is never lower.
    // Under mVendoredListOrder both compactions below place their FIRST survivor last and the
    // incidence writes the pivot before the rest, which yields AMD_2's order without a second
    // pass over anything. Holding one entry back costs a register; doing it afterwards costs a
    // rotate per list per reached vertex, which is a whole extra walk of the structure and was
    // measured at about 50 percent of AMD3's run time before this was folded in. AMD_2 spends
    // three assignments on it for the same reason, letting the list's start shift rather than
    // moving a list. See experiments/ordering/AMD3.md, ledger entry 5 and iteration 10.
    forEachMember(pivot, [&](std::int32_t u) {
        std::int32_t*     source        = mSource.data() + mSourcePtr[u];
        // The two counters are one-dimensional COUNTS, positions in a list bounded by deg(u) and
        // so by n, where std::size_t is for a position into an n x n object. They take the type of
        // what they count, so they move with the array: `std::uint32_t`, and the cast that used to
        // stand here is gone. It was there only because the array was wider than the loop. The
        // dimensional rule in experiments/ordering/REPORT.md asks for this everywhere and it was
        // taken here first because this is the hottest loop in the ordering, Instruments put
        // 277 ms of an 8.53 s run on the incidence loop's header.
        //
        // `heldVertex` is a VERTEX and stays std::int32_t, carrying NIL.
        const std::uint32_t adjacencySize = mAdjacencySize[u];
        std::uint32_t       kept          = 0;
        std::int32_t        heldVertex    = NIL;        // the first survivor, appended last
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (v == pivot) continue;                  // no longer a variable
            if (mMark[v] == inClique) continue;        // both ends inside the new clique
            if (mLiveMerges && mEliminated[v] != 0) continue;   // numbered by a prepass, gone for good
            if (amdOrder && heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mAdjacencySize[u] = kept;                      // A[u] - C[pivot] - {pivot}

        // The read cursor is a hoisted POINTER, not source[adjacencySize + i]. The base and the
        // offset were being added per iteration on a loop whose body is one compare and one
        // conditional store, and Instruments put 277 ms of an 8.53 s run on this loop's header
        // against 100 ms on its body: the largest single line in the ordering. reachableSet has
        // hoisted the same pointer for the same reason since it was written.
        //
        // What CANNOT be hoisted away, and is the rest of that header: the read and the write are
        // into the same buffer, so every conditional store orders the next load behind it and
        // nothing crosses it. Amd.cpp compacts in place too, `Iw[pn++] = e` while reading `Iw[p]`,
        // and pays the same; it is the price of not allocating a scratch, which the note above
        // explains is the right trade here.
        const std::int32_t* incidence     = source + adjacencySize;
        const std::uint32_t incidenceSize = mIncidenceSize[u];
        std::uint32_t       write         = kept;      // follows `kept`; see the note above
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
            if (mMark[c] != absorbed) source[write++] = c;
        }
        source[write++]   = pivot;                     // u joins the new clique, id = pivot
        mIncidenceSize[u] = write - kept;
        // [c1, ..., ck, pivot] to [pivot, c2, ..., ck, c1], which is a SWAP of the two boundary
        // entries: only they move and c2..ck stay put. The pivot cannot be written first, which
        // is what the loop above would otherwise allow: the write cursor starts at `kept` and the
        // read at the original adjacencySize, and those are equal whenever nothing was pruned from
        // A[u], so an extra write before the reads finish clobbers an unread entry. AMD_2 makes
        // its three assignments after both compactions for exactly this reason.
        if (amdOrder && write - kept > 1) std::swap(source[kept], source[write - 1]);
    });

    return finishElimination(pivot);
}

// The same prune, with the driver's first scan folded into the two loops. The header carries the
// argument for why this is one call and why it cannot be folded further; the loops below are the
// plain ones with an accumulation added on each survivor, and nothing else differs.
// The same prune again, with the driver's first scan folded in under `Amd.cpp`'s tagged-W
// encoding. This is the overload above with `outside`, `mark` and `tag` replaced by one array and
// one tag, plus the hash key's adjacency half, which a driver that folds its scan in here has no
// other walk of A[u] to accumulate. Everything the header says about why the fold is sound holds
// unchanged: the scan runs over the untrimmed C[p], and the weights over A[u] are read before mass
// elimination could move any of them.
const std::vector<std::int32_t>& QuotientGraphB::finishElimination(std::int32_t pivot) {
    // Under mLateMassElimination the merge is the caller's, run after it has absorbed, so this
    // hands back an empty list and C[pivot] stays reach(pivot) exactly. See the setter.
    if (mLateMassElimination) {
        mMerged.clear();
    } else {
        massEliminate(pivot);
    }

    mAdjacencySize[pivot] = 0;
    mIncidenceSize[pivot] = 0;
    mEliminated[pivot]    = 1;
    return mMerged;
}

// Mass elimination. u is indistinguishable from the pivot when the two had the same closed
// neighborhood before the step, equivalently when everything u can still reach now lies inside
// the new clique, and eliminating it next then creates no fill at all. The test is a cheap
// sufficient condition for that: nothing explicit left and no clique but the new one. It is
// conservative, and deliberately so; the exact test costs a reachability query per candidate. See
// the mass-elimination section of experiments/ordering/README.md.
//
// It runs from finishElimination by default and from the driver under mLateMassElimination, and
// the body is the same either way: what moves is when the question is asked, since aggressive
// absorption is what makes this cheap test agree with the true one. experiments/ordering/AMD3.md, entry 3.
const std::vector<std::int32_t>& QuotientGraphB::massEliminate(std::int32_t pivot) {
    std::vector<std::int32_t>& merged = mMerged;   // scratch, kept for its capacity
    merged.clear();
    // Walks C[pivot], which is still the full reach: the trim below is this function's own and
    // happens after the loop.
    forEachMember(pivot, [&](std::int32_t u) {
        // Under mVendoredListOrder the new clique sits at the FRONT of I[u] rather than the back,
        // so the single remaining entry is at the head of the incidence run either way: with A[u]
        // empty the run starts with I[u], and with one element there is only one position. The
        // test therefore needs no branch on the flag.
        if (mAdjacencySize[u] == 0 && mIncidenceSize[u] == 1 &&
            mSource[mSourcePtr[u]] == pivot) {         // A[u] empty, so I[u] starts at the run
            mIncidenceSize[u] = 0;
            mEliminated[u]    = 1;
            merged.push_back(u);
        }
    });
    if (!merged.empty()) {                             // C[pivot] - merged, one compaction pass
        ++mTag;
        for (std::int32_t u : merged) mMark[u] = mTag;
        // The compaction follows the chain on BOTH cursors. Only removals happen, so the write
        // cursor trails the read one within a segment; when either meets a link it follows it,
        // and a link is structural and stays where it is. The two can be in different segments,
        // which is why each carries its own position.
        std::size_t   rp = mSourcePtr[pivot], wp = mSourcePtr[pivot];
        std::uint32_t left = mCliqueSize[pivot], kept = 0;
        while (left != 0) {
            const std::int32_t v = mSource[rp];
            if (v < 0) { rp = mSourcePtr[-v - 1]; continue; }
            --left; ++rp;
            if (mMark[v] == mTag) continue;             // merged away, drop it
            while (mSource[wp] < 0 && wp != rp) wp = mSourcePtr[-mSource[wp] - 1];
            mSource[wp++] = v;
            ++kept;
        }
        mCliqueSize[pivot] = kept;

        for (std::int32_t u : merged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            // The weighted clique size follows the clique. `cliqueWeight()` promises the weighted
            // size of C[pivot] AS IT NOW STANDS, and a merged vertex has just left it, so the
            // decrement belongs here rather than in the caller. AMD_2 spells the same line
            // `degme -= nvi` inside its own mass elimination.
            //
            // WITHOUT IT the drivers that mass-eliminate inside the eliminator, Amd1, Amd2 and
            // Amd2B, read the UNTRIMMED size where they had been computing the trimmed one, which
            // is a bound too large per vertex the merge took. That is the same shape as ledger
            // entry 7 and it was caught by `prototype and production agree` in
            // experiments/ordering, with `make amdorder` and all 283 assertions passing: Amd3
            // mass-eliminates late, so its own first read is legitimately of the untrimmed clique
            // and every check that watches Amd3 stayed green.
            mCliqueWeight -= mWeight[u];
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }
    }
    return merged;
}

void QuotientGraphB::merge(std::int32_t u, std::int32_t v) {
    mLiveMerges = true;                            // the reachable set must filter from now on
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    mAdjacencySize[v] = 0;
    mIncidenceSize[v] = 0;
    mEliminated[v]    = 1;
}

// genmmd's mmdnum numbering: each pivot first, then the members of its supervariable by
// ASCENDING VERTEX INDEX rather than by the order they were merged. The members of a
// supervariable are indistinguishable by construction, so this cannot change the fill or the
// elimination forest; it changes only which permutation comes out, and it is here so that a
// comparison against the vendored routine is an equality test rather than a judgement.
//
// Two passes and ONE scratch array. The first walks the pivots, giving each one its base slot,
// and threads its chain marking every member with the root it belongs to. The second scans the
// vertices ASCENDING and drops each member into its root's next slot, so ascending order falls
// out of the scan and nothing is sorted, which is how mmdnum gets it too.
//
// `slot` carries two meanings, told apart by sign, so one array does the work of two. For a ROOT
// it holds the next free position after that root, a non-negative index; for a MEMBER it holds
// `-(root + 1)`, always negative since root is non-negative. A vertex is a root or a member and
// never both, so the two never collide. That is the same trick genmmd plays on `perm`, which
// holds a number for a numbered vertex and a negated parent for a merged one.
//
// The obvious version, written first, allocated four arrays of size n and made six passes; it
// cost 244 ms of a 4.94 s profile where mmdint and mmdnum together cost 116 ms, which was the
// whole reason to come back to it. See experiments/ordering/mmd3.py, ledger entry 6.
std::vector<std::int32_t> QuotientGraphB::orderAscending(
        const std::vector<std::int32_t>& pivots) const {
    const std::size_t n = size();
    std::vector<std::int32_t> order(n);
    std::vector<std::int32_t> slot(n, 0);

    std::size_t pos = 0;
    for (std::int32_t pivot : pivots) {
        // The members first, so marking them cannot overwrite the root's own cursor: the chain
        // starts AT the pivot, and the pivot is a root rather than a member of itself.
        for (std::int32_t u = mSuperNext[pivot]; u != NIL; u = mSuperNext[u]) slot[u] = -(pivot + 1);
        order[pos]  = pivot;
        slot[pivot] = static_cast<std::int32_t>(pos) + 1;   // where its first member goes
        pos += mWeight[pivot];                              // the whole supervariable's room
    }

    for (std::size_t v = 0; v < n; ++v) {                   // ascending, so the members are too
        const std::int32_t s = slot[v];
        if (s >= 0) continue;                               // a root, already placed
        const std::int32_t root = -s - 1;
        order[static_cast<std::size_t>(slot[root]++)] = static_cast<std::int32_t>(v);
    }
    return order;
}

}  // anonymous namespace


std::vector<std::int32_t> orderMmd3B(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    std::int32_t delta) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraphB qg(colPtr, rowIdx);
    // The fourth walk, and the deepest: it fixes the order of C[pivot], hence the content order of
    // every list built from it. The other three are below. All four mirror genmmd holding a linked
    // list pushed at the head and read from the head, `list[nb] = h; h = nb`, where we append to a
    // vector. Same sets, same cost, different winner among equals, and minimum degree is settled by
    // exactly that. See experiments/ordering/mmd3.py.
    qg.setReverseIncidence(true);
    std::vector<std::int32_t> pivots;
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;

    // The filed value, not the degree. mmdint files a degree-0 vertex under 1, `if(dg==0)dg=1`,
    // and the refresh below files under degree + 1, so from here on this array holds what MMD
    // compares and files by rather than the degree itself.
    std::vector<std::uint32_t> degrees(size);
    BucketsB buckets(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        degrees[u] = std::max<std::uint32_t>(qg.adjacencySize(u), 1);
        buckets.file(degrees[u], u);
    }
    std::uint32_t minDegree = *std::min_element(degrees.begin(), degrees.end());

    // Withheld from the buckets rather than refiled, which is genmmd's bwd[nd] = -maxint. Not the
    // same as merged: the vertex is still live and still reachable, it simply cannot be the
    // minimum before the vertex that outmatched it. An elimination that reaches it puts it back.
    // A byte per vertex, not std::vector<bool>; see the note on BucketsB::mFiled. Constructed once
    // per ordering, which is where that specialization costs.
    std::vector<std::uint8_t> outmatched(size, 0);

    std::vector<std::int32_t> touchedRound(size, NIL);
    std::vector<std::int32_t> touched, batch, elementMembers, q2h, qxh, refreshed;
    std::int32_t numRounds = 0;

    // The driver's own membership scratch. Two levels, as mmdupd has: one tag survives a whole
    // element and says "already counted in dg0", the other is fresh per vertex, so one q2h vertex
    // cannot hide a neighbor from the next.
    std::vector<std::int32_t> mark(size, NIL);
    std::int32_t              tag = 0;

    // ---- the prepass ------------------------------------------------------------
    // Bucket 1 holds the isolated and the degree-1 vertices together, by the convention above.
    // Number them and leave the bucket empty. Nothing is eliminated in the quotient-graph sense.
    for (std::int32_t u = buckets.head(1); u != NIL; u = buckets.next(u)) touched.push_back(u);
    for (std::int32_t u : touched) {
        buckets.unfile(degrees[u], u);
        qg.number(u);
        pivots.push_back(u);
        ++numEliminated;
    }
    if (size > 2) minDegree = 2;                // head[1] is empty now, and mdeg starts at 2

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;

        // ---- one batch, no degree refreshed inside it ---------------------------
        std::uint32_t batchLimit = minDegree;
        if (delta > 0)
            batchLimit = std::min(minDegree + static_cast<std::uint32_t>(delta),
                                  static_cast<std::uint32_t>(size) - 1);

        batch.clear();
        touched.clear();
        while (true) {
            if (buckets.empty(minDegree)) {
                if (minDegree >= batchLimit) break;
                ++minDegree;
                continue;
            }
            const std::int32_t pivot = buckets.head(minDegree);
            buckets.unfile(degrees[pivot], pivot);

            const std::vector<std::int32_t>& merged = qg.eliminate(pivot);
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
            for (std::int32_t u : merged) {
                buckets.unfile(degrees[u], u);
                degrees[u] = 0;
            }
            degrees[pivot] = 0;

            qg.forEachMember(pivot, [&](std::int32_t u) {
                outmatched[u] = 0;                  // mmdelm's bwd[rn] = 0, back in the running
                buckets.unfile(degrees[u], u);      // evict, with a stale degree
                if (touchedRound[u] != numRounds) {
                    touchedRound[u] = numRounds;
                    touched.push_back(u);
                }
            });

            if (numEliminated >= size) break;       // genmmd: nothing left to update
            if (delta < 0) break;
        }

        // ---- one refresh, walked element by element -----------------------------
        refreshed.clear();
        // The driver's element list, genmmd's `list[mn] = ehead; ehead = mn`, so the LAST pivot of
        // a batch is the FIRST element refreshed.
        for (auto ee = batch.rbegin(); ee != batch.rend(); ++ee) {
            const std::int32_t element = *ee;
            elementMembers.clear();
            qg.forEachMember(element, [&](std::int32_t v) {
                if (!qg.eliminated(v)) elementMembers.push_back(v);
            });

            ++tag;
            const std::int32_t elementTag = tag;    // dg0's members, marked once for the element
            for (std::int32_t v : elementMembers) mark[v] = elementTag;
            std::uint32_t dg0 = 0;
            for (std::int32_t v : elementMembers) dg0 += qg.weight(v);

            // reach(u) has |A[u]| + |I[u]| sources once the new element is counted, so one other
            // source means everything u reaches is in this element plus that one place. dg0
            // already counts the element, and the other source is walked directly, so no union is
            // formed at all.
            q2h.clear();
            qxh.clear();
            for (std::int32_t u : elementMembers) {
                if (buckets.filed(u) || outmatched[u] != 0) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? q2h : qxh).push_back(u);
            }

            // mmdupd's q2h list, `list[nb] = q2h; q2h = nb`.
            for (auto uu = q2h.rbegin(); uu != q2h.rend(); ++uu) {
                const std::int32_t u = *uu;
                if (qg.eliminated(u) || outmatched[u] != 0) continue;   // by an earlier q2h vertex
                ++tag;
                const std::int32_t vertexTag = tag;
                // dg0 is kept WHOLE and u's own weight subtracted at the end, which is
                // genmmd's `dg - qsize[en] + 1` and not the same as subtracting it now. The
                // walk below can MERGE a vertex into u, and genmmd's merge does
                // `qsize[en] += qsize[nd]` in that same walk, so the weight it subtracts is the
                // one AFTER the merge. Subtracting first files a supervariable one bucket too
                // high per merged vertex, so it is not picked as early as its size has earned.
                // See experiments/ordering/mmd3.py, ledger entry 5.
                std::uint32_t degree = dg0;

                // Not hoisted, deliberately. A q2h vertex has adjacencySize + incidenceSize == 2
                // by the test that put it on this list, so these two loops run over at most two
                // elements between them and a length loaded up front is overhead rather than a
                // saving. Hoist where a loop is long; leave it where the loop is short or exits
                // early. Measured both ways.
                const std::int32_t* adjacency = qg.adjacency(u);
                for (std::uint32_t a = 0; a < qg.adjacencySize(u); ++a) {
                    const std::int32_t v = adjacency[a];
                    if (qg.eliminated(v) || mark[v] == vertexTag) continue;
                    if (mark[v] == elementTag) continue;           // already counted in dg0
                    mark[v] = vertexTag;
                    degree += qg.weight(v);
                }
                const std::int32_t* incidence = qg.incidence(u);
                for (std::uint32_t i = 0; i < qg.incidenceSize(u); ++i) {
                    const std::int32_t c = incidence[i];
                    if (c == element) continue;
                    qg.forEachMember(c, [&](std::int32_t v) {
                        if (v == u || qg.eliminated(v) || mark[v] == vertexTag) return;
                        if (mark[v] == elementTag) {
                            // v is in the new element and in this same other source, so it sees
                            // at least what u sees.
                            if (buckets.filed(v) || outmatched[v] != 0) return;
                            if (qg.adjacencySize(v) + qg.incidenceSize(v) - 1 == 1) {
                                qg.merge(u, v);      // identical reach: u absorbs it
                                ++numEliminated;
                            } else {
                                outmatched[v] = 1;      // reaches more, so never minimal first
                            }
                            return;
                        }
                        mark[v] = vertexTag;
                        degree += qg.weight(v);
                    });
                }

                degrees[u] = std::max<std::uint32_t>(degree - qg.weight(u) + 1, 1);
                buckets.file(degrees[u], u);
                refreshed.push_back(u);
            }

            // mmdupd's qxh list, the same stack.
            for (auto uu = qxh.rbegin(); uu != qxh.rend(); ++uu) {
                const std::int32_t u = *uu;                 // the full union, as md5 computes it
                if (qg.eliminated(u) || outmatched[u] != 0) continue;
                const std::uint32_t degree = qg.reachableWeight(u); // reach excludes u already
                degrees[u] = std::max<std::uint32_t>(degree + 1, 1);
                buckets.file(degrees[u], u);
                refreshed.push_back(u);
            }
        }

        for (std::int32_t u : refreshed) minDegree = std::min(minDegree, degrees[u]);
        ++numRounds;
    }

    return qg.orderAscending(pivots);   // genmmd's mmdnum. See the ledger, entry 6.
}


}  // namespace Oblio
