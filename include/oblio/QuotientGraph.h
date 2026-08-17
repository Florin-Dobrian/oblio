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
    // GENMMD'S `bwd` ENCODING. `mPrev[u]` carries four facts at once, told apart by sign, which
    // is what lets three arrays be one. `private/Mmd.cpp` lines 81 to 84 and 98 are the original.
    //
    //     mPrev[u] >  0            filed; the vertex before it is mPrev[u] - 1
    //     mPrev[u] <  0            filed AT THE HEAD of bucket -mPrev[u] - 1
    //     mPrev[u] == UNFILED      not in any list
    //     mPrev[u] == OUTMATCHED   withheld from the buckets; see outmatch()
    //
    //
    // THE HEAD CASE IS WHAT DELETES THE DEGREE ARRAY. `unfile` used to be given the degree,
    // because a vertex at the head of its list left no record of which list that was; it now
    // reads -mPrev[u] and takes no degree at all. That is exactly why genmmd carries no degree
    // array: the mmd drivers here kept one solely to answer that question, and the amd drivers
    // keep theirs because they READ the value in the bound, which is a different use.
    //
    // BOTH SIDES ARE SHIFTED BY ONE AND BOTH SHIFTS ARE LOAD BEARING. A predecessor is stored as
    // `u + 1` because vertex 0 would otherwise read as UNFILED, and a head as `-(degree + 1)`
    // because DEGREE 0 IS REACHABLE and `-0` would read as UNFILED too. genmmd needs neither: its
    // ids are 1-based, and `mmdint` files a degree-0 vertex under 1, `if(dg==0)dg=1`, so no key
    // of zero ever reaches its buckets. Ours do: Mmd1 and all five amd drivers file
    // `adjacencySize(u)` raw, which is 0 for an isolated vertex. Dropping the second shift builds
    // and passes every 2D case, then corrupts a bucket list on a 3D grid at 6 a side.
    static constexpr std::int32_t UNFILED    = 0;
    static constexpr std::int32_t OUTMATCHED = -2147483647 - 1;   // INT32_MIN; no degree reaches it

    explicit Buckets(std::size_t size)
        : mHead(size + 1, NIL), mNext(size, NIL), mPrev(size, UNFILED) {}

    // buckets[degree].add(u), at the head. The head is the only O(1) end of a singly reachable
    // list, so the winner among equal degrees is whatever was filed last rather than the lowest
    // index. That is the vendored convention and it is why an ordering differs from an exact
    // scan's in its ties.
    void file(std::uint32_t degree, std::int32_t u) {
        mNext[u] = mHead[degree];
        mPrev[u] = -static_cast<std::int32_t>(degree) - 1;   // head of `degree`; see the encoding  // at the head, and this is its bucket
        if (mHead[degree] != NIL) mPrev[mHead[degree]] = u + 1;
        mHead[degree] = u;
    }

    // buckets[degree].discard(u). Idempotent, which matters during a batch: a vertex evicted
    // early can be merged away by a later pivot in the same round, and unfiling it twice must
    // not splice a list it is no longer in.
    void unfile(std::int32_t u) {
        const std::int32_t prev = mPrev[u];
        if (prev == UNFILED || prev == OUTMATCHED) return;   // not in a list; see the encoding
        if (prev > 0) mNext[prev - 1]   = mNext[u];
        else          mHead[-prev - 1] = mNext[u];           // u headed bucket -prev - 1
        if (mNext[u] != NIL) mPrev[mNext[u]] = prev;
        mNext[u] = NIL;
        mPrev[u] = UNFILED;
    }

    // Withhold u from the buckets without filing it anywhere, which is genmmd's
    // `bwd[nd] = -maxint`. It stays live and reachable; it simply cannot be the minimum before
    // the vertex that outmatched it, so it is not a candidate until an elimination reaches it and
    // `restore` puts it back. mmdelm spells that restore `bwd[rn] = 0`, the same store that
    // unfiles, which is why the two sit together at every call site here.
    void outmatch(std::int32_t u)         { unfile(u); mPrev[u] = OUTMATCHED; }
    void restore(std::int32_t u)          { if (mPrev[u] == OUTMATCHED) mPrev[u] = UNFILED; }
    bool outmatched(std::int32_t u) const { return mPrev[u] == OUTMATCHED; }

    // Move u to the bucket for newDegree, carrying the cached degree with it. The three steps go
    // together, which is why this is one call: the bucket a vertex sits in is read from its
    // degree, so writing the degree first would erase it from the wrong list. A vertex whose
    // degree did not change is removed and reinserted into the same list, which is harmless.
    void refile(std::vector<std::uint32_t>& degrees, std::int32_t u, std::uint32_t newDegree) {
        unfile(u);
        degrees[u] = newDegree;
        file(newDegree, u);
    }

    // The next vertex in the same bucket, and whether u is filed at all. MMD reads both: it
    // walks a whole bucket in the prepass, and its refresh asks whether a vertex it reached has
    // already been dealt with this round.
    // AMD_2'S TWO BORROWED SLOTS. A vertex that has been unfiled is in no degree list, so its
    // predecessor and successor links carry nothing, and AMD_2 uses exactly that: `Last[i]` holds
    // the hash key and `Next[i]` the hash chain for the whole middle of an elimination step, which
    // is why it allocates neither a key array nor an Hhead.
    //
    // LEGAL ONLY BETWEEN unfile() AND file(). A key stored here is an arbitrary int32 and can look
    // like any of the encodings mPrev carries, so `unfile()`, `filed()` and `outmatched()` MUST
    // NOT be called on a vertex holding one; they would splice a list on garbage. Only the amd
    // path uses these, and only after taking every member of C[pivot] out of the lists. No mmd
    // driver may call them: mmd leaves its candidates filed and asks `filed()` about them.
    void         setKey(std::int32_t u, std::int32_t k)   { mPrev[u] = k; }
    std::int32_t key(std::int32_t u) const                { return mPrev[u]; }
    void         setChain(std::int32_t u, std::int32_t v) { mNext[u] = v; }
    std::int32_t chain(std::int32_t u) const              { return mNext[u]; }

    std::int32_t next(std::int32_t u) const      { return mNext[u]; }
    bool         filed(std::int32_t u) const     { return mPrev[u] != UNFILED
                                                          && mPrev[u] != OUTMATCHED; }

    std::int32_t head(std::uint32_t degree) const  { return mHead[degree]; }
    bool         empty(std::uint32_t degree) const { return mHead[degree] == NIL; }

private:
    std::vector<std::int32_t> mHead;   // mHead[d], the first live vertex of degree d
    std::vector<std::int32_t> mNext;   // mNext[u], toward the tail
    std::vector<std::int32_t> mPrev;   // mPrev[u], toward the head
    // A byte per vertex, not std::vector<bool>. That specialization packs one bit per element, and
    // what it costs here is CONSTRUCTION rather than access: a vector of n false is built through
    // the bit-reference machinery word by word, where a vector of n zero bytes is a memset. Both
    // this flag and the graph's own arrays are built once per ordering, so a workload that
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
    // THE BUCKETS, OR NOT. A driver that wants AMD_2's hash arrangement passes them: the prune
    // then takes every member of C[pivot] out of the degree lists and parks the adjacency half of
    // the hash key in the predecessor link it has just freed. A driver that refiles inside its own
    // bound pass, which is Amd2 and Amd2B, cannot have either: its links are still degree links
    // when the hash runs. Those pass NULL and build their key in a pass of their own.
    //
    // Null therefore means "leave the degree lists alone and store no key". It does not change
    // what the scan computes, only where the by-products go.
    Buckets*                          buckets;
    std::vector<std::int32_t>&        w;             // per clique, Amd.cpp's tagged W
    // Amd.cpp's `Degree`, which serves a LIVE vertex's degree and a DEAD one's clique weight from
    // one array, the two being disjoint because a clique id is the id of the pivot that formed it.
    // The scan reads only the clique half.
    const std::vector<std::uint32_t>& degree;
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
    // THE MARK'S WIDTH IS A CONSTRUCTOR ARGUMENT, and it is a sizing decision rather than a
    // behaviour flag, which is why it is here and not a setter beside the other three.
    //
    // With `cliqueMarks` off, `mMark` is n and indexes vertices alone, which is all any mmd driver
    // needs: genmmd stamps vertices and never cliques, and its refresh reaches a clique's members
    // through the arena rather than testing the clique's own identity. With it on, `mMark` is 2n,
    // vertices low and cliques at `cliqueBase() + c`, which is what the amd branch needs because
    // supervariable detection tests I[u] == I[v], set equality over CLIQUE IDS.
    //
    // The two halves cannot share one space. A clique id IS the id of the dead pivot that formed
    // it, so stamping a clique in the vertex half would write a live tag over a slot holding GONE,
    // and a walk of an older clique still listing that pivot as a member would read it as live.
    // Neither genmmd nor AMD_2 shares one array between the two kinds.
    //
    // It is an argument rather than always-2n so that the mmd drivers do not allocate and zero an
    // n int32 they never read: the constructor cost of this class is not incidental, and the note
    // beside Buckets' flag byte records a 12 percent ordering-time saving that was almost entirely
    // construction.
    QuotientGraph(const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  bool cliqueMarks = false);

    std::size_t size() const           { return mRun.size(); }

    // HOW MANY ENTRIES THE CLIQUE ARENA HOLDS, for the drivers that report it. This is a SIZE and
    // not a capacity: what was written, not what was reserved. Our arena only grows, so it is also
    // the peak; a scheme that COMPACTS would have a final size below its high-water mark, so if
    // Mmd3B or Amd3B ever answer the same question they must answer it with the peak or the column
    // will be comparing three different things.
    std::size_t arenaEntries() const   { return mCliqueArena.size(); }
    // GONE, which is genmmd's `marker[v] = maxint` exactly: one value reserved above every tag
    // makes the stamp array answer "is v dead" on the load it was making anyway, so no array is
    // spent on liveness at any walk site. AMD_2 does the same with `W[e] = 0` for a dead clique.
    //
    // NOT the weight. `mWeight[v] != 0` is a PARTIAL flag on both sides: `number()` leaves a
    // prepass vertex at weight one deliberately, so its neighbors' degrees still count it, and
    // genmmd's prepass leaves `qsize` at one for the same reason. genmmd uses `qsize[nd] != 0`
    // only inside element walks, where a prepass vertex cannot appear because the mark has
    // already kept it out of every clique. Used as the universal test it lets a numbered vertex
    // back into a reachable set: 201 entries for 200 vertices on a random mmd2 pattern, tried and
    // reverted on 2026-08-08.
    //
    // Safe only because clique ids no longer share this array; see beginElimination.
    static constexpr std::int32_t GONE = 2147483647;   // INT32_MAX, above every reachable tag

    bool eliminated(std::int32_t u) const { return mMark[u] == GONE; }

    // A driver may stamp into this same array rather than allocating one of its own, which is
    // what genmmd's `marker` is: `mmdelm` stamps it at level `tag` and `mmdupd` at level
    // `mt = tag + md0`, one array and one counter serving both. One counter is what makes it
    // safe, since two tags drawn from it can never be equal.
    // The stride between the two halves of `mMark`, and the one place a COUNT becomes an offset
    // in the INDEX space. Meaningful only when the graph was built with cliqueMarks.
    std::int32_t cliqueBase() const              { return static_cast<std::int32_t>(mRun.size()); }

    std::int32_t advanceTag()                    { return ++mTag; }
    std::int32_t mark(std::int32_t v) const      { return mMark[v]; }
    void setMark(std::int32_t v, std::int32_t t) { mMark[v] = t; }

    // The members of clique c, which after eliminate(p) is the pattern of p's column of L: the
    // vertices the pivot reached, less those it absorbed. A pointer and a length, as the
    // adjacency is, for the same reason.
    //
    // **Valid until the next eliminate.** The members live in one arena that grows as cliques are
    // formed, so a reallocation moves them; the offsets are indices and survive it, but a pointer
    // taken beforehand does not. Read a clique at the moment of use, which is the same rule the
    // numeric factor's blocks live by.
    const std::int32_t* clique(std::int32_t c) const {
        return mCliqueArena.data() + mRun[c].sourcePtr;
    }
    std::uint32_t cliqueSize(std::int32_t c) const { return mRun[c].adjacencySize; }

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
        return mSource.data() + mRun[u].sourcePtr;
    }
    std::uint32_t adjacencySize(std::int32_t u) const { return mRun[u].adjacencySize; }

    const std::int32_t* incidence(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr + mRun[u].adjacencySize;
    }
    std::uint32_t incidenceSize(std::int32_t u) const { return mRun[u].incidenceSize; }

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
    void beginElimination(std::int32_t pivot, std::int32_t& inClique);

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
    // ONE OBJECT PER VERTEX, NOT THREE ARRAYS. These three numbers are never useful apart: any
    // walk of u needs where its run starts and at least one of the two lengths, and the prune
    // needs all three. Held as three arrays they sit at three unrelated addresses, so a walk of
    // C[pivot], whose members are scattered vertex ids, pulls in three cache lines per member and
    // uses 4 or 8 bytes of each. Held together they are one line.
    //
    // Measured before the change, on a 100x100 grid: 9148 + 5611 + 5687 = 20446 D1 read misses on
    // the prune's three preamble lines, 15.8 per cent of the ordering's total, against genmmd's
    // 4070 for the same three facts. It gets them from `xadj[rn]` and `xadj[rn+1]`, ADJACENT
    // entries of one array, plus a field of `fwd[rn]` that the same loop reads anyway to unfile
    // the vertex. So this is not a scheme it lacks; it is a scheme we had split apart.
    //
    // EXACTLY 16 BYTES, four to a 64-byte line, and the layout only pays while that holds. The
    // position is `std::size_t` because it offsets into an arena bounded by nnz(A); the two
    // lengths are `std::uint32_t` because they are one dimensional and bounded by n. That is the
    // integer rule producing 8 + 4 + 4 with no padding, rather than a size chosen to fit.
    struct VertexRun {
        std::size_t   sourcePtr;       // where u's run starts in mSource, fixed at construction
        std::uint32_t adjacencySize;   // A[u]'s length, from the run's start
        std::uint32_t incidenceSize;   // I[u]'s length, immediately behind A[u]
    };
    std::vector<VertexRun> mRun;

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
    // A CLIQUE'S DESCRIPTOR IS THE DEAD PIVOT'S OWN RUN, 2026-08-15. `mRun[c].sourcePtr` is where
    // C[c] starts in the arena and `mRun[c].adjacencySize` is how much of it is still live. Two
    // arrays of size n, `mCliquePtr` and `mCliqueSize`, are gone.
    //
    // WHY IT IS SOUND. A clique id IS the id of the pivot that formed it, and that vertex's own
    // A[c] and I[c] are read for the last time inside `beginElimination` and never again, so the
    // three words describing its run are free from that moment. This is what `AMD_2` does with
    // `Pe` and `Len`: an element takes over the slots of the variable it came from, which is why
    // it allocates nothing for elements at all. The one crossing is still there, the arena's new
    // length less the block's start written as a count, and it is still in `beginElimination`.
    //
    // WHAT IT WAS WORTH, and it was not the storage. Before this, every clique visit in a walk
    // probed two separate n-arrays at a dead pivot's id, scattered across the whole id space, once
    // per element of every I[u]. Now one 16-byte run gives both, on a line the walk is already
    // touching. Measured on alpamayo with Amd3B: the amd branch's 2D ratio against the vendored
    // routine had been RISING with n, 1.11x at 32 a side to 1.49x at 400, and it went FLAT at
    // 1.38x from 100 a side up. A differential had already shown the two codes doing the same
    // visits per pivot at every size, so the growth was cost per visit, and this was most of it.
    //
    // THE ARENA ITSELF DOES NOT MOVE. Blocks still live in their own storage in elimination
    // order; only the two words describing them moved. `Mmd3B` priced genmmd's placement, cliques
    // in the pivot's dead segment, and it lost. See docs/DESIGN_DECISIONS.md.
    // APPENDED BY push_back, NOT written through a cursor, and that was measured. A cursor with
    // an explicit used count removes a capacity test and a size update per clique member, about
    // 46000 of each per ordering, which is what genmmd's `adjncy[rl] = nb; rl++` costs nothing
    // for. Built on 2026-08-15 and it came out 109085 instructions and 52880 reads WORSE, because
    // the vector's own length can then no longer be the arena's length: `reserve` does not touch
    // memory while `resize` value-initializes, so the constructor zeroes nnz(A) entries and every
    // growth zeroes its new region, which costs more than the tests it removes.
    // SEPARATE ALLOCATIONS ARE LOAD BEARING, 2026-08-16, and that was found by accident. `AMD_1`
    // carves Pe, Nv, Head, Elen, Degree and W out of ONE block at offsets that are exact multiples
    // of n. On a square grid n is m^2, so a power-of-two side gives a power-of-two n and those six
    // arrays land in the same cache sets at every index. Cachegrind on the vendored routine: the
    // instruction count and the data-read count per vertex are FLAT across 400, 512 and 800 a side
    // to a tenth of a percent, while D1 read misses per vertex go 15.3, 40.3, 17.0. Padding the six
    // apart by one cache line removes 56 percent of the misses at 512 and none at 400, with
    // byte-identical permutations.
    //
    // Every array in this class is its own vector, so nothing here has that property. Which is
    // worth stating where a later reader might be tempted to consolidate them for locality: the
    // consolidation is what creates the hazard. See docs/DESIGN_DECISIONS.md (2026-08-16, later).
    std::vector<std::int32_t>  mCliqueArena;  // every C[c] ever formed, end to end

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

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

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
