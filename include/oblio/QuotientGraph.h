#pragma once

// QuotientGraph.h - the representation Oblio's own minimum-degree orderings run on, and the
// degree buckets they pick from. Shared by Mmd1 and Amd1, which differ only in their drivers:
// one batches eliminations to make the degree refresh rare, the other bounds the degree to make
// each refresh cheap. Everything below that fork is here.
//
// The idea, in one line: an elimination does not create fill edges, it creates a CLIQUE, and a
// clique of d vertices is a d-clique list rather than d(d-1)/2 edges. So the neighbor relation
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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace Oblio {

// Written by every driver that tracks it; read by test_order. Defined with the bodies below,
// `inline` so that every unit including this header shares the one object.
extern std::size_t gPeakCliqueMembers;

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
    // A byte per vertex, not std::vector<bool>. That specialization packs one bit per clique, and
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
    // ONE MARK PER VERTEX, and no longer an argument. `mMark` was sized n or 2n on a constructor
    // flag: the wide form existed so amd's supervariable detection could stamp CLIQUE ids at
    // `cliqueBase() + c`, set equality over I[u] and I[v] needing somewhere to put them. Detection
    // now stamps into the driver's own tagged `W`, which is what AMD_2 does, so nothing asks for
    // the second half and `cliqueBase()` has gone with it. Removed 2026-08-17.
    QuotientGraph(const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const           { return mRun.size(); }

    // HOW MANY ENTRIES THE CLIQUE ARENA HOLDS, for the drivers that report it. This is a SIZE and
    // not a capacity: what was written, not what was reserved. Our arena only grows, so it is also
    // the peak; a scheme that COMPACTS would have a final size below its high-water mark, so if
    // Mmd3B or Amd3B ever answer the same question they must answer it with the peak or the column
    // will be comparing three different things.
    std::size_t arenaEntries() const   { return mCliqueArena.size(); }

    // MEMBERS, NOT ENTRIES, AND THAT IS THE WHOLE DISTINCTION FROM `arenaEntries` ABOVE. An entry
    // is a SLOT the arena has handed out and never takes back; a member is a vertex actually in a
    // live clique at this instant. Both are counts of `std::int32_t` and both compare against
    // nnz(A), so the words have to do the separating: `arenaEntries` is CUMULATIVE and grows by
    // every clique ever built, while these two rise and fall.
    //
    // THE PEAK IS THE NUMBER A DYNAMIC C WOULD NEED, holding each clique in its own allocation and
    // freeing it on death: the most it would ever ask for at one instant. It is the PAYLOAD figure
    // only. A real implementation also pays a header per clique, allocator rounding, and capacity
    // above size, none of which is counted here.
    //
    // A CLIQUE HAS THREE EVENTS AND THE COUNTER SEES ALL OF THEM: it is BORN once, in
    // `beginElimination`; it CONTRACTS, keeping its identity and losing members, at mass
    // elimination and again at supervariable detection; and it DIES, absorbed into a new clique or
    // absorbed aggressively once its external degree reaches zero. mmd has one of each; amd has
    // one birth, two contractions and two deaths. The full inventory, and why a contraction is not
    // a death, is in experiments/ordering/README.md (2026-08-18).
    //
    // IT IS A LIVENESS QUESTION, NOT A PLACEMENT ONE, which is why it is exact under every layout
    // in this tree and not only under this one: add on birth, subtract on both, keep the maximum.
    // The maximum is taken at BIRTH ALONE, deaths and shrinks being the only other events and
    // neither able to raise the total. A compacting layout's own high water mark is a different
    // and less useful number, being the array's rather than the live set's.
    std::size_t numPeakCliqueMembers() const { return mNumPeakCliqueMembers; }
    std::size_t numLiveCliqueMembers() const { return mNumLiveCliqueMembers; }

    // THE COUNTER CHECKED AGAINST A RECOMPUTATION, debug builds only. Births and deaths are spread
    // over four sites and nothing else in the suite would notice if they stopped balancing; the
    // fourth site, `massEliminate` shortening the pivot's own clique, was missed on the first
    // attempt and this is what found it.
    //
    // IT LIVES HERE RATHER THAN IN THE DRIVERS because only this class knows which vertices ever
    // formed a clique. Summing over a driver's pivot list looks equivalent and is not: `Mmd3`
    // pushes prepass vertices onto that list, and for those `cliqueSize` still reports A[p]'s
    // length. That was the second wrong version of this check.
    bool cliqueCountBalances() const;
    // GONE, which is genmmd's `marker[v] = maxint` exactly: one value reserved above every tag
    // makes the stamp array answer "is v dead" on the load it was making anyway, so no array is
    // spent on liveness at any walk site. AMD_2 does the same with `W[e] = 0` for a dead clique.
    //
    // NOT the weight. `mWeight[v] != 0` is a PARTIAL flag on both sides: `number()` leaves a
    // prepass vertex at weight one deliberately, so its neighbors' degrees still count it, and
    // genmmd's prepass leaves `qsize` at one for the same reason. genmmd uses `qsize[nd] != 0`
    // only inside clique walks, where a prepass vertex cannot appear because the mark has
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

    // WRITABLE, for the restore pass alone, which trims the clique as it walks it. See
    // `trimClique` in the .cpp for why that pass and not a pass of its own. Same lifetime rule as
    // the const overload above.
    std::int32_t* clique(std::int32_t c) { return mCliqueArena.data() + mRun[c].sourcePtr; }
    void trimClique(std::int32_t pivot, std::uint32_t kept);

    // THE ONE PLACE A CLIQUE DIES, and the reason it is one place is the counter above. Death has
    // three causes here: absorbed into the new clique, absorbed aggressively once its external
    // degree reaches zero, and merged away with its owner. All three used to write
    // `mRun[c].adjacencySize = 0` where they stood, which is correct and uncountable.
    void killClique(std::int32_t c);

    // |C[pivot]| WEIGHTED, which every amd driver needs and all four used to compute for
    // themselves in a pass of their own, one scattered weight load per member per pivot. It is
    // accumulated in `beginElimination`'s stamping walk instead, which has the member in hand
    // already, so the pass is gone and nothing is added. AMD_2 does the same, `degme += nvi`
    // inside the loops that build the clique, at its lines 1492 and 1636.
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
    // The MAGNITUDE, so that a driver never has to know whether the sign encoding is in force. It
    // is not, at this commit, and the cast is a widening of a value that is always positive.
    std::uint32_t weight(std::int32_t u) const {
        const std::int32_t w = mWeight[u];
        return static_cast<std::uint32_t>(w < 0 ? -w : w);
    }


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
    // clique's outside count. Run afterwards, that is a second and third visit to every clique,
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

    // SET u ASIDE, taking it out of the elimination without numbering it. `AMD_2`'s dense-row
    // rule: `Nv [i] = 0 ; Elen [i] = EMPTY ; nel++ ; Pe [i] = EMPTY`, a variable that is neither
    // eliminated nor available, kept out of every reachable set and every list by its ZERO WEIGHT
    // alone, which is the same mechanism a merged vertex leaves by. The caller owns where it
    // lands in the permutation; `AMD_2` appends the set at the end.
    void setAside(std::int32_t u);

    // Walk I[u] from the back in reachableSet(), matching genmmd's clique stack. A tie-break
    // convention and nothing else: it changes which permutation comes out, never which sets are
    // computed. See the member's note.
    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    // Lay the lists out the way AMD_2 does, which is the opposite of genmmd's on both counts.
    // Two conventions under one switch because they are one fact, that AMD's lists run the other
    // way round, and are only ever wanted together. Like the flag above, this changes which
    // permutation comes out and never which sets are computed. Used by Amd3 alone.
    //
    //   reachableSet   walks the CLIQUES before the explicit adjacency, since AMD_2's
    //                  `for (knt1 = 1; knt1 <= elenme + 1; knt1++)` takes the cliques of me and
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
    void beginElimination(std::int32_t pivot);

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
    // `Pe` and `Len`: a clique takes over the slots of the variable it came from, which is why
    // it allocates nothing for cliques at all. The one crossing is still there, the arena's new
    // length less the block's start written as a count, and it is still in `beginElimination`.
    //
    // WHAT IT WAS WORTH, and it was not the storage. Before this, every clique visit in a walk
    // probed two separate n-arrays at a dead pivot's id, scattered across the whole id space, once
    // per clique of every I[u]. Now one 16-byte run gives both, on a line the walk is already
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

    // See `numPeakCliqueMembers`. Maintained by `killClique`, by `trimClique`, and by the one
    // place a clique is born, which is why all three are funnelled rather than written where they
    // happen.
    std::size_t mNumLiveCliqueMembers = 0;
    std::size_t mNumPeakCliqueMembers = 0;
#ifndef NDEBUG
    std::vector<std::int32_t> mCliqueOwners;   // every vertex that ever formed one; see the check
#endif

    // The supervariable a vertex stands for, as a chain rather than a list per vertex. A list
    // meant one allocation per vertex before anything had happened, n of them for a structure
    // that is usually a singleton, plus a growth every time one absorbed another. The chain is
    // three arrays allocated once: the next member, the last one (so an absorption appends in
    // O(1) and the members keep their order, which the emitted permutation depends on), and the
    // count, which a chain no longer gives away for free.
    // `mWeight` is a one dimensional size, and the bound on an accumulation of weights is not the
    // term count but DISJOINTNESS: the weights partition the original vertices, so a sum over a set
    // of distinct vertices is at most n however many terms it has. That covers every accumulation
    // of weights in the ordering, in `reachableWeight`, in the clique weight below, and in the
    // drivers' bound and refresh passes. An accumulation over an unbounded number of one
    // dimensional terms would otherwise need a wider type, and this is the exception to that.
    //
    // AND IT IS `std::int32_t` RATHER THAN `std::uint32_t`, 2026-08-17, WHICH IS THE RULE DERIVING
    // AND NOT AN EXCEPTION TO IT. A one dimensional size is unsigned because it has nothing to
    // stand in for; this one is about to have. `AMD_2`'s `Nv` carries three facts in one field,
    // positive is the weight, negative means already taken into the clique being built, zero means
    // dead, so one load answers what two arrays answer here. The sign is a SPARE BIT rather than a
    // sentinel: a weight is bounded by n and n is capped at `MAX_IDX`, so no representable value is
    // given up, and the magnitude still carries the number. See docs/CODING_RULES.md, which states
    // the four conditions under which a size may go signed, and docs/DESIGN_DECISIONS.md
    // (2026-08-17).
    //
    // AT THIS COMMIT NOTHING NEGATES IT. The type is landed on its own so that the encoding, which
    // changes what several hot loops read, arrives against a tree where the width is already
    // settled and the digest has already said so.
    std::vector<std::int32_t>  mSuperNext;
    std::vector<std::int32_t>  mSuperLast;
    std::vector<std::int32_t>  mWeight;
    std::uint32_t              mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

    // Which end of I[u] reachableSet() walks from. genmmd threads its clique list through an
    // integer array and pushes at the head, `list[nb] = el; el = nb`, then reads from the head, so
    // the clique seen LAST is expanded FIRST; we hold a vector and append. Same set either way and
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

    // WHETHER ANY VERTEX HAS BEEN NUMBERED WITHOUT BEING ELIMINATED, which is what `number()` does
    // and only an mmd prepass does. False for a whole amd run.
    //
    // IT GUARDS A LOAD RATHER THAN A DECISION. Since the sign of the weight took over membership
    // and liveness, `mMark[v] != GONE` survives in the adjacency walk and in both prunes for one
    // case alone: a prepass vertex keeps a positive weight and sits in every neighbor's adjacency,
    // so a positive weight does not mean live THERE. Amd cannot produce such a vertex, so without
    // this flag its three drivers pay a second scattered load per adjacency clique for a
    // condition that can never hold. `mHasNumbered &&` short circuits before the load, and the
    // branch is perfectly predicted, being constant for the whole run.
    //
    // A FLAG AND NOT A CONSTRUCTOR ARGUMENT, deliberately: it is a fact about what has happened,
    // discovered by `number()` being called, not a mode a caller selects. Nothing has to be kept
    // in step, and a driver that gains a prepass later gets the correct behavior for free.
    bool mHasNumbered = false;

    std::vector<std::int32_t> mMark;     // membership scratch, read against mTag
    std::int32_t              mTag = 0;
};

// ------------------------------------------------------------------------------------------------
// THE BODIES ARE HERE, NOT IN A SOURCE FILE. A driver that calls out of its own translation unit
// reloads the arena's and the run array's bases around every call and spills its pivot loop's
// registers. Every ordering driver is in its own unit with its graph, so that they can be compared
// with EACH OTHER; it is faster where it ships for the same reason. Every out-of-class definition
// below is `inline`, so the units that include this header share one copy. See
// docs/CODING_RULES.md for the rule and the mechanics it needs.
// ------------------------------------------------------------------------------------------------


// PEAK LIVE CLIQUE MEMBERS OF THE LAST ORDERING TO RUN, written by every driver that tracks it and
// read by tests/test_order.cpp. One symbol rather than one per driver, because the whole use is to
// compare two drivers back to back: run one, read this, run the other, read it again.
//
// A GLOBAL RATHER THAN A RETURN VALUE, and deliberately. The figure is a cross-check between
// implementations rather than a result anyone orders a matrix to obtain, so it does not belong in
// the public ordering signature; `gAmd3BCompactions` is here for the same reason. Not thread safe,
// and it does not need to be: nothing writes it outside a test.
inline std::size_t gPeakCliqueMembers = 0;

inline QuotientGraph::QuotientGraph(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mRun(colPtr.empty() ? 0 : colPtr.size() - 1),
      mMark(mRun.size(), NIL) {
    const std::int32_t size = static_cast<std::int32_t>(mRun.size());

    // One pass to place each column's run, dropping the diagonal. The runs are laid out in column
    // order and never move afterwards, so the offsets are written once here and only the lengths
    // change. The run is u's whole allowance: A[u] fills it now and I[u] grows into the room A[u]
    // gives up, never past mRun[u + 1].sourcePtr.
    mSource.reserve(colPtr.empty() ? 0 : colPtr.back());
    for (std::int32_t aj = 0; aj < size; ++aj) {
        // The run's start is the arena's length before this column is appended, which is what the
        // separate `mSourcePtr` array used to hold at index aj. Its one extra entry, the n-th, was
        // read only by the next iteration of this loop and is now the cursor itself.
        mRun[aj].sourcePtr = mSource.size();
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSource.push_back(rowIdx[cp]);
        // THE ONE CROSSING. A difference of two positions is a count, so this is the single place
        // in the class where a two-dimensional quantity is written into a one-dimensional one. It
        // is bounded by deg(aj) and so by n, which the SparseMatrix constructor has already capped
        // at MAX_IDX, but the cast is written rather than left implicit because that bound is an
        // argument and not something the types say.
        mRun[aj].adjacencySize = static_cast<std::uint32_t>(mSource.size() - mRun[aj].sourcePtr);
    }

    // The clique arena is reserved for the same size as the source pool, and that one line was
    // worth 11 percent of AMD3's run time.
    //
    // It is APPEND ONLY: a new clique is written at the end and a dead one's block is left where
    // it lies, so the arena grows to the sum of |C[p]| over the whole elimination, which is a
    // little larger than nnz(A) on a 2D grid, 108705 against 97440 at 140 a side. Unreserved, a
    // vector reaches that by doubling from nothing: 18 reallocations, 131071 entries copied per
    // ordering, and the last few blocks are large enough that the allocator serves them from mmap
    // and every one faults its pages in on first touch. Instruments put 926 ms of a 8.54 s run
    // inside `__builtin_operator_new` beneath reachableSet, which is where the insert below is
    // inlined; the whole of it was growth.
    //
    // nnz(A) is a starting size rather than a bound. It leaves at most one doubling on a 2D grid
    // and it is the same scale as the input, which is the honest thing to say about it; a problem
    // whose fill is many times its input will still grow, now from a large base and so amortized.
    // Amd.cpp does not have this problem at all: its Iw is one pool sized once and COMPACTED IN
    // PLACE when it fills, which is what AMD_NCMPA counts and which is 1 for a whole 140x140 run.
    // Reclaiming the dead blocks the same way is the real fix and is not this.
    mCliqueArena.reserve(colPtr.empty() ? 0 : colPtr.back());

    // Every vertex begins as a supervariable of one: a chain holding itself, so next is NIL and
    // last is the vertex. Mass elimination splices these together and order() walks them.
    mSuperNext.assign(size, NIL);
    mSuperLast.resize(size);
    mWeight.assign(size, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

// APPENDS to `reached` and does not clear it, which is what lets beginElimination point it straight
// at the clique arena instead of at a scratch. The returning overload below clears for itself.
inline void QuotientGraph::reachableSet(std::int32_t u, std::vector<std::int32_t>& reached) {
    // reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}
    //
    // The mark array is the set: mMark[v] == mTag is the membership test, one comparison, and
    // mMark[v] = mTag is the insertion, one store. So the union costs one pass per source rather
    // than a hash per member, and the result comes out in the order the sources were walked.
    //
    // The buffer is the caller's so that a caller in a loop can keep one, which is the whole of
    // what the returning overload costs.
    // Eliminated vertices are skipped rather than purged. Mass elimination cannot leave one
    // behind, since a vertex it merges belongs to the pivot's clique alone and is removed from
    // it; a live merge can, since the vertex it folds away is left where it lies at weight zero
    // and every clique that named it still does.
    // Liveness and membership come off ONE load, `mMark[v] < mTag`, which is genmmd's
    // `marker[nd] < tag`. `Amd.cpp` reaches the same place from the other side, testing
    // `nvi = Nv[i]` before using an entry, so its liveness and its value are one load of an array
    // it needs anyway. See GONE in the header for the encoding and for why the WEIGHT cannot
    // carry it, which is the substitution tried and reverted on 2026-08-08.
    //
    // The premise that attempt got wrong is what the encoding has to survive, so it is kept: the
    // prune removes the pivot from A[u] for every u in C[pivot], and that is NOT every list that
    // can still be walked. A vertex `number()` numbered in the prepass and a vertex `merge()`
    // folded away are both left exactly where they lie. GONE reaches them because it is WRITTEN
    // at every death site rather than inferred from a value.
    //
    // It matters because this is the hottest line in the ordering. Instruments put the clique-walk
    // copy of it at 235 ms of an 8.37 s AMD3 run at 140 a side, two random loads into two arrays
    // where Amd.cpp does one.
    // THE SIGN OF THE WEIGHT IS THE MEMBERSHIP MARK, 2026-08-17, which is `AMD_2`'s `Nv`. The
    // negation IS the insertion: `nv > 0` is "not yet emitted" and `mWeight[v] = -nv` is the emit,
    // where this was `mMark[v] < mTag` and `mMark[v] = mTag`. One load per clique member instead of
    // a mark load and a weight load in the caller's summing pass.
    //
    // C[pivot] IS LEFT NEGATED AND massEliminate PUTS IT BACK. That is the contract, and it is why
    // the negation can be afforded at all: the restore rides in a walk of the same set that already
    // exists. A caller that sets late mass elimination MUST call massEliminate; Amd3 does.
    //
    // THE ADJACENCY LOOPS STILL ASK mMark, the clique loops do not, and the asymmetry is exact
    // rather than cautious. `number()` numbers a prepass vertex, leaves its weight at one so its
    // neighbors' degrees still count it, and leaves it in the adjacency of every one of them, so a
    // positive weight does not mean live there. It cannot appear in a CLIQUE: the prepass completes
    // before the first elimination, so no clique existed when it was numbered, and every clique
    // since is built from a reach that skipped it. genmmd reaches the same arrangement, using
    // `qsize != 0` inside clique walks only and `marker` everywhere else.
    ++mTag;
    mWeight[u] = -mWeight[u];              // never its own neighbor
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique. Measured at 300 ms of AMD1's
    // 6.31 s on alpamayo, in the one accessor that reads as though it were free.
    const VertexRun&    run           = mRun[u];          // one fetch; see the member
    const std::int32_t* source        = mSource.data() + run.sourcePtr;
    const std::uint32_t adjacencySize = run.adjacencySize;
    const std::uint32_t incidenceSize = run.incidenceSize;
    const std::int32_t* incidence     = source + adjacencySize;
    // Hoisted like `live`, and for the same reason: member loads the compiler cannot prove are
    // unaliased by the stores below. Both branches are per SOURCE rather than per member, so
    // neither sits in the loop that does the work. See the members' notes for what each decides.
    const bool reverse  = mReverseIncidence;
    const bool amdOrder = mVendoredListOrder;

    // Which source is walked first. genmmd expands the variables and then the cliques, which is
    // how the whole md ladder is laid out; AMD_2 takes the cliques first and the supervariables
    // only on its last pass. Same set either way, and the order decides C[pivot]'s content order,
    // hence which of two equal-degree candidates a later iteration finds first.
    if (amdOrder) {
        for (std::uint32_t ii = 0; ii < incidenceSize; ++ii) {
            const std::uint32_t i           = reverse ? incidenceSize - 1 - ii : ii;
            const std::int32_t  c           = incidence[i];
            const std::int32_t* members     = mCliqueArena.data() + mRun[c].sourcePtr;
            const std::uint32_t membersSize = mRun[c].adjacencySize;
            for (std::uint32_t k = 0; k < membersSize; ++k) {
                const std::int32_t v  = members[k];
                const std::int32_t nv = mWeight[v];
                if (nv > 0) { mWeight[v] = -nv; reached.push_back(v); }
            }
        }
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v  = source[k];
            const std::int32_t nv = mWeight[v];
            if (nv > 0 && !(mHasNumbered && mMark[v] == GONE)) { mWeight[v] = -nv; reached.push_back(v); }
        }
        return;
    }

    for (std::uint32_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v  = source[k];
        const std::int32_t nv = mWeight[v];
        if (nv > 0 && !(mHasNumbered && mMark[v] == GONE)) { mWeight[v] = -nv; reached.push_back(v); }
    }
    for (std::uint32_t ii = 0; ii < incidenceSize; ++ii) {
        const std::uint32_t i           = reverse ? incidenceSize - 1 - ii : ii;
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mCliqueArena.data() + mRun[c].sourcePtr;
        const std::uint32_t membersSize = mRun[c].adjacencySize;
        for (std::uint32_t k = 0; k < membersSize; ++k) {
            const std::int32_t v  = members[k];
            const std::int32_t nv = mWeight[v];
            if (nv > 0) { mWeight[v] = -nv; reached.push_back(v); }
        }
    }
}

inline std::vector<std::int32_t> QuotientGraph::reachableSet(std::int32_t u) {
    std::vector<std::int32_t> reached;   // empty, so the appending overload needs no clear
    reachableSet(u, reached);
    // AND IT UNDOES THE NEGATION, which the appending overload deliberately leaves for
    // massEliminate. This form is a convenience with no caller inside the elimination, so a reader
    // reaching for it should get a query rather than a half-finished elimination.
    mWeight[u] = -mWeight[u];
    for (std::int32_t v : reached) mWeight[v] = -mWeight[v];
    return reached;
}

inline std::uint32_t QuotientGraph::reachableSize(std::int32_t u) {
    // The same two passes as reachableSet, counting rather than collecting.
    ++mTag;
    std::uint32_t reached = 0;   // DISTINCT vertices, the mark seeing to that, so at most n
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique. Measured at 300 ms of AMD1's
    // 6.31 s on alpamayo, in the one accessor that reads as though it were free.
    const VertexRun&    run           = mRun[u];          // one fetch; see the member
    const std::int32_t* source        = mSource.data() + run.sourcePtr;
    const std::uint32_t adjacencySize = run.adjacencySize;
    const std::uint32_t incidenceSize = run.incidenceSize;
    for (std::uint32_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v = source[k];
        if (mMark[v] != GONE) { mMark[v] = mTag; ++reached; }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::uint32_t i = 0; i < incidenceSize; ++i) {
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mCliqueArena.data() + mRun[c].sourcePtr;
        const std::uint32_t membersSize = mRun[c].adjacencySize;
        for (std::uint32_t k = 0; k < membersSize; ++k) {
            const std::int32_t v = members[k];
            if (mMark[v] < mTag) { mMark[v] = mTag; ++reached; }
        }
    }
    return reached;
}

inline std::uint32_t QuotientGraph::reachableWeight(std::int32_t u) {
    ++mTag;
    std::uint32_t reached = 0;   // a sum over DISTINCT vertices, so bounded by n; see the header
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique. Measured at 300 ms of AMD1's
    // 6.31 s on alpamayo, in the one accessor that reads as though it were free.
    const VertexRun&    run           = mRun[u];          // one fetch; see the member
    const std::int32_t* source        = mSource.data() + run.sourcePtr;
    const std::uint32_t adjacencySize = run.adjacencySize;
    const std::uint32_t incidenceSize = run.incidenceSize;
    for (std::uint32_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v = source[k];
        if (mMark[v] != GONE) { mMark[v] = mTag; reached += static_cast<std::uint32_t>(mWeight[v]); }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::uint32_t i = 0; i < incidenceSize; ++i) {
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mCliqueArena.data() + mRun[c].sourcePtr;
        const std::uint32_t membersSize = mRun[c].adjacencySize;
        for (std::uint32_t k = 0; k < membersSize; ++k) {
            const std::int32_t v = members[k];
            if (mMark[v] < mTag) { mMark[v] = mTag; reached += static_cast<std::uint32_t>(mWeight[v]); }
        }
    }
    return reached;
}

inline void QuotientGraph::number(std::int32_t u) {
    // A numbered vertex lingers in every list that named it, deliberately: its neighbors keep
    // degrees that still count it. GONE is what stops the walks following it back in.
    //
    // AND THE FLAG IS WHAT TELLS THE WALKS TO ASK. See mHasNumbered: this is the only thing that
    // sets it, so a run that never calls this function never pays for the test.
    mHasNumbered = true;
    mMark[u]     = GONE;
}

inline void QuotientGraph::setAside(std::int32_t u) {
    // ZERO WEIGHT IS THE WHOLE MECHANISM. `reachableSet` takes a vertex on `nv > 0` and the prune
    // keeps one on the same test, so a zero-weight vertex is unreachable and is dropped from every
    // list the first time that list is rewritten. GONE additionally stops `eliminated` reporting
    // it live, which is what this class asks rather than the weight.
    //
    // ITS NEIGHBORS KEEP DEGREES THAT STILL COUNT IT, exactly as after `number`, and `AMD_2` does
    // not correct them either: a degree is a bound and one that is too large only delays a pivot.
    mWeight[u] = 0;
    mMark[u]   = GONE;
}

inline void QuotientGraph::beginElimination(std::int32_t pivot) {
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
    // the arena, `mCliqueArena.data() + mRun[c].sourcePtr`, while appending the reach to that same
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
    // The remedy is to make the arena unable to move rather than to re-fetch per clique, which
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
    if (mCliqueArena.capacity() - mCliqueArena.size() < size())
        mCliqueArena.reserve(std::max(2 * mCliqueArena.capacity(), mCliqueArena.size() + size()));

    // THE DESCRIPTOR IS HELD IN LOCALS UNTIL BOTH READERS OF THE PIVOT'S RUN ARE PAST, which is
    // the whole subtlety of storing a clique in the run of the vertex that formed it. The slots it
    // will occupy are still being read: `reachableSet` walks A[pivot] and I[pivot] through them,
    // and the absorbed-clique loop below finds I[pivot] the same way. Writing either early leaves
    // a walk reading the arena through an offset into mSource. See the header.
    const std::size_t cliqueStart = mCliqueArena.size();
    reachableSet(pivot, mCliqueArena);          // appends; see its note
    // THE SECOND AND LAST CROSSING. The arena's new length less this block's start is a member
    // count, so a two-dimensional quantity is written into a one-dimensional one, exactly as the
    // constructor does for the source runs. Bounded by n, a reach having at most n entries, which
    // is the same bound the reserve above relies on.
    const std::uint32_t cliqueLen = static_cast<std::uint32_t>(mCliqueArena.size() - cliqueStart);

    // Taken AFTER the append, since that is what can move the arena. With the reserve above it
    // cannot have moved, and this stays as it is regardless: it costs nothing and it is the shape
    // that remains correct if the reserve is ever revised.
    const std::int32_t* reached     = mCliqueArena.data() + cliqueStart;
    const std::uint32_t reachedSize = cliqueLen;

    // The absorbed cliques are I[pivot], read where they lie. Nothing below writes the pivot's run
    // (the prune rewrites the runs of C[pivot]'s members, and the pivot is not one of them), so no
    // copy and no scratch is needed to keep them alive across the passes that follow.
    const std::int32_t* absorbedCliques = mSource.data() + mRun[pivot].sourcePtr + mRun[pivot].adjacencySize;
    const std::uint32_t absorbedSize    = mRun[pivot].incidenceSize;
    for (std::uint32_t i = 0; i < absorbedSize; ++i)
        killClique(absorbedCliques[i]);              // dead, its block left behind

    // Both readers of the pivot's own run are past, so the run becomes the clique's descriptor.
    mRun[pivot].sourcePtr     = cliqueStart;
    mRun[pivot].adjacencySize = cliqueLen;

    // A CLIQUE IS BORN HERE AND NOWHERE ELSE, which is what makes the peak countable. See
    // `numPeakCliqueMembers`.
    mNumLiveCliqueMembers += cliqueLen;
    mNumPeakCliqueMembers  = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);
#ifndef NDEBUG
    mCliqueOwners.push_back(pivot);
#endif

    // Stamp the new clique. THE ABSORBED CLIQUES ARE NOT STAMPED: the loop above has just set
    // `mRun[c].adjacencySize = 0` for every one of them, and a dead clique is exactly a clique of
    // size zero, so the prune's incidence compaction asks `cliqueSize(c) != 0` and needs no tag
    // and no second pass.
    //
    // That matters beyond the pass it deletes. Clique ids and vertex ids share this array, so
    // stamping a clique wrote a live tag over the slot of the VERTEX that formed it, which is
    // dead. Harmless while the mark carried only "seen this step"; fatal once it also carries
    // "dead", since a dead pivot can still sit as a member of an older clique that is still
    // alive, and a walk of that clique would read the borrowed tag and take it for a live vertex.
    // Neither genmmd nor AMD_2 shares one stamp array between the two kinds.
    //
    // THE STAMPING PASS STAYS, and this is a measured decision rather than an oversight. It looks
    // redundant: reachableSet has just written `mMark[v] = mTag` on every member as it emitted
    // it, so `inClique` could BE that tag and this walk could go, which is what genmmd does and
    // what Mmd3B does. Built on 2026-08-15, with the weighted size accumulated in reachableSet's
    // four emit sites so that nothing was lost, and it measured WORSE: 74000 more instructions
    // and 142000 more data reads on a 100x100 grid.
    //
    // The reason is worth keeping, because it is not obvious and it applies to any fold of this
    // shape. This walk runs over the members ACTUALLY EMITTED; the emit sites run over every
    // candidate EXAMINED, which is more, since a vertex reached through two sources is examined
    // twice and emitted once, and dead ones are examined and never emitted. So moving a per
    // member cost into the emit moves it onto the busier loop. Reverted.
    // NO STAMPING PASS AND NO TAG, 2026-08-17. Membership was written by the walk, in the sign of
    // the weight, so this only sums. The `inClique` out-parameter went with the stamp: both prunes
    // read the sign now and neither looked at it.
    // The weighted size of the new clique is accumulated HERE rather than in a pass of its own in
    // each driver: this loop has the member loaded already. See cliqueWeight().
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t k = 0; k < reachedSize; ++k)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[k]]);
    mCliqueWeight = cliqueWeight;
}

inline const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot) {
    beginElimination(pivot);
    // C[pivot] IS the reach here: beginElimination wrote it there and nothing has trimmed it yet
    // (massEliminate does, and runs after). So the prune walks the arena block rather than a copy.
    const std::int32_t* reached     = mCliqueArena.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;
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
    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        // ONE FETCH OF THE RUN, not three. The three numbers share a 16-byte object, so the
        // reference below brings all of them in on one line; read as `mRun[u].sourcePtr` and
        // friends they would still be three subscripts of the same object and the compiler would
        // still load once, but naming it once is what makes that visible to a reader.
        const VertexRun&  run           = mRun[u];
        std::int32_t*     source        = mSource.data() + run.sourcePtr;
        // The two counters are one-dimensional COUNTS, positions in a list bounded by deg(u) and
        // so by n, where std::size_t is for a position into an n x n object. They take the type of
        // what they count, so they move with the array: `std::uint32_t`, and the cast that used to
        // stand here is gone. It was there only because the array was wider than the loop. The
        // dimensional rule in experiments/ordering/REPORT.md asks for this everywhere and it was
        // taken here first because this is the hottest loop in the ordering, Instruments put
        // 277 ms of an 8.53 s run on the incidence loop's header.
        //
        // `heldVertex` is a VERTEX and stays std::int32_t, carrying NIL.
        const std::uint32_t adjacencySize = run.adjacencySize;
        std::uint32_t       kept          = 0;
        std::int32_t        heldVertex    = NIL;        // the first survivor, appended last
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            // ONE LOAD, THREE QUESTIONS, which is Amd.cpp's `nvj = Nv [j] ; if (nvj > 0)`. A
            // negative weight is a member of the new clique, the pivot included, so the explicit
            // `v == pivot` test goes with the membership test; a zero is a vertex a live merge
            // folded away. The FOURTH question, whether v was numbered by a prepass, still needs
            // mMark, and dropping it is not an option: massEliminate reads `adjacencySize == 0`,
            // so a numbered leftover in A[u] would suppress a merge that used to fire.
            if (mWeight[v] <= 0) continue;             // in the new clique, the pivot, or merged
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by a prepass; see above
            if (amdOrder && heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mRun[u].adjacencySize = kept;                      // A[u] - C[pivot] - {pivot}

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
        const std::uint32_t incidenceSize = run.incidenceSize;
        std::uint32_t       write         = kept;      // follows `kept`; see the note above
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
            if (mRun[c].adjacencySize != 0) source[write++] = c;   // dead is size zero; see above
        }
        source[write++]   = pivot;                     // u joins the new clique, id = pivot
        mRun[u].incidenceSize = write - kept;
        // [c1, ..., ck, pivot] to [pivot, c2, ..., ck, c1], which is a SWAP of the two boundary
        // entries: only they move and c2..ck stay put. The pivot cannot be written first, which
        // is what the loop above would otherwise allow: the write cursor starts at `kept` and the
        // read at the original adjacencySize, and those are equal whenever nothing was pruned from
        // A[u], so an extra write before the reads finish clobbers an unread entry. AMD_2 makes
        // its three assignments after both compactions for exactly this reason.
        if (amdOrder && write - kept > 1) std::swap(source[kept], source[write - 1]);
    }

    return finishElimination(pivot);
}

inline const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot,
                                                                 TaggedScan& scan) {
    // A CLIQUE DIES TWO WAYS AND THE TAGGED W MUST LEARN ABOUT BOTH. Aggressive absorption zeroes
    // `w[c]` in the driver; elimination-time absorption, which is I[pivot], happens inside
    // `beginElimination` and used to be recorded only as a zero size. So the prune could test
    // neither alone and tested the size, an array it reads for nothing else, once per incidence
    // clique on every reached vertex. `AMD_2` writes both deaths into W, `Pe[e] = FLIP(me)` with
    // `W[e] = 0`, and its scan tests `we != 0` off the load it already needs for the value.
    //
    // Done HERE rather than inside beginElimination, which is the only reason it is a separate
    // walk. beginElimination is shared with the mmd drivers and they have no W; giving it a `w`
    // parameter would push an amd concept through mmd's path for nothing. I[pivot] is still intact
    // at this point, so this walks the same short list beginElimination is about to read.
    {
        const std::int32_t* absorbed = mSource.data() + mRun[pivot].sourcePtr
                                                      + mRun[pivot].adjacencySize;
        const std::uint32_t count    = mRun[pivot].incidenceSize;
        for (std::uint32_t i = 0; i < count; ++i) scan.w[absorbed[i]] = 0;
    }

    beginElimination(pivot);
    const std::int32_t* reached     = mCliqueArena.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;

    // EVERY MEMBER OF C[pivot] LEAVES THE DEGREE LISTS HERE, which is AMD_2's "remove variable i
    // from degree list" inside CONSTRUCT NEW CLIQUE. Not bookkeeping moved earlier for its own
    // sake: it is what frees each member's mPrev and mNext, so the hash key and the hash chain can
    // live there for the rest of the step and no hashNext array is needed. The driver files them
    // again in its bound pass, where the new degree is known, so no vertex is out of the lists
    // across a pivot selection.
    //
    // THIS IS WHY THE FOLD IS ON THE TAGGED PATH ONLY. mmd leaves its candidates filed and asks
    // `filed(u)` and `outmatched(u)` about exactly these vertices in its refresh; unfiling them
    // here would change what those tests answer.
    // Only where the driver asked for it; see TaggedScan. Hoisted out of the loop, one test per
    // pivot rather than per member.
    if (scan.buckets != nullptr)
        for (std::uint32_t ri = 0; ri < reachedSize; ++ri) scan.buckets->unfile(reached[ri]);

    const bool          amdOrder    = mVendoredListOrder;
    const std::int32_t  wflg        = scan.wflg;          // hoisted, as the flags are

    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u       = reached[ri];
        // NEGATED, BECAUSE u IS A MEMBER OF C[pivot] AND SO READS NEGATIVE. This is Amd.cpp's
        // `nvi = -Nv [i]` under CONSTRUCT NEW CLIQUE, and the sign is the whole of it: `wnvi`
        // below is `wflg - nvi` and comes out wrong by twice the weight if the magnitude is not
        // taken. The cast that used to sit here is gone, `mWeight` being signed since 2026-08-17;
        // its own comment called it a signedness cast rather than a narrowing one, which was the
        // code saying the type was wrong. `wnvi` must still be able to go negative, which is
        // Amd.cpp's convention and the reason `w` is signed.
        const std::int32_t nvi     = -mWeight[u];
        const std::int32_t wnvi    = wflg - nvi;          // Amd.cpp's wnvi, and signed for it
        // ONE FETCH OF THE RUN, not three. The three numbers share a 16-byte object, so the
        // reference below brings all of them in on one line; read as `mRun[u].sourcePtr` and
        // friends they would still be three subscripts of the same object and the compiler would
        // still load once, but naming it once is what makes that visible to a reader.
        const VertexRun&  run           = mRun[u];
        std::int32_t*     source        = mSource.data() + run.sourcePtr;
        const std::uint32_t adjacencySize = run.adjacencySize;
        std::uint32_t       kept          = 0;
        std::uint32_t       explicitPart  = 0;            // a weight sum, not a count of positions
        // THE HASH KEY IS AMD_2'S EXACTLY. Four things differed and all four are here and in the
        // driver's bound pass. Amd.cpp: `hval = 0 ... hval += e ... hval += j ... hval % n`.
        //   NO `+ 1` ON A TERM; Amd.cpp adds the id itself.
        //   THE PIVOT IS NOT IN THE KEY. Its clique heads every I[u] this step and is shared by
        //     every member of C[p], so it cannot discriminate; Amd.cpp adds `me` to the list after
        //     the key is formed, which says the same thing by placement.
        //   THE MODULUS IS n, not n + 1.
        //   ONE REDUCTION AT THE END, not one per term. Amd.cpp accumulates in an unsigned Int and
        //     lets it WRAP at 2^32, which is what its own comment beside `hval % n` is about, so
        //     the accumulator is uint32 and the overflow is deliberate. Reducing per term gives a
        //     DIFFERENT key, not a cheaper one.
        std::uint32_t       key           = 0;            // wraps, like Amd.cpp's UInt hval
        std::int32_t        heldVertex    = NIL;          // see the plain prune above
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (mWeight[v] <= 0) continue;             // see the plain prune above
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by a prepass
            explicitPart += static_cast<std::uint32_t>(mWeight[v]);
            key += static_cast<std::uint32_t>(v);         // no + 1, no reduction; see above
            if (amdOrder && heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mRun[u].adjacencySize       = kept;
        // THE ADJACENCY HALF GOES INTO w[u], NOT INTO AN ARRAY OF ITS OWN. `w` is indexed by
        // CLIQUE id and a clique id is a dead pivot's id, so for a LIVE vertex the slot carries
        // nothing: u is in C[pivot] and therefore alive, and no clique is named after it until it
        // is eliminated, by which time this value is long consumed. The tagged W thus answers a
        // FOURTH question on top of Amd.cpp's three, and the array that used to carry this one is
        // gone. It cannot collide with the clique writes below: those are indexed by c drawn from
        // I[u], every one of which is a dead pivot, and u is live.
        //
        // The driver's obligation is one store: the slot goes back to alive-and-unseen once the
        // bound has been read. See src/Amd3.cpp.
        scan.w[u]               = static_cast<std::int32_t>(explicitPart);
        // Through the int32 slot the key rides in, bit pattern preserved and read back as uint32
        // in the driver's bound pass. The slot is the vertex's degree-list predecessor, free
        // because every member of C[pivot] was unfiled above. See Buckets.
        // The key is accumulated either way, one add per surviving neighbour, and STORED only
        // where the driver asked for the bucket arrangement. Accumulating unconditionally keeps
        // the inner loop branch-free; the cost to a driver that computes its own key is that one
        // add, against a test per clique if it were guarded.
        if (scan.buckets != nullptr)
            scan.buckets->setKey(u, static_cast<std::int32_t>(key));   // the ADJACENCY half

        const std::int32_t* incidence     = source + adjacencySize;
        const std::uint32_t incidenceSize = run.incidenceSize;
        std::uint32_t       write         = kept;         // follows `kept`; see the plain prune
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
            // Amd.cpp's four lines, transcribed. A clique seen earlier in this step already holds
            // the running value above the tag; one seen for the first time starts from |C[c]| and
            // is listed once; one already absorbed reads ZERO and is dropped from the list here.
            //
            // ONE LOAD, TWO QUESTIONS. `we == 0` is dead, which this used to ask of `mCliqueSize`,
            // an array it read for nothing else, and the value is wanted anyway two lines down.
            // The other half of the change is the w-zeroing at the top of this function.
            std::int32_t we = scan.w[c];
            if (we == 0) continue;                        // absorbed and gone; Amd.cpp's W == 0
            source[write++] = c;
            if (c == pivot) continue;                     // the new clique subtracts from nothing
            if (we >= wflg) {
                we -= nvi;
            } else {
                // A SIGNEDNESS cast, as `nvi` above: `cliqueDegree` is already 32 bits and `wnvi`
                // is negative in the early eliminations.
                we = static_cast<std::int32_t>(scan.degree[c]) + wnvi;
                scan.touchedCliques.push_back(c);
            }
            scan.w[c] = we;
        }
        source[write++]   = pivot;
        mRun[u].incidenceSize = write - kept;
        if (amdOrder && write - kept > 1) std::swap(source[kept], source[write - 1]);
    }

    return finishElimination(pivot);
}

inline const std::vector<std::int32_t>& QuotientGraph::finishElimination(std::int32_t pivot) {
    // Under mLateMassElimination the merge is the caller's, run after it has absorbed, so this
    // hands back an empty list and C[pivot] stays reach(pivot) exactly. See the setter.
    if (mLateMassElimination) {
        mMerged.clear();
    } else {
        massEliminate(pivot);
    }

    // ONLY THE INCIDENCE HALF IS CLEARED. `adjacencySize` is no longer the pivot's A[pivot]: it
    // now holds |C[pivot]|, the clique this elimination just built, and zeroing it here would
    // destroy it. This line cleared both while both were dead; one of them has a second job now.
    mRun[pivot].incidenceSize = 0;
    mMark[pivot]          = GONE;
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
inline const std::vector<std::int32_t>& QuotientGraph::massEliminate(std::int32_t pivot) {
    std::vector<std::int32_t>& merged = mMerged;   // scratch, kept for its capacity
    merged.clear();
    // Walks C[pivot], which is still the full reach: the trim below is this function's own and
    // happens after the loop.
    const std::int32_t* reached     = mCliqueArena.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;
    // THE SIGNS COME BACK HERE, IN A PASS THAT ALREADY EXISTS, and that is what makes the encoding
    // in reachableSet affordable. Amd puts its restore in a bound pass; there is no bound pass on
    // the mmd side, and this walk over C[pivot] is the only other traversal of that same set. The
    // pivot goes first, since the merge below adds into it and both operands must be magnitudes.
    //
    // SO EVERY PATH THROUGH AN ELIMINATION MUST REACH THIS FUNCTION. finishElimination calls it
    // unless mLateMassElimination is set, and the one driver that sets it, Amd3, calls it itself
    // after aggressive absorption. Nothing between the two reads a weight directly, and `weight()`
    // returns the magnitude in any case.
    mWeight[pivot] = -mWeight[pivot];
    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        mWeight[u] = -mWeight[u];                          // live again, and positive
        // Under mVendoredListOrder the new clique sits at the FRONT of I[u] rather than the back,
        // so the single remaining entry is at the head of the incidence run either way: with A[u]
        // empty the run starts with I[u], and with one clique there is only one position. The
        // test therefore needs no branch on the flag.
        if (mRun[u].adjacencySize == 0 && mRun[u].incidenceSize == 1 &&
            mSource[mRun[u].sourcePtr] == pivot) {         // A[u] empty, so I[u] starts at the run
            mRun[u].incidenceSize = 0;
            mMark[u]          = GONE;
            merged.push_back(u);
        }
    }
    if (!merged.empty()) {
        // THE MERGE HAPPENS FIRST, so the compaction can read the ZERO WEIGHT it leaves rather
        // than a stamp of its own. The old order was the reverse and needed a tag pass over
        // `merged` plus a mark read per member; the weight says the same thing and the
        // supervariable bookkeeping had to write it anyway.
        //
        // NO OTHER MEMBER OF C[pivot] CAN READ ZERO, which is what makes the test exact: a vertex a
        // live merge folded away is left at weight zero but is also stamped GONE, so no reach ever
        // emits it into a clique again.
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
            mCliqueWeight -= static_cast<std::uint32_t>(mWeight[u]);
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }

        std::int32_t*     members     = mCliqueArena.data() + mRun[pivot].sourcePtr;
        const std::uint32_t membersSize = mRun[pivot].adjacencySize;
        std::uint32_t       kept        = 0;
        for (std::uint32_t k = 0; k < membersSize; ++k)
            if (mWeight[members[k]] != 0) members[kept++] = members[k];
        trimClique(pivot, kept);         // a shrink is a partial death; see numPeakCliqueMembers
    }
    return merged;
}

// THE SAME COMPACTION, ONE PASS LATER AND FOR THE OTHER REASON. `massEliminate` above drops the
// members it merged into the pivot; this drops the members supervariable DETECTION absorbed into
// each other, which happens after it. `AMD_2` needs one place for both because its detection sits
// inside the scan, so its RESTORE DEGREE LISTS pass sees every casualty at once and writes the
// survivors back with `Iw [p++] = i`, then `Len [me] = p - pme1`. Ours needs two because detection
// is a pass of its own.
//
// WITHOUT THIS THE ABSORBED STAY IN THE CLIQUE FOR THE REST OF THE RUN. No permutation moves,
// every later walk skipping them on the `nv > 0` test, but they are visited, once per walk of that
// clique for as long as it lives, which on a grid is many.
//
// NO CURSOR TO PULL BACK, unlike `AMD_2`'s `if (elenme != 0) pfree = p` and unlike `Amd3B`'s. The
// arena's length is the vector's own and the clique need not be its last block, so the trimmed
// tail is left as a hole like every other dead entry here. This is the arena's whole bargain:
// nothing is reclaimed, so nothing has to be reclaimable.
inline void QuotientGraph::trimClique(std::int32_t pivot, std::uint32_t kept) {
    mNumLiveCliqueMembers -= mRun[pivot].adjacencySize - kept;   // the trimmed tail is not live
    mRun[pivot].adjacencySize = kept;
}

// A DEAD CLIQUE IS EXACTLY A CLIQUE OF SIZE ZERO, which every reader here already relies on: the
// prune's incidence compaction asks `cliqueSize(c) != 0` and needs no tag and no second pass. So
// death is one store, and the only reason it is a function is that the peak has to see it.
//
// ONLY FOR A VERTEX THAT ACTUALLY FORMED ONE. A clique is born in `beginElimination` and nowhere
// else, so `adjacencySize` means a clique's length only for a vertex that has been a pivot; for
// any other it is still A[v]'s length. `merge` therefore does NOT call this, and says so.
inline bool QuotientGraph::cliqueCountBalances() const {
#ifdef NDEBUG
    return true;
#else
    std::size_t live = 0;
    for (std::int32_t c : mCliqueOwners) live += mRun[c].adjacencySize;
    return live == mNumLiveCliqueMembers;
#endif
}

inline void QuotientGraph::killClique(std::int32_t c) {
    mNumLiveCliqueMembers -= mRun[c].adjacencySize;
    mRun[c].adjacencySize = 0;
}

inline void QuotientGraph::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    // NOT `killClique`. v is a live supervariable being absorbed, and it never formed a clique, so
    // this length is A[v]'s and not a clique's; feeding it to the counter would corrupt the peak.
    mRun[v].adjacencySize = 0;
    mRun[v].incidenceSize = 0;
    mMark[v]          = GONE;
}

inline void QuotientGraph::absorb(const std::vector<std::int32_t>& cliques,
                           const std::int32_t* vertices, std::uint32_t vertexCount) {
    if (cliques.empty()) return;

    for (std::int32_t c : cliques) killClique(c);   // dead; the prune reads the size

    for (std::uint32_t k = 0; k < vertexCount; ++k) {  // I[u] - dead, compacted in place
        const std::int32_t u         = vertices[k];
        std::int32_t*      incidence = mSource.data() + mRun[u].sourcePtr + mRun[u].adjacencySize;
        const std::uint32_t size     = mRun[u].incidenceSize;

        // THE ROTATION HAS TO BE REDONE WHEN ITS SUBJECT DIES, 2026-08-18, and this is the whole
        // of the correction. Both codes end I[u] as `[pivot][entries 1..k][entry 0]`, the first
        // entry moved to the back. `AMD_2` performs that rotation inside scan 2, where an
        // aggressively absorbed clique has ALREADY been dropped, so `entry 0` is the first
        // SURVIVOR of this step's absorption. Our prune rotates first and absorbs afterwards, so
        // when the entry it parked at the back is one this pass then removes, the two lists come
        // out in different orders.
        //
        // It costs nothing when it does not apply: the test is one read of the last slot, and the
        // extra rotation runs only for a vertex whose parked entry actually died. Where the parked
        // entry survives, it is still the first survivor and the two agree already.
        //
        // WHAT IT COST BEFORE, and why a six-entry difference was worth this. On GHS_indef/aug2d,
        // n = 29008, two hash buckets of 1314 came out with two entries transposed, which changed
        // which vertex of a supervariable absorbed the others: six positions of the permutation,
        // identical fill, and one of twelve matrices in benchmarks/matrices that differed from
        // `AMD_2` for this reason. Turning aggressive absorption off in both made the orders
        // identical, which is what identified it.
        const bool parkedDied = size > 0 && mRun[incidence[size - 1]].adjacencySize == 0;

        std::uint32_t kept = 0;
        for (std::uint32_t i = 0; i < size; ++i)
            if (mRun[incidence[i]].adjacencySize != 0) incidence[kept++] = incidence[i];
        mRun[u].incidenceSize = kept;

        // Entry 0 is the pivot's own new clique, which the prune put at the front and which is
        // never absorbed, so the rotation runs over positions 1 onward.
        if (parkedDied && kept > 2)
            std::rotate(incidence + 1, incidence + 2, incidence + kept);
    }
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
inline std::vector<std::int32_t> QuotientGraph::orderAscending(
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
        pos += static_cast<std::uint32_t>(mWeight[pivot]);  // the whole supervariable's room
    }

    for (std::size_t v = 0; v < n; ++v) {                   // ascending, so the members are too
        const std::int32_t s = slot[v];
        if (s >= 0) continue;                               // a root, already placed
        const std::int32_t root = -s - 1;
        order[static_cast<std::size_t>(slot[root]++)] = static_cast<std::int32_t>(v);
    }
    return order;
}

inline std::vector<std::int32_t>
QuotientGraph::order(const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order;
    order.reserve(size());
    for (std::int32_t pivot : pivots)
        for (std::int32_t u = pivot; u != NIL; u = mSuperNext[u]) order.push_back(u);
    return order;
}


} // namespace Oblio
