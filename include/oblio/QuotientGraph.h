#pragma once

// QuotientGraph.h - the representation Oblio's own minimum-degree orderings run on, and the degree
// buckets they pick from. Shared by every driver on this layout, mmd and amd alike.
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
// have. Section 5.3 of archive/sparse_factorization.md carries the derivation.
//
// A[u] and I[u] are the two kinds of SOURCE the union is taken over, and their number falls
// monotonically, so the two lists share one block sized once from u's column of A and never grown.
// C[c] has no such bound and gets a store of its own.
//
// A clique id is the pivot that created it, so cliques index into the vertex space and need no id
// space of their own. A dead clique is an empty member list.
//
// No sets anywhere. Membership is a mark array stamped with a monotone tag, so a query is one
// comparison and nothing is allocated. Every list edit is a compaction in place.

#include "oblio/Types.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
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
    // ALL THREE ARE SIZE n. The heads are indexed by DEGREE and the links by VERTEX, two ranges of
    // the same extent: every branch files at the true degree, which reaches n - 1. A caller that
    // reads a bucket must bound its own index; the mmd prepass does.

    // FOUR FACTS IN `mPrev[u]`, told apart by sign, which is what lets three arrays be three
    // rather than four:
    //
    //     mPrev[u] in [0, n)       filed; the vertex before it is mPrev[u]
    //     mPrev[u] <= -2           filed AT THE HEAD of bucket -mPrev[u] - 2
    //     mPrev[u] == UNFILED      not in any list                            (-1)
    //     mPrev[u] == OUTMATCHED   withheld from the buckets; see outmatch()  (INT32_MAX)
    //
    // THE HEAD CASE IS WHAT DELETES A DEGREE ARRAY. `unfile` reads the bucket back out of the
    // slot and takes no degree argument, so no caller has to keep one to unfile with. The amd
    // drivers keep theirs for a different reason: their bound READS the value.
    //
    // A PREDECESSOR IS STORED RAW AND THE HEAD FORM IS SHIFTED BY TWO. The shift is two and not
    // one because DEGREE 0 IS REACHABLE, an isolated vertex filing at 0, and -1 is taken by
    // UNFILED.
    //
    // BOTH SENTINELS ARE EQUALITY-ONLY. Nothing here orders `mPrev` against either, so they sit at
    // the two values a predecessor and a head cannot take, leaving the whole range below -1 to the
    // head form.
    //
    // IT FITS EXACTLY AT THE LARGEST ADMISSIBLE n, with nothing to spare at either end: a head
    // encodes a degree reaching n - 1, so the head form runs to -(n + 1) and lands on the type's
    // minimum at n == MAX_IDX, while a predecessor reaches n - 1 and leaves INT32_MAX free above
    // it. Filing at anything above the true degree does not fit.
    //
    // A DEGREE IS ONE DIMENSIONAL and so `std::uint32_t` in every signature here. The links stay
    // `std::int32_t`, carrying NIL and UNFILED; only the key does not.
    static constexpr std::int32_t UNFILED    = -1;
    static constexpr std::int32_t OUTMATCHED = std::numeric_limits<std::int32_t>::max();

    explicit Buckets(std::size_t size)
        : mHead(size, NIL), mNext(size, NIL), mPrev(size, UNFILED) {}

    // buckets[degree].add(u), at the head. The head is the only O(1) end of a singly reachable
    // list, so the winner among equal degrees is whatever was filed last rather than the lowest
    // index. That is the vendored convention and it is why an ordering differs from an exact
    // scan's in its ties.
    void file(std::uint32_t degree, std::int32_t u) {
        mNext[u] = mHead[degree];
        mPrev[u] = -static_cast<std::int32_t>(degree) - 2;   // head of `degree`; see the encoding
        if (mHead[degree] != NIL) mPrev[mHead[degree]] = u;
        mHead[degree] = u;
    }

    // buckets[degree].discard(u). Idempotent, which matters during a batch: a vertex evicted
    // early can be merged away by a later pivot in the same round, and unfiling it twice must
    // not splice a list it is no longer in.
    void unfile(std::int32_t u) {
        const std::int32_t prev = mPrev[u];
        if (prev == UNFILED || prev == OUTMATCHED) return;   // not in a list; see the encoding
        if (prev >= 0) mNext[prev]      = mNext[u];
        else           mHead[-prev - 2] = mNext[u];          // u headed bucket -prev - 2
        if (mNext[u] != NIL) mPrev[mNext[u]] = prev;
        mNext[u] = NIL;
        mPrev[u] = UNFILED;
    }

    // Withhold u from the buckets without filing it anywhere, which is genmmd's `bwd[nd] = -maxint`
    // at the other end of the range. It stays live and reachable; it simply cannot be the minimum
    // before the vertex that outmatched it, so it is not a candidate until an elimination reaches
    // it and `restore` puts it back. mmdelm spells that restore `bwd[rn] = 0`, the same store that
    // unfiles, which is why the two sit together at every call site here.
    void outmatch(std::int32_t u)         { unfile(u); mPrev[u] = OUTMATCHED; }
    void restore(std::int32_t u)          { if (mPrev[u] == OUTMATCHED) mPrev[u] = UNFILED; }
    bool outmatched(std::int32_t u) const { return mPrev[u] == OUTMATCHED; }

    // The next vertex in the same bucket, and whether u is filed at all. Only mmd reads either.
    // A VERTEX OUT OF EVERY LIST HAS BOTH LINKS FREE, and the amd branch parks the hash key in one
    // and the hash chain in the other for the middle of an elimination step, which is why it needs
    // no key array and no hash-head links of its own.
    //
    // LEGAL ONLY BETWEEN unfile() AND file(). A key is an arbitrary int32 and can look like any of
    // the encodings mPrev carries, so `unfile()`, `filed()` and `outmatched()` MUST NOT be called
    // on a vertex holding one; they would splice a list on garbage. No mmd driver may call these:
    // mmd leaves its candidates filed and asks `filed()` about them.
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
    // A byte per vertex, not std::vector<bool>. That specialization packs one bit per entry, and
    // what it costs here is CONSTRUCTION: a vector of n false is built through the bit-reference
    // machinery word by word where a vector of n zero bytes is a memset, and this is built once
    // per ordering.
};


// The same, for a driver carrying a TAGGED array instead of a value array and a separate
// seen-this-step mark. One array holds three facts: `w[c] == 0` is absorbed, `0 < w[c] < wflg` is
// alive but stale, and `w[c] >= wflg` is seen this step with `w[c] - wflg` the value. So there is
// no mark to carry and no clearing pass.
//
// `key` carries the ADJACENCY HALF of the hash key alone, and the reason is a phase boundary
// rather than a preference: aggressive absorption runs between this prune and the driver's bound
// pass and COMPACTS I[u] in place, so the list the other half must sum over does not exist yet
// here.
//
// REDUCED AS IT ACCUMULATES, modulo the driver's bucket count, which is why it is an int32. The
// key is only ever used modulo that number, so reducing early cannot change which bucket a vertex
// lands in, and the running value stays narrow enough to ride in an array the driver already has.
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
// the eliminator, and refreshes whatever the elimination reached.
class QuotientGraph {
public:
    // Built straight from A's pattern. A is stored with both triangles, so a column's rows are
    // already its neighbors, and the only conversion is dropping the diagonal. Oblio's input
    // assumptions hold by construction, which is why nothing here symmetrizes, deduplicates or
    // sorts.
    //
    // ONE MARK PER VERTEX. Supervariable detection stamps into the driver's own tagged array, so
    // nothing asks for a second half.
    QuotientGraph(const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const { return mSize; }

    // EVERY MEMBER EVER PUT INTO A CLIQUE. Third of three member counters differing only in WHEN:
    // born is cumulative, live rises and falls, peak is live's running maximum.
    //
    // IT READS THE STORE'S LENGTH RATHER THAN A COUNTER, and the two agree only because a
    // contraction leaves the members it drops where they lie and nothing is ever reclaimed. A
    // reclaiming allocator here needs a counter of its own at the one site that bears a clique.
    std::size_t numBornCliqueMembers() const { return mCliqueSrc.size(); }

    // LIVE AND PEAK, the other two thirds of the counter above. Live is the current population and
    // peak its running maximum, taken at BIRTH ALONE since nothing else raises it.
    //
    // A CLIQUE HAS THREE EVENTS AND THE COUNTER SEES ALL OF THEM: born once in `beginElimination`;
    // CONTRACTED, keeping its identity and losing members; and DEAD, absorbed into a new clique.
    // mmd has one of each, amd one birth, two contractions and two deaths.
    //
    // IT IS A LIVENESS QUESTION AND NOT A PLACEMENT ONE, which is why it is exact under every
    // layout here: add on birth, subtract on both others, keep the maximum.
    std::size_t numPeakCliqueMembers() const { return mNumPeakCliqueMembers; }
    std::size_t numLiveCliqueMembers() const { return mNumLiveCliqueMembers; }

    // THE COUNTER CHECKED AGAINST A RECOMPUTATION, debug builds only. Births and deaths are spread
    // over four sites and nothing else in the suite would notice if they stopped balancing; the
    // fourth site, `massEliminate` shortening the pivot's own clique, was missed on the first
    // attempt and this is what found it.
    //
    // IT LIVES HERE RATHER THAN IN THE DRIVERS because only this class knows which vertices ever
    // formed a clique. Summing over a driver's pivot list looks equivalent and is not: `MmdFlat`
    // pushes prepass vertices onto that list, and for those `cliqueSize` still reports A[p]'s
    // length. That was the second wrong version of this check.
    bool cliqueCountBalances() const;
    // GONE: one value reserved above every tag, so the stamp array answers "is v dead" on the load
    // it was making anyway and no array is spent on liveness at any walk site.
    //
    // NOT the weight. `mWeight[v] != 0` is a PARTIAL flag: `number()` leaves a prepass vertex at
    // weight one deliberately, so its neighbors' degrees still count it. Used as the universal
    // test it lets a numbered vertex back into a reachable set.
    static constexpr std::int32_t GONE =
        std::numeric_limits<std::int32_t>::max();   // above every reachable tag

    // ZERO WEIGHT IS THE DEAD STATE on the amd branch, where the three ways a vertex leaves the
    // graph all end in one, so that branch allocates no marks at all.
    //
    // NOT EQUIVALENT TO THE MMD TEST, and the difference is two cases. `number` leaves a prepass
    // vertex live at weight one, which is mmd's alone; and a PIVOT is retired with its weight
    // intact, that weight being the supervariable's and `order` needing it. **So this predicate
    // must never be asked about a pivot.** The amd driver does not: a pivot is unfiled when chosen
    // and never revisited.
    // ALLOCATED ON DEMAND, and the amd branch never calls this. Call once, before any elimination.
    void enableMarks();

    // AND A TAG ON THE MMD BRANCH, because a zero weight is not available to it: `number` leaves
    // a prepass vertex live at weight one so that its neighbors' degrees still count it.
    bool eliminatedMmd(std::int32_t u) const { return mMark[u] == GONE; }

    bool eliminatedAmd(std::int32_t u) const { return mWeight[u] == 0; }

    // A driver may stamp into this same array rather than allocating one of its own, which is
    // what genmmd's `marker` is: `mmdelm` stamps it at level `tag` and `mmdupd` at level
    // `mt = tag + md0`, one array and one counter serving both. One counter is what makes it
    // safe, since two tags drawn from it can never be equal.

    std::int32_t advanceTag()                        { return ++mTag; }
    std::int32_t mark(std::int32_t u) const          { return mMark[u]; }
    void setMark(std::int32_t u, std::int32_t tag)   { mMark[u] = tag; }

    // The members of clique c, which after eliminating p is the pattern of p's column of L: the
    // vertices the pivot reached, less those it absorbed. A pointer and a length, as the
    // adjacency is, for the same reason.
    //
    // **Valid until the next elimination.** The members live in one arena that grows as cliques are
    // formed, so a reallocation moves them; the offsets are indices and survive it, but a pointer
    // taken beforehand does not. Read a clique at the moment of use, which is the same rule the
    // numeric factor's blocks live by.
    const std::int32_t* clique(std::int32_t c) const {
        return mCliqueSrc.data() + mSegment[c].srcPtr;
    }
    std::uint32_t cliqueSize(std::int32_t c) const { return mSegment[c].adjacencySize; }

    // WRITABLE, for the pass that trims the clique as it walks it. See
    // `trimClique` in the .cpp for why that pass and not a pass of its own. Same lifetime rule as
    // the const overload above.
    std::int32_t* clique(std::int32_t c) { return mCliqueSrc.data() + mSegment[c].srcPtr; }
    void trimClique(std::int32_t pivot, std::uint32_t kept);

    // THE ONE PLACE A CLIQUE DIES, and the reason it is one place is the counter above. Death has
    // three causes here: absorbed into the new clique, absorbed aggressively once its external
    // degree reaches zero, and merged away with its owner. All three go through here so the
    // counter sees them.
    void killClique(std::int32_t c);

    // |C[pivot]| WEIGHTED, accumulated in `beginElimination`'s stamping walk, which has the member
    // in hand already.
    //
    // VALID UNTIL THE NEXT ELIMINATION, and a scalar for that reason: read immediately after the
    // eliminator and never afterwards. Same contract as the pointer `clique()` returns.
    //
    // Over the UNTRIMMED clique. A driver that mass-eliminates late and needs the trimmed figure
    // recomputes it.
    std::uint32_t cliqueWeight() const { return mCliqueWeight; }

    // The two halves of the neighbor relation, read by an approximate degree, which decomposes
    // reach(u) rather than forming it. An exact degree has no use for either, uniting them through
    // reachableSet instead.
    //
    // A pointer and a length rather than a container, because that is what the storage is: the two
    // lists share one run, A[u] first and I[u] immediately behind it, which is why the incidence
    // lookup reads the adjacency's length.
    const std::int32_t* adjacencyMmd(std::int32_t u) const {
        return mAdjIncSrc.data() + mSegment[u].srcPtr;
    }
    const std::int32_t* adjacencyAmd(std::int32_t u) const {
        return mAdjIncSrc.data() + mSegment[u].srcPtr;
    }
    std::uint32_t adjacencySize(std::int32_t u) const { return mSegment[u].adjacencySize; }

    const std::int32_t* incidenceMmd(std::int32_t u) const {
        return mAdjIncSrc.data() + mSegment[u].srcPtr + mSegment[u].adjacencySize;
    }
    const std::int32_t* incidenceAmd(std::int32_t u) const {
        return mAdjIncSrc.data() + mSegment[u].srcPtr + mSegment[u].adjacencySize;
    }
    std::uint32_t incidenceSize(std::int32_t u) const { return mSegment[u].incidenceSize; }

    // How many original vertices u stands for. One until mass elimination merges into it, so a
    // degree that counts vertices has to count this rather than entries.
    //
    // Held in a flat array, and this is the condition under which the prototypes said it would
    // have to be: its members are a chain rather than a list, so the size is no longer free to
    // read. The array is here for the chain and not for locality.
    // The MAGNITUDE, so that a driver never has to know whether the sign encoding is in force. It
    // is not, at this commit, and the cast is a widening of a value that is always positive.
    std::uint32_t weight(std::int32_t u) const {
        const std::int32_t w = mWeight[u];
        return static_cast<std::uint32_t>(w < 0 ? -w : w);
    }

    // RESTORING A MEMBER'S SIGN, for a driver that runs mass elimination LATE and so owns the
    // restore; see massEliminate. `AMD_2` is that case, restoring under RESTORE DEGREE LISTS in
    // the pass that computes final degrees and refiles, so the store rides on a load that pass
    // makes anyway.
    

    void formReachableSetMmd(std::int32_t u, std::vector<std::int32_t>& reachableSet);
    // reach(u), as above. Not const: the mark array and its tag are scratch, and threading them
    // through every call site is what the prototypes do only because their display functions
    // needed to borrow them.
    void formReachableSetAmd(std::int32_t u, std::vector<std::int32_t>& reachableSet);

    // The same set, weighted: the number of ORIGINAL vertices reach(u) stands for. Once a branch
    // merges into live vertices, a reached vertex can stand for several, and a degree that counts
    // vertices has to count those rather than entries.
    std::uint32_t reachableSetWeight(std::int32_t u);

    // Eliminate the pivot: turn it into a clique, absorb the cliques it belonged to, prune the
    // edges the new clique implies, and merge in whatever it makes indistinguishable. Returns the
    // vertices merged into the pivot's supervariable.
    //
    // **The returned reference is valid until the next elimination**, being a scratch buffer whose
    // capacity survives from pivot to pivot. A caller that needs the list to outlive the next call
    // copies it.
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
    // The new clique is C[p] and gets no name of its own. None of the three differences builds a
    // set; each is a stamp of the subtrahend and one compaction pass over the minuend.
    const std::vector<std::int32_t>& eliminateMmd(std::int32_t pivot);

    // The same elimination with an approximate-degree driver's first scan folded in. That driver
    // walks exactly the lists the prune has just walked, so folding leaves A[u] visited once and
    // I[u] twice.
    //
    // IT CANNOT BE FOLDED FURTHER: a clique's outside count is complete only once every member of
    // C[p] has been seen, so the sum the bound needs is a second pass by construction.
    //
    // TWO THINGS MAKE IT SOUND AND NEITHER SURVIVES REORDERING THE PHASES. The scan runs over the
    // UNTRIMMED C[p], and the mass-eliminated contribute nothing because the merge test requires
    // I[u] == {pivot} and the scan skips the pivot. And the weights summed over A[u] are read
    // before mass elimination could change any of them, the prune having already removed every
    // member of C[p] from A[u].

    // The same fusion for the tagged-W encoding, which is the one `AmdFlat` carries. Everything the
    // overload above says about why the fold is sound applies here unchanged; only the three-facts
    // array differs. It also fills the hash key's adjacency half, which that overload has no need
    // of because its callers keep a separate walk of A[u].
    const std::vector<std::int32_t>& eliminateAmd(std::int32_t pivot, TaggedScan& scan);


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

    // Walk I[u] from the back in the reachable-set walk, matching genmmd's clique stack. A
    // tie-break
    // convention and nothing else: it changes which permutation comes out, never which sets are
    // computed. See the member's note.
    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    // Stop the eliminator at the prune, leaving mass elimination to the caller. AMD_2 makes the
    // same test in its scan 2, AFTER aggressive absorption has dropped every clique lying inside
    // the new one, and says why in its own comment: with aggressive absorption, `deg == 0` is
    // identical to the structural test. Asking first, which is what the eliminator does by
    // default,
    // asks it of an I[u] that still holds cliques about to be removed, so the cheap test declines
    // vertices AMD merges.
    //
    // With this on, eliminateAmd returns an EMPTY merged list and C[pivot] is reach(pivot) exactly,
    // and the caller must call massEliminate() once it has absorbed. Used by AmdFlat alone. See
    // experiments/ordering/AMD3.md, ledger entry 3.
    void setLateMassElimination(bool on) { mLateMassElimination = on; }

    // The half eliminateAmd no longer does under the flag above: fold into the pivot's
    // supervariable every member of C[pivot] that the new clique now accounts for entirely, and
    // trim C[pivot] of them. Returns the merged vertices, from a member scratch as the
    // eliminator does.
    //
    //     merged   = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
    //     C[pivot] = C[pivot] - merged
    //
    // Calling it without the flag is a caller error and not guarded: eliminateAmd will already have
    // merged, so this would find nothing and cost a pass.
    const std::vector<std::int32_t>& massEliminate(std::int32_t pivot);

    // Kill these cliques and take them out of the incidence lists of these vertices. A clique
    // dies when it is found to lie wholly inside a newer one, which is aggressive absorption, and
    // the vertices to purge are the newer clique's members, since those are the only lists that
    // can still name it.
    // `vertexCount` is `|C[p]|`, so one dimensional; every caller passes `cliqueSize(pivot)`.
    void absorbAggressively(const std::vector<std::int32_t>& cliques,
                            const std::int32_t* vertices, std::uint32_t vertexCount);

    // Expand a pivot sequence into an elimination order over the original vertices. A pivot
    // stands for its whole supervariable, whose members are eliminated consecutively, so this is
    // where a supervariable of size w becomes w columns.
    std::vector<std::int32_t> order(const std::vector<std::int32_t>& pivots) const;

    // The same permutation with each supervariable's members in ASCENDING VERTEX INDEX rather
    // than merge order. Indistinguishable members, so the fill and the forest are unchanged;
    // only the permutation is, and it is genmmd's. Used by MmdFlat alone. See the .cpp.
    std::vector<std::int32_t> orderAscending(const std::vector<std::int32_t>& pivots) const;

private:
    void beginEliminationMmd(std::int32_t pivot);
    // THE THREE STEPS OF AN ELIMINATION, private because between them the graph is half
    // eliminated: the clique is written and stamped but the reached vertices still name the pivot
    // as a variable. Nothing outside may observe that state, which is why the seam is private
    // calls rather than a public begin and end. The two wrappers above are the only callers, and
    // their value is that the ORDER of the three lives here where nothing in a driver could
    // enforce it.
    //
    // ONLY THE AMD BEGIN TAKES THE SCAN, and it can because the begin is already split by branch:
    // the mmd driver reaches `beginEliminationMmd` and never sees a `TaggedScan`. That was not
    // true while there was one `beginElimination`.
    void beginEliminationAmd(std::int32_t pivot, TaggedScan& scan);

    void pruneMmd(std::int32_t pivot);
    void pruneAmd(std::int32_t pivot, TaggedScan& scan);

    const std::vector<std::int32_t>& finishElimination(std::int32_t pivot);

    // Dimensions.
    std::size_t mSize = 0;   // number of vertices

    // A[u] and I[u] share one run, and one array holds every run end to end. Their number is
    // CONSERVED: an elimination that reaches u replaces at least one source with the new clique
    // and never manufactures one, so
    //
    //     |A[u]| + |I[u]| <= the number of off-diagonal entries in u's column of A
    //
    // holds for the whole run and u's block is sized once from the pattern and never grown.
    // Section 5.3 of archive/sparse_factorization.md carries the argument.
    //
    // THE ORDER WITHIN THE RUN IS FORCED. The prune compacts A[u] then I[u], and A[u] shrinks by
    // at least what I[u] gains, so the incidence write always trails the read cursor. The other
    // order would overwrite an adjacency entry not yet read.
    //
    // C[c] has no such bound and gets its own store below.
    std::vector<std::int32_t> mAdjIncSrc;   // every A[u] then I[u], run after run
    // ONE OBJECT PER VERTEX, NOT THREE ARRAYS. The three numbers are never useful apart: any walk
    // of u needs where its run starts and at least one length. Held together they are one cache
    // line; held apart, a walk of C[pivot] pulls three per member.
    //
    // EXACTLY 16 BYTES, four to a line, and the layout only pays while that holds. The position is
    // `std::size_t`, offsetting into a store bounded by nnz(A); the two lengths are
    // `std::uint32_t`, one dimensional and bounded by n.
    //
    // TWO LIVES PER FIELD, decided by whether this id has formed a clique. The turn happens in
    // exactly one place, `beginElimination`, which rewrites the first two fields once both readers
    // of their vertex meaning are past.
    struct Segment {
        // WHERE THE BLOCK STARTS, AND IN WHICH ARRAY. A live vertex's segment is a position in
        // `mAdjIncSrc`, fixed at construction; a clique's is a position in `mCliqueSrc`, written
        // when the clique is born. Nothing in the field says which, so `clique(c)` is the only
        // accessor that adds the clique store's base, and it is reached only for an id the caller
        // already knows is a clique.
        std::size_t   srcPtr;
        // A[u]'S LENGTH for a vertex, THE MEMBER COUNT for a clique. ZERO MEANS A DIFFERENT THING
        // IN EACH, and both readings are live: on a vertex it is an empty adjacency, which is what
        // `massEliminate`'s first conjunct tests, and on a clique it is DEAD, written by
        // `killClique` and read by the prune's incidence compaction. Do not fold the two.
        std::uint32_t adjacencySize;
        // I[u]'S LENGTH, immediately behind A[u]. Zeroed when the clique is born and never read
        // again, a clique having no incidence part.
        std::uint32_t incidenceSize;
    };
    std::vector<Segment> mSegment;

    // Every C[c] ever formed, end to end, appended when a clique is born and never reclaimed: an
    // absorbed clique leaves a hole. The total is bounded by nnz(L).
    //
    // POSITIONS AND NOT POINTERS is what lets this store grow. A reallocation leaves every
    // `mSegment[c].srcPtr` valid where a pointer would dangle.
    //
    // A CLIQUE'S DESCRIPTOR IS ITS DEAD PIVOT'S RUN: `mSegment[c].srcPtr` is where C[c] starts
    // here and `.adjacencySize` is how much of it is still live. Sound because a clique's id is
    // the id of the pivot that formed it, and that vertex's own A[c] and I[c] are read for the
    // last time inside `beginElimination`.
    //
    // DO NOT CARVE THIS AND ITS NEIGHBORS OUT OF ONE BLOCK. Each array here is its own vector.
    // Consolidating them at offsets that are multiples of n makes them share cache sets whenever
    // n is a power of two.
    std::vector<std::int32_t> mCliqueSrc;   // every C[c] ever formed, end to end

    // See `numPeakCliqueMembers`. Maintained by `killClique`, by `trimClique`, and by the one
    // place a clique is born, which is why all three are funnelled rather than written where they
    // happen.
    std::size_t mNumLiveCliqueMembers = 0;
    std::size_t mNumPeakCliqueMembers = 0;
#ifndef NDEBUG
    std::vector<std::int32_t> mCliqueOwners;   // every vertex that ever formed one; see the check
#endif

    // The supervariable a vertex stands for, as a chain rather than a list per vertex:
    // `mSuperNext` links the members and `mSuperLast` names the tail, so a merge is O(1).
    std::vector<std::int32_t> mSuperNext;
    std::vector<std::int32_t> mSuperLast;
    std::vector<std::int32_t> mWeight;
    std::uint32_t              mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

    // Which end of I[u] the reachable-set walk starts from. genmmd threads its clique list
    // through an
    // integer array and pushes at the head, `list[nb] = el; el = nb`, then reads from the head, so
    // the clique seen LAST is expanded FIRST; we hold a vector and append. Same set either way and
    // the same cost, but the order decides C[pivot]'s order, hence which of two equal-degree
    // candidates a later iteration finds first, and minimum degree is settled by exactly that.
    // Off by default, so every existing driver is unaffected; MmdFlat turns it on. See
    // experiments/ordering/mmd3.py, where the same four walks are reversed together.
    bool mReverseIncidence = false;

    // Whether the eliminator stops at the prune and leaves mass elimination to the caller. Off by
    // default. See the setter, and massEliminate() for the half it hands over.
    bool mLateMassElimination = false;

    // WHETHER ANY VERTEX HAS BEEN NUMBERED WITHOUT BEING ELIMINATED, which only the mmd prepass
    // does. While it is false the walks skip their GONE test.
    void markGone(std::int32_t u) { if (!mMark.empty()) mMark[u] = GONE; }

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
// the public ordering signature; `gAmdCompactions` is here for the same reason. Not thread safe,
// and it does not need to be: nothing writes it outside a test.
inline std::size_t gPeakCliqueMembers = 0;

// ALLOCATED ON DEMAND, and the amd branch never calls this. `NIL` rather than zero as the initial
// stamp, since zero is a tag a walk can reach and this array's whole job is to answer "have I seen
// v this step" against `mTag`, which starts there.
inline void QuotientGraph::enableMarks() {
    mMark.assign(mSize, NIL);
}

inline QuotientGraph::QuotientGraph(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mSize(colPtr.empty() ? 0 : colPtr.size() - 1), mSegment(mSize) {
    const std::int32_t size = static_cast<std::int32_t>(mSize);

    // One pass to place each column's run, dropping the diagonal. The runs are laid out in column
    // order and never move afterwards, so the offsets are written once here and only the lengths
    // change. The run is u's whole allowance: A[u] fills it now and I[u] grows into the room A[u]
    // gives up, never past mSegment[u + 1].srcPtr.
    mAdjIncSrc.reserve(colPtr.empty() ? 0 : colPtr.back());
    for (std::int32_t aj = 0; aj < size; ++aj) {
        // The run's start is the store's length before this column is appended.
        mSegment[aj].srcPtr = mAdjIncSrc.size();
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mAdjIncSrc.push_back(rowIdx[cp]);
        // THE ONE CROSSING. A difference of two positions is a count, so this is the single place
        // in the class where a two-dimensional quantity is written into a one-dimensional one. It
        // is bounded by deg(aj) and so by n, which the SparseMatrix constructor has already capped
        // at MAX_IDX, but the cast is written rather than left implicit because that bound is an
        // argument and not something the types say.
        mSegment[aj].adjacencySize =
            static_cast<std::uint32_t>(mAdjIncSrc.size() - mSegment[aj].srcPtr);
    }

    // The clique store is reserved from the source pool's size so that it cannot reallocate under
    // a walk; see beginElimination.
    mCliqueSrc.reserve(colPtr.empty() ? 0 : colPtr.back());

    // Every vertex begins as a supervariable of one: a chain holding itself, so next is NIL and
    // last is the vertex. Mass elimination splices these together and order() walks them.
    mSuperNext.assign(size, NIL);
    mSuperLast.resize(size);
    mWeight.assign(size, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

// APPENDS to `reached` and does not clear it, which is what lets beginElimination point it straight
// at the clique arena instead of at a scratch. The returning overload below clears for itself.
inline void QuotientGraph::formReachableSetMmd(std::int32_t u,
                                               std::vector<std::int32_t>& reachableSet) {
    // reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}
    //
    // The mark array is the set: `mMark[v] == mTag` is the membership test and `mMark[v] = mTag`
    // the insertion, so the union costs one pass per source and comes out in walk order. The
    // buffer is the caller's, so a caller in a loop can keep one.
    //
    // ELIMINATED VERTICES ARE SKIPPED RATHER THAN PURGED, since a live merge leaves the vertex it
    // folds away where it lies and every clique that named it still does. `mMark[v] < mTag`
    // answers liveness and membership from one load, which is why GONE must sort above every tag,
    // and GONE is WRITTEN at every death site rather than inferred from a value.
    //
    // THE ADJACENCY LOOPS ASK mMark AND THE CLIQUE LOOPS DO NOT, and the asymmetry is exact.
    // `number()` leaves a prepass vertex at weight one and in the adjacency of every neighbor, so
    // a positive weight does not mean live there. It cannot appear in a CLIQUE: the prepass
    // completes before the first elimination and every clique since is built from a reach that
    // skipped it.
    mWeight[u] = -mWeight[u];              // never its own neighbor
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mWeight, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique.
    const Segment&      segment       = mSegment[u];          // one fetch; see the member
    const std::int32_t* source        = mAdjIncSrc.data() + segment.srcPtr;
    const std::uint32_t adjacencySize = segment.adjacencySize;
    const std::uint32_t incidenceSize = segment.incidenceSize;
    const std::int32_t* incidence     = source + adjacencySize;
    // Hoisted like `live`, and for the same reason: member loads the compiler cannot prove are
    // unaliased by the stores below. Both branches are per SOURCE rather than per member, so
    // neither sits in the loop that does the work. See the members' notes for what each decides.
    const bool reverse  = mReverseIncidence;

    // Which source is walked first. genmmd expands the variables and then the cliques, which is
    // how the whole md ladder is laid out; AMD_2 takes the cliques first and the supervariables
    // only on its last pass. Same set either way, and the order decides C[pivot]'s content order,
    // hence which of two equal-degree candidates a later iteration finds first.

    for (std::uint32_t vk = 0; vk < adjacencySize; ++vk) {
        const std::int32_t v  = source[vk];
        const std::int32_t vWeight = mWeight[v];
        if (vWeight > 0 && !(mHasNumbered && mMark[v] == GONE)) {
            mWeight[v] = -vWeight; reachableSet.push_back(v);
        }
    }
    for (std::uint32_t ck = 0; ck < incidenceSize; ++ck) {
        const std::uint32_t t           = reverse ? incidenceSize - 1 - ck : ck;
        const std::int32_t  c           = incidence[t];
        const std::int32_t* members     = mCliqueSrc.data() + mSegment[c].srcPtr;
        const std::uint32_t membersSize = mSegment[c].adjacencySize;
        for (std::uint32_t vk = 0; vk < membersSize; ++vk) {
            const std::int32_t v  = members[vk];
            const std::int32_t vWeight = mWeight[v];
            if (vWeight > 0) { mWeight[v] = -vWeight; reachableSet.push_back(v); }
        }
    }
}

// APPENDS to `reached` and does not clear it, which is what lets beginElimination point it straight
// at the clique arena instead of at a scratch. The returning overload below clears for itself.
inline void QuotientGraph::formReachableSetAmd(std::int32_t u,
                                               std::vector<std::int32_t>& reachableSet) {
    // reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}
    //
    // The mark array is the set: `mMark[v] == mTag` is the membership test and `mMark[v] = mTag`
    // the insertion, so the union costs one pass per source and comes out in walk order. The
    // buffer is the caller's, so a caller in a loop can keep one.
    //
    // ELIMINATED VERTICES ARE SKIPPED RATHER THAN PURGED, since a live merge leaves the vertex it
    // folds away where it lies and every clique that named it still does. `mMark[v] < mTag`
    // answers liveness and membership from one load, which is why GONE must sort above every tag,
    // and GONE is WRITTEN at every death site rather than inferred from a value.
    //
    // THE ADJACENCY LOOPS ASK mMark AND THE CLIQUE LOOPS DO NOT, and the asymmetry is exact.
    // `number()` leaves a prepass vertex at weight one and in the adjacency of every neighbor, so
    // a positive weight does not mean live there. It cannot appear in a CLIQUE: the prepass
    // completes before the first elimination and every clique since is built from a reach that
    // skipped it.
    mWeight[u] = -mWeight[u];              // never its own neighbor
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mWeight, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique.
    const Segment&      segment       = mSegment[u];          // one fetch; see the member
    const std::int32_t* source        = mAdjIncSrc.data() + segment.srcPtr;
    const std::uint32_t adjacencySize = segment.adjacencySize;
    const std::uint32_t incidenceSize = segment.incidenceSize;
    const std::int32_t* incidence     = source + adjacencySize;
    // Hoisted like `live`, and for the same reason: member loads the compiler cannot prove are
    // unaliased by the stores below. Both branches are per SOURCE rather than per member, so
    // neither sits in the loop that does the work. See the members' notes for what each decides.
    const bool reverse  = mReverseIncidence;

    // Which source is walked first. genmmd expands the variables and then the cliques, which is
    // how the whole md ladder is laid out; AMD_2 takes the cliques first and the supervariables
    // only on its last pass. Same set either way, and the order decides C[pivot]'s content order,
    // hence which of two equal-degree candidates a later iteration finds first.
    for (std::uint32_t ck = 0; ck < incidenceSize; ++ck) {
        const std::uint32_t t           = reverse ? incidenceSize - 1 - ck : ck;
        const std::int32_t  c           = incidence[t];
        const std::int32_t* members     = mCliqueSrc.data() + mSegment[c].srcPtr;
        const std::uint32_t membersSize = mSegment[c].adjacencySize;
        for (std::uint32_t vk = 0; vk < membersSize; ++vk) {
            const std::int32_t v  = members[vk];
            const std::int32_t vWeight = mWeight[v];
            if (vWeight > 0) { mWeight[v] = -vWeight; reachableSet.push_back(v); }
        }
    }
    for (std::uint32_t vk = 0; vk < adjacencySize; ++vk) {
        const std::int32_t v  = source[vk];
        const std::int32_t vWeight = mWeight[v];
        if (vWeight > 0) { mWeight[v] = -vWeight; reachableSet.push_back(v); }
    }
}



inline std::uint32_t QuotientGraph::reachableSetWeight(std::int32_t u) {
    // A sum over DISTINCT vertices, so bounded by n; see the header.
    std::uint32_t totalWeight = 0;
    ++mTag;
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique.
    const Segment&      segment       = mSegment[u];          // one fetch; see the member
    const std::int32_t* source        = mAdjIncSrc.data() + segment.srcPtr;
    const std::uint32_t adjacencySize = segment.adjacencySize;
    const std::uint32_t incidenceSize = segment.incidenceSize;
    for (std::uint32_t vk = 0; vk < adjacencySize; ++vk) {
        const std::int32_t v = source[vk];
        if (mMark[v] != GONE) {
            mMark[v] = mTag; totalWeight += static_cast<std::uint32_t>(mWeight[v]);
        }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::uint32_t ck = 0; ck < incidenceSize; ++ck) {
        const std::int32_t  c           = incidence[ck];
        const std::int32_t* members     = mCliqueSrc.data() + mSegment[c].srcPtr;
        const std::uint32_t membersSize = mSegment[c].adjacencySize;
        for (std::uint32_t vk = 0; vk < membersSize; ++vk) {
            const std::int32_t v = members[vk];
            if (mMark[v] < mTag) {
                mMark[v] = mTag; totalWeight += static_cast<std::uint32_t>(mWeight[v]);
            }
        }
    }
    return totalWeight;
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
    markGone(u);
}

inline void QuotientGraph::beginEliminationMmd(std::int32_t pivot) {
    // The reach is written STRAIGHT INTO THE ARENA. C[pivot] is the reach, so the block the walk
    // fills is the clique's own block and no scratch is needed.
    //
    // THE ARENA MUST NOT MOVE WHILE THAT WALK RUNS, which is what the reserve below is for and is
    // not an optimization. The walk reads each clique's members through a pointer into the arena
    // while appending to that same arena, so a growth would leave every such pointer dangling. A
    // reach is at most `mSize` entries, so room for one is room for the whole walk.
    if (mCliqueSrc.capacity() - mCliqueSrc.size() < mSize)
        mCliqueSrc.reserve(std::max(2 * mCliqueSrc.capacity(), mCliqueSrc.size() + mSize));

    // THE DESCRIPTOR IS HELD IN LOCALS UNTIL BOTH READERS OF THE PIVOT'S RUN ARE PAST, which is
    // the whole subtlety of storing a clique in the run of the vertex that formed it. The slots it
    // will occupy are still being read: `reachableSet` walks A[pivot] and I[pivot] through them,
    // and the absorbed-clique loop below finds I[pivot] the same way. Writing either early leaves
    // a walk reading the arena through an offset into mAdjIncSrc. See the header.
    const std::size_t cliqueStart = mCliqueSrc.size();
    formReachableSetMmd(pivot, mCliqueSrc);            // appends; see its note
    // THE SECOND AND LAST CROSSING. The arena's new length less this block's start is a member
    // count, so a two-dimensional quantity is written into a one-dimensional one, exactly as the
    // constructor does for the source runs. Bounded by n, a reach having at most n entries, which
    // is the same bound the reserve above relies on.
    const std::uint32_t cliqueLen = static_cast<std::uint32_t>(mCliqueSrc.size() - cliqueStart);

    // Taken AFTER the append, since that is what can move the arena. With the reserve above it
    // cannot have moved, and this stays as it is regardless: it costs nothing and it is the shape
    // that remains correct if the reserve is ever revised.
    const std::int32_t* reached     = mCliqueSrc.data() + cliqueStart;
    const std::uint32_t reachedSize = cliqueLen;

    // The absorbed cliques are I[pivot], read where they lie. Nothing below writes the pivot's run
    // (the prune rewrites the runs of C[pivot]'s members, and the pivot is not one of them), so no
    // copy and no scratch is needed to keep them alive across the passes that follow.
    const std::int32_t* absorbedCliques =
        mAdjIncSrc.data() + mSegment[pivot].srcPtr + mSegment[pivot].adjacencySize;
    const std::uint32_t absorbedSize    = mSegment[pivot].incidenceSize;
    for (std::uint32_t ck = 0; ck < absorbedSize; ++ck)
        killClique(absorbedCliques[ck]);              // dead, its block left behind

    // Both readers of the pivot's own run are past, so the run becomes the clique's descriptor.
    mSegment[pivot].srcPtr     = cliqueStart;
    mSegment[pivot].adjacencySize = cliqueLen;

    // A CLIQUE IS BORN HERE AND NOWHERE ELSE, which is what makes the peak countable. See
    // `numPeakCliqueMembers`.
    mNumLiveCliqueMembers += cliqueLen;
    mNumPeakCliqueMembers  = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);
#ifndef NDEBUG
    mCliqueOwners.push_back(pivot);
#endif

    // Stamp the new clique. THE ABSORBED CLIQUES ARE NOT STAMPED: their members are already in
    // the reach, and stamping them would make a later walk skip entries it must still see.
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t vk = 0; vk < reachedSize; ++vk)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[vk]]);
    mCliqueWeight = cliqueWeight;
}

inline void QuotientGraph::beginEliminationAmd(std::int32_t pivot, TaggedScan& scan) {
    // The reach is written STRAIGHT INTO THE ARENA. C[pivot] is the reach, so the block the walk
    // fills is the clique's own block and no scratch is needed.
    //
    // THE ARENA MUST NOT MOVE WHILE THAT WALK RUNS, which is what the reserve below is for and is
    // not an optimization. The walk reads each clique's members through a pointer into the arena
    // while appending to that same arena, so a growth would leave every such pointer dangling. A
    // reach is at most `mSize` entries, so room for one is room for the whole walk.
    if (mCliqueSrc.capacity() - mCliqueSrc.size() < mSize)
        mCliqueSrc.reserve(std::max(2 * mCliqueSrc.capacity(), mCliqueSrc.size() + mSize));

    // THE DESCRIPTOR IS HELD IN LOCALS UNTIL BOTH READERS OF THE PIVOT'S RUN ARE PAST, which is
    // the whole subtlety of storing a clique in the run of the vertex that formed it. The slots it
    // will occupy are still being read: `reachableSet` walks A[pivot] and I[pivot] through them,
    // and the absorbed-clique loop below finds I[pivot] the same way. Writing either early leaves
    // a walk reading the arena through an offset into mAdjIncSrc. See the header.
    const std::size_t cliqueStart = mCliqueSrc.size();
    formReachableSetAmd(pivot, mCliqueSrc);            // appends; see its note
    // THE SECOND AND LAST CROSSING. The arena's new length less this block's start is a member
    // count, so a two-dimensional quantity is written into a one-dimensional one, exactly as the
    // constructor does for the source runs. Bounded by n, a reach having at most n entries, which
    // is the same bound the reserve above relies on.
    const std::uint32_t cliqueLen = static_cast<std::uint32_t>(mCliqueSrc.size() - cliqueStart);

    // Taken AFTER the append, since that is what can move the arena. With the reserve above it
    // cannot have moved, and this stays as it is regardless: it costs nothing and it is the shape
    // that remains correct if the reserve is ever revised.
    const std::int32_t* reached     = mCliqueSrc.data() + cliqueStart;
    const std::uint32_t reachedSize = cliqueLen;

    // The absorbed cliques are I[pivot], read where they lie. Nothing below writes the pivot's run
    // (the prune rewrites the runs of C[pivot]'s members, and the pivot is not one of them), so no
    // copy and no scratch is needed to keep them alive across the passes that follow.
    //
    // A CLIQUE DIES TWO WAYS AND THE TAGGED W MUST LEARN ABOUT BOTH. Aggressive absorption zeroes
    // `w[c]` in the driver; elimination-time absorption is this list. `AMD_2` writes both deaths
    // into W, `Pe[e] = FLIP(me)` with `W[e] = 0`, and its scan then tests `we != 0` off the load
    // it already needs for the value. The store rides on the walk that kills the clique, so one
    // read of the entry records both facts.
    const std::int32_t* absorbedCliques =
        mAdjIncSrc.data() + mSegment[pivot].srcPtr + mSegment[pivot].adjacencySize;
    const std::uint32_t absorbedSize    = mSegment[pivot].incidenceSize;
    for (std::uint32_t ck = 0; ck < absorbedSize; ++ck) {
        const std::int32_t c = absorbedCliques[ck];
        killClique(c);                               // dead, its block left behind
        scan.w[c] = 0;                               // and dead to the scan; see above
    }

    // Both readers of the pivot's own run are past, so the run becomes the clique's descriptor.
    mSegment[pivot].srcPtr     = cliqueStart;
    mSegment[pivot].adjacencySize = cliqueLen;

    // A CLIQUE IS BORN HERE AND NOWHERE ELSE, which is what makes the peak countable. See
    // `numPeakCliqueMembers`.
    mNumLiveCliqueMembers += cliqueLen;
    mNumPeakCliqueMembers  = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);
#ifndef NDEBUG
    mCliqueOwners.push_back(pivot);
#endif

    // Stamp the new clique. THE ABSORBED CLIQUES ARE NOT STAMPED: their members are already in
    // the reach, and stamping them would make a later walk skip entries it must still see.
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t vk = 0; vk < reachedSize; ++vk)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[vk]]);
    mCliqueWeight = cliqueWeight;
}

inline void QuotientGraph::pruneMmd(std::int32_t pivot) {
    // C[pivot] IS the reach here: beginElimination wrote it there and nothing has trimmed it yet
    // (massEliminate does, and runs after). So the prune walks the arena block rather than a copy.
    const std::int32_t* reached     = mCliqueSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t reachedSize = mSegment[pivot].adjacencySize;

    // Both lists are compacted in place rather than rebuilt into a scratch: every pass here only
    // ever removes, so survivors are written over entries already read and nothing is allocated.
    //
    // Both lists live in u's one run, the adjacency first and the incidence into whatever the
    // adjacency has given up. The pivot is APPENDED and the two boundary entries swapped
    // afterwards; it cannot be written first, because the write cursor starts at `kept` and the
    // read at the original adjacency length, and those are equal whenever nothing was pruned.
    for (std::uint32_t uk = 0; uk < reachedSize; ++uk) {
        const std::int32_t u = reached[uk];
        // ONE FETCH OF THE RUN, not three. The three numbers share a 16-byte object, so the
        // reference below brings all of them in on one line; read as `mSegment[u].srcPtr` and
        // friends they would still be three subscripts of the same object and the compiler would
        // still load once, but naming it once is what makes that visible to a reader.
        const Segment&      segment       = mSegment[u];
        std::int32_t*       source        = mAdjIncSrc.data() + segment.srcPtr;
        // The two counters are one-dimensional COUNTS, positions in a list bounded by deg(u) and
        // so by n, where std::size_t is for a position into an n x n object. They take the type of
        // what they count, so they move with the array and no cast is needed here.
        const std::uint32_t adjacencySize = segment.adjacencySize;
        std::uint32_t       kept          = 0;
        for (std::uint32_t vk = 0; vk < adjacencySize; ++vk) {
            const std::int32_t v = source[vk];
            // ONE LOAD, THREE QUESTIONS, which is Amd.cpp's `nvj = Nv [j] ; if (nvj > 0)`. A
            // negative weight is a member of the new clique, the pivot included, so the explicit
            // `v == pivot` test goes with the membership test; a zero is a vertex a live merge
            // folded away. The FOURTH question, whether v was numbered by a prepass, still needs
            // mMark, and dropping it is not an option: massEliminate reads `adjacencySize == 0`,
            // so a numbered leftover in A[u] would suppress a merge.
            if (mWeight[v] <= 0) continue;             // in the new clique, the pivot, or merged
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by a prepass; see above
            source[kept++] = v;
        }
        mSegment[u].adjacencySize = kept;                      // A[u] - C[pivot] - {pivot}

    // The read cursor is a hoisted POINTER rather than source[adjacencySize + i]. The read and
    // the write are into the same buffer, so every conditional store orders the next load behind
    // it; that is the price of compacting in place rather than into a scratch.
        const std::int32_t* incidence     = source + adjacencySize;
        const std::uint32_t incidenceSize = segment.incidenceSize;
        std::uint32_t       write         = kept;      // follows `kept`; see the note above
        for (std::uint32_t ck = 0; ck < incidenceSize; ++ck) {
            const std::int32_t c = incidence[ck];
            if (mSegment[c].adjacencySize != 0) source[write++] = c;  // dead is size zero; above
        }
        source[write++]   = pivot;                     // u joins the new clique, id = pivot
        mSegment[u].incidenceSize = write - kept;
        // [c1, ..., ck, pivot] to [pivot, c2, ..., ck, c1], which is a SWAP of the two boundary
        // entries: only they move and c2..ck stay put. The pivot cannot be written first, which
        // is what the loop above would otherwise allow: the write cursor starts at `kept` and the
        // read at the original adjacencySize, and those are equal whenever nothing was pruned from
        // A[u], so an extra write before the reads finish clobbers an unread entry. AMD_2 makes
        // its three assignments after both compactions for exactly this reason.
    }
    // THE SIGNS COME BACK HERE, at the end of the prune, which is the LAST READER of them. The
    // walk that built C[pivot] marked membership by negating a weight and the loop above is what
    // consumed that mark, `mWeight[v] <= 0` being "in the new clique, the pivot, or merged".
    // Nothing downstream asks: absorption never touches a weight, mass elimination's merge test is
    // structural, and every other reader takes a magnitude through `weight()`.
    mWeight[pivot] = -mWeight[pivot];
    for (std::uint32_t vk = 0; vk < mSegment[pivot].adjacencySize; ++vk)
        mWeight[mCliqueSrc[mSegment[pivot].srcPtr + vk]] =
            -mWeight[mCliqueSrc[mSegment[pivot].srcPtr + vk]];
}

inline void QuotientGraph::pruneAmd(std::int32_t pivot, TaggedScan& scan) {
    const std::int32_t* reached     = mCliqueSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t reachedSize = mSegment[pivot].adjacencySize;

    // EVERY MEMBER OF C[pivot] LEAVES THE DEGREE LISTS HERE, which is what frees each member's
    // mPrev and mNext for the hash key and the hash chain. A driver that refiles inside its own
    // bound pass passes no buckets, its links still being degree links when the hash runs.
    if (scan.buckets != nullptr)
        for (std::uint32_t uk = 0; uk < reachedSize; ++uk) scan.buckets->unfile(reached[uk]);

    const std::int32_t  wflg        = scan.wflg;          // hoisted, as the flags are

    for (std::uint32_t uk = 0; uk < reachedSize; ++uk) {
        const std::int32_t u       = reached[uk];
        // NEGATED, BECAUSE u IS A MEMBER OF C[pivot] AND SO READS NEGATIVE. This is Amd.cpp's
        // `nvi = -Nv [i]` under CONSTRUCT NEW CLIQUE, and the sign is the whole of it: `wnvi`
        // below is `wflg - nvi` and comes out wrong by twice the weight if the magnitude is not
        // taken. No cast is needed, `mWeight` being signed;
        // its own comment called it a signedness cast rather than a narrowing one, which was the
        // code saying the type was wrong. `wnvi` must still be able to go negative, which is
        // Amd.cpp's convention and the reason `w` is signed.
        const std::int32_t nvi     = -mWeight[u];
        const std::int32_t wnvi    = wflg - nvi;          // Amd.cpp's wnvi, and signed for it
        // ONE FETCH OF THE RUN, not three. The three numbers share a 16-byte object, so the
        // reference below brings all of them in on one line; read as `mSegment[u].srcPtr` and
        // friends they would still be three subscripts of the same object and the compiler would
        // still load once, but naming it once is what makes that visible to a reader.
        const Segment&      segment       = mSegment[u];
        std::int32_t*       source        = mAdjIncSrc.data() + segment.srcPtr;
        const std::uint32_t adjacencySize = segment.adjacencySize;
        std::uint32_t       kept          = 0;
        std::uint32_t       explicitPart  = 0;            // a weight sum, not a count of positions
    // The key is a SUM over the pruned A[u] and the final I[u], reduced modulo the driver's
    // bucket count. It WRAPS in uint32 deliberately and is reduced once, not per term. The
    // modulus must not divide any stride applied to a term, which is why no term carries one.
        std::uint32_t       key           = 0;            // wraps, like Amd.cpp's UInt hval
        std::int32_t        heldVertex    = NIL;          // see the plain prune above
        for (std::uint32_t vk = 0; vk < adjacencySize; ++vk) {
            const std::int32_t v = source[vk];
            if (mWeight[v] <= 0) continue;             // see the plain prune above
            explicitPart += static_cast<std::uint32_t>(mWeight[v]);
            key += static_cast<std::uint32_t>(v);         // no + 1, no reduction; see above
            if (heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mSegment[u].adjacencySize       = kept;
        // THE ADJACENCY HALF GOES INTO w[u], NOT INTO AN ARRAY OF ITS OWN. `w` is indexed by
        // CLIQUE id and a clique id is a dead pivot's id, so for a LIVE vertex the slot carries
        // nothing: u is in C[pivot] and therefore alive, and no clique is named after it until it
        // is eliminated, by which time this value is long consumed. The tagged W thus answers a
        // FOURTH question on top of the other three, and the array that carried this one is
        // gone. It cannot collide with the clique writes below: those are indexed by c drawn from
        // I[u], every one of which is a dead pivot, and u is live.
        //
        // The driver's obligation is one store: the slot goes back to alive-and-unseen once the
        // bound has been read. See src/AmdFlat.cpp.
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
        const std::uint32_t incidenceSize = segment.incidenceSize;
        std::uint32_t       write         = kept;         // follows `kept`; see the plain prune
        for (std::uint32_t ck = 0; ck < incidenceSize; ++ck) {
            const std::int32_t c = incidence[ck];
            // Amd.cpp's four lines, transcribed. A clique seen earlier in this step already holds
            // the running value above the tag; one seen for the first time starts from |C[c]| and
            // is listed once; one already absorbed reads ZERO and is dropped from the list here.
            //
            // ONE LOAD, TWO QUESTIONS. `we == 0` is dead, which would otherwise be a
            // `mCliqueSize` probe,
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
        mSegment[u].incidenceSize = write - kept;
        if (write - kept > 1) std::swap(source[kept], source[write - 1]);
    }
    // THE SIGNS COME BACK HERE, at the end of the prune, which is the LAST READER of them. The
    // walk that built C[pivot] marked membership by negating a weight and the loop above is what
    // consumed that mark, `mWeight[v] <= 0` being "in the new clique, the pivot, or merged".
    // Nothing downstream asks: absorption never touches a weight, mass elimination's merge test is
    // structural, and every other reader takes a magnitude through `weight()`.
    mWeight[pivot] = -mWeight[pivot];
    for (std::uint32_t vk = 0; vk < mSegment[pivot].adjacencySize; ++vk)
        mWeight[mCliqueSrc[mSegment[pivot].srcPtr + vk]] =
            -mWeight[mCliqueSrc[mSegment[pivot].srcPtr + vk]];
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
    mSegment[pivot].incidenceSize = 0;
    markGone(pivot);
    return mMerged;
}

inline const std::vector<std::int32_t>& QuotientGraph::eliminateMmd(std::int32_t pivot) {
    beginEliminationMmd(pivot);
    pruneMmd(pivot);
    return finishElimination(pivot);
}

inline const std::vector<std::int32_t>& QuotientGraph::eliminateAmd(std::int32_t pivot,
                                                                    TaggedScan& scan) {
    beginEliminationAmd(pivot, scan);
    pruneAmd(pivot, scan);
    return finishElimination(pivot);
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
    const std::int32_t* reached     = mCliqueSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t reachedSize = mSegment[pivot].adjacencySize;
    for (std::uint32_t uk = 0; uk < reachedSize; ++uk) {
        const std::int32_t u = reached[uk];
        // In the amd prune the new clique goes to the FRONT of I[u] rather than the back, so the
        // single remaining entry is at the head of the incidence run either way: with A[u] empty
        // the run starts with I[u], and with one clique there is only one position. The test
        // therefore serves both branches unsuffixed.
        if (mSegment[u].adjacencySize == 0 && mSegment[u].incidenceSize == 1 &&
            mAdjIncSrc[mSegment[u].srcPtr] == pivot) {   // A[u] empty, so I[u] starts at the run
            mSegment[u].incidenceSize = 0;
            markGone(u);
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
    // The weighted clique size follows the clique, so `cliqueWeight()` stays true across the
    // merge. Magnitudes only: the sign of a weight is membership of the clique being built.
            mCliqueWeight -= weight(u);
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }

        std::int32_t*     members     = mCliqueSrc.data() + mSegment[pivot].srcPtr;
        const std::uint32_t membersSize = mSegment[pivot].adjacencySize;
        std::uint32_t       kept        = 0;
        for (std::uint32_t uk = 0; uk < membersSize; ++uk)
            if (mWeight[members[uk]] != 0) members[kept++] = members[uk];
        trimClique(pivot, kept);         // a shrink is a partial death; see numPeakCliqueMembers
    }
    return merged;
}

    // Write the survivors over the front of the clique and shorten it. What falls off the end is
    // what supervariable detection absorbed. The trimmed tail is left as a hole: this store never
    // reclaims.
inline void QuotientGraph::trimClique(std::int32_t pivot, std::uint32_t kept) {
    mNumLiveCliqueMembers -= mSegment[pivot].adjacencySize - kept;   // the trimmed tail is not live
    mSegment[pivot].adjacencySize = kept;
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
    for (std::int32_t c : mCliqueOwners) live += mSegment[c].adjacencySize;
    return live == mNumLiveCliqueMembers;
#endif
}

inline void QuotientGraph::killClique(std::int32_t c) {
    mNumLiveCliqueMembers -= mSegment[c].adjacencySize;
    mSegment[c].adjacencySize = 0;
}

inline void QuotientGraph::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    // NOT `killClique`. v is a live supervariable being absorbed, and it never formed a clique, so
    // this length is A[v]'s and not a clique's; feeding it to the counter would corrupt the peak.
    mSegment[v].adjacencySize = 0;
    mSegment[v].incidenceSize = 0;
    markGone(v);
}

inline void QuotientGraph::absorbAggressively(const std::vector<std::int32_t>& cliques,
                                              const std::int32_t*  vertices,
                                              std::uint32_t        vertexCount) {
    if (cliques.empty()) return;

    for (std::int32_t c : cliques) killClique(c);   // dead; the prune reads the size

    for (std::uint32_t uk = 0; uk < vertexCount; ++uk) {  // I[u] - dead, compacted in place
        const std::int32_t u         = vertices[uk];
        const Segment&      segment   = mSegment[u];
        std::int32_t*       incidence = mAdjIncSrc.data() + segment.srcPtr + segment.adjacencySize;
        const std::uint32_t size      = segment.incidenceSize;

    // THE ROTATION HAS TO BE REDONE WHEN ITS SUBJECT DIES. The prune parks one entry at the back
    // of I[u] and the new clique at the front. If the parked entry is one of the cliques absorbed
    // here, dropping it leaves the front entry where the rotation put it, so it is reapplied.
        const bool parkedDied = size > 0 && mSegment[incidence[size - 1]].adjacencySize == 0;

        std::uint32_t kept = 0;
        for (std::uint32_t ck = 0; ck < size; ++ck)
            if (mSegment[incidence[ck]].adjacencySize != 0) incidence[kept++] = incidence[ck];
        mSegment[u].incidenceSize = kept;

        // Entry 0 is the pivot's own new clique, which the prune put at the front and which is
        // never absorbed, so the rotation runs over positions 1 onward.
        if (parkedDied && kept > 2)
            std::rotate(incidence + 1, incidence + 2, incidence + kept);
    }
}

    // Each pivot first, then the members of its supervariable in ASCENDING VERTEX INDEX rather
    // than merge order. Same fill and same forest as `order`, the members being
    // indistinguishable; only the permutation differs.
inline std::vector<std::int32_t> QuotientGraph::orderAscending(
        const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order(mSize);
    std::vector<std::int32_t> cursor(mSize, 0);

    std::uint32_t k = 0;
    for (std::int32_t pivot : pivots) {
        order[k]      = pivot;
        cursor[pivot] = static_cast<std::int32_t>(k + 1);   // where its first member goes
        k += static_cast<std::uint32_t>(mWeight[pivot]);    // the whole supervariable's room
        // The chain starts at mSuperNext[pivot], the pivot not being a member of itself, so the
        // stamp never lands on the cursor just written.
        for (std::int32_t u = mSuperNext[pivot]; u != NIL; u = mSuperNext[u])
            cursor[u] = -(pivot + 1);
    }

    // Ascending, so the members are too. `-(x + 1)` is its own inverse, so the decode below is
    // the encode above, and its sign is the test: a non-negative cursor is a position, so u is a
    // pivot and was placed already.
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(mSize); ++u) {
        const std::int32_t pivot = -(cursor[u] + 1);
        if (pivot < 0) continue;
        order[cursor[pivot]++] = u;
    }
    return order;
}

inline std::vector<std::int32_t>
QuotientGraph::order(const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order;
    order.reserve(mSize);
    for (std::int32_t pivot : pivots)
        for (std::int32_t u = pivot; u != NIL; u = mSuperNext[u]) order.push_back(u);
    return order;
}


} // namespace Oblio
