#pragma once

// QuotientGraphFlat.h - the representation Oblio's own minimum-degree orderings run on, over a
// CLIQUE ARENA: a store that only appends, in elimination order, and never reclaims. Shared by
// every driver on this layout, mmd and amd alike.
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

#include "oblio/Buckets.h"   // TaggedScan, which four of the amd entry points take
#include "oblio/Types.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <cassert>
#include <vector>

namespace Oblio {


// The quotient graph itself: the three lists above, the liveness flags, and the supervariable
// members that mass elimination grows. A driver owns one of these, picks a pivot, calls
// the eliminator, and refreshes whatever the elimination reached.
class QuotientGraphFlat {
public:
    // Built straight from A's pattern. A is stored with both triangles, so a column's rows are
    // already its neighbors, and the only conversion is dropping the diagonal. Oblio's input
    // assumptions hold by construction, which is why nothing here symmetrizes, deduplicates or
    // sorts.
    //
    // ONE MARK PER VERTEX. Supervariable detection stamps into the driver's own tagged array, so
    // nothing asks for a second half.
    QuotientGraphFlat(const std::vector<std::size_t>&  colPtr,
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
    // one array at two levels: the eliminator stamps at level `tag` and the refresh at level
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
    // reach(u) rather than forming it. An exact degree has no use for either, uniting them in one
    // pass through `reachableSetWeight`.
    //
    // A pointer and a length rather than a container, because that is what the storage is: the two
    // lists share one run, and WHICH HALF COMES FIRST IS PER BRANCH. mmd stores A[u] first with
    // I[u] behind it; amd stores I[u] first with A[u] behind it, which is what puts the new clique
    // at index 0 of the whole run and lets hash detection walk the two halves as ONE SPAN. That is
    // the only reason the two pairs exist; the lengths are shared, a length saying nothing about
    // where its half sits.
    const std::int32_t* adjacencyMmd(std::int32_t u) const {
        return mAdjIncSrc.data() + mSegment[u].srcPtr;
    }
    const std::int32_t* adjacencyAmd(std::int32_t u) const {
        return mAdjIncSrc.data() + mSegment[u].srcPtr + mSegment[u].incidenceSize;
    }
    std::uint32_t adjacencySize(std::int32_t u) const { return mSegment[u].adjacencySize; }

    const std::int32_t* incidenceMmd(std::int32_t u) const {
        return mAdjIncSrc.data() + mSegment[u].srcPtr + mSegment[u].adjacencySize;
    }
    const std::int32_t* incidenceAmd(std::int32_t u) const {
        return mAdjIncSrc.data() + mSegment[u].srcPtr;
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
    // restore; see massEliminate. A late mass elimination is that case, restoring in
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
    // it is what the mmd prepass does with the degree-0 and degree-1 vertices. The lists are left
    // alone; every walk skips a numbered vertex from here on.
    void number(std::int32_t u);

    // SET u ASIDE, taking it out of the elimination without numbering it. This is the dense-row
    // rule's state: a variable that is neither eliminated nor available, kept out of every
    // reachable set and every list by its ZERO WEIGHT alone, which is the same mechanism a merged
    // vertex leaves by. The caller owns where it
    // lands in the permutation; the driver appends the set at the end.
    void setAside(std::int32_t u);

    // Walk I[u] from the back in the reachable-set walk. A
    // tie-break
    // convention and nothing else: it changes which permutation comes out, never which sets are
    // computed. See the member's note.
    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    // Stop the eliminator at the prune, leaving mass elimination to the caller. A late run makes
    // the same test AFTER aggressive absorption has dropped every clique lying inside the new
    // one, which is what makes the cheap structural test agree with the true one. Asking first,
    // which is what the eliminator does by default, asks it of an I[u] that still holds cliques
    // about to be removed, so the cheap test declines vertices that should merge.
    //
    // With this on, eliminateAmd returns an EMPTY merged list and C[pivot] is reach(pivot) exactly,
    // and the caller must call massEliminate() once it has absorbed. Used by AmdFlat alone. See
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

    // THE PERMUTATION, and the pair differs in one thing only. Both expand a pivot sequence into an
    // elimination order over the original vertices: a pivot stands for its whole supervariable,
    // whose members are eliminated consecutively, so this is where a supervariable of size w
    // becomes w columns. What differs is the member order WITHIN a supervariable.
    // `orderAsMerged` emits the chain as it stands, which is merge order; `orderAscending`
    // re-emits it in ascending vertex index. The members are
    // indistinguishable by construction, so the fill and the forest are the same either way and
    // only the permutation moves.
    //
    // WHICH ONE A DRIVER CALLS IS THE CALLER'S CHOICE and not a property of either function: both
    // read the same three arrays and neither can tell which branch built them. The mmd drivers
    // take `orderAscending`, which is their oracle's numbering. The amd drivers take
    // `orderAsMerged`, theirs being a raw order taken before the postorder, emitted in the order
    // each pivot accumulated its membership. So the choice records where each oracle is read
    // rather than anything about the branches.
    // mmd today
    std::vector<std::int32_t> orderAscending(const std::vector<std::int32_t>& pivots) const;
    // amd today
    std::vector<std::int32_t> orderAsMerged(const std::vector<std::int32_t>& pivots) const;

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
    std::size_t mSize;   // number of vertices; the constructor's init list is the only writer

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
    std::uint32_t             mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

    // Which end of I[u] the reachable-set walk starts from. The mmd branch expands the clique seen
    // LAST first, so with a list held in order the walk runs BACKWARD. Same set either way and the
    // same cost, but the order decides C[pivot]'s order, hence which of two equal-degree candidates
    // a later iteration finds first, and minimum degree is settled by exactly that.
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
// below is `inline`, so the units that include this header share one copy, and so is any
// variable defined here. No anonymous namespace, which would give each unit its own copy.
// ------------------------------------------------------------------------------------------------


// PEAK LIVE CLIQUE MEMBERS OF THE LAST ORDERING TO RUN, written by every driver that tracks it and
// read by tests/test_order.cpp. One symbol rather than one per driver, because the whole use is to
// compare two drivers back to back: run one, read this, run the other, read it again.
//
// ALLOCATED ON DEMAND, and the amd branch never calls this. `NIL` rather than zero as the initial
// stamp, since zero is a tag a walk can reach and this array's whole job is to answer "have I seen
// v this step" against `mTag`, which starts there.
inline void QuotientGraphFlat::enableMarks() {
    mMark.assign(mSize, NIL);
}

inline QuotientGraphFlat::QuotientGraphFlat(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mSize(colPtr.empty() ? 0 : colPtr.size() - 1), mSegment(mSize) {

    // Initialized from the matrix, diagonals excluded.
    //
    // THE RESERVE IS nnz AND THE STORE ENDS AT nnz - n, the diagonal each column drops. It is an
    // upper bound rather than the exact size because nothing here checks that every column carries
    // one; the store never grows after this loop, so the difference is slack for the object's life.
    mAdjIncSrc.reserve(colPtr.empty() ? 0 : colPtr.back());
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(mSize); ++aj) {
        // The run's start is the store's length before this column is appended.
        mSegment[aj].srcPtr = mAdjIncSrc.size();
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mAdjIncSrc.push_back(rowIdx[cp]);
        mSegment[aj].adjacencySize =
            static_cast<std::uint32_t>(mAdjIncSrc.size() - mSegment[aj].srcPtr);
    }

    // A STARTING GUESS AND NOT A BOUND. The clique store holds every C[c] ever formed, end to end,
    // and nothing here is ever reclaimed, so it passes nnz on any matrix with fill. What keeps it
    // from reallocating under a walk is the per-elimination guard in `beginElimination`, which
    // re-checks every time; this only means the guard usually finds the room already there.
    //
    // AT LEAST n, WHICH IS WHAT LETS THAT GUARD JUST DOUBLE. `capacity() >= mSize` starts here and
    // holds forever, `reserve` never shrinking, and it is what makes one doubling enough for a
    // whole reach: the guard fires only when the free room is under mSize, and then
    // `2 * capacity() - size() >= capacity() >= mSize`. Without it a first nnz below n would leave
    // doubling short of a reach, which `SparseMatrix`'s diagonal per column rules out but this
    // constructor does not check.
    mCliqueSrc.reserve(std::max(colPtr.empty() ? 0 : colPtr.back(), mSize));

    // Every vertex begins as a supervariable of one: a chain holding itself, so next is NIL and
    // last is the vertex. Mass elimination splices these together and order() walks them.
    mSuperNext.assign(mSize, NIL);
    mSuperLast.resize(mSize);
    mWeight.assign(mSize, 1);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(mSize); ++u) mSuperLast[u] = u;
}

// APPENDS to `reached` and does not clear it, which is what lets beginElimination point it straight
// at the clique arena instead of at a scratch. The returning overload below clears for itself.
inline void QuotientGraphFlat::formReachableSetMmd(std::int32_t u,
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
    const Segment&      uSegment       = mSegment[u];          // one fetch; see the member
    const std::int32_t* uAdjacency     = mAdjIncSrc.data() + uSegment.srcPtr;
    const std::uint32_t uAdjacencySize = uSegment.adjacencySize;
    const std::uint32_t uIncidenceSize = uSegment.incidenceSize;
    const std::int32_t* uIncidence     = uAdjacency + uAdjacencySize;
    // Hoisted like `live`, and for the same reason: member loads the compiler cannot prove are
    // unaliased by the stores below. Both branches are per SOURCE rather than per member, so
    // neither sits in the loop that does the work. See the members' notes for what each decides.
    const bool reverse  = mReverseIncidence;

    // Which source is walked first. The mmd branch expands the variables and then the cliques,
    // which is how the whole md ladder is laid out; the amd branch takes the cliques first and the
    // supervariables only on its last pass. Same set either way, and the order decides C[pivot]'s
    // content order, hence which of two equal-degree candidates a later iteration finds first.

    for (std::uint32_t vk = 0; vk < uAdjacencySize; ++vk) {
        const std::int32_t v = uAdjacency[vk];
        const std::int32_t vWeight = mWeight[v];
        if (vWeight > 0 && !(mHasNumbered && mMark[v] == GONE)) {
            mWeight[v] = -vWeight; reachableSet.push_back(v);
        }
    }
    for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
        const std::uint32_t t           = reverse ? uIncidenceSize - 1 - ck : ck;
        const std::int32_t  c           = uIncidence[t];
        const std::int32_t* cClique     = mCliqueSrc.data() + mSegment[c].srcPtr;
        const std::uint32_t cCliqueSize = mSegment[c].adjacencySize;
        for (std::uint32_t vk = 0; vk < cCliqueSize; ++vk) {
            const std::int32_t v = cClique[vk];
            const std::int32_t vWeight = mWeight[v];
            if (vWeight > 0) { mWeight[v] = -vWeight; reachableSet.push_back(v); }
        }
    }
}

// APPENDS to `reached` and does not clear it, which is what lets beginElimination point it straight
// at the clique arena instead of at a scratch. The returning overload below clears for itself.
inline void QuotientGraphFlat::formReachableSetAmd(std::int32_t u,
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
    const Segment&      uSegment       = mSegment[u];          // one fetch; see the member
    const std::int32_t* uIncidence     = mAdjIncSrc.data() + uSegment.srcPtr;
    const std::uint32_t uIncidenceSize = uSegment.incidenceSize;
    const std::uint32_t uAdjacencySize = uSegment.adjacencySize;
    const std::int32_t* uAdjacency     = uIncidence + uIncidenceSize;
    // Hoisted like `live`, and for the same reason: member loads the compiler cannot prove are
    // unaliased by the stores below. Both branches are per SOURCE rather than per member, so
    // neither sits in the loop that does the work. See the members' notes for what each decides.
    const bool reverse  = mReverseIncidence;

    // Which source is walked first. The mmd branch expands the variables and then the cliques,
    // which is how the whole md ladder is laid out; the amd branch takes the cliques first and the
    // supervariables only on its last pass. Same set either way, and the order decides C[pivot]'s
    // content order, hence which of two equal-degree candidates a later iteration finds first.
    for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
        const std::uint32_t t           = reverse ? uIncidenceSize - 1 - ck : ck;
        const std::int32_t  c           = uIncidence[t];
        const std::int32_t* cClique     = mCliqueSrc.data() + mSegment[c].srcPtr;
        const std::uint32_t cCliqueSize = mSegment[c].adjacencySize;
        for (std::uint32_t vk = 0; vk < cCliqueSize; ++vk) {
            const std::int32_t v = cClique[vk];
            const std::int32_t vWeight = mWeight[v];
            if (vWeight > 0) { mWeight[v] = -vWeight; reachableSet.push_back(v); }
        }
    }
    for (std::uint32_t vk = 0; vk < uAdjacencySize; ++vk) {
        const std::int32_t v = uAdjacency[vk];
        const std::int32_t vWeight = mWeight[v];
        if (vWeight > 0) { mWeight[v] = -vWeight; reachableSet.push_back(v); }
    }
}

inline std::uint32_t QuotientGraphFlat::reachableSetWeight(std::int32_t u) {
    // A sum over DISTINCT vertices, so bounded by n; see the header.
    std::uint32_t totalWeight = 0;
    ++mTag;
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique.
    const Segment&      uSegment       = mSegment[u];          // one fetch; see the member
    const std::int32_t* uAdjacency     = mAdjIncSrc.data() + uSegment.srcPtr;
    const std::uint32_t uAdjacencySize = uSegment.adjacencySize;
    const std::uint32_t uIncidenceSize = uSegment.incidenceSize;
    for (std::uint32_t vk = 0; vk < uAdjacencySize; ++vk) {
        const std::int32_t v = uAdjacency[vk];
        if (mMark[v] != GONE) {
            mMark[v] = mTag; totalWeight += static_cast<std::uint32_t>(mWeight[v]);
        }
    }
    const std::int32_t* uIncidence = uAdjacency + uAdjacencySize;
    for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
        const std::int32_t  c           = uIncidence[ck];
        const std::int32_t* cClique     = mCliqueSrc.data() + mSegment[c].srcPtr;
        const std::uint32_t cCliqueSize = mSegment[c].adjacencySize;
        for (std::uint32_t vk = 0; vk < cCliqueSize; ++vk) {
            const std::int32_t v = cClique[vk];
            if (mMark[v] < mTag) {   // includes mMark[v] != GONE, GONE sorting above every tag
                mMark[v] = mTag; totalWeight += static_cast<std::uint32_t>(mWeight[v]);
            }
        }
    }
    return totalWeight;
}

inline void QuotientGraphFlat::number(std::int32_t u) {
    // A numbered vertex lingers in every list that named it, deliberately: its neighbors keep
    // degrees that still count it. GONE is what stops the walks following it back in.
    //
    // AND THE FLAG IS WHAT TELLS THE WALKS TO ASK. See mHasNumbered: this is the only thing that
    // sets it, so a run that never calls this function never pays for the test.
    mHasNumbered = true;
    mMark[u]     = GONE;
}

inline void QuotientGraphFlat::setAside(std::int32_t u) {
    // ZERO WEIGHT IS THE WHOLE MECHANISM. `formReachableSet*` takes a vertex on `nv > 0` and the
    // prune keeps one on the same test, so a zero-weight vertex is unreachable and is dropped from
    // every list the first time that list is rewritten. GONE additionally stops `eliminated`
    // reporting it live, which is what this class asks rather than the weight.
    //
    // ITS NEIGHBORS KEEP DEGREES THAT STILL COUNT IT, exactly as after `number`, and the amd branch
    // does not correct them either: a degree is a bound and one that is too large only delays a
    // pivot.
    mWeight[u] = 0;
    markGone(u);
}

inline void QuotientGraphFlat::beginEliminationMmd(std::int32_t pivot) {
    // The reach is written STRAIGHT INTO THE ARENA. C[pivot] is the reach, so the block the walk
    // fills is the clique's own block and no scratch is needed.
    //
    // THE ARENA MUST NOT MOVE WHILE THAT WALK RUNS, which is what the reserve below is for and is
    // not an optimization. The walk reads each clique's members through a pointer into the arena
    // while appending to that same arena, so a growth would leave every such pointer dangling. A
    // reach is at most `mSize` entries, so room for one is room for the whole walk.
    //
    // ONE DOUBLING IS ALWAYS ENOUGH, on the constructor's `capacity() >= mSize`. See it there.
    if (mCliqueSrc.capacity() - mCliqueSrc.size() < mSize)
        mCliqueSrc.reserve(2 * mCliqueSrc.capacity());

    // THE DESCRIPTOR IS HELD IN LOCALS UNTIL BOTH READERS OF THE PIVOT'S RUN ARE PAST, which is
    // the whole subtlety of storing a clique in the run of the vertex that formed it. The slots it
    // will occupy are still being read: `formReachableSetMmd` walks A[pivot] and I[pivot] through
    // them, and the absorbed-clique loop below finds I[pivot] the same way. Writing either early
    // leaves a walk reading the arena through an offset into mAdjIncSrc. See the header.
    // C[p] = reach(p) (absorb into C[p])
    const std::size_t newCliquePtr = mCliqueSrc.size();
    formReachableSetMmd(pivot, mCliqueSrc);            // appends; see its note
    const std::uint32_t newCliqueSize =
        static_cast<std::uint32_t>(mCliqueSrc.size() - newCliquePtr);

    // reached(pivot) = C[pivot], one block under two readings, which is why the walk needs no
    // scratch. Taken AFTER the append, since that is what can move the arena. With the reserve
    // above it cannot have moved, and this stays as it is regardless: it costs nothing and it is
    // the shape that remains correct if the reserve is ever revised.
    const std::int32_t* newClique = mCliqueSrc.data() + newCliquePtr;

    // The absorbed cliques are I[pivot], read where they lie. Nothing below writes the pivot's run
    // (the prune rewrites the runs of C[pivot]'s members, and the pivot is not one of them), so no
    // copy and no scratch is needed to keep them alive across the passes that follow.
    const std::int32_t* pivotIncidence =
        mAdjIncSrc.data() + mSegment[pivot].srcPtr + mSegment[pivot].adjacencySize;
    const std::uint32_t pivotIncidenceSize = mSegment[pivot].incidenceSize;
    // C = C - I[p] (reclaim I[p])
    for (std::uint32_t ck = 0; ck < pivotIncidenceSize; ++ck)
        killClique(pivotIncidence[ck]);               // dead, its block left behind

    // Both readers of the pivot's own run are past, so the run becomes the clique's descriptor.
    mSegment[pivot].srcPtr        = newCliquePtr;
    mSegment[pivot].adjacencySize = newCliqueSize;

    // A CLIQUE IS BORN HERE AND NOWHERE ELSE, which is what makes the peak countable. See
    // `numPeakCliqueMembers`.
    mNumLiveCliqueMembers += newCliqueSize;
    mNumPeakCliqueMembers = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);
#ifndef NDEBUG
    mCliqueOwners.push_back(pivot);
#endif

    // Stamp the new clique. THE ABSORBED CLIQUES ARE NOT STAMPED: their members are already in
    // the reach, and stamping them would make a later walk skip entries it must still see.
    std::uint32_t newCliqueWeight = 0;
    for (std::uint32_t vk = 0; vk < newCliqueSize; ++vk)
        newCliqueWeight += static_cast<std::uint32_t>(-mWeight[newClique[vk]]);
    mCliqueWeight = newCliqueWeight;
}

inline void QuotientGraphFlat::beginEliminationAmd(std::int32_t pivot, TaggedScan& scan) {
    // The reach is written STRAIGHT INTO THE ARENA. C[pivot] is the reach, so the block the walk
    // fills is the clique's own block and no scratch is needed.
    //
    // THE ARENA MUST NOT MOVE WHILE THAT WALK RUNS, which is what the reserve below is for and is
    // not an optimization. The walk reads each clique's members through a pointer into the arena
    // while appending to that same arena, so a growth would leave every such pointer dangling. A
    // reach is at most `mSize` entries, so room for one is room for the whole walk.
    //
    // ONE DOUBLING IS ALWAYS ENOUGH, on the constructor's `capacity() >= mSize`. See it there.
    if (mCliqueSrc.capacity() - mCliqueSrc.size() < mSize)
        mCliqueSrc.reserve(2 * mCliqueSrc.capacity());

    // THE DESCRIPTOR IS HELD IN LOCALS UNTIL BOTH READERS OF THE PIVOT'S RUN ARE PAST, which is
    // the whole subtlety of storing a clique in the run of the vertex that formed it. The slots it
    // will occupy are still being read: `formReachableSetAmd` walks A[pivot] and I[pivot] through
    // them, and the absorbed-clique loop below finds I[pivot] the same way. Writing either early
    // leaves a walk reading the arena through an offset into mAdjIncSrc. See the header.
    // C[p] = reach(p) (absorb into C[p])
    const std::size_t newCliquePtr = mCliqueSrc.size();
    formReachableSetAmd(pivot, mCliqueSrc);            // appends; see its note
    const std::uint32_t newCliqueSize =
        static_cast<std::uint32_t>(mCliqueSrc.size() - newCliquePtr);

    // reached(pivot) = C[pivot], one block under two readings, which is why the walk needs no
    // scratch. Taken AFTER the append, since that is what can move the arena. With the reserve
    // above it cannot have moved, and this stays as it is regardless: it costs nothing and it is
    // the shape that remains correct if the reserve is ever revised.
    const std::int32_t* newClique = mCliqueSrc.data() + newCliquePtr;

    // The absorbed cliques are I[pivot], read where they lie. Nothing below writes the pivot's run
    // (the prune rewrites the runs of C[pivot]'s members, and the pivot is not one of them), so no
    // copy and no scratch is needed to keep them alive across the passes that follow.
    //
    // A CLIQUE DIES TWO WAYS AND THE TAGGED W MUST LEARN ABOUT BOTH. Aggressive absorption zeroes
    // `work[c]` in the driver; elimination-time absorption is this list. The amd branch writes both
    // deaths into W, `Pe[e] = FLIP(me)` with `W[e] = 0`, and its scan then tests `we != 0` off the
    // load it already needs for the value. The store rides on the walk that kills the clique, so
    // one read of the entry records both facts.
    const std::int32_t* pivotIncidence     = mAdjIncSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t pivotIncidenceSize = mSegment[pivot].incidenceSize;
    // C = C - I[p] (reclaim I[p])
    for (std::uint32_t ck = 0; ck < pivotIncidenceSize; ++ck) {
        const std::int32_t c = pivotIncidence[ck];
        killClique(c);                               // dead, its block left behind
        scan.work[c] = 0;                               // and dead to the scan; see above
    }

    // Both readers of the pivot's own run are past, so the run becomes the clique's descriptor.
    mSegment[pivot].srcPtr        = newCliquePtr;
    mSegment[pivot].adjacencySize = newCliqueSize;

    // A CLIQUE IS BORN HERE AND NOWHERE ELSE, which is what makes the peak countable. See
    // `numPeakCliqueMembers`.
    mNumLiveCliqueMembers += newCliqueSize;
    mNumPeakCliqueMembers = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);
#ifndef NDEBUG
    mCliqueOwners.push_back(pivot);
#endif

    // Stamp the new clique. THE ABSORBED CLIQUES ARE NOT STAMPED: their members are already in
    // the reach, and stamping them would make a later walk skip entries it must still see.
    std::uint32_t newCliqueWeight = 0;
    for (std::uint32_t vk = 0; vk < newCliqueSize; ++vk)
        newCliqueWeight += static_cast<std::uint32_t>(-mWeight[newClique[vk]]);
    mCliqueWeight = newCliqueWeight;
}

// THE PRUNE, mmd. For each member u of C[p], with the survivors written in their original relative
// order:
//
//     in     A[u]   v1 v2 ... vs      dropped: in C[p], the pivot itself, merged, or numbered
//            I[u]   c1 c2 ... ct      dropped: absorbed cliques
//
//     out    A[u] = [ v1 v2 ... vs ]
//            I[u] = [ c1 c2 ... ct p ]        the pivot LAST, nothing rotated
//
// The new clique is read FIRST by walking I[u] REVERSED; see setReverseIncidence.
inline void QuotientGraphFlat::pruneMmd(std::int32_t pivot) {
    const std::int32_t* newClique     = mCliqueSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t newCliqueSize = mSegment[pivot].adjacencySize;

    // ONE CURSOR FOR THE WHOLE RUN. A[u] and I[u] are compacted in place, left to right, and the
    // second compaction picks the cursor up where the first left it. The pivot is APPENDED at the
    // end.
    for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
        const std::int32_t  u              = newClique[uk];
        const Segment&      uSegment       = mSegment[u];
        std::int32_t*       uAdjacency     = mAdjIncSrc.data() + uSegment.srcPtr;
        const std::uint32_t uAdjacencySize = uSegment.adjacencySize;
        const std::int32_t* uIncidence     = uAdjacency + uAdjacencySize;
        const std::uint32_t uIncidenceSize = uSegment.incidenceSize;
        std::uint32_t       cursor         = 0;
        // A[u] = A[u] - C[p] - {p} (prune)
        for (std::uint32_t vk = 0; vk < uAdjacencySize; ++vk) {
            const std::int32_t v = uAdjacency[vk];
            if (mWeight[v] <= 0) continue;             // in the new clique, the pivot, or merged
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by a prepass
            uAdjacency[cursor++] = v;
        }
        const std::uint32_t keptAdjacencySize = cursor;   // the boundary; it does not move again
        mSegment[u].adjacencySize = keptAdjacencySize;

        // I[u] = ( I[u] - I[p] ) | {p}
        // I[u] - I[p] (reclaim I[p])
        for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
            const std::int32_t c = uIncidence[ck];
            if (mSegment[c].adjacencySize != 0) uAdjacency[cursor++] = c;  // dead is size zero
        }
        // | {p} (absorb into C[p])
        uAdjacency[cursor++] = pivot;                  // u joins the new clique, id = pivot
        mSegment[u].incidenceSize = cursor - keptAdjacencySize;
    }
    mWeight[pivot] = -mWeight[pivot];
    {
        const std::size_t   base = mSegment[pivot].srcPtr;
        const std::uint32_t len  = mSegment[pivot].adjacencySize;
        for (std::uint32_t vk = 0; vk < len; ++vk)
            mWeight[mCliqueSrc[base + vk]] = -mWeight[mCliqueSrc[base + vk]];
    }
}

// THE PRUNE, amd. The same two set operations as the plain prune, with the hash key, the bound's
// first term and the first scan riding on the walks, STORED AS [ I , A ] where mmd stores [ A, I ],
// and with BOTH PARTS ROTATED:
//
//     in     A[u]   v1 v2 ... vs      dropped: in C[p], the pivot itself, merged, or numbered
//            I[u]   c1 c2 ... ct      dropped: absorbed cliques
//
//     out    A[u] = [ v2 ... vs v1 ]           first survivor to the END
//            I[u] = [ p c2 ... ct c1 ]         the pivot FIRST, first survivor to the END
//
// THE PIVOT'S POSITION IS LOAD-BEARING: I[u] coming first, the pivot is at index 0 of the WHOLE
// RUN, so `A[u] | ( I[u] - {p} )` is ONE CONTIGUOUS SPAN from index 1 and hash detection walks it
// with a single loop. This branch walks I[u] forward rather than reversed.
//
// NEITHER ROTATION IS A CHOICE. The run is compacted forward from its own start, so the pivot has
// no free slot at the front and one spare at the end, and inserting at the front in O(1) displaces
// two entries rather than shifting the run:
//
//     c1 must leave position 0, and the only slot still inside I[u] once the boundary moves right
//     is the old first adjacency slot; that slot holds v1, so v1 goes to the spare at the end.
//
// So c1 rotating to the end of I[u] is the price of the O(1) insertion, and v1 rotating to the end
// of A[u] is what that displaces. Both orders are also the ones the amd oracle produces, which
// matters because raw order is compared against it entry for entry.
//
// UNDER THE mmd LAYOUT IT WOULD BE CHEAPER, two moves and no adjacency rotation, the boundary not
// moving when I[u] grows rightward. That is the trade this branch makes: three moves here against
// the single detection span the layout buys.
inline void QuotientGraphFlat::pruneAmd(std::int32_t pivot, TaggedScan& scan) {
    const std::int32_t* newClique     = mCliqueSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t newCliqueSize = mSegment[pivot].adjacencySize;

    // EVERY MEMBER OF C[pivot] LEAVES THE DEGREE LISTS HERE, which is what frees each member's
    // mPrev and mNext for the hash key and the hash chain. A driver that refiles inside its own
    // bound pass passes no buckets, its links still being degree links when the hash runs.
    if (scan.buckets != nullptr)
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) scan.buckets->unfile(newClique[uk]);

    const std::int32_t  workTag        = scan.workTag;          // hoisted, as the flags are

    for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
        const std::int32_t  u                = newClique[uk];
        // NEGATED, BECAUSE u IS A MEMBER OF C[pivot] AND SO READS NEGATIVE. Taking the magnitude is
        // the whole of it: `firstSeenBase` below is `workTag - uWeight` and comes out wrong by
        // twice the weight otherwise. No cast, `mWeight` being signed.
        //
        // `firstSeenBase` SEEDS A CLIQUE'S FIRST SIGHTING. The loop maintains, per clique c that
        // C[p] touches, `work[c] = workTag + |C[c]| - (weights of the members seen in c so far)`,
        // and does it in two branches that both subtract this member's weight:
        //
        //     seen already      cWork -= uWeight
        //     first sighting    cWork = degree[c] + workTag - uWeight
        //
        // so this is the constant part of the second, hoisted to once per member rather than once
        // per first-seen clique. It MUST be able to go negative, the tag being small in the early
        // eliminations while the weight is not, and that is the reason `work` is signed.
        const std::int32_t  uWeight          = -mWeight[u];
        const std::int32_t  firstSeenBase    = workTag - uWeight;   // signed, deliberately
        const Segment&      uSegment         = mSegment[u];
        std::int32_t*       source           = mAdjIncSrc.data() + uSegment.srcPtr;
        const std::uint32_t uIncidenceSize   = uSegment.incidenceSize;
        const std::uint32_t uAdjacencySize   = uSegment.adjacencySize;
        const std::int32_t* uAdjacency       = source + uIncidenceSize;
        std::uint32_t       cursor           = 0;
        std::uint32_t       uAdjacencyWeight = 0;    // |A[u] - C[p]|, weighted
    // The key is a SUM over the pruned A[u] and the final I[u], reduced modulo the driver's
    // bucket count. It WRAPS in uint32 deliberately and is reduced once, not per term. The
    // modulus must not divide any stride applied to a term, which is why no term carries one.
        std::uint32_t       uHashKey         = 0;    // wraps, deliberately
        // I[u] = ( I[u] - I[p] ) | {p}
        // I[u] - I[p] (reclaim I[p]), with the first scan riding on it
        for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
            const std::int32_t c = source[ck];
            // A clique seen earlier in this step already holds the running value above the tag; one
            // seen for the first time starts from |C[c]| and is listed once; one already absorbed
            // reads ZERO and is dropped from the list here.
            //
            // ONE LOAD, TWO QUESTIONS. `cWork == 0` is dead, which would otherwise be a
            // `mCliqueSize` probe, an array it read for nothing else, and the value is wanted
            // anyway two lines down. The other half of the change is the work-zeroing above.
            std::int32_t cWork = scan.work[c];
            if (cWork == 0) continue;                        // absorbed and gone
            source[cursor++] = c;
            if (c == pivot) continue;                     // the new clique subtracts from nothing
            if (cWork >= workTag) {
                cWork -= uWeight;
            } else {
                // A SIGNEDNESS cast: the clique degree is already 32 bits, and `firstSeenBase` is
                // negative in the early eliminations.
                cWork = static_cast<std::int32_t>(scan.degree[c]) + firstSeenBase;
                scan.touchedCliques.push_back(c);
            }
            scan.work[c] = cWork;
        }
        const std::uint32_t keptIncidenceSize = cursor;   // the boundary; it does not move again

        // A[u] = A[u] - C[p] - {p} (prune), with the key and the bound's first term riding on it
        for (std::uint32_t vk = 0; vk < uAdjacencySize; ++vk) {
            const std::int32_t v       = uAdjacency[vk];
            const std::int32_t vWeight = mWeight[v];
            if (vWeight <= 0) continue;                // see the plain prune above
            uAdjacencyWeight += static_cast<std::uint32_t>(vWeight);
            uHashKey += static_cast<std::uint32_t>(v);     // no + 1, no reduction; see above
            source[cursor++] = v;
        }
        // THE ADJACENCY HALF GOES INTO work[u], NOT INTO AN ARRAY OF ITS OWN. `work` is indexed by
        // CLIQUE id and a clique id is a dead pivot's id, so for a LIVE vertex the slot carries
        // nothing: u is in C[pivot] and therefore alive, and no clique is named after it until it
        // is eliminated, by which time this value is long consumed. It cannot collide with the
        // clique writes above: those are indexed by c drawn from I[u], every one of which is a
        // dead pivot, and u is live. The driver's obligation is one store: the slot goes back to
        // alive-and-unseen once the bound has been read.
        scan.work[u] = static_cast<std::int32_t>(uAdjacencyWeight);
        // Through the int32 slot the key rides in, read back as uint32 in the driver's bound pass.
        // The slot is the vertex's degree-list predecessor, free because every member of C[pivot]
        // was unfiled above. The key is accumulated either way, one add per surviving neighbour,
        // and STORED only where the driver asked for the bucket arrangement.
        if (scan.buckets != nullptr)
            scan.buckets->setHashKey(u, uHashKey);                             // the ADJACENCY half

        // | {p} (absorb into C[p]), BY THREE MOVES AND NOT BY AN APPEND. The run was compacted
        // forward from its own start, so the pivot has no free slot at the front and one spare at
        // the end: the first adjacency entry goes to that spare, the first incidence entry takes
        // the slot it vacated, and the pivot takes position 0. Both parts come out ROTATED, which
        // is the order the amd oracle produces.
        //
        // THE SPARE ALWAYS EXISTS. u is in C[pivot], so either u is in A[pivot] and A[u] drops the
        // pivot, or u lies in a clique of I[pivot] and I[u] drops that clique. Either way the run
        // loses at least one entry, which is what the assert states.
        assert(cursor < uIncidenceSize + uAdjacencySize &&
               "the rotation has no slot: nothing was dropped from u's run");
        source[cursor]            = source[keptIncidenceSize];
        source[keptIncidenceSize] = source[0];
        source[0]                 = pivot;
        mSegment[u].incidenceSize = keptIncidenceSize + 1;
        mSegment[u].adjacencySize = cursor - keptIncidenceSize;
    }
    mWeight[pivot] = -mWeight[pivot];
    {
        const std::size_t   base = mSegment[pivot].srcPtr;
        const std::uint32_t len  = mSegment[pivot].adjacencySize;
        for (std::uint32_t vk = 0; vk < len; ++vk)
            mWeight[mCliqueSrc[base + vk]] = -mWeight[mCliqueSrc[base + vk]];
    }
}

inline const std::vector<std::int32_t>& QuotientGraphFlat::finishElimination(std::int32_t pivot) {
    // Under mLateMassElimination the merge is the caller's, run after it has absorbed, so this
    // hands back an empty list and C[pivot] stays reach(pivot) exactly. See the setter.
    if (mLateMassElimination) {   // amd (the driver runs it later)
        mMerged.clear();
    } else {                      // mmd
        massEliminate(pivot);
    }

    // ONLY THE INCIDENCE HALF IS CLEARED. `adjacencySize` is no longer the pivot's A[pivot]: it
    // now holds |C[pivot]|, the clique this elimination just built, and zeroing it here would
    // destroy it. This line cleared both while both were dead; one of them has a second job now.
    mSegment[pivot].incidenceSize = 0;
    markGone(pivot);
    return mMerged;
}

inline const std::vector<std::int32_t>& QuotientGraphFlat::eliminateMmd(std::int32_t pivot) {
    beginEliminationMmd(pivot);
    pruneMmd(pivot);
    return finishElimination(pivot);
}

inline const std::vector<std::int32_t>& QuotientGraphFlat::eliminateAmd(std::int32_t pivot,
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
// absorption is what makes this cheap test agree with the true one.
inline const std::vector<std::int32_t>& QuotientGraphFlat::massEliminate(std::int32_t pivot) {
    mMerged.clear();   // a member scratch, kept for its capacity
    // Walks C[pivot], which is still the full reach: the trim below is this function's own and
    // happens after the loop.
    std::int32_t*       newClique     = mCliqueSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t newCliqueSize = mSegment[pivot].adjacencySize;
    for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
        const std::int32_t u = newClique[uk];
        // Two vertices are indistinguishable when
        //
        //     reach(u) | {u} == reach(v) | {v}
        //
        // and mass elimination is that with v FIXED TO THE PIVOT, which makes it a test rather
        // than a search: the counterpart is known before looking. THE TEST BELOW IS THAT FACT
        // READ AFTER p HAS GONE, so it is not the same equation in the same graph. Before the
        // elimination p is still in u's reach, and the condition is
        //
        //     reach(u) | {u} == reach(p) | {p} == C[p] | {p}
        //
        // Afterwards p is out of every reach and C[p] excludes it, so the same fact reads
        //
        //     reach(u) == C[p] - {u}
        //
        // and the three conjuncts below are the definition of reach read backwards,
        //
        //     reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}
        //
        //     adjacencySize == 0                 A[u] = {}, killing the first term
        //     incidenceSize == 1                 the union is over exactly one clique
        //     mAdjIncSrc[srcPtr] == pivot        and that clique is C[p]
        //
        // Free, because all three are descriptor fields the prune has just written for this
        // vertex. Nothing is computed, hashed or compared against another vertex.
        //
        // In the amd prune the new clique goes to the FRONT of I[u] rather than the back, so the
        // single remaining entry is at the head of the incidence run either way: with A[u] empty
        // the run starts with I[u], and with one clique there is only one position. The test
        // therefore serves both branches unsuffixed.
        if (mSegment[u].adjacencySize == 0 && mSegment[u].incidenceSize == 1 &&
            mAdjIncSrc[mSegment[u].srcPtr] == pivot) {   // A[u] empty, so I[u] starts at the run
            mSegment[u].incidenceSize = 0;
            markGone(u);
            mMerged.push_back(u);
        }
    }
    if (!mMerged.empty()) {
        // THE MERGE HAPPENS FIRST, so the compaction can read the ZERO WEIGHT it leaves rather
        // than a stamp of its own. The old order was the reverse and needed a tag pass over
        // `mMerged` plus a mark read per member; the weight says the same thing and the
        // supervariable bookkeeping had to write it anyway.
        //
        // NO OTHER MEMBER OF C[pivot] CAN READ ZERO, which is what makes the test exact: a vertex a
        // live merge folded away is left at weight zero but is also stamped GONE, so no reach ever
        // emits it into a clique again.
        for (std::int32_t u : mMerged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            // The weighted clique size follows the clique, so `cliqueWeight()` stays true across
            // the merge.
            mCliqueWeight -= static_cast<std::uint32_t>(mWeight[u]);
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }

        std::uint32_t cursor = 0;
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk)
            if (mWeight[newClique[uk]] != 0) newClique[cursor++] = newClique[uk];
        trimClique(pivot, cursor);       // a shrink is a partial death; see numPeakCliqueMembers
    }
    return mMerged;
}

    // Write the survivors over the front of the clique and shorten it. What falls off the end is
    // what supervariable detection absorbed. The trimmed tail is left as a hole: this store never
    // reclaims.
inline void QuotientGraphFlat::trimClique(std::int32_t pivot, std::uint32_t kept) {
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
inline bool QuotientGraphFlat::cliqueCountBalances() const {
#ifdef NDEBUG
    return true;
#else
    std::size_t live = 0;
    for (std::int32_t c : mCliqueOwners) live += mSegment[c].adjacencySize;
    return live == mNumLiveCliqueMembers;
#endif
}

inline void QuotientGraphFlat::killClique(std::int32_t c) {
    mNumLiveCliqueMembers -= mSegment[c].adjacencySize;
    mSegment[c].adjacencySize = 0;
}

inline void QuotientGraphFlat::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u]                += mWeight[v];
    mWeight[v]                = 0;

    // v is a live supervariable being absorbed and never formed a clique, so this length is
    // A[v]'s and not a clique's.
    mSegment[v].adjacencySize = 0;
    mSegment[v].incidenceSize = 0;
    markGone(v);          // mmd only; the amd branch reads the zero weight above and has no array
}

inline void QuotientGraphFlat::absorbAggressively(const std::vector<std::int32_t>& cliques,
                                              const std::int32_t*  vertices,
                                              std::uint32_t        vertexCount) {
    if (cliques.empty()) return;

    for (std::int32_t c : cliques) killClique(c);   // dead; the prune reads the size

    for (std::uint32_t uk = 0; uk < vertexCount; ++uk) {  // I[u] - dead, compacted in place
        const std::int32_t u         = vertices[uk];
        const Segment&      segment   = mSegment[u];
        std::int32_t*       incidence = mAdjIncSrc.data() + segment.srcPtr;
        const std::uint32_t size      = segment.incidenceSize;

    // THE ROTATION HAS TO BE REDONE WHEN ITS SUBJECT DIES. The prune parks one entry at the back
    // of I[u] and the new clique at the front. If the parked entry is one of the cliques absorbed
    // here, dropping it leaves the front entry where the rotation put it, so it is reapplied.
        const bool parkedDied = size > 0 && mSegment[incidence[size - 1]].adjacencySize == 0;

        // COMPACTED RIGHTWARD, so the dropped entries fall off the FRONT and `srcPtr` advances.
        // A[u] sits BEHIND I[u] here and must not move; compacting leftward would open a hole
        // between the two and cost a slide of the adjacency part. The run's END does not move,
        // one entry leaving the front for every entry dropped.
        std::uint32_t write = size;
        for (std::uint32_t ck = size; ck-- > 0;)
            if (mSegment[incidence[ck]].adjacencySize != 0) incidence[--write] = incidence[ck];
        const std::uint32_t kept = size - write;
        mSegment[u].srcPtr += write;
        mSegment[u].incidenceSize = kept;

        // Entry 0 is the pivot's own new clique, which the prune put at the front and which is
        // never absorbed, so the rotation runs over positions 1 onward.
        if (parkedDied && kept > 2)
            std::rotate(incidence + write + 1, incidence + write + 2, incidence + size);
    }
}

inline std::vector<std::int32_t> QuotientGraphFlat::orderAscending(
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
QuotientGraphFlat::orderAsMerged(const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order;
    order.reserve(mSize);
    for (std::int32_t pivot : pivots)
        for (std::int32_t u = pivot; u != NIL; u = mSuperNext[u]) order.push_back(u);
    return order;
}


} // namespace Oblio
