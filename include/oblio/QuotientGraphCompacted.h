#pragma once

// QuotientGraphCompacted.h - the quotient graph on a POOLED CLIQUE LAYOUT: one workspace holding a
// vertex's two lists AND every live clique, with a free cursor and a compaction, where
// `QuotientGraph` keeps C[c] in a separate store that only grows.
//
// SAME GRAPH, SAME ALGORITHMS, DIFFERENT STORAGE. Every idea in QuotientGraph.h holds here without
// change. Read that file first; this one comments only what the storage makes different, which is
// where the cliques live and what has to happen when the pool fills.
//
// WHAT THE POOL BUYS. It is sized ONCE from the pattern and never grows: if A fits, the ordering
// completes. A store that grows toward nnz(L) cannot promise that, nnz(L) depending on the
// ordering being computed.
//
// WHAT IT COSTS is a compactor, and that is what makes this class harder than the flat one.
// Absorbed cliques and consumed list prefixes leave dead space, the cursor reaches the end, and
// everything live slides down. **So every walk here is written in POSITIONS AND CURSORS rather
// than pointers and counters**: a compaction moves blocks within a pool that never reallocates, so
// a position survives it and a pointer does not, and a truncation shortens the list being read, so
// a count of what is LEFT survives and a count of what has been done does not.
//
// ONE CLASS, TWO BRANCHES, SPLIT BY SUFFIX. Where the two branches agree there is one method;
// where they disagree there are two, named `...Amd` and `...Mmd`. The count of splits is not a
// constraint: contorting shared code to avoid a suffix is worse than the suffix.
//
// THE THREE THINGS THAT GENUINELY DIFFER, and they are order or lifetime rather than storage:
//
//   - WHERE THE NEW CLIQUE LANDS IN I[u], first or last. That decides the order members enter
//     C[pivot], hence bucket order, hence which of several equal-degree vertices is picked, so
//     neither is free to move. `pruneAmd` and `pruneMmd`.
//   - WHICH HALF OF THE RUN COMES FIRST, and which way the incidence list is walked. The layout is
//     these accessors and the walks and NOTHING ELSE: the pool, the cursor, the elbow room, the
//     compactor's sweep, the descriptors and the counters are all blind to it.
//     `adjacencyAmd`/`adjacencyMmd` with their incidence twins, and the `reachableSet` pair.
//   - HOW A DEAD VERTEX IS RECOGNIZED. Amd reads a zero weight and costs no array at all. Mmd
//     cannot: `number()` leaves a prepass vertex LIVE at weight one and in every list that names
//     it, so only a tag can hide it. `eliminatedAmd` and `eliminatedMmd`, and `mMark` exists for
//     the second alone.
//
// THE MARK ARRAY IS ALLOCATED ON DEMAND, by `enableMarks`, so the amd driver pays a pointer rather
// than n int32. **`mHasNumbered` IS THEREFORE LOAD BEARING** and not merely an optimization: it is
// the short circuit that keeps an empty array safe in the shared bodies carrying the
// `mHasNumbered && mMark[v] == GONE` guard.

#include "oblio/QuotientGraph.h"   // Buckets and TaggedScan, which are shared verbatim
#include "oblio/Types.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace Oblio {

class QuotientGraphCompacted {
public:
    // NO cliqueMarks ARGUMENT, unlike QuotientGraph's. That flag sized a mark array at n or 2n;
    // here the array is absent unless a driver asks for it, and no driver wants clique marks.
    QuotientGraphCompacted(const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const { return mSize; }

    // HOW MANY TIMES THE POOL WAS COMPACTED. This is the class's
    // whole storage figure, and there is deliberately no second one: the pool is sized at
    // construction and never grows, so what the cliques cost is not a quantity that varies. The
    // flat class reports one because its clique store does.
    std::size_t numCompactions() const { return mCompactions; }

    // A BLOCK HEAD IN THE POOL, during compaction only: FLIPPED - e for owner e, so every
    // head is at most FLIPPED and every real entry, a vertex id, is above it.
    //
    // AND THE BOTTOM OF THE RANGE IS EXACTLY INT32_MIN, which is worth stating because nothing
    // about the value says so. An owner runs to n - 1 and n is capped at MAX_IDX, so the encoded
    // heads occupy the whole of `[INT32_MIN, FLIPPED]` and the last one lands on the type's
    // minimum with nothing to spare. It fits at every admissible n and would not at one more.
    // `-2` rather than `-1` is what leaves NIL out of the head range.
    static constexpr std::int32_t FLIPPED = -2;

    // GONE: one value above every reachable tag, so the stamp array
    // answers "is v dead" on the load it was making anyway. Read by the mmd branch alone.
    static constexpr std::int32_t GONE = std::numeric_limits<std::int32_t>::max();

    // ------------------------------------------------------------------ per-branch: liveness

    // AND A TAG ON THE MMD BRANCH, because a zero weight is not available to it: `number()` leaves
    // a prepass vertex live at weight one so that its neighbors' degrees still count it.
    bool eliminatedMmd(std::int32_t u) const { return mMark[u] == GONE; }

    // ZERO WEIGHT IS THE DEAD STATE on the amd branch. The
    // three ways a vertex leaves the graph there all end in a zero weight, so no array is spent.
    bool eliminatedAmd(std::int32_t u) const { return mWeight[u] == 0; }

    // ------------------------------------------------------------------ per-branch: the layout

    // ADJACENCY FIRST, INCIDENCE BEHIND IT, which is the mmd order and the flat class's. The
    // compactor needs no flip here: mmd walks A[u] forward and I[u] backward, so the consumed part
    // is a prefix in one phase and a suffix in the other and the remainder is contiguous in both.
    const std::int32_t* adjacencyMmd(std::int32_t u) const {
        return mSrc.data() + mSegment[u].srcPtr;
    }
    const std::int32_t* incidenceMmd(std::int32_t u) const {
        return mSrc.data() + mSegment[u].srcPtr + mSegment[u].adjacencySize;
    }

    // INCIDENCE FIRST, ADJACENCY BEHIND IT, which is the amd order. The flip is what lets the
    // new clique go in by its three-move rotation and what makes the part a cliques-first walk has
    // consumed a PREFIX, which is what the mid-walk compactor needs.
    const std::int32_t* incidenceAmd(std::int32_t u) const {
        return mSrc.data() + mSegment[u].srcPtr;
    }
    const std::int32_t* adjacencyAmd(std::int32_t u) const {
        return mSrc.data() + mSegment[u].srcPtr + mSegment[u].incidenceSize;
    }

    // THE LENGTHS ARE NOT SPLIT. A length is a count and says nothing about where the half sits.
    std::uint32_t adjacencySize(std::int32_t u) const { return mSegment[u].adjacencySize; }
    std::uint32_t incidenceSize(std::int32_t u) const { return mSegment[u].incidenceSize; }

    // ------------------------------------------------------------------ cliques

    // A clique has no incidence part, so its block starts at the run under either layout and this
    // needs no suffix. Spelled without the offset because a clique is one list rather than two.
    const std::int32_t* clique(std::int32_t c) const {
        return mSrc.data() + mSegment[c].srcPtr;
    }
    std::uint32_t cliqueSize(std::int32_t c) const { return mSegment[c].adjacencySize; }

    // WRITABLE, for the pass that trims the clique as it walks it. See `trimClique`.
    std::int32_t* clique(std::int32_t c) { return mSrc.data() + mSegment[c].srcPtr; }
    void trimClique(std::int32_t pivot, std::uint32_t kept);

    // THE ONE PLACE A CLIQUE DIES, so that the counter sees every death. Absorbed into the new
    // clique, absorbed aggressively once its external degree reaches zero, or merged away with its
    // owner: all three go through here rather than storing the zero where they stand.
    void killClique(std::int32_t c);

    std::uint32_t cliqueWeight() const { return mCliqueWeight; }

    // ------------------------------------------------------------------ weights

    // The MAGNITUDE, for a driver, which does not want to know about the sign. The hot loops
    // inside this class read mWeight directly and test the sign, which is the whole point of it.
    std::uint32_t weight(std::int32_t u) const {
        const std::int32_t w = mWeight[u];
        return static_cast<std::uint32_t>(w < 0 ? -w : w);
    }

    // A ROW THE DENSE RULE SET ASIDE, amd only. Its weight goes to zero and the row is then
    // absent from every reachable set and every prune, which a zero weight already achieves here.
    void setAside(std::int32_t u) { mWeight[u] = 0; }

    // ------------------------------------------------------------------ the mark array, mmd only

    // ALLOCATED ON DEMAND. The amd branch never calls this and so carries an empty vector rather
    // than n int32. Call once, before any elimination.
    void enableMarks();

    std::int32_t advanceTag()                        { return ++mTag; }
    std::int32_t mark(std::int32_t u) const          { return mMark[u]; }
    void setMark(std::int32_t u, std::int32_t tag)   { mMark[u] = tag; }

    // A VERTEX NUMBERED BY THE MMD PREPASS, which is not eliminated in the quotient-graph sense:
    // it keeps its weight and its place in every list that names it, and only the tag hides it.
    void number(std::int32_t u);

    // ------------------------------------------------------------------ elimination

    void beginEliminationMmd(std::int32_t pivot);
    // THE THREE STEPS OF AN ELIMINATION, and the pair of wrappers over them. A driver calls one
    // wrapper; the steps stay public because the wrapper is a sequence rather than a replacement,
    // and its value is that the ORDER lives here where nothing in a driver could enforce it.
    //
    // TWO WRAPPERS AND NOT ONE, which is why they carry the branch in the name: amd's prune
    // computes the degree bound and the hash key in the same pass and mmd's does neither, so one
    // signature would have to carry a parameter one branch ignores.
    //
    // THE FIRST STEP SPLITS AND IT IS THE COSTLIEST OF THE SUFFIXES, because the placement rule it
    // holds is genuinely common: build in place when the incidence list is empty, otherwise at the
    // free cursor, compact when the cursor runs out. What cannot be shared is the two lines that
    // NAME A WALK, and a walk cannot be selected without knowing the branch. So the tail and the
    // capture are factored into private helpers and the six lines of the rule read twice, each
    // against its own reference.
    void beginEliminationAmd(std::int32_t pivot, TaggedScan& scan);

    void pruneMmd(std::int32_t pivot);
    void pruneAmd(std::int32_t pivot, TaggedScan& scan);

    // ONE METHOD FOR BOTH BRANCHES. `markGone` inside it does nothing when the mark array is
    // absent, which is the amd branch, so the two halves this was split into differed in a store
    // that was already a no-op on one of them.
    const std::vector<std::int32_t>& finishElimination(std::int32_t pivot);

    // AND THE THREE AS ONE CALL, overloaded on the scan exactly as `QuotientGraph`'s pair is: the
    // mmd prune takes nothing and the amd one takes a `TaggedScan`, so the argument selects the
    // branch. The value of the wrapper is that the ORDER of the three steps lives here rather
    // than in each driver, where nothing could enforce it. The three remain public and remain
    // suffixed; this is a sequence rather than a replacement.
    const std::vector<std::int32_t>& eliminateMmd(std::int32_t pivot);
    const std::vector<std::int32_t>& eliminateAmd(std::int32_t pivot, TaggedScan& scan);

    // MASS ELIMINATION, one method for both branches. It reads no sign: the merge test is
    // structural, and the prune restored every weight before this runs.
    const std::vector<std::int32_t>& massEliminate(std::int32_t pivot);

    // AGGRESSIVE ABSORPTION, amd only: a clique whose external degree has reached zero is dead,
    // and its members drop it from their incidence lists.
    void absorbAggressively(const std::vector<std::int32_t>& cliques,
                            const std::int32_t* vertices, std::uint32_t vertexCount);

    // ONE SUPERVARIABLE ABSORBS ANOTHER. `v` is a live vertex that never formed a clique, so its
    // length is A[v]'s and NOT a clique's, which is why this does not call `killClique`.
    void merge(std::int32_t u, std::int32_t v);

    // Whether mass elimination is the driver's rather than this class's. Amd sets it: its prune
    // needs the untrimmed clique for the degree bound, so the trim happens later and in the driver.
    void setLateMassElimination(bool on) { mLateMassElimination = on; }

    // Whether I[u] is walked backward, popping the incidence like a stack. Set by the mmd
    // driver; the amd walk has its own method and does not read it.
    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    // ------------------------------------------------------------------ the reachable set

    std::uint32_t formReachableSetMmd(std::int32_t u);
    // THE WALK THAT BUILDS THE NEW CLIQUE, one per branch, differing in the order of the two
    // halves and the direction of the incidence part. Both write at the free cursor, both may
    // COMPACT MID-WALK, and both return the number of members reached, the start having moved.
    std::uint32_t formReachableSetAmd(std::int32_t u);

    std::uint32_t formReachableSetInPlaceMmd(std::int32_t u);
    // AND THE IN-PLACE FORM: a pivot with no cliques has a reach
    // that is a SUBSET of A[pivot], so it is compacted where it stands and the pool is untouched.
    // The mmd form additionally rejects a prepass-numbered vertex, which a weight cannot exclude.
    std::uint32_t formReachableSetInPlaceAmd(std::int32_t u);

    // THE DEGREE OF A REACHABLE SET WITHOUT BUILDING IT, mmd only, called per vertex in the
    // refresh. It stamps rather than negating, so it leaves the weights untouched.
    std::uint32_t reachableSetWeight(std::int32_t u);

    // ------------------------------------------------------------------ counters and output

    // PEAK LIVE CLIQUE MEMBERS. A MEMBER is a vertex in a live clique at this instant, where an
    // ENTRY is a pool slot; the two differ here and not in the flat class, because this layout
    // reclaims. IT IS A PROPERTY OF THE ALGORITHM AND NOT OF THE LAYOUT, which is what makes it
    // comparable against `AmdFlat`'s and `MmdFlat`'s: two drivers agreeing on a permutation can
    // still be caught doing different work, and have been. Checked in tests/test_order.cpp.
    std::size_t numPeakCliqueMembers() const { return mNumPeakCliqueMembers; }
    std::size_t numLiveCliqueMembers() const { return mNumLiveCliqueMembers; }

    // THE COUNTER CHECKED AGAINST A RECOMPUTATION, debug builds only. Five sites move the count: a
    // birth in `beginElimination`, a death in `killClique`, the shrink in mass elimination, the
    // trim in `trimClique`, and the mid-walk truncation inside a walk. A funnel is a claim about
    // the CALL GRAPH and holds only while every writer goes through it, which twice it did not;
    // the way to make such a claim checkable is a recomputation from independent state, and this
    // is it.
    //
    // COMPILED OUT under NDEBUG, where the assertion that calls it is compiled out too.
#ifndef NDEBUG
    bool cliqueCountBalances() const;
#endif

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
    std::vector<std::int32_t> orderAscending(const std::vector<std::int32_t>& pivots) const;
    // amd today
    std::vector<std::int32_t> orderAsMerged(const std::vector<std::int32_t>& pivots) const;

private:
    // COMPACTION, and it takes the half-built clique's start BY
    // REFERENCE. A compaction during a build finds the partial at the top of the pool with no
    // descriptor of its own, so the sweep stops just below it and the partial is then moved down
    // as a block and its start rewritten. Called with `cliqueStart == mFree` outside a build,
    // where the partial is empty and both halves degenerate to the plain sweep.
    void compact(std::size_t& cliqueStart);

    // THE TWO HALVES OF `beginElimination` THAT ARE COMMON. `captureAbsorbed` copies I[pivot]
    // before a walk truncates it; `bearClique` kills those cliques, repoints the pivot's run at
    // the new clique, records the birth and sums the clique's weight. Between them sits the
    // placement rule and the walk, which is the only part that knows the branch.
    void captureAbsorbed(const std::int32_t* incidence, std::uint32_t count);
    void bearClique(std::int32_t pivot, std::size_t cliqueStart, std::uint32_t cliqueLen);

    // GONE IF THE ARRAY EXISTS, a no-op if it does not, which is what lets the three bodies that
    // retire a vertex be shared across the branches. The amd branch never enables marks.
    void markGone(std::int32_t u) { if (!mMark.empty()) mMark[u] = GONE; }

    // Dimensions.
    std::size_t mSize;   // number of vertices; the constructor's init list is the only writer

    // THE POOL. Every vertex's segment and every live clique, end to end, sized once at
    // construction from the off-diagonal count and never grown.
    std::vector<std::int32_t> mSrc;

    // TWO LIVES PER FIELD, DECIDED BY WHETHER THIS ID HAS FORMED A CLIQUE, and both live in the
    // one pool here, where the flat class keeps them in two arrays. The turn happens in exactly
    // one place, `bearClique`. Every other way a vertex leaves the graph, a hash merge, mass
    // elimination or the dense-row rule, forms no clique and leaves the fields meaning what they
    // meant.
    struct Segment {
        // WHERE THE BLOCK STARTS, AND IT MOVES. Three things move it: `compact` sliding every
        // live block down, a walk truncating the pivot's segment to what it has not consumed, and
        // `bearClique` repointing the pivot at its new clique. Anything holding a POSITION
        // survives all three by re-reading this field; a pointer does not.
        std::size_t   srcPtr;
        // A[u]'S LENGTH for a vertex, THE MEMBER COUNT for a clique. ZERO MEANS A DIFFERENT THING
        // IN EACH, and both readings are live: on a vertex it is an empty adjacency, which is what
        // `massEliminate`'s first conjunct tests, and on a clique it is DEAD, written by
        // `killClique` and read by the prune's incidence compaction. Do not fold the two.
        std::uint32_t adjacencySize;
        // I[u]'S LENGTH. Zeroed when the clique is born and never read again, a clique having no
        // incidence part.
        std::uint32_t incidenceSize;
    };
    std::vector<Segment> mSegment;

    std::size_t mFree = 0;             // the free cursor
    std::size_t mCompactions = 0;      // how often the pool was compacted

    // Whether the clique now being eliminated was built in the pivot's own run. Read by mass
    // elimination, which can only give space back to the cursor in the other case.
    bool mBuiltInPlace = false;

    // I[pivot], CAPTURED BEFORE THE WALK CONSUMES IT. A walk truncates that list as it reads it,
    // so after a mid-walk compaction the run is short or empty and the cliques it named could
    // never be killed from it: each would keep a non-zero length, stay a live block for the
    // compactor to copy, and never be subtracted from the live count.
    std::vector<std::int32_t> mAbsorbed;

    std::size_t mNumLiveCliqueMembers = 0;   // see numPeakCliqueMembers
    std::size_t mNumPeakCliqueMembers = 0;
#ifndef NDEBUG
    // Every vertex that ever formed a clique, which is what the recomputation sums over. It cannot
    // be taken from a driver's pivot list: those also carry vertices that never formed one, the
    // dense rows on one branch and the prepass-numbered on the other, and for those
    // `adjacencySize` is still A[p]'s length.
    std::vector<std::int32_t> mCliqueOwners;
#endif

    // The supervariable a vertex stands for, as a chain rather than a list per vertex.
    std::vector<std::int32_t> mSuperNext;
    std::vector<std::int32_t> mSuperLast;

    // THE WEIGHT, CARRYING AN ENCODING IN ITS SIGN, which is why it is signed:
    //
    //     mWeight[v] >  0    live, and not yet taken into the clique being built; the weight
    //     mWeight[v] <  0    live, and taken into it this step; the weight is -mWeight[v]
    //     mWeight[v] == 0    dead, by a hash merge or by mass elimination, or set aside as dense
    //
    // ONE LOAD ANSWERS THREE QUESTIONS in the two hottest loops here. The ZERO is a true sentinel
    // and holds on the amd branch alone, where the only death sites are `merge` and mass
    // elimination; on the mmd branch `number()` leaves a prepass vertex live at weight one, which
    // is why that branch keeps `mMark`. Four conditions have to hold for a signed
    // one dimensional size has to meet.
    std::vector<std::int32_t> mWeight;

    std::uint32_t mCliqueWeight = 0;   // per pivot, not per vertex; see cliqueWeight()

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

    bool mLateMassElimination = false;
    bool mReverseIncidence    = false;

    // MMD ONLY, AND ABSENT UNLESS `enableMarks` IS CALLED. `mHasNumbered` is load bearing rather
    // than an optimization: it is the short circuit that keeps an empty `mMark` safe in the shared
    // bodies carrying the `mHasNumbered && mMark[v] == GONE` guard.
    bool                      mHasNumbered = false;
    std::vector<std::int32_t> mMark;
    std::int32_t              mTag = 0;
};

// ------------------------------------------------------------------------------------------------
// THE BODIES ARE HERE, NOT IN A SOURCE FILE, AND THAT IS DELIBERATE. A class compiled into its own
// translation unit is an opaque call from a driver's pivot loop: the pool's and the run array's
// bases are reloaded around every call and the loop's registers are spilled. In the driver's own
// unit the compiler inlines it and keeps them. Every ordering driver is arranged this way, so that
// they can be compared with EACH OTHER; it says nothing about a comparison against an oracle,
// which differs from us in ways nobody has enumerated. Every out-of-class definition below
// is `inline` so that several drivers may include this and the linker folds the copies.
//
// No anonymous namespace, which would give each unit its own copy.
// ------------------------------------------------------------------------------------------------


// A NAMED NAMESPACE RATHER THAN AN ANONYMOUS ONE, because this file is included by every driver
// that uses the class. An anonymous namespace would give each of them its own copy of what follows,
// and an inline member calling into a copy that differs per unit is an ODR violation waiting for a
// linker that folds differently. `inline` on both is what makes one definition serve all of them.
namespace detail {

#ifdef OBLIO_PAD_ORDERING
// The pool's own step; see `padded` below for why the amounts are staggered.
inline constexpr std::size_t kOrderingPad = 5 * (16384 / sizeof(std::int32_t));
#else
inline constexpr std::size_t kOrderingPad = 0;
#endif

// STAGGERED PADDING FOR THE DATA-PLACEMENT EXPERIMENT, off unless `OBLIO_PAD_ORDERING` is defined,
// and it changes nothing that is computed: each size-n vector is given extra CAPACITY it never
// uses, so every `size()` and every subscript is exactly as it was and only the addresses move.
//
// WHY STAGGERED RATHER THAN A UNIFORM PAGE. The failure mode being tested for is same-sized arrays
// landing in the same cache sets, and a large allocation is page aligned and rounded up to whole
// pages, so adding the SAME amount to each leaves them the same size and therefore still
// congruent. Different amounts is what perturbs the arrangement, and the step must be a whole page
// or the allocator returns the same addresses and the knob does nothing.
template <class T>
inline void padded(std::vector<T>& v, std::size_t size, std::size_t which) {
#ifdef OBLIO_PAD_ORDERING
    constexpr std::size_t page = 16384 / sizeof(T);   // one page of whatever this holds
    v.reserve(size + (which + 1) * page);
#else
    (void)which;
#endif
    v.resize(size);
}

} // namespace detail

inline QuotientGraphCompacted::QuotientGraphCompacted(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mSize(colPtr.empty() ? 0 : colPtr.size() - 1) {
    detail::padded(mSegment, mSize, 0);

    // ELBOW ROOM. The pool holds `offDiagonalCount + offDiagonalCount/5 + n` entries. The figure
    // has to match the amd oracle's exactly rather than merely be of the same shape: the headroom
    // decides how often the compactor runs, so a differential on compaction counts measures the
    // headroom unless the two agree. It is computed from the OFF-DIAGONAL count, an earlier version
    // having counted WITH the diagonal and so run about 1.2n large.
    //
    // Reserved from an upper bound and then sized down, rather than counted in a pass of its own:
    // `nnz` is at least the off-diagonal count, so the reserve cannot reallocate, and the exact
    // figure is known once the runs are laid out. The vector is SIZED rather than left at capacity
    // because the space past the runs is written by cliques and read by nothing until it is, and a
    // cursor into a vector's unused capacity would be undefined.
    const std::size_t nnz = colPtr.empty() ? 0 : colPtr.back();
    mSrc.reserve(nnz + nnz / 5 + mSize + detail::kOrderingPad);
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(mSize); ++aj) {
        mSegment[aj].srcPtr = mSrc.size();
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSrc.push_back(rowIdx[cp]);
        mSegment[aj].adjacencySize = static_cast<std::uint32_t>(mSrc.size() - mSegment[aj].srcPtr);
    }

    // The runs are laid out; everything past them is free space the cliques will use.
    const std::size_t offDiagonalCount = mSrc.size();
    mFree = offDiagonalCount;
    mSrc.resize(offDiagonalCount + offDiagonalCount / 5 + mSize);

    detail::padded(mSuperNext, mSize, 1);
    detail::padded(mSuperLast, mSize, 2);
    detail::padded(mWeight,    mSize, 3);
    std::fill(mSuperNext.begin(), mSuperNext.end(), NIL);
    std::fill(mWeight.begin(),    mWeight.end(),    1);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(mSize); ++u) mSuperLast[u] = u;
}

// WRITES AT `mFree` RATHER THAN APPENDING, which is the storage change, and TESTS FOR ROOM PER
// ENTRY, so a compaction can fire inside this walk and the walk then resumes; the positions and
// cursors below are what make resuming possible. The reach lands exactly where the clique is to
// live, so there is no copy from a scratch into place.
inline std::uint32_t QuotientGraphCompacted::formReachableSetMmd(std::int32_t u) {
    ++mTag;
    // THE PIVOT IS ALREADY NEGATED, by `beginEliminationMmd`, which is where the amd branch has
    // always done it and where the merge moved this branch's. It was here and in the in-place walk
    // before, and leaving it in both places negates twice, which reads as a positive weight and so
    // as a vertex not yet taken. Caught on a 3 by 3 grid by `cliqueCountBalances`, in the first
    // build of the merged class.
    // THE TRUNCATION BELOW IS WRITTEN FOR THE REVERSED WALK, which is the only one this file runs.
    // Reversed, the consumed part of I[u] is a suffix and the survivors are the prefix the run
    // already points at, so a truncation is a length; forward it would be the other way round and
    // `srcPtr` would have to move with it. The flag is the driver's, so the assumption is
    // checked rather than assumed.
    const bool reverse = mReverseIncidence;
    assert(reverse && "formReachableSetMmd's truncation assumes the reversed incidence walk");

    // POSITIONS AND CURSORS, and this is what makes resuming possible at all. The compactor moves
    // every live block, so a pointer taken before it runs points at whatever landed there
    // afterwards; a POSITION survives, the pool being sized once at construction so that a
    // compaction invalidates a block's OFFSET and never the array's base. And a cursor counting
    // what is LEFT survives a truncation, where a counter into the original length would not.
    // A pointer would not survive one, which is the whole reason for keeping positions and
    // restoring them from the descriptors the compaction has just rewritten.
    std::int32_t* const pool = mSrc.data();             // stable: the pool never reallocates
    std::uint32_t       reached     = 0;                // counted; the start moves under us
    std::size_t         cliqueStart = mFree;            // the clique's start, and it moves too

    // A[u] FORWARD, so the consumed part is a PREFIX and what remains is one contiguous block
    // running to the end of the run. `p` is where the unconsumed part starts and `remaining` how
    // much of it is left, which is the pair a truncation writes back.
    std::size_t   p         = mSegment[u].srcPtr;
    std::uint32_t remaining = mSegment[u].adjacencySize;
    while (remaining-- > 0) {
        const std::int32_t v  = pool[p++];
        const std::int32_t vWeight = mWeight[v];
        if (vWeight <= 0 || (mHasNumbered && mMark[v] == GONE)) continue;

        if (mFree >= mSrc.size()) {
            // TRUNCATE TO WHAT IS LEFT OF A[u]: the descriptor moves up and the length drops by
            // knt1`. The consumed part is a prefix, so the survivors are still one block: A[u]'s
            // tail followed by the whole of I[u], which is what `srcPtr` and the two lengths
            // now describe. It is also what makes the space bound hold, the half-built clique
            // being a subset of exactly the prefixes dropped here.
            mSegment[u].srcPtr     = p;
            mSegment[u].adjacencySize = remaining;
            compact(cliqueStart);
            p = mSegment[u].srcPtr;                         // re-read after a compaction
        }

        mWeight[v] = -vWeight;
        pool[mFree++] = v;
        ++reached;
    }

    // I[u] REVERSED, the incidence being pushed and popped like a stack. The
    // consumed part is therefore a SUFFIX, so what remains is the PREFIX the run already points at
    // and a truncation is a length rather than a move. `left` is that length.
    //
    // THE BASE IS RE-READ RATHER THAN HOISTED, for the same reason the walk uses positions: a
    // compaction moves the run. `srcPtr + adjacencySize` is where I[u] begins under every state
    // this loop can be in, including after the first truncation below has set the length to zero.
    std::size_t   incidenceBase = mSegment[u].srcPtr + mSegment[u].adjacencySize;
    const std::uint32_t incidenceLen = mSegment[u].incidenceSize;
    std::uint32_t left = incidenceLen;
    while (left-- > 0) {
        const std::int32_t c  = pool[incidenceBase + (reverse ? left : incidenceLen - 1 - left)];
        std::size_t        pj = mSegment[c].srcPtr;
        std::uint32_t      ln = mSegment[c].adjacencySize;
        while (ln-- > 0) {
            const std::int32_t v  = pool[pj++];
            const std::int32_t vWeight = mWeight[v];
            if (vWeight <= 0) continue;

            if (mFree >= mSrc.size()) {
                // A[u] IS SPENT BY NOW, so the run becomes I[u]'s surviving prefix alone and the
                // consumed adjacency goes with the rest of the dead space. Dropping it is not
                // tidiness: the space bound rests on the half-built clique being covered by what
                // the walk has already released.
                mSegment[u].srcPtr     = incidenceBase;
                mSegment[u].adjacencySize = 0;
                mSegment[u].incidenceSize = left;

                // TRUNCATING A CLIQUE IS A CONTRACTION AND THE COUNTER HAS TO SEE IT. `c` is a
                // live clique being consumed into the new one; dropping its consumed prefix
                // shortens it, and the `killClique` in `beginElimination` would then subtract only
                // what remained, leaving the difference live for the rest of the run. This is the
                // second of the two defects `AmdCompacted` introduced with its own compactor.
                mNumLiveCliqueMembers -= mSegment[c].adjacencySize - ln;
                mSegment[c].srcPtr      = pj;
                mSegment[c].adjacencySize  = ln;

                compact(cliqueStart);
                // AND RESUME FROM THE DESCRIPTORS, both of which the compaction has just rewritten.
                incidenceBase = mSegment[u].srcPtr + mSegment[u].adjacencySize;
                pj            = mSegment[c].srcPtr;
            }

            mWeight[v] = -vWeight;
            pool[mFree++] = v;
            ++reached;
        }
    }
    return reached;
}

// WRITES AT `mFree` RATHER THAN APPENDING, which is the storage change, and TESTS FOR ROOM PER
// ENTRY, the free cursor being tested against the pool's end. The walk resumes
// afterwards, which is what the two steps below it were for.
//
// INCIDENCE THEN ADJACENCY, one loop over the two halves:
// the cliques of me on the first `elenme` passes and the supervariables on the last. It is now
// also the physical order of the run, so the walk reads the two parts in the order they lie.
//
// POSITIONS RATHER THAN POINTERS, and this is what makes resuming possible at all. The compactor
// moves every live block, so a pointer taken before it runs points at whatever landed there
// afterwards. A POSITION survives, because the pool itself never reallocates: `mSrc` is sized
// once at construction and the compactor only slides data within it, so what a compaction
// invalidates is a block's OFFSET and never the array's base. Hoisting that base costs the same as
// holding a pointer, which is the whole reason for keeping positions rather than pointers as
// indices and restoring them from the descriptors the compaction has just rewritten.
//
// CURSORS RATHER THAN COUNTERS, for the same reason again. A compaction TRUNCATES the two lists
// being read, dropping the part already consumed, so after it a list starts at a new base with a
// new length and a counter into the old one means nothing. `p` walks the pivot's run and `pj` the
// clique being read, and both are re-read from their descriptors afterwards.
inline std::uint32_t QuotientGraphCompacted::formReachableSetAmd(std::int32_t u) {
    std::int32_t* const pool = mSrc.data();               // stable: the pool never reallocates
    std::uint32_t reached      = 0;                       // counted; the start moves under us
    std::size_t   cliqueStart  = mFree;                   // the clique's start, and it moves too

    // The pivot's own run, walked incidence part first. `remaining` is what is left of the part
    // being consumed, so a truncation is `srcPtr = p` with the length set from it.
    std::size_t   p         = mSegment[u].srcPtr;
    std::uint32_t remaining = mSegment[u].incidenceSize;
    bool          onCliques = true;

    for (;;) {
        std::size_t   pj;                                 // where the members being read start
        std::uint32_t ln;                                 // how many are left of them
        std::int32_t  c = NIL;                            // the clique, or NIL on the last pass
        if (onCliques) {
            if (remaining == 0) {                         // the incidence part is spent
                onCliques = false;
                remaining = mSegment[u].adjacencySize;
                continue;
            }
            c  = pool[p++];
            --remaining;
            pj = mSegment[c].srcPtr;
            ln = mSegment[c].adjacencySize;
        } else {
            if (remaining == 0) break;                    // and the adjacency part with it
            pj = p;                                       // the supervariables ARE the run's tail
            ln = remaining;
            remaining = 0;
        }

        // ONE LOAD PER MEMBER, and the store that takes a vertex into the clique is the same
        // slot. `> 0` is live-and-not-yet-taken, which is dedup and liveness together.
        while (ln-- > 0) {
            const std::int32_t v  = pool[pj++];
            const std::int32_t vWeight = mWeight[v];
            if (vWeight <= 0) continue;

            if (mFree >= mSrc.size()) {
                // TRUNCATE BOTH LISTS TO WHAT IS LEFT, moving each descriptor forward and
                // dropping its length by what was consumed. This is the step the flipped run
                // exists for: the consumed part is a PREFIX, so what
                // remains is still one contiguous block and the descriptor can simply be moved
                // forward. It is also what makes the space bound hold, the half-built clique
                // being a subset of exactly the prefixes dropped here.
                if (onCliques) {
                    mSegment[u].srcPtr     = p;
                    mSegment[u].incidenceSize = remaining;
                } else {
                    mSegment[u].srcPtr     = pj;
                    mSegment[u].incidenceSize = 0;
                    mSegment[u].adjacencySize = ln;
                }
                // TRUNCATING A CLIQUE IS A CONTRACTION AND THE COUNTER HAS TO SEE IT. `c` is a
                // live clique being consumed into the new one; dropping its consumed prefix
                // shortens it, and `killClique` below will then subtract only what remained. Not
                // telling the counter here left the difference live forever, which is what the
                // `AmdCompacted pC differs` check caught on `JGD_Trefethen/Trefethen_2000`, +31 of
                // 11091. The peak is unaffected either way, being taken only at a birth, and by
                // then the clique is dead in full: the two subtractions sum to its original length.
                //
                // The pivot's own run needs no such line. `u` is a live vertex and the length
                // being shortened is A[u]'s, which this counter never held.
                if (c != NIL) {
                    mNumLiveCliqueMembers -= mSegment[c].adjacencySize - ln;
                    mSegment[c].srcPtr      = pj;
                    mSegment[c].adjacencySize  = ln;
                }
                compact(cliqueStart);
                // AND RESUME FROM THE DESCRIPTORS, both of which the compaction has just rewritten.
                p  = mSegment[u].srcPtr;
                pj = (c != NIL) ? mSegment[c].srcPtr : p;
                if (!onCliques) p = pj;
            }

            mWeight[v] = -vWeight;
            pool[mFree++] = v;
            ++reached;
        }
        if (!onCliques) p = pj;
    }
    return reached;
}

// CONSTRUCTS THE CLIQUE IN THE PIVOT'S OWN RUN, which is the no-clique branch and
// the common case rather than a corner: a pivot with no cliques has a reach that is a SUBSET of
// A[pivot], so it fits where A[pivot] already is and the pool is not touched.
//
// IT IS AN IN-PLACE COMPACTION and safe for the reason every such loop here is: the write cursor
// starts at the read cursor and only ever falls behind it, a vertex being written only when it was
// read.
//
// THE GONE TEST SURVIVES HERE, unlike in AmdCompacted, and for the reason it survives everywhere on
// this branch: `number()` leaves a prepass vertex at weight one and in every neighbour's adjacency,
// so a positive weight does not mean live in an ADJACENCY list. This walk reads nothing else.
inline std::uint32_t QuotientGraphCompacted::formReachableSetInPlaceMmd(std::int32_t u) {
    // THE PIVOT IS ALREADY NEGATED, by `beginEliminationMmd`. It has to be negated somewhere and
    // not merely to keep u out of its own list, which has no self loop: the prune restores the
    // sign of every member AND of the pivot, so a pivot left positive comes out of that restore
    // negative and `orderAscending` then reads a supervariable of negative size. Negating twice,
    // which is what hoisting it without removing it here would do, is the same bug reversed.
    std::int32_t*       source = mSrc.data() + mSegment[u].srcPtr;
    const std::uint32_t count = mSegment[u].adjacencySize;
    std::uint32_t       kept  = 0;
    for (std::uint32_t vk = 0; vk < count; ++vk) {
        const std::int32_t v  = source[vk];
        const std::int32_t vWeight = mWeight[v];
        if (vWeight > 0 && !(mHasNumbered && mMark[v] == GONE)) {
            mWeight[v] = -vWeight; source[kept++] = v;
        }
    }
    return kept;
}

// CONSTRUCTS THE CLIQUE IN THE PIVOT'S OWN RUN, which is the no-clique branch and
// the common case rather than a corner: a pivot with no cliques has a reach that is a SUBSET of
// A[pivot], so it fits where A[pivot] already is and the pool is not touched at all.
//
// IT IS AN IN-PLACE COMPACTION and safe for the reason every such loop here is: the write cursor
// starts at the read cursor and only ever falls behind it, since a vertex is written only when it
// was read.
//
// WHY IT MATTERS TWICE. It keeps two thirds of cliques out of the pool, so the cursor advances far
// more slowly and the compactor runs far less; and it leaves the clique exactly where the pivot's
// adjacency was, which is where the vertices that will read it next are. This file existed for a
// year without it, and the figures it produced were the price of a DIFFERENT layout.
inline std::uint32_t QuotientGraphCompacted::formReachableSetInPlaceAmd(std::int32_t u) {
    std::int32_t*       source = mSrc.data() + mSegment[u].srcPtr;
    const std::uint32_t count = mSegment[u].adjacencySize;
    std::uint32_t       kept  = 0;
    for (std::uint32_t vk = 0; vk < count; ++vk) {
        const std::int32_t v  = source[vk];
        const std::int32_t vWeight = mWeight[v];
        if (vWeight > 0) { mWeight[v] = -vWeight; source[kept++] = v; }
    }
    return kept;
}

inline std::uint32_t QuotientGraphCompacted::reachableSetWeight(std::int32_t u) {
    // A sum over DISTINCT vertices, so bounded by n; see the header.
    std::uint32_t totalWeight = 0;
    ++mTag;
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique.
    const Segment&      uSegment       = mSegment[u];          // one fetch; see the member
    const std::int32_t* uAdjacency     = mSrc.data() + uSegment.srcPtr;
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
        const std::int32_t* cClique     = mSrc.data() + mSegment[c].srcPtr;
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

// THE SAME RECLAIM, ONE PASS LATER AND FOR THE OTHER REASON. `massEliminate` above drops the
// members it merged into the pivot; this drops the members supervariable detection absorbed into
// each other, which happens after it. Doing both in one place needs detection inside the scan, so
// that a single restoring pass sees every casualty at once and writes the survivors back. Ours
// needs two passes because detection is a pass of its own.
//
// WITHOUT THIS THE ABSORBED STAY IN THE CLIQUE FOR THE REST OF THE RUN. No permutation moves, every
// later walk skipping them on the `nv > 0` test, but they are visited, and the space behind them
// is never given back, so the pool fills faster and the compactor runs sooner.
inline void QuotientGraphCompacted::trimClique(std::int32_t pivot, std::uint32_t kept) {
    const std::uint32_t was = mSegment[pivot].adjacencySize;
    if (kept == was) return;
    mNumLiveCliqueMembers -= was - kept;
    mSegment[pivot].adjacencySize = kept;
    if (!mBuiltInPlace) {
        assert(mSegment[pivot].srcPtr + was == mFree &&
               "the pivot's clique is no longer the last block in the pool");
        mFree = mSegment[pivot].srcPtr + kept;
    }
}

inline void QuotientGraphCompacted::killClique(std::int32_t c) {
    mNumLiveCliqueMembers -= mSegment[c].adjacencySize;
    mSegment[c].adjacencySize = 0;
}

// ALLOCATED ON DEMAND, and the amd branch never calls this. `NIL` rather than zero as the initial
// stamp, since zero is a tag a walk can reach and this array's whole job is to answer "have I seen
// v this step" against `mTag`, which starts there.
inline void QuotientGraphCompacted::enableMarks() {
    detail::padded(mMark, mSize, 4);
    std::fill(mMark.begin(), mMark.end(), NIL);
}

inline void QuotientGraphCompacted::number(std::int32_t u) {
    mHasNumbered = true;
    mMark[u]     = GONE;
}

// The same rule on the mmd branch. The walk order is the difference; see the header.
inline void QuotientGraphCompacted::beginEliminationMmd(std::int32_t pivot) {
    captureAbsorbed(incidenceMmd(pivot), mSegment[pivot].incidenceSize);

    mWeight[pivot] = -mWeight[pivot];

    const bool inPlace = mSegment[pivot].incidenceSize == 0;
    mBuiltInPlace = inPlace;

    std::size_t   cliqueStart;
    std::uint32_t cliqueLen;
    if (inPlace) {
        cliqueStart = mSegment[pivot].srcPtr;
        cliqueLen   = formReachableSetInPlaceMmd(pivot);
    } else {
        cliqueLen   = formReachableSetMmd(pivot);
        cliqueStart = mFree - cliqueLen;
    }
    bearClique(pivot, cliqueStart, cliqueLen);
}

// THE PLACEMENT RULE, WRITTEN TWICE. Everything about building a clique is common except the two
// lines that NAME A WALK, and a walk cannot be selected without knowing the branch. So the capture
// before it and the bookkeeping after it are shared helpers and only the rule itself reads twice.
//
// TWO WAYS TO BUILD, and which applies is whether I[pivot] is empty. With no cliques in it the
// reach is a SUBSET of A[pivot] and is compacted where it stands, so the pool is not touched at
// all; otherwise it is assembled at the free cursor. Most eliminations take the first path, which
// is why the branch is worth having and not merely faithful.
//
// NO RESERVATION BEFORE THE WALK. The walks test per entry and compact from inside themselves,
// the free cursor being tested against the pool's end. Reserving room for a worst-case
// reach of n beforehand, which is what a walk holding pointers had to do, compacted where the
// would not.
inline void QuotientGraphCompacted::beginEliminationAmd(std::int32_t pivot, TaggedScan& scan) {
    captureAbsorbed(incidenceAmd(pivot), mSegment[pivot].incidenceSize);

    // The pivot's own weight is negated, and it is what keeps it out of its own clique: the
    // walks take a vertex only when its weight reads positive.
    mWeight[pivot] = -mWeight[pivot];

    const bool inPlace = mSegment[pivot].incidenceSize == 0;
    mBuiltInPlace = inPlace;

    std::size_t   cliqueStart;
    std::uint32_t cliqueLen;
    if (inPlace) {
        cliqueStart = mSegment[pivot].srcPtr;
        cliqueLen   = formReachableSetInPlaceAmd(pivot);
    } else {
        // AFTER the walk, not before it. The cursor is where the clique ends, so its start is a
        // subtraction, and a compaction inside the walk moves both together and leaves this right.
        cliqueLen   = formReachableSetAmd(pivot);
        cliqueStart = mFree - cliqueLen;
    }
    bearClique(pivot, cliqueStart, cliqueLen);

    // AND DEAD TO THE SCAN, from the copy `captureAbsorbed` took. A second loop rather than one,
    // `bearClique` being shared with the mmd branch, which has no scan to write into. The moment is
    // the same one the plain class kills at: after the walk that read their member lists.
    for (std::int32_t c : mAbsorbed) scan.work[c] = 0;

    // A CLIQUE HAS NO INCIDENCE LIST. The mmd branch clears this in `finishEliminationMmd`
    // instead, its mass elimination reading the pivot's run in between.
    mSegment[pivot].incidenceSize = 0;
}

// THE PRUNE, mmd's. One walk of each vertex of C[pivot] drops the dead from A[u], keeps the
// live cliques of I[u], and APPENDS the new clique at the back, which is the mmd convention and
// the reverse of amd's. No degree bound and no hash key: this branch refreshes degrees in a
// pass of its own.
inline void QuotientGraphCompacted::pruneMmd(std::int32_t pivot) {
    const std::int32_t* newClique     = mSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t newCliqueSize = mSegment[pivot].adjacencySize;

    for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
        const std::int32_t  u              = newClique[uk];
        const Segment&      uSegment       = mSegment[u];
        std::int32_t*       uAdjacency     = mSrc.data() + uSegment.srcPtr;
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
        // THE NEW CLIQUE GOES AT THE BACK, which is the mmd convention and the reverse of
        // amd's. It is what the reversed incidence walk then reads first, and it decides the
        // order members enter C[pivot], hence bucket order, hence which of several equal-degree
        // vertices is picked. The held-vertex rotation that put it at the FRONT for amd went
        // | {p} (absorb into C[p])
        uAdjacency[cursor++] = pivot;                  // u joins the new clique, id = pivot
        mSegment[u].incidenceSize = cursor - keptAdjacencySize;
    }

    // THE SIGNS COME BACK HERE, at the end of the prune, which is the LAST READER of them. The
    // walk that built C[pivot] marked membership by negating a weight and the loop above is what
    // consumed that mark. Nothing downstream asks: absorption never touches a weight, mass
    // elimination's merge test is structural, and every other reader takes a magnitude through
    // `weight()`.
    mWeight[pivot] = -mWeight[pivot];
    {
        const std::size_t   base = mSegment[pivot].srcPtr;
        const std::uint32_t len  = mSegment[pivot].adjacencySize;
        for (std::uint32_t vk = 0; vk < len; ++vk)
            mWeight[mSrc[base + vk]] = -mWeight[mSrc[base + vk]];
    }
}

// THE PRUNE, amd's, with the first scan fused in. One walk of each vertex of C[pivot]
// rewrites both halves of its run, accumulates the degree bound's explicit part, builds the
// adjacency half of the hash key, and leaves the incidence half in the tagged `work`.
inline void QuotientGraphCompacted::pruneAmd(std::int32_t pivot, TaggedScan& scan) {
    const std::int32_t* newClique     = mSrc.data() + mSegment[pivot].srcPtr;
    const std::uint32_t newCliqueSize = mSegment[pivot].adjacencySize;

    if (scan.buckets != nullptr)
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) scan.buckets->unfile(newClique[uk]);

    const std::int32_t  workTag        = scan.workTag;          // hoisted, as the flag is

    for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk) {
        const std::int32_t  u                = newClique[uk];
        const std::int32_t  uWeight          = -mWeight[u];
        const std::int32_t  firstSeenBase    = workTag - uWeight;   // signed, deliberately
        const Segment&      uSegment         = mSegment[u];
        std::int32_t*       source           = mSrc.data() + uSegment.srcPtr;
        const std::uint32_t uIncidenceSize   = uSegment.incidenceSize;
        const std::uint32_t uAdjacencySize   = uSegment.adjacencySize;
        const std::int32_t* uAdjacency       = source + uIncidenceSize;
        std::uint32_t       cursor           = 0;
        std::uint32_t       uAdjacencyWeight = 0;    // |A[u] - C[p]|, weighted
        std::uint32_t       uHashKey         = 0;    // wraps, deliberately

        // THE INCIDENCE PART, AT THE FRONT OF THE RUN. Compacted where it lies: the cursor cursor
        // starts at the read cursor and advances only when an entry is kept, so it never passes it.
        for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
            const std::int32_t c = source[ck];
            std::int32_t cWork = scan.work[c];
            if (cWork == 0) continue;                        // absorbed and gone
            source[cursor++] = c;
            if (c == pivot) continue;                     // the new clique subtracts from nothing
            if (cWork >= workTag) {
                cWork -= uWeight;
            } else {
                cWork = static_cast<std::int32_t>(scan.degree[c]) + firstSeenBase;
                scan.touchedCliques.push_back(c);
            }
            scan.work[c] = cWork;
        }
        const std::uint32_t keptIncidenceSize = cursor;        // the boundary

        // THE ADJACENCY PART, BEHIND IT, and the cursor cursor carries on through the boundary. It
        // is still behind the read cursor, having entered this loop at most `uIncidenceSize`.
        for (std::uint32_t vk = 0; vk < uAdjacencySize; ++vk) {
            // ONE LOAD, THREE QUESTIONS, the weight answering all of them. A
            // negative weight is a member of the new clique, including the pivot itself, so the
            // explicit `v == pivot` test is gone with the rest; a zero is absorbed or merged away.
            const std::int32_t v       = uAdjacency[vk];
            const std::int32_t vWeight = mWeight[v];
            if (vWeight <= 0) continue;                // see the plain prune above
            uAdjacencyWeight += static_cast<std::uint32_t>(vWeight);
            uHashKey += static_cast<std::uint32_t>(v);     // no + 1, no reduction; see above
            source[cursor++] = v;
        }
        scan.work[u] = static_cast<std::int32_t>(uAdjacencyWeight);
        if (scan.buckets != nullptr)
            scan.buckets->setHashKey(u, uHashKey);                             // the ADJACENCY half

        // THE NEW CLIQUE GOES IN BY ROTATION, WHICH IS THREE MOVES:
        //
        //     the first adjacency entry moves to the spare slot at the end
        //     the first incidence entry moves into the slot that vacated, the boundary
        //     the pivot takes position 0
        //
        // There is no free slot to insert into, so a rotation is what puts the pivot at the front
        // without a shift, and it is why each part reads with its first entry last. The old code
        // reproduced the same two orders under the reversed layout by holding the first uAdjacency
        // survivor back and swapping the pivot into place afterwards; with the parts in this
        // order the rotation says it directly and takes a test out of the uAdjacency loop.
        //
        // ROOM IS GUARANTEED BY ONE DROPPED ENTRY, and there is always at least one: u is in
        // C[pivot], so either the pivot was in A[u] and the walk above dropped it, or a clique of
        // I[u] was absorbed into the new one and was dropped as `cWork == 0`. Both self-assign
        // harmlessly when a part comes out empty.
        assert(cursor < uIncidenceSize + uAdjacencySize &&
               "the rotation has no slot: nothing was dropped from u's run");
        source[cursor]            = source[keptIncidenceSize];
        source[keptIncidenceSize] = source[0];
        source[0]                 = pivot;
        mSegment[u].incidenceSize = keptIncidenceSize + 1;
        mSegment[u].adjacencySize = cursor - keptIncidenceSize;
    }

    // THE SIGNS COME BACK HERE, at the end of the prune, which is the LAST READER of them. The
    // walk that built C[pivot] marked membership by negating a weight and the loop above is what
    // consumed that mark. Nothing downstream asks: absorption never touches a weight, mass
    // elimination's merge test is structural, and every other reader takes a magnitude through
    // `weight()`.
    mWeight[pivot] = -mWeight[pivot];
    {
        const std::size_t   base = mSegment[pivot].srcPtr;
        const std::uint32_t len  = mSegment[pivot].adjacencySize;
        for (std::uint32_t vk = 0; vk < len; ++vk)
            mWeight[mSrc[base + vk]] = -mWeight[mSrc[base + vk]];
    }
}

// MASS ELIMINATION IS THE DRIVER'S ON THIS BRANCH, so this normally only clears the scratch: the
// amd prune needs the untrimmed clique for its degree bound, so the trim happens later and in the
// driver. The flag is still tested rather than assumed, the class not being the place that decides.
inline const std::vector<std::int32_t>&
QuotientGraphCompacted::finishElimination(std::int32_t pivot) {
    if (mLateMassElimination) {   // amd (the driver runs it later)
        mMerged.clear();
    } else {                      // mmd
        massEliminate(pivot);
    }

    mSegment[pivot].incidenceSize = 0;
    markGone(pivot);
    return mMerged;
}

inline const std::vector<std::int32_t>& QuotientGraphCompacted::eliminateMmd(std::int32_t pivot) {
    beginEliminationMmd(pivot);
    pruneMmd(pivot);
    return finishElimination(pivot);
}

inline const std::vector<std::int32_t>& QuotientGraphCompacted::eliminateAmd(std::int32_t pivot,
                                                                          TaggedScan& scan) {
    beginEliminationAmd(pivot, scan);
    pruneAmd(pivot, scan);
    return finishElimination(pivot);
}

// MASS ELIMINATION, which eliminates non-principal supervariables. A vertex whose reach is
// exactly the new clique is indistinguishable from the pivot and is absorbed into it.
//
// THE AMD FORM DOES NOT RESTORE THE NEGATED WEIGHTS, its driver doing that in a pass it already
// makes over C[pivot]. The mmd form must, `reachableSetWeight` being called per vertex in the
// refresh.
inline const std::vector<std::int32_t>&
QuotientGraphCompacted::massEliminate(std::int32_t pivot) {
    mMerged.clear();   // a member scratch, kept for its capacity
    std::int32_t*       newClique     = mSrc.data() + mSegment[pivot].srcPtr;
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
        //     mSrc[srcPtr] == pivot              and that clique is C[p]
        //
        // Free, because all three are descriptor fields the prune has just written for this
        // vertex. Nothing is computed, hashed or compared against another vertex.
        //
        // In the amd prune the new clique goes to the FRONT of I[u] rather than the back, so the
        // single remaining entry is at the head of the incidence run either way: with A[u] empty
        // the run starts with I[u], and with one clique there is only one position. The test
        // therefore serves both branches unsuffixed.
        if (mSegment[u].adjacencySize == 0 && mSegment[u].incidenceSize == 1 &&
            mSrc[mSegment[u].srcPtr] == pivot) {         // A[u] empty, so I[u] starts at the run
            mSegment[u].incidenceSize = 0;
            markGone(u);
            mMerged.push_back(u);
        }
    }
    if (!mMerged.empty()) {
        for (std::int32_t u : mMerged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            mCliqueWeight -= static_cast<std::uint32_t>(mWeight[u]);
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }

        std::uint32_t cursor = 0;
        for (std::uint32_t uk = 0; uk < newCliqueSize; ++uk)
            if (mWeight[newClique[uk]] != 0) newClique[cursor++] = newClique[uk];
        trimClique(pivot, cursor);   // a shrink is a partial death, and in the pooled
                                     // case it hands the space back to mFree
    }
    return mMerged;
}

inline void QuotientGraphCompacted::absorbAggressively(
    const std::vector<std::int32_t>& cliques,
    const std::int32_t* vertices, std::uint32_t vertexCount) {
    if (cliques.empty()) return;

    // dead; the compaction reads the size
    for (std::int32_t c : cliques) killClique(c);

    // I[u] LESS THE DEAD, COMPACTED RIGHTWARD SO THE ADJACENCY NEVER MOVES. The incidence part is
    // at the front of the run and the adjacency starts at `srcPtr + incidenceSize`, so
    // compacting to the LEFT would drag A[u] down behind it, a copy per surviving vertex per
    // absorbing step. Compacting to the RIGHT instead leaves the boundary where it is: survivors
    // end flush against A[u], `srcPtr` advances by however many were dropped, and the freed
    // space falls off the FRONT of the run rather than the back. The adjacency is not touched at
    // all, and I[u] is the short part after a prune, holding the new clique plus whatever survived.
    //
    // THE PASS ITSELF IS STILL A DIVERGENCE. The reference absorbs inside its second scan, in the
    // walk that rewrites the list, so an absorbed clique is simply not copied. We decide absorption
    // after the whole prune, the prune having the first scan fused into it, and so have to revisit
    // the list. What is removed here is the slide, not the visit.
    for (std::uint32_t uk = 0; uk < vertexCount; ++uk) {
        const std::int32_t u         = vertices[uk];
        std::int32_t*      incidence = mSrc.data() + mSegment[u].srcPtr;
        const std::uint32_t size     = mSegment[u].incidenceSize;

        // THE ROTATION HAS TO BE REDONE WHEN ITS SUBJECT DIES. The reference rotates inside that
        // scan, with this step's absorbed cliques already dropped, so the entry it parks at the
        // back is the first SURVIVOR; our prune rotates first and absorbs here, so when the parked
        // entry is one this pass removes the two lists differ. The same correction has to exist in
        // both files or their columns stop being comparable.
        const bool parkedDied = size > 0 && mSegment[incidence[size - 1]].adjacencySize == 0;

        std::uint32_t write = size;
        for (std::uint32_t i = size; i-- > 0;)
            if (mSegment[incidence[i]].adjacencySize != 0) incidence[--write] = incidence[i];
        mSegment[u].srcPtr += write;                     // the dropped entries fall off the front
        mSegment[u].incidenceSize = size - write;

        // Entry 0 of what survives is the pivot's own new clique, never absorbed, so the rotation
        // runs over the entries behind it.
        const std::uint32_t kept = size - write;
        if (parkedDied && kept > 2)
            std::rotate(incidence + write + 1, incidence + write + 2, incidence + size);
    }
}

inline void QuotientGraphCompacted::merge(std::int32_t u, std::int32_t v) {
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

#ifndef NDEBUG
inline bool QuotientGraphCompacted::cliqueCountBalances() const {
    std::size_t live = 0;
    for (std::int32_t c : mCliqueOwners) live += mSegment[c].adjacencySize;
    return live == mNumLiveCliqueMembers;
}
#endif

inline std::vector<std::int32_t> QuotientGraphCompacted::orderAscending(
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

inline std::vector<std::int32_t> QuotientGraphCompacted::orderAsMerged(
        const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order;
    order.reserve(mSize);
    for (std::int32_t pivot : pivots)
        for (std::int32_t u = pivot; u != NIL; u = mSuperNext[u]) order.push_back(u);
    return order;
}

// THE COMPACTION. The pool fills because absorbed cliques and shortened runs
// leave dead space behind, and nothing reclaims it until this does. Everything live is compacted to
// the front and every descriptor is repointed.
//
// THE TRICK IS WHAT MAKES THIS ONE PASS RATHER THAN A SORT. A block's owner
// cannot be recovered by looking at the pool, so before compacting, the FIRST ENTRY of every live
// block is replaced by FLIP(owner) and the entry it displaced is parked in that owner's own
// `srcPtr`, which is about to be overwritten anyway. The pool can then be scanned from the
// front: a negative entry is a block head, its owner is FLIP of it, and its length is the owner's
// descriptor.
//
// A LIVE BLOCK IS ONE WITH A NONZERO LENGTH, and that covers both kinds without a test for which:
// a live vertex has I[e] then A[e], a dead pivot has C[e], and anything absorbed or merged had its
// lengths zeroed when it died.
//
// AND IT CARRIES THE HALF-BUILT CLIQUE, which is why it takes its start by reference. When the
// compactor runs during a clique build the partial sits at the top of the pool with no descriptor
// of its own, so the sweep stops just below it and the partial is
// then moved down explicitly and its start rewritten. Called with `cliqueStart == mFree` outside a
// build, where the partial is empty and both halves degenerate to the plain sweep.
inline void QuotientGraphCompacted::compact(std::size_t& cliqueStart) {
    ++mCompactions;
    const std::int32_t n = static_cast<std::int32_t>(mSize);

    for (std::int32_t e = 0; e < n; ++e) {
        const std::size_t len = static_cast<std::size_t>(mSegment[e].adjacencySize)
                              + mSegment[e].incidenceSize;
        if (len == 0) continue;
        const std::size_t at = mSegment[e].srcPtr;
        mSegment[e].srcPtr = static_cast<std::size_t>(
                                static_cast<std::uint32_t>(mSrc[at]));      // park the displaced
        mSrc[at] = FLIPPED - e;                                             // and mark the head
    }

    std::size_t src = 0, dst = 0;
    while (src < cliqueStart) {                             // up to the half-built clique
        const std::int32_t head = mSrc[src++];
        if (head > FLIPPED) continue;                       // dead space, skipped
        const std::int32_t e = FLIPPED - head;
        mSrc[dst] = static_cast<std::int32_t>(
                           static_cast<std::uint32_t>(mSegment[e].srcPtr));     // restore it
        mSegment[e].srcPtr = dst;
        ++dst;
        const std::size_t rest = static_cast<std::size_t>(mSegment[e].adjacencySize)
                               + mSegment[e].incidenceSize - 1;
        for (std::size_t k = 0; k < rest; ++k) mSrc[dst + k] = mSrc[src + k];
        src += rest;
        dst += rest;
    }

    // THE PARTIAL CLIQUE FOLLOWS, moved as a block. It has no descriptor to repoint, so its start
    // is the caller's variable and is returned through it.
    const std::size_t moved = dst;
    for (std::size_t k = cliqueStart; k < mFree; ++k) mSrc[dst++] = mSrc[k];
    cliqueStart = moved;
    mFree = dst;
}

// I[pivot] COPIED BEFORE ANY WALK TOUCHES IT. A walk truncates that list as it consumes it, so
// after a mid-walk compaction the run is short or empty and the cliques it named could never be
// killed from it: each would keep a non-zero length, stay a live block for the compactor to copy,
// and never be subtracted from the live count. That defect existed in `AmdCompacted` for a day,
// moved no permutation, and was found by a benchmark comparing peaks across drivers.
inline void QuotientGraphCompacted::captureAbsorbed(const std::int32_t* incidence,
                                                    std::uint32_t count) {
    mAbsorbed.clear();
    for (std::uint32_t ck = 0; ck < count; ++ck) mAbsorbed.push_back(incidence[ck]);
}

// THE ABSORBED CLIQUES DIE HERE: after the walk has read them and BEFORE the new clique is born.
// Both halves of that placement are forced. After the walk, because it needs their member lists.
// Before the birth, because `numPeakCliqueMembers` is a running maximum and the absorbed cliques
// and the one absorbing them are never live at the same instant; killing afterwards would raise
// the peak by their combined size on every step that absorbs anything.
inline void QuotientGraphCompacted::bearClique(std::int32_t pivot, std::size_t cliqueStart,
                                        std::uint32_t cliqueLen) {
    for (std::int32_t c : mAbsorbed) killClique(c);   // dead, their blocks left behind

    // Both readers of the pivot's own run are past, so the run becomes the clique's descriptor.
    mSegment[pivot].srcPtr     = cliqueStart;
    mSegment[pivot].adjacencySize = cliqueLen;

    // A CLIQUE IS BORN HERE AND NOWHERE ELSE, which is what makes the peak countable, and the
    // maximum is taken here alone since nothing else raises the total.
    mNumLiveCliqueMembers += cliqueLen;
    mNumPeakCliqueMembers  = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);
#ifndef NDEBUG
    mCliqueOwners.push_back(pivot);
#endif

    // NO STAMPING PASS AND NO TAG. Membership was written by the walk, in the sign of the weight,
    // so this only sums.
    const std::int32_t* reached = mSrc.data() + cliqueStart;
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t vk = 0; vk < cliqueLen; ++vk)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[vk]]);
    mCliqueWeight = cliqueWeight;
}


} // namespace Oblio
