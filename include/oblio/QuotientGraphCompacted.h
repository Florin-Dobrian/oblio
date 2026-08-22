#pragma once

// QuotientGraphCompacted.h - the quotient graph on AMD_2'S CLIQUE LAYOUT: one pooled workspace
// holding a vertex's two lists AND every live clique, with a free cursor and a compaction,
// where QuotientGraph keeps C[c] in a separate append-only arena that only grows.
//
// SAME GRAPH, SAME ALGORITHMS, DIFFERENT STORAGE. Every idea in QuotientGraph.h holds here without
// change: the reachable set, the two kinds of source, a clique id being the pivot that formed it,
// a dead clique being an empty member list. Read that file first; this one comments only what the
// storage makes different, which is where the cliques live and what has to happen when the pool
// fills.
//
// WHAT THE POOL BUYS, and why the machinery is worth writing. The pool is sized ONCE from the
// pattern, `nzaat + nzaat/5 + n`, and never grows: if A fits, the ordering completes. Our arena
// cannot promise that, since nnz(L) depends on the ordering being computed and so cannot be
// bounded before the run. That is the argument in docs/DESIGN_DECISIONS.md (2026-08-16) and it is
// why both vendored routines carry reclamation machinery.
//
// WHAT IT COSTS is a compactor, and the compactor is what makes this class harder than the flat
// one. Absorbed cliques and consumed list prefixes leave dead space, the cursor eventually reaches
// the end, and everything live has to slide down. So every walk here is written in POSITIONS and
// CURSORS rather than pointers and counters: a compaction moves blocks within a pool that never
// reallocates, so a position survives it and a pointer does not, and a truncation shortens the
// list being read, so a count of what is LEFT survives and a count of what has been done does not.
//
// ONE CLASS, TWO BRANCHES, AND THE SPLIT IS BY SUFFIX. `Amd3B` and `Mmd3C` are the two drivers.
// Where the vendored codes AGREE there is one method; where they disagree there are two, named
// `...Amd` and `...Mmd`, so that each body can be read against `AMD_2` or `mmdelm` directly rather
// than against a flag that means one thing on one branch and another on the other. The count of
// splits is not a constraint: contorting shared code to avoid a suffix is worse than the suffix.
// See docs/DESIGN_DECISIONS.md (2026-08-19).
//
// THE FOUR THINGS THAT GENUINELY DIFFER, and they are all order or lifetime rather than storage:
//
//   - WHERE THE NEW CLIQUE LANDS IN I[u]. `AMD_2` writes it first, genmmd appends it last. That
//     decides the order members enter C[pivot], hence bucket order, hence which of several
//     equal-degree vertices is picked, so neither is free to move: each driver reproduces its
//     reference's permutation exactly or it is not a differential. `pruneAmd` and `pruneMmd`.
//   - WHICH HALF OF THE RUN COMES FIRST, and which way the incidence list is walked. `AMD_2` lays
//     the run out I[u] then A[u] and walks the cliques forward; genmmd lays it A[u] then I[u] and
//     walks the cliques BACKWARD, `mmdelm` popping a stack. `adjacencyAmd`/`adjacencyMmd` with
//     their incidence twins, and `reachableSetAmd`/`reachableSetMmd`.
//
//     The layout is these accessors and the walks and NOTHING ELSE. The pool, the cursor, the
//     elbow room, the compactor's sweep, the descriptors and the counters are all blind to it.
//   - HOW A DEAD VERTEX IS RECOGNIZED. Amd reads a zero weight, which is `AMD_2`'s `Nv [i] == 0`
//     and costs no array at all. Mmd cannot: `number()` leaves a prepass vertex LIVE at weight one
//     and in every list that names it, so only a tag can hide it. `eliminatedAmd` and
//     `eliminatedMmd`, and `mMark` exists for the second alone.
//   - WHEN THE NEGATED WEIGHTS ARE RESTORED. A walk marks membership by flipping a weight's sign,
//     and amd calls the walk ONCE PER PIVOT and restores in a pass it already makes, while mmd
//     calls `reachableWeight` PER VERTEX in its degree refresh and so needs each call to clean up
//     after itself. `massEliminateAmd` and `massEliminateMmd`.
//
// THE MARK ARRAY IS ALLOCATED ON DEMAND, by `enableMarks`, so the amd driver pays a pointer rather
// than n int32 and its footprint claim stays true. `mHasNumbered` IS THEREFORE LOAD BEARING and
// not merely an optimization: it is the short circuit that keeps an empty array safe in the shared
// bodies carrying the `mHasNumbered && mMark[v] == GONE` guard.

#include "oblio/QuotientGraph.h"   // Buckets and TaggedScan, which are shared verbatim
#include "oblio/Types.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

class QuotientGraphCompacted {
public:
    // NO cliqueMarks ARGUMENT, unlike QuotientGraph's. That flag sized a mark array at n or 2n;
    // here the array is absent unless a driver asks for it, and no driver wants clique marks.
    QuotientGraphCompacted(const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const { return mRun.size(); }

    // THE POOL'S SIZE, in place of the flat class's arena entries. Not the same quantity and not
    // meant to be: an arena grows monotonically and this is fixed at construction, so what this
    // reports is the elbow room the layout needed rather than the space the cliques consumed.
    // `compactions()` is the number that says whether the room was enough.
    std::size_t arenaEntries() const { return mSource.size(); }
    std::size_t compactions() const  { return mCompactions; }

    // A BLOCK HEAD IN THE POOL, during compaction only: FLIPPED - e for owner e, so every
    // head is at most FLIPPED and every real entry, a vertex id, is above it. AMD_2's FLIP.
    static constexpr std::int32_t FLIPPED = -2;

    // GONE, genmmd's `marker[v] = maxint`: one value above every reachable tag, so the stamp array
    // answers "is v dead" on the load it was making anyway. Read by the mmd branch alone.
    static constexpr std::int32_t GONE = 2147483647;   // INT32_MAX

    // ------------------------------------------------------------------ per-branch: liveness

    // ZERO WEIGHT IS THE DEAD STATE on the amd branch, which is `Amd.cpp`'s `Nv [i] == 0`. The
    // three ways a vertex leaves the graph there all end in a zero weight, so no array is spent.
    bool eliminatedAmd(std::int32_t u) const { return mWeight[u] == 0; }

    // AND A TAG ON THE MMD BRANCH, because a zero weight is not available to it: `number()` leaves
    // a prepass vertex live at weight one so that its neighbors' degrees still count it.
    bool eliminatedMmd(std::int32_t u) const { return mMark[u] == GONE; }

    // ------------------------------------------------------------------ per-branch: the layout

    // INCIDENCE FIRST, ADJACENCY BEHIND IT, which is `AMD_2`'s order. The flip is what lets the
    // new clique go in by its three-move rotation and what makes the part a cliques-first walk has
    // consumed a PREFIX, which is what the mid-walk compactor needs.
    const std::int32_t* incidenceAmd(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr;
    }
    const std::int32_t* adjacencyAmd(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr + mRun[u].incidenceSize;
    }

    // ADJACENCY FIRST, INCIDENCE BEHIND IT, which is genmmd's order and the flat class's. The
    // compactor needs no flip here: mmd walks A[u] forward and I[u] backward, so the consumed part
    // is a prefix in one phase and a suffix in the other and the remainder is contiguous in both.
    const std::int32_t* adjacencyMmd(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr;
    }
    const std::int32_t* incidenceMmd(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr + mRun[u].adjacencySize;
    }

    // THE LENGTHS ARE NOT SPLIT. A length is a count and says nothing about where the half sits.
    std::uint32_t adjacencySize(std::int32_t u) const { return mRun[u].adjacencySize; }
    std::uint32_t incidenceSize(std::int32_t u) const { return mRun[u].incidenceSize; }

    // ------------------------------------------------------------------ cliques

    // A clique has no incidence part, so its block starts at the run under either layout and this
    // needs no suffix. Spelled without the offset because a clique is one list rather than two.
    const std::int32_t* clique(std::int32_t c) const {
        return mSource.data() + mRun[c].sourcePtr;
    }
    std::uint32_t cliqueSize(std::int32_t c) const { return mRun[c].adjacencySize; }

    // WRITABLE, for a restore pass that trims the clique as it walks it. See `trimClique`.
    std::int32_t* clique(std::int32_t c) { return mSource.data() + mRun[c].sourcePtr; }
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

    // RESTORING A MEMBER'S SIGN, amd only, called from the pass that pass already makes over
    // C[pivot]. The mmd branch restores inside `massEliminateMmd` instead; see the header note.
    std::uint32_t restoreWeight(std::int32_t u) {
        const std::int32_t w = -mWeight[u];
        mWeight[u] = w;
        return static_cast<std::uint32_t>(w);
    }
    void restorePivotWeight(std::int32_t pivot) { mWeight[pivot] = -mWeight[pivot]; }

    // A ROW THE DENSE RULE SET ASIDE, amd only. `AMD_2` writes `Nv [i] = 0` and the row is then
    // absent from every reachable set and every prune, which a zero weight already achieves here.
    void setAside(std::int32_t u) { mWeight[u] = 0; }

    // ------------------------------------------------------------------ the mark array, mmd only

    // ALLOCATED ON DEMAND. The amd branch never calls this and so carries an empty vector rather
    // than n int32. Call once, before any elimination.
    void enableMarks();

    std::int32_t advanceTag()                    { return ++mTag; }
    std::int32_t mark(std::int32_t v) const      { return mMark[v]; }
    void setMark(std::int32_t v, std::int32_t t) { mMark[v] = t; }

    // A VERTEX NUMBERED BY THE MMD PREPASS, which is not eliminated in the quotient-graph sense:
    // it keeps its weight and its place in every list that names it, and only the tag hides it.
    void number(std::int32_t u);

    // ------------------------------------------------------------------ elimination

    // THE THREE STEPS OF AN ELIMINATION, called in this order by both drivers. There is no
    // `eliminate` wrapper: amd's prune computes the degree bound and the hash key in the same pass
    // and mmd's does neither, so one signature would have to carry a parameter one branch ignores.
    // Spelling the sequence at the call site costs two lines and invents nothing.
    //
    // THE FIRST STEP SPLITS AND IT IS THE COSTLIEST OF THE SUFFIXES, because the placement rule it
    // holds is genuinely common: build in place when the incidence list is empty, otherwise at the
    // free cursor, compact when the cursor runs out. What cannot be shared is the two lines that
    // NAME A WALK, and a walk cannot be selected without knowing the branch. So the tail and the
    // capture are factored into private helpers and the six lines of the rule read twice, each
    // against its own reference, which is the trade docs/NEXT.md asked for.
    void beginEliminationAmd(std::int32_t pivot);
    void beginEliminationMmd(std::int32_t pivot);

    void pruneAmd(std::int32_t pivot, TaggedScan& scan);
    void pruneMmd(std::int32_t pivot);

    // ONE METHOD FOR BOTH BRANCHES. `markGone` inside it does nothing when the mark array is
    // absent, which is the amd branch, so the two halves this was split into differed in a store
    // that was already a no-op on one of them.
    const std::vector<std::int32_t>& finishElimination(std::int32_t pivot);

    // AND THE THREE AS ONE CALL, overloaded on the scan exactly as `QuotientGraph`'s pair is: the
    // mmd prune takes nothing and the amd one takes a `TaggedScan`, so the argument selects the
    // branch. The value of the wrapper is that the ORDER of the three steps lives here rather
    // than in each driver, where nothing could enforce it. The three remain public and remain
    // suffixed; this is a sequence rather than a replacement.
    const std::vector<std::int32_t>& eliminate(std::int32_t pivot);
    const std::vector<std::int32_t>& eliminate(std::int32_t pivot, TaggedScan& scan);

    // MASS ELIMINATION, one method for both branches. Whoever runs it restores the negated
    // weights: this method when the eliminator runs it eagerly, the driver when it runs it late.
    // See the body, and `restoreWeight`.
    const std::vector<std::int32_t>& massEliminate(std::int32_t pivot);

    // AGGRESSIVE ABSORPTION, amd only: a clique whose external degree has reached zero is dead,
    // and its members drop it from their incidence lists.
    void absorb(const std::vector<std::int32_t>& cliques,
                const std::int32_t* vertices, std::uint32_t vertexCount);

    // ONE SUPERVARIABLE ABSORBS ANOTHER. `v` is a live vertex that never formed a clique, so its
    // length is A[v]'s and NOT a clique's, which is why this does not call `killClique`.
    void merge(std::int32_t u, std::int32_t v);

    // Whether mass elimination is the driver's rather than this class's. Amd sets it: its prune
    // needs the untrimmed clique for the degree bound, so the trim happens later and in the driver.
    void setLateMassElimination(bool on) { mLateMassElimination = on; }

    // Whether I[u] is walked backward, which is `mmdelm` popping its clique stack. Set by the mmd
    // driver; the amd walk has its own method and does not read it.
    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    // ------------------------------------------------------------------ the reachable set

    // THE WALK THAT BUILDS THE NEW CLIQUE, one per branch, differing in the order of the two
    // halves and the direction of the incidence part. Both write at the free cursor, both may
    // COMPACT MID-WALK, and both return the number of members reached, the start having moved.
    std::uint32_t reachableSetAmd(std::int32_t u);
    std::uint32_t reachableSetMmd(std::int32_t u);

    // AND THE IN-PLACE FORM, `AMD_2`'s `elenme == 0` branch: a pivot with no cliques has a reach
    // that is a SUBSET of A[pivot], so it is compacted where it stands and the pool is untouched.
    // The mmd form additionally rejects a prepass-numbered vertex, which a weight cannot exclude.
    std::uint32_t reachableSetInPlaceAmd(std::int32_t u);
    std::uint32_t reachableSetInPlaceMmd(std::int32_t u);

    // THE DEGREE OF A REACHABLE SET WITHOUT BUILDING IT, mmd only, called per vertex in the
    // refresh. It stamps rather than negating, so it leaves the weights untouched.
    std::uint32_t reachableWeight(std::int32_t u);

    // ------------------------------------------------------------------ counters and output

    // PEAK LIVE CLIQUE MEMBERS. A MEMBER is a vertex in a live clique at this instant, where an
    // ENTRY is a pool slot; the two differ here and not in the flat class, because this layout
    // reclaims. IT IS A PROPERTY OF THE ALGORITHM AND NOT OF THE LAYOUT, which is what makes it
    // comparable against `Amd3`'s and `Mmd3`'s: two drivers agreeing on a permutation can still be
    // caught doing different work, and have been. Checked in tests/test_order.cpp.
    std::size_t numPeakCliqueMembers() const { return mNumPeakCliqueMembers; }
    std::size_t numLiveCliqueMembers() const { return mNumLiveCliqueMembers; }

    // THE COUNTER CHECKED AGAINST A RECOMPUTATION, debug builds only. Five sites move the count: a
    // birth in `beginElimination`, a death in `killClique`, the shrink in mass elimination, the
    // trim in `trimClique`, and the mid-walk truncation inside a walk. A funnel is a claim about
    // the CALL GRAPH and holds only while every writer goes through it, which twice it did not;
    // the way to make such a claim checkable is a recomputation from independent state, and this
    // is it. See docs/DESIGN_DECISIONS.md (2026-08-19).
    //
    // COMPILED OUT under NDEBUG, where the assertion that calls it is compiled out too.
#ifndef NDEBUG
    bool cliqueCountBalances() const;
#endif

    // THE PERMUTATION, one per branch, and the difference is not the layout but what the driver's
    // pivot list holds. `order` expands each pivot's supervariable chain; `orderAscending` does the
    // same and additionally emits a prepass-numbered vertex as the run of one it stands for.
    std::vector<std::int32_t> order(const std::vector<std::int32_t>& pivots) const;
    std::vector<std::int32_t> orderAscending(const std::vector<std::int32_t>& pivots) const;

private:
    // COMPACTION, `AMD_2`'s `garbage_collection`, and it takes the half-built clique's start BY
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
    void markGone(std::int32_t v) { if (!mMark.empty()) mMark[v] = GONE; }

    // THE POOL. Every vertex's run and every live clique, end to end, sized once at construction
    // to `nzaat + nzaat/5 + n` and never grown. `mFree` is `AMD_2`'s `pfree`.
    std::vector<std::int32_t> mSource;

    struct VertexRun {
        // WHERE u'S RUN STARTS, AND IT MOVES. Three things move it: `compact` sliding every
        // live block down, a walk truncating the pivot's run to what it has not consumed, and
        // `beginElimination` repointing the pivot at its new clique. Anything holding a POSITION
        // survives all three by re-reading this field; a pointer does not.
        std::size_t   sourcePtr;
        std::uint32_t adjacencySize;   // A[u]'s length, or a clique's member count
        std::uint32_t incidenceSize;   // I[u]'s length; zero for a clique
    };
    std::vector<VertexRun> mRun;

    std::size_t mFree = 0;             // AMD_2's pfree
    std::size_t mCompactions = 0;      // AMD_2's Info[AMD_NCMPA]

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

    // THE WEIGHT, CARRYING `AMD_2`'s `Nv` ENCODING IN ITS SIGN, which is why it is signed:
    //
    //     mWeight[v] >  0    live, and not yet taken into the clique being built; the weight
    //     mWeight[v] <  0    live, and taken into it this step; the weight is -mWeight[v]
    //     mWeight[v] == 0    dead, by a hash merge or by mass elimination, or set aside as dense
    //
    // ONE LOAD ANSWERS THREE QUESTIONS in the two hottest loops here. The ZERO is a true sentinel
    // and holds on the amd branch alone, where the only death sites are `merge` and mass
    // elimination; on the mmd branch `number()` leaves a prepass vertex live at weight one, which
    // is why that branch keeps `mMark`. See docs/CODING_RULES.md for the four conditions a signed
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
// they can be compared with EACH OTHER; it says nothing about a comparison against a vendored
// routine, which differs from us in ways nobody has enumerated. Every out-of-class definition below
// is `inline` so that several drivers may include this and the linker folds the copies.
//
// The measurement that decided it, and what it cost in compile time, is in
// docs/DESIGN_DECISIONS.md, and `docs/CODING_RULES.md` carries the rule and the four mechanics it
// needs.
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
// congruent. Different amounts is what actually perturbs the arrangement. The step is a whole page
// on Apple Silicon, 16 KiB, because an intervention has to be large enough to intervene: the
// 2026-08-17 investigation spent two rounds on other hypotheses after a sixteen-int pad moved
// nothing, which was not evidence of anything, the allocator having returned the same addresses.
//
// See docs/DESIGN_DECISIONS.md (2026-08-17), where this mechanism cost `Mmd3C` 28 per cent at
// exactly 200 a side and nothing at 199 or 201.
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
                             const std::vector<std::int32_t>& rowIdx) {
    detail::padded(mRun, colPtr.empty() ? 0 : colPtr.size() - 1, 0);
    const std::int32_t size = static_cast<std::int32_t>(mRun.size());

    // ELBOW ROOM, `AMD_2`'S EXACTLY. It sizes `slen = nzaat + nzaat/5 + 7n` and hands
    // `iwlen = slen - 6n` to the main loop, so the pool is `nzaat + nzaat/5 + n` where `nzaat` is
    // the OFF-DIAGONAL entry count, `sum (Len [0..n-1])`. Ours holds the same lists, so the figure
    // transfers, and it has to be the same figure rather than the same shape: the headroom is what
    // decides how often the compactor runs, so a differential against `AMD_2` on compaction counts
    // measures the headroom unless the two agree. Until 2026-08-18 ours was computed from `nnz`
    // WITH the diagonal and so ran about 1.2n large.
    //
    // Reserved from an upper bound and then sized down, rather than counted in a pass of its own:
    // `nnz >= nzaat`, so the reserve cannot reallocate, and `nzaat` is known once the runs are
    // laid out. The vector is SIZED rather than left at capacity because the space past the runs
    // is written by cliques and read by nothing until it is, and a cursor into a vector's unused
    // capacity would be undefined.
    const std::size_t nnz = colPtr.empty() ? 0 : colPtr.back();
    mSource.reserve(nnz + nnz / 5 + mRun.size() + detail::kOrderingPad);
    for (std::int32_t aj = 0; aj < size; ++aj) {
        mRun[aj].sourcePtr = mSource.size();
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSource.push_back(rowIdx[cp]);
        mRun[aj].adjacencySize = static_cast<std::uint32_t>(mSource.size() - mRun[aj].sourcePtr);
    }

    // The runs are laid out; everything past them is free space the cliques will use.
    const std::size_t nzaat = mSource.size();
    mFree = nzaat;
    mSource.resize(nzaat + nzaat / 5 + mRun.size());

    detail::padded(mSuperNext, static_cast<std::size_t>(size), 1);
    detail::padded(mSuperLast, static_cast<std::size_t>(size), 2);
    detail::padded(mWeight,    static_cast<std::size_t>(size), 3);
    std::fill(mSuperNext.begin(), mSuperNext.end(), NIL);
    std::fill(mWeight.begin(),    mWeight.end(),    1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

// THE SAME RECLAIM, ONE PASS LATER AND FOR THE OTHER REASON. `massEliminate` above drops the
// members it merged into the pivot; this drops the members supervariable detection absorbed into
// each other, which happens after it. `AMD_2` does both in one place because its detection sits
// inside the scan, so its RESTORE DEGREE LISTS pass sees every casualty at once and writes the
// survivors back with `Iw [p++] = i`, then `Len [me] = p - pme1` and `if (elenme != 0) pfree = p`.
// Ours needs two because detection is a pass of its own.
//
// WITHOUT THIS THE ABSORBED STAY IN THE CLIQUE FOR THE REST OF THE RUN. No permutation moves, every
// later walk skipping them on the `nv > 0` test, but they are visited, and the space behind them
// is never given back, so the pool fills faster and the compactor runs sooner than `AMD_2`'s.
inline void QuotientGraphCompacted::trimClique(std::int32_t pivot, std::uint32_t kept) {
    const std::uint32_t was = mRun[pivot].adjacencySize;
    if (kept == was) return;
    mNumLiveCliqueMembers -= was - kept;
    mRun[pivot].adjacencySize = kept;
    if (!mBuiltInPlace) {
        assert(mRun[pivot].sourcePtr + was == mFree &&
               "the pivot's clique is no longer the last block in the pool");
        mFree = mRun[pivot].sourcePtr + kept;
    }
}

inline void QuotientGraphCompacted::killClique(std::int32_t c) {
    mNumLiveCliqueMembers -= mRun[c].adjacencySize;
    mRun[c].adjacencySize = 0;
}

// ALLOCATED ON DEMAND, and the amd branch never calls this. `NIL` rather than zero as the initial
// stamp, since zero is a tag a walk can reach and this array's whole job is to answer "have I seen
// v this step" against `mTag`, which starts there.
inline void QuotientGraphCompacted::enableMarks() {
    detail::padded(mMark, mRun.size(), 4);
    std::fill(mMark.begin(), mMark.end(), NIL);
}

inline void QuotientGraphCompacted::number(std::int32_t u) {
    mHasNumbered = true;
    mMark[u]     = GONE;
}

// THE PLACEMENT RULE, WRITTEN TWICE. Everything about building a clique is common except the two
// lines that NAME A WALK, and a walk cannot be selected without knowing the branch. So the capture
// before it and the bookkeeping after it are shared helpers and only the rule itself reads twice.
//
// TWO WAYS TO BUILD, and which applies is `AMD_2`'s `elenme == 0`. With no cliques in I[pivot] the
// reach is a SUBSET of A[pivot] and is compacted where it stands, so the pool is not touched at
// all; otherwise it is assembled at the free cursor. Measured on grids, 62 to 68 per cent of
// eliminations take the first path and the share rises with n, which is why the branch is worth
// having and not merely faithful.
//
// NO RESERVATION BEFORE THE WALK. The walks test per entry and compact from inside themselves,
// which is `AMD_2`'s `if (pfree >= iwlen) garbage_collection`. Reserving room for a worst-case
// reach of n beforehand, which is what a walk holding pointers had to do, compacted where `AMD_2`
// would not.
inline void QuotientGraphCompacted::beginEliminationAmd(std::int32_t pivot) {
    captureAbsorbed(incidenceAmd(pivot), mRun[pivot].incidenceSize);

    // `Nv [me] = -nvpiv` in Amd.cpp, and it is what keeps the pivot out of its own clique: the
    // walks take a vertex only when its weight reads positive.
    mWeight[pivot] = -mWeight[pivot];

    const bool inPlace = mRun[pivot].incidenceSize == 0;
    mBuiltInPlace = inPlace;

    std::size_t   cliqueStart;
    std::uint32_t cliqueLen;
    if (inPlace) {
        cliqueStart = mRun[pivot].sourcePtr;
        cliqueLen   = reachableSetInPlaceAmd(pivot);
    } else {
        // AFTER the walk, not before it. The cursor is where the clique ends, so its start is a
        // subtraction, and a compaction inside the walk moves both together and leaves this right.
        cliqueLen   = reachableSetAmd(pivot);
        cliqueStart = mFree - cliqueLen;
    }
    bearClique(pivot, cliqueStart, cliqueLen);

    // A CLIQUE HAS NO INCIDENCE LIST. The mmd branch clears this in `finishEliminationMmd`
    // instead, its mass elimination reading the pivot's run in between.
    mRun[pivot].incidenceSize = 0;
}

// The same rule against genmmd. The walk order is the difference and the note is in the header.
inline void QuotientGraphCompacted::beginEliminationMmd(std::int32_t pivot) {
    captureAbsorbed(incidenceMmd(pivot), mRun[pivot].incidenceSize);

    mWeight[pivot] = -mWeight[pivot];

    const bool inPlace = mRun[pivot].incidenceSize == 0;
    mBuiltInPlace = inPlace;

    std::size_t   cliqueStart;
    std::uint32_t cliqueLen;
    if (inPlace) {
        cliqueStart = mRun[pivot].sourcePtr;
        cliqueLen   = reachableSetInPlaceMmd(pivot);
    } else {
        cliqueLen   = reachableSetMmd(pivot);
        cliqueStart = mFree - cliqueLen;
    }
    bearClique(pivot, cliqueStart, cliqueLen);
}

// THE PRUNE, `AMD_2`'s scan 2 with its scan 1 fused in. One walk of each vertex of C[pivot]
// rewrites both halves of its run, accumulates the degree bound's explicit part, builds the
// adjacency half of the hash key, and leaves the incidence half in the tagged `w`.
//
// THE ABSORBED CLIQUES' `w` IS ZEROED HERE, from the copy `beginEliminationAmd` took, because the
// compaction below reads `w == 0` as "absorbed and gone", which is `Amd.cpp`'s `W [e] == 0`. It
// used to ride in the capture pass itself, which is the same store one call earlier and nothing
// between the two reads `w`.
inline void QuotientGraphCompacted::pruneAmd(std::int32_t pivot, TaggedScan& scan) {
    for (std::int32_t c : mAbsorbed) scan.w[c] = 0;

    const std::int32_t* reached     = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;

    if (scan.buckets != nullptr)
        for (std::uint32_t ri = 0; ri < reachedSize; ++ri) scan.buckets->unfile(reached[ri]);

    const std::int32_t  wflg        = scan.wflg;          // hoisted, as the flag is

    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u       = reached[ri];
        const std::int32_t nvi     = -mWeight[u];   // u is in C[pivot], so its weight is negative
        const std::int32_t wnvi    = wflg - nvi;          // Amd.cpp's wnvi, and signed for it
        const VertexRun&  run           = mRun[u];
        std::int32_t*     source        = mSource.data() + run.sourcePtr;
        const std::uint32_t incidenceSize = run.incidenceSize;
        const std::uint32_t adjacencySize = run.adjacencySize;
        const std::int32_t* adjacency     = source + incidenceSize;

        // THE INCIDENCE PART, AT THE FRONT OF THE RUN. Compacted where it lies: the write cursor
        // starts at the read cursor and advances only when an entry is kept, so it never passes it.
        std::uint32_t write = 0;
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = source[i];
            std::int32_t we = scan.w[c];
            if (we == 0) continue;                        // absorbed and gone; Amd.cpp's W == 0
            source[write++] = c;
            if (c == pivot) continue;                     // the new clique subtracts from nothing
            if (we >= wflg) {
                we -= nvi;
            } else {
                we = static_cast<std::int32_t>(scan.degree[c]) + wnvi;
                scan.touchedCliques.push_back(c);
            }
            scan.w[c] = we;
        }
        const std::uint32_t incidenceKept = write;        // Amd.cpp's p3, the boundary

        // THE ADJACENCY PART, BEHIND IT, and the write cursor carries on through the boundary. It
        // is still behind the read cursor, having entered this loop at most `incidenceSize`.
        std::uint32_t explicitPart = 0;                   // a weight sum, not a count of positions
        std::uint32_t key          = 0;                   // wraps, like Amd.cpp's UInt hval
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            // ONE LOAD, THREE QUESTIONS, which is Amd.cpp's `nvj = Nv [j] ; if (nvj > 0)`. A
            // negative weight is a member of the new clique, including the pivot itself, so the
            // explicit `v == pivot` test is gone with the rest; a zero is absorbed or merged away.
            const std::int32_t v  = adjacency[k];
            const std::int32_t nv = mWeight[v];
            if (nv <= 0) continue;
            explicitPart += static_cast<std::uint32_t>(nv);
            key += static_cast<std::uint32_t>(v);         // no + 1, no reduction; see above
            source[write++] = v;
        }
        scan.w[u] = static_cast<std::int32_t>(explicitPart);
        if (scan.buckets != nullptr)
            scan.buckets->setKey(u, static_cast<std::int32_t>(key));   // the ADJACENCY half

        // THE NEW CLIQUE GOES IN BY ROTATION, WHICH IS `AMD_2`'S THREE MOVES EXACTLY:
        //
        //     Iw [pn] = Iw [p3] ;   the first adjacency entry moves to the end
        //     Iw [p3] = Iw [p1] ;   the first incidence entry moves to the boundary
        //     Iw [p1] = me ;        and the pivot takes the front
        //
        // There is no free slot to insert into, so a rotation is what puts the pivot at the front
        // without a shift, and it is why each part reads with its first entry last. The old code
        // reproduced the same two orders under the reversed layout by holding the first adjacency
        // survivor back and swapping the pivot into place afterwards; with the parts in `AMD_2`'s
        // order the rotation says it directly and takes a test out of the adjacency loop.
        //
        // ROOM IS GUARANTEED BY ONE DROPPED ENTRY, and there is always at least one: u is in
        // C[pivot], so either the pivot was in A[u] and the walk above dropped it, or a clique of
        // I[u] was absorbed into the new one and was dropped as `we == 0`. Both self-assign
        // harmlessly when a part comes out empty.
        assert(write < incidenceSize + adjacencySize &&
               "the rotation has no slot: nothing was dropped from u's run");
        source[write]         = source[incidenceKept];
        source[incidenceKept] = source[0];
        source[0]             = pivot;
        mRun[u].incidenceSize = incidenceKept + 1;
        mRun[u].adjacencySize = write - incidenceKept;
    }

}

// THE PRUNE, `mmdelm`'s. One walk of each vertex of C[pivot] drops the dead from A[u], keeps the
// live cliques of I[u], and APPENDS the new clique at the back, which is genmmd's convention and
// the reverse of `AMD_2`'s. No degree bound and no hash key: this branch refreshes degrees in a
// pass of its own.
inline void QuotientGraphCompacted::pruneMmd(std::int32_t pivot) {
    const std::int32_t* reached     = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;

    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        const VertexRun&  run           = mRun[u];
        std::int32_t*     source        = mSource.data() + run.sourcePtr;
        const std::uint32_t adjacencySize = run.adjacencySize;
        std::uint32_t       kept          = 0;
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (mWeight[v] <= 0) continue;             // in the new clique, the pivot, or merged
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by a prepass; see above
            source[kept++] = v;
        }
        mRun[u].adjacencySize = kept;                      // A[u] - C[pivot] - {pivot}

        const std::int32_t* incidence     = source + adjacencySize;
        const std::uint32_t incidenceSize = run.incidenceSize;
        std::uint32_t       write         = kept;      // follows `kept`; see the note above
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
            if (mRun[c].adjacencySize != 0) source[write++] = c;   // dead is size zero; see above
        }
        // THE NEW CLIQUE GOES AT THE BACK, which is genmmd's convention and the reverse of
        // `AMD_2`'s. It is what the reversed incidence walk then reads first, and it decides the
        // order members enter C[pivot], hence bucket order, hence which of several equal-degree
        // vertices is picked. The held-vertex rotation that put it at the FRONT for `AMD_2` went
        // with `mVendoredListOrder` on 2026-08-19; see that member's note.
        source[write++]   = pivot;                     // u joins the new clique, id = pivot
        mRun[u].incidenceSize = write - kept;
    }

}

// MASS ELIMINATION IS THE DRIVER'S ON THIS BRANCH, so this normally only clears the scratch: the
// amd prune needs the untrimmed clique for its degree bound, so the trim happens later and in the
// driver. The flag is still tested rather than assumed, the class not being the place that decides.
inline const std::vector<std::int32_t>&
QuotientGraphCompacted::finishElimination(std::int32_t pivot) {
    if (mLateMassElimination) {
        mMerged.clear();
    } else {
        massEliminate(pivot);
    }

    mRun[pivot].incidenceSize = 0;
    markGone(pivot);
    return mMerged;
}

inline const std::vector<std::int32_t>& QuotientGraphCompacted::eliminate(std::int32_t pivot) {
    beginEliminationMmd(pivot);
    pruneMmd(pivot);
    return finishElimination(pivot);
}

inline const std::vector<std::int32_t>& QuotientGraphCompacted::eliminate(std::int32_t pivot,
                                                                          TaggedScan& scan) {
    beginEliminationAmd(pivot);
    pruneAmd(pivot, scan);
    return finishElimination(pivot);
}

// MASS ELIMINATION, `AMD_2`'s "eliminate non principal supervariables". A vertex whose reach is
// exactly the new clique is indistinguishable from the pivot and is absorbed into it.
//
// THE AMD FORM DOES NOT RESTORE THE NEGATED WEIGHTS, its driver doing that in a pass it already
// makes over C[pivot]. The mmd form must, `reachableWeight` being called per vertex in the refresh.
inline const std::vector<std::int32_t>&
QuotientGraphCompacted::massEliminate(std::int32_t pivot) {
    std::vector<std::int32_t>& merged = mMerged;   // scratch, kept for its capacity
    merged.clear();
    const std::int32_t* reached     = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;
    // THE SIGNS COME BACK HERE UNLESS THE DRIVER IS DOING IT. Mass elimination is the last reader
    // of the negated form, so whoever runs it restores: this method when the eliminator runs it
    // eagerly, the driver when it runs it late and has a refile pass to ride the store on.
    // `AMD_2` is the late case, mass-eliminating with the weights still negative, `nvi = -Nv [i]`,
    // and restoring afterwards under RESTORE DEGREE LISTS. See restoreWeight.
    const bool restore = !mLateMassElimination;
    if (restore) mWeight[pivot] = -mWeight[pivot];
    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        if (restore) mWeight[u] = -mWeight[u];             // live again, and positive
        if (mRun[u].adjacencySize == 0 && mRun[u].incidenceSize == 1 &&
            mSource[mRun[u].sourcePtr] == pivot) {         // A[u] empty, so I[u] starts at the run
            mRun[u].incidenceSize = 0;
            markGone(u);
            merged.push_back(u);
        }
    }
    if (!merged.empty()) {
        for (std::int32_t u : merged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            mCliqueWeight -= weight(u);            // magnitude: the sign depends on `restore`
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }

        std::int32_t*     members     = mSource.data() + mRun[pivot].sourcePtr;
        const std::uint32_t membersSize = mRun[pivot].adjacencySize;
        std::uint32_t       kept        = 0;
        for (std::uint32_t k = 0; k < membersSize; ++k)
            if (mWeight[members[k]] != 0) members[kept++] = members[k];
        mNumLiveCliqueMembers -= membersSize - kept;   // a shrink is a partial death
        mRun[pivot].adjacencySize = kept;

        // THE SPACE THE COMPACTION FREED GOES BACK TO THE CURSOR, which is `AMD_2`'s
        // `if (elenme != 0) pfree = p`. Only in the pooled case; the in-place case never took any.
        // Safe because the clique is still the last block in the pool: nothing writes to the pool
        // between beginElimination building it and this function trimming it. Asserted rather than
        // assumed, since if that stops being true this goes wrong silently.
        if (!mBuiltInPlace) {
            assert(mRun[pivot].sourcePtr + membersSize == mFree &&
                   "the pivot's clique is no longer the last block in the pool");
            mFree = mRun[pivot].sourcePtr + kept;
        }
    }
    return merged;
}

inline void QuotientGraphCompacted::absorb(const std::vector<std::int32_t>& cliques,
                           const std::int32_t* vertices, std::uint32_t vertexCount) {
    if (cliques.empty()) return;

    // dead; the compaction reads the size
    for (std::int32_t c : cliques) killClique(c);

    // I[u] LESS THE DEAD, COMPACTED RIGHTWARD SO THE ADJACENCY NEVER MOVES. The incidence part is
    // at the front of the run and the adjacency starts at `sourcePtr + incidenceSize`, so
    // compacting to the LEFT would drag A[u] down behind it, a copy per surviving vertex per
    // absorbing step. Compacting to the RIGHT instead leaves the boundary where it is: survivors
    // end flush against A[u], `sourcePtr` advances by however many were dropped, and the freed
    // space falls off the FRONT of the run rather than the back. The adjacency is not touched at
    // all, and I[u] is the short part after a prune, holding the new clique plus whatever survived.
    //
    // THE PASS ITSELF IS STILL A DIVERGENCE. `AMD_2` absorbs inside scan 2, in the same walk that
    // rewrites the list, so an absorbed clique is simply not copied. We decide absorption after
    // the whole prune, because the prune has scan 1 fused into it, and so have to revisit the list.
    // What is removed here is the slide, not the visit.
    for (std::uint32_t k = 0; k < vertexCount; ++k) {
        const std::int32_t u         = vertices[k];
        std::int32_t*      incidence = mSource.data() + mRun[u].sourcePtr;
        const std::uint32_t size     = mRun[u].incidenceSize;

        // THE ROTATION HAS TO BE REDONE WHEN ITS SUBJECT DIES. `AMD_2` rotates inside scan 2, with
        // this step's absorbed cliques already dropped, so the entry it parks at the back is the
        // first SURVIVOR; our prune rotates first and absorbs here, so when the parked entry is
        // one this pass removes the two lists differ. See the shared class for the measurement
        // that found it, and note the same correction has to exist in both files or their columns
        // stop being comparable.
        const bool parkedDied = size > 0 && mRun[incidence[size - 1]].adjacencySize == 0;

        std::uint32_t write = size;
        for (std::uint32_t i = size; i-- > 0;)
            if (mRun[incidence[i]].adjacencySize != 0) incidence[--write] = incidence[i];
        mRun[u].sourcePtr    += write;                  // the dropped entries fall off the front
        mRun[u].incidenceSize = size - write;

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
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    // NOT `killClique`. v is a live supervariable being absorbed and never formed a clique, so
    // this length is A[v]'s, not a clique's, and feeding it to the counter would corrupt the peak.
    mRun[v].adjacencySize = 0;
    mRun[v].incidenceSize = 0;
    markGone(v);          // mmd only; the amd branch reads the zero weight above and has no array
}

// WRITES AT `mFree` RATHER THAN APPENDING, which is the storage change. The caller still ensures
// room for a whole reach before the walk starts, so the cursor cannot run past the pool inside it;
// that test is the last divergence from `AMD_2`, which tests per entry instead, and closing it is
// what the positions below are for.
//
// INCIDENCE THEN ADJACENCY, which is `AMD_2`'s `for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++)`:
// the cliques of me on the first `elenme` passes and the supervariables on the last. It is now
// also the physical order of the run, so the walk reads the two parts in the order they lie.
//
// POSITIONS RATHER THAN POINTERS, and this is the whole of what the change buys. The compactor
// moves every live block, so a pointer taken before it runs points at whatever landed there
// afterwards, and a walk holding one cannot resume. A POSITION survives, because the pool itself
// never reallocates: `mSource` is sized once at construction and the compactor only slides data
// within it, so what a compaction invalidates is a block's OFFSET and never the array's base.
// Hoisting that base costs the same as holding a pointer, and `AMD_2` walks the same way for the
// same reason, keeping `p` and `pj` as indices and restoring them from `Pe [me]` and `Pe [e]`.
//
// THE OFFSETS ARE THEREFORE RE-READ FROM THE DESCRIPTORS at the top of each clique rather than
// carried, which is what makes them re-derivable at all. `mFree` is likewise read from the member
// and never cached across the walk.
// WRITES AT `mFree` RATHER THAN APPENDING, which is the storage change, and TESTS FOR ROOM PER
// ENTRY exactly as `AMD_2` does with `if (pfree >= iwlen) garbage_collection`. The walk resumes
// afterwards, which is what the two steps below it were for.
//
// INCIDENCE THEN ADJACENCY, which is `AMD_2`'s `for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++)`:
// the cliques of me on the first `elenme` passes and the supervariables on the last. It is now
// also the physical order of the run, so the walk reads the two parts in the order they lie.
//
// POSITIONS RATHER THAN POINTERS, and this is what makes resuming possible at all. The compactor
// moves every live block, so a pointer taken before it runs points at whatever landed there
// afterwards. A POSITION survives, because the pool itself never reallocates: `mSource` is sized
// once at construction and the compactor only slides data within it, so what a compaction
// invalidates is a block's OFFSET and never the array's base. Hoisting that base costs the same as
// holding a pointer. `AMD_2` walks the same way for the same reason, keeping `p` and `pj` as
// indices and restoring them from `Pe [me]` and `Pe [e]`.
//
// CURSORS RATHER THAN COUNTERS, for the same reason again. A compaction TRUNCATES the two lists
// being read, dropping the part already consumed, so after it a list starts at a new base with a
// new length and a counter into the old one means nothing. `p` walks the pivot's run and `pj` the
// clique being read, and both are re-read from their descriptors afterwards.
inline std::uint32_t QuotientGraphCompacted::reachableSetAmd(std::int32_t u) {
    std::int32_t* const pool = mSource.data();            // stable: the pool never reallocates
    std::uint32_t reached      = 0;                       // counted; the start moves under us
    std::size_t   cliqueStart  = mFree;                   // AMD_2's pme1, and it moves too

    // The pivot's own run, walked incidence part first. `remaining` is what is left of the part
    // being consumed, so a truncation is `sourcePtr = p` with the length set from it.
    std::size_t   p         = mRun[u].sourcePtr;
    std::uint32_t remaining = mRun[u].incidenceSize;
    bool          onCliques = true;

    for (;;) {
        std::size_t   pj;                                 // where the members being read start
        std::uint32_t ln;                                 // how many are left of them
        std::int32_t  c = NIL;                            // the clique, or NIL on the last pass
        if (onCliques) {
            if (remaining == 0) {                         // the incidence part is spent
                onCliques = false;
                remaining = mRun[u].adjacencySize;
                continue;
            }
            c  = pool[p++];
            --remaining;
            pj = mRun[c].sourcePtr;
            ln = mRun[c].adjacencySize;
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
            const std::int32_t nv = mWeight[v];
            if (nv <= 0) continue;

            if (mFree >= mSource.size()) {
                // TRUNCATE BOTH LISTS TO WHAT IS LEFT, which is `AMD_2`'s
                // `Pe [me] = p ; Len [me] -= knt1` and `Pe [e] = pj ; Len [e] = ln - knt2`. This
                // is the step the flipped run exists for: the consumed part is a PREFIX, so what
                // remains is still one contiguous block and the descriptor can simply be moved
                // forward. It is also what makes the space bound hold, the half-built clique
                // being a subset of exactly the prefixes dropped here.
                if (onCliques) {
                    mRun[u].sourcePtr     = p;
                    mRun[u].incidenceSize = remaining;
                } else {
                    mRun[u].sourcePtr     = pj;
                    mRun[u].incidenceSize = 0;
                    mRun[u].adjacencySize = ln;
                }
                // TRUNCATING A CLIQUE IS A CONTRACTION AND THE COUNTER HAS TO SEE IT. `c` is a
                // live clique being consumed into the new one; dropping its consumed prefix
                // shortens it, and `killClique` below will then subtract only what remained. Not
                // telling the counter here left the difference live forever, which is what the
                // `AMD3B pC differs` check caught on `JGD_Trefethen/Trefethen_2000`, +31 of 11091.
                // The peak is unaffected either way, being taken only at a birth, and by then the
                // clique is dead in full: the two subtractions sum to its original length.
                //
                // The pivot's own run needs no such line. `u` is a live vertex and the length
                // being shortened is A[u]'s, which this counter never held.
                if (c != NIL) {
                    mNumLiveCliqueMembers -= mRun[c].adjacencySize - ln;
                    mRun[c].sourcePtr      = pj;
                    mRun[c].adjacencySize  = ln;
                }
                compact(cliqueStart);
                // AND RESUME FROM THE DESCRIPTORS, which is `AMD_2`'s `pj = Pe [e] ; p = Pe [me]`.
                p  = mRun[u].sourcePtr;
                pj = (c != NIL) ? mRun[c].sourcePtr : p;
                if (!onCliques) p = pj;
            }

            mWeight[v] = -nv;
            pool[mFree++] = v;
            ++reached;
        }
        if (!onCliques) p = pj;
    }
    return reached;
}

// WRITES AT `mFree` RATHER THAN APPENDING, which is the storage change. The caller has already
// made room for a whole reach, so the cursor cannot run out inside the walk; and the reach lands
// exactly where the clique is to live, so there is no copy from a scratch into place. Production
// appends to its arena and needs a capacity check per pivot for the same guarantee.
inline std::uint32_t QuotientGraphCompacted::reachableSetMmd(std::int32_t u) {
    ++mTag;
    // THE PIVOT IS ALREADY NEGATED, by `beginEliminationMmd`, which is where the amd branch has
    // always done it and where the merge moved this branch's. It was here and in the in-place walk
    // before, and leaving it in both places negates twice, which reads as a positive weight and so
    // as a vertex not yet taken. Caught on a 3 by 3 grid by `cliqueCountBalances`, in the first
    // build of the merged class.
    // THE TRUNCATION BELOW IS WRITTEN FOR THE REVERSED WALK, which is the only one this file runs.
    // Reversed, the consumed part of I[u] is a suffix and the survivors are the prefix the run
    // already points at, so a truncation is a length; forward it would be the other way round and
    // `sourcePtr` would have to move with it. The flag is the driver's, so the assumption is
    // checked rather than assumed.
    const bool reverse = mReverseIncidence;
    assert(reverse && "reachableSetMmd's truncation assumes the reversed incidence walk");

    // POSITIONS AND CURSORS, and this is what makes resuming possible at all. The compactor moves
    // every live block, so a pointer taken before it runs points at whatever landed there
    // afterwards; a POSITION survives, the pool being sized once at construction so that a
    // compaction invalidates a block's OFFSET and never the array's base. And a cursor counting
    // what is LEFT survives a truncation, where a counter into the original length would not.
    // `AMD_2` walks the same way for the same reason, keeping `p` and `pj` as indices and
    // restoring them from `Pe [me]` and `Pe [e]`.
    std::int32_t* const pool = mSource.data();          // stable: the pool never reallocates
    std::uint32_t       reached     = 0;                // counted; the start moves under us
    std::size_t         cliqueStart = mFree;            // AMD_2's pme1, and it moves too

    // A[u] FORWARD, so the consumed part is a PREFIX and what remains is one contiguous block
    // running to the end of the run. `p` is where the unconsumed part starts and `remaining` how
    // much of it is left, which is the pair a truncation writes back.
    std::size_t   p         = mRun[u].sourcePtr;
    std::uint32_t remaining = mRun[u].adjacencySize;
    while (remaining-- > 0) {
        const std::int32_t v  = pool[p++];
        const std::int32_t nv = mWeight[v];
        if (nv <= 0 || (mHasNumbered && mMark[v] == GONE)) continue;

        if (mFree >= mSource.size()) {
            // TRUNCATE TO WHAT IS LEFT OF A[u], which is `AMD_2`'s `Pe [me] = p ; Len [me] -=
            // knt1`. The consumed part is a prefix, so the survivors are still one block: A[u]'s
            // tail followed by the whole of I[u], which is what `sourcePtr` and the two lengths
            // now describe. It is also what makes the space bound hold, the half-built clique
            // being a subset of exactly the prefixes dropped here.
            mRun[u].sourcePtr     = p;
            mRun[u].adjacencySize = remaining;
            compact(cliqueStart);
            p = mRun[u].sourcePtr;                      // AMD_2's `p = Pe [me]`
        }

        mWeight[v] = -nv;
        pool[mFree++] = v;
        ++reached;
    }

    // I[u] REVERSED, which is `mmdelm` pushing cliques onto a stack and popping them. The
    // consumed part is therefore a SUFFIX, so what remains is the PREFIX the run already points at
    // and a truncation is a length rather than a move. `left` is that length.
    //
    // THE BASE IS RE-READ RATHER THAN HOISTED, for the same reason the walk uses positions: a
    // compaction moves the run. `sourcePtr + adjacencySize` is where I[u] begins under every state
    // this loop can be in, including after the first truncation below has set the length to zero.
    std::size_t   incidenceBase = mRun[u].sourcePtr + mRun[u].adjacencySize;
    const std::uint32_t incidenceLen = mRun[u].incidenceSize;
    std::uint32_t left = incidenceLen;
    while (left-- > 0) {
        const std::int32_t c  = pool[incidenceBase + (reverse ? left : incidenceLen - 1 - left)];
        std::size_t        pj = mRun[c].sourcePtr;
        std::uint32_t      ln = mRun[c].adjacencySize;
        while (ln-- > 0) {
            const std::int32_t v  = pool[pj++];
            const std::int32_t nv = mWeight[v];
            if (nv <= 0) continue;

            if (mFree >= mSource.size()) {
                // A[u] IS SPENT BY NOW, so the run becomes I[u]'s surviving prefix alone and the
                // consumed adjacency goes with the rest of the dead space. Dropping it is not
                // tidiness: the space bound rests on the half-built clique being covered by what
                // the walk has already released.
                mRun[u].sourcePtr     = incidenceBase;
                mRun[u].adjacencySize = 0;
                mRun[u].incidenceSize = left;

                // TRUNCATING A CLIQUE IS A CONTRACTION AND THE COUNTER HAS TO SEE IT. `c` is a
                // live clique being consumed into the new one; dropping its consumed prefix
                // shortens it, and the `killClique` in `beginElimination` would then subtract only
                // what remained, leaving the difference live for the rest of the run. This is the
                // second of the two defects `Amd3B` introduced with its own compactor.
                mNumLiveCliqueMembers -= mRun[c].adjacencySize - ln;
                mRun[c].sourcePtr      = pj;
                mRun[c].adjacencySize  = ln;

                compact(cliqueStart);
                // AND RESUME FROM THE DESCRIPTORS, which is `AMD_2`'s `pj = Pe [e] ; p = Pe [me]`.
                incidenceBase = mRun[u].sourcePtr + mRun[u].adjacencySize;
                pj            = mRun[c].sourcePtr;
            }

            mWeight[v] = -nv;
            pool[mFree++] = v;
            ++reached;
        }
    }
    return reached;
}

// CONSTRUCTS THE CLIQUE IN THE PIVOT'S OWN RUN, which is `AMD_2`'s `if (elenme == 0)` branch and
// the common case rather than a corner: a pivot with no cliques has a reach that is a SUBSET of
// A[pivot], so it fits where A[pivot] already is and the pool is not touched at all. Measured on
// grids, 62 to 68 percent of eliminations qualify, and the share rises with n.
//
// IT IS AN IN-PLACE COMPACTION and safe for the reason every such loop here is: the write cursor
// starts at the read cursor and only ever falls behind it, since a vertex is written only when it
// was read. `AMD_2` spells it `Iw [++pme2] = i` with `pme2 = pme1 - 1`.
//
// WHY IT MATTERS TWICE. It keeps two thirds of cliques out of the pool, so the cursor advances far
// more slowly and the compactor runs far less; and it leaves the clique exactly where the pivot's
// adjacency was, which is where the vertices that will read it next are. This file existed for a
// year without it, and the figures it produced were the price of a DIFFERENT layout.
inline std::uint32_t QuotientGraphCompacted::reachableSetInPlaceAmd(std::int32_t u) {
    std::int32_t*       run   = mSource.data() + mRun[u].sourcePtr;
    const std::uint32_t count = mRun[u].adjacencySize;
    std::uint32_t       kept  = 0;
    for (std::uint32_t k = 0; k < count; ++k) {
        const std::int32_t v  = run[k];
        const std::int32_t nv = mWeight[v];
        if (nv > 0) { mWeight[v] = -nv; run[kept++] = v; }
    }
    return kept;
}

// CONSTRUCTS THE CLIQUE IN THE PIVOT'S OWN RUN, which is `AMD_2`'s `if (elenme == 0)` branch and
// the common case rather than a corner: a pivot with no cliques has a reach that is a SUBSET of
// A[pivot], so it fits where A[pivot] already is and the pool is not touched. Measured on grids,
// 62 to 68 percent of mmd's eliminations qualify, and the share rises with n.
//
// IT IS AN IN-PLACE COMPACTION and safe for the reason every such loop here is: the write cursor
// starts at the read cursor and only ever falls behind it, a vertex being written only when it was
// read. `AMD_2` spells it `Iw [++pme2] = i` with `pme2 = pme1 - 1`.
//
// THE GONE TEST SURVIVES HERE, unlike in Amd3B, and for the reason it survives everywhere on this
// branch: `number()` leaves a prepass vertex at weight one and in every neighbour's adjacency, so
// a positive weight does not mean live in an ADJACENCY list. This walk reads nothing else.
inline std::uint32_t QuotientGraphCompacted::reachableSetInPlaceMmd(std::int32_t u) {
    // THE PIVOT IS ALREADY NEGATED, by `beginEliminationMmd`. It has to be negated somewhere and
    // not merely to keep u out of its own list, which has no self loop: `massEliminateMmd` restores
    // the sign of every member AND of the pivot, so a pivot left positive comes out of that restore
    // negative and `orderAscending` then reads a supervariable of negative size. That was found by
    // ASan on a 3 by 3 grid when this walk was written; negating twice, which is what hoisting it
    // without removing it here would do, is the same bug from the other side.
    std::int32_t*       run   = mSource.data() + mRun[u].sourcePtr;
    const std::uint32_t count = mRun[u].adjacencySize;
    std::uint32_t       kept  = 0;
    for (std::uint32_t k = 0; k < count; ++k) {
        const std::int32_t v  = run[k];
        const std::int32_t nv = mWeight[v];
        if (nv > 0 && !(mHasNumbered && mMark[v] == GONE)) { mWeight[v] = -nv; run[kept++] = v; }
    }
    return kept;
}

inline std::uint32_t QuotientGraphCompacted::reachableWeight(std::int32_t u) {
    ++mTag;
    std::uint32_t reached = 0;   // a sum over DISTINCT vertices, so bounded by n; see the header
    mMark[u] = mTag;
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
        const std::int32_t* members     = mSource.data() + mRun[c].sourcePtr;
        const std::uint32_t membersSize = mRun[c].adjacencySize;
        for (std::uint32_t k = 0; k < membersSize; ++k) {
            const std::int32_t v = members[k];
            if (mMark[v] < mTag) { mMark[v] = mTag; reached += static_cast<std::uint32_t>(mWeight[v]); }
        }
    }
    return reached;
}

#ifndef NDEBUG
inline bool QuotientGraphCompacted::cliqueCountBalances() const {
    std::size_t live = 0;
    for (std::int32_t c : mCliqueOwners) live += mRun[c].adjacencySize;
    return live == mNumLiveCliqueMembers;
}
#endif

inline std::vector<std::int32_t> QuotientGraphCompacted::order(
        const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order;
    order.reserve(size());
    for (std::int32_t pivot : pivots)
        for (std::int32_t u = pivot; u != NIL; u = mSuperNext[u]) order.push_back(u);
    return order;
}

inline std::vector<std::int32_t> QuotientGraphCompacted::orderAscending(
        const std::vector<std::int32_t>& pivots) const {
    const std::size_t n = size();
    std::vector<std::int32_t> order(n);
    std::vector<std::int32_t> slot(n, 0);

    std::size_t pos = 0;
    for (std::int32_t pivot : pivots) {
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

// AMD_2'S COMPACTION, ported. The pool fills because absorbed cliques and shortened runs
// leave dead space behind, and nothing reclaims it until this does. Everything live is compacted to
// the front and every descriptor is repointed.
//
// THE TRICK IS AMD_2'S AND IT IS WHAT MAKES THIS ONE PASS RATHER THAN A SORT. A block's owner
// cannot be recovered by looking at the pool, so before compacting, the FIRST ENTRY of every live
// block is replaced by FLIP(owner) and the entry it displaced is parked in that owner's own
// `sourcePtr`, which is about to be overwritten anyway. The pool can then be scanned from the
// front: a negative entry is a block head, its owner is FLIP of it, and its length is the owner's
// descriptor. Amd.cpp does exactly this with `Pe[j] = Iw[pn]` and `Iw[pn] = FLIP(j)`.
//
// A LIVE BLOCK IS ONE WITH A NONZERO LENGTH, and that covers both kinds without a test for which:
// a live vertex has I[e] then A[e], a dead pivot has C[e], and anything absorbed or merged had its
// lengths zeroed when it died.
//
// AND IT CARRIES THE HALF-BUILT CLIQUE, which is why it takes its start by reference. When the
// compactor runs during a clique build the partial sits at the top of the pool with no descriptor
// of its own, so the sweep stops just below it, `AMD_2`'s `pend = pme1 - 1`, and the partial is
// then moved down explicitly and its start rewritten. Called with `cliqueStart == mFree` outside a
// build, where the partial is empty and both halves degenerate to the plain sweep.
inline void QuotientGraphCompacted::compact(std::size_t& cliqueStart) {
    ++mCompactions;
    const std::int32_t n = static_cast<std::int32_t>(mRun.size());

    for (std::int32_t e = 0; e < n; ++e) {
        const std::size_t len = static_cast<std::size_t>(mRun[e].adjacencySize)
                              + mRun[e].incidenceSize;
        if (len == 0) continue;
        const std::size_t at = mRun[e].sourcePtr;
        mRun[e].sourcePtr = static_cast<std::size_t>(
                                static_cast<std::uint32_t>(mSource[at]));   // park the displaced
        mSource[at] = FLIPPED - e;                                          // and mark the head
    }

    std::size_t src = 0, dst = 0;
    while (src < cliqueStart) {                             // AMD_2's `while (psrc <= pend)`
        const std::int32_t head = mSource[src++];
        if (head > FLIPPED) continue;                       // dead space, skipped
        const std::int32_t e = FLIPPED - head;
        mSource[dst] = static_cast<std::int32_t>(
                           static_cast<std::uint32_t>(mRun[e].sourcePtr));  // restore it
        mRun[e].sourcePtr = dst;
        ++dst;
        const std::size_t rest = static_cast<std::size_t>(mRun[e].adjacencySize)
                               + mRun[e].incidenceSize - 1;
        for (std::size_t k = 0; k < rest; ++k) mSource[dst + k] = mSource[src + k];
        src += rest;
        dst += rest;
    }

    // THE PARTIAL CLIQUE FOLLOWS, moved as a block. It has no descriptor to repoint, so its start
    // is the caller's variable and is returned through it.
    const std::size_t moved = dst;
    for (std::size_t k = cliqueStart; k < mFree; ++k) mSource[dst++] = mSource[k];
    cliqueStart = moved;
    mFree = dst;
}

// I[pivot] COPIED BEFORE ANY WALK TOUCHES IT. A walk truncates that list as it consumes it, so
// after a mid-walk compaction the run is short or empty and the cliques it named could never be
// killed from it: each would keep a non-zero length, stay a live block for the compactor to copy,
// and never be subtracted from the live count. That defect existed in `Amd3B` for a day, moved no
// permutation, and was found by a benchmark comparing peaks across drivers.
inline void QuotientGraphCompacted::captureAbsorbed(const std::int32_t* incidence,
                                                    std::uint32_t count) {
    mAbsorbed.clear();
    for (std::uint32_t i = 0; i < count; ++i) mAbsorbed.push_back(incidence[i]);
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
    mRun[pivot].sourcePtr     = cliqueStart;
    mRun[pivot].adjacencySize = cliqueLen;

    // A CLIQUE IS BORN HERE AND NOWHERE ELSE, which is what makes the peak countable, and the
    // maximum is taken here alone since nothing else raises the total.
    mNumLiveCliqueMembers += cliqueLen;
    mNumPeakCliqueMembers  = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);
#ifndef NDEBUG
    mCliqueOwners.push_back(pivot);
#endif

    // NO STAMPING PASS AND NO TAG. Membership was written by the walk, in the sign of the weight,
    // so this only sums.
    const std::int32_t* reached = mSource.data() + cliqueStart;
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t k = 0; k < cliqueLen; ++k)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[k]]);
    mCliqueWeight = cliqueWeight;
}


} // namespace Oblio
