#include "oblio/Amd3B.h"

#include "oblio/Types.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// Amd3B.cpp - Amd3 on AMD_2'S CLIQUE STORAGE: one pool with a free cursor and a garbage collection.
//
// WHAT IT IS FOR, AND IT IS TWO THINGS. Both are permanent; this file is not an experiment awaiting
// a verdict and its old stop condition is withdrawn. It is the amd counterpart of Mmd3B, which says
// the same of itself.
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
// different existing pass serving the same purpose. The detection stamp landed in `Amd3` and in
// `Amd2`, both of which already carried the tagged `w`. `Amd1` needed nothing, having no
// supervariable detection at all. SO THIS FILE'S TIME COLUMN IS NOW THE STORAGE PRICE ALONE, as
// Mmd3B's already was, and it reads about 6 to 8 percent below AMD3 up to 256 a side and 15 percent
// below it from 400 up, reproduced across two runs.
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
// SECOND, IT IS THE PREDICTABLE-SPACE VERSION OF AMD3. From a conversation with Alex Pothen: given
// a machine you know whether A fits, but you cannot know whether L fits, nnz(L) depending on the
// ordering being computed. So a method that stays within `O(n + m)` carries a guarantee no amount
// of speed substitutes for: IF THE INPUT FITS, THE ANSWER IS REACHABLE. The garbage collection
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
// ITS OBLIGATION. It must return Amd3's permutation exactly, which is `AMD_2`'s raw order, so its
// nnz(L) column in benchmarks/ordering must equal AMDraw's on every row and its fill column carries
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
//      collection invalidates a block's OFFSET and never the array's base; a position therefore
//      survives one and a pointer does not. Same cost, and it is how `AMD_2` walks.
//   3  THE COLLECTOR RUNS PER ENTRY and carries the half-built clique, so `beginElimination` no
//      longer reserves room for a worst-case reach of n before the walk. That reservation was the
//      last divergence and it was one-directional: we collected where `AMD_2` would not.
//
// TWO THINGS WENT WITH THEM, both flags this layout cannot serve two values of and neither ever
// set here: the list order and the reverse incidence walk. See the notes where each setter was.
// One thing was added: a slide in `absorb`, the adjacency having to follow the incidence part down
// when that part shrinks. That is a pass `AMD_2` does not make at all, and the flip priced it
// rather than caused it. See the note there.
//
// VERIFIED BY VARYING THE HEADROOM, which is a better check than the digest alone. At the shipped
// reserve the collector never runs on grids, so the mid-walk path would go untested; cut to
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
// THE COPY IS MECHANICAL EXCEPT IN THE STORAGE LAYER. QuotientGraphCompacted and
// BucketsCompacted are
// include/oblio/QuotientGraph.h and src/QuotientGraph.cpp with the whole-line comments removed and
// the members no amd driver calls pruned; the driver is src/Amd3.cpp. Every design note for these
// types lives in those files and is authoritative there. The pool, the collector and the run order
// are hand written and have no counterpart to be regenerated from, so tmp/make_amd3b.py produces
// the rest and the storage is applied on top. Only the DIFFERENCES are commented here.

namespace Oblio {

std::size_t gAmd3BCompactions = 0;   // read by tmp/ probes; see the note at the return below
extern std::size_t gPeakCliqueMembers;   // defined in src/QuotientGraph.cpp; see it there

namespace {

class BucketsCompacted {
public:
    static constexpr std::int32_t UNFILED    = 0;
    static constexpr std::int32_t OUTMATCHED = -2147483647 - 1;   // INT32_MIN; no degree reaches it

    explicit BucketsCompacted(std::size_t size)
        : mHead(size + 1, NIL), mNext(size, NIL), mPrev(size, UNFILED) {}

    void file(std::uint32_t degree, std::int32_t u) {
        mNext[u] = mHead[degree];
        mPrev[u] = -static_cast<std::int32_t>(degree) - 1;   // head of `degree`; see the encoding  // at the head, and this is its bucket
        if (mHead[degree] != NIL) mPrev[mHead[degree]] = u + 1;
        mHead[degree] = u;
    }

    void unfile(std::int32_t u) {
        const std::int32_t prev = mPrev[u];
        if (prev == UNFILED || prev == OUTMATCHED) return;   // not in a list; see the encoding
        if (prev > 0) mNext[prev - 1]   = mNext[u];
        else          mHead[-prev - 1] = mNext[u];           // u headed bucket -prev - 1
        if (mNext[u] != NIL) mPrev[mNext[u]] = prev;
        mNext[u] = NIL;
        mPrev[u] = UNFILED;
    }

    void outmatch(std::int32_t u)         { unfile(u); mPrev[u] = OUTMATCHED; }
    void restore(std::int32_t u)          { if (mPrev[u] == OUTMATCHED) mPrev[u] = UNFILED; }
    bool outmatched(std::int32_t u) const { return mPrev[u] == OUTMATCHED; }

    void refile(std::vector<std::uint32_t>& degrees, std::int32_t u, std::uint32_t newDegree) {
        unfile(u);
        degrees[u] = newDegree;
        file(newDegree, u);
    }

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
};

struct TaggedScanCompacted {
    BucketsCompacted*                          buckets;
    std::vector<std::int32_t>&        w;             // per clique, Amd.cpp's tagged W
    const std::vector<std::uint32_t>& degree;
    std::vector<std::int32_t>&        touchedCliques;// the cliques this step reached, once each
    std::int32_t                      wflg;          // the tag for this elimination
    std::int32_t                      modulus;       // the driver's bucket count, n + 1
};

class QuotientGraphCompacted {
public:
    // NO cliqueMarks ARGUMENT. Production takes one to size its mark array at n or 2n; this file
    // has no mark array, so there was nothing for the flag to decide and it survived only as a
    // member the constructor had to cast to void. Removed 2026-08-17.
    QuotientGraphCompacted(const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const           { return mRun.size(); }

    // A block head in the pool, during garbage collection only: FLIPPED - e for owner e, so every
    // head is at most FLIPPED and every real entry, a vertex id, is above it.
    //
    // GONE IS GONE, and this is the only reserved value left. It marked an eliminated vertex in the
    // mark array; with no mark array there is nothing for a reserved tag to sit above, and the
    // question it answered is now a zero weight.
    static constexpr std::int32_t FLIPPED = -2;

    // ZERO WEIGHT IS THE DEAD STATE, which is Amd.cpp's `Nv [i] == 0` and removes the last
    // mMark read from the hash detection loops. The three ways a vertex leaves the graph all end
    // in a zero weight: `merge` zeroes the absorbed one, `massEliminate` zeroes the merged ones.
    //
    // A PIVOT IS THE ONE CASE THIS DOES NOT COVER, and it does not need to: an eliminated pivot
    // keeps its supervariable weight, but no list this driver walks can still contain it. If w is
    // in A[u] then u is in reach(w), so u was in C[w] and its A[u] was pruned when w was
    // eliminated, with w dropped there. Clique member lists are trimmed the same way. The prune's
    // own `nv <= 0` test has relied on this since the sign fold and 584 permutations agree.
    bool eliminated(std::int32_t u) const { return mWeight[u] == 0; }


    // NO MARK ARRAY, AND THAT IS FINISHED RATHER THAN PENDING. Two comments stood here until
    // 2026-08-17 saying the vertex half survived for the hash exact test and that retiring mMark
    // outright was "the next thing to try". Both were written before the stamp moved into `w` a
    // few hundred lines below, which did exactly that; the member declaration had said so all
    // along. Recorded because a comment describing an intermediate state reads as current, and a
    // reader would have concluded this file still has an array it does not.


    // A clique has no incidence part, so its block starts at the run and `adjacency` would give
    // the same address. Spelled without the offset because a clique is one list rather than two.
    const std::int32_t* clique(std::int32_t c) const {
        return mSource.data() + mRun[c].sourcePtr;
    }
    std::uint32_t cliqueSize(std::int32_t c) const { return mRun[c].adjacencySize; }

    // WRITABLE, for the restore pass alone, which trims the clique as it walks it. See
    // `trimClique` for why that pass and not a pass of its own.
    std::int32_t* clique(std::int32_t c) { return mSource.data() + mRun[c].sourcePtr; }
    void trimClique(std::int32_t pivot, std::uint32_t kept);

    // THE ONE PLACE A CLIQUE DIES, so that the counter below can see every death. Three causes:
    // absorbed into the new clique, absorbed aggressively once its external degree reaches zero,
    // and merged away with its owner. All three used to store the zero where they stood.
    void killClique(std::int32_t c);

    // PEAK LIVE CLIQUE MEMBERS, mirroring QuotientGraph's. A MEMBER is a vertex in a live clique
    // at this instant, where an entry is a pool slot; the two differ here because this layout
    // reclaims and the flat one does not.
    //
    // IT IS A PROPERTY OF THE ALGORITHM, NOT OF THE LAYOUT, which is the point of having it in
    // both classes. `Amd3` and `Amd3B` compute the same permutation, so they form the same cliques
    // and merge the same vertices at the same moments, and these two numbers MUST be equal. The
    // digest says the outputs agree; this says the work behind them agreed too, which is a
    // stronger statement and one nothing else in the suite makes. Checked in tests/test_order.cpp.
    std::size_t numPeakCliqueMembers() const { return mNumPeakCliqueMembers; }

    std::uint32_t cliqueWeight() const { return mCliqueWeight; }

    // INCIDENCE FIRST, ADJACENCY BEHIND IT, which is `AMD_2`'s order and the reverse of
    // production's. See the note on `VertexRun` for why the order is load bearing here and a
    // matter of indifference there.
    const std::int32_t* incidence(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr;
    }
    std::uint32_t incidenceSize(std::int32_t u) const { return mRun[u].incidenceSize; }

    const std::int32_t* adjacency(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr + mRun[u].incidenceSize;
    }
    std::uint32_t adjacencySize(std::int32_t u) const { return mRun[u].adjacencySize; }

    // The MAGNITUDE, for the driver, which does not want to know about the sign. The hot loops
    // inside this class read mWeight directly and test the sign, which is the whole point.
    std::uint32_t weight(std::int32_t u) const {
        const std::int32_t w = mWeight[u];
        return static_cast<std::uint32_t>(w < 0 ? -w : w);
    }
    // Undo the negation ON ONE VERTEX and hand back its weight, which is Amd.cpp's
    // `nvi = -Nv [i] ; ... Nv [i] = nvi` under RESTORE DEGREE LISTS. It is a single vertex rather
    // than a pass because the driver's bound pass already walks C[pivot] and reads exactly this
    // weight: doing it there costs a store on a line already loaded, where a pass of its own cost
    // a whole extra traversal per pivot. Amd.cpp pays neither, and this is why.
    std::uint32_t restoreWeight(std::int32_t u) {
        std::int32_t w = mWeight[u];
        if (w < 0) { w = -w; mWeight[u] = w; }
        return static_cast<std::uint32_t>(w);
    }
    // The pivot's own, which no walk over C[pivot] reaches: it is not a member of its own clique.
    void restorePivotWeight(std::int32_t pivot) {
        if (mWeight[pivot] < 0) mWeight[pivot] = -mWeight[pivot];
    }

    std::uint32_t reachableSet(std::int32_t u);          // writes at mFree, returns the length
    std::uint32_t reachableSetInPlace(std::int32_t u);   // rewrites A[u] as C[u]; see its note




    const std::vector<std::int32_t>& eliminate(std::int32_t pivot, TaggedScanCompacted& scan);

    void merge(std::int32_t u, std::int32_t v);

    // SET u ASIDE, taking it out of the elimination without numbering it. `AMD_2`'s dense-row
    // rule; see the shared class for the whole of it. Here a zero weight is the dead state
    // outright, so this is the one store and there is no mark to write.
    void setAside(std::int32_t u) { mWeight[u] = 0; }


    // NO setReverseIncidence. Walking I[u] from the back is genmmd's convention and cannot
    // coexist with front truncation, which drops what the walk has already consumed. Never set
    // here, so nothing is lost; the same reasoning retired the list-order flag above.

    // NO setVendoredListOrder. Production carries the flag because its run holds the adjacency
    // first and can serve either convention. This layout serves ONE: with the incidence part at
    // the front there is no way to append the new clique behind it without either a free slot or a
    // shift, which is exactly why `AMD_2` inserts by rotation. The flag was always true here, so
    // nothing is lost, and leaving a false branch that the layout cannot express would be worse
    // than not having one.

    void setLateMassElimination(bool on) { mLateMassElimination = on; }

    const std::vector<std::int32_t>& massEliminate(std::int32_t pivot);

    void absorb(const std::vector<std::int32_t>& cliques,
                const std::int32_t* vertices, std::uint32_t vertexCount);

    std::vector<std::int32_t> order(const std::vector<std::int32_t>& pivots) const;


private:
    void beginElimination(std::int32_t pivot);

    const std::vector<std::int32_t>& finishElimination(std::int32_t pivot);

    std::vector<std::int32_t>  mSource;         // every I[u] then A[u], run after run

    // `AMD_2`'S DESCRIPTOR, WITH ITS TWO LENGTHS AND ITS ORDER. `Pe [i]` is the start, `Elen [i]`
    // the number of elements at the FRONT, `Len [i]` the total, so the adjacency length is the
    // difference. We carry both lengths instead, which is the same descriptor with the subtraction
    // taken out.
    //
    // THE ORDER IS LOAD BEARING HERE AND IS NOT IN PRODUCTION, which lays the run out the other
    // way round. `AMD_2` collects garbage in the MIDDLE of building an element and resumes, which
    // it can do only because the part it has already consumed is a PREFIX: `Pe [me] = p ;
    // Len [me] -= knt1` advances the start and shortens the length, and what is left is still one
    // contiguous block. With the adjacency in front, the consumed part of I[u] sits in the middle
    // of the run and there is no start to advance. So the order is what makes the collector
    // portable, and production, which never truncates anything, is free either way.
    // THE START MOVES, IN THREE PLACES, and none of them is construction: the collector repoints
    // every block, a truncation mid-build advances it past what has been consumed, and `absorb`
    // advances it past the cliques it drops from the front. A position held across any of those
    // is stale, which is why every walk re-reads it from here.
    struct VertexRun {
        std::size_t   sourcePtr;       // where u's run starts in mSource; moves, see below
        std::uint32_t incidenceSize;   // I[u]'s length, from the run's start; AMD_2's Elen
        std::uint32_t adjacencySize;   // A[u]'s length, immediately behind I[u]
    };
    std::vector<VertexRun> mRun;

    // ONE POOL, WHICH IS THE WHOLE POINT OF THIS FILE. Production keeps C[c] in a separate
    // append-only arena; `AMD_2` keeps element lists and variable lists in one `Iw` with elbow
    // room, builds a new element at a free cursor, leaves absorbed space dead, and COMPACTS when
    // the cursor reaches the end. `mFree` is that cursor and `garbageCollect` is that compaction.
    //
    // A block's descriptor is `mRun[e]` for both kinds, as in production after the descriptor
    // fold: a live vertex's run is A[e] then I[e], and a dead pivot's is C[e]. With one pool that
    // stops being two meanings of one field and becomes one.
    std::size_t mFree = 0;                    // AMD_2's pfree
    // Whether the clique now being eliminated was built in the pivot's own run. Read by
    // massEliminate, which can only give space back to the cursor in the other case.
    bool        mBuiltInPlace = false;

    std::size_t mNumLiveCliqueMembers = 0;   // see numPeakCliqueMembers
    std::size_t mNumPeakCliqueMembers = 0;
    std::size_t mCompactions = 0;             // AMD_2's Info[AMD_NCMPA]
public:
    std::size_t compactions() const { return mCompactions; }
private:

    void garbageCollect(std::size_t& cliqueStart);

    std::vector<std::int32_t>  mSuperNext;
    std::vector<std::int32_t>  mSuperLast;
    // AMD_2'S `Nv`, AND THE SIGN IS THE MEMBERSHIP MARK. This is the change this file is now
    // about. Amd.cpp's inner loops read ONE array per element:
    //
    //     mWeight[v] >  0    live, principal, not yet taken into the new clique; the weight
    //     mWeight[v] <  0    taken into the new clique this step; the weight is -mWeight[v]
    //     mWeight[v] == 0    absorbed, by a hash merge or by mass elimination
    //
    // so `nvi = Nv[i]; if (nvi > 0)` answers "is it dead", "is it already inside the new clique"
    // and "what does it weigh" off a single load. Production asks the first two of `mMark[v]` and
    // the third of `mWeight[v]`, which is TWO scattered loads per element in the two hottest loops
    // in the ordering, and `mMark` is an array AMD_2 does not have at all.
    //
    // The negation is applied as C[pivot] is built and undone in the bound pass, which is exactly
    // where Amd.cpp does it: `Nv[i] = -nvi` in CONSTRUCT NEW ELEMENT and `Nv[i] = nvi` under
    // RESTORE DEGREE LISTS. Between those two points every member of C[pivot] reads negative, and
    // `weight()` returns the magnitude so the driver does not have to care.
    std::vector<std::int32_t> mWeight;
    std::uint32_t              mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away


    bool mLateMassElimination = false;

    // NO MARK ARRAY. Every question it used to answer is now read off the weight or the tagged
    // W: liveness and clique membership from the sign, absorbed from the zero, and supervariable
    // detection's stamp from `w`, which is what AMD_2 does with `W [Iw [p]] = wflg`. That is the
    // last array this file owned that AMD_2 does not, and it was 2n int32 wide.
};

QuotientGraphCompacted::QuotientGraphCompacted(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mRun(colPtr.empty() ? 0 : colPtr.size() - 1) {
    const std::int32_t size = static_cast<std::int32_t>(mRun.size());

    // ELBOW ROOM, `AMD_2`'S EXACTLY. It sizes `slen = nzaat + nzaat/5 + 7n` and hands
    // `iwlen = slen - 6n` to the main loop, so the pool is `nzaat + nzaat/5 + n` where `nzaat` is
    // the OFF-DIAGONAL entry count, `sum (Len [0..n-1])`. Ours holds the same lists, so the figure
    // transfers, and it has to be the same figure rather than the same shape: the headroom is what
    // decides how often the collector runs, so a differential against `AMD_2` on compaction counts
    // measures the headroom unless the two agree. Until 2026-08-18 ours was computed from `nnz`
    // WITH the diagonal and so ran about 1.2n large.
    //
    // Reserved from an upper bound and then sized down, rather than counted in a pass of its own:
    // `nnz >= nzaat`, so the reserve cannot reallocate, and `nzaat` is known once the runs are
    // laid out. The vector is SIZED rather than left at capacity because the space past the runs
    // is written by cliques and read by nothing until it is, and a cursor into a vector's unused
    // capacity would be undefined.
    const std::size_t nnz = colPtr.empty() ? 0 : colPtr.back();
    mSource.reserve(nnz + nnz / 5 + mRun.size());
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

    mSuperNext.assign(size, NIL);
    mSuperLast.resize(size);
    mWeight.assign(size, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

// CONSTRUCTS THE CLIQUE IN THE PIVOT'S OWN RUN, which is `AMD_2`'s `if (elenme == 0)` branch and
// the common case rather than a corner: a pivot with no elements has a reach that is a SUBSET of
// A[pivot], so it fits where A[pivot] already is and the pool is not touched at all. Measured on
// grids, 62 to 68 percent of eliminations qualify, and the share rises with n.
//
// IT IS AN IN-PLACE COMPACTION and safe for the reason every such loop here is: the write cursor
// starts at the read cursor and only ever falls behind it, since a vertex is written only when it
// was read. `AMD_2` spells it `Iw [++pme2] = i` with `pme2 = pme1 - 1`.
//
// WHY IT MATTERS TWICE. It keeps two thirds of cliques out of the pool, so the cursor advances far
// more slowly and the collector runs far less; and it leaves the clique exactly where the pivot's
// adjacency was, which is where the vertices that will read it next are. This file existed for a
// year without it, and the figures it produced were the price of a DIFFERENT layout.
std::uint32_t QuotientGraphCompacted::reachableSetInPlace(std::int32_t u) {
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

// WRITES AT `mFree` RATHER THAN APPENDING, which is the storage change. The caller still ensures
// room for a whole reach before the walk starts, so the cursor cannot run past the pool inside it;
// that test is the last divergence from `AMD_2`, which tests per entry instead, and closing it is
// what the positions below are for.
//
// INCIDENCE THEN ADJACENCY, which is `AMD_2`'s `for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++)`:
// the elements of me on the first `elenme` passes and the supervariables on the last. It is now
// also the physical order of the run, so the walk reads the two parts in the order they lie.
//
// POSITIONS RATHER THAN POINTERS, and this is the whole of what the change buys. The collector
// moves every live block, so a pointer taken before it runs points at whatever landed there
// afterwards, and a walk holding one cannot resume. A POSITION survives, because the pool itself
// never reallocates: `mSource` is sized once at construction and the collector only slides data
// within it, so what a collection invalidates is a block's OFFSET and never the array's base.
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
// the elements of me on the first `elenme` passes and the supervariables on the last. It is now
// also the physical order of the run, so the walk reads the two parts in the order they lie.
//
// POSITIONS RATHER THAN POINTERS, and this is what makes resuming possible at all. The collector
// moves every live block, so a pointer taken before it runs points at whatever landed there
// afterwards. A POSITION survives, because the pool itself never reallocates: `mSource` is sized
// once at construction and the collector only slides data within it, so what a collection
// invalidates is a block's OFFSET and never the array's base. Hoisting that base costs the same as
// holding a pointer. `AMD_2` walks the same way for the same reason, keeping `p` and `pj` as
// indices and restoring them from `Pe [me]` and `Pe [e]`.
//
// CURSORS RATHER THAN COUNTERS, for the same reason again. A collection TRUNCATES the two lists
// being read, dropping the part already consumed, so after it a list starts at a new base with a
// new length and a counter into the old one means nothing. `p` walks the pivot's run and `pj` the
// clique being read, and both are re-read from their descriptors afterwards.
std::uint32_t QuotientGraphCompacted::reachableSet(std::int32_t u) {
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
                if (c != NIL) {
                    mRun[c].sourcePtr     = pj;
                    mRun[c].adjacencySize = ln;
                }
                garbageCollect(cliqueStart);
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





// AMD_2'S GARBAGE COLLECTION, ported. The pool fills because absorbed cliques and shortened runs
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
// collector runs during a clique build the partial sits at the top of the pool with no descriptor
// of its own, so the sweep stops just below it, `AMD_2`'s `pend = pme1 - 1`, and the partial is
// then moved down explicitly and its start rewritten. Called with `cliqueStart == mFree` outside a
// build, where the partial is empty and both halves degenerate to the plain sweep.
void QuotientGraphCompacted::garbageCollect(std::size_t& cliqueStart) {
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

void QuotientGraphCompacted::beginElimination(std::int32_t pivot) {
    // `Nv [me] = -nvpiv` in Amd.cpp, and it is what keeps the pivot out of its own clique: the
    // walks below take a vertex only when its weight reads positive. The old code needed an
    // explicit `mMark[u] = mTag` for the same purpose.
    mWeight[pivot] = -mWeight[pivot];

    // TWO WAYS TO BUILD THE CLIQUE, and which one applies is `AMD_2`'s `elenme == 0`. With no
    // elements the reach is a subset of A[pivot] and is compacted where it stands; otherwise it is
    // assembled at the free cursor. The pool is touched only in the second case, which is why the
    // branch is worth having and not merely faithful.
    const bool inPlace = mRun[pivot].incidenceSize == 0;
    mBuiltInPlace = inPlace;

    std::size_t   cliqueStart;
    std::uint32_t cliqueLen;
    if (inPlace) {
        cliqueStart = mRun[pivot].sourcePtr;
        cliqueLen   = reachableSetInPlace(pivot);
    } else {
        // NO RESERVATION BEFORE THE WALK. Until 2026-08-18 this tested `mSource.size() - mFree <
        // size()`, room for a worst-case reach of n, because the walk could not survive a
        // collection and so had to be sure of never needing one. `AMD_2` tests per entry inside
        // the walk instead, and now so does `reachableSet`. That was the last divergence in the
        // storage, and it was one-directional: we collected where `AMD_2` would not, never the
        // reverse.
        cliqueLen   = reachableSet(pivot);
        // AFTER the walk, not before it. The cursor is where the clique ends, so its start is a
        // subtraction, and a collection inside the walk moves both together and leaves this right.
        cliqueStart = mFree - cliqueLen;
    }

    const std::int32_t* reached     = mSource.data() + cliqueStart;
    const std::uint32_t reachedSize = cliqueLen;

    // The absorbed cliques die only after the walk has read them. In the in-place case there are
    // none, the incidence list being empty, so this loop is skipped along with the branch above.
    const std::int32_t* absorbedCliques = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t absorbedSize    = mRun[pivot].incidenceSize;
    for (std::uint32_t i = 0; i < absorbedSize; ++i)
        killClique(absorbedCliques[i]);               // dead, its block left behind

    mRun[pivot].sourcePtr     = cliqueStart;
    mRun[pivot].adjacencySize = cliqueLen;
    mRun[pivot].incidenceSize = 0;                    // a clique has no incidence list

    // A CLIQUE IS BORN HERE AND NOWHERE ELSE, which is what makes the peak countable, and the
    // maximum is taken here alone since nothing else raises the total.
    mNumLiveCliqueMembers += cliqueLen;
    mNumPeakCliqueMembers  = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);

    // NO STAMPING PASS, AND NO TAG. Membership was written by the walk, in the sign of the weight,
    // so this only sums. The `inClique` out-parameter and the two `++mTag` that fed it went on
    // 2026-08-17: the prune reads the sign, nothing read the tag, and production's signature is not
    // a reason for a private copy to carry a parameter it never uses.
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t k = 0; k < reachedSize; ++k)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[k]]);
    mCliqueWeight = cliqueWeight;
}


const std::vector<std::int32_t>& QuotientGraphCompacted::eliminate(std::int32_t pivot,
                                                                   TaggedScanCompacted& scan) {
    {
        const std::int32_t* absorbed = mSource.data() + mRun[pivot].sourcePtr;
        const std::uint32_t count    = mRun[pivot].incidenceSize;
        for (std::uint32_t i = 0; i < count; ++i) scan.w[absorbed[i]] = 0;
    }

    beginElimination(pivot);
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

    return finishElimination(pivot);
}

const std::vector<std::int32_t>& QuotientGraphCompacted::finishElimination(std::int32_t pivot) {
    if (mLateMassElimination) {
        mMerged.clear();
    } else {
        massEliminate(pivot);
    }

    mRun[pivot].incidenceSize = 0;
    return mMerged;
}

const std::vector<std::int32_t>& QuotientGraphCompacted::massEliminate(std::int32_t pivot) {
    std::vector<std::int32_t>& merged = mMerged;   // scratch, kept for its capacity
    merged.clear();
    const std::int32_t* reached     = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;
    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        if (mRun[u].adjacencySize == 0 && mRun[u].incidenceSize == 1 &&
            mSource[mRun[u].sourcePtr] == pivot) {         // I[u] starts at the run, always now
            mRun[u].incidenceSize = 0;
            merged.push_back(u);
        }
    }
    if (!merged.empty()) {
        // THE MERGE HAPPENS FIRST, so that the compaction can read the ZERO WEIGHT it leaves
        // rather than a stamp of its own. The old order was the reverse and needed a tag pass over
        // `merged` plus a mark array read per member; the weight says the same thing and the
        // supervariable bookkeeping had to write it anyway. This was the last reader of mMark's
        // VERTEX half, and with it gone the array is clique ids alone.
        for (std::int32_t u : merged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            mCliqueWeight -= static_cast<std::uint32_t>(-mWeight[u]);
            mWeight[pivot] += mWeight[u];              // both negative; magnitudes add
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
        // `if (elenme != 0) pfree = p`: "element was not constructed in place: deallocate part of
        // it since newly nonprincipal variables may have been removed". Only in the pooled case,
        // and the in-place case has nothing to give back since it never took any.
        //
        // SAFE BECAUSE THE CLIQUE IS STILL THE LAST BLOCK IN THE POOL. Nothing is written to the
        // pool between beginElimination building it and this function trimming it: the prune
        // rewrites vertex runs in place, and absorption only zeroes lengths. If that ever stops
        // being true this becomes wrong silently, so the condition is asserted rather than assumed.
        if (!mBuiltInPlace) {
            assert(mRun[pivot].sourcePtr + membersSize == mFree &&
                   "the pivot's clique is no longer the last block in the pool");
            mFree = mRun[pivot].sourcePtr + kept;
        }
    }
    return merged;
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
// is never given back, so the pool fills faster and the collector runs sooner than `AMD_2`'s.
void QuotientGraphCompacted::trimClique(std::int32_t pivot, std::uint32_t kept) {
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

void QuotientGraphCompacted::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    // NOT `killClique`. v is a live supervariable being absorbed and never formed a clique, so
    // this length is A[v]'s, not a clique's, and feeding it to the counter would corrupt the peak.
    mRun[v].adjacencySize = 0;
    mRun[v].incidenceSize = 0;
}

void QuotientGraphCompacted::killClique(std::int32_t c) {
    mNumLiveCliqueMembers -= mRun[c].adjacencySize;
    mRun[c].adjacencySize = 0;
}

void QuotientGraphCompacted::absorb(const std::vector<std::int32_t>& cliques,
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
    // rewrites the list, so an absorbed element is simply not copied. We decide absorption after
    // the whole prune, because the prune has scan 1 fused into it, and so have to revisit the list.
    // What is removed here is the slide, not the visit.
    for (std::uint32_t k = 0; k < vertexCount; ++k) {
        const std::int32_t u         = vertices[k];
        std::int32_t*      incidence = mSource.data() + mRun[u].sourcePtr;
        const std::uint32_t size     = mRun[u].incidenceSize;

        // THE ROTATION HAS TO BE REDONE WHEN ITS SUBJECT DIES. `AMD_2` rotates inside scan 2, with
        // this step's absorbed elements already dropped, so the entry it parks at the back is the
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


std::vector<std::int32_t> QuotientGraphCompacted::order(
        const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order;
    order.reserve(size());
    for (std::int32_t pivot : pivots)
        for (std::int32_t u = pivot; u != NIL; u = mSuperNext[u]) order.push_back(u);
    return order;
}

} // namespace

std::vector<std::int32_t> orderAmd3B(const std::vector<std::size_t>&  colPtr,
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

    // THE EMPTY-ROW PREPASS, riding in the filing loop as `AMD_2` does and as Amd3 now does. See
    // src/Amd3.cpp for why it exists and what it fixes; the two must produce the same permutation
    // and `test_order` asserts it.
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
    BucketsCompacted buckets(size);
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


    std::vector<std::int32_t> hashHead(size + 1, NIL);

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
        TaggedScanCompacted scan{&buckets, w, degrees, touchedCliques, wflg,
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
            const std::int32_t* incidence     = qg.incidence(u);
            const std::uint32_t incidenceSize = qg.incidenceSize(u);

            std::size_t deg = static_cast<std::size_t>(w[u]);   // the adjacency half
            std::uint32_t key = static_cast<std::uint32_t>(buckets.key(u));
            for (std::uint32_t i = 0; i < incidenceSize; ++i) {
                const std::int32_t c = incidence[i];
                if (c != pivot) key += static_cast<std::uint32_t>(c);   // me is not in the key
                if (c != pivot) deg += static_cast<std::size_t>(w[c] - wflg);
            }

            w[u] = static_cast<std::int32_t>(std::min<std::size_t>(deg, degrees[u]));

            if (!qg.eliminated(u)) {
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
            if (qg.eliminated(seed)) continue;
            const std::int32_t hash = buckets.key(seed);
            const std::int32_t headOfBucket = hashHead[hash];
            if (headOfBucket == NIL) continue;      // an earlier member already emptied it
            hashHead[hash] = NIL;

            for (std::int32_t u = headOfBucket; u != NIL && buckets.chain(u) != NIL;
                 u = buckets.chain(u)) {
                if (qg.eliminated(u)) continue;

                // THE STAMP GOES INTO `w`, WHICH RETIRES mMark. Amd.cpp does exactly this,
                // `W [Iw [p]] = wflg` over the whole of i's list, variables and elements alike,
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
                const std::int32_t* runU       = qg.incidence(u);        // the run's first entry
                const std::uint32_t incidenceU = qg.incidenceSize(u);
                const std::uint32_t adjacencyU = qg.adjacencySize(u);
                const std::uint32_t runSizeU   = incidenceU + adjacencyU;
                for (std::uint32_t a = 1; a < runSizeU; ++a) w[runU[a]] = other;

                for (std::int32_t v = buckets.chain(u); v != NIL; v = buckets.chain(v)) {
                    if (qg.eliminated(v)) continue;

                    // THE LENGTHS REJECT BEFORE THE LIST IS TOUCHED, which is Amd.cpp's
                    // `ok = (Len [j] == ln) && (Elen [j] == eln)`. Two compares throw out most
                    // candidates for nothing, where counting live entries as we walk cannot decide
                    // until the walk is over.
                    if (qg.incidenceSize(v) != incidenceU) continue;
                    if (qg.adjacencySize(v) != adjacencyU) continue;

                    bool                same = true;
                    const std::int32_t* runV = qg.incidence(v);
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
            if (qg.eliminated(u)) continue;         // absorbed by the hash a moment ago
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

    // How often the pool actually needed collecting. `AMD_2` reports the same figure as
    // Info[AMD_NCMPA] and its complexity bound assumes it stays constant; this is where that
    // assumption can be checked for OUR storage rather than for its.
    gAmd3BCompactions  = qg.compactions();
    gPeakCliqueMembers = qg.numPeakCliqueMembers();   // see src/QuotientGraph.cpp
    // THE ROWS THE DENSE RULE SET ASIDE GO LAST, in index order, which is where `AMD_2`'s output
    // assembly puts them. They were collected in an ascending pass, so appending the vector is
    // that order; each stands only for itself, having been set aside before it could absorb
    // anything, so `order` expands a chain of one.
    pivots.insert(pivots.end(), denseRows.begin(), denseRows.end());
    return qg.order(pivots);
}

} // namespace Oblio
