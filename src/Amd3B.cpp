#include "oblio/Amd3B.h"

#include "oblio/Types.h"

#include <algorithm>
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
// IT CARRIES A PRIVATE COPY OF QuotientGraph, named QuotientGraphA, so the storage can be changed
// without touching the class every driver shares. THE COST IS EVERY SHARED FOLD LANDING TWICE, and
// it is accepted on purpose; `make digest` catches a copy that has stopped reproducing its original
// in half a second.
//
// ITS OBLIGATION. It must return Amd3's permutation exactly, which is `AMD_2`'s raw order, so its
// nnz(L) column in benchmarks/ordering must equal AMDraw's on every row and its fill column carries
// nothing. That obligation is what makes it an instrument rather than a second ordering, and
// `make digest` checks it across every driver in half a second.
//
// THE COPY IS MECHANICAL. QuotientGraphA and BucketsA are include/oblio/QuotientGraph.h and
// src/QuotientGraph.cpp with the whole-line comments removed and the members no amd driver calls
// pruned; the driver is src/Amd3.cpp. Every design note for these types lives in those files and is
// authoritative there. Only the DIFFERENCES are commented here.

namespace Oblio {

std::size_t gAmd3BCompactions = 0;   // read by tmp/ probes; see the note at the return below

namespace {

class BucketsA {
public:
    static constexpr std::int32_t UNFILED    = 0;
    static constexpr std::int32_t OUTMATCHED = -2147483647 - 1;   // INT32_MIN; no degree reaches it

    explicit BucketsA(std::size_t size)
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

struct TaggedScanA {
    BucketsA*                          buckets;
    std::vector<std::int32_t>&        w;             // per clique, Amd.cpp's tagged W
    const std::vector<std::uint32_t>& degree;
    std::vector<std::int32_t>&        touchedCliques;// the cliques this step reached, once each
    std::int32_t                      wflg;          // the tag for this elimination
    std::int32_t                      modulus;       // the driver's bucket count, n + 1
};

class QuotientGraphA {
public:
    // NO cliqueMarks ARGUMENT. Production takes one to size its mark array at n or 2n; this file
    // has no mark array, so there was nothing for the flag to decide and it survived only as a
    // member the constructor had to cast to void. Removed 2026-08-17.
    QuotientGraphA(const std::vector<std::size_t>&  colPtr,
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


    const std::int32_t* clique(std::int32_t c) const {
        return mSource.data() + mRun[c].sourcePtr;
    }
    std::uint32_t cliqueSize(std::int32_t c) const { return mRun[c].adjacencySize; }

    std::uint32_t cliqueWeight() const { return mCliqueWeight; }

    const std::int32_t* adjacency(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr;
    }
    std::uint32_t adjacencySize(std::int32_t u) const { return mRun[u].adjacencySize; }

    const std::int32_t* incidence(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr + mRun[u].adjacencySize;
    }
    std::uint32_t incidenceSize(std::int32_t u) const { return mRun[u].incidenceSize; }

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

    std::uint32_t reachableSet(std::int32_t u);   // writes at mFree, returns the length




    const std::vector<std::int32_t>& eliminate(std::int32_t pivot, TaggedScanA& scan);

    void merge(std::int32_t u, std::int32_t v);


    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    void setVendoredListOrder(bool on) { mVendoredListOrder = on; }

    void setLateMassElimination(bool on) { mLateMassElimination = on; }

    const std::vector<std::int32_t>& massEliminate(std::int32_t pivot);

    void absorb(const std::vector<std::int32_t>& cliques,
                const std::int32_t* vertices, std::uint32_t vertexCount);

    std::vector<std::int32_t> order(const std::vector<std::int32_t>& pivots) const;


private:
    void beginElimination(std::int32_t pivot);

    const std::vector<std::int32_t>& finishElimination(std::int32_t pivot);

    std::vector<std::int32_t>  mSource;         // every A[u] then I[u], run after run
    struct VertexRun {
        std::size_t   sourcePtr;       // where u's run starts in mSource, fixed at construction
        std::uint32_t adjacencySize;   // A[u]'s length, from the run's start
        std::uint32_t incidenceSize;   // I[u]'s length, immediately behind A[u]
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
    std::size_t mCompactions = 0;             // AMD_2's Info[AMD_NCMPA]
public:
    std::size_t compactions() const { return mCompactions; }
private:

    void garbageCollect();

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

    bool mReverseIncidence = false;

    bool mVendoredListOrder = false;

    bool mLateMassElimination = false;

    // NO MARK ARRAY. Every question it used to answer is now read off the weight or the tagged
    // W: liveness and clique membership from the sign, absorbed from the zero, and supervariable
    // detection's stamp from `w`, which is what AMD_2 does with `W [Iw [p]] = wflg`. That is the
    // last array this file owned that AMD_2 does not, and it was 2n int32 wide.
};

QuotientGraphA::QuotientGraphA(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mRun(colPtr.empty() ? 0 : colPtr.size() - 1) {
    const std::int32_t size = static_cast<std::int32_t>(mRun.size());

    // ELBOW ROOM, `AMD_2`'s. It sizes `slen = nzaat + nzaat/5 + 7n` and hands `iwlen = slen - 6n`
    // to the main loop, so the pool is the pattern plus a fifth plus n. Ours holds the same lists,
    // so the same shape applies. The vector is SIZED rather than reserved: the space past the
    // runs is written by cliques and read by nothing until it is, and a cursor into a vector's
    // unused capacity would be undefined.
    const std::size_t nnz = colPtr.empty() ? 0 : colPtr.back();
    mSource.reserve(nnz + nnz / 5 + mRun.size() + 1);
    for (std::int32_t aj = 0; aj < size; ++aj) {
        mRun[aj].sourcePtr = mSource.size();
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSource.push_back(rowIdx[cp]);
        mRun[aj].adjacencySize = static_cast<std::uint32_t>(mSource.size() - mRun[aj].sourcePtr);
    }

    // The runs are laid out; everything past them is free space the cliques will use.
    mFree = mSource.size();
    mSource.resize(mSource.capacity());

    mSuperNext.assign(size, NIL);
    mSuperLast.resize(size);
    mWeight.assign(size, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

// WRITES AT `mFree` RATHER THAN APPENDING, which is the storage change. The caller has already
// ensured room for a whole reach, so the cursor cannot run past the pool inside this walk, and the
// pool never reallocates, so the pointers taken into it here stay valid for the whole walk. That
// second property is what the production version needs its reserve for and gets here for free.
std::uint32_t QuotientGraphA::reachableSet(std::int32_t u) {
    const std::size_t start = mFree;
    const VertexRun&    run           = mRun[u];          // one fetch; see the member
    const std::int32_t* source        = mSource.data() + run.sourcePtr;
    const std::uint32_t adjacencySize = run.adjacencySize;
    const std::uint32_t incidenceSize = run.incidenceSize;
    const std::int32_t* incidence     = source + adjacencySize;
    const bool reverse  = mReverseIncidence;
    const bool amdOrder = mVendoredListOrder;

    if (amdOrder) {
        for (std::uint32_t ii = 0; ii < incidenceSize; ++ii) {
            const std::uint32_t i           = reverse ? incidenceSize - 1 - ii : ii;
            const std::int32_t  c           = incidence[i];
            const std::int32_t* members     = mSource.data() + mRun[c].sourcePtr;
            const std::uint32_t membersSize = mRun[c].adjacencySize;
            // ONE LOAD PER MEMBER, and the store that takes a vertex into the clique is the same
            // slot. `> 0` is live-and-not-yet-taken, which is dedup and liveness together.
            for (std::uint32_t k = 0; k < membersSize; ++k) {
                const std::int32_t v  = members[k];
                const std::int32_t nv = mWeight[v];
                if (nv > 0) {
                    mWeight[v] = -nv;
                    mSource[mFree++] = v;
                }
            }
        }
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v  = source[k];
            const std::int32_t nv = mWeight[v];
            if (nv > 0) {
                mWeight[v] = -nv;
                mSource[mFree++] = v;
            }
        }
        return static_cast<std::uint32_t>(mFree - start);
    }

    // The adjacency-first order, on the same encoding. This driver never takes it, setting the
    // vendored list order, but it is folded rather than left reading a mark array that no longer
    // exists.
    for (std::uint32_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v  = source[k];
        const std::int32_t nv = mWeight[v];
        if (nv > 0) { mWeight[v] = -nv; mSource[mFree++] = v; }
    }
    for (std::uint32_t ii = 0; ii < incidenceSize; ++ii) {
        const std::uint32_t i           = reverse ? incidenceSize - 1 - ii : ii;
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mSource.data() + mRun[c].sourcePtr;
        const std::uint32_t membersSize = mRun[c].adjacencySize;
        for (std::uint32_t k = 0; k < membersSize; ++k) {
            const std::int32_t v  = members[k];
            const std::int32_t nv = mWeight[v];
            if (nv > 0) { mWeight[v] = -nv; mSource[mFree++] = v; }
        }
    }
    return static_cast<std::uint32_t>(mFree - start);
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
// a live vertex has A[e] then I[e], a dead pivot has C[e], and anything absorbed or merged had its
// lengths zeroed when it died.
void QuotientGraphA::garbageCollect() {
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
    while (src < mFree) {
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
    mFree = dst;
}

void QuotientGraphA::beginElimination(std::int32_t pivot) {
    // ROOM FOR A WHOLE REACH BEFORE THE WALK STARTS, which is what makes the cursor safe inside
    // it. A reach has at most n entries, so room for n is room for any of them. `AMD_2` tests the
    // same thing per entry, `if (pfree >= iwlen) garbage_collection`, and can because it writes one
    // entry at a time from a loop it can resume; ours writes from a walk with pointers taken into
    // the pool, so the test has to come first and be for the worst case.
    if (mSource.size() - mFree < size()) garbageCollect();

    // `Nv [me] = -nvpiv` in Amd.cpp, and it is what keeps the pivot out of its own clique: the
    // walk below takes a vertex only when its weight reads positive. The old code needed an
    // explicit `mMark[u] = mTag` for the same purpose.
    mWeight[pivot] = -mWeight[pivot];

    const std::size_t   cliqueStart = mFree;
    const std::uint32_t cliqueLen   = reachableSet(pivot);

    const std::int32_t* reached     = mSource.data() + cliqueStart;
    const std::uint32_t reachedSize = cliqueLen;

    const std::int32_t* absorbedCliques = mSource.data() + mRun[pivot].sourcePtr + mRun[pivot].adjacencySize;
    const std::uint32_t absorbedSize    = mRun[pivot].incidenceSize;
    for (std::uint32_t i = 0; i < absorbedSize; ++i)
        mRun[absorbedCliques[i]].adjacencySize = 0;   // dead, its block left behind

    mRun[pivot].sourcePtr     = cliqueStart;
    mRun[pivot].adjacencySize = cliqueLen;

    // NO STAMPING PASS, AND NO TAG. Membership was written by the walk, in the sign of the weight,
    // so this only sums. The `inClique` out-parameter and the two `++mTag` that fed it went on
    // 2026-08-17: the prune reads the sign, nothing read the tag, and production's signature is not
    // a reason for a private copy to carry a parameter it never uses.
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t k = 0; k < reachedSize; ++k)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[k]]);
    mCliqueWeight = cliqueWeight;
}


const std::vector<std::int32_t>& QuotientGraphA::eliminate(std::int32_t pivot, TaggedScanA& scan) {
    {
        const std::int32_t* absorbed = mSource.data() + mRun[pivot].sourcePtr
                                                      + mRun[pivot].adjacencySize;
        const std::uint32_t count    = mRun[pivot].incidenceSize;
        for (std::uint32_t i = 0; i < count; ++i) scan.w[absorbed[i]] = 0;
    }

    beginElimination(pivot);
    const std::int32_t* reached     = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;

    if (scan.buckets != nullptr)
        for (std::uint32_t ri = 0; ri < reachedSize; ++ri) scan.buckets->unfile(reached[ri]);

    const bool          amdOrder    = mVendoredListOrder;
    const std::int32_t  wflg        = scan.wflg;          // hoisted, as the flags are

    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u       = reached[ri];
        const std::int32_t nvi     = -mWeight[u];   // u is in C[pivot], so its weight is negative
        const std::int32_t wnvi    = wflg - nvi;          // Amd.cpp's wnvi, and signed for it
        const VertexRun&  run           = mRun[u];
        std::int32_t*     source        = mSource.data() + run.sourcePtr;
        const std::uint32_t adjacencySize = run.adjacencySize;
        std::uint32_t       kept          = 0;
        std::uint32_t       explicitPart  = 0;            // a weight sum, not a count of positions
        std::uint32_t       key           = 0;            // wraps, like Amd.cpp's UInt hval
        std::int32_t        heldVertex    = NIL;          // see the plain prune above
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            // ONE LOAD, THREE QUESTIONS, which is Amd.cpp's `nvj = Nv [j] ; if (nvj > 0)`. A
            // negative weight is a member of the new clique, including the pivot itself, so the
            // explicit `v == pivot` test is gone with the rest; a zero is absorbed or merged away.
            const std::int32_t v  = source[k];
            const std::int32_t nv = mWeight[v];
            if (nv <= 0) continue;
            explicitPart += static_cast<std::uint32_t>(nv);
            key += static_cast<std::uint32_t>(v);         // no + 1, no reduction; see above
            if (amdOrder && heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mRun[u].adjacencySize       = kept;
        scan.w[u]               = static_cast<std::int32_t>(explicitPart);
        if (scan.buckets != nullptr)
            scan.buckets->setKey(u, static_cast<std::int32_t>(key));   // the ADJACENCY half

        const std::int32_t* incidence     = source + adjacencySize;
        const std::uint32_t incidenceSize = run.incidenceSize;
        std::uint32_t       write         = kept;         // follows `kept`; see the plain prune
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
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
        source[write++]   = pivot;
        mRun[u].incidenceSize = write - kept;
        if (amdOrder && write - kept > 1) std::swap(source[kept], source[write - 1]);
    }

    return finishElimination(pivot);
}

const std::vector<std::int32_t>& QuotientGraphA::finishElimination(std::int32_t pivot) {
    if (mLateMassElimination) {
        mMerged.clear();
    } else {
        massEliminate(pivot);
    }

    mRun[pivot].incidenceSize = 0;
    return mMerged;
}

const std::vector<std::int32_t>& QuotientGraphA::massEliminate(std::int32_t pivot) {
    std::vector<std::int32_t>& merged = mMerged;   // scratch, kept for its capacity
    merged.clear();
    const std::int32_t* reached     = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;
    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        if (mRun[u].adjacencySize == 0 && mRun[u].incidenceSize == 1 &&
            mSource[mRun[u].sourcePtr] == pivot) {         // A[u] empty, so I[u] starts at the run
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
        mRun[pivot].adjacencySize = kept;
    }
    return merged;
}

void QuotientGraphA::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    mRun[v].adjacencySize = 0;
    mRun[v].incidenceSize = 0;
}

void QuotientGraphA::absorb(const std::vector<std::int32_t>& cliques,
                           const std::int32_t* vertices, std::uint32_t vertexCount) {
    if (cliques.empty()) return;

    for (std::int32_t c : cliques) mRun[c].adjacencySize = 0;   // dead; the compaction reads the size

    for (std::uint32_t k = 0; k < vertexCount; ++k) {  // I[u] - dead, compacted in place
        const std::int32_t u         = vertices[k];
        std::int32_t*      incidence = mSource.data() + mRun[u].sourcePtr + mRun[u].adjacencySize;
        const std::uint32_t size     = mRun[u].incidenceSize;
        std::uint32_t       kept     = 0;
        for (std::uint32_t i = 0; i < size; ++i)
            if (mRun[incidence[i]].adjacencySize != 0) incidence[kept++] = incidence[i];
        mRun[u].incidenceSize = kept;
    }
}


std::vector<std::int32_t> QuotientGraphA::order(const std::vector<std::int32_t>& pivots) const {
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

    QuotientGraphA qg(colPtr, rowIdx);

    qg.setVendoredListOrder(true);      // cliques before adjacency; the new clique at the front
    qg.setLateMassElimination(true);    // and mass elimination becomes this driver's, below

    std::vector<std::int32_t> pivots;               // the order over supervariables
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;                // a counter, not a scan of the flags

    std::uint32_t numLive = static_cast<std::uint32_t>(size);

    std::vector<std::uint32_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    BucketsA buckets(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        buckets.file(degrees[u], u);
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
        TaggedScanA scan{&buckets, w, degrees, touchedCliques, wflg,
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
                // IT HAS TO INTERLEAVE WITH THE TAG PROTOCOL, not clobber it. `stamp` starts above
                // the values this step's scan wrote, `wflg + lemax`, and rises by one per
                // candidate; `wflg` is then set past all of them at the end of the step. So next
                // step every stamped entry reads BELOW wflg, which the prune's `we >= wflg` test
                // treats as alive-and-unseen and rebuilds from the clique degree. That is the same
                // reading Amd.cpp relies on, and it is why stamping is safe here at all.
                //
                // The zeros survive: only entries of a live vertex's list are stamped, and a dead
                // clique is not in one.
                const std::int32_t other = ++stamp;
                std::uint32_t sizeU = 0;       // list entries, so at most deg(u)
                const std::int32_t* adjacencyU = qg.adjacency(u);
                for (std::uint32_t a = 0; a < qg.adjacencySize(u); ++a) {
                    const std::int32_t x = adjacencyU[a];
                    if (!qg.eliminated(x)) { w[x] = other; ++sizeU; }
                }
                const std::int32_t* incidenceU = qg.incidence(u);
                for (std::uint32_t i = 1; i < qg.incidenceSize(u); ++i) {
                    w[incidenceU[i]] = other;
                    ++sizeU;
                }

                for (std::int32_t v = buckets.chain(u); v != NIL; v = buckets.chain(v)) {
                    if (qg.eliminated(v)) continue;

                    std::uint32_t sizeV = 0;   // list entries, so at most deg(v)
                    bool        same  = true;
                    const std::int32_t* adjacencyV = qg.adjacency(v);
                    for (std::uint32_t a = 0; a < qg.adjacencySize(v) && same; ++a) {
                        const std::int32_t x = adjacencyV[a];
                        if (qg.eliminated(x)) continue;
                        ++sizeV;
                        if (w[x] != other) same = false;
                    }
                    if (same) {
                        const std::int32_t* incidenceV = qg.incidence(v);
                        for (std::uint32_t i = 1; i < qg.incidenceSize(v) && same; ++i) {
                            ++sizeV;
                            if (w[incidenceV[i]] != other) same = false;
                        }
                    }
                    if (!same || sizeU != sizeV) continue;

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
        }

        // PAST EVERY STAMP THIS STEP LAID DOWN, not merely past the scan's values. Amd.cpp
        // advances wflg through detection for the same reason: a stamp must read as alive-unseen
        // next step, which means strictly below the new tag.
        stamp = std::max(stamp, wflg + lemax);
        wflg  = stamp + 1;
    }

    // How often the pool actually needed collecting. `AMD_2` reports the same figure as
    // Info[AMD_NCMPA] and its complexity bound assumes it stays constant; this is where that
    // assumption can be checked for OUR storage rather than for its.
    gAmd3BCompactions = qg.compactions();
    return qg.order(pivots);
}

} // namespace Oblio
