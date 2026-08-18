#include "oblio/Mmd3C.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

// Mmd3C.cpp - Mmd3 on AMD_2'S CLIQUE LAYOUT: one pooled workspace with a free cursor and a garbage
// collection, where production keeps a separate append-only arena in elimination order.
//
// IT IS A CELL OF THE LAYOUT MATRIX, and the matrix is why it exists. Two algorithms and three
// clique layouts give six cells; with only the diagonal filled, every reading confounds the layout
// with the algorithm that happens to use it. B is a driver on its OWN branch's vendored layout and
// C is a driver on the OTHER branch's, so this is the mmd counterpart of Amd3B and pairs with it
// down a column. See docs/DESIGN_DECISIONS.md (2026-08-16, the layout matrix, and 2026-08-17).
//
// WHAT IT IS FOR, precisely. With encoding held equal, genmmd's dead-segment scheme COSTS Mmd3
// about 8 percent on squares and 29 on cubes, and AMD_2's pool EARNS Amd3 about 9 percent, rising
// to 15 at large n. Each of those is a single reading on a single algorithm, so neither can be told
// apart from something about how that algorithm walks. This file is the second reading of the pool:
// if mmd on it also comes in below 1.0, the pool wins on both algorithms and our arena is the thing
// to change.
//
// A PREVIOUS Mmd3C HELD THIS NAME AND WAS SOMETHING ELSE. It was mmd on the PRODUCTION layout, a
// transitional vehicle for working the amd array folds out on the mmd side without disturbing the
// class six drivers run. It did that; the folds are in QuotientGraph and the file was replaced.
// Nothing of it survives here but the name.
//
// ITS OBLIGATION. It must return Mmd3's permutation exactly, which is genmmd's. `make digest` in
// benchmarks/ordering hashes every driver's permutation over 73 grids and names which one moved;
// `make mmdorder` in experiments/ordering is what says correct. And it must stay ENCODING-IDENTICAL
// to Mmd3: a fold that lands in QuotientGraph lands here too, or its time column stops being about
// storage and starts being about whatever drifted. That is the same obligation Mmd3B carries, and
// Mmd3B was found to have broken it on 2026-08-17.
//
// THE COPY IS MECHANICAL WHERE IT CAN BE. QuotientGraphC and BucketsC are include/oblio/
// QuotientGraph.h and src/QuotientGraph.cpp with whole-line comments removed and the members no mmd
// driver calls pruned; the driver is src/Mmd3.cpp. tmp/make_mmd3c.py generated it, and regenerating
// from CURRENT sources beats restoring a stale copy, which reintroduces a divergence instead of
// starting from zero. THE STORAGE LAYER IS NOT MECHANICAL: it is a hand port from src/Amd3B.cpp,
// which is the only other file that has it. Design notes for the shared types live in the files
// above and are authoritative there; only the DIFFERENCES are commented here.
//
// AT THIS COMMIT THE STORAGE HAS NOT LANDED YET and this is a verbatim copy of Mmd3, which is
// deliberate: a vehicle that starts identical proves the extraction before it carries anything.

namespace Oblio {

// HOW OFTEN THE POOL RAN OUT, read by tmp/ probes as gAmd3BCompactions is on the amd side. It is
// the one number that says whether `AMD_2`'s elbow room, the pattern plus a fifth plus n, is the
// right size for MMD's cliques: the two branches fill the pool at different rates, amd absorbing
// cliques aggressively where mmd does not, so the sizing is the part of this port that is not a
// copy. Non-zero and growing would mean the room is too small and the compactions are being paid
// for repeatedly.
std::size_t gMmd3CCompactions = 0;

namespace {

class BucketsC {
public:
    static constexpr std::int32_t UNFILED    = 0;
    static constexpr std::int32_t OUTMATCHED = -2147483647 - 1;   // INT32_MIN; no degree reaches it

    explicit BucketsC(std::size_t size)
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

struct TaggedScanC {
    BucketsC*                          buckets;
    std::vector<std::int32_t>&        w;             // per clique, Amd.cpp's tagged W
    const std::vector<std::uint32_t>& degree;
    std::vector<std::int32_t>&        touchedCliques;// the cliques this step reached, once each
    std::int32_t                      wflg;          // the tag for this elimination
    std::int32_t                      modulus;       // the driver's bucket count, n + 1
};

class QuotientGraphC {
public:
    QuotientGraphC(const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const           { return mRun.size(); }

    // THE POOL'S SIZE, in place of production's arena entries. Not the same quantity and not
    // meant to be: production's grows monotonically and this one is fixed at construction, so what
    // this reports is the elbow room the layout needed rather than the space the cliques consumed.
    // `compactions()` is the number that says whether the room was enough.
    std::size_t arenaEntries() const   { return mSource.size(); }
    // A BLOCK HEAD IN THE POOL, during garbage collection only: FLIPPED - e for owner e, so every
    // head is at most FLIPPED and every real entry, a vertex id, is above it. AMD_2's FLIP.
    static constexpr std::int32_t FLIPPED = -2;

    static constexpr std::int32_t GONE = 2147483647;   // INT32_MAX, above every reachable tag

    bool eliminated(std::int32_t u) const { return mMark[u] == GONE; }

    std::int32_t advanceTag()                    { return ++mTag; }
    std::int32_t mark(std::int32_t v) const      { return mMark[v]; }
    void setMark(std::int32_t v, std::int32_t t) { mMark[v] = t; }

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

    std::uint32_t weight(std::int32_t u) const {
        const std::int32_t w = mWeight[u];
        return static_cast<std::uint32_t>(w < 0 ? -w : w);
    }

    void reachableSet(std::int32_t u);                   // writes at mFree
    std::uint32_t reachableSetInPlace(std::int32_t u);   // rewrites A[u] as C[u]; see its note


    std::uint32_t reachableWeight(std::int32_t u);

    const std::vector<std::int32_t>& eliminate(std::int32_t pivot);


    void merge(std::int32_t u, std::int32_t v);

    void number(std::int32_t u);

    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    void setVendoredListOrder(bool on) { mVendoredListOrder = on; }

    void setLateMassElimination(bool on) { mLateMassElimination = on; }

    const std::vector<std::int32_t>& massEliminate(std::int32_t pivot);



    std::vector<std::int32_t> orderAscending(const std::vector<std::int32_t>& pivots) const;

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
    // append-only arena that only grows; `AMD_2` keeps variable lists and clique lists in one
    // workspace with elbow room, builds a new clique at a free cursor, leaves absorbed space dead,
    // and COMPACTS when the cursor reaches the end. `mFree` is that cursor, `garbageCollect` that
    // compaction, and `mSource` now holds both kinds of block.
    //
    // A block's descriptor is `mRun[e]` for either kind, which is what makes one pool possible at
    // all: a live vertex's run is A[e] then I[e], and a dead pivot's is C[e].
    std::size_t mFree        = 0;             // AMD_2's pfree
    // Whether the clique now being eliminated was built in the pivot's own run. Read by
    // massEliminate, which can only give space back to the cursor in the other case.
    bool        mBuiltInPlace = false;
    std::size_t mCompactions = 0;             // AMD_2's Info[AMD_NCMPA]
public:
    std::size_t compactions() const { return mCompactions; }
private:
    void garbageCollect();

    std::vector<std::int32_t>  mSuperNext;
    std::vector<std::int32_t>  mSuperLast;
    std::vector<std::int32_t>  mWeight;
    std::uint32_t              mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

    bool mReverseIncidence = false;

    bool mVendoredListOrder = false;

    bool mLateMassElimination = false;

    bool mHasNumbered = false;

    std::vector<std::int32_t> mMark;     // membership scratch, read against mTag
    std::int32_t              mTag = 0;
};

QuotientGraphC::QuotientGraphC(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mRun(colPtr.empty() ? 0 : colPtr.size() - 1),
      mMark(mRun.size(), NIL) {
    const std::int32_t size = static_cast<std::int32_t>(mRun.size());

    // ELBOW ROOM, `AMD_2`'s. It sizes `slen = nzaat + nzaat/5 + 7n` and hands `iwlen = slen - 6n`
    // to the main loop, so the pool is the pattern plus a fifth plus n. Ours holds the same lists,
    // so the same shape applies. Reserved BEFORE the runs are laid down, so nothing reallocates
    // under a `sourcePtr` already taken.
    const std::size_t nnz = colPtr.empty() ? 0 : colPtr.back();
    mSource.reserve(nnz + nnz / 5 + mRun.size() + 1);
    for (std::int32_t aj = 0; aj < size; ++aj) {
        mRun[aj].sourcePtr = mSource.size();
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSource.push_back(rowIdx[cp]);
        mRun[aj].adjacencySize = static_cast<std::uint32_t>(mSource.size() - mRun[aj].sourcePtr);
    }

    // The runs are laid out; everything past them is free space the cliques will use. The vector is
    // SIZED rather than left at its length: a cursor into a vector's unused capacity would be
    // undefined, and this space is written by cliques and read by nothing until it is.
    mFree = mSource.size();
    mSource.resize(mSource.capacity());

    mSuperNext.assign(size, NIL);
    mSuperLast.resize(size);
    mWeight.assign(size, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

// COMPACTION, `AMD_2`'s `garbage_collection`. Absorbed blocks are left dead where they lie, so the
// cursor eventually reaches the end of the pool and everything still live has to be slid down.
//
// THE TRICK, and it is what makes this possible without a second array: each live block's FIRST
// ENTRY is temporarily replaced by `FLIPPED - e`, naming its owner, and the displaced entry is
// parked in `mRun[e].sourcePtr`, which is about to be overwritten anyway. The sweep then reads the
// pool linearly and can tell a block head from dead space by the value alone, every head being at
// most FLIPPED and every real entry a vertex id above it. Amd.cpp does exactly this with FLIP.
//
// A DEAD BLOCK IS ONE WITH NO LENGTH. Both kinds are covered by one test: a vertex merged away has
// both list lengths zeroed, and an absorbed clique has its adjacencySize zeroed by
// beginElimination.
void QuotientGraphC::garbageCollect() {
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

// CONSTRUCTS THE CLIQUE IN THE PIVOT'S OWN RUN, which is `AMD_2`'s `if (elenme == 0)` branch and
// the common case rather than a corner: a pivot with no elements has a reach that is a SUBSET of
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
std::uint32_t QuotientGraphC::reachableSetInPlace(std::int32_t u) {
    // THE PIVOT IS NEGATED HERE, as the appending walk does at its own head. It is not needed to
    // keep u out of its own list, which has no self loop, but massEliminate restores the SIGN of
    // every member and of the pivot, so a pivot left positive comes out of the restore negative
    // and `orderAscending` then reads a supervariable of negative size. Found by ASan on a 3 by 3
    // grid; the write went four billion entries past the permutation.
    mWeight[u] = -mWeight[u];
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

// WRITES AT `mFree` RATHER THAN APPENDING, which is the storage change. The caller has already
// made room for a whole reach, so the cursor cannot run out inside the walk; and the reach lands
// exactly where the clique is to live, so there is no copy from a scratch into place. Production
// appends to its arena and needs a capacity check per pivot for the same guarantee.
void QuotientGraphC::reachableSet(std::int32_t u) {
    ++mTag;
    mWeight[u] = -mWeight[u];              // never its own neighbor
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
            for (std::uint32_t k = 0; k < membersSize; ++k) {
                const std::int32_t v  = members[k];
                const std::int32_t nv = mWeight[v];
                if (nv > 0) { mWeight[v] = -nv; mSource[mFree++] = v; }
            }
        }
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v  = source[k];
            const std::int32_t nv = mWeight[v];
            if (nv > 0 && !(mHasNumbered && mMark[v] == GONE)) { mWeight[v] = -nv; mSource[mFree++] = v; }
        }
        return;
    }

    for (std::uint32_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v  = source[k];
        const std::int32_t nv = mWeight[v];
        if (nv > 0 && !(mHasNumbered && mMark[v] == GONE)) { mWeight[v] = -nv; mSource[mFree++] = v; }
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
}



std::uint32_t QuotientGraphC::reachableWeight(std::int32_t u) {
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

void QuotientGraphC::number(std::int32_t u) {
    mHasNumbered = true;
    mMark[u]     = GONE;
}

void QuotientGraphC::beginElimination(std::int32_t pivot) {
    // TWO WAYS TO BUILD THE CLIQUE, and which one applies is `AMD_2`'s `elenme == 0`. With no
    // elements the reach is a subset of A[pivot] and is compacted where it stands; otherwise it is
    // assembled at the free cursor. The pool is touched only in the second case.
    const bool inPlace = mRun[pivot].incidenceSize == 0;
    mBuiltInPlace = inPlace;

    std::size_t   cliqueStart;
    std::uint32_t cliqueLen;
    if (inPlace) {
        cliqueStart = mRun[pivot].sourcePtr;
        cliqueLen   = reachableSetInPlace(pivot);
    } else {
        // ROOM FOR A WHOLE REACH BEFORE THE WALK STARTS, which is what makes the cursor safe
        // inside it. A reach has at most n entries, so room for n is room for any of them.
        // `AMD_2` tests the same thing per entry, `if (pfree >= iwlen) garbage_collection`, and
        // can because it writes one entry at a time from a loop it can resume; ours writes from a
        // walk holding pointers into the pool, so the test has to come first and be for the worst
        // case. THIS IS THE ONE KNOWN DIVERGENCE LEFT in the storage: we can collect where AMD_2
        // would not, never the reverse.
        if (mSource.size() - mFree < size()) garbageCollect();
        cliqueStart = mFree;
        reachableSet(pivot);                    // writes at mFree; see its note
        cliqueLen   = static_cast<std::uint32_t>(mFree - cliqueStart);
    }

    const std::int32_t* reached     = mSource.data() + cliqueStart;
    const std::uint32_t reachedSize = cliqueLen;

    const std::int32_t* absorbedCliques = mSource.data() + mRun[pivot].sourcePtr + mRun[pivot].adjacencySize;
    const std::uint32_t absorbedSize    = mRun[pivot].incidenceSize;
    for (std::uint32_t i = 0; i < absorbedSize; ++i)
        mRun[absorbedCliques[i]].adjacencySize = 0;   // dead, its block left behind

    mRun[pivot].sourcePtr     = cliqueStart;
    mRun[pivot].adjacencySize = cliqueLen;

    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t k = 0; k < reachedSize; ++k)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[k]]);
    mCliqueWeight = cliqueWeight;
}

const std::vector<std::int32_t>& QuotientGraphC::eliminate(std::int32_t pivot) {
    beginElimination(pivot);
    const std::int32_t* reached     = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;
    const bool amdOrder = mVendoredListOrder;      // hoisted, as the other flags are

    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        const VertexRun&  run           = mRun[u];
        std::int32_t*     source        = mSource.data() + run.sourcePtr;
        const std::uint32_t adjacencySize = run.adjacencySize;
        std::uint32_t       kept          = 0;
        std::int32_t        heldVertex    = NIL;        // the first survivor, appended last
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (mWeight[v] <= 0) continue;             // in the new clique, the pivot, or merged
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by a prepass; see above
            if (amdOrder && heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mRun[u].adjacencySize = kept;                      // A[u] - C[pivot] - {pivot}

        const std::int32_t* incidence     = source + adjacencySize;
        const std::uint32_t incidenceSize = run.incidenceSize;
        std::uint32_t       write         = kept;      // follows `kept`; see the note above
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
            if (mRun[c].adjacencySize != 0) source[write++] = c;   // dead is size zero; see above
        }
        source[write++]   = pivot;                     // u joins the new clique, id = pivot
        mRun[u].incidenceSize = write - kept;
        if (amdOrder && write - kept > 1) std::swap(source[kept], source[write - 1]);
    }

    return finishElimination(pivot);
}


const std::vector<std::int32_t>& QuotientGraphC::finishElimination(std::int32_t pivot) {
    if (mLateMassElimination) {
        mMerged.clear();
    } else {
        massEliminate(pivot);
    }

    mRun[pivot].incidenceSize = 0;
    mMark[pivot]          = GONE;
    return mMerged;
}

const std::vector<std::int32_t>& QuotientGraphC::massEliminate(std::int32_t pivot) {
    std::vector<std::int32_t>& merged = mMerged;   // scratch, kept for its capacity
    merged.clear();
    const std::int32_t* reached     = mSource.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;
    mWeight[pivot] = -mWeight[pivot];
    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        mWeight[u] = -mWeight[u];                          // live again, and positive
        if (mRun[u].adjacencySize == 0 && mRun[u].incidenceSize == 1 &&
            mSource[mRun[u].sourcePtr] == pivot) {         // A[u] empty, so I[u] starts at the run
            mRun[u].incidenceSize = 0;
            mMark[u]          = GONE;
            merged.push_back(u);
        }
    }
    if (!merged.empty()) {
        for (std::int32_t u : merged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            mCliqueWeight -= static_cast<std::uint32_t>(mWeight[u]);
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }

        std::int32_t*     members     = mSource.data() + mRun[pivot].sourcePtr;
        const std::uint32_t membersSize = mRun[pivot].adjacencySize;
        std::uint32_t       kept        = 0;
        for (std::uint32_t k = 0; k < membersSize; ++k)
            if (mWeight[members[k]] != 0) members[kept++] = members[k];
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

void QuotientGraphC::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    mRun[v].adjacencySize = 0;
    mRun[v].incidenceSize = 0;
    mMark[v]          = GONE;
}


std::vector<std::int32_t> QuotientGraphC::orderAscending(
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

} // namespace

namespace {

std::vector<std::int32_t> orderMmd3CImpl(const std::vector<std::size_t>&  colPtr,
                                        const std::vector<std::int32_t>& rowIdx,
                                        std::int32_t delta,
                                        std::size_t* arenaEntries) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraphC qg(colPtr, rowIdx);
    qg.setReverseIncidence(true);
    std::vector<std::int32_t> pivots;
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;

    BucketsC buckets(size);
    std::uint32_t minDegree = static_cast<std::uint32_t>(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u) {
        const std::uint32_t degree = std::max<std::uint32_t>(qg.adjacencySize(u), 1);
        buckets.file(degree, u);
        minDegree = std::min(minDegree, degree);
    }

    std::vector<std::int32_t> batch, elementMembers, q2h, qxh;

    for (std::int32_t u = buckets.head(1); u != NIL; ) {
        const std::int32_t next = buckets.next(u);   // before the unfile invalidates it
        buckets.unfile(u);
        qg.number(u);
        pivots.push_back(u);
        ++numEliminated;
        u = next;
    }
    if (size > 2) minDegree = 2;                // head[1] is empty now, and mdeg starts at 2

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;

        std::uint32_t batchLimit = minDegree;
        if (delta > 0)
            batchLimit = std::min(minDegree + static_cast<std::uint32_t>(delta),
                                  static_cast<std::uint32_t>(size) - 1);

        batch.clear();
        while (true) {
            if (buckets.empty(minDegree)) {
                if (minDegree >= batchLimit) break;
                ++minDegree;
                continue;
            }
            const std::int32_t pivot = buckets.head(minDegree);
            buckets.unfile(pivot);

            const std::vector<std::int32_t>& merged = qg.eliminate(pivot);
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
            for (std::int32_t u : merged) buckets.unfile(u);

            const std::int32_t* clique     = qg.clique(pivot);
            const std::uint32_t cliqueSize = qg.cliqueSize(pivot);
            for (std::uint32_t k = 0; k < cliqueSize; ++k) {
                const std::int32_t u = clique[k];
                buckets.unfile(u);                  // evict; mmdelm's bwd[rn] = 0 is both this
                buckets.restore(u);                 //   and putting a withheld vertex back
            }

            if (numEliminated >= size) break;       // genmmd: nothing left to update
            if (delta < 0) break;
        }

        for (auto ee = batch.rbegin(); ee != batch.rend(); ++ee) {
            const std::int32_t element = *ee;
            const std::int32_t* members     = qg.clique(element);
            const std::uint32_t membersSize = qg.cliqueSize(element);

            elementMembers.clear();
            for (std::uint32_t k = 0; k < membersSize; ++k)
                if (!qg.eliminated(members[k])) elementMembers.push_back(members[k]);

            const std::int32_t elementTag = qg.advanceTag();   // marked once for the element
            for (std::int32_t v : elementMembers) qg.setMark(v, elementTag);
            std::uint32_t dg0 = 0;
            for (std::int32_t v : elementMembers) dg0 += qg.weight(v);

            q2h.clear();
            qxh.clear();
            for (std::int32_t u : elementMembers) {
                if (buckets.filed(u) || buckets.outmatched(u)) continue;   // done, or withheld
                const std::uint32_t otherSources = qg.adjacencySize(u) + qg.incidenceSize(u) - 1;
                (otherSources == 1 ? q2h : qxh).push_back(u);
            }

            for (auto uu = q2h.rbegin(); uu != q2h.rend(); ++uu) {
                const std::int32_t u = *uu;
                if (qg.eliminated(u) || buckets.outmatched(u)) continue;   // by an earlier q2h vertex
                const std::int32_t vertexTag = qg.advanceTag();
                std::uint32_t degree = dg0;

                const std::int32_t* adjacency = qg.adjacency(u);
                for (std::uint32_t a = 0; a < qg.adjacencySize(u); ++a) {
                    const std::int32_t v = adjacency[a];
                    const std::int32_t m = qg.mark(v);
                    if (m >= vertexTag) continue;                  // seen this pass, or dead
                    if (m == elementTag) continue;                 // already counted in dg0
                    qg.setMark(v, vertexTag);
                    degree += qg.weight(v);
                }
                const std::int32_t* incidence = qg.incidence(u);
                for (std::uint32_t i = 0; i < qg.incidenceSize(u); ++i) {
                    const std::int32_t c = incidence[i];
                    if (c == element) continue;
                    const std::int32_t* other     = qg.clique(c);
                    const std::uint32_t otherSize = qg.cliqueSize(c);
                    for (std::uint32_t k = 0; k < otherSize; ++k) {
                        const std::int32_t v = other[k];
                        const std::int32_t m = qg.mark(v);
                        if (v == u || m >= vertexTag) continue;    // seen this pass, or dead
                        if (m == elementTag) {
                            if (buckets.filed(v) || buckets.outmatched(v)) continue;
                            if (qg.adjacencySize(v) + qg.incidenceSize(v) - 1 == 1) {
                                qg.merge(u, v);      // identical reach: u absorbs it
                                ++numEliminated;
                            } else {
                                buckets.outmatch(v);    // reaches more, so never minimal first
                            }
                            continue;
                        }
                        qg.setMark(v, vertexTag);
                        degree += qg.weight(v);
                    }
                }

                const std::uint32_t filed = std::max<std::uint32_t>(degree - qg.weight(u) + 1, 1);
                buckets.file(filed, u);
                minDegree = std::min(minDegree, filed);
            }

            for (auto uu = qxh.rbegin(); uu != qxh.rend(); ++uu) {
                const std::int32_t u = *uu;                 // the full union, as md5 computes it
                if (qg.eliminated(u) || buckets.outmatched(u)) continue;
                const std::uint32_t degree = qg.reachableWeight(u); // reach excludes u already
                const std::uint32_t filed = std::max<std::uint32_t>(degree + 1, 1);
                buckets.file(filed, u);
                minDegree = std::min(minDegree, filed);
            }
        }

    }

    if (arenaEntries != nullptr) *arenaEntries = qg.arenaEntries();
    gMmd3CCompactions = qg.compactions();
    return qg.orderAscending(pivots);   // genmmd's mmdnum. See the ledger, entry 6.
}

} // namespace

std::vector<std::int32_t> orderMmd3C(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    std::int32_t delta) {
    return orderMmd3CImpl(colPtr, rowIdx, delta, nullptr);
}

std::vector<std::int32_t> orderMmd3C(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    std::int32_t delta,
                                    std::size_t& arenaEntries) {
    return orderMmd3CImpl(colPtr, rowIdx, delta, &arenaEntries);
}

} // namespace Oblio
