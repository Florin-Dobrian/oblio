#include "oblio/Mmd3C.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

// Mmd3C.cpp - Mmd3 on the PRODUCTION clique layout, carrying a private quotient graph so that the
// amd branch's array folds can be worked out on the mmd side without disturbing the class six
// drivers share.
//
// IT IS TRANSITIONAL, WHICH IS THE ONE THING TO KNOW BEFORE READING IT. Mmd3B and Amd3B are
// permanent, each being the predictable-space version of its own ordering as well as an alignment
// vehicle. This file is neither. It exists to carry one transition, on the precedent of the AMD3C
// that held a re-schedule until it was proved and then went. When the transition is done it is
// expected to be REPLACED by an Mmd3C on AMD_2's clique layout, which is the cell the layout matrix
// in docs/DESIGN_DECISIONS.md (2026-08-16) reserves under this name. Same name, different content,
// and the second is a permanent file where this one is not.
//
// WHAT IT IS FOR, precisely. Amd3B carries five array folds found by a differential against AMD_2,
// three of which live in the shared quotient graph. Whether those three can be taken on the mmd
// branch is an open question with a known obstacle, and it is not a question to open by editing a
// class every driver runs. See experiments/ordering/README.md, "What each branch has to remember,
// and why an encoding does not transfer", for what is being asked and why the answer is not the
// amd one.
//
// AT THIS COMMIT IT IS A VERBATIM COPY and computes exactly what Mmd3 computes, which is genmmd's
// permutation. That is deliberate: a vehicle that starts identical proves the extraction before it
// carries anything, and it gives every other column in the benchmark an error bar measured under
// identical conditions, which docs/NEXT.md asks for and which costs one column.
//
// ITS OBLIGATION, for as long as it exists. It must return Mmd3's permutation exactly. `make
// digest` in benchmarks/ordering hashes every driver's permutation over 73 grids and names which
// one moved, which is what makes carrying a fourth copy of this class tolerable; `make mmdorder` in
// experiments/ordering is what says correct.
//
// THE COPY IS MECHANICAL. QuotientGraphC and BucketsC are include/oblio/QuotientGraph.h and
// src/QuotientGraph.cpp with the whole-line comments removed and the members no mmd driver calls
// pruned; the driver is src/Mmd3.cpp. Every design note for these types lives in those files and is
// authoritative there. Only the DIFFERENCES are commented here, and at this commit there are none.
// tmp/make_mmd3c.py is what generated it; regenerate from current sources rather than restoring a
// stale copy, which would reintroduce a divergence instead of starting from zero.

// THE PADDING PROBE, 2026-08-17, and it is temporary. `MMD3C` reads about 1.33x `MMD3` at exactly
// 200 a side and 0.96x and 0.88x at 199 and 201, on alpamayo, while computing the same permutation
// from the same algorithm. A spike at one n with flat neighbors cannot be algorithmic; it is
// addressing. This shifts every size-n allocation in this file by `kPad` elements and changes
// nothing else, which is the same intervention that settled `AMD_2`'s self-aliasing on 2026-08-16.
//
// In benchmarks/ordering, two runs. The baseline is `make clean && make scale-mmd-2d`. The padded
// one passes the whole flag set, since CXXFLAGS is assigned with `?=` and naming it on the command
// line replaces it rather than appending:
//
//   make clean
//   make CXXFLAGS="-std=c++17 -O3 -DNDEBUG -Wall -Wextra -I../../include -DOBLIO_BLAS_UNDERSCORE -DOBLIO_MMD3C_PAD=16" scale-mmd-2d
//
// IT IS SELF-CHECKING: only addresses move, so the permutation must come out identical and every
// nnzL column must be unchanged. If one moves, the probe is wrong and the reading is void.
// `make digest` says so in half a second.
//
// IT ANSWERED, 2026-08-17: DATA PLACEMENT. At `kPad = 1024` the 200^2 ratio falls from 1.28 to
// 0.99 while 199 and 201 hold at 0.95 and 0.97, permutations byte-identical. Same mechanism as
// AMD_2's self-aliasing; see docs/DESIGN_DECISIONS.md (2026-08-17).
//
// THE VALUE HAS TO BE A PAGE, NOT A CACHE LINE, and getting that wrong cost two rounds. The first
// attempt used sixteen ints, one line, and moved nothing, which read as evidence against data
// placement and was evidence of nothing: at 200 a side these arrays are 160,000 bytes, 40 pages,
// and 160,064 rounds to the same 40 pages, so the allocator very probably returned the same
// addresses and no intervention took place. A perturbation that does not perturb is inconclusive
// rather than negative.
#ifndef OBLIO_MMD3C_PAD
#define OBLIO_MMD3C_PAD 0
#endif

namespace Oblio {
namespace {

constexpr std::uint32_t kPad = OBLIO_MMD3C_PAD;   // extra elements per size-n vector; see above

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
    // NO cliqueMarks ARGUMENT. Production takes one to size mMark at n or 2n, the 2n form being
    // what amd's supervariable detection needs to stamp clique ids. No mmd driver stamps a clique,
    // so the true branch was unreachable in this file and the parameter decided nothing. Removed
    // 2026-08-17; mMark is n here, always.
    QuotientGraphC(const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const           { return mRun.size(); }

    std::size_t arenaEntries() const   { return mCliqueArena.size(); }
    static constexpr std::int32_t GONE = 2147483647;   // INT32_MAX, above every reachable tag

    bool eliminated(std::int32_t u) const { return mMark[u] == GONE; }


    std::int32_t advanceTag()                    { return ++mTag; }
    std::int32_t mark(std::int32_t v) const      { return mMark[v]; }
    void setMark(std::int32_t v, std::int32_t t) { mMark[v] = t; }

    const std::int32_t* clique(std::int32_t c) const {
        return mCliqueArena.data() + mRun[c].sourcePtr;
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

    // THE MAGNITUDE, for the driver, which does not want to know about the sign. The walks inside
    // this class read mWeight directly and test the sign, which is the whole point.
    std::uint32_t weight(std::int32_t u) const {
        const std::int32_t w = mWeight[u];
        return static_cast<std::uint32_t>(w < 0 ? -w : w);
    }

    void reachableSet(std::int32_t u, std::vector<std::int32_t>& reached);


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

    std::vector<std::int32_t>  mCliqueArena;  // every C[c] ever formed, end to end

    std::vector<std::int32_t>  mSuperNext;
    std::vector<std::int32_t>  mSuperLast;
    // SIGNED, WHICH IS AMD_2'S `Nv` AND THE POINT OF THIS FILE. A one dimensional size is normally
    // `std::uint32_t` because it has nothing to stand in for; this one has, so the rule derives the
    // other answer rather than being excepted. See docs/CODING_RULES.md and
    // docs/DESIGN_DECISIONS.md (2026-08-17).
    //
    //     mWeight[v] >  0    live, and not yet taken into the clique being built; the weight
    //     mWeight[v] <  0    taken into it this step, the pivot included; the weight is -mWeight[v]
    //     mWeight[v] == 0    dead, by a live merge or by mass elimination
    //
    // A FOURTH STATE FOR `number()` IS NOT AVAILABLE, which is what keeps mMark alive in this file.
    // INT32_MIN is unreachable by the encoding and looks like a free slot for "numbered", but
    // `orderAscending` reads `mWeight[pivot]` for every pivot and a prepass vertex IS one, so the
    // weight cannot be spent. GONE in mMark stays the answer for that case.
    //
    // No range is lost: a weight is bounded by n and n is capped at MAX_IDX, so the sign bit was
    // never reachable.
    std::vector<std::int32_t> mWeight;
    std::uint32_t              mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

    bool mReverseIncidence = false;

    bool mVendoredListOrder = false;

    bool mLateMassElimination = false;

    // MIRRORS QuotientGraph::mHasNumbered, and it is here for FAITHFULNESS rather than for speed.
    // In this file it is always true by the time anything reads it, mmd having a prepass, so it
    // saves nothing; what it buys is that a fold landing in the shared class lands here too. The
    // moment this copy stops matching, its time column stops being about storage and starts being
    // about whatever drifted. Same obligation Mmd3B carries.
    bool mHasNumbered = false;

    std::vector<std::int32_t> mMark;     // membership scratch, read against mTag
    std::int32_t              mTag = 0;
};

QuotientGraphC::QuotientGraphC(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mRun(colPtr.empty() ? 0 : colPtr.size() - 1),
      mMark(mRun.size() + kPad, NIL) {
    const std::int32_t size = static_cast<std::int32_t>(mRun.size());

    mSource.reserve((colPtr.empty() ? 0 : colPtr.back()) + kPad);
    for (std::int32_t aj = 0; aj < size; ++aj) {
        mRun[aj].sourcePtr = mSource.size();
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSource.push_back(rowIdx[cp]);
        mRun[aj].adjacencySize = static_cast<std::uint32_t>(mSource.size() - mRun[aj].sourcePtr);
    }

    mCliqueArena.reserve((colPtr.empty() ? 0 : colPtr.back()) + kPad);

    // Padded, not resized: the extra elements are never indexed, and `size()` reads mRun rather
    // than any of these, so nothing downstream can see the difference.
    mSuperNext.assign(static_cast<std::uint32_t>(size) + kPad, NIL);
    mSuperLast.resize(static_cast<std::uint32_t>(size) + kPad);
    mWeight.assign(static_cast<std::uint32_t>(size) + kPad, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

void QuotientGraphC::reachableSet(std::int32_t u, std::vector<std::int32_t>& reached) {
    // THE SIGN OF THE WEIGHT IS THE MEMBERSHIP MARK, which is `AMD_2`'s `Nv` and the fold this file
    // exists to try on the mmd side. The negation IS the set insertion: `nv > 0` is "not yet
    // emitted" and `mWeight[v] = -nv` is the emit, exactly as `mMark[v] < mTag` and
    // `mMark[v] = mTag` were. It is undone in massEliminate, which walks this same set.
    //
    // AND IT SERVES IN ONE OF THE TWO LOOPS ONLY. See the two comments below; the split is the
    // whole finding, and it is genmmd's own arrangement.
    // NO ++mTag HERE. This walk stopped stamping mMark when the sign took over, so the increment
    // fed nothing and only burned tag values. mTag itself STAYS, unlike in Amd3B: reachableWeight
    // still marks against it and the driver still draws from it through advanceTag.
    mWeight[u] = -mWeight[u];              // never its own neighbor; was mMark[u] = mTag
    const VertexRun&    run           = mRun[u];          // one fetch; see the member
    const std::int32_t* source        = mSource.data() + run.sourcePtr;
    const std::uint32_t adjacencySize = run.adjacencySize;
    const std::uint32_t incidenceSize = run.incidenceSize;
    const std::int32_t* incidence     = source + adjacencySize;
    const bool reverse  = mReverseIncidence;
    const bool amdOrder = mVendoredListOrder;

    // ONE LOAD PER CLIQUE MEMBER, and this is where the fold pays. A NUMBERED VERTEX CANNOT BE IN A
    // CLIQUE: the prepass runs to completion before the first elimination, so no clique existed
    // when `number()` was called, and every clique formed since comes from a reach that already
    // skipped it. So `nv > 0` is exactly "live and not yet taken" here, zero covering a vertex a
    // live merge folded away and the sign covering one already inside the new clique.
    //
    // THE ADJACENCY HALF CANNOT HAVE IT, which is the other side of the same fact. `number()`
    // leaves a prepass vertex at weight one deliberately, so that its neighbors' degrees still
    // count it, and leaves it in the adjacency of every one of them. A positive weight therefore
    // does not mean live there, and GONE is what says otherwise. Two arrays, unavoidably.
    //
    // genmmd reaches the same arrangement from the other side: `qsize != 0` inside element walks
    // only, `marker` everywhere else.
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



// DELIBERATELY NOT FOLDED, and this is the one place the encoding was declined rather than blocked.
// The sign could carry liveness and dedup in the clique loop below exactly as it does in
// reachableSet. What it cannot do is clean up after itself: `++mTag` retires every mark in the
// array at once and costs nothing, where every negation has to be undone one at a time, and this
// function is called PER REFRESHED VERTEX with no later traversal of its result to hide a restore
// in. reachableSet is called once per pivot and massEliminate walks its result anyway, which is why
// the fold is affordable there and not here. Folding it needs a scratch list of what was negated,
// which is a pass this function exists to avoid; see the header's note on why it counts without
// materializing the set.
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
        const std::int32_t* members     = mCliqueArena.data() + mRun[c].sourcePtr;
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
    if (mCliqueArena.capacity() - mCliqueArena.size() < size())
        mCliqueArena.reserve(std::max(2 * mCliqueArena.capacity(), mCliqueArena.size() + size()));

    const std::size_t cliqueStart = mCliqueArena.size();
    reachableSet(pivot, mCliqueArena);          // appends; see its note
    const std::uint32_t cliqueLen = static_cast<std::uint32_t>(mCliqueArena.size() - cliqueStart);

    const std::int32_t* reached     = mCliqueArena.data() + cliqueStart;
    const std::uint32_t reachedSize = cliqueLen;

    const std::int32_t* absorbedCliques = mSource.data() + mRun[pivot].sourcePtr + mRun[pivot].adjacencySize;
    const std::uint32_t absorbedSize    = mRun[pivot].incidenceSize;
    for (std::uint32_t i = 0; i < absorbedSize; ++i)
        mRun[absorbedCliques[i]].adjacencySize = 0;   // dead, its block left behind

    mRun[pivot].sourcePtr     = cliqueStart;
    mRun[pivot].adjacencySize = cliqueLen;

    // NO STAMPING PASS, AND NO TAG. Membership was written by the walk, in the sign of the weight,
    // so this only sums. The `inClique` out-parameter and the `++mTag` feeding it went on
    // 2026-08-17, the prune having started reading the sign; matching production's signature is not
    // a reason for a private copy to carry a parameter it never uses.
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t k = 0; k < reachedSize; ++k)
        cliqueWeight += static_cast<std::uint32_t>(-mWeight[reached[k]]);
    mCliqueWeight = cliqueWeight;
}

const std::vector<std::int32_t>& QuotientGraphC::eliminate(std::int32_t pivot) {
    beginElimination(pivot);
    const std::int32_t* reached     = mCliqueArena.data() + mRun[pivot].sourcePtr;
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
            // ONE LOAD, THREE QUESTIONS. A negative weight is a member of the new clique, the
            // pivot included, so the explicit `v == pivot` test goes with the membership test; a
            // zero is a vertex a live merge folded away. The FOURTH question, whether v was
            // numbered by the prepass, still needs mMark: `number()` leaves such a vertex at
            // weight one and in this very list. Dropping that test is not an option, since
            // mass elimination reads `adjacencySize == 0` and a numbered leftover would suppress
            // a merge that used to fire.
            const std::int32_t v = source[k];
            if (mWeight[v] <= 0) continue;             // in the new clique, the pivot, or merged
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by the prepass; see above
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
        // NOT REACHED FROM Mmd3, which never sets the flag, and correct if it ever is. The restore
        // rides in massEliminate's walk over C[pivot]; with that walk deferred to the caller, the
        // signs would otherwise still be negative when the caller next reads a weight.
        const std::int32_t* reached     = mCliqueArena.data() + mRun[pivot].sourcePtr;
        const std::uint32_t reachedSize = mRun[pivot].adjacencySize;
        mWeight[pivot] = -mWeight[pivot];
        for (std::uint32_t ri = 0; ri < reachedSize; ++ri)
            mWeight[reached[ri]] = -mWeight[reached[ri]];
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
    const std::int32_t* reached     = mCliqueArena.data() + mRun[pivot].sourcePtr;
    const std::uint32_t reachedSize = mRun[pivot].adjacencySize;

    // THE SIGNS COME BACK HERE, IN A PASS THAT ALREADY EXISTS, which is what makes the encoding
    // affordable on this branch at all. Amd puts the restore in its bound pass; mmd has no bound
    // pass, and this loop over C[pivot] is the only other traversal of that same set. The pivot
    // goes first, since the merge below adds into it and both operands must be magnitudes by then.
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
        // THE MERGE HAPPENS FIRST, so that the compaction can read the ZERO WEIGHT it leaves
        // rather than a stamp of its own. The old order was the reverse and needed a tag pass over
        // `merged` plus a mark read per member; the weight says the same thing and the
        // supervariable bookkeeping had to write it anyway.
        //
        // NO OTHER MEMBER OF C[pivot] CAN READ ZERO, which is what makes the test exact: a vertex a
        // live merge folded away is left at weight zero but is also stamped GONE, so no reach ever
        // emits it into a clique again.
        for (std::int32_t u : merged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            mCliqueWeight -= static_cast<std::uint32_t>(mWeight[u]);
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }

        std::int32_t*     members     = mCliqueArena.data() + mRun[pivot].sourcePtr;
        const std::uint32_t membersSize = mRun[pivot].adjacencySize;
        std::uint32_t       kept        = 0;
        for (std::uint32_t k = 0; k < membersSize; ++k)
            if (mWeight[members[k]] != 0) members[kept++] = members[k];
        mRun[pivot].adjacencySize = kept;
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
