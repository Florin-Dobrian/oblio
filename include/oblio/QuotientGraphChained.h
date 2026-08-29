#pragma once

// QuotientGraphChained.h - the quotient graph on a CHAINED CLIQUE STORE. Same algorithm and same
// encoding as `QuotientGraph`; the difference is where a clique's members live, and pricing that
// difference is the whole reason this class exists.
//
// A clique lives in the DEAD SEGMENT of the pivot that formed it, chaining into further dead
// segments when one will not hold it, so no clique store is allocated at all. What that costs is a
// link test on every read of every clique, forever, and no headroom reduces it, because chaining
// exists precisely to need none.
//
// ITS OBLIGATION IS TO STAY ENCODING-IDENTICAL TO `QuotientGraph`, so that the only difference
// between them is storage. A fold that lands there lands here, or the comparison quietly stops
// being about storage.

#include "oblio/QuotientGraph.h"   // Buckets, which is shared verbatim
#include "oblio/Types.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <functional>
#include <vector>

namespace Oblio {

class QuotientGraphChained {
public:
    // Built straight from A's pattern. A is stored with both triangles, so a column's rows are
    // already its neighbors, and the only conversion is dropping the diagonal, which is a self
    // loop and says nothing about fill. Oblio's other input assumptions (sorted, unique, both
    // triangles, a structurally present diagonal) hold by construction, which is why nothing here
    // symmetrizes, deduplicates or sorts. See the pass-5 discussion in
    // experiments/ordering/README.md.
    QuotientGraphChained(const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx);

    std::size_t size() const { return mSize; }
    // GONE: one value reserved above every tag, so the stamp array answers "is v dead" on the load
    // it was making anyway and no array is spent on liveness at any walk site.
    //
    // NOT the weight. `mWeight[v] != 0` is a PARTIAL flag: the prepass leaves a numbered vertex at
    // weight one deliberately, so its neighbors' degrees still count it. A zero weight is tested
    // inside clique walks only, where a numbered vertex cannot appear.
    static constexpr std::int32_t GONE =
        std::numeric_limits<std::int32_t>::max();   // above every reachable tag

    bool eliminated(std::int32_t u) const { return mMark[u] == GONE; }

    // The driver stamps into this same array rather than allocating a second one: the eliminator
    // stamps at one level and the refresh at another, one array and one counter for both. One
    // counter is what makes it safe, so two
    // tags can never collide, and a comparison against a captured tag means what it says.
    std::int32_t advanceTag()                          { return ++mTag; }
    std::int32_t mark(std::int32_t u) const            { return mMark[u]; }
    void setMark(std::int32_t u, std::int32_t tag)     { mMark[u] = tag; }

    // The members of clique c, which after eliminating p is the pattern of p's column of L: the
    // vertices the pivot reached, less those it absorbed. A pointer and a length, as the
    // adjacency is, for the same reason.
    //
    // **Valid until the next elimination.** The members live in one arena that grows as cliques are
    // formed, so a reallocation moves them; the offsets are indices and survive it, but a pointer
    // taken beforehand does not. Read a clique at the moment of use, which is the same rule the
    // numeric factor's blocks live by.
    const std::int32_t* clique(std::int32_t c) const {
        return mSrc.data() + mSegment[c].srcPtr;
    }
    // NO CLIQUE SIZE ARRAY. A clique ends at a TERMINATOR value, or at the end of its last segment
    // when the members fill it exactly, and the walk carries a stop condition for each. The
    // terminator has to be a value no member and no link can take, so zero will not do: ids are
    // 0-based and vertex 0 is a real member. A link is `-(c + 1)` and so
    // lies in [-n, -1], which leaves the bottom of the range.
    static constexpr std::int32_t TERMINATOR = std::numeric_limits<std::int32_t>::lowest();

    // Every read of a clique goes through here, because a clique is not one flat run. The walk
    // runs to the end of the current segment, follows any link it meets, and stops at the
    // terminator. TWO STOP CONDITIONS AND BOTH ARE NEEDED: the terminator ends a clique that left
    // room, the segment end ends one whose members fill it exactly and so had nowhere to put one.
    //
    // Templated on the body so it inlines: these are the hottest loops in the ordering.
    template <class F>
    void forEachMember(std::int32_t c, F f) const {
        std::size_t p   = mSegment[c].srcPtr;
        std::size_t end = mSegment[c + 1].srcPtr;
        while (p != end) {
            const std::int32_t v = mSrc[p];
            if (v >= 0) { f(v); ++p; continue; }
            if (v == TERMINATOR) return;
            const std::int32_t d = -v - 1;             // a link, to clique d's segment
            p   = mSegment[d].srcPtr;
            end = mSegment[d + 1].srcPtr;
        }
    }

    // |C[pivot]| WEIGHTED, accumulated in the walk that builds the clique, which has the member in
    // hand already.
    //
    // VALID UNTIL THE NEXT ELIMINATION, and a scalar for that reason: read immediately after the
    // eliminator and never afterwards.
    //
    // Over the UNTRIMMED clique. A driver that mass-eliminates late and needs the trimmed figure
    // recomputes it.
    std::uint32_t cliqueWeight() const { return mCliqueWeight; }

    // PEAK LIVE CLIQUE MEMBERS, for the cross-driver check in tests/test_order.cpp and in
    // benchmarks/matrices. The three mmd drivers return the same permutation, so they form the
    // same cliques and lose the same members at the same moments; this figure MUST be equal across
    // the three however differently they store them.
    //
    // IT IS THE NOTIONAL COUNT, NOT THIS FILE'S STORED ONE, and that is deliberate. `MmdFlat`
    // drops the mass-eliminated from C[pivot] and this file does not, so the size tracked here is
    // one the storage does not have. Reading it as a description of the chained store is a
    // mistake.
    //
    // IT COSTS AN ARRAY that the flat drivers do not pay, this file keeping no clique length at
    // all. INSTRUMENTATION rather than mechanism, and present in release because the benchmark
    // that reads it builds with NDEBUG.
    std::size_t numPeakCliqueMembers() const { return mNumPeakCliqueMembers; }

    // The two halves of the neighbor relation, read by an approximate degree, which decomposes
    // reach(u) rather than forming it.
    //
    // A pointer and a length rather than a container, because that is what the storage is: the two
    // lists share one run, A[u] first and I[u] immediately behind it.
    const std::int32_t* adjacencyMmd(std::int32_t u) const {
        return mSrc.data() + mSegment[u].srcPtr;
    }
    std::uint32_t adjacencySize(std::int32_t u) const { return mSegment[u].adjacencySize; }

    const std::int32_t* incidenceMmd(std::int32_t u) const {
        return mSrc.data() + mSegment[u].srcPtr + mSegment[u].adjacencySize;
    }
    std::uint32_t incidenceSize(std::int32_t u) const { return mSegment[u].incidenceSize; }

    // How many original vertices u stands for. One until mass elimination merges into it, so a
    // degree that counts vertices has to count this rather than entries.
    //
    // Held in a flat array, and this is the condition under which the prototypes said it would
    // have to be: its members are a chain rather than a list, so the size is no longer free to
    // read. The array is here for the chain and not for locality.
    // The MAGNITUDE, for the driver, which does not want to know about the sign. The walks inside
    // this class read mWeight directly and test it, which is the whole point.
    std::uint32_t weight(std::int32_t u) const {
        const std::int32_t w = mWeight[u];
        return static_cast<std::uint32_t>(w < 0 ? -w : w);
    }


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
    const std::vector<std::int32_t>& eliminateMmd(std::int32_t pivot);




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

    // Walk I[u] from the back in the reachable-set walk. A
    // tie-break
    // convention and nothing else: it changes which permutation comes out, never which sets are
    // computed. See the member's note.
    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    // Lay the lists out cliques-first, which is the opposite of this file's default on both
    // counts. Two conventions under one switch because they are one fact, and like the flag above
    // it changes which permutation comes out and never which sets are computed.
    //
    //   reachableSet   walks the CLIQUES before the explicit adjacency
    //   the prune      puts the new clique at the FRONT of I[u] rather than appending, with the
    //                  displaced entries ROTATED rather than shifted
    void setAmdListOrder(bool on) { mAmdListOrder = on; }

    // Stop the eliminator at the prune, leaving mass elimination to the caller. With this on the
    // eliminator returns an EMPTY merged list and C[pivot] is reach(pivot) exactly, and the caller
    // MUST call massEliminate once it has absorbed.
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

    // Expand a pivot sequence into an elimination order over the original vertices. A pivot stands
    // for its whole supervariable, whose members are eliminated consecutively, so this is where a
    // supervariable of size w becomes w columns, with the members in ASCENDING VERTEX INDEX rather
    // than merge order. Indistinguishable members, so only the permutation turns on that; the fill
    // and the forest do not.
    //
    // NO `orderAsMerged` HERE, where the other two classes carry the pair. That one emits merge
    // order and the amd drivers are its only callers, and no amd driver uses this store.
    std::vector<std::int32_t> orderAscending(const std::vector<std::int32_t>& pivots) const;

private:
    // The head and tail of an elimination, shared by the two overloads above and private because
    // between them the graph is half eliminated: the clique is written and stamped but the reached
    // vertices still name the pivot as a variable. Nothing outside may observe that state, which is
    // why the seam is two private calls rather than a public begin and end.
    void beginElimination(std::int32_t pivot, std::int32_t& absorbed);

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
    // at least what I[u] gains, so the incidence write always trails the read cursor.
    std::vector<std::int32_t> mSrc;   // every A[u] then I[u], run after run

    // ONE OBJECT PER VERTEX, NOT THREE ARRAYS. The three numbers are never useful apart: any walk
    // of u needs where its run starts and at least one length. The shared `QuotientGraph` carries
    // the same struct; read its member for the reasoning.
    struct Segment {
        // WHERE u'S SEGMENT STARTS in `mSrc`, fixed at construction and never moved. It is also
        // where the clique formed by u begins, once u has eliminated.
        std::size_t   srcPtr;
        // A[u]'S LENGTH, and zero means the adjacency is EMPTY, which `massEliminate`'s first
        // conjunct tests. It is zeroed at elimination beside `incidenceSize`, so on a clique it
        // carries nothing and is not the member count the other two classes keep here.
        std::uint32_t adjacencySize;
        // I[u]'S LENGTH, immediately behind A[u], and zeroed at elimination with the other.
        std::uint32_t incidenceSize;
    };
    std::vector<Segment> mSegment;

    // A clique's members live in the dead segment of the pivot that formed it, chained into
    // further dead segments when one does not hold them. No store of its own, which is the whole
    // of what this layout buys and what it is here to price.
    //
    // NO CLIQUE LENGTH ANYWHERE. A walk ends at a terminator or at the segment's end, which is why
    // every read goes through the walker above.
    std::vector<std::int32_t> mAbsorbed;   // I[pivot], copied before its segment is overwritten

    // The supervariable a vertex stands for, as a chain rather than a list per vertex:
    // `mSuperNext` links the members and `mSuperLast` names the tail, so a merge is O(1).
    std::vector<std::int32_t> mSuperNext;
    std::vector<std::int32_t> mSuperLast;
    // SIGNED, MIRRORING QuotientGraph. A one dimensional size is normally unsigned
    // because it has nothing to stand in for; this one has. Positive is the
    // weight, negative means already taken into the clique being built, zero means dead, so one
    // load answers what two arrays answered. No range is lost, a weight being bounded by n.
    //
    // IT IS HERE BECAUSE THIS FILE'S OBLIGATION IS TO STAY ENCODING-IDENTICAL TO `MmdFlat` rather
    // than because it pays on its own. Without it the time column stops being the price of the
    // storage and becomes that plus an encoding difference, which is what this file exists not to
    // be.
    std::vector<std::int32_t> mWeight;
    // Mirrors QuotientGraph::mHasNumbered; always true here once the prepass has run.
    bool                      mHasNumbered = false;
    std::uint32_t              mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

    // Which end of I[u] the reachable-set walk starts from. The mmd branch expands the clique seen
    // LAST first, so with a list held in order the walk runs BACKWARD. Same set either way and the
    // same cost, but the order decides C[pivot]'s order, hence which of two equal-degree candidates
    // a later iteration finds first, and minimum degree is settled by exactly that.
    // Off by default, so every existing driver is unaffected; MmdFlat turns it on. See
    // experiments/ordering/mmd3.py, where the same four walks are reversed together.
    bool mReverseIncidence = false;

    // The amd list conventions, off by default so the other drivers are untouched. Read in
    // the reachable-set walk and in the prune, hoisted at both sites for the same reason
    // mReverseIncidence
    // is: a member load the compiler cannot prove is unaliased by the stores in the loop. See the
    // setter for what each half does and why they are one flag.
    bool mAmdListOrder = false;

    // Whether the eliminator stops at the prune and leaves mass elimination to the caller. Off by
    // default. See the setter, and massEliminate() for the half it hands over.
    bool mLateMassElimination = false;

    // TWO STAMP SPACES IN ONE ARRAY: vertices at [v], cliques at [mSize + c]. Clique ids ARE
    // vertex ids, so one space cannot hold both once the vertex half carries GONE: stamping a
    // clique would write a live tag over the slot of the dead pivot that formed it, and a walk of
    // an older clique still holding that pivot as a member would take it for live. Cliques and
    // vertices therefore need marks that cannot be confused, and this is the cheap way to get
    // them. It is what lets the dead-clique test be a stamp again rather than a size, which is what
    // retires mCliqueSize.
    // PEAK LIVE CLIQUE MEMBERS, and this array is the price of having it here. See
    // `numPeakCliqueMembers`.
    std::vector<std::uint32_t> mCliqueLiveMembers;
    std::size_t mNumLiveCliqueMembers = 0;
    std::size_t mNumPeakCliqueMembers = 0;

    std::vector<std::int32_t> mMark;     // membership scratch, read against mTag
    std::int32_t              mTag = 0;
};



inline QuotientGraphChained::QuotientGraphChained(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mSize(colPtr.empty() ? 0 : colPtr.size() - 1),
      mSegment(mSize + 1),
      mCliqueLiveMembers(mSize, 0),
      mMark(2 * mSize, NIL) {

    // One pass to place each column's run, dropping the diagonal. The runs are laid out in column
    // order and never move afterwards, so the offsets are written once here and only the lengths
    // change. The run is u's whole allowance: A[u] fills it now and I[u] grows into the room A[u]
    // gives up, never past mSegment[u + 1].srcPtr.
    mSrc.reserve(colPtr.empty() ? 0 : colPtr.back());
    mSegment[0].srcPtr = 0;
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(mSize); ++aj) {
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSrc.push_back(rowIdx[cp]);
        mSegment[aj + 1].srcPtr = mSrc.size();
        // THE ONE CROSSING. A difference of two positions is a count, so this is the single place
        // in the class where a two-dimensional quantity is written into a one-dimensional one. It
        // is bounded by deg(aj) and so by n, which the SparseMatrix constructor has already capped
        // at MAX_IDX, but the cast is written rather than left implicit because that bound is an
        // argument and not something the types say.
        mSegment[aj].adjacencySize = static_cast<std::uint32_t>(mSegment[aj + 1].srcPtr
                                                            - mSegment[aj].srcPtr);
    }

    // There is no second arena to reserve: cliques live in mSrc, in the segments their
    // pivots vacate. That is the change this file exists to measure.

    // Every vertex begins as a supervariable of one: a chain holding itself, so next is NIL and
    // last is the vertex. Mass elimination splices these together and order() walks them.
    mSuperNext.assign(mSize, NIL);
    mSuperLast.resize(mSize);
    mWeight.assign(mSize, 1);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(mSize); ++u) mSuperLast[u] = u;
}

inline std::uint32_t QuotientGraphChained::reachableSetWeight(std::int32_t u) {
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
    // NO `cClique` POINTER HERE, where the other two classes take one. A clique's members are a
    // chain of segments joined by links, so there is no contiguous block to point at and no
    // member count to hoist; `forEachMember` is what walks it.
    for (std::uint32_t ck = 0; ck < uIncidenceSize; ++ck) {
        const std::int32_t c = uIncidence[ck];
        forEachMember(c, [&](std::int32_t v) {
            if (mMark[v] < mTag) {   // includes mMark[v] != GONE, GONE sorting above every tag
                mMark[v] = mTag; totalWeight += static_cast<std::uint32_t>(mWeight[v]);
            }
        });
    }
    return totalWeight;
}

inline void QuotientGraphChained::number(std::int32_t u) {
    mHasNumbered = true;     // see the member: what makes the GONE test worth asking
    mMark[u]     = GONE;     // a numbered vertex lingers in lists; GONE is what filters it
}

inline void QuotientGraphChained::beginElimination(std::int32_t pivot, std::int32_t& absorbed) {
    // The reach is written STRAIGHT INTO THE SEGMENT. C[pivot] is the reach, so the block the walk
    // fills is the clique's own block and no scratch is needed.
    //
    // A CLIQUE LIVES IN ITS OWN PIVOT'S DEAD SEGMENT, which is what this layout exists to price:
    // the pivot's run is finished with the moment the walk is past it, so the members go there and
    // no second store is allocated at all. A clique too long for one segment CHAINS into the next
    // dead one through a negative link, and a terminator ends it where it leaves room.
    const bool reverse = mReverseIncidence;

    const std::size_t   base = mSegment[pivot].srcPtr;
    const std::uint32_t adjN = mSegment[pivot].adjacencySize;
    const std::uint32_t incN = mSegment[pivot].incidenceSize;

    // I[pivot] must be copied out before the segment holding it is overwritten. A small vector kept
    // for its capacity costs no allocation after the first few pivots and keeps the walk order
    // explicit.
    mAbsorbed.clear();
    for (std::uint32_t ck = 0; ck < incN; ++ck)
        mAbsorbed.push_back(mSrc[base + adjN + (reverse ? incN - 1 - ck : ck)]);

    // THE SIGN OF THE WEIGHT IS THE MEMBERSHIP MARK, mirroring QuotientGraph. The negation IS the
    // insertion, and it is undone in massEliminate, which walks this same set. The pivot is
    // negated so the walk cannot take it into its own clique.
    //
    // The adjacency loop keeps the GONE test, guarded: `number()` leaves a prepass vertex at
    // weight one and in every neighbour's adjacency, so a positive weight does not mean live
    // there. A clique cannot hold one, the prepass completing before the first elimination, so
    // the clique loop asks nothing else.
    ++mTag;
    mWeight[pivot] = -mWeight[pivot];          // never its own neighbor

    std::size_t   rl    = base;                        // write cursor
    std::size_t   rm    = mSegment[pivot + 1].srcPtr - 1;   // last entry of the segment filled

    // One member written, following a link first if the segment is full. The link can only be
    // there to be followed: it is written before the walk that can reach it.
    const auto emit = [&](std::int32_t v) {
        while (rl >= rm) {
            const std::int32_t l = mSrc[rm];            // -(c + 1); see the link encoding above
            rl = mSegment[-l - 1].srcPtr;
            rm = mSegment[-l].srcPtr - 1;
        }
        mSrc[rl++] = v;
    };
    // THE LIVE NEIGHBORS FIRST, AND WITHOUT THE BOUND CHECK. They are a subset of A[pivot], which
    // was read from this same segment, so the cursor cannot pass the reader and cannot leave the
    // segment: at most adjN entries are written into a segment holding adjN + incN. The first loop
    // is unchecked for exactly this reason, and checking here is not merely wasteful, it is wrong:
    // with incN == 0 the cursor legitimately lands on the last entry, and a check would read it as
    // a link that was never written.
    std::uint32_t born = 0;                     // the clique's size; see numPeakCliqueMembers
    for (std::uint32_t vk = 0; vk < adjN; ++vk) {
        const std::int32_t v  = mSrc[base + vk];
        const std::int32_t vWeight = mWeight[v];
        if (vWeight > 0 && !(mHasNumbered && mMark[v] == GONE)) {
            mWeight[v] = -vWeight;
            mSrc[rl++] = v;
            ++born;
        }
    }

    // Then the cliques, each writing its continuation before it can be reached. This loop is the
    // one that can outgrow the segment, and it is the only one that follows a link.
    for (const std::int32_t c : mAbsorbed) {
        // c DIES HERE, absorbed into the clique being built, and this is where its members leave
        // the live count. The size is read rather than derived: walking c to count it would double
        // the cost of the only loop that matters, and nothing else in this file carries a clique's
        // length. See numPeakCliqueMembers.
        mNumLiveCliqueMembers -= mCliqueLiveMembers[c];
        mCliqueLiveMembers[c]  = 0;                 // so a second sighting subtracts nothing

        mSrc[rm] = -(c + 1);                        // the continuation, before it can be read
        forEachMember(c, [&](std::int32_t v) {
            const std::int32_t vWeight = mWeight[v];        // one load; see the note above the walk
            if (vWeight > 0) { mWeight[v] = -vWeight; emit(v); ++born; }
        });
    }

    // THE TERMINATOR. The cursor stops one
    // short of the segment end whenever a link was needed there, so there is room; the one case
    // with no room is a clique whose members fill its last segment exactly, and the walk's second
    // stop condition covers that.
    // Against `rm`, the CURRENT segment's last entry, not the pivot's: the cursor has followed
    // every link the emit needed and is wherever that left it. Comparing against the pivot's own
    // segment end instead leaves a chain
    // unterminated and the walk runs off into whatever the next segment holds, which on a 2x2
    // grid is an immediate hang.
    if (rl <= rm) mSrc[rl] = TERMINATOR;

    // The absorbed cliques die only now: they were being READ until the loop above finished, and
    // their segments now hold part of the new clique, so nothing can be written into them to say
    // so. The stamp below is what says it.

    // NO STAMPING PASS AND NO `inClique`. Membership is written by the walk, in the sign of the
    // weight, so there is nothing to stamp and nothing for a tag to say. Only `absorbed` survives,
    // and it marks the CLIQUE half of mMark, which this file alone keeps.
    //
    // The pivot reads negative from the negation above, so it reads as a member of its own clique
    // and the prune's `mWeight[v] <= 0` drops it with the rest.
    mCliqueLiveMembers[pivot] = born;
    mNumLiveCliqueMembers    += born;
    mNumPeakCliqueMembers     = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);

    ++mTag;
    absorbed = mTag;
    const std::size_t cliqueBase = mSize;
    for (const std::int32_t c : mAbsorbed) mMark[cliqueBase + static_cast<std::size_t>(c)] = absorbed;
}

inline const std::vector<std::int32_t>& QuotientGraphChained::eliminateMmd(std::int32_t pivot) {
    std::int32_t absorbed = NIL;
    beginElimination(pivot, absorbed);
    const std::size_t cliqueBase = mSize;
    // C[pivot] IS the reach here: beginElimination wrote it there and nothing has trimmed it yet
    // (massEliminate does, and runs after). So the prune walks the arena block rather than a copy.
    const bool amdOrder = mAmdListOrder;      // hoisted, as the other flags are

    // Both lists are compacted in place rather than rebuilt into a scratch: every pass here only
    // ever removes, so survivors are written over entries already read and nothing is allocated.
    //
    // Both lists live in u's one run, the adjacency first and the incidence into whatever the
    // adjacency has given up. The pivot is APPENDED and the two boundary entries swapped
    // afterwards; it cannot be written first, because the write cursor starts at `kept` and the
    // read at the original adjacency length, and those are equal whenever nothing was pruned.
    forEachMember(pivot, [&](std::int32_t u) {
        std::int32_t*     source        = mSrc.data() + mSegment[u].srcPtr;
        // The two counters are one-dimensional COUNTS, positions in a list bounded by deg(u) and
        // so by n, where std::size_t is for a position into an n x n object. They take the type of
        // what they count, so they move with the array and no cast is needed here.
        // `heldVertex` is a VERTEX and stays std::int32_t, carrying NIL.
        const std::uint32_t adjacencySize = mSegment[u].adjacencySize;
        std::uint32_t       kept          = 0;
        std::int32_t        heldVertex    = NIL;        // the first survivor, appended last
        for (std::uint32_t vk = 0; vk < adjacencySize; ++vk) {
            // ONE LOAD, THREE QUESTIONS. A negative weight is a member of the new clique, the
            // pivot included, so the explicit pivot test goes with the membership test; a zero is
            // a vertex a live merge folded away. The fourth question, whether v was numbered by
            // the prepass, still needs mMark.
            const std::int32_t v = source[vk];
            if (mWeight[v] <= 0) continue;             // in the new clique, the pivot, or merged
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by a prepass
            if (amdOrder && heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mSegment[u].adjacencySize = kept;                      // A[u] - C[pivot] - {pivot}

    // The read cursor is a hoisted POINTER rather than source[adjacencySize + i]. The read and
    // the write are into the same buffer, so every conditional store orders the next load behind
    // it; that is the price of compacting in place rather than into a scratch.
        const std::int32_t* incidence     = source + adjacencySize;
        const std::uint32_t incidenceSize = mSegment[u].incidenceSize;
        std::uint32_t       write         = kept;      // follows `kept`; see the note above
        for (std::uint32_t ck = 0; ck < incidenceSize; ++ck) {
            const std::int32_t c = incidence[ck];
            if (mMark[cliqueBase + static_cast<std::size_t>(c)] != absorbed)
                source[write++] = c;
        }
        source[write++]   = pivot;                     // u joins the new clique, id = pivot
        mSegment[u].incidenceSize = write - kept;
        // [c1, ..., ck, pivot] to [pivot, c2, ..., ck, c1], which is a SWAP of the two boundary
        // entries: only they move and c2..ck stay put. The pivot cannot be written first, which
        // is what the loop above would otherwise allow: the write cursor starts at `kept` and the
        // read at the original adjacencySize, and those are equal whenever nothing was pruned from
        // A[u], so an extra write before the reads finish clobbers an unread entry. The three
        // assignments come after both compactions for exactly this reason.
        if (amdOrder && write - kept > 1) std::swap(source[kept], source[write - 1]);
    });

    return finishElimination(pivot);
}

// The same prune, with the driver's first scan folded into the two loops. The header carries the
// argument for why this is one call and why it cannot be folded further; the loops below are the
// plain ones with an accumulation added on each survivor, and nothing else differs.
// The same prune again, with the driver's first scan folded in under the tagged-workspace
// encoding. This is the overload above with `outside`, `mark` and `tag` replaced by one array and
// one tag, plus the hash key's adjacency half, which a driver that folds its scan in here has no
// other walk of A[u] to accumulate. Everything the header says about why the fold is sound holds
// unchanged: the scan runs over the untrimmed C[p], and the weights over A[u] are read before mass
// elimination could move any of them.
inline const std::vector<std::int32_t>&
QuotientGraphChained::finishElimination(std::int32_t pivot) {
    // Under mLateMassElimination the merge is the caller's, run after it has absorbed, so this
    // hands back an empty list and C[pivot] stays reach(pivot) exactly. See the setter.
    if (mLateMassElimination) {
        mMerged.clear();
    } else {
        massEliminate(pivot);
    }

    mSegment[pivot].adjacencySize = 0;
    mSegment[pivot].incidenceSize = 0;
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
// absorption is what makes this cheap test agree with the true one.
inline const std::vector<std::int32_t>& QuotientGraphChained::massEliminate(std::int32_t pivot) {
    mMerged.clear();   // a member scratch, kept for its capacity
    // THE SIGNS COME BACK HERE, IN A PASS THAT ALREADY EXISTS, mirroring QuotientGraph. The walk
    // below is the only other traversal of C[pivot], so the restore rides in it and costs no pass.
    // The pivot goes first, since the merge at the end adds into it and both operands must be
    // magnitudes by then. Every path through an elimination reaches this function: no mmd driver
    // sets late mass elimination.
    mWeight[pivot] = -mWeight[pivot];
    // Walks C[pivot], which is the full reach and STAYS the full reach: see the note below on why
    // nothing here shortens it.
    forEachMember(pivot, [&](std::int32_t u) {
        mWeight[u] = -mWeight[u];                          // live again, and positive
        // Under mAmdListOrder the new clique sits at the FRONT of I[u] rather than the back,
        // so the single remaining entry is at the head of the incidence run either way: with A[u]
        // empty the run starts with I[u], and with one clique there is only one position. The
        // test therefore needs no branch on the flag.
        if (mSegment[u].adjacencySize == 0 && mSegment[u].incidenceSize == 1 &&
            mSrc[mSegment[u].srcPtr] == pivot) {         // A[u] empty, so I[u] starts at the run
            mSegment[u].incidenceSize = 0;
            mMark[u]          = GONE;
            mMerged.push_back(u);
        }
    });
    // NO COMPACTION OF C[pivot]. A clique is placed once and never shortened: mass elimination
    // zeroes the mMerged vertex's weight and leaves it in the list, and every later reader skips it
    // on that.
    //
    // EVERY READER ALREADY SKIPS THE DEAD, which is what makes this safe: `reachableSet` tests
    // `nv > 0`, the reach count tests the mark with GONE outranking any tag, the driver's clique
    // walk tests `eliminated`, and its pair test rejects on the tag. The one loop with no test is
    // the eviction over C[pivot], and the unfile-and-restore pair is idempotent, so a mMerged
    // vertex is evicted harmlessly.
    //
    // THE COST IS SPACE. A clique keeps its dead members for as long as it lives, so live clique
    // storage here is strictly above what the compacting classes report for the same ordering.
    mCliqueLiveMembers[pivot] -= static_cast<std::uint32_t>(mMerged.size());
    mNumLiveCliqueMembers     -= mMerged.size();

    if (!mMerged.empty()) {
        for (std::int32_t u : mMerged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
    // The weighted clique size follows the clique, so `cliqueWeight()` stays true across the
    // merge. Magnitudes only: the sign of a weight is membership of the clique being built.
    //
    // THE MERGED LEAVE THE LIVE COUNT even though they do NOT leave this file's storage, a mMerged
    // vertex keeping its place and being skipped on a zero weight. Tracking the notional size is
    // what makes the figure comparable with `MmdFlat`, which does drop them.
            mCliqueWeight -= static_cast<std::uint32_t>(mWeight[u]);
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }
    }
    return mMerged;
}

inline void QuotientGraphChained::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u]                += mWeight[v];
    mWeight[v]                = 0;

    mSegment[v].adjacencySize = 0;
    mSegment[v].incidenceSize = 0;
    mMark[v]                  = GONE;
}

inline std::vector<std::int32_t> QuotientGraphChained::orderAscending(
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

} // namespace Oblio
