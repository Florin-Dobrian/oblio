#include <algorithm>
#include "oblio/QuotientGraph.h"

#include <utility>

namespace Oblio {

QuotientGraph::QuotientGraph(const std::vector<std::size_t>&  colPtr,
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
void QuotientGraph::reachableSet(std::int32_t u, std::vector<std::int32_t>& reached) {
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
    // `qsize != 0` inside element walks only and `marker` everywhere else.
    ++mTag;
    mWeight[u] = -mWeight[u];              // never its own neighbor
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per element. Measured at 300 ms of AMD1's
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

    // Which source is walked first. genmmd expands the variables and then the elements, which is
    // how the whole md ladder is laid out; AMD_2 takes the elements first and the supervariables
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

std::vector<std::int32_t> QuotientGraph::reachableSet(std::int32_t u) {
    std::vector<std::int32_t> reached;   // empty, so the appending overload needs no clear
    reachableSet(u, reached);
    // AND IT UNDOES THE NEGATION, which the appending overload deliberately leaves for
    // massEliminate. This form is a convenience with no caller inside the elimination, so a reader
    // reaching for it should get a query rather than a half-finished elimination.
    mWeight[u] = -mWeight[u];
    for (std::int32_t v : reached) mWeight[v] = -mWeight[v];
    return reached;
}

std::uint32_t QuotientGraph::reachableSize(std::int32_t u) {
    // The same two passes as reachableSet, counting rather than collecting.
    ++mTag;
    std::uint32_t reached = 0;   // DISTINCT vertices, the mark seeing to that, so at most n
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per element. Measured at 300 ms of AMD1's
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

std::uint32_t QuotientGraph::reachableWeight(std::int32_t u) {
    ++mTag;
    std::uint32_t reached = 0;   // a sum over DISTINCT vertices, so bounded by n; see the header
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per element. Measured at 300 ms of AMD1's
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

void QuotientGraph::number(std::int32_t u) {
    // A numbered vertex lingers in every list that named it, deliberately: its neighbors keep
    // degrees that still count it. GONE is what stops the walks following it back in.
    //
    // AND THE FLAG IS WHAT TELLS THE WALKS TO ASK. See mHasNumbered: this is the only thing that
    // sets it, so a run that never calls this function never pays for the test.
    mHasNumbered = true;
    mMark[u]     = GONE;
}

void QuotientGraph::beginElimination(std::int32_t pivot) {
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
    // The remedy is to make the arena unable to move rather than to re-fetch per element, which
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
        mRun[absorbedCliques[i]].adjacencySize = 0;   // dead, its block left behind

    // Both readers of the pivot's own run are past, so the run becomes the clique's descriptor.
    mRun[pivot].sourcePtr     = cliqueStart;
    mRun[pivot].adjacencySize = cliqueLen;

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

const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot) {
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

const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot, TaggedScan& scan) {
    // A CLIQUE DIES TWO WAYS AND THE TAGGED W MUST LEARN ABOUT BOTH. Aggressive absorption zeroes
    // `w[c]` in the driver; elimination-time absorption, which is I[pivot], happens inside
    // `beginElimination` and used to be recorded only as a zero size. So the prune could test
    // neither alone and tested the size, an array it reads for nothing else, once per incidence
    // element on every reached vertex. `AMD_2` writes both deaths into W, `Pe[e] = FLIP(me)` with
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
    // from degree list" inside CONSTRUCT NEW ELEMENT. Not bookkeeping moved earlier for its own
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
        // `nvi = -Nv [i]` under CONSTRUCT NEW ELEMENT, and the sign is the whole of it: `wnvi`
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
        // add, against a test per element if it were guarded.
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

const std::vector<std::int32_t>& QuotientGraph::finishElimination(std::int32_t pivot) {
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
const std::vector<std::int32_t>& QuotientGraph::massEliminate(std::int32_t pivot) {
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
        // empty the run starts with I[u], and with one element there is only one position. The
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
        mRun[pivot].adjacencySize = kept;
    }
    return merged;
}

void QuotientGraph::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    mRun[v].adjacencySize = 0;
    mRun[v].incidenceSize = 0;
    mMark[v]          = GONE;
}

void QuotientGraph::absorb(const std::vector<std::int32_t>& cliques,
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
std::vector<std::int32_t> QuotientGraph::orderAscending(
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

std::vector<std::int32_t> QuotientGraph::order(const std::vector<std::int32_t>& pivots) const {
    std::vector<std::int32_t> order;
    order.reserve(size());
    for (std::int32_t pivot : pivots)
        for (std::int32_t u = pivot; u != NIL; u = mSuperNext[u]) order.push_back(u);
    return order;
}

} // namespace Oblio
