#include <algorithm>
#include "oblio/QuotientGraph.h"

#include <utility>

namespace Oblio {

QuotientGraph::QuotientGraph(const std::vector<std::size_t>&  colPtr,
                             const std::vector<std::int32_t>& rowIdx)
    : mRun(colPtr.empty() ? 0 : colPtr.size() - 1),
      mCliquePtr(mRun.size(), 0),
      mCliqueSize(mRun.size(), 0),
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
    ++mTag;
    mMark[u] = mTag;                       // never its own neighbor
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
            const std::int32_t* members     = mCliqueArena.data() + mCliquePtr[c];
            const std::uint32_t membersSize = mCliqueSize[c];
            for (std::uint32_t k = 0; k < membersSize; ++k) {
                const std::int32_t v = members[k];
                if (mMark[v] < mTag) {
                    mMark[v] = mTag;
                    reached.push_back(v);
                }
            }
        }
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (mMark[v] < mTag) {
                mMark[v] = mTag;
                reached.push_back(v);
            }
        }
        return;
    }

    for (std::uint32_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v = source[k];
        if (mMark[v] != GONE) {
            mMark[v] = mTag;
            reached.push_back(v);
        }
    }
    for (std::uint32_t ii = 0; ii < incidenceSize; ++ii) {
        const std::uint32_t i           = reverse ? incidenceSize - 1 - ii : ii;
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mCliqueArena.data() + mCliquePtr[c];
        const std::uint32_t membersSize = mCliqueSize[c];
        for (std::uint32_t k = 0; k < membersSize; ++k) {
            const std::int32_t v = members[k];
            if (mMark[v] < mTag) {
                mMark[v] = mTag;
                reached.push_back(v);
            }
        }
    }
}

std::vector<std::int32_t> QuotientGraph::reachableSet(std::int32_t u) {
    std::vector<std::int32_t> reached;   // empty, so the appending overload needs no clear
    reachableSet(u, reached);
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
        const std::int32_t* members     = mCliqueArena.data() + mCliquePtr[c];
        const std::uint32_t membersSize = mCliqueSize[c];
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
        if (mMark[v] != GONE) { mMark[v] = mTag; reached += mWeight[v]; }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::uint32_t i = 0; i < incidenceSize; ++i) {
        const std::int32_t  c           = incidence[i];
        const std::int32_t* members     = mCliqueArena.data() + mCliquePtr[c];
        const std::uint32_t membersSize = mCliqueSize[c];
        for (std::uint32_t k = 0; k < membersSize; ++k) {
            const std::int32_t v = members[k];
            if (mMark[v] < mTag) { mMark[v] = mTag; reached += mWeight[v]; }
        }
    }
    return reached;
}

void QuotientGraph::number(std::int32_t u) {
    // A numbered vertex lingers in every list that named it, deliberately: its neighbors keep
    // degrees that still count it. GONE is what stops the walks following it back in.
    mMark[u] = GONE;
}

void QuotientGraph::beginElimination(std::int32_t pivot, std::int32_t& inClique) {
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
    // the arena, `mCliqueArena.data() + mCliquePtr[c]`, while appending the reach to that same
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

    mCliquePtr[pivot] = mCliqueArena.size();
    reachableSet(pivot, mCliqueArena);          // appends; see its note
    // THE SECOND AND LAST CROSSING. The arena's new length less this block's start is a member
    // count, so a two-dimensional quantity is written into a one-dimensional one, exactly as the
    // constructor does for the source runs. Bounded by n, a reach having at most n entries, which
    // is the same bound the reserve above relies on.
    mCliqueSize[pivot] = static_cast<std::uint32_t>(mCliqueArena.size() - mCliquePtr[pivot]);

    // Taken AFTER the append, since that is what can move the arena. With the reserve above it
    // cannot have moved, and this stays as it is regardless: it costs nothing and it is the shape
    // that remains correct if the reserve is ever revised.
    const std::int32_t* reached     = mCliqueArena.data() + mCliquePtr[pivot];
    const std::uint32_t reachedSize = mCliqueSize[pivot];

    // The absorbed cliques are I[pivot], read where they lie. Nothing below writes the pivot's run
    // (the prune rewrites the runs of C[pivot]'s members, and the pivot is not one of them), so no
    // copy and no scratch is needed to keep them alive across the passes that follow.
    const std::int32_t* absorbedCliques = mSource.data() + mRun[pivot].sourcePtr + mRun[pivot].adjacencySize;
    const std::uint32_t absorbedSize    = mRun[pivot].incidenceSize;
    for (std::uint32_t i = 0; i < absorbedSize; ++i)
        mCliqueSize[absorbedCliques[i]] = 0;    // dead, its block left behind

    // Stamp the new clique. THE ABSORBED CLIQUES ARE NOT STAMPED: the loop above has just set
    // `mCliqueSize[c] = 0` for every one of them, and a dead clique is exactly a clique of size
    // zero, so the prune's incidence compaction asks `mCliqueSize[c] != 0` and needs no tag and
    // no second pass.
    //
    // That matters beyond the pass it deletes. Clique ids and vertex ids share this array, so
    // stamping a clique wrote a live tag over the slot of the VERTEX that formed it, which is
    // dead. Harmless while the mark carried only "seen this step"; fatal once it also carries
    // "dead", since a dead pivot can still sit as a member of an older clique that is still
    // alive, and a walk of that clique would read the borrowed tag and take it for a live vertex.
    // Neither genmmd nor AMD_2 shares one stamp array between the two kinds.
    ++mTag;
    inClique = mTag;
    // The weighted size of the new clique is accumulated HERE rather than in a pass of its own in
    // each driver: this loop has the member loaded already, so the weight is one more read off a
    // line the stamp is touching anyway. See cliqueWeight().
    std::uint32_t cliqueWeight = 0;
    for (std::uint32_t k = 0; k < reachedSize; ++k) {
        const std::int32_t v = reached[k];
        mMark[v] = inClique;
        cliqueWeight += mWeight[v];
    }
    mCliqueWeight = cliqueWeight;
}

const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot) {
    std::int32_t inClique = NIL;
    beginElimination(pivot, inClique);
    // C[pivot] IS the reach here: beginElimination wrote it there and nothing has trimmed it yet
    // (massEliminate does, and runs after). So the prune walks the arena block rather than a copy.
    const std::int32_t* reached     = mCliqueArena.data() + mCliquePtr[pivot];
    const std::uint32_t reachedSize = mCliqueSize[pivot];
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
            if (v == pivot) continue;                  // no longer a variable
            if (mMark[v] == inClique) continue;        // both ends inside the new clique
            if (mMark[v] == GONE) continue;            // dead for good; see GONE   // numbered by a prepass, gone for good
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
            if (mCliqueSize[c] != 0) source[write++] = c;   // dead is size zero; see above
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

// The same prune, with the driver's first scan folded into the two loops. The header carries the
// argument for why this is one call and why it cannot be folded further; the loops below are the
// plain ones with an accumulation added on each survivor, and nothing else differs.
const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot,
                                                         ApproximateScan& scan) {
    std::int32_t inClique = NIL;
    beginElimination(pivot, inClique);
    // C[pivot] IS the reach here: beginElimination wrote it there and nothing has trimmed it yet
    // (massEliminate does, and runs after). So the prune walks the arena block rather than a copy.
    const std::int32_t* reached     = mCliqueArena.data() + mCliquePtr[pivot];
    const std::uint32_t reachedSize = mCliqueSize[pivot];
    const bool amdOrder = mVendoredListOrder;      // hoisted, as the other flags are

    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
        const std::uint32_t weightU     = mWeight[u];
        // ONE FETCH OF THE RUN, not three. The three numbers share a 16-byte object, so the
        // reference below brings all of them in on one line; read as `mRun[u].sourcePtr` and
        // friends they would still be three subscripts of the same object and the compiler would
        // still load once, but naming it once is what makes that visible to a reader.
        const VertexRun&  run           = mRun[u];
        std::int32_t*     source        = mSource.data() + run.sourcePtr;
        const std::uint32_t adjacencySize = run.adjacencySize;
        std::uint32_t       kept          = 0;
        std::uint32_t       explicitPart  = 0;         // a weight sum, not a count of positions
        std::int32_t        heldVertex    = NIL;       // see the plain prune above
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (v == pivot) continue;
            if (mMark[v] == inClique) continue;
            if (mMark[v] == GONE) continue;            // dead for good; see GONE
            explicitPart += mWeight[v];                // the bound's explicit term, in this visit
            if (amdOrder && heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mRun[u].adjacencySize       = kept;
        scan.explicitPart[u]    = explicitPart;

        const std::int32_t* incidence     = source + adjacencySize;   // hoisted; see the plain prune
        const std::uint32_t incidenceSize = run.incidenceSize;
        std::uint32_t       write         = kept;      // follows `kept`; see the plain prune
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
            if (mCliqueSize[c] == 0) continue;             // dead is size zero; see above
            source[write++] = c;
            if (c == pivot) continue;                  // the new clique subtracts from nothing
            if (scan.mark[c] != scan.tag) {            // first sighting: start from |C[c]|
                scan.mark[c] = scan.tag;
                scan.touchedCliques.push_back(c);
                scan.outside[c] = scan.cliqueDegree[c] - weightU;
            } else {
                scan.outside[c] -= weightU;
            }
        }
        source[write++]   = pivot;
        mRun[u].incidenceSize = write - kept;
        if (amdOrder && write - kept > 1) std::swap(source[kept], source[write - 1]);
    }

    return finishElimination(pivot);
}

// The same prune again, with the driver's first scan folded in under `Amd.cpp`'s tagged-W
// encoding. This is the overload above with `outside`, `mark` and `tag` replaced by one array and
// one tag, plus the hash key's adjacency half, which a driver that folds its scan in here has no
// other walk of A[u] to accumulate. Everything the header says about why the fold is sound holds
// unchanged: the scan runs over the untrimmed C[p], and the weights over A[u] are read before mass
// elimination could move any of them.
const std::vector<std::int32_t>& QuotientGraph::eliminate(std::int32_t pivot, TaggedScan& scan) {
    std::int32_t inClique = NIL;
    beginElimination(pivot, inClique);
    const std::int32_t* reached     = mCliqueArena.data() + mCliquePtr[pivot];
    const std::uint32_t reachedSize = mCliqueSize[pivot];
    const bool          amdOrder    = mVendoredListOrder;
    const std::int32_t  wflg        = scan.wflg;          // hoisted, as the flags are
    const std::int32_t  modulus     = scan.modulus;

    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u       = reached[ri];
        // A SIGNEDNESS cast, not a narrowing one: `mWeight` is already 32 bits and `wnvi` below
        // must be able to go negative, which is Amd.cpp's convention and the reason `w` is signed.
        const std::int32_t nvi     = static_cast<std::int32_t>(mWeight[u]);
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
        std::int32_t        key           = 0;            // reduced as it goes; see the header
        std::int32_t        heldVertex    = NIL;          // see the plain prune above
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            const std::int32_t v = source[k];
            if (v == pivot) continue;
            if (mMark[v] == inClique) continue;
            if (mMark[v] == GONE) continue;            // dead for good; see GONE
            explicitPart += mWeight[v];
            key = (key + v + 1) % modulus;
            if (amdOrder && heldVertex == NIL) { heldVertex = v; continue; }
            source[kept++] = v;
        }
        if (heldVertex != NIL) source[kept++] = heldVertex;
        mRun[u].adjacencySize       = kept;
        scan.explicitPart[u]    = explicitPart;
        scan.key[u]             = key;                    // the ADJACENCY half alone; see the header

        const std::int32_t* incidence     = source + adjacencySize;
        const std::uint32_t incidenceSize = run.incidenceSize;
        std::uint32_t       write         = kept;         // follows `kept`; see the plain prune
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
            if (mCliqueSize[c] == 0) continue;             // dead is size zero; see above
            source[write++] = c;
            if (c == pivot) continue;                     // the new clique subtracts from nothing
            // Amd.cpp's four lines, transcribed. A clique seen earlier in this step already holds
            // the running value above the tag; one seen for the first time starts from |C[c]| and
            // is listed once; one already absorbed reads zero and is left alone.
            std::int32_t we = scan.w[c];
            if (we >= wflg) {
                we -= nvi;
            } else if (we != 0) {
                // A SIGNEDNESS cast, as `nvi` above: `cliqueDegree` is already 32 bits and `wnvi`
                // is negative in the early eliminations.
                we = static_cast<std::int32_t>(scan.cliqueDegree[c]) + wnvi;
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

    mRun[pivot].adjacencySize = 0;
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
    const std::int32_t* reached     = mCliqueArena.data() + mCliquePtr[pivot];
    const std::uint32_t reachedSize = mCliqueSize[pivot];
    for (std::uint32_t ri = 0; ri < reachedSize; ++ri) {
        const std::int32_t u = reached[ri];
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
    if (!merged.empty()) {                             // C[pivot] - merged, one compaction pass
        ++mTag;
        for (std::int32_t u : merged) mMark[u] = mTag;
        std::int32_t*     members     = mCliqueArena.data() + mCliquePtr[pivot];
        const std::uint32_t membersSize = mCliqueSize[pivot];
        std::uint32_t       kept        = 0;
        for (std::uint32_t k = 0; k < membersSize; ++k)
            if (mMark[members[k]] != mTag) members[kept++] = members[k];
        mCliqueSize[pivot] = kept;

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
            mCliqueWeight -= mWeight[u];
            mWeight[pivot] += mWeight[u];
            mWeight[u] = 0;
        }
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

    for (std::int32_t c : cliques) mCliqueSize[c] = 0;   // dead; the compaction reads the size

    for (std::uint32_t k = 0; k < vertexCount; ++k) {  // I[u] - dead, compacted in place
        const std::int32_t u         = vertices[k];
        std::int32_t*      incidence = mSource.data() + mRun[u].sourcePtr + mRun[u].adjacencySize;
        const std::uint32_t size     = mRun[u].incidenceSize;
        std::uint32_t       kept     = 0;
        for (std::uint32_t i = 0; i < size; ++i)
            if (mCliqueSize[incidence[i]] != 0) incidence[kept++] = incidence[i];
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
        pos += mWeight[pivot];                              // the whole supervariable's room
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
