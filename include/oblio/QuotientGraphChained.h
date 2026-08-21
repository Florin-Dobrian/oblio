#pragma once

// QuotientGraphChained.h - the quotient graph on GENMMD'S CLIQUE LAYOUT: every clique lives in the
// DEAD SEGMENT of the run belonging to the vertex that formed it, and a clique too large for that
// run is CHAINED across the dead segments of several, linked head to tail.
//
// SAME GRAPH, SAME ALGORITHMS, DIFFERENT STORAGE. Every idea in QuotientGraph.h holds here without
// change. Read that file first; this one comments only what the storage makes different, which is
// where a clique's members live and what a walk over them has to follow.
//
// WHAT THE SCHEME BUYS is that no space is ever allocated for cliques at all: a run's dead segment
// is space a vertex has already given up, so the workspace is exactly the pattern and never grows.
// That is a stronger bound than `AMD_2`'s pool, which reserves a fifth extra, and stronger again
// than our arena, which cannot be bounded before the run.
//
// WHAT IT COSTS IS A POINTER CHASE ON THE HOTTEST PATH. A clique's members are not contiguous, so
// every walk over one follows links between dead segments, and the reachable set walks cliques.
// The price is measured and it is large; see docs/DESIGN_DECISIONS.md.
//
// ONE DRIVER, `Mmd3B`, AND NO SUFFIXES. This layout is genmmd's alone and no amd driver uses it,
// so unlike QuotientGraphCompacted.h nothing here is split two ways.
//
// THE BODIES ARE HERE, NOT IN A SOURCE FILE. A driver that calls out of its own translation unit
// reloads the run array's base around every call and spills its pivot loop's registers. Every
// ordering driver is in its own unit with its graph, so that they can be compared with EACH OTHER.
// Every out-of-class definition below is `inline`. See docs/CODING_RULES.md for the rule and the
// mechanics it needs.

#include "oblio/QuotientGraph.h"   // Buckets, which is shared verbatim
#include "oblio/Types.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
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

    std::size_t size() const           { return (mRun.size() - 1); }
    // GONE, and it is genmmd's `marker[v] = maxint` exactly. One value reserved above every tag
    // makes the stamp array answer "is v dead" on the load it was making anyway, which is why
    // genmmd spends no array on liveness at any of its walk sites.
    //
    // NOT the weight: `mWeight[v] != 0` is a PARTIAL flag on both sides. `number()` leaves a
    // prepass vertex at weight one deliberately, so its neighbors' degrees still count it, and
    // genmmd's prepass leaves `qsize` at one for the same reason; it uses `qsize[nd] != 0` only
    // inside clique walks, where a prepass vertex cannot appear because the mark kept it out of
    // every clique. Used as the universal test it lets a numbered vertex back into a reachable
    // set, which is 201 entries for 200 vertices on a random mmd2 pattern.
    //
    // Safe only because clique ids no longer share this array; see beginElimination.
    static constexpr std::int32_t GONE = 2147483647;   // INT32_MAX, above every reachable tag

    bool eliminated(std::int32_t u) const { return mMark[u] == GONE; }

    // The driver stamps into this same array rather than allocating a second one, which is what
    // genmmd's `marker` does: `mmdelm` stamps at level `tag` and `mmdupd` at level
    // `mt = tag + md0`, one array and one counter for both. Ours is one counter likewise, so two
    // tags can never collide, and a comparison against a captured tag means what it says.
    std::int32_t advanceTag()                          { return ++mTag; }
    std::int32_t mark(std::int32_t v) const            { return mMark[v]; }
    void setMark(std::int32_t v, std::int32_t t)       { mMark[v] = t; }

    // The members of clique c, which after eliminate(p) is the pattern of p's column of L: the
    // vertices the pivot reached, less those it absorbed. A pointer and a length, as the
    // adjacency is, for the same reason.
    //
    // **Valid until the next eliminate.** The members live in one arena that grows as cliques are
    // formed, so a reallocation moves them; the offsets are indices and survive it, but a pointer
    // taken beforehand does not. Read a clique at the moment of use, which is the same rule the
    // numeric factor's blocks live by.
    const std::int32_t* clique(std::int32_t c) const {
        return mSource.data() + mRun[c].sourcePtr;
    }
    // NO CLIQUE SIZE ARRAY. A clique ends where genmmd's ends: at a TERMINATOR value, or at the
    // end of its last segment when the members fill it exactly. `mmdelm` writes
    // `if (rl <= rm) adjncy[rl] = 0` and its walk carries both `j <= jt` and `if (nd == 0) break`,
    // for precisely these two cases; ours does the same with one difference. genmmd's ids are
    // 1-based so zero is free to mean "end"; ours are 0-based and vertex 0 is a real member, so
    // the terminator has to be a value no member and no link can take. A link is `-(c + 1)` and so
    // lies in [-n, -1], which leaves the bottom of the range.
    static constexpr std::int32_t TERMINATOR = -2147483647 - 1;   // INT32_MIN

    // Every read of a clique goes through here, because a clique is no longer one flat run. The
    // walk runs to the end of the current segment, follows any link it meets, and stops at the
    // terminator. Two stop conditions rather than one, and both are needed: the terminator ends a
    // clique that left room, the segment end ends one whose members fill it exactly and so had
    // nowhere to put a terminator. `mmdelm`'s `n400` loop carries the same pair.
    //
    // What this replaces is a counter decremented per member, which was a second loop-carried
    // dependency beside the cursor and cost a size-n array to feed.
    //
    // Templated on the body so it inlines: these are the hottest loops in the ordering, and a
    // std::function or a virtual would cost more than the whole scheme saves.
    template <class F>
    void forEachMember(std::int32_t c, F f) const {
        std::size_t p   = mRun[c].sourcePtr;
        std::size_t end = mRun[c + 1].sourcePtr;
        while (p != end) {
            const std::int32_t v = mSource[p];
            if (v >= 0) { f(v); ++p; continue; }
            if (v == TERMINATOR) return;
            const std::int32_t d = -v - 1;             // a link, to clique d's segment
            p   = mRun[d].sourcePtr;
            end = mRun[d + 1].sourcePtr;
        }
    }

    // |C[pivot]| WEIGHTED, which every amd driver needs and all four used to compute for
    // themselves in a pass of their own, one scattered weight load per member per pivot. It is
    // accumulated in `beginElimination`'s stamping walk instead, which has the member in hand
    // already, so the pass is gone and nothing is added. AMD_2 does the same, `degme += nvi`
    // inside the loops that build the clique, at its lines 1492 and 1636.
    //
    // MEASURED AT ZERO, and that is the point of saying so here. CPU Counters before and after,
    // same session: useful cycles unchanged within half a percent in both families, AMD3 still
    // 1.61x the vendored routine on cubes and 1.56x in 2D. The 25 visits per pivot this removes on
    // cubes were contiguous walks over hot arrays and cost near nothing. It is kept as a faithful
    // port that removes two passes and a stale-value hazard, NOT as a speed fix, and it should not
    // be cited as one. `benchmarks/ordering/README.md` (2026-08-09) has the numbers.
    //
    // VALID UNTIL THE NEXT ELIMINATION, and a scalar rather than an array for that reason: it is
    // read immediately after `eliminate` and never afterwards. Same contract as the pointer
    // `clique()` returns, and stated here because a scalar carrying per-pivot state is the kind
    // of thing that goes stale silently.
    //
    // Over the UNTRIMMED clique, which is what the drivers want: mass elimination runs after the
    // absorption in Amd3, and the value that pass needs is the one before any trimming. Amd3
    // recomputes it afterwards over the trimmed clique, which is ledger entry 7 and stays.
    std::uint32_t cliqueWeight() const { return mCliqueWeight; }

    // PEAK LIVE CLIQUE MEMBERS, for the cross-driver check in tests/test_order.cpp and in
    // benchmarks/matrices. `Mmd3`, `Mmd3B` and `Mmd3C` return the same permutation, so they form
    // the same cliques and lose the same members at the same moments; this figure MUST be equal
    // across the three however differently they store them. The digest says the outputs agree,
    // this says the work behind them agreed too, and on the amd side the same check found two
    // defects in a day.
    //
    // IT IS THE NOTIONAL COUNT, NOT THIS FILE'S STORED ONE, which is the odd part and is
    // deliberate. `Mmd3` drops the mass-eliminated from C[pivot]; this file does not, `mmdelm`
    // leaving them in place and skipping them on `qsize != 0`. So the size tracked here is one
    // the storage does not have: it is what the flat drivers hold, maintained so the comparison
    // is possible at all. Reading it as a description of the chained store is a mistake.
    //
    // AND IT COSTS AN ARRAY, which the flat drivers do not pay: they read a clique's size from
    // `mRun[c].adjacencySize`, a descriptor they keep anyway, and this file keeps no clique length
    // at all. INSTRUMENTATION rather than mechanism, and it is present in release because the
    // benchmark that reads it builds with NDEBUG. It should not be counted against this file's
    // argument about how many arrays a vertex's state lives in.
    std::size_t numPeakCliqueMembers() const { return mNumPeakCliqueMembers; }

    // The two halves of the neighbor relation, read by an approximate degree, which decomposes
    // reach(u) rather than forming it: A[u] contributes its own term and each clique in I[u]
    // contributes one number computed for the clique rather than for u. An exact degree has no
    // use for either, since it unites them through reachableSet instead.
    //
    // Both are a pointer and a length rather than a container, because that is what their storage
    // is: the two lists share one run of one flat array. A[u] starts at the run and I[u] starts
    // immediately behind it, which is why the incidence lookup reads the adjacency's length. See
    // the note on the members below.
    const std::int32_t* adjacency(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr;
    }
    std::uint32_t adjacencySize(std::int32_t u) const { return mRun[u].adjacencySize; }

    const std::int32_t* incidence(std::int32_t u) const {
        return mSource.data() + mRun[u].sourcePtr + mRun[u].adjacencySize;
    }
    std::uint32_t incidenceSize(std::int32_t u) const { return mRun[u].incidenceSize; }

    // How many original vertices u stands for. One until mass elimination merges into it, so a
    // degree that counts vertices has to count this rather than entries.
    //
    // Held in a flat array, and this is the condition under which the prototypes said it would
    // have to be: its members are a chain rather than a list, so the size is no longer free to
    // read. Caching it purely for locality was measured earlier and made no difference at all,
    // 2.77x against 2.76x, so the array is here for the chain and not for the cache.
    // The MAGNITUDE, for the driver, which does not want to know about the sign. The walks inside
    // this class read mWeight directly and test it, which is the whole point.
    std::uint32_t weight(std::int32_t u) const {
        const std::int32_t w = mWeight[u];
        return static_cast<std::uint32_t>(w < 0 ? -w : w);
    }


    // reach(u), as above. Not const: the mark array and its tag are scratch, and threading them
    // through every call site is what the prototypes do only because their display functions
    // needed to borrow them.
    void reachableSet(std::int32_t u, std::vector<std::int32_t>& reached);
    std::vector<std::int32_t> reachableSet(std::int32_t u);

    // |reach(u)|, counted without materializing it. A maintained degree wants the size and never
    // the set, so the same two passes run with a counter where reachableSet has a push_back, and
    // an exact refresh stops allocating once per refreshed vertex. That is the dominant
    // allocation in a branch that keeps its degrees exact, and none at all in one that bounds
    // them, which is most of why the two differ in speed by more than they differ in work.
    std::uint32_t reachableSize(std::int32_t u);

    // The same set, weighted: the number of ORIGINAL vertices reach(u) stands for. Once a branch
    // merges into live vertices, a reached vertex can stand for several, and a degree that counts
    // vertices has to count those rather than entries.
    std::uint32_t reachableWeight(std::int32_t u);

    // Eliminate the pivot: turn it into a clique, absorb the cliques it belonged to, prune the
    // edges the new clique now implies, and merge in whatever it makes indistinguishable.
    // Returns the vertices merged into the pivot's supervariable, which the caller needs in
    // order to take them out of its own bookkeeping.
    //
    // **The returned reference is valid until the next eliminate**, being a scratch buffer whose
    // capacity survives from pivot to pivot. Returning by value cost one allocation per
    // elimination that merged anything, about 1700 per ordering at 140x140. A caller that needs
    // the list to outlive the next call copies it.
    //
    // This is now the ONLY scratch here. The reachable set used to have one too and no longer
    // does: the walk writes it straight into the clique arena, since C[pivot] is the reach and
    // there was never a reason for it to exist anywhere else first. See beginElimination.
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
    // The new clique is C[p] and gets no name of its own, which is what makes the first line
    // read as what an elimination is: the pivot stops being a vertex with a reachable set and
    // becomes a clique holding that same set. Not one of the three differences builds a set;
    // each is a stamp of the subtrahend and one compaction pass over the minuend.
    const std::vector<std::int32_t>& eliminate(std::int32_t pivot);




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
    // it is what genmmd's prepass does with the degree-0 and degree-1 vertices. The lists are left
    // alone; every walk skips a numbered vertex from here on.
    void number(std::int32_t u);

    // Walk I[u] from the back in reachableSet(), matching genmmd's clique stack. A tie-break
    // convention and nothing else: it changes which permutation comes out, never which sets are
    // computed. See the member's note.
    void setReverseIncidence(bool on) { mReverseIncidence = on; }

    // Lay the lists out the way AMD_2 does, which is the opposite of genmmd's on both counts.
    // Two conventions under one switch because they are one fact, that AMD's lists run the other
    // way round, and are only ever wanted together. Like the flag above, this changes which
    // permutation comes out and never which sets are computed. Used by Amd3 alone.
    //
    //   reachableSet   walks the CLIQUES before the explicit adjacency, since AMD_2's
    //                  `for (knt1 = 1; knt1 <= elenme + 1; knt1++)` takes the cliques of me and
    //                  the supervariables only on its last pass. genmmd is the other way round,
    //                  which is why the md ladder is laid out as it is.
    //   the prune      puts the new clique at the FRONT of I[u] rather than appending, with the
    //                  displaced entries ROTATED rather than shifted, which is AMD_2's
    //                  `Iw[pn] = Iw[p3]; Iw[p3] = Iw[p1]; Iw[p1] = me`. Our two lists share one
    //                  run exactly as its do, so the same three moves apply unchanged.
    //
    // See experiments/ordering/AMD3.md, ledger entries 2 and 5.
    void setVendoredListOrder(bool on) { mVendoredListOrder = on; }

    // Stop eliminate() at the prune, leaving mass elimination to the caller. AMD_2 makes the same
    // test in its scan 2, AFTER aggressive absorption has dropped every clique lying inside the
    // new one, and says why in its own comment: with aggressive absorption, `deg == 0` is
    // identical to the structural test. Asking first, which is what eliminate() does by default,
    // asks it of an I[u] that still holds cliques about to be removed, so the cheap test declines
    // vertices AMD merges.
    //
    // With this on, eliminate() returns an EMPTY merged list and C[pivot] is reach(pivot) exactly,
    // and the caller must call massEliminate() once it has absorbed. Used by Amd3 alone. See
    // experiments/ordering/AMD3.md, ledger entry 3.
    void setLateMassElimination(bool on) { mLateMassElimination = on; }

    // The half eliminate() no longer does under the flag above: fold into the pivot's supervariable
    // every member of C[pivot] that the new clique now accounts for entirely, and trim C[pivot] of
    // them. Returns the merged vertices, from a member scratch as eliminate() does.
    //
    //     merged   = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
    //     C[pivot] = C[pivot] - merged
    //
    // Calling it without the flag is a caller error and not guarded: eliminate() will already have
    // merged, so this would find nothing and cost a pass.
    const std::vector<std::int32_t>& massEliminate(std::int32_t pivot);

    // Kill these cliques and take them out of the incidence lists of these vertices. A clique
    // dies when it is found to lie wholly inside a newer one, which is aggressive absorption, and
    // the vertices to purge are the newer clique's members, since those are the only lists that
    // can still name it.
    // `vertexCount` is `|C[p]|`, so one dimensional; every caller passes `cliqueSize(pivot)`.
    void absorb(const std::vector<std::int32_t>& cliques,
                const std::int32_t* vertices, std::uint32_t vertexCount);

    // Expand a pivot sequence into an elimination order over the original vertices. A pivot
    // stands for its whole supervariable, whose members are eliminated consecutively, so this is
    // where a supervariable of size w becomes w columns.
    std::vector<std::int32_t> order(const std::vector<std::int32_t>& pivots) const;

    // The same permutation with each supervariable's members in ASCENDING VERTEX INDEX rather
    // than merge order. Indistinguishable members, so the fill and the forest are unchanged;
    // only the permutation is, and it is genmmd's. Used by Mmd3 alone. See the .cpp.
    std::vector<std::int32_t> orderAscending(const std::vector<std::int32_t>& pivots) const;

private:
    // The head and tail of an elimination, shared by the two overloads above and private because
    // between them the graph is half eliminated: the clique is written and stamped but the reached
    // vertices still name the pivot as a variable. Nothing outside may observe that state, which is
    // why the seam is two private calls rather than a public begin and end.
    void beginElimination(std::int32_t pivot, std::int32_t& absorbed);

    const std::vector<std::int32_t>& finishElimination(std::int32_t pivot);

    // A[u] and I[u] share one run, and one array holds every run end to end. The two lists are
    // the SOURCES of reach(u), one per explicit neighbor and one per clique, and their number is
    // conserved: each elimination that reaches u replaces at least one source with the new
    // clique, the explicit neighbor `pivot` where u was reached through A[pivot], or an absorbed
    // clique where it was reached through one. So
    //
    //     |A[u]| + |I[u]| <= the number of off-diagonal entries in u's column of A
    //
    // holds for the whole run, and u's block can be sized once from the pattern and never grown.
    // Section 5.3 of archive/sparse_factorization.md carries the argument; 5.15 records that both
    // vendored routines rely on it and that neither states it.
    //
    // The order within the run is forced rather than chosen. The prune compacts A[u] and then
    // I[u], and A[u] shrinks by at least what I[u] gains, so writing the incidence behind the
    // adjacency always trails the read cursor. The other order would have the incidence overwrite
    // an adjacency entry it had not read yet.
    //
    // C[c] gets none of this and needs its own arena below: a clique's members are its pivot's
    // reach, which nothing about the pivot's own column bounds.
    // ONE DIMENSIONAL SIZES ARE `std::uint32_t`, POSITIONS ARE `std::size_t`. `sourcePtr` offsets
    // into the pool and is bounded by nnz(A), so it is two dimensional and stays wide; the two
    // lengths are bounded by deg(u) and so by n, which the constructor caps at `MAX_IDX`. The two
    // kinds meet at exactly one place, the constructor's crossing, where a difference of positions
    // is written as a count, and that is where the cast belongs. `incidenceSize` has no such
    // crossing at all: it is only ever written from a cursor or from another length.
    //
    // The CONSERVATION LEMMA is why one uint32 covers both. Their sum is bounded by u's column of
    // A for the whole run, so neither can reach n on its own where the other is nonzero, and
    // `adjacencySize(u) + incidenceSize(u)` cannot overflow either.
    std::vector<std::int32_t>  mSource;         // every A[u] then I[u], run after run

    // ONE OBJECT PER VERTEX, NOT THREE ARRAYS, and the shared QuotientGraph carries the same
    // struct with the same reasoning; read its member for the argument and the measurement. The
    // one difference is the LENGTH. This scheme addresses a clique by the pivot that formed it
    // and needs the END of that pivot's segment, `mRun[c + 1].sourcePtr`, so the array has n + 1
    // entries where the shared class has n. The extra entry's two lengths are never read.
    struct VertexRun {
        std::size_t   sourcePtr;       // where u's run starts in mSource, fixed at construction
        std::uint32_t adjacencySize;   // A[u]'s length, from the run's start
        std::uint32_t incidenceSize;   // I[u]'s length, immediately behind A[u]
    };
    std::vector<VertexRun> mRun;

    // C[c] is flat too, and for a second property rather than the adjacency's. Its members are
    // not known in advance, so its block cannot be placed at construction; but they are known
    // exactly when the clique is formed, at the moment its pivot is eliminated, and the set only
    // ever shrinks afterwards. So a bump allocator suffices: take the next block, write the
    // members, and that block is the clique's for as long as it lives.
    //
    // Nothing is reclaimed. An absorbed clique leaves a hole, which is what the vendored pool
    // answers with a compaction pass, and we do not need the answer: the arena's total is the sum
    // of every clique ever formed, which is bounded by nnz(L) and is a few megabytes at the sizes
    // we run. The offsets being indices rather than pointers is what lets the arena simply grow
    // instead, since a reallocation leaves every offset valid.
    // THE CHANGE UNDER TEST. There is no clique arena. C[c] lives in the segment of the vertex
    // that formed it, mSource[mRun[c].sourcePtr ...], which is dead the moment c is eliminated, so a
    // clique costs no storage and is addressed by its pivot. When the reach outgrows that segment
    // the list CONTINUES in the segment of a clique the pivot has just absorbed, and the last
    // entry of the segment being filled holds a LINK to it.
    //
    // A link is -(c + 1) and not genmmd's -c: ids here start at 0, so -c cannot distinguish
    // clique 0 from member 0. The encoding maps [0, 2^31-1] onto [-2^31, -1] exactly, no value
    // spare and none overlapping, which is what makes a sign test the whole discriminator.
    //
    // genmmd ends a list with a 0 entry, which works only because its ids start at 1. Ours is a
    // reserved value, TERMINATOR, for the same job. THIS FILE KEEPS NO CLIQUE LENGTH: it did, in
    // `mCliqueSize`, and the stamp scheme retired it; see the note on mMark below. The paragraph
    // that stood here said the opposite and had been stale since.
    std::vector<std::int32_t>  mAbsorbed;     // I[pivot], copied out before its segment is overwritten

    // The supervariable a vertex stands for, as a chain rather than a list per vertex. A list
    // meant one allocation per vertex before anything had happened, n of them for a structure
    // that is usually a singleton, plus a growth every time one absorbed another. The chain is
    // three arrays allocated once: the next member, the last one (so an absorption appends in
    // O(1) and the members keep their order, which the emitted permutation depends on), and the
    // count, which a chain no longer gives away for free.
    // `mWeight` is one dimensional and so `std::uint32_t`, and the bound is not the term count but
    // DISJOINTNESS: the weights partition the original vertices, so a sum of them over a set of
    // distinct vertices is at most n however many terms it has. That covers every accumulation of
    // weights in the ordering, in `reachableWeight`, in the clique weight below, and in the
    // drivers' bound and refresh passes. An accumulation over an unbounded number of one
    // dimensional terms would otherwise need a wider type, and this is the exception to that.
    std::vector<std::int32_t>  mSuperNext;
    std::vector<std::int32_t>  mSuperLast;
    // SIGNED, MIRRORING QuotientGraph, 2026-08-17. A one dimensional size is normally unsigned
    // because it has nothing to stand in for; this one has. `AMD_2`'s `Nv`: positive is the
    // weight, negative means already taken into the clique being built, zero means dead, so one
    // load answers what two arrays answered. No range is lost, a weight being bounded by n.
    //
    // IT IS HERE BECAUSE THIS FILE'S OBLIGATION IS TO STAY ENCODING-IDENTICAL TO Mmd3, not because
    // it pays on its own: it measured a wash on the mmd side. Without it the time column stops
    // being the price of genmmd's storage and becomes that plus four encoding folds, which is
    // exactly what this file exists not to be. See docs/CODING_RULES.md for when a size may go
    // signed and docs/DESIGN_DECISIONS.md (2026-08-17).
    std::vector<std::int32_t> mWeight;
    // Mirrors QuotientGraph::mHasNumbered; always true here once the prepass has run.
    bool                      mHasNumbered = false;
    std::uint32_t              mCliqueWeight = 0;  // see cliqueWeight(); per-pivot, not per-vertex

    std::vector<std::int32_t> mMerged;   // scratch for the vertices an elimination merges away

    // Which end of I[u] reachableSet() walks from. genmmd threads its clique list through an
    // integer array and pushes at the head, `list[nb] = el; el = nb`, then reads from the head, so
    // the clique seen LAST is expanded FIRST; we hold a vector and append. Same set either way and
    // the same cost, but the order decides C[pivot]'s order, hence which of two equal-degree
    // candidates a later iteration finds first, and minimum degree is settled by exactly that.
    // Off by default, so every existing driver is unaffected; Mmd3 turns it on. See
    // experiments/ordering/mmd3.py, where the same four walks are reversed together.
    bool mReverseIncidence = false;

    // AMD_2's list conventions, off by default so the other five drivers are untouched. Read in
    // reachableSet() and in the prune, hoisted at both sites for the same reason mReverseIncidence
    // is: a member load the compiler cannot prove is unaliased by the stores in the loop. See the
    // setter for what each half does and why they are one flag.
    bool mVendoredListOrder = false;

    // Whether eliminate() stops at the prune and leaves mass elimination to the caller. Off by
    // default. See the setter, and massEliminate() for the half it hands over.
    bool mLateMassElimination = false;

    // TWO STAMP SPACES IN ONE ARRAY: vertices at [v], cliques at [size() + c]. Clique ids ARE
    // vertex ids, so one space cannot hold both once the vertex half carries GONE: stamping a
    // clique would write a live tag over the slot of the dead pivot that formed it, and a walk of
    // an older clique still holding that pivot as a member would take it for live. Neither genmmd
    // nor AMD_2 shares one array between the two kinds, and this is the cheap way to stop doing
    // so. It is what lets the dead-clique test be a stamp again rather than a size, which is what
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
    : mRun(colPtr.empty() ? 1 : colPtr.size()),
      mCliqueLiveMembers(mRun.size() - 1, 0),
      mMark(2 * (mRun.size() - 1), NIL) {
    const std::int32_t size = static_cast<std::int32_t>((mRun.size() - 1));

    // One pass to place each column's run, dropping the diagonal. The runs are laid out in column
    // order and never move afterwards, so the offsets are written once here and only the lengths
    // change. The run is u's whole allowance: A[u] fills it now and I[u] grows into the room A[u]
    // gives up, never past mRun[u + 1].sourcePtr.
    mSource.reserve(colPtr.empty() ? 0 : colPtr.back());
    mRun[0].sourcePtr = 0;
    for (std::int32_t aj = 0; aj < size; ++aj) {
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) mSource.push_back(rowIdx[cp]);
        mRun[aj + 1].sourcePtr = mSource.size();
        // THE ONE CROSSING. A difference of two positions is a count, so this is the single place
        // in the class where a two-dimensional quantity is written into a one-dimensional one. It
        // is bounded by deg(aj) and so by n, which the SparseMatrix constructor has already capped
        // at MAX_IDX, but the cast is written rather than left implicit because that bound is an
        // argument and not something the types say.
        mRun[aj].adjacencySize = static_cast<std::uint32_t>(mRun[aj + 1].sourcePtr
                                                            - mRun[aj].sourcePtr);
    }

    // There is no second arena to reserve: cliques live in mSource, in the segments their
    // pivots vacate. That is the change this file exists to measure.

    // Every vertex begins as a supervariable of one: a chain holding itself, so next is NIL and
    // last is the vertex. Mass elimination splices these together and order() walks them.
    mSuperNext.assign(size, NIL);
    mSuperLast.resize(size);
    mWeight.assign(size, 1);
    for (std::int32_t u = 0; u < size; ++u) mSuperLast[u] = u;
}

inline std::uint32_t QuotientGraphChained::reachableWeight(std::int32_t u) {
    ++mTag;
    std::uint32_t reached = 0;   // a sum over DISTINCT vertices, so bounded by n; see the header
    mMark[u] = mTag;
    // The bounds are hoisted, all of them. Each is a load from a member vector, and the bodies
    // below store through mMark, which the compiler cannot prove does not alias the sizes, so a
    // bound left in the condition is re-loaded once per clique. Measured at 300 ms of AMD1's
    // 6.31 s on alpamayo, in the one accessor that reads as though it were free.
    const std::int32_t* source        = mSource.data() + mRun[u].sourcePtr;
    const std::uint32_t adjacencySize = mRun[u].adjacencySize;
    const std::uint32_t incidenceSize = mRun[u].incidenceSize;
    for (std::uint32_t k = 0; k < adjacencySize; ++k) {
        const std::int32_t v = source[k];
        if (mMark[v] != GONE) { mMark[v] = mTag; reached += static_cast<std::uint32_t>(mWeight[v]); }
    }
    const std::int32_t* incidence = source + adjacencySize;
    for (std::uint32_t i = 0; i < incidenceSize; ++i) {
        const std::int32_t c = incidence[i];
        forEachMember(c, [&](std::int32_t v) {
            if (mMark[v] < mTag) { mMark[v] = mTag; reached += static_cast<std::uint32_t>(mWeight[v]); }
        });
    }
    return reached;
}

inline void QuotientGraphChained::number(std::int32_t u) {
    mHasNumbered = true;     // see the member: what makes the GONE test worth asking
    mMark[u]     = GONE;     // a numbered vertex lingers in lists; GONE is what filters it
}

inline void QuotientGraphChained::beginElimination(std::int32_t pivot, std::int32_t& absorbed) {
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
    // The remedy is to make the arena unable to move rather than to re-fetch per clique, which
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
    // C[pivot] IS BUILT IN THE PIVOT'S OWN SEGMENT, which is dead from this moment: A[pivot] and
    // I[pivot] are read here and never again. When the reach outgrows it the list continues in the
    // segment of a clique the pivot absorbs, linked by a negative entry. This is mmdelm, ported.
    //
    // THE ONE INVARIANT THAT MAKES IT SAFE, and it is worth stating because it is not obvious: the
    // link written at the current segment's last entry always names the clique ABOUT TO BE WALKED.
    // So the write cursor can only ever enter a segment we are reading, never one still to be
    // read, and within that segment writes trail reads because every entry written was read first.
    const bool reverse = mReverseIncidence;

    const std::size_t   base = mRun[pivot].sourcePtr;
    const std::uint32_t adjN = mRun[pivot].adjacencySize;
    const std::uint32_t incN = mRun[pivot].incidenceSize;

    // I[pivot] must be copied out before the segment holding it is overwritten. genmmd threads the
    // same list through its `list[]` array; a small vector kept for its capacity costs no
    // allocation after the first few pivots and keeps the walk order explicit.
    mAbsorbed.clear();
    for (std::uint32_t ii = 0; ii < incN; ++ii)
        mAbsorbed.push_back(mSource[base + adjN + (reverse ? incN - 1 - ii : ii)]);

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
    std::size_t   rm    = mRun[pivot + 1].sourcePtr - 1;   // last entry of the segment being filled

    // One member written, following a link first if the segment is full. The link can only be
    // there to be followed: it is written before the walk that can reach it.
    const auto emit = [&](std::int32_t v) {
        while (rl >= rm) {
            const std::int32_t l = mSource[rm];         // -(c + 1); see the link encoding above
            rl = mRun[-l - 1].sourcePtr;
            rm = mRun[-l].sourcePtr - 1;
        }
        mSource[rl++] = v;
    };
    // THE LIVE NEIGHBORS FIRST, AND WITHOUT THE BOUND CHECK. They are a subset of A[pivot], which
    // was read from this same segment, so the cursor cannot pass the reader and cannot leave the
    // segment: at most adjN entries are written into a segment holding adjN + incN. genmmd's
    // first loop is unchecked for exactly this reason, and checking here is not merely wasteful,
    // it is wrong -- with incN == 0 the cursor legitimately lands on the last entry, and a check
    // would read it as a link that was never written. Found by ASan on a 2 by 2 grid.
    std::uint32_t born = 0;                     // the clique's size; see numPeakCliqueMembers
    for (std::uint32_t k = 0; k < adjN; ++k) {
        const std::int32_t v  = mSource[base + k];
        const std::int32_t nv = mWeight[v];
        if (nv > 0 && !(mHasNumbered && mMark[v] == GONE)) {
            mWeight[v] = -nv;
            mSource[rl++] = v;
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

        mSource[rm] = -(c + 1);                     // the continuation, before it can be read
        forEachMember(c, [&](std::int32_t v) {
            const std::int32_t nv = mWeight[v];        // one load; see the note above the walk
            if (nv > 0) { mWeight[v] = -nv; emit(v); ++born; }
        });
    }

    // THE TERMINATOR, where genmmd writes `if (rl <= rm) adjncy[rl] = 0`. The cursor stops one
    // short of the segment end whenever a link was needed there, so there is room; the one case
    // with no room is a clique whose members fill its last segment exactly, and the walk's second
    // stop condition covers that.
    // Against `rm`, the CURRENT segment's last entry, not the pivot's: the cursor has followed
    // every link the emit needed and is wherever that left it. genmmd writes `if (rl <= rm)` for
    // exactly this reason. Comparing against the pivot's own segment end instead leaves a chain
    // unterminated and the walk runs off into whatever the next segment holds, which on a 2x2
    // grid is an immediate hang.
    if (rl <= rm) mSource[rl] = TERMINATOR;

    // The absorbed cliques die only now: they were being READ until the loop above finished, and
    // their segments now hold part of the new clique, so nothing can be written into them to say
    // so. The stamp below is what says it.

    // NO STAMPING PASS, AND NO `inClique`. Membership was written by the walk, in the sign of the
    // weight, so there is nothing to stamp and nothing for a tag to say; the out-parameter went on
    // 2026-08-17 when the prune started reading the sign. genmmd has no such pass either. Only
    // `absorbed` survives, and that is this file's own scheme rather than a shared one: it marks
    // the CLIQUE half of mMark, which the shared class no longer has at all.
    //
    // What kept the pass alive was the weighted clique size accumulated beside it. That value has
    // NO READER in this file: `cliqueWeight()` exists for the amd drivers on the shared class, and
    // this one carries the mmd driver alone, whose refresh computes `dg0` itself from the clique
    // members. So the accumulation goes with the walk rather than moving into the emit.
    //
    // The pivot reads negative from the negation above, so it reads as a member of its own
    // clique, and the prune's `mWeight[v] <= 0` drops it with the rest.
    // AND THE NEW CLIQUE IS BORN, the maximum taken here alone since nothing else raises the
    // total. See numPeakCliqueMembers.
    mCliqueLiveMembers[pivot] = born;
    mNumLiveCliqueMembers    += born;
    mNumPeakCliqueMembers     = std::max(mNumPeakCliqueMembers, mNumLiveCliqueMembers);

    ++mTag;
    absorbed = mTag;
    const std::size_t cliqueBase = mRun.size() - 1;
    for (const std::int32_t c : mAbsorbed) mMark[cliqueBase + static_cast<std::size_t>(c)] = absorbed;
}

inline const std::vector<std::int32_t>& QuotientGraphChained::eliminate(std::int32_t pivot) {
    std::int32_t absorbed = NIL;
    beginElimination(pivot, absorbed);
    const std::size_t cliqueBase = mRun.size() - 1;
    // C[pivot] IS the reach here: beginElimination wrote it there and nothing has trimmed it yet
    // (massEliminate does, and runs after). So the prune walks the arena block rather than a copy.
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
    forEachMember(pivot, [&](std::int32_t u) {
        std::int32_t*     source        = mSource.data() + mRun[u].sourcePtr;
        // The two counters are one-dimensional COUNTS, positions in a list bounded by deg(u) and
        // so by n, where std::size_t is for a position into an n x n object. They take the type of
        // what they count, so they move with the array: `std::uint32_t`, and the cast that used to
        // stand here is gone. It was there only because the array was wider than the loop. The
        // dimensional rule in experiments/ordering/REPORT.md asks for this everywhere and it was
        // taken here first because this is the hottest loop in the ordering, Instruments put
        // 277 ms of an 8.53 s run on the incidence loop's header.
        //
        // `heldVertex` is a VERTEX and stays std::int32_t, carrying NIL.
        const std::uint32_t adjacencySize = mRun[u].adjacencySize;
        std::uint32_t       kept          = 0;
        std::int32_t        heldVertex    = NIL;        // the first survivor, appended last
        for (std::uint32_t k = 0; k < adjacencySize; ++k) {
            // ONE LOAD, THREE QUESTIONS. A negative weight is a member of the new clique, the
            // pivot included, so the explicit pivot test goes with the membership test; a zero is
            // a vertex a live merge folded away. The fourth question, whether v was numbered by
            // the prepass, still needs mMark.
            const std::int32_t v = source[k];
            if (mWeight[v] <= 0) continue;             // in the new clique, the pivot, or merged
            if (mHasNumbered && mMark[v] == GONE) continue;   // numbered by a prepass
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
        const std::uint32_t incidenceSize = mRun[u].incidenceSize;
        std::uint32_t       write         = kept;      // follows `kept`; see the note above
        for (std::uint32_t i = 0; i < incidenceSize; ++i) {
            const std::int32_t c = incidence[i];
            if (mMark[cliqueBase + static_cast<std::size_t>(c)] != absorbed)
                source[write++] = c;
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
    });

    return finishElimination(pivot);
}

// The same prune, with the driver's first scan folded into the two loops. The header carries the
// argument for why this is one call and why it cannot be folded further; the loops below are the
// plain ones with an accumulation added on each survivor, and nothing else differs.
// The same prune again, with the driver's first scan folded in under `Amd.cpp`'s tagged-W
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
inline const std::vector<std::int32_t>& QuotientGraphChained::massEliminate(std::int32_t pivot) {
    std::vector<std::int32_t>& merged = mMerged;   // scratch, kept for its capacity
    merged.clear();
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
        // Under mVendoredListOrder the new clique sits at the FRONT of I[u] rather than the back,
        // so the single remaining entry is at the head of the incidence run either way: with A[u]
        // empty the run starts with I[u], and with one clique there is only one position. The
        // test therefore needs no branch on the flag.
        if (mRun[u].adjacencySize == 0 && mRun[u].incidenceSize == 1 &&
            mSource[mRun[u].sourcePtr] == pivot) {         // A[u] empty, so I[u] starts at the run
            mRun[u].incidenceSize = 0;
            mMark[u]          = GONE;
            merged.push_back(u);
        }
    });
    // NO COMPACTION OF C[pivot], AND THAT IS GENMMD'S BEHAVIOUR RATHER THAN AN OMISSION,
    // 2026-08-18.
    // `mmdelm`'s n1100 loop mass-eliminates with `qsize[md] += qsize[rn] ; qsize[rn] = 0` and DOES
    // NOT rewrite the list: the merged vertex keeps its place and every later reader skips it on
    // `qsize[nb] != 0`. The same is true of the further absorption inside `mmdupd`. So genmmd
    // places a clique once and never shortens it, where `AMD_2` compacts in its restore pass and
    // hands the tail back with `pfree = p`.
    //
    // WE USED TO COMPACT HERE, following the chain on two cursors, which cost a pass genmmd never
    // pays and saved the skipped entries it pays on every later read. Neither is visible in a
    // permutation, a skipped member and an absent one reading alike, so no check in this tree
    // would have caught it; it was found by reading `mmdelm`. Since this file exists so that a
    // differential against genmmd measures the layout and nothing else, the pass had to go.
    //
    // EVERY READER ALREADY SKIPS THE DEAD, which is what makes the removal a deletion:
    // `reachableSet` tests `nv > 0`, genmmd's own test; the reach count tests the mark, and GONE
    // outranks any tag; the driver's clique walk tests `eliminated`; its pair test rejects on
    // `m >= vertexTag`. The one loop with no test is the eviction over C[pivot], which is
    // `mmdelm`'s `bwd[rn] = 0` and which genmmd likewise runs over the uncompacted list;
    // the unfile-and-restore pair is idempotent, so a merged vertex is evicted harmlessly.
    //
    // THE COST IS SPACE, and it is the space genmmd spends. A clique keeps its dead members for as
    // long as it lives, so live clique storage here is strictly above what the compacting classes
    // report for the same ordering.
    // THE MERGED LEAVE THE LIVE COUNT even though they do NOT leave this file's storage: `mmdelm`
    // keeps them in the list and skips them on `qsize != 0`. Tracking the notional size is what
    // makes the figure comparable with `Mmd3`, which does drop them. See numPeakCliqueMembers.
    mCliqueLiveMembers[pivot] -= static_cast<std::uint32_t>(merged.size());
    mNumLiveCliqueMembers     -= merged.size();

    if (!merged.empty()) {
        for (std::int32_t u : merged) {                // the pivot now stands for them too
            mSuperNext[mSuperLast[pivot]] = u;         // append u's chain, order preserved
            mSuperLast[pivot]             = mSuperLast[u];
            // The weighted clique size follows the clique. `cliqueWeight()` promises the weighted
            // size of the LIVE members of C[pivot], and a merged vertex has just stopped being
            // one, so the decrement belongs here rather than in the caller. It is unaffected by
            // the list no longer being compacted: the entry stays, the weight does not count.
            // AMD_2 spells the same line `degme -= nvi` inside its own mass elimination.
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
    }
    return merged;
}

inline void QuotientGraphChained::merge(std::int32_t u, std::int32_t v) {
    mSuperNext[mSuperLast[u]] = v;                 // append v's chain, order preserved
    mSuperLast[u]             = mSuperLast[v];
    mWeight[u] += mWeight[v];
    mWeight[v] = 0;

    mRun[v].adjacencySize = 0;
    mRun[v].incidenceSize = 0;
    mMark[v]          = GONE;
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
inline std::vector<std::int32_t> QuotientGraphChained::orderAscending(
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

} // namespace Oblio
