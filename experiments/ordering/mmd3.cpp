// mmd3.cpp -- multiple minimum degree, matching genmmd's permutation.
//
// mmd2 has every mechanism genmmd has, and still returns a different permutation.
// This layer adds no mechanism at all. It changes the ORDER FOUR LISTS ARE WALKED IN,
// and that is the whole difference between the two files.
//
// THE ALIGNMENT LEDGER. One row per divergence found, root-caused and closed. It is kept
// here rather than in a note because the code is where it will be read, and it is APPEND
// ONLY: a row is never edited once closed, so the sequence stays a record of what was
// actually wrong rather than a tidy summary written afterwards. The same table is in
// mmd3.py and is summarized in the README's mmd3 section.
//
//   #  what diverged                where here            genmmd                nature
//   -  ---------------------------  --------------------  --------------------  ----------
//   1  clique expansion            mmd3Neighbors, I[u]   mmdelm, the el stack  convention
//   2  two-source walk               the refresh           mmdupd, q2h           convention
//   3  many-source walk              the refresh           mmdupd, qxh           convention
//   4  batch clique order          the driver            genmmd, ehead         convention
//   5  merged weight in a           the refresh, two-source  mmdupd, dg -          DEFECT
//      supervariable's bucket                             qsize[en] + 1
//   6  supervariable member order   the final expansion   mmdnum, the scan      cosmetic
//
// All six closed 2026-08-07. The nature column is the one to read first.
//
// DEFECT means wrong on its own terms, independent of genmmd: entry 5 filed a supervariable
// one bucket too high per vertex merged into it, so it was never picked as early as its size
// had earned, and the code did not do what its own comment said. It was FIXED IN mmd2 AND IN
// PRODUCTION Mmd2 as well, where it had been costing fill since those were written; mmd2's
// gap against genmmd fell from about 20 percent to about 7. mmd1 cannot have it, having no
// two-source path and no live merges.
//
// THE AMD LAYERS COULD, and this paragraph used to say they could not. It argued that AMD
// files at an external degree, which excludes a vertex's own supervariable and so does not
// move when its weight changes. The external degree does not move; the `- weight(u)` term
// inside the amd bound does, and it is the term that decides the bucket. amd2, Amd2 and
// Amd2B carried exactly this defect and it was costing 3 to 9 percent of fill on grids,
// found on 2026-08-08 by aligning amd3 against the vendored AMD, where it is that ledger's
// entry 4, and fixed the same day. amd1 and Amd1B genuinely cannot have it, having no live
// merges at all.
//
// Worth keeping as a lesson about ledgers rather than only as a correction. The claim was
// reasoned from a true premise and never checked against the quantity actually filed. A
// defect ruled out by an argument is ruled out only as far as the argument reaches, and the
// same sentence stood in the README for a month.
//
// CONVENTION means neither side is wrong. A tie-break has no right answer, and 1 to 4 only
// change which of several equal-degree vertices wins. What they cost is not quality but
// COMPARABILITY: while they differed, no measurement against genmmd could separate a
// difference of mechanism from a difference of arbitrary choice.
//
// COSMETIC means it cannot change the answer at all. Entry 6 reorders the members of a
// supervariable, which are indistinguishable by construction, so the fill and the elimination
// forest are identical either way.
//
// 5 is a different kind, and the first real defect rather than a convention. In the two-source
// refresh we subtracted u's own weight from cliqueWeight BEFORE the walk; genmmd keeps cliqueWeight
// whole and subtracts at the end. Those differ exactly when the walk MERGES a vertex into u, since
// genmmd's merge does `qsize[en] += qsize[nd]` in the same walk and the weight it then
// subtracts is the post-merge one. Subtracting first files a supervariable one bucket too
// high per vertex merged into it, so it is not picked as early as its size has earned.
// Closing it took 32x32 from 1.5 percent fill above genmmd to ZERO, and matched 5x5 and 7x7
// outright.
//
// 1 to 4 are one defect found four times, described under THE FOUR LISTS below. They were
// closed together because none can be judged alone: with 2, 3 and 4 the seven examples
// reach four of seven, and only 1 takes them to seven of seven, since 1 fixes the order of
// C[pivot] and therefore the content order of what 2, 3 and 4 walk.
//
// 6 cannot change the fill and does not pretend to. A supervariable's members are
// indistinguishable by construction, so any order among them gives the same factor, which
// is why fill reached genmmd's exactly while the printed permutation still differed. It is
// closed anyway, because an exact permutation makes the comparison an equality test rather
// than a judgement, and that is the instrument the next layer will be aligned with.
//
// WHERE IT STANDS. **mmd3 returns genmmd's permutation EXACTLY**, on all seven
// examples and on every square grid tested from 5 a side to 80, n = 6400, and its
// fill is genmmd's to the digit at every size on the scale ladder where mmd2 ran
// 12 to 25 percent above. Six alignments got it there, four conventions, one
// cosmetic numbering and one real defect, with no mechanism added to mmd2 at all.
//
// (This block used to record the position after entries 1 to 4, seven of seven
// examples but 67 to 78 percent of the pivots on grids, and concluded that one
// more mechanism was unaccounted for. Entry 5 was that mechanism. The paragraph
// was not revisited when it closed, which is worth one line of warning about
// headers written mid-work.)
//
// WHY THAT MATTERS AND IS NOT A DETAIL. Minimum degree is a tie-break algorithm. At
// almost every iteration several vertices share the least degree, and which one is taken
// is decided by whatever the data structure hands over first. Two codes can agree on
// every rule and still part company on the first tie, and from there they are ordering
// different graphs. Matching the permutation is the only way to tell a difference of
// MECHANISM from a difference of ARBITRARY CHOICE; until it holds, a fill comparison
// measures both at once.
//
// THE FOUR LISTS, all the same idiom. genmmd threads a linked list through an integer
// array, pushes at the head, `list[nb] = h; h = nb`, then reads from the head, so the
// entry seen LAST is processed FIRST. We hold a vector and append: same set, opposite
// order, same cost.
//
//   mmd3Neighbors, the I[u] walk   mmdelm's clique stack, list[nb] = el; el = nb
//   the two-source walk              mmdupd's, list[nb] = q2h; q2h = nb
//   the many-source walk             mmdupd's, list[nb] = qxh; qxh = nb
//   the batch walk                 the driver's, list[mn] = ehead; ehead = mn
//
// The first is the deepest and the others cannot be judged without it: it fixes the order
// of C[pivot], the content order of every list built downstream. The other three alone
// reach four of seven examples; all four reach seven of seven.
//
// WHY THIS IS A LAYER RATHER THAN AN EDIT TO mmd2. mmd2 is a different ordering and
// a correct one: it has every mechanism genmmd has, and only its arbitrary choices
// differ. Folding these six alignments into it would delete that ordering rather
// than add one, and would leave nothing to measure the conventions against. So the
// two sit side by side, mmd2 carrying our tie-break and mmd3 carrying genmmd's, and
// `Ordering::MMD2` remains selectable.
//
// (This block used to say the opposite reason: that grids still parted company in
// the last third of the run, so one mechanism was unaccounted for and this file was
// the anchor it would be found against. That was true when it was written and entry
// 5 closed it the same day. Stale until 2026-08-08.)
//
// NOT A SORT. genmmd never sorts and neither does this. Sorting the clique members
// reaches five of seven and reproduces the effect without the mechanism, at the cost of a
// sort per clique in a routine designed around not having one.
//
// Everything below this header is mmd2, unchanged.

// PASS 1, THE PREPASS. genmmd numbers every vertex in the degree-1 list before the
// main loop starts, marks each marker[mn] = maxint, and never refreshes a
// neighbor. Two things travel with it. mmdint files a degree-0 vertex under degree
// 1, `if(dg==0)dg=1`, so isolated and degree-1 vertices are numbered together, and
// the bucket a vertex sits in stops being its true degree. And head[1] = 0
// afterwards, with the main loop starting at mdeg = 2.
//
// A prepass vertex is NOT eliminated in the quotient-graph sense: no clique is
// formed, nothing is pruned, and its neighbors keep degrees that still count it.
// It is simply numbered and then skipped, which marker[mn] = maxint does there and
// `eliminated[u]` does here. The neighbor query skips such vertices, and the prune
// loop drops them when it compacts, which is what mmdelm's `marker[nb] < tag` test
// does.
//
// PASS 2, THE TWO-SOURCE SPLIT. mmdupd does not walk a flat list of reached vertices. It
// walks the CLIQUES created this iteration, `el = list[el]`, and for each one it
// computes cliqueWeight once, the weighted size of that clique, then visits the clique's
// members. A member is classified by what it has left BESIDES the new clique:
// mmdelm stashes fwd[rn] = nq + 1 where nq counts the survivors of the compaction,
// which here is A[u].size() + I[u].size() - 1. nq == 1 puts the vertex on mmdupd's q2h
// chain, anything else on its qxh. Ours are twoSourceQueue and manySourceQueue.
//
// The two-source case is answered without a union. Everything the vertex reaches is
// either inside the clique, already counted in cliqueWeight, or comes from that one other
// source, so the walk adds only what the other source contributes. Two mark levels
// keep that straight, as mmdupd's mt and *tag do: cliqueTag says "already in
// cliqueWeight" and survives the whole clique, while vertexTag is fresh per vertex, so one
// two-source vertex cannot hide a neighbor from the next.
//
// The degrees come out identical either way, which is the check on this pass. What
// does move is the ORDER of the filing, since the refresh is now clique by clique
// with twoSourceQueue before manySourceQueue, and filing order decides what a bucket
// holds.
//
// A vertex reached by two pivots in the same iteration is refreshed once: mmdupd skips
// it on the second visit with `if(bwd[en]!=0) goto n2200`, since a refiled vertex
// has a bucket again. Here that is the `filed` flag.
//
// PASSES 3 AND 4, THE PAIRWISE MERGE AND OUTMATCHED MARKING. Both live in one branch
// of the two-source walk, reached when a member of the one other source is ALSO a member
// of the new clique:
//
//   else if(bwd[nd]==0){
//       if(fwd[nd]==2){qsize[en]+=qsize[nd];qsize[nd]=0;marker[nd]=maxint;
//                      fwd[nd]=-en;bwd[nd]=-maxint;}
//       else if(bwd[nd]==0)bwd[nd]=-maxint;}
//
// MERGE. If nd is two-source too, its only other source is that same clique, so en and nd
// reach exactly the same vertices and are indistinguishable. en absorbs nd. This is
// the first merge in the whole sequence that folds a vertex into a LIVE one, which
// is why the weight array returns: from here a candidate can stand for several
// original vertices, and every degree has to count them rather than count entries.
// It is also what makes MMD's supervariables coarser than md3's, whose test only
// ever compares a vertex against the pivot.
//
// OUTMATCHED. If nd is not two-source it has other sources besides these two, so its reach
// contains en's. It can never be the minimum before en, and MMD withdraws it from
// the degree lists rather than refiling it: bwd[nd] = -maxint. It is not merged and
// not eliminated, just held out until something reaches it again, at which point
// mmdelm restores it with bwd[rn] = 0. Here that is the `outmatched` flag, cleared
// in mmd3Eliminate for every vertex the new clique reaches.
//
// PASS 5, THE FILING CONVENTION. mmdupd does not file a vertex under its degree.
// It files under `dg = dg - qsize[en] + 1`, floored at 1, where dg was the weighted
// reach INCLUDING en's own members. So the bucket index is the external degree plus
// one, and the floor catches the case where a vertex reaches nothing outside itself.
//
// mmdint, meanwhile, files at the plain degree, `dg = xadj[nd+1] - xadj[nd]`, with
// only the zero case lifted to 1. So MMD runs on two scales: the initial buckets
// hold degrees, every refiled bucket holds degree + 1. That is genuine, not a
// misreading, and it tilts the pivot choice slightly against refreshed vertices,
// which sit one bucket higher than an untouched vertex of the same reach.
//
// From here degrees[] holds the FILED value, which is what the picker compares and
// what minDegree tracks. The nnz(L) accounting does not use it: that sums weights
// over the live members of C[pivot] and is unaffected.
//
// PASS 6, THE COUNTERS. Two small things in genmmd's main loop.
//
// ncsub, `*ncsub += mdeg + qsize[mn] - 2`, accumulated per pivot. It is the
// statistic genmmd returns alongside the permutation, an estimate of the subscript
// storage the factor will need, and it is computed from values the loop already
// has. Not checked against the vendored number yet: mmd_order takes it as a local
// and drops it, so one temporary line in the wrapper prints it, the same move
// tools/hook_amd.py makes on the AMD side. Unchecked, not uncheckable.
//
// The early termination, `if((num+qsize[mn])>neqns)goto n1000`, checked after the
// pivot is numbered and before it is eliminated. When the last supervariable is
// reached there is nothing left to update, so genmmd skips the elimination and goes
// straight to the numbering. Ours is the same test on numEliminatedVertices.
//
// Build:  g++ -std=c++17 -O3 mmd3.cpp -o mmd3_cpp  (or: make)
// Run:    ./mmd3_cpp
//         ./mmd3_cpp 3      just the third example

#include <cstdlib>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Plain vectors, UNSORTED, and a vector indexed by clique id, not std::map. A set
// costs O(log d) per membership test and per insertion; keeping a vector sorted
// costs a merge per union. Neither is needed: membership comes from a MARK array
// stamped with a tag, so "is v in the new clique" is one comparison, and every
// pass is linear in what it touches. That is what the vendored codes and Oblio's
// own SymFactorEngine do. See the README section on complexity.
//
// Types follow Oblio's rule, and there are THREE of them. An INDEX names a vertex
// or a clique and is a std::int32_t, signed only because NIL has to share a type
// with the values it stands in for. A ONE DIMENSIONAL SIZE is bounded by n and is
// a std::uint32_t: nothing to stand in for, so no sentinel and no sign bit spent.
// A TWO DIMENSIONAL size or position is bounded by nnz and is a std::size_t. An
// entity loop therefore has a signedness cast where its int32 counter meets a
// uint32 bound.
constexpr std::int32_t NIL = -1;

// The mark array is a set and the tag names it, so a tag must never repeat: a
// repeat makes a stale stamp read as a match, which is wrong silently. The tag
// only ever climbs, so the ceiling is where it has to be swept back. Half the
// positive range of std::int32_t, which is a pragmatic choice and not a derived
// one: nothing here stores anything but a tag, so the true ceiling is the type's
// own maximum, and the room left over is against a later layer wanting some of it.
constexpr std::int32_t TAG_CEILING = (1 << 30) - 1;

// Above this n, nothing is printed from inside the run: no initial state, no
// per-iteration trace. That output is for reading a small example by eye, and at any
// size worth calling large it is O(n) lines of O(n) each, so it is unreadable and slow
// to produce. What still prints at every size is the end of the run, the counters and
// the order, since each is O(1) lines and that is what the twin comparison comes down
// to. To watch a larger run, raise this.
constexpr std::uint32_t SHOW_THRESHOLD = 32;

// TWO GRAPHS, NOT ONE ALIAS. Both hold one list of int32 per vertex and differ only in what the
// entries MEAN: A[u] holds vertices and I[u] holds clique ids. An alias made them one type, so
// nothing but a variable name said which was which and nothing stopped one being passed for the
// other. Two classes cost the duplication below and buy a compiler that knows the difference.
//
// `mSize` is what makes them classes rather than pairs of vectors: it is the id space, a one
// dimensional size, so `std::uint32_t`. Holding it here is what keeps `n` out of `std::size_t` for
// the whole layer, and every vector length in the file is then bounded by it. Oblio's
// QuotientGraph owns both roles in one class and will get the same treatment iteratively.
class AdjacencyGraph {
public:
    explicit AdjacencyGraph(std::uint32_t size) : mSize(size), mAdjacency(size) {}
    // For the examples at the bottom, which are written as brace lists of neighbor lists.
    AdjacencyGraph(std::initializer_list<std::vector<std::int32_t>> rows)
        : mSize(static_cast<std::uint32_t>(rows.size())), mAdjacency(rows) {}

    std::uint32_t size() const { return mSize; }
    const std::vector<std::int32_t>& operator[](std::int32_t u) const { return mAdjacency[u]; }
    std::vector<std::int32_t>&       operator[](std::int32_t u)       { return mAdjacency[u]; }

private:
    std::uint32_t                          mSize;
    std::vector<std::vector<std::int32_t>> mAdjacency;
};

class IncidenceGraph {
public:
    explicit IncidenceGraph(std::uint32_t size) : mSize(size), mIncidence(size) {}

    std::uint32_t size() const { return mSize; }
    const std::vector<std::int32_t>& operator[](std::int32_t u) const { return mIncidence[u]; }
    std::vector<std::int32_t>&       operator[](std::int32_t u)       { return mIncidence[u]; }

private:
    std::uint32_t                          mSize;
    std::vector<std::vector<std::int32_t>> mIncidence;
};


// C[c] holds the members of clique c, and cliqueLive[c] says whether c exists.
// A clique id is the pivot that created it, so the id space is the vertex space.
class Cliques {
public:
    explicit Cliques(std::uint32_t size) : mSize(size), mMembers(size), mLive(size, false) {}

    // THE ACCESSORS FOLLOW THE MEMBERS, one for one and in the same order, so the four here read
    // as the four in the private section below. `size` and `live` exist at all because the m
    // prefix made the reads from outside visible: the show pass and the storage count had been
    // walking the fields directly.
    //
    // ONE ID SPACE, c in [0, size()), of which numLive() are live at any moment. `size` is the id
    // space and not the population, which is the sense every other size() in the tree has.
    std::uint32_t size() const    { return mSize; }
    std::uint32_t numLive() const { return mNumLive; }

    // Both overloads, so C[c] reads the same whether the reference is const or not. There was an
    // at() here instead of the const one until 2026-08-12, which compiled but borrowed a name the
    // standard library uses for the BOUNDS-CHECKED subscript. This at() checked nothing, so the
    // name promised something it did not do and the two spellings looked like a choice.
    const std::vector<std::int32_t>& operator[](std::int32_t c) const { return mMembers[c]; }
    std::vector<std::int32_t>&       operator[](std::int32_t c)       { return mMembers[c]; }

    bool live(std::int32_t c) const { return mLive[c]; }

    void create(std::int32_t c, std::vector<std::int32_t> members) {
        if (!mLive[c]) ++mNumLive;
        mLive[c]    = true;
        mMembers[c] = std::move(members);
    }
    void erase(std::int32_t c) {
        if (mLive[c]) --mNumLive;
        mLive[c] = false;
        mMembers[c].clear();
    }

private:
    // Declaration order is initialization order, so the two scalars come first and the vectors
    // they size come after. `mSize` is set by the caller and carries no default; `mNumLive` is
    // zero by the type, a fresh Cliques having nothing live, and says so where it is declared
    // rather than in a constructor that a second constructor could forget.
    std::uint32_t                          mSize;          // the id space, and both vectors' length
    std::uint32_t                          mNumLive = 0;   // how many of those ids are live now
    std::vector<std::vector<std::int32_t>> mMembers;
    std::vector<bool>                      mLive;
};

// The degree buckets, as the vendored codes hold them: one doubly linked list per
// degree, threaded through arrays of size n. Push, pop and splice are all O(1),
// which an ordered container cannot give. MMD spells these fwd/bwd and AMD
// Next/Last. The Python twin mirrors the same sequence with a list whose position
// 0 is the head, so both pick the same pivot.
//
// Set view: buckets[d] is the set of live vertices whose current degree is d, and
// the two operations here are add and discard. md5 has a third, refile, which is
// the two together with the degree written between them; from this layer up the
// eviction splits them, so there is nothing left for it to do and it is gone. A
// linked list gives both in O(1) and gives the head in O(1) too, which is
// everything the picker asks of it. What it does not give is a minimum, which is
// why minDegree walks. A sorted container would hand over the minimum directly and
// charge a log on every file, and files outnumber picks.
class Buckets {
public:
    explicit Buckets(std::uint32_t size)
        : mSize(size), mHead(size, NIL), mNext(size, NIL), mPrev(size, NIL), mFiled(size, false) {}

    // The accessors follow the members, one for one and in the same order. `size` and `next`
    // exist for the display walk, which is the only thing outside this class that needs to see a
    // chain; production's Buckets has no display and so has neither. `filed` is this layer's: a
    // batch WITHHOLDS vertices rather than refiling them, so the driver asks whether one is in
    // a bucket at all. Production derives the same answer from mPrev == UNFILED and keeps no
    // flag array. A DEGREE is one dimensional
    // and so std::uint32_t; the links are std::int32_t, carrying NIL.
    std::uint32_t size() const                 { return mSize; }
    std::int32_t  head(std::uint32_t d) const  { return mHead[d]; }
    bool          empty(std::uint32_t d) const { return mHead[d] == NIL; }
    std::int32_t  next(std::int32_t u) const   { return mNext[u]; }
    bool          filed(std::int32_t u) const  { return mFiled[u]; }

    void file(std::uint32_t d, std::int32_t u) {        // buckets[d].add(u), at the head
        mNext[u] = mHead[d];
        mPrev[u] = NIL;
        if (mHead[d] != NIL) mPrev[mHead[d]] = u;
        mHead[d] = u;
        mFiled[u] = true;
    }
    void unfile(std::uint32_t d, std::int32_t u) {      // buckets[d].discard(u)
        if (!mFiled[u]) return;                         // idempotent, as set.discard was
        if (mPrev[u] != NIL) mNext[mPrev[u]] = mNext[u];
        else mHead[d] = mNext[u];
        if (mNext[u] != NIL) mPrev[mNext[u]] = mPrev[u];
        mNext[u] = NIL;
        mPrev[u] = NIL;
        mFiled[u] = false;
    }

private:
    // INDEXED BY one thing, HOLDING another, and it is what they hold that sets the type. mHead is
    // subscripted by a DEGREE, which is one dimensional and so std::uint32_t in every signature
    // above; what it stores is a VERTEX or NIL, so its elements are std::int32_t. mNext and mPrev
    // are subscripted by a vertex and store one too.
    std::uint32_t             mSize;    // the id space, and every vector's length
    std::vector<std::int32_t> mHead;    // [degree] -> first live vertex of that degree, or NIL
    std::vector<std::int32_t> mNext;    // [vertex] -> next vertex toward the tail, or NIL
    std::vector<std::int32_t> mPrev;    // [vertex] -> previous vertex toward the head, or NIL
    std::vector<bool>         mFiled;   // [vertex] -> whether it is in a bucket at all
};

// I[u] cliques that contain u
// C[c] vertices that c contains

std::vector<std::int32_t> mmd3Neighbors(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
                                       const std::vector<bool>& eliminated,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u);


// Print a quotient graph: adjacency, incidence, cliques, in the order the
// structure holds them.
void mmd3Show(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
             const std::vector<std::uint32_t>& degrees, const std::string& title = "",
             const std::vector<bool>* eliminated = nullptr) {
    const std::uint32_t n = A.size();
    int width = static_cast<int>(std::to_string(n > 0 ? n - 1 : 0).size());
    std::vector<std::int32_t> liveVertices;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        if (eliminated == nullptr || !(*eliminated)[u]) liveVertices.push_back(u);
    std::size_t numLiveEdges = 0;
    for (std::int32_t u : liveVertices) numLiveEdges += A[u].size();
    numLiveEdges /= 2;
    std::size_t numLiveIncidences = 0;
    for (std::int32_t u : liveVertices) numLiveIncidences += I[u].size();
    std::uint32_t numLiveCliques = C.numLive();
    if (!title.empty()) std::cout << title << "\n";
    std::ostringstream liveVerticesText;
    if (eliminated == nullptr) liveVerticesText << n;
    else liveVerticesText << liveVertices.size() << " of " << n;
    std::cout << "num live vertices = " << liveVerticesText.str()
              << ", num live edges = " << numLiveEdges
              << ", num live cliques = " << numLiveCliques
              << ", storage = " << 2 * numLiveEdges << " + " << 2 * numLiveIncidences
              << " = " << 2 * (numLiveEdges + numLiveIncidences) << "\n";
    for (std::int32_t u : liveVertices) {
        std::ostringstream adjacencyText;
        bool first = true;
        for (std::int32_t v : A[u]) {
            adjacencyText << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::ostringstream incidenceText;
        first = true;
        for (std::int32_t c : I[u]) {
            incidenceText << (first ? "" : " ") << "c" << c;
            first = false;
        }
        std::cout << "  " << std::setw(width) << u << ": {" << adjacencyText.str()
                  << "} {" << incidenceText.str() << "} degree " << degrees[u] << "\n";
    }
    for (std::int32_t c = 0; c < static_cast<std::int32_t>(n); ++c) {
        if (!C.live(c)) continue;
        std::ostringstream cliqueMembersText;
        bool first = true;
        for (std::int32_t u : C[c]) {
            cliqueMembersText << (first ? "" : " ") << std::setw(width) << u;
            first = false;
        }
        std::cout << "  c" << c << ": {" << cliqueMembersText.str() << "}\n";
    }
    std::cout << "\n";
}

// Print the state arrays: degrees, buckets, min degree, members, eliminated,
// and the order so far.
void mmd3ShowState(const std::vector<std::uint32_t>& degrees, const Buckets& buckets,
                  std::uint32_t minDegree,
                  const std::vector<std::vector<std::int32_t>>& superMembers,
                  const std::vector<bool>& eliminated,
                  const std::vector<std::int32_t>& pivots, const std::string& title = "") {
    const std::uint32_t n = static_cast<std::uint32_t>(superMembers.size());
    int width = static_cast<int>(std::to_string(n > 0 ? n - 1 : 0).size());
    if (!title.empty()) std::cout << title << "\n";
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        std::string status;
        if (!eliminated[u]) status = "live";
        else if (!superMembers[u].empty()) status = "done";
        else status = "merged";
        std::ostringstream superMemberList;
        bool first = true;
        for (std::int32_t v : superMembers[u]) {
            superMemberList << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::cout << "  " << std::setw(width) << u << ": members ["
                  << superMemberList.str() << "] " << status << "\n";
    }
    std::ostringstream degreesText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        degreesText << (u == 0 ? "" : " ") << std::setw(width) << degrees[u];
    std::ostringstream superMembersText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        superMembersText << (u == 0 ? "" : " ") << "[";
        bool firstMember = true;
        for (std::int32_t v : superMembers[u]) {
            superMembersText << (firstMember ? "" : " ") << v;
            firstMember = false;
        }
        superMembersText << "]";
    }
    std::ostringstream eliminatedText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        eliminatedText << (u == 0 ? "" : " ") << std::setw(width) << (eliminated[u] ? 1 : 0);
    std::ostringstream bucketsText;
    bool firstBucket = true;
    for (std::uint32_t d = 0; d < buckets.size(); ++d) {
        if (buckets.empty(d)) continue;
        bucketsText << (firstBucket ? "" : "  ") << d << ": [";
        bool firstMember = true;
        for (std::int32_t v = buckets.head(d); v != NIL; v = buckets.next(v)) {
            bucketsText << (firstMember ? "" : " ") << v;
            firstMember = false;
        }
        bucketsText << "]";
        firstBucket = false;
    }
    std::cout << "  degrees: [" << degreesText.str() << "]\n";
    std::cout << "  buckets: " << (firstBucket ? "all empty" : bucketsText.str()) << "\n";
    std::cout << "  min degree: " << minDegree << "\n";
    std::cout << "  members: [" << superMembersText.str() << "]\n";
    std::cout << "  eliminated: [" << eliminatedText.str() << "]\n";
    std::cout << "  pivots: [";
    for (std::uint32_t k = 0; k < pivots.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << pivots[k];
    std::cout << "]\n";
    std::cout << "  order: [";
    bool firstOrder = true;
    for (std::int32_t pivot : pivots)
        for (std::int32_t u : superMembers[pivot]) {
            std::cout << (firstOrder ? "" : ", ") << u;
            firstOrder = false;
        }
    std::cout << "]\n\n";
}

// Entries actually stored. Each edge costs two, one per endpoint in A. Each
// incidence costs two as well, the clique id in I and the member in C. Watch
// the total fall monotonically; the naive graph's only rises.
std::size_t mmd3Storage(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C) {
    std::size_t total = 0;   // TWO DIMENSIONAL, a count of entries, so it stays wide
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(A.size()); ++u) total += A[u].size();
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(I.size()); ++u) total += I[u].size();
    for (std::int32_t c = 0; c < static_cast<std::int32_t>(C.size()); ++c)
        if (C.live(c)) total += C[c].size();
    return total;
}

// The neighbors of live vertex u: its explicit adjacency A[u] together with the
// members of every clique that contains u, minus u itself, which the cliques
// always carry. This is George and Liu's reachable set, and it is what the
// elimination graph would hold explicitly.
std::vector<std::int32_t> mmd3Neighbors(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
                                       const std::vector<bool>& eliminated,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u) {
    // In set terms this is one line, and it is worth keeping in view because the
    // code below is that line with the set taken away:
    //
    //     reach(u) = ( A[u] | C[c] for every c in I[u] ) - {u}
    //
    // The mark array IS the set. mark[v] == tag is the membership test, one
    // comparison; mark[v] = tag is the insertion, one store. So the union costs one
    // pass per source rather than a hash per member, and nothing is allocated.
    //
    // One pass per source, with the mark array doing the deduplication, so the
    // cost is linear in what is touched. Nothing is sorted: the order is the order
    // the sources were walked in. A vertex numbered by the prepass is skipped,
    // which is what marker[mn] = maxint does in genmmd.
    ++tag;
    std::vector<std::int32_t> neighbors;
    mark[u] = tag;                          // never its own neighbor
    for (std::int32_t v : A[u])
        if (!eliminated[v]) { mark[v] = tag; neighbors.push_back(v); }
    // THE ONE CHANGE FROM mmd2, AND IT IS NOT AN ALGORITHM. genmmd holds this list as a
    // linked list pushed at the head, `list[nb] = el; el = nb`, then reads from the head,
    // so the entry seen LAST is processed FIRST. We hold a vector and append, which is the
    // same set in the opposite order. Same set, same cost; what differs is which of two
    // equal candidates wins, and minimum degree is decided by exactly that. Walking
    // backwards restores genmmd's order without paying for a head insertion.
    //
    // This one is mmdelm's CLIQUE stack, and it is the deepest of the four: it fixes the
    // order of C[pivot], hence the content order of every list built from it.
    for (auto c = I[u].rbegin(); c != I[u].rend(); ++c)
        for (std::int32_t v : C[*c])
            if (mark[v] != tag && !eliminated[v]) { mark[v] = tag; neighbors.push_back(v); }
    return neighbors;
}

// The weighted degree of u: its neighbors counted in ORIGINAL vertices, since a
// neighbor may stand for several.
//
// Set view: sum of |superMembers[v]| over v in reach(u), which is |reach(u)| once
// every supervariable is expanded back to the vertices it stands for.
//
// The count of a supervariable is
// superMembers[v].size(), which is O(1), so no weight array is kept: md3 through
// mmd1 have none for the same reason. One becomes necessary only when the members
// are chains over a flat array, where a size stops being free.
std::uint32_t mmd3Degree(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
                       const std::vector<bool>& eliminated,
                       const std::vector<std::vector<std::int32_t>>& superMembers,
                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                       std::int32_t u) {
    std::uint32_t degree = 0;
    for (std::int32_t v : mmd3Neighbors(A, I, C, eliminated, mark, tag, u))
        degree += superMembers[v].size();
    return degree;
}

// Turn the pivot into a clique, then merge in every member it makes
// indistinguishable. Identical to md5Eliminate: this layer changes how often
// degrees are refreshed, not what an elimination does.
//
// Returns (neighbors, absorbedCliques, prunedEdges, mergedVertices): as in md2,
// plus the vertices folded into the pivot by mass elimination. The last three
// are reported for display; only neighbors is used by the caller.//
// Set view of the whole function, in the order the code does it:
//
//     C[pivot] = reach(pivot)                    absorb into C[pivot]
//     C        = C - I[pivot]                    reclaim I[pivot]
//     for u in C[pivot]:
//         A[u] = A[u] - C[pivot] - {pivot}       prune
//         I[u] = ( I[u] - I[pivot] ) | {pivot}   absorb into C[pivot], reclaim I[pivot]
//
// The new clique is C[pivot] and gets no name of its own, so the first line reads
// as what an elimination IS: the pivot stops being a vertex with a reachable set
// and becomes a clique holding that same set. The last line is the first two
// written on the I side, since u is in C[c] exactly when c is in I[u].
//
// Three set differences, and not one of them builds a set. Each is a single stamp
// of the subtrahend followed by one compaction pass over the minuend, which turns
// |A[u]| * |C[pivot]| comparisons into |A[u]| + |C[pivot]|.
//
// Mass elimination adds two more lines, and breaks the identity in the first one:
// from here C[pivot] is reach(pivot) minus what the pivot absorbed.
//
//     merged   = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
//     C[pivot] = C[pivot] - merged
std::tuple<std::vector<std::int32_t>, std::vector<std::int32_t>,
           std::vector<std::pair<std::int32_t, std::int32_t>>, std::vector<std::int32_t>>
mmd3Eliminate(AdjacencyGraph& A, IncidenceGraph& I, Cliques& C, std::vector<bool>& eliminated,
              std::vector<bool>& outmatched,
              std::vector<std::int32_t>& mark, std::int32_t& tag, std::int32_t pivot) {
    const std::vector<std::int32_t> neighbors = mmd3Neighbors(A, I, C, eliminated, mark, tag, pivot);
    const std::vector<std::int32_t> absorbedCliques = I[pivot];
    for (std::int32_t c : absorbedCliques)
        C.erase(c);
    C.create(pivot, neighbors);     // becomes the column pattern of the pivot

    // THREE TAGS IN THIS ELIMINATOR, and it is worth naming what each is about and what it
    // labels, since the mark array is shared and a tag is the only thing saying which set a stamp
    // belongs to. Two are md2's, about cliques, differing in WHICH SIDE of a clique they name:
    //
    //     pivotCliqueTag       about the PIVOT'S clique, labels its MEMBERS,  so stamps VERTICES
    //     absorbedCliquesTag   about the ABSORBED cliques, labels their IDS,  so stamps CLIQUE IDS
    //
    // Each side is what one loop below needs: pruning A[u] asks whether a VERTEX is in C[pivot],
    // pruning I[u] asks whether a CLIQUE ID is one of the absorbed. Each set is built in one pass
    // and then queried for free, and both loops are compactions in place.
    //
    // THE THIRD IS MASS ELIMINATION'S, further down, and it is the first tag in the family on the
    // MEMBER side alongside pivotCliqueTag: it labels the merged vertices so C[pivot] can be
    // compacted against them. Two consequences. It is CONDITIONAL, fired only when something
    // merged, so this layer's eliminate advance is a bound, at most 4, where md2's is exactly 3.
    // And md2's two tags could have shared one value, their sides being disjoint; here they could
    // not, because two of the three now stamp vertices.
    ++tag;
    const std::int32_t pivotCliqueTag = tag;
    for (std::int32_t v : neighbors) mark[v] = pivotCliqueTag;
    ++tag;
    const std::int32_t absorbedCliquesTag = tag;
    for (std::int32_t c : absorbedCliques) mark[c] = absorbedCliquesTag;

    std::vector<std::pair<std::int32_t, std::int32_t>> prunedEdges;

    // ONE SCRATCH BUFFER FOR EVERY COMPACTION, here and in the mass elimination below. Each use
    // is a filter into it followed by a swap, so after the swap it holds the list that was just
    // replaced; clear() then empties it while keeping that capacity, and the next fill reuses the
    // allocation. Named buffers would read no better and would allocate more. The uses are
    // labeled, since the same name means a different list a few lines apart.
    std::vector<std::int32_t> kept;
    for (std::int32_t u : neighbors) {
        outmatched[u] = false;      // mmdelm's bwd[rn] = 0: back in the running
        kept.clear();                            // KEPT IS ADJACENCY here: A[u] - C[pivot] - {pivot}
        for (std::int32_t v : A[u]) {
            if (v == pivot) continue;            // the pivot is no longer a variable
            if (eliminated[v]) continue;         // numbered by the prepass, gone for good
            if (mark[v] == pivotCliqueTag) {          // both ends inside the new clique
                if (u < v) prunedEdges.push_back({u, v});
                continue;                        // implicit now: drop the explicit copy
            }
            kept.push_back(v);
        }
        A[u].swap(kept);                         // what survives is A[u] - C[pivot] - {pivot}

        kept.clear();                            // KEPT IS INCIDENCE here: I[u] - I[pivot], + pivot
        for (std::int32_t c : I[u])
            if (mark[c] != absorbedCliquesTag) kept.push_back(c);
        kept.push_back(pivot);                   // u joins the new clique, id = pivot
        I[u].swap(kept);
    }

    // Mass elimination. u is INDISTINGUISHABLE from the pivot when the two have
    // the same closed neighborhood, mmd3Neighbors(u) | {u} == mmd3Neighbors(pivot)
    // | {pivot}, as it stood before the iteration. Equivalently, now that the clique is
    // formed, when everything u can still reach lies inside it. The test below is
    // a cheap sufficient condition for that: nothing explicit left and no clique
    // but the new one means u sees exactly what the pivot sees, so eliminating it
    // next would cost no fill. Fold it into the pivot now and strip it from the
    // cliques, since it is no longer a vertex.
    // merged = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
    std::vector<std::int32_t> mergedVertices;
    for (std::int32_t u : neighbors) {
        if (A[u].empty() && I[u].size() == 1 && I[u][0] == pivot) {
            I[u].clear();
            eliminated[u] = true;
            mergedVertices.push_back(u);
        }
    }
    if (!mergedVertices.empty()) {           // C[pivot] - merged, one compaction pass
        // THE MERGED TAG: about the vertices this step absorbed into the pivot, labeling them
        // directly. On the MEMBER side, like pivotCliqueTag, and the only conditional advance in
        // the file.
        ++tag;
        for (std::int32_t u : mergedVertices) mark[u] = tag;
        kept.clear();                        // KEPT IS CLIQUE MEMBERS here: C[pivot] - merged
        for (std::int32_t v : C[pivot])
            if (mark[v] != tag) kept.push_back(v);
        C[pivot].swap(kept);
    }

    A[pivot].clear();
    I[pivot].clear();
    eliminated[pivot] = true;
    return {neighbors, absorbedCliques, prunedEdges, mergedVertices};
}

// Multiple elimination: a batch of independent pivots per degree refresh.
//
// delta widens the batch to vertices within delta of the minimum degree, which
// buys still fewer refreshes for a real concession, since those vertices are not
// minimal. delta = 0 keeps the batch to true minima. A negative delta takes one
// pivot per iteration, which is md5's behavior reached through this code path.
// delta is signed: negative means one pivot per iteration. It is compared against a
// degree and its useful range stops at n - 1, so it is an index-like quantity by
// Oblio's rule, a std::int32_t rather than a count.
std::vector<std::int32_t> mmd3MinimumDegree(const AdjacencyGraph& G, std::int32_t delta = 0) {
    const std::uint32_t n = G.size();

    // An empty graph has no prepass to run and no bucket 1 to read, and production
    // returns here too, `if (size == 0) return std::vector<std::int32_t>()`. The
    // summary is written out rather than derived because at n = 0 every quantity in
    // it is zero by inspection; the cost is that a new counter has to be added in two
    // places, which is why there is exactly one line of it.
    if (n == 0) {
        std::cout << "n = 0, nnz(L) = 0 against nnz(tril A) = 0, fill = 0\n";
        std::cout << "iterations: 0\n";
        std::cout << "eliminations: 0\n";
        std::cout << "sum of |C[p]|: 0\n";
        std::cout << "degree computations: 0, degree updates: 0, bucket probes: 0, "
                     "prepass: 0, pair merges: 0, outmatched: 0, ncsub: 0\n";
        std::cout << "tag sweeps: 0\n";
        std::cout << "order: []\n";
        return std::vector<std::int32_t>();
    }

    std::size_t nnzTrilA = 0;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    AdjacencyGraph A = G;                                  // explicit vertex neighbors
    IncidenceGraph I(n);                                   // cliques that contain each vertex
    Cliques        C(n);      // clique id -> member list
    std::vector<std::int32_t> mark(n, NIL);       // scratch for membership, with tag
    std::int32_t tag = 0;
    // Calls to the eliminate procedure, one per pivot. Not the count of vertices
    // removed: a pivot can carry mass-merged vertices out with it, and from mmd1 up
    // an iteration batches several eliminations before one degree update pass. The three
    // counts coincide only where both of those are absent.
    std::uint32_t numEliminations = 0;
    // Summed over the eliminations, |C[p]| being the new clique AFTER the trim, so
    // in supernodal terms the update rather than the front. It is the raw reach of
    // the eliminations, undeduplicated: where a layer deduplicates, the degree
    // update count comes out below this, and the gap is what the batching saved.
    // In md2 it is nnz(L) - n, there being no mass elimination to shrink a clique.
    std::size_t numCliqueEntries = 0;
    std::uint32_t numIterations = 0;                    // batches, the metric this layer adds
    std::vector<std::vector<std::int32_t>> superMembers(n);   // for the expansion
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        superMembers[u].push_back(u);
    std::vector<bool> eliminated(n, false);
    std::vector<bool> outmatched(n, false);       // withheld from the buckets, not merged
    std::vector<std::int32_t> pivots;             // the order over supervariables
    std::uint32_t numEliminatedVertices = 0;                // a counter, not a scan of eliminated
    std::size_t nnzL = 0;

    std::vector<std::uint32_t> degrees(n);          // a degree counts, so it measures
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) degrees[u] = A[u].size();
    // Only the updates are counted. The total, including the initial pass over all
    // n vertices, is that plus n, so the report derives it rather than keeping a
    // second counter that could drift from this one.
    std::size_t numDegreeUpdates = 0;
    // Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    // is here as the witness that the guard is inert rather than as a statistic.
    std::size_t numTagSweeps = 0;

    // mmdint files a degree-0 vertex under degree 1, `if(dg==0)dg=1`, so the
    // bucket a vertex sits in is max(degree, 1) rather than its degree. From here
    // degrees[] holds that filed value, which is what MMD compares and files by.
    //
    // THE FLOOR IS HERE TO MATCH THE VENDORED OUTPUT, not to find the prepass
    // vertices. Taking bucket 0 and then bucket 1 finds the same vertices, and the
    // fill comes out identical; what changes is the ORDER within the prepass, since
    // the floor puts both degrees on ONE list where they interleave by insertion and
    // separate buckets group them by degree. Measured: same prepass set and same
    // nnz(L) on all 300 random graphs tried, different permutation on 212 of them.
    // So this is a tie-break of the same kind as mmd3's four, and dropping it would
    // be a fifth alignment defect that no fill check could see.
    // n + 1 and not n, which is how production sizes it: `mHead(size + 1, NIL)` beside
    // `mNext(size, NIL)`, because a head is indexed by a DEGREE and a link by a VERTEX.
    // The two index spaces are not the same size, and the floor above is what makes the
    // difference bite: it files a degree-0 vertex under 1, so at n = 1 index 1 has to
    // exist and holds the only vertex there is. Bucket 0 goes unused from here on, the
    // floor having taken it out of the range.
    //
    // This layer's Buckets sizes all four of its arrays from the one argument, so the
    // three vertex-indexed ones come out a slot longer than they need. Production does
    // not, and the honest repair is there rather than here; the cost of the slack is one
    // int32 per array on a structure built once per ordering.
    Buckets buckets(n + 1);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        degrees[u] = std::max<std::uint32_t>(degrees[u], 1);
        buckets.file(degrees[u], u);
    }
    std::uint32_t minDegree = n > 0 ? *std::min_element(degrees.begin(), degrees.end()) : 0;
    std::size_t numBucketProbes = 0;
    std::size_t ncsub = 0;                        // genmmd's subscript estimate
    std::size_t pairMerges = 0;                   // two-source merges, coarser supervariables
    std::size_t outmatchedCount = 0;              // withheld rather than refiled

    // NOT PRODUCTION: display only. The trace is what makes these files teachable and
    // is the whole reason they exist; nothing downstream reads it.
    if (n <= SHOW_THRESHOLD) {
        mmd3Show(A, I, C, degrees, "start: every edge explicit, no clique yet", &eliminated);
        mmd3ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
    }

    // ---- the PREPASS -------------------------------------------------------
    // Number everything in bucket 1, which after the floor above holds the
    // isolated and the degree-1 vertices together, and then leave the bucket
    // empty. Nothing is eliminated in the quotient-graph sense: no clique is
    // formed, nothing is pruned, and the neighbors keep degrees that still count
    // these vertices. That staleness is the point, and it is what genmmd does.
    std::vector<std::int32_t> prepassVertices;
    for (std::int32_t v = buckets.head(1); v != NIL; v = buckets.next(v))
        prepassVertices.push_back(v);
    for (std::int32_t u : prepassVertices) {
        buckets.unfile(degrees[u], u);
        std::uint32_t externalDegree = 0;
        for (std::int32_t v : A[u]) if (!eliminated[v]) ++externalDegree;
        nnzL += externalDegree + 1;
        eliminated[u] = true;
        pivots.push_back(u);
        ++numEliminatedVertices;
    }
    if (!prepassVertices.empty()) {
        std::ostringstream prepassText;
        for (std::uint32_t k = 0; k < prepassVertices.size(); ++k)
            prepassText << (k == 0 ? "" : ", ") << prepassVertices[k];
        std::ostringstream title;
        title << "prepass: numbered " << prepassVertices.size() << ": " << prepassText.str();
        // NOT PRODUCTION: display only. The trace is what makes these files teachable and
        // is the whole reason they exist; nothing downstream reads it.
        if (n <= SHOW_THRESHOLD) {
            mmd3Show(A, I, C, degrees, title.str(), &eliminated);
            mmd3ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
        }
    }
    if (n > 2) minDegree = 2;                  // head[1] = 0, and mdeg starts at 2

    while (numEliminatedVertices < n) {
        while (buckets.empty(minDegree)) {      // walk up to the first live bucket
            ++minDegree;
            ++numBucketProbes;
        }
        ++numBucketProbes;

        // ---- one BATCH, no degree refreshed inside it ----------------------
        // Take pivots from buckets [minDegree, minDegree + delta]. Eviction is
        // what keeps them independent: eliminating a pivot pulls every vertex it
        // reached out of the buckets, so whatever is still filed was not reached,
        // hence is not adjacent to anything taken this iteration.
        //
        // Set view of the invariant the eviction maintains, where reached is the
        // union of C[p] over the pivots taken so far:
        //
        //     filed = live - reached,  so  batch & reached == {}
        //
        // No set is built for either side. Membership in filed is the filed flag,
        // and nothing else is needed: unlike mmd1 this layer carries no evicted
        // list, the refresh below re-deriving its vertices from the cliques.
        // Clamped: a degree is at most n - 1, so a wider window would walk the
        // bucket array off its end.
        std::size_t batchLimit = minDegree;      // delta > 0 here, so no narrowing
        if (delta > 0)
            batchLimit = std::min(minDegree + delta, static_cast<std::uint32_t>(n) - 1);
        std::vector<std::int32_t> batch;
        while (true) {
            if (buckets.empty(minDegree)) {     // this degree is drained
                if (minDegree >= batchLimit) break;
                ++minDegree;
                ++numBucketProbes;
                continue;
            }
            std::int32_t pivot = buckets.head(minDegree);
            std::uint32_t degree = degrees[pivot];
            buckets.unfile(degree, pivot);

            // Sweep the tag back before it can wrap. Two sites in this layer, one
            // before each region that advances the tag, and each placed where nothing
            // in mark is live. This one is INSIDE the batch loop rather than before
            // it, since a batch takes several pivots and each calls the eliminator.
            // Safe between eliminations because the eviction that follows stamps
            // filed, which is a separate array. Not inside
            // mmd3Eliminate, which holds three stamps live in turn: pivotCliqueTag and
            // absorbedCliquesTag across the prune loop, then the merged set across the
            // C[pivot] compaction. Never observed to fire.
            if (tag > TAG_CEILING) {
                std::fill(mark.begin(), mark.end(), NIL);
                tag = 0;
                ++numTagSweeps;
            }
            auto [neighbors, absorbedCliques, prunedEdges, mergedVertices] =
                mmd3Eliminate(A, I, C, eliminated, outmatched, mark, tag, pivot);
            ++numEliminations;
            numCliqueEntries += C[pivot].size();
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminatedVertices += 1 + mergedVertices.size();
            for (std::int32_t u : mergedVertices) {   // the pivot now stands for them too
                superMembers[pivot].insert(superMembers[pivot].end(),
                                           superMembers[u].begin(), superMembers[u].end());
                superMembers[u].clear();
                buckets.unfile(degrees[u], u);
                degrees[u] = 0;
            }
            degrees[pivot] = 0;

            for (std::int32_t u : C[pivot]) {   // EVICT, with a stale degree
                buckets.unfile(degrees[u], u);
            }

            // A supervariable of size w is w consecutive columns of L. Its
            // external degree is what remains of the clique after the merges,
            // since a merged vertex joins the supervariable instead of
            // neighboring it. The first column holds ext + w - 1 entries below
            // its diagonal, the next ext + w - 2, down to ext, and each column
            // contributes its own diagonal.
            ncsub += degree + superMembers[pivot].size() - 2;   // genmmd's *ncsub
            std::uint32_t superSize = superMembers[pivot].size();
            std::uint32_t externalDegree = 0;
            for (std::int32_t v : C[pivot])
                if (!eliminated[v]) externalDegree += superMembers[v].size();
            // ONE FACTOR WIDENED: a product of two one-dimensional quantities is TWO dimensional, so
            // it is formed in std::size_t. Widening cannot be done after the multiply the way
            // narrowing is done after a subtraction; one operand is enough, the other promoting to
            // meet it. Both factors are bounded by n, so the product reaches n^2.
            nnzL += static_cast<std::size_t>(superSize) * externalDegree
                  + static_cast<std::size_t>(superSize) * (superSize - 1) / 2 + superSize;

            // NOT PRODUCTION: display only, and silent above the threshold.
            if (n <= SHOW_THRESHOLD) {
                std::ostringstream absorbedCliquesText;
                if (absorbedCliques.empty()) {
                    absorbedCliquesText << "none";
                } else {
                    bool first = true;
                    for (std::int32_t c : absorbedCliques) {
                        absorbedCliquesText << (first ? "" : ", ") << "c" << c;
                        first = false;
                    }
                }
                std::ostringstream prunedEdgesText;
                if (prunedEdges.empty()) {
                    prunedEdgesText << "none";
                } else {
                    bool first = true;
                    for (auto [u, v] : prunedEdges) {
                        prunedEdgesText << (first ? "" : ", ") << u << "-" << v;
                        first = false;
                    }
                }
                std::ostringstream mergedVerticesText;
                if (mergedVertices.empty()) {
                    mergedVerticesText << "none";
                } else {
                    bool first = true;
                    for (std::int32_t u : mergedVertices) {
                        mergedVerticesText << (first ? "" : ", ") << u;
                        first = false;
                    }
                }
                std::ostringstream evictedText;
                if (C[pivot].empty()) {
                    evictedText << "none";
                } else {
                    bool first = true;
                    for (std::int32_t u : C[pivot]) {
                        evictedText << (first ? "" : ", ") << u;
                        first = false;
                    }
                }
                std::cout << "iteration " << numIterations << ": eliminate " << pivot << " (degree "
                          << degree << ", size " << superSize << ", external degree "
                          << externalDegree << "), absorbed cliques: "
                          << absorbedCliquesText.str() << ", pruned edges: "
                          << prunedEdgesText.str() << ", merged vertices: "
                          << mergedVerticesText.str() << ", evicted: " << evictedText.str() << "\n";
            }
            if (numEliminatedVertices >= n) break;      // genmmd's num + qsize[mn] > neqns:
                                                // nothing left to update
            if (delta < 0) break;               // one pivot per iteration, as md5 does
        }

        // ---- one REFRESH, walked CLIQUE BY CLIQUE ------------------------
        // mmdupd walks the cliques this iteration created, not the vertices it
        // reached, and computes cliqueWeight once per clique: the size of that clique,
        // which every member of it reaches in full. A member with exactly one other
        // source goes on the twoSourceQueue and is answered from cliqueWeight plus that source;
        // everything else goes on manySourceQueue and pays for the full union. Same degrees,
        // different work, and a different filing order.
        //
        // This is also why no evicted list is carried. mmd1 accumulates one during
        // the batch and walks it here; walking cliques re-derives the same vertices
        // from C[clique], deduplicating with the filed flag, so the list and the
        // second stamp array it needs both go. genmmd makes the same trade, chaining
        // its new cliques in `list` and building no vertex set at all.
        std::vector<std::int32_t> refreshedVertices;
        std::vector<std::int32_t> cliqueMembers, twoSourceQueue, manySourceQueue;
        // The second site, before the refresh, and OUTSIDE the clique loop rather
        // than inside it. cliqueTag is stamped once per clique and read all the
        // way through that clique's two-source walk, where it decides both the merge
        // and the outmatched case, with vertexTag fresh per vertex nested inside it.
        // Two levels live at once, which is mmdupd's mt against its tag, so a sweep
        // within an clique erases marks about to be read.
        if (tag > TAG_CEILING) {
            std::fill(mark.begin(), mark.end(), NIL);
            tag = 0;
            ++numTagSweeps;
        }
        // The driver's clique list, `list[mn] = ehead; ehead = mn`, so the LAST pivot of a
        // batch is the FIRST clique refreshed. See mmd3Neighbors.
        for (auto cliqueIt = batch.rbegin(); cliqueIt != batch.rend(); ++cliqueIt) {
            const std::int32_t clique = *cliqueIt;
            cliqueMembers.clear();
            for (std::int32_t u : C[clique])
                if (!eliminated[u]) cliqueMembers.push_back(u);
            ++tag;                              // cliqueWeight's members, marked once
            const std::int32_t cliqueTag = tag;
            for (std::int32_t u : cliqueMembers) mark[u] = cliqueTag;
            std::uint32_t cliqueWeight = 0;
            for (std::int32_t u : cliqueMembers) cliqueWeight += superMembers[u].size();

            // Set view of the split. reach(u) has |A[u]| + |I[u]| sources once the
            // new clique is counted, so |A[u]| + |I[u]| - 1 == 1 says everything u
            // reaches lies in this clique plus ONE other source. That is the case a
            // union is not needed for: cliqueWeight already counts the clique, and the one
            // other source is walked directly.
            twoSourceQueue.clear();
            manySourceQueue.clear();
            for (std::int32_t u : cliqueMembers) {
                if (buckets.filed(u) || outmatched[u]) continue;  // done, or withheld
                std::uint32_t otherSources = A[u].size() + I[u].size() - 1;
                (otherSources == 1 ? twoSourceQueue : manySourceQueue).push_back(u);
            }

            // mmdupd's q2h chain, head-pushed there. See mmd3Neighbors.
            for (auto uit = twoSourceQueue.rbegin(); uit != twoSourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;
                if (eliminated[u] || outmatched[u]) continue;   // merged or withheld
                                                               // by an earlier two-source
                // Everything u reaches is in the clique or comes from its one
                // other source. cliqueWeight counts the clique, minus u itself. Two mark
                // levels, as mmdupd has: cliqueTag says "already in cliqueWeight" and
                // survives the whole clique, while vertexTag is fresh per vertex,
                // so one two-source vertex cannot hide a neighbor from the next.
                ++tag;
                const std::int32_t vertexTag = tag;
                // cliqueWeight is kept WHOLE here and u's own weight subtracted at the end, which
                // is genmmd's `dg - qsize[en] + 1` and not the same as subtracting it now. The
                // walk below can MERGE a vertex into u, and genmmd's merge does
                // `qsize[en] += qsize[nd]` in that same walk, so the weight it subtracts is
                // the one AFTER the merge. Subtracting first files a supervariable one bucket
                // too high per merged vertex, so it is not picked as early as its size has
                // earned. See the ledger, entry 5.
                std::uint32_t degree = cliqueWeight;
                for (std::int32_t v : A[u]) {
                    if (eliminated[v] || mark[v] == vertexTag) continue;
                    if (mark[v] == cliqueTag) continue;        // already in cliqueWeight
                    mark[v] = vertexTag;
                    degree += superMembers[v].size();
                }
                for (std::int32_t c : I[u]) {
                    if (c == clique) continue;
                    for (std::int32_t v : C[c]) {
                        if (v == u || eliminated[v] || mark[v] == vertexTag) continue;
                        if (mark[v] == cliqueTag) {
                            // v is in the new clique AND in this same other
                            // source, so it sees at least what u sees.
                            if (buckets.filed(v) || outmatched[v]) continue;
                            if (A[v].size() + I[v].size() - 1 == 1) {
                                // v is two-source too, so its only other source is this
                                // one: identical reach, and u absorbs it. Set view:
                                // reach(u) | {u} == reach(v) | {v}, decided without
                                // forming either side, because both sets are pinned
                                // to the same two sources.
                                superMembers[u].insert(superMembers[u].end(),
                                                       superMembers[v].begin(),
                                                       superMembers[v].end());
                                superMembers[v].clear();
                                eliminated[v] = true;
                                ++numEliminatedVertices;
                                ++pairMerges;
                            } else {
                                // v reaches more than u does, so it can never be
                                // the minimum first. Withhold it from the buckets.
                                // Set view: reach(u) <= reach(v), a containment
                                // rather than an equality, so v is withheld and not
                                // merged.
                                outmatched[v] = true;
                                ++outmatchedCount;
                            }
                            continue;
                        }
                        mark[v] = vertexTag;
                        degree += superMembers[v].size();
                    }
                }
                degrees[u] = std::max<std::uint32_t>(
                    degree - static_cast<std::uint32_t>(superMembers[u].size()) + 1, 1);
                buckets.file(degrees[u], u);
                refreshedVertices.push_back(u);
            }

            // mmdupd's qxh chain, the same shape. See mmd3Neighbors.
            for (auto uit = manySourceQueue.rbegin(); uit != manySourceQueue.rend(); ++uit) {
                const std::int32_t u = *uit;
                if (eliminated[u] || outmatched[u]) continue;
                std::uint32_t degree = mmd3Degree(A, I, C, eliminated, superMembers, mark, tag, u);
                degrees[u] = std::max<std::uint32_t>(degree + 1, 1);  // dg - qsize + 1
                buckets.file(degrees[u], u);
                refreshedVertices.push_back(u);
            }
        }
        numDegreeUpdates += refreshedVertices.size();
        for (std::int32_t u : refreshedVertices) minDegree = std::min(minDegree, degrees[u]);
        ++numIterations;

        // NOT PRODUCTION: display only, and silent above the threshold. Built INSIDE
        // the guard, as the per-elimination line above is, so a run above the
        // threshold formats nothing.
        if (n <= SHOW_THRESHOLD) {
            std::ostringstream batchText;
            for (std::uint32_t k = 0; k < batch.size(); ++k)
                batchText << (k == 0 ? "" : ", ") << batch[k];
            std::ostringstream refreshedVerticesText;
            if (refreshedVertices.empty()) {
                refreshedVerticesText << "none";
            } else {
                bool first = true;
                for (std::int32_t u : refreshedVertices) {
                    refreshedVerticesText << (first ? "" : ", ") << u;
                    first = false;
                }
            }
            std::ostringstream title;
            title << "iteration " << (numIterations - 1) << " done: batch of " << batch.size() << ": "
                  << batchText.str() << ", refreshed vertices: " << refreshedVerticesText.str();
            mmd3Show(A, I, C, degrees, title.str(), &eliminated);
            mmd3ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
        }
    }

    // mmdnum, and the last thing that differed from genmmd. A supervariable's members are
    // indistinguishable by construction, so their order among themselves cannot change the
    // fill; it does change the permutation, which is why every grid matched on FILL and on
    // the PIVOT SEQUENCE while the printed order still diverged. genmmd numbers them root
    // first, then by ASCENDING VERTEX INDEX, and gets that from a single ascending scan
    // rather than a sort. See the ledger, entry 6.
    std::vector<std::int32_t> rootOf(n, NIL);
    for (std::int32_t pivot : pivots)
        for (std::int32_t u : superMembers[pivot]) rootOf[u] = pivot;
    std::vector<std::vector<std::int32_t>> membersOf(n);   // ascending, by the scan
    for (std::int32_t v = 0; v < static_cast<std::int32_t>(n); ++v)
        if (rootOf[v] != NIL && rootOf[v] != v) membersOf[rootOf[v]].push_back(v);

    std::vector<std::int32_t> order;
    for (std::int32_t pivot : pivots) {
        order.push_back(pivot);
        for (std::int32_t u : membersOf[pivot]) order.push_back(u);
    }
    std::cout << "n = " << n << ", nnz(L) = " << nnzL
              << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "iterations: " << numIterations << "\n";
    std::cout << "eliminations: " << numEliminations << "\n";
    std::cout << "sum of |C[p]|: " << numCliqueEntries << "\n";
    std::cout << "degree computations: " << (numDegreeUpdates + n)
              << ", degree updates: " << numDegreeUpdates
              << ", bucket probes: " << numBucketProbes
              << ", prepass: " << prepassVertices.size()
              << ", pair merges: " << pairMerges
              << ", outmatched: " << outmatchedCount
              << ", ncsub: " << ncsub << "\n";
    std::cout << "tag sweeps: " << numTagSweeps << "\n";
    std::cout << "order: [";
    for (std::uint32_t k = 0; k < order.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << order[k];
    std::cout << "]\n";
    return order;
}


// A square grid graph, four-neighbor, for running the counters at a size the seven examples
// cannot reach. It is here rather than among them because it is not an example: nothing about it
// illustrates a mechanism, and above SHOW_THRESHOLD its trace is not printed at all.
//
// It must match the Python twin's grid_graph exactly, vertex for vertex, or `make test` would be
// diffing two different problems.
static AdjacencyGraph gridGraph(int side) {
    const int n = side * side;
    AdjacencyGraph graph(static_cast<std::uint32_t>(n));
    for (int r = 0; r < side; ++r)
        for (int c = 0; c < side; ++c) {
            const int u = r * side + c;
            if (r > 0)        graph[u].push_back(u - side);
            if (c > 0)        graph[u].push_back(u - 1);
            if (c + 1 < side) graph[u].push_back(u + 1);
            if (r + 1 < side) graph[u].push_back(u + side);
        }
    return graph;
}

void run(const std::string& name, const AdjacencyGraph& G) {
    std::cout << "=== " << name << " ===\n";
    mmd3MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Grid mode, spelled the same way in every layer and in both twins, so `make test` can diff
    // at a size the seven examples cannot reach. Nothing is filtered: the run is silent above
    // SHOW_THRESHOLD and prints its closing lines as always.
    //
    //   ./mmd3_cpp grid 22
    if (argc > 2 && std::string(argv[1]) == "grid") {
        const int side = std::atoi(argv[2]);
        std::cout << "=== grid " << side << "x" << side << " (n = " << side * side << ") ===\n";
        mmd3MinimumDegree(gridGraph(side));
        return 0;
    }

    // The same three graphs as md1 and md2.
    //
    //   graph1, a 4-cycle: eliminating any vertex forces its two neighbors
    //   together, so it is the smallest graph that fills (one fill edge).
    //
    //      0---1          edges: 0-1 1-2 2-3 3-0
    //      |   |
    //      3---2
    //
    //   graph2, uneven degrees so the picker actually chooses; it fills twice.
    //
    //        0            edges: 0-1 0-2 1-3 2-4
    //       / \                  3-4 3-5 4-5
    //      1   2
    //      |   |
    //      3---4
    //       \ /
    //        5
    //
    //   graph3, twelve vertices: a path 0-1-...-11 with eight extra edges. Big
    //   enough that cliques grow past two members, which is where the quotient
    //   graph starts to pay, and its elimination order is not the identity.
    //
    //      edges: 0-1 0-3 0-8 1-2 1-6 1-8 2-3 2-5 3-4 4-5
    //             5-6 5-9 6-7 6-10 7-8 8-9 9-10 10-11
    AdjacencyGraph graph1 = {
        {1, 3}, {0, 2}, {1, 3}, {0, 2},
    };
    AdjacencyGraph graph2 = {
        {1, 2}, {0, 3}, {0, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4},
    };
    AdjacencyGraph graph3 = {
        {1, 3, 8},        // 0
        {0, 2, 6, 8},     // 1
        {1, 3, 5},        // 2
        {0, 2, 4},        // 3
        {3, 5},           // 4
        {2, 4, 6, 9},     // 5
        {1, 5, 7, 10},    // 6
        {6, 8},           // 7
        {0, 1, 7, 9},     // 8
        {5, 8, 10},       // 9
        {6, 9, 11},       // 10
        {10},             // 11
    };

    // graph4, eight vertices and fourteen edges. Denser than the others, and here
    // for one specific reason: it is the smallest graph we could find on which
    // AMD's degree BOUND is ever loose. The bound overcounts only when a vertex
    // belongs to two cliques that overlap outside the new one, which needs enough
    // eliminations to have made several cliques and enough fill for them to
    // intersect. Every connected graph on five or six vertices is exact (checked
    // exhaustively), and so are graph1 to graph3, so without this one the amd
    // trace would never show the approximation approximating. The other layers use
    // it as an ordinary denser test.
    //
    //   edges: 0-2 0-3 0-4 0-7 1-3 1-4 1-6 1-7 2-3 2-5 3-6 3-7 4-5 5-6
    AdjacencyGraph graph4 = {
        {2, 3, 4, 7},     // 0
        {3, 4, 6, 7},     // 1
        {0, 3, 5},        // 2
        {0, 1, 2, 6, 7},  // 3
        {0, 1, 5},        // 4
        {2, 4, 6},        // 5
        {1, 3, 5},        // 6
        {0, 1, 3},        // 7
    };

    // graph5, five vertices and four edges, two paths joined at 4: 2-1-4-0-3.
    // Small and fill free, and here for one reason: it is the smallest graph on
    // which mmd3's merge test declines a genuine supervariable. At the iteration whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test mmd3Neighbors(A, I, C, u) contained in
    // C[pivot] would merge it. See the README section on mass elimination.
    //
    //   edges: 0-3 0-4 1-2 1-4
    AdjacencyGraph graph5 = {
        {3, 4},           // 0
        {2, 4},           // 1
        {1},              // 2
        {0},              // 3
        {0, 1},           // 4
    };

    // graph6, six vertices and eight edges. Here because one small graph carries
    // three things at once. Its supervariable {0, 4} is a supernode but NOT a
    // fundamental one: the elimination forest is 2 -> 1 -> 4 and 3 -> 0 -> 4, so
    // 4 already has 1 as a child when 0 merges into it. The merge happens at iteration
    // 2 of 5, so the run continues afterwards and the selection degree, 3 over
    // {2, 3, 4}, differs from the external degree, 2 over {2, 3}, with the
    // difference being the size of what merged. And superMembers ends with a hole
    // in the middle, slot 4 empty between two used ones, while no pivot equals
    // its own iteration number. See the README sections on mass elimination and on
    // external degree.
    //
    //   edges: 0-2 0-3 0-4 1-3 2-3 2-4 2-5 3-4
    AdjacencyGraph graph6 = {
        {2, 3, 4},        // 0
        {3},              // 1
        {0, 3, 4, 5},     // 2
        {0, 1, 2, 4},     // 3
        {0, 2, 3},        // 4
        {2},              // 5
    };

    // graph7, five vertices and six edges. The pairwise case: at the iteration whose
    // pivot is 0 and whose clique is {2, 4}, vertices 2 and 4 are
    // indistinguishable FROM EACH OTHER, both reaching the same closed
    // neighborhood, yet neither is absorbable into the pivot, since each still
    // reaches 3 from outside the clique. No test framed against the pivot finds
    // them, and the exact test does not help either: both orders are 1 0 (2 3 4).
    // Catching such pairs needs a comparison between candidates, which is what
    // amd's hashing does. See the README section on detecting supervariables
    // against each other.
    //
    //   edges: 0-1 0-2 0-4 1-4 2-3 2-4 3-4
    AdjacencyGraph graph7 = {
        {1, 2, 4},        // 0
        {0, 4},           // 1
        {0, 3, 4},        // 2
        {2, 4},           // 3
        {0, 1, 2, 3},     // 4
    };

    std::vector<std::pair<std::string, AdjacencyGraph>> examples = {
        {"graph1", graph1}, {"graph2", graph2},
        {"graph3", graph3}, {"graph4", graph4},
        {"graph5", graph5}, {"graph6", graph6},
        {"graph7", graph7},
    };

    // All of them by default. To run just one, pass its number: ./mmd3_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
