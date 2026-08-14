// mmd1.cpp -- multiple minimum degree.
//
// md5 finished the cheap wins. It has the quotient graph, supervariables,
// maintained degrees and degree buckets, and it returns exactly the ordering md1
// returns, only far faster. Everything left costs something.
//
// This is the first layer that changes the ANSWER. Section 5.11 of
// archive/sparse_factorization.md.
//
// The idea, from Liu (1985), is the M in MMD. Refreshing degrees is the expensive
// iteration, so do it less often: eliminate a whole INDEPENDENT SET of least-degree
// vertices before refreshing anything. Non-adjacent pivots cannot disturb each
// other's degrees, so every pivot in a batch is still a true minimum-degree
// vertex when it is taken.
//
// We never search for the independent set. It falls out of the bookkeeping:
// eliminating a pivot EVICTS every vertex it reached from the degree buckets, so
// whatever is still sitting in the bucket was not reached, hence is non-adjacent
// to everything already taken this iteration.
//
// WHAT THIS GIVES UP, and it is not what one would guess. The pivots are exact,
// but the vertices the batch evicted are invisible for the rest of the iteration, so
// the choice is made among the untouched remainder rather than among all
// candidates. The batch does not pick a worse vertex, it picks a different vertex
// OF THE SAME DEGREE. Minimum degree is famously sensitive to tie-breaks, so the
// fill moves by a fraction of a percent, in either direction.
//
// WHAT IS HERE, AND WHAT MMD2 ADDS. This file is the idea alone. Everything else
// genmmd does is deliberately left to mmd1's successor, which completes it:
//
//   - the PREPASS that numbers degree 0 and 1 vertices before the main loop,
//     leaving their neighbors' degrees stale (genmmd, the loop over head[1])
//   - mmdupd's q2h path. mmdelm stashes each reached vertex's pruned adjacency
//     count as fwd[rn] = nq+1, and mmdupd routes the nq==1 cases into a separate
//     list where it merges indistinguishable PAIRS. The merge test here catches
//     only vertices indistinguishable from the pivot, so MMD's supervariables are
//     at least as coarse as ours and sometimes coarser.
//   - OUTMATCHED marking, bwd[nd] = -maxint, which takes a vertex out of the
//     degree lists without merging it.
//   - the filing convention: MMD files at `dg - qsize[en] + 1` floored at 1, so
//     its least bucket is 1 where ours is 0, and it never uses bucket 0. Plus the
//     ncsub subscript statistic.
//
// NO WEIGHT ARRAY, for the same reason md3 through md5 have none: mass elimination
// merges only into the PIVOT, which is eliminated in the same call, so no live
// vertex ever stands for more than one original vertex, and a supervariable's size
// is superMembers[pivot].size() whenever it is wanted. mmd2 needs one, because its
// q2h merge folds a vertex into a LIVE one.
//
// TIE-BREAKS. Our buckets are index-ordered, *buckets[minDegree].begin(), which is
// md5's convention and the reason md1 through md5 agree. MMD's degree lists are
// linked chains prepended at head[dg], so its bucket is a stack and the winner is
// whatever was pushed last, which after construction is the highest-numbered
// vertex of that degree. There is no quality claim behind it: prepending is the
// cheap end of a linked list. We keep our convention and the orderings differ in
// ties; see the README.
//
// The tag/marker machinery with its maxint overflow reset is not modeled at all.
// It exists because the marks live in reusable integer arrays; ours are the mark
// and evictedMark arrays, which do the same job with an explicit tag.
//
//
// COMPLEXITY, AND ONE PLACE THE PYTHON PAYS MORE THAN THE C++. The goal is the
// same asymptotic cost as the vendored routines, without their coding style. Two
// things were wrong and are fixed: the driver loop counts eliminations rather than
// scanning `eliminated` (O(n) per iteration before, O(1) now), and the mass elimination
// block strips a merged vertex from C[pivot] alone rather than from every clique,
// which is sound because I[u] was {pivot}. On a 20 by 20 grid those two cost 14800
// and 4247 elementary iterations before, against 34 and 47 after, with the real
// neighbor work at 26408.
//
// The containers are flat: A and I are sorted vectors, C is indexed by clique id,
// membership comes from a mark array with a tag, and a bucket is a linked list,
// head[d] with next and prev over n, so filing, unfiling and popping are O(1).
// With that this file performs the same operations at the same cost as the
// vendored genmmd; what it does not yet have is genmmd's remaining features,
// which are mmd2's business.
//
// The Python mirrors the buckets with a list whose position 0 is the head, so both
// twins hold the same sequence and pick the same pivot, and it pays O(bucket) for
// insert and remove where the C++ splices in O(1). That is the one place the
// Python is asymptotically worse than its twin; A, I and C stay sets there,
// because set algebra is what makes the layer readable.
//
// Build:  g++ -std=c++17 -O3 mmd1.cpp -o mmd1_cpp  (or: make)
// Run:    ./mmd1_cpp
//         ./mmd1_cpp 3      just the third example

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
    // chain; production's Buckets has no display and so has neither. A DEGREE is one dimensional
    // and so std::uint32_t; the links are std::int32_t, carrying NIL.
    std::uint32_t size() const                 { return mSize; }
    std::int32_t  head(std::uint32_t d) const  { return mHead[d]; }
    bool          empty(std::uint32_t d) const { return mHead[d] == NIL; }
    std::int32_t  next(std::int32_t u) const   { return mNext[u]; }

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

std::vector<std::int32_t> mmd1Neighbors(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u);


// Print a quotient graph: adjacency, incidence, cliques, in the order the
// structure holds them.
void mmd1Show(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
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
void mmd1ShowState(const std::vector<std::uint32_t>& degrees, const Buckets& buckets,
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
std::size_t mmd1Storage(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C) {
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
std::vector<std::int32_t> mmd1Neighbors(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
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
    // the sources were walked in.
    //
    // A REACH TAG, ABOUT VERTEX u, LABELLING reach(u) together with u. Not about any clique: the
    // cliques in I[u] are read here as SOURCES of members, never stamped as ids. Consumed before
    // this function returns, unlike the eliminator's two, which stay live across its whole prune
    // loop; that is why the sweep guard may sit before this call and not before those.
    ++tag;
    std::vector<std::int32_t> neighbors;
    mark[u] = tag;                          // never its own neighbor
    for (std::int32_t v : A[u]) { mark[v] = tag; neighbors.push_back(v); }
    for (std::int32_t c : I[u])
        for (std::int32_t v : C[c])
            if (mark[v] != tag) { mark[v] = tag; neighbors.push_back(v); }
    return neighbors;
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
mmd1Eliminate(AdjacencyGraph& A, IncidenceGraph& I, Cliques& C, std::vector<bool>& eliminated,
             std::vector<std::int32_t>& mark, std::int32_t& tag, std::int32_t pivot) {
    const std::vector<std::int32_t> neighbors = mmd1Neighbors(A, I, C, mark, tag, pivot);
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
        kept.clear();                            // KEPT IS ADJACENCY here: A[u] - C[pivot] - {pivot}
        for (std::int32_t v : A[u]) {
            if (v == pivot) continue;            // the pivot is no longer a variable
            if (mark[v] == pivotCliqueTag) {     // both ends inside the new clique
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
    // the same closed neighborhood, mmd1Neighbors(u) | {u} == mmd1Neighbors(pivot)
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
std::vector<std::int32_t> mmd1MinimumDegree(const AdjacencyGraph& G, std::int32_t delta = 0) {
    const std::uint32_t n = G.size();
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

    Buckets buckets(n);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        buckets.file(degrees[u], u);
    std::uint32_t minDegree = n > 0 ? *std::min_element(degrees.begin(), degrees.end()) : 0;
    std::size_t numBucketProbes = 0;
    std::vector<std::int32_t> evictedMark(n, NIL);       // the iteration u was last evicted in

    // NOT PRODUCTION: display only. The trace is what makes these files teachable and
    // is the whole reason they exist; nothing downstream reads it.
    if (n <= SHOW_THRESHOLD) {
        mmd1Show(A, I, C, degrees, "start: every edge explicit, no clique yet", &eliminated);
        mmd1ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
    }
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
        // and evictedMark is the same idea one level up: it stamps the iteration a
        // vertex was evicted in, so the refresh set is accumulated without a set
        // and without a sort. The stamp is numIterations, which is monotone and
        // already maintained, so there is no separate tag to bump and the array
        // never needs clearing.
        // Clamped: a degree is at most n - 1, so a wider window would walk the
        // bucket array off its end.
        std::size_t batchLimit = minDegree;      // delta > 0 here, so no narrowing
        if (delta > 0)
            batchLimit = std::min(minDegree + delta, static_cast<std::uint32_t>(n) - 1);
        std::vector<std::int32_t> batch;
        std::vector<std::int32_t> evicted;    // first-eviction order, no set and no sort
        while (true) {
            if (buckets.empty(minDegree)) {     // this degree is drained
                if (minDegree >= batchLimit) break;
                ++minDegree;
                ++numBucketProbes;
                continue;
            }
            std::int32_t pivot = buckets.head(minDegree);

            // Sweep the tag back before it can wrap. Two sites in this layer, one
            // before each region that advances the tag, and each placed where nothing
            // in mark is live. This one is INSIDE the batch loop rather than before
            // it, since a batch takes several pivots and each calls the eliminator.
            // Safe between eliminations because the eviction that follows stamps
            // evictedMark and filed, which are separate arrays. Not inside
            // mmd1Eliminate, which holds three stamps live in turn: pivotCliqueTag and
            // absorbedCliquesTag across the prune loop, then the merged set across the
            // C[pivot] compaction. Never observed to fire.
            if (tag > TAG_CEILING) {
                std::fill(mark.begin(), mark.end(), NIL);
                tag = 0;
                ++numTagSweeps;
            }
            auto [neighbors, absorbedCliques, prunedEdges, mergedVertices] =
                mmd1Eliminate(A, I, C, eliminated, mark, tag, pivot);
            ++numEliminations;
            numCliqueEntries += C[pivot].size();
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminatedVertices += 1 + mergedVertices.size();
            for (std::int32_t u : mergedVertices) {   // the pivot now stands for them too
                superMembers[pivot].insert(superMembers[pivot].end(),
                                           superMembers[u].begin(), superMembers[u].end());
                superMembers[u].clear();
            }

            buckets.unfile(degrees[pivot], pivot);    // the pivot has left
            degrees[pivot] = 0;
            for (std::int32_t u : mergedVertices) {   // and so have the merged vertices
                buckets.unfile(degrees[u], u);
                degrees[u] = 0;
            }

            for (std::int32_t u : C[pivot]) {   // EVICT, with a stale degree
                buckets.unfile(degrees[u], u);
                if (evictedMark[u] != static_cast<std::int32_t>(numIterations)) {
                    evictedMark[u] = static_cast<std::int32_t>(numIterations);
                    evicted.push_back(u);       // a marker, so O(1) per eviction
                }
            }

            // A supervariable of size w is w consecutive columns of L. Its
            // external degree is what remains of the clique after the merges,
            // since a merged vertex joins the supervariable instead of
            // neighboring it. The first column holds ext + w - 1 entries below
            // its diagonal, the next ext + w - 2, down to ext, and each column
            // contributes its own diagonal.
            std::uint32_t superSize = superMembers[pivot].size();
            std::uint32_t externalDegree = C[pivot].size();
            // ONE FACTOR WIDENED: a product of two one-dimensional quantities is TWO dimensional, so
            // it is formed in std::size_t. Widening cannot be done after the multiply the way
            // narrowing is done after a subtraction; one operand is enough, the other promoting to
            // meet it. Both factors are bounded by n, so the product reaches n^2.
            nnzL += static_cast<std::size_t>(superSize) * externalDegree
                  + static_cast<std::size_t>(superSize) * (superSize - 1) / 2 + superSize;

            // NOT PRODUCTION: display only, and silent above the threshold. Everything
            // the line needs is built inside the guard, so a run above it formats
            // nothing. The degree printed is neighbors.size(), the reach the eliminator
            // found, which for a PIVOT equals the cached degrees[pivot] it was picked
            // at: a pivot comes off a bucket, so it was not evicted this iteration, so
            // nothing since its last refresh changed a source of its reach.
            if (n <= SHOW_THRESHOLD) {
                const std::uint32_t degree = neighbors.size();
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
            if (delta < 0) break;               // one pivot per iteration, as md5 does
        }

        // ---- one REFRESH, for everything the batch reached -----------------
        std::vector<std::int32_t> refreshedVertices;
        for (std::int32_t u : evicted) if (!eliminated[u]) refreshedVertices.push_back(u);
        // The second site, before the degree update pass. Safe here because the
        // batch's stamps are all spent, and because every mmd1Neighbors call stamps
        // what it reads in the same call.
        if (tag > TAG_CEILING) {
            std::fill(mark.begin(), mark.end(), NIL);
            tag = 0;
            ++numTagSweeps;
        }
        for (std::int32_t u : refreshedVertices) {
            degrees[u] = mmd1Neighbors(A, I, C, mark, tag, u).size();
            buckets.file(degrees[u], u);
        }
        numDegreeUpdates += refreshedVertices.size();
        for (std::int32_t u : refreshedVertices) minDegree = std::min(minDegree, degrees[u]);
        ++numIterations;

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
        // NOT PRODUCTION: display only. The trace is what makes these files teachable and
        // is the whole reason they exist; nothing downstream reads it.
        if (n <= SHOW_THRESHOLD) {
            mmd1Show(A, I, C, degrees, title.str(), &eliminated);
            mmd1ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
        }
    }

    std::vector<std::int32_t> order;
    for (std::int32_t pivot : pivots)
        for (std::int32_t u : superMembers[pivot]) order.push_back(u);
    std::cout << "n = " << n << ", nnz(L) = " << nnzL
              << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "iterations: " << numIterations << "\n";
    std::cout << "eliminations: " << numEliminations << "\n";
    std::cout << "sum of |C[p]|: " << numCliqueEntries << "\n";
    std::cout << "degree computations: " << (numDegreeUpdates + n)
              << ", degree updates: " << numDegreeUpdates
              << ", bucket probes: " << numBucketProbes << "\n";
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
    mmd1MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Grid mode, spelled the same way in every layer and in both twins, so `make test` can diff
    // at a size the seven examples cannot reach. Nothing is filtered: the run is silent above
    // SHOW_THRESHOLD and prints its closing lines as always.
    //
    //   ./mmd1_cpp grid 22
    if (argc > 2 && std::string(argv[1]) == "grid") {
        const int side = std::atoi(argv[2]);
        std::cout << "=== grid " << side << "x" << side << " (n = " << side * side << ") ===\n";
        mmd1MinimumDegree(gridGraph(side));
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
    // belongs to two elements that overlap outside the new one, which needs enough
    // eliminations to have made several elements and enough fill for them to
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
    // which mmd1's merge test declines a genuine supervariable. At the iteration whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test mmd1Neighbors(A, I, C, u) contained in
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

    // All of them by default. To run just one, pass its number: ./mmd1_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
