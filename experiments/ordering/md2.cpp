// md2.cpp -- minimum degree, iteration 2: the quotient graph.
//
// Same ordering as md1, computed WITHOUT ever storing fill. When a vertex is
// eliminated it becomes a CLIQUE on the vertices it would have joined. A clique is
// fully described by its vertex list, so every edge inside it is implicit, and
// that cuts twice:
//
//   - the fill edges are never added, and
//   - the edges ALREADY present between two members are now redundant, so they
//     are pruned from the explicit adjacency.
//
// So an elimination adds nothing and removes something. Each A[u] only ever
// shrinks, which is why this representation never needs more room than the
// original graph. Section 5.3 of archive/sparse_factorization.md.
//
// A live vertex u is stored as A[u], its remaining explicit vertex neighbors, and
// I[u], the ids of the cliques that contain u; C[c] holds the members of clique c,
// so an incidence is stored twice, once from each side, just as an edge is. The
// true neighborhood of u is the union of the two, formed only when asked.
//
// Naming: the literature calls the cliques ELEMENTS and writes A_i and E_i for
// what we call A[u] and I[u]. They are cliques; we name them for what they are.
//
// The order and the per-iteration degrees match md1 exactly: same algorithm, cheaper
// storage. What this layer does NOT yet fix is that the degree is still a full
// union every time it is asked; a cheap degree is a later layer.
//
// Build:  g++ -std=c++17 -O3 md2.cpp -o md2_cpp  (or: make)
// Run:    ./md2_cpp
//         ./md2_cpp 3      just the third example

#include <algorithm>
#include <cstdlib>
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
// QuotientGraphFlat owns both roles in one class and will get the same treatment iteratively.
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

// I[u] cliques that contain u
// C[c] vertices that c contains

std::vector<std::int32_t> md2Neighbors(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u);

// Print a quotient graph: adjacency, incidence, cliques, in the order the
// structure holds them.
void md2Show(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
             std::vector<std::int32_t>& mark, std::int32_t& tag,
             const std::string& title = "",
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
        std::uint32_t degree = md2Neighbors(A, I, C, mark, tag, u).size();
        std::cout << "  " << std::setw(width) << u << ": {" << adjacencyText.str()
                  << "} {" << incidenceText.str() << "} degree " << degree << "\n";
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

// Entries actually stored. Each edge costs two, one per endpoint in A. Each
// incidence costs two as well, the clique id in I and the member in C. Watch
// the total fall monotonically; the naive graph's only rises.
std::size_t md2Storage(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C) {
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
std::vector<std::int32_t> md2Neighbors(const AdjacencyGraph& A, const IncidenceGraph& I, const Cliques& C,
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

// Turn the pivot into a clique.
//
// Returns (neighbors, absorbedCliques, prunedEdges): the pivot's neighbor set,
// which becomes the clique and the pattern of its column of L; the cliques that
// the new one swallows; and the explicit edges the new clique makes redundant.
// The last two are reported for display; only neighbors is used by the caller.//
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
std::tuple<std::vector<std::int32_t>, std::vector<std::int32_t>,
           std::vector<std::pair<std::int32_t, std::int32_t>>>
md2Eliminate(AdjacencyGraph& A, IncidenceGraph& I, Cliques& C, std::vector<bool>& eliminated,
             std::vector<std::int32_t>& mark, std::int32_t& tag, std::int32_t pivot) {
    const std::vector<std::int32_t> neighbors = md2Neighbors(A, I, C, mark, tag, pivot);
    const std::vector<std::int32_t> absorbedCliques = I[pivot];
    for (std::int32_t c : absorbedCliques)
        C.erase(c);
    C.create(pivot, neighbors);     // becomes the column pattern of the pivot

    // TWO TAGS, AND BOTH ARE ABOUT CLIQUES. They differ in WHICH SIDE of a clique they name,
    // which is the quotient graph's two-sided representation showing up in the mark array:
    //
    //     pivotCliqueTag       about the PIVOT'S clique, labels its MEMBERS,  so stamps VERTICES
    //     absorbedCliquesTag   about the ABSORBED cliques, labels their IDS,  so stamps CLIQUE IDS
    //
    // Each side is what one loop below needs: pruning A[u] asks whether a VERTEX is in C[pivot],
    // pruning I[u] asks whether a CLIQUE ID is one of the absorbed. Each set is built in one pass
    // and then queried for free, and both loops are compactions in place.
    //
    // ONE TAG WOULD ALSO WORK, and it is worth knowing why two are used. The two tags stamp
    // different sides and the two loops query different sides, so they cannot meet: A[u] holds
    // only vertices, I[u] only clique ids. With a single value each test would still name its own
    // set. Two make that true of the tag alone rather than of the lists as well, and the cost of
    // the second is one advance per elimination, the difference between 3 and 2 in the
    // tag-overflow table's eliminate column. Parked until md3, whose third stamped set, the
    // merged vertices, is on the MEMBER side alongside C[pivot] and so is the first case where
    // two tags would share a side.
    ++tag;
    const std::int32_t pivotCliqueTag = tag;
    for (std::int32_t v : neighbors) mark[v] = pivotCliqueTag;
    ++tag;
    const std::int32_t absorbedCliquesTag = tag;
    for (std::int32_t c : absorbedCliques) mark[c] = absorbedCliquesTag;

    std::vector<std::pair<std::int32_t, std::int32_t>> prunedEdges;

    // ONE SCRATCH BUFFER FOR BOTH COMPACTIONS, and for every neighbor. Each use is a filter into
    // it followed by a swap, so after the swap it holds the list that was just replaced; clear()
    // then empties it while keeping that capacity, and the next fill reuses the allocation. Two
    // named buffers would read no better and would allocate twice. The two uses are labeled
    // below, since the same name means a different list four lines apart.
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

    A[pivot].clear();
    I[pivot].clear();
    eliminated[pivot] = true;
    return {neighbors, absorbedCliques, prunedEdges};
}

// Same heuristic as md1, on the quotient graph. No fill is ever stored.
std::vector<std::int32_t> md2MinimumDegree(const AdjacencyGraph& G) {
    const std::uint32_t n = G.size();
    std::size_t nnzTrilA = 0;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    AdjacencyGraph A = G;                   // explicit vertex neighbors
    IncidenceGraph I(n);                    // cliques each vertex belongs to
    Cliques        C(n);                    // clique id -> member list
    std::vector<std::int32_t> mark(n, NIL);    // scratch for membership, with tag
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
    // Passes of the outer loop, each one a batch of eliminations followed by one
    // degree update pass. Here the batch is always a single elimination, so this
    // equals numEliminations; from mmd1 up the two come apart.
    std::uint32_t numIterations = 0;
    std::size_t numDegreeComputations = 0;
    // Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    // is here as the witness that the guard is inert rather than as a statistic.
    std::size_t numTagSweeps = 0;
    std::vector<bool> eliminated(n, false);
    std::vector<std::int32_t> order;
    std::size_t degreeSum = 0;

    // NOT PRODUCTION: display only. The trace is what makes these files teachable and
    // is the whole reason they exist; nothing downstream reads it.
    if (n <= SHOW_THRESHOLD) {
        md2Show(A, I, C, mark, tag, "start: every edge explicit, no clique yet", &eliminated);
    }
    for (std::uint32_t iteration = 0; iteration < n; ++iteration) {
        ++numIterations;
        // Sweep the tag back before it can wrap. Two sites in this layer, one before
        // each region that advances the tag, and each placed where nothing in mark is
        // live. Here the region is the pivot search, which calls md2Neighbors once per
        // live vertex; every call stamps what it reads in the same call, so there is
        // nothing to erase. Never observed to fire.
        if (tag > TAG_CEILING) {
            std::fill(mark.begin(), mark.end(), NIL);
            tag = 0;
            ++numTagSweeps;
        }
        std::int32_t pivot = NIL;
        std::uint32_t minDegree = 0;   // a degree, so one dimensional
        for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
            if (eliminated[u]) continue;
            ++numDegreeComputations;
            std::uint32_t degree = md2Neighbors(A, I, C, mark, tag, u).size();
            if (pivot == NIL || degree < minDegree) { pivot = u; minDegree = degree; }
        }
        // The second site. Not inside md2Eliminate: it holds pivotCliqueTag and
        // absorbedCliquesTag live across the whole prune loop, so a sweep in there would
        // erase marks about to be read.
        if (tag > TAG_CEILING) {
            std::fill(mark.begin(), mark.end(), NIL);
            tag = 0;
            ++numTagSweeps;
        }
        auto [neighbors, absorbedCliques, prunedEdges] =
            md2Eliminate(A, I, C, eliminated, mark, tag, pivot);
        ++numEliminations;
        numCliqueEntries += C[pivot].size();
        std::uint32_t degree = neighbors.size();
        order.push_back(pivot);
        degreeSum += degree;

        // NOT PRODUCTION: display only, and silent above the threshold. Built INSIDE
        // the guard, so a run above the threshold formats nothing: these are per
        // elimination, and on a grid that is work for a line nobody prints.
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
            std::ostringstream title;
            title << "iteration " << iteration << ": eliminate " << pivot << " (degree " << degree
                  << "), absorbed cliques: " << absorbedCliquesText.str()
                  << ", pruned edges: " << prunedEdgesText.str();
            md2Show(A, I, C, mark, tag, title.str(), &eliminated);
        }
    }

    std::size_t nnzL = degreeSum + n;
    std::cout << "n = " << n << ", nnz(L) = " << nnzL
              << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "iterations: " << numIterations << "\n";
    std::cout << "eliminations: " << numEliminations << "\n";
    std::cout << "sum of |C[p]|: " << numCliqueEntries << "\n";
    std::cout << "degree computations: " << numDegreeComputations << "\n";
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
    md2MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Grid mode, spelled the same way in every layer and in both twins, so `make test` can diff
    // at a size the seven examples cannot reach. Nothing is filtered: the run is silent above
    // SHOW_THRESHOLD and prints its closing lines as always.
    //
    //   ./md2_cpp grid 22
    if (argc > 2 && std::string(argv[1]) == "grid") {
        const int side = std::atoi(argv[2]);
        std::cout << "=== grid " << side << "x" << side << " (n = " << side * side << ") ===\n";
        md2MinimumDegree(gridGraph(side));
        return 0;
    }

    // Two examples.
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
    // which md3's merge test declines a genuine supervariable. At the iteration whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test md3Neighbors(A, I, C, u) contained in
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
    // difference being the weight that merged. And superMembers ends with a hole
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

    // All of them by default. To run just one, pass its number: ./md2_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (std::size_t number = 1; number <= examples.size(); ++number) {
        if (selected != 0 && number != static_cast<std::size_t>(selected)) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
