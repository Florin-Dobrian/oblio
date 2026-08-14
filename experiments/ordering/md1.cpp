// md1.cpp -- minimum degree, iteration 1: the smallest version.
//
// Naive minimum degree, nothing else. Eliminate the vertex of least degree,
// make its neighbors a clique, repeat. The new edges are FILL: the whole point
// of the ordering is to keep them few. Section 5.1 of
// archive/sparse_factorization.md as code. We build on it later.
//
// It names each fill edge as it is created, so the ordering can be seen earning
// (or wasting) its keep, iteration by iteration.
//
// Build:  g++ -std=c++17 -O3 md1.cpp -o md1_cpp  (or: make)
// Run:    ./md1_cpp
//         ./md1_cpp 3      just the third example

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// The adjacency of a vertex is a plain vector, UNSORTED. A set costs O(log d) per
// membership test and per insertion; a sorted vector costs a merge to keep the
// order. Neither is needed: membership comes from a mark array stamped with a
// tag, which answers in one comparison and keeps every pass linear in what it
// touches, and that is what the vendored codes and Oblio's own SymFactorEngine
// do. See the README section on complexity.
//
// Types follow Oblio's rule, since this code is meant to grow into the ordering
// engine, and there are THREE of them. An INDEX names a vertex and is a
// std::int32_t, signed only because NIL has to share a type with the values it
// stands in for. A ONE DIMENSIONAL SIZE is bounded by n and is a std::uint32_t:
// nothing to stand in for, so no sentinel and no sign bit spent. A TWO
// DIMENSIONAL size or position is bounded by nnz and is a std::size_t. A vector
// of indices is std::vector<std::int32_t>; an entity loop therefore has a
// signedness cast where its int32 counter meets a uint32 bound.
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

// ONE GRAPH HERE, and it is a class rather than an alias. `mSize` is the id space, a one
// dimensional size, so `std::uint32_t`, and holding it here is what keeps `n` out of
// `std::size_t` for the whole layer: every vector length in the file is then bounded by it.
// md2 onward add an IncidenceGraph beside this one, which is the point at which having two
// TYPES rather than two variables of one alias starts to matter.
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

// Print a graph: adjacency lists, in the order the structure holds them.
void md1Show(const AdjacencyGraph& A, const std::string& title = "",
             const std::vector<bool>* eliminated = nullptr) {
    const std::uint32_t n = A.size();
    int width = static_cast<int>(std::to_string(n > 0 ? n - 1 : 0).size());  // setw field
    std::vector<std::int32_t> liveVertices;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        if (eliminated == nullptr || !(*eliminated)[u]) liveVertices.push_back(u);
    std::size_t numLiveEdges = 0;
    for (std::int32_t u : liveVertices) numLiveEdges += A[u].size();
    numLiveEdges /= 2;
    if (!title.empty()) std::cout << title << "\n";
    std::ostringstream liveVerticesText;
    if (eliminated == nullptr) liveVerticesText << n;
    else liveVerticesText << liveVertices.size() << " of " << n;
    std::cout << "num live vertices = " << liveVerticesText.str()
              << ", num live edges = " << numLiveEdges
              << ", storage = " << 2 * numLiveEdges << "\n";
    for (std::int32_t u : liveVertices) {
        std::ostringstream adjacencyText;
        bool first = true;
        for (std::int32_t v : A[u]) {
            adjacencyText << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::cout << "  " << std::setw(width) << u << ": {" << adjacencyText.str()
                  << "} degree " << A[u].size() << "\n";
    }
    std::cout << "\n";
}

// What the graph currently costs: one entry per edge endpoint. Compare with md2,
// where the same number falls monotonically. Here fill pushes it back up.
std::size_t md1Storage(const AdjacencyGraph& A) {
    std::size_t total = 0;   // TWO DIMENSIONAL, a count of entries, so it stays wide
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(A.size()); ++u) total += A[u].size();
    return total;
}

// Make the pivot's neighbors a clique, then remove the pivot.
//
// Returns (neighbors, fillEdges): the pivot's adjacency at elimination, which is
// the pattern of its column of L, and the fill edges created among those
// neighbors, pairs that were not already adjacent.
//
// In set terms this is the elimination game itself, three lines:
//
//     for u in A[pivot]:
//         fill(u) = A[pivot] - A[u] - {u}        what was not already there
//         A[u]    = ( A[u] | fill(u) ) - {pivot}
//     A[pivot] = {}
//
// The loop below is that difference without a set: stamp A[u], then walk A[pivot]
// and keep whatever is unstamped. Two passes over vectors rather than a hash per
// element, and A[pivot] is captured before the loop because u is inside it.
std::pair<std::vector<std::int32_t>, std::vector<std::pair<std::int32_t, std::int32_t>>>
md1Eliminate(AdjacencyGraph& A, std::vector<std::int32_t>& mark, std::int32_t& tag,
             std::vector<bool>& eliminated, std::int32_t pivot) {
    const std::vector<std::int32_t> neighbors = A[pivot];
    std::vector<std::pair<std::int32_t, std::int32_t>> fillEdges;

    // Nothing is sorted. For each neighbor u, stamp what u already sees, then walk
    // the clique once and append what is missing. One pass per neighbor, O(d) each,
    // against the O(d log d) a merge would cost to keep an order nobody needs.
    //
    // ONE TAG PER NEIGHBOR, EACH ABOUT ONE NEIGHBOR u AND LABELLING WHAT u ALREADY SEES: A[u],
    // with u and the pivot stamped alongside so they fail the test below and never become fill.
    // A tag is about a VERTEX here, never about a set of vertices shared by several of them,
    // because there is no such set in this layer: the fill is pairwise and each neighbor's
    // missing edges are its own. So the set a stamp belongs to changes every time round the
    // loop, and the eliminator's advance is |A[pivot]|, at most n - 1, which is the whole of
    // this layer's cost in the tag-overflow table.
    //
    // md2 is where that collapses, and the reason is what its tags are ABOUT rather than how
    // many there are. The quotient graph gives the fill a name, the clique, so its eliminator
    // stamps sets that belong to the PIVOT rather than to each neighbor: the pivot's reach, the
    // members of the clique that reach becomes, and the ids of the cliques it absorbs. Three,
    // whatever the pivot's degree, against one per neighbor here.
    for (std::int32_t u : neighbors) {
        ++tag;
        for (std::int32_t v : A[u]) mark[v] = tag;
        mark[u] = tag;                      // never adjacent to itself
        mark[pivot] = tag;                  // the pivot is leaving anyway
        for (std::int32_t v : neighbors) {
            if (mark[v] != tag) {           // A[pivot] - A[u] - {u, pivot}
                mark[v] = tag;
                A[u].push_back(v);
                if (u < v) fillEdges.push_back({u, v});
            }
        }
    }

    std::vector<std::int32_t> kept;
    for (std::int32_t u : neighbors) {       // drop the pivot, compacting in place
        kept.clear();
        for (std::int32_t v : A[u])
            if (v != pivot) kept.push_back(v);
        A[u].swap(kept);
    }
    A[pivot].clear();
    eliminated[pivot] = true;
    return {neighbors, fillEdges};
}

// Eliminate the least-degree vertex each iteration, naming the fill it makes.
std::vector<std::int32_t> md1MinimumDegree(const AdjacencyGraph& G) {
    const std::uint32_t n = G.size();
    std::size_t nnzTrilA = 0;                      // before we mutate it
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    AdjacencyGraph A = G;
    std::vector<std::int32_t> mark(n, NIL);   // scratch for membership, stamped with tag
    std::int32_t tag = 0;
    // Calls to the eliminate procedure, one per pivot. Not the count of vertices
    // removed: a pivot can carry mass-merged vertices out with it, and from mmd1 up
    // an iteration batches several eliminations before one degree update pass. The three
    // counts coincide only where both of those are absent.
    std::uint32_t numEliminations = 0;
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
    std::size_t totalFill = 0;
    std::size_t degreeSum = 0;     // sum of pivot degrees == sum of column counts of L

    // NOT PRODUCTION: display only. The trace is what makes these files teachable and
    // is the whole reason they exist; nothing downstream reads it.
    if (n <= SHOW_THRESHOLD) {
        md1Show(A, "start: every edge explicit, no fill yet", &eliminated);
    }
    for (std::uint32_t iteration = 0; iteration < n; ++iteration) {
        ++numIterations;
        std::int32_t pivot = NIL;
        // The scan asks every live vertex for its degree, so the count is the live
        // count summed over iterations, n(n+1)/2 here since exactly one vertex leaves
        // per iteration. Two things keep it from being comparable with the layers
        // above. There is no initial build to charge for, degrees being computed here
        // and nowhere else, so this starts at 0 where md4 and md5 start at n. And a
        // degree computation is A[u].size() rather than a union over A[u] and the
        // cliques in I[u], so it is the same count of a much cheaper operation.
        for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
            if (eliminated[u]) continue;
            ++numDegreeComputations;
            if (pivot == NIL || A[u].size() < A[pivot].size()) pivot = u;
        }
        // Sweep the tag back before it can wrap. Here because nothing in mark is
        // live between eliminations: every pass inside md1Eliminate stamps what it
        // reads in the same pass, so there is nothing to erase. One elimination
        // advances the tag once per neighbor of the pivot, at most n, which is the
        // room the ceiling has to leave and does. Never observed to fire.
        if (tag > TAG_CEILING) {
            std::fill(mark.begin(), mark.end(), NIL);
            tag = 0;
            ++numTagSweeps;
        }
        auto [neighbors, fillEdges] = md1Eliminate(A, mark, tag, eliminated, pivot);
        ++numEliminations;
        std::uint32_t degree = neighbors.size();
        order.push_back(pivot);
        totalFill += fillEdges.size();
        degreeSum += degree;

        // NOT PRODUCTION: display only, and silent above the threshold. Built INSIDE
        // the guard, so a run above the threshold formats nothing: these are per
        // elimination, and on a grid that is work for a line nobody prints.
        if (n <= SHOW_THRESHOLD) {
            std::ostringstream fillEdgesText;
            if (fillEdges.empty()) {
                fillEdgesText << "none";
            } else {
                bool first = true;
                for (auto [u, v] : fillEdges) {
                    fillEdgesText << (first ? "" : ", ") << u << "-" << v;
                    first = false;
                }
            }
            std::ostringstream title;
            title << "iteration " << iteration << ": eliminate " << pivot << " (degree " << degree
                  << "), fill edges: " << fillEdgesText.str()
                  << ", fill so far: " << totalFill;
            md1Show(A, title.str(), &eliminated);
        }
    }

    // The degree of a pivot at elimination is the count of its column of L, so
    // the degrees already computed give nnz(L) with no extra work (Section 5.1).
    std::size_t nnzL = degreeSum + n;
    std::cout << "n = " << n << ", nnz(L) = " << nnzL
              << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA)
              << " (fill edges counted: " << totalFill << ")\n";
    std::cout << "iterations: " << numIterations << "\n";
    std::cout << "eliminations: " << numEliminations << "\n";
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
    md1MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Grid mode, spelled the same way in every layer and in both twins, so `make test` can diff
    // at a size the seven examples cannot reach. Nothing is filtered: the run is silent above
    // SHOW_THRESHOLD and prints its closing lines as always.
    //
    //   ./md1_cpp grid 22
    if (argc > 2 && std::string(argv[1]) == "grid") {
        const int side = std::atoi(argv[2]);
        std::cout << "=== grid " << side << "x" << side << " (n = " << side * side << ") ===\n";
        md1MinimumDegree(gridGraph(side));
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
        {1, 3},   // 0
        {0, 2},   // 1
        {1, 3},   // 2
        {0, 2},   // 3
    };
    AdjacencyGraph graph2 = {
        {1, 2},      // 0
        {0, 3},      // 1
        {0, 4},      // 2
        {1, 4, 5},   // 3
        {2, 3, 5},   // 4
        {3, 4},      // 5
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

    // All of them by default. To run just one, pass its number: ./md1_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (std::size_t number = 1; number <= examples.size(); ++number) {
        if (selected != 0 && number != static_cast<std::size_t>(selected)) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
