#pragma once

// graphs.h -- the test graphs, defined once for every driver in this folder.
//
// WHY THIS EXISTS. The seven examples were written out twice, in vendored.cpp and in
// production.cpp, and a third copy was about to go into amdorder.cpp. Three copies of a literal
// is not a style complaint here: the drivers answer different questions ABOUT THE SAME GRAPHS,
// and `make test` compares their outputs line for line, so a graph that drifts in one copy makes
// two drivers disagree for a reason that is not the code. That failure would read exactly like a
// defect in an ordering. One definition removes the possibility rather than the symptom.
//
// It also means a graph added for one driver is available to the others, which is the half that
// matters going forward: the acceptance test wants shapes the seven examples cannot reach, and
// nothing should have to be transcribed to get them there.
//
// THE REPRESENTATION IS THE PROTOTYPES'. `Graph` is an adjacency list per vertex, ascending, with
// no diagonal, which is what every layer in this folder takes and what the Python twins print.
// vendored.cpp used to hold its copy as `std::vector<std::set<int>>`; a set iterates ascending
// and these vectors are written ascending, so the two produce identical CSC and the change is a
// representation only. Each driver keeps its OWN conversion to CSC, and that is deliberate: they
// need different ones. Ours takes a full-symmetric pattern with the diagonal present, which is
// what a SparseMatrix holds; the vendored routines take the off-diagonal pattern alone.
//
// WHAT EACH EXAMPLE IS FOR is documented in the prototype that introduced it and in the README's
// "The test graphs" section, which is the place to read rather than here.

#include <cstdint>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace OrderingExperiment {

// One adjacency list per vertex, ascending, no diagonal.
using Graph = std::vector<std::vector<std::int32_t>>;

// The seven examples, in the order the drivers number them: `./production_cpp amd3 3` is the
// third of these.
inline const std::vector<std::pair<std::string, Graph>>& exampleGraphs() {
    static const std::vector<std::pair<std::string, Graph>> examples = {
        {"graph1", Graph{
            {1, 3}, {0, 2}, {1, 3}, {0, 2},
        }},
        {"graph2", Graph{
            {1, 2}, {0, 3}, {0, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4},
        }},
        {"graph3", Graph{
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
        }},
        {"graph4", Graph{
            {2, 3, 4, 7},     // 0
            {3, 4, 6, 7},     // 1
            {0, 3, 5},        // 2
            {0, 1, 2, 6, 7},  // 3
            {0, 1, 5},        // 4
            {2, 4, 6},        // 5
            {1, 3, 5},        // 6
            {0, 1, 3},        // 7
        }},
        {"graph5", Graph{
            {3, 4},           // 0
            {2, 4},           // 1
            {1},              // 2
            {0},              // 3
            {0, 1},           // 4
        }},
        {"graph6", Graph{
            {2, 3, 4},        // 0
            {3},              // 1
            {0, 3, 4, 5},     // 2
            {0, 1, 2, 4},     // 3
            {0, 2, 3},        // 4
            {2},              // 5
        }},
        {"graph7", Graph{
            {1, 2, 4},        // 0
            {0, 4},           // 1
            {0, 3, 4},        // 2
            {2, 4},           // 3
            {0, 1, 2, 3},     // 4
        }},
    };
    return examples;
}

// A square five-point grid, the same one every prototype builds.
//
// A grid is not an eighth example: nothing about it illustrates a mechanism and its trace is
// unreadable. It is here because the seven examples are at most twelve vertices, which is too
// small to fire a mechanism that needs real structure, and two defects in amd2 lived behind
// exactly that gap. See the amd2 section of the README.
//
// The push order is ascending by construction, `u - side` then `u - 1` then `u + 1` then
// `u + side`, and that is load bearing rather than incidental. A column's row indices must be
// ascending: it is the CSC precondition SparseMatrix states, and it is also a TIE-BREAK INPUT,
// since the order within a column decides the content order of C[pivot] and so which of several
// equal-degree vertices a bucket hands over first. A builder that emitted an unsorted column
// would make two drivers diverge on the first elimination for a reason that is not the code.
inline Graph gridGraph(int side) {
    const int n = side * side;
    Graph graph(n);
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

// A cubic seven-point grid. A DIFFERENT SHAPE and not a larger size: it fills faster, makes
// larger cliques, and mass-eliminates far more often than a square grid does, which is why a
// defect can live behind a 2D-only check at every size and surface here at 16 a side.
//
// The push order is `z-1, y-1, x-1, x+1, y+1, x+side`-ish and is NOT ascending, so the lists are
// sorted before returning. This is the one that caught us out: written in the natural
// `x-1, x+1, y-1, y+1, z-1, z+1` order the far corner comes out `62, 59, 47`, descending, and a
// driver that consumes the pattern as given then parts company with one that sorts it on the
// FIRST elimination. See the note on ascending order above gridGraph; the 2D builder is
// ascending by construction and so never exposed it.
inline Graph grid3dGraph(int side) {
    const int n = side * side * side;
    Graph graph(n);
    const auto id = [side](int x, int y, int z) { return (z * side + y) * side + x; };
    for (int z = 0; z < side; ++z)
        for (int y = 0; y < side; ++y)
            for (int x = 0; x < side; ++x) {
                const int u = id(x, y, z);
                if (z > 0)        graph[u].push_back(id(x, y, z - 1));
                if (y > 0)        graph[u].push_back(id(x, y - 1, z));
                if (x > 0)        graph[u].push_back(id(x - 1, y, z));
                if (x + 1 < side) graph[u].push_back(id(x + 1, y, z));
                if (y + 1 < side) graph[u].push_back(id(x, y + 1, z));
                if (z + 1 < side) graph[u].push_back(id(x, y, z + 1));
            }
    return graph;   // ascending by construction: z-1 < y-1 < x-1 < x+1 < y+1 < z+1
}

// An irregular pattern AT SIZE, which neither grid family gives. Each vertex draws `degree`
// partners and every edge is symmetrized, so the mean degree comes out near twice that and the
// spread is wide, which is what exercises a mechanism whose cost depends on how uneven the
// structure is.
//
// Deterministic across platforms, and that is a requirement rather than a nicety: this is an
// equality test, so a pattern that differed between two machines would make the check
// unreproducible. `std::mt19937`'s output sequence is fixed by the standard while the
// distribution templates are implementation defined, so the partner is taken from the raw
// engine and no `<random>` distribution appears. The same discipline, for the same reason, as
// the tier-1 matrices in notes/TESTING_SPECIFICATION.md.
inline Graph randomGraph(int n, int degree, unsigned seed) {
    std::mt19937 engine(seed);
    std::vector<std::set<std::int32_t>> neighbors(n);
    for (int u = 0; u < n; ++u)
        for (int k = 0; k < degree; ++k) {
            const int v = static_cast<int>(engine() % static_cast<unsigned>(n));
            if (v == u) continue;                       // no self loop: it is not an edge
            neighbors[u].insert(v);
            neighbors[v].insert(u);                     // symmetric, as a graph must be
        }
    Graph graph(n);
    for (int u = 0; u < n; ++u)
        graph[u].assign(neighbors[u].begin(), neighbors[u].end());
    return graph;   // a std::set iterates ascending, which is the order required above
}

} // namespace OrderingExperiment
