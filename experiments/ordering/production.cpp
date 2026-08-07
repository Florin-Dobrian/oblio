// production.cpp -- run Oblio's own ordering drivers on the prototypes' graphs.
//
// The mirror image of vendored.cpp. That one runs the routines we are reading against and
// prints their permutations in our format; this one runs the routines we pulled out of these
// prototypes, so that "production still says what the prototype says" is a diff rather than a
// judgement. `make test` compares the order lines of each pair.
//
// It links ../../src directly rather than a copy, deliberately and unlike vendored.cpp: a copy
// is right for code that is not ours to edit and wrong for code we are actively changing, since
// the whole point here is to notice when the two come apart.
//
// The graphs are given as adjacency lists, as the prototypes give them, and then built into a
// full-symmetric CSC with the diagonal present, which is what a SparseMatrix holds and what
// OrderEngine hands down. So the CSC overload is what runs, and buildGraph's dropping of the
// diagonal is under the same check as everything else.
//
// Build:  make production   (the sources are the two under ../../src, plus this file)
// Run:    ./production_cpp mmd1
//         ./production_cpp amd1
//         ./production_cpp mmd1 3      just the third example

#include "oblio/Amd1.h"
#include "oblio/Amd2.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/Mmd1.h"
#include "oblio/Mmd2.h"
#include "oblio/Mmd3.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"
#include "oblio/SymFactor.h"
#include "oblio/SymFactorEngine.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// The prototypes' own shape for a graph. The production side no longer has a type for it:
// QuotientGraph is built from a matrix pattern directly, its adjacency being one flat array
// rather than a list per vertex, so the conversion below is the whole of the difference.
using Graph = std::vector<std::vector<std::int32_t>>;

// Adjacency lists to full-symmetric CSC with a structurally present diagonal, which is the
// input Oblio's own matrices satisfy by construction.
static void toCsc(const Graph& graph,
                  std::vector<std::size_t>& colPtr, std::vector<std::int32_t>& rowIdx) {
    const std::int32_t size = static_cast<std::int32_t>(graph.size());
    colPtr.assign(graph.size() + 1, 0);
    rowIdx.clear();
    for (std::int32_t aj = 0; aj < size; ++aj) {
        bool diagonalWritten = false;
        for (std::int32_t ai : graph[aj]) {
            if (!diagonalWritten && ai > aj) { rowIdx.push_back(aj); diagonalWritten = true; }
            rowIdx.push_back(ai);
        }
        if (!diagonalWritten) rowIdx.push_back(aj);
        colPtr[aj + 1] = rowIdx.size();
    }
}

static void printOrder(const std::vector<std::int32_t>& order) {
    std::cout << "order: [";
    for (std::size_t k = 0; k < order.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << order[k];
    std::cout << "]\n";
}

// The same square five-point grid the prototypes build, transcribed from gridGraph in each of
// them. It has to match exactly, since the check below compares permutations and any difference
// in the numbering would make them differ for a reason that is not the code.
//
// A grid is not an eighth example: nothing about it illustrates a mechanism. It is here because
// the seven examples are at most twelve vertices, which is too small to fire a mechanism that
// needs real structure, and two defects in amd2 lived behind exactly that gap. See the amd2
// section of the README.
static Graph gridGraph(int side) {
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

static void run(const std::string& layer, const std::string& name, const Graph& graph) {
    std::cout << "=== " << name << " ===\n";
    std::vector<std::size_t>  colPtr;
    std::vector<std::int32_t> rowIdx;
    toCsc(graph, colPtr, rowIdx);
    if (layer == "mmd1") printOrder(Oblio::orderMmd1(colPtr, rowIdx));
    if (layer == "amd1") printOrder(Oblio::orderAmd1(colPtr, rowIdx));
    if (layer == "mmd2") printOrder(Oblio::orderMmd2(colPtr, rowIdx));
    if (layer == "mmd3") printOrder(Oblio::orderMmd3(colPtr, rowIdx));
    if (layer == "amd2") printOrder(Oblio::orderAmd2(colPtr, rowIdx));
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // The prototypes' seven graphs, copied as vendored.cpp copies them. What each one is for is
    // documented in the prototype that introduced it and in the README.
    Graph graph1 = {
        {1, 3}, {0, 2}, {1, 3}, {0, 2},
    };
    Graph graph2 = {
        {1, 2}, {0, 3}, {0, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4},
    };
    Graph graph3 = {
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
    Graph graph4 = {
        {2, 3, 4, 7},     // 0
        {3, 4, 6, 7},     // 1
        {0, 3, 5},        // 2
        {0, 1, 2, 6, 7},  // 3
        {0, 1, 5},        // 4
        {2, 4, 6},        // 5
        {1, 3, 5},        // 6
        {0, 1, 3},        // 7
    };
    Graph graph5 = {
        {3, 4},           // 0
        {2, 4},           // 1
        {1},              // 2
        {0},              // 3
        {0, 1},           // 4
    };
    Graph graph6 = {
        {2, 3, 4},        // 0
        {3},              // 1
        {0, 3, 4, 5},     // 2
        {0, 1, 2, 4},     // 3
        {0, 2, 3},        // 4
        {2},              // 5
    };
    Graph graph7 = {
        {1, 2, 4},        // 0
        {0, 4},           // 1
        {0, 3, 4},        // 2
        {2, 4},           // 3
        {0, 1, 2, 3},     // 4
    };

    std::vector<std::pair<std::string, Graph>> examples = {
        {"graph1", graph1}, {"graph2", graph2},
        {"graph3", graph3}, {"graph4", graph4},
        {"graph5", graph5}, {"graph6", graph6},
        {"graph7", graph7},
    };

    const std::string layer = (argc > 1) ? argv[1] : "mmd1";

    // Grid mode, matching the prototypes' own:  ./production_cpp amd2 grid 20
    if (argc > 3 && std::string(argv[2]) == "grid") {
        const int side = std::atoi(argv[3]);
        run(layer, "grid " + std::to_string(side) + "x" + std::to_string(side), gridGraph(side));
        return 0;
    }

    const int selected      = (argc > 2) ? std::atoi(argv[2]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(layer, examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
