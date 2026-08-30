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
// Run:    ./production_cpp mmd3
//         ./production_cpp amd3
//         ./production_cpp mmd3 3      just the third example

#include "oblio/AmdFlat.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/MmdFlat.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"
#include "oblio/SymFactor.h"
#include "oblio/SymFactorEngine.h"

#include "graphs.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// The prototypes' own shape for a graph, from graphs.h. The production side has no type for it:
// QuotientGraphFlat is built from a matrix pattern directly, its adjacency being one flat array
// rather than a list per vertex, so the conversion below is the whole of the difference.
using OrderingExperiment::Graph;
using OrderingExperiment::gridGraph;

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

static void run(const std::string& layer, const std::string& name, const Graph& graph) {
    std::cout << "=== " << name << " ===\n";
    std::vector<std::size_t>  colPtr;
    std::vector<std::int32_t> rowIdx;
    toCsc(graph, colPtr, rowIdx);
    if (layer == "mmd3") printOrder(Oblio::orderMmdFlat(colPtr, rowIdx).order());
    if (layer == "amd3") printOrder(Oblio::orderAmdFlat(colPtr, rowIdx).order());
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // The seven examples come from graphs.h, shared with vendored.cpp and amdorder.cpp.
    const auto& examples = OrderingExperiment::exampleGraphs();

    const std::string layer = (argc > 1) ? argv[1] : "mmd3";

    // Grid mode, matching the prototypes' own:  ./production_cpp amd3 grid 20
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
