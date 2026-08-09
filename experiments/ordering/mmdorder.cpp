// mmdorder.cpp -- production Mmd3 against genmmd's elimination order.
//
// THE ACCEPTANCE TEST FOR THE MMD ALIGNMENT, and one of a pair with amdorder.cpp. `Mmd3` adds no
// mechanism to `Mmd2`: it exists to return genmmd's permutation, so that a comparison against the
// vendored routine is an EQUALITY TEST rather than a fill number somebody has to judge. This is
// what makes that claim checkable.
//
// THE TWO ARE NAMED FOR THE BRANCH, NOT FOR THE MECHANISM, and the mechanisms are not alike:
//
//   genmmd emits the order directly.  `mmd_order` returns `perm`, the order it eliminates in, and
//                                     there is no postorder anywhere in the routine. So the
//                                     vendored output vector IS the object to compare: no hook, no
//                                     generated copy, no anchors to assert, no Control array. This
//                                     file is the whole of the check.
//
//   AMD gets it through a HOOK.       `amd_order` returns a vector AMD_postorder has relabeled, so
//                                     the raw order has to be reconstructed upstream of it by a
//                                     generated, hooked copy of the vendored source. That is the
//                                     whole of hook_amd.py, and none of it is needed here.
//
// The difference is genmmd's rather than ours, and it is written down here and in amdorder.cpp
// rather than encoded in the target names, which say only which branch each checks.
//
// AND THERE IS NO KNOB TO SET WRONG HERE, which is the other half of the asymmetry. amdorder.cpp
// has to switch off dense-row removal, the one mechanism amd3 lacks, and getting that wrong cost a
// day: the threshold is computed by a double-to-Int conversion that overflows for a large enough
// alpha, so the same line meant different things on two platforms. genmmd has no dense handling
// and no Control array, so nothing is passed and nothing can be misread. Its one tunable is the
// `maxint` ceiling our wrapper supplies for the marker sweep, and `experiments/ordering/REPORT.md`
// records that sweep firing zero times at every size we run, so it is not a hidden term in this
// comparison either.
//
// `make aligned` runs this and amdorder.cpp together. See experiments/ordering/README.md, "mmd3,
// and the alignment ledger".

#include "oblio/Mmd3.h"

#include "graphs.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// The vendored routine, from private/Mmd.cpp. Linked rather than copied or generated: it emits
// what we need without alteration, which is exactly why this file has no generator beside it.
void mmd_order(int n, const int colPtr[], const int rowIdx[], int perm[], int invp[]);

// The graphs come from graphs.h, shared with vendored.cpp, production.cpp and amdorder.cpp, so a
// graph added for one driver is available to the others and none of them can drift.
using OrderingExperiment::Graph;
using OrderingExperiment::grid3dGraph;
using OrderingExperiment::gridGraph;
using OrderingExperiment::randomGraph;

// A graph, two ways: ours takes the diagonal, genmmd takes a pattern without it. Both are the
// same conversions amdorder.cpp makes, and they stay in each driver rather than moving into
// graphs.h, since what a driver feeds its routine is the driver's own business.
static void toCsc(const Graph& graph, std::vector<std::size_t>& colPtr,
                  std::vector<std::int32_t>& rowIdx) {
    colPtr.assign(graph.size() + 1, 0);
    rowIdx.clear();
    for (std::size_t j = 0; j < graph.size(); ++j) {
        for (std::int32_t i : graph[j]) if (i < static_cast<std::int32_t>(j)) rowIdx.push_back(i);
        rowIdx.push_back(static_cast<std::int32_t>(j));
        for (std::int32_t i : graph[j]) if (i > static_cast<std::int32_t>(j)) rowIdx.push_back(i);
        colPtr[j + 1] = rowIdx.size();
    }
}

static void toCscNoDiagonal(const Graph& graph,
                            std::vector<int>& colPtr, std::vector<int>& rowIdx) {
    colPtr.assign(graph.size() + 1, 0);
    rowIdx.clear();
    for (std::size_t j = 0; j < graph.size(); ++j) {
        for (std::int32_t i : graph[j]) rowIdx.push_back(i);
        colPtr[j + 1] = static_cast<int>(rowIdx.size());
    }
}

static bool check(const std::string& name, const Graph& graph) {
    std::vector<std::size_t>  colPtr;
    std::vector<std::int32_t> rowIdx;
    toCsc(graph, colPtr, rowIdx);
    const std::vector<std::int32_t> ours = Oblio::orderMmd3(colPtr, rowIdx);

    std::vector<int> ap, ai;
    toCscNoDiagonal(graph, ap, ai);
    const int n = static_cast<int>(graph.size());
    std::vector<int> perm(n), invp(n);
    mmd_order(n, ap.data(), ai.data(), perm.data(), invp.data());

    if (ours.size() != perm.size()) {
        std::printf("  %-22s SIZE MISMATCH: ours %zu, genmmd %zu\n",
                    name.c_str(), ours.size(), perm.size());
        return false;
    }
    for (std::size_t k = 0; k < ours.size(); ++k) {
        if (ours[k] != perm[k]) {
            std::printf("  %-22s DIFFER at %zu: ours %d, genmmd %d\n",
                        name.c_str(), k, ours[k], perm[k]);
            return false;
        }
    }
    std::printf("  %-22s n = %6d   elimination order matches\n", name.c_str(), n);
    return true;
}

int main(int argc, char** argv) {
    std::printf("production Mmd3 against genmmd's elimination order\n");
    std::printf("  (the permutation itself: genmmd emits the order it eliminates in and\n");
    std::printf("   does no postorder, so nothing has to be reconstructed)\n\n");

    std::vector<std::pair<std::string, Graph>> cases;

    // Sizes on the command line run 2D grids only, for bisecting a failure.
    if (argc > 1) {
        for (int a = 1; a < argc; ++a) {
            const int side = std::atoi(argv[a]);
            cases.emplace_back("grid " + std::to_string(side) + "x" + std::to_string(side),
                               gridGraph(side));
        }
    } else {
        // The same FOUR SHAPES amdorder.cpp runs, and for the same reason: widening one grid from
        // 4 to 140 exercises scale and never mechanism, so a defect that needs a structure grids
        // do not produce passes every size of it. On the amd branch that was not hypothetical, a
        // 3D grid at 16 a side finding a defect eleven 2D sizes had missed. Nothing on this branch
        // is known to need the wider shapes, which is not a reason to run less: what a check
        // cannot reach it cannot check, and this one costs a fraction of a second.
        for (const auto& example : OrderingExperiment::exampleGraphs())
            cases.emplace_back(example.first, example.second);

        for (int side : {4, 5, 8, 10, 16, 20, 32, 40, 64, 100, 140})
            cases.emplace_back("grid " + std::to_string(side) + "x" + std::to_string(side),
                               gridGraph(side));

        for (int side : {2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24})
            cases.emplace_back("grid3d " + std::to_string(side) + "^3", grid3dGraph(side));

        for (int degree : {3, 6, 12})
            for (unsigned seed : {1u, 2u, 3u})
                cases.emplace_back("random d" + std::to_string(degree) + " s" +
                                       std::to_string(seed),
                                   randomGraph(2000, degree, seed));
    }

    int failed = 0;
    for (auto& c : cases) if (!check(c.first, c.second)) ++failed;

    std::printf("\n%zu cases, %d failed\n", cases.size(), failed);
    return failed == 0 ? 0 : 1;
}
