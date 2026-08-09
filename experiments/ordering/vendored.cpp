// vendored.cpp -- run the vendored MMD and AMD on our test graphs.
//
// Nothing here is ours to change: vendored/vendored_mmd.cpp and
// vendored/vendored_amd.cpp are copies of src/Mmd.cpp and src/Amd.cpp, and this
// file only feeds them the same seven graphs the prototypes use and prints their
// permutations in our format. It exists so that "no feature missing" is a diff
// rather than a judgement.
//
// The lowercase names are deliberate. Oblio capitalizes source files, but macOS
// is case-insensitive by default, so Amd.cpp and amd.cpp are one path: dropping
// the vendored copy next to the prototype silently overwrote it once already.
//
// Build:  g++ -std=c++17 -O3 vendored.cpp vendored/vendored_mmd.cpp \
//             vendored/vendored_amd.cpp -o vendored_cpp
// Run:    ./vendored_cpp
//         ./vendored_cpp 3      just the third example

#include "graphs.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using OrderingExperiment::Graph;

void mmd_order(int n, const int colPtr[], const int rowIdx[], int perm[], int invp[]);
extern "C" int amd_order(int32_t n, const int32_t Ap[], const int32_t Ai[],
                         int32_t P[], double Control[], double Info[]);

// Off-diagonal-only symmetric CSC, which is what both routines expect. The graphs used to be
// held here as `std::vector<std::set<int>>` and are now the shared `Graph`; a set iterates
// ascending and those lists are ascending, so this produces the identical pattern.
static void toCsc(const Graph& graph,
                  std::vector<int>& colPtr, std::vector<int>& rowIdx) {
    int n = static_cast<int>(graph.size());
    colPtr.assign(n + 1, 0);
    rowIdx.clear();
    for (int j = 0; j < n; ++j) {
        colPtr[j] = static_cast<int>(rowIdx.size());
        for (int i : graph[j]) rowIdx.push_back(i);
    }
    colPtr[n] = static_cast<int>(rowIdx.size());
}

static void printOrder(const std::string& label, const std::vector<int>& order) {
    std::cout << label << ": [";
    for (std::size_t k = 0; k < order.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << order[k];
    std::cout << "]\n";
}

static void run(const std::string& name, const Graph& graph) {
    int n = static_cast<int>(graph.size());
    std::cout << "=== " << name << " ===\n";
    std::vector<int> colPtr, rowIdx;
    toCsc(graph, colPtr, rowIdx);

    std::vector<int> perm(n), invp(n);
    mmd_order(n, colPtr.data(), rowIdx.data(), perm.data(), invp.data());
    printOrder("mmd order", perm);

    std::vector<int32_t> ap(colPtr.begin(), colPtr.end()), ai(rowIdx.begin(), rowIdx.end());
    std::vector<int32_t> p(n);
    double control[5] = {0, 0, 0, 0, 0}, info[20] = {0};
    int status = amd_order(n, ap.data(), ai.data(), p.data(), nullptr, info);
    std::vector<int> amdOrder(p.begin(), p.end());
    printOrder("amd order", amdOrder);
    if (status != 0) std::cout << "amd status " << status << "\n";
    (void)control;
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // The seven examples come from graphs.h, which vendored.cpp, production.cpp and raworder.cpp
    // all read, so a graph cannot drift between the drivers whose outputs `make test` compares.
    const auto& examples = OrderingExperiment::exampleGraphs();
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
