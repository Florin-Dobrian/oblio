// vendored.cpp -- run the vendored MMD and AMD on our test graphs.
//
// Nothing here is ours to change: vendored/Mmd.cpp and vendored/Amd.cpp are
// copies of src/Mmd.cpp and src/Amd.cpp, and this file only feeds them the same
// seven graphs the prototypes use and prints their permutations in our format.
// It exists so that "no feature missing" is a diff rather than a judgement.
//
// Build:  g++ -std=c++17 -O3 vendored.cpp vendored/Mmd.cpp vendored/Amd.cpp -o vendored_cpp
// Run:    ./vendored_cpp
//         ./vendored_cpp 3      just the third example

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

void mmd_order(int n, const int colPtr[], const int rowIdx[], int perm[], int invp[]);
extern "C" int amd_order(int32_t n, const int32_t Ap[], const int32_t Ai[],
                         int32_t P[], double Control[], double Info[]);

// Off-diagonal-only symmetric CSC, which is what both routines expect.
static void toCsc(const std::vector<std::set<int>>& graph,
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

static void run(const std::string& name, const std::vector<std::set<int>>& graph) {
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
    std::vector<std::set<int>> graph1 = {
        {1, 3},
        {0, 2},
        {1, 3},
        {0, 2},
    };
    std::vector<std::set<int>> graph2 = {
        {1, 2},
        {0, 3},
        {0, 4},
        {1, 4, 5},
        {2, 3, 5},
        {3, 4},
    };
    std::vector<std::set<int>> graph3 = {
        {1, 3, 8},
        {0, 2, 6, 8},
        {1, 3, 5},
        {0, 2, 4},
        {3, 5},
        {2, 4, 6, 9},
        {1, 5, 7, 10},
        {6, 8},
        {0, 1, 7, 9},
        {5, 8, 10},
        {6, 9, 11},
        {10},
    };
    std::vector<std::set<int>> graph4 = {
        {2, 3, 4, 7},
        {3, 4, 6, 7},
        {0, 3, 5},
        {0, 1, 2, 6, 7},
        {0, 1, 5},
        {2, 4, 6},
        {1, 3, 5},
        {0, 1, 3},
    };
    std::vector<std::set<int>> graph5 = {
        {3, 4},
        {2, 4},
        {1},
        {0},
        {0, 1},
    };
    std::vector<std::set<int>> graph6 = {
        {2, 3, 4},
        {3},
        {0, 3, 4, 5},
        {0, 1, 2, 4},
        {0, 2, 3},
        {2},
    };
    std::vector<std::set<int>> graph7 = {
        {1, 2, 4},
        {0, 4},
        {0, 3, 4},
        {2, 4},
        {0, 1, 2, 3},
    };
    std::vector<std::pair<std::string, std::vector<std::set<int>>>> examples = {
        {"graph1", graph1}, {"graph2", graph2},
        {"graph3", graph3}, {"graph4", graph4},
        {"graph5", graph5}, {"graph6", graph6},
        {"graph7", graph7},
    };
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
