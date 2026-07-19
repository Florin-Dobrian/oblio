// md1.cpp -- minimum degree, step 1: the smallest version.
//
// Naive minimum degree, nothing else. Eliminate the vertex of least degree,
// make its neighbors a clique, repeat. The new edges are FILL: the whole point
// of the ordering is to keep them few. Section 5.1 of
// archive/sparse_factorization.md as code. We build on it later.
//
// It names each fill edge as it is created, so the ordering can be seen earning
// (or wasting) its keep, step by step.
//
// Build:  g++ -std=c++17 -O3 md1.cpp -o md1_cpp  (or: make)
// Run:    ./md1_cpp

#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

using Graph = std::vector<std::set<int>>;

// Make the pivot's neighbors a clique, then remove the pivot. Returns the fill
// edges created: pairs of neighbors that were not already adjacent.
std::vector<std::pair<int, int>> eliminate(Graph& graph, int pivot) {
    const std::set<int> neighbors = graph[pivot];
    std::vector<std::pair<int, int>> fill;
    for (int u : neighbors)
        for (int w : neighbors)
            if (u < w && graph[u].count(w) == 0) {   // a genuinely new edge
                graph[u].insert(w);
                graph[w].insert(u);
                fill.push_back({u, w});
            }
    for (int u : neighbors)
        graph[u].erase(pivot);
    graph[pivot].clear();
    return fill;
}

// What the graph currently costs: one entry per edge endpoint. Compare with md3,
// where the same number falls monotonically. Here fill pushes it back up.
std::size_t storage(const Graph& graph, const std::vector<bool>& eliminated) {
    std::size_t total = 0;
    for (int v = 0; v < static_cast<int>(graph.size()); ++v)
        if (!eliminated[v]) total += graph[v].size();
    return total;
}

// Print the vertices still live and their neighbors.
void showState(const Graph& graph, const std::vector<bool>& eliminated) {
    std::cout << "         graph:   ";
    for (int v = 0; v < static_cast<int>(graph.size()); ++v) {
        if (eliminated[v]) continue;
        std::cout << v << ":{";
        bool first = true;
        for (int u : graph[v]) { std::cout << (first ? "" : ",") << u; first = false; }
        std::cout << "} ";
    }
    std::cout << "  storage " << storage(graph, eliminated) << "\n";
}

// Eliminate the least-degree vertex each step, naming the fill it makes.
std::vector<int> minimumDegree(Graph& graph) {
    int n = static_cast<int>(graph.size());
    std::vector<bool> eliminated(n, false);
    std::vector<int> order;
    int totalFill = 0;
    int degreeSum = 0;             // sum of pivot degrees == column counts of L
    std::size_t nnzTrilA = 0;      // measured before we mutate the graph
    for (int v = 0; v < n; ++v) nnzTrilA += graph[v].size();
    nnzTrilA = nnzTrilA / 2 + n;

    std::cout << "start: every edge explicit, no fill yet\n";
    showState(graph, eliminated);
    for (int step = 0; step < n; ++step) {
        int pivot = -1;
        for (int v = 0; v < n; ++v) {
            if (eliminated[v]) continue;
            if (pivot == -1 || graph[v].size() < graph[pivot].size())
                pivot = v;
        }
        int degree = static_cast<int>(graph[pivot].size());
        std::vector<std::pair<int, int>> fill = eliminate(graph, pivot);
        eliminated[pivot] = true;
        order.push_back(pivot);
        totalFill += static_cast<int>(fill.size());
        degreeSum += degree;

        std::cout << "step " << step << ": eliminate " << pivot
                  << " (degree " << degree << ")\n";
        std::cout << "         fill: ";
        if (fill.empty()) {
            std::cout << "none";
        } else {
            bool first = true;
            for (auto [u, w] : fill) { std::cout << (first ? "" : ", ") << u << "-" << w; first = false; }
        }
        std::cout << "   (fill so far " << totalFill << ")\n";
        showState(graph, eliminated);
    }

    // The degree of a pivot at elimination is the count of its column of L, so
    // the degrees already computed give nnz(L) with no extra work (Section 5.1).
    std::size_t nnzL = static_cast<std::size_t>(degreeSum) + n;
    std::cout << "nnz(L) = " << nnzL << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA)
              << "   (edges counted: " << totalFill << ")\n";
    return order;
}

void run(const std::string& name, Graph graph) {
    std::cout << "=== " << name << " ===\n";
    std::vector<int> order = minimumDegree(graph);
    std::cout << "order:";
    for (int v : order) std::cout << " " << v;
    std::cout << "\n\n";
}

int main() {
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
    Graph graph1 = {
        {1, 3},   // 0
        {0, 2},   // 1
        {1, 3},   // 2
        {0, 2},   // 3
    };
    Graph graph2 = {
        {1, 2},      // 0
        {0, 3},      // 1
        {0, 4},      // 2
        {1, 4, 5},   // 3
        {2, 3, 5},   // 4
        {3, 4},      // 5
    };

    std::vector<std::set<int>> graph3 = {
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
    std::vector<std::set<int>> graph4 = {
        {2, 3, 4, 7},     // 0
        {3, 4, 6, 7},     // 1
        {0, 3, 5},        // 2
        {0, 1, 2, 6, 7},  // 3
        {0, 1, 5},        // 4
        {2, 4, 6},        // 5
        {1, 3, 5},        // 6
        {0, 1, 3},        // 7
    };

    run("graph1", graph1);
    run("graph2", graph2);
    run("graph3", graph3);
    run("graph4", graph4);
    return 0;
}
