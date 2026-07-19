// minimum_degree.cpp -- the minimum degree ordering, written for legibility.
//
// Layer 1 of the legible MMD study. This is the NAIVE version: it materializes
// fill, the one thing a real code refuses to do (Section 5.1 of
// archive/sparse_factorization.md). That is deliberate. It is the vertex
// elimination game transcribed line for line, so it is the ground truth that
// the quotient graph (layer 2) and everything after it optimize *without
// changing what is computed*. Slow, up to O(n^3), and completely readable.
//
// Scope: correct, not equivalent. It produces a valid minimum degree ordering,
// but does not reproduce the vendored genmmd's exact tie-breaking. The vendored
// src/Mmd.cpp is the oracle for a later, equivalent version; this one is for
// understanding the mechanism.
//
// Build:  g++ -std=c++17 -O3 minimum_degree.cpp -o minimum_degree  (or: make)
// Run:    ./minimum_degree

#include <cstddef>
#include <iostream>
#include <set>
#include <vector>

// The elimination graph: one neighbor set per vertex. A std::set keeps each
// neighborhood sorted and duplicate-free, which is all we need, and it reads
// plainly. Slow next to a flat array, but flexibility first.
using Graph = std::vector<std::set<int>>;

// Build the graph from a symmetric off-diagonal CSC pattern. We take CSC so the
// input matches the vendored orderer (and the rest of the tree); we convert to
// the friendly Graph immediately and never look at CSC again.
Graph buildGraph(int n, const std::vector<std::size_t>& colPtr,
                 const std::vector<int>& rowIdx) {
    Graph graph(n);
    for (int column = 0; column < n; ++column)
        for (std::size_t position = colPtr[column]; position < colPtr[column + 1]; ++position) {
            int row = rowIdx[position];
            if (row != column) {
                graph[column].insert(row);
                graph[row].insert(column);  // keep it symmetric
            }
        }
    return graph;
}

// Eliminate one vertex: make its neighbors a clique, then detach it. Making the
// neighbors mutually adjacent is where fill is created, and here we actually
// store it. This function *is* the definition of a pivot's effect on the graph.
void eliminate(Graph& graph, int pivot) {
    const std::set<int> neighbors = graph[pivot];  // copy: we mutate the graph below
    for (int u : neighbors)
        for (int w : neighbors)
            if (u != w)
                graph[u].insert(w);  // the fill edge, materialized
    for (int u : neighbors)
        graph[u].erase(pivot);
    graph[pivot].clear();
}

// The whole algorithm: repeatedly eliminate the live vertex of least degree.
// Returns the elimination order, order[k] being the vertex eliminated k-th,
// which read as a sequence is the new-to-old permutation.
std::vector<int> minimumDegree(int n, const std::vector<std::size_t>& colPtr,
                               const std::vector<int>& rowIdx) {
    Graph graph = buildGraph(n, colPtr, rowIdx);
    std::vector<bool> eliminated(n, false);
    std::vector<int> order;
    order.reserve(n);

    for (int step = 0; step < n; ++step) {
        // PICK the live vertex of least degree. Linear scan: the naive part.
        int pivot = -1;
        std::size_t leastDegree = 0;
        for (int vertex = 0; vertex < n; ++vertex) {
            if (eliminated[vertex]) continue;
            std::size_t degree = graph[vertex].size();
            if (pivot == -1 || degree < leastDegree) {
                pivot = vertex;
                leastDegree = degree;
            }
        }

        eliminate(graph, pivot);
        eliminated[pivot] = true;
        order.push_back(pivot);
    }
    return order;
}

// Count the fill an ordering produces, for illustration. Replays the
// elimination on a fresh copy of the graph and tallies each new edge once.
int countFill(Graph graph, const std::vector<int>& order) {
    int fill = 0;
    for (int pivot : order) {
        const std::set<int> neighbors = graph[pivot];
        for (int u : neighbors)
            for (int w : neighbors)
                if (u < w && graph[u].count(w) == 0) {
                    graph[u].insert(w);
                    graph[w].insert(u);
                    ++fill;
                }
        for (int u : neighbors) graph[u].erase(pivot);
        graph[pivot].clear();
    }
    return fill;
}

int main() {
    // Arrowhead: vertex 0 is a hub joined to 1..n-1; every other vertex touches
    // only the hub. This is Example 1 of the doc: eliminate the hub first and
    // its neighbors become a full clique (much fill); eliminate it last and
    // there is none. Minimum degree finds the good order on its own.
    const int n = 6;
    std::vector<std::vector<int>> adjacency(n);
    for (int leaf = 1; leaf < n; ++leaf) {
        adjacency[0].push_back(leaf);
        adjacency[leaf].push_back(0);
    }

    std::vector<std::size_t> colPtr(n + 1, 0);
    std::vector<int> rowIdx;
    for (int vertex = 0; vertex < n; ++vertex) {
        colPtr[vertex + 1] = colPtr[vertex] + adjacency[vertex].size();
        for (int neighbor : adjacency[vertex]) rowIdx.push_back(neighbor);
    }

    std::vector<int> mdOrder = minimumDegree(n, colPtr, rowIdx);
    std::vector<int> naturalOrder(n);
    for (int vertex = 0; vertex < n; ++vertex) naturalOrder[vertex] = vertex;  // hub first

    std::cout << "minimum degree order:      ";
    for (int vertex : mdOrder) std::cout << vertex << " ";
    std::cout << "  fill = " << countFill(buildGraph(n, colPtr, rowIdx), mdOrder) << "\n";

    std::cout << "natural order (hub first): ";
    for (int vertex : naturalOrder) std::cout << vertex << " ";
    std::cout << "  fill = " << countFill(buildGraph(n, colPtr, rowIdx), naturalOrder) << "\n";

    return 0;
}
