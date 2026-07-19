// md3.cpp -- minimum degree, step 3: the quotient graph.
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
// So an elimination adds nothing and removes something. Each A[i] only ever
// shrinks, which is why this representation never needs more room than the
// original graph. Section 5.2 of archive/sparse_factorization.md.
//
// A live variable i is stored as A[i], its remaining explicit variable neighbors,
// and C[i], the cliques it belongs to. Its true neighborhood is the union of the
// two, formed only when asked.
//
// Naming: the literature calls these objects ELEMENTS and writes E_i for C[i].
// They are cliques; we name them for what they are.
//
// The order and the per-step degrees match md1 exactly: same algorithm, cheaper
// storage. What this layer does NOT yet fix is that the degree is still a full
// union every time it is asked; a cheap degree is a later layer.
//
// Build:  g++ -std=c++17 -O3 md3.cpp -o md3_cpp  (or: make)
// Run:    ./md3_cpp

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using Cliques = std::map<int, std::set<int>>;

// The true neighbors of live variable i: its explicit neighbors A[i] together
// with the members of every element it belongs to, minus itself.
std::set<int> reachable(const std::vector<std::set<int>>& A,
                        const std::vector<std::set<int>>& C,
                        const Cliques& cliques, int i) {
    std::set<int> neighbors = A[i];
    for (int c : C[i])
        neighbors.insert(cliques.at(c).begin(), cliques.at(c).end());
    neighbors.erase(i);
    return neighbors;
}

// What the quotient graph currently costs: explicit endpoints plus element
// members. Watch it fall. The naive graph's edge count only rises.
std::size_t storage(const std::vector<std::set<int>>& A, const Cliques& cliques) {
    std::size_t total = 0;
    for (const auto& a : A) total += a.size();
    for (const auto& [e, members] : cliques) total += members.size();
    return total;
}

// Turn the pivot into an element. Fills neighbors (the pivot's reachable set,
// stored once as an element; in that role the doc calls it the element's pattern
// L_pivot, also the nonzero pattern of column pivot of L), absorbed (cliques the
// pivot belonged to), and pruned (explicit edges the new element now implies).
void eliminate(std::vector<std::set<int>>& A, std::vector<std::set<int>>& C,
               Cliques& cliques, std::vector<bool>& eliminated, int pivot,
               std::set<int>& neighbors, std::set<int>& absorbed,
               std::vector<std::pair<int, int>>& pruned) {
    neighbors = reachable(A, C, cliques, pivot);
    absorbed = C[pivot];
    for (int c : absorbed) cliques.erase(c);
    cliques[pivot] = neighbors;      // becomes L_pivot: the element's pattern

    pruned.clear();
    for (int i : neighbors) {
        std::set<int> redundant;         // both ends inside the new clique
        std::set_intersection(A[i].begin(), A[i].end(), neighbors.begin(), neighbors.end(),
                              std::inserter(redundant, redundant.begin()));
        for (int j : redundant)
            if (i < j) pruned.push_back({i, j});
        for (int j : redundant) A[i].erase(j);   // implicit now: delete the copy
        A[i].erase(pivot);                        // the pivot is no longer a variable
        for (int e : absorbed) C[i].erase(e);     // its absorbed cliques are gone
        C[i].insert(pivot);                       // i belongs to the new element
    }

    A[pivot].clear();
    C[pivot].clear();
    eliminated[pivot] = true;
}

std::string braces(const std::set<int>& s) {
    std::string out = "{";
    bool first = true;
    for (int x : s) { out += (first ? "" : ", ") + std::to_string(x); first = false; }
    return out + "}";
}

// A, C and cliques on their own lines; then the three results of the eliminate
// call, in the order it returns them. `hasResult` false for the start-of-run dump.
void showState(const std::vector<std::set<int>>& A, const std::vector<std::set<int>>& C,
               const Cliques& cliques, bool hasResult,
               const std::set<int>& neighbors, const std::set<int>& absorbed,
               const std::vector<std::pair<int, int>>& pruned) {
    int n = static_cast<int>(A.size());

    std::cout << "         A       = {";
    for (int v = 0; v < n; ++v)
        std::cout << (v ? ", " : "") << v << ": " << braces(A[v]);
    std::cout << "}\n";

    std::cout << "         C       = {";
    for (int v = 0; v < n; ++v) {
        std::cout << (v ? ", " : "") << v << ": [";
        bool f = true;
        for (int x : C[v]) { std::cout << (f ? "" : ", ") << "c" << x; f = false; }
        std::cout << "]";
    }
    std::cout << "}\n";

    std::cout << "         cliques = {";
    bool f = true;
    for (const auto& [x, members] : cliques) {
        std::cout << (f ? "" : ", ") << "c" << x << ": " << braces(members);
        f = false;
    }
    std::cout << "}   storage " << storage(A, cliques) << "\n";

    if (!hasResult) {
        std::cout << "         neighbors = none, absorbed = none, pruned = none"
                     "   (nothing eliminated yet)\n";
        return;
    }
    std::cout << "         neighbors = " << braces(neighbors) << ", absorbed = ";
    if (absorbed.empty()) std::cout << "none";
    else { bool g = true; for (int x : absorbed) { std::cout << (g ? "" : ", ") << "c" << x; g = false; } }
    std::cout << ", pruned = ";
    if (pruned.empty()) std::cout << "none";
    else { bool g = true; for (auto [u, w] : pruned) { std::cout << (g ? "" : ", ") << u << "-" << w; g = false; } }
    std::cout << "\n";
}

std::vector<int> minimumDegree(const std::vector<std::set<int>>& graph) {
    int n = static_cast<int>(graph.size());
    std::size_t nnzTrilA = 0;               // measured before we mutate
    for (const auto& row : graph) nnzTrilA += row.size();
    nnzTrilA = nnzTrilA / 2 + graph.size();
    std::vector<std::set<int>> A = graph;   // explicit variable neighbors
    std::vector<std::set<int>> C(n);        // cliques each variable belongs to
    Cliques cliques;                      // element id -> its live members
    std::vector<bool> eliminated(n, false);
    std::vector<int> order;
    std::size_t degreeSum = 0;              // sum of pivot degrees == column counts of L

    std::cout << "start: no cliques yet, every neighbor explicit\n";
    showState(A, C, cliques, false, {}, {}, {});
    for (int step = 0; step < n; ++step) {
        int pivot = -1;
        std::size_t leastDegree = 0;
        for (int v = 0; v < n; ++v) {
            if (eliminated[v]) continue;
            std::size_t d = reachable(A, C, cliques, v).size();
            if (pivot == -1 || d < leastDegree) { pivot = v; leastDegree = d; }
        }
        std::size_t degree = reachable(A, C, cliques, pivot).size();

        std::set<int> neighbors, absorbed;
        std::vector<std::pair<int, int>> pruned;
        eliminate(A, C, cliques, eliminated, pivot, neighbors, absorbed, pruned);
        order.push_back(pivot);
        degreeSum += degree;

        std::cout << "step " << step << ": eliminate " << pivot
                  << " (degree " << degree << ")\n";
        showState(A, C, cliques, true, neighbors, absorbed, pruned);
    }

    // The degree of a pivot at elimination is the count of its column of L, so
    // the degrees already computed give nnz(L) with no extra work (Section 5.1).
    std::size_t nnzL = degreeSum + n;
    std::cout << "nnz(L) = " << nnzL << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    return order;
}

void run(const std::string& name, const std::vector<std::set<int>>& graph) {
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
    std::vector<std::set<int>> graph1 = {
        {1, 3}, {0, 2}, {1, 3}, {0, 2},
    };
    std::vector<std::set<int>> graph2 = {
        {1, 2}, {0, 3}, {0, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4},
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
