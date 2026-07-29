// md2.cpp -- minimum degree, step 2: the quotient graph.
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
// So an elimination adds nothing and removes something. Each A[u] only ever
// shrinks, which is why this representation never needs more room than the
// original graph. Section 5.3 of archive/sparse_factorization.md.
//
// A live vertex u is stored as A[u], its remaining explicit vertex neighbors, and
// I[u], the ids of the cliques that contain u; C[c] holds the members of clique c,
// so an incidence is stored twice, once from each side, just as an edge is. The
// true neighborhood of u is the union of the two, formed only when asked.
//
// Naming: the literature calls the cliques ELEMENTS and writes A_i and E_i for
// what we call A[u] and I[u]. They are cliques; we name them for what they are.
//
// The order and the per-step degrees match md1 exactly: same algorithm, cheaper
// storage. What this layer does NOT yet fix is that the degree is still a full
// union every time it is asked; a cheap degree is a later layer.
//
// Build:  g++ -std=c++17 -O3 md2.cpp -o md2_cpp  (or: make)
// Run:    ./md2_cpp

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using Graph = std::vector<std::set<int>>;
using Cliques = std::map<int, std::set<int>>;

// I[u] cliques that contain u
// C[c] vertices that c contains

std::set<int> md2Neighbors(const Graph& A, const Graph& I, const Cliques& C, int u);

// Print a quotient graph: adjacency sets, incidence sets, cliques.
void md2Show(const Graph& A, const Graph& I, const Cliques& C,
             const std::string& title = "",
             const std::vector<bool>* eliminated = nullptr) {
    int n = static_cast<int>(A.size());
    int width = static_cast<int>(std::to_string(std::max(n - 1, 0)).size());
    std::vector<int> aliveVertices;
    for (int u = 0; u < n; ++u)
        if (eliminated == nullptr || !(*eliminated)[u]) aliveVertices.push_back(u);
    std::size_t numAliveEdges = 0;
    for (int u : aliveVertices) numAliveEdges += A[u].size();
    numAliveEdges /= 2;
    std::size_t numAliveIncidences = 0;
    for (int u : aliveVertices) numAliveIncidences += I[u].size();
    std::size_t numAliveCliques = C.size();
    if (!title.empty()) std::cout << title << "\n";
    std::ostringstream aliveVerticesText;
    if (eliminated == nullptr) aliveVerticesText << n;
    else aliveVerticesText << aliveVertices.size() << " of " << n;
    std::cout << "num alive vertices = " << aliveVerticesText.str()
              << ", num alive edges = " << numAliveEdges
              << ", num alive cliques = " << numAliveCliques
              << ", storage = " << 2 * numAliveEdges << " + " << 2 * numAliveIncidences
              << " = " << 2 * (numAliveEdges + numAliveIncidences) << "\n";
    for (int u : aliveVertices) {
        std::ostringstream adjacencyText;
        bool first = true;
        for (int v : A[u]) {
            adjacencyText << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::ostringstream incidenceText;
        first = true;
        for (int c : I[u]) {
            incidenceText << (first ? "" : " ") << "c" << c;
            first = false;
        }
        std::size_t degree = md2Neighbors(A, I, C, u).size();
        std::cout << "  " << std::setw(width) << u << ": {" << adjacencyText.str()
                  << "} {" << incidenceText.str() << "} degree " << degree << "\n";
    }
    for (const auto& [c, cliqueMembers] : C) {
        std::ostringstream cliqueMembersText;
        bool first = true;
        for (int u : cliqueMembers) {
            cliqueMembersText << (first ? "" : " ") << std::setw(width) << u;
            first = false;
        }
        std::cout << "  c" << c << ": {" << cliqueMembersText.str() << "}\n";
    }
    std::cout << "\n";
}

// Entries actually stored. Each edge costs two, one per endpoint in A. Each
// incidence costs two as well, the clique id in I and the member in C. Watch
// the total fall monotonically; the naive graph's only rises.
std::size_t md2Storage(const Graph& A, const Graph& I, const Cliques& C) {
    std::size_t total = 0;
    for (const std::set<int>& adjacency : A) total += adjacency.size();
    for (const std::set<int>& incidence : I) total += incidence.size();
    for (const auto& [c, cliqueMembers] : C) { (void)c; total += cliqueMembers.size(); }
    return total;
}

// The neighbors of live vertex u: its explicit adjacency A[u] together with the
// members of every clique that contains u, minus u itself, which the cliques
// always carry. This is George and Liu's reachable set, and it is what the
// elimination graph would hold explicitly.
std::set<int> md2Neighbors(const Graph& A, const Graph& I, const Cliques& C, int u) {
    std::set<int> neighbors = A[u];
    for (int c : I[u]) {
        const std::set<int>& cliqueMembers = C.at(c);
        neighbors.insert(cliqueMembers.begin(), cliqueMembers.end());
    }
    neighbors.erase(u);
    return neighbors;
}

// Turn the pivot into a clique.
//
// Returns (neighbors, absorbedCliques, prunedEdges): the pivot's neighbor set,
// which becomes the clique and the pattern of its column of L; the cliques that
// the new one swallows; and the explicit edges the new clique makes redundant.
// The last two are reported for display; only neighbors is used by the caller.
std::tuple<std::set<int>, std::set<int>, std::vector<std::pair<int, int>>> md2Eliminate(
        Graph& A, Graph& I, Cliques& C, std::vector<bool>& eliminated, int pivot) {
    const std::set<int> neighbors = md2Neighbors(A, I, C, pivot);
    const std::set<int> absorbedCliques = I[pivot];
    for (int c : absorbedCliques)
        C.erase(c);
    C[pivot] = neighbors;           // becomes L_pivot, the column pattern

    std::vector<std::pair<int, int>> prunedEdges;
    for (int u : neighbors) {
        std::set<int> redundant;    // both ends inside the new clique
        for (int v : A[u])
            if (neighbors.count(v) != 0) redundant.insert(v);
        for (int v : redundant)
            if (u < v) prunedEdges.push_back({u, v});
        for (int v : redundant)
            A[u].erase(v);          // implicit now: delete the explicit copy
        A[u].erase(pivot);          // the pivot is no longer a variable
        for (int c : absorbedCliques)
            I[u].erase(c);          // its absorbed cliques are gone
        I[u].insert(pivot);         // u joins the new clique, whose id is the pivot
    }

    A[pivot].clear();
    I[pivot].clear();
    eliminated[pivot] = true;
    return {neighbors, absorbedCliques, prunedEdges};
}

// Same heuristic as md1, on the quotient graph. No fill is ever stored.
std::vector<int> md2MinimumDegree(const Graph& G) {
    int n = static_cast<int>(G.size());
    std::size_t nnzTrilA = 0;
    for (int u = 0; u < n; ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    Graph A = G;                            // explicit variable neighbors
    Graph I(n);                             // cliques each variable belongs to
    Cliques C;                              // clique id -> member set
    std::vector<bool> eliminated(n, false);
    std::vector<int> order;
    int degreeSum = 0;

    md2Show(A, I, C, "start: every edge explicit, no clique yet", &eliminated);
    for (int step = 0; step < n; ++step) {
        int pivot = -1;
        std::size_t best = 0;
        for (int u = 0; u < n; ++u) {
            if (eliminated[u]) continue;
            std::size_t degree = md2Neighbors(A, I, C, u).size();
            if (pivot == -1 || degree < best) { pivot = u; best = degree; }
        }
        auto [neighbors, absorbedCliques, prunedEdges] =
            md2Eliminate(A, I, C, eliminated, pivot);
        int degree = static_cast<int>(neighbors.size());
        order.push_back(pivot);
        degreeSum += degree;

        std::ostringstream absorbedCliquesText;
        if (absorbedCliques.empty()) {
            absorbedCliquesText << "none";
        } else {
            bool first = true;
            for (int c : absorbedCliques) {
                absorbedCliquesText << (first ? "" : ", ") << "c" << c;
                first = false;
            }
        }
        std::ostringstream prunedEdgesText;
        if (prunedEdges.empty()) {
            prunedEdgesText << "none";
        } else {
            bool first = true;
            for (auto [u, v] : prunedEdges) {
                prunedEdgesText << (first ? "" : ", ") << u << "-" << v;
                first = false;
            }
        }
        std::ostringstream title;
        title << "step " << step << ": eliminate " << pivot << " (degree " << degree
              << "), absorbed cliques: " << absorbedCliquesText.str()
              << ", pruned edges: " << prunedEdgesText.str();
        md2Show(A, I, C, title.str(), &eliminated);
    }

    std::size_t nnzL = static_cast<std::size_t>(degreeSum) + n;
    std::cout << "nnz(L) = " << nnzL << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "order: [";
    for (std::size_t k = 0; k < order.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << order[k];
    std::cout << "]\n";
    return order;
}

void run(const std::string& name, const Graph& G) {
    std::cout << "=== " << name << " ===\n";
    md2MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
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

    // graph5, five vertices and four edges, two paths joined at 4: 2-1-4-0-3.
    // Small and fill free, and here for one reason: it is the smallest graph on
    // which md3's merge test declines a genuine supervariable. At the step whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test md3Neighbors(A, I, C, u) contained in
    // C[pivot] would merge it. See the README section on mass elimination.
    //
    //   edges: 0-3 0-4 1-2 1-4
    std::vector<std::set<int>> graph5 = {
        {3, 4},           // 0
        {2, 4},           // 1
        {1},              // 2
        {0},              // 3
        {0, 1},           // 4
    };

    // graph6, six vertices and eight edges. Here because one small graph carries
    // three things at once. Its supervariable {0, 4} is a supernode but NOT a
    // fundamental one: the elimination forest is 2 -> 1 -> 4 and 3 -> 0 -> 4, so
    // 4 already has 1 as a child when 0 merges into it. The merge happens at step
    // 2 of 5, so the run continues afterwards and the selection degree, 3 over
    // {2, 3, 4}, differs from the external degree, 2 over {2, 3}, with the
    // difference being the weight that merged. And superMembers ends with a hole
    // in the middle, slot 4 empty between two used ones, while no pivot equals
    // its own step number. See the README sections on mass elimination and on
    // external degree.
    //
    //   edges: 0-2 0-3 0-4 1-3 2-3 2-4 2-5 3-4
    std::vector<std::set<int>> graph6 = {
        {2, 3, 4},        // 0
        {3},              // 1
        {0, 3, 4, 5},     // 2
        {0, 1, 2, 4},     // 3
        {0, 2, 3},        // 4
        {2},              // 5
    };

    // graph7, five vertices and six edges. The pairwise case: at the step whose
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
    std::vector<std::set<int>> graph7 = {
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

    // All of them by default. To run just one, pass its number: ./md2_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
