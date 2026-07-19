// md5.cpp -- minimum degree, step 5: maintained degrees.
//
// Every version so far recomputed a reachable set for EVERY live vertex at EVERY
// step, just to find the smallest, then threw all but one away. On a 3D grid that
// is roughly ten times the necessary work, and the ratio grows with n.
// Section 5.7 of archive/sparse_factorization.md.
//
// The waste is easy to see once stated: eliminating a pivot can only change the
// degrees of the vertices it REACHED. Every other vertex has the same A, the same
// cliques and the same weights as before, so its degree is still whatever it was.
//
// So we keep a degrees[] array and refresh only the reached set. The picker then
// scans cached integers instead of building set unions. This is the first half of
// what both MMD and AMD do before their ideas diverge; md6 adds the second half,
// degree buckets, which removes the scan itself.
//
// Why refreshing the clique alone is enough, in three parts:
//
//   - PRUNING and clique membership change only for members of the new clique
//   - ABSORPTION deletes cliques the pivot belonged to, and every member of such
//     a clique is reachable from the pivot, hence in the new clique
//   - MERGING removes a vertex i, but i merged only because everything it could
//     see lay inside the new clique, so nobody outside sees i disappear
//
// Build:  g++ -std=c++17 -O3 md5.cpp -o md5_cpp  (or: make)
// Run:    ./md5_cpp

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
// with the members of every clique it belongs to, minus itself.
std::set<int> reachable(const std::vector<std::set<int>>& A,
                        const std::vector<std::set<int>>& C,
                        const Cliques& cliques, int i) {
    std::set<int> neighbors = A[i];
    for (int c : C[i])
        neighbors.insert(cliques.at(c).begin(), cliques.at(c).end());
    neighbors.erase(i);
    return neighbors;
}

// External degree, WEIGHTED. Called only on refresh now, not on every pick.
int degreeOf(const std::vector<std::set<int>>& A, const std::vector<std::set<int>>& C,
             const Cliques& cliques, const std::vector<int>& weight, int i) {
    int total = 0;
    for (int j : reachable(A, C, cliques, i)) total += weight[j];
    return total;
}

std::size_t storage(const std::vector<std::set<int>>& A, const Cliques& cliques) {
    std::size_t total = 0;
    for (const auto& a : A) total += a.size();
    for (const auto& [c, members] : cliques) total += members.size();
    return total;
}

// Turn the pivot into a clique, then absorb every member it makes
// indistinguishable.
void eliminate(std::vector<std::set<int>>& A, std::vector<std::set<int>>& C,
               Cliques& cliques, std::vector<int>& weight,
               std::vector<bool>& eliminated, int pivot,
               std::set<int>& neighbors, std::set<int>& absorbed,
               std::vector<std::pair<int, int>>& pruned, std::vector<int>& merged) {
    neighbors = reachable(A, C, cliques, pivot);
    absorbed = C[pivot];
    for (int c : absorbed) cliques.erase(c);
    cliques[pivot] = neighbors;

    pruned.clear();
    for (int i : neighbors) {
        std::set<int> redundant;
        std::set_intersection(A[i].begin(), A[i].end(), neighbors.begin(), neighbors.end(),
                              std::inserter(redundant, redundant.begin()));
        for (int j : redundant)
            if (i < j) pruned.push_back({i, j});
        for (int j : redundant) A[i].erase(j);
        A[i].erase(pivot);
        for (int c : absorbed) C[i].erase(c);
        C[i].insert(pivot);
    }

    // MASS ELIMINATION. i is indistinguishable from the pivot when everything it
    // still sees lies inside the new clique: nothing explicit left, and no other
    // clique to reach through. Eliminating it next would create no fill, so we
    // merge it into the pivot now and let the weight remember how many there are.
    merged.clear();
    for (int i : neighbors) {
        if (A[i].empty() && C[i].size() == 1 && *C[i].begin() == pivot) {
            weight[pivot] += weight[i];
            weight[i] = 0;
            C[i].clear();
            eliminated[i] = true;
            merged.push_back(i);
        }
    }
    for (int i : merged) {                 // a merged vertex is no longer a vertex
        for (auto& [c, members] : cliques) members.erase(i);
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

void showState(const std::vector<std::set<int>>& A, const std::vector<std::set<int>>& C,
               const Cliques& cliques, const std::vector<int>& weight,
               const std::vector<int>& degrees, const std::vector<bool>& eliminated,
               bool hasResult,
               const std::set<int>& neighbors, const std::set<int>& absorbed,
               const std::vector<std::pair<int, int>>& pruned,
               const std::vector<int>& merged, const std::vector<int>& refreshed) {
    int n = static_cast<int>(A.size());

    std::cout << "         A       = {";
    for (int v = 0; v < n; ++v) std::cout << (v ? ", " : "") << v << ": " << braces(A[v]);
    std::cout << "}\n         C       = {";
    for (int v = 0; v < n; ++v) {
        std::cout << (v ? ", " : "") << v << ": [";
        bool f = true;
        for (int x : C[v]) { std::cout << (f ? "" : ", ") << "c" << x; f = false; }
        std::cout << "]";
    }
    std::cout << "}\n         cliques = {";
    bool f = true;
    for (const auto& [x, members] : cliques) {
        std::cout << (f ? "" : ", ") << "c" << x << ": " << braces(members); f = false;
    }
    std::cout << "}\n         weights = {";
    for (int v = 0; v < n; ++v) std::cout << (v ? ", " : "") << v << ": " << weight[v];
    std::cout << "}   storage " << storage(A, cliques) << "\n         degrees = {";
    for (int v = 0; v < n; ++v) {
        std::cout << (v ? ", " : "") << v << ": ";
        if (eliminated[v]) std::cout << "-"; else std::cout << degrees[v];
    }
    std::cout << "}\n";

    if (!hasResult) {
        std::cout << "         neighbors = none, absorbed = none, pruned = none, "
                     "merged = none, refreshed = all (initial)\n";
        return;
    }
    std::cout << "         neighbors = " << braces(neighbors) << ", absorbed = ";
    if (absorbed.empty()) std::cout << "none";
    else { bool g = true; for (int x : absorbed) { std::cout << (g ? "" : ", ") << "c" << x; g = false; } }
    std::cout << ", pruned = ";
    if (pruned.empty()) std::cout << "none";
    else { bool g = true; for (auto [u, w] : pruned) { std::cout << (g ? "" : ", ") << u << "-" << w; g = false; } }
    std::cout << ", merged = ";
    if (merged.empty()) std::cout << "none";
    else { bool g = true; for (int x : merged) { std::cout << (g ? "" : ", ") << x; g = false; } }
    std::cout << "\n         refreshed = ";
    if (refreshed.empty()) std::cout << "none";
    else { bool g = true; for (int x : refreshed) { std::cout << (g ? "" : ", ") << x; g = false; } }
    std::cout << "\n";
}

std::vector<int> minimumDegree(const std::vector<std::set<int>>& graph) {
    int n = static_cast<int>(graph.size());
    std::size_t nnzTrilA = 0;
    for (const auto& row : graph) nnzTrilA += row.size();
    nnzTrilA = nnzTrilA / 2 + n;

    std::vector<std::set<int>> A = graph;
    std::vector<std::set<int>> C(n);
    Cliques cliques;
    std::vector<int> weight(n, 1);
    std::vector<std::vector<int>> members(n);
    for (int v = 0; v < n; ++v) members[v].push_back(v);
    std::vector<bool> eliminated(n, false);
    std::vector<int> pivots;
    std::size_t nnzL = 0;

    // The cache. Built once here, then touched only where it can be wrong.
    std::vector<int> degrees(n);
    for (int v = 0; v < n; ++v) degrees[v] = static_cast<int>(A[v].size());
    std::size_t refreshes = n;          // count them: this is the whole point

    std::cout << "start: no cliques yet, every neighbor explicit, degrees from A\n";
    showState(A, C, cliques, weight, degrees, eliminated, false, {}, {}, {}, {}, {});

    int step = 0;
    while (std::find(eliminated.begin(), eliminated.end(), false) != eliminated.end()) {
        // The pick is now a scan over cached integers, not over set unions.
        int pivot = -1;
        for (int v = 0; v < n; ++v) {
            if (eliminated[v]) continue;
            if (pivot == -1 || degrees[v] < degrees[pivot]) pivot = v;
        }
        int d = degrees[pivot];
        int w = weight[pivot];

        std::set<int> neighbors, absorbed;
        std::vector<std::pair<int, int>> pruned;
        std::vector<int> merged;
        eliminate(A, C, cliques, weight, eliminated, pivot,
                  neighbors, absorbed, pruned, merged);
        for (int i : merged)
            members[pivot].insert(members[pivot].end(), members[i].begin(), members[i].end());
        pivots.push_back(pivot);

        // REFRESH, and only here. The surviving members of the new clique are
        // exactly the vertices whose degree can have changed.
        std::vector<int> refreshed;
        for (int i : neighbors)
            if (!eliminated[i]) refreshed.push_back(i);
        for (int i : refreshed)
            degrees[i] = degreeOf(A, C, cliques, weight, i);
        refreshes += refreshed.size();

        // A supervariable of weight wp is wp consecutive columns of L. Its
        // EXTERNAL degree is what remains of the clique after the merges, since
        // a merged vertex joins the supervariable instead of neighboring it.
        int wp = weight[pivot], ext = 0;
        for (int j : cliques[pivot]) ext += weight[j];
        nnzL += static_cast<std::size_t>(wp) * ext + wp * (wp - 1) / 2 + wp;

        std::cout << "step " << step << ": eliminate " << pivot << " (degree " << d
                  << ", weight " << w << " -> " << wp << ")\n";
        showState(A, C, cliques, weight, degrees, eliminated, true,
                  neighbors, absorbed, pruned, merged, refreshed);
        ++step;
    }

    std::vector<int> order;
    for (int p : pivots)
        order.insert(order.end(), members[p].begin(), members[p].end());
    std::cout << "supervariable pivots =";
    for (int p : pivots) std::cout << " " << p;
    std::cout << "\nnnz(L) = " << nnzL << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "degree computations: " << refreshes
              << "   (md4 would do one per live vertex per step)\n";
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
    // The same three graphs as md1 and md3.
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
