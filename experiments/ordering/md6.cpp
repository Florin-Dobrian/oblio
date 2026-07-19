// md6.cpp -- minimum degree, step 6: degree buckets.
//
// md5 stopped recomputing degrees that could not have changed. What it left in
// place is the scan: the picker still walks every live vertex to find the
// smallest cached degree, O(n) per step, now over integers rather than set
// unions. Cheap, but still the only remaining O(n) per pivot.
//
// The fix is to file each supervariable in a bucket indexed by its degree, so the
// minimum can be found by walking UP from the last known minimum rather than
// looking at everything. Section 5.9 of archive/sparse_factorization.md describes
// this as common ground: both MMD and AMD do it, neither invented it.
//
// Two things make the walk cheap:
//
//   - mdeg, a LOWER BOUND on the current minimum degree. Each search starts at
//     mdeg rather than at 0, and every bucket below it is known to be empty.
//   - a vertex whose degree changes must be pulled out of the middle of its old
//     bucket, so buckets need O(1) removal. Here that is a std::set; the vendored
//     codes use doubly linked lists (MMD's fwd/bwd, AMD's Next/Last) because
//     Fortran and C of that era had no such container.
//
// Keeping mdeg correct is the whole of the difficulty, and it is a lower bound
// rather than the true minimum on purpose: it may lag, and the walk fixes it.
// What it must never do is overshoot, since a bucket below mdeg is never examined
// and a vertex sitting there would never be chosen.
//
// Build:  g++ -std=c++17 -O3 md6.cpp -o md6_cpp  (or: make)
// Run:    ./md6_cpp

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
               const std::vector<std::set<int>>& buckets, int mdeg,
               bool hasResult,
               const std::set<int>& neighbors, const std::set<int>& absorbed,
               const std::vector<std::pair<int, int>>& pruned,
               const std::vector<int>& merged, const std::vector<int>& refreshed,
               int walked) {
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
    std::cout << "}\n         buckets = {";
    { bool g = true;
      for (std::size_t d = 0; d < buckets.size(); ++d)
          if (!buckets[d].empty()) {
              std::cout << (g ? "" : ", ") << d << ": " << braces(buckets[d]); g = false;
          }
    }
    std::cout << "}   mdeg = " << mdeg << "\n";

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
    std::cout << "   (walked " << walked << " bucket(s) to find the pivot)\n";
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

    std::vector<int> degrees(n);
    for (int v = 0; v < n; ++v) degrees[v] = static_cast<int>(A[v].size());

    // Vertices filed by degree, so the pick is a walk rather than a scan.
    std::vector<std::set<int>> buckets(n + 1);
    for (int v = 0; v < n; ++v) buckets[degrees[v]].insert(v);
    int mdeg = n ? *std::min_element(degrees.begin(), degrees.end()) : 0;
    std::size_t probes = 0;             // bucket slots examined, the metric

    // Move i to the bucket for its new degree. The bound may only ever fall.
    auto refile = [&](int i, int newDegree) {
        buckets[degrees[i]].erase(i);
        degrees[i] = newDegree;
        buckets[newDegree].insert(i);
        if (newDegree < mdeg) mdeg = newDegree;
    };

    std::cout << "start: no cliques yet, degrees from A, vertices filed by degree\n";
    showState(A, C, cliques, weight, degrees, eliminated, buckets, mdeg,
              false, {}, {}, {}, {}, {}, 0);

    int step = 0;
    while (std::find(eliminated.begin(), eliminated.end(), false) != eliminated.end()) {
        // PICK: walk up from mdeg to the first non-empty bucket. No scan.
        int walked = 0;
        while (buckets[mdeg].empty()) { ++mdeg; ++walked; }
        probes += walked + 1;
        int pivot = *buckets[mdeg].begin();   // lowest index, matching md5's tie-break
        int d = degrees[pivot];
        int w = weight[pivot];
        buckets[d].erase(pivot);

        std::set<int> neighbors, absorbed;
        std::vector<std::pair<int, int>> pruned;
        std::vector<int> merged;
        eliminate(A, C, cliques, weight, eliminated, pivot,
                  neighbors, absorbed, pruned, merged);
        for (int i : merged) {
            members[pivot].insert(members[pivot].end(), members[i].begin(), members[i].end());
            buckets[degrees[i]].erase(i);     // merged vertices leave the buckets
        }
        pivots.push_back(pivot);

        // REFRESH, and only here. The surviving members of the new clique are
        // exactly the vertices whose degree can have changed.
        std::vector<int> refreshed;
        for (int i : neighbors)
            if (!eliminated[i]) refreshed.push_back(i);
        for (int i : refreshed)
            refile(i, degreeOf(A, C, cliques, weight, i));

        // A supervariable of weight wp is wp consecutive columns of L. Its
        // EXTERNAL degree is what remains of the clique after the merges, since
        // a merged vertex joins the supervariable instead of neighboring it.
        int wp = weight[pivot], ext = 0;
        for (int j : cliques[pivot]) ext += weight[j];
        nnzL += static_cast<std::size_t>(wp) * ext + wp * (wp - 1) / 2 + wp;

        std::cout << "step " << step << ": eliminate " << pivot << " (degree " << d
                  << ", weight " << w << " -> " << wp << ")\n";
        showState(A, C, cliques, weight, degrees, eliminated, buckets, mdeg, true,
                  neighbors, absorbed, pruned, merged, refreshed, walked);
        ++step;
    }

    std::vector<int> order;
    for (int p : pivots)
        order.insert(order.end(), members[p].begin(), members[p].end());
    std::cout << "supervariable pivots =";
    for (int p : pivots) std::cout << " " << p;
    std::cout << "\nnnz(L) = " << nnzL << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "bucket probes: " << probes
              << "   (md5 would scan every live vertex per step)\n";
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
