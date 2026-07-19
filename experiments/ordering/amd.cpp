// amd.cpp -- approximate minimum degree.
//
// The other fork from md6. Section 5.13 of archive/sparse_factorization.md.
//
// md6 has the quotient graph, supervariables, maintained degrees and buckets,
// and returns exactly md1's ordering. What is left costing anything is the
// refresh itself, which for each reached vertex i unites the members of every
// element in C[i] and counts the result. That union is the expensive object.
//
// MMD made the refresh RARE (mmd.cpp). AMD makes each one CHEAP, and the two are
// the same answer reached from opposite ends: do the expensive thing less.
//
// THE BOUND. Rather than uniting the elements, sum their separate contributions:
//
//   degree[i] <= min( n - k,                                 nothing exceeds what remains
//                     degree_old[i] + |L \ i|,               it can only grow by the new element
//                     |A_i \ L| + |L \ i| + sum |L_e \ L| )  over e in C_i, e != the new element
//
// The third line overcounts, because two elements may overlap outside L and the
// overlap is counted twice. So it is an UPPER BOUND, not the degree.
//
// WHY THAT IS FAST, which is the entire point and is easy to miss. The quantity
// |L_e \ L| depends only on the element e, not on the vertex i, so it is computed
// ONCE PER ELEMENT and then read by every vertex whose element list contains e.
// The exact degree costs, per vertex, a walk over the members of all its
// elements. The bound costs, per vertex, one addition per element. The two
// counters below report each, and the gap widens with the size of the elements,
// which is to say with the amount of fill, which is to say exactly where it
// matters.
//
// TWO MORE MECHANISMS, both beyond md6 and neither about the degree:
//
//   - AGGRESSIVE ABSORPTION. If |L_e \ L| == 0, element e lies entirely inside
//     the new element, so it is dead and can be absorbed at once. Ordinary
//     absorption only kills the elements the PIVOT touched; this kills elements
//     that any reached vertex touched. Cheap, since |L_e \ L| is already known.
//   - HASH SUPERVARIABLE DETECTION. Our mass-elimination test only merges
//     vertices indistinguishable from the PIVOT. Two vertices can be
//     indistinguishable from EACH OTHER without either being absorbable into the
//     pivot. We hash (A_i, C_i), compare only within a hash bucket, and merge on
//     an exact match. MMD reaches the same vertices through mmdupd's q2h list;
//     the goal is shared and the mechanism is not.
//
// WHAT IS GIVEN UP. Unlike MMD, whose pivots are exact and whose ordering moves
// only through tie-breaks, AMD can genuinely pick the wrong vertex: an
// overcounted bound can hide the true minimum. So the quality loss here is of a
// different kind, and this is the first prototype whose pivot is not guaranteed
// to be a minimum-degree vertex at all.
//
// HOW THIS DIFFERS FROM THE VENDORED Amd.cpp. One algorithmic difference, and it
// shows in the output; see section 5.13.
//
//   - THE UPDATE RUNS IN ONE PASS HERE, IN TWO THERE. Amd.cpp computes a bound
//     that excludes the new element, mass-eliminating and shrinking degme as it
//     goes, and then a SECOND loop adds the final degme to every survivor. We fold
//     both into one loop, so a vertex handled early sees a larger degme than one
//     handled late, where the vendored code gives them all the same value. Since
//     degme only shrinks, our early vertices get looser bounds. On eleven test
//     graphs this changes the ordering on four, all grids, with fill moving a few
//     tenths of a percent either way. One pass is easier to read; two is faithful.
//
// And several missing features, none of them about the degree: dense-row handling
// (amd_preprocess, AMD_DENSE), which changes orderings on matrices that have such
// rows; the postorder (amd_postorder), which leaves fill alone but changes the
// permutation; amd_aat, which forms A + A' so unsymmetric input can be ordered on
// its symmetric pattern; amd_valid; the Control and Info arrays, where we hardcode
// the choices; and the workspace compression, which is a memory strategy for flat
// arrays with no counterpart in a prototype built on sets.
//
// Build:  g++ -std=c++17 -O3 amd.cpp -o amd_cpp  (or: make)
// Run:    ./amd_cpp

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <string>
#include <vector>

using Cliques = std::map<int, std::set<int>>;

std::set<int> reachable(const std::vector<std::set<int>>& A,
                        const std::vector<std::set<int>>& C,
                        const Cliques& cliques, int i) {
    std::set<int> neighbors = A[i];
    for (int c : C[i])
        neighbors.insert(cliques.at(c).begin(), cliques.at(c).end());
    neighbors.erase(i);
    return neighbors;
}

// Kept only to measure the bound against the truth.
int exactDegree(const std::vector<std::set<int>>& A, const std::vector<std::set<int>>& C,
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

std::vector<int> amd(const std::vector<std::set<int>>& graph, bool aggressive = true) {
    int n = static_cast<int>(graph.size());
    std::size_t nnzTrilA = 0;
    for (const auto& neigh : graph) nnzTrilA += neigh.size();
    nnzTrilA = nnzTrilA / 2 + n;

    std::vector<std::set<int>> A = graph;
    std::vector<std::set<int>> C(n);
    Cliques cliques;
    std::vector<int> weight(n, 1);
    std::vector<std::vector<int>> members(n);
    for (int v = 0; v < n; ++v) members[v] = {v};
    std::vector<bool> eliminated(n, false);
    std::vector<int> pivots;
    std::size_t nnzL = 0;

    std::vector<int> degrees(n);
    for (int v = 0; v < n; ++v) degrees[v] = static_cast<int>(A[v].size());
    std::vector<std::set<int>> buckets(n + 1);
    for (int v = 0; v < n; ++v) buckets[degrees[v]].insert(v);
    int mdeg = n ? *std::min_element(degrees.begin(), degrees.end()) : 0;

    std::size_t memberVisits = 0;   // what an EXACT degree would cost
    std::size_t elementReads = 0;   // what the BOUND costs
    std::size_t overcounts = 0, boundChecks = 0;
    int eliminatedCount = 0;

    auto show = [&](const std::string& label) {
        std::cout << "  A       = {";
        for (int v = 0; v < n; ++v) {
            std::cout << (v ? ", " : "") << v << ": [";
            bool g = true;
            for (int x : A[v]) { std::cout << (g ? "" : ", ") << x; g = false; }
            std::cout << "]";
        }
        std::cout << "}\n  C       = {";
        for (int v = 0; v < n; ++v) {
            std::cout << (v ? ", " : "") << v << ": [";
            bool g = true;
            for (int x : C[v]) { std::cout << (g ? "" : ", ") << "c" << x; g = false; }
            std::cout << "]";
        }
        std::cout << "}\n  cliques = {";
        bool first = true;
        for (const auto& [c, mem] : cliques) {
            std::cout << (first ? "" : ", ") << "c" << c << ": [";
            first = false;
            bool g = true;
            for (int x : mem) { std::cout << (g ? "" : ", ") << x; g = false; }
            std::cout << "]";
        }
        std::cout << "}\n  weights = {";
        for (int v = 0; v < n; ++v)
            std::cout << (v ? ", " : "") << v << ": " << weight[v];
        std::cout << "}   storage " << storage(A, cliques) << "\n  bounds  = {";
        for (int v = 0; v < n; ++v) {
            std::cout << (v ? ", " : "") << v << ": ";
            if (eliminated[v]) std::cout << "-"; else std::cout << degrees[v];
        }
        std::cout << "}   " << label << "\n  exact   = {";
        for (int v = 0; v < n; ++v) {
            std::cout << (v ? ", " : "") << v << ": ";
            if (eliminated[v]) std::cout << "-";
            else std::cout << exactDegree(A, C, cliques, weight, v);
        }
        std::cout << "}\n  buckets = {";
        first = true;
        for (int d = 0; d <= n; ++d) {
            if (buckets[d].empty()) continue;
            std::cout << (first ? "" : ", ") << d << ": [";
            first = false;
            bool g = true;
            for (int x : buckets[d]) { std::cout << (g ? "" : ", ") << x; g = false; }
            std::cout << "]";
        }
        std::cout << "}   mdeg = " << mdeg << "\n";
    };

    std::cout << "start: no cliques yet, degrees exact (no elements to approximate over)\n";
    show("(exact at the start)");

    while (eliminatedCount < n) {
        while (mdeg <= n && buckets[mdeg].empty()) ++mdeg;
        if (mdeg > n) break;
        int pivot = *buckets[mdeg].begin();
        int dPivot = degrees[pivot];
        buckets[dPivot].erase(pivot);

        // ---- form the new element ----
        std::set<int> L = reachable(A, C, cliques, pivot);
        std::set<int> absorbed = C[pivot];
        for (int c : absorbed) cliques.erase(c);
        cliques[pivot] = L;
        int degme = 0;
        for (int j : L) degme += weight[j];

        for (int i : L) {
            std::set<int> keep;
            for (int x : A[i]) if (!L.count(x) && x != pivot) keep.insert(x);
            A[i] = keep;
            for (int c : absorbed) C[i].erase(c);
            C[i].insert(pivot);
        }
        A[pivot].clear();
        C[pivot].clear();
        eliminated[pivot] = true;
        eliminatedCount += weight[pivot];

        // ---- |L_e \ L| ONCE PER ELEMENT, and aggressive absorption ----
        std::set<int> touched;
        for (int i : L)
            for (int c : C[i]) if (c != pivot) touched.insert(c);
        std::map<int, int> we;
        std::vector<int> dead;
        for (int e : touched) {
            int outside = 0;
            for (int j : cliques[e]) if (!L.count(j)) outside += weight[j];
            memberVisits += cliques[e].size();   // the exact degree pays this PER VERTEX
            we[e] = outside;
            if (aggressive && outside == 0) dead.push_back(e);
        }
        for (int e : dead) {
            cliques.erase(e);
            for (int i = 0; i < n; ++i) C[i].erase(e);
        }

        // ---- the bound, plus mass elimination ----
        int nleft = n - eliminatedCount;
        std::vector<int> merged;
        for (int i : L) {
            if (eliminated[i]) continue;
            int outsideA = 0;
            for (int j : A[i]) outsideA += weight[j];
            std::set<int> elems;
            for (int c : C[i]) if (c != pivot) elems.insert(c);
            elementReads += elems.size();        // the bound pays only this
            int bound = outsideA;
            for (int e : elems) if (we.count(e)) bound += we[e];

            if (A[i].empty() && elems.empty()) {
                // MASS ELIMINATION: nothing outside the new clique
                weight[pivot] += weight[i];
                eliminatedCount += weight[i];
                degme -= weight[i];   // i joins the pivot, so it is no longer
                                      // part of the new element
                weight[i] = 0;
                C[i].clear();
                eliminated[i] = true;
                merged.push_back(i);
                continue;
            }

            bound += degme - weight[i];
            bound = std::min(bound, std::min(nleft - weight[i],
                                             degrees[i] + degme - weight[i]));
            int trueDegree = exactDegree(A, C, cliques, weight, i);
            ++boundChecks;
            if (bound > trueDegree) ++overcounts;
            buckets[degrees[i]].erase(i);
            degrees[i] = bound;
            buckets[bound].insert(i);
            if (bound < mdeg) mdeg = bound;
        }

        for (int i : merged) {
            members[pivot].insert(members[pivot].end(),
                                  members[i].begin(), members[i].end());
            buckets[degrees[i]].erase(i);
            cliques[pivot].erase(i);
            for (auto& [c, mem] : cliques) mem.erase(i);
        }

        // ---- HASH SUPERVARIABLE DETECTION ----
        // Vertices indistinguishable from EACH OTHER, which the pivot test above
        // cannot see. Hash first so the exact comparison runs only within a
        // bucket; the hash is a filter, never the decision.
        std::vector<int> survivors;
        for (int i : L) if (!eliminated[i]) survivors.push_back(i);
        std::map<std::pair<std::size_t, std::size_t>, std::vector<int>> byHash;
        for (int i : survivors) {
            std::size_t h1 = 0, h2 = 0;
            for (int x : A[i]) h1 = h1 * 1000003u + static_cast<std::size_t>(x) + 1;
            for (int x : C[i]) h2 = h2 * 1000003u + static_cast<std::size_t>(x) + 1;
            byHash[{h1, h2}].push_back(i);
        }
        std::vector<std::pair<int, int>> pairs;
        for (auto& [h, group] : byHash) {
            (void)h;
            if (group.size() < 2) continue;
            for (std::size_t x = 0; x < group.size(); ++x) {
                int i = group[x];
                if (eliminated[i]) continue;
                for (std::size_t y = x + 1; y < group.size(); ++y) {
                    int j = group[y];
                    if (eliminated[j]) continue;
                    std::set<int> ai = A[i], aj = A[j];
                    ai.erase(j);
                    aj.erase(i);
                    if (ai == aj && C[i] == C[j]) {
                        weight[i] += weight[j];
                        weight[j] = 0;
                        members[i].insert(members[i].end(),
                                          members[j].begin(), members[j].end());
                        buckets[degrees[j]].erase(j);
                        A[j].clear();
                        C[j].clear();
                        eliminated[j] = true;
                        // j is absorbed into the live supervariable i, not
                        // eliminated: its weight is counted when i is chosen
                        for (auto& [c, mem] : cliques) mem.erase(j);
                        for (int v = 0; v < n; ++v) A[v].erase(j);
                        pairs.push_back({i, j});
                    }
                }
            }
        }

        pivots.push_back(pivot);
        int wp = weight[pivot], ext = 0;
        for (int j : cliques[pivot]) ext += weight[j];
        nnzL += static_cast<std::size_t>(wp) * ext + wp * (wp - 1) / 2 + wp;

        std::cout << "eliminate " << pivot << " (bound " << dPivot
                  << ", weight -> " << wp << ")   merged = ";
        if (merged.empty()) std::cout << "none";
        else { bool g = true; for (int x : merged) { std::cout << (g ? "" : ", ") << x; g = false; } }
        std::cout << "   absorbed = ";
        if (dead.empty()) std::cout << "none";
        else { std::sort(dead.begin(), dead.end()); bool g = true;
               for (int x : dead) { std::cout << (g ? "" : ", ") << "c" << x; g = false; } }
        std::cout << "   hash-merged = ";
        if (pairs.empty()) std::cout << "none";
        else { bool g = true; for (auto& [i, j] : pairs) {
                   std::cout << (g ? "" : ", ") << i << "<-" << j; g = false; } }
        std::cout << "\n";
        show("(bounds after)");
    }

    std::vector<int> order;
    for (int p : pivots)
        order.insert(order.end(), members[p].begin(), members[p].end());

    std::cout << "nnz(L) = " << nnzL << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "element-member visits an exact degree would need: "
              << memberVisits << "\n";
    std::cout << "element reads the bound needed:                   "
              << elementReads << "\n";
    std::cout << "bound was loose " << overcounts << " times out of "
              << boundChecks << "\n";
    return order;
}

int main() {
    std::vector<std::pair<std::string, std::vector<std::set<int>>>> tests = {
        {"graph1", {{1, 3}, {0, 2}, {1, 3}, {0, 2}}},
        {"graph2", {{1, 2}, {0, 3}, {0, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4}}},
        {"graph3", {{1, 3, 8}, {0, 2, 6, 8}, {1, 3, 5}, {0, 2, 4}, {3, 5},
                    {2, 4, 6, 9}, {1, 5, 7, 10}, {6, 8}, {0, 1, 7, 9},
                    {5, 8, 10}, {6, 9, 11}, {10}}},
        // graph4, eight vertices and fourteen edges: the smallest graph we could
        // find on which the BOUND is ever loose. The bound overcounts only when a
        // vertex belongs to two elements that overlap outside the new one, which
        // needs enough eliminations to have made several elements and enough fill
        // for them to intersect. Every connected graph on five or six vertices is
        // exact (checked exhaustively), and so are graph1 to graph3, so without
        // this one the trace would never show the approximation approximating.
        //   edges: 0-2 0-3 0-4 0-7 1-3 1-4 1-6 1-7 2-3 2-5 3-6 3-7 4-5 5-6
        {"graph4", {{2, 3, 4, 7}, {3, 4, 6, 7}, {0, 3, 5}, {0, 1, 2, 6, 7},
                    {0, 1, 5}, {2, 4, 6}, {1, 3, 5}, {0, 1, 3}}},
    };
    for (const auto& [name, graph] : tests) {
        std::cout << "=== " << name << " ===\n";
        std::vector<int> order = amd(graph);
        std::cout << "order:";
        for (int v : order) std::cout << " " << v;
        std::cout << "\n\n";
    }
    return 0;
}
