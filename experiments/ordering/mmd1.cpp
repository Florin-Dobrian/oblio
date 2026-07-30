// mmd1.cpp -- multiple minimum degree.
//
// md5 finished the cheap wins. It has the quotient graph, supervariables,
// maintained degrees and degree buckets, and it returns exactly the ordering md1
// returns, only far faster. Everything left costs something.
//
// This is the first layer that changes the ANSWER. Section 5.11 of
// archive/sparse_factorization.md.
//
// The idea, from Liu (1985), is the M in MMD. Refreshing degrees is the expensive
// step, so do it less often: eliminate a whole INDEPENDENT SET of least-degree
// vertices before refreshing anything. Non-adjacent pivots cannot disturb each
// other's degrees, so every pivot in a batch is still a true minimum-degree
// vertex when it is taken.
//
// We never search for the independent set. It falls out of the bookkeeping:
// eliminating a pivot EVICTS every vertex it reached from the degree buckets, so
// whatever is still sitting in the bucket was not reached, hence is non-adjacent
// to everything already taken this round.
//
// WHAT THIS GIVES UP, and it is not what one would guess. The pivots are exact,
// but the vertices the batch evicted are invisible for the rest of the round, so
// the choice is made among the untouched remainder rather than among all
// candidates. The batch does not pick a worse vertex, it picks a different vertex
// OF THE SAME DEGREE. Minimum degree is famously sensitive to tie-breaks, so the
// fill moves by a fraction of a percent, in either direction.
//
// WHAT IS HERE, AND WHAT MMD2 ADDS. This file is the idea alone. Everything else
// genmmd does is deliberately left to mmd1's successor, which completes it:
//
//   - the PREPASS that numbers degree 0 and 1 vertices before the main loop,
//     leaving their neighbors' degrees stale (genmmd, the loop over head[1])
//   - mmdupd's q2h path. mmdelm stashes each reached vertex's pruned adjacency
//     count as fwd[rn] = nq+1, and mmdupd routes the nq==1 cases into a separate
//     list where it merges indistinguishable PAIRS. The merge test here catches
//     only vertices indistinguishable from the pivot, so MMD's supervariables are
//     at least as coarse as ours and sometimes coarser.
//   - OUTMATCHED marking, bwd[nd] = -maxint, which takes a vertex out of the
//     degree lists without merging it.
//   - the filing convention: MMD files at `dg - qsize[en] + 1` floored at 1, so
//     its least bucket is 1 where ours is 0, and it never uses bucket 0. Plus the
//     ncsub subscript statistic.
//
// NO WEIGHT ARRAY, for the same reason md3 through md5 have none: mass elimination
// merges only into the PIVOT, which is eliminated in the same call, so no live
// vertex ever stands for more than one original vertex, and a supervariable's size
// is superMembers[pivot].size() whenever it is wanted. mmd2 needs one, because its
// q2h merge folds a vertex into a LIVE one.
//
// TIE-BREAKS. Our buckets are index-ordered, *buckets[minDegree].begin(), which is
// md5's convention and the reason md1 through md5 agree. MMD's degree lists are
// linked chains prepended at head[dg], so its bucket is a stack and the winner is
// whatever was pushed last, which after construction is the highest-numbered
// vertex of that degree. There is no quality claim behind it: prepending is the
// cheap end of a linked list. We keep our convention and the orderings differ in
// ties; see the README.
//
// The tag/marker machinery with its maxint overflow reset is not modeled at all.
// It exists because the marks live in reusable integer arrays; we use std::set.
//
//
// COMPLEXITY, AND ONE PLACE THE PYTHON PAYS MORE THAN THE C++. The goal is the
// same asymptotic cost as the vendored routines, without their coding style. Two
// things were wrong and are fixed: the driver loop counts eliminations rather than
// scanning `eliminated` (O(n) per step before, O(1) now), and the mass elimination
// block strips a merged vertex from C[pivot] alone rather than from every clique,
// which is sound because I[u] was {pivot}. On a 20 by 20 grid those two cost 14800
// and 4247 elementary steps before, against 34 and 47 after, with the real
// neighbor work at 26408.
//
// What remains is min(buckets[min_degree]), which is O(bucket size) because a
// Python set is unordered. The C++ twin does not pay it: std::set is ordered, so
// *buckets[minDegree].begin() is O(1) and matches the vendored head[dg] in cost
// while keeping our index-ordered tie-break. Closing the gap in Python would need
// a heap per bucket with lazy deletion, plus a membership set to skip stale
// entries, which is more machinery than a prototype should carry. It is the one
// documented place where the Python is asymptotically worse than the C++.
//
// Build:  g++ -std=c++17 -O3 mmd1.cpp -o mmd1_cpp  (or: make)
// Run:    ./mmd1_cpp
//         ./mmd1_cpp 3      just the third example

#include <algorithm>
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

std::set<int> mmd1Neighbors(const Graph& A, const Graph& I, const Cliques& C, int u);

// Print a quotient graph with supervariables: adjacency, incidence, cliques.
void mmd1Show(const Graph& A, const Graph& I, const Cliques& C,
             const std::vector<int>& degrees, const std::string& title = "",
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
        std::cout << "  " << std::setw(width) << u << ": {" << adjacencyText.str()
                  << "} {" << incidenceText.str() << "} degree " << degrees[u] << "\n";
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

// Print the state arrays: members, eliminated, and the order so far.
void mmd1ShowState(const std::vector<int>& degrees,
                  const std::vector<std::set<int>>& buckets, int minDegree,
                  const std::vector<std::vector<int>>& superMembers,
                  const std::vector<bool>& eliminated,
                  const std::vector<int>& pivots, const std::string& title = "") {
    int n = static_cast<int>(superMembers.size());
    int width = static_cast<int>(std::to_string(std::max(n - 1, 0)).size());
    if (!title.empty()) std::cout << title << "\n";
    for (int u = 0; u < n; ++u) {
        std::string status;
        if (!eliminated[u]) status = "live";
        else if (!superMembers[u].empty()) status = "done";
        else status = "merged";
        std::ostringstream superMemberList;
        bool first = true;
        for (int v : superMembers[u]) {
            superMemberList << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::cout << "  " << std::setw(width) << u << ": members ["
                  << superMemberList.str() << "] " << status << "\n";
    }
    std::ostringstream degreesText;
    for (int u = 0; u < n; ++u)
        degreesText << (u == 0 ? "" : " ") << std::setw(width) << degrees[u];
    std::ostringstream superMembersText;
    for (int u = 0; u < n; ++u) {
        superMembersText << (u == 0 ? "" : " ") << "[";
        bool firstMember = true;
        for (int v : superMembers[u]) {
            superMembersText << (firstMember ? "" : " ") << v;
            firstMember = false;
        }
        superMembersText << "]";
    }
    std::ostringstream eliminatedText;
    for (int u = 0; u < n; ++u)
        eliminatedText << (u == 0 ? "" : " ") << std::setw(width) << (eliminated[u] ? 1 : 0);
    std::ostringstream bucketsText;
    bool firstBucket = true;
    for (std::size_t d = 0; d < buckets.size(); ++d) {
        if (buckets[d].empty()) continue;
        bucketsText << (firstBucket ? "" : "  ") << d << ": {";
        bool firstMember = true;
        for (int v : buckets[d]) {
            bucketsText << (firstMember ? "" : " ") << v;
            firstMember = false;
        }
        bucketsText << "}";
        firstBucket = false;
    }
    std::cout << "  degrees: [" << degreesText.str() << "]\n";
    std::cout << "  buckets: " << (firstBucket ? "all empty" : bucketsText.str()) << "\n";
    std::cout << "  min degree: " << minDegree << "\n";
    std::cout << "  members: [" << superMembersText.str() << "]\n";
    std::cout << "  eliminated: [" << eliminatedText.str() << "]\n";
    std::cout << "  pivots: [";
    for (std::size_t k = 0; k < pivots.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << pivots[k];
    std::cout << "]\n";
    std::cout << "  order: [";
    bool firstOrder = true;
    for (int pivot : pivots)
        for (int u : superMembers[pivot]) {
            std::cout << (firstOrder ? "" : ", ") << u;
            firstOrder = false;
        }
    std::cout << "]\n\n";
}

// Entries actually stored, as in md5. Batching changes when degrees are
// refreshed, not what the quotient graph holds.
std::size_t mmd1Storage(const Graph& A, const Graph& I, const Cliques& C) {
    std::size_t total = 0;
    for (const std::set<int>& adjacency : A) total += adjacency.size();
    for (const std::set<int>& incidence : I) total += incidence.size();
    for (const auto& [c, cliqueMembers] : C) { (void)c; total += cliqueMembers.size(); }
    return total;
}

// The neighbors of live vertex u, exactly as in md2: its explicit adjacency A[u]
// together with the members of every clique that contains u, minus u.
std::set<int> mmd1Neighbors(const Graph& A, const Graph& I, const Cliques& C, int u) {
    std::set<int> neighbors = A[u];
    for (int c : I[u]) {
        const std::set<int>& cliqueMembers = C.at(c);
        neighbors.insert(cliqueMembers.begin(), cliqueMembers.end());
    }
    neighbors.erase(u);
    return neighbors;
}

// Turn the pivot into a clique, then merge in every member it makes
// indistinguishable. Identical to md5Eliminate: this layer changes how often
// degrees are refreshed, not what an elimination does.
//
// Returns (neighbors, absorbedCliques, prunedEdges, mergedVertices): as in md2,
// plus the vertices folded into the pivot by mass elimination. The last three
// are reported for display; only neighbors is used by the caller.
std::tuple<std::set<int>, std::set<int>, std::vector<std::pair<int, int>>,
           std::vector<int>> mmd1Eliminate(
        Graph& A, Graph& I, Cliques& C, std::vector<bool>& eliminated, int pivot) {
    const std::set<int> neighbors = mmd1Neighbors(A, I, C, pivot);
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

    // Mass elimination. u is INDISTINGUISHABLE from the pivot when the two have
    // the same closed neighborhood, mmd1Neighbors(u) | {u} == mmd1Neighbors(pivot)
    // | {pivot}, as it stood before the step. Equivalently, now that the clique is
    // formed, when everything u can still reach lies inside it. The test below is
    // a cheap sufficient condition for that: nothing explicit left and no clique
    // but the new one means u sees exactly what the pivot sees, so eliminating it
    // next would cost no fill. Fold it into the pivot now and strip it from the
    // cliques, since it is no longer a vertex.
    std::vector<int> mergedVertices;
    for (int u : neighbors) {
        if (A[u].empty() && I[u].size() == 1 && *I[u].begin() == pivot) {
            I[u].clear();
            eliminated[u] = true;
            mergedVertices.push_back(u);
        }
    }
    for (int u : mergedVertices)
        C[pivot].erase(u);      // I[u] was {pivot}, so no other clique holds u

    A[pivot].clear();
    I[pivot].clear();
    eliminated[pivot] = true;
    return {neighbors, absorbedCliques, prunedEdges, mergedVertices};
}

// Move u from the bucket for its old degree to the one for newDegree. Removal
// from the middle of a bucket must be O(1), which is why a bucket is a set here;
// the vendored codes use doubly linked lists for the same reason.
void mmd1Refile(std::vector<std::set<int>>& buckets, std::vector<int>& degrees,
               int u, int newDegree) {
    buckets[degrees[u]].erase(u);
    degrees[u] = newDegree;
    buckets[newDegree].insert(u);
}

// Multiple elimination: a batch of independent pivots per degree refresh.
//
// delta widens the batch to vertices within delta of the minimum degree, which
// buys still fewer refreshes for a real concession, since those vertices are not
// minimal. delta = 0 keeps the batch to true minima. A negative delta takes one
// pivot per round, which is md5's behavior reached through this code path.
std::vector<int> mmd1MinimumDegree(const Graph& G, int delta = 0) {
    int n = static_cast<int>(G.size());
    std::size_t nnzTrilA = 0;
    for (int u = 0; u < n; ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    Graph A = G;                                  // explicit vertex neighbors
    Graph I(n);                                   // cliques that contain each vertex
    Cliques C;                                    // clique id -> member set
    std::vector<std::vector<int>> superMembers(n);   // which ones, for the expansion
    for (int u = 0; u < n; ++u) superMembers[u].push_back(u);
    std::vector<bool> eliminated(n, false);
    std::vector<int> pivots;                      // the order over supervariables
    int numEliminated = 0;                        // a counter, not a scan of eliminated
    std::size_t nnzL = 0;

    std::vector<int> degrees(n);
    for (int u = 0; u < n; ++u) degrees[u] = static_cast<int>(A[u].size());
    int numDegreeComputations = n;

    std::vector<std::set<int>> buckets(n);       // buckets[d] holds the live degree-d
    for (int u = 0; u < n; ++u) buckets[degrees[u]].insert(u);
    int minDegree = n > 0 ? *std::min_element(degrees.begin(), degrees.end()) : 0;
    int numBucketProbes = 0;
    int numRounds = 0;                           // batches, the metric this layer adds

    mmd1Show(A, I, C, degrees, "start: every edge explicit, no clique yet", &eliminated);
    mmd1ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
    while (numEliminated < n) {
        while (buckets[minDegree].empty()) {     // walk up to the first live bucket
            ++minDegree;
            ++numBucketProbes;
        }
        ++numBucketProbes;

        // ---- one BATCH, no degree refreshed inside it ----------------------
        // Take pivots from buckets [minDegree, minDegree + delta]. Eviction is
        // what keeps them independent: eliminating a pivot pulls every vertex it
        // reached out of the buckets, so whatever is still filed was not reached,
        // hence is not adjacent to anything taken this round.
        int batchLimit = minDegree + delta;
        std::vector<int> batch;
        std::set<int> touched;
        while (true) {
            if (buckets[minDegree].empty()) {    // this degree is drained
                if (minDegree >= batchLimit) break;
                ++minDegree;
                ++numBucketProbes;
                continue;
            }
            int pivot = *buckets[minDegree].begin();
            int degree = degrees[pivot];
            buckets[degree].erase(pivot);

            auto [neighbors, absorbedCliques, prunedEdges, mergedVertices] =
                mmd1Eliminate(A, I, C, eliminated, pivot);
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + static_cast<int>(mergedVertices.size());
            for (int u : mergedVertices) {        // the pivot now stands for them too
                superMembers[pivot].insert(superMembers[pivot].end(),
                                           superMembers[u].begin(), superMembers[u].end());
                superMembers[u].clear();
                buckets[degrees[u]].erase(u);
                degrees[u] = 0;
            }
            degrees[pivot] = 0;

            for (int u : C[pivot]) {              // EVICT, with a stale degree
                buckets[degrees[u]].erase(u);
                touched.insert(u);
            }

            int superSize = static_cast<int>(superMembers[pivot].size());
            int externalDegree = static_cast<int>(C[pivot].size());
            nnzL += static_cast<std::size_t>(superSize) * externalDegree
                    + static_cast<std::size_t>(superSize) * (superSize - 1) / 2
                    + superSize;

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
            std::ostringstream mergedVerticesText;
            if (mergedVertices.empty()) {
                mergedVerticesText << "none";
            } else {
                bool first = true;
                for (int u : mergedVertices) {
                    mergedVerticesText << (first ? "" : ", ") << u;
                    first = false;
                }
            }
            std::ostringstream evictedText;
            if (C[pivot].empty()) {
                evictedText << "none";
            } else {
                bool first = true;
                for (int u : C[pivot]) {
                    evictedText << (first ? "" : ", ") << u;
                    first = false;
                }
            }
            std::cout << "round " << numRounds << ": eliminate " << pivot << " (degree "
                      << degree << ", size " << superSize << ", external degree "
                      << externalDegree << "), absorbed cliques: "
                      << absorbedCliquesText.str() << ", pruned edges: "
                      << prunedEdgesText.str() << ", merged vertices: "
                      << mergedVerticesText.str() << ", evicted: " << evictedText.str() << "\n";
            if (delta < 0) break;                 // one pivot per round, as md5 does
        }

        // ---- one REFRESH, for everything the batch reached -----------------
        std::vector<int> refreshedVertices;
        for (int u : touched) if (!eliminated[u]) refreshedVertices.push_back(u);
        for (int u : refreshedVertices) {
            degrees[u] = static_cast<int>(mmd1Neighbors(A, I, C, u).size());
            buckets[degrees[u]].insert(u);
        }
        numDegreeComputations += static_cast<int>(refreshedVertices.size());
        for (int u : refreshedVertices) minDegree = std::min(minDegree, degrees[u]);
        ++numRounds;

        std::ostringstream batchText;
        for (std::size_t k = 0; k < batch.size(); ++k)
            batchText << (k == 0 ? "" : ", ") << batch[k];
        std::ostringstream refreshedVerticesText;
        if (refreshedVertices.empty()) {
            refreshedVerticesText << "none";
        } else {
            bool first = true;
            for (int u : refreshedVertices) {
                refreshedVerticesText << (first ? "" : ", ") << u;
                first = false;
            }
        }
        std::ostringstream title;
        title << "round " << (numRounds - 1) << " done: batch of " << batch.size() << ": "
              << batchText.str() << ", refreshed: " << refreshedVerticesText.str();
        mmd1Show(A, I, C, degrees, title.str(), &eliminated);
        mmd1ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
    }

    std::vector<int> order;
    for (int pivot : pivots)
        for (int u : superMembers[pivot]) order.push_back(u);
    std::cout << "nnz(L) = " << nnzL << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "degree computations: " << numDegreeComputations
              << ", bucket probes: " << numBucketProbes
              << ", rounds: " << numRounds << "\n";
    std::cout << "order: [";
    for (std::size_t k = 0; k < order.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << order[k];
    std::cout << "]\n";
    return order;
}

void run(const std::string& name, const Graph& G) {
    std::cout << "=== " << name << " ===\n";
    mmd1MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // The same three graphs as md1 and md2.
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
    // which mmd1's merge test declines a genuine supervariable. At the step whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test mmd1Neighbors(A, I, C, u) contained in
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
    // difference being the size of what merged. And superMembers ends with a hole
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

    // All of them by default. To run just one, pass its number: ./mmd1_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
