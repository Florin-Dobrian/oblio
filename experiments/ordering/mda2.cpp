// mda2.cpp -- minimum degree, iteration 2a: the quotient graph, with the degree bounded.
//
// md2 with one call site changed: the picker asks for a BOUND on the degree instead
// of the degree. Everything else, the quotient graph, the elimination, the pruning,
// the absorption, is md2's and untouched.
//
// This is a fork, not a rung. md1 has no cliques, so there is nothing to
// approximate; md2 is the FIRST layer where the choice exists, because the degree
// has become a union over a structure rather than a length.
//
//     reach(u)  = ( A[u] | C[c] for every c in I[u] ) - {u}
//     degree(u) = |reach(u)|                                   a union
//     bound(u)  = |A[u]| + sum over c in I[u] of ( |C[c]| - 1 ) additions only
//
// Why the minus one, and why it is exact: u is in C[c] for every c in I[u], that
// being what an incidence means, and reach excludes u. Why A[u] contributes its
// whole length: it cannot overlap any of the cliques, because joining c pruned C[c]
// out of A[u] and A[u] only ever shrinks after. So the ONLY source of overcounting
// is one vertex lying in two cliques of I[u].
//
// This is the PIVOT-FREE bound, and it is not the tighter one. The bound used from
// mdam2 onward, and by amd1, is stated against the new clique:
//
//     PIVOT-FREE            sum over c in I[u]              of ( |C[c]| - 1 )
//     AGAINST C[pivot]      sum over c in I[u], c != pivot, of |C[c] - C[pivot]|
//
// and the second is strictly better, since u is one member of C[pivot]:
//
//     |C[c] - C[pivot]|   <=   |C[c] - {u}|   =   |C[c]| - 1
//
// mda2 cannot use it. That form is stated against C[pivot], so it applies only to a
// vertex IN C[pivot], and this picker produces a number for every live vertex at
// every iteration, touched or not. What makes the tighter form available is MAINTAINED
// degrees, which narrow the refresh set to exactly C[pivot]'s members: that is
// mdam2. See the README section "Zooming in on md2", and 5.13 of
// archive/sparse_factorization.md.
//
// The trace prints bound and degree side by side, so the looseness is visible
// rather than argued. Computing the exact degree costs the very union the bound
// exists to avoid; it is instrumentation and nothing else.
//
// What follows is md2's own description of the quotient graph, unchanged.
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
// The order can differ from md2's, and that is the result rather than a defect:
// where the bound is loose the picker may prefer a vertex the exact degree would
// not have chosen.
//
// Build:  g++ -std=c++17 -O3 mda2.cpp -o mda2_cpp  (or: make)
// Run:    ./mda2_cpp
//         ./mda2_cpp 3      just the third example

#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Plain vectors, UNSORTED, and a vector indexed by clique id, not std::map. A set
// costs O(log d) per membership test and per insertion; keeping a vector sorted
// costs a merge per union. Neither is needed: membership comes from a MARK array
// stamped with a tag, so "is v in the new clique" is one comparison, and every
// pass is linear in what it touches. That is what the vendored codes and Oblio's
// own SymFactorEngine do. See the README section on complexity.
//
// Types follow Oblio's rule: an INDEX names a vertex or a clique and is a
// std::int32_t, with NIL for "none"; a POSITION locates something inside a vector
// and is a std::size_t.
constexpr std::int32_t NIL = -1;

// The mark array is a set and the tag names it, so a tag must never repeat: a
// repeat makes a stale stamp read as a match, which is wrong silently. The tag
// only ever climbs, so the ceiling is where it has to be swept back. Half the
// positive range of std::int32_t, which is a pragmatic choice and not a derived
// one: nothing here stores anything but a tag, so the true ceiling is the type's
// own maximum, and the room left over is against a later layer wanting some of it.
constexpr std::int32_t TAG_CEILING = (1 << 30) - 1;

// Above this n, nothing is printed from inside the run: no initial state, no
// per-iteration trace. That output is for reading a small example by eye, and at any
// size worth calling large it is O(n) lines of O(n) each, so it is unreadable and slow
// to produce. What still prints at every size is the end of the run, the counters and
// the order, since each is O(1) lines and that is what the twin comparison comes down
// to. To watch a larger run, raise this.
constexpr std::size_t SHOW_THRESHOLD = 32;

using Graph = std::vector<std::vector<std::int32_t>>;

// C[c] holds the members of clique c, and cliqueLive[c] says whether c exists.
// A clique id is the pivot that created it, so the id space is the vertex space.
struct Cliques {
    std::vector<std::vector<std::int32_t>> members;
    std::vector<bool> live;
    std::size_t count = 0;

    explicit Cliques(std::size_t n) : members(n), live(n, false) {}
    const std::vector<std::int32_t>& at(std::int32_t c) const { return members[c]; }
    std::vector<std::int32_t>& operator[](std::int32_t c) { return members[c]; }
    void create(std::int32_t c, std::vector<std::int32_t> m) {
        if (!live[c]) ++count;
        live[c] = true;
        members[c] = std::move(m);
    }
    void erase(std::int32_t c) {
        if (live[c]) --count;
        live[c] = false;
        members[c].clear();
    }
    std::size_t size() const { return count; }
};

// I[u] cliques that contain u
// C[c] vertices that c contains

std::size_t mda2Bound(const Graph& A, const Graph& I, const Cliques& C, std::int32_t u);

std::vector<std::int32_t> mda2Neighbors(const Graph& A, const Graph& I, const Cliques& C,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u);

// Print a quotient graph: adjacency, incidence, cliques, in the order the
// structure holds them.
void mda2Show(const Graph& A, const Graph& I, const Cliques& C,
             std::vector<std::int32_t>& mark, std::int32_t& tag,
             const std::string& title = "",
             const std::vector<bool>* eliminated = nullptr) {
    const std::size_t n = A.size();
    int width = static_cast<int>(std::to_string(n > 0 ? n - 1 : 0).size());
    std::vector<std::int32_t> aliveVertices;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        if (eliminated == nullptr || !(*eliminated)[u]) aliveVertices.push_back(u);
    std::size_t numAliveEdges = 0;
    for (std::int32_t u : aliveVertices) numAliveEdges += A[u].size();
    numAliveEdges /= 2;
    std::size_t numAliveIncidences = 0;
    for (std::int32_t u : aliveVertices) numAliveIncidences += I[u].size();
    std::size_t numAliveCliques = C.size();
    std::size_t loose = 0;
    if (!title.empty()) std::cout << title << "\n";
    std::ostringstream aliveVerticesText;
    if (eliminated == nullptr) aliveVerticesText << n;
    else aliveVerticesText << aliveVertices.size() << " of " << n;
    std::cout << "num alive vertices = " << aliveVerticesText.str()
              << ", num alive edges = " << numAliveEdges
              << ", num alive cliques = " << numAliveCliques
              << ", storage = " << 2 * numAliveEdges << " + " << 2 * numAliveIncidences
              << " = " << 2 * (numAliveEdges + numAliveIncidences) << "\n";
    for (std::int32_t u : aliveVertices) {
        std::ostringstream adjacencyText;
        bool first = true;
        for (std::int32_t v : A[u]) {
            adjacencyText << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::ostringstream incidenceText;
        first = true;
        for (std::int32_t c : I[u]) {
            incidenceText << (first ? "" : " ") << "c" << c;
            first = false;
        }
        std::size_t degree = mda2Neighbors(A, I, C, mark, tag, u).size();
        std::size_t bound = mda2Bound(A, I, C, u);
        if (bound > degree) ++loose;
        std::cout << "  " << std::setw(width) << u << ": {" << adjacencyText.str()
                  << "} {" << incidenceText.str() << "} bound " << bound
                  << " degree " << degree << "\n";
    }
    for (std::int32_t c = 0; c < static_cast<std::int32_t>(n); ++c) {
        if (!C.live[c]) continue;
        std::ostringstream cliqueMembersText;
        bool first = true;
        for (std::int32_t u : C.at(c)) {
            cliqueMembersText << (first ? "" : " ") << std::setw(width) << u;
            first = false;
        }
        std::cout << "  c" << c << ": {" << cliqueMembersText.str() << "}\n";
    }
    if (loose)
        std::cout << "  bound loose on " << loose << " of " << aliveVertices.size()
                  << " live vertices\n";
    std::cout << "\n";
}

// An upper bound on |reach(u)|, computed without forming reach(u).
//
//     bound(u) = |A[u]|                               -> A[u].size()
//              + sum over c in I[u] of ( |C[c]| - 1 )  -> sum of C[c].size() - 1
//
// TWO terms, where the bound against C[pivot] has three: there is no designated
// clique to separate out here, so where u belongs to the new clique it is simply one
// more c in the sum. And nothing is precomputed, since |C[c]| - 1 is a length read
// straight off C[c]; this form needs no array where the other needs one.
std::size_t mda2Bound(const Graph& A, const Graph& I, const Cliques& C, std::int32_t u) {
    std::size_t bound = A[u].size();
    for (std::int32_t c : I[u]) bound += C.at(c).size() - 1;
    return bound;
}

// Entries actually stored. Each edge costs two, one per endpoint in A. Each
// incidence costs two as well, the clique id in I and the member in C. Watch
// the total fall monotonically; the naive graph's only rises.
std::size_t mda2Storage(const Graph& A, const Graph& I, const Cliques& C) {
    std::size_t total = 0;
    for (const std::vector<std::int32_t>& adjacency : A) total += adjacency.size();
    for (const std::vector<std::int32_t>& incidence : I) total += incidence.size();
    for (std::size_t c = 0; c < C.members.size(); ++c)
        if (C.live[c]) total += C.members[c].size();
    return total;
}

// The neighbors of live vertex u: its explicit adjacency A[u] together with the
// members of every clique that contains u, minus u itself, which the cliques
// always carry. This is George and Liu's reachable set, and it is what the
// elimination graph would hold explicitly.
std::vector<std::int32_t> mda2Neighbors(const Graph& A, const Graph& I, const Cliques& C,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u) {
    // In set terms this is one line, and it is worth keeping in view because the
    // code below is that line with the set taken away:
    //
    //     reach(u) = ( A[u] | C[c] for every c in I[u] ) - {u}
    //
    // The mark array IS the set. mark[v] == tag is the membership test, one
    // comparison; mark[v] = tag is the insertion, one store. So the union costs one
    // pass per source rather than a hash per member, and nothing is allocated.
    //
    // One pass per source, with the mark array doing the deduplication, so the
    // cost is linear in what is touched. Nothing is sorted: the order is the order
    // the sources were walked in.
    ++tag;
    std::vector<std::int32_t> neighbors;
    mark[u] = tag;                          // never its own neighbor
    for (std::int32_t v : A[u]) { mark[v] = tag; neighbors.push_back(v); }
    for (std::int32_t c : I[u])
        for (std::int32_t v : C.at(c))
            if (mark[v] != tag) { mark[v] = tag; neighbors.push_back(v); }
    return neighbors;
}

// Turn the pivot into a clique.
//
// Returns (neighbors, absorbedCliques, prunedEdges): the pivot's neighbor set,
// which becomes the clique and the pattern of its column of L; the cliques that
// the new one swallows; and the explicit edges the new clique makes redundant.
// The last two are reported for display; only neighbors is used by the caller.//
// Set view of the whole function, in the order the code does it:
//
//     C[pivot] = reach(pivot)                    absorb into C[pivot]
//     C        = C - I[pivot]                    reclaim I[pivot]
//     for u in C[pivot]:
//         A[u] = A[u] - C[pivot] - {pivot}       prune
//         I[u] = ( I[u] - I[pivot] ) | {pivot}   absorb into C[pivot], reclaim I[pivot]
//
// The new clique is C[pivot] and gets no name of its own, so the first line reads
// as what an elimination IS: the pivot stops being a vertex with a reachable set
// and becomes a clique holding that same set. The last line is the first two
// written on the I side, since u is in C[c] exactly when c is in I[u].
//
// Three set differences, and not one of them builds a set. Each is a single stamp
// of the subtrahend followed by one compaction pass over the minuend, which turns
// |A[u]| * |C[pivot]| comparisons into |A[u]| + |C[pivot]|.
//
// This file bounds the degree and this function is untouched by that, which is
// worth saying here rather than only in the header. The new clique IS the
// reachable set, so forming it needs the members and not a count and there is
// nothing to approximate. The approximation buys nothing at this call and
// everything at the picker's, which is the asymmetry the whole bounded branch
// rests on: one union per pivot either way, and the exact branch pays another per
// vertex it refreshes.
std::tuple<std::vector<std::int32_t>, std::vector<std::int32_t>,
           std::vector<std::pair<std::int32_t, std::int32_t>>>
mda2Eliminate(Graph& A, Graph& I, Cliques& C, std::vector<bool>& eliminated,
             std::vector<std::int32_t>& mark, std::int32_t& tag, std::int32_t pivot) {
    const std::vector<std::int32_t> neighbors = mda2Neighbors(A, I, C, mark, tag, pivot);
    const std::vector<std::int32_t> absorbedCliques = I[pivot];
    for (std::int32_t c : absorbedCliques)
        C.erase(c);
    C.create(pivot, neighbors);     // becomes the column pattern of the pivot

    // Stamp the new clique once, and the absorbed cliques once. Membership is then
    // a comparison, and both loops below are compactions in place. cliqueTag is the
    // set C[pivot] and absorbedTag is the set I[pivot], each built in one pass and
    // then queried for free.
    ++tag;
    const std::int32_t cliqueTag = tag;
    for (std::int32_t v : neighbors) mark[v] = cliqueTag;
    ++tag;
    const std::int32_t absorbedTag = tag;
    for (std::int32_t c : absorbedCliques) mark[c] = absorbedTag;

    std::vector<std::pair<std::int32_t, std::int32_t>> prunedEdges;
    std::vector<std::int32_t> kept;
    for (std::int32_t u : neighbors) {
        kept.clear();
        for (std::int32_t v : A[u]) {
            if (v == pivot) continue;            // the pivot is no longer a variable
            if (mark[v] == cliqueTag) {          // both ends inside the new clique
                if (u < v) prunedEdges.push_back({u, v});
                continue;                        // implicit now: drop the explicit copy
            }
            kept.push_back(v);
        }
        A[u].swap(kept);                         // what survives is A[u] - C[pivot] - {pivot}

        kept.clear();                            // I[u] loses the absorbed cliques
        for (std::int32_t c : I[u])
            if (mark[c] != absorbedTag) kept.push_back(c);
        kept.push_back(pivot);                   // u joins the new clique, id = pivot
        I[u].swap(kept);
    }

    A[pivot].clear();
    I[pivot].clear();
    eliminated[pivot] = true;
    return {neighbors, absorbedCliques, prunedEdges};
}

// Same heuristic as md1, on the quotient graph. No fill is ever stored.
std::vector<std::int32_t> mda2MinimumDegree(const Graph& G) {
    const std::size_t n = G.size();
    std::size_t nnzTrilA = 0;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    Graph A = G;                            // explicit vertex neighbors
    Graph I(n);                             // cliques each vertex belongs to
    Cliques C(n);   // clique id -> member list
    std::vector<std::int32_t> mark(n, NIL);    // scratch for membership, with tag
    std::int32_t tag = 0;
    // Calls to the eliminate procedure, one per pivot. Not the count of vertices
    // removed: a pivot can carry mass-merged vertices out with it, and from mmd1 up
    // an iteration batches several eliminations before one degree update pass. The three
    // counts coincide only where both of those are absent.
    std::size_t numEliminations = 0;
    // Summed over the eliminations, |C[p]| being the new clique AFTER the trim, so
    // in supernodal terms the update rather than the front. It is the raw reach of
    // the eliminations, undeduplicated: where a layer deduplicates, the degree
    // update count comes out below this, and the gap is what the batching saved.
    // In md2 it is nnz(L) - n, there being no mass elimination to shrink a clique.
    std::size_t numCliqueEntries = 0;
    // Passes of the outer loop, each one a batch of eliminations followed by one
    // degree update pass. Here the batch is always a single elimination, so this
    // equals numEliminations; from mmd1 up the two come apart.
    std::size_t numIterations = 0;
    // Every bound this layer computes. No split into a build and updates, because
    // nothing is maintained: the picker recomputes each candidate's bound from
    // scratch on every iteration, as md2 and md3 do with exact degrees.
    std::size_t numBoundComputations = 0;
    // Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    // is here as the witness that the guard is inert rather than as a statistic.
    std::size_t numTagSweeps = 0;
    std::vector<bool> eliminated(n, false);
    std::vector<std::int32_t> order;
    std::size_t degreeSum = 0;
    // NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    std::size_t loosePicks = 0;

    // NOT PRODUCTION: display only. The trace is what makes these files teachable and
    // is the whole reason they exist; nothing downstream reads it.
    if (n <= SHOW_THRESHOLD) {
        mda2Show(A, I, C, mark, tag, "start: every edge explicit, no clique yet", &eliminated);
    }
    for (std::int32_t iteration = 0; iteration < static_cast<std::int32_t>(n); ++iteration) {
        ++numIterations;
        std::int32_t pivot = NIL;          // O(n) scan, O(live) additions, no unions
        std::size_t best = 0;
        for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
            if (eliminated[u]) continue;
            ++numBoundComputations;
            std::size_t candidateBound = mda2Bound(A, I, C, u);
            if (pivot == NIL || candidateBound < best) { pivot = u; best = candidateBound; }
        }
        // Sweep the tag back before it can wrap. One site in this layer, unlike its
        // three neighbors in the square: mda2Bound reads lengths and takes no mark
        // or tag, and the recomputing column has no refresh, so the elimination is
        // the only region that spends a tag. Not inside mda2Eliminate, which holds
        // cliqueTag and absorbedTag live across the whole prune loop. Never
        // observed to fire.
        if (tag > TAG_CEILING) {
            std::fill(mark.begin(), mark.end(), NIL);
            tag = 0;
            ++numTagSweeps;
        }
        auto [neighbors, absorbedCliques, prunedEdges] =
            mda2Eliminate(A, I, C, eliminated, mark, tag, pivot);
        ++numEliminations;
        numCliqueEntries += C[pivot].size();
        std::size_t degree = neighbors.size();
        // NOT PRODUCTION: instrumentation, counting how often the bound was loose.
        // No union is needed for it: the eliminator has just formed reach(pivot) as
        // the new clique, so degree IS the pivot's exact degree, free of charge.
        if (best > degree) ++loosePicks;
        order.push_back(pivot);
        degreeSum += degree;

        std::ostringstream absorbedCliquesText;
        if (absorbedCliques.empty()) {
            absorbedCliquesText << "none";
        } else {
            bool first = true;
            for (std::int32_t c : absorbedCliques) {
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
        title << "iteration " << iteration << ": eliminate " << pivot << " (bound " << best
              << ", degree " << degree
              << "), absorbed cliques: " << absorbedCliquesText.str()
              << ", pruned edges: " << prunedEdgesText.str();
        // NOT PRODUCTION: display only. The trace is what makes these files teachable and
        // is the whole reason they exist; nothing downstream reads it.
        if (n <= SHOW_THRESHOLD) {
            mda2Show(A, I, C, mark, tag, title.str(), &eliminated);
        }
    }

    std::size_t nnzL = degreeSum + n;
    std::cout << "n = " << n << ", nnz(L) = " << nnzL
              << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "iterations: " << numIterations << "\n";
    std::cout << "eliminations: " << numEliminations << "\n";
    std::cout << "sum of |C[p]|: " << numCliqueEntries << "\n";
    // NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    std::cout << "loose picks = " << loosePicks << " of " << n << "\n";
    std::cout << "bound computations: " << numBoundComputations << "\n";
    std::cout << "tag sweeps: " << numTagSweeps << "\n";
    std::cout << "order: [";
    for (std::size_t k = 0; k < order.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << order[k];
    std::cout << "]\n";
    return order;
}

// A square grid graph, four-neighbor, for running the counters at a size the seven examples
// cannot reach. It is here rather than among them because it is not an example: nothing about it
// illustrates a mechanism, and above SHOW_THRESHOLD its trace is not printed at all.
//
// It must match the Python twin's grid_graph exactly, vertex for vertex, or `make test` would be
// diffing two different problems.
static Graph gridGraph(int side) {
    const int n = side * side;
    Graph graph(n);
    for (int r = 0; r < side; ++r)
        for (int c = 0; c < side; ++c) {
            const int u = r * side + c;
            if (r > 0)        graph[u].push_back(u - side);
            if (c > 0)        graph[u].push_back(u - 1);
            if (c + 1 < side) graph[u].push_back(u + 1);
            if (r + 1 < side) graph[u].push_back(u + side);
        }
    return graph;
}

void run(const std::string& name, const Graph& G) {
    std::cout << "=== " << name << " ===\n";
    mda2MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Grid mode, spelled the same way in every layer and in both twins, so `make test` can diff
    // at a size the seven examples cannot reach. Nothing is filtered: the run is silent above
    // SHOW_THRESHOLD and prints its closing lines as always.
    //
    //   ./mda2_cpp grid 22
    if (argc > 2 && std::string(argv[1]) == "grid") {
        const int side = std::atoi(argv[2]);
        std::cout << "=== grid " << side << "x" << side << " (n = " << side * side << ") ===\n";
        mda2MinimumDegree(gridGraph(side));
        return 0;
    }

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
        {1, 3}, {0, 2}, {1, 3}, {0, 2},
    };
    Graph graph2 = {
        {1, 2}, {0, 3}, {0, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4},
    };
    Graph graph3 = {
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
    Graph graph4 = {
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
    // which md3's merge test declines a genuine supervariable. At the iteration whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test md3Neighbors(A, I, C, u) contained in
    // C[pivot] would merge it. See the README section on mass elimination.
    //
    //   edges: 0-3 0-4 1-2 1-4
    Graph graph5 = {
        {3, 4},           // 0
        {2, 4},           // 1
        {1},              // 2
        {0},              // 3
        {0, 1},           // 4
    };

    // graph6, six vertices and eight edges. Here because one small graph carries
    // three things at once. Its supervariable {0, 4} is a supernode but NOT a
    // fundamental one: the elimination forest is 2 -> 1 -> 4 and 3 -> 0 -> 4, so
    // 4 already has 1 as a child when 0 merges into it. The merge happens at iteration
    // 2 of 5, so the run continues afterwards and the selection degree, 3 over
    // {2, 3, 4}, differs from the external degree, 2 over {2, 3}, with the
    // difference being the weight that merged. And superMembers ends with a hole
    // in the middle, slot 4 empty between two used ones, while no pivot equals
    // its own iteration number. See the README sections on mass elimination and on
    // external degree.
    //
    //   edges: 0-2 0-3 0-4 1-3 2-3 2-4 2-5 3-4
    Graph graph6 = {
        {2, 3, 4},        // 0
        {3},              // 1
        {0, 3, 4, 5},     // 2
        {0, 1, 2, 4},     // 3
        {0, 2, 3},        // 4
        {2},              // 5
    };

    // graph7, five vertices and six edges. The pairwise case: at the iteration whose
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
    Graph graph7 = {
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

    // All of them by default. To run just one, pass its number: ./mda2_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (std::size_t number = 1; number <= examples.size(); ++number) {
        if (selected != 0 && number != static_cast<std::size_t>(selected)) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
