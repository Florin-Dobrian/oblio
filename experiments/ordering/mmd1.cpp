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
// It exists because the marks live in reusable integer arrays; ours are the mark
// and touchedRound arrays, which do the same job with an explicit tag.
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
// The containers are flat: A and I are sorted vectors, C is indexed by clique id,
// membership comes from a mark array with a tag, and a bucket is a linked list,
// head[d] with next and prev over n, so filing, unfiling and popping are O(1).
// With that this file performs the same operations at the same cost as the
// vendored genmmd; what it does not yet have is genmmd's remaining features,
// which are mmd2's business.
//
// The Python mirrors the buckets with a list whose position 0 is the head, so both
// twins hold the same sequence and pick the same pivot, and it pays O(bucket) for
// insert and remove where the C++ splices in O(1). That is the one place the
// Python is asymptotically worse than its twin; A, I and C stay sets there,
// because set algebra is what makes the layer readable.
//
// Build:  g++ -std=c++17 -O3 mmd1.cpp -o mmd1_cpp  (or: make)
// Run:    ./mmd1_cpp
//         ./mmd1_cpp 3      just the third example

#include <cstdlib>
#include <algorithm>
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

std::vector<std::int32_t> mmd1Neighbors(const Graph& A, const Graph& I, const Cliques& C,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u);

// The degree buckets, as the vendored codes hold them: one doubly linked list per
// degree, threaded through arrays of size n. Push, pop and splice are all O(1),
// which an ordered container cannot give. MMD spells these fwd/bwd and AMD
// Next/Last. The Python twin mirrors the same sequence with a list whose position
// 0 is the head, so both pick the same pivot.
struct Buckets {
    std::vector<std::int32_t> head;   // head[d], the first live vertex of degree d
    std::vector<std::int32_t> next;   // next[u], toward the tail
    std::vector<std::int32_t> prev;   // prev[u], toward the head
    std::vector<bool> filed;          // whether u is in a bucket at all

    explicit Buckets(std::size_t n)
        : head(n, NIL), next(n, NIL), prev(n, NIL), filed(n, false) {}

    void file(std::size_t d, std::int32_t u) {          // buckets[d].add(u), at the head
        next[u] = head[d];
        prev[u] = NIL;
        if (head[d] != NIL) prev[head[d]] = u;
        head[d] = u;
        filed[u] = true;
    }
    void unfile(std::size_t d, std::int32_t u) {        // buckets[d].discard(u)
        if (!filed[u]) return;                          // idempotent, as set.discard was
        if (prev[u] != NIL) next[prev[u]] = next[u];
        else head[d] = next[u];
        if (next[u] != NIL) prev[next[u]] = prev[u];
        next[u] = NIL;
        prev[u] = NIL;
        filed[u] = false;
    }
    bool empty(std::size_t d) const { return head[d] == NIL; }
};

// Print a quotient graph: adjacency, incidence, cliques, in the order the
// structure holds them.
void mmd1Show(const Graph& A, const Graph& I, const Cliques& C,
             const std::vector<std::size_t>& degrees, const std::string& title = "",
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
        std::cout << "  " << std::setw(width) << u << ": {" << adjacencyText.str()
                  << "} {" << incidenceText.str() << "} degree " << degrees[u] << "\n";
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
    std::cout << "\n";
}

// Print the state arrays: degrees, buckets, min degree, members, eliminated,
// and the order so far.
void mmd1ShowState(const std::vector<std::size_t>& degrees, const Buckets& buckets,
                  std::size_t minDegree,
                  const std::vector<std::vector<std::int32_t>>& superMembers,
                  const std::vector<bool>& eliminated,
                  const std::vector<std::int32_t>& pivots, const std::string& title = "") {
    const std::size_t n = superMembers.size();
    int width = static_cast<int>(std::to_string(n > 0 ? n - 1 : 0).size());
    if (!title.empty()) std::cout << title << "\n";
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        std::string status;
        if (!eliminated[u]) status = "live";
        else if (!superMembers[u].empty()) status = "done";
        else status = "merged";
        std::ostringstream superMemberList;
        bool first = true;
        for (std::int32_t v : superMembers[u]) {
            superMemberList << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::cout << "  " << std::setw(width) << u << ": members ["
                  << superMemberList.str() << "] " << status << "\n";
    }
    std::ostringstream degreesText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        degreesText << (u == 0 ? "" : " ") << std::setw(width) << degrees[u];
    std::ostringstream superMembersText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        superMembersText << (u == 0 ? "" : " ") << "[";
        bool firstMember = true;
        for (std::int32_t v : superMembers[u]) {
            superMembersText << (firstMember ? "" : " ") << v;
            firstMember = false;
        }
        superMembersText << "]";
    }
    std::ostringstream eliminatedText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        eliminatedText << (u == 0 ? "" : " ") << std::setw(width) << (eliminated[u] ? 1 : 0);
    std::ostringstream bucketsText;
    bool firstBucket = true;
    for (std::size_t d = 0; d < buckets.head.size(); ++d) {
        if (buckets.empty(d)) continue;
        bucketsText << (firstBucket ? "" : "  ") << d << ": [";
        bool firstMember = true;
        for (std::int32_t v = buckets.head[d]; v != NIL; v = buckets.next[v]) {
            bucketsText << (firstMember ? "" : " ") << v;
            firstMember = false;
        }
        bucketsText << "]";
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
    for (std::int32_t pivot : pivots)
        for (std::int32_t u : superMembers[pivot]) {
            std::cout << (firstOrder ? "" : ", ") << u;
            firstOrder = false;
        }
    std::cout << "]\n\n";
}

// Entries actually stored. Each edge costs two, one per endpoint in A. Each
// incidence costs two as well, the clique id in I and the member in C. Watch
// the total fall monotonically; the naive graph's only rises.
std::size_t mmd1Storage(const Graph& A, const Graph& I, const Cliques& C) {
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
std::vector<std::int32_t> mmd1Neighbors(const Graph& A, const Graph& I, const Cliques& C,
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
    for (std::int32_t v : A[u])
        if (mark[v] != tag) { mark[v] = tag; neighbors.push_back(v); }
    for (std::int32_t c : I[u])
        for (std::int32_t v : C.at(c))
            if (mark[v] != tag) { mark[v] = tag; neighbors.push_back(v); }
    return neighbors;
}

// Turn the pivot into a clique, then merge in every member it makes
// indistinguishable. Identical to md5Eliminate: this layer changes how often
// degrees are refreshed, not what an elimination does.
//
// Returns (neighbors, absorbedCliques, prunedEdges, mergedVertices): as in md2,
// plus the vertices folded into the pivot by mass elimination. The last three
// are reported for display; only neighbors is used by the caller.//
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
// Mass elimination adds two more lines, and breaks the identity in the first one:
// from here C[pivot] is reach(pivot) minus what the pivot absorbed.
//
//     merged   = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
//     C[pivot] = C[pivot] - merged
std::tuple<std::vector<std::int32_t>, std::vector<std::int32_t>,
           std::vector<std::pair<std::int32_t, std::int32_t>>, std::vector<std::int32_t>>
mmd1Eliminate(Graph& A, Graph& I, Cliques& C, std::vector<bool>& eliminated,
             std::vector<std::int32_t>& mark, std::int32_t& tag, std::int32_t pivot) {
    const std::vector<std::int32_t> neighbors = mmd1Neighbors(A, I, C, mark, tag, pivot);
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

    // Mass elimination. u is INDISTINGUISHABLE from the pivot when the two have
    // the same closed neighborhood, mmd1Neighbors(u) | {u} == mmd1Neighbors(pivot)
    // | {pivot}, as it stood before the step. Equivalently, now that the clique is
    // formed, when everything u can still reach lies inside it. The test below is
    // a cheap sufficient condition for that: nothing explicit left and no clique
    // but the new one means u sees exactly what the pivot sees, so eliminating it
    // next would cost no fill. Fold it into the pivot now and strip it from the
    // cliques, since it is no longer a vertex.
    // merged = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
    std::vector<std::int32_t> mergedVertices;
    for (std::int32_t u : neighbors) {
        if (A[u].empty() && I[u].size() == 1 && I[u][0] == pivot) {
            I[u].clear();
            eliminated[u] = true;
            mergedVertices.push_back(u);
        }
    }
    if (!mergedVertices.empty()) {           // C[pivot] - merged, one compaction pass
        ++tag;
        for (std::int32_t u : mergedVertices) mark[u] = tag;
        kept.clear();
        for (std::int32_t v : C[pivot])
            if (mark[v] != tag) kept.push_back(v);
        C[pivot].swap(kept);
    }

    A[pivot].clear();
    I[pivot].clear();
    eliminated[pivot] = true;
    return {neighbors, absorbedCliques, prunedEdges, mergedVertices};
}

// Move u from the bucket for its old degree to the one for newDegree. Filing
// pushes at the head, which is the O(1) end of the list.
//
// Set view: buckets[d] is the set of live vertices whose current degree is d, and
// filing, unfiling and refiling are add, discard and move between two of them. A
// linked list gives all three in O(1) and gives the head in O(1) too, which is
// everything the picker asks of it. What it does not give is a minimum, which is
// why minDegree walks. A sorted container would hand over the minimum directly and
// charge a log on every file, and files outnumber picks.
void mmd1Refile(Buckets& buckets, std::vector<std::size_t>& degrees,
               std::int32_t u, std::size_t newDegree) {
    buckets.unfile(degrees[u], u);
    degrees[u] = newDegree;
    buckets.file(newDegree, u);
}

// Multiple elimination: a batch of independent pivots per degree refresh.
//
// delta widens the batch to vertices within delta of the minimum degree, which
// buys still fewer refreshes for a real concession, since those vertices are not
// minimal. delta = 0 keeps the batch to true minima. A negative delta takes one
// pivot per round, which is md5's behavior reached through this code path.
// delta is signed: negative means one pivot per round. It is compared against a
// degree and its useful range stops at n - 1, so it is an index-like quantity by
// Oblio's rule, a std::int32_t rather than a count.
std::vector<std::int32_t> mmd1MinimumDegree(const Graph& G, std::int32_t delta = 0) {
    const std::size_t n = G.size();
    std::size_t nnzTrilA = 0;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    Graph A = G;                                  // explicit vertex neighbors
    Graph I(n);                                   // cliques that contain each vertex
    Cliques C(n);      // clique id -> member list
    std::vector<std::int32_t> mark(n, NIL);       // scratch for membership, with tag
    std::int32_t tag = 0;
    std::vector<std::vector<std::int32_t>> superMembers(n);   // for the expansion
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        superMembers[u].push_back(u);
    std::vector<bool> eliminated(n, false);
    std::vector<std::int32_t> pivots;             // the order over supervariables
    std::size_t numEliminated = 0;                // a counter, not a scan of eliminated
    std::size_t nnzL = 0;

    std::vector<std::size_t> degrees(n);          // a degree counts, so it measures
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) degrees[u] = A[u].size();
    std::size_t numDegreeComputations = n;

    Buckets buckets(n);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        buckets.file(degrees[u], u);
    std::size_t minDegree = n > 0 ? *std::min_element(degrees.begin(), degrees.end()) : 0;
    std::size_t numBucketProbes = 0;
    std::size_t numRounds = 0;                    // batches, the metric this layer adds
    std::vector<std::int32_t> touchedRound(n, NIL);  // the round u was last evicted in

    // NOT PRODUCTION: display only. The trace is what makes these files teachable and
    // is the whole reason they exist; nothing downstream reads it.
    mmd1Show(A, I, C, degrees, "start: every edge explicit, no clique yet", &eliminated);
    mmd1ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
    while (numEliminated < n) {
        while (buckets.empty(minDegree)) {      // walk up to the first live bucket
            ++minDegree;
            ++numBucketProbes;
        }
        ++numBucketProbes;

        // ---- one BATCH, no degree refreshed inside it ----------------------
        // Take pivots from buckets [minDegree, minDegree + delta]. Eviction is
        // what keeps them independent: eliminating a pivot pulls every vertex it
        // reached out of the buckets, so whatever is still filed was not reached,
        // hence is not adjacent to anything taken this round.
        //
        // Set view of the invariant the eviction maintains, where reached is the
        // union of C[p] over the pivots taken so far:
        //
        //     filed = live - reached,  so  batch & reached == {}
        //
        // No set is built for either side. Membership in filed is the filed flag,
        // and touchedRound is the same idea one level up: it stamps the round a
        // vertex was evicted in, so the refresh set is accumulated without a set
        // and without a sort.
        // Clamped: a degree is at most n - 1, so a wider window would walk the
        // bucket array off its end.
        std::size_t batchLimit = minDegree;      // delta > 0 here, so no narrowing
        if (delta > 0)
            batchLimit = std::min(minDegree + delta, n - 1);
        std::vector<std::int32_t> batch;
        std::vector<std::int32_t> touched;    // first-touch order, no set and no sort
        while (true) {
            if (buckets.empty(minDegree)) {     // this degree is drained
                if (minDegree >= batchLimit) break;
                ++minDegree;
                ++numBucketProbes;
                continue;
            }
            std::int32_t pivot = buckets.head[minDegree];
            std::size_t degree = degrees[pivot];
            buckets.unfile(degree, pivot);

            auto [neighbors, absorbedCliques, prunedEdges, mergedVertices] =
                mmd1Eliminate(A, I, C, eliminated, mark, tag, pivot);
            batch.push_back(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + mergedVertices.size();
            for (std::int32_t u : mergedVertices) {   // the pivot now stands for them too
                superMembers[pivot].insert(superMembers[pivot].end(),
                                           superMembers[u].begin(), superMembers[u].end());
                superMembers[u].clear();
                buckets.unfile(degrees[u], u);
                degrees[u] = 0;
            }
            degrees[pivot] = 0;

            for (std::int32_t u : C[pivot]) {   // EVICT, with a stale degree
                buckets.unfile(degrees[u], u);
                if (touchedRound[u] != static_cast<std::int32_t>(numRounds)) {
                    touchedRound[u] = static_cast<std::int32_t>(numRounds);
                    touched.push_back(u);       // a marker, so O(1) per eviction
                }
            }

            // A supervariable of size w is w consecutive columns of L. Its
            // external degree is what remains of the clique after the merges,
            // since a merged vertex joins the supervariable instead of
            // neighboring it. The first column holds ext + w - 1 entries below
            // its diagonal, the next ext + w - 2, down to ext, and each column
            // contributes its own diagonal.
            std::size_t superSize = superMembers[pivot].size();
            std::size_t externalDegree = C[pivot].size();
            nnzL += superSize * externalDegree + superSize * (superSize - 1) / 2 + superSize;

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
            std::ostringstream mergedVerticesText;
            if (mergedVertices.empty()) {
                mergedVerticesText << "none";
            } else {
                bool first = true;
                for (std::int32_t u : mergedVertices) {
                    mergedVerticesText << (first ? "" : ", ") << u;
                    first = false;
                }
            }
            std::ostringstream evictedText;
            if (C[pivot].empty()) {
                evictedText << "none";
            } else {
                bool first = true;
                for (std::int32_t u : C[pivot]) {
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
            if (delta < 0) break;               // one pivot per round, as md5 does
        }

        // ---- one REFRESH, for everything the batch reached -----------------
        std::vector<std::int32_t> refreshedVertices;
        for (std::int32_t u : touched) if (!eliminated[u]) refreshedVertices.push_back(u);
        for (std::int32_t u : refreshedVertices) {
            degrees[u] = mmd1Neighbors(A, I, C, mark, tag, u).size();
            buckets.file(degrees[u], u);
        }
        numDegreeComputations += refreshedVertices.size();
        for (std::int32_t u : refreshedVertices) minDegree = std::min(minDegree, degrees[u]);
        ++numRounds;

        std::ostringstream batchText;
        for (std::size_t k = 0; k < batch.size(); ++k)
            batchText << (k == 0 ? "" : ", ") << batch[k];
        std::ostringstream refreshedVerticesText;
        if (refreshedVertices.empty()) {
            refreshedVerticesText << "none";
        } else {
            bool first = true;
            for (std::int32_t u : refreshedVertices) {
                refreshedVerticesText << (first ? "" : ", ") << u;
                first = false;
            }
        }
        std::ostringstream title;
        title << "round " << (numRounds - 1) << " done: batch of " << batch.size() << ": "
              << batchText.str() << ", refreshed: " << refreshedVerticesText.str();
        // NOT PRODUCTION: display only. The trace is what makes these files teachable and
        // is the whole reason they exist; nothing downstream reads it.
        mmd1Show(A, I, C, degrees, title.str(), &eliminated);
        mmd1ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
    }

    std::vector<std::int32_t> order;
    for (std::int32_t pivot : pivots)
        for (std::int32_t u : superMembers[pivot]) order.push_back(u);
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

// A square grid graph, four-neighbor, for running the counters at a size the seven examples
// cannot reach. It is here rather than among them because it is not an example: nothing about it
// illustrates a mechanism, and its trace is far too long to read. The grid mode below discards the
// trace and prints only the closing counter lines, which is what a comparison between layers wants.
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

// Keep the trace's summary lines and discard everything else, as it is written rather than
// afterwards. A grid trace is far too large to hold: every step prints the whole quotient graph,
// so at n = 10000 the captured text runs to gigabytes and the process dies holding it. This
// filters line by line instead, keeping only what the whitelist names, so the memory is one line.
class CounterSink : public std::streambuf {
public:
    explicit CounterSink(std::vector<std::string> keys) : mKeys(std::move(keys)) {}

protected:
    int overflow(int ch) override {
        if (ch == traits_type::eof()) return traits_type::not_eof(ch);
        const char c = static_cast<char>(ch);
        if (c != '\n') { mLine.push_back(c); return ch; }
        for (const std::string& key : mKeys)
            if (mLine.rfind(key, 0) == 0) { mKept.push_back(mLine); break; }
        mLine.clear();
        return ch;
    }

public:
    const std::vector<std::string>& kept() const { return mKept; }

private:
    std::vector<std::string> mKeys;
    std::vector<std::string> mKept;
    std::string              mLine;
};
void run(const std::string& name, const Graph& G) {
    std::cout << "=== " << name << " ===\n";
    mmd1MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Grid mode: one square grid, the trace discarded, the counters kept.
    //
    //   ./mmd1_cpp grid 22
    if (argc > 2 && std::string(argv[1]) == "grid") {
        const int side = std::atoi(argv[2]);
        std::cout << "=== grid " << side << "x" << side << " (n = " << side * side << ") ===\n";
        CounterSink sink({"order:", "nnz(L)", "degree computations"});
        std::streambuf* saved = std::cout.rdbuf(&sink);
        mmd1MinimumDegree(gridGraph(side));
        std::cout.rdbuf(saved);
        for (const std::string& line : sink.kept()) std::cout << line << "\n";
        return 0;
    }

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
    // which mmd1's merge test declines a genuine supervariable. At the step whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test mmd1Neighbors(A, I, C, u) contained in
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
    // 4 already has 1 as a child when 0 merges into it. The merge happens at step
    // 2 of 5, so the run continues afterwards and the selection degree, 3 over
    // {2, 3, 4}, differs from the external degree, 2 over {2, 3}, with the
    // difference being the size of what merged. And superMembers ends with a hole
    // in the middle, slot 4 empty between two used ones, while no pivot equals
    // its own step number. See the README sections on mass elimination and on
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

    // All of them by default. To run just one, pass its number: ./mmd1_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
