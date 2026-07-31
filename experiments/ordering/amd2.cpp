// amd2.cpp -- approximate minimum degree, complete.
//
// amd1 has the idea: the degree bound, computed once per clique and read once per
// vertex. What it does not have is the rest of what amd_1 and amd_2 do, and this
// file adds it, one pass at a time. All seven are in. Section 5.13 of
// archive/sparse_factorization.md, plus the vendored routine itself in
// vendored/vendored_amd.cpp.
//
// The list:
//
//   1. AGGRESSIVE ABSORPTION, which kills a clique the moment the bound work
//      shows it lies inside the new one                                   [done]
//   2. HASH SUPERVARIABLE DETECTION, which finds vertices indistinguishable
//      from EACH OTHER rather than from the pivot                         [done]
//   3. the two-pass degree update, scan 1 obtaining |C[c] - C[p]| by
//      subtraction from a maintained clique degree                        [done]
//   4. dense row and column detection by the alpha ratio, held out and
//      placed last                                                        [done]
//   5. amd_aat and amd_preprocess, forming the pattern of A + A'          [done]
//   6. amd_postorder, so the output is a postorder of the assembly tree    [done]
//   7. amd_valid and the Control/Info interface                           [done]
//
// PASSES 1 AND 2, THE TWO MECHANISMS THAT RIDE ALONG WITH THE BOUND. Neither is
// about the degree, and both are cheap only because the bound's work has already
// been done.
//
// AGGRESSIVE ABSORPTION. |C[c] - C[p]| has just been computed for every clique c
// that the new one touched. If it is zero, C[c] lies entirely inside C[p], so that
// clique is dead and can be absorbed at once. Ordinary absorption only kills the
// cliques the PIVOT touched; this kills cliques that any reached vertex touched,
// and it costs nothing extra because the quantity was needed anyway.
//
// HASH SUPERVARIABLE DETECTION. Mass elimination merges a vertex into the pivot.
// Two vertices can be indistinguishable from EACH OTHER without either being
// absorbable into the pivot, and no pivot test can see that. AMD hashes (A[u],
// I[u]),
// compares only within a hash bucket, and merges on an exact match. The hash is a
// filter, never the decision. MMD reaches the same vertices through mmdupd's q2h
// list, so the goal is shared and the mechanism is not.
//
// Build:  g++ -std=c++17 -O3 amd2.cpp -o amd2_cpp  (or: make)
// Run:    ./amd2_cpp
//         ./amd2_cpp 3      just the third example

#include <cstdlib>
#include <algorithm>
#include <cmath>
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

std::vector<std::int32_t> amd2Neighbors(const Graph& A, const Graph& I, const Cliques& C,
                                       const std::vector<bool>& eliminated,
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
void amd2Show(const Graph& A, const Graph& I, const Cliques& C,
              const std::vector<std::size_t>& degrees,
              const std::vector<std::size_t>& exact, const std::string& title = "",
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
                  << "} {" << incidenceText.str() << "} bound " << degrees[u]
                  << " exact " << exact[u] << "\n";
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
void amd2ShowState(const std::vector<std::size_t>& degrees, const Buckets& buckets,
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
std::size_t amd2Storage(const Graph& A, const Graph& I, const Cliques& C) {
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
// The three statuses amd_valid can return, and the two Control knobs.
const int AMD_OK = 0;
const int AMD_OK_BUT_JUMBLED = 1;
const int AMD_INVALID = -2;

// What amd_valid checks before anything else touches the input: the column
// pointers start at zero and ascend, and every row index is in range. It is not a
// yes or no. Unsorted columns and duplicate entries return a THIRD answer,
// AMD_OK_BUT_JUMBLED, which is not an error but a request: it is what tells
// amd_order to run amd_preprocess rather than use the pattern directly.
//
// So pass 7 and pass 5 are one mechanism seen from two sides. The validity check
// decides whether the conditioning pass is needed, and the conditioning pass exists
// because the check tolerates what it tolerates.
//
// For Oblio the whole of this is a constructor invariant. A SparseMatrix is sorted,
// duplicate free and in range before anything is asked of it, so there is nothing
// here to decide at ordering time.
int amd2Valid(std::size_t n, const std::vector<std::int32_t>& Ap,
              const std::vector<std::int32_t>& Ai) {
    if (Ap.size() != n + 1 || Ap[0] != 0 || Ap[n] < 0) return AMD_INVALID;
    if (Ai.size() != static_cast<std::size_t>(Ap[n])) return AMD_INVALID;
    int result = AMD_OK;
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(n); ++j) {
        if (Ap[j] > Ap[j + 1]) return AMD_INVALID;      // pointers must ascend
        std::int32_t last = NIL;
        for (std::int32_t p = Ap[j]; p < Ap[j + 1]; ++p) {
            const std::int32_t i = Ai[p];
            if (i < 0 || i >= static_cast<std::int32_t>(n)) return AMD_INVALID;
            if (i <= last) result = AMD_OK_BUT_JUMBLED;  // unsorted, or a duplicate
            last = i;
        }
    }
    return result;
}

// R, the row form of the pattern of A with duplicates removed.
//
// Amd.cpp calls this when the input may be unsorted or hold duplicates, since
// A + A' can be formed from R without either problem. R is the pattern of A
// transposed, so R + R' is A + A' and nothing is lost by working from it.
//
// The diagonal is NOT dropped here; amd2Aat deals with it, because it has to count
// it anyway.
//
// Set view: R[i] = { j : A(i,j) != 0 }, and flag[i] == j is the membership test
// that makes the deduplication one comparison rather than a search.
void amd2Preprocess(std::size_t n, const std::vector<std::int32_t>& Ap,
                    const std::vector<std::int32_t>& Ai,
                    std::vector<std::int32_t>& Rp, std::vector<std::int32_t>& Ri) {
    std::vector<std::int32_t> counts(n, 0), flag(n, NIL);
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(n); ++j)
        for (std::int32_t p = Ap[j]; p < Ap[j + 1]; ++p) {
            const std::int32_t i = Ai[p];
            if (flag[i] != j) {                 // i has not appeared in column j yet
                ++counts[i];
                flag[i] = j;
            }
        }
    Rp.assign(n + 1, 0);
    for (std::size_t i = 0; i < n; ++i) Rp[i + 1] = Rp[i] + counts[i];
    std::vector<std::int32_t> position(Rp.begin(), Rp.begin() + n);
    std::fill(flag.begin(), flag.end(), NIL);
    Ri.assign(static_cast<std::size_t>(Rp[n]), 0);
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(n); ++j)
        for (std::int32_t p = Ap[j]; p < Ap[j + 1]; ++p) {
            const std::int32_t i = Ai[p];
            if (flag[i] != j) {
                Ri[position[i]++] = j;
                flag[i] = j;
            }
        }
}

// The pattern of A + A' with the diagonal dropped, as adjacency lists, plus the
// statistics Amd.cpp reports about the input. (Ap, Ai) is the column form and
// (Rp, Ri) its row form, both free of duplicates, as amd2Preprocess leaves them.
//
// The ordering is defined on a symmetric structure and a general matrix is not
// one, so this is where an arbitrary pattern becomes a graph. Two things go: the
// diagonal, which is a self loop and says nothing about fill, and the distinction
// between A(i,j) and A(j,i), since either one forces the same elimination.
//
// The symmetry reported is the vendored definition, with B the strictly triangular
// parts of A:
//
//     sym = nnz(B & B') / nnz(B),   or 1 when nnz(B) is zero
//
// Amd.cpp computes it with a two-pointer scan that walks the two triangles
// together, which needs sorted columns and saves a pass. Ours asks the row form
// whether the transposed entry exists, one stamp and one comparison, which is what
// the rest of this file does and works on unsorted input. That is a deviation from
// the vendored code, and it is stated rather than hidden.
struct AatInfo {
    std::size_t nz;
    std::size_t numDiagonal;
    std::size_t numBoth;
    double symmetry;
    std::size_t nzaat;
};

Graph amd2Aat(std::size_t n, const std::vector<std::int32_t>& Ap,
              const std::vector<std::int32_t>& Ai,
              const std::vector<std::int32_t>& Rp,
              const std::vector<std::int32_t>& Ri, AatInfo& info) {
    const std::size_t nz = static_cast<std::size_t>(Ap[n]);
    std::vector<std::int32_t> rowMark(n, NIL);  // rowMark[k] == j: A(j,k) present
    std::size_t numDiagonal = 0;
    std::size_t numBoth = 0;
    Graph A(n);
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(n); ++j) {
        for (std::int32_t p = Rp[j]; p < Rp[j + 1]; ++p)   // row j, the transposed
            rowMark[Ri[p]] = j;
        for (std::int32_t p = Ap[j]; p < Ap[j + 1]; ++p) {
            const std::int32_t i = Ai[p];
            if (i == j) {                       // a self loop, dropped
                ++numDiagonal;
                continue;
            }
            if (i > j && rowMark[i] == j) ++numBoth;   // below and above both there
            A[j].push_back(i);                  // both directions, deduplicated below
            A[i].push_back(j);
        }
    }

    std::vector<std::int32_t> stamp(n, NIL);    // one pass per list, as everywhere
    std::vector<std::int32_t> kept;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        kept.clear();
        for (std::int32_t v : A[u])
            if (stamp[v] != u) {
                stamp[v] = u;
                kept.push_back(v);
            }
        A[u].swap(kept);
    }

    std::size_t nzaat = 0;
    for (std::size_t u = 0; u < n; ++u) nzaat += A[u].size();
    info = {nz, numDiagonal, numBoth,
            (nz == numDiagonal) ? 1.0
                                : (2.0 * static_cast<double>(numBoth))
                                      / static_cast<double>(nz - numDiagonal),
            nzaat};
    return A;
}


std::vector<std::int32_t> amd2Neighbors(const Graph& A, const Graph& I, const Cliques& C,
                                       const std::vector<bool>& eliminated,
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
        if (mark[v] != tag && !eliminated[v]) { mark[v] = tag; neighbors.push_back(v); }
    for (std::int32_t c : I[u])
        for (std::int32_t v : C.at(c))
            if (mark[v] != tag && !eliminated[v]) { mark[v] = tag; neighbors.push_back(v); }
    return neighbors;
}

// The degree md5 would have computed: the union of A[u] with the members of every
// clique in I[u], counted in original vertices. Kept only so the trace can show
// the bound beside the truth and count how often the bound is loose.
//
// Set view: sum of |superMembers[v]| over v in reach(u). It is the union the bound
// exists to avoid, so this function is instrumentation and nothing more.
std::size_t amd2ExactDegree(const Graph& A, const Graph& I, const Cliques& C,
                            const std::vector<bool>& eliminated,
                            const std::vector<std::vector<std::int32_t>>& superMembers,
                            std::vector<std::int32_t>& mark, std::int32_t& tag,
                            std::int32_t u) {
    std::size_t degree = 0;
    for (std::int32_t v : amd2Neighbors(A, I, C, eliminated, mark, tag, u))
        degree += superMembers[v].size();
    return degree;
}

// Turn the pivot into a clique, then merge in every member it makes
// indistinguishable. Identical to md5Eliminate: this layer changes how a degree is
// estimated afterwards, not what an elimination does.
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
amd2Eliminate(Graph& A, Graph& I, Cliques& C, std::vector<bool>& eliminated,
             std::vector<std::int32_t>& mark, std::int32_t& tag, std::int32_t pivot) {
    const std::vector<std::int32_t> neighbors =
        amd2Neighbors(A, I, C, eliminated, mark, tag, pivot);
    const std::vector<std::int32_t> absorbedCliques = I[pivot];
    for (std::int32_t c : absorbedCliques)
        C.erase(c);
    C.create(pivot, neighbors);     // becomes L_pivot, the column pattern

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
    // the same closed neighborhood, amd2Neighbors(u) | {u} == amd2Neighbors(pivot)
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
void amd2Refile(Buckets& buckets, std::vector<std::size_t>& degrees,
               std::int32_t u, std::size_t newDegree) {
    buckets.unfile(degrees[u], u);
    degrees[u] = newDegree;
    buckets.file(newDegree, u);
}

// amd1 plus the mechanisms that ride along with the bound: aggressive absorption
// and hash supervariable detection.
// amd1 plus the passes of amd_1 and amd_2, one at a time. See the header.
//
// alpha and aggressive are amd_order's Control array. alpha is the dense row and
// column ratio, AMD_DEFAULT_DENSE, with a negative value removing only completely
// dense rows. aggressive switches aggressive absorption off, which is the only
// other thing the vendored routine lets a caller change.
std::vector<std::int32_t> amd2MinimumDegree(const Graph& G, double alpha = 10.0,
                                            bool aggressive = true) {
    const std::size_t n = G.size();
    std::size_t nnzTrilA = 0;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    Graph A = G;                                  // explicit vertex neighbors
    Graph I(n);                                   // cliques that contain each vertex
    Cliques C(n);      // clique id -> member list
    // Twice n, because the hash test stamps cliques at c + n so a clique and a
    // vertex of the same index cannot be confused. Everything else uses the first n.
    std::vector<std::int32_t> mark(2 * n, NIL);       // scratch for membership, with tag
    std::int32_t tag = 0;
    std::vector<std::vector<std::int32_t>> superMembers(n);   // for the expansion
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        superMembers[u].push_back(u);
    std::vector<bool> eliminated(n, false);
    std::vector<std::int32_t> pivots;             // the order over supervariables
    std::size_t numEliminated = 0;
    // Live ORIGINAL vertices, which is not n - numEliminated. numEliminated counts
    // what has left the SELECTION, and a hash merge folds v into a LIVE u, so v
    // stops being selectable while the vertices it stands for are still live inside
    // u. The first cap of the bound needs the second reading, so it gets its own
    // counter: only an elimination and a dense removal reduce it.
    std::size_t numLive = n;                // a counter, not a scan of eliminated
    std::size_t nnzL = 0;

    // The cache, and the count of degree computations, which is what this layer
    // exists to reduce. Built once, then touched only where it can be wrong.
    // The cache, as in md5, except that from the first elimination it holds a
    // BOUND rather than a degree. exact[] is carried alongside for the trace only.
    std::vector<std::size_t> degrees(n);          // a degree counts, so it measures
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) degrees[u] = A[u].size();
    std::vector<std::size_t> exact = degrees;
    std::size_t numDegreeComputations = n;
    std::size_t numMemberVisits = 0;              // what an exact refresh would cost
    std::size_t numCliqueReads = 0;               // clique reads in the bound itself
    std::size_t numIncidenceReads = 0;            // incidence entries scan 1 walks
    std::size_t numBoundChecks = 0;
    std::size_t numLooseBounds = 0;
    std::size_t numAbsorbed = 0;                  // cliques killed aggressively
    // |C[c]| for every live clique, weighted. Exact, not an estimate, and the
    // invariants that keep it so are worth stating because they are not obvious:
    // a live clique never holds an eliminated vertex, since eliminating v absorbs
    // every clique in I[v]; mass elimination only ever removes a vertex whose I[u]
    // is {pivot}, so no other clique is touched; and a hash merge folds v into u
    // where I[u] == I[v], so every clique holding v holds u and its weighted size
    // does not move.
    std::vector<std::size_t> cliqueDegree(n, 0);
    Graph hashBucket(n + 1);                      // Amd.cpp's Head[hval], reused
    // The assembly tree, for the postorder. parent[c] is the clique that absorbed
    // c, and frontSize[e] is the front e would form, its pivots plus what it
    // reaches. Amd.cpp keeps the same two in Pe and Elen.
    std::vector<std::int32_t> parent(n, NIL);
    std::vector<std::size_t> frontSize(n, 0);
    std::vector<bool> isElement(n, false);
    std::size_t numHashMerges = 0;                // pairs found by the hash
    std::size_t numDivides = 0;                   // AMD_NDIV, one per off-diagonal
    std::size_t numMultsubsLdl = 0;               // AMD_NMULTSUBS_LDL
    std::size_t numMultsubsLu = 0;                // AMD_NMULTSUBS_LU
    std::size_t frontMax = 0;                     // AMD_DMAX, the largest front

    // The buckets, and minDegree, a LOWER BOUND on the current minimum degree.
    // The search starts at minDegree rather than at 0, so it never looks at
    // buckets known to be empty. The bound may lag, and the walk corrects it; what
    // it must never do is overshoot, since a vertex below it would never be seen.
    //
    // n buckets is exactly right. A live vertex counts only live neighbors, so its
    // degree is at most n - 1, and the walk stops at the first non-empty bucket,
    // which exists while anything is live.
    Buckets buckets(n);

    // ---- DENSE ROWS AND COLUMNS, before anything is filed -------------------
    // A row whose INITIAL degree exceeds the threshold is taken out of the graph
    // and placed last. One dense row touches nearly everything, so it inflates the
    // degree of nearly everything, and the ordering spends its effort avoiding a
    // vertex it cannot avoid. Removing it is cheaper than ordering around it.
    //
    // The threshold is the vendored one, dense = alpha * sqrt(n), floored at 16 and
    // capped at n, with a negative alpha meaning n - 2, which removes only rows
    // that are completely dense. It is computed from the INITIAL degrees and never
    // revisited: a vertex that becomes dense later is not caught.
    std::size_t denseThreshold = (alpha < 0)
        ? (n >= 2 ? n - 2 : 0)
        : static_cast<std::size_t>(alpha * std::sqrt(static_cast<double>(n)));
    denseThreshold = std::max<std::size_t>(16, denseThreshold);
    denseThreshold = std::min<std::size_t>(n, denseThreshold);
    std::vector<std::int32_t> denseVertices;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        if (degrees[u] > denseThreshold) denseVertices.push_back(u);

    // Amd.cpp sets Nv[i] = 0 and Pe[i] = EMPTY, which leaves the vertex in other
    // lists while contributing zero weight to every degree. Clearing its own lists
    // and purging it from the rest is the same statement in our representation,
    // since a weight of zero and an absence are indistinguishable to every count
    // this file makes.
    for (std::int32_t u : denseVertices) {
        A[u].clear();
        I[u].clear();
        superMembers[u].clear();               // weight 0, as Nv[i] = 0
        eliminated[u] = true;
        ++numEliminated;
        --numLive;                             // out of the graph, not merged
    }
    if (!denseVertices.empty()) {
        std::vector<std::int32_t> keptAdjacency;
        for (std::int32_t w = 0; w < static_cast<std::int32_t>(n); ++w) {
            if (eliminated[w]) continue;
            keptAdjacency.clear();
            for (std::int32_t v : A[w])
                if (!eliminated[v]) keptAdjacency.push_back(v);
            A[w].swap(keptAdjacency);
            degrees[w] = A[w].size();
            exact[w] = degrees[w];
        }
    }

    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        if (!eliminated[u]) buckets.file(degrees[u], u);
    std::size_t minDegree = 0;
    bool haveMin = false;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        if (eliminated[u]) continue;
        if (!haveMin || degrees[u] < minDegree) { minDegree = degrees[u]; haveMin = true; }
    }
    std::size_t numBucketProbes = 0;

    amd2Show(A, I, C, degrees, exact,
             "start: every edge explicit, no clique yet, degrees exact", &eliminated);
    amd2ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
    int step = 0;
    while (numEliminated < n) {
        while (buckets.empty(minDegree)) {      // walk up to the first live bucket
            ++minDegree;
            ++numBucketProbes;
        }
        ++numBucketProbes;
        std::int32_t pivot = buckets.head[minDegree];   // whatever was filed last
        auto [neighbors, absorbedCliques, prunedEdges, mergedVertices] =
            amd2Eliminate(A, I, C, eliminated, mark, tag, pivot);
        isElement[pivot] = true;
        for (std::int32_t c : absorbedCliques) parent[c] = pivot;   // its children
        std::size_t degree = neighbors.size();
        pivots.push_back(pivot);
        numEliminated += 1 + mergedVertices.size();
        for (std::int32_t u : mergedVertices) {   // the pivot now stands for them too
            superMembers[pivot].insert(superMembers[pivot].end(),
                                       superMembers[u].begin(), superMembers[u].end());
            superMembers[u].clear();
        }
        numLive -= superMembers[pivot].size();  // every original the pivot stands for

        buckets.unfile(degrees[pivot], pivot);  // the pivot has left the graph
        degrees[pivot] = 0;
        for (std::int32_t u : mergedVertices) { // and so have the merged vertices
            buckets.unfile(degrees[u], u);
            degrees[u] = 0;
        }

        // ---- the BOUND, in place of md5's exact refresh --------------------
        // C[pivot] is the new clique. Everything it reached needs a new degree, and
        // the bound replaces the union with a sum of separate contributions.
        //
        // Set view of the three quantities, none of which is built as a set:
        //
        //     pivotClique = C[pivot]
        //     degme       = |C[pivot]|             weighted, original vertices
        //     outside[c]  = |C[c] - C[pivot]|      ONE value per CLIQUE
        //
        // The last line is the whole idea. |C[c] - C[pivot]| depends on c and not
        // on the vertex asking, so it is computed once and read many times, where
        // the exact degree recomputes a union per vertex. mark[v] == inClique is
        // the membership test for C[pivot]; a second tag makes touchedCliques a set
        // too, so a clique is listed once however many vertices reach it.
        const std::vector<std::int32_t> pivotClique = C[pivot];
        ++tag;
        const std::int32_t inClique = tag;      // membership of C[pivot], one test
        for (std::int32_t v : pivotClique) mark[v] = inClique;
        std::size_t degme = 0;
        for (std::int32_t v : pivotClique) degme += superMembers[v].size();
        cliqueDegree[pivot] = degme;            // what scan 1 subtracts from

        // ---- SCAN 1, |C[c] - C[pivot]| ONCE PER CLIQUE ---------------------
        // This is the whole reason the bound is cheap: the quantity depends on c
        // alone, so every vertex whose incidence list holds c reads it rather than
        // recomputing it.
        //
        // And it is obtained by SUBTRACTION, never by looking at C[c] at all:
        //
        //     |C[c] - C[pivot]| = |C[c]| - sum of weight(u) over u in C[c] & C[pivot]
        //
        // cliqueDegree[c] supplies the first term, and the members of C[pivot]
        // supply the second, since c is in I[u] exactly when u is in C[c]. So the
        // scan walks the INCIDENCE lists of the new clique's members and pays
        // sum |I[u]|, where amd1 walked the member lists of every touched clique
        // and paid sum |C[c]|. Amd.cpp's scan 1 does the same thing at
        // `we = Degree[e] + wnvi` then `we -= nvi`.
        std::vector<std::int32_t> touchedCliques;
        std::vector<std::size_t> outside(n, 0);
        ++tag;
        const std::int32_t seenClique = tag;
        for (std::int32_t u : pivotClique) {
            const std::size_t weightU = superMembers[u].size();
            for (std::int32_t c : I[u]) {
                if (c == pivot) continue;
                if (mark[c] != seenClique) {    // first sighting: start from |C[c]|
                    mark[c] = seenClique;
                    touchedCliques.push_back(c);
                    outside[c] = cliqueDegree[c] - weightU;
                } else {                        // every later member just subtracts
                    outside[c] -= weightU;
                }
                ++numIncidenceReads;
                numMemberVisits += C[c].size();
            }
        }

        // AGGRESSIVE ABSORPTION. Set view: dead = { c : C[c] <= C[pivot] }, the
        // containment decided by the count already computed for the bound, since
        // |C[c] - C[pivot]| == 0 IS C[c] <= C[pivot]. Then I[u] = I[u] - dead for
        // every u in C[pivot], one stamp and one compaction pass, the same shape as
        // the absorption in the eliminator. Ordinary absorption kills only what the
        // pivot touched; this kills what any reached vertex touched, and it is free
        // because the quantity was computed for the bound anyway.
        std::vector<std::int32_t> deadCliques;
        if (aggressive)
            for (std::int32_t c : touchedCliques)
                if (outside[c] == 0) deadCliques.push_back(c);
        if (!deadCliques.empty()) {
            ++tag;
            const std::int32_t deadTag = tag;
            for (std::int32_t c : deadCliques) { C.erase(c); mark[c] = deadTag; }
            std::vector<std::int32_t> keptCliques;
            for (std::int32_t u : pivotClique) {
                keptCliques.clear();
                for (std::int32_t c : I[u])
                    if (mark[c] != deadTag) keptCliques.push_back(c);
                I[u].swap(keptCliques);         // I[u] - dead
            }
            for (std::int32_t c : deadCliques) parent[c] = pivot;  // same tree
            numAbsorbed += deadCliques.size();
        }

        const std::size_t numLeft = numLive;
        const std::vector<std::int32_t>& refreshedVertices = pivotClique;
        for (std::int32_t u : refreshedVertices) {
            // bound = |A[u]| + |C[pivot] - {u}| + sum |C[c] - C[pivot]| over the
            // cliques in I[u] - {pivot}, against the exact
            // |( A[u] | C[c] for c in I[u] ) - {u}|. The bound replaces the union
            // by a sum, so an overlap outside C[pivot] is counted once per clique
            // that holds it, which is exactly where it overcounts.
            std::size_t explicitPart = 0;
            for (std::int32_t v : A[u]) explicitPart += superMembers[v].size();
            std::size_t bound = explicitPart + degme - superMembers[u].size();
            for (std::int32_t c : I[u]) {
                if (c == pivot) continue;
                bound += outside[c];
                ++numCliqueReads;               // what the bound pays instead
            }
            bound = std::min(bound, numLeft - superMembers[u].size());
            bound = std::min(bound, degrees[u] + degme - superMembers[u].size());
            exact[u] = amd2ExactDegree(A, I, C, eliminated, superMembers, mark, tag, u);
            ++numBoundChecks;
            if (bound > exact[u]) ++numLooseBounds;
            amd2Refile(buckets, degrees, u, bound);
        }
        // HASH SUPERVARIABLE DETECTION. Vertices indistinguishable from EACH
        // OTHER, which the pivot test cannot see. Hash first so the exact
        // comparison runs only within a bucket; the hash is a filter, never the
        // decision.
        // The buckets are an array indexed by the hash value, allocated once and
        // cleared only where it was used, which is Amd.cpp's Head[hval]. A map keyed
        // by the hash would cost a log per insertion and a node per group, for a
        // quantity that is already an index into 0 .. n.
        std::vector<std::size_t> usedKeys;
        for (std::int32_t u : pivotClique) {
            if (eliminated[u]) continue;
            // The hash stands for the PAIR of sets (A[u], I[u]), so equal sets
            // always collide and unequal ones usually do not. It is a filter and
            // never the decision: a collision costs a comparison, not a wrong merge.
            //
            // A SUM, because addition has no order and neither do the sets. Sorting
            // to build a key would be a log factor for nothing; Amd.cpp sums the
            // indices and reduces modulo n for the same reason.
            std::size_t key = 0;
            for (std::int32_t v : A[u])
                if (!eliminated[v]) key += static_cast<std::size_t>(v) + 1;
            for (std::int32_t c : I[u])
                key += (static_cast<std::size_t>(c) + 1) * (n + 1);
            const std::size_t k = key % (n + 1);
            if (hashBucket[k].empty()) usedKeys.push_back(k);
            hashBucket[k].push_back(u);
        }
        std::vector<std::pair<std::int32_t, std::int32_t>> hashPairs;
        for (std::size_t k : usedKeys) {
            std::vector<std::int32_t>& group = hashBucket[k];
            if (group.size() < 2) continue;
            for (std::size_t x = 0; x < group.size(); ++x) {
                std::int32_t u = group[x];
                if (eliminated[u]) continue;
                for (std::size_t y = x + 1; y < group.size(); ++y) {
                    std::int32_t v = group[y];
                    if (eliminated[v]) continue;
                    // The exact test, which the hash only filters for:
                    //     A[u] - {v} == A[v] - {u}  and  I[u] == I[v]
                    // Decided by stamping one side and counting matches on the
                    // other, as every other membership test in this file is, so it
                    // costs one pass and no sort.
                    ++tag;
                    const std::int32_t other = tag;
                    std::size_t sizeV = 0;
                    for (std::int32_t w : A[v])
                        if (w != u && !eliminated[w]) { mark[w] = other; ++sizeV; }
                    for (std::int32_t c : I[v]) {        // stamped past the vertices
                        mark[c + static_cast<std::int32_t>(n)] = other;
                        ++sizeV;
                    }
                    std::size_t sizeU = 0;
                    bool same = true;
                    for (std::int32_t w : A[u]) {
                        if (w == v || eliminated[w]) continue;
                        ++sizeU;
                        if (mark[w] != other) { same = false; break; }
                    }
                    if (same)
                        for (std::int32_t c : I[u]) {
                            ++sizeU;
                            if (mark[c + static_cast<std::int32_t>(n)] != other) {
                                same = false;
                                break;
                            }
                        }
                    if (!same || sizeU != sizeV) continue;
                    // v is folded into u and left exactly where it lies, with a
                    // weight of zero, which is Amd.cpp's Nv[v] = 0 and the same move
                    // the dense prepass makes. Removing it from every clique and
                    // every adjacency would cost a pass over the whole structure per
                    // merge, against O(1) for this.
                    //
                    // Nothing is lost by leaving it. The merge required
                    // A[u] - {v} == A[v] - {u}, so every list holding v holds u as
                    // well: v is redundant wherever it appears, never the only way
                    // to reach anything. The walks skip it.
                    superMembers[u].insert(superMembers[u].end(),
                                           superMembers[v].begin(), superMembers[v].end());
                    superMembers[v].clear();
                    buckets.unfile(degrees[v], v);
                    A[v].clear();
                    I[v].clear();
                    eliminated[v] = true;
                    ++numEliminated;
                    hashPairs.push_back({u, v});
                    ++numHashMerges;
                }
            }
        }

        for (std::size_t k : usedKeys) hashBucket[k].clear();  // only what was used

        numDegreeComputations += refreshedVertices.size();
        for (std::int32_t u : refreshedVertices)
            if (!eliminated[u]) minDegree = std::min(minDegree, degrees[u]);

        // A supervariable of size w is w consecutive columns of L. Its external
        // degree is what remains of the clique after the merges, since a merged
        // vertex joins the supervariable instead of neighboring it, and every
        // member left there is a live vertex standing for itself alone. The first
        // column then holds ext + w - 1 entries below its diagonal, the next
        // ext + w - 2, down to ext, and each column contributes its own diagonal.
        std::size_t superSize = superMembers[pivot].size();
        std::size_t externalDegree = 0;
        for (std::int32_t v : C[pivot])
            if (!eliminated[v]) externalDegree += superMembers[v].size();
        // The dense rows were taken out but they still sit below every column of
        // L, so each one adds an entry to each. Amd.cpp does the same at
        // r = degme + ndense, which makes nnz(L) an upper bound rather than a count
        // once anything is dense: a dense row is assumed nonzero everywhere.
        const std::size_t reachSize = externalDegree + denseVertices.size();
        nnzL += superSize * reachSize + superSize * (superSize - 1) / 2 + superSize;
        frontSize[pivot] = superSize + externalDegree;   // Amd.cpp: nvpiv + degme
        // Amd.cpp's Info, with f the pivots of this front and r what it reaches.
        {
            const std::size_t f = superSize;
            const std::size_t r = reachSize;
            const std::size_t lnzMe = f * r + (f - 1) * f / 2;
            numDivides += lnzMe;
            const std::size_t multsubs = f * r * r + r * (f - 1) * f
                                       + (f - 1) * f * (2 * f - 1) / 6;
            numMultsubsLu += multsubs;
            numMultsubsLdl += (multsubs + lnzMe) / 2;
            frontMax = std::max(frontMax, f + r);
        }

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
        title << "step " << step << ": eliminate " << pivot << " (degree " << degree
              << ", size " << superSize << ", external degree " << externalDegree
              << "), absorbed cliques: " << absorbedCliquesText.str()
              << ", pruned edges: " << prunedEdgesText.str()
              << ", merged vertices: " << mergedVerticesText.str()
              << ", refreshed: " << refreshedVerticesText.str();
        amd2Show(A, I, C, degrees, exact, title.str(), &eliminated);
        amd2ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
        ++step;
    }

    // ---- THE POSTORDER, in place of raw elimination order -------------------
    // The elements form the assembly tree, and any postorder of it gives the same
    // factor: a node is numbered after all its descendants either way, so no fill
    // moves. What a postorder buys is locality. Children finish before their parent
    // starts, so the update from a child is consumed while it is still warm, and a
    // supernode's columns come out contiguous instead of scattered.
    //
    // Two details from Amd.cpp, both about which child goes first. The child lists
    // are built by walking the elements DOWNWARD, so a list comes out ascending.
    // Then the BIGGEST child by front size is moved to the end, so the largest
    // subtree is traversed last and the stack of pending updates stays small.
    std::vector<std::int32_t> child(n, NIL), sibling(n, NIL);
    for (std::int32_t e = static_cast<std::int32_t>(n) - 1; e >= 0; --e)
        if (isElement[e] && parent[e] != NIL) {   // downward, so lists end ascending
            sibling[e] = child[parent[e]];
            child[parent[e]] = e;
        }
    for (std::int32_t e = 0; e < static_cast<std::int32_t>(n); ++e) {
        if (!isElement[e] || child[e] == NIL) continue;
        std::int32_t biggest = NIL, biggestPrevious = NIL, previous = NIL;
        std::size_t largest = 0;
        bool haveBiggest = false;
        for (std::int32_t f = child[e]; f != NIL; f = sibling[f]) {
            if (!haveBiggest || frontSize[f] >= largest) {  // the LAST maximal one
                largest = frontSize[f];
                biggestPrevious = previous;
                biggest = f;
                haveBiggest = true;
            }
            previous = f;
        }
        if (sibling[biggest] != NIL) {            // already last means nothing to do
            if (biggestPrevious == NIL) child[e] = sibling[biggest];
            else sibling[biggestPrevious] = sibling[biggest];
            sibling[biggest] = NIL;
            sibling[previous] = biggest;
        }
    }

    std::vector<std::int32_t> postorder, stack;
    for (std::int32_t root = 0; root < static_cast<std::int32_t>(n); ++root) {
        if (!isElement[root] || parent[root] != NIL) continue;   // roots in order
        stack.push_back(root);
        while (!stack.empty()) {                  // explicit stack, no recursion
            const std::int32_t e = stack.back();
            if (child[e] != NIL) {
                const std::int32_t f = child[e];
                child[e] = sibling[f];            // each child pushed exactly once
                stack.push_back(f);
            } else {
                postorder.push_back(e);           // all descendants done, number it
                stack.pop_back();
            }
        }
    }

    std::vector<std::int32_t> order;
    for (std::int32_t e : postorder)
        for (std::int32_t u : superMembers[e]) order.push_back(u);
    // The dense rows go last, in index order, and Amd.cpp counts the block they
    // form as completely full, which it is in the worst case and usually is not.
    if (!denseVertices.empty()) {
        const std::size_t numDense = denseVertices.size();
        order.insert(order.end(), denseVertices.begin(), denseVertices.end());
        nnzL += numDense * (numDense - 1) / 2 + numDense;
        const std::size_t lnzDense = numDense * (numDense - 1) / 2;  // counted full
        numDivides += lnzDense;
        const std::size_t denseMultsubs =
            (numDense - 1) * numDense * (2 * numDense - 1) / 6;
        numMultsubsLu += denseMultsubs;
        numMultsubsLdl += (denseMultsubs + lnzDense) / 2;
        frontMax = std::max(frontMax, numDense);
    }

    std::cout << "nnz(L) = " << nnzL << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "degree computations: " << numDegreeComputations
              << ", bucket probes: " << numBucketProbes << "\n";
    std::cout << "clique-member visits an exact degree would need: "
              << numMemberVisits << "\n";
    std::cout << "clique reads the bound needed:                    "
              << numCliqueReads << "\n";
    std::cout << "incidence entries scan 1 walked:                  "
              << numIncidenceReads << "\n";
    std::cout << "bound was loose " << numLooseBounds << " times out of "
              << numBoundChecks << "\n";
    std::cout << "aggressively absorbed: " << numAbsorbed
              << ", hash merges: " << numHashMerges << "\n";
    std::cout << "dense threshold: " << denseThreshold << ", dense rows removed: "
              << denseVertices.size() << "\n";
    // The rest of Amd.cpp's Info array, which is a factorization cost PREDICTION
    // and so belongs to a symbolic phase rather than to an ordering. It is here for
    // the record: nnz(L) above is AMD_LNZ, computed from the same expression, and
    // the three below come from the same f and r with no extra structure.
    //
    // Two fields are deliberately missing. AMD_MEMORY and AMD_NCMPA report the peak
    // workspace and the number of garbage collections in the vendored flat pool,
    // and this file has no pool to compact, so there is nothing to report.
    std::cout << "predicted: divides " << numDivides << ", multiply-subtracts LDL "
              << numMultsubsLdl << ", LU " << numMultsubsLu << ", largest front "
              << frontMax << "\n";
    std::cout << "order: [";
    for (std::size_t k = 0; k < order.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << order[k];
    std::cout << "]\n";
    return order;
}

void run(const std::string& name, const Graph& G) {
    std::cout << "=== " << name << " ===\n";
    amd2MinimumDegree(G);
    std::cout << "\n";
}

// Order a general sparse matrix given in column form, which is what amd_order
// takes. The pattern may be unsorted, may hold duplicates, may carry a diagonal
// and need not be symmetric. Two preprocess passes give the deduplicated row and
// column forms, amd2Aat turns them into a graph, and the ordering runs on that.
//
// alpha and aggressive are the whole of amd_order's Control array, AMD_DENSE and
// AMD_AGGRESSIVE. Everything else it takes is the matrix.
std::vector<std::int32_t> amd2OrderMatrix(std::size_t n,
                                          const std::vector<std::int32_t>& Ap,
                                          const std::vector<std::int32_t>& Ai,
                                          double alpha = 10.0,
                                          bool aggressive = true) {
    const int status = amd2Valid(n, Ap, Ai);
    if (status == AMD_INVALID) {
        std::cout << "status: invalid input\n";
        return {};
    }
    std::cout << "status: " << (status ? "ok but jumbled" : "ok") << "\n";
    std::vector<std::int32_t> Rp, Ri, Cp, Ci;
    amd2Preprocess(n, Ap, Ai, Rp, Ri);          // row form of A
    amd2Preprocess(n, Rp, Ri, Cp, Ci);          // and back, so the column form is clean
    AatInfo info;
    Graph A = amd2Aat(n, Cp, Ci, Rp, Ri, info);
    // The one sort, at construction, exactly as the Graph path sorts its input:
    // outside every loop, and only so the two twins hold the same order.
    for (std::size_t u = 0; u < n; ++u) std::sort(A[u].begin(), A[u].end());
    std::cout << "n = " << n << ", nz = " << info.nz << ", symmetry = "
              << std::fixed << std::setprecision(3) << info.symmetry
              << std::defaultfloat << ", nz diagonal = " << info.numDiagonal
              << ", nz(A+A\') = " << info.nzaat << "\n";
    return amd2MinimumDegree(A, alpha, aggressive);
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
    // belongs to two cliques that overlap outside the new one, which needs enough
    // eliminations to have made several cliques and enough fill for them to
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
    // which amd2's merge test declines a genuine supervariable. At the step whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test amd2Neighbors(A, I, C, u) contained in
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

    // matrix1, the one example given as a MATRIX rather than as a graph, so the
    // input path has something to chew on. Six by six in column form, and
    // deliberately awful: unsymmetric, with a diagonal, with duplicate entries, and
    // with one column whose rows are out of order. amd_order accepts exactly this,
    // and amd2Preprocess and amd2Aat are what turn it into a graph.
    //
    //   column 0: rows 0 1 3          column 3: rows 3 0 2   (unsorted)
    //   column 1: rows 1 2 2          column 4: rows 4 5
    //   column 2: rows 0 2 5          column 5: rows 1 5
    const std::size_t matrix1N = 6;
    const std::vector<std::int32_t> matrix1Ap = {0, 3, 6, 9, 12, 14, 16};
    const std::vector<std::int32_t> matrix1Ai = {0, 1, 3,
                                                 1, 2, 2,
                                                 0, 2, 5,
                                                 3, 0, 2,
                                                 4, 5,
                                                 1, 5};

    // All of them by default. To run just one, pass its number: ./amd2_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    if (selected == 0 || selected == static_cast<int>(examples.size()) + 1) {
        std::cout << "=== matrix1 ===\n";
        amd2OrderMatrix(matrix1N, matrix1Ap, matrix1Ai);
        std::cout << "\n";
    }
    return 0;
}
