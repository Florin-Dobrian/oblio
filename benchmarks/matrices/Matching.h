#pragma once

// Matching.h -- maximum bipartite matching, and the structural rank it gives us.
//
// WHY THIS IS HERE. A matrix is SINGULAR if a column carries no nonzero value, which is the test
// this folder started with, and that test finds only the isolated vertex. The complete structural
// statement is a maximum matching between rows and columns: a permutation contributing to the
// determinant needs a nonzero in every position, which is exactly a perfect matching in the
// bipartite graph of the pattern. Where none exists, every term of the determinant expansion
// vanishes and the matrix is singular for almost every choice of values. The size of the maximum
// matching is the STRUCTURAL RANK, and n minus it is what the driver labels.
//
// A NOTE ON WHICH MATCHING, because the obvious reading is wrong and cost us a wrong proposal.
// The bipartite matching is not the graph's own matching. It was proposed here that a connected
// component of odd size with no diagonal must be singular, since its vertices cannot be paired
// off, and that is false: the permutation need not be an involution. The triangle refutes it,
//
//     A = [[0,1,1],[1,0,1],[1,1,0]]     det 2, eigenvalues -1, -1, 2
//
// which is odd, has no diagonal, and is nonsingular, its determinant coming from the 3-cycle.
// General (non-bipartite) matching answers the pairing question and is the wrong tool here; the
// bipartite graph of rows against columns is bipartite by construction whatever A's symmetry.
//
// WHERE THE CODE CAME FROM. The matching itself is taken from the `combinatorial-suite` package,
// file `hk_it_csr.cpp`: an iterative Hopcroft-Karp over CSR adjacency, O(E sqrt(V)), with the
// lean single-level BFS and a stack-based DFS carrying an edge index so no phase rescans. It is
// carried here with THREE MECHANICAL CHANGES and no change to any algorithm:
//
//   1. the `main` and its file-reading harness are dropped, this being a header rather than a
//      program, and the folder having drivers of its own;
//   2. the six definitions are marked `inline` and the two static ones follow, a header being
//      included rather than compiled;
//   3. everything sits in `namespace Matching`, which also keeps its `NIL` off `Oblio::NIL`.
//
// It therefore keeps its own conventions, bare `size_t` and `int32_t` where the house rule says
// `std::`-qualified, and its own naming. That is deliberate: it is imported code with a life of
// its own upstream, and the cost of re-syncing it is what a local cleanup would buy against.
//
// WHAT IS NOT BUILT ON IT, and what the rest of Dulmage-Mendelsohn is actually for. DM has two
// halves and only the first is about singularity. Its FINE decomposition comes straight off the
// maximum matching and splits rows and columns into an underdetermined part, a square part and an
// overdetermined part, both outer parts being empty exactly when a perfect matching exists. Its
// COARSE decomposition then takes the strongly connected components of the digraph the matching
// induces on the square part, which gives the block triangular form: that half is REDUCIBILITY,
// and it is the part that is distinctively DM rather than distinctively matching.
//
// Neither is built here, and the second would buy little if it were. Our matrices are SYMMETRIC,
// and a symmetric permutation cannot produce a nontrivial triangular form, so the block
// triangular form of a symmetric matrix is block DIAGONAL and its blocks are the connected
// components of the graph, which one depth-first pass gives. DM proper earns its keep on
// unsymmetric matrices, where the blocks really are triangular and the order among them matters.

#include "oblio/SparseMatrix.h"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Matching {

constexpr int32_t NIL = -1;

/* ---------- Input: BipartiteGraph ---------- */

struct BipartiteGraph {
    size_t sNumVtxs;
    size_t tNumVtxs;
    size_t numEdgs;
    std::vector<size_t>  sIdx;
    std::vector<size_t>  tIdx;
    std::vector<int32_t> sAdj;
    std::vector<int32_t> tAdj;
};

/* ---------- Output: BipartiteMatching ---------- */

struct BipartiteMatching {
    size_t sNumVtxs;
    size_t tNumVtxs;
    size_t numEdgs;
    std::vector<int32_t> sMate;
    std::vector<int32_t> tMate;
};

/* ---------- State: HKIState ---------- */

struct HKIState {
    std::vector<int32_t> sLevel;       // length sNumVtxs+1; sLevel[sNumVtxs] is the NIL sentinel
    std::vector<int32_t> sIdx;       // length sNumVtxs; relative offset within s's adjacency range, persistent within a phase
    std::vector<int32_t> sPrcbStk;       // length sNumVtxs; DFS stack (s-vertices)
    std::vector<int32_t> tPrcbStk;       // length sNumVtxs; t-vertex chosen at each depth
    int32_t stkTop;
};

/* ---------- BipartiteGraph construction ---------- */

inline BipartiteGraph buildBipartiteGraph(size_t sNumVtxs, size_t tNumVtxs,
                                   const std::vector<std::pair<int32_t,int32_t>>& edges) {
    std::vector<std::vector<int32_t>> sTmp(sNumVtxs);
    std::vector<std::vector<int32_t>> tTmp(tNumVtxs);
    for (auto& e : edges) {
        int32_t s = e.first, t = e.second;
        if (s >= 0 && static_cast<size_t>(s) < sNumVtxs && t >= 0 && static_cast<size_t>(t) < tNumVtxs) {
            sTmp[s].push_back(t);
            tTmp[t].push_back(s);
        }
    }
    for (size_t s = 0; s < sNumVtxs; s++) {
        std::sort(sTmp[s].begin(), sTmp[s].end());
        sTmp[s].erase(std::unique(sTmp[s].begin(), sTmp[s].end()), sTmp[s].end());
    }
    for (size_t t = 0; t < tNumVtxs; t++) {
        std::sort(tTmp[t].begin(), tTmp[t].end());
        tTmp[t].erase(std::unique(tTmp[t].begin(), tTmp[t].end()), tTmp[t].end());
    }

    BipartiteGraph graph;
    graph.sNumVtxs = sNumVtxs;
    graph.tNumVtxs = tNumVtxs;

    graph.sIdx.assign(sNumVtxs + 1, 0);
    for (size_t s = 0; s < sNumVtxs; s++)
        graph.sIdx[s + 1] = graph.sIdx[s] + sTmp[s].size();
    graph.sAdj.resize(graph.sIdx[sNumVtxs]);
    for (size_t s = 0; s < sNumVtxs; s++)
        std::copy(sTmp[s].begin(), sTmp[s].end(), graph.sAdj.begin() + graph.sIdx[s]);

    graph.tIdx.assign(tNumVtxs + 1, 0);
    for (size_t t = 0; t < tNumVtxs; t++)
        graph.tIdx[t + 1] = graph.tIdx[t] + tTmp[t].size();
    graph.tAdj.resize(graph.tIdx[tNumVtxs]);
    for (size_t t = 0; t < tNumVtxs; t++)
        std::copy(tTmp[t].begin(), tTmp[t].end(), graph.tAdj.begin() + graph.tIdx[t]);

    graph.numEdgs = graph.sIdx[sNumVtxs];
    return graph;
}

/* ---------- BipartiteMatching construction ---------- */

inline BipartiteMatching emptyBipartiteMatching(const BipartiteGraph& graph) {
    BipartiteMatching matching;
    matching.sNumVtxs = graph.sNumVtxs;
    matching.tNumVtxs = graph.tNumVtxs;
    matching.numEdgs = 0;
    matching.sMate.assign(graph.sNumVtxs, NIL);
    matching.tMate.assign(graph.tNumVtxs, NIL);
    return matching;
}

/* ---------- Greedy initial matching: simple ---------- */

inline int32_t greedyInit(const BipartiteGraph& graph, BipartiteMatching& matching) {
    int32_t numEdgs = 0;
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] != NIL) continue;
        size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
        for (size_t k = sBegin; k < sEnd; k++) {
            int32_t t = graph.sAdj[k];
            if (matching.tMate[t] == NIL) {
                matching.sMate[s] = t;
                matching.tMate[t] = static_cast<int32_t>(s);
                numEdgs++;
                break;
            }
        }
    }
    matching.numEdgs += numEdgs;
    return numEdgs;
}

/* ---------- Greedy initial matching: min-degree ---------- */

inline int32_t greedyInitMd(const BipartiteGraph& graph, BipartiteMatching& matching) {
    int32_t numEdgs = 0;
    std::vector<int32_t> deg(graph.tNumVtxs, 0);
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
        for (size_t k = sBegin; k < sEnd; k++) deg[graph.sAdj[k]]++;
    }
    std::vector<int32_t> sOrder(graph.sNumVtxs);
    for (size_t s = 0; s < graph.sNumVtxs; s++) sOrder[s] = static_cast<int32_t>(s);
    /* Sort s-vertices in increasing order of degree, breaking ties by vertex label. */
    std::sort(sOrder.begin(), sOrder.end(), [&](int32_t s1, int32_t s2){
        size_t s1Deg = graph.sIdx[s1 + 1] - graph.sIdx[s1];
        size_t s2Deg = graph.sIdx[s2 + 1] - graph.sIdx[s2];
        return s1Deg < s2Deg || (s1Deg == s2Deg && s1 < s2);
    });
    for (int32_t s : sOrder) {
        if (matching.sMate[s] != NIL) continue;
        int32_t best = NIL, bestDeg = INT_MAX;
        size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
        for (size_t k = sBegin; k < sEnd; k++) {
            int32_t t = graph.sAdj[k];
            if (matching.tMate[t] == NIL && deg[t] < bestDeg) {
                best = t;
                bestDeg = deg[t];
            }
        }
        if (best != NIL) {
            matching.sMate[s] = best;
            matching.tMate[best] = s;
            numEdgs++;
        }
    }
    matching.numEdgs += numEdgs;
    return numEdgs;
}

/* ---------- HK BFS ---------- */

inline bool bfs(const BipartiteGraph& graph, const BipartiteMatching& matching, HKIState& state) {
    std::vector<int32_t> sPrcbQue(graph.sNumVtxs);
    int32_t queHead = 0, queTail = 0;

    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] == NIL) { state.sLevel[s] = 0; sPrcbQue[queTail++] = static_cast<int32_t>(s); }
        else state.sLevel[s] = INT_MAX;
    }
    state.sLevel[graph.sNumVtxs] = INT_MAX;  /* NIL sentinel */

    while (queHead < queTail) {
        int32_t s = sPrcbQue[queHead++];
        if (state.sLevel[s] < state.sLevel[graph.sNumVtxs]) {
            size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
            for (size_t k = sBegin; k < sEnd; k++) {
                int32_t t = graph.sAdj[k];
                int32_t ss = (matching.tMate[t] == NIL) ? static_cast<int32_t>(graph.sNumVtxs) : matching.tMate[t];
                if (state.sLevel[ss] == INT_MAX) {
                    state.sLevel[ss] = state.sLevel[s] + 1;
                    if (matching.tMate[t] != NIL) sPrcbQue[queTail++] = matching.tMate[t];
                }
            }
        }
    }
    return state.sLevel[graph.sNumVtxs] != INT_MAX;
}

/*
 * DFS: iterative with edge index.
 *
 * state.sIdx[s] is an offset WITHIN s's adjacency range [graph.sIdx[s], graph.sIdx[s+1]).
 * So the "current candidate edge" is graph.sAdj[graph.sIdx[s] + state.sIdx[s]].
 */
inline bool dfs(int32_t sFirst, const BipartiteGraph& graph,
                BipartiteMatching& matching, HKIState& state) {
    state.stkTop = 0;
    state.sPrcbStk[state.stkTop++] = sFirst;

    while (state.stkTop > 0) {
        int32_t s = state.sPrcbStk[state.stkTop - 1];
        size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
        int32_t sNumEdgs = static_cast<int32_t>(sEnd - sBegin);

        bool pushed = false;
        while (state.sIdx[s] < sNumEdgs) {
            int32_t t = graph.sAdj[sBegin + state.sIdx[s]];
            int32_t ss = (matching.tMate[t] == NIL) ? static_cast<int32_t>(graph.sNumVtxs) : matching.tMate[t];
            if (state.sLevel[ss] != state.sLevel[s] + 1) {
                state.sIdx[s]++;
                continue;
            }

            state.tPrcbStk[state.stkTop - 1] = t;
            state.sIdx[s]++;

            if (matching.tMate[t] == NIL) {
                /* Found augmenting path — augment all the way back */
                for (int32_t k = state.stkTop - 1; k >= 0; k--) {
                    matching.tMate[state.tPrcbStk[k]] = state.sPrcbStk[k];
                    matching.sMate[state.sPrcbStk[k]] = state.tPrcbStk[k];
                }
                return true;
            }

            state.sPrcbStk[state.stkTop++] = matching.tMate[t];
            pushed = true;
            break;
        }

        if (!pushed) {
            state.sLevel[s] = INT_MAX;
            state.stkTop--;
        }
    }
    return false;
}

/* ---------- Top-level Hopcroft-Karp Iterative ---------- */

inline int32_t hkIterativeMcm(const BipartiteGraph& graph, BipartiteMatching& matching) {
    HKIState state;
    state.sLevel.assign(graph.sNumVtxs + 1, 0);
    state.sIdx.assign(graph.sNumVtxs, 0);
    state.sPrcbStk.assign(graph.sNumVtxs, 0);
    state.tPrcbStk.assign(graph.sNumVtxs, 0);
    state.stkTop = 0;

    int32_t numPhases = 0;
    int32_t newEdgs = 0;
    while (bfs(graph, matching, state)) {
        numPhases++;
        for (size_t s = 0; s < graph.sNumVtxs; s++) state.sIdx[s] = 0;
        for (size_t s = 0; s < graph.sNumVtxs; s++) {
            if (matching.sMate[s] == NIL && dfs(static_cast<int32_t>(s), graph, matching, state)) {
                newEdgs++;
            }
        }
    }
    matching.numEdgs += newEdgs;
    return numPhases;
}

/* ---------- Validation ---------- */

inline void validateBipartiteMatching(const BipartiteGraph& graph, const BipartiteMatching& matching) {
    int32_t errors = 0;
    int32_t matchedS = 0, matchedT = 0;

    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] != NIL) {
            matchedS++;
            int32_t t = matching.sMate[s];
            if (t < 0 || static_cast<size_t>(t) >= graph.tNumVtxs) {
                fprintf(stderr, "ERROR: sMate[%zu] = %d out of range\n", s, t);
                errors++;
            } else if (matching.tMate[t] != static_cast<int32_t>(s)) {
                fprintf(stderr, "ERROR: sMate[%zu]=%d but tMate[%d]=%d\n", s, t, t, matching.tMate[t]);
                errors++;
            } else {
                size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
                if (!std::binary_search(graph.sAdj.begin() + sBegin,
                                        graph.sAdj.begin() + sEnd, t)) {
                    fprintf(stderr, "ERROR: edge (%zu,%d) not in graph\n", s, t);
                    errors++;
                }
            }
        }
    }
    for (size_t t = 0; t < graph.tNumVtxs; t++) {
        if (matching.tMate[t] != NIL) matchedT++;
    }

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %zu\n", matching.numEdgs);
    printf("S matched: %d, T matched: %d\n", matchedS, matchedT);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

} // namespace Matching

namespace MatrixBenchmark {

// The adapter, and the whole of what this folder calls.
//
// TWO THINGS IT DOES, and the second is the one that makes it correct.
//
// It bypasses `buildBipartiteGraph`, which takes an edge list and sorts, deduplicates and
// rebuilds it into CSR. A `SparseMatrix` already IS that: `colPtr` and `rowIdx`, sorted ascending
// within a column with no duplicates, which the reader guarantees. So the s-side arrays are the
// matrix's own, and building an edge list first would cost half a million pairs on the larger
// matrices here for nothing.
//
// And IT FILTERS ON THE VALUE, not on the pattern. Our conversion inserts a structural zero on
// every diagonal that nothing landed on, so the pattern carries entries whose value is zero.
// Handing that pattern to the matching would match every one of them and report full structural
// rank on a matrix that is nothing but isolated vertices. That filter is the difference between
// 1424 and 1589 on netscience.
//
// The t-side arrays are left empty. `bfs` and `dfs` read `sIdx`, `sAdj` and `tMate` and never
// touch `tIdx` or `tAdj`, so building them would be work with no reader. `validateBipartiteMatching`
// does read `sAdj`, and works.
inline Matching::BipartiteGraph bipartiteGraph(const Oblio::SparseMatrix<double>& A) {
    Matching::BipartiteGraph graph;
    graph.sNumVtxs = A.size();   // s is a column of A
    graph.tNumVtxs = A.size();   // t is a row of A

    graph.sIdx.assign(A.size() + 1, 0);
    graph.sAdj.reserve(A.nnz());

    for (std::size_t j = 0; j < A.size(); ++j) {
        for (std::size_t cp = A.colPtr()[j]; cp < A.colPtr()[j + 1]; ++cp)
            if (A.val()[cp] != 0.0)
                graph.sAdj.push_back(A.rowIdx()[cp]);
        graph.sIdx[j + 1] = graph.sAdj.size();
    }

    graph.numEdgs = graph.sAdj.size();
    return graph;
}

// The structural rank: the size of a maximum matching between columns and rows. Equal to n when
// the matrix is structurally nonsingular, and below it by the structural deficiency otherwise.
//
// The min-degree greedy runs first because it is nearly free and it is worth several phases: on
// netscience it starts the matching at 1318 of the final 1424 and takes Hopcroft-Karp from six
// phases to two.
inline std::size_t structuralRank(const Oblio::SparseMatrix<double>& A) {
    const Matching::BipartiteGraph graph = bipartiteGraph(A);
    Matching::BipartiteMatching matching = Matching::emptyBipartiteMatching(graph);
    Matching::greedyInitMd(graph, matching);
    Matching::hkIterativeMcm(graph, matching);
    return matching.numEdgs;
}

} // namespace MatrixBenchmark
