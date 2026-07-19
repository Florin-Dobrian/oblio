#pragma once
// Shared test utilities.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/Types.h"   // NIL

#include <vector>
#include <set>
#include <utility>
#include <cstddef>
#include <cstdint>

namespace OblioTest {

// Structural symmetry: every stored (row i, col j) has a matching (row j, col i).
// All of Oblio's structural code assumes A is symmetric, and A is stored fully — so
// test inputs must satisfy this. There is no runtime guardrail in the library yet;
// this catches asymmetric *test* matrices early instead of letting them silently
// produce a wrong ordering or elimination forest.
template<class Val>
bool isStructurallySymmetric(const Oblio::SparseMatrix<Val>& A) {
    const auto& colPtr = A.colPtr();
    const auto& rowIdx = A.rowIdx();
    const std::size_t size = A.size();
    std::set<std::pair<std::size_t, std::size_t>> entries;
    for (std::size_t j = 0; j < size; ++j)
        for (std::size_t p = colPtr[j]; p < colPtr[j + 1]; ++p)
            entries.insert({static_cast<std::size_t>(rowIdx[p]), j});
    for (const auto& e : entries)
        if (entries.find({e.second, e.first}) == entries.end())
            return false;
    return true;
}

// Dense oracle for the pattern of L, independent of any of Oblio's structural code.
// A is permuted into factor order, then Cholesky fill is simulated densely:
// eliminating column j makes every pair of its subdiagonal rows adjacent, which is
// the definition of fill. Returns, for each factor column, its nonzero rows in
// increasing order (the diagonal first, then the subdiagonal rows).
//
// This is deliberately the naive O(n^3) formulation. It shares no code path with
// the elimination forest or the symbolic factorization, so agreement between the
// two is real evidence, not a tautology.
template<class Val>
std::vector<std::vector<std::int32_t>>
denseFactorPattern(const Oblio::SparseMatrix<Val>& A, const Oblio::Permutation& p) {
    const std::size_t size = A.size();
    const auto& colPtr = A.colPtr();
    const auto& rowIdx = A.rowIdx();
    const auto& oldToNew = p.oldToNew();

    // Lower triangle of the permuted A.
    std::vector<std::vector<bool>> L(size, std::vector<bool>(size, false));
    for (std::size_t ac = 0; ac < size; ++ac) {
        const std::size_t lc = static_cast<std::size_t>(oldToNew[ac]);
        for (std::size_t t = colPtr[ac]; t < colPtr[ac + 1]; ++t) {
            const std::size_t lr =
                static_cast<std::size_t>(oldToNew[static_cast<std::size_t>(rowIdx[t])]);
            if (lr >= lc)
                L[lr][lc] = true;
        }
    }

    // Fill: eliminating column j connects its subdiagonal rows pairwise.
    for (std::size_t j = 0; j < size; ++j) {
        std::vector<std::size_t> sub;
        for (std::size_t i = j + 1; i < size; ++i)
            if (L[i][j])
                sub.push_back(i);
        for (std::size_t a = 0; a < sub.size(); ++a)
            for (std::size_t b = a + 1; b < sub.size(); ++b)
                L[sub[b]][sub[a]] = true;
    }

    std::vector<std::vector<std::int32_t>> pattern(size);
    for (std::size_t j = 0; j < size; ++j)
        for (std::size_t i = j; i < size; ++i)
            if (L[i][j])
                pattern[j].push_back(static_cast<std::int32_t>(i));
    return pattern;
}

// Independent oracle for the fundamental supernodes, derived from the dense factor
// pattern alone. Recomputes the nodal elimination tree from that pattern (the parent of
// column j is its first subdiagonal nonzero), then merges column j into its parent k
// exactly when k has one child and the two columns share a sparsity pattern, which for
// the nodal case means |pattern(j)| == |pattern(k)| + 1.
//
// Returns the column -> supernode map, labeled by the same rule the engine uses: scan
// columns in increasing order, a column joining its child keeps that child's label, a
// column starting a supernode takes the next free label. Sharing no code path with the
// engine, this checks the merge decisions and the labelling independently.
template<class Val>
std::vector<std::int32_t>
fundamentalSupernodes(const Oblio::SparseMatrix<Val>& A, const Oblio::Permutation& p,
                      std::size_t& snodeSizeOut) {
    const auto pattern = denseFactorPattern(A, p);
    const std::size_t size = A.size();

    // The nodal elimination tree, and how many children each column has.
    std::vector<std::int32_t> parent(size, Oblio::NIL);
    std::vector<std::size_t>  numChildren(size, 0);
    for (std::size_t j = 0; j < size; ++j) {
        // pattern[j][0] is the diagonal; the next entry, if any, is the parent.
        if (pattern[j].size() > 1) {
            parent[j] = pattern[j][1];
            ++numChildren[static_cast<std::size_t>(parent[j])];
        }
    }

    // Merge k into its child j when j is k's only child and their patterns agree.
    std::vector<std::int32_t> map(size, Oblio::NIL);
    std::size_t snodeSize = 0;
    for (std::size_t k = 0; k < size; ++k) {
        std::int32_t onlyChild = Oblio::NIL;
        if (numChildren[k] == 1) {
            for (std::size_t j = 0; j < k; ++j)
                if (parent[j] == static_cast<std::int32_t>(k)) {
                    onlyChild = static_cast<std::int32_t>(j);
                    break;
                }
        }
        const bool merge = (onlyChild != Oblio::NIL)
            && (pattern[static_cast<std::size_t>(onlyChild)].size() == pattern[k].size() + 1);

        if (merge)
            map[k] = map[static_cast<std::size_t>(onlyChild)];
        else
            map[k] = static_cast<std::int32_t>(snodeSize++);
    }
    snodeSizeOut = snodeSize;
    return map;
}

} // namespace OblioTest
