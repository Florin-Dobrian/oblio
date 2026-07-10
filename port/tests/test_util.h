#pragma once
// Shared test utilities.

#include "oblio/SparseMatrix.h"

#include <vector>
#include <set>
#include <utility>
#include <cstddef>

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
    const std::size_t n = A.numCols();
    std::set<std::pair<std::size_t, std::size_t>> entries;
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t p = colPtr[j]; p < colPtr[j + 1]; ++p)
            entries.insert({static_cast<std::size_t>(rowIdx[p]), j});
    for (const auto& e : entries)
        if (entries.find({e.second, e.first}) == entries.end())
            return false;
    return true;
}

} // namespace OblioTest
