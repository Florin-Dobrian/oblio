#pragma once

// SymFactEngine.h - computes the symbolic factorization of a SparseMatrix under a
// given Permutation, from an already-computed ElmForest. Copies the forest
// attributes the factor needs, then computes each supernode's index set as a
// union: the sparsity patterns of the supernode's own front columns, plus the
// update indices of its children. Returns true on success.
//
// The recurrence runs over supernodes in increasing order, which is a topological
// order (a parent's index is larger than its children's), so a supernode's
// children are complete by the time it is reached.
//
// Reads the matrix, permutation and forest via their public accessors and fills
// the symbolic factorization via friend access.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/ElmForest.h"
#include "oblio/SymFact.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class SymFactEngine {
public:
    SymFactEngine() = default;

    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p,
                 const ElmForest& f, SymFact& s) const;

private:
    // The front indices of each supernode, gathered contiguously. The map gives
    // each column's supernode, but a supernode's columns are generally scattered
    // through it, and the union recurrence needs them contiguous to be efficient.
    // Counting sort by supernode: walking columns in increasing order leaves each
    // supernode's front indices sorted in increasing order.
    void gatherFrontIdx(std::size_t size, std::size_t supSize,
                        const std::vector<std::int32_t>& idxToSupIdx,
                        const std::vector<std::size_t>&  frontSize,
                        std::vector<std::size_t>&  frontIdxPtr,
                        std::vector<std::int32_t>& frontIdx) const;

    // Sort each supernode's index set in increasing order, by a transposition
    // bucket sort: build the transpose (the supernodes each index appears in),
    // then read it back walking indices in increasing order. Linear in the number
    // of factor indices. The indices need not be sorted for correctness (an entry
    // is reached by indirection anyway), but numeric factorization assembles the
    // original matrix values assuming sorted front indices.
    void sortIdx(std::size_t size, std::size_t supSize,
                 const std::vector<std::size_t>& idxPtr,
                 std::vector<std::int32_t>& idx) const;
};

extern template bool SymFactEngine::compute(const SparseMatrix<double>&, const Permutation&, const ElmForest&, SymFact&) const;
extern template bool SymFactEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, const ElmForest&, SymFact&) const;

} // namespace Oblio
