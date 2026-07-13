#pragma once

// SymFactorEngine.h - computes the symbolic factorization of a SparseMatrix under a
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
#include "oblio/SymFactor.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class SymFactorEngine {
public:
    SymFactorEngine() = default;

    // Compute the symbolic factorization of A under p, from the forest f.
    //
    // The factor's pattern depends on A's pattern alone, never on its values, so the
    // implementation is the non-templated overload, taking colPtr and rowIdx, compiled
    // once. The templated overload is an adapter over it. The pattern overload is public:
    // a caller holding a graph with no numbers attached can compute the symbolic
    // factorization without inventing a scalar type to satisfy the signature.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p,
                 const ElmForest& f, SymFactor& s) const;

    bool compute(const std::vector<std::size_t>&  colPtr,
                 const std::vector<std::int32_t>& rowIdx,
                 const Permutation& p, const ElmForest& f, SymFactor& s) const;

private:
    // The front indices of each supernode, gathered contiguously. The map gives
    // each column's supernode, but a supernode's columns are generally scattered
    // through it, and the union recurrence needs them contiguous to be efficient.
    // Counting sort by supernode: walking columns in increasing order leaves each
    // supernode's front indices sorted in increasing order.
    //
    //   reads:   s.mSize, s.mSupSize, s.mIdxToSupIdx, s.mFrontSize
    // The two outputs are scratch for the union recurrence, not fields of s, so they stay
    // out-parameters.
    void gatherFrontalIndices(const SymFactor& s,
                        std::vector<std::size_t>&  frontSupPtr,
                        std::vector<std::int32_t>& frontRowIdx) const;

    // Sort each supernode's index set in increasing order, by a transposition
    // bucket sort: build the transpose (the supernodes each index appears in),
    // then read it back walking indices in increasing order. Linear in the number
    // of factor indices. The indices need not be sorted for correctness (an entry
    // is reached by indirection anyway), but numeric factorization assembles the
    // original matrix values assuming sorted front indices.
    //
    //   reads:   s.mSize, s.mSupSize, s.mSupPtr
    //   writes:  s.mRowIdx  (in place, permuted into sorted order)
    void sortIndices(SymFactor& s) const;
};

extern template bool SymFactorEngine::compute(const SparseMatrix<double>&, const Permutation&, const ElmForest&, SymFactor&) const;
extern template bool SymFactorEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, const ElmForest&, SymFactor&) const;

} // namespace Oblio
