#pragma once

// ElmForestEngine.h — computes the elimination forest of a SparseMatrix under a
// given Permutation. Minimal: produces the parent links (etree). Returns true on
// success. Reads the matrix/permutation via their public accessors and fills the
// forest via friend access.
//
// A is stored full-symmetric, so the etree reads each column's neighbours
// directly (the diagonal self-skips) — no expansion step.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/ElmForest.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class ElmForestEngine {
public:
    ElmForestEngine() = default;

    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p, ElmForest& f) const;

private:
    // Parent links (etree) via Liu's path-compression construction, on the
    // full-symmetric structure, in the permuted (factor) order.
    void computeParent(std::size_t n,
                       const std::vector<std::size_t>&  colPtr,
                       const std::vector<std::int32_t>& rowIdx,
                       const std::vector<std::int32_t>& oldToNew,
                       const std::vector<std::int32_t>& newToOld,
                       std::vector<std::int32_t>& parent) const;
};

extern template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
extern template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
