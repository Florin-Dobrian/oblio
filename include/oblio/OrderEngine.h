#pragma once

// OrderEngine.h — computes a fill-reducing permutation of a SparseMatrix.
//
//   Natural — identity (no reordering)
//   MMD     — Multiple Minimum Degree (Liu/Sparspak, via 0.9)
//   AMD     — Approximate Minimum Degree (SuiteSparse 3.3.4, BSD-3)
//
// The MMD/AMD algorithms are vendored, self-contained codes operating on raw int
// CSC arrays (src/Mmd.cpp, src/Amd.cpp); this engine is the seam that reads the
// matrix structure and fills the Permutation. Returns true on success.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {

enum class OrderMethod { Natural, MMD, AMD };

class OrderEngine {
public:
    OrderEngine() = default;
    explicit OrderEngine(OrderMethod method) : mMethod(method) {}

    void        setMethod(OrderMethod method) { mMethod = method; }
    OrderMethod method() const                { return mMethod; }

    // Order A into p.
    //
    // An ordering is a pure graph operation: AMD and MMD read the sparsity pattern and
    // would not know what to do with a value. So the implementation is the non-templated
    // overload, taking colPtr and rowIdx, and it is compiled once. The templated overload
    // is an adapter over it, for the common case of holding a matrix. The pattern overload
    // is public: a caller holding a graph with no numbers attached can order it without
    // inventing a scalar type to satisfy the signature.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, Permutation& p) const;

    bool compute(const std::vector<std::size_t>&  colPtr,
                 const std::vector<std::int32_t>& rowIdx,
                 Permutation& p) const;

private:
    OrderMethod mMethod = OrderMethod::MMD;

    bool orderNatural(std::size_t size, Permutation& p) const;
    bool orderMMD(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& p) const;
    bool orderAMD(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& p) const;
};

extern template bool OrderEngine::compute(const SparseMatrix<double>&, Permutation&) const;
extern template bool OrderEngine::compute(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
