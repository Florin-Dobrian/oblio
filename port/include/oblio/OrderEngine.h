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

    // Order A into p. Reads only A's structure (colPtr/rowIdx); values unused.
    template<class Val>
    bool order(const SparseMatrix<Val>& A, Permutation& p) const;

private:
    OrderMethod mMethod = OrderMethod::MMD;

    bool orderNatural(std::size_t n, Permutation& p) const;
    bool orderMMD(std::size_t n,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& p) const;
    bool orderAMD(std::size_t n,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& p) const;
};

extern template bool OrderEngine::order(const SparseMatrix<double>&, Permutation&) const;
extern template bool OrderEngine::order(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
