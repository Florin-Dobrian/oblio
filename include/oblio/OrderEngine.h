#pragma once

// OrderEngine.h, computes a fill-reducing permutation of a SparseMatrix.
//
//   Natural, identity (no reordering)
//   MMD,     Multiple Minimum Degree (Liu/Sparspak, via 0.9)
//   MMD1,    Multiple Minimum Degree, Oblio's own, the batch alone (src/Mmd1.cpp)
//   MMD2,    MMD1 plus the prepass, q2h refresh, pair merging and outmatched marking
//   AMD,     Approximate Minimum Degree (SuiteSparse 3.3.4, BSD-3)
//   AMD1,    Approximate Minimum Degree, Oblio's own, the bound alone (src/Amd1.cpp)
//   AMD2,    AMD1 plus aggressive absorption and hash detection (src/Amd2.cpp)
//
// Two lineages sit behind those names. MMD and AMD are vendored, self-contained codes
// operating on raw int CSC arrays (src/Mmd.cpp, src/Amd.cpp). MMD1 and AMD1 are ours,
// built over the shared quotient graph in src/QuotientGraph.cpp, and each carries the
// base algorithm without its vendored counterpart's later refinements, so they order
// differently and are not drop-in replacements. Either way this engine is the seam that
// reads the matrix structure and fills the Permutation. Returns true on success.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {

enum class OrderMethod { Natural, MMD, MMD1, MMD2, AMD, AMD1, AMD2 };

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
    bool orderMMD1(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& p) const;
    bool orderMMD2(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& p) const;
    bool orderAMD(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& p) const;
    bool orderAMD1(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& p) const;
    bool orderAMD2(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& p) const;
};

extern template bool OrderEngine::compute(const SparseMatrix<double>&, Permutation&) const;
extern template bool OrderEngine::compute(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
