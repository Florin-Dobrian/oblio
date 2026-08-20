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
//   AMD3,    AMD2 with the vendored routine's list order, to match its permutation (src/Amd3.cpp)
// A trailing digit means a different ordering: AMD2 has mechanisms AMD1 lacks, so their
// permutations and their fill differ and both are correct. Every one of the nine above is a
// supported ordering a caller may ask for.
//
// THE B LAYERS ARE NOT IN THIS ENUM, DELIBERATELY. Two remain and are PERMANENT, `orderMmd3B` and
// `orderAmd3B`, each declared in its own header and reached as a free function: they are their
// originals on the two VENDORED clique storage schemes, both of which keep the ordering inside
// `O(n + m)` where ours does not. `orderAmd1B` and `orderAmd2B` were two more until 2026-08-16,
// when the schedule they carried won and moved into their originals. A
// B is not an ordering but the SAME ordering computed differently, so it must return exactly its
// original's permutation and a difference is a defect in one of them. That makes it a measuring
// instrument, and an enumerator would put a benchmark's oracle into the library's public enum
// and into every switch over it, which is a cost paid by every reader of every switch for
// something no caller should choose. They are built, tested and benchmarked; they are simply not
// offered. See docs/NEXT.md, item 2c.
//
// Two lineages sit behind those names. MMD and AMD are vendored, self-contained codes
// operating on raw int CSC arrays (src/Mmd.cpp, src/Amd.cpp). MMD1 and AMD1 are ours,
// built over the shared quotient graph in include/oblio/QuotientGraph.h, and each carries the
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


enum class Ordering { Natural, MMD, MMD1, MMD2, MMD3, AMD, AMD1, AMD2, AMD3 };

class OrderEngine {
public:
    OrderEngine() = default;
    explicit OrderEngine(Ordering ordering) : mOrdering(ordering) {}

    void     setOrdering(Ordering ordering) { mOrdering = ordering; }
    Ordering ordering() const               { return mOrdering; }

    // Order A into P.
    //
    // An ordering is a pure graph operation: AMD and MMD read the sparsity pattern and
    // would not know what to do with a value. So the implementation is the non-templated
    // overload, taking colPtr and rowIdx, and it is compiled once. The templated overload
    // is an adapter over it, for the common case of holding a matrix. The pattern overload
    // is public: a caller holding a graph with no numbers attached can order it without
    // inventing a scalar type to satisfy the signature.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, Permutation& P) const;

    bool compute(const std::vector<std::size_t>&  colPtr,
                 const std::vector<std::int32_t>& rowIdx,
                 Permutation& P) const;

private:
    // Our MMD, not the vendored one, because the vendored orderings are optional: a default has to
    // be an ordering that is always there. MMD3 since 2026-08-07, and the reason is not that it
    // measured best. It returns genmmd's permutation EXACTLY, on every example and every square
    // grid tested, so its behavior is whatever thirty years of use have established for a
    // reference implementation, where MMD2's tie-break is ours and has been exercised on grids
    // alone. On those grids MMD2 is in fact very slightly better at 32 a side and a few percent
    // worse above it; a default is a bet on the cases nobody has run, and reproducing the
    // reference is the better bet. See experiments/ordering's mmd3 section.
    Ordering mOrdering = Ordering::MMD3;

    bool orderNatural(std::size_t size, Permutation& P) const;
    // The two vendored orderings, declared only when private/ supplies their sources. Everything
    // else here is ours and is always present.
#ifdef OBLIO_VENDORED_ORDERINGS
    bool orderMMD(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& P) const;
    bool orderAMD(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& P) const;
#endif
    bool orderMMD1(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderMMD2(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderMMD3(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderAMD3(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderAMD1(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderAMD2(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
};

extern template bool OrderEngine::compute(const SparseMatrix<double>&, Permutation&) const;
extern template bool OrderEngine::compute(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
