#pragma once

// OrderEngine.h, computes a fill-reducing permutation of a SparseMatrix.
//
//   Natural, identity (no reordering)
//   MMD,     Multiple Minimum Degree (Liu/Sparspak, via 0.9)
//   MmdFlat,    Multiple Minimum Degree, Oblio's own, matching MMD's permutation (src/MmdFlat.cpp)
//   AMD,     Approximate Minimum Degree (SuiteSparse 3.3.4, BSD-3)
//   AmdFlat, Approximate Minimum Degree, Oblio's own, matching AMD's permutation (src/AmdFlat.cpp)
//
// FIVE ENUMERATORS, TWO OF THEM OURS, and each of ours returns its reference's exact permutation
// on every matrix the benchmarks cover. So a caller choosing MmdFlat over MMD, or AmdFlat over AMD,
// is choosing an implementation and not an ordering: the fill is identical by construction.
//
// THE EARLIER LADDER LAYERS WERE RETIRED ON 2026-08-21 and are in retired/. `MMD1`, `MMD2`, `AMD1`
// and `AMD2` were ours too, each carrying the base algorithm without its reference's later
// refinements, so they ordered differently and were not drop-in replacements. They were the
// ladder that got us to MmdFlat and AmdFlat, and the same ladder lives as working C++ and Python
// twins in experiments/ordering/. See retired/README.md.
//
// THE B AND C LAYERS ARE NOT IN THIS ENUM, DELIBERATELY. Three exist and are PERMANENT,
// `orderMmd3B`, `orderMmd3C` and `orderAmd3B`, each declared in its own header and reached as a
// free function: they are their originals on a different CLIQUE STORAGE scheme, the chained one
// genmmd uses and the compacted one `AMD_2` uses, both of which keep the ordering inside
// `O(n + m)` where ours does not. A B or a C is not an ordering but the SAME ordering computed
// differently, so it must return exactly its original's permutation and a difference is a defect
// in one of them. That makes it a measuring instrument, and an enumerator would put a benchmark's
// oracle into the library's public enum and into every switch over it, which is a cost paid by
// every reader of every switch for something no caller should choose. They are built, tested and
// benchmarked; they are simply not offered.
//
// Two lineages sit behind these names. MMD and AMD are vendored, self-contained codes operating on
// raw int CSC arrays (src/Mmd.cpp, src/Amd.cpp). MmdFlat and AmdFlat are ours, built over the
// shared quotient graph in include/oblio/QuotientGraph.h. Either way this engine is the seam that reads
// the matrix structure and fills the Permutation. Returns true on success.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {


enum class Ordering { Natural, MmdVendored, MmdFlat, AmdVendored, AmdFlat };

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
    // be an ordering that is always there. MmdFlat since 2026-08-07, and the reason is not that it
    // measured best. It returns genmmd's permutation EXACTLY, on every example, every square grid
    // and all 246 matrices in benchmarks/matrices, so its behavior is whatever thirty years of use
    // have established for a reference implementation. The alternative then was an earlier ladder
    // layer whose tie-break was ours and had been exercised on grids alone; a default is a bet on
    // the cases nobody has run, and reproducing the reference is the better bet. See
    // experiments/ordering's mmd3 section.
    Ordering mOrdering = Ordering::MmdFlat;

    bool orderNatural(std::size_t size, Permutation& P) const;
    // The two vendored orderings, declared only when private/ supplies their sources. Everything
    // else here is ours and is always present.
#ifdef OBLIO_VENDORED_ORDERINGS
    bool orderMmdVendored(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& P) const;
    bool orderAmdVendored(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& P) const;
#endif
    bool orderMmdFlat(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderAmdFlat(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
};

extern template bool OrderEngine::compute(const SparseMatrix<double>&, Permutation&) const;
extern template bool OrderEngine::compute(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
