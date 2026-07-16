#pragma once

// SparseMatrix.h — sparse matrix in compressed sparse column (CSC) form.
//
// Minimal on purpose: enough structure for an ordering engine to compute a
// Permutation. Stored as three flat, contiguous vectors:
//   colPtr — length size()+1; column j occupies rowIdx/val[colPtr[j] .. colPtr[j+1]-1]
//   rowIdx — size nnz(); row indices, sorted ascending within each column
//   val    — size nnz(); the corresponding values
// A symmetric matrix is stored FULLY (both triangles), matching Oblio 0.9/10.12:
// each column holds its complete neighbour list plus the diagonal. Full storage
// lets the structural phases (ordering, elimination forest, symbolic) read each
// column's neighbours directly, with no lower->full expansion, and is the natural
// substrate for a future unsymmetric extension (factor the symmetrized structure).
// Values are carried but unused by the structural phases.
//
// Index types (see the "index types" design decision):
//   colPtr — OFFSETS into rowIdx/val -> std::size_t (never negative, may exceed 2^31).
//   rowIdx — row IDS -> std::int32_t (IDs may carry a -1/NIL sentinel elsewhere and
//            must match the graph/ordering convention). Loop counters are std::size_t;
//            an ID is cast to std::size_t only where it subscripts an array.
//
// Construction here is the single basic path: hand it arrays already in CSC form.
// Other builders (e.g. from COO triplets, with sorting / duplicate merging / zero-
// diagonal insertion) are deliberately not included yet — 0.9 is the oracle for those.

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class OrderEngine;

template<class Val>
class SparseMatrix {
public:
    SparseMatrix() = default;

    // Take arrays already in CSC form (moved in). colPtr has size+1 entries;
    // rowIdx and val each have colPtr[size] entries.
    SparseMatrix(std::size_t size,
                 std::vector<std::size_t>  colPtr,
                 std::vector<std::int32_t> rowIdx,
                 std::vector<Val>          val);

    std::size_t size() const;   // matrix dimension (number of rows / columns)
    std::size_t nnz()  const;   // number of stored entries

    // Structural read access. All callers (engines, tests, users) read A through
    // these; A is input and has no writer, so there is no friend. Bulk traversals
    // bind the returned container once and loop over it (vectorizes; no per-element
    // call), and a BLAS call takes .data() directly.
    const std::vector<std::size_t>&  colPtr() const;
    const std::vector<std::int32_t>& rowIdx() const;
    const std::vector<Val>&          val()    const;

private:
    std::size_t               mSize = 0;
    std::vector<std::size_t>  mColPtr;   // size mSize + 1 (offsets)
    std::vector<std::int32_t> mRowIdx;   // size nnz (row IDs)
    std::vector<Val>          mVal;      // size nnz
    std::size_t               mNnz = 0;  // number of stored entries (== mRowIdx.size()); nnz() returns this
};

extern template class SparseMatrix<double>;
extern template class SparseMatrix<std::complex<double>>;

} // namespace Oblio
