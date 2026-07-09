#pragma once

// SparseMatrix.h — sparse matrix in compressed sparse column (CSC) form.
//
// Minimal on purpose: enough structure for an ordering engine to compute a
// Permutation. Stored as three flat, contiguous vectors:
//   colPtr — size numCols()+1; column j occupies rowIdx/val[colPtr[j] .. colPtr[j+1]-1]
//   rowIdx — size nnz(); row indices, sorted ascending within each column
//   val    — size nnz(); the corresponding values
// For the symmetric case this is the lower triangle, diagonal included. Values are
// carried but unused by the structural (ordering/symbolic) phases.
//
// Construction here is the single basic path: hand it arrays already in CSC form.
// Other builders (e.g. from COO triplets, with sorting / duplicate merging / zero-
// diagonal insertion) are deliberately not included yet — 0.9 is the oracle for those.

#include <vector>
#include <complex>
#include <cstddef>

namespace Oblio {

class OrderEngine;

template<class Val>
class SparseMatrix {
public:
    SparseMatrix() = default;

    // Take arrays already in CSC form (moved in). colPtr has numCols+1 entries;
    // rowIdx and val each have colPtr[numCols] entries.
    SparseMatrix(std::size_t numCols,
                 std::vector<std::size_t> colPtr,
                 std::vector<std::size_t> rowIdx,
                 std::vector<Val>         val);

    std::size_t numCols() const;   // matrix dimension (number of columns)
    std::size_t nnz()     const;   // number of stored entries

    // Structural read access for non-hot-path callers. Engines that traverse the
    // structure on hot paths use friend access to the members directly.
    const std::vector<std::size_t>& colPtr() const;
    const std::vector<std::size_t>& rowIdx() const;
    const std::vector<Val>&         val()    const;

private:
    std::size_t              mNumCols = 0;
    std::vector<std::size_t> mColPtr;   // size mNumCols + 1
    std::vector<std::size_t> mRowIdx;   // size nnz
    std::vector<Val>         mVal;      // size nnz

    friend class OrderEngine;   // reads structure directly; add engines as needed
};

extern template class SparseMatrix<double>;
extern template class SparseMatrix<std::complex<double>>;

} // namespace Oblio
