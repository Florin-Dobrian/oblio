#pragma once

// MultiplyEngine.h - one sparse matvec, over either storage.
//
// The question: can a single multiply function serve both a CSC matrix and a
// vector-of-vectors matrix?
//
// The answer here is yes, and not by templating the algorithm. There is exactly **one**
// compiled multiply. It takes no matrix at all. It takes three arrays:
//
//     rowIdxPtr[j]   where column j's row indices start
//     valPtr[j]      where column j's values start
//     len[j]         how many entries column j has
//
// Both storages can fill those. CSC points into its single contiguous buffer
// (&mRowIdx[mColPtr[j]]); VV reads each inner vector's data(). Once filled, the arrays are
// indistinguishable, and multiply() cannot tell which class produced them, because there is
// nothing left to tell apart.
//
// That is the whole demonstration. The storage question and the algorithm question are
// separable: the layout decides where the pointers come from, and nothing else.
//
// Two extractors (columnPointers, one overload per storage) and one algorithm. Plus a
// hand-written CSC multiply as the honest baseline: it skips the pointer arrays entirely and
// walks colPtr directly, which is what any sane CSC code does. If the pointer-array version
// matches it, the indirection costs nothing.
//
// Val is fixed to double. This experiment is about storage, not scalar type.

#include "SparseMatrixCsc.h"
#include "SparseMatrixVv.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace StorageOptions {

class MultiplyEngine {
public:
    // The one algorithm. y = A*x, column by column, scattering into y.
    //
    // It sees pointers and lengths. No matrix, no template parameter, no virtual call, no
    // knowledge of any storage class. One function, compiled once.
    void multiply(std::size_t size,
                  const std::int32_t* const* rowIdxPtr,
                  const double* const*       valPtr,
                  const std::size_t*         len,
                  const double* x, double* y) const;

    // Extract the column pointers. One overload per storage; both produce the same three
    // arrays, and after this call the storage is out of the picture.
    void columnPointers(const SparseMatrixCsc& A,
                        std::vector<const std::int32_t*>& rowIdxPtr,
                        std::vector<const double*>&       valPtr,
                        std::vector<std::size_t>&         len) const;

    void columnPointers(const SparseMatrixVv& A,
                        std::vector<const std::int32_t*>& rowIdxPtr,
                        std::vector<const double*>&       valPtr,
                        std::vector<std::size_t>&         len) const;

    // The baseline. Hand-written CSC: no pointer arrays, walks colPtr directly. This is what
    // one would write if CSC were the only storage, and it is what the general version must
    // match to earn its keep.
    void multiplyCsc(const SparseMatrixCsc& A, const double* x, double* y) const;
};

} // namespace StorageOptions
