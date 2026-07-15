#include "MultiplyEngine.h"

namespace StorageOptions {

// The one multiply source. y += A*x, column by column, scattering into y. It reads each column
// through the storage's own lookups (rowIdxPtr / valPtr / colLen) and names no member, buffer, or
// layout, so this single body serves both storages; the explicit instantiations below specialize
// it per storage. Direct access: the pointer and length are read at the moment of use, held across
// nothing, so a growing dynamic storage has no cached pointer to dangle.
template <class Matrix>
void MultiplyEngine::multiply(const Matrix& A, const double* x, double* y) const {
    const std::size_t size = A.size();
    for (std::size_t j = 0; j < size; ++j) {
        const std::int32_t* rowIdx = A.rowIdxPtr(j);
        const double*       val    = A.valPtr(j);
        const std::size_t   n      = A.colLen(j);
        const double        xj     = x[j];

        for (std::size_t p = 0; p < n; ++p)
            y[rowIdx[p]] += val[p] * xj;
    }
}

// One line per storage, and adding a storage is one more line. This is the whole cost of the
// declaration-in-header / definition-here pattern.
template void MultiplyEngine::multiply<SparseMatrixStatic>(
    const SparseMatrixStatic&, const double*, double*) const;
template void MultiplyEngine::multiply<SparseMatrixDynamic>(
    const SparseMatrixDynamic&, const double*, double*) const;

// The baseline: the static matrix walked directly, no lookup call and no pointer array. It reaches
// the raw buffers through friendship, and that rawness is the point: it is the reference the
// templated multiply must match, so that "reaching a column through the lookup costs nothing" is a
// measured claim rather than a hoped-for one.
void MultiplyEngine::multiplyStatic(const SparseMatrixStatic& A, const double* x, double* y) const {
    const std::size_t   size   = A.mSize;
    const std::size_t*  colPtr = A.mColPtr.data();
    const std::int32_t* rowIdx = A.mRowIdx.data();
    const double*       val    = A.mVal.data();

    for (std::size_t j = 0; j < size; ++j) {
        const double xj = x[j];
        for (std::size_t p = colPtr[j]; p < colPtr[j + 1]; ++p)
            y[rowIdx[p]] += val[p] * xj;
    }
}

} // namespace StorageOptions
