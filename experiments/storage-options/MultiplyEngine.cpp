#include "MultiplyEngine.h"

namespace StorageOptions {

// The one multiply source. y = A*x, a pure multiply: it overwrites y and never reads the old value
// (BLAS's beta = 0 case). The loop is column-outer, so each y[i] is touched by several columns and
// cannot be written in one shot; y is therefore zeroed once up front and then accumulated into,
// which a column-outer sweep must do regardless. It reads each column through the storage's own
// accessors (rowIdx / val / colSize) and names no member, buffer, or layout, so this single body
// serves both storages; the explicit instantiations below specialize it per storage. Direct access:
// the pointer and size are read at the moment of use, held across nothing, so a growing dynamic
// storage has no cached pointer to dangle.
template <class Matrix>
void MultiplyEngine::multiply(const Matrix& A, const double* x, double* y) const {
    const std::size_t size = A.size();
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(size); ++i)
        y[i] = 0.0;
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(size); ++j) {
        const std::int32_t* rowIdx  = A.rowIdx(j);
        const double*       val     = A.val(j);
        const std::size_t   colSize = A.colSize(j);
        const double        xj      = x[j];

        // cp is the column position: it scans this column's rowIdx/val, whatever backs them (a slice
        // of a flat buffer here for static, an inner vector for dynamic). Not a colPtr index.
        for (std::size_t cp = 0; cp < colSize; ++cp)
            y[rowIdx[cp]] += val[cp] * xj;
    }
}

// One line per storage, and adding a storage is one more line. This is the whole cost of the
// declaration-in-header / definition-here pattern.
template void MultiplyEngine::multiply<SparseMatrixStatic>(
    const SparseMatrixStatic&, const double*, double*) const;
template void MultiplyEngine::multiply<SparseMatrixDynamic>(
    const SparseMatrixDynamic&, const double*, double*) const;

// The baseline: the static matrix walked directly, no per-column accessor call and no pointer array.
// It reads the raw CSC buffers through the class's public colPtr() / rowIdx() / val(), the same
// arrays the main-code matrix exposes, so it needs no friendship. That rawness is the point: it is
// the reference the templated multiply must match, so that "reaching a column through the accessor
// costs nothing" is a measured claim rather than a hoped-for one. Same pure-multiply contract: y is
// overwritten, zeroed once then accumulated, since it too is column-outer.
void MultiplyEngine::multiplyStatic(const SparseMatrixStatic& A, const double* x, double* y) const {
    const std::size_t                size   = A.size();
    const std::vector<std::size_t>&  colPtr = A.colPtr();
    const std::vector<std::int32_t>& rowIdx = A.rowIdx();
    const std::vector<double>&       val    = A.val();

    for (std::int32_t i = 0; i < static_cast<std::int32_t>(size); ++i)
        y[i] = 0.0;
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(size); ++j) {
        const double xj = x[j];
        for (std::size_t cp = colPtr[j]; cp < colPtr[j + 1]; ++cp)
            y[rowIdx[cp]] += val[cp] * xj;
    }
}

} // namespace StorageOptions
