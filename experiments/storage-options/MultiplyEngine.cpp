#include "MultiplyEngine.h"

namespace StorageOptions {

// The one algorithm. Note what is absent: any mention of SparseMatrixStatic, of SparseMatrixDynamic,
// of colPtr, of an inner vector. It walks pointers.
void MultiplyEngine::multiply(std::size_t size,
                              const std::int32_t* const* rowIdxPtr,
                              const double* const*       valPtr,
                              const std::size_t*         len,
                              const double* x, double* y) const {
    for (std::size_t j = 0; j < size; ++j) {
        const std::int32_t* rowIdx = rowIdxPtr[j];
        const double*       val    = valPtr[j];
        const std::size_t   n      = len[j];
        const double        xj     = x[j];

        for (std::size_t p = 0; p < n; ++p)
            y[rowIdx[p]] += val[p] * xj;
    }
}

// CSC: colPtr holds indices, so take the address of the element they index. The pointers all
// land in one buffer, at increasing addresses, which is exactly why CSC streams well.
void MultiplyEngine::columnPointers(const SparseMatrixStatic& A,
                                    std::vector<const std::int32_t*>& rowIdxPtr,
                                    std::vector<const double*>&       valPtr,
                                    std::vector<std::size_t>&         len) const {
    const std::size_t size = A.mSize;
    rowIdxPtr.resize(size);
    valPtr.resize(size);
    len.resize(size);

    for (std::size_t j = 0; j < size; ++j) {
        const std::size_t from = A.mColPtr[j];
        rowIdxPtr[j] = A.mRowIdx.data() + from;
        valPtr[j]    = A.mVal.data() + from;
        len[j]       = A.mColPtr[j + 1] - from;
    }
}

// VV: the pointers are already there, one per inner vector. They land wherever the allocator
// put each column, which is the whole difference, and the only difference.
void MultiplyEngine::columnPointers(const SparseMatrixDynamic& A,
                                    std::vector<const std::int32_t*>& rowIdxPtr,
                                    std::vector<const double*>&       valPtr,
                                    std::vector<std::size_t>&         len) const {
    const std::size_t size = A.mSize;
    rowIdxPtr.resize(size);
    valPtr.resize(size);
    len.resize(size);

    for (std::size_t j = 0; j < size; ++j) {
        rowIdxPtr[j] = A.mRowIdx[j].data();
        valPtr[j]    = A.mVal[j].data();
        len[j]       = A.mRowIdx[j].size();
    }
}

// The baseline: the static matrix with no pointer arrays at all.
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
