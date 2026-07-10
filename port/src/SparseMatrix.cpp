#include "oblio/SparseMatrix.h"

#include "oblio/Types.h"   // MAX_IDX

#include <stdexcept>   // std::length_error
#include <utility>     // std::move

namespace Oblio {

template<class Val>
SparseMatrix<Val>::SparseMatrix(std::size_t numCols,
                                std::vector<std::size_t>  colPtr,
                                std::vector<std::int32_t> rowIdx,
                                std::vector<Val>          val)
    : mNumCols(numCols),
      mColPtr(std::move(colPtr)),
      mRowIdx(std::move(rowIdx)),
      mVal(std::move(val)) {
    // Row indices are std::int32_t, so the dimension and the number of stored
    // entries must both fit in that range (indices name rows; nnz is cast to int
    // at the AMD/MMD boundary). Guard at construction so an over-range input fails
    // loudly here rather than silently wrapping at the first narrowing cast.
    if (mNumCols > MAX_IDX || mRowIdx.size() > MAX_IDX)
        throw std::length_error(
            "SparseMatrix: dimension or nnz exceeds the std::int32_t index range");
}

template<class Val>
std::size_t SparseMatrix<Val>::numCols() const {
    return mNumCols;
}

template<class Val>
std::size_t SparseMatrix<Val>::nnz() const {
    return mColPtr.empty() ? 0 : mColPtr[mNumCols];
}

template<class Val>
const std::vector<std::size_t>& SparseMatrix<Val>::colPtr() const {
    return mColPtr;
}

template<class Val>
const std::vector<std::int32_t>& SparseMatrix<Val>::rowIdx() const {
    return mRowIdx;
}

template<class Val>
const std::vector<Val>& SparseMatrix<Val>::val() const {
    return mVal;
}

template class SparseMatrix<double>;
template class SparseMatrix<std::complex<double>>;

} // namespace Oblio
