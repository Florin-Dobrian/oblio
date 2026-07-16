#include "oblio/SparseMatrix.h"

#include "oblio/Types.h"   // checkIndexRange

#include <utility>     // std::move

namespace Oblio {

template<class Val>
SparseMatrix<Val>::SparseMatrix(std::size_t size,
                                std::vector<std::size_t>  colPtr,
                                std::vector<std::int32_t> rowIdx,
                                std::vector<Val>          val)
    : mSize(size),
      mColPtr(std::move(colPtr)),
      mRowIdx(std::move(rowIdx)),
      mVal(std::move(val)),
      mNnz(mRowIdx.size()) {
    // Row indices are std::int32_t, so the dimension and the number of stored entries must both fit
    // in that range (indices name rows; nnz is cast to int at the AMD/MMD boundary). The vectors are
    // moved in, not allocated here, so there is nothing to precede: guard after the moves. Two calls,
    // so the message names whichever overflowed.
    checkIndexRange(mSize, "SparseMatrix size");
    checkIndexRange(mNnz, "SparseMatrix nnz");
}

template<class Val>
std::size_t SparseMatrix<Val>::size() const {
    return mSize;
}

template<class Val>
std::size_t SparseMatrix<Val>::nnz() const {
    return mNnz;
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
