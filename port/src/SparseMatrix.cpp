#include "oblio/SparseMatrix.h"

#include <utility>   // std::move

namespace Oblio {

template<class Val>
SparseMatrix<Val>::SparseMatrix(std::size_t numCols,
                                std::vector<std::size_t> colPtr,
                                std::vector<std::size_t> rowIdx,
                                std::vector<Val>         val)
    : mNumCols(numCols),
      mColPtr(std::move(colPtr)),
      mRowIdx(std::move(rowIdx)),
      mVal(std::move(val)) {}

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
const std::vector<std::size_t>& SparseMatrix<Val>::rowIdx() const {
    return mRowIdx;
}

template<class Val>
const std::vector<Val>& SparseMatrix<Val>::val() const {
    return mVal;
}

template class SparseMatrix<double>;
template class SparseMatrix<std::complex<double>>;

} // namespace Oblio
