#include "Matrix.h"

namespace Oblio {

template<class Val>
Matrix<Val>::Matrix()
    : mRows(0), mCols(0) {}

template<class Val>
Matrix<Val>::Matrix(std::size_t rows, std::size_t cols, const std::vector<Val>& vals)
    : mRows(rows), mCols(cols), mVals(vals) {
    assert(vals.size() == rows * cols);
}

template<class Val>
Val Matrix<Val>::operator()(std::size_t i, std::size_t j) const {
    assert(i < mRows && j < mCols);
    return mVals[i * mCols + j];
}

template<class Val>
Val& Matrix<Val>::operator()(std::size_t i, std::size_t j) {
    assert(i < mRows && j < mCols);
    return mVals[i * mCols + j];
}

template<class Val>
std::size_t Matrix<Val>::rows() const { return mRows; }

template<class Val>
std::size_t Matrix<Val>::cols() const { return mCols; }

template class Matrix<double>;
template class Matrix<std::complex<double>>;

} // namespace Oblio
