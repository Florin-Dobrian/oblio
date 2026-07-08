// Matrix_ext.cpp — Case 3: explicit instantiation + extern template
//
// Full template implementation lives here, not in the header. This translation
// unit is compiled once. The two explicit instantiations at the bottom cause
// the compiler to emit object code for double and complex<double>. All other
// translation units pick up those symbols at link time via Matrix_ext.o.

#include "Matrix_ext.h"
#include <cassert>

namespace Oblio {

template<class Val>
Matrix<Val>::Matrix()
    : mRows(0), mCols(0) {}

template<class Val>
Matrix<Val>::Matrix(std::size_t rows, std::size_t cols,
                    const std::vector<Val>& vals)
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

// ── Explicit instantiations ───────────────────────────────────────────────
// Compiled once. All other .cpp files link against these.
// To add float support: append one line below and one extern template in
// Matrix_ext.h. No other files need to change.

template class Matrix<double>;
template class Matrix<std::complex<double>>;

} // namespace Oblio
