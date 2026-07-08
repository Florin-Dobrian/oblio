// Matrix_exp.cc — Case 2: explicit instantiation (forcing only)
//
// Full template implementation lives here, not in the header. This translation
// unit is compiled once. The two explicit instantiations at the bottom force the
// compiler to emit object code for double and complex<double>. Other translation
// units see only declarations in Matrix_exp.h, so they cannot implicitly
// instantiate — they pick up these symbols at link time via Matrix_exp.o.

#include "Matrix_exp.h"
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
// Compiled once. All other .cc files link against these.
// To add float support: append one line below. Nothing else changes — unlike
// Case 3, there is no extern template line in the header to keep in sync.

template class Matrix<double>;
template class Matrix<std::complex<double>>;

} // namespace Oblio
