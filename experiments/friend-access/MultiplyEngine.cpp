#include "MultiplyEngine.h"

namespace Oblio {

// API path: element access through the public operators.
template<class Val>
Vector<Val> MultiplyEngine::multiplyByApi(const Matrix<Val>& A,
                                          const Vector<Val>& x) const {
    Vector<Val> y(A.rows());
    for (std::size_t i = 0; i < A.rows(); ++i) {
        Val sum{0};
        for (std::size_t j = 0; j < A.cols(); ++j)
            sum += A(i, j) * x[j];
        y[i] = sum;
    }
    return y;
}

// Fast path: friend access to contiguous storage. Pointers fetched once; the inner
// loop is plain contiguous memory with no calls or bounds checks, so it vectorizes.
template<class Val>
Vector<Val> MultiplyEngine::multiplyDirectly(const Matrix<Val>& A,
                                             const Vector<Val>& x) const {
    const std::size_t rows = A.mRows;
    const std::size_t cols = A.mCols;
    Vector<Val> y(rows);

    const Val* a  = A.mVals.data();
    const Val* xp = x.mVals.data();
    Val*       yp = y.mVals.data();

    for (std::size_t i = 0; i < rows; ++i) {
        const Val* arow = a + i * cols;
        Val sum{0};
        for (std::size_t j = 0; j < cols; ++j)
            sum += arow[j] * xp[j];
        yp[i] = sum;
    }
    return y;
}

template Vector<double>
    MultiplyEngine::multiplyByApi(const Matrix<double>&, const Vector<double>&) const;
template Vector<std::complex<double>>
    MultiplyEngine::multiplyByApi(const Matrix<std::complex<double>>&,
                                  const Vector<std::complex<double>>&) const;
template Vector<double>
    MultiplyEngine::multiplyDirectly(const Matrix<double>&, const Vector<double>&) const;
template Vector<std::complex<double>>
    MultiplyEngine::multiplyDirectly(const Matrix<std::complex<double>>&,
                                     const Vector<std::complex<double>>&) const;

} // namespace Oblio
