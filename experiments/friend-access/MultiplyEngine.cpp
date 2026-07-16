#include "MultiplyEngine.h"

#include <cstdint>   // std::int32_t

// Fortran BLAS gemv. Declared here directly (underscore ABI, as used by Accelerate
// on macOS and reference/OpenBLAS on Linux). The real tree routes these through a
// BlasLapack layer with the OBLIO_BLAS_UNDERSCORE portability macro; kept inline here
// so the example is self-contained. std::complex<double> is layout-compatible with
// Fortran complex*16, so its pointers pass through directly.
extern "C" {
void dgemv_(const char* trans, const int* m, const int* n,
            const double* alpha, const double* a, const int* lda,
            const double* x, const int* incx,
            const double* beta, double* y, const int* incy);
void zgemv_(const char* trans, const int* m, const int* n,
            const std::complex<double>* alpha, const std::complex<double>* a, const int* lda,
            const std::complex<double>* x, const int* incx,
            const std::complex<double>* beta, std::complex<double>* y, const int* incy);
}

namespace Oblio {
namespace {

// y = (A^T) x, where A is blasM x blasN column-major with lda = blasM. Overloaded on
// scalar type so the caller is type-generic.
void blasGemvTrans(int blasM, int blasN, const double* a, const double* x, double* y) {
    const char trans = 'T';
    const int  one   = 1;
    const double alpha = 1.0, beta = 0.0;
    dgemv_(&trans, &blasM, &blasN, &alpha, a, &blasM, x, &one, &beta, y, &one);
}
void blasGemvTrans(int blasM, int blasN, const std::complex<double>* a,
                   const std::complex<double>* x, std::complex<double>* y) {
    const char trans = 'T';   // 'T' = transpose, NOT 'C' — no conjugation
    const int  one   = 1;
    const std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
    zgemv_(&trans, &blasM, &blasN, &alpha, a, &blasM, x, &one, &beta, y, &one);
}

} // namespace

// API path: element access through the public operators.
template<class Val>
Vector<Val> MultiplyEngine::multiplyByApi(const Matrix<Val>& A,
                                          const Vector<Val>& x) const {
    Vector<Val> y(A.rows());
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(A.rows()); ++i) {
        Val sum{0};
        for (std::int32_t j = 0; j < static_cast<std::int32_t>(A.cols()); ++j)
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

    for (std::int32_t i = 0; i < static_cast<std::int32_t>(rows); ++i) {
        const Val* arow = a + i * cols;
        Val sum{0};
        for (std::int32_t j = 0; j < static_cast<std::int32_t>(cols); ++j)
            sum += arow[j] * xp[j];
        yp[i] = sum;
    }
    return y;
}

// BLAS path: same friend access to the raw block, then gemv. Our matrix is row-major
// (rows x cols); the same buffer is column-major A^T (cols x rows). gemv with TRANS='T'
// on that computes (A^T)^T x = A x — so BLAS M = cols, N = rows, lda = cols.
// (The real solver stores dense blocks column-major to feed BLAS with no transpose;
// the transpose here is only because this example reuses the row-major Matrix.)
template<class Val>
Vector<Val> MultiplyEngine::multiplyWithBlas(const Matrix<Val>& A,
                                             const Vector<Val>& x) const {
    const int rows = static_cast<int>(A.mRows);
    const int cols = static_cast<int>(A.mCols);
    Vector<Val> y(static_cast<std::size_t>(rows));
    if (rows == 0 || cols == 0) return y;

    const Val* a  = A.mVals.data();
    const Val* xp = x.mVals.data();
    Val*       yp = y.mVals.data();

    blasGemvTrans(cols, rows, a, xp, yp);   // BLAS M=cols, N=rows, lda=cols
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
template Vector<double>
    MultiplyEngine::multiplyWithBlas(const Matrix<double>&, const Vector<double>&) const;
template Vector<std::complex<double>>
    MultiplyEngine::multiplyWithBlas(const Matrix<std::complex<double>>&,
                                     const Vector<std::complex<double>>&) const;

} // namespace Oblio
