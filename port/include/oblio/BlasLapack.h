#pragma once

// BlasLapack.h - the BLAS and LAPACK routines the numeric factorization uses, wrapped.
//
// Two jobs. The first is mechanical: BLAS is a Fortran interface, so every argument is passed by
// address, and a call site without a wrapper reads
//
//     char uplo = 'L', trans = 'N';  int n = ..., k = ...;  double alpha = -1, beta = 1;
//     dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
//
// The inline overloads below hide that, taking values and dispatching on the scalar type. 0.9
// does the same and it is right.
//
// The second job is the one 0.9 gets wrong, and it is the reason this header is shaped the way
// it is.
//
// **The wrappers name operations, not routines.** 0.9 wraps BLAS routine by routine (SYRK, GEMM,
// TRSM, POTRF) and leaves the caller to pick 'T' versus 'C' and SYRK versus HERK. The convention
// therefore leaks into the engine, and 0.9's engine gets it wrong: its *complex* Cholesky calls
// SYRK, TRSM('T') and GEMM('N','T'), which are the complex-symmetric pattern, while its POTRF
// maps to zpotrf_, which is Hermitian. There is no HERK anywhere in 0.9's BlasLapack. For
// Hermitian A = LL^H the rank-k update must be L21 L21^H and the off-diagonal solve must be
// against L11^H; a plain transpose is correct only when L11 is real. Almost certainly never
// exercised on a genuinely complex Hermitian matrix.
//
// So `herk` here means **"A times A-conjugate-transpose"**, whatever that means for the scalar
// type: dsyrk_ for double (transpose and conjugate-transpose coincide over the reals), zherk_ for
// complex. And `Blas<Val>::conjTrans` is the character that means conjugate-transpose: 'T' for
// double, 'C' for complex. The Cholesky kernel is then one piece of code, correct for both, and
// 0.9's bug is unwriteable, because the engine never names SYRK or HERK and so cannot pick the
// wrong one.
//
// `syrk` and a literal 'T' remain available for LDL^T, where a *plain* transpose is what the
// mathematics asks for (complex LDL^T is the complex-symmetric case). There the algorithm asks
// for them deliberately rather than inheriting them by accident. See the factorization-space
// entry in DESIGN_DECISIONS for why symmetry is determined by (Val, factorization type) and is
// never a separate parameter.

#include <complex>

namespace Oblio {

// ---------------------------------------------------------------------------------------------
// The Fortran symbols. OBLIO_BLAS_UNDERSCORE selects the trailing-underscore convention, which
// is what Accelerate and the reference BLAS use.
// ---------------------------------------------------------------------------------------------

#ifdef OBLIO_BLAS_UNDERSCORE
#define OBLIO_BLAS(name) name##_
#else
#define OBLIO_BLAS(name) name
#endif

using Cplx = std::complex<double>;

extern "C" {

// Cholesky of a symmetric (Hermitian) positive definite matrix.
void OBLIO_BLAS(dpotrf)(const char* uplo, const int* n,
                        double* a, const int* lda, int* info);
void OBLIO_BLAS(zpotrf)(const char* uplo, const int* n,
                        Cplx* a, const int* lda, int* info);

// Triangular solve with multiple right-hand sides.
void OBLIO_BLAS(dtrsm)(const char* side, const char* uplo, const char* transa, const char* diag,
                       const int* m, const int* n,
                       const double* alpha, const double* a, const int* lda,
                       double* b, const int* ldb);
void OBLIO_BLAS(ztrsm)(const char* side, const char* uplo, const char* transa, const char* diag,
                       const int* m, const int* n,
                       const Cplx* alpha, const Cplx* a, const int* lda,
                       Cplx* b, const int* ldb);

// Symmetric rank-k update: C := alpha A A^T + beta C.
void OBLIO_BLAS(dsyrk)(const char* uplo, const char* trans,
                       const int* n, const int* k,
                       const double* alpha, const double* a, const int* lda,
                       const double* beta, double* c, const int* ldc);
void OBLIO_BLAS(zsyrk)(const char* uplo, const char* trans,
                       const int* n, const int* k,
                       const Cplx* alpha, const Cplx* a, const int* lda,
                       const Cplx* beta, Cplx* c, const int* ldc);

// Hermitian rank-k update: C := alpha A A^H + beta C. Note alpha and beta are *real*, which is
// forced: the diagonal of a Hermitian matrix is real, so a complex scale would break the
// symmetry. That is why the two herk() overloads below can share one signature.
void OBLIO_BLAS(zherk)(const char* uplo, const char* trans,
                       const int* n, const int* k,
                       const double* alpha, const Cplx* a, const int* lda,
                       const double* beta, Cplx* c, const int* ldc);

// General matrix multiply: C := alpha op(A) op(B) + beta C.
void OBLIO_BLAS(dgemm)(const char* transa, const char* transb,
                       const int* m, const int* n, const int* k,
                       const double* alpha, const double* a, const int* lda,
                       const double* b, const int* ldb,
                       const double* beta, double* c, const int* ldc);
void OBLIO_BLAS(zgemm)(const char* transa, const char* transb,
                       const int* m, const int* n, const int* k,
                       const Cplx* alpha, const Cplx* a, const int* lda,
                       const Cplx* b, const int* ldb,
                       const Cplx* beta, Cplx* c, const int* ldc);

}   // extern "C"

// ---------------------------------------------------------------------------------------------
// The scalar-type traits. One member, and it is the whole of what the scalar type decides.
// ---------------------------------------------------------------------------------------------

template<class Val>
struct Blas;

template<>
struct Blas<double> {
    // Conjugate-transpose. Over the reals conjugation is the identity, so this is the plain
    // transpose, and 'T' and 'C' would in fact behave identically in dtrsm/dgemm. We write 'T'
    // because it is what the operation *is* here, not because 'C' would fail.
    static constexpr char conjTrans = 'T';
};

template<>
struct Blas<Cplx> {
    static constexpr char conjTrans = 'C';
};

// ---------------------------------------------------------------------------------------------
// The operations.
// ---------------------------------------------------------------------------------------------

// Cholesky of the leading n x n block: A = L L^H (which for real is L L^T). On return the lower
// triangle holds L. info > 0 means the leading minor of that order is not positive definite,
// which is how a bad pivot is reported.
inline void potrf(char uplo, int n, double* a, int lda, int* info) {
    OBLIO_BLAS(dpotrf)(&uplo, &n, a, &lda, info);
}
inline void potrf(char uplo, int n, Cplx* a, int lda, int* info) {
    OBLIO_BLAS(zpotrf)(&uplo, &n, a, &lda, info);
}

// Triangular solve. B := alpha op(A)^-1 B, or B := alpha B op(A)^-1 for side == 'R'.
inline void trsm(char side, char uplo, char transa, char diag,
                 int m, int n, double alpha, const double* a, int lda,
                 double* b, int ldb) {
    OBLIO_BLAS(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
inline void trsm(char side, char uplo, char transa, char diag,
                 int m, int n, Cplx alpha, const Cplx* a, int lda,
                 Cplx* b, int ldb) {
    OBLIO_BLAS(ztrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

// Rank-k update with a CONJUGATE transpose: C := alpha A A^H + beta C, taking the lower triangle.
//
// This is the operation Cholesky needs, in both scalar types. For double, A^H is A^T and dsyrk_
// is the routine; for complex it is genuinely conjugate and zherk_ is. The engine asks for the
// operation and never learns which routine ran, which is precisely the point.
//
// alpha and beta are real in both overloads. For zherk_ that is forced by the mathematics (see
// above); for dsyrk_ it is simply what the routine takes. The signatures therefore agree, and one
// call site serves both.
inline void herk(char uplo, char trans, int n, int k,
                 double alpha, const double* a, int lda,
                 double beta, double* c, int ldc) {
    OBLIO_BLAS(dsyrk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}
inline void herk(char uplo, char trans, int n, int k,
                 double alpha, const Cplx* a, int lda,
                 double beta, Cplx* c, int ldc) {
    OBLIO_BLAS(zherk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

// Rank-k update with a PLAIN transpose: C := alpha A A^T + beta C, taking the lower triangle.
//
// Not used by Cholesky, and deliberately so. It is what LDL^T wants, where the complex case is
// complex-symmetric and a plain transpose is the correct operation. Kept here so that when LDL
// arrives it asks for this explicitly, rather than Cholesky reaching for it by accident, which is
// exactly the mistake 0.9 makes.
inline void syrk(char uplo, char trans, int n, int k,
                 double alpha, const double* a, int lda,
                 double beta, double* c, int ldc) {
    OBLIO_BLAS(dsyrk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}
inline void syrk(char uplo, char trans, int n, int k,
                 Cplx alpha, const Cplx* a, int lda,
                 Cplx beta, Cplx* c, int ldc) {
    OBLIO_BLAS(zsyrk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

// General matrix multiply. The transpose characters are the caller's, so a Cholesky call site
// passes Blas<Val>::conjTrans and an LDL one passes 'T'.
inline void gemm(char transa, char transb, int m, int n, int k,
                 double alpha, const double* a, int lda,
                 const double* b, int ldb,
                 double beta, double* c, int ldc) {
    OBLIO_BLAS(dgemm)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}
inline void gemm(char transa, char transb, int m, int n, int k,
                 Cplx alpha, const Cplx* a, int lda,
                 const Cplx* b, int ldb,
                 Cplx beta, Cplx* c, int ldc) {
    OBLIO_BLAS(zgemm)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

}   // namespace Oblio
