#include "Gemm.h"

// The Fortran BLAS prototype, spelled as the main tree spells it. See
// include/oblio/BlasLapack.h; OBLIO_BLAS_UNDERSCORE selects the trailing-underscore
// symbol convention, which Accelerate and the reference BLAS both use.
extern "C" {
#ifdef OBLIO_BLAS_UNDERSCORE
void dgemm_(const char* transa, const char* transb, const int* m, const int* n,
            const int* k, const double* alpha, const double* a, const int* lda,
            const double* b, const int* ldb, const double* beta, double* c,
            const int* ldc);
#define OBLIO_DGEMM dgemm_
#else
void dgemm(const char* transa, const char* transb, const int* m, const int* n,
           const int* k, const double* alpha, const double* a, const int* lda,
           const double* b, const int* ldb, const double* beta, double* c,
           const int* ldc);
#define OBLIO_DGEMM dgemm
#endif
}

namespace Oblio {

// Column major, so the innermost loop walks a column of C and a column of A with
// unit stride and holds one entry of B fixed. The natural ordering for this layout,
// and nothing beyond it: no blocking, no unrolling, no intrinsics.
void gemmHand(std::size_t n, const double* a, const double* b, double* c)
{
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t k = 0; k < n; ++k) {
            const double bkj = b[k + j * n];
            const double* ak = a + k * n;
            double* cj = c + j * n;
            for (std::size_t i = 0; i < n; ++i) cj[i] += ak[i] * bkj;
        }
    }
}

void gemmBlas(std::size_t n, const double* a, const double* b, double* c)
{
    // The one narrowing cast, size_t to the BLAS int, in the direction that loses
    // information. Safe here because n is a test parameter of a few thousand at most.
    const int m = static_cast<int>(n);
    const double alpha = 1.0;
    const double beta = 1.0;  // accumulate, matching gemmHand's C := C + A B
    OBLIO_DGEMM("N", "N", &m, &m, &m, &alpha, a, &m, b, &m, &beta, c, &m);
}

}  // namespace Oblio
