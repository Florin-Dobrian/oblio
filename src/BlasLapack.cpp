#include "oblio/BlasLapack.h"

#include <cmath>

namespace Oblio {

namespace {

// The pieces the three kernels need, per scalar type. Kept local: nothing outside wants them.

double absVal(double v)      { return std::abs(v); }
double absVal(Cplx v)        { return std::abs(v); }

// Conjugate, or not. `hermitian` decides, and for double it is a no-op either way, which is the
// whole point: real symmetric and real Hermitian are the same case, so the branch costs nothing
// and says nothing.
double maybeConj(double v, bool)          { return v; }
Cplx   maybeConj(Cplx v, bool hermitian)  { return hermitian ? std::conj(v) : v; }

// A Hermitian factorization's D is real. Floating point will not make it exactly so, and the
// imaginary part it leaves is noise, not information. Discard it, which is what LAPACK's zpotrf
// does with the diagonal too.
double forceReal(double v, bool)          { return v; }
Cplx   forceReal(Cplx v, bool hermitian)  { return hermitian ? Cplx(v.real(), 0.0) : v; }

} // namespace

template<class Val>
void formUpper(int n, int k, const Val* l, int ldl, Val* u, int ldu,
               const Val* d, int ldd, bool hermitian) {
    // For each column i of L (there are k of them), and each row j (there are n):
    //
    //     U[i][j] = D[i] * conj?(L[j][i])
    //
    // U is L's transpose, scaled row-wise by D. The pointer walk is 0.9's: lp steps down L's
    // columns (by ldl), up steps across U's rows (by 1), dp walks D's diagonal (by ldd + 1).
    for (int i = 0, lp = 0, up = 0, dp = 0; i < k; ++i, lp += ldl, ++up, dp += ldd + 1)
        for (int j = 0, lq = lp, uq = up; j < n; ++j, ++lq, uq += ldu)
            u[uq] = d[dp] * maybeConj(l[lq], hermitian);
}

template<class Val>
void gemmLower(int n, int k, const Val* l, int ldl, const Val* u, int ldu, Val* a, int lda) {
    if (n == 1) {
        gemm('N', 'N', 1, 1, k, Val(-1), l, ldl, u, ldu, Val(1), a, lda);
        return;
    }

    const int n1 = n / 2;
    const int n2 = n - n1;

    const Val* l11 = l;
    const Val* l21 = l11 + n1;
    const Val* u11 = u;
    const Val* u12 = u11 + ldu * n1;
    Val*       a11 = a;
    Val*       a21 = a11 + n1;
    Val*       a22 = a21 + lda * n1;   // the (n1, n1) element: down n1 and across n1

    // The two diagonal blocks are symmetric, so recurse. The off-diagonal one is not, so a plain
    // GEMM fills it whole.
    gemmLower(n1, k, l11, ldl, u11, ldu, a11, lda);
    gemm('N', 'N', n2, n1, k, Val(-1), l21, ldl, u11, ldu, Val(1), a21, lda);
    gemmLower(n2, k, l21, ldl, u12, ldu, a22, lda);
}

template<class Val>
int ldl(int n, Val* a, int lda, double perturbation, int* numPerturbations, bool hermitian) {
    if (n == 1) {
        // The pivot. Force it real first if Hermitian, since that is what D is.
        *a = forceReal(*a, hermitian);

        // A static factorization cannot pivot, so a pivot too small to divide by has no remedy
        // but replacement. Count it: the caller is entitled to know we factored a different
        // matrix.
        if (absVal(*a) < std::abs(perturbation)) {
            *a = Val(perturbation);
            ++(*numPerturbations);
        }
        return 0;
    }

    const int n1 = n / 2;
    const int n2 = n - n1;

    Val* a11 = a;
    Val* a12 = a11 + lda * n1;
    Val* a21 = a11 + n1;
    Val* a22 = a21 + lda * n1;

    const int info = ldl(n1, a11, lda, perturbation, numPerturbations, hermitian);
    if (info > 0)
        return info;

    // L21 = A21 U11^-1, where U11 = D11 L11^T sits in a11's upper triangle. The solve is against
    // the upper, untransposed, which is exactly what having U stored buys us.
    trsm('R', 'U', 'N', 'N', n2, n1, Val(1), a11, lda, a21, lda);

    // U12 = D11 L21^T, into a12. Needed twice: it is the upper block of the factor, and it is the
    // right-hand factor of the Schur complement below.
    formUpper(n2, n1, a21, lda, a12, lda, a11, lda, hermitian);

    // A22 -= L21 U12. Symmetric, so only the lower triangle is filled.
    gemmLower(n2, n1, a21, lda, a12, lda, a22, lda);

    return ldl(n2, a22, lda, perturbation, numPerturbations, hermitian);
}

// Explicit instantiation. The headers declare; this file defines and instantiates, once. Adding a
// scalar type is one line per template.
template int  ldl<double>(int, double*, int, double, int*, bool);
template int  ldl<Cplx>(int, Cplx*, int, double, int*, bool);

template void formUpper<double>(int, int, const double*, int, double*, int,
                                const double*, int, bool);
template void formUpper<Cplx>(int, int, const Cplx*, int, Cplx*, int,
                              const Cplx*, int, bool);

template void gemmLower<double>(int, int, const double*, int, const double*, int, double*, int);
template void gemmLower<Cplx>(int, int, const Cplx*, int, const Cplx*, int, Cplx*, int);

} // namespace Oblio
