#include "oblio/BlasLapack.h"

#include "oblio/Types.h"

#include <cmath>

namespace Oblio {

namespace {

// The pieces the three kernels need, per scalar type. Kept local: nothing outside wants them.

double absVal(double v)      { return std::abs(v); }
double absVal(Cplx v)        { return std::abs(v); }

// maybeConjugate and forceReal live in Types.h, beside the `hermitian` predicate that drives them.
// A Hermitian factorization's D is real; floating point will not make it exactly so, and the
// imaginary part left behind is noise rather than information, which is why forceReal discards it,
// as LAPACK's zpotrf does with its diagonal.

} // namespace

template<class Val>
void formStaticUpper(int n, int k, const Val* l, int ldl, Val* u, int ldu,
                     const Val* d, int ldd, bool hermitian) {
    // For each column i of L (there are k of them), and each row j (there are n):
    //
    //     U[i][j] = D[i] * conj?(L[j][i])
    //
    // U is L's transpose, scaled row-wise by D. D is a plain diagonal, one entry per column. The
    // pointer walk is 0.9's: lp steps down L's columns (by ldl), up steps across U's rows (by 1),
    // dp walks D's diagonal (by ldd + 1).
    for (int i = 0, lp = 0, up = 0, dp = 0; i < k; ++i, lp += ldl, ++up, dp += ldd + 1)
        for (int j = 0, lq = lp, uq = up; j < n; ++j, ++lq, uq += ldu)
            u[uq] = d[dp] * maybeConjugate(l[lq], hermitian);
}

template<class Val>
void formDynamicUpper(int n, int k, const Val* l, int ldl, Val* u, int ldu,
                      const Val* d, int ldd, const std::int32_t* pivotType,
                      const std::int32_t* nodeIdx, bool hermitian) {
    // The same U := D conj?(L^T) as formStaticUpper, but D is block-diagonal: dynamic LDL pivots in
    // 1x1 and 2x2 blocks, marked by pivotType (indexed by the global node nodeIdx[c], as LAPACK's ipiv
    // marks its own). A 1x1 column scales like the static case, one d per column. A 2x2 couples its
    // two columns through the four entries of the pivot, which no single-d-per-column walk expresses.
    //
    // Positions are column-major: `at` into L and D (leading dimension ldl == ldd, both into the
    // supernode's val), `atU` into U (leading dimension ldu).
    const auto at  = [ldl](int r, int c) { return static_cast<std::ptrdiff_t>(c) * ldl + r; };
    const auto atU = [ldu](int r, int c) { return static_cast<std::ptrdiff_t>(c) * ldu + r; };
    (void)ldd;   // ldd == ldl here: D and L share the supernode's leading dimension.

    for (int c = 0; c < k; ) {
        if (pivotType[nodeIdx[c]] == 1) {
            const Val dc = d[at(c, c)];
            for (int j = 0; j < n; ++j)
                u[atU(c, j)] = dc * maybeConjugate(l[at(j, c)], hermitian);
            ++c;
        } else {                                   // a 2x2 pivot: two columns solved together
            const int c1 = c;
            const int c2 = c + 1;

            const Val d11 = d[at(c1, c1)];
            const Val d12 = d[at(c1, c2)];   // the upper part, written back by the pivot
            const Val d21 = d[at(c2, c1)];   // the lower part, never overwritten
            const Val d22 = d[at(c2, c2)];

            for (int j = 0; j < n; ++j) {
                const Val l1 = maybeConjugate(l[at(j, c1)], hermitian);
                const Val l2 = maybeConjugate(l[at(j, c2)], hermitian);
                u[atU(c1, j)] = d11 * l1 + d12 * l2;
                u[atU(c2, j)] = d21 * l1 + d22 * l2;
            }
            c += 2;
        }
    }
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
    formStaticUpper(n2, n1, a21, lda, a12, lda, a11, lda, hermitian);

    // A22 -= L21 U12. Symmetric, so only the lower triangle is filled.
    gemmLower(n2, n1, a21, lda, a12, lda, a22, lda);

    return ldl(n2, a22, lda, perturbation, numPerturbations, hermitian);
}

// Explicit instantiation. The headers declare; this file defines and instantiates, once. Adding a
// scalar type is one line per template.
template int  ldl<double>(int, double*, int, double, int*, bool);
template int  ldl<Cplx>(int, Cplx*, int, double, int*, bool);

template void formStaticUpper<double>(int, int, const double*, int, double*, int,
                                      const double*, int, bool);
template void formStaticUpper<Cplx>(int, int, const Cplx*, int, Cplx*, int,
                                    const Cplx*, int, bool);

template void formDynamicUpper<double>(int, int, const double*, int, double*, int,
                                       const double*, int, const std::int32_t*,
                                       const std::int32_t*, bool);
template void formDynamicUpper<Cplx>(int, int, const Cplx*, int, Cplx*, int,
                                     const Cplx*, int, const std::int32_t*,
                                     const std::int32_t*, bool);

template void gemmLower<double>(int, int, const double*, int, const double*, int, double*, int);
template void gemmLower<Cplx>(int, int, const Cplx*, int, const Cplx*, int, Cplx*, int);

} // namespace Oblio
