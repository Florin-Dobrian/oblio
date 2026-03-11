// test_complex_extended.cc — complex extended test suite for oblio.
//
// T1: Hermitian 2D Laplacian (Cholesky) — Natural, MMD, AMD
// T2: Symmetric 2D Laplacian (StaticLDL, DynamicLDL) — Natural, MMD, AMD
// T3: n=1 complex
// T4: Block diagonal complex
// T5: Near-singular complex (perturbation)

#include "oblio/OblioEngine.h"
#include "oblio/Vector.h"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

using namespace Oblio;
using C = std::complex<double>;

static int g_pass = 0, g_fail = 0;

// Hermitian residual: upper entry = conj(lower).
static double residualH(const Matrix<C>& A, const Vector<C>& x,
                         const Vector<C>& b) {
    Size n = A.mSize;
    std::vector<C> ax(n, C{0});
    for (Size j = 0; j < n; j++)
        for (Size p = A.mColPtr[j]; p < A.mColPtr[j+1]; p++) {
            Size i = A.mRowIdx[p];
            ax[i] += A.mVal[p] * x[j];
            if (i != j) ax[j] += std::conj(A.mVal[p]) * x[i];
        }
    double r = 0, nb = 0;
    for (Size i = 0; i < n; i++) {
        C d = ax[i] - b[i]; r += std::norm(d);
        nb += std::norm(b[i]);
    }
    return std::sqrt(r / (nb > 0 ? nb : 1));
}

// Symmetric residual: upper entry = lower (no conj).
static double residualS(const Matrix<C>& A, const Vector<C>& x,
                         const Vector<C>& b) {
    Size n = A.mSize;
    std::vector<C> ax(n, C{0});
    for (Size j = 0; j < n; j++)
        for (Size p = A.mColPtr[j]; p < A.mColPtr[j+1]; p++) {
            Size i = A.mRowIdx[p];
            ax[i] += A.mVal[p] * x[j];
            if (i != j) ax[j] += A.mVal[p] * x[i];
        }
    double r = 0, nb = 0;
    for (Size i = 0; i < n; i++) {
        C d = ax[i] - b[i]; r += std::norm(d);
        nb += std::norm(b[i]);
    }
    return std::sqrt(r / (nb > 0 ? nb : 1));
}

static void check(const char* label, double res, double tol = 1e-10) {
    bool ok = (res < tol);
    printf("  %-55s res=%.2e %s\n", label, res, ok ? "PASS" : "FAIL");
    if (ok) g_pass++; else g_fail++;
}

static const char* algName(FactorAlg a) {
    switch (a) {
        case FactorAlg::eLeftLooking:  return "LL";
        case FactorAlg::eRightLooking: return "RL";
        case FactorAlg::eMultifrontal: return "MF";
    }
    return "?";
}
static const char* typeName(FactorType t) {
    switch (t) {
        case FactorType::eCholesky:   return "CC";
        case FactorType::eStaticLDL:  return "SL";
        case FactorType::eDynamicLDL: return "DL";
    }
    return "?";
}
static const char* ordName(OrderAlg o) {
    switch (o) {
        case OrderAlg::eNatural: return "Nat";
        case OrderAlg::eMMD:     return "MMD";
        case OrderAlg::eAMD:     return "AMD";
    }
    return "?";
}

// Build Hermitian 2D Laplacian: diag=5 (real), off-diag = -1+0.2i in lower tri.
// Upper = conj = -1-0.2i.  HPD by diagonal dominance (5 > 4*|−1+0.2i| ≈ 4.08).
static Matrix<C> buildHermitianLap(int m) {
    int n = m * m;
    std::vector<Size> rows, cols;
    std::vector<C> vals;
    C offval(-1.0, 0.2);  // lower triangle value
    for (int j = 0; j < n; j++) {
        int r = j / m, c = j % m;
        rows.push_back(j); cols.push_back(j); vals.push_back(C{5.0, 0.0});
        // Lower triangle entries (row > col)
        if (c > 0)   { rows.push_back(j); cols.push_back(j-1); vals.push_back(offval); }
        if (r > 0)   { rows.push_back(j); cols.push_back(j-m); vals.push_back(offval); }
        if (c < m-1) { rows.push_back(j); cols.push_back(j+1); vals.push_back(offval); }
        if (r < m-1) { rows.push_back(j); cols.push_back(j+m); vals.push_back(offval); }
    }
    return Matrix<C>::fromCOO(n, rows, cols, vals);
}

// Build complex symmetric 2D Laplacian: diag=4+0.1i, off-diag = -1+0.3i.
// Symmetric (not Hermitian): A^T = A.
static Matrix<C> buildSymmetricLap(int m) {
    int n = m * m;
    std::vector<Size> rows, cols;
    std::vector<C> vals;
    C diagval(4.0, 0.1);
    C offval(-1.0, 0.3);
    for (int j = 0; j < n; j++) {
        int r = j / m, c = j % m;
        rows.push_back(j); cols.push_back(j); vals.push_back(diagval);
        if (c > 0)   { rows.push_back(j); cols.push_back(j-1); vals.push_back(offval); }
        if (r > 0)   { rows.push_back(j); cols.push_back(j-m); vals.push_back(offval); }
        if (c < m-1) { rows.push_back(j); cols.push_back(j+1); vals.push_back(offval); }
        if (r < m-1) { rows.push_back(j); cols.push_back(j+m); vals.push_back(offval); }
    }
    return Matrix<C>::fromCOO(n, rows, cols, vals);
}

int main() {
    FactorAlg algs[] = {FactorAlg::eLeftLooking, FactorAlg::eRightLooking,
                        FactorAlg::eMultifrontal};
    OrderAlg orders[] = {OrderAlg::eNatural, OrderAlg::eMMD, OrderAlg::eAMD};

    // ---- T1: Hermitian 2D Laplacian, Cholesky ----
    printf("=== T1: Hermitian 2D Laplacian (Cholesky) ===\n");
    {
        int sizes[] = {4, 5, 10};
        for (int si = 0; si < 3; si++) {
            int m = sizes[si], n = m * m;
            auto A = buildHermitianLap(m);
            Vector<C> b(n, C{1.0, 0.5}), x;

            for (int oi = 0; oi < 3; oi++)
                for (int ai = 0; ai < 3; ai++) {
                    OblioEngine<C> eng;
                    eng.setOrderAlg(orders[oi]);
                    eng.setFactorAlg(algs[ai]);
                    eng.setFactorType(FactorType::eCholesky);
                    Err e = eng.analyzeAndFactor(A);
                    char label[80];
                    snprintf(label, sizeof(label), "%dx%d(n=%d) %s %s CC",
                             m, m, n, ordName(orders[oi]), algName(algs[ai]));
                    if (e != Err::eNone) {
                        printf("  %-55s ERR=%d FAIL\n", label, (int)e);
                        g_fail++; continue;
                    }
                    eng.solve(b, x);
                    check(label, residualH(A, x, b));
                }
        }
    }

    // ---- T2: Symmetric 2D Laplacian, StaticLDL + DynamicLDL ----
    printf("\n=== T2: Symmetric 2D Laplacian (LDL) ===\n");
    {
        FactorType ldlTypes[] = {FactorType::eStaticLDL, FactorType::eDynamicLDL};
        const char* ldlNames[] = {"SL", "DL"};
        int sizes[] = {4, 5, 10};
        for (int si = 0; si < 3; si++) {
            int m = sizes[si], n = m * m;
            auto A = buildSymmetricLap(m);
            Vector<C> b(n, C{1.0, 0.5}), x;

            for (int oi = 0; oi < 3; oi++)
                for (int ai = 0; ai < 3; ai++)
                    for (int ti = 0; ti < 2; ti++) {
                        OblioEngine<C> eng;
                        eng.setOrderAlg(orders[oi]);
                        eng.setFactorAlg(algs[ai]);
                        eng.setFactorType(ldlTypes[ti]);
                        Err e = eng.analyzeAndFactor(A);
                        char label[80];
                        snprintf(label, sizeof(label), "%dx%d(n=%d) %s %s %s",
                                 m, m, n, ordName(orders[oi]), algName(algs[ai]),
                                 ldlNames[ti]);
                        if (e != Err::eNone) {
                            printf("  %-55s ERR=%d FAIL\n", label, (int)e);
                            g_fail++; continue;
                        }
                        eng.solve(b, x);
                        check(label, residualS(A, x, b));
                    }
        }
    }

    // ---- T3: n=1 complex ----
    printf("\n=== T3: n=1 complex ===\n");
    {
        std::vector<Size> r = {0}, c = {0};
        std::vector<C> v = {C{7.0, 0.0}};
        auto A = Matrix<C>::fromCOO(1, r, c, v);
        Vector<C> b(1, C{3.5, 1.0}), x;

        FactorType types[] = {FactorType::eCholesky, FactorType::eStaticLDL,
                              FactorType::eDynamicLDL};
        const char* tn[] = {"CC", "SL", "DL"};
        for (int ai = 0; ai < 3; ai++)
            for (int ti = 0; ti < 3; ti++) {
                OblioEngine<C> eng;
                eng.setOrderAlg(OrderAlg::eNatural);
                eng.setFactorAlg(algs[ai]);
                eng.setFactorType(types[ti]);
                eng.analyzeAndFactor(A);
                eng.solve(b, x);
                char label[80];
                snprintf(label, sizeof(label), "n=1 %s %s", algName(algs[ai]), tn[ti]);
                // For n=1, residual is trivial
                double res = std::abs(A.mVal[0] * x[0] - b[0]) / std::abs(b[0]);
                check(label, res);
            }
    }

    // ---- T4: Block diagonal complex ----
    printf("\n=== T4: Block diagonal complex (2 x 4-tridiag, n=8) ===\n");
    {
        int n = 8;
        std::vector<Size> rows, cols;
        std::vector<C> vals;
        // Block 1: Hermitian, diag=4, off=-1+0.2i
        C off1(-1.0, 0.2);
        for (int j = 0; j < 4; j++) {
            rows.push_back(j); cols.push_back(j); vals.push_back(C{4.0, 0.0});
            if (j > 0) { rows.push_back(j); cols.push_back(j-1); vals.push_back(off1); }
            if (j < 3) { rows.push_back(j); cols.push_back(j+1); vals.push_back(off1); }
        }
        // Block 2: Hermitian, diag=6, off=-1+0.2i
        for (int j = 4; j < 8; j++) {
            rows.push_back(j); cols.push_back(j); vals.push_back(C{6.0, 0.0});
            if (j > 4) { rows.push_back(j); cols.push_back(j-1); vals.push_back(off1); }
            if (j < 7) { rows.push_back(j); cols.push_back(j+1); vals.push_back(off1); }
        }
        auto A = Matrix<C>::fromCOO(n, rows, cols, vals);
        Vector<C> b(n, C{1.0, 0.5}), x;

        // Cholesky (Hermitian)
        for (int oi = 0; oi < 3; oi++)
            for (int ai = 0; ai < 3; ai++) {
                OblioEngine<C> eng;
                eng.setOrderAlg(orders[oi]);
                eng.setFactorAlg(algs[ai]);
                eng.setFactorType(FactorType::eCholesky);
                Err e = eng.analyzeAndFactor(A);
                char label[80];
                snprintf(label, sizeof(label), "BlkDiag %s %s CC",
                         ordName(orders[oi]), algName(algs[ai]));
                if (e != Err::eNone) {
                    printf("  %-55s ERR=%d FAIL\n", label, (int)e);
                    g_fail++; continue;
                }
                eng.solve(b, x);
                check(label, residualH(A, x, b));
            }
    }

    // ---- T5: Near-singular complex (perturbation) ----
    printf("\n=== T5: Near-singular complex (StaticLDL perturbation) ===\n");
    {
        int n = 3;
        std::vector<Size> r = {0, 1, 2};
        std::vector<Size> c = {0, 1, 2};
        std::vector<C> v = {C{4.0, 0.0}, C{1e-20, 0.0}, C{4.0, 0.0}};
        auto A = Matrix<C>::fromCOO(n, r, c, v);

        for (int ai = 0; ai < 3; ai++) {
            OblioEngine<C> eng;
            eng.setOrderAlg(OrderAlg::eNatural);
            eng.setFactorAlg(algs[ai]);
            eng.setFactorType(FactorType::eStaticLDL);
            Err e = eng.analyzeAndFactor(A);
            Size np = eng.nPert();
            char label[80];
            snprintf(label, sizeof(label), "NearSing %s SLDL nPert=%zu",
                     algName(algs[ai]), np);
            if (e != Err::eNone) {
                printf("  %-55s ERR=%d FAIL\n", label, (int)e);
                g_fail++; continue;
            }
            bool ok = (np > 0);
            printf("  %-55s %s\n", label, ok ? "PASS" : "FAIL");
            if (ok) g_pass++; else g_fail++;
        }
    }

    printf("\n========================================\n");
    printf("Total: %d/%d passed\n", g_pass, g_pass + g_fail);
    return g_fail > 0 ? 1 : 0;
}
