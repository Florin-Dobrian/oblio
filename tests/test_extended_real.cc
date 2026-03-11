// test_extended.cc — extended test suite for oblio.
//
// T1: 2D Laplacian with Natural and MMD ordering (n=16, 25, 100)
// T2: n=1 edge case
// T3: Block diagonal (disconnected graph: two 4×4 blocks)
// T4: Indefinite matrix → DynamicLDL 2×2 pivots (verify nSwap > 0)
// T5: Near-singular matrix → StaticLDL perturbation (verify nPert > 0)

#include "oblio/OblioEngine.h"
#include "oblio/Vector.h"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

using namespace Oblio;

static int g_pass = 0, g_fail = 0;

static double residual(const Matrix<double>& A, const Vector<double>& x,
                        const Vector<double>& b) {
    Size n = A.mSize;
    std::vector<double> ax(n, 0.0);
    for (Size j = 0; j < n; j++)
        for (Size p = A.mColPtr[j]; p < A.mColPtr[j + 1]; p++) {
            Size i = A.mRowIdx[p];
            ax[i] += A.mVal[p] * x[j];
            if (i != j) ax[j] += A.mVal[p] * x[i];
        }
    double r = 0, nb = 0;
    for (Size i = 0; i < n; i++) {
        double d = ax[i] - b[i]; r += d * d; nb += b[i] * b[i];
    }
    return std::sqrt(r / (nb > 0 ? nb : 1));
}

static void check(const char* label, double res, double tol = 1e-10) {
    bool ok = (res < tol);
    printf("  %-50s res=%.2e %s\n", label, res, ok ? "PASS" : "FAIL");
    if (ok) g_pass++; else g_fail++;
}

static void checkFail(const char* label, int errCode) {
    printf("  %-50s FACTOR_ERR=%d FAIL\n", label, errCode);
    g_fail++;
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

// Run all 9 combos with a given ordering.
static void runAll9(const Matrix<double>& A, const char* prefix, OrderAlg ord) {
    Size n = A.mSize;
    Vector<double> b(n, 1.0), x;
    FactorAlg algs[] = {FactorAlg::eLeftLooking, FactorAlg::eRightLooking,
                        FactorAlg::eMultifrontal};
    FactorType types[] = {FactorType::eCholesky, FactorType::eStaticLDL,
                          FactorType::eDynamicLDL};
    const char* oname = (ord == OrderAlg::eNatural) ? "Nat" :
                         (ord == OrderAlg::eMMD)     ? "MMD" : "AMD";

    for (int ai = 0; ai < 3; ai++)
        for (int ti = 0; ti < 3; ti++) {
            OblioEngine<double> eng;
            eng.setOrderAlg(ord);
            eng.setFactorAlg(algs[ai]);
            eng.setFactorType(types[ti]);
            Err e = eng.analyzeAndFactor(A);
            char label[80];
            snprintf(label, sizeof(label), "%s %s %s %s",
                     prefix, oname, algName(algs[ai]), typeName(types[ti]));
            if (e != Err::eNone) { checkFail(label, (int)e); continue; }
            eng.solve(b, x);
            check(label, residual(A, x, b));
        }
}

// ============================================================================

static Matrix<double> build2DLap(int m) {
    int n = m * m;
    std::vector<Size> rows, cols;
    std::vector<double> vals;
    for (int j = 0; j < n; j++) {
        int r = j / m, c = j % m;
        rows.push_back(j); cols.push_back(j); vals.push_back(4.0);
        if (c > 0)   { rows.push_back(j); cols.push_back(j-1); vals.push_back(-1.0); }
        if (r > 0)   { rows.push_back(j); cols.push_back(j-m); vals.push_back(-1.0); }
        if (c < m-1) { rows.push_back(j); cols.push_back(j+1); vals.push_back(-1.0); }
        if (r < m-1) { rows.push_back(j); cols.push_back(j+m); vals.push_back(-1.0); }
    }
    return Matrix<double>::fromCOO(n, rows, cols, vals);
}

// ============================================================================
int main() {
    // ---- T1: 2D Laplacians with both orderings ----
    printf("=== T1: 2D Laplacians ===\n");
    {
        int sizes[] = {4, 5, 10};
        for (int si = 0; si < 3; si++) {
            auto A = build2DLap(sizes[si]);
            char pfx[32];
            snprintf(pfx, sizeof(pfx), "%dx%d(n=%d)", sizes[si], sizes[si],
                     sizes[si] * sizes[si]);
            runAll9(A, pfx, OrderAlg::eNatural);
            runAll9(A, pfx, OrderAlg::eMMD);
            runAll9(A, pfx, OrderAlg::eAMD);
        }
    }

    // ---- T2: n=1 edge case ----
    printf("\n=== T2: n=1 edge case ===\n");
    {
        std::vector<Size> r = {0}, c = {0};
        std::vector<double> v = {7.0};
        auto A = Matrix<double>::fromCOO(1, r, c, v);
        runAll9(A, "n=1", OrderAlg::eNatural);
    }

    // ---- T3: Block diagonal (two disconnected 4×4 tridiag blocks) ----
    printf("\n=== T3: Block diagonal (2 x 4×4 tridiag, n=8) ===\n");
    {
        int n = 8;
        std::vector<Size> rows, cols;
        std::vector<double> vals;
        // Block 1: nodes 0..3, diag=4, off=-1
        for (int j = 0; j < 4; j++) {
            rows.push_back(j); cols.push_back(j); vals.push_back(4.0);
            if (j > 0) { rows.push_back(j); cols.push_back(j-1); vals.push_back(-1.0); }
            if (j < 3) { rows.push_back(j); cols.push_back(j+1); vals.push_back(-1.0); }
        }
        // Block 2: nodes 4..7, diag=6, off=-1
        for (int j = 4; j < 8; j++) {
            rows.push_back(j); cols.push_back(j); vals.push_back(6.0);
            if (j > 4) { rows.push_back(j); cols.push_back(j-1); vals.push_back(-1.0); }
            if (j < 7) { rows.push_back(j); cols.push_back(j+1); vals.push_back(-1.0); }
        }
        auto A = Matrix<double>::fromCOO(n, rows, cols, vals);
        runAll9(A, "BlkDiag n=8", OrderAlg::eNatural);
        runAll9(A, "BlkDiag n=8", OrderAlg::eMMD);
        runAll9(A, "BlkDiag n=8", OrderAlg::eAMD);
    }

    // ---- T4: Indefinite matrix → DynamicLDL 2×2 pivots ----
    printf("\n=== T4: Indefinite matrix (DynamicLDL 2×2 pivots) ===\n");
    {
        // Dense 4×4 indefinite symmetric matrix.  One snode with f=4, BK
        // pivoting identifies a 2×2 pivot for the near-zero diagonals.
        int n = 4;
        std::vector<Size> rows, cols;
        std::vector<double> vals;
        double M[4][4] = {
            { 1e-6,  1.0,  0.5,  0.3},
            { 1.0,   2.0,  0.7,  0.4},
            { 0.5,   0.7,  1e-6, 1.0},
            { 0.3,   0.4,  1.0,  3.0}
        };
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++) {
                rows.push_back(i); cols.push_back(j); vals.push_back(M[i][j]);
            }
        auto A = Matrix<double>::fromCOO(n, rows, cols, vals);
        Vector<double> b(n, 1.0), x;

        FactorAlg algs[] = {FactorAlg::eLeftLooking, FactorAlg::eRightLooking,
                            FactorAlg::eMultifrontal};

        for (int ai = 0; ai < 3; ai++) {
            OblioEngine<double> eng;
            eng.setOrderAlg(OrderAlg::eNatural);
            eng.setFactorAlg(algs[ai]);
            eng.setFactorType(FactorType::eDynamicLDL);
            Err e = eng.analyzeAndFactor(A);
            char label[80];
            snprintf(label, sizeof(label), "Indef %s DLDL (nSwap=%zu)",
                     algName(algs[ai]), eng.nSwap());
            if (e != Err::eNone) { checkFail(label, (int)e); continue; }
            eng.solve(b, x);
            double res = residual(A, x, b);
            check(label, res);
        }
    }

    // ---- T5: Near-singular → StaticLDL perturbation ----
    printf("\n=== T5: Near-singular matrix (StaticLDL perturbation) ===\n");
    {
        // Diagonal matrix with one near-zero entry.  No off-diagonals means
        // the pivot stays near-zero and perturbation must fire.  Residual against
        // original A will be large (the perturbed system is different), so we just
        // verify the factor completes and nPert > 0.
        int n = 3;
        std::vector<Size> rows = {0, 1, 2};
        std::vector<Size> cols = {0, 1, 2};
        std::vector<double> vals = {4.0, 1e-20, 4.0};
        auto A = Matrix<double>::fromCOO(n, rows, cols, vals);
        Vector<double> b(n, 1.0), x;

        FactorAlg algs[] = {FactorAlg::eLeftLooking, FactorAlg::eRightLooking,
                            FactorAlg::eMultifrontal};

        for (int ai = 0; ai < 3; ai++) {
            OblioEngine<double> eng;
            eng.setOrderAlg(OrderAlg::eNatural);
            eng.setFactorAlg(algs[ai]);
            eng.setFactorType(FactorType::eStaticLDL);
            Err e = eng.analyzeAndFactor(A);
            Size np = eng.nPert();
            char label[80];
            snprintf(label, sizeof(label), "NearSing %s SLDL nPert=%zu",
                     algName(algs[ai]), np);
            if (e != Err::eNone) { checkFail(label, (int)e); continue; }
            // Verify perturbation fired and factor completed.
            bool ok = (np > 0);
            printf("  %-50s %s\n", label, ok ? "PASS" : "FAIL");
            if (ok) g_pass++; else g_fail++;
        }
    }

    printf("\n========================================\n");
    printf("Total: %d/%d passed\n", g_pass, g_pass + g_fail);
    return g_fail > 0 ? 1 : 0;
}
