// test_smoke.cc — smoke test for all 9 factorization combinations.
// Uses a 5x5 tridiagonal SPD matrix: diag=4, off=-1.
// Tests both single-RHS (Vector) and multi-RHS (DenseMatrix) solves.

#include "oblio/OblioEngine.h"
#include "oblio/Vector.h"
#include "oblio/DenseMatrix.h"
#include <cmath>
#include <cstdio>
#include <vector>

using namespace Oblio;

static const char* algName(FactorAlg a) {
    switch(a) {
        case FactorAlg::eLeftLooking:  return "LeftLooking";
        case FactorAlg::eRightLooking: return "RightLooking";
        case FactorAlg::eMultifrontal: return "Multifrontal";
    } return "?";
}
static const char* typeName(FactorType t) {
    switch(t) {
        case FactorType::eCholesky:   return "Cholesky";
        case FactorType::eStaticLDL:  return "StaticLDL";
        case FactorType::eDynamicLDL: return "DynamicLDL";
    } return "?";
}

// Build n×n tridiagonal: diag=4, off=-1, lower triangle only.
static Matrix<double> buildTridiag(Size n) {
    std::vector<Size> rows, cols;
    std::vector<double> vals;
    for (Size i = 0; i < n; ++i) {
        rows.push_back(i); cols.push_back(i); vals.push_back(4.0);
        if (i+1 < n) { rows.push_back(i+1); cols.push_back(i); vals.push_back(-1.0); }
    }
    return Matrix<double>::fromCOO(n, rows, cols, vals);
}

// Compute ||Ax - b|| for symmetric CSC A.
static double residual(const Matrix<double>& A, const double* x, const double* b, Size n) {
    std::vector<double> ax(n, 0.0);
    for (Size j = 0; j < n; ++j)
        for (Size p = A.mColPtr[j]; p < A.mColPtr[j+1]; ++p) {
            Size i = A.mRowIdx[p];
            ax[i] += A.mVal[p] * x[j];
            if (i != j) ax[j] += A.mVal[p] * x[i];
        }
    double r = 0.0;
    for (Size i = 0; i < n; ++i) { double d = ax[i]-b[i]; r += d*d; }
    return std::sqrt(r);
}

int main() {
    const Size n = 5;
    auto A = buildTridiag(n);

    FactorAlg  algs[]  = { FactorAlg::eLeftLooking, FactorAlg::eRightLooking, FactorAlg::eMultifrontal };
    FactorType types[] = { FactorType::eCholesky,   FactorType::eStaticLDL,   FactorType::eDynamicLDL  };

    int pass = 0, total = 0;

    // ---- Single-RHS (Vector) ------------------------------------------------
    printf("=== Single RHS (Vector) ===\n");
    for (auto alg : algs) {
        for (auto typ : types) {
            ++total;
            OblioEngine<double> eng;
            eng.setOrderAlg(OrderAlg::eMMD);
            eng.setFactorAlg(alg);
            eng.setFactorType(typ);

            Err e = eng.analyzeAndFactor(A);
            if (e != Err::eNone) {
                printf("FAIL [%s x %s]: factor error %d\n", algName(alg), typeName(typ), (int)e);
                continue;
            }

            Vector<double> b(n, 1.0), x;
            e = eng.solve(b, x);
            if (e != Err::eNone) {
                printf("FAIL [%s x %s]: solve error %d\n", algName(alg), typeName(typ), (int)e);
                continue;
            }

            double res = residual(A, x.data(), b.data(), n);
            if (res < 1e-10) {
                printf("PASS [%s x %s]: residual=%.2e\n", algName(alg), typeName(typ), res);
                ++pass;
            } else {
                printf("FAIL [%s x %s]: residual=%.2e\n", algName(alg), typeName(typ), res);
            }
        }
    }

    // ---- Multi-RHS (DenseMatrix) --------------------------------------------
    printf("\n=== Multi-RHS (DenseMatrix, nRHS=3) ===\n");
    // Three RHS: all-ones, all-twos, identity-like (b[i]=i+1)
    const Size nRHS = 3;
    DenseMatrix<double> B(n, nRHS);
    for (Size i = 0; i < n; ++i) {
        B(i,0) = 1.0;
        B(i,1) = 2.0;
        B(i,2) = (double)(i+1);
    }

    for (auto alg : algs) {
        for (auto typ : types) {
            ++total;
            OblioEngine<double> eng;
            eng.setOrderAlg(OrderAlg::eMMD);
            eng.setFactorAlg(alg);
            eng.setFactorType(typ);

            Err e = eng.analyzeAndFactor(A);
            if (e != Err::eNone) {
                printf("FAIL [%s x %s]: factor error %d\n", algName(alg), typeName(typ), (int)e);
                continue;
            }

            DenseMatrix<double> X;
            e = eng.solve(B, X);
            if (e != Err::eNone) {
                printf("FAIL [%s x %s]: solve error %d\n", algName(alg), typeName(typ), (int)e);
                continue;
            }

            // Check each column
            double maxRes = 0.0;
            for (Size r = 0; r < nRHS; ++r) {
                std::vector<double> bcol(n), xcol(n);
                for (Size i = 0; i < n; ++i) { bcol[i]=B(i,r); xcol[i]=X(i,r); }
                double res = residual(A, xcol.data(), bcol.data(), n);
                if (res > maxRes) maxRes = res;
            }
            if (maxRes < 1e-10) {
                printf("PASS [%s x %s]: max_residual=%.2e\n", algName(alg), typeName(typ), maxRes);
                ++pass;
            } else {
                printf("FAIL [%s x %s]: max_residual=%.2e\n", algName(alg), typeName(typ), maxRes);
            }
        }
    }

    printf("\n%d/%d passed\n", pass, total);
    return (pass == total) ? 0 : 1;
}
