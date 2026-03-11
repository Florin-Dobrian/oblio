// test_complex.cc — verify complex<double> factorization and solve.
//
// Two classes of test matrices:
//   1. Hermitian (A = A^H): diagonal real, off-diag conj-symmetric.
//      Used for Cholesky (zpotrf assumes A = L*L^H).
//   2. Complex symmetric (A = A^T): diagonal real, off-diag symmetric.
//      Used for LDL variants (no Hermitian assumption).
//
// Residual uses the appropriate symmetry for each case.
//
// Tests all 9 combinations for both Vector and DenseMatrix (2 RHS) solves.

#include "oblio/OblioEngine.h"
#include "oblio/Vector.h"
#include "oblio/DenseMatrix.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

using namespace Oblio;
using C = std::complex<double>;

static const char* algName(FactorAlg a) {
    switch(a) {
        case FactorAlg::eLeftLooking:  return "LeftLooking ";
        case FactorAlg::eRightLooking: return "RightLooking";
        case FactorAlg::eMultifrontal: return "Multifrontal";
    } return "?";
}
static const char* typeName(FactorType t) {
    switch(t) {
        case FactorType::eCholesky:   return "Cholesky  ";
        case FactorType::eStaticLDL:  return "StaticLDL ";
        case FactorType::eDynamicLDL: return "DynamicLDL";
    } return "?";
}

// Build n×n tridiagonal (lower triangle only).
// diag = diagVal, sub-diag = offVal.
// For Hermitian tests: offVal in lower triangle, conj(offVal) in upper.
// For symmetric tests: offVal both.
// The matrix we store is always lower-triangle CSC regardless.
static Matrix<C> buildTridiag(Size n, C diagVal, C offVal) {
    std::vector<Size> rows, cols;
    std::vector<C> vals;
    for (Size i = 0; i < n; ++i) {
        rows.push_back(i); cols.push_back(i); vals.push_back(diagVal);
        if (i+1 < n) {
            rows.push_back(i+1); cols.push_back(i); vals.push_back(offVal);
        }
    }
    return Matrix<C>::fromCOO(n, rows, cols, vals);
}

// Residual for Hermitian matrix: A[i,j] = conj(A[j,i]) for i<j.
// Lower-triangle CSC stores A[i,j] for i>=j; upper entry = conj.
static double residualH(const Matrix<C>& A, const C* x, const C* b, Size n) {
    std::vector<C> ax(n, C{0});
    for (Size j = 0; j < n; ++j)
        for (Size p = A.mColPtr[j]; p < A.mColPtr[j+1]; ++p) {
            Size i = A.mRowIdx[p];
            ax[i] += A.mVal[p] * x[j];
            if (i != j) ax[j] += std::conj(A.mVal[p]) * x[i]; // Hermitian upper = conj
        }
    double r = 0.0;
    for (Size i = 0; i < n; ++i) { C d = ax[i]-b[i]; r += std::norm(d); }
    return std::sqrt(r);
}

// Residual for complex symmetric matrix: A[i,j] = A[j,i] (no conj).
static double residualS(const Matrix<C>& A, const C* x, const C* b, Size n) {
    std::vector<C> ax(n, C{0});
    for (Size j = 0; j < n; ++j)
        for (Size p = A.mColPtr[j]; p < A.mColPtr[j+1]; ++p) {
            Size i = A.mRowIdx[p];
            ax[i] += A.mVal[p] * x[j];
            if (i != j) ax[j] += A.mVal[p] * x[i]; // symmetric upper = same value
        }
    double r = 0.0;
    for (Size i = 0; i < n; ++i) { C d = ax[i]-b[i]; r += std::norm(d); }
    return std::sqrt(r);
}

int main() {
    const Size n = 6;
    int pass = 0, total = 0;

    // ---- Cholesky: Hermitian matrix ----------------------------------------
    // diag=8, off=-1+0.5i in lower triangle (so upper = conj = -1-0.5i).
    // A is Hermitian positive definite (diag dominant: 8 > 2*sqrt(1^2+0.25^2) ≈ 2.06).
    auto AH = buildTridiag(n, C{8,0}, C{-1,0.5});

    printf("=== Cholesky (Hermitian matrix, single RHS) ===\n");
    for (auto alg : {FactorAlg::eLeftLooking, FactorAlg::eRightLooking, FactorAlg::eMultifrontal}) {
        ++total;
        OblioEngine<C> eng;
        eng.setOrderAlg(OrderAlg::eMMD);
        eng.setFactorAlg(alg);
        eng.setFactorType(FactorType::eCholesky);
        Err e = eng.analyzeAndFactor(AH);
        if (e != Err::eNone) {
            printf("FAIL [%s x Cholesky  ]: factor error %d\n", algName(alg), (int)e);
            continue;
        }
        Vector<C> b(n), x;
        for (Size i = 0; i < n; ++i) b[i] = C{1.0, 0.5*(double)i};
        eng.solve(b, x);
        double res = residualH(AH, x.data(), b.data(), n);
        bool ok = res < 1e-10;
        printf("%s [%s x Cholesky  ]: residual=%.2e\n", ok?"PASS":"FAIL", algName(alg), res);
        if (ok) ++pass;
    }

    printf("\n=== LDL variants (complex symmetric matrix, single RHS) ===\n");
    // For LDL: use a complex symmetric diagonally-dominant matrix.
    // diag=8+2i (complex diagonal — symmetric, not Hermitian), off=-1+0.5i.
    // Diagonal dominance: |8+2i|=8.25 > 2*|{-1,0.5}|≈2.24. StaticLDL/DynamicLDL ok.
    auto AS = buildTridiag(n, C{8,2}, C{-1,0.5});

    FactorType ldlTypes[] = { FactorType::eStaticLDL, FactorType::eDynamicLDL };
    for (auto alg : {FactorAlg::eLeftLooking, FactorAlg::eRightLooking, FactorAlg::eMultifrontal}) {
        for (auto typ : ldlTypes) {
            ++total;
            OblioEngine<C> eng;
            eng.setOrderAlg(OrderAlg::eMMD);
            eng.setFactorAlg(alg);
            eng.setFactorType(typ);
            Err e = eng.analyzeAndFactor(AS);
            if (e != Err::eNone) {
                printf("FAIL [%s x %s]: factor error %d\n", algName(alg), typeName(typ), (int)e);
                continue;
            }
            Vector<C> b(n), x;
            for (Size i = 0; i < n; ++i) b[i] = C{1.0, 0.5*(double)i};
            eng.solve(b, x);
            double res = residualS(AS, x.data(), b.data(), n);
            bool ok = res < 1e-10;
            printf("%s [%s x %s]: residual=%.2e\n", ok?"PASS":"FAIL", algName(alg), typeName(typ), res);
            if (ok) ++pass;
        }
    }

    // ---- Multi-RHS (DenseMatrix) -------------------------------------------
    printf("\n=== Cholesky (Hermitian, DenseMatrix nRHS=2) ===\n");
    DenseMatrix<C> BH(n, 2);
    for (Size i = 0; i < n; ++i) {
        BH(i,0) = C{1.0, 0.5*(double)i};
        BH(i,1) = C{(double)(i+1), -(double)i};
    }
    for (auto alg : {FactorAlg::eLeftLooking, FactorAlg::eRightLooking, FactorAlg::eMultifrontal}) {
        ++total;
        OblioEngine<C> eng;
        eng.setOrderAlg(OrderAlg::eMMD);
        eng.setFactorAlg(alg);
        eng.setFactorType(FactorType::eCholesky);
        Err e = eng.analyzeAndFactor(AH);
        if (e != Err::eNone) { printf("FAIL [%s x Cholesky  ]: factor error %d\n", algName(alg),(int)e); continue; }
        DenseMatrix<C> X;
        eng.solve(BH, X);
        double maxRes = 0;
        for (Size r = 0; r < 2; ++r) {
            std::vector<C> bcol(n), xcol(n);
            for (Size i=0;i<n;++i){bcol[i]=BH(i,r);xcol[i]=X(i,r);}
            double res=residualH(AH,xcol.data(),bcol.data(),n);
            if(res>maxRes) maxRes=res;
        }
        bool ok=maxRes<1e-10;
        printf("%s [%s x Cholesky  ]: max_residual=%.2e\n",ok?"PASS":"FAIL",algName(alg),maxRes);
        if(ok) ++pass;
    }

    printf("\n=== LDL variants (complex symmetric, DenseMatrix nRHS=2) ===\n");
    DenseMatrix<C> BS(n, 2);
    for (Size i = 0; i < n; ++i) {
        BS(i,0) = C{1.0, 0.5*(double)i};
        BS(i,1) = C{(double)(i+1), -(double)i};
    }
    for (auto alg : {FactorAlg::eLeftLooking, FactorAlg::eRightLooking, FactorAlg::eMultifrontal}) {
        for (auto typ : ldlTypes) {
            ++total;
            OblioEngine<C> eng;
            eng.setOrderAlg(OrderAlg::eMMD);
            eng.setFactorAlg(alg);
            eng.setFactorType(typ);
            Err e = eng.analyzeAndFactor(AS);
            if (e != Err::eNone) { printf("FAIL [%s x %s]: factor error %d\n",algName(alg),typeName(typ),(int)e); continue; }
            DenseMatrix<C> X;
            eng.solve(BS, X);
            double maxRes=0;
            for (Size r=0;r<2;++r){
                std::vector<C> bcol(n),xcol(n);
                for(Size i=0;i<n;++i){bcol[i]=BS(i,r);xcol[i]=X(i,r);}
                double res=residualS(AS,xcol.data(),bcol.data(),n);
                if(res>maxRes) maxRes=res;
            }
            bool ok=maxRes<1e-10;
            printf("%s [%s x %s]: max_residual=%.2e\n",ok?"PASS":"FAIL",algName(alg),typeName(typ),maxRes);
            if(ok) ++pass;
        }
    }

    printf("\n%d/%d passed\n", pass, total);
    return (pass == total) ? 0 : 1;
}
