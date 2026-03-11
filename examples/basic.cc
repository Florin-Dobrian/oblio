// examples/basic.cc
// Minimal usage example: build a 4x4 tridiagonal matrix, factor it, solve.
//
// Compile (from repo root, after cmake build):
//   g++ -std=c++17 -I include examples/basic.cc -L build -loblio -lblas -llapack -o basic

#include "oblio/OblioEngine.h"
#include <cstdio>
#include <vector>

int main() {
    using namespace Oblio;

    // Build 4x4 tridiagonal: diag=4, off=-1 (lower triangle only, COO format).
    //
    //   [ 4 -1  0  0 ]
    //   [-1  4 -1  0 ]
    //   [ 0 -1  4 -1 ]
    //   [ 0  0 -1  4 ]

    const Size n = 4;
    std::vector<Size> rows = {0, 1, 1, 2, 2, 3, 3};
    std::vector<Size> cols = {0, 0, 1, 1, 2, 2, 3};
    std::vector<double> vals = {4, -1, 4, -1, 4, -1, 4};

    auto A = Matrix<double>::fromCOO(n, rows, cols, vals);

    // Right-hand side: all ones.
    std::vector<double> b(n, 1.0);

    // Set up solver: MMD ordering, multifrontal, Cholesky.
    OblioEngine<double> eng;
    eng.setOrderAlg(OrderAlg::eMMD);
    eng.setFactorAlg(FactorAlg::eMultifrontal);
    eng.setFactorType(FactorType::eCholesky);

    // Analyze + factor.
    Err e = eng.analyzeAndFactor(A);
    if (e != Err::eNone) {
        printf("Factor failed: %d\n", (int)e);
        return 1;
    }

    // Solve A*x = b.
    std::vector<double> x;
    e = eng.solve(b, x);
    if (e != Err::eNone) {
        printf("Solve failed: %d\n", (int)e);
        return 1;
    }

    // Print solution.
    printf("Solution x:\n");
    for (Size i = 0; i < n; ++i)
        printf("  x[%zu] = %.10f\n", i, x[i]);

    // Compute residual ||Ax - b||.
    double res = 0.0;
    for (Size j = 0; j < n; ++j) {
        double ax = 0.0;
        for (Size p = A.mColPtr[j]; p < A.mColPtr[j+1]; ++p) {
            Size i = A.mRowIdx[p];
            ax += A.mVal[p] * x[i];
            if (i != j) ax += A.mVal[p] * x[j]; // symmetric off-diagonal
        }
        double r = ax - b[j];
        res += r * r;
    }
    printf("Residual ||Ax-b|| = %.3e\n", std::sqrt(res));

    return 0;
}
