// examples/basic.cpp
// Minimal usage example: build a 4x4 tridiagonal matrix, factor it, solve, check the residual.
//
// This is the whole pipeline behind one object. For the same thing wired by hand, one engine at a
// time, see examples/pipeline.cpp.
//
// Compile (macOS, from repo root):
//   g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude examples/basic.cpp src/*.cpp -framework Accelerate -o basic
// Linux: replace `-framework Accelerate` with `-llapack -lblas`.

#include "oblio/DirectSolver.h"

#include <cstdio>
#include <vector>

using namespace Oblio;

int main() {
    // The 4x4 tridiagonal, diagonal 4, off-diagonal -1, stored full (both triangles) in CSC:
    //
    //   [ 4 -1  0  0 ]
    //   [-1  4 -1  0 ]
    //   [ 0 -1  4 -1 ]
    //   [ 0  0 -1  4 ]
    //
    const std::size_t n = 4;
    const std::vector<std::size_t>  colPtr = {0, 2, 5, 8, 10};
    const std::vector<std::int32_t> rowIdx = {0, 1,  0, 1, 2,  1, 2, 3,  2, 3};
    const std::vector<double>       val    = {4, -1, -1, 4, -1, -1, 4, -1, -1, 4};
    const SparseMatrix<double> A(n, colPtr, rowIdx, val);

    // Right-hand side: all ones.
    Vector<double> b(n);
    for (std::size_t i = 0; i < n; ++i)
        b[i] = 1.0;

    DirectSolver<double> solver(Factorization::Cholesky, Traversal::LeftLooking);
    solver.setOrderMethod(OrderMethod::AMD);

    // The three phases. Split like this because they have different lifetimes: analyze depends only
    // on the pattern, factor on the values, solve on the right-hand side.
    if (!solver.analyze(A)) { printf("analysis failed\n"); return 1; }
    if (!solver.factor(A))  { printf("factorization failed\n"); return 1; }

    Vector<double> x(n);
    if (!solver.solve(b, x)) { printf("solve failed\n"); return 1; }

    printf("Solution x:\n");
    for (std::size_t i = 0; i < n; ++i)
        printf("  x[%zu] = %.10f\n", i, x[i]);

    printf("\nRelative residual ||Ax - b|| / ||b|| = %.3e\n", solver.relativeResidual(A, b, x));

    // A second right-hand side reuses the factorization: no analysis, no factorization, just the
    // triangular solves. This is the reason the phases are separate.
    Vector<double> b2(n);
    for (std::size_t i = 0; i < n; ++i)
        b2[i] = static_cast<double>(i + 1);

    Vector<double> x2(n);
    if (!solver.solve(b2, x2)) { printf("second solve failed\n"); return 1; }

    printf("\nSecond right-hand side, same factorization:\n");
    for (std::size_t i = 0; i < n; ++i)
        printf("  x[%zu] = %.10f\n", i, x2[i]);
    printf("\nRelative residual ||Ax - b|| / ||b|| = %.3e\n", solver.relativeResidual(A, b2, x2));

    return 0;
}
