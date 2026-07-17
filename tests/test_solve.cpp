// test_solve.cpp - solve A x = b, and check the residual.
//
// **This is the first test that exercises the pipeline as a pipeline.** Everything before it checks
// one phase against an oracle: the forest against a recomputation, the symbolic factor against a
// dense pattern, the numeric factor against a dense Cholesky or by reconstruction. Each says "this
// phase computed what it should". None says the phases *compose*.
//
//     || A x - b ||  /  || b ||
//
// says that, in one number, through ordering, elimination forest, symbolic factorization, numeric
// factorization, triangular solve, and sparse matvec. Six phases. If any of them is subtly wrong,
// or if any two disagree about a convention (an ordering, a conjugate, an index base), this is
// where it shows.
//
// The right-hand side is manufactured from a known solution: pick x, form b = A x, solve, and see
// whether x comes back. That way the test needs no reference solver.

#include "oblio/ElmForestEngine.h"
#include "oblio/MultiplyEngine.h"
#include "oblio/NumFactorEngine.h"
#include "oblio/OrderEngine.h"
#include "oblio/SolveEngine.h"
#include "oblio/SymFactorEngine.h"

#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

using namespace Oblio;
using Cplx = std::complex<double>;

namespace {

int pass = 0, fail = 0;
void ck(bool ok, const std::string& what) {
    if (ok) { ++pass; std::cout << "  PASS  " << what << "\n"; }
    else    { ++fail; std::cout << "  FAIL  " << what << "\n"; }
}

double maybeConj(double v, bool)          { return v; }
Cplx   maybeConj(Cplx v, bool hermitian)  { return hermitian ? std::conj(v) : v; }

// A random matrix, Hermitian or complex-symmetric as asked, diagonally dominant so that it is
// definite (for Cholesky) and its pivots stay away from zero (for LDL).
template<class Val>
std::vector<std::vector<Val>> randomMatrix(std::size_t n, std::mt19937& rng, bool hermitian) {
    std::vector<std::vector<Val>> A(n, std::vector<Val>(n, Val(0)));
    std::uniform_real_distribution<double> u(-1.0, 1.0);

    auto entry = [&]() -> Val {
        if constexpr (std::is_same_v<Val, double>) return Val(u(rng));
        else                                       return Val(u(rng), u(rng));
    };

    for (std::size_t i = 1; i < n; ++i) {
        const Val v = entry();
        A[i][i - 1] = v;
        A[i - 1][i] = maybeConj(v, hermitian);
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 2; j < n; ++j)
            if (rng() % 100 < 25) {
                const Val v = entry();
                A[j][i] = v;
                A[i][j] = maybeConj(v, hermitian);
            }
    for (std::size_t i = 0; i < n; ++i) {
        double off = 0;
        for (std::size_t j = 0; j < n; ++j)
            if (i != j) off += std::abs(A[i][j]);
        A[i][i] = Val(off + 1.0);
    }
    return A;
}

template<class Val>
SparseMatrix<Val> toSparse(const std::vector<std::vector<Val>>& A) {
    const std::size_t n = A.size();
    std::vector<std::size_t>  colPtr(n + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<Val>          val;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i)
            if (std::abs(A[i][j]) != 0.0) {
                rowIdx.push_back(static_cast<std::int32_t>(i));
                val.push_back(A[i][j]);
            }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<Val>(n, colPtr, rowIdx, val);
}

// The whole pipeline, and the one number that judges it.
//
// Manufacture b from a known x, solve, and measure ||A x - b|| / ||b||. Note this checks the
// *residual*, not the distance to the manufactured x: they differ by the conditioning of A, and
// the residual is the honest thing to require of a direct solver.
template<class Val, class Factor = NumFactorStatic<Val>>
double solveResidual(const std::vector<std::vector<Val>>& dense, std::mt19937& rng,
                     Factorization factorization, Traversal traversal, int& failures) {
    const std::size_t n = dense.size();
    const SparseMatrix<Val> A = toSparse(dense);

    OrderEngine ord(OrderMethod::AMD);
    Permutation p;
    if (!ord.compute(A, p)) { ++failures; return -1; }

    ElmForest f;
    ElmForestEngine fe;
    if (!fe.compute(A, p, f)) { ++failures; return -1; }

    SymFactor s;
    SymFactorEngine se;
    if (!se.compute(A, p, f, s)) { ++failures; return -1; }

    Factor nf;
    NumFactorEngine ne(factorization, traversal);
    if (!ne.compute(A, p, s, nf)) { ++failures; return -1; }

    // A known solution, and the right-hand side it implies.
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    Vector<Val> xTrue(n);
    for (std::size_t i = 0; i < n; ++i) {
        if constexpr (std::is_same_v<Val, double>) xTrue[i] = Val(u(rng));
        else                                       xTrue[i] = Val(u(rng), u(rng));
    }

    MultiplyEngine mul;
    Vector<Val> b(n);
    if (!mul.compute(A, xTrue, b)) { ++failures; return -1; }

    SolveEngine sol;
    Vector<Val> x(n);
    if (!sol.compute(p, nf, b, x)) { ++failures; return -1; }

    Vector<Val> r(n);
    if (!mul.residual(A, x, b, r)) { ++failures; return -1; }

    const double normB = b.norm();
    return normB > 0 ? r.norm() / normB : r.norm();
}

// The whole pipeline, swept over every factorization and both traversals, into one storage. Run
// once per factor type: the point of the storage split is that the solve does not care which
// factor it is handed, so the residual must reach machine precision either way.
template<template<class> class FactorT>
void pipelineSweep(std::mt19937& rng, double tol, const std::string& storage) {
    int failures = 0;
    double wRC = 0, wRT = 0, wCC = 0, wCT = 0, wCH = 0;

    for (int trial = 0; trial < 40; ++trial) {
        const std::size_t n = 5 + rng() % 15;

        const auto herm  = randomMatrix<double>(n, rng, true);    // real: symmetric = Hermitian
        const auto cherm = randomMatrix<Cplx>(n, rng, true);      // complex Hermitian
        const auto csym  = randomMatrix<Cplx>(n, rng, false);     // complex symmetric

        for (Traversal tr : {Traversal::LeftLooking, Traversal::RightLooking}) {
            wRC = std::max(wRC, solveResidual<double, FactorT<double>>(herm,  rng, Factorization::Cholesky,   tr, failures));
            wRT = std::max(wRT, solveResidual<double, FactorT<double>>(herm,  rng, Factorization::StaticLDLT, tr, failures));
            wCC = std::max(wCC, solveResidual<Cplx,   FactorT<Cplx>>  (cherm, rng, Factorization::Cholesky,   tr, failures));
            wCH = std::max(wCH, solveResidual<Cplx,   FactorT<Cplx>>  (cherm, rng, Factorization::StaticLDLH, tr, failures));
            wCT = std::max(wCT, solveResidual<Cplx,   FactorT<Cplx>>  (csym,  rng, Factorization::StaticLDLT, tr, failures));
        }
    }

    ck(failures == 0, storage + " x40 x10  : every system solved");
    ck(wRC < tol, storage + " real    Cholesky : ||Ax - b|| / ||b|| at machine precision");
    ck(wRT < tol, storage + " real    LDLT     : ||Ax - b|| / ||b|| at machine precision");

    // The complex Hermitian ones are where a missing conjugate would show, and 10.12's solve omits
    // it: its backward pass applies L^-T where a Hermitian factor needs L^-H.
    ck(wCC < tol, storage + " complex Cholesky : ||Ax - b|| / ||b|| at machine precision (Hermitian)");
    ck(wCH < tol, storage + " complex LDLH     : ||Ax - b|| / ||b|| at machine precision (Hermitian)");
    ck(wCT < tol, storage + " complex LDLT     : ||Ax - b|| / ||b|| at machine precision (symmetric)");
}

} // namespace

int main() {
    std::mt19937 rng(20260713);
    const double tol = 1e-10;

    // MultiplyEngine on its own: y = A x, against a hand-computed answer.
    {
        std::vector<std::vector<double>> A = {{2, 1, 0}, {1, 3, 1}, {0, 1, 4}};
        const SparseMatrix<double> S = toSparse(A);
        Vector<double> x(3), y(3);
        x[0] = 1; x[1] = 2; x[2] = 3;
        MultiplyEngine mul;
        const bool ok = mul.compute(S, x, y);
        // A x = (2+2, 1+6+3, 2+12) = (4, 10, 14)
        ck(ok && std::abs(y[0] - 4) < 1e-12 && std::abs(y[1] - 10) < 1e-12
              && std::abs(y[2] - 14) < 1e-12,
           "multiply 3x3        : y = A x matches by hand");
    }

    // The pipeline, every factorization and both traversals, run into each factor storage in turn.
    pipelineSweep<NumFactorStatic> (rng, tol, "static ");
    pipelineSweep<NumFactorDynamic>(rng, tol, "dynamic");

    // A grid, which is where the structure is real: deep forest, genuine fill, an ordering that
    // matters. The residual has to survive all of it.
    {
        const std::size_t g = 10, n = g * g;
        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        for (std::size_t y = 0; y < g; ++y)
            for (std::size_t x = 0; x < g; ++x) {
                const std::size_t j = y * g + x;
                A[j][j] = 4.5;
                if (x + 1 < g) { A[j][j + 1] = -1.0; A[j + 1][j] = -1.0; }
                if (y + 1 < g) { A[j][j + g] = -1.0; A[j + g][j] = -1.0; }
            }
        int failures = 0;
        const double rcS = solveResidual<double, NumFactorStatic<double>> (A, rng, Factorization::Cholesky,   Traversal::LeftLooking,  failures);
        const double rtS = solveResidual<double, NumFactorStatic<double>> (A, rng, Factorization::StaticLDLT, Traversal::RightLooking, failures);
        const double rcD = solveResidual<double, NumFactorDynamic<double>>(A, rng, Factorization::Cholesky,   Traversal::LeftLooking,  failures);
        const double rtD = solveResidual<double, NumFactorDynamic<double>>(A, rng, Factorization::StaticLDLT, Traversal::RightLooking, failures);
        ck(failures == 0 && rcS < tol && rtS < tol && rcD < tol && rtD < tol,
           "10x10 grid          : ||Ax - b|| / ||b|| at machine precision, static and dynamic");
    }

    std::cout << "\nSolve tests: " << pass << "/" << (pass + fail) << " passed\n";
    return fail == 0 ? 0 : 1;
}
