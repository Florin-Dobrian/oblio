// test_multiply.cpp — friend-access experiment.
//
// (1) Correctness: multiplyByApi and multiplyDirectly must agree with each other and
//     with hand-computed results, for real and complex.
// (2) Timing: on a large dense matrix, show the direct (friend) path is faster than
//     the API path. Indicative only — a single timed run, not a rigorous benchmark.
//
// Build & run:  make test   (produces ./test_multiply_cpp)

#include "MultiplyEngine.h"

#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <vector>

using Oblio::Matrix;
using Oblio::Vector;
using Oblio::MultiplyEngine;

namespace {

int gPass = 0;
int gFail = 0;

template<class Val>
bool close(Val a, Val b) { return std::abs(a - b) < 1e-9; }

template<class Val>
Vector<Val> makeVector(const std::vector<Val>& v) {
    Vector<Val> x(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) x[i] = v[i];
    return x;
}

template<class Val>
void check(const char* name,
           std::size_t rows, std::size_t cols,
           const std::vector<Val>& aVals,
           const std::vector<Val>& xVals,
           const std::vector<Val>& yExpect) {
    Matrix<Val> A(rows, cols, aVals);
    Vector<Val> x = makeVector(xVals);

    MultiplyEngine mul;
    Vector<Val> yApi  = mul.multiplyByApi(A, x);
    Vector<Val> yDir  = mul.multiplyDirectly(A, x);
    Vector<Val> yBlas = mul.multiplyWithBlas(A, x);

    bool ok = (yApi.size() == yExpect.size())
           && (yDir.size() == yExpect.size())
           && (yBlas.size() == yExpect.size());
    for (std::size_t i = 0; ok && i < yExpect.size(); ++i)
        ok = close(yApi[i], yExpect[i])
          && close(yDir[i], yExpect[i])
          && close(yBlas[i], yExpect[i]);

    std::cout << "  " << name << (ok ? "  PASS" : "  FAIL") << std::endl;
    if (ok) ++gPass; else ++gFail;
}

} // namespace

int main() {
    using cplx = std::complex<double>;

    std::cout << "=== friend-access matvec: correctness ===" << std::endl;

    // real 3x2: A=[[1,2],[3,4],[5,6]], x=[7,8] -> [23,53,83]
    check<double>("real 3x2   ", 3, 2, {1,2, 3,4, 5,6}, {7,8}, {23,53,83});
    // real 1x1
    check<double>("real 1x1   ", 1, 1, {5}, {3}, {15});
    // complex 2x2 (straight product, no conjugation)
    check<cplx>("complex 2x2", 2, 2,
        { cplx(1,1), cplx(2,0), cplx(0,1), cplx(1,-1) },
        { cplx(1,0), cplx(0,1) },
        { cplx(1,3), cplx(1,2) });

    std::cout << "Correctness: " << gPass << "/" << (gPass + gFail)
              << " passed\n" << std::endl;

    // --- timing (double, n x n dense) — indicative, not a rigorous benchmark ---
    const std::size_t n = 2000;
    std::cout << "=== timing (double, " << n << "x" << n << " dense) ===" << std::endl;

    std::vector<double> av(n * n), xv(n);
    for (std::size_t k = 0; k < n * n; ++k) av[k] = static_cast<double>((k % 7) + 1);
    for (std::size_t k = 0; k < n; ++k)     xv[k] = static_cast<double>((k % 5) + 1);
    Matrix<double> A(n, n, av);
    Vector<double> x = makeVector(xv);
    MultiplyEngine mul;

    auto timed = [&](const char* label, auto fn) {
        auto t0 = std::chrono::steady_clock::now();
        Vector<double> y = fn();
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        volatile double sink = y[0]; (void)sink;   // keep the result live
        std::cout << "  " << label << ": " << ms << " ms" << std::endl;
        return ms;
    };

    double tApi  = timed("multiplyByApi   ", [&] { return mul.multiplyByApi(A, x); });
    double tDir  = timed("multiplyDirectly", [&] { return mul.multiplyDirectly(A, x); });
    double tBlas = timed("multiplyWithBlas", [&] { return mul.multiplyWithBlas(A, x); });
    if (tDir  > 0.0) std::cout << "  speedup (api / direct): " << (tApi / tDir)  << "x" << std::endl;
    if (tBlas > 0.0) std::cout << "  speedup (api / blas)  : " << (tApi / tBlas) << "x" << std::endl;
    if (tBlas > 0.0) std::cout << "  speedup (direct / blas): " << (tDir / tBlas) << "x" << std::endl;

    return gFail == 0 ? 0 : 1;
}
