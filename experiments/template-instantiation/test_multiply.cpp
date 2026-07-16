// test_multiply.cpp - the Matrix / Vector / MultiplyEngine template-instantiation
// example (dense mat-vec, exercised across all three instantiation strategies).
//
// ONE source, THREE build configurations selected by macro (TI = template
// instantiation). The point of the example is that the variants differ only in HOW
// the code is compiled, never in WHAT it computes, so all builds must produce
// identical results.
//
//   -DOBLIO_TI_IMPLICIT        : implicit - body in header. Nothing to link.
//   -DOBLIO_TI_PLAIN_EXPLICIT  : plain explicit - bodies in .cpp, header signatures
//                                only. Link the three _PlainExplicit.cpp.
//   -DOBLIO_TI_GUARDED_EXPLICIT: guarded explicit - plain explicit + extern template
//                                guard. Link the three _GuardedExplicit.cpp.
//
// Build & run (from the example folder), e.g.:
//   implicit:        g++ -std=c++17 -DOBLIO_TI_IMPLICIT test_multiply.cpp -o test_multiply_implicit_cpp
//   plain explicit:  g++ -std=c++17 -DOBLIO_TI_PLAIN_EXPLICIT test_multiply.cpp
//                    MatrixPlainExplicit.cpp VectorPlainExplicit.cpp MultiplyEnginePlainExplicit.cpp
//                    -o test_multiply_plain_explicit_cpp
//   guarded explicit: g++ -std=c++17 -DOBLIO_TI_GUARDED_EXPLICIT test_multiply.cpp
//                    MatrixGuardedExplicit.cpp VectorGuardedExplicit.cpp MultiplyEngineGuardedExplicit.cpp
//                    -o test_multiply_guarded_explicit_cpp
// Or just: make test   (builds and runs all three)
//
// Bonus check: build the plain explicit or guarded explicit variant WITHOUT its
// .cpp files. Both must fail at link with undefined references, proof the bodies
// live only in the .o files.

#include <complex>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#if defined(OBLIO_TI_IMPLICIT)
#  include "MultiplyEngineImplicit.h"
#  define OBLIO_TI_VARIANT "implicit"
#elif defined(OBLIO_TI_PLAIN_EXPLICIT)
#  include "MultiplyEnginePlainExplicit.h"
#  define OBLIO_TI_VARIANT "plain explicit"
#elif defined(OBLIO_TI_GUARDED_EXPLICIT)
#  include "MultiplyEngineGuardedExplicit.h"
#  define OBLIO_TI_VARIANT "guarded explicit"
#else
#  error "Define OBLIO_TI_IMPLICIT, OBLIO_TI_PLAIN_EXPLICIT, or OBLIO_TI_GUARDED_EXPLICIT to select the example variant."
#endif

using Oblio::Matrix;
using Oblio::Vector;
using Oblio::MultiplyEngine;

namespace {

int gPass = 0;
int gFail = 0;

template<class Val>
bool close(Val a, Val b) {
    return std::abs(a - b) < 1e-12;
}

template<class Val>
Vector<Val> makeVector(const std::vector<Val>& v) {
    Vector<Val> x(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) x[i] = v[i];
    return x;
}

// Build A (row-major, rows*cols), multiply by x (cols), compare to yExpect (rows).
template<class Val>
void checkMultiply(const char* name,
                   std::size_t rows, std::size_t cols,
                   const std::vector<Val>& aVals,
                   const std::vector<Val>& xVals,
                   const std::vector<Val>& yExpect) {
    Matrix<Val> A(rows, cols, aVals);
    Vector<Val> x = makeVector(xVals);

    MultiplyEngine mul;
    Vector<Val> y = mul.Multiply(A, x);

    bool ok = (y.size() == yExpect.size());
    for (std::size_t i = 0; ok && i < yExpect.size(); ++i)
        ok = close(y[i], yExpect[i]);

    std::cout << "  [" OBLIO_TI_VARIANT "] " << name
              << (ok ? "  PASS" : "  FAIL") << std::endl;
    if (ok) ++gPass; else ++gFail;
}

} // namespace

int main() {
    using cplx = std::complex<double>;

    std::cout << "=== MultiplyEngine example (" OBLIO_TI_VARIANT ") ===" << std::endl;

    // Real, non-square, distinct entries, catches row/col and index slips.
    // A = [[1,2],[3,4],[5,6]] (3x2), x = [7,8]  ->  y = [23,53,83]
    checkMultiply<double>("real 3x2   ", 3, 2,
        {1,2, 3,4, 5,6}, {7,8}, {23,53,83});

    // Real 1x1 sanity.  [5] * [3] = [15]
    checkMultiply<double>("real 1x1   ", 1, 1, {5}, {3}, {15});

    // Complex 2x2 with nonzero imaginary parts. Verifies a STRAIGHT product (no
    // conjugation): a conj-transpose bug would change these values.
    // A = [[(1,1),(2,0)],[(0,1),(1,-1)]], x = [(1,0),(0,1)] -> y = [(1,3),(1,2)]
    checkMultiply<cplx>("complex 2x2", 2, 2,
        { cplx(1,1), cplx(2,0), cplx(0,1), cplx(1,-1) },
        { cplx(1,0), cplx(0,1) },
        { cplx(1,3), cplx(1,2) });

    std::cout << "Total (" OBLIO_TI_VARIANT "): "
              << gPass << "/" << (gPass + gFail) << " passed" << std::endl;

    return gFail == 0 ? 0 : 1;
}
