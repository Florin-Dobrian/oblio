// test_multiply.cc — exercises the Matrix / Vector / MultiplyEngine spike.
//
// ONE source, THREE build configurations selected by macro. The point of the
// spike is that the variants differ only in HOW the code is compiled, never in
// WHAT it computes — so all builds must produce identical results.
//
//   -DOBLIO_SPIKE_TPL : Case 1 — header-only template inclusion. Nothing to link.
//   -DOBLIO_SPIKE_EXP : Case 2 — explicit instantiation, forcing only. Link _exp.cc.
//   -DOBLIO_SPIKE_EXT : Case 3 — explicit instantiation + extern template. Link _ext.cc.
//
// Build & run (from archive/), e.g.:
//   tpl:  g++ -std=c++17 -DOBLIO_SPIKE_TPL test_multiply.cc -o test_multiply_tpl
//   exp:  g++ -std=c++17 -DOBLIO_SPIKE_EXP test_multiply.cc
//         Matrix_exp.cc Vector_exp.cc MultiplyEngine_exp.cc -o test_multiply_exp
//   ext:  g++ -std=c++17 -DOBLIO_SPIKE_EXT test_multiply.cc
//         Matrix_ext.cc Vector_ext.cc MultiplyEngine_ext.cc -o test_multiply_ext
// Or just: make test   (builds and runs all three)
//
// Bonus check: build exp or ext WITHOUT their .cc files. Both must fail at link
// with undefined references — proof the bodies live only in the .o files.

#include <complex>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#if defined(OBLIO_SPIKE_TPL)
#  include "MultiplyEngine_tpl.h"
#  define OBLIO_SPIKE_VARIANT "tpl"
#elif defined(OBLIO_SPIKE_EXP)
#  include "MultiplyEngine_exp.h"
#  define OBLIO_SPIKE_VARIANT "exp"
#elif defined(OBLIO_SPIKE_EXT)
#  include "MultiplyEngine_ext.h"
#  define OBLIO_SPIKE_VARIANT "ext"
#else
#  error "Define OBLIO_SPIKE_TPL, OBLIO_SPIKE_EXP, or OBLIO_SPIKE_EXT to select the spike variant."
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

    std::cout << "  [" OBLIO_SPIKE_VARIANT "] " << name
              << (ok ? "  PASS" : "  FAIL") << std::endl;
    if (ok) ++gPass; else ++gFail;
}

} // namespace

int main() {
    using cplx = std::complex<double>;

    std::cout << "=== MultiplyEngine spike (" OBLIO_SPIKE_VARIANT ") ===" << std::endl;

    // Real, non-square, distinct entries — catches row/col and index slips.
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

    std::cout << "Total (" OBLIO_SPIKE_VARIANT "): "
              << gPass << "/" << (gPass + gFail) << " passed" << std::endl;

    return gFail == 0 ? 0 : 1;
}
