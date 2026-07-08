// MultiplyEngine_exp.cpp — Case 2: explicit instantiation (forcing only)
//
// The Multiply function body is defined here and compiled once per scalar type.
// Other translation units that call Multiply<double>() or Multiply<complex<double>>()
// see only the declaration in the header, so they cannot implicitly instantiate —
// they resolve the symbol from this object file at link time.

#include "MultiplyEngine_exp.h"
#include <cassert>

namespace Oblio {

template<class Val>
Vector<Val> MultiplyEngine::Multiply(const Matrix<Val>& A,
                                     const Vector<Val>& x) const {
    assert(A.cols() == x.size());
    Vector<Val> y(A.rows());
    for (std::size_t i = 0; i < A.rows(); ++i) {
        Val sum{0};
        for (std::size_t j = 0; j < A.cols(); ++j) {
            sum += A(i, j) * x[j];
        }
        y[i] = sum;
    }
    return y;
}

// ── Explicit instantiations ───────────────────────────────────────────────
// Each line forces a complete, linkable definition of Multiply for that scalar
// type. This is the only place it is compiled. Adding float support = one more
// line here; nothing changes in the header.

template Vector<double>
    MultiplyEngine::Multiply(const Matrix<double>&,
                             const Vector<double>&) const;

template Vector<std::complex<double>>
    MultiplyEngine::Multiply(const Matrix<std::complex<double>>&,
                             const Vector<std::complex<double>>&) const;

} // namespace Oblio
