// MultiplyEnginePlainExplicit.cpp — Plain explicit: bodies in .cpp, header signatures only
//
// The Multiply function body is defined here and compiled once per scalar type.
// Other translation units that call Multiply<double>() or Multiply<complex<double>>()
// see only the declaration in the header, so they cannot implicitly instantiate —
// they resolve the symbol from this object file at link time.

#include "MultiplyEnginePlainExplicit.h"
#include <cassert>
#include <cstdint>   // std::int32_t

namespace Oblio {

template<class Val>
Vector<Val> MultiplyEngine::Multiply(const Matrix<Val>& A,
                                     const Vector<Val>& x) const {
    assert(A.cols() == x.size());
    Vector<Val> y(A.rows());
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(A.rows()); ++i) {
        Val sum{0};
        for (std::int32_t j = 0; j < static_cast<std::int32_t>(A.cols()); ++j) {
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
