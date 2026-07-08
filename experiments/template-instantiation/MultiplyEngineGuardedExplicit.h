#pragma once

// MultiplyEngineGuardedExplicit.h — Guarded explicit: plain explicit + extern template guard
//
// MultiplyEngine itself is a non-templated class. Only its Multiply member
// function is a template, so extern template applies at the function level
// rather than the class level. The function body lives in MultiplyEngineGuardedExplicit.cpp.
//
// Note: this is the pattern that applies to most of Oblio's engine classes —
// the class itself does not need to be templated, only the methods that touch
// Val-typed data.

#include "MatrixGuardedExplicit.h"
#include "VectorGuardedExplicit.h"

namespace Oblio {

class MultiplyEngine {
public:
    MultiplyEngine() = default;

    template<class Val>
    Vector<Val> Multiply(const Matrix<Val>& A, const Vector<Val>& x) const;

private:
    MultiplyEngine(const MultiplyEngine&) = delete;
    MultiplyEngine& operator=(const MultiplyEngine&) = delete;
};

// Suppress implicit instantiation of the two specialisations we provide.
extern template Vector<double>
    MultiplyEngine::Multiply(const Matrix<double>&,
                             const Vector<double>&) const;

extern template Vector<std::complex<double>>
    MultiplyEngine::Multiply(const Matrix<std::complex<double>>&,
                             const Vector<std::complex<double>>&) const;

} // namespace Oblio
