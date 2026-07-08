#pragma once

// MultiplyEngine.h — computes y = A*x three ways, to contrast access patterns:
//
//   multiplyByApi    — through the public operator()/operator[]. Readable, but each
//                      element is a non-inlined, cross-TU call, so the compiler can't
//                      vectorize across the calls. (Bounds-check asserts add little;
//                      the cost is the calls, not the checks — see the README timing.)
//   multiplyDirectly — `friend` access to the contiguous storage: fetch the raw
//                      pointers once and walk them. No per-element calls, so the inner
//                      loop vectorizes. Hand-written fast path.
//   multiplyWithBlas — same `friend` access to get the raw block, then hand it to BLAS
//                      gemv (dgemv_/zgemv_). This is what `friend` ultimately enables,
//                      and the real solver's fast path.
//
// All three compute the same result.

#include "Matrix.h"
#include "Vector.h"

namespace Oblio {

class MultiplyEngine {
public:
    MultiplyEngine() = default;

    template<class Val>
    Vector<Val> multiplyByApi(const Matrix<Val>& A, const Vector<Val>& x) const;

    template<class Val>
    Vector<Val> multiplyDirectly(const Matrix<Val>& A, const Vector<Val>& x) const;

    template<class Val>
    Vector<Val> multiplyWithBlas(const Matrix<Val>& A, const Vector<Val>& x) const;

private:
    MultiplyEngine(const MultiplyEngine&) = delete;
    MultiplyEngine& operator=(const MultiplyEngine&) = delete;
};

extern template Vector<double>
    MultiplyEngine::multiplyByApi(const Matrix<double>&, const Vector<double>&) const;
extern template Vector<std::complex<double>>
    MultiplyEngine::multiplyByApi(const Matrix<std::complex<double>>&,
                                  const Vector<std::complex<double>>&) const;
extern template Vector<double>
    MultiplyEngine::multiplyDirectly(const Matrix<double>&, const Vector<double>&) const;
extern template Vector<std::complex<double>>
    MultiplyEngine::multiplyDirectly(const Matrix<std::complex<double>>&,
                                     const Vector<std::complex<double>>&) const;
extern template Vector<double>
    MultiplyEngine::multiplyWithBlas(const Matrix<double>&, const Vector<double>&) const;
extern template Vector<std::complex<double>>
    MultiplyEngine::multiplyWithBlas(const Matrix<std::complex<double>>&,
                                     const Vector<std::complex<double>>&) const;

} // namespace Oblio
