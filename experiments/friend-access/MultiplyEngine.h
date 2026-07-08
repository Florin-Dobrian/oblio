#pragma once

// MultiplyEngine.h — computes y = A*x two ways, to contrast access patterns:
//
//   multiplyByApi    — through the public operator()/operator[]. Readable, but each
//                      element is a (non-inlined, cross-TU) call with a bounds-check
//                      assert, and the compiler can't vectorize across the calls.
//   multiplyDirectly — `friend` access to the contiguous storage: fetch the raw
//                      pointers once and walk them. No per-element calls, no asserts,
//                      and the inner loop vectorizes. This is the performance path.
//
// Both compute the same result. (No BLAS yet — a later step hands the raw block to
// dgemv/dgemm, which is faster still; friend access is what enables that handoff.)

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

} // namespace Oblio
