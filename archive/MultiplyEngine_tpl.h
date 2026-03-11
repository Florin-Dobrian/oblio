#pragma once

// MultiplyEngine_tpl.h — header-only template style
//
// The Multiply member function template is defined entirely here. Any .cc
// that includes this header and calls Multiply<double>() or
// Multiply<complex<double>>() will compile its own copy of the function body.
// With N translation units and 2 scalar types, that is up to 2*N compilations
// of the same code.

#include "Matrix_tpl.h"
#include "Vector_tpl.h"

namespace Oblio {

class MultiplyEngine {
public:
    MultiplyEngine() = default;

    // Computes y = A * x.
    // Works for any Val type for which operator* and operator+= are defined,
    // including double and std::complex<double>.
    template<class Val>
    Vector<Val> Multiply(const Matrix<Val>& A, const Vector<Val>& x) const {
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

private:
    MultiplyEngine(const MultiplyEngine&) = delete;
    MultiplyEngine& operator=(const MultiplyEngine&) = delete;
};

} // namespace Oblio
