#pragma once

// Vector.h - a dense vector of scalars, one column.
//
// The right-hand side of a solve, and its solution. Nothing more: no arithmetic, no BLAS, no
// cleverness. The engines do the work.
//
// **One column, for now.** 0.9 has two vector classes, and the split is principled: with a single
// right-hand side there is no level-3 BLAS to be had, so its solve is scalar and works directly on
// the vector with indirect indexing. With *many* right-hand sides a supernode's rows become a
// dense block, and the solve becomes TRSM and GEMM, which is worth the gather and scatter that
// packing demands. 0.9's `SingleVector` and `MultipleVector` are those two cases. 10.12 kept only
// the first, and so do we.
//
// A multi-column version is a real thing to add later, not a hypothetical: it is where the solve
// phase gets its level-3 BLAS. But it is a performance path, not a correctness one, and the single
// column is what a residual check needs.

#include "oblio/Types.h"

#include <complex>
#include <cstddef>
#include <vector>

namespace Oblio {

class MultiplyEngine;
class SolveEngine;

template<class Val>
class Vector {
public:
    Vector() = default;

    // Defined in Vector.cpp, not here: both guard the size (via checkIndexRange) and so can throw,
    // and keeping them out of the header keeps that exception path out of the translation units that
    // compile the hot multiply and solve kernels (an in-header throw was measured to perturb such a
    // loop's codegen). The default constructor and the inline accessors below stay in the header.
    explicit Vector(std::size_t size);

    // Takes the values and derives the size from them. This used to be Vector(size, val), taking
    // both and checking neither against the other, so Vector(10, std::vector<double>(3)) built an
    // object reporting size 10 over three values. Deriving makes that state unrepresentable rather
    // than merely detected, which needs no test to stay true. Nothing in the tree called the
    // two-argument form.
    explicit Vector(std::vector<Val> val);

    // **Values change in place; the SIZE changes by assignment.** Unlike SparseMatrix, this
    // class is writable, so `x[i] = v` and the non-const `val()` are how a caller fills a
    // vector, and that is what the solve and multiply paths do. Assignment is for the other
    // case:
    //
    //   x = Vector<double>(n);                  // n zeros, whatever x held before
    //   x = Vector<double>(std::move(vals));    // size derived from the values
    //
    // The copy and move operations are the compiler's, since none is declared here, so the old
    // buffer is released and a copy is independent of its source. Assignment runs the
    // constructor, so the size guard applies to the new object as it did to the first.

    std::size_t size() const { return mSize; }

    const std::vector<Val>& val() const { return mVal; }
    std::vector<Val>&       val()       { return mVal; }

    Val  operator[](std::size_t i) const { return mVal[i]; }
    Val& operator[](std::size_t i)       { return mVal[i]; }

    // The two-norm. Used for residuals, where the ratio is what matters, so the square root is
    // worth taking rather than leaving the caller to.
    double norm() const {
        double s = 0;
        for (const Val& v : mVal)
            s += std::norm(v);   // |v|^2, for both double and complex
        return std::sqrt(s);
    }

private:
    std::size_t      mSize = 0;
    std::vector<Val> mVal;

    friend class MultiplyEngine;
    friend class SolveEngine;
};

extern template class Vector<double>;
extern template class Vector<std::complex<double>>;

} // namespace Oblio
