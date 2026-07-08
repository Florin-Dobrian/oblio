#pragma once

// Matrix.h — dense row-major matrix (guarded-explicit instantiation).
//
// Public API: operator()(i,j) for convenient, bounds-checked element access.
// MultiplyEngine is a `friend`, so the performance path can reach the contiguous
// storage (mVals) directly — no per-element calls, no bounds checks — while the
// public API stays available for readable, non-hot-path use. Both coexist by design.

#include <vector>
#include <complex>
#include <cstddef>
#include <cassert>

namespace Oblio {

class MultiplyEngine;  // befriended below

template<class Val>
class Matrix {
public:
    Matrix();
    Matrix(std::size_t rows, std::size_t cols, const std::vector<Val>& vals);

    Val  operator()(std::size_t i, std::size_t j) const;
    Val& operator()(std::size_t i, std::size_t j);

    std::size_t rows() const;
    std::size_t cols() const;

private:
    std::size_t      mRows;
    std::size_t      mCols;
    std::vector<Val> mVals;

    // Grants MultiplyEngine direct access to mRows/mCols/mVals for the fast path.
    friend class MultiplyEngine;
};

extern template class Matrix<double>;
extern template class Matrix<std::complex<double>>;

} // namespace Oblio
