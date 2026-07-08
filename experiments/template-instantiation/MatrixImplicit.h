#pragma once

// MatrixImplicit.h — Implicit: template inclusion (body in header)
//
// The full implementation lives here in the header. Every translation unit
// that includes this file will implicitly instantiate the template for every
// Val type it uses. For a two-type library (double, complex<double>) compiled
// across many .cpp files, this means the same code is compiled repeatedly.

#include <vector>
#include <complex>
#include <cassert>
#include <cstddef>

namespace Oblio {

template<class Val>
class Matrix {
public:
    // Constructs an empty 0x0 matrix.
    Matrix()
        : mRows(0), mCols(0) {}

    // Constructs a rows x cols matrix from a flat row-major value array.
    // The values vector must have exactly rows * cols entries.
    Matrix(std::size_t rows, std::size_t cols, const std::vector<Val>& vals)
        : mRows(rows), mCols(cols), mVals(vals) {
        assert(vals.size() == rows * cols);
    }

    // Element access (row-major storage).
    Val operator()(std::size_t i, std::size_t j) const {
        assert(i < mRows && j < mCols);
        return mVals[i * mCols + j];
    }

    Val& operator()(std::size_t i, std::size_t j) {
        assert(i < mRows && j < mCols);
        return mVals[i * mCols + j];
    }

    std::size_t rows() const { return mRows; }
    std::size_t cols() const { return mCols; }

private:
    std::size_t      mRows;
    std::size_t      mCols;
    std::vector<Val> mVals;
};

} // namespace Oblio
