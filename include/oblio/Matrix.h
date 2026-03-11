#pragma once
#include "oblio/Types.h"
#include <vector>
#include <complex>
#include <stdexcept>

namespace Oblio {

// Sparse symmetric matrix stored as lower triangle in CSC format.
// mColPtr[j]..mColPtr[j+1]-1 are the row indices in column j.
// Row indices are sorted ascending, diagonal included.
template<class Val>
struct Matrix {
    Size              mSize   = 0;
    bool              mIsSymm = false;
    std::vector<Size> mColPtr;   // size+1 entries
    std::vector<Size> mRowIdx;   // nnz entries (lower tri, diag included)
    std::vector<Val>  mVal;      // nnz entries

    Matrix() = default;

    // Build from COO lower triangle (including diagonal). 0-based indices.
    static Matrix fromCOO(Size n,
                          const std::vector<Size>& rows,
                          const std::vector<Size>& cols,
                          const std::vector<Val>&  vals);

    Size getSize()    const { return mSize; }
    bool isValid()    const { return mIsSymm && mSize > 0; }
    Size nnz()        const { return mColPtr.empty() ? 0 : mColPtr[mSize]; }
    Size numOffDiag() const;
};

// Explicit instantiation declarations — definitions in Matrix.cc
extern template struct Matrix<double>;
extern template struct Matrix<std::complex<double>>;

} // namespace Oblio
