#pragma once
#include "oblio/Types.h"
#include <vector>
#include <complex>

namespace Oblio {

// DenseMatrix<Val> — column-major dense matrix (rows × cols), owns its storage.
// Element (i,j) is at data()[j*rows + i], matching LAPACK/BLAS convention.
// The leading dimension (ld) always equals rows.

template<class Val>
class DenseMatrix {
public:
    // ---- construction ----
    DenseMatrix();
    DenseMatrix(Size rows, Size cols);              // zero-initialised
    DenseMatrix(Size rows, Size cols, Val fill);    // filled with fill
    DenseMatrix(Size rows, Size cols, const Val* colMajorData); // copy from raw
    DenseMatrix(const DenseMatrix&);
    DenseMatrix(DenseMatrix&&) noexcept;
    DenseMatrix& operator=(const DenseMatrix&);
    DenseMatrix& operator=(DenseMatrix&&) noexcept;
    ~DenseMatrix() = default;

    // ---- dimensions ----
    Size rows()  const;
    Size cols()  const;
    Size ld()    const;   // leading dimension = rows
    bool empty() const;

    // ---- data access ----
    Val*        data();
    const Val*  data()                           const;

    // (i,j): element at row i, column j
    Val&        operator()(Size i, Size j);
    const Val&  operator()(Size i, Size j)       const;

    // raw pointer to column j (useful for BLAS calls)
    Val*        col(Size j);
    const Val*  col(Size j)                      const;

    // ---- mutators ----
    void        resize(Size rows, Size cols);     // data lost, re-zeroed
    void        setZero();
    void        fill(Val v);

private:
    Size             mRows;
    Size             mCols;
    std::vector<Val> mData;  // column-major, length mRows*mCols
};

// Explicit instantiation declarations — definitions in DenseMatrix.cc
extern template class DenseMatrix<double>;
extern template class DenseMatrix<std::complex<double>>;

} // namespace Oblio
