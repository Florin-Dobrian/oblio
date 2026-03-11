#include "oblio/DenseMatrix.h"
#include <algorithm>
#include <complex>

namespace Oblio {

template<class Val>
DenseMatrix<Val>::DenseMatrix()
    : mRows(0), mCols(0) {}

template<class Val>
DenseMatrix<Val>::DenseMatrix(Size rows, Size cols)
    : mRows(rows), mCols(cols), mData(rows * cols, Val{0}) {}

template<class Val>
DenseMatrix<Val>::DenseMatrix(Size rows, Size cols, Val fill)
    : mRows(rows), mCols(cols), mData(rows * cols, fill) {}

template<class Val>
DenseMatrix<Val>::DenseMatrix(Size rows, Size cols, const Val* src)
    : mRows(rows), mCols(cols), mData(src, src + rows * cols) {}

template<class Val>
DenseMatrix<Val>::DenseMatrix(const DenseMatrix& o)
    : mRows(o.mRows), mCols(o.mCols), mData(o.mData) {}

template<class Val>
DenseMatrix<Val>::DenseMatrix(DenseMatrix&& o) noexcept
    : mRows(o.mRows), mCols(o.mCols), mData(std::move(o.mData)) {
    o.mRows = 0; o.mCols = 0;
}

template<class Val>
DenseMatrix<Val>& DenseMatrix<Val>::operator=(const DenseMatrix& o) {
    if (this != &o) {
        mRows = o.mRows; mCols = o.mCols; mData = o.mData;
    }
    return *this;
}

template<class Val>
DenseMatrix<Val>& DenseMatrix<Val>::operator=(DenseMatrix&& o) noexcept {
    if (this != &o) {
        mRows = o.mRows; mCols = o.mCols;
        mData = std::move(o.mData);
        o.mRows = 0; o.mCols = 0;
    }
    return *this;
}

template<class Val>
Size DenseMatrix<Val>::rows()  const { return mRows; }

template<class Val>
Size DenseMatrix<Val>::cols()  const { return mCols; }

template<class Val>
Size DenseMatrix<Val>::ld()    const { return mRows; }

template<class Val>
bool DenseMatrix<Val>::empty() const { return mRows == 0 || mCols == 0; }

template<class Val>
Val* DenseMatrix<Val>::data() { return mData.data(); }

template<class Val>
const Val* DenseMatrix<Val>::data() const { return mData.data(); }

template<class Val>
Val& DenseMatrix<Val>::operator()(Size i, Size j) {
    return mData[j * mRows + i];
}

template<class Val>
const Val& DenseMatrix<Val>::operator()(Size i, Size j) const {
    return mData[j * mRows + i];
}

template<class Val>
Val* DenseMatrix<Val>::col(Size j) {
    return mData.data() + j * mRows;
}

template<class Val>
const Val* DenseMatrix<Val>::col(Size j) const {
    return mData.data() + j * mRows;
}

template<class Val>
void DenseMatrix<Val>::resize(Size rows, Size cols) {
    mRows = rows; mCols = cols;
    mData.assign(rows * cols, Val{0});
}

template<class Val>
void DenseMatrix<Val>::setZero() {
    std::fill(mData.begin(), mData.end(), Val{0});
}

template<class Val>
void DenseMatrix<Val>::fill(Val v) {
    std::fill(mData.begin(), mData.end(), v);
}

// ---- explicit instantiations ------------------------------------------------
template class DenseMatrix<double>;
template class DenseMatrix<std::complex<double>>;

} // namespace Oblio
