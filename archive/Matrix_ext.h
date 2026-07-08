#pragma once

// Matrix_exp.h — explicit instantiation style
//
// The header contains only the class declaration. The implementation lives in
// Matrix_exp.cc, which is compiled exactly once. extern template tells every
// other translation unit not to instantiate these specialisations themselves;
// they will be resolved at link time from Matrix_exp.o.
//
// Adding a new scalar type (e.g. float) requires:
//   1. One new extern template line here.
//   2. One new explicit instantiation line in Matrix_exp.cc.
//   Nothing else changes.

#include <vector>
#include <complex>
#include <cstddef>

namespace Oblio {

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
};

// Suppress implicit instantiation in all other translation units.
// The definitions are provided by Matrix_exp.cc.
extern template class Matrix<double>;
extern template class Matrix<std::complex<double>>;

} // namespace Oblio
