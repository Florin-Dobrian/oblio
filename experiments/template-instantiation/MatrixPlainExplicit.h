#pragma once

// MatrixPlainExplicit.h - Plain explicit: bodies in .cpp, header signatures only
//
// The header contains only the class declaration; the implementation lives in
// MatrixPlainExplicit.cpp. Because no member bodies are visible here, a translation unit
// that includes this header cannot implicitly instantiate Matrix<Val>, it emits
// undefined references and resolves them at link time against the explicit
// instantiations forced in MatrixPlainExplicit.cpp. No `extern template` is needed: the
// build win comes from the definitions being absent from the header, not from
// suppressing instantiation.
//
// (The _GuardedExplicit variant keeps this exact layout and adds `extern template`
// lines. With a declaration-only header those lines suppress nothing, they are
// documentation, taking effect only if definitions were later kept in the header.)
//
// Adding a new scalar type (e.g. float) = one new explicit instantiation line in
// MatrixPlainExplicit.cpp. Nothing changes here.

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

} // namespace Oblio
