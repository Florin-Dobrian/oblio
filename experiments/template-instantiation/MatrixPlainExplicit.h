#pragma once

// MatrixPlainExplicit.h - Plain explicit: bodies in .cpp, header signatures only
//
// The header contains the class declaration and one body; the rest of the implementation
// lives in MatrixPlainExplicit.cpp. A translation unit that includes this header cannot
// implicitly instantiate the members it cannot see, so it emits undefined references and
// resolves them at link time against the explicit instantiations forced in
// MatrixPlainExplicit.cpp. No `extern template` is needed for those: the build win comes
// from the definitions being absent from the header, not from suppressing instantiation.
//
// **The one body is the defaulted default constructor**, and it is worth knowing that this
// does NOT make this variant differ from its guarded twin. A defaulted member is defined at
// its declaration, so the obvious reading is that it can be instantiated per translation
// unit here, where nothing suppresses it, and not in _GuardedExplicit, whose `extern
// template` lines would act on it. Measured, that is false: a defaulted default constructor
// over scalars and std::vector members is trivial, so no out-of-line function is emitted
// for it at all and there is no symbol either to suppress or to link. Linking a program
// that default-constructs a Matrix without this .cpp leaves the same undefined references
// under both variants, and Matrix<double>::Matrix() is in neither list. See the README.
//
// It is written this way deliberately, because the real tree does exactly this (see
// include/oblio/Vector.h and CLAUDE.md's inline-under-the-guard exception), and a study of
// how Oblio instantiates its templates should model the arrangement Oblio uses.
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
    // The one exception to this header's declaration-only rule, and a deliberate one: a
    // defaulted default constructor is a body, so it can be instantiated here rather than
    // linked from the .cpp. See the note above on what that costs this variant.
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, const std::vector<Val>& vals);

    Val  operator()(std::size_t i, std::size_t j) const;
    Val& operator()(std::size_t i, std::size_t j);

    std::size_t rows() const;
    std::size_t cols() const;

private:
    std::size_t      mRows = 0;
    std::size_t      mCols = 0;
    std::vector<Val> mVals;
};

} // namespace Oblio
