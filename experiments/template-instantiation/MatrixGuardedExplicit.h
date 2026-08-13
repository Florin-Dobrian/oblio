#pragma once

// MatrixGuardedExplicit.h - Guarded explicit: plain explicit + extern template guard
//
// The header contains the class declaration and one body, the defaulted default
// constructor. The rest of the implementation lives in MatrixGuardedExplicit.cpp, which is
// compiled exactly once. extern template tells every other translation unit not to
// instantiate these specialisations themselves; they will be resolved at link time from
// MatrixGuardedExplicit.o.
//
// **This is the variant the real tree uses**, and the defaulted constructor is why: a body
// in the header is permitted precisely because the class stays explicitly instantiated with
// the guard present. That is CLAUDE.md's exception to definitions-in-cpp, and
// include/oblio/Vector.h is the same arrangement.
//
// **What the guard does NOT do here, measured rather than assumed.** The tempting reading
// is that the defaulted constructor finally gives extern template something to suppress,
// where a member with no visible body offers it nothing. It does not, and not for the reason
// one might guess: the constructor is NOT trivial, Matrix holding a std::vector whose own
// default constructor is non-trivial. It is that the definition is visible in this header, so
// a translation unit needing the constructor instantiates its own weak copy and the linker
// folds them, and extern template does not stop that. Compiled at -O0, where nothing is
// inlined away, this variant and the unguarded one produce byte-identical symbol tables, both
// carrying Matrix<double>::Matrix() as a weak symbol. Linking a program that
// default-constructs a Matrix without this .cpp leaves the same undefined references either
// way, and Matrix<double>::Matrix() appears in neither list. So the guard remains
// documentation here, exactly as it was before the body arrived. The rule is what matters,
// not the mechanism: inline under the guard is a choice, inline without it is the bug. See the
// README.
//
// Adding a new scalar type (e.g. float) requires:
//   1. One new extern template line here.
//   2. One new explicit instantiation line in MatrixGuardedExplicit.cpp.
//   Nothing else changes.

#include <vector>
#include <complex>
#include <cstddef>

namespace Oblio {

template<class Val>
class Matrix {
public:
    Matrix() = default;   // the one body in this header; see the note above on the guard
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

// Suppress implicit instantiation in all other translation units.
// The definitions are provided by MatrixGuardedExplicit.cpp.
extern template class Matrix<double>;
extern template class Matrix<std::complex<double>>;

} // namespace Oblio
