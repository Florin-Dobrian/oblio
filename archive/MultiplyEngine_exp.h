#pragma once

// MultiplyEngine_exp.h — Case 2: explicit instantiation (forcing only)
//
// MultiplyEngine itself is a non-templated class; only its Multiply member
// function is a template. The function body lives in MultiplyEngine_exp.cc, not
// here — so a translation unit that calls Multiply<Val> sees only the declaration
// and cannot implicitly instantiate it; it links the specialisation forced in the
// .cc. No `extern template` is needed, because there is no visible body to
// suppress. (Case 3 — _ext — adds extern-template function declarations here;
// with a body-less header they are belt-and-suspenders.)
//
// This is the pattern for most of Oblio's engine classes — the class need not be
// templated, only the methods that touch Val-typed data.

#include "Matrix_exp.h"
#include "Vector_exp.h"

namespace Oblio {

class MultiplyEngine {
public:
    MultiplyEngine() = default;

    template<class Val>
    Vector<Val> Multiply(const Matrix<Val>& A, const Vector<Val>& x) const;

private:
    MultiplyEngine(const MultiplyEngine&) = delete;
    MultiplyEngine& operator=(const MultiplyEngine&) = delete;
};

} // namespace Oblio
