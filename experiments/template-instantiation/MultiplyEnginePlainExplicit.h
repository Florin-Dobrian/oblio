#pragma once

// MultiplyEnginePlainExplicit.h — Plain explicit: bodies in .cpp, header signatures only
//
// MultiplyEngine itself is a non-templated class; only its Multiply member
// function is a template. The function body lives in MultiplyEnginePlainExplicit.cpp, not
// here — so a translation unit that calls Multiply<Val> sees only the declaration
// and cannot implicitly instantiate it; it links the specialisation forced in the
// .cpp. No `extern template` is needed, because there is no visible body to
// suppress. (The _GuardedExplicit variant adds extern-template function declarations here;
// with a body-less header they suppress nothing — documentation only.)
//
// This is the pattern for most of Oblio's engine classes — the class need not be
// templated, only the methods that touch Val-typed data.

#include "MatrixPlainExplicit.h"
#include "VectorPlainExplicit.h"

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
