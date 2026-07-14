#pragma once

// SolveEngine.h - solve A x = b, given a factorization of A.
//
// Three passes, and which of them run depends on the factorization:
//
//   Cholesky   A = L L^H       forward:  L y = b        backward: L^H x = y
//   LDL        A = L D L^H     forward:  L y = b        diagonal: D z = y        backward: L^H x = z
//
// For Cholesky, L's diagonal holds the factor's own diagonal, so the forward pass divides by it.
// For LDL, L is *unit* lower triangular (the 1s implicit), the diagonal holds D, and the division
// is a separate pass. That is the whole difference, and it falls straight out of the storage.
//
// **The conjugate matters in the backward pass, and this is where 10.12 repeats its bug.** Its
// backward solve is
//
//     y[col] -= y[row] * val[...];        // no conjugate: this is L^-T, not L^-H
//
// which is right for a complex *symmetric* factor and wrong for a Hermitian one. Ours conjugates
// when the factorization does (Cholesky and LDLH), and does not when it does not (LDLT). Same
// `hermitian()` predicate the factorization uses, so it is one rule stated once, not two rules
// that must be kept in step.
//
// **Scalar, not BLAS, and deliberately.** With a single right-hand side there is no level-3 BLAS
// to be had: a supernode's contribution is a triangular solve against one vector, which is BLAS-2
// at best, and reaching it would need the supernode's rows gathered into a dense block and
// scattered back afterwards. 0.9 makes exactly this call: its `SingleVector` solve is scalar and
// works directly on the vector through indirect indexing, while its `MultipleVector` solve gathers,
// calls TRSM and GEMM, and scatters. Many right-hand sides make a supernode a *matrix* operation,
// and only then does the packing pay. That path is worth adding; it is not needed to be correct.

#include "oblio/NumFactorStatic.h"
#include "oblio/Permutation.h"
#include "oblio/Types.h"
#include "oblio/Vector.h"

#include <complex>

namespace Oblio {

class SolveEngine {
public:
    SolveEngine() = default;

    // Solve in the factor's own ordering: y := A' ^-1 y, where A' = P A P^T. In place.
    //
    // Rarely what a caller wants directly, but it is the whole of the arithmetic, and the
    // permuting overload below is a thin wrapper on it.
    template<class Val>
    bool compute(const NumFactorStatic<Val>& f, Vector<Val>& y) const;

    // Solve A x = b, in A's ordering. Permutes b into the factor's ordering, solves, and permutes
    // back. This is the one to call.
    template<class Val>
    bool compute(const Permutation& p, const NumFactorStatic<Val>& f,
                 const Vector<Val>& b, Vector<Val>& x) const;

private:
    // L y = b. Descending the supernodes, which is a topological order, so a supernode's own
    // columns are solved before anything below them is updated.
    template<class Val>
    void forward(const NumFactorStatic<Val>& f, Vector<Val>& y) const;

    // D z = y. LDL only. Trivial while pivots are 1x1; a 2x2 block solve when dynamic LDL brings
    // them.
    template<class Val>
    void diagonal(const NumFactorStatic<Val>& f, Vector<Val>& y) const;

    // L^H x = z (or L^T, for LDLT). Ascending, the mirror of the forward pass.
    template<class Val>
    void backward(const NumFactorStatic<Val>& f, Vector<Val>& y) const;

    // Supernode kk's dense block. The same seam NumFactorEngine has, and for the same reason: one
    // overload per factor storage, fetched at the moment of use and never hoisted. The solve reads
    // a *finished* factor, so nothing grows under it and the hazard does not arise here; the shape
    // is kept identical so the dynamic factor slots in without a second solve.
    template<class Val>
    static const Val* blockOf(const NumFactorStatic<Val>& f, std::int32_t kk) {
        return f.mVal.data() + f.mValPtr[kk];
    }

    // Whether the backward pass conjugates. The same rule the factorization uses, stated once.
    static bool hermitian(Factorization factorization) {
        return factorization == Factorization::Cholesky
            || factorization == Factorization::StaticLDLH
            || factorization == Factorization::DynamicLDLH;
    }

    // Whether the diagonal is a separate pass. Cholesky folds it into L; LDL does not.
    static bool separateDiagonal(Factorization factorization) {
        return factorization != Factorization::Cholesky;
    }
};

} // namespace Oblio
