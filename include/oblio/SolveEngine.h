#pragma once

// SolveEngine.h - solve A x = b, given a factorization of A.
//
// Three passes, and which of them run depends on the factorization:
//
//   Cholesky   A = C C^H       forward:  C y = b        backward: C^H x = y
//   LDL        A = L D L^H     forward:  L y = b        diagonal: D z = y        backward: L^H x = z
//
// For Cholesky, C's diagonal holds the factor's own diagonal, so the forward pass divides by it.
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

#include "oblio/NumFactorDynamic.h"
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
    template<class Val, class Factor>
    bool compute(const Factor& nf, Vector<Val>& y) const;

    // Solve A x = b, in A's ordering. Permutes b into the factor's ordering, solves, and permutes
    // back. This is the one to call.
    template<class Val, class Factor>
    bool compute(const Permutation& p, const Factor& nf,
                 const Vector<Val>& b, Vector<Val>& x) const;

private:
    // The three passes come in pairs, split on the *pivoting* axis, exactly as NumFactorEngine's
    // kernels are. The static three serve Cholesky and static LDL and are templated on the factor,
    // since either storage holds them; the dynamic three name NumFactorDynamic outright, because
    // dynamic pivoting requires that storage. Same rule, same reason, stated once by
    // dynamicPivoting() in Types.h.
    //
    // Two things separate a dynamic pass from its static twin, and both come from delayed columns:
    //
    //   The stride.  A block's leading dimension is frontSize + numberOfDelayedColumns +
    //                updateSize. The delayed columns kept their rows when shrinkEntry dropped their
    //                columns, and those rows are genuine rows of L that the solve must walk.
    //   The 2x2s.    Where pivotType says a column opens a 2x2, the entry just below its diagonal
    //                holds D's off-diagonal rather than an entry of L, so the triangular passes
    //                step over it and the diagonal pass solves the pair together.

    // L y = b. Descending the supernodes, which is a topological order, so a supernode's own
    // columns are solved before anything below them is updated.
    template<class Val, class Factor>
    void forwardStatic(const Factor& nf, Vector<Val>& y) const;
    template<class Val>
    void forwardDynamic(const NumFactorDynamic<Val>& nf, Vector<Val>& y) const;

    // D z = y. LDL only. One division per column while every pivot is 1x1; a 2x2 block solve where
    // dynamic pivoting chose a pair.
    template<class Val, class Factor>
    void diagonalStatic(const Factor& nf, Vector<Val>& y) const;
    template<class Val>
    void diagonalDynamic(const NumFactorDynamic<Val>& nf, Vector<Val>& y) const;

    // L^H x = z (or L^T, for LDLT). Ascending, the mirror of the forward pass.
    template<class Val, class Factor>
    void backwardStatic(const Factor& nf, Vector<Val>& y) const;
    template<class Val>
    void backwardDynamic(const NumFactorDynamic<Val>& nf, Vector<Val>& y) const;
};

} // namespace Oblio
