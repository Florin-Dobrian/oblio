#pragma once
#include "oblio/Types.h"
#include "oblio/Factors.h"
#include "oblio/Vector.h"
#include "oblio/DenseMatrix.h"
#include <complex>

namespace Oblio {

// SolveEngine<Val> — triangular solves for LL^T (Cholesky) and LDL^T variants.
//
// Single RHS:   solve(lu, y)   — y is overwritten in place (b → x)
// Multiple RHS: solve(lu, Y)   — each column of Y is solved; batched BLAS per
//                                supernode for efficiency.
//
// Three passes: Forward (L^-1 b), Diagonal (D^-1 y), Backward (L^-T y).
// The scope can be restricted to any subset of the three passes.

template<class Val>
class SolveEngine {
public:
    SolveEngine();
    void setScope(SolveScope s);

    // Single RHS — y is overwritten (b on entry, x on exit).
    Err solve(const Factors<Val>& lu, Vector<Val>& y) const;

    // Multiple RHS — each column of Y is solved independently using
    // per-supernode batched BLAS (dtrsm / dgemm on all nRHS columns at once).
    Err solve(const Factors<Val>& lu, DenseMatrix<Val>& Y) const;

private:
    SolveScope mScope;

    // ---- single-RHS passes ----
    void fwdCC   (const Factors<Val>& lu, Vector<Val>& y) const;
    void fwdSLDL (const Factors<Val>& lu, Vector<Val>& y) const;
    void fwdDLDL (const Factors<Val>& lu, Vector<Val>& y) const;
    void diagSLDL(const Factors<Val>& lu, Vector<Val>& y) const;
    void diagDLDL(const Factors<Val>& lu, Vector<Val>& y) const;
    void bwdCC   (const Factors<Val>& lu, Vector<Val>& y) const;
    void bwdSLDL (const Factors<Val>& lu, Vector<Val>& y) const;
    void bwdDLDL (const Factors<Val>& lu, Vector<Val>& y) const;

    // ---- multi-RHS passes ----
    // Each pass operates on all nRHS columns simultaneously.
    // A local gather buffer (nI × nRHS) is used per supernode so that
    // BLAS sees a contiguous column-major matrix.
    void fwdCC_m   (const Factors<Val>& lu, DenseMatrix<Val>& Y) const;
    void fwdSLDL_m (const Factors<Val>& lu, DenseMatrix<Val>& Y) const;
    void fwdDLDL_m (const Factors<Val>& lu, DenseMatrix<Val>& Y) const;
    void diagSLDL_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const;
    void diagDLDL_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const;
    void bwdCC_m   (const Factors<Val>& lu, DenseMatrix<Val>& Y) const;
    void bwdSLDL_m (const Factors<Val>& lu, DenseMatrix<Val>& Y) const;
    void bwdDLDL_m (const Factors<Val>& lu, DenseMatrix<Val>& Y) const;
};

// Explicit instantiation declarations — definitions in SolveEngine.cc
extern template class SolveEngine<double>;
extern template class SolveEngine<std::complex<double>>;

} // namespace Oblio
