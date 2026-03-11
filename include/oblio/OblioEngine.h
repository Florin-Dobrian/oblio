#pragma once
#include "oblio/Types.h"
#include "oblio/Matrix.h"
#include "oblio/Permutation.h"
#include "oblio/Symbolic.h"
#include "oblio/Factors.h"
#include "oblio/OrderEngine.h"
#include "oblio/SymbolicEngine.h"
#include "oblio/FactorEngine.h"
#include "oblio/SolveEngine.h"
#include "oblio/Vector.h"
#include "oblio/DenseMatrix.h"
#include <complex>

namespace Oblio {

// OblioEngine<Val> — top-level driver.
//   eng.analyze(A)         — order + symbolic analysis
//   eng.factor(A)          — numerical factorization
//   eng.solve(b, x)        — triangular solve
//   eng.analyzeAndFactor(A)— convenience wrapper
template<class Val>
class OblioEngine {
public:
    OblioEngine();

    // ---- configuration ----
    void setOrderAlg     (OrderAlg   a);
    void setFactorAlg    (FactorAlg  a);
    void setFactorType   (FactorType t);
    void setSolveScope   (SolveScope s);
    void setPerturbation (RVal p);
    void setPivotThreshold(RVal a);

    // ---- phases ----
    Err analyze         (const Matrix<Val>& A);
    Err factor          (const Matrix<Val>& A);
    // Single RHS: x = A^-1 b.  x may alias b (safe to pass same object).
    Err solve           (const Vector<Val>& b,      Vector<Val>& x)      const;
    // Multiple RHS: each column of X = A^-1 (corresponding column of B).
    // X is resized to match B if needed.
    Err solve           (const DenseMatrix<Val>& B, DenseMatrix<Val>& X) const;
    Err analyzeAndFactor(const Matrix<Val>& A);

    // ---- accessors ----
    const Permutation&  perm()     const;
    const Symbolic&     symbolic() const;
    const Factors<Val>& factors()  const;
    Size nPert()  const;
    Size nSwap()  const;
    Size nDelay() const;
    Size rank()   const;

private:
    OrderEngine       mOrd;
    FactorEngine<Val> mFac;
    SolveEngine<Val>  mSol;
    Permutation       mPerm;
    Symbolic          mSym;
    Factors<Val>      mLU;
    bool              mReady;
};

// Explicit instantiation declarations — definitions in OblioEngine.cc
extern template class OblioEngine<double>;
extern template class OblioEngine<std::complex<double>>;

} // namespace Oblio
