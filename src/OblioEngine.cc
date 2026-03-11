#include "oblio/OblioEngine.h"
#include <complex>

namespace Oblio {

template<class Val>
OblioEngine<Val>::OblioEngine()
    : mReady(false)
{
    mOrd.setAlg(OrderAlg::eMMD);
    mFac.setAlg(FactorAlg::eMultifrontal);
    mFac.setType(FactorType::eStaticLDL);
    mSol.setScope(SolveScope::eFull);
}

template<class Val>
void OblioEngine<Val>::setOrderAlg(OrderAlg a)      { mOrd.setAlg(a);   }

template<class Val>
void OblioEngine<Val>::setFactorAlg(FactorAlg a)    { mFac.setAlg(a);   }

template<class Val>
void OblioEngine<Val>::setFactorType(FactorType t)  { mFac.setType(t);  }

template<class Val>
void OblioEngine<Val>::setSolveScope(SolveScope s)  { mSol.setScope(s); }

template<class Val>
void OblioEngine<Val>::setPerturbation(RVal p)      { mFac.setPert(p);  }

template<class Val>
void OblioEngine<Val>::setPivotThreshold(RVal a)    { mFac.setAlpha(a); }

template<class Val>
Err OblioEngine<Val>::analyze(const Matrix<Val>& A) {
    mReady = false;
    Err e = mOrd.order(A, &mPerm);
    if (e != Err::eNone) return e;
    SymbolicEngine se;
    e = se.compute(A, mPerm, &mSym);
    if (e != Err::eNone) return e;
    mReady = true;
    return Err::eNone;
}

template<class Val>
Err OblioEngine<Val>::factor(const Matrix<Val>& A) {
    if (!mReady) return Err::eInvArg;
    return mFac.factor(A, mPerm, mSym, mLU);
}

template<class Val>
Err OblioEngine<Val>::solve(const Vector<Val>& b, Vector<Val>& x) const {
    Size n = mSym.mSize;
    if (b.size() != n) return Err::eDimMismatch;
    // Permute b into the factored ordering, solve, permute solution back.
    Vector<Val> yPerm(n);
    for (Size i = 0; i < n; ++i) yPerm[mPerm.oldToNew(i)] = b[i];
    Err e = mSol.solve(mLU, yPerm);
    if (e != Err::eNone) return e;
    x.resize(n);
    for (Size i = 0; i < n; ++i) x[mPerm.newToOld(i)] = yPerm[i];
    return Err::eNone;
}

template<class Val>
Err OblioEngine<Val>::solve(const DenseMatrix<Val>& B, DenseMatrix<Val>& X) const {
    Size n = mSym.mSize;
    if (B.rows() != n) return Err::eDimMismatch;
    Size nRHS = B.cols();
    // Permute all columns of B into the factored ordering.
    DenseMatrix<Val> YPerm(n, nRHS);
    for (Size r = 0; r < nRHS; ++r)
        for (Size i = 0; i < n; ++i)
            YPerm(mPerm.oldToNew(i), r) = B(i, r);
    Err e = mSol.solve(mLU, YPerm);
    if (e != Err::eNone) return e;
    // Permute solution back.
    X.resize(n, nRHS);
    for (Size r = 0; r < nRHS; ++r)
        for (Size i = 0; i < n; ++i)
            X(mPerm.newToOld(i), r) = YPerm(i, r);
    return Err::eNone;
}

template<class Val>
Err OblioEngine<Val>::analyzeAndFactor(const Matrix<Val>& A) {
    Err e = analyze(A);
    if (e != Err::eNone) return e;
    return factor(A);
}

template<class Val>
const Permutation& OblioEngine<Val>::perm()     const { return mPerm; }

template<class Val>
const Symbolic& OblioEngine<Val>::symbolic()    const { return mSym;  }

template<class Val>
const Factors<Val>& OblioEngine<Val>::factors() const { return mLU;   }

template<class Val>
Size OblioEngine<Val>::nPert()  const { return mFac.nPert();  }

template<class Val>
Size OblioEngine<Val>::nSwap()  const { return mFac.nSwap();  }

template<class Val>
Size OblioEngine<Val>::nDelay() const { return mFac.nDelay(); }

template<class Val>
Size OblioEngine<Val>::rank()   const { return mFac.rank();   }

// ── Explicit instantiations ──────────────────────────────────────────────────
template class OblioEngine<double>;
template class OblioEngine<std::complex<double>>;

} // namespace Oblio
