#include "oblio/DirectSolver.h"

#include "oblio/MultiplyEngine.h"
#include "oblio/NumFactorEngine.h"
#include "oblio/SolveEngine.h"
#include "oblio/SymFactorEngine.h"

#include <complex>

namespace Oblio {

namespace {

// Dynamic LDL delays a pivot it cannot take, which grows a front. A flat buffer cannot grow, so
// these two need per-supernode storage; the rest are cheaper in the flat factor.
bool needsDynamicStorage(Factorization factorization) {
    return factorization == Factorization::DynamicLDLT
        || factorization == Factorization::DynamicLDLH;
}

} // namespace

template<class Val>
void DirectSolver<Val>::setOrderMethod(OrderMethod method) {
    mOrderMethod = method;
    mAnalyzed    = false;   // the ordering is the first thing the analysis computes
    mFactored    = false;
}

template<class Val>
void DirectSolver<Val>::setFactorization(Factorization factorization) {
    mFactorization = factorization;
    mFactored      = false;   // the analysis survives: it depends on the pattern, not the method
}

template<class Val>
void DirectSolver<Val>::setTraversal(Traversal traversal) {
    mTraversal = traversal;
    mFactored  = false;
}

template<class Val>
bool DirectSolver<Val>::analyze(const SparseMatrix<Val>& A) {
    mAnalyzed = false;
    mFactored = false;

    const OrderEngine ord(mOrderMethod);
    if (!ord.compute(A, mPermutation))
        return false;

    const ElmForestEngine fe;
    if (!fe.compute(A, mPermutation, mForest))
        return false;

    const SymFactorEngine se;
    if (!se.compute(A, mPermutation, mForest, mSymFactor))
        return false;

    mSize     = A.size();
    mAnalyzed = true;
    return true;
}

template<class Val>
bool DirectSolver<Val>::factor(const SparseMatrix<Val>& A) {
    mFactored = false;
    if (!mAnalyzed || A.size() != mSize)
        return false;

    NumFactorEngine ne(mFactorization, mTraversal);
    ne.setPivotThreshold(mPivotThreshold);

    mUsesDynamicStorage = needsDynamicStorage(mFactorization);
    mFactored = mUsesDynamicStorage ? ne.compute(A, mPermutation, mSymFactor, mDynamicFactor)
                                    : ne.compute(A, mPermutation, mSymFactor, mStaticFactor);
    return mFactored;
}

template<class Val>
bool DirectSolver<Val>::solve(const Vector<Val>& b, Vector<Val>& x) const {
    if (!mFactored || b.size() != mSize)
        return false;

    const SolveEngine sol;
    return mUsesDynamicStorage ? sol.compute(mPermutation, mDynamicFactor, b, x)
                               : sol.compute(mPermutation, mStaticFactor, b, x);
}

template<class Val>
bool DirectSolver<Val>::compute(const SparseMatrix<Val>& A, const Vector<Val>& b, Vector<Val>& x) {
    return analyze(A) && factor(A) && solve(b, x);
}

template<class Val>
std::size_t DirectSolver<Val>::numPerturbations() const {
    return mUsesDynamicStorage ? mDynamicFactor.numPerturbations()
                               : mStaticFactor.numPerturbations();
}

template<class Val>
double DirectSolver<Val>::relativeResidual(const SparseMatrix<Val>& A, const Vector<Val>& b,
                                           const Vector<Val>& x) const {
    const MultiplyEngine mul;
    Vector<Val>          r(b.size());
    if (!mul.residual(A, x, b, r))
        return -1;

    const double normB = b.norm();
    return normB > 0 ? r.norm() / normB : r.norm();
}

template class DirectSolver<double>;
template class DirectSolver<std::complex<double>>;

} // namespace Oblio
