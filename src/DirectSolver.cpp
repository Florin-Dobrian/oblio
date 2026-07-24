#include "oblio/DirectSolver.h"

#include "oblio/MultiplyEngine.h"
#include "oblio/NumFactorEngine.h"
#include "oblio/SolveEngine.h"
#include "oblio/SymFactorEngine.h"

#include <complex>

namespace Oblio {

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
void DirectSolver<Val>::setSupernodes(Supernodes supernodes) {
    mSupernodes = supernodes;
    mAnalyzed   = false;   // the forest is built from this, and the analysis is the forest
    mFactored   = false;
}

template<class Val>
void DirectSolver<Val>::setAmalgamation(std::optional<std::size_t> amalgamation) {
    mAmalgamation = amalgamation;
    mAnalyzed     = false;   // as above: merging supernodes changes the forest
    mFactored     = false;
}

template<class Val>
void DirectSolver<Val>::setTraversal(Traversal traversal) {
    // The analysis survives a change between left- and right-looking, which read the same forest,
    // but not one into or out of multifrontal: that reorders the children and relabels the
    // supernodes, so the forest itself differs. See ElmForestEngine::sortForOptimalMultifrontal.
    const bool wasMultifrontal = mTraversal  == Traversal::Multifrontal;
    const bool isMultifrontal  = traversal   == Traversal::Multifrontal;

    mTraversal = traversal;
    mFactored  = false;
    if (wasMultifrontal != isMultifrontal)
        mAnalyzed = false;
}

template<class Val>
bool DirectSolver<Val>::analyze(const SparseMatrix<Val>& A) {
    mAnalyzed = false;
    mFactored = false;

    const OrderEngine ord(mOrderMethod);
    if (!ord.compute(A, mPermutation))
        return false;

    // The multifrontal child ordering is part of the analysis, not the factorization, so the
    // traversal has to be known here. setTraversal invalidates the analysis when it changes this.
    const ElmForestEngine fe(mSupernodes, mAmalgamation, mTraversal == Traversal::Multifrontal);
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

    // Dynamic pivoting delays a column, which expands a front, which a flat buffer cannot do.
    mUsesDynamicStorage = dynamicPivoting(mFactorization);
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
