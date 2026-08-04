#include "oblio/DirectSolver.h"

#include "oblio/MultiplyEngine.h"
#include "oblio/NumFactorEngine.h"
#include "oblio/SolveEngine.h"
#include "oblio/SymFactorEngine.h"

#include <algorithm>
#include <type_traits>
#include <complex>

namespace Oblio {

template<class Val>
void DirectSolver<Val>::setOrdering(Ordering ordering) {
    mOrdering = ordering;
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

    const OrderEngine ordEng(mOrdering);
    if (!ordEng.compute(A, mPermutation))
        return false;

    // The multifrontal child ordering is part of the analysis, not the factorization, so the
    // traversal has to be known here. setTraversal invalidates the analysis when it changes this.
    const ElmForestEngine efEng(mSupernodes, mAmalgamation, mTraversal == Traversal::Multifrontal);
    if (!efEng.compute(A, mPermutation, mElmForest))
        return false;

    const SymFactorEngine sfEng;
    if (!sfEng.compute(A, mPermutation, mElmForest, mSymFactor))
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

    NumFactorEngine nfEng(mFactorization, mTraversal);
    nfEng.setPivotThreshold(mPivotThreshold);

    // Dynamic pivoting delays a column, which expands a front, which a flat buffer cannot do.
    mUsesDynamicStorage = dynamicPivoting(mFactorization);
    mFactored = mUsesDynamicStorage ? nfEng.compute(A, mPermutation, mSymFactor, mNumFactorDynamic)
                                    : nfEng.compute(A, mPermutation, mSymFactor, mNumFactorStatic);
    return mFactored;
}

template<class Val>
bool DirectSolver<Val>::solve(const Vector<Val>& b, Vector<Val>& x) const {
    if (!mFactored || b.size() != mSize)
        return false;

    const SolveEngine solEng;
    return mUsesDynamicStorage ? solEng.compute(mPermutation, mNumFactorDynamic, b, x)
                               : solEng.compute(mPermutation, mNumFactorStatic, b, x);
}

template<class Val>
bool DirectSolver<Val>::compute(const SparseMatrix<Val>& A, const Vector<Val>& b, Vector<Val>& x) {
    return analyze(A) && factor(A) && solve(b, x);
}

template<class Val>
std::size_t DirectSolver<Val>::numPerturbations() const {
    return mUsesDynamicStorage ? mNumFactorDynamic.numPerturbations()
                               : mNumFactorStatic.numPerturbations();
}

// The three below describe the pivot search, which exists only under dynamic LDL. A static
// factorization pivots on the diagonal in symbolic order and chooses nothing, so it delays no
// column and takes no 2x2 block, and reporting zero for it is the accurate answer rather than a
// placeholder.

template<class Val>
std::size_t DirectSolver<Val>::numDelayedColumns() const {
    if (!mUsesDynamicStorage) return 0;
    std::size_t delayed = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(mNumFactorDynamic.snodeSize()); ++kk)
        delayed += mNumFactorDynamic.delaySize(kk);
    return delayed;
}

// pivotType marks a 1x1 pivot with 1 and the two halves of a 2x2 with 2 and 3, so counting the 2s
// counts blocks rather than columns.
template<class Val>
std::size_t DirectSolver<Val>::numPivots1x1() const {
    if (!mUsesDynamicStorage) return 0;
    const std::vector<std::int32_t>& pivotType = mNumFactorDynamic.pivotType();
    return static_cast<std::size_t>(std::count(pivotType.begin(), pivotType.end(), 1));
}

template<class Val>
std::size_t DirectSolver<Val>::numPivots2x2() const {
    if (!mUsesDynamicStorage) return 0;
    const std::vector<std::int32_t>& pivotType = mNumFactorDynamic.pivotType();
    return static_cast<std::size_t>(std::count(pivotType.begin(), pivotType.end(), 2));
}

// Sylvester's law of inertia, read off D. See the header for what this is for and where it stops
// being reliable. Two walks, because the two storages hold D differently: statically pivoted, D is
// the block's diagonal and nothing else; dynamically pivoted, pivotType says which columns pair
// into a 2x2 and the block's four entries sit where the factorization left them, at the same
// offsets diagonalDynamic reads them from in SolveEngine.
//
// Cholesky needs no special case. It stores C, whose diagonal is positive wherever it succeeded,
// and it succeeds only on a positive definite matrix, so reading signs gives all positive and that
// is the right answer.
template<class Val>
bool DirectSolver<Val>::inertia(Inertia& inertia) const {
    inertia = Inertia{};
    if (!mFactored)
        return false;

    // A complex-symmetric LDL^T has complex eigenvalues, so there are no signs to count. Over the
    // reals the two transposes coincide and every factorization qualifies.
    if (!std::is_same<Val, double>::value && !hermitian(mFactorization))
        return false;

    // One sign, from a quantity that is real whatever Val is: D's diagonal is real for every
    // factorization that reaches here.
    const auto tally = [&inertia](double d) {
        if      (d > 0) ++inertia.positive;
        else if (d < 0) ++inertia.negative;
        else            ++inertia.zero;
    };

    if (!mUsesDynamicStorage) {
        for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(mNumFactorStatic.snodeSize()); ++jj) {
            const std::size_t frontSize  = mNumFactorStatic.frontSize(jj);
            const std::size_t numNodeIdx = frontSize + mNumFactorStatic.updateSize(jj);
            const Val*        val        = mNumFactorStatic.val(jj);
            for (std::size_t j = 0; j < frontSize; ++j)
                tally(std::real(val[j * numNodeIdx + j]));
        }
        return true;
    }

    const std::vector<std::int32_t>& pivotType = mNumFactorDynamic.pivotType();
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(mNumFactorDynamic.snodeSize()); ++jj) {
        const std::size_t   frontSize  = mNumFactorDynamic.frontSize(jj);
        const std::size_t   numNodeIdx = frontSize + mNumFactorDynamic.delaySize(jj)
                                                   + mNumFactorDynamic.updateSize(jj);
        const std::int32_t* nodeIdx    = mNumFactorDynamic.nodeIdx(jj);
        const Val*          val        = mNumFactorDynamic.val(jj);

        const auto at = [numNodeIdx](std::size_t r, std::size_t c) { return c * numNodeIdx + r; };

        for (std::size_t j = 0; j < frontSize; ) {
            if (pivotType[nodeIdx[j]] == 1) {
                tally(std::real(val[at(j, j)]));
                ++j;
                continue;
            }

            // A 2x2 block. Its two eigenvalues have opposite signs exactly when the determinant is
            // negative, and share the trace's sign when it is positive, which needs no square root
            // and no eigenvalue.
            const double d11 = std::real(val[at(j, j)]);
            const double d22 = std::real(val[at(j + 1, j + 1)]);
            const double d21 = std::real(val[at(j + 1, j)]);
            const double d12 = std::real(val[at(j, j + 1)]);
            const double det = d11 * d22 - d21 * d12;
            const double trace = d11 + d22;

            if (det < 0)      { ++inertia.positive; ++inertia.negative; }
            else if (det > 0) { if (trace > 0) inertia.positive += 2; else inertia.negative += 2; }
            else              { ++inertia.zero; tally(trace); }
            j += 2;
        }
    }
    return true;
}

template<class Val>
double DirectSolver<Val>::relativeResidual(const SparseMatrix<Val>& A, const Vector<Val>& b,
                                           const Vector<Val>& x) const {
    const MultiplyEngine mulEng;
    Vector<Val>          r(b.size());
    if (!mulEng.residual(A, x, b, r))
        return -1;

    const double normB = b.norm();
    return normB > 0 ? r.norm() / normB : r.norm();
}

template class DirectSolver<double>;
template class DirectSolver<std::complex<double>>;

} // namespace Oblio
