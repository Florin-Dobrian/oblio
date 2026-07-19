#ifndef OBLIO_DIRECT_SOLVER_H
#define OBLIO_DIRECT_SOLVER_H

// DirectSolver.h - the whole pipeline behind one object.
//
// Everything else in Oblio is an engine: a stateless worker that turns inputs into one output
// object. OrderEngine produces a Permutation, SymFactorEngine a SymFactor, NumFactorEngine a
// numeric factor. Wiring them by hand is what examples/pipeline.cpp does, and it is the right way
// to see the seams, but it is not how a solver should be used.
//
// This class is not an engine, which is why it is not named one. It holds state: the permutation,
// the forest, the symbolic factor and the numeric factor, all of it living between calls. That
// state is the point. The three phases have different costs and different lifetimes:
//
//   analyze(A)   ordering, elimination forest, symbolic factorization. Depends only on the
//                *pattern* of A, so a sequence of matrices sharing a pattern analyzes once.
//   factor(A)    the numeric factorization. Depends on the values, so it reruns when they change,
//                reusing the analysis.
//   solve(b, x)  one triangular solve pair. Cheap, and rerun per right-hand side.
//
// The name says direct, the class of method that factors the matrix outright, as opposed to the
// iterative solvers that approach a solution by repeated multiplication. If those ever arrive here
// they sit beside this class rather than inside it.
//
// The storage of the numeric factor is chosen for the caller. Dynamic LDL delays pivots it cannot
// take, growing a front, so it requires the per-supernode storage of NumFactorDynamic; everything
// else takes the cheaper flat NumFactorStatic. See NumFactorDynamic.h for why the two exist.

#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/NumFactorDynamic.h"
#include "oblio/NumFactorStatic.h"
#include "oblio/OrderEngine.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"
#include "oblio/SymFactor.h"
#include "oblio/Types.h"
#include "oblio/Vector.h"

#include <complex>
#include <cstddef>

namespace Oblio {

template<class Val>
class DirectSolver {
public:
    // Every setting is available two ways: in the constructor, in pipeline order, and through a
    // setter afterwards. The defaults are each engine's own default, so a DirectSolver built with
    // no arguments behaves as the engines it drives would.
    explicit DirectSolver(OrderMethod   method        = OrderMethod::MMD,
                          Factorization factorization = Factorization::Cholesky,
                          Traversal     traversal     = Traversal::LeftLooking)
        : mOrderMethod(method), mFactorization(factorization), mTraversal(traversal) {}

    // Configuration. Changing any of these invalidates what has been computed so far: the ordering
    // choice changes the analysis, the factorization and traversal change the numeric factor.
    void          setOrderMethod(OrderMethod method);
    OrderMethod   orderMethod() const   { return mOrderMethod; }
    void          setFactorization(Factorization factorization);
    Factorization factorization() const { return mFactorization; }
    void          setTraversal(Traversal traversal);
    Traversal     traversal() const     { return mTraversal; }

    // Dynamic LDL only: how large a candidate pivot must be relative to its column before it is
    // accepted in place rather than delayed. See NumFactorEngine.
    void   setPivotThreshold(double threshold) { mPivotThreshold = threshold; }
    double pivotThreshold() const              { return mPivotThreshold; }

    // The three phases, in order. factor requires a prior analyze, solve a prior factor, and each
    // returns false rather than pretending otherwise. A matrix passed to factor must have the same
    // size as the one analyzed, since the analysis describes its pattern.
    bool analyze(const SparseMatrix<Val>& A);
    bool factor(const SparseMatrix<Val>& A);
    bool solve(const Vector<Val>& b, Vector<Val>& x) const;

    // All three at once, for the caller with one matrix and one right-hand side.
    bool compute(const SparseMatrix<Val>& A, const Vector<Val>& b, Vector<Val>& x);

    bool analyzed() const { return mAnalyzed; }
    bool factored() const { return mFactored; }

    // The intermediate results, for a caller who wants to inspect the ordering or the fill rather
    // than only the solution.
    const Permutation& permutation() const { return mPermutation; }
    const ElmForest&   forest() const      { return mForest; }
    const SymFactor&   symFactor() const   { return mSymFactor; }

    // Whether the numeric factor went into per-supernode storage, and how many pivots the
    // factorization had to replace (static LDL perturbs rather than failing).
    bool        usesDynamicStorage() const { return mUsesDynamicStorage; }
    std::size_t numPerturbations() const;

    // ||A x - b|| / ||b||, the one number that says whether the pipeline worked. Recomputes A x, so
    // it costs a multiplication; it is a convenience, not part of solving.
    double relativeResidual(const SparseMatrix<Val>& A, const Vector<Val>& b,
                            const Vector<Val>& x) const;

private:
    // Set by the constructor, which is where their defaults live, so they are not repeated here.
    OrderMethod   mOrderMethod;
    Factorization mFactorization;
    Traversal     mTraversal;

    double mPivotThreshold = 0.1;   // not a constructor argument: dynamic LDL tuning, as on the engine

    Permutation mPermutation;
    ElmForest   mForest;
    SymFactor   mSymFactor;

    // Both storages are declared, one is filled. An unused vector costs nothing, and this keeps the
    // factor a concrete member rather than something reached through a pointer or a variant.
    NumFactorStatic<Val>  mStaticFactor;
    NumFactorDynamic<Val> mDynamicFactor;

    std::size_t mSize               = 0;
    bool        mUsesDynamicStorage = false;
    bool        mAnalyzed           = false;
    bool        mFactored           = false;
};

extern template class DirectSolver<double>;
extern template class DirectSolver<std::complex<double>>;

} // namespace Oblio

#endif // OBLIO_DIRECT_SOLVER_H
