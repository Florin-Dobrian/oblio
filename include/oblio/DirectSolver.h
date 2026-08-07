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
// take, expanding a front, so it requires the per-supernode storage of NumFactorDynamic; everything
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
#include <optional>

namespace Oblio {

// The inertia of a symmetric matrix: how many of its eigenvalues are positive, negative and zero.
// The three sum to the matrix order. It is a scalar-type-free fact about A, not about the factor,
// which is why it is a plain struct here rather than a member of either.
struct Inertia {
    std::size_t positive = 0;
    std::size_t negative = 0;
    std::size_t zero     = 0;
};

template<class Val>
class DirectSolver {
public:
    // Every setting is available two ways: in the constructor, in pipeline order, and through a
    // setter afterwards. The defaults are each engine's own default, so a DirectSolver built with
    // no arguments behaves as the engines it drives would.
    explicit DirectSolver(Ordering      ordering      = Ordering::MMD3,
                          Factorization factorization = Factorization::Cholesky,
                          Traversal     traversal     = Traversal::LeftLooking,
                          Supernodes    supernodes    = Supernodes::Fundamental,
                          std::optional<std::size_t> amalgamation = std::nullopt)
        : mOrdering(ordering), mFactorization(factorization), mTraversal(traversal),
          mSupernodes(supernodes), mAmalgamation(amalgamation) {}

    // Configuration. Every setting is available in the constructor, in pipeline order, and through
    // a setter here. Changing one invalidates what it must and no more: the ordering and the two
    // forest settings change the analysis, the factorization changes only the numeric factor, and
    // the traversal changes the analysis only across the multifrontal boundary, since that is where
    // the forest differs.
    void          setOrdering(Ordering ordering);
    Ordering      ordering() const      { return mOrdering; }
    void          setFactorization(Factorization factorization);
    Factorization factorization() const { return mFactorization; }
    void          setTraversal(Traversal traversal);
    Traversal     traversal() const     { return mTraversal; }

    // Dynamic LDL only: how large a candidate pivot must be relative to its column before it is
    // accepted in place rather than delayed. See NumFactorEngine.
    void   setPivotThreshold(double threshold) { mPivotThreshold = threshold; }
    double pivotThreshold() const              { return mPivotThreshold; }

    // The two forest settings. Both change the forest itself, so both invalidate a completed
    // analysis, as setOrdering does. Fundamental supernodes and no amalgamation are the
    // defaults, matching ElmForestEngine; Nodal exists to make an unamalgamated,
    // one-column-per-supernode factorization reachable for comparison.
    void       setSupernodes(Supernodes supernodes);
    Supernodes supernodes() const { return mSupernodes; }

    // Absent means no amalgamation, which is the default. A value is the number of explicitly
    // stored zeros a merge may buy per column. See ElmForestEngine.
    void                       setAmalgamation(std::optional<std::size_t> amalgamation);
    std::optional<std::size_t> amalgamation() const { return mAmalgamation; }

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
    const ElmForest&   elmForest() const   { return mElmForest; }
    const SymFactor&   symFactor() const   { return mSymFactor; }

    // Whether the numeric factor went into per-supernode storage, and how many pivots the
    // factorization had to replace (static LDL perturbs rather than failing).
    bool        usesDynamicStorage() const { return mUsesDynamicStorage; }
    std::size_t numPerturbations() const;

    // What the pivot search did, which only dynamic LDL has: a statically pivoted factorization
    // takes the diagonal in the order the symbolic factorization fixed and makes no choices, so all
    // three are zero for it rather than describing it badly.
    //
    //   numDelayedColumns  columns a supernode could not pivot on and passed to its parent, summed
    //                      over supernodes. Each one widens the parent's front, so this is what
    //                      dynamic pivoting costs in fill and work
    //   numPivots1x1       single-column pivots accepted
    //   numPivots2x2       two-column blocks accepted, taken when no single column was acceptable
    //
    // The columns are partitioned, so numPivots1x1 + 2 * numPivots2x2 == size() once factored.
    // All three are counted on demand by a scan of the factor, not maintained during it.
    std::size_t numDelayedColumns() const;
    std::size_t numPivots1x1() const;
    std::size_t numPivots2x2() const;

    // The numeric factor's three sizes, forwarded to whichever storage is live, in the same order
    // and with the same meanings as on the four classes that define them: numNodeIdx counts index
    // entries, numVal the values allocated, nnz the entries of L. Zero before factor(), which is
    // the true count for a factor that does not exist yet rather than a sentinel.
    //
    // **These are what the factorization did, not what the analysis predicted.** For a statically
    // pivoted factor the two agree, nothing having moved. For a dynamically pivoted one they do
    // not: a delayed column widens its parent's front, so these exceed symFactor()'s by exactly
    // what the delays cost. Comparing the two is how that cost is read.
    std::size_t numNodeIdx() const;
    std::size_t numVal() const;
    std::size_t nnz() const;

    // How many eigenvalues of A are positive, negative and zero. Not computed but *read*, from the
    // signs of D, which is Sylvester's law of inertia: A = L D L^H is a congruence, a congruence
    // preserves those signs, so counting them in D counts them in A without forming an eigenvalue.
    // A 1x1 pivot contributes the sign of its diagonal; a 2x2 block contributes one of each sign
    // when its determinant is negative, and two of its trace's sign when positive.
    //
    // This is what an indefinite solver is expected to report, and it is why one is used to count
    // eigenvalues in an interval (factor A - sigma I and read the negatives) or to check
    // second-order optimality on a saddle-point system.
    //
    // Returns false rather than guessing when it cannot answer: before a factorization, and for a
    // complex-symmetric LDL^T, whose eigenvalues are complex and so have no signs to count. Every
    // other case is covered, Cholesky's being trivially all positive since it factors nothing else.
    //
    // Two caveats on the answer itself. It is exact for a nonsingular A, and for a singular one a
    // zero eigenvalue lands on whichever side rounding puts it, so the split is unreliable exactly
    // where it would be most interesting. And a static LDL that perturbed a pivot
    // (numPerturbations() > 0) reports the inertia of the matrix it actually factored, which is the
    // perturbed one.
    bool inertia(Inertia& inertia) const;

    // ||A x - b|| / ||b||, the one number that says whether the pipeline worked. Recomputes A x, so
    // it costs a multiplication; it is a convenience, not part of solving.
    double relativeResidual(const SparseMatrix<Val>& A, const Vector<Val>& b,
                            const Vector<Val>& x) const;

private:
    // Set by the constructor, which is where their defaults live, so they are not repeated here.
    Ordering                   mOrdering;
    Factorization              mFactorization;
    Traversal                  mTraversal;
    Supernodes                 mSupernodes;
    std::optional<std::size_t> mAmalgamation;

    double mPivotThreshold = 0.1;   // not a constructor argument: dynamic LDL tuning, as on the engine

    Permutation mPermutation;
    ElmForest   mElmForest;
    SymFactor   mSymFactor;

    // Both storages are declared, one is filled. An unused vector costs nothing, and this keeps the
    // factor a concrete member rather than something reached through a pointer or a variant.
    NumFactorStatic<Val>  mNumFactorStatic;
    NumFactorDynamic<Val> mNumFactorDynamic;

    std::size_t mSize               = 0;
    bool        mUsesDynamicStorage = false;
    bool        mAnalyzed           = false;
    bool        mFactored           = false;
};

extern template class DirectSolver<double>;
extern template class DirectSolver<std::complex<double>>;

} // namespace Oblio

#endif // OBLIO_DIRECT_SOLVER_H
