#pragma once

// NumFactorEngine.h - computes the numeric factorization of a SparseMatrix, under a given
// Permutation, from an already-computed SymFactor.
//
// Two settings, and between them they determine everything:
//
//   Factorization   what is computed (Cholesky first; static and dynamic LDL to follow)
//   Traversal       how the supernodes are walked (left-looking, right-looking; multifrontal
//                   to follow)
//
// Symmetry is *not* a setting. The scalar type and the factorization fix it: real is symmetric,
// complex Cholesky is Hermitian, complex LDL is complex-symmetric. Complex symmetric Cholesky is
// rejected, permanently, because positive definiteness is meaningless for it. See the
// factorization-space entry in DESIGN_DECISIONS, and BlasLapack.h, where that decision becomes
// code: the engine asks for the *operation* (`herk`, meaning A times A-conjugate-transpose) and
// never names the routine, so it cannot pick the wrong one, which is the mistake 0.9 makes.
//
// The kernels take (Val* block, rows, cols, ld) and never see a factor object. That is not
// fastidiousness: experiments/storage-options established that a plain array of pointers is
// enough to let one compiled algorithm serve both a flat and a per-supernode storage, so writing
// the kernels this way is what lets NumFactorDynamic reach them later without a line changing.
//
// Reads the matrix, permutation and symbolic factor through their public accessors; writes the
// numeric factor through friendship.

#include "oblio/NumFactorDynamic.h"
#include "oblio/NumFactorStatic.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"
#include "oblio/SymFactor.h"
#include "oblio/Types.h"
#include "oblio/UpdateBlock.h"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

class NumFactorEngine {
public:
    NumFactorEngine() = default;
    NumFactorEngine(Factorization factorization, Traversal traversal)
        : mFactorization(factorization), mTraversal(traversal) {}

    void          setFactorization(Factorization factorization) { mFactorization = factorization; }
    Factorization factorization() const                         { return mFactorization; }

    void      setTraversal(Traversal traversal) { mTraversal = traversal; }
    Traversal traversal() const                 { return mTraversal; }

    // The threshold for accepting a 1x1 (or 2x2) pivot in dynamic LDL: a candidate pivots in place
    // only if its magnitude is at least this fraction of the largest off-diagonal in its column,
    // otherwise it is delayed to an ancestor. 0.9's default, in (0, 1]; smaller delays fewer
    // columns (less fill) at the cost of stability.
    void   setPivotThreshold(double t) { mPivotThreshold = t; }
    double pivotThreshold() const      { return mPivotThreshold; }

    // The perturbation threshold, for LDL only.
    //
    // A static factorization does not pivot, so a pivot smaller than this has no remedy: it is
    // replaced by this value, and counted. The count comes back on the factor
    // (NumFactorStatic::numPerturbations), because it is a property of what was computed, not of
    // the engine that computed it. Cholesky ignores this: it fails on a bad pivot rather than
    // perturbing, which is what positive definiteness entitles it to do.
    void   setPerturbation(double perturbation) { mPerturbation = perturbation; }
    double perturbation() const                 { return mPerturbation; }

    // The public entry point. Returns false if the matrix is not positive definite (Cholesky
    // found a non-positive pivot), or if the matrix does not match the structure the symbolic
    // factorization predicted.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& sf,
                 NumFactorStatic<Val>& nf) const;

    // The same entry point, producing the dynamic factor. The static factorizations run into
    // per-supernode storage unchanged; the dynamic-LDL cases are the one thing this overload will
    // gain that its sibling cannot (a flat buffer cannot grow a front).
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& sf,
                 NumFactorDynamic<Val>& nf) const;

private:
    Factorization mFactorization = Factorization::Cholesky;
    Traversal     mTraversal     = Traversal::LeftLooking;
    double        mPivotThreshold = 0.1;
    double        mPerturbation  = 1e-14;

    // Copy the structure from SymFactor and allocate the value blocks, zeroed. Every traversal
    // starts here, so it is shared.
    //
    // The blocks are rectangles: indexSize rows by frontSize columns, column-major. The strictly
    // upper triangle of the front is allocated and left zero, which is what lets the whole block
    // go to a level-3 BLAS call with a leading dimension.
    template<class Val>
    void setSymFactor(const SymFactor& sf, NumFactorStatic<Val>& nf) const;

    // The dynamic twin: same structure copied, but the value blocks are one vector per supernode
    // rather than offsets into one buffer, so a front can later grow without moving its neighbors.
    template<class Val>
    void setSymFactor(const SymFactor& sf, NumFactorDynamic<Val>& nf) const;

    // A supernode's rows, in *its own* local coordinates: gblToLcl[globalRow] = local position.
    // Set before working on a supernode, cleared after, so the array is allocated once and reused
    // across all of them. The same discipline as symbolic's mark array, and for the same reason:
    // clearing only what was written keeps this out of the O(n) per supernode it would otherwise
    // cost.
    void setGlobalToLocal(std::size_t numIdx, const std::int32_t* rowIdx,
                          std::vector<std::int32_t>& gblToLcl) const;
    void clearGlobalToLocal(std::size_t numIdx, const std::int32_t* rowIdx,
                            std::vector<std::int32_t>& gblToLcl) const;

    // Scatter A's original values into a supernode's block. The only place A is read.
    //
    // `numberOfDelayedColumns` is where A's own columns begin. It is zero for every static
    // factorization; under dynamic pivoting a front is grown at the left by the columns its
    // children delayed into it, and those columns hold no entry of A, so the scatter starts past
    // them. The count is therefore a property of the *destination*, not of A.
    //
    // Assigns rather than accumulates, and so does assembleDelayed below: the block was zeroed and
    // each of these positions has exactly one writer. Only assembleUpdate accumulates, because only
    // it lands where A's values and other descendants' updates already sit. (0.9 writes `+=` here,
    // which differs only if A stores a duplicate entry, and A is assumed valid.)
    //
    // Returns false if A holds an entry the symbolic structure does not predict, which means the
    // two do not describe the same matrix. That is a caller error, not a numeric failure.
    template<class Val>
    bool assembleFromA(const SparseMatrix<Val>& A, const Permutation& p,
                       const std::vector<std::int32_t>& gblToLcl,
                       std::size_t numberOfDelayedColumns,
                       std::size_t frontSize, std::size_t numIdx,
                       const std::int32_t* rowIdx, Val* block) const;

    // Scatter an UpdateBlock into the ancestor it was formed for. Adds, does not overwrite: an
    // ancestor collects updates from many descendants.
    template<class Val>
    void assembleUpdate(const std::vector<std::int32_t>& gblToLcl,
                        const UpdateBlock<Val>& t,
                        std::size_t numIdx, Val* block) const;

    // Factor one supernode's block in place, and form the update it owes an ancestor. **These two
    // dispatch on the factorization type**, which is what keeps the traversals below identical
    // across Cholesky and LDL: left-looking and right-looking decide *when* to factor and *when*
    // to update, never *how*.
    //
    // Cholesky:
    //   factor  POTRF on the front, then TRSM on the update rows: L21 = A21 (L11^H)^-1.
    //   update  HERK on the square part, GEMM on the rectangle below. One rank-k call, because
    //           L21 L21^H is exactly what HERK computes.
    //
    // LDL:
    //   factor  our own kernel (LAPACK has no unpivoted LDL), then TRSM against the *upper*
    //           triangle, untransposed, because the block holds U = D L^T there.
    //   update  T -= L21 D L21^T, which no BLAS routine computes: the D in the middle rules out
    //           a rank-k call. So U := D L21^T into a scratch, then gemmLower on the symmetric
    //           part and GEMM on the rectangle. The scratch is the price of the D.
    //
    // factorStaticSupernode returns false only for Cholesky, on a non-positive pivot. LDL cannot
    // fail:
    // it perturbs instead, and reports how often.
    template<class Val>
    bool factorStaticSupernode(std::size_t frontSize, std::size_t numIdx, Val* block,
                         std::size_t& numPerturbations) const;

    template<class Val>
    void updateStaticSupernode(std::size_t frontSize, std::size_t numIdx, const Val* block,
                         std::size_t offset, UpdateBlock<Val>& t) const;

    // The traversals, named for the *pivoting*, which is the axis the Factorization enum names:
    // static pivoting is Cholesky and static LDL, dynamic pivoting is dynamic LDL. Left-looking and
    // right-looking are the same arithmetic in opposite directions:
    //
    //   Left-looking:  for each supernode, PULL every update owed to it, then factor.
    //                  Needs a list per supernode of who still owes it, and a position per
    //                  supernode tracking how far it has got through its ancestors.
    //   Right-looking: for each supernode, factor, then PUSH its update to every ancestor.
    //                  Needs no lists at all.
    //
    // **Static pivoting runs in either storage; dynamic pivoting requires the dynamic one**, since
    // delaying a column grows a front. So these two are templated on the factor, while the dynamic
    // traversal below names NumFactorDynamic outright. See dynamicPivoting() in Types.h.
    //
    // These two do not know which factorization they are running. They call
    // factorStaticSupernode and updateStaticSupernode, and those branch on mFactorization to
    // choose Cholesky or LDL.
    template<class Val, class Factor>
    bool factorStaticLeftLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                 const SymFactor& sf, Factor& nf) const;
    template<class Val, class Factor>
    bool factorStaticRightLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                  const SymFactor& sf, Factor& nf) const;

    // Dynamic LDL, left-looking: the full traversal. Grows a front by whatever its children could
    // not pivot, assembles A past those columns, folds each child's delayed columns into it before
    // shrinking them away, then factors, which may delay in turn. Real only so far, serving both
    // `DynamicLDLT` and `DynamicLDLH`, which are the same computation over the reals. Complex is
    // where the two part company; right-looking and multifrontal come after.
    // The two pivot eliminations, applied once a selection loop has accepted one. Shared by the
    // kernel's two passes, which differ only in *selection*: once a pivot is accepted the arithmetic
    // is the same, so it lives in one place. See the note above their definitions.
    template<class Val>
    void applyPivot1x1(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t j_,
                       std::int32_t k1_, std::int32_t k1, std::int32_t jjFrontSize,
                       std::int32_t rows, std::vector<std::int32_t>& gblToLcl) const;

    template<class Val>
    void applyPivot2x2(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t j_,
                       std::int32_t k1_, std::int32_t k2_, std::int32_t k1, std::int32_t k2,
                       std::int32_t jjFrontSize, std::int32_t rows,
                       std::vector<std::int32_t>& gblToLcl) const;

    template<class Val>
    bool factorDynamicLeftLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                  const SymFactor& sf, NumFactorDynamic<Val>& nf) const;

    // Dynamic LDL, right-looking. The same factorization by the same two kernels: 0.9's
    // factorDynamicLDL_ and updateDynamicLDL_ are byte-identical between its left- and
    // right-looking engines, so only the driver differs.
    //
    // Where it differs is growth. This traversal assembles A into every front at the start and
    // pushes each supernode's update into its ancestors as it goes, so a front that grows already
    // holds values and must keep them: it calls extendEntry where left-looking calls resetEntry.
    // That is the whole of the difference, and it is why extendEntry existed unported until now.
    template<class Val>
    bool factorDynamicRightLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                   const SymFactor& sf, NumFactorDynamic<Val>& nf) const;

    // Factor one supernode's dense front in place with threshold pivoting, delaying the columns it
    // cannot pivot to an ancestor. Records pivotType (1 / 2,3), numberOfDelayedColumns, and reduces
    // frontSize by the number delayed. The dynamic counterpart of factorStaticSupernode, and
    // unlike it this one takes the factor rather than a raw block: pivoting swaps columns and
    // records their kind, so it cannot be storage-blind. Ported from 0.9 factorDynamicLDL_
    // (updateSize == 0 pass).
    template<class Val>
    bool factorDynamicSupernode(NumFactorDynamic<Val>& nf, std::int32_t jj,
                                std::vector<std::int32_t>& gblToLcl) const;

    // Form the update supernode jj owes one ancestor, taking jj's rows from `offset` down. The
    // dynamic counterpart of updateStaticSupernode, and it differs from it in exactly two places.
    //
    // The leading dimension is the block's *height*, which still counts the delayed columns:
    // shrinkEntry dropped their columns and kept their rows, so a delayed row is a genuine row of
    // L and the stride steps over it.
    //
    // And D is no longer diagonal. The static twin hands the whole D L21^T scratch to formUpper in
    // one call; here a 2x2 pivot makes D block-diagonal, so the scratch is built by walking the
    // front and consulting pivotType to know whether to step one column or two. That walk is why
    // this takes the factor rather than a block: the pivot kind is stored per *global* node, so it
    // needs the index set to be read at all. Ported from 0.9 updateDynamicLDL_.
    template<class Val>
    void updateDynamicSupernode(const NumFactorDynamic<Val>& nf, std::int32_t jj,
                                std::size_t offset, UpdateBlock<Val>& t) const;

    // Fold supernode jj's delayed columns into kk, its parent, which has already been grown to
    // hold them. The third assemble, and the only one dynamic pivoting adds: A's values and a
    // descendant's update both land in a block the symbolic factorization predicted, while these
    // columns are the part it did not.
    //
    // jj's block still carries them at this point. factorDynamicSupernode has decremented
    // frontSize[jj] and set numberOfDelayedColumns[jj], so the columns are the run just past the
    // new front, and shrinkEntry has not yet reclaimed them. Which fixes the order: assemble, then
    // shrink.
    template<class Val>
    void assembleDelayed(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t kk,
                         const std::vector<std::int32_t>& gblToLcl) const;
};

} // namespace Oblio
