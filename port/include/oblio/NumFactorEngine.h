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

    // The public entry point. Returns false if the matrix is not positive definite (Cholesky
    // found a non-positive pivot), or if the matrix does not match the structure the symbolic
    // factorization predicted.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& s,
                 NumFactorStatic<Val>& f) const;

private:
    Factorization mFactorization = Factorization::Cholesky;
    Traversal     mTraversal     = Traversal::LeftLooking;

    // Copy the structure from SymFactor and allocate the value blocks, zeroed. Every traversal
    // starts here, so it is shared.
    //
    // The blocks are rectangles: indexSize rows by frontSize columns, column-major. The strictly
    // upper triangle of the front is allocated and left zero, which is what lets the whole block
    // go to a level-3 BLAS call with a leading dimension.
    template<class Val>
    void setSymFactor(const SymFactor& s, NumFactorStatic<Val>& f) const;

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
    // Returns false if A holds an entry the symbolic structure does not predict, which means the
    // two do not describe the same matrix. That is a caller error, not a numeric failure.
    template<class Val>
    bool assembleFromA(const SparseMatrix<Val>& A, const Permutation& p,
                       const std::vector<std::int32_t>& gblToLcl,
                       std::size_t frontSize, std::size_t numIdx,
                       const std::int32_t* rowIdx, Val* block) const;

    // Scatter an UpdateBlock into the ancestor it was formed for. Adds, does not overwrite: an
    // ancestor collects updates from many descendants.
    template<class Val>
    void assembleUpdate(const std::vector<std::int32_t>& gblToLcl,
                        const UpdateBlock<Val>& t,
                        std::size_t numIdx, Val* block) const;

    // Factor one supernode's block in place.
    //
    //   POTRF on the front (the diagonal block), giving L11.
    //   TRSM  on the update rows, giving L21 = A21 (L11^H)^-1.
    //
    // The conjugate-transpose character comes from Blas<Val>, so this is correct for real and for
    // complex Hermitian without a branch. Returns false on a non-positive pivot.
    template<class Val>
    bool factorSupernode(std::size_t frontSize, std::size_t numIdx, Val* block) const;

    // Form the update that a factored supernode contributes to one ancestor, into t.
    //
    //   HERK on the square part:      t -= L21' L21'^H
    //   GEMM on the rectangle below:  t -= L21'' L21'^H
    //
    // where L21' is the `width` rows of the supernode's update block that land in the ancestor,
    // and L21'' the rows below them. Two calls rather than one because the diagonal part is
    // symmetric and HERK is the routine that knows it.
    template<class Val>
    void updateSupernode(std::size_t frontSize, std::size_t numIdx, const Val* block,
                         std::size_t offset, UpdateBlock<Val>& t) const;

    // The two traversals. Same arithmetic, opposite direction.
    //
    //   Left-looking:  for each supernode, PULL every update owed to it, then factor.
    //                  Needs a list per supernode of who still owes it, and a position per
    //                  supernode tracking how far it has got through its ancestors.
    //   Right-looking: for each supernode, factor, then PUSH its update to every ancestor.
    //                  Needs no lists at all.
    template<class Val>
    bool factorLeftLooking(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& s,
                           NumFactorStatic<Val>& f) const;
    template<class Val>
    bool factorRightLooking(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& s,
                            NumFactorStatic<Val>& f) const;
};

} // namespace Oblio
