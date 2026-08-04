#include "oblio/SolveEngine.h"

#include "oblio/NumFactorDynamic.h"
#include "oblio/NumFactorStatic.h"

#include <cmath>
#include <complex>
#include <type_traits>

namespace Oblio {

// maybeConjugate lives in Types.h, beside the `hermitian` predicate that drives it, because the
// dense kernels, the dynamic pivot code and this solve all need the same thing.

// =================================================================================================
// The supernodal block. Supernode jj holds a dense column-major rectangle:
//
//     number of rows    = frontSize + updateSize
//     number of columns = frontSize
//
// Columns are stored full-height with nothing between them, so the stride from one to the next is
// the index-set size: entry (i, j) is at val[j * numNodeIdx + i], and the global row it stands
// for is nodeIdx[i].
//
// For Cholesky the diagonal holds C's own diagonal. For LDL it holds D, and L is *unit* lower
// triangular, its 1s implicit. That single difference is the whole of what separates the two
// solves, and it is why the diagonal is a separate pass for one and not the other.
// =================================================================================================

template<class Val, class Factor>
void SolveEngine::forwardStatic(const Factor& nf, Vector<Val>& y) const {
    const bool withSeparateDiagonal = separateDiagonal(nf.factorization());

    // Ascending supernode order is a topological order, so a supernode's columns are finished
    // before any supernode below it needs them.
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(nf.snodeSize()); ++jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
        const Val*          jjVal        = nf.val(jj);

        for (std::int32_t j = 0; j < static_cast<std::int32_t>(jjFrontSize); ++j) {
            const std::int32_t lj = jjNodeIdx[j];   // the global column
            const std::size_t  cp = j * jjNumNodeIdx;

            // Divide by the diagonal, unless L is unit (LDL), where the diagonal holds D and is
            // dealt with in its own pass.
            if (!withSeparateDiagonal)
                y.mVal[lj] /= jjVal[cp + j];

            // Scatter the column's contribution down. Note the rows run to the *end of the index
            // set*, not the end of the front: a supernode's update rows are exactly the rows of L
            // below its own columns, and they belong to supernodes not yet reached.
            const Val yjVal = y.mVal[lj];
            for (std::int32_t i = j + 1; i < static_cast<std::int32_t>(jjNumNodeIdx); ++i)
                y.mVal[jjNodeIdx[i]] -= jjVal[cp + i] * yjVal;
        }
    }
}

template<class Val, class Factor>
void SolveEngine::diagonalStatic(const Factor& nf, Vector<Val>& y) const {
    // D z = y. LDL only, and every pivot here is 1x1 by construction: a static factorization does
    // not pivot, so it never forms a 2x2 block. That is diagonalDynamic's business.
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(nf.snodeSize()); ++jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
        const Val*          jjVal        = nf.val(jj);

        for (std::int32_t j = 0; j < static_cast<std::int32_t>(jjFrontSize); ++j)
            y.mVal[jjNodeIdx[j]] /= jjVal[j * jjNumNodeIdx + j];
    }
}

template<class Val, class Factor>
void SolveEngine::backwardStatic(const Factor& nf, Vector<Val>& y) const {
    const bool withSeparateDiagonal = separateDiagonal(nf.factorization());
    const bool withHermitian        = hermitian(nf.factorization());

    // Descending, the mirror of the forward pass: a supernode's columns are solved only once
    // everything below them is known.
    for (std::int32_t jj = static_cast<std::int32_t>(nf.snodeSize()) - 1; jj >= 0; --jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
        const Val*          jjVal        = nf.val(jj);

        for (std::int32_t j = static_cast<std::int32_t>(jjFrontSize) - 1; j >= 0; --j) {
            const std::int32_t lj = jjNodeIdx[j];
            const std::size_t  cp = j * jjNumNodeIdx;

            // Gather the contributions from below. **The conjugate is the point.** This pass
            // applies L^H, not L^T, whenever the factorization is Hermitian (Cholesky, LDLH), and
            // L^T when it is not (LDLT). 10.12 omits it, which is right for its complex-symmetric
            // LDL and wrong for its Cholesky, and nothing at the call site reveals that.
            Val acc = y.mVal[lj];
            for (std::int32_t i = j + 1; i < static_cast<std::int32_t>(jjNumNodeIdx); ++i)
                acc -= maybeConjugate(jjVal[cp + i], withHermitian) * y.mVal[jjNodeIdx[i]];

            // And divide, unless L is unit, in which case D was dealt with already. The divide
            // branch is Cholesky-only, so the diagonal is real and C^H's conjugation is the
            // identity: no maybeConjugate here, unlike the off-diagonal above.
            y.mVal[lj] = withSeparateDiagonal ? acc : acc / jjVal[cp + j];
        }
    }
}


// =================================================================================================
// The dynamic three. Same three solves, against a factor whose fronts moved while it was computed.
//
// A delayed column left its row behind: contractVal reclaimed its column and kept every row, so the
// leading dimension is jjFrontSize + delaySize + updateSize while the columns to solve
// are only the jjFrontSize of them. And a 2x2 pivot puts D's off-diagonal in the slot immediately
// below a diagonal, where L's first sub-diagonal entry would otherwise be, so the triangular passes
// step over it and the diagonal pass takes the pair together.
//
// Three conservation facts make these solves need no knowledge of what delayed where. The supernode
// count is fixed: delaying moves columns between fronts, never creates or destroys a supernode, so
// the loop runs over the same snodeSize it always would. Each front's row count is conserved
// (jjFrontSize + delaySize + updateSize, the height note in NumFactorEngine). And the total column
// count is conserved: every original column ends up as a surviving front column of exactly one
// supernode, the one where it finally pivoted, so summing jjFrontSize over all supernodes is still n.
// The delay history is therefore invisible here. The solve walks every supernode's *final* front,
// bounded by its settled jjFrontSize, and those finals partition the same n columns the symbolic
// fronts did. A supernode that gained columns does more work, one emptied by delay does none (its
// column loop runs zero times), and the sum is unchanged.
//
// pivotType is per *global* node: 1 for a 1x1, 2 for the first column of a 2x2, 3 for its second.
// Only 2 changes a loop bound, which is why the tests below read `!= 2` rather than enumerating.
// =================================================================================================

template<class Val>
void SolveEngine::forwardDynamic(const NumFactorDynamic<Val>& nf, Vector<Val>& y) const {
    const std::vector<std::int32_t>& pivotType = nf.pivotType();

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(nf.snodeSize()); ++jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.delaySize(jj) + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
        const Val*          jjVal        = nf.val(jj);

        for (std::int32_t j = 0; j < static_cast<std::int32_t>(jjFrontSize); ++j) {
            const std::int32_t  lj = jjNodeIdx[j];
            const std::size_t   cp = j * jjNumNodeIdx;

            // L is unit, so no division. Scatter down, starting one row lower where this column
            // opens a 2x2: that row holds D's off-diagonal, not an entry of L.
            const Val yjVal = y.mVal[lj];
            for (std::int32_t i = (pivotType[lj] != 2 ? j + 1 : j + 2); i < static_cast<std::int32_t>(jjNumNodeIdx); ++i)
                y.mVal[jjNodeIdx[i]] -= jjVal[cp + i] * yjVal;
        }
    }
}

template<class Val>
void SolveEngine::diagonalDynamic(const NumFactorDynamic<Val>& nf, Vector<Val>& y) const {
    const std::vector<std::int32_t>& pivotType = nf.pivotType();

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(nf.snodeSize()); ++jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.delaySize(jj) + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
        const Val*          jjVal        = nf.val(jj);

        const auto at = [jjNumNodeIdx](std::int32_t r, std::int32_t c) {
            return c * jjNumNodeIdx + r;
        };

        for (std::int32_t j = 0; j < static_cast<std::int32_t>(jjFrontSize); ) {
            const std::int32_t lj = jjNodeIdx[j];

            if (pivotType[lj] == 1) {
                // A zero pivot is left alone rather than treated as a failure. It reaches a 1x1
                // only by delay: a column that never found an acceptable pivot rides up the forest
                // and, at a root, pass 1 is forced to accept it as a 1x1 with whatever sits on the
                // diagonal (see the pivoting section of sparse_factorization.md). A zero there is
                // *expected* to mean A is singular, but that is conjectured, not proven: the clean
                // result holds for dense Bunch-Kaufman, whose forcing path is an empty column
                // (which is a zero row of A, hence singular); pass 1's forcing path is threshold
                // exhaustion instead, a different predicate, and the singularity argument has not
                // been shown to cross that gap. The divide is skipped either way to avoid inf/nan;
                // y[lj] keeps the value the forward sweep gave it, which is a valid choice of the
                // free null-space component if A is indeed singular here.
                //
                // What this does NOT distinguish is the two cases a zero pivot splits into. With
                // y[lj] == 0 the row reads 0 = 0: consistent, and leaving y[lj] alone is correct.
                // With y[lj] != 0 it reads 0 = nonzero: the system is inconsistent and has no
                // solution, and that is passed over silently here. 0.9 has exactly this test with a
                // SET_ERROR(InconsistentSystem) in the else branch, written and commented out; we
                // keep its live behavior. The check is disabled for a reason: at solve time y[lj]
                // has already been touched by the forward sweep, so a floating-point != 0 on it is
                // an unreliable inconsistency test, missing real cases and firing on rounded ones.
                // The honest diagnostic is the post-solve residual ||Ax - b|| / ||b||, which the
                // pipeline already computes; a flag here, if ever wanted, belongs there, not in
                // this per-pivot branch.
                if (jjVal[at(j, j)] != Val(0))
                    y.mVal[lj] /= jjVal[at(j, j)];
                ++j;
                continue;
            }

            // A 2x2 block, solved by an explicit LU with partial pivoting. D's four entries sit
            // where the factorization left them: the lower two are the original matrix entries,
            // untouched because the elimination starts two rows down, and the upper one was
            // written back when the pivot was accepted.
            //
            // Why LU and not a symmetric 2x2 solve, when the block is symmetric: a symmetric 2x2
            // can carry a zero diagonal (the [[0,1],[1,0]] shape), and a symmetric solve would try
            // to pivot on it and divide by zero. The row swap below (a11 = the larger of the two
            // first-column entries) moves a nonzero into the pivot instead. The asymmetry is the
            // safety: solving this "symmetrically" would reintroduce the division-by-zero it
            // avoids.
            //
            // Three divisions follow (l21 = a21/a11, x2 = y2/u22, x1 = .../u11), none guarded, and
            // whether that is safe depends on which pass accepted the block. Note u22 = det/a11 and
            // u11 = a11, so all three are nonzero exactly when a11 != 0 and det != 0.
            //   Pass 2 (has update rows): accepted on the threshold test |det| >= u * maxmax with
            //   u > 0, so det is bounded away from zero, and a11 is the larger first-column
            //   entry of a block whose off-diagonal is nonzero, so a11 != 0. Both divisions safe.
            //   Pass 1 (dense front): accepted on max1 == max2, the off-diagonal magnitudes alone,
            //   *without reading the determinant*. a11 != 0 still holds (a zero column exits as a
            //   1x1 at the max1 == 0 guard before any 2x2 forms, so a formed block has nonzero
            //   off-diagonals). But det is NOT bounded: a singular block with equal off-diagonals
            //   (e.g. [[d, d],[d, d]]) passes max1 == max2, reaches here with u22 = det/a11 = 0,
            //   and x2 = y2/u22 divides by zero, unguarded. This is precisely what dense
            //   Bunch-Kaufman structurally forbids: its case-4 2x2 forces det < 0 via the alpha
            //   constant, so BK never forms a singular 2x2 and singularity surfaces only as a zero
            //   1x1. Pass 1's max1 == max2 has no such guarantee, so it can construct the block BK
            //   cannot. See 7.8 of sparse_factorization.md. d == d exactly is measure-zero in
            //   floating point, so a random matrix will not hit it, but a structured one (a
            //   repeated KKT block, an assembled duplicate) could. Reachability is unverified; the
            //   fix, if it is reachable, is a determinant test on pass 1's 2x2 in the factor kernel
            //   (matching pass 2's, which is safe), not a guard here. Tracked in TODO alongside the
            //   pass-1 zero-pivot question, which is the same
            //   boundary seen from the 1x1 side.
            const std::int32_t j1 = lj;
            const std::int32_t j2 = jjNodeIdx[j + 1];

            // **A x = b**, with A = D and one right-hand side, solved by explicit LU with partial
            // pivoting. A column system to begin with, so A itself is what gets factored. factor2x2
            // in NumFactorEngine has the row system x A = b instead, which turns into a column
            // system in either A^T or A^H; it takes A^T in both cases, and see there for why. That
            // makes it these same four lines with a12 and a21 exchanged. The two are written to
            // look alike on purpose, and that exchange is the whole of what distinguishes them.
            Val  a11, a12, a21, a22;
            bool swapped;
            if (std::abs(jjVal[at(j, j)]) >= std::abs(jjVal[at(j + 1, j)])) {
                a11 = jjVal[at(j, j)];      a12 = jjVal[at(j, j + 1)];
                a21 = jjVal[at(j + 1, j)];  a22 = jjVal[at(j + 1, j + 1)];  swapped = false;
            } else {                                // pivot the rows: the second one is larger
                a11 = jjVal[at(j + 1, j)];  a12 = jjVal[at(j + 1, j + 1)];
                a21 = jjVal[at(j, j)];      a22 = jjVal[at(j, j + 1)];      swapped = true;
            }
            const Val l21 = a21 / a11;
            const Val u11 = a11;
            const Val u12 = a12;
            const Val u22 = a22 - l21 * u12;

            const Val t1 = y.mVal[j1];
            const Val t2 = y.mVal[j2];
            // the swap permutes the equations, not the unknowns
            const Val b1 = swapped ? t2 : t1;
            const Val b2 = swapped ? t1 : t2;

            const Val x2 = (b2 - l21 * b1) / u22;
            const Val x1 = (b1 - u12 * x2) / u11;

            y.mVal[j1] = x1;
            y.mVal[j2] = x2;

            j += 2;
        }
    }
}

template<class Val>
void SolveEngine::backwardDynamic(const NumFactorDynamic<Val>& nf, Vector<Val>& y) const {
    const bool withHermitian = hermitian(nf.factorization());

    const std::vector<std::int32_t>& pivotType = nf.pivotType();

    for (std::int32_t jj = static_cast<std::int32_t>(nf.snodeSize()) - 1; jj >= 0; --jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.delaySize(jj) + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
        const Val*          jjVal        = nf.val(jj);

        for (std::int32_t i = static_cast<std::int32_t>(jjFrontSize) - 1; i >= 0; --i) {
            const std::int32_t li = jjNodeIdx[i];
            const std::size_t  cp = i * jjNumNodeIdx;

            // Gather from below, skipping the 2x2's own off-diagonal exactly as forwardDynamic
            // does. The conjugate is the same rule as in backwardStatic: L^H where the
            // factorization is Hermitian, L^T where it is not. Identity for DynamicLDLT, which is
            // all that runs today, and correct in advance for DynamicLDLH.
            Val acc = y.mVal[li];
            for (std::int32_t j = (pivotType[li] != 2 ? i + 1 : i + 2); j < static_cast<std::int32_t>(jjNumNodeIdx); ++j)
                acc -= maybeConjugate(jjVal[cp + j], withHermitian) * y.mVal[jjNodeIdx[j]];

            y.mVal[li] = acc;   // L is unit: no division, D was dealt with in its own pass
        }
    }
}

template<class Val, class Factor>
bool SolveEngine::compute(const Factor& nf, Vector<Val>& y) const {
    if (y.size() != nf.size())
        return false;

    // Dynamic pivoting requires the dynamic storage, so for a static factor this branch is not
    // merely never taken, it is impossible, and `if constexpr` says so: the dynamic passes are not
    // instantiated for NumFactorStatic at all. Same rule NumFactorEngine follows, and the reason it
    // is a compile-time question rather than a runtime one is dynamicPivoting() in Types.h.
    if constexpr (std::is_same_v<Factor, NumFactorDynamic<Val>>) {
        if (dynamicPivoting(nf.factorization())) {
            forwardDynamic(nf, y);
            diagonalDynamic(nf, y);   // LDL always, and dynamic LDL is the only dynamic kind
            backwardDynamic(nf, y);
            return true;
        }
    }

    forwardStatic(nf, y);
    if (separateDiagonal(nf.factorization()))
        diagonalStatic(nf, y);
    backwardStatic(nf, y);
    return true;
}

template<class Val, class Factor>
bool SolveEngine::compute(const Permutation& P, const Factor& nf,
                          const Vector<Val>& b, Vector<Val>& x) const {
    const std::size_t size = nf.size();
    if (b.size() != size || P.size() != size)
        return false;

    const std::vector<std::int32_t>& oldToNew = P.oldToNew();
    const std::vector<std::int32_t>& newToOld = P.newToOld();

    // The factor is of P A P^T, so the right-hand side must be permuted into its ordering, and the
    // answer permuted back. Row lk of the permuted system is row newToOld[lk] of the original.
    Vector<Val> y(size);
    for (std::int32_t ak = 0; ak < static_cast<std::int32_t>(size); ++ak)
        y.mVal[oldToNew[ak]] = b.mVal[ak];

    if (!compute(nf, y))
        return false;

    x.mVal.assign(size, Val(0));
    x.mSize = size;
    for (std::int32_t lk = 0; lk < static_cast<std::int32_t>(size); ++lk)
        x.mVal[newToOld[lk]] = y.mVal[lk];

    return true;
}

template bool SolveEngine::compute(const NumFactorStatic<double>&, Vector<double>&) const;
template bool SolveEngine::compute(const NumFactorStatic<std::complex<double>>&,
                                   Vector<std::complex<double>>&) const;
template bool SolveEngine::compute(const NumFactorDynamic<double>&, Vector<double>&) const;
template bool SolveEngine::compute(const NumFactorDynamic<std::complex<double>>&,
                                   Vector<std::complex<double>>&) const;
template bool SolveEngine::compute(const Permutation&, const NumFactorStatic<double>&,
                                   const Vector<double>&, Vector<double>&) const;
template bool SolveEngine::compute(const Permutation&,
                                   const NumFactorStatic<std::complex<double>>&,
                                   const Vector<std::complex<double>>&,
                                   Vector<std::complex<double>>&) const;
template bool SolveEngine::compute(const Permutation&, const NumFactorDynamic<double>&,
                                   const Vector<double>&, Vector<double>&) const;
template bool SolveEngine::compute(const Permutation&,
                                   const NumFactorDynamic<std::complex<double>>&,
                                   const Vector<std::complex<double>>&,
                                   Vector<std::complex<double>>&) const;

} // namespace Oblio
