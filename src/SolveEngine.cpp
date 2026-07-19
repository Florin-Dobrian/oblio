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
        const std::size_t   frontSize  = nf.frontSize(jj);
        const std::size_t   updateSize = nf.updateSize(jj);
        const std::size_t   numNodeIdx = frontSize + updateSize;
        const std::int32_t* nodeIdx    = nf.nodeIdx(jj);
        const Val*          val        = nf.val(jj);

        for (std::int32_t j = 0; j < static_cast<std::int32_t>(frontSize); ++j) {
            const std::int32_t lj = nodeIdx[j];   // the global column
            const std::size_t  cp = static_cast<std::size_t>(j) * numNodeIdx;

            // Divide by the diagonal, unless L is unit (LDL), where the diagonal holds D and is
            // dealt with in its own pass.
            if (!withSeparateDiagonal)
                y.mVal[lj] /= val[cp + j];

            // Scatter the column's contribution down. Note the rows run to the *end of the index
            // set*, not the end of the front: a supernode's update rows are exactly the rows of L
            // below its own columns, and they belong to supernodes not yet reached.
            const Val yj = y.mVal[lj];
            for (std::int32_t i = j + 1; i < static_cast<std::int32_t>(numNodeIdx); ++i)
                y.mVal[nodeIdx[i]] -= val[cp + i] * yj;
        }
    }
}

template<class Val, class Factor>
void SolveEngine::diagonalStatic(const Factor& nf, Vector<Val>& y) const {
    // D z = y. LDL only, and every pivot here is 1x1 by construction: a static factorization does
    // not pivot, so it never forms a 2x2 block. That is diagonalDynamic's business.
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(nf.snodeSize()); ++jj) {
        const std::size_t   frontSize  = nf.frontSize(jj);
        const std::size_t   updateSize = nf.updateSize(jj);
        const std::size_t   numNodeIdx = frontSize + updateSize;
        const std::int32_t* nodeIdx    = nf.nodeIdx(jj);
        const Val*          val        = nf.val(jj);

        for (std::int32_t j = 0; j < static_cast<std::int32_t>(frontSize); ++j)
            y.mVal[nodeIdx[j]] /= val[j * numNodeIdx + j];
    }
}

template<class Val, class Factor>
void SolveEngine::backwardStatic(const Factor& nf, Vector<Val>& y) const {
    const bool withSeparateDiagonal = separateDiagonal(nf.factorization());
    const bool withHermitian        = hermitian(nf.factorization());

    // Descending, the mirror of the forward pass: a supernode's columns are solved only once
    // everything below them is known.
    for (std::int32_t jj = static_cast<std::int32_t>(nf.snodeSize()) - 1; jj >= 0; --jj) {
        const std::size_t   frontSize  = nf.frontSize(jj);
        const std::size_t   updateSize = nf.updateSize(jj);
        const std::size_t   numNodeIdx = frontSize + updateSize;
        const std::int32_t* nodeIdx    = nf.nodeIdx(jj);
        const Val*          val        = nf.val(jj);

        for (std::int32_t j = static_cast<std::int32_t>(frontSize) - 1; j >= 0; --j) {
            const std::int32_t lj = nodeIdx[j];
            const std::size_t  cp = static_cast<std::size_t>(j) * numNodeIdx;

            // Gather the contributions from below. **The conjugate is the point.** This pass
            // applies L^H, not L^T, whenever the factorization is Hermitian (Cholesky, LDLH), and
            // L^T when it is not (LDLT). 10.12 omits it, which is right for its complex-symmetric
            // LDL and wrong for its Cholesky, and nothing at the call site reveals that.
            Val acc = y.mVal[lj];
            for (std::int32_t i = j + 1; i < static_cast<std::int32_t>(numNodeIdx); ++i)
                acc -= maybeConjugate(val[cp + i], withHermitian) * y.mVal[nodeIdx[i]];

            // And divide, unless L is unit, in which case D was dealt with already. The divide
            // branch is Cholesky-only, so the diagonal is real and C^H's conjugation is the
            // identity: no maybeConjugate here, unlike the off-diagonal above.
            y.mVal[lj] = withSeparateDiagonal ? acc : acc / val[cp + j];
        }
    }
}


// =================================================================================================
// The dynamic three. Same three solves, against a factor whose fronts moved while it was computed.
//
// A delayed column left its row behind: shrinkEntry reclaimed its column and kept every row, so the
// leading dimension is frontSize + numberOfDelayedColumns + updateSize while the columns to solve
// are only the frontSize of them. And a 2x2 pivot puts D's off-diagonal in the slot immediately
// below a diagonal, where L's first sub-diagonal entry would otherwise be, so the triangular passes
// step over it and the diagonal pass takes the pair together.
//
// pivotType is per *global* node: 1 for a 1x1, 2 for the first column of a 2x2, 3 for its second.
// Only 2 changes a loop bound, which is why the tests below read `!= 2` rather than enumerating.
// =================================================================================================

template<class Val>
void SolveEngine::forwardDynamic(const NumFactorDynamic<Val>& nf, Vector<Val>& y) const {
    const std::vector<std::int32_t>& pivotType = nf.pivotType();

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(nf.snodeSize()); ++jj) {
        const std::int32_t  frontSize  = static_cast<std::int32_t>(nf.frontSize(jj));
        const std::int32_t  rows       = frontSize + nf.numberOfDelayedColumns(jj)
                                       + static_cast<std::int32_t>(nf.updateSize(jj));
        const std::int32_t* nodeIdx    = nf.nodeIdx(jj);
        const Val*          val        = nf.val(jj);

        for (std::int32_t j = 0; j < frontSize; ++j) {
            const std::int32_t  lj = nodeIdx[j];
            const std::ptrdiff_t cp = static_cast<std::ptrdiff_t>(j) * rows;

            // L is unit, so no division. Scatter down, starting one row lower where this column
            // opens a 2x2: that row holds D's off-diagonal, not an entry of L.
            const Val yj = y.mVal[lj];
            for (std::int32_t i = (pivotType[lj] != 2 ? j + 1 : j + 2); i < rows; ++i)
                y.mVal[nodeIdx[i]] -= val[cp + i] * yj;
        }
    }
}

template<class Val>
void SolveEngine::diagonalDynamic(const NumFactorDynamic<Val>& nf, Vector<Val>& y) const {
    const std::vector<std::int32_t>& pivotType = nf.pivotType();

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(nf.snodeSize()); ++jj) {
        const std::int32_t  frontSize = static_cast<std::int32_t>(nf.frontSize(jj));
        const std::int32_t  rows      = frontSize + nf.numberOfDelayedColumns(jj)
                                      + static_cast<std::int32_t>(nf.updateSize(jj));
        const std::int32_t* nodeIdx   = nf.nodeIdx(jj);
        const Val*          val       = nf.val(jj);

        const auto at = [rows](std::int32_t r, std::int32_t c) {
            return static_cast<std::ptrdiff_t>(c) * rows + static_cast<std::ptrdiff_t>(r);
        };

        for (std::int32_t j = 0; j < frontSize; ) {
            const std::int32_t lj = nodeIdx[j];

            if (pivotType[lj] == 1) {
                // A zero pivot is left alone rather than treated as a failure. It means the pivot
                // search found nothing to eliminate against, and the system is singular or
                // inconsistent; 0.9 has the inconsistency error written and commented out, and we
                // keep its live behavior.
                if (val[at(j, j)] != Val(0))
                    y.mVal[lj] /= val[at(j, j)];
                ++j;
                continue;
            }

            // A 2x2 block, solved by an explicit LU with partial pivoting. D's four entries sit
            // where the factorization left them: the lower two are the original matrix entries,
            // untouched because the elimination starts two rows down, and the upper one was
            // written back when the pivot was accepted.
            const std::int32_t j1 = lj;
            const std::int32_t j2 = nodeIdx[j + 1];

            Val a11, a12, a21, a22, b1, b2;
            if (std::abs(val[at(j, j)]) >= std::abs(val[at(j + 1, j)])) {
                a11 = val[at(j, j)];         a12 = val[at(j, j + 1)];
                a21 = val[at(j + 1, j)];     a22 = val[at(j + 1, j + 1)];
                b1  = y.mVal[j1];            b2  = y.mVal[j2];
            } else {                          // pivot the rows: the second is the larger
                a11 = val[at(j + 1, j)];     a12 = val[at(j + 1, j + 1)];
                a21 = val[at(j, j)];         a22 = val[at(j, j + 1)];
                b1  = y.mVal[j2];            b2  = y.mVal[j1];
            }

            const Val l21 = a21 / a11;
            const Val u11 = a11;
            const Val u12 = a12;
            const Val u22 = a22 - l21 * u12;

            const Val y1 = b1;
            const Val y2 = b2 - l21 * y1;

            const Val x2 = y2 / u22;
            const Val x1 = (y1 - u12 * x2) / u11;

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
        const std::int32_t  frontSize = static_cast<std::int32_t>(nf.frontSize(jj));
        const std::int32_t  rows      = frontSize + nf.numberOfDelayedColumns(jj)
                                      + static_cast<std::int32_t>(nf.updateSize(jj));
        const std::int32_t* nodeIdx   = nf.nodeIdx(jj);
        const Val*          val       = nf.val(jj);

        for (std::int32_t i = frontSize - 1; i >= 0; --i) {
            const std::int32_t   li = nodeIdx[i];
            const std::ptrdiff_t cp = static_cast<std::ptrdiff_t>(i) * rows;

            // Gather from below, skipping the 2x2's own off-diagonal exactly as forwardDynamic
            // does. The conjugate is the same rule as in backwardStatic: L^H where the
            // factorization is Hermitian, L^T where it is not. Identity for DynamicLDLT, which is
            // all that runs today, and correct in advance for DynamicLDLH.
            Val acc = y.mVal[li];
            for (std::int32_t j = (pivotType[li] != 2 ? i + 1 : i + 2); j < rows; ++j)
                acc -= maybeConjugate(val[cp + j], withHermitian) * y.mVal[nodeIdx[j]];

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
bool SolveEngine::compute(const Permutation& p, const Factor& nf,
                          const Vector<Val>& b, Vector<Val>& x) const {
    const std::size_t size = nf.size();
    if (b.size() != size || p.size() != size)
        return false;

    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

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
