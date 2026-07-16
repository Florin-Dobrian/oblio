#include "oblio/SolveEngine.h"

#include <complex>

namespace Oblio {

namespace {

// Conjugate, or not, according to the factorization. For double it is a no-op in both branches,
// which is the honest expression of "real symmetric and real Hermitian are the same case".
double maybeConj(double v, bool)                              { return v; }
std::complex<double> maybeConj(std::complex<double> v, bool h) { return h ? std::conj(v) : v; }

} // namespace

// =================================================================================================
// The block. Supernode kk holds a dense column-major rectangle:
//
//     rows    = frontSize + updateSize   = the index-set size
//     columns = frontSize
//     ld      = the index-set size
//
// so entry (localRow, localCol) is at block[localCol * numIdx + localRow], and the global row it
// stands for is rowIdx[localRow].
//
// For Cholesky the diagonal holds L's own diagonal. For LDL it holds D, and L is *unit* lower
// triangular, its 1s implicit. That single difference is the whole of what separates the two
// solves, and it is why the diagonal is a separate pass for one and not the other.
// =================================================================================================

template<class Val>
void SolveEngine::forward(const NumFactorStatic<Val>& f, Vector<Val>& y) const {
    const bool unitDiagonal = separateDiagonal(f.mFactorization);   // LDL: L has implicit 1s

    // Ascending supernode order is a topological order, so a supernode's columns are finished
    // before any supernode below it needs them.
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(f.mSupSize); ++kk) {
        const std::size_t   frontSize = f.mFrontSize[kk];
        const std::size_t   numIdx    = frontSize + f.mUpdateSize[kk];
        const std::int32_t* rowIdx    = f.mRowIdx.data() + f.mSupPtr[kk];
        const Val*          block     = f.blockPtr(kk);

        for (std::size_t lclCol = 0; lclCol < frontSize; ++lclCol) {
            const std::int32_t lk = rowIdx[lclCol];               // the global column
            const Val*         col = block + lclCol * numIdx;

            // Divide by the diagonal, unless L is unit (LDL), where the diagonal holds D and is
            // dealt with in its own pass.
            if (!unitDiagonal)
                y.mVal[lk] /= col[lclCol];

            // Scatter the column's contribution down. Note the rows run to the *end of the index
            // set*, not the end of the front: a supernode's update rows are exactly the rows of L
            // below its own columns, and they belong to supernodes not yet reached.
            const Val yk = y.mVal[lk];
            for (std::size_t lclRow = lclCol + 1; lclRow < numIdx; ++lclRow)
                y.mVal[rowIdx[lclRow]] -= col[lclRow] * yk;
        }
    }
}

template<class Val>
void SolveEngine::diagonal(const NumFactorStatic<Val>& f, Vector<Val>& y) const {
    // D z = y. LDL only, and while every pivot is 1x1 this is one division per column. Dynamic LDL
    // will bring 2x2 pivots, and then a pair of columns is solved together against a 2x2 block.
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(f.mSupSize); ++kk) {
        const std::size_t   frontSize = f.mFrontSize[kk];
        const std::size_t   numIdx    = frontSize + f.mUpdateSize[kk];
        const std::int32_t* rowIdx    = f.mRowIdx.data() + f.mSupPtr[kk];
        const Val*          block     = f.blockPtr(kk);

        for (std::size_t lclCol = 0; lclCol < frontSize; ++lclCol)
            y.mVal[rowIdx[lclCol]] /= block[lclCol * numIdx + lclCol];
    }
}

template<class Val>
void SolveEngine::backward(const NumFactorStatic<Val>& f, Vector<Val>& y) const {
    const bool unitDiagonal = separateDiagonal(f.mFactorization);
    const bool conj         = hermitian(f.mFactorization);

    // Descending, the mirror of the forward pass: a supernode's columns are solved only once
    // everything below them is known.
    for (std::int32_t kk = static_cast<std::int32_t>(f.mSupSize) - 1; kk >= 0; --kk) {
        const std::size_t   frontSize = f.mFrontSize[kk];
        const std::size_t   numIdx    = frontSize + f.mUpdateSize[kk];
        const std::int32_t* rowIdx    = f.mRowIdx.data() + f.mSupPtr[kk];
        const Val*          block     = f.blockPtr(kk);

        for (std::size_t s = frontSize; s > 0; --s) {
            const std::size_t  lclCol = s - 1;
            const std::int32_t lk     = rowIdx[lclCol];
            const Val*         col    = block + lclCol * numIdx;

            // Gather the contributions from below. **The conjugate is the point.** This pass
            // applies L^H, not L^T, whenever the factorization is Hermitian (Cholesky, LDLH), and
            // L^T when it is not (LDLT). 10.12 omits it, which is right for its complex-symmetric
            // LDL and wrong for its Cholesky, and nothing at the call site reveals that.
            Val acc = y.mVal[lk];
            for (std::size_t lclRow = lclCol + 1; lclRow < numIdx; ++lclRow)
                acc -= maybeConj(col[lclRow], conj) * y.mVal[rowIdx[lclRow]];

            // And divide, unless L is unit, in which case D was dealt with already.
            y.mVal[lk] = unitDiagonal ? acc : acc / maybeConj(col[lclCol], conj);
        }
    }
}

template<class Val>
bool SolveEngine::compute(const NumFactorStatic<Val>& f, Vector<Val>& y) const {
    if (y.size() != f.size())
        return false;

    forward(f, y);
    if (separateDiagonal(f.mFactorization))
        diagonal(f, y);
    backward(f, y);
    return true;
}

template<class Val>
bool SolveEngine::compute(const Permutation& p, const NumFactorStatic<Val>& f,
                          const Vector<Val>& b, Vector<Val>& x) const {
    const std::size_t size = f.size();
    if (b.size() != size || p.size() != size)
        return false;

    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // The factor is of P A P^T, so the right-hand side must be permuted into its ordering, and the
    // answer permuted back. Row li of the permuted system is row newToOld[li] of the original.
    Vector<Val> y(size);
    for (std::int32_t ai = 0; ai < static_cast<std::int32_t>(size); ++ai)
        y.mVal[oldToNew[ai]] = b.mVal[ai];

    if (!compute(f, y))
        return false;

    x.mVal.assign(size, Val(0));
    x.mSize = size;
    for (std::int32_t li = 0; li < static_cast<std::int32_t>(size); ++li)
        x.mVal[newToOld[li]] = y.mVal[li];

    return true;
}

template bool SolveEngine::compute(const NumFactorStatic<double>&, Vector<double>&) const;
template bool SolveEngine::compute(const NumFactorStatic<std::complex<double>>&,
                                   Vector<std::complex<double>>&) const;
template bool SolveEngine::compute(const Permutation&, const NumFactorStatic<double>&,
                                   const Vector<double>&, Vector<double>&) const;
template bool SolveEngine::compute(const Permutation&,
                                   const NumFactorStatic<std::complex<double>>&,
                                   const Vector<std::complex<double>>&,
                                   Vector<std::complex<double>>&) const;

} // namespace Oblio
