#include "oblio/SolveEngine.h"

#include "oblio/NumFactorDynamic.h"
#include "oblio/NumFactorStatic.h"

#include <complex>

namespace Oblio {

namespace {

// Conjugate, or not, according to the factorization. For double it is a no-op in both branches,
// which is the honest expression of "real symmetric and real Hermitian are the same case".
double maybeConjugate(double val, bool)                                       { return val; }
std::complex<double> maybeConjugate(std::complex<double> val, bool hermitian) { return hermitian ? std::conj(val) : val; }

} // namespace

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
void SolveEngine::forward(const Factor& nf, Vector<Val>& y) const {
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
void SolveEngine::diagonal(const Factor& nf, Vector<Val>& y) const {
    // D z = y. LDL only, and while every pivot is 1x1 this is one division per column. Dynamic LDL
    // will bring 2x2 pivots, and then a pair of columns is solved together against a 2x2 block.
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
void SolveEngine::backward(const Factor& nf, Vector<Val>& y) const {
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

template<class Val, class Factor>
bool SolveEngine::compute(const Factor& nf, Vector<Val>& y) const {
    if (y.size() != nf.size())
        return false;

    forward(nf, y);
    if (separateDiagonal(nf.factorization()))
        diagonal(nf, y);
    backward(nf, y);
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
