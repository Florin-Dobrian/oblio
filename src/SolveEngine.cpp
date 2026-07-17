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
// so entry (localRow, localCol) is at val[localCol * numNodeIdx + localRow], and the global row it
// stands for is nodeIdx[localRow].
//
// For Cholesky the diagonal holds L's own diagonal. For LDL it holds D, and L is *unit* lower
// triangular, its 1s implicit. That single difference is the whole of what separates the two
// solves, and it is why the diagonal is a separate pass for one and not the other.
// =================================================================================================

template<class Val>
void SolveEngine::forward(const NumFactorStatic<Val>& nf, Vector<Val>& y) const {
    const bool unitDiagonal = separateDiagonal(nf.factorization());   // LDL: L has implicit 1s

    // Ascending supernode order is a topological order, so a supernode's columns are finished
    // before any supernode below it needs them.
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(nf.snodeSize()); ++kk) {
        const std::size_t   frontSize  = nf.frontSize()[kk];
        const std::size_t   updateSize = nf.updateSize()[kk];
        const std::size_t   numNodeIdx = frontSize + updateSize;
        const std::int32_t* nodeIdx    = nf.nodeIdxPtr(kk);
        const Val*          val        = nf.valPtr(kk);

        for (std::int32_t j = 0; j < static_cast<std::int32_t>(frontSize); ++j) {
            const std::int32_t lj = nodeIdx[j];               // the global column
            const Val*         col = val + j * numNodeIdx;

            // Divide by the diagonal, unless L is unit (LDL), where the diagonal holds D and is
            // dealt with in its own pass.
            if (!unitDiagonal)
                y.mVal[lj] /= col[j];

            // Scatter the column's contribution down. Note the rows run to the *end of the index
            // set*, not the end of the front: a supernode's update rows are exactly the rows of L
            // below its own columns, and they belong to supernodes not yet reached.
            const Val yj = y.mVal[lj];
            for (std::int32_t i = j + 1; i < static_cast<std::int32_t>(numNodeIdx); ++i)
                y.mVal[nodeIdx[i]] -= col[i] * yj;
        }
    }
}

template<class Val>
void SolveEngine::diagonal(const NumFactorStatic<Val>& nf, Vector<Val>& y) const {
    // D z = y. LDL only, and while every pivot is 1x1 this is one division per column. Dynamic LDL
    // will bring 2x2 pivots, and then a pair of columns is solved together against a 2x2 block.
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(nf.snodeSize()); ++kk) {
        const std::size_t   frontSize  = nf.frontSize()[kk];
        const std::size_t   updateSize = nf.updateSize()[kk];
        const std::size_t   numNodeIdx = frontSize + updateSize;
        const std::int32_t* nodeIdx    = nf.nodeIdxPtr(kk);
        const Val*          val        = nf.valPtr(kk);

        for (std::int32_t j = 0; j < static_cast<std::int32_t>(frontSize); ++j)
            y.mVal[nodeIdx[j]] /= val[j * numNodeIdx + j];
    }
}

template<class Val>
void SolveEngine::backward(const NumFactorStatic<Val>& nf, Vector<Val>& y) const {
    const bool unitDiagonal = separateDiagonal(nf.factorization());
    const bool conj         = hermitian(nf.factorization());

    // Descending, the mirror of the forward pass: a supernode's columns are solved only once
    // everything below them is known.
    for (std::int32_t kk = static_cast<std::int32_t>(nf.snodeSize()) - 1; kk >= 0; --kk) {
        const std::size_t   frontSize  = nf.frontSize()[kk];
        const std::size_t   updateSize = nf.updateSize()[kk];
        const std::size_t   numNodeIdx = frontSize + updateSize;
        const std::int32_t* nodeIdx    = nf.nodeIdxPtr(kk);
        const Val*          val        = nf.valPtr(kk);

        for (std::int32_t j = static_cast<std::int32_t>(frontSize) - 1; j >= 0; --j) {
            const std::int32_t lj  = nodeIdx[j];
            const Val*         col = val + j * numNodeIdx;

            // Gather the contributions from below. **The conjugate is the point.** This pass
            // applies L^H, not L^T, whenever the factorization is Hermitian (Cholesky, LDLH), and
            // L^T when it is not (LDLT). 10.12 omits it, which is right for its complex-symmetric
            // LDL and wrong for its Cholesky, and nothing at the call site reveals that.
            Val acc = y.mVal[lj];
            for (std::int32_t i = j + 1; i < static_cast<std::int32_t>(numNodeIdx); ++i)
                acc -= maybeConj(col[i], conj) * y.mVal[nodeIdx[i]];

            // And divide, unless L is unit, in which case D was dealt with already.
            y.mVal[lj] = unitDiagonal ? acc : acc / maybeConj(col[j], conj);
        }
    }
}

template<class Val>
bool SolveEngine::compute(const NumFactorStatic<Val>& nf, Vector<Val>& y) const {
    if (y.size() != nf.size())
        return false;

    forward(nf, y);
    if (separateDiagonal(nf.factorization()))
        diagonal(nf, y);
    backward(nf, y);
    return true;
}

template<class Val>
bool SolveEngine::compute(const Permutation& p, const NumFactorStatic<Val>& nf,
                          const Vector<Val>& b, Vector<Val>& x) const {
    const std::size_t size = nf.size();
    if (b.size() != size || p.size() != size)
        return false;

    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // The factor is of P A P^T, so the right-hand side must be permuted into its ordering, and the
    // answer permuted back. Row li of the permuted system is row newToOld[li] of the original.
    Vector<Val> y(size);
    for (std::int32_t ai = 0; ai < static_cast<std::int32_t>(size); ++ai)
        y.mVal[oldToNew[ai]] = b.mVal[ai];

    if (!compute(nf, y))
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
