#include "oblio/MultiplyEngine.h"

namespace Oblio {

template<class Val>
bool MultiplyEngine::compute(const SparseMatrix<Val>& A, const Vector<Val>& x,
                             Vector<Val>& y) const {
    const std::size_t size = A.size();
    if (x.size() != size)
        return false;

    y.mVal.assign(size, Val(0));
    y.mSize = size;

    const std::vector<std::size_t>&  colPtr = A.colPtr();
    const std::vector<std::int32_t>& rowIdx = A.rowIdx();
    const std::vector<Val>&          val    = A.val();

    // Column by column, scattering into y. A is stored full-symmetric, so every entry appears
    // once in its own column and nothing has to be reflected.
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(size); ++aj) {
        const Val xj = x.mVal[aj];
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            y.mVal[rowIdx[cp]] += val[cp] * xj;
    }
    return true;
}

template<class Val>
bool MultiplyEngine::residual(const SparseMatrix<Val>& A, const Vector<Val>& x,
                              const Vector<Val>& b, Vector<Val>& r) const {
    if (b.size() != A.size())
        return false;
    if (!compute(A, x, r))
        return false;
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(r.mSize); ++i)
        r.mVal[i] -= b.mVal[i];
    return true;
}

template bool MultiplyEngine::compute(const SparseMatrix<double>&, const Vector<double>&,
                                      Vector<double>&) const;
template bool MultiplyEngine::compute(const SparseMatrix<std::complex<double>>&,
                                      const Vector<std::complex<double>>&,
                                      Vector<std::complex<double>>&) const;
template bool MultiplyEngine::residual(const SparseMatrix<double>&, const Vector<double>&,
                                       const Vector<double>&, Vector<double>&) const;
template bool MultiplyEngine::residual(const SparseMatrix<std::complex<double>>&,
                                       const Vector<std::complex<double>>&,
                                       const Vector<std::complex<double>>&,
                                       Vector<std::complex<double>>&) const;

} // namespace Oblio
