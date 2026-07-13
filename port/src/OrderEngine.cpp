#include "oblio/OrderEngine.h"

#include <vector>

// Entry points of the vendored ordering codes (raw int CSC arrays).
extern void mmd_order(int n, const int colPtr[], const int rowIdx[],
                      int perm[], int invp[]);                         // src/Mmd.cpp
extern "C" int amd_order(int n, const int Ap[], const int Ai[],
                         int P[], double Control[], double Info[]);    // src/Amd.cpp

namespace Oblio {

// Adapter: an ordering needs only the sparsity pattern, so the matrix overload pulls it
// out and forwards. The implementation below is free of Val and compiled once.
template<class Val>
bool OrderEngine::compute(const SparseMatrix<Val>& A, Permutation& p) const {
    return compute(A.colPtr(), A.rowIdx(), p);
}

bool OrderEngine::compute(const std::vector<std::size_t>&  colPtr,
                          const std::vector<std::int32_t>& rowIdx,
                          Permutation& p) const {
    if (colPtr.empty())
        return false;
    const std::size_t size = colPtr.size() - 1;

    if (mMethod == OrderMethod::Natural)
        return orderNatural(size, p);

    // Non-natural: size the maps; the algorithms fill them by index below.
    p.mOldToNew.assign(size, 0);
    p.mNewToOld.assign(size, 0);

    if (mMethod == OrderMethod::AMD)
        // A is full-symmetric; AMD ignores the diagonal and symmetrizes internally,
        // so its structure can be passed directly.
        return orderAMD(size, colPtr, rowIdx, p);

    // A is stored full-symmetric; MMD wants the off-diagonal structure only.
    // Strip the diagonal (no expansion needed, A already holds both triangles). Columns are
    // indices, so aj is an int32_t and the comparison against rowIdx[cp] needs no cast; cp is
    // a position into A's arrays.
    std::vector<std::size_t> colPtrOff(size + 1, 0);
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(size); ++aj)
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) colPtrOff[aj + 1]++;
    for (std::size_t j = 0; j < size; ++j) colPtrOff[j + 1] += colPtrOff[j];
    std::vector<std::int32_t> rowIdxOff(colPtrOff[size]);
    std::vector<std::size_t> cur(colPtrOff.begin(), colPtrOff.end());
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(size); ++aj)
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) rowIdxOff[cur[aj]++] = rowIdx[cp];
    return orderMMD(size, colPtrOff, rowIdxOff, p);
}

bool OrderEngine::orderNatural(std::size_t size, Permutation& p) const {
    p.mOldToNew.assign(size, 0);
    p.mNewToOld.assign(size, 0);
    p.setIdentity();
    return true;
}

bool OrderEngine::orderMMD(std::size_t size,
                           const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx,
                           Permutation& p) const {
    if (size == 0) return true;   // maps already sized to 0 by compute()

    // Crossing into the vendored C API, which is int-based throughout. These casts are not the
    // index/position crossings of our own type rules; they are the boundary of a foreign
    // interface, and the arrays below exist only to feed it.
    const int N   = static_cast<int>(size);
    const int nnz = static_cast<int>(rowIdx.size());

    std::vector<int> cp(N + 1), ri(nnz);
    for (int j = 0; j <= N; ++j) cp[j] = static_cast<int>(colPtr[j]);
    for (int k = 0; k < nnz; ++k) ri[k] = static_cast<int>(rowIdx[k]);

    std::vector<int> perm(N), invp(N);
    mmd_order(N, cp.data(), ri.data(), perm.data(), invp.data());

    for (int j = 0; j < N; ++j) {
        p.mOldToNew[j] = static_cast<std::int32_t>(invp[j]);
        p.mNewToOld[j] = static_cast<std::int32_t>(perm[j]);
    }
    return true;
}

bool OrderEngine::orderAMD(std::size_t size,
                           const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx,
                           Permutation& p) const {
    if (size == 0) return true;   // maps already sized to 0 by compute()
    const int N  = static_cast<int>(size);
    const int nz = static_cast<int>(colPtr[size]);

    std::vector<int> Ap(N + 1), Ai(nz);
    for (int j = 0; j <= N; ++j) Ap[j] = static_cast<int>(colPtr[j]);
    for (int k = 0; k < nz; ++k) Ai[k] = static_cast<int>(rowIdx[k]);

    std::vector<int> P(N);
    const int status = amd_order(N, Ap.data(), Ai.data(), P.data(), nullptr, nullptr);
    if (status < 0) return false;

    for (int k = 0; k < N; ++k) {
        p.mNewToOld[k]    = static_cast<std::int32_t>(P[k]);
        p.mOldToNew[P[k]] = static_cast<std::int32_t>(k);
    }
    return true;
}

template bool OrderEngine::compute(const SparseMatrix<double>&, Permutation&) const;
template bool OrderEngine::compute(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
