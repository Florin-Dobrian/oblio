#include "oblio/OrderEngine.h"

#include <vector>

// Entry points of the vendored ordering codes (raw int CSC arrays).
extern void mmd_order(int n, const int colPtr[], const int rowIdx[],
                      int perm[], int invp[]);                         // src/Mmd.cpp
extern "C" int amd_order(int n, const int Ap[], const int Ai[],
                         int P[], double Control[], double Info[]);    // src/Amd.cpp

namespace Oblio {

template<class Val>
bool OrderEngine::order(const SparseMatrix<Val>& A, Permutation& p) const {
    const std::size_t size = A.size();

    if (mMethod == OrderMethod::Natural)
        return orderNatural(size, p);

    // Read A's structure through the public API, bound once (no per-element call).
    const std::vector<std::size_t>& colPtrA = A.colPtr();
    const std::vector<std::int32_t>& rowIdxA = A.rowIdx();

    // Non-natural: size the maps; the algorithms fill them by index below.
    p.mOldToNew.assign(size, 0);
    p.mNewToOld.assign(size, 0);

    if (mMethod == OrderMethod::AMD)
        // A is full-symmetric; AMD ignores the diagonal and symmetrizes internally,
        // so its structure can be passed directly.
        return orderAMD(size, colPtrA, rowIdxA, p);

    // A is stored full-symmetric; MMD wants the off-diagonal structure only.
    // Strip the diagonal (no expansion needed — A already holds both triangles).
    std::vector<std::size_t> colPtr(size + 1, 0);
    for (std::size_t j = 0; j < size; ++j)
        for (std::size_t sp = colPtrA[j]; sp < colPtrA[j + 1]; ++sp)
            if (static_cast<std::size_t>(rowIdxA[sp]) != j) colPtr[j + 1]++;
    for (std::size_t j = 0; j < size; ++j) colPtr[j + 1] += colPtr[j];
    std::vector<std::int32_t> rowIdx(colPtr[size]);
    std::vector<std::size_t> cur(colPtr.begin(), colPtr.end());
    for (std::size_t j = 0; j < size; ++j)
        for (std::size_t sp = colPtrA[j]; sp < colPtrA[j + 1]; ++sp)
            if (static_cast<std::size_t>(rowIdxA[sp]) != j) rowIdx[cur[j]++] = rowIdxA[sp];
    return orderMMD(size, colPtr, rowIdx, p);
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
    if (size == 0) return true;   // maps already sized to 0 by order()
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
    if (size == 0) return true;   // maps already sized to 0 by order()
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
        p.mOldToNew[static_cast<std::size_t>(P[k])] = static_cast<std::int32_t>(k);
    }
    return true;
}

template bool OrderEngine::order(const SparseMatrix<double>&, Permutation&) const;
template bool OrderEngine::order(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
