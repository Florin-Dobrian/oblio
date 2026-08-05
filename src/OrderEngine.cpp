#include "oblio/OrderEngine.h"
#include "oblio/Mmd1.h"
#include "oblio/Mmd2.h"
#include "oblio/Amd1.h"
#include "oblio/Amd1B.h"
#include "oblio/Amd2.h"
#include "oblio/Amd2B.h"

#include <vector>

// Entry points of the vendored ordering codes (raw int CSC arrays), declared only when those
// sources are being compiled. They live in private/, which is not part of the published tree, so
// the two orderings they implement are optional; see the Makefile and CMakeLists.txt for the
// detection, and docs/DESIGN_DECISIONS.md for why.
#ifdef OBLIO_VENDORED_ORDERINGS
extern void mmd_order(int n, const int colPtr[], const int rowIdx[],
                      int perm[], int invp[]);                         // private/Mmd.cpp
extern "C" int amd_order(int n, const int Ap[], const int Ai[],
                         int P[], double Control[], double Info[]);    // private/Amd.cpp
#endif

namespace Oblio {

// Adapter: an ordering needs only the sparsity pattern, so the matrix overload pulls it
// out and forwards. The implementation below is free of Val and compiled once.
template<class Val>
bool OrderEngine::compute(const SparseMatrix<Val>& A, Permutation& P) const {
    return compute(A.colPtr(), A.rowIdx(), P);
}

// Dispatch, and nothing else. A switch naming every enumerator with no default is the house
// rule, and here it is load-bearing rather than stylistic: this was an if-chain ending in an
// unguarded fall-through to MMD, so a new method would have become MMD silently and produced a
// valid permutation that nothing would have questioned. Each order* below sizes the maps it
// fills, so no step happens before the method is known.
bool OrderEngine::compute(const std::vector<std::size_t>&  colPtr,
                          const std::vector<std::int32_t>& rowIdx,
                          Permutation& P) const {
    if (colPtr.empty())
        return false;
    const std::size_t size = colPtr.size() - 1;

    switch (mOrdering) {
        case Ordering::Natural: return orderNatural(size, P);
#ifdef OBLIO_VENDORED_ORDERINGS
        case Ordering::MMD:     return orderMMD(size, colPtr, rowIdx, P);
#else
        case Ordering::MMD:     return false;   // vendored, and private/ is not present
#endif
        case Ordering::MMD1:    return orderMMD1(size, colPtr, rowIdx, P);
        case Ordering::MMD2:    return orderMMD2(size, colPtr, rowIdx, P);
#ifdef OBLIO_VENDORED_ORDERINGS
        case Ordering::AMD:     return orderAMD(size, colPtr, rowIdx, P);
#else
        case Ordering::AMD:     return false;   // vendored, and private/ is not present
#endif
        case Ordering::AMD1:    return orderAMD1(size, colPtr, rowIdx, P);
        case Ordering::AMD2:    return orderAMD2(size, colPtr, rowIdx, P);
        case Ordering::AMD1B:   return orderAMD1B(size, colPtr, rowIdx, P);
        case Ordering::AMD2B:   return orderAMD2B(size, colPtr, rowIdx, P);
    }
    return false;   // unreachable: every enumerator is named above, which is what -Wall checks
}

bool OrderEngine::orderNatural(std::size_t size, Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    P.setIdentity();
    return true;
}

#ifdef OBLIO_VENDORED_ORDERINGS
bool OrderEngine::orderMMD(std::size_t size,
                           const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx,
                           Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    if (size == 0) return true;

    // A is stored full-symmetric; the vendored MMD wants the off-diagonal structure only, so
    // strip the diagonal (no expansion needed, A already holds both triangles). This is MMD's
    // requirement alone: AMD symmetrizes and drops the diagonal itself, and MMD1 builds its own
    // adjacency lists skipping i == j. Columns are indices, so aj is an int32_t and the
    // comparison against rowIdx[cp] needs no cast; cp is a position into A's arrays.
    std::vector<std::size_t> colPtrOff(size + 1, 0);
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(size); ++aj)
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) colPtrOff[aj + 1]++;
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(size); ++aj)
        colPtrOff[aj + 1] += colPtrOff[aj];
    std::vector<std::int32_t> rowIdxOff(colPtrOff[size]);
    std::vector<std::size_t> cur(colPtrOff.begin(), colPtrOff.end());
    for (std::int32_t aj = 0; aj < static_cast<std::int32_t>(size); ++aj)
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] != aj) rowIdxOff[cur[aj]++] = rowIdx[cp];

    // Crossing into the vendored C API, which is int-based throughout. These casts are not the
    // index/position crossings of our own type rules; they are the boundary of a foreign
    // interface, and the arrays below exist only to feed it.
    const int N   = static_cast<int>(size);
    const int nnz = static_cast<int>(rowIdxOff.size());

    std::vector<int> cp(N + 1), ri(nnz);
    for (int j = 0; j <= N; ++j) cp[j] = static_cast<int>(colPtrOff[j]);
    for (int k = 0; k < nnz; ++k) ri[k] = static_cast<int>(rowIdxOff[k]);

    std::vector<int> perm(N), invp(N);
    mmd_order(N, cp.data(), ri.data(), perm.data(), invp.data());

    for (int j = 0; j < N; ++j) {
        P.mOldToNew[j] = static_cast<std::int32_t>(invp[j]);
        P.mNewToOld[j] = static_cast<std::int32_t>(perm[j]);
    }
    return true;
}
#endif

// Ours, over the shared quotient graph. The driver returns an elimination order, exactly as the
// vendored AMD returns its P, and the two maps are written from it here: the friend grant on
// Permutation belongs to the engine that fills it, so the algorithms below stay free functions
// with no access of their own and can be called without a Permutation in hand.
bool OrderEngine::orderMMD1(std::size_t size,
                            const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    if (size == 0) return true;

    const std::vector<std::int32_t> order = orderMmd1(colPtr, rowIdx);
    if (order.size() != size) return false;

    for (std::size_t k = 0; k < size; ++k) {
        P.mNewToOld[k]        = order[k];
        P.mOldToNew[order[k]] = static_cast<std::int32_t>(k);
    }
    return true;
}

// Ours, the completed MMD, written into the maps the same way.
bool OrderEngine::orderMMD2(std::size_t size,
                            const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    if (size == 0) return true;

    const std::vector<std::int32_t> order = orderMmd2(colPtr, rowIdx);
    if (order.size() != size) return false;

    for (std::size_t k = 0; k < size; ++k) {
        P.mNewToOld[k]        = order[k];
        P.mOldToNew[order[k]] = static_cast<std::int32_t>(k);
    }
    return true;
}

// A is full-symmetric; AMD ignores the diagonal and symmetrizes internally, so its structure
// can be passed straight through.
#ifdef OBLIO_VENDORED_ORDERINGS
bool OrderEngine::orderAMD(std::size_t size,
                           const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx,
                           Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    if (size == 0) return true;
    const int N  = static_cast<int>(size);
    const int nz = static_cast<int>(colPtr[size]);

    std::vector<int> Ap(N + 1), Ai(nz);
    for (int j = 0; j <= N; ++j) Ap[j] = static_cast<int>(colPtr[j]);
    for (int k = 0; k < nz; ++k) Ai[k] = static_cast<int>(rowIdx[k]);

    std::vector<int> perm(N);
    const int status = amd_order(N, Ap.data(), Ai.data(), perm.data(), nullptr, nullptr);
    if (status < 0) return false;

    for (int k = 0; k < N; ++k) {
        P.mNewToOld[k]       = static_cast<std::int32_t>(perm[k]);
        P.mOldToNew[perm[k]] = static_cast<std::int32_t>(k);
    }
    return true;
}
#endif

// Ours, over the same quotient graph as MMD1, and written into the maps the same way.
bool OrderEngine::orderAMD1(std::size_t size,
                            const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    if (size == 0) return true;

    const std::vector<std::int32_t> order = orderAmd1(colPtr, rowIdx);
    if (order.size() != size) return false;

    for (std::size_t k = 0; k < size; ++k) {
        P.mNewToOld[k]        = order[k];
        P.mOldToNew[order[k]] = static_cast<std::int32_t>(k);
    }
    return true;
}

bool OrderEngine::orderAMD2(std::size_t size,
                            const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    if (size == 0) return true;

    const std::vector<std::int32_t> order = orderAmd2(colPtr, rowIdx);
    if (order.size() != size) return false;

    for (std::size_t k = 0; k < size; ++k) {
        P.mNewToOld[k]        = order[k];
        P.mOldToNew[order[k]] = static_cast<std::int32_t>(k);
    }
    return true;
}

bool OrderEngine::orderAMD1B(std::size_t size,
                            const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    if (size == 0) return true;

    const std::vector<std::int32_t> order = orderAmd1B(colPtr, rowIdx);
    if (order.size() != size) return false;

    for (std::size_t k = 0; k < size; ++k) {
        P.mNewToOld[k]        = order[k];
        P.mOldToNew[order[k]] = static_cast<std::int32_t>(k);
    }
    return true;
}

bool OrderEngine::orderAMD2B(std::size_t size,
                            const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            Permutation& P) const {
    P.mOldToNew.assign(size, 0);
    P.mNewToOld.assign(size, 0);
    if (size == 0) return true;

    const std::vector<std::int32_t> order = orderAmd2B(colPtr, rowIdx);
    if (order.size() != size) return false;

    for (std::size_t k = 0; k < size; ++k) {
        P.mNewToOld[k]        = order[k];
        P.mOldToNew[order[k]] = static_cast<std::int32_t>(k);
    }
    return true;
}

template bool OrderEngine::compute(const SparseMatrix<double>&, Permutation&) const;
template bool OrderEngine::compute(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
