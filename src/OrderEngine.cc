#include "oblio/OrderEngine.h"
#include <vector>

// MMD ordering — src/mmd/Mmd.cc
extern void mmd_order(int n, const int colPtr[], const int rowIdx[],
                      int perm[], int invp[]);

// AMD ordering — src/amd/Amd.cc
extern "C" int amd_order(int n, const int Ap[], const int Ai[],
                         int P[], double Control[], double Info[]);

namespace Oblio {

Err OrderEngine::rNatural(Size n, Permutation* p) const {
    p->setIdentity(); return Err::eNone;
}

Err OrderEngine::rMMD(Size n,
                       const std::vector<Size>& colPtr,
                       const std::vector<Size>& rowIdx,
                       Permutation* p) const {
    if (!n) { p->setIdentity(); return Err::eNone; }
    int N = (int)n;
    int nnz = (int)rowIdx.size();

    // Convert Size → int for mmd_order.
    std::vector<int> cp(N + 1), ri(nnz);
    for (int j = 0; j <= N; j++) cp[j] = (int)colPtr[j];
    for (int k = 0; k < nnz; k++) ri[k] = (int)rowIdx[k];

    std::vector<int> perm(N), invp(N);
    mmd_order(N, cp.data(), ri.data(), perm.data(), invp.data());

    Size* o2n = p->oldToNewData();
    Size* n2o = p->newToOldData();
    for (int j = 0; j < N; j++) {
        o2n[j] = (Size)invp[j];
        n2o[j] = (Size)perm[j];
    }
    return Err::eNone;
}

Err OrderEngine::rAMD(Size n,
                       const std::vector<Size>& colPtr,
                       const std::vector<Size>& rowIdx,
                       Permutation* p) const {
    if (!n) { p->setIdentity(); return Err::eNone; }
    int N = (int)n;
    int nz = (int)colPtr[n];

    // Convert Size → int32_t for amd_order (0-based, same as our CSC).
    std::vector<int> Ap(N + 1), Ai(nz);
    for (int j = 0; j <= N; j++) Ap[j] = (int)colPtr[j];
    for (int k = 0; k < nz; k++) Ai[k] = (int)rowIdx[k];

    std::vector<int> P(N);
    int status = amd_order(N, Ap.data(), Ai.data(), P.data(), nullptr, nullptr);
    if (status < 0) return Err::eInvArg;

    Size* o2n = p->oldToNewData();
    Size* n2o = p->newToOldData();
    for (int k = 0; k < N; k++) {
        n2o[k] = (Size)P[k];
        o2n[P[k]] = (Size)k;
    }
    return Err::eNone;
}

} // namespace Oblio
