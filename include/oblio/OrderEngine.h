#pragma once
#include "oblio/Types.h"
#include "oblio/Matrix.h"
#include "oblio/Permutation.h"
#include <vector>

namespace Oblio {

// ============================================================================
// OrderEngine — computes fill-reducing permutation.
// Supports eNatural, eMMD (Multiple Minimum Degree), eAMD (Approximate MD).
// MMD: src/mmd/Mmd.cc (Liu/Sparspak via oblio 0.9)
// AMD: src/amd/Amd.cc (SuiteSparse AMD 3.3.4, Davis/Amestoy/Duff, BSD-3-clause)
// ============================================================================

class OrderEngine {
public:
    OrderEngine() : mAlg(OrderAlg::eMMD) {}
    void setAlg(OrderAlg a) { mAlg=a; }

    template<class Val>
    Err order(const Matrix<Val>& A, Permutation* p) const;

private:
    OrderAlg mAlg;
    Err rNatural(Size n, Permutation* p) const;
    Err rMMD(Size n, const std::vector<Size>& colPtr,
             const std::vector<Size>& rowIdx, Permutation* p) const;
    Err rAMD(Size n, const std::vector<Size>& colPtr,
             const std::vector<Size>& rowIdx, Permutation* p) const;
};

template<class Val>
Err OrderEngine::order(const Matrix<Val>& A, Permutation* p) const {
    Size n = A.mSize;
    p->resize(n);
    if (mAlg == OrderAlg::eNatural) return rNatural(n, p);

    // AMD can take the raw lower-triangle CSC directly (it computes A+A'
    // internally).  Pass it without modification.
    if (mAlg == OrderAlg::eAMD)
        return rAMD(n, A.mColPtr, A.mRowIdx, p);

    // MMD needs the full-symmetric off-diagonal-only CSC.
    std::vector<Size> colPtr(n+1, 0);
    for (Size j=0;j<n;++j)
        for (Size sp=A.mColPtr[j];sp<A.mColPtr[j+1];++sp) {
            Size i=A.mRowIdx[sp];
            if (i!=j) { colPtr[j+1]++; colPtr[i+1]++; }
        }
    for (Size j=0;j<n;++j) colPtr[j+1]+=colPtr[j];
    std::vector<Size> rowIdx(colPtr[n]);
    std::vector<Size> cur(colPtr.begin(),colPtr.end());
    for (Size j=0;j<n;++j)
        for (Size sp=A.mColPtr[j];sp<A.mColPtr[j+1];++sp) {
            Size i=A.mRowIdx[sp];
            if (i!=j) { rowIdx[cur[j]++]=i; rowIdx[cur[i]++]=j; }
        }
    return rMMD(n, colPtr, rowIdx, p);
}

} // namespace Oblio

// Implementation in OrderEngine.cc
