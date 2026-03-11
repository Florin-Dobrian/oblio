#pragma once
#include "oblio/Types.h"
#include "oblio/BlasLapack.h"
#include "oblio/Matrix.h"
#include "oblio/Permutation.h"
#include "oblio/Symbolic.h"
#include "oblio/Factors.h"
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <cassert>

namespace Oblio {

template<class Val>
class FactorEngine {
public:
    FactorEngine()
      : mAlg(FactorAlg::eMultifrontal), mType(FactorType::eCholesky),
        mPert(1e-14), mAlpha(0.1),
        mNPert(0), mNSwap(0), mNDelay(0), mN1x1(0), mN2x2(0), mRank(0) {}

    void setAlg (FactorAlg  a) { mAlg  = a; }
    void setType(FactorType t) { mType = t; }
    void setPert(RVal p)       { mPert = p; }
    void setAlpha(RVal a)      { mAlpha= a; }

    // Statistics after factorization.
    Size nPert()  const { return mNPert;  }
    Size nSwap()  const { return mNSwap;  }
    Size nDelay() const { return mNDelay; }
    Size n1x1()   const { return mN1x1;   }
    Size n2x2()   const { return mN2x2;   }
    Size rank()   const { return mRank;   }

    Err factor(const Matrix<Val>& A, const Permutation& p,
               const Symbolic& s, Factors<Val>& lu) const;

private:
    FactorAlg  mAlg;
    FactorType mType;
    RVal       mPert, mAlpha;
    mutable Size mNPert,mNSwap,mNDelay,mN1x1,mN2x2,mRank;

    // Precomputed permuted lower-triangle matrix (in new ordering).
    // Built once per factor() call; used by asmOrig.
    mutable std::vector<Size> mPColPtr, mPRowIdx;
    mutable std::vector<Val>  mPVal;

    // Build permuted lower-triangle CSC from A and p.
    void buildPermuted(const Matrix<Val>& A, const Permutation& p,
                       bool hermitian=false) const {
        Size n=A.mSize;
        const Size* o2n=p.oldToNewData(); const Size* n2o=p.newToOldData();
        mPColPtr.assign(n+1,0);
        // Count entries per column in new ordering.
        for(Size jOld=0;jOld<n;++jOld){
            for(Size sp=A.mColPtr[jOld];sp<A.mColPtr[jOld+1];++sp){
                Size iOld=A.mRowIdx[sp];
                Size iNew=o2n[iOld], jNew=o2n[jOld];
                Size col=(iNew<=jNew)?iNew:jNew; // min = column in lower tri new order
                mPColPtr[col+1]++;
            }
        }
        for(Size j=0;j<n;++j) mPColPtr[j+1]+=mPColPtr[j];
        Size nnz=mPColPtr[n];
        mPRowIdx.resize(nnz); mPVal.resize(nnz);
        std::vector<Size> cur(mPColPtr.begin(),mPColPtr.end());
        for(Size jOld=0;jOld<n;++jOld){
            for(Size sp=A.mColPtr[jOld];sp<A.mColPtr[jOld+1];++sp){
                Size iOld=A.mRowIdx[sp]; Val v=A.mVal[sp];
                Size iNew=o2n[iOld], jNew=o2n[jOld];
                bool flipped = (iNew < jNew); // lower entry → upper in new ordering
                Size col=flipped?iNew:jNew;
                Size row=flipped?jNew:iNew;
                // For Hermitian matrices (Cholesky): A[col,row] = conj(A[row,col]).
                // Conjugate when a lower-tri entry is flipped to upper in new ordering.
                // For symmetric (LDL): no conjugation — A[col,row] = A[row,col].
                Size pos=cur[col]++;
                mPRowIdx[pos]=row;
                mPVal[pos]=(hermitian&&flipped)?conjv(v):v;
            }
        }
    }

    FactorEngine(const FactorEngine&)=delete;
    FactorEngine& operator=(const FactorEngine&)=delete;

    // ---- shared helpers ----
    void bldG2L(const Factors<Val>& lu, Size jj, std::vector<Size>& g2l) const;
    void clrG2L(const Factors<Val>& lu, Size jj, std::vector<Size>& g2l) const;
    bool asmOrig(const Matrix<Val>& A, const Permutation& p,
                 Factors<Val>& lu, Size jj, Size nDel,
                 const std::vector<Size>& g2l) const;
    void asmTemp(const std::vector<Size>& tIdx, const std::vector<Val>& tVal,
                 Size tH, Size tW,
                 Factors<Val>& lu, Size kk,
                 const std::vector<Size>& g2l) const;

    // ---- local factor kernels ----
    Err factCC (Factors<Val>& lu, Size jj) const;
    Err factSLDL(Factors<Val>& lu, Size jj) const;
    Err factDLDL(Factors<Val>& lu, Size jj, std::vector<Size>& g2l) const;

    // ---- Schur complement (into Temporary arrays) ----
    void updCC (const Factors<Val>& lu, Size jj, Size off,
                std::vector<Val>& tv, Size tH, Size tW) const;
    void updSLDL(const Factors<Val>& lu, Size jj, Size off,
                 std::vector<Val>& tv, Size tH, Size tW) const;
    void updDLDL(const Factors<Val>& lu, Size jj, Size off,
                 std::vector<Val>& tv, Size tH, Size tW) const;

    // ---- traversal loops ----
    Err lftLook(const Matrix<Val>&, const Permutation&, const Symbolic&, Factors<Val>&) const;
    Err rgtLook(const Matrix<Val>&, const Permutation&, const Symbolic&, Factors<Val>&) const;
    Err multiFnt(const Matrix<Val>&, const Permutation&, const Symbolic&, Factors<Val>&) const;

    // ---- multifrontal specific ----
    void mfAsmUpdt(const UpdateStack<Val>& u, Factors<Val>& lu,
                   Size jj, Size kk, const std::vector<Size>& g2l) const;
    void mfUpdCC (const Factors<Val>& lu, UpdateStack<Val>& u, Size jj) const;
    void mfUpdSLDL(const Factors<Val>& lu, UpdateStack<Val>& u, Size jj) const;
    void mfUpdDLDL(const Factors<Val>& lu, UpdateStack<Val>& u, Size jj) const;

    static RVal av(Val v){ return std::abs(v); }
};


// Explicit instantiation declarations — definitions in FactorEngine.cc
extern template class FactorEngine<double>;
extern template class FactorEngine<std::complex<double>>;

} // namespace Oblio
