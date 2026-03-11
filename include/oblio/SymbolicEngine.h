#pragma once
// SymbolicEngine.h — builds supernodal symbolic factorization.
// Uses Liu's etree algorithm on the PERMUTED full-symmetric CSC.

#include "oblio/Types.h"
#include "oblio/Matrix.h"
#include "oblio/Permutation.h"
#include "oblio/Symbolic.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cassert>

namespace Oblio {

class SymbolicEngine {
public:
    SymbolicEngine() {}
    template<class Val>
    Err compute(const Matrix<Val>& A, const Permutation& p, Symbolic* s) const;
private:
    // Build full-symmetric CSC in new ordering from lower-tri CSC in old ordering.
    template<class Val>
    void buildFullSymCSC(const Matrix<Val>& A, const Permutation& p,
                         std::vector<Size>& colPtr,
                         std::vector<Size>& rowIdx,
                         std::vector<Size>& /*val unused for symbolic*/) const;

    void computeEtree(Size n,
                      const std::vector<Size>& colPtr,
                      const std::vector<Size>& rowIdx,
                      std::vector<Size>& prnt) const;

    void computeColCounts(Size n,
                          const std::vector<Size>& colPtr,
                          const std::vector<Size>& rowIdx,
                          const std::vector<Size>& prnt,
                          std::vector<Size>& colCount) const;

    void amalgamate(Size n, const std::vector<Size>& prnt,
                    const std::vector<Size>& colCount,
                    std::vector<Size>& snodeOf) const;

    void buildChildLinks(Size n, const std::vector<Size>& prnt,
                         std::vector<Size>& fstChld,
                         std::vector<Size>& nxtSblg) const;

    // Compute row pattern of column j in L using the etree and original structure.
    void computeRowPattern(Size n, Size j,
                           const std::vector<Size>& colPtr,
                           const std::vector<Size>& rowIdx,
                           const std::vector<Size>& prnt,
                           std::vector<Size>& rowPattern,
                           std::vector<Size>& mark) const;
};

template<class Val>
void SymbolicEngine::buildFullSymCSC(const Matrix<Val>& A, const Permutation& p,
                                      std::vector<Size>& cp,
                                      std::vector<Size>& ri,
                                      std::vector<Size>& /*unused*/) const {
    Size n=A.mSize;
    const Size* o2n=p.oldToNewData();
    const Size* n2o=p.newToOldData();
    // Count off-diagonal entries per column in new ordering (both directions).
    cp.assign(n+1,0);
    for(Size jOld=0;jOld<n;++jOld){
        for(Size sp=A.mColPtr[jOld];sp<A.mColPtr[jOld+1];++sp){
            Size iOld=A.mRowIdx[sp];
            if(iOld==jOld) continue; // skip diagonal
            Size iNew=o2n[iOld], jNew=o2n[jOld];
            cp[jNew+1]++; // entry in column jNew
            cp[iNew+1]++; // symmetric entry in column iNew
        }
    }
    for(Size j=0;j<n;++j) cp[j+1]+=cp[j];
    ri.resize(cp[n]);
    std::vector<Size> cur(cp.begin(),cp.end());
    for(Size jOld=0;jOld<n;++jOld){
        for(Size sp=A.mColPtr[jOld];sp<A.mColPtr[jOld+1];++sp){
            Size iOld=A.mRowIdx[sp];
            if(iOld==jOld) continue;
            Size iNew=o2n[iOld], jNew=o2n[jOld];
            ri[cur[jNew]++]=iNew;
            ri[cur[iNew]++]=jNew;
        }
    }
}

inline void SymbolicEngine::computeEtree(Size n,
                                          const std::vector<Size>& colPtr,
                                          const std::vector<Size>& rowIdx,
                                          std::vector<Size>& prnt) const {
    // Liu's column elimination tree — matches oblio 0.9 EliminationForestEngine.
    // For each column k, process rows j < k.  Find root h of subtree containing
    // j; if h is not yet linked to k, set parent(h) = k.
    prnt.assign(n, cNullIdx);
    std::vector<Size> anc(n, cNullIdx);
    for (Size k = 0; k < n; ++k) {
        for (Size sp = colPtr[k]; sp < colPtr[k + 1]; ++sp) {
            Size j = rowIdx[sp];
            if (j >= k) continue;          // only j < k (upper triangle)
            Size h = j;
            while (anc[h] != cNullIdx && anc[h] != k) {
                Size t = anc[h]; anc[h] = k; h = t;
            }
            if (anc[h] == cNullIdx) { anc[h] = k; prnt[h] = k; }
        }
    }
}

// Helper to safely read prnt (avoid shadowing)

inline void SymbolicEngine::computeColCounts(Size n,
                                              const std::vector<Size>& colPtr,
                                              const std::vector<Size>& rowIdx,
                                              const std::vector<Size>& prnt,
                                              std::vector<Size>& colCount) const {
    // Matches oblio 0.9 EliminationForestEngine column size computation.
    // For each column k, walk from every j < k (that appears in A's column k)
    // up to k via the etree parent links, counting unmarked nodes on the path.
    colCount.assign(n, 0);
    std::vector<Size> color(n, cNullIdx);
    for (Size k = 0; k < n; ++k) {
        colCount[k]++;           // diagonal
        color[k] = k;
        for (Size sp = colPtr[k]; sp < colPtr[k + 1]; ++sp) {
            Size j = rowIdx[sp];
            if (j >= k) continue;   // only j < k
            Size h = j;
            while (color[h] != k) {
                colCount[h]++;
                color[h] = k;
                h = prnt[h];
            }
        }
    }
}

inline void SymbolicEngine::amalgamate(Size n,
                                        const std::vector<Size>& prnt,
                                        const std::vector<Size>& colCount,
                                        std::vector<Size>& snodeOf) const {
    std::vector<Size> nch(n,0);
    for(Size j=0;j<n;++j) if(prnt[j]!=cNullIdx) nch[prnt[j]]++;
    snodeOf.resize(n);
    for(Size j=0;j<n;++j) snodeOf[j]=j;
    for(Size j=0;j<n;++j){
        Size pj=prnt[j]; if(pj==cNullIdx) continue;
        if(nch[pj]==1 && colCount[j]==colCount[pj]+1)
            snodeOf[pj]=snodeOf[j];
    }
    // iterative path compression
    for(Size j=0;j<n;++j){
        Size r=j; while(snodeOf[r]!=r) r=snodeOf[r];
        Size k=j; while(snodeOf[k]!=r){Size nk=snodeOf[k];snodeOf[k]=r;k=nk;}
    }
}

inline void SymbolicEngine::buildChildLinks(Size n,
                                             const std::vector<Size>& prnt,
                                             std::vector<Size>& fstChld,
                                             std::vector<Size>& nxtSblg) const {
    fstChld.assign(n,cNullIdx); nxtSblg.assign(n,cNullIdx);
    for(Size j=n;j-->0;){
        if(prnt[j]==cNullIdx) continue;
        nxtSblg[j]=fstChld[prnt[j]]; fstChld[prnt[j]]=j;
    }
}

inline void SymbolicEngine::computeRowPattern(Size n, Size j,
                                               const std::vector<Size>& colPtr,
                                               const std::vector<Size>& rowIdx,
                                               const std::vector<Size>& prnt,
                                               std::vector<Size>& rows,
                                               std::vector<Size>& mark) const {
    // Row pattern of column j in L = {j} union {rows of A in col j with i>j}
    //   union (row patterns of snodes updating j via etree)
    rows.clear(); rows.push_back(j); mark[j]=j;
    for(Size sp=colPtr[j];sp<colPtr[j+1];++sp){
        Size i=rowIdx[sp]; if(i<=j) continue;
        // climb etree from i, collecting rows
        Size r=i;
        while(mark[r]!=j){
            mark[r]=j; rows.push_back(r);
            if(prnt[r]==cNullIdx) break;
            r=prnt[r]; if(r<=j) break;
        }
    }
    std::sort(rows.begin(),rows.end());
}

template<class Val>
Err SymbolicEngine::compute(const Matrix<Val>& A, const Permutation& p, Symbolic* s) const {
    Size n=A.mSize;
    if(!A.isValid())          return Err::eInvMat;
    if(!p.isValid())          return Err::eInvArg;
    if(A.mSize!=p.getSize()) return Err::eDimMismatch;

    // 1. Build full-symmetric CSC of A in new ordering (off-diagonal only).
    std::vector<Size> colPtr, rowIdx, dummy;
    buildFullSymCSC(A, p, colPtr, rowIdx, dummy);

    // 2. Elimination tree on permuted structure.
    std::vector<Size> prnt;
    computeEtree(n, colPtr, rowIdx, prnt);

    // 3. Column counts.
    std::vector<Size> colCount;
    computeColCounts(n, colPtr, rowIdx, prnt, colCount);

    // 4. Amalgamate into fundamental supernodes.
    std::vector<Size> snodeOf;
    amalgamate(n, prnt, colCount, snodeOf);

    // 5. Assign supernode IDs.
    std::vector<Size> repToId(n,cNullIdx);
    Size nSnd=0;
    for(Size j=0;j<n;++j) if(snodeOf[j]==j) repToId[j]=nSnd++;
    std::vector<Size> idxToSnode(n);
    for(Size j=0;j<n;++j) idxToSnode[j]=repToId[snodeOf[j]];

    // 6. Supernode parent links.
    std::vector<Size> sndPrnt(nSnd,cNullIdx);
    for(Size j=0;j<n;++j){
        if(prnt[j]==cNullIdx) continue;
        Size sj=idxToSnode[j], sp2=idxToSnode[prnt[j]];
        if(sj!=sp2) sndPrnt[sj]=sp2;
    }
    std::vector<Size> fstChld, nxtSblg;
    buildChildLinks(nSnd, sndPrnt, fstChld, nxtSblg);

    // 7. Build index sets bottom-up (matching oblio 0.9 SymbolicEngine::run_).
    //    For each snode kk:
    //      a) Add all rows >= first_front_col from the original matrix columns
    //      b) Add all update rows inherited from children
    //    Then sort indices; front entries come first (sorted), then update rows.
    std::vector<Size> frntSz(nSnd,0), updtSz(nSnd,0);
    for(Size j=0;j<n;++j) frntSz[idxToSnode[j]]++;

    // Build front-index arrays per snode (sorted column indices).
    std::vector<std::vector<Size>> frontCols(nSnd);
    for(Size j=0;j<n;++j) frontCols[idxToSnode[j]].push_back(j);

    std::vector<std::vector<Size>> idxVec(nSnd);
    std::vector<Size> color(n, cNullIdx);

    for(Size kk=0;kk<nSnd;++kk){
        auto& idx=idxVec[kk];
        // a) Original matrix entries for each front column of kk.
        for(Size k : frontCols[kk]){
            for(Size sp=colPtr[k];sp<colPtr[k+1];++sp){
                Size i=rowIdx[sp];
                if(i<k) continue;          // only rows >= k
                if(color[i]==kk) continue;  // already added
                color[i]=kk; idx.push_back(i);
            }
            // Also add k itself (diagonal) if not yet added.
            if(color[k]!=kk){ color[k]=kk; idx.push_back(k); }
        }
        // b) Inherit update rows from children.
        for(Size jj=fstChld[kk];jj!=cNullIdx;jj=nxtSblg[jj]){
            for(Size r : idxVec[jj]){
                if(idxToSnode[r]==jj) continue;  // skip jj's own front columns
                if(color[r]==kk) continue;        // already added
                color[r]=kk; idx.push_back(r);
            }
        }
        // Sort: front columns (< first update row) come first.
        std::sort(idx.begin(),idx.end());
        // Count update size.
        Size f=0; for(Size r:idx) if(idxToSnode[r]==kk) f++;
        updtSz[kk]=idx.size()-f;
        assert(frntSz[kk]==f);
    }

    // 8. Populate Symbolic.
    s->mSize=n; s->mNumSnodes=nSnd;
    s->mPrntVec=sndPrnt; s->mFstChldVec=fstChld; s->mNxtSblgVec=nxtSblg;
    s->mIdxToSnodeVec=idxToSnode;
    s->mFrntSzVec=frntSz; s->mUpdtSzVec=updtSz;
    s->mIdxVecVec=idxVec;

    s->mNumTrees=0; s->mFstRoot=cNullIdx; s->mLstRoot=cNullIdx;
    for(Size jj=0;jj<nSnd;++jj)
        if(sndPrnt[jj]==cNullIdx){s->mNumTrees++;s->mLstRoot=jj;if(s->mFstRoot==cNullIdx)s->mFstRoot=jj;}

    Size ti=0,tv=0;
    for(Size jj=0;jj<nSnd;++jj){ti+=idxVec[jj].size();tv+=frntSz[jj]*idxVec[jj].size();}
    s->mNumIdxs=ti; s->mNumVals=tv;
    return Err::eNone;
}

} // namespace Oblio
