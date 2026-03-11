#include "oblio/FactorEngine.h"
#include <complex>

namespace Oblio {

// ============================================================================
// Implementation
// ============================================================================

template<class Val>
Err FactorEngine<Val>::factor(const Matrix<Val>& A, const Permutation& p,
                               const Symbolic& s, Factors<Val>& lu) const {
    mNPert=mNSwap=mNDelay=mN1x1=mN2x2=0; mRank=A.mSize;
    if (!A.isValid())          return Err::eInvMat;
    if (!p.isValid())          return Err::eInvArg;
    if (A.mSize!=p.getSize()) return Err::eDimMismatch;
    if (A.mSize!=s.mSize)     return Err::eDimMismatch;
    switch(mAlg){
        case FactorAlg::eLeftLooking:  return lftLook(A,p,s,lu);
        case FactorAlg::eRightLooking: return rgtLook(A,p,s,lu);
        case FactorAlg::eMultifrontal: return multiFnt(A,p,s,lu);
    }
    return Err::eInvArg;
}

// ---- bldG2L / clrG2L -------------------------------------------------------
template<class Val>
void FactorEngine<Val>::bldG2L(const Factors<Val>& lu, Size jj,
                                std::vector<Size>& g2l) const {
    const Size* ix=lu.idx(jj); Size nI=lu.nIdx(jj);
    for(Size k=0;k<nI;++k) g2l[ix[k]]=k;
}
template<class Val>
void FactorEngine<Val>::clrG2L(const Factors<Val>& lu, Size jj,
                                std::vector<Size>& g2l) const {
    const Size* ix=lu.idx(jj); Size nI=lu.nIdx(jj);
    for(Size k=0;k<nI;++k) g2l[ix[k]]=cNullIdx;
}

// ---- asmOrig ---------------------------------------------------------------
// Scatter entries of the PERMUTED matrix AP into front matrix of snode jj.
// AP is given in CSC (lower triangle in new ordering):
//   apColPtr[j] .. apColPtr[j+1]-1 → row indices apRowIdx[], values apVal[].
// apColPtr[j] indexes column j in NEW ordering (already permuted).
template<class Val>
bool FactorEngine<Val>::asmOrig(const Matrix<Val>& A, const Permutation& p,
                                 Factors<Val>& lu, Size jj, Size nDel,
                                 const std::vector<Size>& g2l) const {
    Val* jjVal=lu.val(jj);
    Size jjF=lu.frntSz(jj), nI=lu.nIdx(jj);
    const Size* jjIdx=lu.idx(jj);

    // Use the precomputed permuted lower-triangle colPtr/rowIdx/val arrays.
    // These are stored in mPColPtr, mPRowIdx, mPVal (set in lftLook/rgtLook/multiFnt).
    const std::vector<Size>& cp=mPColPtr;
    const std::vector<Size>& ri=mPRowIdx;
    const std::vector<Val>&  rv=mPVal;

    for(Size j_=nDel;j_<jjF;++j_){
        Size jNew=jjIdx[j_];
        for(Size sp=cp[jNew];sp<cp[jNew+1];++sp){
            Size iNew=ri[sp];
            Size i_=g2l[iNew]; if(i_==cNullIdx) continue;
            jjVal[j_*nI+i_]+=rv[sp];
        }
    }
    return true;
}

// ---- asmTemp ---------------------------------------------------------------
// Scatter Temporary tile tVal (height tH, width tW, col-major) into snode kk.
template<class Val>
void FactorEngine<Val>::asmTemp(const std::vector<Size>& tIdx,
                                 const std::vector<Val>& tVal,
                                 Size tH, Size tW,
                                 Factors<Val>& lu, Size kk,
                                 const std::vector<Size>& g2l) const {
    Val* kkVal=lu.val(kk); Size kkNI=lu.nIdx(kk);
    for(Size j_=0;j_<tW;++j_){
        Size gj=tIdx[j_]; Size dj=g2l[gj]; if(dj==cNullIdx)continue;
        for(Size i_=j_;i_<tH;++i_){
            Size di=g2l[tIdx[i_]]; if(di==cNullIdx)continue;
            kkVal[dj*kkNI+di]+=tVal[j_*tH+i_];
        }
    }
}

// ============================================================================
// Local factor kernels
// ============================================================================

template<class Val>
Err FactorEngine<Val>::factCC(Factors<Val>& lu, Size jj) const {
    Val* v=lu.val(jj); Size f=lu.frntSz(jj), nI=lu.nIdx(jj);
    int info=BT<Val>::potrf('L',(int)f,v,(int)nI);
    if(info>0) return Err::eInvPivot;
    mN1x1+=f;
    if(nI>f) { Val one{1}; BT<Val>::trsm('R','L',conjTrans<Val>(),'N',(int)(nI-f),(int)f,one,v,(int)nI,v+f,(int)nI); }
    return Err::eNone;
}

template<class Val>
Err FactorEngine<Val>::factSLDL(Factors<Val>& lu, Size jj) const {
    Val* v=lu.val(jj); Size f=lu.frntSz(jj), nI=lu.nIdx(jj);
    RVal eps=mPert; if(std::abs(eps)<1e-16)eps=1e-16;
    int np=0;
    int info=oblioPotrfLDL<Val>((int)f,v,(int)nI,eps,&np);
    mNPert+=(Size)np;
    if(info>0) return Err::eInvPivot;
    mN1x1+=f;
    if(nI>f) { Val one{1}; BT<Val>::trsm('R','U','N','N',(int)(nI-f),(int)f,one,v,(int)nI,v+f,(int)nI); }
    return Err::eNone;
}

// Bunch-Kaufman dynamic pivoting — ported from 0.9 factorDynamicLDL_.
template<class Val>
Err FactorEngine<Val>::factDLDL(Factors<Val>& lu, Size jj,
                                  std::vector<Size>& g2l) const {
    Size* ix=lu.idx(jj); Val* v=lu.val(jj);
    Size f=lu.frntSz(jj), nI=lu.nIdx(jj);

    std::list<Size> pending;
    for(Size j_=0;j_<f;++j_) pending.push_back(ix[j_]);

    Size j_cur=0;
    while(!pending.empty()){
        bool found=false;
        Size nTry=pending.size();
        while(nTry-->0){
            Size k1=pending.front(); pending.pop_front();
            Size k1_=g2l[k1];
            Size pk1k1=k1_*nI+k1_;
            // max off-diagonal in column k1_
            RVal max1=0; Size maxI_=cNullIdx;
            for(Size i_=j_cur;i_<nI;++i_){
                if(i_==k1_) continue;
                RVal a=av(v[k1_*nI+i_]);
                if(a>max1){max1=a;maxI_=i_;}
            }
            RVal d1=av(v[pk1k1]);
            // 1x1 pivot?
            if(max1==0 || (d1>0 && d1>=mAlpha*max1)){
                found=true;
                if(j_cur!=k1_){lu.swapCols(jj,j_cur,k1_,g2l);++mNSwap;}
                Val diag=v[j_cur*nI+j_cur];
                if(av(diag)==0) --mRank;
                if(av(diag)>0){
                    for(Size i_=j_cur+1;i_<nI;++i_) v[j_cur*nI+i_]/=diag;
                }
                for(Size k_=j_cur+1;k_<f;++k_)
                    v[k_*nI+j_cur]=diag*v[j_cur*nI+k_];
                for(Size k_=j_cur+1;k_<f;++k_)
                    for(Size i_=k_;i_<nI;++i_)
                        v[k_*nI+i_]-=v[j_cur*nI+i_]*v[k_*nI+j_cur];
                lu.appendPivot(jj,1); ++mN1x1; ++j_cur; break;
            }
            // try 2x2
            if(!pending.empty() && maxI_!=cNullIdx){
                Size k2_=maxI_; Size k2=ix[k2_];
                Size pk2k2=k2_*nI+k2_;
                RVal max2=0;
                for(Size i_=j_cur;i_<nI;++i_){
                    if(i_==k2_) continue;
                    RVal a=av(v[k2_*nI+i_]); if(a>max2) max2=a;
                }
                Val d11=v[k1_*nI+k1_], d22=v[pk2k2];
                Val d12=(k2_>k1_)?v[k1_*nI+k2_]:v[k2_*nI+k1_];
                Val det=d11*d22-d12*d12;
                RVal mm=std::max(av(d22)*max1+max1*max2, av(d11)*max2+max1*max1);
                bool acc=(max1==max2&&max1!=0)||(av(det)>0&&av(det)>=mAlpha*mm);
                if(acc){
                    found=true; pending.remove(k2);
                    Size j1_=j_cur, j2_=j_cur+1;
                    if(k1_<k2_){
                        if(j1_!=k1_){lu.swapCols(jj,j1_,k1_,g2l);++mNSwap;}
                        if(j2_!=k2_){lu.swapCols(jj,j2_,k2_,g2l);++mNSwap;}
                    }else if(j1_==k2_&&j2_==k1_){
                        lu.swapCols(jj,j1_,j2_,g2l);++mNSwap;
                    }else{
                        if(j2_!=k2_){lu.swapCols(jj,j2_,k2_,g2l);++mNSwap;}
                        lu.swapCols(jj,j1_,k1_,g2l);++mNSwap;
                    }
                    // reread after swaps
                    d11=v[j1_*nI+j1_]; d22=v[j2_*nI+j2_];
                    d12=v[j1_*nI+j2_]; det=d11*d22-d12*d12;
                    v[j2_*nI+j1_]=d12; // store in upper for solve
                    for(Size i_=j1_+2;i_<nI;++i_){
                        Val t1=v[j1_*nI+i_], t2=v[j2_*nI+i_];
                        v[j1_*nI+i_]=(t1*d22-t2*d12)/det;
                        v[j2_*nI+i_]=(t2*d11-t1*d12)/det;
                    }
                    for(Size k_=j1_+2;k_<f;++k_){
                        Val l1=v[j1_*nI+k_], l2=v[j2_*nI+k_];
                        v[k_*nI+j1_]=d11*l1+d12*l2;
                        v[k_*nI+j2_]=d12*l1+d22*l2;
                    }
                    for(Size k_=j1_+2;k_<f;++k_)
                        for(Size i_=k_;i_<nI;++i_)
                            v[k_*nI+i_]-=v[j1_*nI+i_]*v[k_*nI+j1_];
                    for(Size k_=j2_+1;k_<f;++k_)
                        for(Size i_=k_;i_<nI;++i_)
                            v[k_*nI+i_]-=v[j2_*nI+i_]*v[k_*nI+j2_];
                    lu.appendPivot(jj,2); lu.appendPivot(jj,3);
                    ++mN2x2; j_cur+=2; break;
                }
            }
            pending.push_back(k1);
        }
        if(!found) break;
    }
    Size nDel=pending.size(); mNDelay+=nDel;
    lu.setNDelayed(jj,nDel); lu.shrinkFront(jj,nDel);
    return Err::eNone;
}

// ============================================================================
// Schur complement (into tVal, column-major, tH x tW)
// off = first row index into the snode's column (>= frontSz)
// ============================================================================

template<class Val>
void FactorEngine<Val>::updCC(const Factors<Val>& lu, Size jj, Size off,
                               std::vector<Val>& tv, Size tH, Size tW) const {
    const Val* v=lu.val(jj); Size f=lu.frntSz(jj), nI=lu.nIdx(jj);
    // t -= L21*L21^H  (lower tri for Hermitian/symmetric part, full for rectangular)
    Val m1{-1},p1{1};
    BT<Val>::rankUpdate('L',(int)tW,(int)f,m1,v+off,(int)nI,p1,tv.data(),(int)tH);
    if(tW<tH)
        BT<Val>::gemm('N',conjTrans<Val>(),(int)(tH-tW),(int)tW,(int)f,
                      m1,v+off+tW,(int)nI, v+off,(int)nI,
                      p1,tv.data()+tW,(int)tH);
}

template<class Val>
void FactorEngine<Val>::updSLDL(const Factors<Val>& lu, Size jj, Size off,
                                 std::vector<Val>& tv, Size tH, Size tW) const {
    const Val* v=lu.val(jj); Size f=lu.frntSz(jj), nI=lu.nIdx(jj);
    std::vector<Val> u2(f*tW);
    oblioComputeU<Val>((int)tW,(int)f, v+off,(int)nI, u2.data(),(int)f, v,(int)nI);
    oblioGemm<Val>((int)tW,(int)f, v+off,(int)nI, u2.data(),(int)f, tv.data(),(int)tH);
    if(tW<tH){ Val m1{-1},p1{1};
        BT<Val>::gemm('N','N',(int)(tH-tW),(int)tW,(int)f,
                      m1,v+off+tW,(int)nI, u2.data(),(int)f,
                      p1,tv.data()+tW,(int)tH); }
}

template<class Val>
void FactorEngine<Val>::updDLDL(const Factors<Val>& lu, Size jj, Size off,
                                 std::vector<Val>& tv, Size tH, Size tW) const {
    const Val* v=lu.val(jj); const int* pt=lu.piv(jj);
    const Size* ix=lu.idx(jj);
    Size f=lu.frntSz(jj), nI=lu.nIdx(jj);
    if(!f) return;
    std::vector<Val> u2(f*tW);
    Size j_=0,pjj=0;
    while(j_<f){
        if(pt[j_]==1){
            Val d=v[pjj];
            for(Size k_=off,uk=j_;k_<off+tW;++k_,uk+=f)
                u2[uk]=d*v[pjj+k_-j_];
            j_++; pjj+=nI+1;
        }else{
            Val d11=v[pjj],d22=v[pjj+nI+1],d12=v[pjj+nI];
            for(Size k_=off;k_<off+tW;++k_){
                Size uk1=j_+(k_-off)*f, uk2=(j_+1)+(k_-off)*f;
                Val l1=v[pjj+k_-j_], l2=v[pjj+nI+1+k_-(j_+1)];
                u2[uk1]=d11*l1+d12*l2; u2[uk2]=d12*l1+d22*l2;
            }
            j_+=2; pjj+=2*nI+2;
        }
    }
    oblioGemm<Val>((int)tW,(int)f, v+off,(int)nI, u2.data(),(int)f, tv.data(),(int)tH);
    if(tW<tH){ Val m1{-1},p1{1};
        BT<Val>::gemm('N','N',(int)(tH-tW),(int)tW,(int)f,
                      m1,v+off+tW,(int)nI, u2.data(),(int)f,
                      p1,tv.data()+tW,(int)tH); }
}

// ============================================================================
// Left-looking traversal
// ============================================================================
template<class Val>
Err FactorEngine<Val>::lftLook(const Matrix<Val>& A, const Permutation& p,
                                const Symbolic& s, Factors<Val>& lu) const {
    lu.allocate(s); lu.zero();
    Size n=A.mSize, nSnd=s.mNumSnodes;
    lu.setFactorType(mType);
    buildPermuted(A,p, mType==FactorType::eCholesky);
    std::vector<Size> g2l(n,cNullIdx);
    std::vector<Size> pp(nSnd,0);
    std::vector<std::list<Size>> upd(nSnd);
    const Size* sMap=s.snodeMapData();
    const Size* sFst=s.fstChldData();
    const Size* sNxt=s.nxtSblgData();
    const Size* sPrnt=s.prntData();

    bool dynLDL=(mType==FactorType::eDynamicLDL);

    // For static types: pre-scatter A.
    if(!dynLDL){
        for(Size jj=0;jj<nSnd;++jj){
            bldG2L(lu,jj,g2l); asmOrig(A,p,lu,jj,0,g2l); clrG2L(lu,jj,g2l);
        }
    }

    for(Size kk=0;kk<nSnd;++kk){
        Size nDel=0;
        if(dynLDL){
            for(Size jj=sFst[kk];jj!=cNullIdx;jj=sNxt[jj]) nDel+=lu.nDelayed(jj);
            if(nDel){ lu.extendFront(kk,nDel,s.mFstChldVec,s.mNxtSblgVec); lu.reallocVal(kk); lu.zero(kk); }
            bldG2L(lu,kk,g2l);
            asmOrig(A,p,lu,kk,nDel,g2l);
        } else { bldG2L(lu,kk,g2l); }

        // process all updaters
        while(!upd[kk].empty()){
            Size jj=upd[kk].front(); upd[kk].pop_front();
            const Size* jjIdx=lu.idx(jj);
            Size jjNI=lu.nIdx(jj);
            // for DynamicLDL, assemble delayed cols if jj is child of kk
            if(dynLDL && sPrnt[jj]==kk){ lu.assembleDelayed(jj,kk,g2l); lu.shrinkFront(jj,lu.nDelayed(jj)); }
            Size rr=pp[jj];
            Size tH=jjNI-rr, tW=0, tmpRr=rr;
            do{++tmpRr;++tW;}while(tmpRr<jjNI && sMap[jjIdx[tmpRr]]==kk);
            std::vector<Size> tIdx(tH); for(Size i=0;i<tH;++i) tIdx[i]=jjIdx[rr+i];
            std::vector<Val>  tVal(tH*tW, Val{0});
            switch(mType){
                case FactorType::eCholesky:   updCC  (lu,jj,rr,tVal,tH,tW); break;
                case FactorType::eStaticLDL:  updSLDL(lu,jj,rr,tVal,tH,tW); break;
                case FactorType::eDynamicLDL: updDLDL(lu,jj,rr,tVal,tH,tW); break;
            }
            pp[jj]=tmpRr;
            asmTemp(tIdx,tVal,tH,tW,lu,kk,g2l);
            jjNI=lu.nIdx(jj); // may have changed (shrinkFront)
            if(pp[jj]<jjNI){ Size nx=sMap[lu.idx(jj)[pp[jj]]]; upd[nx].push_back(jj); }
        }

        Err e=Err::eNone;
        switch(mType){
            case FactorType::eCholesky:   e=factCC  (lu,kk); break;
            case FactorType::eStaticLDL:  e=factSLDL(lu,kk); break;
            case FactorType::eDynamicLDL: e=factDLDL(lu,kk,g2l); break;
        }
        if(e!=Err::eNone) return e;

        Size kkF=lu.frntSz(kk), kkND=lu.nDelayed(kk);
        pp[kk]+=(dynLDL?(kkF+kkND):kkF);
        Size kkNI=lu.nIdx(kk);
        if(pp[kk]<kkNI){ Size nx=sMap[lu.idx(kk)[pp[kk]]]; upd[nx].push_back(kk); }
        if(dynLDL&&kkND>0&&sPrnt[kk]==cNullIdx) return Err::eInvPivot;
        clrG2L(lu,kk,g2l);
    }
    lu.symmetrize(); return Err::eNone;
}

// ============================================================================
// Right-looking traversal
// ============================================================================
template<class Val>
Err FactorEngine<Val>::rgtLook(const Matrix<Val>& A, const Permutation& p,
                                const Symbolic& s, Factors<Val>& lu) const {
    lu.allocate(s); lu.zero();
    Size n=A.mSize, nSnd=s.mNumSnodes;
    lu.setFactorType(mType);
    buildPermuted(A,p, mType==FactorType::eCholesky);
    std::vector<Size> g2l(n,cNullIdx);
    std::vector<bool> init(nSnd,false);
    const Size* sMap=s.snodeMapData();

    for(Size jj=0;jj<nSnd;++jj){
        bldG2L(lu,jj,g2l);
        if(!init[jj]){ asmOrig(A,p,lu,jj,0,g2l); init[jj]=true; }
        Err e=Err::eNone;
        switch(mType){
            case FactorType::eCholesky:   e=factCC  (lu,jj); break;
            case FactorType::eStaticLDL:  e=factSLDL(lu,jj); break;
            case FactorType::eDynamicLDL: e=factDLDL(lu,jj,g2l); break;
        }
        if(e!=Err::eNone) return e;

        const Size* jjIdx=lu.idx(jj); Size jjF=lu.frntSz(jj), jjNI=lu.nIdx(jj);
        Size pp=jjF;
        while(pp<jjNI){
            Size kk=sMap[jjIdx[pp]]; Size rr=pp, tW=0;
            do{++rr;++tW;}while(rr<jjNI&&sMap[jjIdx[rr]]==kk);
            Size tH=jjNI-pp;
            std::vector<Size> tIdx(tH); for(Size i=0;i<tH;++i) tIdx[i]=jjIdx[pp+i];
            std::vector<Val>  tVal(tH*tW,Val{0});
            switch(mType){
                case FactorType::eCholesky:   updCC  (lu,jj,pp,tVal,tH,tW); break;
                case FactorType::eStaticLDL:  updSLDL(lu,jj,pp,tVal,tH,tW); break;
                case FactorType::eDynamicLDL: updDLDL(lu,jj,pp,tVal,tH,tW); break;
            }
            pp=rr;
            bldG2L(lu,kk,g2l);
            if(!init[kk]){ asmOrig(A,p,lu,kk,0,g2l); init[kk]=true; }
            asmTemp(tIdx,tVal,tH,tW,lu,kk,g2l);
            clrG2L(lu,kk,g2l);
        }
        clrG2L(lu,jj,g2l);
    }
    lu.symmetrize(); return Err::eNone;
}

// ============================================================================
// Multifrontal traversal
// ============================================================================

// Assemble update matrix u[jj] into front lu[kk].
template<class Val>
void FactorEngine<Val>::mfAsmUpdt(const UpdateStack<Val>& u, Factors<Val>& lu,
                                   Size jj, Size kk,
                                   const std::vector<Size>& g2l) const {
    const Size* uIdx=u.idx(jj); const Val* uVal=u.val(jj);
    Val* kv=lu.val(kk); Size kNI=lu.nIdx(kk);
    Size uSz=u.sz(jj);
    for(Size sj_=0;sj_<uSz;++sj_){
        Size gj=uIdx[sj_]; Size dj=g2l[gj]; if(dj==cNullIdx)continue;
        Size dp0j=dj*kNI;
        for(Size si_=sj_;si_<uSz;++si_){
            Size di=g2l[uIdx[si_]]; if(di==cNullIdx)continue;
            kv[dp0j+di]+=uVal[sj_*uSz+si_];
        }
    }
    // Also assemble delayed columns from jj factor into kk front.
    lu.assembleDelayed(jj,kk,g2l);
}

template<class Val>
void FactorEngine<Val>::mfUpdCC(const Factors<Val>& lu,
                                 UpdateStack<Val>& u, Size jj) const {
    const Val* v=lu.val(jj); Size f=lu.frntSz(jj), nI=lu.nIdx(jj);
    Val* uv=u.val(jj); Size us=u.sz(jj);
    Val m1{-1},p1{1};
    BT<Val>::rankUpdate('L',(int)us,(int)f,m1,v+f,(int)nI,p1,uv,(int)us);
}
template<class Val>
void FactorEngine<Val>::mfUpdSLDL(const Factors<Val>& lu,
                                   UpdateStack<Val>& u, Size jj) const {
    const Val* v=lu.val(jj); Size f=lu.frntSz(jj), nI=lu.nIdx(jj);
    Val* uv=u.val(jj); Size us=u.sz(jj);
    std::vector<Val> u2(f*us);
    oblioComputeU<Val>((int)us,(int)f, v+f,(int)nI, u2.data(),(int)f, v,(int)nI);
    oblioGemm<Val>((int)us,(int)f, v+f,(int)nI, u2.data(),(int)f, uv,(int)us);
}
template<class Val>
void FactorEngine<Val>::mfUpdDLDL(const Factors<Val>& lu,
                                   UpdateStack<Val>& u, Size jj) const {
    const Val* v=lu.val(jj); const int* pt=lu.piv(jj);
    const Size* ix=lu.idx(jj);
    Size f=lu.frntSz(jj), nI=lu.nIdx(jj);
    Val* uv=u.val(jj); Size us=u.sz(jj);
    if(!f||!us) return;
    std::vector<Val> u2(f*us,Val{0});
    Size j_=0,pjj=0;
    while(j_<f){
        if(pt[j_]==1){
            Val d=v[pjj];
            for(Size k_=0;k_<us;++k_) u2[j_+k_*f]=d*v[pjj+f-j_+k_];
            j_++;pjj+=nI+1;
        }else{
            Val d11=v[pjj],d22=v[pjj+nI+1],d12=v[pjj+nI];
            for(Size k_=0;k_<us;++k_){
                Val l1=v[pjj+f-j_+k_], l2=v[pjj+nI+1+f-(j_+1)+k_];
                u2[j_+k_*f]=d11*l1+d12*l2;
                u2[(j_+1)+k_*f]=d12*l1+d22*l2;
            }
            j_+=2;pjj+=2*nI+2;
        }
    }
    oblioGemm<Val>((int)us,(int)f, v+f,(int)nI, u2.data(),(int)f, uv,(int)us);
}

template<class Val>
Err FactorEngine<Val>::multiFnt(const Matrix<Val>& A, const Permutation& p,
                                 const Symbolic& s, Factors<Val>& lu) const {
    lu.allocate(s); lu.zero();
    Size n=A.mSize, nSnd=s.mNumSnodes;
    lu.setFactorType(mType);
    buildPermuted(A,p, mType==FactorType::eCholesky);
    std::vector<Size> g2l(n,cNullIdx);
    const Size* sPrnt=s.prntData();
    const Size* sFst=s.fstChldData();
    const Size* sNxt=s.nxtSblgData();
    UpdateStack<Val> us(nSnd);
    bool dynLDL=(mType==FactorType::eDynamicLDL);

    for(Size kk=0;kk<nSnd;++kk){
        Size nDel=0;
        if(dynLDL){
            for(Size jj=sFst[kk];jj!=cNullIdx;jj=sNxt[jj]) nDel+=lu.nDelayed(jj);
            if(nDel){ lu.extendFront(kk,nDel,s.mFstChldVec,s.mNxtSblgVec); lu.reallocVal(kk); lu.zero(kk); }
        }
        bldG2L(lu,kk,g2l);
        asmOrig(A,p,lu,kk,nDel,g2l);

        // Allocate us[kk] BEFORE child assembly so bypass entries can land there.
        Size kkUS=lu.updtSz(kk);
        if(kkUS>0){
            us.alloc(kk,kkUS);
            us.setIdx(kk,lu.idx(kk)+lu.frntSz(kk),kkUS);
            us.zero(kk);
        }

        Size kF=lu.frntSz(kk);

        // Assemble children's update matrices into lu[kk] and us[kk].
        // Following the 0.9 reference pattern:
        //   - entries whose column maps to kk's pivot region (dj < kF) → lu[kk]
        //   - entries whose column maps to kk's update region (dj >= kF) → us[kk]
        for(Size jj=sFst[kk];jj!=cNullIdx;jj=sNxt[jj]){
            Size jjUS=us.sz(jj);
            const Size* uIdx=us.idx(jj); const Val* uVal=us.val(jj);
            Val* kv=lu.val(kk); Size kNI=lu.nIdx(kk);
            Val* ukv=(kkUS>0)?us.val(kk):nullptr;
            for(Size sj_=0;sj_<jjUS;++sj_){
                Size gj=uIdx[sj_]; Size dj=g2l[gj]; if(dj==cNullIdx)continue;
                if(dj<kF){
                    // Column in kk's pivot front → factor storage.
                    for(Size si_=sj_;si_<jjUS;++si_){
                        Size di=g2l[uIdx[si_]]; if(di==cNullIdx)continue;
                        kv[dj*kNI+di]+=uVal[sj_*jjUS+si_];
                    }
                }else{
                    // Column in kk's update region → update stack.
                    Size udj=dj-kF;
                    for(Size si_=sj_;si_<jjUS;++si_){
                        Size di=g2l[uIdx[si_]]; if(di==cNullIdx)continue;
                        Size udi=di-kF;
                        ukv[udj*kkUS+udi]+=uVal[sj_*jjUS+si_];
                    }
                }
            }
            // Also assemble delayed columns from jj factor (DynamicLDL).
            if(dynLDL) lu.assembleDelayed(jj,kk,g2l);
            if(dynLDL) lu.shrinkFront(jj,lu.nDelayed(jj));
            us.discard(jj);
        }

        // Factor the frontal matrix.
        Err e=Err::eNone;
        switch(mType){
            case FactorType::eCholesky:   e=factCC  (lu,kk); break;
            case FactorType::eStaticLDL:  e=factSLDL(lu,kk); break;
            case FactorType::eDynamicLDL: e=factDLDL(lu,kk,g2l); break;
        }
        if(e!=Err::eNone) return e;
        if(dynLDL&&lu.nDelayed(kk)>0&&sPrnt[kk]==cNullIdx) return Err::eInvPivot;

        // Compute Schur complement and ADD to us[kk] (which may already
        // contain bypass residuals from children).  All mfUpd* functions
        // use beta=1 so they accumulate rather than overwrite.
        if(kkUS>0){
            switch(mType){
                case FactorType::eCholesky:   mfUpdCC  (lu,us,kk); break;
                case FactorType::eStaticLDL:  mfUpdSLDL(lu,us,kk); break;
                case FactorType::eDynamicLDL: mfUpdDLDL(lu,us,kk); break;
            }
        }
        clrG2L(lu,kk,g2l);
    }
    lu.symmetrize(); return Err::eNone;
}

// ── Explicit instantiations ──────────────────────────────────────────────────
template class FactorEngine<double>;
template class FactorEngine<std::complex<double>>;

} // namespace Oblio
