#include "oblio/SolveEngine.h"
#include "oblio/BlasLapack.h"
#include <vector>
#include <utility>
#include <complex>

namespace Oblio {

template<class Val>
SolveEngine<Val>::SolveEngine() : mScope(SolveScope::eFull) {}

template<class Val>
void SolveEngine<Val>::setScope(SolveScope s) { mScope = s; }

// ============================================================================
// Top-level dispatch
// ============================================================================

template<class Val>
Err SolveEngine<Val>::solve(const Factors<Val>& lu, Vector<Val>& y) const {
    if (y.size() != lu.getSize()) return Err::eDimMismatch;
    FactorType ft = lu.getFactorType();
    bool doFwd  = (mScope == SolveScope::eFull || mScope == SolveScope::eLower);
    bool doDiag = (mScope == SolveScope::eFull || mScope == SolveScope::eDiagonal);
    bool doBwd  = (mScope == SolveScope::eFull || mScope == SolveScope::eUpper);
    if (doFwd) {
        switch (ft) {
            case FactorType::eCholesky:   fwdCC   (lu, y); break;
            case FactorType::eStaticLDL:  fwdSLDL (lu, y); break;
            case FactorType::eDynamicLDL: fwdDLDL (lu, y); break;
        }
    }
    if (doDiag && ft != FactorType::eCholesky) {
        switch (ft) {
            case FactorType::eStaticLDL:  diagSLDL(lu, y); break;
            case FactorType::eDynamicLDL: diagDLDL(lu, y); break;
            default: break;
        }
    }
    if (doBwd) {
        switch (ft) {
            case FactorType::eCholesky:   bwdCC   (lu, y); break;
            case FactorType::eStaticLDL:  bwdSLDL (lu, y); break;
            case FactorType::eDynamicLDL: bwdDLDL (lu, y); break;
        }
    }
    return Err::eNone;
}

template<class Val>
Err SolveEngine<Val>::solve(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    if (Y.rows() != lu.getSize()) return Err::eDimMismatch;
    if (Y.empty()) return Err::eNone;
    FactorType ft = lu.getFactorType();
    bool doFwd  = (mScope == SolveScope::eFull || mScope == SolveScope::eLower);
    bool doDiag = (mScope == SolveScope::eFull || mScope == SolveScope::eDiagonal);
    bool doBwd  = (mScope == SolveScope::eFull || mScope == SolveScope::eUpper);
    if (doFwd) {
        switch (ft) {
            case FactorType::eCholesky:   fwdCC_m   (lu, Y); break;
            case FactorType::eStaticLDL:  fwdSLDL_m (lu, Y); break;
            case FactorType::eDynamicLDL: fwdDLDL_m (lu, Y); break;
        }
    }
    if (doDiag && ft != FactorType::eCholesky) {
        switch (ft) {
            case FactorType::eStaticLDL:  diagSLDL_m(lu, Y); break;
            case FactorType::eDynamicLDL: diagDLDL_m(lu, Y); break;
            default: break;
        }
    }
    if (doBwd) {
        switch (ft) {
            case FactorType::eCholesky:   bwdCC_m   (lu, Y); break;
            case FactorType::eStaticLDL:  bwdSLDL_m (lu, Y); break;
            case FactorType::eDynamicLDL: bwdDLDL_m (lu, Y); break;
        }
    }
    return Err::eNone;
}

// ============================================================================
// Single-RHS passes
// ============================================================================

template<class Val>
void SolveEngine<Val>::fwdCC(const Factors<Val>& lu, Vector<Val>& y) const {
    Size nSnd = lu.getNumSnodes();
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (Size j_ = 0, pjj = 0; j_ < f; ++j_, pjj += nI + 1) {
            y[ix[j_]] /= v[pjj];
            for (Size i_ = j_+1, pij = pjj+1; i_ < nI; ++i_, ++pij)
                y[ix[i_]] -= v[pij] * y[ix[j_]];
        }
    }
}

template<class Val>
void SolveEngine<Val>::fwdSLDL(const Factors<Val>& lu, Vector<Val>& y) const {
    Size nSnd = lu.getNumSnodes();
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (Size j_ = 0, pjj = 0; j_ < f; ++j_, pjj += nI + 1)
            for (Size i_ = j_+1, pij = pjj+1; i_ < nI; ++i_, ++pij)
                y[ix[i_]] -= v[pij] * y[ix[j_]];
    }
}

template<class Val>
void SolveEngine<Val>::fwdDLDL(const Factors<Val>& lu, Vector<Val>& y) const {
    Size nSnd = lu.getNumSnodes();
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        const int* pt  = lu.piv(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (Size j_ = 0, pjj = 0; j_ < f; ) {
            if (pt[j_] == 1) {
                for (Size i_ = j_+1, pij = pjj+1; i_ < nI; ++i_, ++pij)
                    y[ix[i_]] -= v[pij] * y[ix[j_]];
                ++j_; pjj += nI + 1;
            } else {
                for (Size i_ = j_+2, pij = pjj+2; i_ < nI; ++i_, ++pij)
                    y[ix[i_]] -= v[pij] * y[ix[j_]];
                for (Size i_ = j_+2, pij = pjj+nI+2; i_ < nI; ++i_, ++pij)
                    y[ix[i_]] -= v[pij] * y[ix[j_+1]];
                j_ += 2; pjj += 2*nI + 2;
            }
        }
    }
}

template<class Val>
void SolveEngine<Val>::diagSLDL(const Factors<Val>& lu, Vector<Val>& y) const {
    Size nSnd = lu.getNumSnodes();
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (Size j_ = 0, pjj = 0; j_ < f; ++j_, pjj += nI + 1) {
            Val d = v[pjj];
            if (d != Val{0}) y[ix[j_]] /= d;
        }
    }
}

template<class Val>
void SolveEngine<Val>::diagDLDL(const Factors<Val>& lu, Vector<Val>& y) const {
    Size nSnd = lu.getNumSnodes();
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        const int* pt  = lu.piv(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (Size j_ = 0, pjj = 0; j_ < f; ) {
            if (pt[j_] == 1) {
                Val d = v[pjj];
                if (d != Val{0}) y[ix[j_]] /= d;
                ++j_; pjj += nI + 1;
            } else {
                Val a = v[pjj], b = v[pjj+nI], c = v[pjj+nI+1];
                Val det = a*c - b*b;
                if (std::abs(det) > 0) {
                    Val y1 = y[ix[j_]], y2 = y[ix[j_+1]];
                    y[ix[j_  ]] = (c*y1 - b*y2) / det;
                    y[ix[j_+1]] = (a*y2 - b*y1) / det;
                }
                j_ += 2; pjj += 2*nI + 2;
            }
        }
    }
}

template<class Val>
void SolveEngine<Val>::bwdCC(const Factors<Val>& lu, Vector<Val>& y) const {
    Size nSnd = lu.getNumSnodes();
    for (Size jj = nSnd; jj-- > 0; ) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (int i_ = (int)f-1; i_ >= 0; --i_) {
            int pii = i_*(int)nI + i_;
            for (Size j_ = i_+1, pij = pii+1; j_ < nI; ++j_, ++pij)
                y[ix[i_]] -= conjv(v[pij]) * y[ix[j_]];
            y[ix[i_]] /= conjv(v[pii]);
        }
    }
}

template<class Val>
void SolveEngine<Val>::bwdSLDL(const Factors<Val>& lu, Vector<Val>& y) const {
    Size nSnd = lu.getNumSnodes();
    for (Size jj = nSnd; jj-- > 0; ) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (int i_ = (int)f-1; i_ >= 0; --i_) {
            int pii = i_*(int)nI + i_;
            for (Size j_ = i_+1, pij = pii+1; j_ < nI; ++j_, ++pij)
                y[ix[i_]] -= v[pij] * y[ix[j_]];
        }
    }
}

template<class Val>
void SolveEngine<Val>::bwdDLDL(const Factors<Val>& lu, Vector<Val>& y) const {
    Size nSnd = lu.getNumSnodes();
    for (Size jj = nSnd; jj-- > 0; ) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        const int* pt  = lu.piv(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        std::vector<std::pair<Size,Size>> pivots;
        for (Size j_ = 0, pjj = 0; j_ < f; ) {
            pivots.push_back({j_, pjj});
            if (pt[j_] == 1) { ++j_; pjj += nI+1; }
            else              { j_ += 2; pjj += 2*nI+2; }
        }
        for (Size pi = pivots.size(); pi-- > 0; ) {
            Size j_ = pivots[pi].first, pjj = pivots[pi].second;
            if (pt[j_] == 1) {
                for (Size k_ = j_+1, pkj = pjj+1; k_ < nI; ++k_, ++pkj)
                    y[ix[j_]] -= v[pkj] * y[ix[k_]];
            } else {
                for (Size k_ = j_+2, pk = pjj+2;    k_ < nI; ++k_, ++pk)
                    y[ix[j_  ]] -= v[pk] * y[ix[k_]];
                for (Size k_ = j_+2, pk = pjj+nI+2; k_ < nI; ++k_, ++pk)
                    y[ix[j_+1]] -= v[pk] * y[ix[k_]];
            }
        }
    }
}

// ============================================================================
// Multi-RHS helpers
// ============================================================================

template<class Val>
static void gather(const Size* ix, Size nI, const DenseMatrix<Val>& Y,
                   std::vector<Val>& G) {
    Size nRHS = Y.cols();
    G.resize(nI * nRHS);
    for (Size r = 0; r < nRHS; ++r)
        for (Size k = 0; k < nI; ++k)
            G[r*nI + k] = Y(ix[k], r);
}

template<class Val>
static void scatter(const Size* ix, Size nI, const std::vector<Val>& G,
                    DenseMatrix<Val>& Y) {
    Size nRHS = Y.cols();
    for (Size r = 0; r < nRHS; ++r)
        for (Size k = 0; k < nI; ++k)
            Y(ix[k], r) = G[r*nI + k];
}

// ============================================================================
// Multi-RHS passes
// ============================================================================

template<class Val>
void SolveEngine<Val>::fwdCC_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    Size nSnd = lu.getNumSnodes(); int nRHS = toBI(Y.cols());
    std::vector<Val> G;
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        gather(ix, nI, Y, G);
        int inf = toBI(nI), iff = toBI(f);
        Val one{1}, mone{-1};
        BT<Val>::trsm('L','L','N','N', iff, nRHS, one, v, inf, G.data(), inf);
        if (nI > f)
            BT<Val>::gemm('N','N', toBI(nI-f), nRHS, iff, mone,
                          v+f, inf, G.data(), inf, one, G.data()+f, inf);
        scatter(ix, nI, G, Y);
    }
}

template<class Val>
void SolveEngine<Val>::fwdSLDL_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    Size nSnd = lu.getNumSnodes(); int nRHS = toBI(Y.cols());
    std::vector<Val> G;
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        gather(ix, nI, Y, G);
        int inf = toBI(nI), iff = toBI(f);
        Val one{1}, mone{-1};
        BT<Val>::trsm('L','L','N','U', iff, nRHS, one, v, inf, G.data(), inf);
        if (nI > f)
            BT<Val>::gemm('N','N', toBI(nI-f), nRHS, iff, mone,
                          v+f, inf, G.data(), inf, one, G.data()+f, inf);
        scatter(ix, nI, G, Y);
    }
}

template<class Val>
void SolveEngine<Val>::fwdDLDL_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    Size nSnd = lu.getNumSnodes(); Size nRHS = Y.cols();
    std::vector<Val> G;
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        const int* pt  = lu.piv(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        gather(ix, nI, Y, G);
        for (Size j_ = 0, pjj = 0; j_ < f; ) {
            if (pt[j_] == 1) {
                for (Size i_ = j_+1; i_ < nI; ++i_)
                    for (Size r = 0; r < nRHS; ++r)
                        G[r*nI+i_] -= v[pjj+(i_-j_)] * G[r*nI+j_];
                ++j_; pjj += nI+1;
            } else {
                for (Size i_ = j_+2; i_ < nI; ++i_)
                    for (Size r = 0; r < nRHS; ++r) {
                        G[r*nI+i_] -= v[pjj+(i_-j_)]         * G[r*nI+j_];
                        G[r*nI+i_] -= v[pjj+nI+1+(i_-j_-1)] * G[r*nI+j_+1];
                    }
                j_ += 2; pjj += 2*nI+2;
            }
        }
        scatter(ix, nI, G, Y);
    }
}

template<class Val>
void SolveEngine<Val>::diagSLDL_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    Size nSnd = lu.getNumSnodes(); Size nRHS = Y.cols();
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (Size j_ = 0, pjj = 0; j_ < f; ++j_, pjj += nI+1) {
            Val d = v[pjj]; if (d == Val{0}) continue;
            Size row = ix[j_];
            for (Size r = 0; r < nRHS; ++r) Y(row,r) /= d;
        }
    }
}

template<class Val>
void SolveEngine<Val>::diagDLDL_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    Size nSnd = lu.getNumSnodes(); Size nRHS = Y.cols();
    for (Size jj = 0; jj < nSnd; ++jj) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        const int* pt  = lu.piv(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        for (Size j_ = 0, pjj = 0; j_ < f; ) {
            if (pt[j_] == 1) {
                Val d = v[pjj];
                if (d != Val{0}) { Size row=ix[j_]; for(Size r=0;r<nRHS;++r) Y(row,r)/=d; }
                ++j_; pjj += nI+1;
            } else {
                Val a=v[pjj], b=v[pjj+nI], c=v[pjj+nI+1], det=a*c-b*b;
                if (std::abs(det) > 0) {
                    Size r1=ix[j_], r2=ix[j_+1];
                    for (Size r=0; r<nRHS; ++r) {
                        Val y1=Y(r1,r), y2=Y(r2,r);
                        Y(r1,r)=(c*y1-b*y2)/det; Y(r2,r)=(a*y2-b*y1)/det;
                    }
                }
                j_ += 2; pjj += 2*nI+2;
            }
        }
    }
}

template<class Val>
void SolveEngine<Val>::bwdCC_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    Size nSnd = lu.getNumSnodes(); int nRHS = toBI(Y.cols());
    std::vector<Val> G;
    for (Size jj = nSnd; jj-- > 0; ) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        gather(ix, nI, Y, G);
        int inf = toBI(nI), iff = toBI(f);
        Val one{1}, mone{-1};
        char ct = conjTrans<Val>();
        if (nI > f)
            BT<Val>::gemm(ct,'N', iff, nRHS, toBI(nI-f), mone,
                          v+f, inf, G.data()+f, inf, one, G.data(), inf);
        BT<Val>::trsm('L','L',ct,'N', iff, nRHS, one, v, inf, G.data(), inf);
        scatter(ix, nI, G, Y);
    }
}

template<class Val>
void SolveEngine<Val>::bwdSLDL_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    Size nSnd = lu.getNumSnodes(); int nRHS = toBI(Y.cols());
    std::vector<Val> G;
    for (Size jj = nSnd; jj-- > 0; ) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        gather(ix, nI, Y, G);
        int inf = toBI(nI), iff = toBI(f);
        Val one{1}, mone{-1};
        if (nI > f)
            BT<Val>::gemm('T','N', iff, nRHS, toBI(nI-f), mone,
                          v+f, inf, G.data()+f, inf, one, G.data(), inf);
        BT<Val>::trsm('L','L','T','U', iff, nRHS, one, v, inf, G.data(), inf);
        scatter(ix, nI, G, Y);
    }
}

template<class Val>
void SolveEngine<Val>::bwdDLDL_m(const Factors<Val>& lu, DenseMatrix<Val>& Y) const {
    Size nSnd = lu.getNumSnodes(); Size nRHS = Y.cols();
    std::vector<Val> G;
    for (Size jj = nSnd; jj-- > 0; ) {
        const Size* ix = lu.idx(jj); const Val* v = lu.val(jj);
        const int* pt  = lu.piv(jj);
        Size f = lu.frntSz(jj), nI = lu.nIdx(jj);
        gather(ix, nI, Y, G);
        std::vector<std::pair<Size,Size>> pivots;
        for (Size j_ = 0, pjj = 0; j_ < f; ) {
            pivots.push_back({j_, pjj});
            if (pt[j_]==1) { ++j_; pjj+=nI+1; } else { j_+=2; pjj+=2*nI+2; }
        }
        for (Size pi = pivots.size(); pi-- > 0; ) {
            Size j_ = pivots[pi].first, pjj = pivots[pi].second;
            if (pt[j_] == 1) {
                for (Size k_ = j_+1; k_ < nI; ++k_)
                    for (Size r = 0; r < nRHS; ++r)
                        G[r*nI+j_] -= v[pjj+(k_-j_)] * G[r*nI+k_];
            } else {
                for (Size k_ = j_+2; k_ < nI; ++k_)
                    for (Size r = 0; r < nRHS; ++r) {
                        G[r*nI+j_  ] -= v[pjj+(k_-j_)]         * G[r*nI+k_];
                        G[r*nI+j_+1] -= v[pjj+nI+1+(k_-j_-1)] * G[r*nI+k_];
                    }
            }
        }
        scatter(ix, nI, G, Y);
    }
}

// ============================================================================
template class SolveEngine<double>;
template class SolveEngine<std::complex<double>>;

} // namespace Oblio
