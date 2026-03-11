#include "oblio/Factors.h"
#include <complex>

namespace Oblio {

// ============================================================================
// Factors<Val>
// ============================================================================

template<class Val>
Factors<Val>::Factors()
    : mSize(0), mNumSnodes(0), mFactorType(FactorType::eCholesky) {}

template<class Val>
void Factors<Val>::allocate(const Symbolic& s) {
    mSize = s.mSize; mNumSnodes = s.mNumSnodes;
    mFrntSz .assign(mNumSnodes, 0);
    mUpdtSz .assign(mNumSnodes, 0);
    mNDelayed.assign(mNumSnodes, 0);
    mIdx    .resize(mNumSnodes);
    mVal    .resize(mNumSnodes);
    mPivot  .resize(mNumSnodes);
    for (Size jj = 0; jj < mNumSnodes; ++jj) {
        mFrntSz[jj] = s.mFrntSzVec[jj];
        mUpdtSz[jj] = s.mUpdtSzVec[jj];
        Size nIdx   = mFrntSz[jj] + mUpdtSz[jj];
        mIdx[jj].assign(s.mIdxVecVec[jj].begin(), s.mIdxVecVec[jj].end());
        mVal[jj].assign(mFrntSz[jj] * nIdx, Val{0});
        mPivot[jj].clear();
    }
}

template<class Val>
void Factors<Val>::zero()        { for (Size jj = 0; jj < mNumSnodes; ++jj) zero(jj); }

template<class Val>
void Factors<Val>::zero(Size jj) { std::fill(mVal[jj].begin(), mVal[jj].end(), Val{0}); }

template<class Val>
void Factors<Val>::symmetrize()  {}   // lower triangle sufficient for solve

template<class Val>
FactorType Factors<Val>::getFactorType()        const { return mFactorType; }

template<class Val>
void Factors<Val>::setFactorType(FactorType ft)       { mFactorType = ft; }

template<class Val>
Size Factors<Val>::getSize()                    const { return mSize; }

template<class Val>
Size Factors<Val>::getNumSnodes()               const { return mNumSnodes; }

template<class Val>
Size Factors<Val>::frntSz(Size jj)              const { return mFrntSz[jj]; }

template<class Val>
Size Factors<Val>::updtSz(Size jj)              const { return mUpdtSz[jj]; }

template<class Val>
Size Factors<Val>::nIdx(Size jj)                const { return mFrntSz[jj] + mUpdtSz[jj]; }

template<class Val>
Size Factors<Val>::nDelayed(Size jj)            const { return mNDelayed[jj]; }

template<class Val>
const Size* Factors<Val>::idx(Size jj) const { return mIdx[jj].data(); }

template<class Val>
Size*       Factors<Val>::idx(Size jj)       { return mIdx[jj].data(); }

template<class Val>
const Val*  Factors<Val>::val(Size jj) const { return mVal[jj].data(); }

template<class Val>
Val*        Factors<Val>::val(Size jj)       { return mVal[jj].data(); }

template<class Val>
const int*  Factors<Val>::piv(Size jj) const { return mPivot[jj].data(); }

template<class Val>
int*        Factors<Val>::piv(Size jj)       { return mPivot[jj].data(); }

template<class Val>
void Factors<Val>::setNDelayed(Size jj, Size n) { mNDelayed[jj] = n; }

template<class Val>
void Factors<Val>::appendPivot(Size jj, int t)  { mPivot[jj].push_back(t); }

template<class Val>
void Factors<Val>::extendFront(Size kk, Size nDelay,
                                const std::vector<Size>& fstChld,
                                const std::vector<Size>& nxtSblg) {
    Size old = mFrntSz[kk] + mUpdtSz[kk];
    mIdx[kk].resize(old + nDelay, cNullIdx);
    for (Size k = old; k-- > 0;) mIdx[kk][k + nDelay] = mIdx[kk][k];
    Size dst = 0;
    for (Size jj = fstChld[kk]; jj != cNullIdx; jj = nxtSblg[jj]) {
        Size f = mFrntSz[jj], nd = mNDelayed[jj];
        for (Size k = f; k < f + nd; ++k, ++dst) mIdx[kk][dst] = mIdx[jj][k];
    }
    mFrntSz[kk] += nDelay;
    mPivot[kk].clear();
}

template<class Val>
void Factors<Val>::reallocVal(Size kk) {
    Size nI = mFrntSz[kk] + mUpdtSz[kk];
    mVal[kk].assign(mFrntSz[kk] * nI, Val{0});
}

template<class Val>
void Factors<Val>::shrinkFront(Size jj, Size nDelay) {
    if (!nDelay) return;
    mFrntSz[jj]  -= nDelay;
    mNDelayed[jj]  = 0;
    mIdx[jj].erase(mIdx[jj].begin(), mIdx[jj].begin() + nDelay);
    Size nI = mFrntSz[jj] + mUpdtSz[jj];
    mVal[jj].resize(mFrntSz[jj] * nI);
}

template<class Val>
void Factors<Val>::swapCols(Size jj, Size a_, Size b_, std::vector<Size>& g2l) {
    if (a_ == b_) return;
    Size* ix = mIdx[jj].data();
    Val*  v  = mVal[jj].data();
    Size  f  = mFrntSz[jj];
    Size  nI = f + mUpdtSz[jj];
    std::swap(ix[a_], ix[b_]);
    g2l[ix[a_]] = a_;  g2l[ix[b_]] = b_;
    for (Size i = 0; i < nI; ++i) std::swap(v[a_ * nI + i], v[b_ * nI + i]);
    for (Size j = 0; j < f;  ++j) std::swap(v[j  * nI + a_], v[j * nI + b_]);
}

template<class Val>
void Factors<Val>::assembleDelayed(Size jj, Size kk, const std::vector<Size>& g2l) {
    const Size* ji = mIdx[jj].data();
    const Val*  jv = mVal[jj].data();
    Val*        kv = mVal[kk].data();
    Size jf   = mFrntSz[jj], jd = mNDelayed[jj];
    Size jnI  = jf + jd + mUpdtSz[jj];
    Size knI  = mFrntSz[kk] + mUpdtSz[kk];
    Size sp0j = jf * jnI;
    for (Size sj_ = jf; sj_ < jf + jd; ++sj_) {
        Size j = ji[sj_], dj = g2l[j];
        if (dj == cNullIdx) { sp0j += jnI; continue; }
        for (Size si_ = sj_; si_ < jnI; ++si_) {
            Size di = g2l[ji[si_]];
            if (di != cNullIdx) kv[dj * knI + di] = jv[sp0j + si_];
        }
        sp0j += jnI;
    }
}

// ============================================================================
// UpdateMatrix<Val>
// ============================================================================

template<class Val>
void UpdateMatrix<Val>::alloc(Size s) {
    sz = s;
    idx.assign(s, cNullIdx);
    val.assign(s * s, Val{0});
}

template<class Val>
void UpdateMatrix<Val>::setIdx(const Size* src, Size n) {
    assert(n == sz);
    std::copy(src, src + n, idx.begin());
}

template<class Val>
void UpdateMatrix<Val>::zero()  { std::fill(val.begin(), val.end(), Val{0}); }

template<class Val>
void UpdateMatrix<Val>::clear() { sz = 0; idx.clear(); val.clear(); }

// ============================================================================
// UpdateStack<Val>
// ============================================================================

template<class Val>
UpdateStack<Val>::UpdateStack(Size n) : um(n) {}

template<class Val>
void UpdateStack<Val>::alloc(Size jj, Size s)               { um[jj].alloc(s); }

template<class Val>
void UpdateStack<Val>::setIdx(Size jj, const Size* src, Size s) { um[jj].setIdx(src, s); }

template<class Val>
void UpdateStack<Val>::zero(Size jj)                         { um[jj].zero(); }

template<class Val>
void UpdateStack<Val>::discard(Size jj)                      { um[jj].clear(); }

template<class Val>
Size UpdateStack<Val>::sz(Size jj)           const { return um[jj].sz; }

template<class Val>
const Size* UpdateStack<Val>::idx(Size jj)   const { return um[jj].idx.data(); }

template<class Val>
const Val*  UpdateStack<Val>::val(Size jj)   const { return um[jj].val.data(); }

template<class Val>
Val*        UpdateStack<Val>::val(Size jj)         { return um[jj].val.data(); }

// ── Explicit instantiations ──────────────────────────────────────────────────
template class Factors<double>;
template class Factors<std::complex<double>>;
template struct UpdateMatrix<double>;
template struct UpdateMatrix<std::complex<double>>;
template struct UpdateStack<double>;
template struct UpdateStack<std::complex<double>>;

} // namespace Oblio
