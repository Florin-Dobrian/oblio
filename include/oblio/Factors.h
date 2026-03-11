#pragma once
#include "oblio/Types.h"
#include "oblio/Symbolic.h"
#include <vector>
#include <complex>
#include <algorithm>
#include <cassert>

namespace Oblio {

template<class Val> class FactorEngine;
template<class Val> class SolveEngine;

// Per-supernode factor data.
// Column-major: shape (frontSz + updtSz) x frontSz, ldim = frontSz + updtSz.
template<class Val>
class Factors {
public:
    Factors();

    void allocate(const Symbolic& s);
    void zero();
    void zero(Size jj);
    void symmetrize();

    // ---- type / size -------------------------------------------------------
    FactorType getFactorType()        const;
    void       setFactorType(FactorType ft);
    Size       getSize()              const;
    Size       getNumSnodes()         const;
    Size       frntSz(Size jj)        const;
    Size       updtSz(Size jj)        const;
    Size       nIdx  (Size jj)        const;
    Size       nDelayed(Size jj)      const;

    // ---- data pointers -----------------------------------------------------
    const Size* idx(Size jj) const;
    Size*       idx(Size jj);
    const Val*  val(Size jj) const;
    Val*        val(Size jj);
    const int*  piv(Size jj) const;
    int*        piv(Size jj);

    // ---- DynamicLDL mutations ----------------------------------------------
    void setNDelayed(Size jj, Size n);
    void appendPivot(Size jj, int t);
    void extendFront(Size kk, Size nDelay,
                     const std::vector<Size>& fstChld,
                     const std::vector<Size>& nxtSblg);
    void reallocVal(Size kk);
    void shrinkFront(Size jj, Size nDelay);
    void swapCols(Size jj, Size a_, Size b_, std::vector<Size>& g2l);
    void assembleDelayed(Size jj, Size kk, const std::vector<Size>& g2l);

private:
    Size       mSize, mNumSnodes;
    FactorType mFactorType;
    std::vector<Size>              mFrntSz, mUpdtSz, mNDelayed;
    std::vector<std::vector<Size>> mIdx;
    std::vector<std::vector<Val>>  mVal;
    std::vector<std::vector<int>>  mPivot;

    Factors(const Factors&)            = delete;
    Factors& operator=(const Factors&) = delete;

    template<class V> friend class FactorEngine;
    template<class V> friend class SolveEngine;
};

// ============================================================================
// Update matrix / stack (multifrontal)
// ============================================================================
template<class Val>
struct UpdateMatrix {
    Size              sz = 0;
    std::vector<Size> idx;
    std::vector<Val>  val;

    void alloc(Size s);
    void setIdx(const Size* src, Size n);
    void zero();
    void clear();
};

template<class Val>
struct UpdateStack {
    explicit UpdateStack(Size n);
    std::vector<UpdateMatrix<Val>> um;

    void         alloc  (Size jj, Size s);
    void         setIdx (Size jj, const Size* src, Size s);
    void         zero   (Size jj);
    void         discard(Size jj);
    Size         sz     (Size jj) const;
    const Size*  idx    (Size jj) const;
    const Val*   val    (Size jj) const;
    Val*         val    (Size jj);
};

// Explicit instantiation declarations — definitions in Factors.cc
extern template class Factors<double>;
extern template class Factors<std::complex<double>>;
extern template struct UpdateMatrix<double>;
extern template struct UpdateMatrix<std::complex<double>>;
extern template struct UpdateStack<double>;
extern template struct UpdateStack<std::complex<double>>;

} // namespace Oblio
