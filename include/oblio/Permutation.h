#pragma once
#include "oblio/Types.h"
#include <vector>

namespace Oblio {

class Permutation {
public:
    Permutation() : mSize(0) {}
    explicit Permutation(Size n) { resize(n); }

    void resize(Size n) {
        mSize = n;
        mOldToNew.assign(n, cNullIdx);
        mNewToOld.assign(n, cNullIdx);
    }
    void setIdentity() {
        for (Size i=0;i<mSize;++i) { mOldToNew[i]=i; mNewToOld[i]=i; }
    }

    Size getSize() const { return mSize; }
    bool isValid() const {
        if (mOldToNew.size()!=mSize||mNewToOld.size()!=mSize) return false;
        for (Size i=0;i<mSize;++i) if(mOldToNew[i]>=mSize||mNewToOld[i]>=mSize) return false;
        return true;
    }

    Size oldToNew(Size i) const { return mOldToNew[i]; }
    Size newToOld(Size j) const { return mNewToOld[j]; }
    const Size* oldToNewData() const { return mOldToNew.data(); }
    const Size* newToOldData() const { return mNewToOld.data(); }
    Size*       oldToNewData()       { return mOldToNew.data(); }
    Size*       newToOldData()       { return mNewToOld.data(); }

private:
    Size              mSize;
    std::vector<Size> mOldToNew;
    std::vector<Size> mNewToOld;
};

} // namespace Oblio
