#pragma once
#include "oblio/Types.h"
#include <vector>

namespace Oblio {

// Supernodal symbolic factorization.
// All arrays are indexed by supernode number jj in [0, mNumSnodes).
// mIdxVecVec[jj] = [front indices | update indices] (global, sorted).
struct Symbolic {
    Size mSize      = 0;
    Size mNumSnodes = 0;
    Size mNumTrees  = 0;
    Size mHeight    = 0;
    Size mFstRoot   = cNullIdx;
    Size mLstRoot   = cNullIdx;
    Size mNumIdxs   = 0;   // total indices across all snodes
    Size mNumVals   = 0;   // total values (= sum frontSz*nIdxs)

    std::vector<Size> mPrntVec;         // parent[jj]; cNullIdx for roots
    std::vector<Size> mFstChldVec;      // first child of jj
    std::vector<Size> mNxtSblgVec;      // next sibling of jj
    std::vector<Size> mIdxToSnodeVec;   // column j -> supernode
    std::vector<Size> mFrntSzVec;       // front columns owned by jj
    std::vector<Size> mUpdtSzVec;       // update (trailing) rows of jj
    std::vector<std::vector<Size>> mIdxVecVec;  // indices per snode

    // Accessors used by FactorEngine / Factors.
    Size getSize()      const { return mSize; }
    Size getNumSnodes() const { return mNumSnodes; }

    Size snodeFrontSize(Size jj) const { return mFrntSzVec[jj]; }
    Size snodeUpdtSize (Size jj) const { return mUpdtSzVec[jj]; }

    const Size* snodeIdxData(Size jj)   const { return mIdxVecVec[jj].data(); }
    const Size* frntSzData()            const { return mFrntSzVec.data(); }
    const Size* updtSzData()            const { return mUpdtSzVec.data(); }
    const Size* prntData()              const { return mPrntVec.data(); }
    const Size* fstChldData()           const { return mFstChldVec.data(); }
    const Size* nxtSblgData()           const { return mNxtSblgVec.data(); }
    const Size* snodeMapData()          const { return mIdxToSnodeVec.data(); }
};

} // namespace Oblio
