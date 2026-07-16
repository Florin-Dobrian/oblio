#pragma once

// NumFactorStatic.h - the numeric factorization of a sparse matrix, statically stored.
//
// Static means the structure is fixed before any arithmetic runs: symbolic factorization has
// already sized every supernode's block, and nothing grows. That covers Cholesky and static LDL.
// Dynamic LDL does grow, by delaying unstable pivots into an ancestor, and gets its own class
// (NumFactorDynamic), because the storage that suits growth does not suit this.
//
// This object is SymFactor plus the values. It copies what it needs from SymFactor rather than
// referring to it, exactly as SymFactor copies from ElmForest: each object is self-contained, so
// SymFactor may be discarded once factorization is done, and the solve reads only this.
//
// Storage is flat, one contiguous buffer with per-supernode offsets:
//
//   supPtr -> rowIdx     the index sets, copied from SymFactor
//   valPtr -> val        the values, computed here
//
// Supernode kk's block is val[valPtr[kk] .. valPtr[kk + 1]), a **dense column-major rectangle**:
//
//   rows    = frontSize(kk) + updateSize(kk)   = the index-set size
//   columns = frontSize(kk)
//   leading dimension = the index-set size
//
// Entry (li, lj) of the block, in local coordinates, is at valPtr[kk] + lj * indexSize + li.
//
// Note the block is a **rectangle, not a trapezoid**. The strictly upper triangle of the front is
// allocated and left as zeros, because BLAS wants a rectangular block with a leading dimension,
// and paying for `frontSize * (frontSize - 1) / 2` zeros per supernode is the price of handing
// the whole thing to a level-3 kernel. 0.9 counts these two quantities separately for exactly
// this reason: numberOfAllocatedEntries is the rectangle, numberOfEntries the trapezoid.
//
// The engine writes this object through friendship; nothing else may.

#include "oblio/Types.h"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

class NumFactorEngine;
class SolveEngine;

template<class Val>
class NumFactorStatic {
public:
    NumFactorStatic() = default;

    std::size_t   size()      const { return mSize; }
    std::size_t   supSize()   const { return mSupSize; }
    Factorization factorization() const { return mFactorization; }

    // Total row indices, and total values. The second is the rectangle count: the sum over
    // supernodes of indexSize * frontSize.
    std::size_t numRowIdx() const { return mNumRowIdx; }
    std::size_t numVal()    const { return mNumVal; }

    // How many pivots the factorization had to replace.
    //
    // Only LDL can perturb, and only because it must: a *static* factorization does not pivot, so
    // a pivot too small to divide by has no remedy but replacement. A nonzero count means we
    // factored a matrix slightly different from the one handed to us, and the caller is entitled
    // to know. Cholesky never perturbs; it fails instead, which is what positive definiteness
    // entitles it to do.
    std::size_t numPerturbations() const { return mNumPerturbations; }

    // Column to supernode.
    const std::vector<std::int32_t>& idxToSupIdx() const { return mIdxToSupIdx; }

    // Per supernode: its own columns, and the rows below them.
    const std::vector<std::size_t>& frontSize()  const { return mFrontSize; }
    const std::vector<std::size_t>& updateSize() const { return mUpdateSize; }

    // The index sets, flat: offsets (length supSize() + 1), then row indices.
    const std::vector<std::size_t>&  supPtr() const { return mSupPtr; }
    const std::vector<std::int32_t>& rowIdx() const { return mRowIdx; }

    // The values, flat: offsets (length supSize() + 1), then the dense blocks.
    const std::vector<std::size_t>& valPtr() const { return mValPtr; }
    const std::vector<Val>&         val()    const { return mVal; }

private:
    // Where supernode kk's dense block lives.
    //
    // **A lookup, not a view**, and that is why it belongs here rather than on the engines. It
    // computes an address inside the existing storage: no allocation, nothing materialized,
    // nothing owned. The layout it knows about (one flat buffer, offsets in mValPtr) is a fact
    // about *this class*, and no consumer should have to restate it.
    //
    // The matrix's per-column lookups (rowIdx / val / colSize in experiments/storage-options)
    // are the same idea over a different storage: a fact about the layout, answered by the class
    // that owns it. NumFactorDynamic supplies its own blockPtr, over its own layout, and the
    // engines cannot tell them apart.
    //
    // **Call it at the moment of use, never hoist it.** In the dynamic factor a delayed pivot
    // grows an ancestor's front, which reallocates its buffer, which dangles every pointer
    // previously taken into it, silently. experiments/storage-options demonstrates the rule
    // (structural mutation invalidates, value mutation does not) and measures the cost of obeying
    // it: one indirection, which is nothing.
    Val*       blockPtr(std::int32_t kk)       { return mVal.data() + mValPtr[kk]; }
    const Val* blockPtr(std::int32_t kk) const { return mVal.data() + mValPtr[kk]; }

    std::size_t   mSize    = 0;
    std::size_t   mSupSize = 0;
    Factorization mFactorization = Factorization::Cholesky;

    // Copied from SymFactor.
    std::vector<std::int32_t> mIdxToSupIdx;
    std::vector<std::size_t>  mFrontSize;
    std::vector<std::size_t>  mUpdateSize;
    std::size_t               mNumRowIdx = 0;
    std::vector<std::size_t>  mSupPtr;
    std::vector<std::int32_t> mRowIdx;

    // Computed here. mVal holds every supernode's dense block, end to end.
    std::size_t              mNumVal = 0;
    std::vector<std::size_t> mValPtr;
    std::vector<Val>         mVal;

    std::size_t mNumPerturbations = 0;

    friend class NumFactorEngine;
    friend class SolveEngine;
};

extern template class NumFactorStatic<double>;
extern template class NumFactorStatic<std::complex<double>>;

} // namespace Oblio
