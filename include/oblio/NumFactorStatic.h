#pragma once

// NumFactorStatic.h - the numeric factorization of a sparse matrix, statically stored.
//
// Static means the structure is fixed before any arithmetic runs: symbolic factorization has
// already sized every supernode's block, and nothing expands. That covers Cholesky and static LDL.
// Dynamic LDL does expand, by delaying unstable pivots into an ancestor, and gets its own class
// (NumFactorDynamic), because the storage that suits expansion does not suit this.
//
// This object is SymFactor plus the values. It copies what it needs from SymFactor rather than
// referring to it, exactly as SymFactor copies from ElmForest: each object is self-contained, so
// SymFactor may be discarded once factorization is done, and the solve reads only this.
//
// Storage is flat, one contiguous buffer with per-supernode offsets:
//
//   snodeNodeIdxPtr -> nodeIdx    the index sets, copied from SymFactor
//   snodeValPtr     -> val        the values, computed here
//
// Supernode kk's block is val[snodeValPtr[kk] .. snodeValPtr[kk + 1]), a **dense column-major rectangle**:
//
//   rows    = frontSize(kk) + updateSize(kk)   = the index-set size
//   columns = frontSize(kk)
//   leading dimension = the index-set size
//
// Entry (li, lj) of the block, in local coordinates, is at snodeValPtr[kk] + lj * indexSize + li.
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

template<class Val>
class NumFactorStatic {
public:
    NumFactorStatic() = default;

    std::size_t   size()      const { return mSize; }
    std::size_t   snodeSize() const { return mSnodeSize; }
    Factorization factorization() const { return mFactorization; }

    // Total node indices, and total values. The second is the rectangle count: the sum over
    // supernodes of indexSize * frontSize.
    std::size_t numNodeIdx() const { return mNumNodeIdx; }
    std::size_t numVal()     const { return mNumVal; }

    // How many pivots the factorization had to replace.
    //
    // Only LDL can perturb, and only because it must: a *static* factorization does not pivot, so
    // a pivot too small to divide by has no remedy but replacement. A nonzero count means we
    // factored a matrix slightly different from the one handed to us, and the caller is entitled
    // to know. Cholesky never perturbs; it fails instead, which is what positive definiteness
    // entitles it to do.
    std::size_t numPerturbations() const { return mNumPerturbations; }

    // Node to supernode: the whole map, and the supernode a single node belongs to. Indexed by a
    // *node*, unlike frontSize and its kin, which are indexed by a supernode.
    const std::vector<std::int32_t>& nodeToSnode()                  const { return mNodeToSnode; }
    std::int32_t                     nodeToSnode(std::int32_t node) const { return mNodeToSnode[node]; }

    // Per supernode: its own columns, and the rows below them.
    const std::vector<std::size_t>& frontSize()  const { return mFrontSize; }
    const std::vector<std::size_t>& updateSize() const { return mUpdateSize; }
    std::size_t frontSize(std::int32_t kk)  const { return mFrontSize[kk]; }
    std::size_t updateSize(std::int32_t kk) const { return mUpdateSize[kk]; }

    // The index sets, flat: offsets (length snodeSize() + 1), then node indices.
    const std::vector<std::size_t>&  snodeNodeIdxPtr() const { return mSnodeNodeIdxPtr; }
    const std::vector<std::int32_t>& nodeIdx()         const { return mNodeIdx; }

    // The values, flat: offsets (length snodeSize() + 1), then the dense blocks.
    const std::vector<std::size_t>& snodeValPtr() const { return mSnodeValPtr; }
    const std::vector<Val>&         val()         const { return mVal; }

    // Where supernode kk's node set and dense block live.
    //
    // **A lookup, not a view**, and that is why it belongs on this class rather than on the
    // engines. It computes an address inside the existing storage: no allocation, nothing
    // materialized, nothing owned. The layout it knows about (flat buffers, offsets in
    // mSnodeNodeIdxPtr and mSnodeValPtr) is a fact about *this class*, and no consumer should have
    // to restate it.
    //
    // The matrix's per-column lookups (rowIdx / val / colSize in experiments/storage-options)
    // are the same idea over a different storage: a fact about the layout, answered by the class
    // that owns it. NumFactorDynamic supplies its own nodeIdx and val, over its own layout,
    // and the engines cannot tell them apart.
    //
    // Reading is public: a consumer that only reads the factor (SolveEngine) reaches these through
    // the const overloads and needs no friendship. The mutable overloads, which hand out a writable
    // block, stay private for the engine that fills them.
    //
    // **Call it at the moment of use, never hoist it.** In the dynamic factor a delayed pivot
    // expands an ancestor's front, which reallocates its buffer, which dangles every pointer
    // previously taken into it, silently. experiments/storage-options demonstrates the rule
    // (structural mutation invalidates, value mutation does not) and measures the cost of obeying
    // it: one indirection, which is nothing.
    const std::int32_t* nodeIdx(std::int32_t kk) const { return mNodeIdx.data() + mSnodeNodeIdxPtr[kk]; }
    const Val*          val(std::int32_t kk)     const { return mVal.data()     + mSnodeValPtr[kk]; }

private:
    // The write path: NumFactorEngine fills each supernode's block through these, reached only
    // through friendship.
    std::int32_t*       nodeIdx(std::int32_t kk) { return mNodeIdx.data() + mSnodeNodeIdxPtr[kk]; }
    Val*                val(std::int32_t kk)     { return mVal.data()     + mSnodeValPtr[kk]; }

    // Also the write path: the engine accumulates the perturbation count through this reference
    // (factorStaticSupernode increments it). The const read overload above is public: the caller is
    // entitled to the count.
    std::size_t& numPerturbations() { return mNumPerturbations; }

    std::size_t   mSize      = 0;
    std::size_t   mSnodeSize = 0;
    Factorization mFactorization = Factorization::Cholesky;

    // Copied from SymFactor.
    std::vector<std::int32_t> mNodeToSnode;
    std::vector<std::size_t>  mFrontSize;
    std::vector<std::size_t>  mUpdateSize;
    std::size_t               mNumNodeIdx = 0;
    std::vector<std::size_t>  mSnodeNodeIdxPtr;
    std::vector<std::int32_t> mNodeIdx;

    // Computed here. mVal holds every supernode's dense block, end to end.
    std::size_t              mNumVal = 0;
    std::vector<std::size_t> mSnodeValPtr;
    std::vector<Val>         mVal;

    std::size_t mNumPerturbations = 0;

    friend class NumFactorEngine;
};

extern template class NumFactorStatic<double>;
extern template class NumFactorStatic<std::complex<double>>;

} // namespace Oblio
