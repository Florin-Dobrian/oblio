#pragma once

// SymFactor.h - the symbolic factorization of a sparse matrix.
//
// The symbolic factorization is the sparsity pattern of the factor L, computed
// per supernode: for every supernode, the set of factor node indices its columns
// touch. That index set splits into the front indices (the supernode's own
// columns) and the update indices (the rows below it, which its children and it
// itself update). Numeric factorization consumes this structure.
//
// SymFactor is self-contained: it copies the forest attributes it needs (the map,
// the parent/child/sibling links, the front and update sizes, the tree
// attributes), so numeric factorization needs SymFactor alone, not the forest.
// Both 0.9 and 10.12 do this.
//
// Storage is flat, and deliberately the same shape as SparseMatrix. The index sets are
// written once and their sizes are known in advance from the forest (front size + update
// size per supernode), so they go in one contiguous array of node indices with per-supernode
// offsets. Exactly as in 0.9; 10.12 uses a vector of vectors here (see the flat-vs-VV
// decision).
//
// The parallel with SparseMatrix is exact in shape, with L's names marking the compression (see
// the A-and-L naming decision): L holds one node index array with per-supernode offsets where A
// holds one row index array with per-column offsets.
//
//   SparseMatrix:  rowIdx[colPtr[j] .. colPtr[j + 1])        the rows of column j
//   SymFactor:     nodeIdx[snodePtr[s] .. snodePtr[s + 1])   the nodes of supernode s
//
// A supernode's nodes are the pattern its columns share: the first frontSize(s) of them are
// its front indices (its own columns), the rest its update indices (the rows below). Within
// a supernode they are sorted in increasing order.
//
// Index roles (see the index-types decision): supernode indices, column indices,
// factor node indices and the links are IDs -> std::int32_t with NIL = -1; sizes,
// counts and offsets are std::size_t.

#include "oblio/Types.h"

#include <vector>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class SymFactorEngine;

class SymFactor {
public:
    SymFactor() = default;

    std::size_t size()       const { return mSize; }         // number of columns
    std::size_t snodeSize()  const { return mSnodeSize; }    // number of supernodes
    std::size_t numTrees()   const { return mNumTrees; }     // number of trees (roots)
    std::size_t height()     const { return mHeight; }       // forest height (max depth + 1)
    // Three sizes of the factor, in the same order on all four classes that can answer. They count
    // different things and the difference grows with front width, so the names are not
    // interchangeable:
    //
    //   numNodeIdx  the index sets, one entry per row a supernode touches
    //   numVal      values *allocated*: the full rectangle, whose front block is stored whole with
    //               its strict upper triangle left zero so BLAS can take it in one call
    //   nnz         entries of L: each supernode's own lower triangle plus its update rectangle
    //
    // In magnitude they nest, numNodeIdx <= nnz <= numVal, with equality only when every
    // supernode is one column. On an 8x8 grid Laplacian under AMD, nnz is 354 and numVal 375;
    // amalgamated at 8 they are 499 and 652. 0.9 spelled the last two numberOfEntries and
    // numberOfAllocatedEntries.
    //
    // **nnz counts entries, not nonzeros.** Amalgamation merges supernodes and pads their index
    // sets with rows that are structurally zero, and those are counted here. The two agree exactly
    // when exactPatterns() is true.
    //
    // Summed on demand rather than maintained, deliberately: the scan is linear in supernodes and
    // measured at under 100 microseconds on a 400x400 grid, about 1e-4 of that factorization, while
    // a maintained counter would need updating wherever the structure changes and fails silently
    // when one site is missed.
    //
    // numNodeIdx is a member here rather than a sum: nodeIdx is one flat array and its length falls
    // out of the prefix sum that builds it. The other two are predictions, no values existing yet.
    std::size_t numNodeIdx() const { return mNumNodeIdx; }
    std::size_t numVal() const {
        std::size_t sum = 0;
        for (std::size_t jj = 0; jj < mSnodeSize; ++jj)
            sum += (mFrontSize[jj] + mUpdateSize[jj]) * mFrontSize[jj];
        return sum;
    }
    std::size_t nnz() const {
        std::size_t sum = 0;
        for (std::size_t jj = 0; jj < mSnodeSize; ++jj)
            sum += mFrontSize[jj] * (mFrontSize[jj] + 1) / 2 + mFrontSize[jj] * mUpdateSize[jj];
        return sum;
    }

    std::int32_t firstRoot() const { return mFirstRoot; }   // first root snodeIdx, or NIL
    std::int32_t lastRoot()  const { return mLastRoot; }    // last root snodeIdx, or NIL

    // Node-to-supernode map (length size()): the whole map, and the supernode a single node
    // belongs to. Indexed by a *node*, unlike the per-supernode accessors.
    const std::vector<std::int32_t>& nodeToSnode()                  const { return mNodeToSnode; }
    std::int32_t                     nodeToSnode(std::int32_t node) const { return mNodeToSnode[node]; }

    // Per-supernode links (length snodeSize()).
    const std::vector<std::int32_t>& parent()      const { return mParent; }
    const std::vector<std::int32_t>& firstChild()  const { return mFirstChild; }
    const std::vector<std::int32_t>& nextSibling() const { return mNextSibling; }

    // Per-supernode sizes (length snodeSize()).
    const std::vector<std::size_t>&  frontSize()   const { return mFrontSize; }
    const std::vector<std::size_t>&  updateSize()  const { return mUpdateSize; }
    std::size_t frontSize(std::int32_t kk)  const { return mFrontSize[kk]; }
    std::size_t updateSize(std::int32_t kk) const { return mUpdateSize[kk]; }

    // The index sets, flat: offsets (length snodeSize() + 1), then node indices
    // (length numNodeIdx()).
    const std::vector<std::size_t>&  snodePtr() const { return mSnodePtr; }
    const std::vector<std::int32_t>& nodeIdx()         const { return mNodeIdx; }

private:
    // Dimensions and tree attributes, copied from the forest.
    std::size_t  mSize      = 0;     // number of columns
    std::size_t  mSnodeSize = 0;     // number of supernodes
    std::size_t  mNumTrees  = 0;     // number of trees (roots)
    std::size_t  mHeight    = 0;     // forest height (max depth + 1)
    std::int32_t mFirstRoot = NIL;   // first root snodeIdx, or NIL if empty
    std::int32_t mLastRoot  = NIL;   // last root snodeIdx, or NIL if empty

    // Column-indexed (length mSize), copied from the forest.
    std::vector<std::int32_t> mNodeToSnode;   // node -> snode

    // Per-supernode links (length mSnodeSize), copied from the forest. Singly linked
    // here: the child list is walked forward, which needs only these three.
    std::vector<std::int32_t> mParent;        // parent snodeIdx, or NIL at a root
    std::vector<std::int32_t> mFirstChild;    // first child snodeIdx, or NIL at a leaf
    std::vector<std::int32_t> mNextSibling;   // next sibling snodeIdx, or NIL at last

    // Per-supernode sizes (length mSnodeSize), copied from the forest. Together they
    // give each supernode's index count, and frontSize splits its index set.
    std::vector<std::size_t>  mFrontSize;     // front indices (columns) in the supernode
    std::vector<std::size_t>  mUpdateSize;    // update indices below the supernode

    // The index sets, flat and computed here.
    std::size_t               mNumNodeIdx = 0;   // total node indices (== mSnodePtr[mSnodeSize])
    std::vector<std::size_t>  mSnodePtr;   // offsets into mNodeIdx (length mSnodeSize + 1)
    std::vector<std::int32_t> mNodeIdx;           // factor node indices (length mNumNodeIdx)

    friend class SymFactorEngine;   // fills the symbolic factorization via the engine
};

} // namespace Oblio
