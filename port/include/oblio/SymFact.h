#pragma once

// SymFact.h - the symbolic factorization of a sparse matrix.
//
// The symbolic factorization is the sparsity pattern of the factor L, computed
// per supernode: for every supernode, the set of factor row indices its columns
// touch. That index set splits into the front indices (the supernode's own
// columns) and the update indices (the rows below it, which its children and it
// itself update). Numeric factorization consumes this structure.
//
// SymFact is self-contained: it copies the forest attributes it needs (the map,
// the parent/child/sibling links, the front and update sizes, the tree
// attributes), so numeric factorization needs SymFact alone, not the forest.
// Both 0.9 and 10.12 do this.
//
// Storage is flat. The index sets are written once and their sizes are known in
// advance from the forest (front size + update size per supernode), so they go in
// one contiguous array mIdx with per-supernode offsets mIdxPtr, exactly as in 0.9
// (10.12 uses a vector of vectors here; see the flat-vs-VV decision). For
// supernode s, its index set is mIdx[mIdxPtr[s] .. mIdxPtr[s + 1]), whose first
// frontSize(s) entries are the front indices and the rest the update indices.
// Within a supernode the indices are sorted in increasing order.
//
// Index roles (see the index-types decision): supernode indices, column indices,
// factor row indices and the links are IDs -> std::int32_t with NIL = -1; sizes,
// counts and offsets are std::size_t.

#include "oblio/Types.h"

#include <vector>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class SymFactEngine;

class SymFact {
public:
    SymFact() = default;

    std::size_t size()     const { return mSize; }       // number of columns
    std::size_t supSize()  const { return mSupSize; }    // number of supernodes
    std::size_t numTrees() const { return mNumTrees; }   // number of trees (roots)
    std::size_t height()   const { return mHeight; }     // forest height (max depth + 1)
    std::size_t numIdx()   const { return mNumIdx; }     // total factor indices

    std::int32_t firstRoot() const { return mFirstRoot; }   // first root supIdx, or NIL
    std::int32_t lastRoot()  const { return mLastRoot; }    // last root supIdx, or NIL

    // Column-indexed map (length size()).
    const std::vector<std::int32_t>& idxToSupIdx() const { return mIdxToSupIdx; }

    // Per-supernode links (length supSize()).
    const std::vector<std::int32_t>& parent()      const { return mParent; }
    const std::vector<std::int32_t>& firstChild()  const { return mFirstChild; }
    const std::vector<std::int32_t>& nextSibling() const { return mNextSibling; }

    // Per-supernode sizes (length supSize()).
    const std::vector<std::size_t>&  frontSize()   const { return mFrontSize; }
    const std::vector<std::size_t>&  updateSize()  const { return mUpdateSize; }

    // The index sets, flat. Offsets (length supSize() + 1), then indices
    // (length numIdx()).
    const std::vector<std::size_t>&  idxPtr() const { return mIdxPtr; }
    const std::vector<std::int32_t>& idx()    const { return mIdx; }

private:
    // Dimensions and tree attributes, copied from the forest.
    std::size_t  mSize      = 0;     // number of columns
    std::size_t  mSupSize   = 0;     // number of supernodes
    std::size_t  mNumTrees  = 0;     // number of trees (roots)
    std::size_t  mHeight    = 0;     // forest height (max depth + 1)
    std::int32_t mFirstRoot = NIL;   // first root supIdx, or NIL if empty
    std::int32_t mLastRoot  = NIL;   // last root supIdx, or NIL if empty

    // Column-indexed (length mSize), copied from the forest.
    std::vector<std::int32_t> mIdxToSupIdx;   // idx -> supIdx

    // Per-supernode links (length mSupSize), copied from the forest. Singly linked
    // here: the child list is walked forward, which needs only these three.
    std::vector<std::int32_t> mParent;        // parent supIdx, or NIL at a root
    std::vector<std::int32_t> mFirstChild;    // first child supIdx, or NIL at a leaf
    std::vector<std::int32_t> mNextSibling;   // next sibling supIdx, or NIL at last

    // Per-supernode sizes (length mSupSize), copied from the forest. Together they
    // give each supernode's index count, and frontSize splits its index set.
    std::vector<std::size_t>  mFrontSize;     // front indices (columns) in the supernode
    std::vector<std::size_t>  mUpdateSize;    // update indices below the supernode

    // The index sets, flat and computed here.
    std::size_t               mNumIdx = 0;    // total factor indices (== mIdxPtr[mSupSize])
    std::vector<std::size_t>  mIdxPtr;        // offsets into mIdx (length mSupSize + 1)
    std::vector<std::int32_t> mIdx;           // factor row indices (length mNumIdx)

    friend class SymFactEngine;   // fills the symbolic factorization via the engine
};

} // namespace Oblio
