#pragma once

// ElmForest.h — the elimination forest of a factored sparse matrix.
//
// "Forest" (not tree) because a reducible/disconnected matrix yields several
// trees. A vertex of the forest is a supernode: a set of factor columns in the
// permuted (factor) order. Supernodes are trivial for now (one column each, so a
// supernode index equals its column index), but the structure is supernodal
// throughout, so it keeps its shape when fundamental supernodes are added later.
// Symbolic factorization consumes this structure.
//
// The forest's attributes: its dimensions (column count, supernode count) and tree
// attributes (number of trees, height, first/last root); the column -> supernode
// map; the doubly-linked per-supernode structure (parent, first/last child,
// next/previous sibling), where the root list is itself a sibling chain; and the
// per-supernode front and update sizes. These are structural attributes of the
// forest, not optional statistics. Only fundamental-supernode compression, an
// engine step that merges columns, is deferred.
//
// Index roles (see the index-types decision): supernode indices, column indices,
// the links, and the roots are IDs -> std::int32_t with NIL = -1; sizes and
// counts are std::size_t. While supernodes are trivial, parent[s] > s still holds
// (a parent is a later factor column).

#include "oblio/Types.h"

#include <vector>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class ElmForestEngine;

class ElmForest {
public:
    ElmForest() = default;

    std::size_t size()     const { return mSize; }       // number of columns
    std::size_t supSize()  const { return mSupSize; }    // number of supernodes
    std::size_t numTrees() const { return mNumTrees; }   // number of trees (roots)
    std::size_t height()   const { return mHeight; }     // forest height (max depth + 1)

    std::int32_t firstRoot() const { return mFirstRoot; }   // first root supIdx, or NIL
    std::int32_t lastRoot()  const { return mLastRoot; }    // last root supIdx, or NIL

    // Column-indexed map (length size()).
    const std::vector<std::int32_t>& idxToSupIdx() const { return mIdxToSupIdx; }

    // Per-supernode links (length supSize()), doubly linked.
    const std::vector<std::int32_t>& parent()          const { return mParent; }
    const std::vector<std::int32_t>& firstChild()      const { return mFirstChild; }
    const std::vector<std::int32_t>& lastChild()       const { return mLastChild; }
    const std::vector<std::int32_t>& nextSibling()     const { return mNextSibling; }
    const std::vector<std::int32_t>& previousSibling() const { return mPreviousSibling; }

    // Per-supernode sizes (length supSize()).
    const std::vector<std::size_t>&  frontSize()   const { return mFrontSize; }
    const std::vector<std::size_t>&  updateSize()  const { return mUpdateSize; }

private:
    // Dimensions and tree attributes.
    std::size_t  mSize      = 0;     // number of columns
    std::size_t  mSupSize   = 0;     // number of supernodes (== mSize while trivial)
    std::size_t  mNumTrees  = 0;     // number of trees (roots)
    std::size_t  mHeight    = 0;     // forest height (max depth + 1)
    std::int32_t mFirstRoot = NIL;   // first root supIdx, or NIL if empty
    std::int32_t mLastRoot  = NIL;   // last root supIdx, or NIL if empty

    // Column-indexed (length mSize): the supernode owning each column.
    // Identity while supernodes are trivial (supIdx == column index).
    std::vector<std::int32_t> mIdxToSupIdx;   // idx -> supIdx

    // Per-supernode links (length mSupSize), doubly linked. IDs, NIL terminator.
    std::vector<std::int32_t> mParent;           // parent supIdx, or NIL at a root
    std::vector<std::int32_t> mFirstChild;       // first child supIdx, or NIL at a leaf
    std::vector<std::int32_t> mLastChild;        // last child supIdx, or NIL at a leaf
    std::vector<std::int32_t> mNextSibling;      // next sibling supIdx, or NIL at last
    std::vector<std::int32_t> mPreviousSibling;  // previous sibling supIdx, or NIL at first

    // Per-supernode sizes (length mSupSize). Counts, so std::size_t.
    std::vector<std::size_t>  mFrontSize;     // front indices (columns) in the supernode
    std::vector<std::size_t>  mUpdateSize;    // update indices below the supernode

    friend class ElmForestEngine;   // fills the forest via the engine
};

} // namespace Oblio
