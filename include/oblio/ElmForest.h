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

    // Do the columns of every supernode share exactly one sparsity pattern?
    //
    // True for a nodal forest (a supernode is a column, trivially) and for fundamental
    // supernodes (that is condition 2 of their definition). True also after amalgamation at
    // threshold zero, which merges only where the merge is free, hence only where the patterns
    // already agree. False as soon as amalgamation stores an explicit zero: the merged columns
    // then have *nearly* identical patterns, and the later ones carry rows the first does not.
    //
    // Consumers need this. Symbolic factorization can read one front column per supernode when
    // it holds, and must read them all when it does not.
    bool exactPatterns() const { return mExactPatterns; }

    std::size_t size()      const { return mSize; }        // number of columns
    std::size_t snodeSize() const { return mSnodeSize; }   // number of supernodes
    std::size_t numTrees()  const { return mNumTrees; }    // number of trees (roots)
    std::size_t height()    const { return mHeight; }      // forest height (max depth + 1)

    std::int32_t firstRoot() const { return mFirstRoot; }   // first root snodeIdx, or NIL
    std::int32_t lastRoot()  const { return mLastRoot; }    // last root snodeIdx, or NIL

    // Node-to-supernode map (length size()).
    const std::vector<std::int32_t>& nodeToSnode() const { return mNodeToSnode; }

    // Per-supernode links (length snodeSize()), doubly linked.
    const std::vector<std::int32_t>& parent()          const { return mParent; }
    const std::vector<std::int32_t>& firstChild()      const { return mFirstChild; }
    const std::vector<std::int32_t>& lastChild()       const { return mLastChild; }
    const std::vector<std::int32_t>& nextSibling()     const { return mNextSibling; }
    const std::vector<std::int32_t>& previousSibling() const { return mPreviousSibling; }

    // Per-supernode sizes (length snodeSize()).
    const std::vector<std::size_t>&  frontSize()   const { return mFrontSize; }
    const std::vector<std::size_t>&  updateSize()  const { return mUpdateSize; }

private:
    // Dimensions and tree attributes.
    std::size_t  mSize      = 0;     // number of columns
    std::size_t  mSnodeSize = 0;     // number of supernodes (== mSize while trivial)
    std::size_t  mNumTrees  = 0;     // number of trees (roots)
    std::size_t  mHeight    = 0;     // forest height (max depth + 1)
    std::int32_t mFirstRoot = NIL;   // first root snodeIdx, or NIL if empty
    std::int32_t mLastRoot  = NIL;   // last root snodeIdx, or NIL if empty

    // Set false only by amalgamation, and only when it actually stores a zero.
    bool mExactPatterns = true;

    // Column-indexed (length mSize): the supernode owning each column.
    // Identity while supernodes are trivial (snodeIdx == column index).
    std::vector<std::int32_t> mNodeToSnode;   // node -> snode

    // Per-supernode links (length mSnodeSize), doubly linked. IDs, NIL terminator.
    std::vector<std::int32_t> mParent;           // parent snodeIdx, or NIL at a root
    std::vector<std::int32_t> mFirstChild;       // first child snodeIdx, or NIL at a leaf
    std::vector<std::int32_t> mLastChild;        // last child snodeIdx, or NIL at a leaf
    std::vector<std::int32_t> mNextSibling;      // next sibling snodeIdx, or NIL at last
    std::vector<std::int32_t> mPreviousSibling;  // previous sibling snodeIdx, or NIL at first

    // Per-supernode sizes (length mSnodeSize). Counts, so std::size_t.
    std::vector<std::size_t>  mFrontSize;     // front indices (columns) in the supernode
    std::vector<std::size_t>  mUpdateSize;    // update indices below the supernode

    friend class ElmForestEngine;   // fills the forest via the engine
};

} // namespace Oblio
