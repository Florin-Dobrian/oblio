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
    // The forest holds no index set and no values, so all three are predictions here, of what the
    // symbolic factorization and then the numeric one will need. That is the point: they are
    // available a phase before anything is allocated.
    std::size_t numNodeIdx() const {
        std::size_t sum = 0;
        for (std::size_t jj = 0; jj < mSnodeSize; ++jj) sum += mFrontSize[jj] + mUpdateSize[jj];
        return sum;
    }
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

    // Per-supernode links (length snodeSize()), doubly linked.
    const std::vector<std::int32_t>& parent()          const { return mParent; }
    const std::vector<std::int32_t>& firstChild()      const { return mFirstChild; }
    const std::vector<std::int32_t>& lastChild()       const { return mLastChild; }
    const std::vector<std::int32_t>& nextSibling()     const { return mNextSibling; }
    const std::vector<std::int32_t>& previousSibling() const { return mPreviousSibling; }

    // Per-supernode sizes (length snodeSize()).
    const std::vector<std::size_t>&  frontSize()   const { return mFrontSize; }
    const std::vector<std::size_t>&  updateSize()  const { return mUpdateSize; }
    std::size_t frontSize(std::int32_t kk)  const { return mFrontSize[kk]; }
    std::size_t updateSize(std::int32_t kk) const { return mUpdateSize[kk]; }

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
