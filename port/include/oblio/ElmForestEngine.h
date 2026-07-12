#pragma once

// ElmForestEngine.h — computes the elimination forest of a SparseMatrix under a
// given Permutation. Builds the forest's attributes: the parent links (etree),
// the (trivial) supernode structure, the doubly-linked child/sibling links with
// tree roots and height, and the front/update sizes. Returns true on success.
// Reads the matrix/permutation via their public accessors and fills the forest
// via friend access.
//
// A is stored full-symmetric, so the etree reads each column's neighbours
// directly (the diagonal self-skips) — no expansion step.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/ElmForest.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {

// Whether to merge the trivial (one column) supernodes into fundamental supernodes.
// Nodal keeps one supernode per column, which is the regime every consumer degenerates
// to gracefully, and is useful for isolating whether a fault lies in compression or
// downstream of it. Both references carry the same switch (0.9's supernodal_, 10.12's
// mCompressFundamental), defaulting to on.
//
// Threshold-based compression, when it lands, is an orthogonal setting rather than a
// third value here: both references apply it on top of either regime.
enum class Supernodes { Nodal, Fundamental };

class ElmForestEngine {
public:
    ElmForestEngine() = default;
    explicit ElmForestEngine(Supernodes supernodes) : mSupernodes(supernodes) {}

    void       setSupernodes(Supernodes supernodes) { mSupernodes = supernodes; }
    Supernodes supernodes() const                   { return mSupernodes; }

    // Compute the elimination forest of A under the permutation p.
    //
    // Two overloads, the same layering used throughout this engine. The forest depends on
    // the sparsity pattern alone, never on the values, so the implementation is the
    // non-templated one: it takes colPtr and rowIdx and is compiled once. The templated
    // overload is an adapter over it, for the common case of holding a matrix.
    //
    // The pattern overload is public, not merely an internal convenience: a caller holding
    // a graph with no numbers attached can compute its elimination forest without inventing
    // a scalar type to satisfy the signature.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p, ElmForest& f) const;

    bool compute(const std::vector<std::size_t>&  colPtr,
                 const std::vector<std::int32_t>& rowIdx,
                 const Permutation& p, ElmForest& f) const;

private:
    Supernodes mSupernodes = Supernodes::Fundamental;
    // Every helper below is free of Val. The adaptation happens once, at the public
    // boundary above, so nothing in here is templated and nothing in here is compiled twice.
    // The matrix appears only as its sparsity pattern, colPtr and rowIdx. The Permutation
    // and the ElmForest are passed whole: neither is templated, the Permutation is read
    // through its public accessors, and the ElmForest is written through friendship.

    // Parent links (etree) via Liu's path-compression construction, on the full-symmetric
    // structure, in the permuted (factor) order. Writes f.mSize and f.mParent.
    void computeParent(const std::vector<std::size_t>&  colPtr,
                       const std::vector<std::int32_t>& rowIdx,
                       const Permutation& p,
                       ElmForest& f) const;

    // Given the parent links, finalize the rest: the doubly-linked child/sibling structure
    // (first/last child, next/previous sibling) plus the tree count and the first/last
    // root. Front-insertion in decreasing order yields children and roots in increasing
    // label order.
    //
    //   reads:   mSupSize, mParent
    //   writes:  mFirstChild, mLastChild, mNextSibling, mPreviousSibling,
    //            mNumTrees, mFirstRoot, mLastRoot
    void finalizeLinks(ElmForest& f) const;

    // Forest height by a roots-down breadth-first traversal: depth[s] =
    // depth[parent[s]] + 1, height = max depth + 1. Reads the doubly-linked child/sibling
    // structure, so it runs after the links are built. Returns the height rather than
    // storing it, so it stays a plain function of the forest, and takes it by const
    // reference to say so.
    std::size_t computeHeight(const ElmForest& f) const;

    // Per-supernode front and update sizes, for a nodal forest. This runs before
    // compression and only before it: the update-size walk indexes by column, which
    // coincides with supernode only while supernodes are trivial. compressFundamental
    // derives the merged sizes itself, so this is never re-run afterwards.
    //
    // Front size counts the columns mapping to each supernode. On the nodal forest the map
    // is the identity, so this is all 1s and a fixed fill (10.12's choice) would do as
    // well. Counting derives the value from the map rather than asserting it, which is the
    // only reason to prefer it; it buys no generality here.
    //
    // Update size is the subdiagonal nonzero count of each L column, by the pruned-row-
    // subtree walk. This is 0.9's columnSize with the diagonal dropped, so it equals
    // columnSize - 1. Depends only on the map and the parent links, not on the
    // child/sibling structure, so it could equally run before the links are built.
    //
    //   reads:   mSize, mSupSize, mParent, mIdxToSupIdx
    //   writes:  mFrontSize, mUpdateSize
    void computeFrontAndUpdateSizes(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    const Permutation& p,
                                    ElmForest& f) const;

    // Merge the columns of a nodal forest into fundamental supernodes: maximal paths in
    // the forest whose columns share one sparsity pattern, each vertex on the path (bar
    // the bottom one) having exactly one child. Fundamental supernodes are unique for a
    // given forest.
    //
    // Precondition: f is nodal, one column per supernode. That is the only useful input,
    // since fundamental supernodes are maximal and a second pass would find nothing to
    // merge.
    //
    // Takes the forest rather than its arrays. The engine is a friend, so passing pieces
    // restricts nothing, and the contract here is not one a signature can carry anyway:
    // this rewrites the forest's shape, and what matters is which parts it leaves valid.
    //
    //   reads:   mParent, mFirstChild, mLastChild, mFrontSize, mUpdateSize
    //   writes:  mSupSize, mIdxToSupIdx, mParent, mFrontSize, mUpdateSize
    //   stales:  mFirstChild, mLastChild, mNextSibling, mPreviousSibling,
    //            mNumTrees, mFirstRoot, mLastRoot   (caller rebuilds with finalizeLinks)
    //   leaves:  mSize, mHeight   (height is computed once, afterwards)
    //
    // The stale set is the point. A const& parameter would have said "not written", which
    // reads as "still valid", and these links are precisely not that: they describe the
    // nodal forest and are wrong for the compressed one until rebuilt.
    void compressFundamental(ElmForest& f) const;
};

extern template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
extern template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
