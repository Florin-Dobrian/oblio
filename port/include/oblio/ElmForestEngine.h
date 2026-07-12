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

    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p, ElmForest& f) const;

private:
    Supernodes mSupernodes = Supernodes::Fundamental;
    // Parent links (etree) via Liu's path-compression construction, on the
    // full-symmetric structure, in the permuted (factor) order.
    void computeParent(std::size_t size,
                       const std::vector<std::size_t>&  colPtr,
                       const std::vector<std::int32_t>& rowIdx,
                       const std::vector<std::int32_t>& oldToNew,
                       const std::vector<std::int32_t>& newToOld,
                       std::vector<std::int32_t>& parent) const;

    // Given parent, finalize the remaining links: the doubly-linked child/sibling
    // structure (first/last child, next/previous sibling) plus the tree count and
    // the first/last root. Front-insertion in decreasing order yields children and
    // roots in increasing label order.
    void finalizeLinks(std::size_t supSize,
                       const std::vector<std::int32_t>& parent,
                       std::vector<std::int32_t>& firstChild,
                       std::vector<std::int32_t>& lastChild,
                       std::vector<std::int32_t>& nextSibling,
                       std::vector<std::int32_t>& previousSibling,
                       std::size_t&  numTrees,
                       std::int32_t& firstRoot,
                       std::int32_t& lastRoot) const;

    // Forest height by a roots-down breadth-first traversal: depth[s] =
    // depth[parent[s]] + 1, height = max depth + 1. Reads the doubly-linked
    // child/sibling structure, so it runs after the links are built.
    std::size_t computeHeight(std::size_t supSize, std::int32_t lastRoot,
                              const std::vector<std::int32_t>& parent,
                              const std::vector<std::int32_t>& lastChild,
                              const std::vector<std::int32_t>& previousSibling) const;

    // Per-supernode front and update sizes. Front size is the number of columns a
    // supernode owns, counted from the map (all 1s while trivial, and correct once
    // compression merges columns, unlike a fixed fill). Update size is the
    // subdiagonal nonzero count of each L column, by the pruned-row-subtree walk;
    // this is 0.9's columnSize with the diagonal dropped, so == columnSize - 1.
    // Depends only on the map and the parent links, not on the child/sibling
    // structure, so it could equally run before the links are built.
    void computeFrontAndUpdateSizes(std::size_t size, std::size_t supSize,
                                    const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    const std::vector<std::int32_t>& oldToNew,
                                    const std::vector<std::int32_t>& newToOld,
                                    const std::vector<std::int32_t>& parent,
                                    const std::vector<std::int32_t>& idxToSupIdx,
                                    std::vector<std::size_t>& frontSize,
                                    std::vector<std::size_t>& updateSize) const;

    // Merge the trivial supernodes into fundamental supernodes: maximal paths in the
    // forest whose columns share one sparsity pattern, each vertex on the path (bar the
    // bottom one) having exactly one child. Fundamental supernodes are unique for a given
    // forest. Rewrites the map, the parent links and the sizes, and shrinks supSize.
    //
    // Reads the child links (the merge test asks whether a supernode has exactly one
    // child) and the sizes, so it runs after both. It leaves the child, sibling and root
    // links describing the old forest, so the caller rebuilds them with finalizeLinks.
    void compressFundamental(std::size_t size, std::size_t& supSize,
                             std::vector<std::int32_t>& idxToSupIdx,
                             std::vector<std::int32_t>& parent,
                             const std::vector<std::int32_t>& firstChild,
                             const std::vector<std::int32_t>& lastChild,
                             std::vector<std::size_t>&  frontSize,
                             std::vector<std::size_t>&  updateSize) const;
};

extern template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
extern template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
