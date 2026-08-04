#pragma once

// ElmForestEngine.h — computes the elimination forest of a SparseMatrix under a
// given Permutation. Builds the forest's attributes: the parent links (etree),
// the (trivial) supernode structure, the doubly-linked child/sibling links with
// tree roots and height, and the front/update sizes. Returns true on success.
// Reads the matrix/permutation via their public accessors and fills the forest
// via friend access.
//
// A is stored full-symmetric, so the etree reads each column's neighbors
// directly (the diagonal self-skips) — no expansion step.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/ElmForest.h"

#include <optional>
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
    // The three settings in pipeline order: which supernodes to form, whether to amalgamate them
    // and how hard, and whether to order the children for the multifrontal stack. Each is also a
    // setter below; the constructor exists so a caller that knows all three up front, as
    // DirectSolver does, can say so in one place and hold the engine const.
    explicit ElmForestEngine(Supernodes                 supernodes,
                             std::optional<std::size_t> threshold            = std::nullopt,
                             bool                       optimizeMultifrontal = false)
        : mSupernodes(supernodes), mThreshold(threshold),
          mOptimizeMultifrontal(optimizeMultifrontal) {}

    void       setSupernodes(Supernodes supernodes) { mSupernodes = supernodes; }
    Supernodes supernodes() const                   { return mSupernodes; }

    // The amalgamation threshold: how many explicitly stored zeros we are willing to buy a
    // wider front with. Amalgamation merges a supernode with several of its children, unlike
    // fundamental compression, which merges only an only-child, and unlike it, it may pay
    // fill to do so. It runs after fundamental compression when both are on.
    //
    // Absent means no amalgamation at all, which is not the same as a threshold of zero: at
    // zero it still merges, but only where the merge is free. See the DESIGN_DECISIONS entry
    // and Section 4.5 of the sparse-factorization notes.
    void setThreshold(std::optional<std::size_t> threshold) { mThreshold = threshold; }
    std::optional<std::size_t> threshold() const            { return mThreshold; }

    // Reorder each supernode's children so that assembling them first to last minimizes the
    // multifrontal update stack's peak. 0.9's `optimizeMultifrontalCore_`, and off by default
    // there as here, because it is not free of consequence:
    //
    // It changes the child order for *every* traversal, not only multifrontal, the forest being
    // shared. Left- and right-looking do not care for the static factorizations, where the
    // assembly is a commutative sum, but dynamic LDL does: a parent's front is expanded by
    // prepending its children's delayed columns in sibling order, so reordering the siblings
    // reorders those columns, which the pivot sequence then reads. The result is a different and
    // equally valid factorization with different delay and pivot counts, not wrong, just different.
    //
    // Setting it runs two steps, not one, as 0.9 does. `sortForOptimalMultifrontal` chooses the
    // child order, and `labelDepthFirst` then relabels the supernodes into a postorder so that the
    // labels follow that order. The second is what realizes the saving: the multifrontal drivers
    // loop over supernodes in increasing label order, so a contribution block lives from its own
    // supernode to its parent *in the numbering*, which the child links alone do not affect. The
    // sort without the relabeling would be inert.
    //
    // **Whoever drives the engines by hand owns this decision.** The forest is built long before a
    // traversal is chosen, so this object cannot know that multifrontal is coming; a caller running
    // OrderEngine, ElmForestEngine and SymFactorEngine themselves must set it when they intend
    // multifrontal. `DirectSolver` does it for them, passing `traversal() == Multifrontal` to the
    // constructor, which is the only place both facts are known at once.
    //
    // Forgetting it is silent: the factorization is correct and the residual unchanged, only the
    // update stack's peak is larger, measured at 16 to 38 percent on grids. Nothing detects that,
    // which is why it is written down here. See TODO for the measurements.
    void setOptimizeMultifrontal(bool optimize) { mOptimizeMultifrontal = optimize; }
    bool optimizeMultifrontal() const           { return mOptimizeMultifrontal; }

    // Compute the elimination forest of A under the permutation P.
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
    bool compute(const SparseMatrix<Val>& A, const Permutation& P, ElmForest& ef) const;

    bool compute(const std::vector<std::size_t>&  colPtr,
                 const std::vector<std::int32_t>& rowIdx,
                 const Permutation& P, ElmForest& ef) const;

private:
    Supernodes                 mSupernodes = Supernodes::Fundamental;
    std::optional<std::size_t> mThreshold;   // absent: no amalgamation
    bool                       mOptimizeMultifrontal = false;   // reorder children for the stack
    // Every helper below is free of Val. The adaptation happens once, at the public
    // boundary above, so nothing in here is templated and nothing in here is compiled twice.
    // The matrix appears only as its sparsity pattern, colPtr and rowIdx. The Permutation
    // and the ElmForest are passed whole: neither is templated, the Permutation is read
    // through its public accessors, and the ElmForest is written through friendship.

    // Parent links (etree) via Liu's path-compression construction, on the full-symmetric
    // structure, in the permuted (factor) order. Writes ef.mSize and ef.mParent.
    void computeParent(const std::vector<std::size_t>&  colPtr,
                       const std::vector<std::int32_t>& rowIdx,
                       const Permutation& P,
                       ElmForest& ef) const;

    // Given the parent links, finalize the rest: the doubly-linked child/sibling structure
    // (first/last child, next/previous sibling) plus the tree count and the first/last
    // root. Front-insertion in decreasing order yields children and roots in increasing
    // label order.
    //
    //   reads:   mSnodeSize, mParent
    //   writes:  mFirstChild, mLastChild, mNextSibling, mPreviousSibling,
    //            mNumTrees, mFirstRoot, mLastRoot
    void finalizeLinks(ElmForest& ef) const;

    // Forest height by a roots-down breadth-first traversal: depth[s] =
    // depth[parent[s]] + 1, height = max depth + 1. Reads the doubly-linked child/sibling
    // structure, so it runs after the links are built. Returns the height rather than
    // storing it, so it stays a plain function of the forest, and takes it by const
    // reference to say so.
    std::size_t computeHeight(const ElmForest& ef) const;

    // The column sizes of a nodal forest: how many nonzeros each column of L holds.
    //
    // Two parts, of very different weight. The front size is trivial, a supernode is a column
    // here and owns exactly itself, so it is 1. The update size is the real work: the
    // subdiagonal nonzero count of each column, computed by the pruned-row-subtree walk
    // without ever building the pattern it counts. We need the counts to size the symbolic
    // factorization's storage before we know what goes in it.
    //
    // Runs before compression and only before it: the walk indexes by column, which coincides
    // with supernode only while the forest is nodal. compressFundamental derives the merged
    // sizes itself, so this is never re-run afterwards.
    //
    // Row-oriented, like computeParent and unlike everything downstream, and for the same
    // reason it reads A's columns above the diagonal rather than its rows: under full storage
    // those are the same entries.
    //
    // Nodal throughout: it reads and writes nothing supernodal, and sizes its outputs by the
    // column count, not the supernode count. The two coincide here, but leaning on that would
    // be leaning on a coincidence.
    //
    //   reads:   mSize, mParent
    //   writes:  mFrontSize, mUpdateSize
    void computeColumnSizes(const std::vector<std::size_t>&  colPtr,
                            const std::vector<std::int32_t>& rowIdx,
                            const Permutation& P,
                            ElmForest& ef) const;

    // Merge the columns of a nodal forest into fundamental supernodes: maximal paths in
    // the forest whose columns share one sparsity pattern, each vertex on the path (bar
    // the bottom one) having exactly one child. Fundamental supernodes are unique for a
    // given forest.
    //
    // Precondition: ef is nodal, one column per supernode. That is the only useful input,
    // since fundamental supernodes are maximal and a second pass would find nothing to
    // merge.
    //
    // Takes the forest rather than its arrays. The engine is a friend, so passing pieces
    // restricts nothing, and the contract here is not one a signature can carry anyway:
    // this rewrites the forest's shape, and what matters is which parts it leaves valid.
    //
    //   reads:   mParent, mFirstChild, mLastChild, mFrontSize, mUpdateSize
    //   writes:  mSnodeSize, mNodeToSnode, mParent, mFrontSize, mUpdateSize
    //   stales:  mFirstChild, mLastChild, mNextSibling, mPreviousSibling,
    //            mNumTrees, mFirstRoot, mLastRoot   (caller rebuilds with finalizeLinks)
    //   leaves:  mSize, mHeight   (height is computed once, afterwards)
    //
    // The stale set is the point. A const& parameter would have said "not written", which
    // reads as "still valid", and these links are precisely not that: they describe the
    // nodal forest and are wrong for the compressed one until rebuilt.
    void compressFundamental(ElmForest& ef) const;

    // Amalgamate: merge a supernode with several of its children, paying explicitly stored
    // zeros for a wider front. The threshold is the budget, in zeros.
    //
    // Unlike compressFundamental this works on any forest, nodal or supernodal, since a nodal
    // forest is one whose supernodes are all trivial. Unlike it, it is also **not canonical**:
    // where two children of a supernode could each be merged for free, only one can be (the
    // first widens the front, which prices out the second), so the algorithm must break ties,
    // and the tie-breaking rule is a convention rather than a theorem. Fundamental supernodes
    // are unique precisely because they refuse to make that choice.
    //
    // Even at threshold zero this merges more than compressFundamental does: a free merge
    // asks only that the child's pattern match, while a fundamental merge additionally
    // demands that the child be an only child.
    //
    //   reads:   mSize, mSnodeSize, mNodeToSnode, mParent, mFirstChild, mNextSibling,
    //            mFrontSize, mUpdateSize
    //   writes:  mSnodeSize, mNodeToSnode, mParent, mFrontSize, mUpdateSize
    //   stales:  the child, sibling and root links, as above
    void compressThreshold(ElmForest& ef, std::size_t threshold) const;

    // The two steps of the multifrontal ordering, run as a pair. They are written apart because
    // they are separate ideas, one choosing an order among siblings and the other making the
    // labels carry it, but they are not independently useful and this is their only caller.
    //
    // The asymmetry is worth knowing. `labelDepthFirst` alone is harmless: it produces a valid
    // postorder, the factor is unchanged, and it simply buys nothing when no sort has chosen an
    // order for it to follow. `sortForOptimalMultifrontal` alone is worse than useless, because it
    // is *inert while looking like it worked*: the links carry the chosen order, the drivers walk
    // labels and never read them, and the stack peak does not move. Measured, on a 3D grid: sort
    // alone leaves the peak at 357398 while the model it optimizes reports 255433.
    //
    //   reads / writes:  whatever the two below do
    void orderForMultifrontal(ElmForest& ef) const;

    // Reorder each supernode's children into the order that minimizes the multifrontal stack
    // peak. Ported from 0.9 `sortForOptimalMultifrontal_`. Never call this without the relabeling
    // below; call `orderForMultifrontal` instead.
    // Runs after both compressions, on the final supernodes, and rewrites the links in place,
    // so unlike the compressions it leaves nothing stale and needs no finalizeLinks after it.
    //
    //   reads:   mSnodeSize, mFirstChild, mNextSibling, mUpdateSize
    //   writes:  mFirstChild, mLastChild, mNextSibling, mPreviousSibling
    void sortForOptimalMultifrontal(ElmForest& ef) const;

    // Relabel the supernodes into a postorder, so that each subtree holds a contiguous run of
    // labels. Ported from 0.9 `labelDepthFirst_`, which 0.9 calls straight after the sort above.
    // Safe on its own, though pointless on its own: see `orderForMultifrontal`.
    // It is the half that realizes the saving: the drivers loop in increasing label order, so the
    // stack peak follows the labels, not the links, and only a postorder makes the two agree.
    //
    //   reads:   mSize, mSnodeSize, mNodeToSnode, mFirstRoot, mFirstChild, mNextSibling
    //   writes:  mNodeToSnode, mParent, mFirstChild, mLastChild, mNextSibling, mPreviousSibling,
    //            mFrontSize, mUpdateSize, mFirstRoot, mLastRoot
    void labelDepthFirst(ElmForest& ef) const;
};

extern template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
extern template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
