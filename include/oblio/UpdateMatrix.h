#pragma once

// UpdateMatrix.h - one supernode's entire Schur complement, for the multifrontal traversal.
//
// When supernode J is factored, eliminating its pivots leaves a dense Schur complement over J's
// update rows. In left- and right-looking that complement is never materialized whole: it is
// formed one ancestor at a time as an UpdateBlock and scattered immediately. Multifrontal instead
// forms it *once*, in full, and carries it up the tree: J's complement is assembled into J's
// parent, whose own complement is assembled into *its* parent, and so on. This object is that
// full complement.
//
// The geometry. Let u be J's update size. The matrix is u by u, square, symmetric (only the lower
// triangle is meaningful), column-major, paired with the global (factor-ordered) node index of
// each of its u rows-and-columns:
//
//     size     = u, J's update size; the block is size by size
//     nodeIdx  = the u global indices, one per row/column
//     val      = size * size values, column-major, lower triangle meaningful
//
// **UpdateBlock versus UpdateMatrix**, the distinction UpdateBlock.h draws from the other side:
//
//   UpdateBlock     J's contribution to ONE ancestor K. Rectangular, one per (J, K) pair, formed
//                   and scattered and discarded at once. What left- and right-looking use.
//   UpdateMatrix    J's ENTIRE Schur complement, u by u, symmetric, one per supernode, assembled
//                   into J's parent. What multifrontal uses.
//
// A slice versus the whole; per-pair versus per-supernode; immediate versus carried.
//
// **The container.** Multifrontal keeps one UpdateMatrix per supernode in a flat
// `std::vector<UpdateMatrix<Val>>`, named `updateMatrix` and sized once to the supernode count. A
// supernode's matrix is `allocate`d when the supernode is reached and `discard`ed once its parent
// has assembled it, and because everything in multifrontal flows child to parent, that lifetime is
// bounded. `discard` releases the storage rather than clearing it, so the peak is the storage live
// at one cut across the tree, not the sum over the whole run.
//
// **It is not called a stack, deliberately.** The literature calls this the update stack and the
// lifetimes are indeed nested, but nothing here is pushed or popped: the slots are independent
// allocations indexed by supernode, freed in child order rather than from a top, and a parent's own
// matrix is allocated *before* its children are freed, so it sits above blocks that then die
// beneath it, which no stack pointer can express. The name would claim a discipline the code does
// not keep. The abstract UpdateStack / concrete UpdateStackDynamic split in 0.9 existed to carry an
// out-of-core hook and storage-size counters; in core, the bare vector suffices and neither is
// needed. See the TODO entry on the stack for what a true LIFO would and would not buy.
//
// allocate and discard are the engine's, through friendship, so each slot's lifetime stays under
// the traversal's control and no other caller can grow or free one.

#include "oblio/Types.h"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

class NumFactorEngine;

template<class Val>
class UpdateMatrix {
public:
    UpdateMatrix() = default;   // empty; the engine allocates it when the supernode is reached

    std::size_t         size()    const { return mSize; }   // rows == cols == the supernode's update size

    // The row indices and the values, read-only. Reads go through these even from the engine,
    // which is a friend: friendship is for writing, not for reaching past the interface to read.
    const std::int32_t* rowIdx()  const { return mRowIdx.data(); }
    const Val*          val()     const { return mVal.data(); }

    // The block's leading dimension, for BLAS. Column-major square, so it is the size.
    std::size_t ld() const { return mSize; }

    // Never copied (a contribution block is heavy and single-owner); movable so the backing
    // vector can be built and, if ever resized, relocate its slots without a deep copy.
    UpdateMatrix(const UpdateMatrix&)            = delete;
    UpdateMatrix& operator=(const UpdateMatrix&) = delete;
    UpdateMatrix(UpdateMatrix&&)                 = default;
    UpdateMatrix& operator=(UpdateMatrix&&)      = default;

private:
    // Size to a `size`-by-`size` block: its index set and its values, the latter zeroed because
    // assembly adds into it. The engine fills the indices afterward through the non-const accessor.
    void allocate(std::size_t size);

    // Release the storage, back to the empty state, once the parent has assembled this block. This
    // frees rather than clears, which is what keeps the peak bounded.
    void discard();

    // The writable views, private, so only the engine reaches them through friendship.
    std::int32_t* rowIdx() { return mRowIdx.data(); }
    Val*          val()    { return mVal.data(); }

    std::size_t               mSize = 0;
    std::vector<std::int32_t> mRowIdx;    // the `size` global row indices, as in UpdateBlock
    std::vector<Val>          mVal;       // size * size, column-major, lower triangle meaningful

    friend class NumFactorEngine;
};

extern template class UpdateMatrix<double>;
extern template class UpdateMatrix<std::complex<double>>;

} // namespace Oblio
