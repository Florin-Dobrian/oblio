#pragma once

// UpdateBlock.h - the update one supernode contributes to one ancestor.
//
// When supernode J is factored, its columns update every supernode above it that shares rows with
// it. That contribution is a dense block, and this is it: `height` rows by `width` columns, in
// column-major order, plus the global row index each of its rows stands for.
//
// The geometry. Take J's index set, and the run of update rows that lie in one particular
// ancestor K. Those rows are the block's **columns** (width). Every row of J's index set from
// there down is the block's **rows** (height). So the block is the lower-trapezoidal slice of
// J's Schur complement that lands in K, held as a rectangle:
//
//     width  = how many of J's update rows belong to K
//     height = how many of J's update rows lie at or below the first of them
//     idx[t] = the global (factor-ordered) row index of local row t
//
// It is formed, scattered into K, and discarded. Nothing keeps it.
//
// **Why not "Temporary"** (0.9's name). That names the *lifetime*, which is the least interesting
// thing about the object and is in any case what a C++ object already manages. UpdateBlock names
// what it holds, and lines up with the vocabulary the rest of the solver already uses: update
// indices, update size, "a child's update indices carry into the parent". This is the numeric
// counterpart of exactly those.
//
// **And it is not the multifrontal update matrix**, which is a different object and will get a
// different name. The distinction matters:
//
//   UpdateBlock     J's contribution to ONE ancestor K. Rectangular. One per (J, K) pair.
//                   Formed, assembled, discarded immediately. What left- and right-looking use.
//
//   UpdateMatrix    J's ENTIRE Schur complement, u by u, symmetric, destined for J's parent.
//                   One per supernode. Pushed on a stack and carried up the tree. What
//                   multifrontal uses.
//
// A slice versus the whole; per-pair versus per-supernode; immediate versus stacked. Calling both
// "the update" would blur precisely the distinction that separates the three traversals.

#include "oblio/Types.h"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

class NumFactorEngine;

template<class Val>
class UpdateBlock {
public:
    UpdateBlock(std::size_t height, std::size_t width)
        : mHeight(height), mWidth(width),
          mRowIdx(height), mVal(height * width, Val(0)) {}

    std::size_t height() const { return mHeight; }
    std::size_t width()  const { return mWidth; }

    // The block's leading dimension, for BLAS. Column-major, so it is the height.
    std::size_t ld() const { return mHeight; }

private:
    std::size_t mHeight = 0;
    std::size_t mWidth  = 0;

    // The global row index of each local row. Length mHeight.
    std::vector<std::int32_t> mRowIdx;

    // The values, column-major. Length mHeight * mWidth, zeroed on construction, because the
    // update accumulates into it.
    std::vector<Val> mVal;

    friend class NumFactorEngine;
};

extern template class UpdateBlock<double>;
extern template class UpdateBlock<std::complex<double>>;

} // namespace Oblio
