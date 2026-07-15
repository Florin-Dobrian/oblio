// SparseMatrixDynamic.cpp - just the constructor.
//
// Same reason as the static sibling: the constructor can throw (it guards the dimension and nnz),
// and an in-header throw was measured to perturb the codegen of the templated multiply compiled in
// the same translation unit. Defining it here confines the exception path to this file. This
// constructor also sums the column sizes once to seed mNnz, so nnz() can be O(1) thereafter.

#include "SparseMatrixDynamic.h"

#include <numeric>     // std::accumulate
#include <stdexcept>   // std::length_error
#include <utility>     // std::move

namespace StorageOptions {

SparseMatrixDynamic::SparseMatrixDynamic(std::size_t size,
                                         std::vector<std::vector<std::int32_t>> rowIdx,
                                         std::vector<std::vector<double>>       val)
    : mSize(size), mRowIdx(std::move(rowIdx)), mVal(std::move(val)),
      mNnz(std::accumulate(mRowIdx.begin(), mRowIdx.end(), std::size_t{0},
                           [](std::size_t sum, const auto& rowIdx) { return sum + rowIdx.size(); })) {
    // mNnz is the sum of the per-column entry counts, seeded once here and maintained by setColumn
    // thereafter, so nnz() stays O(1). The std::size_t{0} seed is load-bearing: it fixes the
    // accumulator type, so the sum cannot silently overflow an int.
    //
    // Same cap as the static sibling and the main-code matrix: dimension and nnz must fit the
    // std::int32_t index range, since nnz narrows to int at the AMD/MMD ordering boundary.
    if (mSize > MAX_IDX || mNnz > MAX_IDX)
        throw std::length_error(
            "SparseMatrixDynamic: dimension or nnz exceeds the std::int32_t index range");
}

} // namespace StorageOptions
