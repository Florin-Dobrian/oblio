// SparseMatrixDynamic.cpp - just the constructor.
//
// Same reason as the static sibling: the constructor can throw (it guards the dimension and nnz),
// and an in-header throw was measured to perturb the codegen of the templated multiply compiled in
// the same translation unit. Defining it here confines the exception path to this file. This
// constructor also sums the column sizes once to seed mNnz, so nnz() can be O(1) thereafter.

#include "SparseMatrixDynamic.h"

#include <stdexcept>   // std::length_error
#include <utility>     // std::move

namespace StorageOptions {

SparseMatrixDynamic::SparseMatrixDynamic(std::size_t size,
                                         std::vector<std::vector<std::int32_t>> rowIdx,
                                         std::vector<std::vector<double>>       val)
    : mSize(size), mRowIdx(std::move(rowIdx)), mVal(std::move(val)) {
    // mNnz is the sum of the per-column entry counts, summed once here and maintained by setColumn
    // thereafter, so nnz() stays O(1). A loop is what an initializer cannot do, so it belongs in
    // the body, and mNnz's own `= 0` at its declaration is the seed the loop accumulates onto,
    // which is why that initializer is live rather than a value waiting to be overwritten. The
    // accumulator's type is the member's, so the sum cannot silently overflow an int.
    for (const std::vector<std::int32_t>& column : mRowIdx)
        mNnz += column.size();

    // Same cap as the static sibling and the main-code matrix: dimension and nnz must fit the
    // std::int32_t index range, since nnz narrows to int at the AMD/MMD ordering boundary.
    if (mSize > MAX_IDX || mNnz > MAX_IDX)
        throw std::length_error(
            "SparseMatrixDynamic: dimension or nnz exceeds the std::int32_t index range");
}

} // namespace StorageOptions
