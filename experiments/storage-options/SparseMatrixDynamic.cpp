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
    // No flat buffer to read nnz from, so sum the column sizes once here and keep the total in mNnz
    // (setColumn maintains it thereafter, so nnz() stays O(1)). Same cap as the static sibling and
    // the main-code matrix: dimension and nnz must fit the std::int32_t index range, since nnz
    // narrows to int at the AMD/MMD ordering boundary.
    std::size_t sum = 0;
    for (const auto& rowIdx : mRowIdx) sum += rowIdx.size();
    mNnz = sum;
    if (mSize > MAX_IDX || mNnz > MAX_IDX)
        throw std::length_error(
            "SparseMatrixDynamic: dimension or nnz exceeds the std::int32_t index range");
}

} // namespace StorageOptions
