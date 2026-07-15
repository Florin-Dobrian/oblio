// SparseMatrixStatic.cpp - just the constructor.
//
// It lives here, not in the header, on purpose. The constructor can throw (it guards the dimension
// and nnz against the std::int32_t index range), and an in-header throw was measured to perturb the
// codegen of the templated multiply compiled in the same translation unit: the hot loop slowed
// noticeably even though its own source was unchanged, and removing the throw restored it. Defining
// the constructor here confines the exception path to this file, out of the units that compile the
// multiply. This is also what the main-code SparseMatrix does (its constructor is in SparseMatrix.cpp).

#include "SparseMatrixStatic.h"

#include <stdexcept>   // std::length_error
#include <utility>     // std::move

namespace StorageOptions {

SparseMatrixStatic::SparseMatrixStatic(std::size_t size,
                                       std::vector<std::size_t>  colPtr,
                                       std::vector<std::int32_t> rowIdx,
                                       std::vector<double>       val)
    : mSize(size), mColPtr(std::move(colPtr)),
      mRowIdx(std::move(rowIdx)), mVal(std::move(val)) {
    // Row indices are std::int32_t, and nnz is cast to int at the AMD/MMD ordering boundary, so both
    // the dimension and the entry count must fit that range. Guard at construction so an over-range
    // input fails loudly here rather than wrapping silently at the first narrowing cast. This is the
    // matrix's own nnz cap; setColumn on the dynamic sibling enforces the same ceiling without
    // throwing (it returns false).
    if (mSize > MAX_IDX || mRowIdx.size() > MAX_IDX)
        throw std::length_error(
            "SparseMatrixStatic: dimension or nnz exceeds the std::int32_t index range");
}

} // namespace StorageOptions
