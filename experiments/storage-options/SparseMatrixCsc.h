#pragma once

// SparseMatrixCsc.h - a sparse matrix in compressed sparse column form.
//
// One contiguous array of row indices, one of values, and an offset array saying where each
// column begins. Column j occupies rowIdx/val[colPtr[j] .. colPtr[j+1]-1].
//
// Note that colPtr holds *indices*, not pointers. But real pointers are one step away:
// &mRowIdx[mColPtr[j]] is where column j's row indices start, and it is a plain
// const std::int32_t*. That is the observation this experiment turns on.
//
// MultiplyEngine is a friend, and reaches the storage directly (the friend-access decision).

#include <cstddef>
#include <cstdint>
#include <vector>

namespace StorageOptions {

class MultiplyEngine;

class SparseMatrixCsc {
public:
    SparseMatrixCsc(std::size_t size,
                    std::vector<std::size_t>  colPtr,
                    std::vector<std::int32_t> rowIdx,
                    std::vector<double>       val)
        : mSize(size), mColPtr(std::move(colPtr)),
          mRowIdx(std::move(rowIdx)), mVal(std::move(val)) {}

    std::size_t size() const { return mSize; }
    std::size_t nnz()  const { return mRowIdx.size(); }

private:
    std::size_t               mSize;
    std::vector<std::size_t>  mColPtr;   // length mSize + 1
    std::vector<std::int32_t> mRowIdx;   // length nnz, one contiguous run
    std::vector<double>       mVal;      // length nnz, one contiguous run

    friend class MultiplyEngine;
};

} // namespace StorageOptions
