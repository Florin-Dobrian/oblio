#pragma once

// SparseMatrixVv.h - a sparse matrix as a vector of vectors.
//
// One std::vector of row indices per column, and one of values. No offset array: each column
// carries its own length, and its own allocation.
//
// This is the layout Oblio would need for dynamic LDL, where a delayed pivot passes columns
// up to an ancestor and that ancestor's front grows at runtime by an amount symbolic never
// predicted. The growth is local, one front grows while its siblings do not, which is what
// a vector of vectors does cheaply and a flat buffer does not.
//
// The pointers are already here: mRowIdx[j].data() is where column j's row indices start,
// and it is a plain const std::int32_t*. The same type CSC produces, from a different place.
//
// MultiplyEngine is a friend, as in the CSC class.

#include <cstddef>
#include <cstdint>
#include <vector>

namespace StorageOptions {

class MultiplyEngine;

class SparseMatrixVv {
public:
    SparseMatrixVv(std::size_t size,
                   std::vector<std::vector<std::int32_t>> rowIdx,
                   std::vector<std::vector<double>>       val)
        : mSize(size), mRowIdx(std::move(rowIdx)), mVal(std::move(val)) {}

    std::size_t size() const { return mSize; }
    std::size_t nnz()  const {
        std::size_t n = 0;
        for (const auto& c : mRowIdx) n += c.size();
        return n;
    }

private:
    std::size_t                            mSize;
    std::vector<std::vector<std::int32_t>> mRowIdx;   // one vector per column
    std::vector<std::vector<double>>       mVal;      // one vector per column

    friend class MultiplyEngine;
};

} // namespace StorageOptions
