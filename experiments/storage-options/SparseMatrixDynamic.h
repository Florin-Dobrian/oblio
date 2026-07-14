#pragma once

// SparseMatrixDynamic.h - a sparse matrix whose structure may change. Stored as a vector of
// vectors, one per column.
//
// **The name is aspirational, and this header says so rather than pretending.** Nothing in this
// experiment mutates: both matrices are built once and read. What is being measured is what the
// *layout* costs, not what the mutation buys. The layout is the one a mutable matrix would need,
// which is why it carries that name.
//
// What would make the name honest is a `setColumn`: replace one column's pattern and values
// without touching its neighbours. That is cheap here (the column owns its buffer) and expensive
// in the flat sibling (every later column shifts). It is not written yet. One step at a time.
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

class SparseMatrixDynamic {
public:
    SparseMatrixDynamic(std::size_t size,
                   std::vector<std::vector<std::int32_t>> rowIdx,
                   std::vector<std::vector<double>>       val)
        : mSize(size), mRowIdx(std::move(rowIdx)), mVal(std::move(val)) {}

    std::size_t size() const { return mSize; }
    std::size_t nnz()  const {
        std::size_t n = 0;
        for (const auto& c : mRowIdx) n += c.size();
        return n;
    }

    // Replace the values, keeping the structure. Cheap here as it is in the static sibling: this
    // is the mutation both layouts support, because it moves nothing.
    bool setValues(const std::vector<std::vector<double>>& val) {
        if (val.size() != mVal.size())
            return false;
        for (std::size_t j = 0; j < val.size(); ++j)
            if (val[j].size() != mRowIdx[j].size())
                return false;
        mVal = val;
        return true;
    }

    // **Replace one column, structure and all. This is what the dynamic layout is for, and the
    // static sibling deliberately does not offer it.**
    //
    // Cheap here, and cheap for a reason that is worth stating: the column owns its own buffer, so
    // giving it a different pattern touches that buffer and nothing else. Its neighbours do not
    // move; no other column is even read. In a flat layout the same operation shifts every later
    // column, which is O(nnz).
    //
    // So the asymmetry between these two classes is not a wart to be tidied away. It **is** the
    // design: an object offers what its storage makes cheap. A common base class, or a matching
    // API, would force one of two lies, either setColumn on the flat matrix (pretending an O(nnz)
    // shift is a column operation), or setColumn on neither (crippling this one to match its
    // sibling's weakness).
    //
    // **This invalidates every pointer previously taken into column j.** The column's buffer is
    // replaced outright, so anything holding `mVal[j].data()` from before is dangling. Callers
    // must re-extract. Nothing enforces this, and nothing can: it is C++.
    //
    // Note the contrast with setValues, which does *not* invalidate: it overwrites the buffers'
    // contents and leaves the buffers where they are. So the rule is exactly "structural mutation
    // invalidates; value mutation does not", which is the same rule the dynamic numeric factor
    // will live by when delayed pivoting grows a front.
    //
    // Returns false if the row indices are out of range or not sorted, or if the two arrays
    // disagree in length. Sorted, because every consumer downstream assumes it and checking here
    // is cheaper than discovering it later.
    bool setColumn(std::size_t j, std::vector<std::int32_t> rowIdx, std::vector<double> val) {
        if (j >= mSize || rowIdx.size() != val.size())
            return false;
        for (std::size_t k = 0; k < rowIdx.size(); ++k) {
            if (rowIdx[k] < 0 || static_cast<std::size_t>(rowIdx[k]) >= mSize)
                return false;
            if (k > 0 && rowIdx[k] <= rowIdx[k - 1])
                return false;   // must be sorted and distinct
        }
        mRowIdx[j] = std::move(rowIdx);
        mVal[j]    = std::move(val);
        return true;
    }
private:
    std::size_t                            mSize;
    std::vector<std::vector<std::int32_t>> mRowIdx;   // one vector per column
    std::vector<std::vector<double>>       mVal;      // one vector per column

    friend class MultiplyEngine;
};

} // namespace StorageOptions
