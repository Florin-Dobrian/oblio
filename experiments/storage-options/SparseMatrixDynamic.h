#pragma once

// SparseMatrixDynamic.h - a sparse matrix whose structure may change. Stored as a vector of
// vectors, one per column.
//
// **Named for the purpose: its structure can change.** What makes the name honest is setColumn,
// which replaces one column's pattern and values without touching its neighbours: cheap here (the
// column owns its buffer), expensive in the flat sibling (every later column shifts). The timed
// multiply does not mutate, both matrices are built once and read there, so the benchmark measures
// what the *layout* costs; mutation is exercised separately in testMutation and testInvalidation.
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
// MultiplyEngine is no longer a friend: it reads this matrix through the public per-column lookups
// (rowIdx / val / colSize), the same ones the static sibling offers, so it needs no access to
// the private storage.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace StorageOptions {

class SparseMatrixDynamic {
public:
    // Defined in SparseMatrixDynamic.cpp, not here, on purpose: it sums the column sizes to seed
    // mNnz and can throw on an over-range dimension/nnz, and an in-header throw was measured to
    // perturb the codegen of the templated multiply compiled in the same translation unit. Keeping
    // the throwing body in the .cpp confines the exception path there, matching the static sibling
    // and the main-code SparseMatrix.
    SparseMatrixDynamic(std::size_t size,
                   std::vector<std::vector<std::int32_t>> rowIdx,
                   std::vector<std::vector<double>>       val);

    std::size_t size() const { return mSize; }
    std::size_t nnz()  const { return mNnz; }

    // Per-column accessors, the same three the static sibling offers, under the same names,
    // answered from this layout instead. Here each pointer is the inner vector's own data() and the
    // size is its size(); in the static sibling they come from one flat buffer via colPtr. Same
    // signature, same meaning, different storage, which is exactly what lets a consumer read either
    // matrix through one storage-blind path.
    const std::int32_t* rowIdx(std::int32_t j)  const { return mRowIdx[j].data(); }
    const double*       val(std::int32_t j)     const { return mVal[j].data(); }
    std::size_t         colSize(std::int32_t j) const { return mRowIdx[j].size(); }

    // Replace one column's values, keeping its structure. Identical in signature to the static
    // sibling, and cheap for the same reason: it overwrites this column's value buffer in place,
    // mVal[j], touching nothing else. Value mutation at column granularity, the numbers change and
    // the pattern does not.
    //
    // **This does not invalidate any pointer.** It overwrites the inner vector's contents and
    // leaves the vector (and its data()) where it is, which is the value-mutation half of the rule
    // setColumn completes: structural mutation invalidates, value mutation does not.
    //
    // Returns false if j is out of range, or if the value count does not match the column's length.
    bool setValues(std::int32_t j, const std::vector<double>& val) {
        if (j < 0 || static_cast<std::size_t>(j) >= mSize)
            return false;
        if (val.size() != mRowIdx[j].size())
            return false;
        std::copy(val.begin(), val.end(), mVal[j].begin());
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
    // Returns false if the row indices are out of range or not sorted, if the two arrays disagree
    // in length, or if the change would push nnz past the index cap. Unlike the constructor, which
    // throws, a mutator reports the ceiling the same way it reports its other errors, by returning
    // false. Sorted, because every consumer downstream assumes it and checking here is cheaper than
    // discovering it later.
    bool setColumn(std::int32_t j, std::vector<std::int32_t> rowIdx, std::vector<double> val) {
        if (j < 0 || static_cast<std::size_t>(j) >= mSize || rowIdx.size() != val.size())
            return false;
        for (std::size_t cp = 0; cp < rowIdx.size(); ++cp) {
            if (rowIdx[cp] < 0 || static_cast<std::size_t>(rowIdx[cp]) >= mSize)
                return false;
            if (cp > 0 && rowIdx[cp] <= rowIdx[cp - 1])
                return false;   // must be sorted and distinct
        }
        // nnz after swapping column j's old size for the new one. mNnz >= the old size, so the
        // subtraction cannot underflow.
        const std::size_t newNnz = mNnz - mRowIdx[j].size() + rowIdx.size();
        if (newNnz > MAX_IDX)
            return false;   // would push nnz past the std::int32_t index range
        mNnz = newNnz;
        mRowIdx[j] = std::move(rowIdx);
        mVal[j]    = std::move(val);
        return true;
    }
private:
    // Largest representable index, INT32_MAX: the ceiling both the dimension and nnz are checked
    // against, since indices are std::int32_t and nnz narrows to int at the ordering boundary.
    static constexpr std::size_t MAX_IDX =
        static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());

    std::size_t                            mSize;
    std::vector<std::vector<std::int32_t>> mRowIdx;    // one vector per column
    std::vector<std::vector<double>>       mVal;       // one vector per column
    std::size_t                            mNnz;       // sum of column sizes; set in ctor, kept by setColumn
};

} // namespace StorageOptions
