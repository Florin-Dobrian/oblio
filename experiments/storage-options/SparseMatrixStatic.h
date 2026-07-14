#pragma once

// SparseMatrixStatic.h - a sparse matrix whose structure is fixed. Stored flat, in compressed
// sparse column form.
//
// **Named for the purpose, not the layout.** The flat layout is a *consequence* of the structure
// being fixed, not the thing being chosen: nothing grows, so nothing needs room to grow in, so one
// contiguous buffer with offsets is the obvious storage. Its sibling here (SparseMatrixDynamic) is
// the same matrix with the opposite premise. Naming them Csc and Vv would name the bytes; static
// and dynamic name the reason, and that is the vocabulary the solver uses for NumFactorStatic and
// NumFactorDynamic too.
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

class SparseMatrixStatic {
public:
    SparseMatrixStatic(std::size_t size,
                    std::vector<std::size_t>  colPtr,
                    std::vector<std::int32_t> rowIdx,
                    std::vector<double>       val)
        : mSize(size), mColPtr(std::move(colPtr)),
          mRowIdx(std::move(rowIdx)), mVal(std::move(val)) {}

    std::size_t size() const { return mSize; }
    std::size_t nnz()  const { return mRowIdx.size(); }

    // Replace the values, keeping the structure. Cheap: nothing moves.
    //
    // This is the mutation the solver actually does most often, a Newton iteration, a time step,
    // same pattern, new numbers, refactorize, and the flat layout is perfectly happy with it.
    //
    // **This does not invalidate any pointer.** The buffer stays where it is; only its contents
    // change. That is worth saying, because its dynamic sibling's setColumn *does* invalidate, and
    // the difference is not cosmetic: the rule is "structural mutation invalidates, value mutation
    // does not", and it holds in both storages.
    //
    // Returns false if the count does not match, since a value array of the wrong length can only
    // mean the caller thinks this is a different matrix.
    bool setValues(std::vector<double> val) {
        if (val.size() != mVal.size())
            return false;
        mVal = std::move(val);
        return true;
    }

    // **There is deliberately no setColumn here, and its absence is the design.**
    //
    // Changing one column's *structure* in a flat layout means shifting every later column: it is
    // O(nnz), not O(column). An API that looks cheap and is secretly linear in the whole matrix is
    // a trap, and the caller who writes it in a loop will not find out until their program crawls.
    //
    // Refusing to offer it is not a limitation. It is telling the truth about the storage. The
    // caller knows whether they are changing one column or rebuilding, and can pick the object
    // that suits: SparseMatrixDynamic has setColumn, because there a column owns its own buffer
    // and the operation really is cheap. To restructure a static matrix, build a new one.
    //
    // The rule, in one line: **an object offers what its storage makes cheap, and nothing else.**
private:
    std::size_t               mSize;
    std::vector<std::size_t>  mColPtr;   // length mSize + 1
    std::vector<std::int32_t> mRowIdx;   // length nnz, one contiguous run
    std::vector<double>       mVal;      // length nnz, one contiguous run

    friend class MultiplyEngine;
};

} // namespace StorageOptions
