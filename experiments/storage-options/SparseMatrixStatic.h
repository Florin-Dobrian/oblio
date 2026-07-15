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
// The raw CSC arrays are public (colPtr / rowIdx / val, exactly as the main-code matrix exposes
// them), so the hand-written baseline reads the flat buffers straight through them. No friendship:
// reading is public, and the storage-blind path reaches columns through the per-column accessors.
// MultiplyEngine is a friend of neither class.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace StorageOptions {

class SparseMatrixStatic {
public:
    // Defined in SparseMatrixStatic.cpp, not here, on purpose: the constructor guards the
    // dimension and nnz and can throw, and an in-header throw was measured to perturb the codegen
    // of the templated multiply compiled in the same translation unit (it slowed the hot loop even
    // though the loop's own source was unchanged). Keeping the throwing body in the .cpp confines
    // the exception path there. This mirrors the main-code SparseMatrix, whose constructor is
    // likewise in its .cpp.
    SparseMatrixStatic(std::size_t size,
                    std::vector<std::size_t>  colPtr,
                    std::vector<std::int32_t> rowIdx,
                    std::vector<double>       val);

    std::size_t size() const { return mSize; }
    std::size_t nnz()  const { return mNnz; }

    // The raw CSC arrays, exposed as the main-code SparseMatrix exposes them. A consumer that knows
    // it holds CSC (the hand-written baseline, or any main-code-style raw walk) reads the flat
    // buffers straight through these and needs no friendship, because reading is public. colPtr
    // appears here, in the static class's own interface, and never in the storage-blind per-column
    // interface below, which both classes share and which has no colPtr because the dynamic sibling
    // has none.
    const std::vector<std::size_t>&  colPtr() const { return mColPtr; }
    const std::vector<std::int32_t>& rowIdx() const { return mRowIdx; }
    const std::vector<double>&       val()    const { return mVal; }

    // Per-column accessors: where column j's row indices and values start, and how many entries it
    // has. Each returns an address (or a size) into the existing storage, O(1), no allocation,
    // nothing owned. This is a fact about the layout, so it belongs to the storage that holds it,
    // and it is the matrix-side twin of the factor's blockPtr. Named for what they return (a
    // column's row indices, its values), not for the array behind them, so the dynamic sibling can
    // offer the same three under the same names. This is the storage-blind interface; colPtr is not
    // in it (see the CSC-specific accessors just above).
    const std::int32_t* rowIdx(std::int32_t j)  const { return mRowIdx.data() + mColPtr[j]; }
    const double*       val(std::int32_t j)     const { return mVal.data()    + mColPtr[j]; }
    std::size_t         colSize(std::int32_t j) const { return mColPtr[j + 1] - mColPtr[j]; }

    // Replace one column's values, keeping its structure. Cheap: the column's values are the
    // contiguous run mVal[colPtr[j] .. colPtr[j+1]-1], so overwriting them touches that run and
    // nothing else, O(colSize), no shift. This is the value half of mutation, the numbers change and
    // the pattern does not, which is what a solver does most often (a Newton iteration, a time step,
    // refactorize). Its identical twin on SparseMatrixDynamic overwrites that column's inner vector.
    //
    // The signature is the same on both classes on purpose: value mutation is cheap at column
    // granularity on either layout, so both offer it the same way. Setting every value is a loop
    // over columns, exactly as reading every column is a loop over the lookups. (A one-shot
    // whole-matrix load would be the bulk form, storage-specific, and we do not need it.)
    //
    // **This does not invalidate any pointer.** The buffer stays where it is; only its contents
    // change. That is the value-mutation half of the rule its dynamic sibling's setColumn completes:
    // structural mutation invalidates, value mutation does not, and it holds in both storages.
    //
    // Returns false if j is out of range, or if the value count does not match the column's length.
    bool setValues(std::int32_t j, const std::vector<double>& val) {
        if (j < 0 || static_cast<std::size_t>(j) >= mSize)
            return false;
        if (val.size() != colSize(j))
            return false;
        std::copy(val.begin(), val.end(), mVal.begin() + mColPtr[j]);
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
    // Largest representable index, INT32_MAX: the ceiling both the dimension and nnz are checked
    // against, since indices are std::int32_t and nnz narrows to int at the ordering boundary.
    static constexpr std::size_t MAX_IDX =
        static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());

    std::size_t               mSize;
    std::vector<std::size_t>  mColPtr;   // length mSize + 1
    std::vector<std::int32_t> mRowIdx;   // length nnz, one contiguous run
    std::vector<double>       mVal;      // length nnz, one contiguous run
    std::size_t               mNnz;      // == mRowIdx.size(); stored for parity with the dynamic class
};

} // namespace StorageOptions
