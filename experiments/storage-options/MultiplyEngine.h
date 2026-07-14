#pragma once

// MultiplyEngine.h - one sparse matvec, over either storage.
//
// The question: can a single multiply function serve both a **static** matrix (fixed structure,
// stored flat, in CSC) and a **dynamic** one (mutable structure, stored as a vector of vectors)?
//
// The classes are named for the purpose, the comments for the layout. The layout is a consequence
// of the purpose, and the solver names them the same way (NumFactorStatic, NumFactorDynamic).
//
// The answer here is yes, and not by templating the algorithm. There is exactly **one**
// compiled multiply. It takes no matrix at all. It takes three arrays:
//
//     rowIdxPtr[j]   where column j's row indices start
//     valPtr[j]      where column j's values start
//     len[j]         how many entries column j has
//
// Both storages can fill those. The static one points into its single contiguous buffer
// (&mRowIdx[mColPtr[j]]); the dynamic one reads each inner vector's data(). Once filled, the arrays are
// indistinguishable, and multiply() cannot tell which class produced them, because there is
// nothing left to tell apart.
//
// That is the whole demonstration. The storage question and the algorithm question are
// separable: the layout decides where the pointers come from, and nothing else.
//
// Two extractors (columnPointers, one overload per storage) and one algorithm. Plus a
// hand-written flat multiply as the honest baseline: it skips the pointer arrays entirely and
// walks colPtr directly, which is what any sane CSC code does. If the pointer-array version
// matches it, the indirection costs nothing.
//
// Val is fixed to double. This experiment is about storage, not scalar type.

#include "SparseMatrixStatic.h"
#include "SparseMatrixDynamic.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace StorageOptions {

class MultiplyEngine {
public:
    // The one algorithm. y = A*x, column by column, scattering into y.
    //
    // It sees pointers and lengths. No matrix, no template parameter, no virtual call, no
    // knowledge of any storage class. One function, compiled once.
    void multiply(std::size_t size,
                  const std::int32_t* const* rowIdxPtr,
                  const double* const*       valPtr,
                  const std::size_t*         len,
                  const double* x, double* y) const;

    // Extract the column pointers. One overload per storage; both produce the same three
    // arrays, and after this call the storage is out of the picture.
    //
    // **The pointers are only valid until the matrix's structure changes.** This is the rule that
    // matters, and nothing enforces it:
    //
    //   setValues   preserves them. The buffers stay put; only their contents change.
    //   setColumn   destroys them. The column's buffer is replaced, so any pointer into it dangles.
    //
    // So the invalidation rule tracks *structural* mutation exactly, and in both storages. Refactor
    // with new numbers and the pointers hold; change a pattern and they do not.
    //
    // This is not a quirk of the experiment. It is the rule the solver's dynamic factor will live
    // by: delayed pivoting grows a front, which reallocates its buffer, which invalidates every
    // pointer previously taken into it. A use-after-free there is silent, and it is exactly the
    // failure this experiment exists to make visible before it costs a day.
    void columnPointers(const SparseMatrixStatic& A,
                        std::vector<const std::int32_t*>& rowIdxPtr,
                        std::vector<const double*>&       valPtr,
                        std::vector<std::size_t>&         len) const;

    void columnPointers(const SparseMatrixDynamic& A,
                        std::vector<const std::int32_t*>& rowIdxPtr,
                        std::vector<const double*>&       valPtr,
                        std::vector<std::size_t>&         len) const;

    // The baseline. Hand-written flat: no pointer arrays, walks colPtr directly. This is what
    // one would write if CSC were the only storage, and it is what the general version must
    // match to earn its keep.
    void multiplyStatic(const SparseMatrixStatic& A, const double* x, double* y) const;
};

} // namespace StorageOptions
