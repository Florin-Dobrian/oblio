#pragma once

// MultiplyEngine.h - one sparse matvec source, over either storage.
//
// The question: can a single multiply serve both a **static** matrix (fixed structure, stored
// flat, in CSC) and a **dynamic** one (mutable structure, stored as a vector of vectors), without
// the engine restating either layout?
//
// The classes are named for the purpose, the comments for the layout. The layout is a consequence
// of the purpose, and the solver names them the same way (NumFactorStatic, NumFactorDynamic).
//
// The answer is yes, through the storage's own per-column lookups. Each matrix carries three:
//
//     rowIdx(j)    where column j's row indices start
//     val(j)       where column j's values start
//     colSize(j)   how many entries column j has
//
// The static one answers from its single contiguous buffer (mRowIdx.data() + mColPtr[j]); the
// dynamic one from each inner vector's data(). A lookup is a fact about the layout, so it lives on
// the storage that holds it, exactly as blockPtr lives on the numeric factor. That is the one place
// the two layouts differ.
//
// multiply is then a template over the matrix. It calls those lookups and never names a member, a
// buffer, or a layout, so its source is written once and serves both storages; the compiler
// specializes it per storage. This is **direct access, the consumer templated on the storage**, and
// it is the same shape the numeric engine uses on the factor. There is deliberately no extractor:
// nothing materializes a pointer-array view, so nothing is owned, nothing goes stale, and no
// consumer carries a columnPointers of its own to restate. That per-engine extractor was the
// repetition this design exists to avoid: a view written once per consumer is a view written many
// times.
//
// Direct is also what the numeric factorization must use regardless (a dynamic factor grows under
// it, so any pointer extracted up front would dangle), so the matvec matches it rather than
// inventing a second pattern. A read-only consumer that later wants the bulk form (one compiled
// kernel over pointer arrays, for the deferred multi-RHS solve) builds those arrays from these same
// lookups, in one shared helper, never a method per engine.
//
// Plus a hand-written flat multiply as the honest baseline: it walks colPtr directly, which is what
// one would write if CSC were the only storage. If the templated version matches it, reaching a
// column through the lookup costs nothing.
//
// Val is fixed to double. This experiment is about storage, not scalar type.

#include "SparseMatrixStatic.h"
#include "SparseMatrixDynamic.h"

#include <cstddef>
#include <cstdint>

namespace StorageOptions {

class MultiplyEngine {
public:
    // The one multiply source. y = A*x, column by column, scattering into y.
    //
    // Templated over the matrix, it reads each column through the storage's own lookups
    // (rowIdx / val / colSize) and knows nothing else about the storage. Direct access: it
    // asks the matrix where a column lives at the moment it needs it and holds no pointer across
    // the loop, so there is nothing to invalidate. This is the same shape the numeric engine uses
    // on the factor, and the reason one source serves both static and dynamic.
    //
    // A pure multiply, y = A*x (BLAS's beta = 0): it overwrites y and does not read the old value,
    // so the caller passes y and does not pre-zero. Because the sweep is column-outer the kernel
    // zeros y once itself, then accumulates, which a column-outer multiply must do anyway.
    //
    // Declared here, defined and explicitly instantiated in the .cpp (for SparseMatrixStatic and
    // SparseMatrixDynamic), the same declaration-in-header / definition-in-cpp discipline the main
    // tree uses for its value templates. The extern template lines below are the in-place reminder
    // of that rule: they suppress an implicit instantiation here, so a body accidentally moved into
    // this header would stop compiling rather than silently instantiate per translation unit.
    template <class Matrix>
    void multiply(const Matrix& A, const double* x, double* y) const;

    // The baseline. Hand-written flat: walks colPtr directly, no lookup call. This is what one
    // would write if CSC were the only storage, and it is what the templated version must match to
    // earn its keep. It reaches the raw buffers through friendship, which is exactly the honest
    // reference the general path is measured against.
    void multiplyStatic(const SparseMatrixStatic& A, const double* x, double* y) const;
};

// The two instantiations that exist, defined in the .cpp. These extern declarations are the
// header-side reminder of the definition-in-cpp rule; the matching `template ...` definitions live
// in MultiplyEngine.cpp.
extern template void MultiplyEngine::multiply<SparseMatrixStatic>(
    const SparseMatrixStatic&, const double*, double*) const;
extern template void MultiplyEngine::multiply<SparseMatrixDynamic>(
    const SparseMatrixDynamic&, const double*, double*) const;

} // namespace StorageOptions
