#pragma once

// Types.h - shared index-type conventions (see the "index types" design decision),
// plus the enums that select what the numeric factorization computes and how.
//
// IDs (vertices/rows/columns/supernodes, which may carry a "none" marker) are
// std::int32_t; offsets/counts/sizes are std::size_t. NIL is the shared "none"
// sentinel for ID arrays (parent/child/sibling links, supernode maps, etc.).
// MAX_IDX is the largest representable ID; construction guards check inputs against
// it so an over-range value fails loudly rather than wrapping at a narrowing cast.

#include <cstddef>
#include <cstdint>
#include <limits>

namespace Oblio {

// "No such ID" sentinel for std::int32_t link/ID arrays (e.g. a tree root's parent).
inline constexpr std::int32_t NIL = -1;

// Largest valid ID, as a std::size_t for comparison against sizes/offsets. Equal to
// INT32_MAX. Bounds both the matrix dimension (indices are std::int32_t) and, for now,
// nnz (the entry count is narrowed to int at the vendored AMD/MMD boundary); those two
// caps happen to coincide at this value today but have different origins.
inline constexpr std::size_t MAX_IDX =
    static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());

// What the numeric factorization computes. **The symmetry is part of this, not a separate
// setting.**
//
// For Cholesky there is nothing to choose: `A = LL^H` always, and over the reals `L^H` is `L^T`,
// so real symmetric and real Hermitian are the same case. But for LDL, a complex matrix may be
// symmetric (`A = A^T`, `D` complex) or Hermitian (`A = A^H`, `D` real), and those are genuinely
// different factorizations. So the transpose is named:
//
//   ...LDLT     A = LDL^T    plain transpose.        Complex: D is complex.
//   ...LDLH     A = LDL^H    conjugate transpose.    Complex: D is real.
//
// Over the reals the two coincide, and either name computes the same thing.
//
// Folding symmetry into this enum, rather than adding a `Symmetry` setting alongside it, is
// deliberate and it is the same principle as BlasLapack's operation-named wrappers: **make the
// wrong thing unwriteable.** A separate flag would let a caller ask for a complex *symmetric*
// Cholesky, which we would then have to reject at runtime. Here there is no such value to name.
// The combination is not forbidden; it does not exist.
//
// (Why it could not exist: positive definiteness means `x* A x > 0`, which requires that quantity
// to be real, which holds for all `x` exactly when `A` is Hermitian. For a complex symmetric `A`,
// `x^T A x` is a complex number and the inequality does not typecheck. Cholesky has no pivoting to
// recover with, because not needing pivoting is the entire point of Cholesky, and that guarantee
// comes from positive definiteness. LAPACK has no complex-symmetric Cholesky for this reason.)
enum class Factorization {
    Cholesky,      // A = LL^H  (LL^T for real). Requires positive definiteness, hence Hermitian.
    StaticLDLT,    // A = LDL^T, pivots fixed by the symbolic structure. Complex: D complex.
    StaticLDLH,    // A = LDL^H, ditto.                                  Complex: D real.
    DynamicLDLT,   // A = LDL^T, with 2x2 pivoting and delayed columns.
    DynamicLDLH    // A = LDL^H, ditto.
};

// How the supernodes are traversed. All three compute the same factor; they differ in when the
// updates are formed and where they are kept.
//
//   LeftLooking   for each supernode: pull every update owed to it, then factor it.
//   RightLooking  for each supernode: factor it, then push its updates to every ancestor.
//   Multifrontal  for each supernode: assemble its children's update matrices from a stack,
//                 factor, and push its own Schur complement back on the stack.
enum class Traversal {
    LeftLooking,
    RightLooking,
    Multifrontal
};

} // namespace Oblio
