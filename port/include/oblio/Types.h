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

// What the numeric factorization computes.
//
// The scalar type and the factorization together determine the symmetry, which is therefore
// never a separate parameter (see the factorization-space entry in DESIGN_DECISIONS):
//
//   real    + anything      symmetric,  conjugate-transpose is plain transpose
//   complex + Cholesky      Hermitian,  A = LL^H
//   complex + LDL           symmetric,  A = LDL^T
//
// Cholesky on a complex *symmetric* matrix is not merely unimplemented, it is meaningless:
// positive definiteness requires x* A x to be real, which holds for all x exactly when A is
// Hermitian. That combination is rejected, permanently.
enum class Factorization {
    Cholesky,     // A = LL^H (LL^T for real). Requires positive definiteness.
    StaticLDL,    // A = LDL^T, pivots fixed by the symbolic structure.
    DynamicLDL    // A = LDL^T, with 2x2 pivoting and delayed columns.
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
