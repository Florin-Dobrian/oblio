#pragma once

// Types.h — shared index-type conventions (see the "index types" design decision).
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

} // namespace Oblio
