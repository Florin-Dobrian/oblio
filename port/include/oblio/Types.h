#pragma once

// Types.h — shared index-type conventions (see the "index types" design decision).
//
// IDs (vertices/rows/columns/supernodes, which may carry a "none" marker) are
// std::int32_t; offsets/counts/sizes are std::size_t. NIL is the shared "none"
// sentinel for ID arrays (parent/child/sibling links, supernode maps, etc.).

#include <cstdint>

namespace Oblio {

// "No such ID" sentinel for std::int32_t link/ID arrays (e.g. a tree root's parent).
inline constexpr std::int32_t NIL = -1;

} // namespace Oblio
