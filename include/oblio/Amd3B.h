#pragma once

// Amd3B.h - Amd3 on AMD_2's clique storage. TEMPORARY; see src/Amd3B.cpp for what is being
// measured and for the stop condition. It returns Amd3's permutation, which is AMD_2's raw order,
// so it is an oracle for itself.
//
// It is otherwise Amd3 exactly. Every design note for that layer is in Amd3.h and is authoritative
// there. The one difference is where cliques live: Amd3 keeps them in a separate append-only arena
// in elimination order, and this keeps them in the same pool as the adjacency and incidence lists,
// with a free cursor and a garbage collection, which is what `AMD_2` does with `Iw`.
//
// It is the amd counterpart of Mmd3B, which prices our arena against genmmd's dead-segment scheme.
// AMD_2's is a third design and had never been compared against.

#include "oblio/Types.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderAmd3B(const std::vector<std::size_t>&  colPtr,
                                     const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
