#pragma once

// Amd3B.h - Amd3 on AMD_2's clique storage: one pool with a free cursor and a garbage collection.
// PERMANENT, and for two reasons. It is the ALIGNMENT VEHICLE for a differential against AMD_2,
// holding cliques the way AMD_2 does so that whatever still differs is either layout or an
// improvement to carry back into our own ladder; and it is the PREDICTABLE-SPACE version of AMD3,
// staying inside `O(n + m)` so the answer is reachable whenever the input fits, which our arena
// cannot promise. See src/Amd3B.cpp and docs/DESIGN_DECISIONS.md (2026-08-16).
//
// It returns Amd3's permutation, which is AMD_2's raw order, so it is an oracle for itself.
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

// How often the pool had to be compacted during the last ordering. `AMD_2` reports the same figure
// as Info[AMD_NCMPA]; `Mmd3C` publishes `gMmd3CCompactions` for its half. See src/Amd3B.cpp.
extern std::size_t gAmd3BCompactions;

// The elimination order, over the original vertices.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderAmd3B(const std::vector<std::size_t>&  colPtr,
                                     const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
