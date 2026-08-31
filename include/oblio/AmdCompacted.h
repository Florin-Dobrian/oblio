#pragma once

// AmdCompacted.h - AmdFlat on THE AMD ORACLE'S CLIQUE STORAGE: one pool with a free cursor and a
// compaction. PERMANENT, and for two reasons. It is the ALIGNMENT VEHICLE for a differential
// against that oracle, holding cliques the way it does so that whatever still differs is either
// layout or an improvement to carry back into our own ladder; and it is the PREDICTABLE-SPACE
// version of AmdFlat, staying inside `O(n + m)` so the answer is reachable whenever the input
// fits, which our arena cannot promise. See src/AmdCompacted.cpp.
//
// It returns AmdFlat's permutation, which is the amd oracle's raw order, so it is an oracle for
// itself.
//
// It is otherwise AmdFlat exactly. Every design note for that layer is in AmdFlat.h and is
// authoritative there. The one difference is where cliques live: AmdFlat keeps them in a separate
// append-only arena in elimination order, and this keeps them in the same pool as the adjacency
// and incidence lists, with a free cursor and a compaction, which is what the oracle does with its
// single work array.
//
// It is the amd counterpart of MmdChained, which prices our arena against the mmd oracle's
// dead-segment scheme. The amd oracle's is a third design and had never been compared against.

#include "oblio/ElmOrder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
ElmOrder orderAmdCompacted(const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
