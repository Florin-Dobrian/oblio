#pragma once

// MmdCompacted.h - MmdFlat on AMD_2'S CLIQUE LAYOUT, one pooled workspace with a free cursor and a
// compaction, over a private copy of the quotient graph.
//
// A CELL OF THE LAYOUT MATRIX. B is a driver on its own branch's vendored layout, C is a driver on
// the other branch's, so this is the mmd counterpart of AmdCompacted and pairs with it down a
// column. It exists because the pool's advantage is currently a single reading on a single algorithm: with
// encoding held equal it earns AmdFlat about 9 percent, rising to 15 at large n, and nothing says
// whether that is the layout or something about how amd walks. See src/MmdCompacted.cpp and
// docs/DESIGN_DECISIONS.md (2026-08-16, the layout matrix).
//
// It returns MmdFlat's permutation, which is genmmd's, and must go on doing so. `make digest` in
// benchmarks/ordering catches a drift; `make mmdorder` in experiments/ordering says correct.
//
// It is otherwise MmdFlat, whose header describes the ordering itself: the prepass, the filing
// convention, the clique-by-clique refresh and its q2h path, and pairwise merging with outmatched
// marking. That description is authoritative there and is not repeated here.

#include "oblio/QuotientGraph.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// How many times the clique pool ran out and had to be compacted. Non-zero says `AMD_2`'s elbow
// room is too small for mmd's cliques, which is the one part of this layout that could not be
// copied across: the branches fill the pool at different rates. Written by orderMmdCompacted.
extern std::size_t gMmdCompactions;

// The elimination order, over the original vertices. delta is as in Mmd1: zero keeps a batch to
// true minima, negative takes one pivot per round.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderMmdCompacted(const std::vector<std::size_t>&  colPtr,
                                     const std::vector<std::int32_t>& rowIdx,
                                     std::int32_t delta = 0);

// The same, reporting how many entries the clique arena ended up holding. An OVERLOAD rather than a
// fourth defaulted parameter, for the reason MmdFlat.h gives at its own pair.
std::vector<std::int32_t> orderMmdCompacted(const std::vector<std::size_t>&  colPtr,
                                     const std::vector<std::int32_t>& rowIdx,
                                     std::int32_t delta,
                                     std::size_t& arenaEntries);

} // namespace Oblio
