#pragma once

// Mmd3C.h - Mmd3 on the PRODUCTION clique layout, over a private copy of the quotient graph.
//
// TRANSITIONAL, and unlike Mmd3B and Amd3B it is not a permanent alternative. Those two are the
// predictable-space versions of their orderings, staying inside `O(n + m)` so that an answer is
// reachable whenever the input fits. This one changes no storage at all. It exists to carry the
// amd branch's array folds onto the mmd side without editing a class six drivers run, and it is
// expected to be REPLACED by an Mmd3C on AMD_2's clique layout once that is done, which is the cell
// the layout matrix reserves under this name. See src/Mmd3C.cpp and
// docs/DESIGN_DECISIONS.md (2026-08-16).
//
// It returns Mmd3's permutation, which is genmmd's, and must go on doing so for as long as it
// exists. `make digest` in benchmarks/ordering is what catches a drift; `make mmdorder` in
// experiments/ordering is what says correct.
//
// It is otherwise Mmd3, whose header carries the description of the ordering itself: the prepass,
// the filing convention, the element-by-element refresh and its q2h path, and pairwise merging with
// outmatched marking. That description is authoritative there and is not repeated here.

#include "oblio/QuotientGraph.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices. delta is as in Mmd1: zero keeps a batch to
// true minima, negative takes one pivot per round.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderMmd3C(const std::vector<std::size_t>&  colPtr,
                                     const std::vector<std::int32_t>& rowIdx,
                                     std::int32_t delta = 0);

// The same, reporting how many entries the clique arena ended up holding. An OVERLOAD rather than a
// fourth defaulted parameter, for the reason Mmd3.h gives at its own pair.
std::vector<std::int32_t> orderMmd3C(const std::vector<std::size_t>&  colPtr,
                                     const std::vector<std::int32_t>& rowIdx,
                                     std::int32_t delta,
                                     std::size_t& arenaEntries);

} // namespace Oblio
