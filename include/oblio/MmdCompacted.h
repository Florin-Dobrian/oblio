#pragma once

// MmdCompacted.h - MmdFlat on THE AMD ORACLE'S CLIQUE LAYOUT, one pooled workspace with a free
// cursor and a compaction, over a private copy of the quotient graph.
//
// A CELL OF THE LAYOUT MATRIX. B is a driver on its own branch's vendored layout, C is a driver on
// the other branch's, so this is the mmd counterpart of AmdCompacted and pairs with it down a
// column. It exists because the pool's advantage is currently a single reading on a single algorithm: with
// encoding held equal it earns AmdFlat about 9 percent, rising to 15 at large n, and nothing says
// whether that is the layout or something about how amd walks. See src/MmdCompacted.cpp and
// notes/DESIGN_DECISIONS.md (2026-08-16, the layout matrix).
//
// It returns MmdFlat's permutation, which is the mmd oracle's, and must go on doing so.
// `make digest` in benchmarks/ordering catches a drift; `make mmdorder` in experiments/ordering
// says correct.
//
// It is otherwise MmdFlat, whose header describes the ordering itself: the prepass, the filing
// convention, the clique-by-clique refresh and its two-source path, and pairwise merging with
// outmatched marking. That description is authoritative there and is not repeated here.

#include "oblio/ElmOrder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices. delta is as in Mmd1: zero keeps a batch to
// true minima, negative takes one pivot per round.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
ElmOrder orderMmdCompacted(const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx,
                           std::int32_t delta = 0);

} // namespace Oblio
