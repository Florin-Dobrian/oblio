#pragma once

// Mmd1.h - multiple minimum degree, Oblio's own implementation of the idea, over the shared
// quotient graph. Section 5.11 of archive/sparse_factorization.md, and the mmd1 layer of
// experiments/ordering, which is the specification this is pulled from.
//
// The idea is Liu's (1985) and it is the M in MMD. Refreshing degrees is the expensive step, so
// do it less often: eliminate a whole independent set of least-degree vertices before refreshing
// anything. Non-adjacent pivots cannot disturb each other's degrees, so every pivot in a batch
// is still a true minimum-degree vertex when it is taken.
//
// The independent set is never searched for. It falls out of the bookkeeping: eliminating a
// pivot evicts every vertex it reached from the degree buckets, so whatever is still filed was
// not reached, hence is not adjacent to anything already taken this round, and draining a bucket
// drains an independent set.
//
// What the batch gives up is not what one would guess. The pivots are exact. What is lost is the
// evicted set, which is invisible for the rest of the round, so the choice is made among the
// untouched remainder: a different vertex of the same degree rather than a worse one. Minimum
// degree is famously sensitive to ties, so the fill moves by a fraction of a percent in either
// direction.
//
// This is the idea alone, not the whole of genmmd. The vendored routine additionally carries the
// prepass that numbers degree-0 and degree-1 vertices, mmdupd's q2h path with its merging of
// vertices indistinguishable from each other rather than from the pivot, outmatched marking, and
// its own filing convention. Those are the mmd2 layer's business and are not here, so this
// ordering is not the vendored one and is not trying to be.

#include "oblio/QuotientGraph.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices.
//
// delta widens the batch to vertices within delta of the minimum degree, which buys still fewer
// refreshes at a real concession, those vertices not being minimal. Zero keeps the batch to true
// minima and is what the vendored driver, Sparspak and Liu's paper all use. Negative takes one
// pivot per round, which is plain minimum degree reached through this path and is how the
// prototype checks its batching against the layer below it.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderMmd1(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    std::int32_t delta = 0);

} // namespace Oblio
