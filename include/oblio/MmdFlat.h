#pragma once

// MmdFlat.h - multiple minimum degree with the mechanisms that ride with the batch, over the shared
// quotient graph. Section 5.11 of archive/sparse_factorization.md, and the mmd2 layer of
// experiments/ordering.
//
// MMD1 is the batch idea alone. This adds the four things the mmd oracle does around it, and on a
// grid they are worth about three times MMD1's run time, the largest measured gap in the ordering
// work.
//
//   THE PREPASS. Everything in bucket 1, which after the filing convention below holds the
//   isolated and the degree-1 vertices together, is numbered outright: no clique is formed,
//   nothing is pruned, and the neighbors keep degrees that still count these vertices. The
//   staleness is deliberate and is what the oracle does.
//
//   THE FILING CONVENTION, READ AND NOT TAKEN. The oracle files a degree-0 vertex under degree 1
//   and a refreshed vertex under its degree PLUS ONE, two scales in one bucket array, so a vertex
//   the refresh has touched is penalised by one against a vertex no pivot has reached. WE FILE AT
//   THE TRUE DEGREE at every site and the bucket a vertex sits in IS its degree. See
//   private/MmdCorrected.cpp, which is the reference this driver matches.
//
//   THE CLIQUE-BY-CLIQUE REFRESH. The refresh walks the cliques the round created rather than the
//   vertices they reached, and computes the refreshed clique's weight once per clique, which every
//   member reaches in full. A member with exactly one other source is answered from that weight
//   plus that source, without forming a union at all. That is the TWO-SOURCE PATH, and on grids it
//   takes 36 to 44 percent of all refreshes.
//
//   PAIRWISE MERGING AND OUTMATCHED MARKING. Walking a two-source vertex's one other source
//   reveals vertices whose reach equals its own, which are merged into it, and vertices whose
//   reach contains its own, which cannot be the minimum before it and are WITHHELD from the
//   buckets rather than refiled. A withheld vertex is restored when a later elimination reaches
//   it.
//
// Degrees are weighted here, unlike MMD1's: pairwise merging makes a LIVE vertex stand for
// several originals, so a degree counts what its neighbors stand for rather than how many entries
// they occupy.

#include "oblio/ElmOrder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices. delta is as in Mmd1: zero keeps a batch to
// true minima, negative takes one pivot per round.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
ElmOrder orderMmdFlat(const std::vector<std::size_t>&  colPtr,
                      const std::vector<std::int32_t>& rowIdx,
                      std::int32_t delta = 0);

} // namespace Oblio
