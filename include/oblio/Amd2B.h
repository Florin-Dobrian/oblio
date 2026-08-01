#pragma once

// Amd2B.h - AMD2's ordering, computed with the eliminator's walk and the driver's first scan fused
// into one pass. Same heuristic, same tie-breaks, same permutation; only the schedule of the work
// differs, which is what the trailing B means. Amd1B is the same transformation applied to Amd1.
//
// **B is not a second ordering.** A trailing digit means a different ordering: Amd2 has aggressive
// absorption and hash supervariable detection where Amd1 has neither, so their permutations and
// their fill legitimately differ and both are correct. A trailing B means the same ordering on a
// different schedule, so Amd2B must return exactly Amd2's permutation and a difference is a defect
// in one of them.
//
// **And that oracle is worth more here than on the Amd1 pair.** Amd2 carries two mechanisms Amd1
// does not, an absorption pass that kills cliques and a hash pass that folds one live vertex into
// another, so it has more places to go quietly wrong. An identity check against it guards all of
// them at once, and it is the only check in the ordering suite that can: every other pair of
// methods here is two different orderings and can be compared on fill alone.
//
// What it changes, and does not. `QuotientGraph::eliminate` takes an `ApproximateScan` and folds
// the clique-degree scan into the prune, so A[u] is visited once and I[u] twice rather than two and
// three times. Aggressive absorption, the hash pass, the bound, the caps, the buckets and mass
// elimination are Amd2's, untouched: absorption runs after the elimination and touches only the
// incidence lists, so the adjacency weights accumulated during the prune stay valid, and the hash
// merges happen after the bound.
//
// **Measured, it is not faster.** Amd1B, the same transformation on the simpler driver, came out
// about five percent slower at 140x140 on alpamayo, and the diagnosis it was built to act on turned
// out to be wrong: AMD1's gap to the vendored routine is 46 percent data stalls and 27 percent
// work, so merging work was never going to reach it. See benchmarks/ordering/README.md. This file
// exists for the oracle and for the pair being complete, not for speed.

#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderAmd2B(const std::vector<std::size_t>&  colPtr,
                                     const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
