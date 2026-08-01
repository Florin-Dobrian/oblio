#pragma once

// Amd1B.h - AMD1's ordering, computed with the eliminator's walk and the driver's first scan
// fused into one pass. Same heuristic, same tie-breaks, same permutation; only the schedule of
// the work differs, which is what the trailing B means.
//
// **B is not a second ordering.** Where Amd1 and Amd2 differ, they differ in what they compute:
// Amd2 has mechanisms Amd1 lacks and produces different permutations and different fill, and both
// are legitimate. Amd1B produces identical output to Amd1 on every input, and if it ever does
// not, one of them is wrong. That makes the pair its own acceptance test, and a stronger one than
// anything else in the benchmark: Amd1 against Amd2 can only be compared on fill, where Amd1
// against Amd1B must agree on the permutation itself.
//
// What it changes. An approximate degree decomposes reach(u) rather than forming it, so its
// driver reads exactly the lists the eliminator has just finished rewriting. Amd1 therefore
// visits each element of A[u] twice and each element of I[u] three times: once to prune, once for
// the clique-degree scan, once for the bound. Amd1B hands the eliminator an `ApproximateScan` and
// gets the first two of those in one visit, leaving the bound as the only second pass, which is
// what `amd_2`'s two scans cost. On a 100x100 grid the visit count falls from 483677 to about
// 216662, the vendored figure.
//
// What it does not change. The bound, the caps, the buckets, the tie-breaks, mass elimination and
// the expansion are Amd1's, untouched. See Amd1.h for the algorithm and section 5.13 of
// archive/sparse_factorization.md for the derivation.
//
// Whether it survives beside Amd1 is a question with a stated answer: when it is
// permutation-identical and faster, it replaces Amd1 and drops the suffix. See docs/TODO.md.

#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderAmd1B(const std::vector<std::size_t>&  colPtr,
                                     const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
