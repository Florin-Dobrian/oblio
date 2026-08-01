#pragma once

// Amd1.h - approximate minimum degree, Oblio's own implementation of the idea, over the shared
// quotient graph. Section 5.13 of archive/sparse_factorization.md, and the amd1 layer of
// experiments/ordering, which is the specification this is pulled from.
//
// The idea is Amestoy, Davis and Duff's. Refreshing degrees is the expensive step, and where MMD
// makes it rare (Mmd1.h), AMD makes each one cheap: it refuses to compute the degree exactly at
// all. The exact external degree of u is the size of a union, one per vertex per step, and a
// union is not decomposable. So decompose reach(u) against the clique just formed and replace the
// union by a sum:
//
//     reach(u) = ( A[u] - C[p] ) | ( C[p] - {u} ) | ( C[c] - C[p]  for c in I[u] - {p} )
//     bound(u) = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]|
//
// with two further caps, that nothing exceeds what is still live and that a degree can only have
// grown by the new clique. Every term is exact except the sum, which counts an overlap between
// two old cliques once per clique holding it, so bound(u) >= degree(u) and the pick can be wrong.
//
// What makes it cheap is one property: |C[c] - C[p]| depends on the clique and not on the vertex
// asking, so it is computed once per clique and read once per vertex, where an exact degree walks
// every clique's members per vertex. The gap grows with clique size, which is to say with fill,
// which is to say exactly where the ordering is expensive.
//
// So the two branches give up different things, and the difference matters more than its size.
// MMD1's pivots are exact and only the tie among equals moves. AMD1's may simply not be minimal,
// because an overcounted bound can hide the true minimum. MMD perturbs; AMD can be wrong. Both
// cost well under a percent of fill in either direction.
//
// This is the bound and nothing else. The vendored routine additionally carries aggressive
// absorption, hash detection of supervariables indistinguishable from each other, its scan-1
// subtraction from a maintained clique degree, dense row and column removal, the input
// conditioning of amd_aat and amd_preprocess, the postorder of the assembly tree, and the
// Control/Info interface. Those are the amd2 layer's business and are not here, so this ordering
// is not the vendored one and is not trying to be.

#include "oblio/QuotientGraph.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderAmd1(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
