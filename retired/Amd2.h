#pragma once

// Amd2.h - approximate minimum degree with the mechanisms that ride with the bound, over the
// shared quotient graph. Section 5.13 of notes/SPARSE_FACTORIZATION.md, and the amd2 layer of
// experiments/ordering.
//
// AMD1 is the bound alone. This adds the two mechanisms that make the vendored routine what it
// is, both of which act on the same reached set the bound has just been computed over:
//
//   AGGRESSIVE ABSORPTION. A clique whose members lie wholly inside the new one can never
//   contribute anything again, and scan 1 already computed the number that says so. Killing it
//   shortens every incidence list that named it, so the saving compounds.
//
//   HASH SUPERVARIABLE DETECTION. Mass elimination finds vertices indistinguishable from the
//   PIVOT. Two vertices can become indistinguishable from EACH OTHER without either being
//   indistinguishable from it, and those merges are invisible to AMD1. The hash is a filter over
//   the reached set and the exact set comparison is the decision.
//
// What it does not carry, and why, the experiment's catalogue having sorted these:
//
//   THE POSTORDER. Amd.cpp returns a postorder of the assembly tree it builds while eliminating.
//   Oblio builds its own forest from the permuted matrix and orders it in ElmForestEngine, so
//   this would be work done twice and then overwritten. Skipping it changes the permutation and
//   not the fill, a postorder being a relabeling of an elimination tree.
//
//   THE INPUT CONDITIONING. amd_aat and amd_preprocess symmetrize, deduplicate and drop the
//   diagonal. Oblio's matrices satisfy all three by construction.
//
//   CONTROL AND INFO. Oblio has no such interface, and nothing downstream reads an ordering's
//   estimate of anything: SymFactorEngine computes the exact structure later.
//
//   DENSE ROW AND COLUMN REMOVAL. This one is a real capability rather than an artifact, a hub
//   row wrecking an approximate degree, and it is left out for want of evidence rather than on
//   principle: every matrix measured so far is a grid. See notes/TODO.md.
//
#include "oblio/QuotientGraph.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderAmd2(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
