#pragma once

// Amd3.h - AMD2 with the vendored routine's list order, over the shared quotient graph. Section
// 5.13 of archive/sparse_factorization.md, and the amd3 layer of experiments/ordering.
//
// **This adds no mechanism.** It carries exactly what AMD2 carries, the approximate degree bound
// with aggressive absorption and hash supervariable detection, and differs from it only in the
// order four things are laid out and one is timed. It is MMD3's counterpart on the AMD branch and
// the digit means the same on both: 3 is the layer aligned to the vendored code.
//
// WHY A SEPARATE ORDERING FOR A TIE-BREAK. Minimum degree is a tie-break algorithm. At almost
// every iteration several vertices share the least degree, and which is taken is decided by
// whatever the data structure hands over first, so two codes can agree on every rule and part
// company at the first tie. While that is true, no comparison against the vendored routine can
// separate a difference of MECHANISM from a difference of ARBITRARY CHOICE, and a fill gap cannot
// be attributed to either. Aligning turns the comparison into an equality test.
//
// WHAT IT RETURNS. `AMD_2`'s permutation exactly, up to the postorder, on the seven example graphs
// and on every square grid tested from 3 a side to 40, member order within each supervariable
// included. Its nnz(L) is 206332 at 100 a side and 474995 at 140, which is what
// benchmarks/ordering/README.md records for the vendored AMD.
//
// THE POSTORDER IS NOT DONE, AND THAT IS THE DESIGN. Amd.cpp ends by relabeling its output as a
// postorder of the assembly tree it built while eliminating. That cannot change the fill, since a
// node is numbered after all its descendants either way, and it cannot change correctness
// downstream, since every traversal needs only a TOPOLOGICAL order, which holds for any
// permutation ElmForestEngine sees. What it buys is a smaller multifrontal stack peak, and
// ElmForest reaches that better from the other side, by Liu's rule on the supernodal tree with
// real front and update sizes, which DirectSolver turns on whenever the traversal is multifrontal.
// So AMD's postorder would be work done twice and then overwritten. See experiments/ordering/AMD3.md.
//
// The five differences from AMD2, three of which are conventions in the shared quotient graph and
// so arrive as flags rather than as code here:
//
//   THE HASH BUCKET is walked backward. Amd.cpp builds it by pushing at the head while its scan 2
//   walks Lme forward, so its chain comes out reversed against C[p] and the pair loop takes the
//   head: the survivor of an indistinguishable pair is the member seen LAST. This one is ours,
//   below.
//
//   THE LIST ORDER, both halves, by QuotientGraph::setVendoredListOrder. The reachable set walks
//   the cliques before the explicit adjacency, and the prune puts the new clique at the front of
//   I[u] with Amd.cpp's rotation rather than appending.
//
//   MASS ELIMINATION runs after aggressive absorption rather than inside the eliminator, by
//   QuotientGraph::setLateMassElimination and a call to massEliminate below. Absorption is what
//   makes the cheap structural test agree with the true one, which Amd.cpp says in its own
//   comment, so asking first declines merges it makes.
//
//   THE OWN-WEIGHT SUBTRACTION happens after supervariable detection, not before it. That one was
//   a DEFECT rather than a convention and was fixed in Amd2 and Amd2B as well, where it had been
//   filing every supervariable one bucket too high per vertex a hash merge absorbed.
//
// Which of the four is which, and what each cost, is experiments/ordering/AMD3.md.
//
// TWO THINGS IT ALSO TAKES FROM Amd.cpp THAT ARE NOT ORDERING AT ALL, and so are not ledger
// entries: the tagged W array, which carries seen, absorbed and value in one place where three
// were used, and the hoisted stamp in supervariable detection, which stamps the outer vertex once
// instead of the inner one per pair. Both are pure re-schedules and were landed as a separate
// driver, AMD3C, until `AMD3C == AMD3` had confirmed them; that driver is gone and this is it.
// experiments/ordering/AMD3.md, iteration 15.

#include "oblio/QuotientGraph.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
std::vector<std::int32_t> orderAmd3(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
