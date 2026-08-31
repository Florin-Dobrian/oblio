#pragma once

// AmdFlat.h - AMD2 with the amd oracle's list order, over the shared quotient graph.
//
// **This adds no mechanism.** It carries exactly what AMD2 carries, the approximate degree bound
// with aggressive absorption and hash supervariable detection, and differs from it only in the
// order four things are laid out and one is timed. It is MmdFlat's counterpart on the AMD branch
// and the digit means the same on both: 3 is the layer aligned to the oracle.
//
// WHY A SEPARATE ORDERING FOR A TIE-BREAK. Minimum degree is a tie-break algorithm. At almost
// every iteration several vertices share the least degree, and which is taken is decided by
// whatever the data structure hands over first, so two codes can agree on every rule and part
// company at the first tie. While that is true, no comparison against the oracle can
// separate a difference of MECHANISM from a difference of ARBITRARY CHOICE, and a fill gap cannot
// be attributed to either. Aligning turns the comparison into an equality test.
//
// WHAT IT RETURNS. The amd oracle's permutation exactly, up to the postorder, on the seven example
// graphs
// and on every square grid tested from 3 a side to 40, member order within each supervariable
// included.
//
// THE POSTORDER IS NOT DONE, AND THAT IS THE DESIGN. The oracle ends by relabeling its output as a
// postorder of the assembly tree it built while eliminating. That cannot change correctness
// downstream, since every traversal needs only a TOPOLOGICAL order, which holds for any
// permutation ElmForestEngine sees.
//
// It NEARLY cannot change the fill either, and the qualifier matters. A
// postorder of the ELIMINATION tree cannot: a node is numbered after all its descendants either
// way, so the tree and its fill are unchanged. The oracle postorders its ASSEMBLY tree, and its own
// header
// says that need not be the precise supernodal elimination tree, because mass elimination under an
// approximate degree merges vertices that were never adjacent. So its relabeling is not guaranteed
// fill-neutral. Measured: it agrees on every square grid to 140 a side and on every cubic grid from
// 7 a side up, and differs by one to three entries at 4^3, 5^3 and 6^3. What follows for the
// figures below is one word: they are the RAW order's fill, which is what this layer computes, and
// they equal the oracle's published fill everywhere the benchmarks report rather than by
// construction. What it buys is a smaller multifrontal stack peak, and
// ElmForest reaches that better from the other side, by Liu's rule on the supernodal tree with
// real front and update sizes, which DirectSolver turns on whenever the traversal is multifrontal.
// So that postorder would be work done twice and then overwritten.
//
// AND IT IS A HEURISTIC ONE, which the oracle says itself and which is the stronger half of the
// argument: mass elimination combined with the approximate degree can merge nodes of lower exact
// degree than the pivot, so a clique need not be a fundamental supernode and its diagonal block
// can carry zeros. Its own header therefore states that the assembly tree "is not guaranteed to be
// the precise supernodal elimination tree" and that its postordering "is not guaranteed to be a
// precise postordering" of it. So it is not that we reach the same result by a better route. Theirs
// is a depth-first tidy of an APPROXIMATE tree; Liu's rule on the exact supernodal tree, which is
// what ElmForest computes, supersedes it outright.
//
// It is defensible where it sits. The oracle is a standalone library returning a permutation to
// callers who may have no symbolic phase at all, and for those a depth-first clustering beats the
// raw elimination order with nothing downstream to do better. But it is not optional, so every
// caller who DOES have a symbolic phase pays for a result they discard. Note also that it is a
// separate function, called after the ordering itself returns: the decomposition is
// there in the source and the profile only shows it fused because everything inlines into
// amd_order.
//
// The five differences from AMD2, three of which are conventions in the shared quotient graph and
// so arrive as flags rather than as code here:
//
//   THE HASH BUCKET is walked backward. The oracle builds it by pushing at the head while its
//   second scan walks the new clique forward, so its chain comes out reversed against C[p] and the
//   pair loop takes the
//   head: the survivor of an indistinguishable pair is the member seen LAST. This one is ours,
//   below.
//
//   THE LIST ORDER, both halves, by the amd list-order flag. The reachable set walks
//   the cliques before the explicit adjacency, and the prune puts the new clique at the front of
//   I[u] with the oracle's rotation rather than appending.
//
//   MASS ELIMINATION runs after aggressive absorption rather than inside the eliminator, by
//   QuotientGraphFlat::setLateMassElimination and a call to massEliminate below. Absorption is what
//   makes the cheap structural test agree with the true one, which the oracle says in its own
//   comment, so asking first declines merges it makes.
//
//   THE OWN-WEIGHT SUBTRACTION happens after supervariable detection, not before it. That one was
//   a DEFECT rather than a convention and was fixed in Amd2 and Amd2B as well, where it had been
//   filing every supervariable one bucket too high per vertex a hash merge absorbed.
//
// THREE THINGS IT ALSO TAKES FROM THE ORACLE THAT ARE NOT ORDERING AT ALL, and so are not ledger
// entries: the tagged mark array, which carries seen, absorbed and value in one place where three
// were used, and the hoisted stamp in supervariable detection, which stamps the outer vertex once
// instead of the inner one per pair. Both are pure re-schedules and were landed as a separate
// driver, AMD3C, until `AMD3C == AmdFlat` had confirmed them; that driver is gone and this is it.
//
// And the third: THE HASH KEY IS ACCUMULATED IN THE BOUND'S WALKS rather than in
// a pass of its own, accumulated inside the oracle's second scan. It removes a
// sweep over C[p] and a second walk of A[u] and I[u], 26.70 of this driver's 149.96 element visits
// per pivot at 140 a side and 66.77 of 352.57 at 26 cubed. Measured on alpamayo at about 4 to 7
// percent at six consecutive square grids from 64 to 400 a side and 5 to 14 percent on cubic grids
// from 12 to 32, with nnz(L) unchanged everywhere.
//
// THE CAVEAT THAT USED TO END THIS PARAGRAPH IS WITHDRAWN. It read that the variant was
// timed as a free function and this driver through OrderEngine, a bias of up to 2.4 percent in the
// variant's favor. Measured directly, with an `AMD3f` column that ran this driver down the free
// function path beside the ordinary one, that seam is ZERO across the whole ladder in both
// families. It was a real reading on a different driver on a different day and it does not
// reproduce; every figure that leaned on it should be re-read without the correction.
//
// It is the FIRST of five attempts at this gap to move anything, and it had been tried and
// reverted once already. The difference is that the earlier version carried the key in an array of
// size n; this one files each vertex into its bucket where the key completes and stores nothing.
// The failure was the array rather than the fusion, which is a distinction the four other null
// results had made it easy to stop looking for.
//
// And the fourth, added later the same day: THE FIRST SCAN IS FOLDED INTO THE PRUNE, through
// QuotientGraphFlat's TaggedScan overload of `eliminate`. I[u] is walked twice per pivot and A[u]
// once, which is the oracle's count exactly, where this driver walked them three times and twice.
// Worth 10 to 16 percent on cubic grids, measured with both codes down the same harness path, and
// 0 to 8 percent in 2D. Over eight runs this layer reads 0.83 to 0.89 ms at 16 cubed where the
// oracle reads 0.74 to 0.86, so the two overlap there, and about 1.2x at 32 a side. A
// ratio per row is a poor summary here: the oracle's column moves more between runs than ours
// does. It carries no arrays of
// its own: the two values crossing from the prune to the bound live in `partial` and `hashNext`,
// which are dead at that point, and a version with two fresh vectors of size n was 12 percent
// SLOWER in 2D.

#include "oblio/ElmOrder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// The elimination order, over the original vertices, with what computing it cost in clique
// storage.
//
// It takes A's pattern, which is all an ordering reads, and builds its own quotient graph from it.
ElmOrder orderAmdFlat(const std::vector<std::size_t>&  colPtr,
                      const std::vector<std::int32_t>& rowIdx);

} // namespace Oblio
