# Dynamic LDL, slice 2 kickoff

2026-07-18. Written at the end of the slice 1 session, so the next session can start
without re-deriving the model. Slice 1 is committed and green (119/119).

## Where we are

Slice 1 ported the pivot core in isolation: `factorDynamicLDL` (0.9's `updateSize == 0`
pass) plus a single-front driver `factorDynamicLeftLooking`, validated by reconstructing
`L D L^T` from the block and `pivotType` against the pivoted, permuted matrix. Both 1x1
and 2x2 pivots fire, and the test asserts that they do. The three growth verbs live on
`NumFactorDynamic` as private, engine-only methods: `extendIndex`, `swap`, `shrinkEntry`.
The dynamic `compute` dispatches `DynamicLDLT` to the new driver for real, left-looking
input only, guarded by `if constexpr (std::is_same_v<Val, double>)`; everything else
still returns false.

The single-front driver refuses anything with `updateSize != 0`, and refuses a supernode
that delays a column, since in slice 1 there is no ancestor to take it. Slice 2 removes
both restrictions.

## The invariant everything rests on

A supernode's block height is `frontSize + numberOfDelayedColumns + updateSize`.

`updateSize` is never rewritten. `factorDynamicLDL` ends by setting
`numberOfDelayedColumns[jj]` to the leftover pivot count and doing `frontSize[jj] -= n`,
so the height is conserved: the delayed columns reclassify from front to delayed rather
than disappearing. This is why `shrinkEntry` is a column truncation that keeps every row.
Any bookkeeping change in slice 2 has to preserve this identity, and it is the first
thing to check when a residual comes out wrong.

## The 0.9 map

All line numbers are `FactorLeftLookingEngineReal-0_9.cc`, 1993 lines total. The complex
twin is the same file with `Real` replaced, so slice 3 can diff the two rather than
re-read from scratch.

The pivot function `factorDynamicLDL_` spans 355 to 942. Its two passes split at the
`if (jjUpdateSize == 0)` test on line 390: pass 1 runs to 656 (ported in slice 1), pass 2
from 657 to about 900, and the shared tail from 900 to 942 sets
`numberOfDelayedColumns` and decrements `frontSize` (also ported). Pass 2 is close to a
copy of pass 1 with the trailing update extended over the update rows, so the port is
mostly a matter of confirming which loop bounds change from `jjFrontSize` to
`jjIndexSize`, not of new logic.

The update kernel `updateDynamicLDL_` spans 1036 to 1130. It builds `D L^T` over the
front into a scratch array, walking `pivotType` to apply 1x1 scalings and 2x2 blocks,
then forms the update with two GEMM calls. Those map onto our existing `gemm` and
`gemmLower` wrappers, so this is the piece that should feel most familiar.

Three `assemble_` overloads matter. The one at 84 takes the matrix and a
`numberOfDelayedColumns` offset, assembling A's original values into a front that has
already been grown. The one at 153 takes two supernode indices, `jj` and `kk`, and folds
`jj`'s delayed columns into its parent `kk`. The one at 210 assembles a temporary into a
front, and is the analogue of our `assembleUpdate`.

The driver is inside `run`, roughly 1527 to 1700 in absolute terms.

## The driver sequence

Per supernode `kk`, in order:

Sum the delayed column counts of the children of `kk`. If that sum is nonzero, call
`extendIndex(kk, n)`, shift `kk`'s existing indices right by `n`, prepend the children's
delayed global indices, then discard the entry block, increase `frontSize[kk]` by `n`,
reallocate the block, and zero it. Left-looking never calls `extendEntry`, it discards
and reallocates, which is why that verb is still unported.

Assemble A's original values into `kk`, passing the delayed count as the offset.

For each updater `jj` of `kk`: if `kk` is `jj`'s parent, assemble `jj`'s delayed columns
into `kk` and then call `shrinkEntry(jj, numberOfDelayedColumns[jj])`. Build the
temporary, call `updateDynamicLDL_`, and assemble the temporary into `kk`.

Factor `kk` with `factorDynamicLDL`. Advance `pp[kk]` by `frontSize + numberOfDelayedColumns`,
not by `frontSize` alone, and enqueue `kk` against its next ancestor.

The ordering matters twice. The delayed columns must be assembled into the parent before
`shrinkEntry` drops them, and the parent must be grown before A is assembled into it,
since the offset assumes the new width.

## Open questions to settle while porting

Whether the owed-list bookkeeping from our static left-looking traversal survives
unchanged once fronts change size mid-traversal, or whether the position array has to be
recomputed rather than advanced. 0.9 advances it, so the answer is probably yes, but it
deserves an explicit check rather than an assumption.

Whether pass 2 can share code with pass 1 rather than being transcribed separately. The
faithful choice is to transcribe first, get it green, and merge afterward if the
difference really is only the loop bounds.

Whether `pivotType` needs to be consulted by `SolveEngine::diagonal` before slice 2 can
be validated end to end. The 2x2 block solve is already anticipated in a comment there,
and a residual test needs it, so this may pull a small amount of solve work into slice 2.

## Validation plan

Slice 2 has no isolation seam the way slice 1 did. The kernels and the driver only prove
out together, so the test is a multi-supernode residual on matrices that actually force
delaying: symmetric indefinite input, several supernodes, and a check that at least one
supernode reports a nonzero `numberOfDelayedColumns` so the machinery is genuinely
exercised rather than trivially skipped. Assert `||Ax - b|| / ||b||` at machine precision,
and keep the slice 1 reconstruction test passing untouched.

When the residual is wrong, the first checks are the height invariant above, then whether
the delayed assembly happened before the shrink, then the `pp` advance.

## Not in scope for slice 2

Complex, both symmetric and Hermitian, is slice 3. Right-looking and multifrontal, and
the `UpdateStackDynamic` that the multifrontal engines befriend, come after. `extendEntry`
stays unported until right-looking needs it.
