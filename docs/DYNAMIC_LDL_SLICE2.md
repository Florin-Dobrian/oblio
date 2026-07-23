# Dynamic LDL, slice 2 kickoff

2026-07-18. Written at the end of the slice 1 session, so the next session can start
without re-deriving the model. Slice 1 is committed and green (119/119).

## Where we are

Slice 1 ported the pivot core in isolation: `factorDynamicSupernode` (0.9's `updateSize == 0`
pass) plus a single-front driver `factorDynamicLeftLooking`, validated by reconstructing
`L D L^T` from the block and `pivotType` against the pivoted, permuted matrix. Both 1x1
and 2x2 pivots fire, and the test asserts that they do. The three expansion and contraction verbs live on
`NumFactorDynamic` as private, engine-only methods: `expandNodeIdx`, `swap`, `contractVal`.
The dynamic `compute` dispatches `DynamicLDLT` to the new driver for real, left-looking
input only, guarded by `if constexpr (std::is_same_v<Val, double>)`; everything else
still returns false.

The single-front driver refuses anything with `updateSize != 0`, and refuses a supernode
that delays a column, since in slice 1 there is no ancestor to take it. Slice 2 removes
both restrictions.

## The invariant everything rests on

A supernode's block height is `frontSize + delaySize + updateSize`.

`updateSize` is never rewritten. `factorDynamicSupernode` ends by setting
`delaySize[jj]` to the leftover pivot count and doing `frontSize[jj] -= n`,
so the height is conserved: the delayed columns reclassify from front to delayed rather
than disappearing. This is why `contractVal` is a column truncation that keeps every row.
Any bookkeeping change in slice 2 has to preserve this identity, and it is the first
thing to check when a residual comes out wrong.

## The 0.9 map

All line numbers are `FactorLeftLookingEngineReal-0_9.cc`, 1993 lines total. The complex
twin is the same file with `Real` replaced, so slice 3 can diff the two rather than
re-read from scratch.

The pivot function `factorDynamicLDL_` spans 355 to 942. Its two passes split at the
`if (jjUpdateSize == 0)` test on line 390: pass 1 runs to 656 (ported in slice 1), pass 2
from 657 to about 900, and the shared tail from 900 to 942 sets
`delaySize` and decrements `frontSize` (also ported). Pass 2 is close to a
copy of pass 1 with the trailing update extended over the update rows, so the port is
mostly a matter of confirming which loop bounds change from `jjFrontSize` to
`jjIndexSize`, not of new logic.

The update kernel `updateDynamicLDL_` spans 1036 to 1130. It builds `D L^T` over the
front into a scratch array, walking `pivotType` to apply 1x1 scalings and 2x2 blocks,
then forms the update with two GEMM calls. Those map onto our existing `gemm` and
`gemmLower` wrappers, so this is the piece that should feel most familiar.

Three `assemble_` overloads matter. The one at 84 takes the matrix and a
`delaySize` offset, assembling A's original values into a front that has
already been expanded. The one at 153 takes two supernode indices, `jj` and `kk`, and assembles
`jj`'s delayed columns into its parent `kk`. The one at 210 assembles a temporary into a
front, and is the analogue of our `assembleUpdate`.

The driver is inside `run`, roughly 1527 to 1700 in absolute terms.

## The driver sequence

Per supernode `kk`, in order:

Sum the delayed column counts of the children of `kk`. If that sum is nonzero, call
`expandNodeIdx(kk, n)`, shift `kk`'s existing indices right by `n`, prepend the children's
delayed global indices, then discard the entry block, increase `frontSize[kk]` by `n`,
reallocate the block, and zero it. Left-looking never calls `expandVal`, it discards
and reallocates, which is why that verb is still unported.

Assemble A's original values into `kk`, passing the delayed count as the offset.

For each updater `jj` of `kk`: if `kk` is `jj`'s parent, assemble `jj`'s delayed columns
into `kk` and then call `contractVal(jj, delaySize[jj])`. Build the
temporary, call `updateDynamicLDL_`, and assemble the temporary into `kk`.

Factor `kk` with `factorDynamicSupernode`. Advance `pp[kk]` by `frontSize + delaySize`,
not by `frontSize` alone, and enqueue `kk` against its next ancestor.

The ordering matters twice. The delayed columns must be assembled into the parent before
`contractVal` drops them, and the parent must be expanded before A is assembled into it,
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
supernode reports a nonzero `delaySize` so the machinery is genuinely
exercised rather than trivially skipped. Assert `||Ax - b|| / ||b||` at machine precision,
and keep the slice 1 reconstruction test passing untouched.

When the residual is wrong, the first checks are the height invariant above, then whether
the delayed assembly happened before the contract, then the `pp` advance.

## Names

The port renamed its traversals onto the pivoting axis after this note was first written:
`factorStaticLeftLooking` and `factorStaticRightLooking` for Cholesky and static LDL,
`factorDynamicLeftLooking` for the dynamic driver, and `factorDynamicSupernode` for the pivot
kernel (the dynamic counterpart of `factorStaticSupernode`). By the same rule the update kernel
ported in this slice should land as `updateDynamicSupernode`, beside the existing
`updateStaticSupernode`. 0.9's
own names keep their trailing underscore here (`factorDynamicLDL_`, `updateDynamicLDL_`) so that
line references stay greppable against the reference sources.

## What the reference covers beyond this slice

0.9 implements dynamic LDL in all three traversals and for both scalar types, so nothing after this
slice needs inventing. Two things about that are worth knowing before opening those files.

**Right-looking is a near twin of left-looking.** `FactorRightLookingEngineReal-0_9.cc` has
`factorDynamicLDL_` at 357 and `updateDynamicLDL_` at 1038, against 355 and 1036 in the
left-looking file, and the two files are the same length to within a few lines. So once this slice
is green, right-looking should be closer to a transcription with the direction flipped than to new
work.

**Multifrontal contains three versions of `factorDynamicLDL_`, and two of them are dead.** They
share one signature, which compiles only because the preprocessor hides all but one:

```
352   #ifdef UNDEF          355   version 1     599   #endif      (dead)
600   //#ifdef UNDEF        603   version 2    1189   //#endif    (LIVE: the guard is commented out)
1190  #ifdef UNDEF         1193   version 3    1508   #endif      (dead)
```

**The live one is the middle version, at 603.** Its guard is commented out, which is what leaves it
in the build; the other two are compiled out and were evidently abandoned. Porting version 1 or 3
would mean porting code the author rejected. `updateDynamicLDL_` at 1591 has no such twin.

That multifrontal needed three attempts is itself a signal. It is the traversal where delayed
columns meet the update stack (`UpdateStackDynamic`, which only the multifrontal engines befriend),
and it is the right one to leave for last.

## Not in scope for slice 2 (as written at the time)

Complex, both symmetric and Hermitian, is slice 3. Right-looking and multifrontal, and
the `UpdateStackDynamic` that the multifrontal engines befriend, come after. `expandVal`
stays unported until right-looking needs it.

*All of that except multifrontal was finished the same day. Right-looking came first, on the rule
that the two traversals travel together, and `expandVal` was exactly what it needed. Complex
followed: `DynamicLDLT` needed only a dispatch guard, since the port was already in 0.9's complex
form, and `DynamicLDLH` was an extension with no reference. The live record is PORTING_LEDGER and
TODO; this file is history.*

## Closed, 2026-07-19

Slice 2 is done and verified end to end, and right-looking came with it rather than after it, on
the rule that the two traversals travel together. So this note is history now; the current record is
PORTING_LEDGER for what is ported, TESTING_SPECIFICATION for what is checked, and TODO for what is
owed.

Two of its predictions were wrong and are worth keeping visible, since the same guesses are
available to make again. Pass 2 is **not** pass 1 with different loop bounds: the arithmetic bodies
are character-identical and the whole difference is in the selection (no forced 1x1, a 2x2 partner
restricted to front columns, and a Bunch-Kaufman determinant test in place of `max1 == max2`). And
the solve needed all three passes rewritten, not just the diagonal, because the leading dimension
changes as well as the pivot handling.

One prediction was right and paid off immediately: `expandVal` was exactly what right-looking
needed and nothing else did.
