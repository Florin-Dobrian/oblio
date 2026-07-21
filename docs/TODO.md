# TODO

Work we intend to do that is not a porting question. The porting ledger tracks whether a 0.9 unit
has been carried over and checked, which is a question with an end; this file holds everything that
outlives it, including things 0.9 never did and that we want anyway.

Two rules keep the split clean. If the item is "0.9 has this and we do not yet", it belongs in
PORTING_LEDGER under Owed. If the item is "neither reference has this, or we want it done
differently", it belongs here. Where an item straddles the two, it lives here and the ledger points
at it, because this is the file meant for browsing.

Each entry says what, why it matters, and what makes it urgent. Order within a section is rough
priority, but nothing here is scheduled.

## Robustness

### Validate the input matrix

Nothing checks that `A` is a well-formed input, and three distinct properties are assumed silently.
Today we generate our own test matrices, so the assumption holds by construction; it stops holding
the moment anyone else calls the solver.

**The diagonal must be structurally present**, even where it is numerically zero. Symbolic
factorization builds a column's index set from A's column structure, so a column with no stored
diagonal never enters its own index set, and the resulting supernode index sets are wrong. This is
not a corner case: a zero diagonal is completely ordinary in indefinite problems, which are exactly
what dynamic LDL exists for. A KKT matrix whose zero block is stored sparsely, dropping the zeros,
hits it immediately. Found on 2026-07-19 while validating dynamic LDL, where it produced both wrong
answers and spurious refusals before the cause was clear.

**No duplicate entries.** `assembleFromA` and `assembleDelay` assign rather than accumulate,
which is correct when each position has one writer and silently keeps only the last value if a
column stores the same row twice. 0.9 accumulates in the matrix case and so tolerates duplicates.
The header already documents CSC with sorted, and by implication unique, row indices, so this is an
unchecked precondition rather than a disagreement.

**Symmetry, and which kind.** Complex Cholesky requires a Hermitian matrix and complex LDL^T
requires a complex-symmetric one, and neither is verified. The failure is silent: `zpotrf` reads
only the lower triangle and assumes the upper is its conjugate, so a complex symmetric matrix is
factored as though it were the Hermitian matrix agreeing with its lower triangle, succeeds, and
returns a plausible wrong answer. The agreed fix is two flags computed on `SparseMatrix` in one
construction pass, `isHermitian()` and `isSymmetric()`, with the engine requiring the first for
Cholesky. For `double` the two coincide; for complex a matrix may be one, the other, or neither.
This entry supersedes the Hermitian item under Owed in PORTING_LEDGER, which is the same fix seen
from the porting side.

The three want one validation pass and one place to report from, so they should be done together
rather than one at a time. Trigger: before anyone outside this effort runs the solver.

## Capability

### Complex Hermitian dynamic LDL: DONE, 2026-07-19

Recorded here as done rather than deleted, because it is the one factorization in the library with
no reference behind it and that is worth being able to find later. 0.9's complex LDL is symmetric
only, so nothing was transcribed; the oracles are the residual in `test_pipeline` and reconstruction
of `L D L^H` on dense fronts.

Complex `DynamicLDLT` turned out to need no kernel change at all, 0.9's complex `factorDynamicLDL_`
differs from its real one in six lines, all declaring the pivot magnitudes real rather than scalar,
which this port had done from the start, so only the dispatch guard widened.

`DynamicLDLH` needed four changes, one place each thanks to the pivot-body merge landing first:
the conjugate in `readPivotBlock2x2`, conjugated `L` where the `D L^H` rows are formed in both
eliminations, and `forceReal` on the diagonal.

**And one that was not in the plan, in `swap`.** Two of its loops move values across the diagonal
and needed conjugating, which was anticipated; the one that was not is the entry *between* the two
swapped positions, `block[at(k_, j_)]`. It is its own reflection under the permutation, so no loop
touches it, and for a symmetric factor it is invariant, which is why 0.9 leaves it alone and why
the omission is invisible until the factorization is Hermitian. Left out, the factor reconstructs
the *conjugate* of the matrix in the affected rows, with no crash, no NaN, and no symptom at all
when the swaps happen to be the identity. Roughly a third of random Hermitian matrices failed.

The lesson worth keeping: **a symmetric-to-Hermitian port is not a search for `conj` opportunities
in the arithmetic.** The arithmetic was nearly all fine. The bug was in a data-movement verb, in an
entry whose symmetry made it a no-op for every case the reference ever ran.

### Multiple right-hand sides### Multiple right-hand sides

`Vector<Val>` carries one right-hand side and the solve is scalar, which is the right call for one
column: there is no level-3 BLAS to be had. Many right-hand sides make a supernode a matrix
operation, and then gathering the rows into a dense block, calling TRSM and GEMM, and scattering
back pays for itself. 0.9 has this as its `MultipleVector` path, so the algorithm is available; what
is missing on our side is `DenseMatrix`, which the ledger lists as a unit and which nothing has
needed yet.

## Structure

Two extractions are available in the dynamic factorization code, recorded here in enough detail to
be picked up cold. Both are refactors of working, tested code: 143 assertions cover this path,
including tier 1 and tier 2 matrices that exercise both pivot selections and both traversals, so
either can be done with a real safety net rather than by eye.

Measured 2026-07-19, code lines with comments stripped: `factorDynamicSupernode` 230 against 0.9's
313, `updateDynamicSupernode` 49 against 61. The port is already about a quarter shorter than the
reference, almost all of it from dropping the error-macro scaffolding and the raw-pointer preamble
rather than from restructuring. What follows is the restructuring that has not been done.

### Merge the duplicated pivot bodies: DONE, 2026-07-19

The eliminations are now `applyPivot1x1` and `applyPivot2x2`, shared by both selection loops, and
the 2x2 block is read once in `readPivotBlock2x2`. 230 code lines became 204 across three functions,
under 147 passing assertions throughout.

The non-goal held: **the selection loops stayed separate.** They are two algorithms rather than one
with flags, differing in three ways, no forced 1x1 in pass 2, a partner scan bounded by the front,
and a Bunch-Kaufman determinant test in place of `max1 == max2`. Merging those behind parameters
would save lines and cost the reader the ability to see what each pass does. That remains the
recommendation if anyone revisits it.

Two things surfaced that the duplication had hidden, which is the usual argument for removing
duplication and is worth recording as evidence for it. **Pass 1 never reads the 2x2 block**: it
decides on magnitudes alone, and the compiler pointed this out with an unused-variable warning the
moment one shared body served both. And `readPivotBlock2x2` is the single place the symmetry of D is
decided, `d12 = d21` being the symmetric statement, which makes complex `LDL^H` a change in one
function rather than four.

### Extract the front expansion shared by both drivers

Not previously recorded. `factorDynamicLeftLooking` and `factorDynamicRightLooking` expand a parent's
index set with about 16 identical lines: extend the index array by the children's total, shift the
existing indices right by that amount, then prepend each child's delayed globals in sibling order.
Verified identical apart from local variable names.

They differ only in the verb that closes the block, `resetVal` for left-looking and `expandVal`
for right-looking, which is the one place the traversals genuinely part company on delayed columns.

So the extraction is clean and the signature is narrow, something like

```
std::int32_t expandFrontIndexSet(NumFactorDynamic<Val>& nf, std::int32_t parent,
                                 const std::vector<std::int32_t>& firstChild,
                                 const std::vector<std::int32_t>& nextSibling) const;
```

returning the number of delayed columns folded in, with each driver calling its own val verb
afterwards (`resetVal` or `expandVal`). That also has a documentation benefit beyond the line count:
it puts the drivers' real difference on one visible line each, instead of at the bottom of a block
that otherwise reads identically in both. (Name uses the current `expand` verb, matching
`expandNodeIdx` / `expandVal`; an earlier draft of this entry said `grow`, before that rename.)

Lower risk and smaller than the pivot merge, and independent of it. Either can be done first.

### The max1 == 0 swap branch: leave duplicated, do not extract

`factorDynamicSupernode` has the same five-line `max1 == 0` branch in both pivot passes (mark pivot
found, `swap` if `j_ != k1_`, set `mPivotType` to 1, advance `j_`, break). It reads as a candidate
for extraction alongside the front expansion, but it is not one, and this note records the decision
so it is not re-proposed. The five lines are welded to two different scan structures: pass 1 does a
two-loop partner scan, pass 2 a one-loop scan, and the branch's `break` and `++j_` act on *that
pass's* loop and cursor. Lifting it would mean either threading "what to do next" signals back out
to the caller, or extracting the whole pivot-accept region, which would bury the
pass-1-versus-pass-2 distinction that is the actual content there. Five visible lines repeated is
the cheaper reading than a helper that has to route control flow. Unlike the front expansion, the
duplication here is the lesser evil.

### Left-looking A-assembly ordering: whether to align its loop skeleton with right-looking

Open design question, deferred for thought, not a bug. The two dynamic drivers currently order the
per-supernode work differently, and the difference is the storage-level shadow of pull-versus-push.
Three orderings are on the table.

Right-looking, as it stands. A separate prepass assembles A into every front first; then the main
loop, per supernode jj, expands jj from its children (`expandVal`, which *preserves* what is already
there), contracts children and assembles their delays, factors jj, and pushes jj's updates up into
ancestors. A is in place before the loop because a push lands on an ancestor not yet reached, and
`assembleFromA` assigns, so it cannot run after a push without overwriting it.

Left-looking, as it stands. No prepass. Per supernode kk, the loop expands kk from its children
(`resetVal`, which *discards*: nothing has been written into kk yet, so the widened block is simply
re-zeroed), assembles A into the now-expanded front at the delayed-column offset, contracts children
and assembles their delays, pulls updates from descendants, and factors kk last. Assembling A after
the expand is what lets the expand discard, and discarding is cheaper than the shift `expandVal`
does: A is placed directly at its final offset rather than written low and moved.

Left-looking, the proposed aligned form. Move `assembleFromA` to the first line of the `kk` loop, so
every driver reads `assemble A, expand, delays, updates-or-factor` in one skeleton. This is not a
pure reorder: assembling A before the expand means the expand can no longer discard, so `resetVal`
must become `expandVal` (preserve and shift), exactly right-looking's verb. The cost is real if
small: left-looking would assemble A into the unexpanded front and then shift those values
during the expand, rather than placing them once at the final offset. It also erases a distinction
the code and `docs/ARCHITECTURE.md` currently draw deliberately: left-looking fills a front from
scratch each turn (`resetVal`), right-looking carries a front's accumulated state through expansions
(`expandVal`). Note this does not remove the prepass asymmetry: left-looking would still have no
prepass, since a front's shape still is not known until its descendants are factored. It removes
only the `resetVal`-versus-`expandVal` asymmetry.

The trade is uniform loop skeleton against a genuine (if small) efficiency and a meaningful
distinction. If taken, the change must switch `resetVal` to `expandVal` in left-looking, re-verify
that A lands identically under assemble-then-shift versus expand-then-assemble-at-offset, and
rewrite the `resetVal`/`expandVal` header comments in `NumFactorDynamic.h` (currently at the two
verbs) and the lazy-assembly paragraph in `ARCHITECTURE.md`, which would no longer be accurate.

### Align, and perhaps unify, the two update kernels

`updateStaticSupernode` and `updateDynamicSupernode` compute the same thing, `block -= L21 D L21^H`,
over the same geometry. Both walk jj's front (`f` columns) against the update rows the ancestor
receives (`L21`, starting at `offset`), form an upper factor `U := D L21^H` into an `f` by `width`
scratch, then close with `gemmLower` for the symmetric square and a plain `gemm` for the rectangle
below. The closing two multiplies are already identical. The geometry matches once the names are
read correctly: static's `f` is the whole front and its stride is `frontSize + updateSize`;
dynamic's `f` is the post-factor front (reduced by the delayed columns, which no longer drive
updates) and its stride is `frontSize + delaySize + updateSize`, the one extra term being the
delayed columns' rows, which a delayed column keeps even though it left as a column.

Near term, align the coding style so the two read as twins. The prologue is where they diverge
cosmetically: static pulls `frontSize`/`numNodeIdx` as named `size_t` locals through accessors and
casts once to `int`; dynamic goes straight to `int` off raw members (`nf.mFrontSize[jj]`,
`nf.mDelaySize[jj]`, `nf.mUpdateSize[jj]`) and is inconsistent with itself, using accessors for
`val`/`nodeIdx` but raw members for the sizes. Bring dynamic to static's form: named `size_t` locals
via accessors, then the cast. Keep the differences that are essential, not cosmetic: dynamic's
`f == 0` early return (every column delayed, no pivot to update with), and its `withHermitian`/`idx`
locals, which the pivot walk needs. Leave the kernel free of the pre/post/expand vocabulary; inside
it jj is simply whatever the caller passed, so bare `frontSize` matches static and does not leak the
expand-block naming.

The one substantive difference is forming `U`. Static calls `formUpper`, which walks D as a pure
diagonal, one `d[dp]` per column (`BlasLapack.cpp`), correct because static's D is diagonal
(identity for Cholesky, a plain diagonal for static LDL). Dynamic inlines its own form-upper loop
because its D is block diagonal: `mPivotType` marks each column 1x1 or 2x2, the 1x1 branch mirrors
`formUpper` line for line, and the 2x2 branch couples two columns through four D entries, which
`formUpper` cannot express. So dynamic's loop is `formUpper` generalized to block-diagonal D.

Longer term, decide whether to truly unify. A pivot-type-aware `formUpper` that handled 2x2 blocks
would let dynamic call out instead of inlining, after which the two kernels would differ only in
where `f` and the stride come from, and could plausibly become one templated kernel. The cost is
coupling the BLAS-adjacent layer to pivot types, and losing the inline pivot walk, which is readable
where it sits. Align first, then take this decision with the symmetric versions side by side.

## Performance

### Measure left-looking against right-looking

The two traversals compute the same factor by the same kernels and differ only in bookkeeping, so
the comparison is clean, and neither side is obviously cheaper. **The reasoning lives in "The life of
an update" in docs/ARCHITECTURE.md**, what the per-pair costs are, why the unit is pairs rather than
supernodes, and what could be done about left-looking's queues if it loses. That belongs there
rather than here: it describes how the thing works, and stays true whether or not this item is ever
picked up.

What is owed here is only the measurement. It wants matrices large enough to separate an allocation
from a linear sweep, and it should report per-pair costs rather than totals, so the two sides are
comparable. If left-looking measures badly, try the flat-array queues described in ARCHITECTURE
before concluding anything about the traversal itself.

Worth doing alongside any other timing work rather than on its own, and worth doing before either
traversal is chosen as a default for anything.

## Testing

### Examples are built but never run

`make` compiles everything in `examples/`, and nothing ever executes it or checks what it prints. So
an example can go stale silently, and one did: `examples/pipeline.cpp` fixed its storage at
`NumFactorStatic` and therefore reported every dynamic LDL cell as "not implemented" for as long as
it took someone to run it by hand and notice. The code was right; the example was lying about it.

An example is documentation that compiles, which is most of its value, but compiling is not the
claim it makes. The cheap fix is to run each example under `make test` and require exit status zero,
which catches a crash but not a wrong number. The honest fix is for examples with deterministic
output to have that output checked. Worth deciding which, rather than drifting into the first
because it is easier.

### Assert the LDL perturbation branch

The branch runs and nothing checks it, which PORTING_LEDGER now records with measurements. The
recipe is known: `StaticLDLT` on `bandIndefinite(40, 3, 0.3, 7)` under Natural ordering perturbs one
pivot and lands at 1.4e-03.

Two assertions, and the second is the one with content: that `numPerturbations()` is nonzero, and
that the factor reconstructs a matrix differing from `A` by about the perturbation. A residual
assertion alone would just record that the answer is worse, which says nothing about whether the
right thing was done.

Worth pairing with a note in the specification that a poor residual from a static factorization on
an indefinite matrix is expected behavior rather than a defect. That confusion has now cost time
twice in one session.

### A specification the suite can be regenerated from

The suite is currently defined by its own source: what is tested is whatever the test files happen
to assert. That is enough to catch regressions and not enough to notice a gap, and it makes the
combinations we deliberately do not support indistinguishable from the ones we forgot.

What we want is a document that states, per combination of ordering, factorization, traversal and
scalar type, whether it is supported and what is asserted about it, plus the catalog of test
matrices and what each is for. Tests are then written from the spec rather than the spec being
recovered from the tests, and a combination that is not yet implemented is a cell that asserts a
refusal rather than a hole.

Being worked on now; this entry stands until the spec exists and the suite matches it.
