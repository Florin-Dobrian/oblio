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
std::int32_t growFrontIndexSet(NumFactorDynamic<Val>& nf, std::int32_t parent,
                               const std::vector<std::int32_t>& firstChild,
                               const std::vector<std::int32_t>& nextSibling) const;
```

returning the number of delayed columns folded in, with each driver calling its own expansion verb
afterwards. That also has a documentation benefit beyond the line count: it puts the drivers'
real difference on one visible line each, instead of at the bottom of a block that otherwise reads
identically in both.

Lower risk and smaller than the pivot merge, and independent of it. Either can be done first.

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
