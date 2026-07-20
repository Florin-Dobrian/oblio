# Testing Specification

What the test suite covers, stated independently of the test sources. The intent is that tests are
written from this document rather than this document being recovered from the tests, so that a
combination we do not support is a stated expectation rather than an absence, and a gap is
visible as a gap.

This first version is a **catalog**: an inventory of the suite as it stands on 2026-07-19, taken
from the tests as they actually run rather than from reading their sources. Reorganizing it into a
specification proper, with the full combination space stated cell by cell, comes next.

**This document and the suite move together**, which is an invariant and is stated as one in
CLAUDE.md, where it is loaded every session. It is repeated here for a reader who arrives at this
file first, but CLAUDE.md is where it lives.

## Running

`make test` builds and runs every suite. Each suite is a standalone binary that prints one line per
assertion and a count, and returns nonzero on any failure. Suites are discovered by wildcard, so a
new `tests/*.cpp` file needs no Makefile change. Totals today: **148 assertions across 8 suites**.

| suite | assertions | what it establishes |
|---|---|---|
| `smoke` | 5 | the tree builds and the basic objects work |
| `test_permutation` | 11 | the index map and its composition |
| `test_order` | 21 | AMD and MMD produce valid permutations |
| `test_forest` | 23 | elimination forest, supernodes, amalgamation |
| `test_symfactor` | 29 | supernodal index sets against a dense oracle |
| `test_numfactor` | 16 | the numeric factor, by oracle and by reconstruction |
| `test_solve` | 14 | the solve step, by residual |
| `test_pipeline` | 29 | whole-pipeline combinations, by residual |

## A note on the word dynamic

It carries two meanings in these names, and both are correct under the rule in WRITING_RULES: on a
data object static and dynamic describe *storage*, on an algorithm they describe *pivoting*. So in
`test_solve` the pairs named `static` and `dynamic` are the two storages, `NumFactorStatic` and
`NumFactorDynamic`, both running statically pivoted factorizations. Dynamic *pivoting* appears in
exactly one place today, the slice 1 assertion in `test_numfactor`.

## Oracles

Four independent checks are in use, and the distinction matters more than the count. An oracle
shares no code path with what it checks, so agreement is evidence rather than restatement.

**Dense fill simulation** (`denseFactorPattern` in `tests/test_util.h`). Permutes A, then simulates
Cholesky fill densely: eliminating a column makes every pair of its subdiagonal rows adjacent. The
naive cubic formulation, deliberately. Checks the elimination forest and the symbolic factorization.

**Dense Cholesky** (`denseCholesky` in `test_numfactor.cpp`). Textbook `A = L L^H`, no sparsity, no
supernodes. Checks the numeric factor for Cholesky.

**Reconstruction.** Multiply the computed factors back out and compare against A, or against the
pivoted and permuted A where pivoting reordered it. Used for LDL in all its variants, where no
independent dense implementation exists.

**Residual.** `||Ax - b|| / ||b||` through the whole pipeline, using `MultiplyEngine` for the
product. The only check that exercises ordering, symbolic, numeric and solve together.

## Catalog

### smoke, 5 assertions

Structural symmetry of a hand-built matrix, its size and its entry count under full storage, and
that the identity permutation is valid and round-trips.

### test_permutation, 11 assertions

`setOldToNew` and `setNewToOld` each adopt a map and build the inverse, and agree on the same
permutation from either direction. Both setters reject a duplicate, an out-of-range value and a
negative one, leaving the object untouched. Composition is checked for six properties: with its own
inverse it gives the identity, it is order sensitive in both directions, it is not commutative,
identity on either side is neutral, and a size mismatch is refused. A random sweep of 500 checks
composition against direct application and the inverse against the identity.

### test_order, 21 assertions

Seven matrices, each checked for structural symmetry and then ordered by AMD and by MMD, with the
result checked for validity as a permutation. Matrices: a 6x6 arrow, tridiagonals at n = 1, 2, 10
and 100, a 5x5 diagonal, and a complex arrow.

The orderings are checked for *validity*, not against 0.9's output, and not for quality. Nothing
asserts that AMD or MMD reduces fill.

### test_forest, 23 assertions

Parent, child and sibling links, roots, height, column sizes, fundamental compression and threshold
amalgamation. Small cases with hand-computable answers: tridiagonals at n = 4 and 6 (a path, with
the last two columns merging), a 6x6 arrow (one supernode), two disjoint blocks (two trees), a star
at n = 4 (a sibling chain), and a grid. A random sweep of 200 checks links, height and sizes in both
the nodal and fundamental regimes, and that amalgamation never increases the supernode count; at
threshold 8 it merges 768 supernodes into 214. Four assertions cover `exactPatterns`, the predicate
distinguishing an index set with no stored zeros from one carrying them.

Amalgamation is greedy and not canonical, so only its tie-break-invariant properties are asserted.

### test_symfactor, 29 assertions

Index sets against the dense oracle. Small cases with counted fill: a path (no fill), a vee (1), a
5-cycle (2), a 3x3 grid (8), two blocks (none, two trees), and a 6x6 arrow ordered hub first (10
fill) and hub last (0), which is the standard demonstration that ordering is the whole game. A
dense 4x4 gives one supernode; a tail case gives supernodes {0} and {1,2,3,4}.

A random sweep of 200 checks four things: nodal and fundamental compression give the same factor,
compression merges (1494 supernodes into 784), natural and AMD orderings both match the oracle, and
amalgamation loses no true nonzero while introducing explicit zeros as extra indices. Both index
paths are asserted to be taken, 6 cases exact and 194 carrying stored zeros, so neither branch is
silently untested.

### test_numfactor, 16 assertions

Cholesky against dense Cholesky: a 3x3, then left- and right-looking for real and complex, then a
10x10 grid under natural and AMD orderings in both traversals. A non-positive-definite 2x2 is
refused rather than mangled. A random sweep of 40 matrices confirms every one factors. One assertion
records that the grid case is structurally non-trivial (forest height 90, AMD cutting supernodes
from 90 to 81 and indices from 954 to 501), so the comparison is not passing on a degenerate case.

LDL by reconstruction, in all three symmetries: real `L D L^T` in both traversals, complex symmetric
`L D L^T` with complex D, and complex Hermitian `L D L^H` with real D. One assertion records that no
perturbation was needed, the inputs being diagonally dominant.

Two assertions concern storage rather than arithmetic. `static into dynamic` runs every static
factorization into `NumFactorDynamic` and checks the result is identical to the flat one, across all
symmetries and both traversals. `dynamic LDL slice 1` is the only dynamic *pivoting* assertion:
`L D L^T` reconstructs a single dense front, with both 1x1 and 2x2 pivots asserted to fire.

### test_solve, 14 assertions

The multiply itself is checked once against a hand-computed product. A random sweep of 40 matrices
with 10 right-hand sides each confirms every system solves, run separately for the two storages.

The residual assertions cover five combinations, each in both traversals, and each for both
storages: real Cholesky, real static LDL^T, complex Cholesky and complex LDL^H against Hermitian
input, and complex LDL^T against complex-symmetric input. A 10x10 grid is checked separately in both
storages. All are ordered by AMD.

### test_pipeline, 29 assertions

Added 2026-07-19, with slice 2 of dynamic LDL. Where `test_numfactor` checks the factor against an
oracle and `test_solve` checks the solve, this suite checks that the phases *compose*, for a given
combination of ordering, factorization, traversal and storage. When a residual here goes wrong, the
focused suites are where to look next.

**Ordering is Natural throughout, deliberately.** A fill-reducing ordering would make these
assertions depend on AMD's tie-breaking, so the matrices are instead built already in a good order,
banded or grid-structured, and the ordering step is asked to do nothing. Numerical behavior is then
the only variable. Real arithmetic only so far.

The tiers describe how hard a matrix is to pivot, not how large it is.

**Tier 0, no pivoting required.** A 6x6 grid Laplacian, 36 columns, numbered row-major. Symmetric
positive definite and diagonally dominant. Ten assertions covering all five factorizations, each in
both implemented traversals: the three statically pivoted ones (Cholesky, StaticLDLT, StaticLDLH)
each in both storages, worst-cased over the traversals; then each dynamic one asserted twice, once
that the residual is right in both directions and once that neither delayed anything nor chose a
2x2.

The second of each pair is the assertion with teeth. A pivot search that delays unnecessarily still
produces a correct answer and a machine-precision residual, so the residual cannot tell "did the
right thing" from "did something pointless but arithmetically sound". Only the count can.

For real input LDL^T and LDL^H are the same computation. Both are run rather than one being assumed
to stand in for the other, and tier 1 asserts the equivalence directly.

**Tier 0 refusals, two assertions.** Combinations that must return false rather than answer:
dynamic pivoting into flat storage, and the multifrontal traversal. These are as much a part of the
specification as the working cells, because a cell that begins returning a plausible wrong answer
instead of a refusal is exactly the failure a port invites.

**Tier 1, mild pivoting.** Two banded indefinite matrices of half-bandwidth 3 with random
off-diagonals and half the diagonal entries zeroed, at n = 40 seed 7 and n = 24 seed 7. Four
assertions, a residual and a count triple for each:

| matrix | delayed | supernodes delaying | 2x2 pivots | residual |
|---|---|---|---|---|
| band n=40, w=3, zeros 0.50, seed 7 | 5 | 5 | 4 | 7.21e-16 |
| band n=24, w=3, zeros 0.50, seed 7 | 3 | 3 | 3 | 1.00e-15 |

Two further assertions run the n = 40 matrix again and require **bit-identical** results, same
counts and the same residual to the last bit:

`DynamicLDLH` against `DynamicLDLT`, which is where the claim that the two transposes coincide over
the reals is actually tested, on an input that pivots, rather than asserted in a comment.

Right-looking against left-looking. These are two different drivers over the same two kernels
(0.9's `factorDynamicLDL_` and `updateDynamicLDL_` are byte-identical between its two engines), and
they expand a front by opposite means: left-looking discards an empty front and rebuilds it, while
right-looking must carry forward the values A and the already-factored descendants left in it.
Agreement on a matrix that delays is what says the second of those is right, and it is the only
assertion that exercises `expandVal` at all.

**The counts are pinned exactly, and the matrices are built to make that meaningful.** `std::mt19937`
has its output sequence fixed by the standard, but the distribution templates do not: their
algorithms are implementation-defined, so the same seed yields different doubles under libstdc++ and
libc++. Assertions that only bound the residual are indifferent to this, since the bound holds for
the whole family, but an exact count is a property of one particular matrix. The test therefore
derives its doubles from the engine's raw output rather than through `uniform_real_distribution`,
which makes each matrix a pure function of its seed on every platform.

Reproducible matrices are necessary for an exact count but not sufficient, and the gap is worth
naming. A pivot is accepted when `|d| >= threshold * max`, and both sides of that come out of BLAS
updates, so a matrix with a pivot near the threshold could be decided differently under different
rounding and the count would then differ between platforms with nothing wrong. **Verified on
2026-07-19 across two standard libraries and two BLAS implementations** (libstdc++ with reference
BLAS, libc++ with Accelerate): the residuals differ, as they should, and both count triples are
identical, so neither matrix is near-tied. That is a property of these two matrices rather than a
guarantee, and a future tier 1 matrix should be checked the same way before its counts are pinned.

Exact counts make these change detectors on purpose. What they detect is a change in *pivoting
behavior*, which nothing else in the suite can see: a pivot search that delays unnecessarily still
produces a correct answer and a residual at machine precision. Legitimate causes for the numbers to
move exist, a different default pivot threshold or a different amalgamation among them, and when one
lands the numbers are re-recorded here and in the test together.

Worth recording, because it cost an hour to learn: **most zero diagonals do not delay.** They fill
in from the Schur complement before they are reached, which is precisely why quasi-definite systems
factor without pivoting. Grid Laplacians with zeros punched into the diagonal, banded KKT blocks,
and hub-and-spoke structures were all tried and all delayed nothing. What does delay is a column
reached while still small relative to its own column, and the banded family above produces those
reliably.

**Complex, five assertions.** All ten (factorization, scalar type) cells are supported as of
2026-07-19, and this section covers the dynamic ones. `StaticLDLT` on a diagonally dominant complex
symmetric band; `DynamicLDLT` on the same family with half the diagonal zeroed, asserted for residual
in both traversals and for the traversals agreeing on 2 delayed columns and 2 two-by-two pivots; and
and `DynamicLDLH` on a genuinely Hermitian
band, conjugate off-diagonals and a real diagonal, likewise asserted for residual in both traversals
and for the traversals agreeing. The counts are bounded rather than pinned, and this is the case
that proved that necessary: the two platforms disagree on them (see the tier 2 note above). The
traversals must still agree with each other, which is the claim worth asserting, since both run on
the same machine.

`DynamicLDLH` is the one factorization in the library with **no reference behind it**: 0.9's complex
LDL is symmetric only. Its oracles are the residual here and reconstruction of `L D L^H` on dense
fronts.

`DynamicLDLT`'s matrices are complex **symmetric**, `A = A^T` with complex diagonal entries, which
is what LDL^T means over the complex field and is not the same thing as Hermitian. `DynamicLDLH`
gets a Hermitian one. Handing either the other's matrix would test nothing.

**Three factorizations, three different matrices, and they cannot be shared.** This cost two failing
assertions to learn. A static factorization cannot pivot, so a zero diagonal is *perturbed* rather
than delayed and the residual is then honestly poor (4.8e-03 observed), the perturbation branch
working as designed, not a defect. So `StaticLDLT` gets the dominant matrix and `DynamicLDLT` gets
the one with zero diagonals. Cholesky is absent on purpose: it needs Hermitian positive definite
input, a third matrix and a different property, and a complex symmetric matrix is not a valid input
for it at all. Complex Cholesky is covered in `test_numfactor` and `test_solve`.

Worth stating as a rule, since it has now caught the suite three times: **a failing assertion on a
new matrix is more often the matrix than the code.** Singular inputs, structurally absent diagonals,
and now a matrix handed to a factorization whose preconditions it does not meet.

**The facade, two assertions.** The tier 0 matrix again, this time through `DirectSolver`, over all
five factorizations and both traversals: that all ten are reached, and that the worst residual is at
machine precision.

Not redundant with the by-hand sweep. The facade owns both factors and chooses between them with
`dynamicPivoting()`, so it can fail to reach a combination that works. That is exactly the defect
`examples/pipeline.cpp` carried: it fixed the storage at `NumFactorStatic`, so every dynamic cell
printed "not implemented" long after it was implemented. Nothing caught it because examples are
built by `make` and never run. The `reached == 10` assertion is the guard against that shape of
error, and it is why the count is asserted separately from the residual: a silently skipped
combination would otherwise leave the worst residual looking perfect.

**Tier 2, heavy pivoting, four assertions.** Two families, and between them they cover the two ways
pivoting gets hard: many delays, and no 1x1 pivot available at all.

**Saddle point**, `[[H, B^T], [B, 0]]` at 30 variables and 12 constraints, with **both** blocks
carrying a zero diagonal. The honest use case for an indefinite solver, and the family that delays
hardest, since a constraint column has no diagonal to pivot on and no update can give it one. Two
assertions: residual under 1e-11 in both traversals, and heavy delaying with the traversals
agreeing. Observed 92 delayed columns across 25 supernodes, 16 1x1 and 13 2x2 pivots, residual
9.43e-14.

Worth knowing before reaching for this family: a *nonzero* `H` diagonal makes it tier 0 again and it
delays nothing at all.

**Zero-diagonal tridiagonal** at n = 12 and n = 24. The extreme case: with nothing on the diagonal
no 1x1 pivot can ever be accepted, so every pivot is a 2x2 and the factorization is essentially
exact, residuals at 9.05e-17 and 1.69e-16. One assertion each, pinning `1x1 == 0` and
`2x2 == n / 2`.

**Even order only.** At odd order this matrix is exactly singular, condition number around 1e16, so
it has no residual to hit and tests nothing. That trap produced a residual of 0.27 that looked like
a defect in dynamic LDL and was not.

**Counts are bounded at tier 2 and pinned at tier 1, which is deliberate and the opposite way round
from what looks natural.** The more pivot decisions a matrix forces, the likelier one of them sits
near enough the acceptance threshold to be decided differently under different rounding; tier 1
makes a handful and was verified identical across two BLAS implementations, tier 2 makes dozens and
that verification would not be worth leaning on.

In the event tier 2 agreed too. On 2026-07-19 the saddle point reported 92 delayed columns, 16 1x1
and 13 2x2 pivots under both libstdc++ with reference BLAS and libc++ with Accelerate, with only the
residual moving (9.43e-14 against 4.61e-14), and both zero-diagonal cases agreed as well.

**A counterexample turned up the same day, and it settles the question.** The complex Hermitian case
in the complex section below does *not* agree across platforms:

| | delayed | supernodes | 1x1 | 2x2 | residual |
|---|---|---|---|---|---|
| x86-64, libstdc++, reference BLAS | 2 | 2 | 28 | 2 | 1.76e-15 |
| Apple Silicon, libc++, Accelerate | 3 | 3 | 26 | 3 | 5.40e-15 |

Same matrix, the generator is deterministic and complex `DynamicLDLT` reports identical counts on
both, and genuinely different pivot decisions.

**Why it happens.** The pivot search asks a yes/no question, `|d| >= threshold * max1`, of quantities
that come out of BLAS updates. Two machines can produce different last bits there for several
reasons at once: a different BLAS implementation, different vectorization and therefore a different
summation order, and fused multiply-add contracting differently. Near the boundary those last bits
decide the answer. One platform accepts a column where the other delays it to its parent, and from
that point the two factorizations differ in structure.

**Why it does not matter for accuracy, which is the part worth being clear about.** Both are valid
factorizations of the same matrix. The threshold is not selecting one canonical answer; it is
refusing pivots small enough to amplify error. Any choice that passes it yields a backward-stable
factorization, so both platforms land at machine precision, 1.76e-15 and 5.40e-15 are both simply
zero. The arithmetic is not fragile; the decision boundary is, and by design.

**The diagnostic, when this recurs.** Residual at machine precision but counts differing across
machines is expected and is not a defect. Residual *degraded* is a defect, wherever the counts land.
And the traversals disagreeing with each other on one machine is always a defect, since they share
hardware, library and BLAS, which is why that comparison, rather than the absolute counts, is what
the complex assertions check.

Pinned counts here would have gone red on a correct build, and the newest, least-referenced change
in the tree would have been the obvious first suspect.

So the rule is not "counts agree in practice, pin them where convenient". It is: **pin a count only
where the platforms have been checked and the matrix makes few enough decisions to make that
meaningful, and bound it everywhere else.** Tier 1 qualifies, at a handful of decisions and verified
on both. Tier 2 and the complex cases do not. The all-2x2 claim is the exception and is pinned,
because it is structural rather than numerical: with nothing on the diagonal there is no 1x1 to
choose, so the count follows from the matrix and not from a comparison rounding could tip.

**Singular matrices are excluded throughout**: they have no residual to hit, and asserting something
weaker about them would look like coverage without being any.

## What the catalog shows

Three gaps stand out, stated as facts rather than as a plan.

**The combination space is still barely sampled.** Three orderings, five factorizations, three
traversals and two scalar types is ninety cells. `test_solve` exercises ten, all under AMD, and
`test_pipeline` adds ten more under Natural. Multifrontal is asserted to refuse and is otherwise
absent. MMD appears only in `test_order`, where its output is checked for validity and never used
to factor anything.

**Refusals are now asserted, but only two of them.** `test_pipeline` covers the two that remain
reachable. The remaining unsupported cells, both complex dynamic factorizations among them, are
still indistinguishable from tests nobody wrote. Those two are unsupported for different reasons and
should not be recorded as one debt: `DynamicLDLT` is unported, `DynamicLDLH` is underived.

**Dynamic pivoting is covered at tiers 0 and 1 only.** Delayed columns crossing a forest, 2x2
pivots in the solve, and the expansion of a parent's front are now exercised, but mildly: five delayed
columns is the most any assertion sees. Tier 2 is where the machinery would be pushed, and it does
not exist yet.

**Nothing tests the pivot threshold.** It defaults to 0.1 and it directly controls how much
delaying happens, so every dynamic assertion silently tests one value of it.
