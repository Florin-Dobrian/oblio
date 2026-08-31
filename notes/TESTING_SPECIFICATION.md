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
new `tests/*.cpp` file needs no Makefile change.

It then runs each example in `examples/` and checks its exit status, printing one line apiece with
their output discarded. That is a weaker check than a suite and is not counted among the assertions
below: it catches an example that crashes, that returns a failure, or that has quietly stopped
being built, and says nothing about whether the numbers it prints are right. The stronger version,
checking deterministic output, is open in notes/TODO.md and is awkward while residuals are in the
output, since those legitimately differ in the last bits across BLAS
implementations. Totals today: **265 assertions across 8 suites** with the vendored orderings
present, **237 without**. Measured 2026-08-24 by building and running both modes.

**RETIRED 2026-08-21: MMD1, MMD2, AMD1 AND AMD2.** Those four drivers are in `retired/`, out of the
build and out of the `Ordering` enum, which is now five values. Sixteen validity assertions in
`test_order` went with them and twelve more across the other suites, taking the totals from 279 and
265 to the figures above. Both were read off a run. See `retired/README.md`.

**CORRECTED 2026-08-17, AND THE SAME FAILURE AS THE ONE RECORDED BELOW.** This file said 283 and
269, `README.md` and `CLAUDE.md` said the same, `test_order.cpp`'s own header comment said a third
thing, and the suite ran **261 and 247**. Retiring `AMD1B` and `AMD2B` on 2026-08-15 removed
fourteen sameness assertions and, because the pair also left the enum, fourteen validity ones;
`Mmd3B` added six back, for a net 22. None of the four files moved with them. The figures above are
the post-correction totals INCLUDING the coverage added the same day, and they were measured by
building and running every suite rather than by arithmetic.

**Why nothing caught it, which is the part worth keeping.** Every figure was internally consistent,
each file agreed with itself, and the suite passes whatever a document claims about how many
assertions it has. This is the second recorded instance, after the 2026-08-08 one below, and both
have the same shape: a test removed or added and a count left behind. The only reliable check is to
run the suites and count, which takes one command.

**The count depends on the build, and that is deliberate.** The vendored MMD, the corrected MMD and
the vendored AMD live in
`private/`, which is not published, and both builds detect rather than require it. Twenty-eight
assertions in `test_order` check those three routines and compile only when they are there, so that
suite reports 73 or 45. It was fourteen and 59 until `MmdCorrected` arrived on 2026-08-23, which
brought a validity assertion and a `MmdFlat == MmdCorrected` assertion on each of the seven
matrices. Nothing else varies: every other suite asserts the same thing either way,
everything they use being ours. The one place the difference shows outside `test_order` is the
ordering sweep in `test_pipeline`, which expects every ordering the build has rather than a fixed
number. That sweep covers three orderings without `private/` and five with it, having dropped
AMD1B and AMD2B on 2026-08-15 and MMD1, MMD2, AMD1 and AMD2 on 2026-08-21 when they left the enum.

| suite | assertions | what it establishes |
|---|---|---|
| `smoke` | 5 | the tree builds and the basic objects work |
| `test_permutation` | 11 | the index map and its composition |
| `test_order` | 73 / 45 | the orderings are valid, and each driver of ours reproduces its branch's reference |
| `test_forest` | 29 | elimination forest, supernodes, amalgamation, multifrontal child order |
| `test_symfactor` | 29 | supernodal index sets against a dense oracle |
| `test_numfactor` | 18 | the numeric factor, by oracle and by reconstruction |
| `test_solve` | 14 | the solve step, by residual |
| `test_pipeline` | 86 | whole-pipeline combinations, by residual |

The counts in this table had drifted from the suite before MMD1 was added, `test_pipeline` reading
48 against 56 on disk and the total reading 153 against 183. They are corrected here rather than
left, since a stale count is exactly the kind of quiet untruth the invariant above exists to
prevent.

## What the suite does NOT check, and how the amd and mmd alignments are verified instead

`test_order` checks that each ordering returns a valid permutation, and that each of the three
non-enum layers reproduces its original. It does NOT check either alignment against its vendored
oracle, and it cannot: `genmmd` and `AMD_2` live in `private/`, which is unpublished, and the
comparison needs a hooked scratch copy of the vendored source rather than a link against it.

That verification lives in `experiments/ordering`, is run by hand, and is described in that folder's
README under "The two alignment checks". In short:

```
experiments/ordering:  make test          63 checks, twins against each other and
                                          prototype against production
                       make aligned       both alignment checks, 38 cases each
                       make amdorder      production Amd3 against the vendored AMD's raw
                                          elimination order
                       make mmdorder      production Mmd3 against genmmd's elimination order
                       make mmdmatrices   the same mmd check on REAL matrices, over data/*.mtx
```

**`make mmdmatrices` is the wider half of the mmd check, added 2026-08-15**, and it reported 243
matched, 0 differed, 3 skipped on its first run over the SuiteSparse matrices already fetched for
the two benchmark reports. It asserts exactly what `make mmdorder` asserts, the permutation entry
for entry against genmmd, on input nothing in `graphs.h` can generate: dense rows, pure diagonals
with every vertex isolated, social and nearest-neighbor graphs, near-identical siblings that test
tie-breaking, and sizes two orders past the grid check. It is not part of `make aligned`, since its
result depends on what has been downloaded into a gitignored directory, and it exits clean when
nothing has been. `experiments/ordering/README.md` carries the account.

**The two are named for the branch and not for the mechanism**, and the mechanisms differ: genmmd
emits the order DIRECTLY, while AMD gets it through a HOOK. `make amdorder` generates a hooked copy
of `private/Amd.cpp` on demand, because `amd_order` returns a vector `AMD_postorder` has already
relabeled and the raw order has to be reconstructed upstream of it; the copy is gitignored and
removed by `make clean`, so it cannot become a stale oracle. `make mmdorder` needs none of that:
`mmd_order` returns the order genmmd eliminates in and there is no postorder anywhere in the
routine, so the vendored output vector is the object to compare. The amd target was `raworder`
until 2026-08-09, and was renamed because a name describing one target's mechanism made a matched
pair read as two unlike things.

**Four shapes each, not one shape at many sizes**: the seven examples, 2D grids from 4 a side to
140, 3D grids from 2 to 24, and nine random patterns at n = 2000. The distinction is load bearing
rather than thorough for its own sake. Widening a square grid exercises scale and never mechanism,
and the 2D-only version of the amd check was green while production `Amd3` carried a stale clique
degree that a 3D grid at 16 a side finds. It also found a use-after-free in the shared
`QuotientGraphFlat` that every ordering had, which no assertion in this suite could see because the
program was reading its own freed memory and getting the right answer back. Both are in
`notes/DESIGN_DECISIONS.md` (2026-08-09).

**The mmd check is newer than the mmd alignment**, which is worth knowing when reading dates. The
alignment was established on 2026-08-07 by a scratch probe that did not survive its session, and
what ran between then and 2026-08-09 was the benchmark's fill column, which `MmdFlat.md` iteration 6
shows is not sufficient: fill was exact at every size while the permutation still diverged at pivot
700 of 1024.

The same widening is available to `make test`, whose prototype-against-production comparison still
runs on 2D grids alone; `graphs.h` now holds the 3D and random builders, so that is a change of one
list rather than of any code.

**The object compared for amd is the RAW elimination order**, not `amd_order`'s output vector. The
vendored routine relabels its result with a postorder that Oblio deliberately does not reproduce, so
the output vectors differ by construction. The fill agrees anyway, 206332 at 100 a side and 474995
at 140, but NOT by construction, which is a correction made on 2026-08-09. A postorder of the
elimination tree cannot change fill; AMD postorders its assembly tree, which its own header says
need not be that tree, so its relabeling is not guaranteed fill-neutral.

**AND IT IS NOT FILL-NEUTRAL AT LARGE SIZES EITHER, corrected 2026-08-24.** This paragraph said the
two agree "on every square grid and on cubic grids from 7 a side up", differing by one to three
entries at 4^3, 5^3 and 6^3 alone. `make scale3d` on alpamayo disagrees at 81 a side:

```
81^3    AmdVendored  628021440       the postordered output
        AMDraw       628021443       the raw elimination order, and ours
```

Three entries, reproduced in two consecutive runs, with `AmdFlat` and `AmdCompacted` matching
`AMDraw` exactly as they must. So the difference is a property of the vendored postorder and not of
ours, it is not confined to tiny cubes, and it is not monotone in size: 65^3 agrees in the same
runs.

The safe statement is the one the paragraph opens with, that the relabeling is not guaranteed
fill-neutral, and the list of sizes where it happens not to be is an observation about the sizes
measured rather than a rule. The acceptance test compares the PERMUTATION, so none of this reaches
it.

Anyone extending the suite should know that the strongest evidence for the two orderings'
correctness is not in it.

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

### test_order, 73 assertions (28 of them the vendored routines', and optional)

Seven matrices, each checked for structural symmetry and then ordered by five enum methods and
checked for validity as a permutation: AmdVendored, MmdVendored, MmdCorrected, MmdFlat and AmdFlat.
Matrices: a 6x6 arrow, tridiagonals at n = 1, 2, 10 and 100, a 5x5 diagonal, and a complex arrow.

**`MmdFlat == MmdCorrected` is asserted on all seven**, which is the mmd branch's oracle in the
suite rather than only in the benchmarks. `MmdCorrected` is genmmd with its degree scale repaired
and is what our three mmd drivers reproduce; `MmdVendored` is the frozen original and nothing is
compared against it. That pair of assertions per matrix is 14 of the 28 that need `private/`.

**Three drivers are reached as FREE FUNCTIONS here and are the strongest oracle in this suite**,
because each must reproduce its arena twin ENTRY FOR ENTRY. Every other pair of orderings here can
only be checked for validity, each being a different ordering whose permutation legitimately
differs; a driver that is its twin computed differently has no such licence, so an identical
permutation is a requirement and any difference is a defect in one of the two.

**THE SAMENESS ASSERTION COMPARES THREE THINGS, and only the first is the permutation.** Each also
compares PEAK LIVE CLIQUE MEMBERS and MEMBERS BORN between the twin and its arena original, both of
which are properties of the ALGORITHM rather than of the layout: two drivers running the same method
form the same cliques and merge the same vertices at the same moments whatever their storage, so
both figures must agree exactly. Two drivers agreeing on the permutation while doing different work
is the failure this catches and the permutation cannot, and it is a shape of defect this tree has
found by hand more than once. A figure of zero means the driver does not track it and that
comparison alone is skipped, which today is `MmdChained` on both counts: chained storage ends a
clique at a terminator and keeps no length to subtract on death.

The permutation is compared against the ENGINE's, so the assertion also covers the enum reaching the
driver it names; the two work figures are compared between the two FREE FUNCTIONS, because the
engine calls one of them and comparing its figures with that driver's would compare a call with
itself.

| driver | is | permanent |
|---|---|---|
| `MmdChained` | MmdFlat on genmmd's dead-segment clique storage | yes |
| `AmdCompacted` | AmdFlat on AMD_2's compacted storage | yes |
| `MmdCompacted` | MmdFlat on AMD_2's compacted storage | yes |

All three joined the `Ordering` enum on 2026-08-21 and this suite still reaches them as free
functions, deliberately: the free function is the entry point the comparison wants, with no
`Permutation` built around it. **None of the three is in `test_pipeline`'s ordering sweep**, which
names five enumerators by hand, so this suite and `make digest` in `benchmarks/ordering` remain the
only things checking them end to end. That is why the coverage here is uniform rather than sampled:
one validity assertion per driver on the arrow, and one sameness assertion per driver on every one
of the seven matrices, 24 in all.

**`AmdCompacted` is the default ordering and no suite runs a whole solve through it**, the sweep
naming `AmdFlat` instead. The two return one permutation, so nothing is unchecked in the ordering;
what is unchecked is the default as a default.

**That uniformity is new on 2026-08-17 and it closed two real gaps.** `MmdChained` ran on five of
the
seven, missing the 5x5 diagonal and the complex arrow, and `AmdCompacted` had no assertion anywhere
in the
suite despite being a permanent maintained alternative. The two matrices it was missing are the ones
`experiments/ordering`'s own graphs cannot produce, since all seven of those are connected and none
is trivial: the diagonal is n isolated vertices where every degree is zero and nothing merges, and
the complex arrow exercises the structural overloads through a second instantiation.

`AMD1B` and `AMD2B` were two more such layers until 2026-08-16, when the fused schedule they carried
measured identical and faster and moved into their originals. Their 14 sameness assertions and 14
validity ones went with them, which is the drift corrected at the top of this file.

**AND THE ENUM LOST FOUR MEMBERS ON 2026-08-21.** `MMD1`, `MMD2`, `AMD1` and `AMD2` were retired to
`retired/`, taking sixteen validity assertions from this suite. What they were and why they went is
in `retired/README.md`; the short version is that the flat quotient graph was serving six drivers
across three list-order conventions, and the third one blocked aligning the two that ship.

The distinction the old suffixes carried is worth restating here because the test depended on it: a
trailing digit meant a different ordering, a trailing B or C the same ordering computed over
a different clique store. A driver is now named for its branch and its store, so the name says it
outright.

The orderings are checked for *validity*, not against 0.9's output, and not for quality. Nothing
asserts that any method reduces fill, ours included: each is a new ordering rather
than a reimplementation of a vendored one, so agreement with MMD or AMD is not a property to
assert. Fill
is measured in `experiments/ordering` instead.

Three of the seven matrices are inputs the prototypes in that experiment have never seen, and they
are the ones worth watching when a driver changes: the 5x5 diagonal is n isolated vertices, so every
degree is zero and nothing ever merges, and the tridiagonals at n = 1 and n = 2 are a single vertex
and a single edge. The experiment's seven graphs are all connected and none is trivial.

### test_forest, 29 assertions

Parent, child and sibling links, roots, height, column sizes, fundamental compression and threshold
amalgamation. Small cases with hand-computable answers: tridiagonals at n = 4 and 6 (a path, with
the last two columns merging), a 6x6 arrow (one supernode), two disjoint blocks (two trees), a star
at n = 4 (a sibling chain), and a grid. A random sweep of 200 checks links, height and sizes in both
the nodal and fundamental regimes, and that amalgamation never increases the supernode count; at
threshold 8 it merges 768 supernodes into 214. Four assertions cover `exactPatterns`, the predicate
distinguishing an index set with no stored zeros from one carrying them.

Amalgamation is greedy and not canonical, so only its tie-break-invariant properties are asserted.

Six assertions cover `setOptimizeMultifrontal`, the child reordering and postorder relabeling
ported from 0.9. The forest is built twice from one 8x8 grid under AMD, once with the option off and
once on, and the two compared. Off is asserted to be the default. On, the forest must be the same
forest, same supernode count, height, tree count, and the same multiset of front and update sizes,
since the pair moves labels and links and nothing else; the links must stay consistent in both
directions, checked by the same `validLinks` the rest of the file uses; every child list must come
out in non-increasing key order, the key being `maxStorage(jj) - updateSize(jj)^2` recomputed
from the sorted forest; and the labels must form a postorder, every subtree holding a contiguous run
ending at its own label, verified by walking each range back up the parent chain. That last one is
the assertion that matters most, because it is the property the multifrontal drivers depend on: they
loop over labels, not links, so the sort buys nothing without it. A grid is used rather than a chain
deliberately: a chain has no sibling choice to make and would assert nothing.

The sixth pins the *correspondence* between the two, that a parent's children increase in label order
as they run first to last in the link order, and it is not implied by the other five. Reverse the
direction in which `labelDepthFirst` pushes children and the links are still correctly sorted, the
labels are still a valid postorder with contiguous subtrees, and every other assertion still passes,
while the order actually realized is the reverse of the one the sort chose. Since reversing an
optimum tends toward the pessimum, that mutation would turn the saving into a loss silently. It was
checked by making exactly that mutation and confirming this assertion alone fails.

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

### test_numfactor, 18 assertions

Cholesky against dense Cholesky: a 3x3, then left-looking, right-looking and multifrontal for real
and complex, then a 10x10 grid under natural and AMD orderings in all three traversals. A
non-positive-definite 2x2 is refused rather than mangled. A random sweep of 40 matrices per
(traversal, scalar type) confirms every one factors. One assertion records that the grid case is
structurally non-trivial (forest height 90, AMD cutting supernodes from 90 to 81 and indices from 954
to 501), so the comparison is not passing on a degenerate case.

The complex multifrontal assertion is the sharpest of the Cholesky set. 0.9's complex Cholesky is
already wrong in the left- and right-looking engines (SYRK where HERK is needed), and its
*multifrontal* engine forms the contribution block with SYRK for both scalar types too, so complex
multifrontal is exactly where that bug lives. The port's operation-named `herk` wrapper resolves to
zherk, so the factor is correct, and this assertion is what proves it against the dense oracle rather
than only by residual.

LDL by reconstruction, in all three symmetries and now all three traversals (left-looking,
right-looking and multifrontal): real `L D L^T`, complex symmetric `L D L^T` with complex D, and
complex Hermitian `L D L^H` with real D. So static LDL multifrontal is checked here directly, factor
against reconstruction, the same bar the other two traversals meet. One assertion records that no
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

### test_pipeline, 86 assertions

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

**Tier 0 multifrontal and the refusal, five assertions.** Multifrontal is now the third working
traversal for every factorization, static and dynamic alike: each is the same factor as left- and
right-looking, reached by a postorder pass that carries each supernode's contribution block up a
stack to its parent, checked here by residual (all at machine precision). The three static ones and
dynamic LDL on an input that needs no pivoting are all checked here; the delayed-column path is the
business of tiers 1 and 2. One combination must still return false rather than answer: dynamic
pivoting into flat storage, which has no per-supernode index storage to grow. That refusal is as
much a part of the specification as the working cells, because a cell that begins returning a
plausible wrong answer instead of a refusal is exactly the failure a port invites.

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

A third traversal assertion runs the n = 40 matrix through multifrontal. This is the assertion that
exercises delayed columns meeting the update stack: the multifrontal driver reaches them by a third
route, carrying each supernode's contribution block up a stack rather than pulling or pushing per
ancestor. Its delay and 2x2 counts are required to match left-looking exactly (the pivoting
decisions are the same), but the residual is checked to tolerance rather than bit-for-bit, since the
assembly sums a front in a different order than the pull and the last bits legitimately differ.

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

**Complex, seven assertions.** All ten (factorization, scalar type) cells are supported as of
2026-07-19, and this section covers the dynamic ones. `StaticLDLT` on a diagonally dominant complex
symmetric band; `DynamicLDLT` on the same family with half the diagonal zeroed, asserted for residual
in both traversals and for the traversals agreeing on 2 delayed columns and 2 two-by-two pivots; and
and `DynamicLDLH` on a genuinely Hermitian
band, conjugate off-diagonals and a real diagonal, likewise asserted for residual in both traversals
and for the traversals agreeing. The counts are bounded rather than pinned, and this is the case
that proved that necessary: the two platforms disagree on them (see the tier 2 note above). The
traversals must still agree with each other, which is the claim worth asserting, since both run on
the same machine.

Each dynamic factorization then gets a multifrontal assertion, and the Hermitian one is the last cell
of the whole matrix to be filled: complex Hermitian dynamic through the update stack, the conjugating
pivot kernel reached by the multifrontal driver, on the factorization with no reference behind it.
Residual is the oracle for both (7.24e-16 symmetric, 1.83e-15 Hermitian).

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

**The facade, fourteen assertions.** The tier 0 matrix again, this time through `DirectSolver`, over
all five factorizations and all three traversals: that all fifteen are reached, and that the worst
residual is at machine precision.

Not redundant with the by-hand sweep. The facade owns both factors and chooses between them with
`dynamicPivoting()`, so it can fail to reach a combination that works. That is exactly the defect
`examples/pipeline.cpp`, since renamed `examples/example_pipeline_real.cpp`, carried: it fixed the
storage at `NumFactorStatic`, so every dynamic cell
printed "not implemented" long after it was implemented. Nothing caught it because examples are
built by `make` and never run. The `reached == 15` assertion is the guard against that shape of
error, and it is why the count is asserted separately from the residual: a silently skipped
combination would otherwise leave the worst residual looking perfect.

**Two assertions cover the ordering methods**, on the same tier 0 matrix and through the same
facade: that all eleven are reached, and that the worst residual over them is at machine precision.
Validity, which `test_order` checks, says a permutation is well formed; it does not say the rest of
the pipeline can use it. This is the assertion that says so, and it is what a new ordering method
has to pass before it counts as working. Fill is not asserted here either.

Five further assertions cover the one thing the facade decides that the by-hand caller does not: the
multifrontal child ordering is computed during `analyze`, not during `factor`, so `DirectSolver` has
to know the traversal by analysis time and passes it to `ElmForestEngine`'s constructor. That makes
`setTraversal` able to invalidate an analysis, which no other setter of its kind does. The assertions
are that `analyze` succeeds and reports itself analyzed; that the analysis *survives* a switch from
left- to right-looking, which read the same forest; that it is *invalidated* switching into
multifrontal, and again switching back out, because the forest itself differs, its children reordered
and its supernodes relabeled; and that a re-analyzed multifrontal solver still solves to machine
precision. The pair of invalidation assertions is the point: asserting only that it invalidates would
pass for a setter that invalidated unconditionally, which would make every left-right switch redo the
ordering for nothing.

Seven more cover the two forest settings the facade now exposes, supernodes and amalgamation,
available in the constructor and through setters like every other setting. Three take the
constructor route: that the defaults are fundamental supernodes and no amalgamation, matching
`ElmForestEngine`; that nodal supernodes are reachable and still solve to machine precision; and that
an amalgamating solver does too, taken with multifrontal so the two forest settings are exercised
together. Neither changes what is computed, only the block structure it is computed in, so the
residual is the right oracle for both.

The other four take the setter route, and what they pin is the invalidation. Both settings change the
forest, and the analysis *is* the forest, so both must clear it, which is asserted for each in turn
and then followed by a full re-analyze, factor and solve to machine precision. That last one matters
because clearing the flag is only half the contract: a solver that invalidates and then cannot
recompute would pass the first two assertions.

What none of these observe is whether the multifrontal ordering was actually *applied*, only that the
code path differs. Nothing exposes the forest or the stack peak through the facade, so the saving measured in
TODO is not pinned by any assertion here. That is the same gap as `maxStorage` being computed and
discarded, recorded in TODO.

**Tier 2, heavy pivoting, five assertions.** Two families, and between them they cover the two ways
pivoting gets hard: many delays, and no 1x1 pivot available at all.

**Saddle point**, `[[H, B^T], [B, 0]]` at 30 variables and 12 constraints, with **both** blocks
carrying a zero diagonal. The honest use case for an indefinite solver, and the family that delays
hardest, since a constraint column has no diagonal to pivot on and no update can give it one. Three
assertions: residual under 1e-11 in both left- and right-looking with those two agreeing on the
counts, then the same matrix through multifrontal, whose residual is checked to the same tolerance
while its counts are only bounded, not matched to left-looking. Dozens of threshold decisions under a
different summation order can tip one, which is a different valid factorization rather than an error,
so the residual is the oracle. Observed 92 delayed columns across 25 supernodes, 16 1x1 and 13 2x2
pivots, residual 9.43e-14, and multifrontal at 9.05e-14.

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

**The facade's pivot statistics and inertia, eleven assertions.** `DirectSolver` does not expose the
numeric factor, so four accessors are its only window onto what the pivot search did:
`numDelayedColumns`, `numPivots1x1`, `numPivots2x2` and `inertia`. They are checked two ways, and
neither way is the accessor checking itself.

The three counts are compared against the by-hand sweep on the tier 1 matrix, which reaches the
factor directly and whose numbers this suite already pins at 5 delayed and 4 two-by-twos. Agreement
therefore says the forwarders report what the factor holds. A fourth assertion pins the partition,
`1x1 + 2 * 2x2 == n`, which is what catches a 2x2 counted in columns rather than in blocks, the
`pivotType` encoding using 2 and 3 for the two halves. A fifth checks that a statically pivoted
factorization reports zero for all three, which is the accurate report rather than a placeholder: it
takes the diagonal in symbolic order and makes no choices.

The inertia is checked against the closed form for a grid Laplacian's eigenvalues,
`4 - 2cos(pi i / (g+1)) - 2cos(pi j / (g+1))`, which shares no code with the factorization. Three
shifts of one 8x8 grid, definite, mildly indefinite and strongly indefinite, with the negative count
matching at every shift, the three components summing to the order, and at least one shift taking a
2x2 so that case runs. One assertion checks that Cholesky, static LDL and dynamic LDL agree on the
definite matrix, which is Sylvester's law across three different factors of one matrix. Two more
check the cases it declines rather than guessing: before a factorization exists, and for a
complex-symmetric `LDL^T`, whose eigenvalues are complex and have no signs to count.

**The three shifts are chosen for conditioning, not rounded.** This grid's eigenvalues hit 2, 4 and
6 exactly by symmetry, so those shifts make `A` exactly singular, and a zero eigenvalue then lands
on whichever side rounding puts it. Shifts of 0, 1 and 3 give condition numbers of about 32, 105
and 40.

**One branch is not exercised and the mutation test says so.** A 2x2 block contributes two
eigenvalues of its trace's sign when its determinant is positive, and that has never been observed:
707 accepted blocks over 400 random indefinite matrices all had a negative determinant, and deleting
the positive case leaves this suite green. At a root it is a guarantee, bounded Bunch-Kaufman
reaching a 2x2 only after both diagonals fail their own tests, which forces `det < 0`. At a non-root
it is not, Figure 3.3's test reading `|det|`, so the branch is defensive rather than dead and a
witness would be worth having.

**The factor's three sizes, seventeen assertions.** `numNodeIdx`, `numVal` and `nnz` are answered
by `ElmForest`, `SymFactor` and both numeric factors, and forwarded by `DirectSolver`. They are
checked two ways, and the difference between the two is the point.

*Against each other.* The four classes describe one factorization from four points in the pipeline,
so they are constrained: the forest and the symbolic factor must agree on all three, a statically
pivoted factor must match what symbolic predicted, a dynamically pivoted one must exceed it, and
its index sets must exceed it by exactly the delayed columns. `numNodeIdx <= nnz <= numVal` holds
in every class.

*Against a recomputation.* Each class is also checked against the three sums recomputed in the test
from the per-supernode sizes it publishes. **This is the check that pins the formulae**, and the
relations above cannot replace it: mutation testing found that a dynamic `nnz` ignoring `delaySize`
entirely still satisfies every inequality, because a dynamic factor's own `frontSize` differs from
the symbolic one once columns have moved, so it comes out larger than predicted anyway. Only the
recomputation catches it.

**The facade's guard is tested through a refusal.** `DirectSolver` reports zero for all three before
`factor()`, which is the count for a factor that does not exist rather than a sentinel. The case
that makes the guard load-bearing is a factorization that *failed*: a Cholesky handed an indefinite
matrix analyzes, refuses, and leaves behind whatever the attempt wrote, and the facade must report
zero rather than that debris. Removing the guard passes every other assertion here and fails that
one.

Five mutations were run against this group and all five fail it: a dynamic `nnz` ignoring
`delaySize`, a `SymFactor` `nnz` dropping the update rectangle, an `ElmForest` `numVal` off by one
per supernode, a facade forwarding the wrong storage, and a facade without the guard.

## What the catalog shows

Three gaps stand out, stated as facts rather than as a plan.

**The combination space is still barely sampled.** Three orderings, five factorizations, three
traversals and two scalar types is ninety cells. `test_solve` exercises ten, all under AMD, and
`test_pipeline` adds ten more under Natural. Multifrontal now works for every factorization, static
and dynamic: Cholesky is checked against the dense-Cholesky oracle directly (real and complex,
natural and AMD, in `test_numfactor`), the static LDLs by reconstruction across all three traversals
(`test_numfactor`), and all of them end to end by residual (`test_pipeline`), the dynamic ones
including the delayed-column path at tiers 1 and 2 and complex Hermitian through the stack. Every
cell of the order-by-factorization-by-traversal matrix now resolves to a residual rather than a
refusal. MMD appears only in `test_order`, where its output is checked for validity and never used to
factor anything.

**One refusal is asserted.** `test_pipeline` covers the one reachable unsupported combination:
dynamic pivoting into flat storage, which has no per-supernode index storage for a delayed column to
grow. Static factorization into dynamic storage through multifrontal is the one remaining cell that
still returns false without being asserted to; it is a convenience path, not a missing capability,
since the same factor is reached through flat storage.

**Dynamic pivoting is now covered at all three tiers and all three traversals.** Delayed columns
crossing a forest, 2x2 pivots in the solve, the expansion of a parent's front, and, with
multifrontal, delayed columns meeting the update stack are all exercised: mildly at tier 1 with the
counts pinned, and hard at tier 2 where the saddle point delays dozens of columns.

**Nothing tests the pivot threshold.** It defaults to 0.1 and it directly controls how much
delaying happens, so every dynamic assertion silently tests one value of it.
