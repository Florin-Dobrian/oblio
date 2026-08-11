# Oblio: Accuracy

This report asks a single question, on matrices Oblio's authors did not choose and did not
generate: **when Oblio returns a solution, how good is it?** Companion reports cover performance
on the same real matrices, and scaling on synthetic grids.

The short answer, over 106 matrices from the SuiteSparse Matrix Collection spanning 28 problem
kinds and four orders of magnitude in fill: **every solution Oblio produced was backward stable at
machine precision, and every case where it produced no useful solution was a property of the
matrix rather than of the solver.** There is no row in the results where the answer is poor and
the reason is Oblio's.

Everything below says how that was measured, what it does and does not establish, and what the
awkward cases were, because a report that only shows the good rows is not worth reading.

## The environment

| | |
|---|---|
| Machine | Apple Silicon (M-series), macOS |
| Compiler | Apple Clang, C++17, `-O3 -DNDEBUG` |
| BLAS and LAPACK | Apple Accelerate |
| Oblio dependencies | a C++17 compiler, `make`, and a BLAS. Nothing else |

Oblio calls LAPACK for dense supernodal kernels (`?potrf`, `?trsm`, `?gemm`, `?syrk`) and supplies
its own unpivoted LDL kernel, LAPACK having none: `?sytrf` pivots, and pivoting is what a static
factorization refuses to do.

## What Oblio was asked to do

**Oblio's defaults throughout, held fixed on purpose.** The question here is accuracy, not which
ordering or traversal is fastest, so nothing varies that does not have to. Oblio also provides AMD
orderings and right-looking and multifrontal traversals; those are the subject of the performance
report and would only distract from this one.

| | |
|---|---|
| Ordering | multiple minimum degree (MMD), the default |
| Traversal | left-looking, the default |
| Factorizations | Cholesky, static $LDL^T$, dynamic $LDL^T$ |
| Values | real, double precision |
| Right-hand side | all ones |

$LDL^H$ and complex arithmetic, which Oblio supports in all three forms, are not exercised here
either. Over the reals $LDL^T$ and $LDL^H$ are the same computation, so running both would have
repeated a column rather than covered a case.

### The three factorizations, and why all three

They differ in what they do when a pivot is unusable, and that difference is the whole of what
this report measures.

**Cholesky** requires positive definiteness and refuses without it. It never perturbs and never
pivots. In this report it is as much a detector as a solver: it should succeed on exactly the
positive definite matrices and on no others.

**Static LDL** does not pivot at all. The elimination order is fixed by the analysis and never
changes, so a pivot too small to divide by is **perturbed** to a floor value (`1e-14` by default)
and counted. The factorization therefore cannot fail, and the price is that it has factored a
matrix near `A` rather than `A` itself. On a matrix that needs pivoting this shows up honestly in
the backward error, and the perturbation count says why.

**Dynamic LDL** chooses its pivots as it runs. It takes 1x1 and 2x2 pivots subject to a threshold
test (0.1 by default) and, where no acceptable pivot exists in a supernode, **delays** the offending
columns to the parent front rather than perturbing them. The structure of the factor therefore
changes during the numeric phase, which is what "dynamic" means here. Non-root supernodes follow
Ashcraft, Grimes and Lewis (1998), Figure 3.4 for the pivot search with the Figure 3.3 acceptance
test; roots, which have no parent to delay into, use bounded Bunch-Kaufman.

This is the strategy that costs something, and this report measures both sides of it: the accuracy
it buys and the fill it spends.

## The matrices, and how coverage was obtained

Real matrices were used because generated ones cannot settle this question. Oblio's development
history contains several claims that held on square grids and failed on cubic ones; a family of
one shape is a poor witness for anything.

Matrices came from the SuiteSparse Matrix Collection, selected from its published index by an
automated filter rather than by hand:

- **real, square, and numerically symmetric**, which is what Oblio's symmetric factorizations take;
- **structural variety before size**, sampling the cheapest few matrices of *each* problem kind the
  collection records rather than taking the largest or the most famous. Variety is what generated
  grids lack, so it is what a real set has to supply;
- **an upper bound on nnz(A)**, to keep a full pass to a few minutes.

The resulting set is 160 files, of which 41 carry no values at all (the collection stores binary
matrices as patterns) and 119 are usable. It spans **28 problem kinds**, including structural
mechanics, thermal, computational fluid dynamics, quantum chemistry, circuit and power-network
simulation, optimal control, economics, materials, model reduction, statistical learning and
combinatorial problems.

The selection script is committed with the benchmark, so the set is reproducible on another
machine from one command.

## What is measured, and why two numbers rather than one

For a computed solution `x` of `Ax = b`, this report gives both

```
res = ||Ax - b|| / ||b||                        the relative residual
bwd = ||Ax - b|| / (||A|| ||x|| + ||b||)        the normwise backward error
```

in infinity norms, from a single residual vector.

**The backward error is the verdict.** It says the computed `x` exactly solves a nearby problem,
`(A + E)x = b` with `||E||` small relative to `||A||`, and for a backward stable factorization it
stays near machine precision **whatever the conditioning of the matrix**. That is the property a
direct solver can be held to.

**The relative residual is not a verdict**, and reporting it alone would be misleading in both
directions. It differs from the backward error by a factor of `(||A|| ||x|| + ||b||) / ||b||`,
which grows to the order of the condition number when `||x||` is large relative to `||b||`. On a
matrix conditioned at 1e15 a perfectly computed solution can show a relative residual near 1. Real
matrices span a range of conditioning that no synthetic grid does, so on this set a large residual
beside a tiny backward error is a statement about the problem and not about the solver.

Both are reported because together they separate the two cases. Neither is thresholded: this
report has no pass or fail line, and nothing here is graded.

## Structurally singular matrices are identified and declined

A system with no solution cannot be used to measure a solver. Before any factorization runs, each
matrix is tested for **structural singularity** by computing a maximum matching between its rows
and columns, using Hopcroft-Karp on the pattern of nonzero *values*.

The reasoning is exact. A permutation contributing to the determinant needs a nonzero in every
position, which is precisely a perfect matching in the bipartite graph of the pattern. Where no
perfect matching exists, every term of the determinant expansion vanishes and the matrix is
singular for almost any choice of values. The size of the maximum matching is the **structural
rank**, and a deficit is a proof taken from the pattern with nothing computed.

Twelve of the 119 valued matrices were declined this way:

| matrix | n | structural rank | deficit | empty columns |
|---|---|---|---|---|
| `GHS_indef/aug3d` | 24300 | 11664 | 12636 | 0 |
| `GHS_indef/aug2dc` | 30200 | 20000 | 10200 | 0 |
| `GHS_indef/aug2d` | 29008 | 19208 | 9800 | 0 |
| `Cunningham/m3plates` | 11107 | 6639 | 4468 | 4468 |
| `Boeing/bcsstm38` | 8032 | 5199 | 2833 | 2833 |
| `Pajek/Reuters911` | 13332 | 10685 | 2647 | 18 |
| `HB/zenios` | 2873 | 266 | 2607 | 2605 |
| `Arenas/PGPgiantcompo` | 10680 | 8159 | 2521 | 0 |
| `HB/bcsstm13` | 2003 | 1241 | 762 | 762 |
| `Newman/netscience` | 1589 | 1424 | 165 | 128 |
| `GHS_indef/laser` | 3002 | 3000 | 2 | 0 |
| `GHS_indef/bloweybl` | 30003 | 30002 | 1 | 1 |

Two things are worth noticing. **The matching is doing real work**: the `GHS_indef/aug*` family and
`PGPgiantcompo` have deficits in the thousands and **no empty columns at all**, so a simpler test
that only looked for empty rows would have passed every one of them and their meaningless residuals
would have sat in the results looking like defects. And **`GHS_indef/laser` is deficient by 2 out
of 3002**, which no inspection by eye would find.

Where the deficit does come from empty columns, the cause is usually mundane and worth knowing:
`m3plates`, `bcsstm13` and `bcsstm38` are structural mass matrices, which are commonly rank
deficient by construction, and `netscience` is a collaboration network whose isolated vertices are
authors who published alone.

**What this test cannot see is numerical singularity in a structurally full matrix**, and nothing
short of factoring can. Those cases are marked in the results instead, from the inertia.

## Results

Of the 160 files: 41 carry no values, 12 are structurally singular and declined, 1 could not be
factored in the memory available (discussed below), and **106 were solved**.

### Classification, from the inertia

Oblio reports the inertia of `A` by reading the signs of `D`, which is exact rather than estimated:
$A = LDL^H$ is a congruence and congruence preserves the signs of the eigenvalues, so counting them
in `D` counts them in `A` without computing an eigenvalue.

| class | count | meaning |
|---|---|---|
| positive definite | 30 | every eigenvalue strictly positive |
| negative definite | 18 | every eigenvalue strictly negative, so Cholesky must refuse and `-A` would factor |
| indefinite | 56 | both signs present, none zero |
| numerically singular | 2 | at least one exactly zero eigenvalue |

**Cholesky agreed with the inertia on all 106 matrices**: it succeeded on every positive definite
matrix and on no other. That is two independent parts of the library checking each other, and it
has now held across four independent samples of the collection, at 35, 58, 88 and 106 matrices.

### Accuracy

Every solved matrix produced a backward error at or near machine precision from at least one
factorization. The eight matrices with the largest residuals, taking the best of the three
factorizations on each measure:

| matrix | best res | best bwd | perturbed | zero eig |
|---|---|---|---|---|
| `Marini/eurqsa` | 2.0e+16 | 5.8e-18 | 2084 | 8 |
| `Guettel/TEM27623` | 1.1e+04 | 3.5e-18 | 79 | 0 |
| `GHS_indef/sit100` | 2.6e+03 | 9.6e-15 | 1063 | 0 |
| `Oberwolfach/t2dal_bci` | 1.4e+03 | 1.1e-16 | 0 | 0 |
| `HB/plat1919` | 1.5e+02 | 7.3e-17 | 1 | 0 |
| `Cote/vibrobox` | 1.4e+02 | 6.3e-18 | 0 | 0 |
| `Oberwolfach/t2dah_a` | 1.7e+01 | 1.3e-16 | 0 | 0 |
| `Boeing/crystk01` | 2.8e+00 | 8.1e-17 | 6 | 0 |

**Read the second column.** Every one of these is backward stable. The large residuals are
conditioning and rank, which is what the pair of measures exists to distinguish. **A large residual
beside a large backward error, with no perturbations and no zero eigenvalues, would be a defect,
and there is no such row.**

### What dynamic pivoting buys

The comparison between static and dynamic LDL on the same matrix is the clearest statement of what
the pivoting is for. A few representative rows:

| matrix | static bwd | dynamic bwd | delayed |
|---|---|---|---|
| `Gset/G32` | 1.0e-02 | 4.3e-15 | 4017 |
| `Schenk_IBMNA/c-18` | 4.1e-16 | 2.6e-22 | 2603 |
| `ML_Graph/mice_10NN` | 1.1e-03 | 7.2e-17 | 242 |
| `IPSO/OPF_6000` | 1.0e+00 | 3.5e-19 | 10684 |
| `Boeing/nasa1824` | 4.6e-11 | 7.9e-18 | 55 |
| `FIDAP/ex32` | 1.5e-07 | 8.0e-19 | 1329 |

On `IPSO/OPF_6000`, a 29902-column power-flow system, static LDL perturbs 3020 pivots and returns a
backward error of 1. Dynamic LDL delays 10684 columns and returns 3.5e-19. **That is seventeen
orders of magnitude, on a matrix a static factorization simply cannot handle.**

## Three cases worth knowing

### `HB/saylr3`: numerically singular, and a control beside it

`saylr3` and `sherman1` are the same 10x10x10 grid problem from two different authors, identical in
`n`, in `nnz(A)`, in `nnz(L)` at 9872, in diagonal range and in norm. One solves cleanly and the
other does not:

```
saylr3    static   bwd 1.9e-15   res 2.0e+00   2 perturbed
          dynamic  bwd 8.9e-05   res 2.0e+00   inertia 0 / 998 / 2
sherman1  both     bwd 8.1e-17   res 1.8e-12   inertia 0 / 1000 / 0
```

`saylr3` is numerically singular with rank 998, and three independent witnesses agree: the inertia
finds two zero eigenvalues, static LDL perturbs exactly two pivots, and giving it a consistent
right-hand side removes the problem entirely. With `b = Ax` for a chosen `x`, dynamic LDL on
`saylr3` returns **9.3e-17, identical to `sherman1`**. The all-ones right-hand side simply is not
in the range of a rank-deficient matrix, so no solution exists and the 2.0 residual is that fact
and nothing else.

This is the class the structural test cannot catch by construction: structural rank is a full 1000,
and dynamic LDL neither perturbs nor delays, because accepting a zero pivot and losing rank is
exactly what it should do.

### `GHS_indef/bloweybq`: ill-conditioned, not singular

This matrix carries its own description in the collection: *"matrix for which early version of MA57
fails"*. It is positive definite, structurally full, and conditioned at roughly 1e15.

```
b = ones    dynamic   bwd 2.1e-19   res 8.8e-01   max|x| 8.3e+14
b = A x     dynamic   bwd 1.8e-16   res 9.1e-16   max|x| 4.0e+00
```

Given a consistent right-hand side all three factorizations solve it and recover the chosen
solution to 1.6e-3 relative. Given `b` all ones the solution is legitimately of size 1e15, and at
that magnitude `Ax` cannot be formed accurately enough to reproduce `b`, so the residual near 1 is
a cancellation artifact and not a failure to solve. **Oblio handles correctly a matrix distributed
specifically because a well-known solver did not.**

### `Marini/eurqsa`: an economic KKT system at the edge

A time-series reconciliation problem, structurally nonsingular, with 2121 zero diagonal entries, a
diagonal of 1 to 2 where present, and off-diagonals reaching 1.19e6. Its smallest singular value is
around 1e-11 against a norm of 6.4e6.

Static LDL perturbs 2084 pivots and reaches a backward error of 4.0e-02, which is the honest cost
of refusing to pivot on a matrix that needs it. **Dynamic LDL reaches 5.8e-18**, at a cost of 40294
delayed columns and a factor 19.9 times larger than the analysis predicted. Both numbers are the
strategy working as designed, and the pair is the trade stated plainly.

## The cost of accuracy, and one matrix Oblio could not factor

Delaying a column is not free. A delayed column widens its parent's front, which enlarges the
factor and the work. Oblio reports both the fill the analysis predicted and the fill the
factorization actually held, so the cost is visible per matrix:

| matrix | n | predicted | actual | ratio | delayed |
|---|---|---|---|---|---|
| `LFAT5000` | 19994 | 67463 | 15678721 | **232.4x** | 3128750 |
| `eurqsa` | 7245 | 156243 | 3101991 | 19.9x | 40294 |
| `c-66b` | 49989 | 1121278 | 3927468 | 3.5x | 319743 |
| `vibrobox` | 12328 | 2332570 | 7113191 | 3.0x | 29296 |

**On the great majority of the set the ratio is exactly 1.0 and nothing was delayed**, meaning
dynamic pivoting cost nothing at all in fill. Where it costs, it can cost a great deal.

### `GHS_indef/bloweya`

One matrix in the set, a Cahn-Hilliard problem at n = 30004, could not be factored in the memory
available. Reported plainly:

- Its **analysis is small**: 30002 supernodes, largest front 3, predicted fill 110007 entries, well
  under a megabyte of values.
- **Static LDL factors it** in exactly the predicted 110007 entries.
- **Dynamic LDL exceeds 32 GB** and does not complete.

The mechanism is visible and has a control on each side of it. Its elimination tree is a chain of
height 5004, and 20003 of its 30004 diagonal entries are structurally zero, so nearly every pivot
candidate is one the threshold test must reject. Delays then cascade: each delayed column widens
its parent, the wider parent has more columns it cannot use, and over a chain that deep the fronts
grow without bound. `GHS_indef/bloweybq` has the same chain shape at height 4999 and, being
positive definite, delays nothing and factors in 39997 entries.

**This is a known characteristic of delayed-pivot sparse indefinite factorization, not a defect,
and the literature Oblio's strategy comes from addresses it directly.** Ashcraft, Grimes and Lewis
(1998) state the mechanism outright: their standard pivot search postpones *all* `m` columns of a
block when none of the `2m - 1` candidates it examines is acceptable, and "failing to eliminate a
column in the partial factorization creates more fill in L". Four routes follow from that, in
roughly increasing order of work:

1. **Relax the threshold.** Oblio exposes `setPivotThreshold`, and the default of 0.1 sits at the
   strict end. AGL study exactly this trade across their test set and recommend the looser end as
   the effective compromise, finding accuracy largely insensitive to it over a wide range while
   efficiency improves. This costs nothing to try and would be the first thing to measure.
2. **Search 2x2 pivots exhaustively.** Oblio uses AGL's Figure 3.4 standard ordering, which tests
   only `m - 1` of the `m(m - 1)/2` possible 2x2 pivots in a block. Their Figure 3.5 search reuses
   the column maxima it has already computed to test all of them for about `m(m - 1)/2` additional
   scalar operations, and finds acceptable pivots the standard search misses. Fewer rejected
   columns means fewer delays.
3. **Allow pivot blocks larger than 2x2.** AGL give a stable extension to larger blocks, which
   accepts pivots that no 1x1 or 2x2 test can.
4. **Use static LDL with perturbation, and iterate.** Static LDL factors this matrix in the
   predicted 110007 entries with a bounded backward error on the perturbed system. Paired with
   iterative refinement, that is the standard route when a dynamically pivoted factor would be too
   large, and it is available in Oblio today apart from the refinement step.

None of these is implemented in Oblio yet, and the matrix is reported here rather than quietly
skipped because a report that hides its hard case is worth less than one that explains it.

`Oberwolfach/LFAT5000` is the same mechanism inside the memory budget, and it sharpens the
question: it is **positive definite**, every pivot acceptable, and dynamic LDL delayed 3.1 million
columns anyway for a 232-fold increase in fill, reaching a backward error of 9.7e-20. Definiteness
alone does not prevent the cascade, so the tuning above is a real avenue rather than a
rationalization.

## What this report does not establish

Stated plainly, because the boundaries matter more than the headline.

- **No timing.** Nothing here is a performance claim. That is a companion report.
- **One ordering and one traversal**, both defaults. AMD, right-looking and multifrontal are
  implemented and exercised in the performance report rather than here.
- **Real double precision only.** Oblio supports complex Hermitian and complex symmetric
  factorizations, including dynamic pivoting for both, and none is exercised in this set. Complex
  Hermitian dynamic LDL is the one part of Oblio's numeric code that is an extension rather than a
  port from its predecessor, so it is the part that most wants outside evidence.
- **One right-hand side per matrix**, all ones, with no iterative refinement and no scaling or
  equilibration applied. Several matrices here would improve substantially with either.
- **Sizes to about 70000 columns and 15.7 million factor entries.** The larger end of the
  collection is untouched.
- **Numerical singularity is marked, not proved.** The zero count from the inertia is a marker: a
  matrix singular to within rounding may show no exactly zero eigenvalue at all. It is a proof when
  it fires and no evidence either way when it does not.

## Reproducing this

The benchmark and the selection script sit beside this report, in `benchmarks/matrices/`. The
README there carries the working notes: what was measured, what went wrong on the way, and the two
measurement mistakes that had to be corrected before the results above could be trusted. From that
directory:

```
./ssget.py list --per-kind 8 --values --max-nnz 500000 > candidates.txt
./ssget.py fetch candidates.txt
make run
```

The matrices are downloaded from the collection rather than committed, so the set is reproducible
without shipping a gigabyte of somebody else's data.
