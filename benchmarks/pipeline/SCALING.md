# Oblio: Scaling

This report asks how a solve's cost **grows**, on two families of generated grid, across every
factorization Oblio has and all three value types. Its companions ask different questions on real
matrices: `ACCURACY.md` how good the answers are, `PERFORMANCE.md` what they cost at a given size.
Both live in `benchmarks/matrices/`; this one is in `benchmarks/pipeline/` beside the driver that
produced it.

Generated grids are used here precisely because real matrices cannot answer a growth question: two
real matrices of the same order can differ in fill by orders of magnitude, so a curve through them
measures the sample rather than the algorithm. A grid ladder holds the structure fixed and moves
one parameter.

Four results, over six square grids to n = 360000 and six cubic grids to n = 110592:

- **Factorization time tracks the fill almost exactly** for Cholesky and static LDL, at an
  exponent of 0.94 to 1.04 against nnz(L) in both families. The implementation is not falling
  behind the arithmetic as problems grow.
- **The two families are different problems, not one problem at different sizes.** Fill grows as
  $n^{1.17}$ in 2D and $n^{1.56}$ in 3D, and every conclusion below splits on that.
- **The multifrontal traversal wins by 2.1x in 2D and by 1.13x in 3D.** The advantage is real in
  both and its size is not transferable between them.
- **Dynamic pivoting costs 1.5x in 2D and 7x in 3D, with nothing to pivot.** No pivot was
  perturbed and no column delayed on any matrix in this report, so that is the price of the
  search alone, and it is the one place where cost outruns fill.

## The environment

| | |
|---|---|
| Machine | Apple M4, 32 GB, macOS 26.6.1 |
| Compiler | Apple Clang (Xcode 26.6), C++17, `-O3 -DNDEBUG` |
| BLAS and LAPACK | Apple Accelerate |
| Timing | best of three after a warm-up, per phase |

Single-threaded throughout. On Apple Silicon the benchmark asks the scheduler for a performance
core before timing anything, since a command-line process may otherwise be parked on an efficiency
core for a whole run.

## What was measured

| | |
|---|---|
| Families | square grids, 100 to 600 a side; cubic grids, 16 to 48 |
| Orderings | multiple minimum degree (MMD) and approximate minimum degree (AMD) |
| Traversals | left-looking, right-looking, multifrontal |
| Factorizations | eight, across three value types |

The eight are every combination that exists:

| values | factorizations |
|---|---|
| real | Cholesky, static $LDL^T$, dynamic $LDL^T$ |
| complex Hermitian | Cholesky, static $LDL^H$, dynamic $LDL^H$ |
| complex symmetric | static $LDL^T$, dynamic $LDL^T$ |

The two absences are not gaps. Over the reals $LDL^H$ **is** $LDL^T$, the conjugate being the
identity, so running both would repeat a column. And Cholesky requires Hermitian positive
definiteness, so it does not apply to a complex symmetric matrix at all.

**This is the only place complex Cholesky can be measured.** The SuiteSparse Matrix Collection has
23 complex square matrices with symmetric pattern, of which exactly one is marked positive
definite, and that one is complex symmetric rather than Hermitian. A generated grid has no such
scarcity.

### The matrices, and why nothing pivots

All three arms share one pattern, so a difference between them is arithmetic and never structure.
Each is **diagonally dominant by 25 percent**: the diagonal is 1.25 times the sum of the
off-diagonal magnitudes in its row.

| arm | diagonal | off-diagonal |
|---|---|---|
| real | 1.25 * degree | -1 |
| complex Hermitian | 1.25 * degree * sqrt(1.25), real | -1 + 0.5i below, conjugate above |
| complex symmetric | the same, real | -1 + 0.5i on both sides |

The margin is deliberate and it is what makes the report interpretable: **no pivot is ever too
small to divide by and no column is ever unusable**, in any arm. Static LDL therefore never
perturbs and dynamic LDL never delays, which is confirmed on every matrix of both ladders and
reported by the run itself. The dynamic-against-static comparison below is then a clean measurement
of the pivot search, with nothing found.

The complex off-diagonal carries a genuine imaginary part rather than a zero one, so the complex
kernels do complex arithmetic rather than multiplying zeros in a complex type.

**Accuracy is not the subject here and was checked anyway.** Every one of the 96 factorizations
reported returned a normwise backward error between 4.9e-16 and 9.1e-15, which is machine
precision throughout, with relative residuals from 2.1e-15 to 9.1e-14. Nothing perturbed, nothing
delayed. `ACCURACY.md` is where that question is asked properly, on matrices chosen to be hard.

## How the problem grows

Fill, from the matrices themselves:

| | square | cubic |
|---|---|---|
| smallest | n = 10000, nnz(L) = 186835 | n = 4096, nnz(L) = 295113 |
| largest | n = 360000, nnz(L) = 12234021 | n = 110592, nnz(L) = 51113859 |
| growth | $n^{1.17}$ | $n^{1.56}$ |

**A cubic grid of 110592 columns fills four times more than a square grid of 360000.** That single
fact is why both ladders exist and why no result below is quoted without saying which family it
came from.

## Time against fill, which is the question that matters

Growth exponents, first rung to last:

**Square grids:**

| phase | against n | against nnz(L) |
|---|---|---|
| ordering | 0.70 | |
| analysis | 0.87 | |
| Cholesky, left-looking | 1.10 | **0.94** |
| static LDL | 1.11 | **0.95** |
| dynamic LDL | 1.21 | 1.04 |

**Cubic grids:**

| phase | against n | against nnz(L) |
|---|---|---|
| ordering | 1.02 | |
| analysis | 1.22 | |
| Cholesky, left-looking | 1.62 | **1.04** |
| static LDL | 1.59 | **1.01** |
| dynamic LDL | 2.05 | **1.31** |

**Read the nnz(L) columns.** Cholesky and static LDL cost very nearly one unit of time per unit of
fill, in both families, across a factor of 65 in problem size and 173 in fill. That is the result a
scaling report exists to establish: **the cost is the arithmetic and not the bookkeeping**, and
nothing in the implementation degrades as the problem grows.

**The analysis grows more slowly than the factorization**, at 0.87 against 1.10 in 2D and 1.22
against 1.62 in 3D. Its share therefore falls with size, which is the opposite of what
`PERFORMANCE.md` finds on real matrices, where the split stays near constant. Grids grow their fill
faster than real matrices of comparable order do.

**Dynamic LDL is the exception and the finding**, at 1.31 against fill in 3D. It is discussed
below.

## The traversal, and why the answer differs by family

Numeric factorization, real Cholesky, milliseconds:

| grid | left-looking | right-looking | multifrontal | MF advantage |
|---|---|---|---|---|
| square 400 | 79.87 | 74.85 | 38.25 | 2.09x |
| square 600 | 195.68 | 188.98 | 92.67 | 2.11x |
| cube 40 | 304.68 | 373.13 | 262.87 | 1.16x |
| cube 48 | 819.24 | 983.13 | 722.12 | 1.13x |

**Multifrontal is the fastest traversal in both families and by very different margins**, slightly
over 2x in 2D and slightly over 1.1x in 3D. `PERFORMANCE.md` measures 1.33x on real positive
definite matrices, which sits between the two.

The mechanism explains the split. The multifrontal traversal assembles each supernode's
contributions into a dense frontal matrix and hands that to the BLAS, so more of the work happens
in `dgemm` and `dsyrk` on contiguous memory. That is worth most where the alternative's indirection
dominates, which is where fronts are small and numerous, and 2D grids are exactly that. In 3D the
fronts are large enough that left-looking is already spending its time inside dense kernels, and
the assembly buys less.

**Right-looking is also family-dependent**: slightly faster than left-looking in 2D, and 20 percent
slower in 3D.

## What complex arithmetic costs

Real against complex Hermitian, same pattern, same margins, left-looking:

| grid | real | complex Hermitian | ratio |
|---|---|---|---|
| square 600, Cholesky | 195.68 | 536.61 | 2.74x |
| square 600, static LDL | 279.52 | 856.47 | 3.06x |
| cube 48, Cholesky | 819.24 | 3690.37 | 4.50x |
| cube 48, static LDL | 1064.44 | 4703.57 | 4.42x |

A complex multiply-add is four real multiplies and two adds, so **4x is the arithmetic
expectation**, and the measurements straddle it: below in 2D and at or slightly above in 3D. A
complex problem does four times the arithmetic on twice the memory, so it is relatively less
memory-bound, and the ratio approaches the flop count as fronts grow and the kernels become
compute-bound. In 2D it never gets there.

**Complex Hermitian and complex symmetric cost the same**, within a percent at every size, which is
the expected result and a useful check: the two arms have identical magnitudes and differ only in
whether the kernel conjugates.

## What dynamic pivoting costs when there is nothing to pivot

Dynamic against static LDL, left-looking, **zero delayed columns everywhere**:

| grid | static | dynamic | ratio |
|---|---|---|---|
| square 100 | 5.17 | 5.42 | 1.05x |
| square 400 | 114.89 | 153.59 | 1.34x |
| square 600 | 279.52 | 413.07 | 1.48x |
| cube 16 | 5.71 | 8.66 | 1.52x |
| cube 32 | 141.22 | 568.53 | 4.03x |
| cube 48 | 1064.44 | 7426.54 | **6.98x** |

**This is the price of the pivot search alone**, since not one column was ever delayed and the
factor is byte-for-byte the size the analysis predicted. In 2D it is a manageable and slowly
growing overhead. In 3D it reaches a factor of seven and is still climbing: dynamic LDL grows as
$nnz(L)^{1.31}$ where static grows as $nnz(L)^{1.01}$.

The mechanism is that the search cost depends on the front, not on the fill. A threshold test
examines candidate pivots within a supernode and takes column maxima to do it, and 3D grids have
far larger fronts than 2D grids at comparable fill. So the search is charged per front-column
against a front that is growing, while the factorization itself is charged per entry.

**This is a real cost and it is a cost worth paying when it is needed.** `ACCURACY.md` records
`IPSO/OPF_6000`, where static LDL returns a backward error of 1 and dynamic returns 3.5e-19, a
difference of seventeen orders of magnitude on a matrix static factorization cannot handle. The
point of measuring it here, on matrices that need none of it, is that **the two are separable**: a
caller who knows the matrix is definite should ask for Cholesky or static LDL and not pay this,
and Oblio lets them.

## The two orderings

AMD relative to MMD, on fill and on analysis time:

| grid | fill | ordering time | analysis time |
|---|---|---|---|
| square 200 | +10.2% | 1.06x | 1.01x |
| square 400 | +16.5% | 1.15x | 1.10x |
| square 600 | +16.5% | 1.15x | 1.06x |
| cube 26 | -1.1% | 0.46x | 0.66x |
| cube 32 | -1.9% | 0.48x | 0.72x |
| cube 48 | +2.3% | 0.44x | 0.79x |

**On square grids MMD fills 10 to 18 percent less; on cubic grids the two are within a few percent
of each other**, in both directions. And **AMD orders cubic grids in less than half the time**,
consistently, while the two cost about the same in 2D.

Both of these are worth setting beside `PERFORMANCE.md`, where the fill difference on 107 real
matrices is one to three percent and AMD is slightly ahead. Grids, and square grids in particular,
are where the fill gap between minimum degree variants lives: nearly every live vertex has the same
degree, so the tie-break decides almost every pick. **A fill claim measured on a square grid should
not be carried to a cubic one, let alone to a real matrix.**

## What this report does not establish

- **Two structures only**, and highly regular ones. Grids answer a growth question precisely
  because they are uniform, and that is also their limitation. `PERFORMANCE.md` covers 107 real
  matrices at a fixed size.
- **Sizes to n = 360000 in 2D and n = 110592 in 3D**, with factors to 51 million entries. The next
  rungs are documented but not run.
- **A single machine, a single BLAS, no parallelism.** Traversal and complex-arithmetic results in
  particular depend on the BLAS, since both turn on how much work reaches dense kernels.
- **Matrices constructed to need no pivoting.** That is what makes the dynamic-against-static
  measurement clean, and it means nothing here says what dynamic pivoting costs when it is doing
  its job.
- **An exponent fitted from two endpoints** rather than a regression over the ladder. The
  intermediate rungs are in the raw output and are consistent with it, but this is a shape and not
  a precise constant.
- **One right-hand side per matrix and no iterative refinement.**

## Reproducing this

The benchmark sits in `benchmarks/pipeline/`, whose README carries the working notes and the two
longer ladders that were priced but not run. From that directory:

```
make scale2d
make scale3d

./pipeline_scaling_cpp 3d 16 20 26     # any sides
```

Nothing is downloaded and nothing is stored: the matrices are generated, so the whole report is
reproducible from the repository alone.
