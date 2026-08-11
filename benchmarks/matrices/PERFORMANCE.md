# Oblio: Performance

This report asks what a solve costs, on matrices Oblio's authors did not choose and did not
generate: **where does the time go, and which choices move it?** Its companions ask the other two
questions: `ACCURACY.md` how good the answers are, on the same kind of input, and
`benchmarks/pipeline/SCALING.md` how the cost grows, on generated grids.

Three results, over 107 positive definite matrices from the SuiteSparse Matrix Collection:

- **The multifrontal traversal is the fastest numeric factorization, by roughly a third.** Against
  the best of the three on each matrix it averages 1.03 where left-looking averages 1.36 and
  right-looking 1.39.
- **The choice between multiple minimum degree and approximate minimum degree barely matters for
  fill**, one to three percent, on matrices with real structure. On generated grids the same
  comparison reaches thirteen percent, which turns out to be a property of grids rather than of
  the orderings.
- **The analysis is 45 percent of a one-shot solve at the median**, against 47 percent for the
  numeric factorization, and that split holds across three orders of magnitude of factor size.

Everything below says how those were measured and what they do not establish.

## The environment

| | |
|---|---|
| Machine | Apple M4, 32 GB, macOS 26.6.1 |
| Compiler | Apple Clang (Xcode 26.6), C++17, `-O3 -DNDEBUG` |
| BLAS and LAPACK | Apple Accelerate |
| Timing | best of three after a warm-up, per phase |

On Apple Silicon the benchmark asks the scheduler for a performance core before timing anything.
A command-line process runs at the default quality-of-service class, which prefers a performance
core but permits the scheduler to park the thread on an efficiency one, and that placement is
**sticky over long stretches** rather than jittering per sample, so a minimum over repeats does not
filter it. In a sibling benchmark this one call turned a scattered null result into a clean one.

## What was measured

| | |
|---|---|
| Factorization | Cholesky, real, double precision |
| Orderings | multiple minimum degree (MMD) and approximate minimum degree (AMD) |
| Traversals | left-looking, right-looking, multifrontal |
| Phases timed | ordering, analysis, numeric factorization, solve |

**Cholesky only, and that is the point rather than a limitation.** Cholesky never pivots, which
buys three things at once:

- the factor's structure is exactly what the analysis predicted, so **nnz(L), the number of
  nonzeros in the factor, is a property of the ordering alone** and is comparable across
  orderings;
- no column is delayed, so no row can surprise the memory budget;
- and the numeric phase does the same arithmetic under every traversal, so a difference between
  left-looking, right-looking and multifrontal is a difference in **how the work is scheduled**
  rather than in how much of it there is.

**The phases are separated so that no cost is counted twice.** `order` is the ordering computation
alone. `analyze` is the whole analysis, so the elimination forest and symbolic factorization cost
the difference between them. The numeric factorization is timed once per traversal, with the
analysis already done and outside the timed region, and the analysis is timed separately for
multifrontal because that traversal sorts each supernode's children and relabels them into a
postorder.

## The matrices

Cholesky requires positive definiteness, so the set is every real symmetric positive definite
matrix the collection marks as such in a size range: 1000 to 100000 columns and at most two
million nonzeros. That is 114 matrices, taken **whole rather than sampled**, since the pool is
small enough that sampling would lose more than it saves.

Of the 114: **107 measured**, one was not in fact positive definite, five were skipped for
predicted fill above the benchmark's cap, and one carried no values.

The set spans 20 problem kinds, among them structural mechanics, computational fluid dynamics,
thermal, model reduction, optimization and materials. Factor sizes run from a thousand entries to
**46 million**, and multifrontal factorization times from under a tenth of a millisecond to 1.1
seconds.

**Accuracy is not the subject here and was checked anyway.** Across all 428 factorizations, four
orderings on each of 107 matrices, the worst normwise backward error was **4.1e-15**, which is
machine precision. Nothing perturbed a pivot and nothing delayed a column, as Cholesky on a
positive definite matrix never does either. Relative residuals range far more widely, up to 1.5e+03
on the worst-conditioned matrix in the set, which is conditioning rather than error and is the
subject of `ACCURACY.md`, where the distinction between the two measures is drawn properly.

## Results

### The traversal is the decision that matters

Geometric mean of the numeric factorization time relative to the best traversal on each matrix, at
a fixed ordering:

| traversal | relative time |
|---|---|
| left-looking | 1.364 |
| right-looking | 1.388 |
| **multifrontal** | **1.028** |

**Multifrontal is close to the best on nearly every matrix and the other two are close to a third
behind.** The result held unchanged when the set grew from 87 matrices to 107.

The reason is where the arithmetic happens. All three traversals do the same operations; the
multifrontal one assembles each supernode's contributions into a dense frontal matrix and hands
that to the BLAS, so more of the work happens inside `dgemm` and `dsyrk` on contiguous memory. The
other two accumulate updates into the factor itself, with more indirection per operation.

This is the one place in the table where a default matters, and it is worth stating plainly: on
this set, **the traversal choice is worth more than the ordering choice**.

### The ordering barely moves the fill

Geometric mean relative to the best ordering on each matrix:

| ordering | nnz(L) | ordering time | analysis | factorization |
|---|---|---|---|---|
| MMD | 1.032 | 1.566 | 1.335 | 1.081 |
| AMD | 1.020 | 1.126 | 1.045 | 1.073 |

**One to three percent of fill separates them.** For comparison, the same two orderings on
generated grid Laplacians differ by up to thirteen percent, and Oblio's own experiments have
measured a sixteen percent fill spread among four tie-breaking rules of a single algorithm. That
gap is a property of grids: nearly every live vertex has the same degree, so the tie-break decides
almost every pick. **Real structure has enough degree variation that the choice stops mattering
much.**

AMD fills slightly less here, which is the opposite of what square grids say, and by a margin small
enough that the honest conclusion is that neither wins.

**The spread lives in the analysis, not the factorization**, 1.045 against 1.335, and the next two
sections explain where it comes from. It is also the phase that disappears when a pattern is
factored more than once, since the analysis can be reused.

### A dense row costs minimum degree a factor of eighty

Two matrices ordered far more slowly than the rest under MMD, and the cause is precise.

**A single dense row or column makes minimum degree quadratic.** A vertex adjacent to everything
appears in every degree update for the whole elimination, and the published AMD algorithm addresses
this directly: rows with more than `max(16, 10 * sqrt(n))` entries are set aside before ordering,
the rest is ordered, and the dense rows are placed last. **MMD has no such step.**

Removing every column above that threshold by hand, in milliseconds:

| matrix | n | dense columns | MMD | AMD |
|---|---|---|---|---|
| bloweybq | 10001 | 1 | 70.7 -> 0.8 | 1.4 -> 1.4 |
| bundle1 | 10581 | 252 | 259.3 -> 11.0 | 7.5 -> 2.5 |

`bloweybq` has **exactly one column of degree 10000 and 9992 columns of degree 5**, and that single
vertex accounts for the entire difference. MMD falls by a factor of eighty when it is removed; AMD
barely moves, having set it aside already.

**Fill is unaffected**, 39996 against 39997 on `bloweybq`, so this is a cost in time and not in
ordering quality, and it is confined to matrices with a genuinely dense row. Two of the 107 were
slow enough under MMD to draw attention; the set was not swept for dense rows beyond that.

### And a weakness no preprocessing fixes

`Mulvey/finan512` is slower under MMD by a factor of thirty and has **no dense column at all**: its
maximum degree is 54 against a threshold of 2734, and removing nothing changes anything. It also
fills 6.5 million entries against AMD's 2.8 million, losing on both axes at once.

Its degree histogram is the explanation: 512 columns of degree 54, 512 of 51, 512 of 22, 1024 of
20, 24064 of 6. A nested block structure with massive degree ties is exactly where exact minimum
degree spends its time, and where AMD's approximate degree does not.

**So AMD's advantage on this set is not only the dense-row step.** Where degrees are heavily tied,
computing them exactly costs more than approximating them, and no amount of preprocessing changes
that.

## Where a solve's time actually goes

Median share of a one-shot solve, at multiple minimum degree and the multifrontal traversal,
over the 107
matrices, and again over the 26 with the largest factors:

| phase | all | largest quarter |
|---|---|---|
| ordering | 21.9% | 21.7% |
| rest of the analysis | 23.5% | 22.5% |
| numeric factorization | 46.7% | 47.9% |
| triangular solve | 6.6% | 7.8% |

Each median is taken over its own column, so the four are not required to sum to 100.

**The analysis is not a rounding error.** Ordering plus the elimination forest and symbolic
factorization is **45 percent of a one-shot solve at the median**, against 47 percent for the
numeric factorization. Anyone who has been told that a sparse direct solve is dominated by its
factorization is being told something that is true asymptotically and not true at these sizes.

**And the split barely moves with size**, which was not what we expected: the largest quarter of
the set looks almost exactly like the whole of it. Fill grows faster than n, so the factorization
share does rise on individual large matrices, but the analysis grows with it and the two stay
comparable across three orders of magnitude of factor size.

**The analysis is also the part that disappears on reuse.** A caller factoring the same pattern
repeatedly, with different values, pays it once. These are the one-shot shares and therefore the
upper bound on what the analysis can cost; Oblio supports reusing an analysis across
factorizations.

### Eight matrices in full

Milliseconds at multiple minimum degree and the multifrontal traversal, spanning the set from a
thousand
factor entries to 46 million:

| matrix | n | nnz(L) | order | analysis | factor | solve |
|---|---|---|---|---|---|---|
| bcsstm08 | 1074 | 1074 | 0.01 | 0.07 | 0.07 | 0.01 |
| 1138_bus | 1138 | 3301 | 0.05 | 0.04 | 0.14 | 0.01 |
| bcsstk14 | 1806 | 112507 | 0.75 | 0.38 | 0.92 | 0.12 |
| fv1 | 9604 | 300080 | 0.51 | 0.88 | 2.07 | 0.31 |
| wathen100 | 30401 | 2139398 | 2.01 | 5.48 | 12.31 | 2.24 |
| apache1 | 80800 | 12509459 | 29.66 | 31.83 | 95.66 | 14.18 |
| cfd1 | 70656 | 33765596 | 50.37 | 86.74 | 293.45 | 40.50 |
| uni_chimera_i2 | 100000 | 46605358 | 227.19 | 109.53 | 1129.30 | 55.89 |

`analysis` is the elimination forest and symbolic factorization, the ordering being its own column,
and the four phases sum to the solve. Matrix names are the collection's with the group prefix
dropped.

As shares of the total:

| matrix | ordering | analysis | factorization | solve |
|---|---|---|---|---|
| bcsstm08 | 6% | 44% | 44% | 6% |
| 1138_bus | 21% | 17% | 58% | 4% |
| bcsstk14 | 35% | 18% | 42% | 6% |
| fv1 | 14% | 23% | 55% | 8% |
| wathen100 | 9% | 25% | 56% | 10% |
| apache1 | 17% | 19% | 56% | 8% |
| cfd1 | 11% | 18% | 62% | 9% |
| uni_chimera_i2 | 15% | 7% | 74% | 4% |

**The factorization share does climb with size**, from 42 to 74 percent across these eight, and
the analysis falls correspondingly. What it does not do is take over: on the largest matrix in the
set, with a 46 million entry factor and a 1.5 second solve, the analysis is still a fifth of the
time.

**`bcsstk14` is the row worth pausing on.** Its ordering alone is 35 percent of the solve,
because multiple minimum degree is expensive on it relative to a factorization of only 112507
entries. The
two sections above explain when that happens.

## What this report does not establish

- **One factorization.** Cholesky only. LDL, static and dynamic, is not timed here; its accuracy
  is the subject of `ACCURACY.md` and its cost includes delayed columns, which Cholesky has none of.
- **Real double precision only**, and one right-hand side per matrix.
- **Positive definite matrices only**, which narrows the structural variety considerably: 20
  problem kinds here against 28 in the accuracy set. Anything said about "real matrices" is said
  about the definite ones.
- **A single machine and a single BLAS.** Traversal results in particular depend on the BLAS, since
  multifrontal's advantage comes from handing dense blocks to it.
- **No parallelism.** Everything is single-threaded.
- **Sizes to 100000 columns and 46 million factor entries.** The large end of the collection is
  untouched, and five matrices were declined for predicted fill above the cap, the largest at 1.2
  billion entries.
- **No claim about asymptotic scaling.** Measuring how cost grows with problem size needs a
  controlled family rather than a heterogeneous set, which is what `SCALING.md` does on generated
  grids.

## Reproducing this

The benchmark and the selection script sit beside this report, in `benchmarks/matrices/`. The
README there carries the working notes, including the two ordering investigations above in full.
From that directory:

```
./ssget.py list --posdef --max-nnz 2000000 > performance_candidates.txt
./ssget.py fetch performance_candidates.txt
make performance
```

The matrices are downloaded from the collection rather than committed, so the set is reproducible
without shipping somebody else's data.
