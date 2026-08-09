# Pipeline Benchmark

Where a solve's time actually goes: ordering, the rest of analyze, factorization in each of the
three traversals, and the solve. Nine orderings, Cholesky, real, on grid Laplacians.

**A benchmark, not an experiment**, on the same terms as `../ordering`: it links `../../src`
directly and is expected to keep compiling as the tree moves, which is why `make` builds it and
only `make run` measures.

```
make          build
make run      build and run, grid sides 32, 64, 100, 140
make clean

./pipeline_timing_cpp 200      any sides
```

## Why this folder exists

`../ordering` measures one phase against itself, which answers how fast an ordering is and how much
it fills. It cannot answer whether either matters, because it never runs the phases the ordering
feeds. This folder does, and the two questions it settles are:

- **what share of a solve the ordering is**, which decides whether optimizing it is worth an
  afternoon
- **the break-even factor count**, which is the number a caller actually faces, since analyze runs
  once per pattern and factor runs per Newton step or per time step

An ordering that saves a millisecond of analysis and costs two per factorization is a loss for
anyone who factors twice. Neither the time column nor the fill column says that on its own.

## What it measures, and how

`order` is `OrderEngine::compute` timed on its own; `analyze` is the whole facade call, so the
forest and symbolic cost is the difference and no phase is double counted. `analyzeMF` is the same
call under multifrontal, which is a separate column because **analyze is not traversal
independent**: multifrontal sorts each supernode's children and relabels them into a postorder,
which left- and right-looking do not, and switching into or out of it invalidates an analysis.

nnz(L) comes from the symbolic factor, a supernode's own triangle plus its update rows, exactly as
in `../ordering`. The computation is duplicated rather than shared, deliberately: a benchmark
should stand alone.

Best of three after a warm-up, `-O3 -DNDEBUG`. Natural is included as a row because on a grid its
fill is enormous, and it is what sets the scale the other eight are read against.

**Grid Laplacians only, and that is a limitation rather than a starting point.** A grid is where
the MMD fill gap is largest, since nearly every live vertex has the same degree and the tie-break
decides almost every pick; `experiments/ordering` measured a 16 percent fill spread across four
filing orders of one algorithm. So a break-even computed here may be a property of grids rather
than of the orderings. Widening the matrices is the first thing this folder needs.

## Results

**SUPERSEDED IN PART, 2026-08-08.** Every AMD2 figure below predates a defect fix: `Amd2` and
`Amd2B` filed a supervariable one bucket too high per vertex a hash merge absorbed, because the
bound subtracts the vertex's own weight before the merge that grows it, where `AMD_2` subtracts
it after supervariable detection. Corrected, AMD2's fill is 11900 at 32 a side, 199386 at 100 and
444191 at 140, against the 12364, 212496 and 487111 recorded here, so it now beats AMD1 at every
size and the vendored AMD at the two larger ones. The figures below are kept as the record of the
run that produced them. `docs/DESIGN_DECISIONS.md` (2026-08-08) and `docs/TODO.md` carry the
finding; nothing about MMD, AMD or AMD1 moved.

**AND SUPERSEDED AGAIN, 2026-08-09.** Ledger entry 8 fixed the amd hash key, which is worth a
factor of two to three on the ordering time of AMD2, AMD2B and AMD3 and moves AMD2's and
AMD2B's fill a second time. Every `order` and `analyze` figure for those three below is
therefore stale, and so is every break-even computed from them: AMD2 and AMD2B broke even at
twenty-seven and seventy factorizations here and will now break even far sooner. This folder
has not been re-run. `benchmarks/ordering/README.md` carries the new ordering times and
`docs/DESIGN_DECISIONS.md` (2026-08-09) the finding.

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate, 2026-08-01.** Milliseconds, best of
three, grid 140x140, n = 19600.

```
ordering    order   analyze  analyzeMF     nnz(L)    factLL    factRL    factMF    solve
Natural      0.00      9.45       9.51    2744139    535.80    749.14    264.20     3.08
MMD          1.30      2.79       2.97     412921      7.27      7.18      3.43     0.45
MMD1         2.94      4.58       4.78     492921      7.35      7.68      3.93     0.53
MMD2         1.68      3.34       3.51     501951      7.46      7.61      3.95     0.53
AMD          1.25      2.65       2.83     474995      7.59      8.85      3.91     0.57
AMD1         1.91      3.43       3.64     455472      7.33      7.52      3.79     0.54
AMD1B        2.05      3.53       3.71     455472      7.25      7.34      3.69     0.49
AMD2         3.09      4.75       4.91     487111      7.51      7.72      3.84     0.52
AMD2B        3.19      4.81       4.97     487111      7.56      7.58      3.79     0.52
```

### Analysis is a large share of a one-shot solve

**Analyze is 27 to 40 percent of one analyze-plus-factor, and ordering is about half of analyze.**
At 140x140 the ordering alone runs 1.25 to 3.19 ms against a 7.3 ms factorization. Against
multifrontal's 3.8 ms factor, the ordering is comparable to the factorization itself.

That is worth stating because it was guessed at, badly, before it was measured: the working
assumption an hour before this table existed was that ordering was well under one percent of a
solve and that optimizing it was close to pointless. It is nearer a quarter.

### But the fill differences do not propagate, and that reverses a claim

```
              nnz(L)     vs AMD    factLL    vs AMD
MMD           412921      -13%      7.27      -4%
AMD1          455472       -4%      7.33      -3%
AMD           474995        --      7.59       --
MMD2          501951       +6%      7.46      -2%
```

**Fill spans 19 percent across the eight non-trivial orderings and factorization time spans 5**,
with no visible correlation: MMD fills 13 percent less than AMD and factors 4 percent faster, while
MMD2 fills 6 percent more and also factors faster.

The claim this replaces, made before the measurement, was that factorization cost grows with the
sum of squared column counts, so 22 percent more fill should mean 40 or 50 percent more work, and
that the fill gap therefore mattered far more than the ordering time. **That is true in the limit
and false in the range these orderings occupy.** Natural proves the limit: 4.8 times AMD's fill and
22 times its factorization time. Among orderings that are all roughly good, the extra fill lands in
supernodes that were already being handed to a dense BLAS, where a few more columns cost almost
nothing.

So the practical reading inverts. **Ordering time is the axis on which these methods actually
differ**, and the fill column, which `../ordering` reports first, is nearly flat in its effect.

### The break-even, which is the number a caller faces

Against the vendored AMD at 140x140, the analysis time each ordering saves and the factorization
time it costs:

```
against       analyze    factor LL     break-even
AMD          saved ms      cost ms factorizations
MMD             -0.14        -0.32            0.4
MMD1            -1.92        -0.24            8.1
MMD2            -0.68        -0.12            5.5
AMD1            -0.78        -0.26            3.0
AMD1B           -0.87        -0.33            2.6
AMD2            -2.10        -0.08           26.6
AMD2B           -2.16        -0.03           70.4
```

Every one of ours analyzes slower and factors imperceptibly faster, so each pays only above some
number of factorizations. **AMD1 and AMD1B break even at about three, MMD2 at five, AMD2 and AMD2B
at twenty-seven and seventy.** A caller who factors once per pattern should use a vendored routine;
one who factors a dozen times per analysis can use ours and not notice.

**And MMD is the ordering to beat, not AMD. ON SQUARE GRIDS, corrected 2026-08-09.** The fill half
of that does not hold on cubic ones: `benchmarks/ordering` now measures both families, and MMD
fills 13 percent BELOW AMD on squares and slightly ABOVE it on cubes, 2869267 against 2836813 at 26
a side. The analysis-time half is untested outside squares, this folder having only the one family.
So the claim below is a claim about square grids, which is what everything in this folder is, and
the first item under "What this folder still needs" is what would settle it. As written:

**MMD is the ordering to beat, not AMD.** It analyzes fastest of the two vendored routines,
fills the least of all nine, and factors fastest. The closest of ours is MMD2, 0.68 ms behind on
analysis and 0.12 ms ahead per factorization.

### Two things this turned up that belong to other parts of the tree

**Multifrontal is about twice as fast as left-looking here**, 3.4 to 3.9 ms against 7.3 to 7.6, on
every ordering and at every size. That is the largest single effect in the table.
`docs/ARCHITECTURE.md` describes the traversal trade and says explicitly that it is not measured;
on grid Laplacians with Cholesky it now is, and multifrontal wins outright. What that section says
about the trade depending on front size is untested by this, since a grid has fat fronts by
construction, and a tree-like matrix with thin fronts is exactly the case predicted to go the other
way.

**Right-looking is anomalous on the vendored AMD**, 8.85 ms against left-looking's 7.59, where the
two agree within noise on all eight other orderings. Reproduced across sizes at 100x100 as well,
4.22 against 3.83. No explanation, and it is small enough that it may be an artifact of one
permutation. Recorded rather than chased.

## What this folder still needs

- **Matrices that are not grids**, which the caveat above makes the first item, and which the
  ordering folder's own open questions also ask for.
- **The other factorizations.** Cholesky only, so nothing here says what LDL or dynamic pivoting
  costs, and dynamic pivoting is where the interesting work happens.
- **Complex.** Real only.
- **A profile target.** The ordering folder has one and this does not, because nothing here has yet
  needed a call tree.
