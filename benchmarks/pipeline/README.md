# Pipeline Benchmark

Where a solve's time actually goes: ordering, the rest of analyze, factorization in each of the
three traversals, and the solve. Nine orderings, Cholesky, real, on grid Laplacians.

**A benchmark, not an experiment**, on the same terms as `../ordering`: it links `../../src`
directly and is expected to keep compiling as the tree moves, which is why `make` builds it and
only `make run2d` and `make run3d` measure.

```
make          build
make run2d    build and run, square grid sides 32, 64, 100, 140
make run3d    the same on CUBIC grids, sides 6, 12, 16, 20, 26
make clean
make help     print this list

./pipeline_timing_cpp 200        any square sides
./pipeline_timing_cpp 3d 26      any cubic sides
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

**AND SUPERSEDED AGAIN, 2026-08-09, and this time the folder WAS re-run.** Ledger entry 8 fixed the
amd hash key, worth a factor of two to three on the ordering time of AMD2, AMD2B and AMD3, and the
mmd entry-5 filing defect had moved MMD2's fill on 2026-08-07 without this folder noticing. The
2026-08-01 table below is therefore stale in both branches and is kept as the record of the run
that produced it. The current numbers are in the section that follows it.

**The break-evens are the figures that moved most, and they moved by an order of magnitude.** AMD2
went from 26.6 factorizations to 2.2 and AMD2B from 70.4 to 3.2, so the sentence this folder was
best known for, that a caller who factors once should use a vendored routine and one who factors a
dozen times can use ours, now holds at two or three rather than at twenty-seven and seventy. Every
one of ours except AMD3 breaks even inside four factorizations.

### The re-run, 2026-08-09

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate.** Milliseconds, best of three, grid
140x140, n = 19600, after ledger entry 8 and after the mmd entry-5 fix. Eleven rows where the run
above has nine, MMD3, AMD3 and AMD2B having landed since.

```
ordering    order   analyze  analyzeMF     nnz(L)    factLL    factRL    factMF    solve
Natural      0.00      9.45       9.61    2744139    529.22    747.63    261.93     3.05
MMD          1.19      2.60       2.77     412921      6.88      6.83      3.42     0.43
MMD1         2.80      4.46       4.62     492921      7.37      7.48      3.94     0.51
MMD2         1.65      3.24       3.48     447712      7.25      7.41      3.71     0.47
MMD3         1.58      3.01       3.17     412921      6.86      6.87      3.44     0.47
AMD          1.32      2.77       2.91     474995      7.61      8.76      3.98     0.57
AMD1         1.89      3.40       3.57     455472      7.28      7.35      3.70     0.52
AMD1B        1.96      3.45       3.55     455472      7.11      7.21      3.77     0.47
AMD2         2.12      3.73       3.92     450190      7.18      7.10      3.80     0.48
AMD3         2.35      3.88       4.07     474995      7.44      7.55      3.89     0.51
AMD2B        2.19      3.79       3.93     450190      7.29      7.44      3.67     0.55

against       analyze    factor LL     break-even
AMD          saved ms      cost ms factorizations
MMD              0.17        -0.73         always, wins on both
MMD1            -1.69        -0.24            7.0
MMD2            -0.47        -0.36            1.3
MMD3            -0.24        -0.76            0.3
AMD1            -0.63        -0.33            1.9
AMD1B           -0.69        -0.50            1.4
AMD2            -0.96        -0.44            2.2
AMD3            -1.11        -0.18            6.3
AMD2B           -1.03        -0.32            3.2
```

**Four things to read off it, and the last is the one that changes what this folder says.**

**Every ordering of ours now breaks even inside seven factorizations and most inside two.** MMD3 at
0.3 is effectively free, and it is the default. Where the 2026-08-01 run had two rows out past
twenty-seven, the worst now is AMD3 at 6.3, and AMD3 is the one row whose cost is deliberate: it
reproduces `AMD_2`'s permutation exactly, so its fill is the vendored routine's rather than ours.

**MMD3 costs 0.24 ms of analysis against the vendored MMD and gives 0.76 ms per factorization
back.** It has the vendored routine's fill by construction and factors as fast, so against AMD it
wins on both axes at every count. That is the strongest position any of ours has ever held here.

**The analysis share is unchanged and the earlier reading stands.** Analyze is 2.6 to 3.9 ms
against a 6.9 to 7.6 ms left-looking factorization, so 27 to 36 percent of one analyze-plus-factor,
with ordering about half of analyze. Optimizing the ordering was worth the days; the fill spread
still does not propagate, 19 percent of fill against 5 of factorization time.

**And the amd branch's cost is now almost entirely AMD3's alignment rather than the extras.** AMD1
at 1.89 ms of ordering, AMD2 at 2.12 and AMD3 at 2.35: the extras cost 12 percent where before the
key fix they cost a factor. What separates AMD3 from AMD2 is the fourth pass and the seven percent
more pairs the vendored list order tests, which is the price of matching `AMD_2` exactly and is
paid for a reason.

**One caution, since this folder has exactly one problem family.** Everything above is square grids.
`benchmarks/ordering` measures cubic ones and the two disagree about our orderings, so no
break-even here should be quoted as a property of the methods rather than of this family. That is
still the first item under "What this folder still needs".

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
time it costs.

**What the number is, since the column has twice been read as a property of the methods.** It is
the count of factorizations per analysis at which our ordering overtakes the vendored one, and it
is a quotient of the two columns beside it. Both of those are differences against AMD and both
carry a sign that is easy to misread: `analyze -0.96` means our ordering takes 0.96 ms LONGER to
analyze, and `factor -0.44` means each factorization is 0.44 ms FASTER. So one analysis and `k`
factorizations cost

```
AMD     2.77 + 7.61 k
AMD2    3.73 + 7.18 k        equal at k = 0.96 / 0.44 = 2.2
```

and a caller who analyzes a pattern once and factors it three times is better off with ours.

**Both halves move independently, which is why a break-even can improve for two unrelated reasons.**
AMD2 went from 26.6 to 2.2 between 2026-08-01 and 2026-08-09, and that is about a factor of two
from the numerator, the ordering itself getting cheaper with ledger entry 8, and about a factor of
five from the denominator, AMD2's fill falling from 487111 to 450190 through the entry-4 filing
defect and the entry-8 tie-break. Neither of those is the other, and quoting the twelvefold
improvement as a speed result would be wrong.

**And the denominator is the least reliable number in this folder**, which bounds what the column
can be used for. 0.44 ms against a 7.2 ms factorization is six percent, and the section above this
one establishes that fill differences among orderings that are all roughly good do not propagate
into factorization time in any orderly way, 19 percent of fill spread giving 5 percent of factor
spread with no visible correlation. A small denominator makes the quotient swing hard: AMD3's 6.3
rests on 0.18 ms, which is 2.4 percent of a factorization and inside what this harness can resolve,
so it is the softest figure in the column and should not be quoted alone.

**It is also traversal specific**, and the tables give only the left-looking one. Under
multifrontal the same AMD2 row is 3.98 ms against 3.80, an advantage of 0.18 rather than 0.44, so
its break-even is nearer 5 than 2. Read every figure here as "if you factor left-looking".

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

**SUPERSEDED 2026-08-09, and by an order of magnitude.** The table above predates ledger entry 8
and the mmd entry-5 fix. Re-run: AMD2 breaks even at 2.2 factorizations rather than 26.6, AMD2B at
3.2 rather than 70.4, MMD2 at 1.3 rather than 5.5, and MMD3, which did not exist when this was
written, at 0.3. The worst of ours is now AMD3 at 6.3, and the sentence above should be read as two
or three factorizations rather than a dozen. The numbers are in the re-run under Results.

**And MMD is the ordering to beat, not AMD. ON SQUARE GRIDS, corrected 2026-08-09.** The fill half
of that does not hold on cubic ones: `benchmarks/ordering` now measures both families, and MMD
fills 13 percent BELOW AMD on squares and slightly ABOVE it on cubes, 2869267 against 2836813 at 26
a side. The analysis-time half is untested outside squares, this folder having only the one family.
So the claim below is a claim about square grids, which is what everything in this folder is, and
the first item under "What this folder still needs" is what would settle it. As written:

**MMD is the ordering to beat, not AMD.** It analyzes fastest of the two vendored routines,
fills the least of all nine, and factors fastest. The closest of ours is MMD2, 0.68 ms behind on
analysis and 0.12 ms ahead per factorization.

**And the closest of ours is MMD3 now, 2026-08-09, which is a different statement.** MMD2 was close
on time and behind on fill; MMD3 has the vendored routine's fill BY CONSTRUCTION, being its
permutation, so on this family it gives up 0.24 ms of analysis and nothing else. Against AMD it
wins on both axes at every count. The gap to beat is therefore 0.24 ms of ordering rather than
anything about the ordering's quality, which is a narrower and more tractable target than this
section was written against.

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

## Cubic grids, 2026-08-10, and what actually changes between the families

`make run3d`, beside `make run2d`, on the same ladders `../ordering` uses so a row here can be read
beside a row there without converting sizes. **Neither folder has a bare `make run` any more**: the
family is always named, because this section is what happens when it is not.

Every table above this section is square, and on square grids **MMD wins on all three terms**:

```
2D, MMD against the vendored AMD        order      nnz(L)     factor LL
32x32                                   0.86x      0.993x       0.938x
64x64                                   0.83x      0.941x       0.957x
100x100                                 1.03x      0.906x       0.935x
140x140                                 0.98x      0.869x       0.933x
```

It orders no slower, fills up to 13 percent less, and factors 5 to 7 percent faster, which is why
`MMD` reads "wins on both, at every count" against AMD at three of the four sizes.

**On cubic grids only ONE of those three reverses, and it is the smallest term.**

```
3D, MMD against the vendored AMD        order      nnz(L)     factor LL
6^3                                     1.33x      1.004x       1.000x
12^3                                    1.93x      0.995x       0.902x
16^3                                    1.97x      1.050x       0.925x
20^3                                    1.89x      0.998x       0.917x
26^3                                    2.24x      1.011x       0.961x
```

**AMD orders about twice as fast**, and that is the whole of its advantage. **Fill is a wash**,
within one percent at every size but 16 cubed. And **MMD still factors 4 to 10 percent faster**,
at the same nnz(L), which is the genuinely surprising row: an ordering that fills the same produces
a factor that runs measurably quicker, presumably through supernode shape rather than count.

So the break-even runs the other way from 2D rather than the ordering flipping: on cubes `MMD` pays
above about 3 to 6 factorizations, `MMD3` similarly. **For a one-shot solve AMD wins on cubes; for
anything factored repeatedly MMD still does.** In 2D MMD wins at every count.

**A correction, and it is the third time the same mistake has been made this week.** This section
first said the mmd branch fills five percent more on cubes and factors eleven percent slower, and
that every mmd method loses on both terms at every count. That was written from ONE sandbox run at
16 cubed. Sixteen is the only cubic size where MMD fills more, and the factor figure was backwards.
**A claim measured at one size of one family is a claim about that size**, which is what
`../ordering` learned on 2026-08-09 and this folder has now learned twice.

**Where `AMD3` sits.** It returns `AMD_2`'s permutation, so its fill is the vendored routine's
exactly, and against it `AMD3` reads "wins on both, at every count" at 12 and 16 cubed and pays
above 1.6 to 3.6 factorizations at 20 and 26. Against `AMD2`, which uses its own tie-break, the
fill comparison is two-sided on cubes as it is on squares: `AMD2` fills less at 12, 20 and 26 a
side and more at 6 and 16, and in 2D it fills 5 percent less at 140. So `AMD2` is generally the
better filler of the two and `AMD3` the faster orderer with a permutation that has decades of use
behind it. Neither dominates.

## What the ordering phase is worth, in proportion, 2026-08-10

Worth stating plainly beside a week of ordering work. With `AMD3`, ordering as a share of
analyze plus one left-looking factorization plus one solve:

```
140x140   order 2.04 of 11.49 ms    17.8%
16^3      order 0.94 of  5.61 ms    16.8%
26^3      order 4.96 of 47.99 ms    10.3%
```

**Between a tenth and a fifth for one solve**, falling as the factorization grows, and a smaller
share again of anything factored more than once, which is the case the break-even column exists
for. Note that `analyze` already contains the ordering, so those are shares of the whole one-shot
cost rather than of something the ordering sits beside.

So a 10 to 16 percent gain on the ordering phase is one or two percent of a single solve, and less
thereafter. That is not an argument against the work: the amd branch went from 3.0x the vendored
routine to overlapping it in three days, the alignment that made it possible is the strongest
correctness oracle in the tree, and at 16 cubed the branch now wins this table outright. It is an
argument about **where to look next**, and this table says the factorization, which is 2 to 7 times
the ordering at every size and 30 times it under a Natural ordering.

## What this folder still needs

- **Matrices that are not grids**, which the caveat above makes the first item, and which the
  ordering folder's own open questions also ask for.
- **The other factorizations.** Cholesky only, so nothing here says what LDL or dynamic pivoting
  costs, and dynamic pivoting is where the interesting work happens.
- **Complex.** Real only.
- **A profile target.** The ordering folder has one and this does not, because nothing here has yet
  needed a call tree.
- **A repeat count sized from a probe.** `bestOfThree` is a fixed three trials where the ordering
  folder sizes its count so every row is measured for about the same wall time and takes the best
  of fifteen to thirty. Three samples is thin against thermal state and page placement, which the
  performance-core request added on 2026-08-10 does not fix: the ordering benchmark still moves 7
  to 11 percent run to run on unchanged code with that call in. The callable interface
  `bestOfThree` already has makes the protocol a drop-in.
