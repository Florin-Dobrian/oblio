# Ordering Benchmark

What each ordering method costs: wall time to produce the permutation, and nnz(L) under it. Six
methods, two lineages: the vendored MMD and AMD against Oblio's own MMD1, MMD2, AMD1 and AMD2.

**A benchmark, not an experiment.** The studies under `experiments/` are frozen: each answers one
question with a measurement and is not maintained afterwards, so nothing may depend on their
contents staying current. This is the opposite. It links `../../src` directly and is expected to
keep compiling as the tree moves, which is why `make` builds it without running it: a benchmark
that silently stopped compiling would be discovered on the day it was wanted.

```
make              build
make run2d        build and run, every method, square grids at 32, 64, 100, 140
make run3d        the same on CUBIC grids at 12, 16, 20, 26
make scale2d      one branch at a time over the ladder to 400 a side, with the gap columns
make scale3d      the same, cubic, to 32 a side
make scale-mmd-2d, scale-amd-2d, scale-mmd-3d, scale-amd-3d    one branch
make clean

./order_timing_cpp 200 280         any square sides
./order_timing_cpp amd 3d 20 26    one branch, one family, chosen sides
```

**Every target names both axes**, the question and the problem family, and the targets were `run`
and `scale` until 2026-08-09. The cubic half arrived then and left one shape named and the other
implied, which is how every figure recorded below came to be 2D without anyone deciding it should
be. A default invisible in the name is the kind that nobody revisits.

**Why `scale` is per branch.** The gap columns are ratios against the vendored routine of the same
lineage, so they only mean anything within one; and ten methods at 400 a side is a long wait for
columns nobody is reading. The AMD list is AMD, AMD1, AMD2, AMD3, with the vendored routine first
because the gap columns take the first entry as their baseline. The B variants are left out: they
are their originals' permutations on a different schedule, so their fill column carries no
information and their time column belongs to the question about the seam rather than to the
question about the branch.

## What it measures, and how

Grid Laplacians in two families. **Square** grids, since fill on one is the case every
minimum-degree paper reports and the one where the methods are known to separate. **Cubic** grids,
because that case flatters us and the other does not: the same code beats the vendored AMD on fill
in 2D and loses to it in 3D, so a conclusion drawn from squares alone is a conclusion about
squares. Ordering time only, not analysis or factorization: `OrderEngine::compute` and nothing
else. Timed as the best of however many repeats fit a fixed wall time, after a warm-up, because a
single cold reading is a reading rather than a result, and `-O3 -DNDEBUG` for the same reason.

**The two families are never rows of one table.** A row is a matrix, and a column read down across
two problem families is exactly the error the second family exists to prevent.

**`AMDraw` is a column and not a method.** `AMD` is `amd_order` as SuiteSparse ships it, which is
what a caller selecting `Ordering::AMD` gets; `AMDraw` is the same routine with an additive hook
reporting the order `AMD_2` would emit if it stopped at the end of its main loop, before
`AMD_postorder` relabels it. It is there because `AMD3` reproduces that raw order and deliberately
does not postorder, so `AMD` is the wrong thing for it to sit against: they are different
permutations and their fill need not agree. It does agree on every size in these tables, and
differs by one to three entries on cubes of 4, 5 and 6 a side, because a postorder of AMD's
ASSEMBLY tree is not fill-neutral. Its own header says that tree need not be the precise supernodal
elimination tree, mass elimination under an approximate degree merging vertices that were never
adjacent. Against `AMDraw` the comparison is exact by construction, which is why the fill gap
columns are taken against it while the time gaps stay against the shipped routine. Its own time is
not reported: the hooked copy carries the hook's bookkeeping, so timing it would measure our
instrumentation. The column needs `private/` and leaves the table with `AMD` in a published build.

nnz(L) comes from the symbolic factor rather than from any ordering's own estimate, summing each
supernode's triangle and its update rows. That is exact, it is computed from the permutation
actually emitted, and it is the same number whichever method produced the permutation, which is
what makes the comparison fair.

## Reading it

Time and fill move independently and both matter. A method that orders twice as fast and fills ten
percent more is not obviously better, and which way that trade falls depends on how many times the
matrix will be factored against one analysis.

**Record the machine, the date and the compiler beside any numbers kept here.** A timing without
them cannot be compared against anything later, and this tree already has one instance of a
measurement being read as a result when it was a reading of something else.

## AMD2, added 2026-07-31

**SUPERSEDED IN PART, 2026-08-08.** Every AMD2 figure below predates a defect fix: `Amd2` and
`Amd2B` filed a supervariable one bucket too high per vertex a hash merge absorbed, because the
bound subtracts the vertex's own weight before the merge that grows it, where `AMD_2` subtracts
it after supervariable detection. Corrected, AMD2's fill is 11900 at 32 a side, 199386 at 100 and
444191 at 140, against the 12364, 212496 and 487111 recorded here, so it now beats AMD1 at every
size and the vendored AMD at the two larger ones. The figures below are kept as the record of the
run that produced them. `docs/DESIGN_DECISIONS.md` (2026-08-08) and `docs/TODO.md` carry the
finding; nothing about MMD, AMD or AMD1 moved.

AMD1 plus aggressive absorption and hash supervariable detection. Measured on the Linux sandbox at
140x140, as ratios to the vendored routines, the machine being too loaded for the milliseconds to
mean anything: **AMD2 about 2.9x AMD, against AMD1's 1.5x**, with fill of 487111 against AMD1's
455472 and the vendored AMD's 474995.

So on grids it is worse on both counts than the branch it extends: 5 percent more fill, which is
the coarser supervariables the hash produces, and roughly twice the time, the hash pass being new
work and its live merges forcing the reachable set to start checking for dead vertices. It does
refresh fewer degrees, which is what it was meant to buy, and on a grid that does not pay for
itself.

That is worth stating plainly rather than burying: **AMD2 is currently the worst of our three on
this test set**, and the test set is one problem family. Whether hash detection earns itself on a
matrix with real supernodal structure is the question the grids cannot answer, and it is now the
strongest argument for widening the benchmark's matrices.

**No longer true on fill, 2026-08-08.** AMD2 is now the BEST of ours on every grid size measured,
11900 at 32 a side against AMD1's 12074, 199386 against 201856 at 100 and 444191 against 455472 at
140, and it beats the vendored AMD at the two larger sizes. The paragraph above was measuring a
defect in the filing rather than the cost of the mechanism: the bound subtracted a vertex's own
weight before the hash merge that grows it, so every supervariable the hash produced was filed one
bucket too high. On time nothing has changed and AMD2 remains the slowest of ours, so the argument
for widening the matrices stands on that half and on the 3D result in `REPORT.md`.

## Where the gap actually is, 2026-07-31

The question this folder was made for: ours order with comparable fill to the vendored routines
while doing less work, and are slower. The method is in `../README.md`; this is what it found.

**Cycles, split into work and efficiency.** 3000 orderings of a 140x140 grid, alpamayo, Instruments
CPU Counters:

```
                cycles           useful    useful cycles   cycles vs vendored   work vs vendored
MMD1            43.9 G           42.19%    18.5 G          3.05x                2.99x
MMD             14.4 G           43.02%     6.2 G
AMD1            29.5 G           49.02%    14.4 G          1.95x                1.68x
AMD             15.1 G           56.85%     8.6 G
```

**MMD1's gap is entirely work.** Its efficiency is the vendored routine's, 42.2 percent against
43.0, so every cycle of the 3x is spent doing three times as much. There is nothing structural to
win: no layout, no allocator, no locality. The mechanisms `genmmd` has and `mmd1` does not are the
whole story, and MMD2 is the answer to it.

**AMD1's gap is 86 percent work and 14 percent stalling**, 1.68x on useful cycles with the rest
showing as instruction processing at 20.5 percent against 15.0.

**Element counts agree, independently.** Instrumented copies of both, 100x100 grid:

```
                          ours (AMD1)        vendored AMD
build the clique               73891      (their Lme pass, not counted)
prune A and I                 179630      fused into their scan 2
degree scan                   304047                  216662
adjacency scan                 29991                   50230
TOTAL                         557568                  266892
pivots                          7512                    7635
```

Their scan 2 prunes the vertex, accumulates the degree and builds the hash in **one visit per
element**. Ours visits three times: once to prune, once for scan 1, once for the bound. That is
where AMD1's extra work is, and it is a restructuring of the driver-eliminator seam rather than
anything to do with memory.

**And the mechanisms are not the difference.** Absorption and hash detection fire identically on a
grid, 1 clique absorbed and 2488 merges, for the vendored routine and for our AMD2 alike. So AMD1's
lack of them costs it almost nothing here, and AMD2's having them buys it almost nothing, which is
why AMD2 measures worse than AMD1 overall.

**Careful with the fusion reading**, because it looks like a contradiction and is not. Fusing those
three visits into one was tried and bought nothing (see the negative results in `../README.md`).
Fusing removes loop setup and a re-fetch of data still in L1; it does not remove the per-element
work, which is the test, the copy and the weight lookup. Reducing visits means doing the work of
three visits once, not running three loops as one.

## MMD2, added 2026-08-01

MMD1 plus the prepass, the filing convention, the element-by-element q2h refresh, pairwise merging
and outmatched marking. The profile had said MMD1's whole gap was work rather than efficiency, so
this is the change that was predicted to pay, and it is the only prediction all day that did.

Linux sandbox, 140x140: **MMD2 6.57 ms against MMD1's 9.36**, so 30 percent off, taking the MMD
branch from about 2.6x the vendored MMD to 1.8x. Fill moves either way with the tie-breaks, 217102
against MMD1's 223806 at 100x100 and 501951 against 492921 at 140x140.

That leaves the MMD branch's remaining gap at roughly 1.8x, which is close to where AMD1 sits, and
the two now have the same character: neither is missing a mechanism, both visit elements more times
than the vendored scan does. The driver restructuring is what addresses that, for both.

## The incidence lists into the adjacency run, 2026-08-01

The measured target from the previous session's profile: `I[u]`'s `push_back` in `eliminate` was
1.12 s of `amd1`'s 7.08 s, the largest single line in the program. The expectation was an arena.

**It needed no arena, because `I[u]` never had to grow.** `A[u]` and `I[u]` are the two kinds of
source `reach(u)` is a union over, and an elimination destroys one for each it creates: where `u`
was reached through `A[pivot]`, the prune drops `pivot` from `A[u]`; where it was reached through a
clique, that clique is absorbed and leaves `I[u]`. So `|A[u]| + |I[u]|` never exceeds `u`'s
original degree, the two lists share one block sized once from the pattern, and the incidence is
written into the room the adjacency has just given up. Both vendored routines do this and neither
says why. Section 5.3 of `archive/sparse_factorization.md` carries the argument, 5.15 the reading
of the two codes.

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate, 2026-08-01.** Ordering time in
milliseconds, best of three, 140x140:

```
                 before     after    vs vendored, before -> after
MMD                1.18      1.18
MMD1               3.46      3.12         2.98x -> 2.64x
MMD2               2.28      1.86         1.97x -> 1.58x
AMD                1.27      1.26
AMD1               2.42      2.12         1.90x -> 1.68x
AMD2               3.81      3.58         3.00x -> 2.84x
```

The two vendored columns move under three percent, which is the noise floor and is what says the
rest is real. Every fill figure is unchanged, digit for digit, on both machines.

**MMD2 at 1.58x is now the closest any of ours has come**, and the MMD branch has passed the AMD
branch for the first time.

**Allocations, which moved much further than the time did.** At 140x140, counted by replacing
global `new`:

```
              before     after
MMD1           31915      2457
MMD2           29499       105
AMD1           32256      2760
AMD2           46351     16855
```

The quotient graph itself now allocates 11 times to build and under 200 over a whole ordering, so
what remains belongs to the drivers and to two specific places: the `merged` list `eliminate`
returns by value, about 1700 of MMD1's and AMD1's remainder, and AMD2's `hashGroup`, which is the
whole of its 16855. Neither is `I[u]` and neither is the clique arena, which is one buffer that
doubles and costs a couple of dozen allocations for the entire run.

**And the gap between those two tables is the result, again.** Removing 92 percent of the
allocations bought 12 percent on AMD1, which is the third time in this folder that allocation
counts have failed to predict time. What was bought here was a line in the innermost loop, not the
mallocs, and the honest reading is that even that came in under the profile's attribution: the
trace credited `push_back` with 15.9 percent alone, before `insert`, `resize` and the 19600 vector
headers in the constructor that went with it. Where the difference went is the next thing to
measure rather than to argue about.

**What this does not touch is AMD's remaining gap**, and it was never going to. The cycle profile
had already said AMD1 is 86 percent work and 14 percent stalling, and the element counts say the
work is the three visits per element against the vendored scan's one. That is the driver
restructuring, and it is unchanged by this.

## Hoisting the loop bounds, 2026-08-01

The Time Profiler trace taken after the change above put `QuotientGraph::incidenceSize(int) const`
fourth in the run, 309 ms weight and 300 ms self, which is absurd for an inlined accessor that
returns one `std::size_t`. The cause was the call sites, written as

```cpp
for (std::size_t i = 0; i < qg.incidenceSize(u); ++i)
```

The loop bodies store through `mark`, `outside` and `touchedCliques`, and nothing tells the
compiler those cannot alias the size vector, so the bound is re-loaded once per element. Hoisting
it to a local before the loop took that entry to 80 ms weight and 76 ms self, and AMD1 from 6.31 s
to 6.09 s.

**But only where the loop is long, which took two measurements to learn.** Hoisting everywhere made
AMD2 2.3 percent worse in a cache simulation and MMD2 3.2 percent worse on alpamayo. Both
regressions were in loops that do not run: AMD2's hash comparison exits on the first mismatch
through `&& same`, and MMD2's q2h walk is entered only by vertices with
`adjacencySize + incidenceSize == 2`, so those loops cover at most two elements between them. A
length loaded before a loop that runs once is overhead. Reverting those eight sites recovered both.

So the rule, and it is in the comments at the sites that decline it: **hoist where the loop is long,
leave it where the loop is short or exits early.**

Worth knowing that `g++` at `-O3` on x86-64 does this hoisting itself, so a Linux measurement of
this change reports nothing at all. Apple Clang on arm64 evidently does not. A negative result on
one toolchain says nothing about the other, which is the same lesson as the gprof entry in
`../README.md` in a different form.

## The boolean flags off `std::vector<bool>`, 2026-08-01

`std::vector<bool>` is the standard's bit-packed specialization: one bit per element, a proxy
rather than a reference, and no `data()`. `QuotientGraph::mEliminated`, `Buckets::mFiled` and
`Mmd2`'s `outmatched` were all built on it, inherited from the prototypes, where it is the right
choice. All three are now `std::vector<std::uint8_t>`, a plain byte each.

**alpamayo, 140x140, best of three, with the two vendored routines as controls:**

```
                 before    after    delta
MMD (control)      1.18     1.15
AMD (control)      1.26     1.22
MMD1               3.05     2.80    -8.2%
MMD2               1.87     1.65   -11.8%
AMD1               2.04     1.81   -11.3%
AMD2               3.50     3.23    -7.7%
```

The controls drifted about 3 percent downward in the same run, so read the branch gains as roughly
5 to 9 percent. Every one of them clears the drift and three clear it threefold.

**And the mechanism is not the one to guess.** The obvious story is the access cost, a load, a
variable shift and a mask against one byte load. That is not what the traces show. Two Time
Profiler runs of MMD2 agreeing with each other:

```
                      before    after    after
QuotientGraph (ctor)   710 ms   188 ms   184 ms
eliminated(int)        208 ms   189 ms   225 ms
Buckets::filed(int)     94 ms    92 ms    89 ms
total                 5.44 s    4.78 s   4.92 s
```

**The graph constructor fell from 710 ms to about 186 ms and the flag reads moved by noise.** A
`std::vector<bool>` of n false is built through the bit-reference machinery word by word; a
`std::vector<std::uint8_t>` of n zeros is a `memset`. `Buckets` and the graph are constructed once
per ordering, 3000 times in this harness, and about 170 microseconds apiece is the difference
between a loop and a `memset` at n = 19600.

That also explains the branch split, which the access story gets wrong: every branch constructs
both arrays once per ordering whatever its access pattern, so the saving is nearly uniform, which
is what the table shows. **A per-construction cost, not a per-read one.**

Two honest caveats. `order_profile` orders 3000 times, which weights construction more heavily than
a workload that orders one matrix once, so the saving is real per ordering but this harness makes
it the headline. And `outmatched`, the third array, measured at nothing: MMD2 moved 1.2 percent
when it changed, the same as two branches it cannot touch. It is kept for consistency, not on
evidence.

**Where this leaves the four**, alpamayo, 140x140, ratios against the same run's own controls:

```
              this morning     now
MMD1                 2.98x    2.46x
MMD2                 1.97x    1.45x
AMD1                 1.90x    1.48x
AMD2                 3.00x    2.57x
```

Both branches we would ship are inside 1.5x. Every fill figure in this folder is unchanged, digit
for digit, through all three changes.

**And the instrument has run out.** Four of the last five decisions landed in the 1 to 4 percent
band against 1.7 to 4 percent of run-to-run drift, which `make run`'s best-of-three cannot resolve.
Below about 5 percent, the Time Profiler line is the instrument, and it is what actually settled
both changes above.

## Where the cycles went, 2026-08-01

CPU Counters on AMD1, 3000 orderings of a 140x140 grid, taken at three points in the day so the two
kinds of gain can be told apart:

```
                     cycles    useful    useful cycles   what changed
2026-07-31          29.50 G    49.02%        14.46 G     (baseline)
after the shared run 25.93 G   43.27%        11.22 G     work removed
after both passes   22.49 G    47.05%        10.58 G     efficiency recovered
AMD (vendored)      15.10 G    56.85%         8.58 G
```

**The two gains are of opposite kinds, and the middle row is why the split is worth taking.** The
shared run cut useful cycles by 22 percent while efficiency fell nearly six points: it deleted
allocator bookkeeping, which is high-IPC work, so what remained was a harder mixture. The hoist and
the byte arrays then cut cycles another 13 percent while useful cycles fell only 6, putting 3.8
points of efficiency back. Removing work and removing stalls look identical in a timing table and
opposite here.

**What it means for the driver restructuring.** Against the vendored routine AMD1 went from 1.95x
cycles and 1.68x useful cycles to **1.49x and 1.23x**. So the remaining gap is now roughly half work
and half stalling, where on 2026-07-31 it was 86 percent work. The restructuring, which removes two
of three visits per element, attacks the work half alone: its ceiling is now about 20 percent of
AMD1 rather than the 40 it looked like this morning. Still the largest single item, no longer the
whole answer, and it should not be budgeted as one.

## The last of the fine-grained allocation, 2026-08-01

Two sites, both inherited from the prototypes, and after them the ordering allocates fewer than a
hundred times per run.

**AMD2's hash buckets** were a `std::vector<std::vector<std::int32_t>>` over n + 1, constructed and
destroyed once per ordering plus one allocation per bucket a step used. They are now `hashHead` over
n + 1 and `hashNext` over n, which is the idiom `Buckets` already uses and which `Amd.cpp` uses for
the same job. **`QuotientGraph::eliminate` returned its `merged` list by value**, one allocation per
elimination that merged anything; it is a member scratch returned by const reference now, the shape
`mReached` already had.

```
allocations at 140x140     this morning     now
MMD1                              31915      70
MMD2                              29499     105
AMD1                              32256      62
AMD2                              46351      56
```

**alpamayo, 140x140, best of three:**

```
              before    after    delta      controls drifted -2%
MMD1            2.83     2.71    -4.2%
MMD2            1.67     1.61    -3.6%
AMD1            1.88     1.86    -1.1%
AMD2            3.27     3.06    -6.4%
```

So about 4 percent for AMD2 and 1.5 for the two MMDs once the drift is taken out, and nothing for
AMD1. That is the expected shape: AMD2 held nearly all the remaining allocations, and AMD1 has now
gained nothing from two consecutive allocation changes, having been the least allocation-bound of
the four throughout.

**One hazard worth carrying forward.** A chain pushed at the head comes out reversed, and the order
within a hash bucket decides which of two indistinguishable vertices absorbs the other, so the first
version of this change moved the permutation on four graphs. Filling the chain in reverse restores
it exactly. **A container swapped for a chain can be a tie-break change wearing a data-structure
change's clothes**, which is the same trap the degree buckets' filing convention represents, and the
only thing that caught it was diffing 76 permutations against the previous tree.

**And the sandbox over-predicted twice today.** Instruction counts under `g++` on x86 said this pass
was worth 8.5 percent to MMD2; it measured about 1.5 on alpamayo. Earlier the same instrument said
the loop-bound hoist was worth nothing, and it was worth 3.5 percent. So a cache simulation on
another machine is useful for deciding **whether a change is worth trying here**, and worth nothing
as a prediction of what it will buy.

**Where the four stand after the day**, ratios against each run's own controls:

```
        this morning     now
MMD1           2.98x    2.40x
MMD2           1.97x    1.42x
AMD1           1.90x    1.50x
AMD2           3.00x    2.47x
```

Fill unchanged, digit for digit, through all five changes.

**Fine-grained allocation is now closed as a subject here.** What remains is coarse: a handful of
whole arrays per ordering, which is what the elimination forest and the symbolic factorization have
always had. Note that closing it does not close *memory*: the `std::vector<bool>` result above was
about initialization rather than allocation, and the graph constructor is still 3 to 4 percent of
MMD2 and mostly first touch of fresh pages, which no reduction in allocation count reaches.

## AMD1B, and the diagnosis it falsified, 2026-08-01

The largest item on the ordering list since 2026-07-31 was this: the vendored `amd_2` prunes a
vertex, accumulates its degree and builds its hash in one visit per element, where AMD1 visits each
element of `A[u]` twice and each element of `I[u]` three times. 557568 element visits against
266892 on a 100x100 grid. It was called the driver restructuring, and it was expected to be worth a
fifth of AMD1.

**AMD1B is that change, built.** `QuotientGraph::eliminate` split into a private
`beginElimination` and `finishElimination` with the prune loop between them, and a second public
overload taking an `ApproximateScan` that folds the driver's first scan into the prune. `A[u]`
falls to one visit and `I[u]` to two, which is what `amd_2` costs; the bound cannot be folded
further, because `outside[c]` is complete only after every member of `C[p]` has been seen, which is
why the vendored routine has two scans as well.

**It is not faster.** alpamayo, 140x140: AMD1B 1.93 ms against AMD1's 1.84, five percent slower,
and one to two percent faster at the three smaller sizes. Fill identical at every size, permutation
identical on 136 graphs.

**So the three-visit diagnosis was wrong**, and it had been cited three times before being tested.
The visits were never the cost: a second pass over a list still in L1 is a load and a compare, and
the per-element work the merge saves is a test and a copy that were cheap to begin with. This is
the same lesson as the loop-fusion negative already in `../README.md`, and the distinction drawn
there, that fusion removes loop setup while merging removes work, turns out not to have been the
thing that mattered.

## Where AMD1's gap actually is, and two hypotheses that failed, 2026-08-01

The first comparative profile in this folder: the vendored AMD and AMD1 through the same binary,
3000 orderings of a 140x140 grid on alpamayo.

```
                    AMD (vendored)        AMD1
cycles                  15.03 G          22.49 G       1.50x
useful                   8.58 G  57.1%   10.58 G  47.1%
instruction processing   2.27 G  15.1%    5.71 G  25.4%   waiting on data
instruction delivery     3.10 G  20.6%    3.85 G  17.1%   waiting on instructions
discarded                1.04 G   6.9%    2.27 G  10.1%   branch mispredicts
```

**Both stall heavily, and on opposite things.** The vendored routine wastes 43 percent of its
cycles and ours 53, so this is a memory-bound problem for everyone. Their dominant waste is the
front end, which is what one 1800-line function does to an instruction cache; ours is the back end,
waiting on data.

Splitting the 7.46 G gap by where it goes:

```
data stalls        +3.44 G    46%
useful work        +2.00 G    27%
branch mispredicts +1.23 G    17%
front end          +0.75 G    10%
```

**Work is the third-largest component, not the first.** AMD1B attacked it and gained nothing, which
is consistent.

**The second hypothesis, and it failed too.** Data stalls being the largest component, the obvious
lever was footprint: `QuotientGraph` holds six `std::size_t` arrays read in the innermost loops,
940 KB at n = 19600 where `int32_t` would need 470. Narrowed as an experiment, cachegrind reported
**D1 misses down 17 percent and last-level misses down 13, with the instruction count flat.** On
alpamayo it measured nothing at all, and if anything slightly slower, with both vendored controls
sitting still. Reverted; `std::size_t` for a position stands, and the convention now has a
measurement behind it rather than only a rule.

**Which retires cache simulation for this question.** Cachegrind was wrong three times in one
afternoon and in three different directions: it reported zero for the loop-bound hoist, which was
worth 3.5 percent; 8.5 percent for an allocation change worth about 1.5; and a large improvement in
the exact quantity being targeted here, worth zero. Instruction counts transferred faithfully all
day. **Cache behavior did not transfer even comparatively**, which is what Apple Silicon's much
larger caches and more aggressive prefetching would predict, and is worth knowing before the next
person reaches for it.

**So AMD1's 1.46x is unexplained.** Two hypotheses with good counter evidence have been built and
falsified. What is left is unexamined rather than unlikely, and none of it is more than a guess:
the branch mispredicts, which more than doubled and are 17 percent of the gap; `Buckets`, which has
never been profiled; and the possibility that a dozen separate arrays cost through the prefetcher
and the TLB rather than through cache misses, which is exactly what cachegrind cannot model and
what CPU Counters reports only in aggregate. **The next step is a new instrument, not a fourth
idea.**

One correction to the headline while we are here. The vendored run's 3.67 s includes `AMD_aat` at
139 ms and `AMD_postorder` at 154 ms, neither of which we do: we take the pattern directly and
`ElmForestEngine` does its own postorder. The comparable part is about 3.3 s, so **the true gap is
nearer 1.6x than 1.46x.**

## AMD2B, and a pattern both B variants share, 2026-08-01

The same fusion applied to AMD2, added for the oracle rather than for speed: AMD2 carries an
absorption pass and a hash pass that AMD1 does not, so it has more places to go quietly wrong, and
`AMD2B == AMD2` guards all of them at once. 136 graphs, zero mismatches, fill identical at every
size.

**alpamayo, milliseconds, best of three:**

```
grid        AMD1    AMD1B          AMD2    AMD2B
 32x32      0.20     0.16  -20%    0.35     0.28  -20%
 64x64      0.56     0.54   -4%    0.93     0.85   -9%
100x100     1.01     0.97   -4%    1.61     1.56   -3%
140x140     1.86     1.94   +4%    3.06     3.10   +1%
```

**Both pairs show the same shape: the fusion helps at small n and stops helping, or hurts, at
large n.** Two independent instances of one pattern is worth more than either alone, and it points
at a mechanism rather than at noise.

**The likely mechanism, and it is uncomfortable.** The fusion is not free of memory: it adds
`explicitPart`, an array of size n that the prune writes and the bound reads, 156 KB at n = 19600.
At small n everything is resident and the saved visits are pure gain; at large n the new array is
another stream competing with the dozen already there. That is consistent with the comparative
profile above, which found AMD1's largest gap component to be data stalls, and it means **the
fusion trades element visits for footprint, which is the wrong direction for the branch's actual
constraint.**

It is a hypothesis, not a finding. What would test it is removing `explicitPart` and accepting one
more walk of `A[u]`, which is a third variant nobody has asked for. Recorded rather than pursued.

**Both B variants are kept**, on the oracle and on the seam being reusable, not on speed. Their
collapse condition, written before either was built, was that a B variant replaces its original
when permutation-identical and faster; neither is faster, so neither fires.

## AMD3, added 2026-08-08, and the one line that was 75 percent of it

AMD2 with the vendored routine's list order, adding no mechanism, so its fill is the vendored AMD's
exactly at every size and the only question it poses is cost.

**alpamayo, `make scale-amd`, ratios against the vendored AMD in the same row:**

```
grid          AMD ms    AMD1     AMD2     AMD3        AMD nnzL   AMD3 nnzL
 32x32          0.07    1.16x    2.06x    2.18x          11900       11900
100x100         0.66    1.52x    2.46x    2.56x         206332      206332
140x140         1.27    1.51x    2.42x    2.55x         474995      474995
400x400        10.60    1.80x    2.58x    2.81x        5663298     5663298
```

**Before entry 6 those AMD3 ratios were 3.54x, 4.41x, 4.07x and 4.63x**, and the whole difference
was one line. The alignment's fifth entry put the new clique at the front of every `I[u]`, which is
what `Amd.cpp` does, and did not also skip that entry in the hash pass's exact test, which
`Amd.cpp` does in the same breath. So a guaranteed match sat at the head of a short-circuiting walk
and every failing pair paid one extra iteration.

**Four hypotheses were built or argued before the right one, and three were wrong**, which is worth
recording because each was plausible and each had a mechanism:

- a `std::rotate` in the prune, `O(|A[u]|)` per reached vertex. Real extra work, folded into the
  compaction, and worth nothing: the sandbox reading that suggested it was noise.
- `partial`, a size-`n` array the entry-4 split carries, argued twice from the `explicitPart`
  precedent in the AMD1B entry above, which it resembles exactly. The trace never implicated it.
- the hash pair count, measured and equal within 7 percent.
- and the one that was right, which no count could see: **the pair loops short-circuit, so what
  they cost depends on WHERE the mismatch is.** Every counter measured what the test COULD cost.
  Counting iterations actually executed gave 2.08 per pair against AMD2's 1.21, exactly one extra,
  which is what a guaranteed match at position zero predicts. 1.08 after.

```
                        AMD2      AMD3 before    AMD3 after
Time Profiler total    8.90 s        14.90 s        9.44 s
orderAmdN SELF         4.76 s        10.46 s        4.76 s
```

Self weight identical to AMD2's to the millisecond. The residual 6 to 9 percent is a fourth pass
AMD2 does not have plus 7 percent more pairs tested, which is about its size.

**The method note, since this folder exists for it.** Instruments' call tree could not find this:
everything is inlined into one symbol, so the tree said only that 5.70 s sat in `orderAmd3`'s self
weight. The SOURCE VIEW named the line, at 6.22 s of a 14.90 s run. When a driver is one inlined
symbol, the call tree bounds the search and only the source annotation ends it. And counting, which
this folder recommends before profiling, was actively misleading here: it is the right instrument
for "are we doing more work" and it cannot see a loop whose cost is decided by an early exit.

## Cubic grids, 2026-08-09, and the 2D advantage does not survive

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate.** `make run3d`, milliseconds, best of a
timed repeat count; nnz(L) in entries. The first cubic figures this folder has ever carried, and
the reason they exist is that every other number here is square.

The default cubic ladder gained a 6 cube after this run, so a rerun has one more row than the
tables below. It is there because 6 a side is the largest size at which `AMD` and `AMDraw`
disagree, 3265 against 3266, and a pair of columns that always matches is a pair nobody reads. The
ratio table in the section after this one comes from `make scale2d` and `make scale3d`, which run a
wider ladder than these.

```
grid3d      n      MMD    MMD1    MMD2    MMD3     AMD    AMD1    AMD2    AMD3
12^3     1728     0.47    1.43    0.46    0.46    0.24    0.28    0.91    0.74
16^3     4096     1.66    4.77    1.70    1.62    0.76    0.95    2.76    2.29
20^3     8000     3.66   11.02    3.65    3.48    1.76    2.24    6.09    5.29
26^3    17576     8.89   34.60    8.18    7.78    4.07    5.76   14.88   12.30

grid3d      n      MMD nnzL   MMD1 nnzL   MMD2 nnzL   MMD3 nnzL
12^3     1728         75674       81415       77504       75674
16^3     4096        295113      316502      303848      295113
20^3     8000        840560      860800      808926      840560
26^3    17576       2869267     3086339     2943504     2869267

grid3d      n      AMD nnzL  AMDraw nnzL   AMD1 nnzL   AMD2 nnzL   AMD3 nnzL
12^3     1728         76038        76038       77240       71848       76038
16^3     4096        281014       281014      310541      287924      281014
20^3     8000        842282       842282      877430      819997      842282
26^3    17576       2836813      2836813     3225460     2907586     2836813
```

**The claim that our tie-break beats AMD's was a 2D result, and it does not survive.** `AMD2`
against the vendored AMD reads `-5.5, +2.5, -2.6, +2.5` percent across the cubes, and `-4.6` at
32 a side from `make scale-amd-3d`. Two-sided, and it does not settle with size. In 2D the same
comparison is a clean monotone `-1.7` to `-6.5` percent, improving with n. So the honest statement
is that our tie-break is **different, not better**, and the square grids were one problem family
behaving consistently rather than a property of the ordering. `MMD2` is the same story with the
same signs, `+2.4, +3.0, -3.8, +2.6, -8.0`.

That closes the question `experiments/ordering/REPORT.md` parks under "is LIFO actually better, or
is genmmd merely good?", at least on the amd side: the 6.5 percent was the whole of the evidence
for it, and it is a 2D artifact.

## What the two families say about where the amd gap is, 2026-08-09

The first cubic run answered the question it was built for, above, and then said something the
square grids could not. It is recorded here rather than in `REPORT.md` because the evidence is a
benchmark table, and it is the strongest hint this folder has produced about the amd constant
factor.

**Split each branch into its base and its extras**, ratios against the vendored routine of the same
lineage, from `make scale2d` and `make scale3d` on alpamayo:

```
                2D, 32 to 400 a side       3D, 12 to 32 a side
MMD1            1.47x  ->  2.85x           3.06x  ->  4.91x        base,   WORSE in 3D
MMD3            1.17x  ->  1.51x           0.96x  ->  0.90x        extras, BETTER, and under 1

AMD1            1.19x  ->  1.76x           1.22x  ->  1.55x        base,   FLAT
AMD3            1.96x  ->  2.54x           3.42x  ->  3.23x        extras, WORSE
```

**`AMD1` costs the same on both families**, about 1.2 to 1.8x, and if anything slightly less in 3D
at comparable n. So the whole of the amd branch's degradation on cubic grids is in the EXTRAS,
aggressive absorption and hash supervariable detection, which take the branch from 1.5x to 3.0x in
3D against 1.5x to 2.4x in 2D.

**That is a hint about where the remaining gap is, and it points away from where we have been
looking.** The parked proposal in `docs/DESIGN_DECISIONS.md`, giving cliques their own mark space
so that liveness folds into the vertex marks, is a change to the SHARED QUOTIENT GRAPH. It would
help `AMD1` exactly as much as `AMD3`, and `AMD1` is the half that is already fine on both
families. The same is true of the locality hypothesis, one `Iw` pool against our two structures:
whatever it costs, it is being paid by a driver whose ratio does not move between families. Neither
is refuted by this, and both would still be worth what they are worth in 2D. But a fixed cost in
the shared class cannot explain a gap that appears only when the extras are switched on and only on
one family.

**And the two branches move in mirror image**, which is `experiments/ordering/REPORT.md`'s finding
2 in a new form: on mmd the base collapses in 3D and the extras rescue it, on amd the base holds
and the extras spoil it.

**One mechanism accounts for both directions, and it is CLIQUE SIZE.** A cubic grid makes far
larger cliques than a square one at comparable n, so any cost proportional to clique MEMBERS grows
with the family while any cost proportional to clique COUNT does not.

- `MMD1` recomputes a reach per refreshed vertex, walking the members of every clique that vertex
  names, so it pays members per vertex and reaches 4.9x.
- `MMD2` and `MMD3` route through the q2h path, which computes per clique and shares the result
  across its members. That is why they survive the family change, and apparently why they overtake:
  genmmd has the same idea and we are executing it faster once the loop dominates.
- `AMD1`'s bound reads one number per clique and never opens it, which is the whole content of
  section 5.13. It cannot care how large a clique is, and the measurement says it does not.

**Which leaves the amd extras, and the hash pass is the suspect.** Its pair loop costs the sum of
SQUARED bucket sizes over `C[p]`, and `C[p]` is much larger in 3D, so the term that is quadratic in
bucket size is exactly the one the family change inflates. Aggressive absorption is one test per
touched clique, already computed for the bound, and cannot be it. The vendored routine carries both
mechanisms and does not degrade, so this is our implementation rather than the idea.

**The measurement that would settle it is a COUNT and not a profile**, which is worth saying since
this folder's own advice is to profile: the question is whether we are doing more work, not whether
we are doing it less efficiently, and a count answers that machine-independently. Pairs tested per
iteration, and the bucket size distribution, 2D against 3D at comparable n. The prototypes already
count hash merges, so pairs tested is one counter beside an existing one. If the pair count per
pivot is flat across families the hypothesis dies and the extras' cost is elsewhere; if it grows
with clique size, the fix is a better filter rather than a faster loop.

**Three further readings, none of them expected.**

- **`MMD3` is FASTER than genmmd in 3D**, 7.78 ms against 8.89 at 26 a side, 0.88x, where it runs
  about 1.2x in 2D. The first time one of ours has beaten a vendored routine on time here, and it
  corroborates the 0.86x `REPORT.md` measured for `MMD2` at 32 a side on a different machine.
- **`AMD3`'s time gap WIDENS in 3D**, 12.30 against 4.07, so about 3.0x where 2D gives 2.3x. The
  parked constant-factor work is worth more on this family than the square figure suggested.
- **"MMD is the ordering to beat" is 2D-only too.** MMD fills 13 percent below AMD on squares and
  slightly ABOVE it on cubes, 2869267 against 2836813 at 26 a side. `benchmarks/pipeline/README.md`
  states that conclusion without qualification and should not.

**`MMD3` and `AMD3` are 0.0 percent at every size**, which is the alignment holding on a family it
was never measured on, and `AMD nnzL == AMDraw nnzL` throughout, which is the postorder being
fill-neutral at every size a table here reports.

**And the fill columns reproduce across machines exactly**, digit for digit against a Linux run of
the same binary. That has always been true of this benchmark's fill and is worth restating for the
new family: fill is a deterministic function of the pattern, so it is a usable regression signal
where timings never are.

## The hash key, and what the section above was actually seeing, 2026-08-09

The tables above put the amd branch's whole cubic degradation in the extras and named the hash
pass's pair loop as the suspect, on the grounds that its cost is the sum of squared bucket sizes
over `C[p]` and `C[p]` grows with the family. They also said the measurement wanted was a COUNT
rather than a profile. That count was taken, and it found a defect: our hash key multiplied its
incidence half by a stride of `size + 1` and then reduced modulo the same number, which annihilates
that half exactly, so the bucket was a function of the adjacency alone. Against the vendored
routine on the same graphs, for the SAME MERGES, we tested 19.0 pairs per pivot at 140 a side
against its 0.333, and 155.3 at 26 cubed against its 0.484.

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate, 2026-08-09.** `make scale3d` and
`make scale2d`, against the rows recorded earlier the same day:

```
grid3d          AMD      AMD1      AMD2      AMD3          AMD3 vs AMD
26^3  before   4.07      5.76     14.88     12.30                3.02x
26^3  after    4.06      5.69      5.45      5.83                1.44x
32^3  after    8.81     13.03     12.33     12.71                1.44x

grid            AMD      AMD1      AMD2      AMD3          AMD3 vs AMD
140x140 before  1.27      1.51x     2.42x     2.55x
140x140 after   1.21      1.80      2.04      2.24                1.81x
400x400 after  10.49     17.47     19.86     21.80                2.08x
```

**The two controls do not move**, the vendored routine and `AMD1` sitting within a percent, which
is the drift and is what says the rest is real. `AMD3` on cubes goes from about 3.0x to 1.44x and
`AMD2` from 2.85x to 1.40x. In 2D the same change is worth about 1.4x, which is the ratio the count
predicted from 155 pairs per pivot against 19: the defect was family dependent because `A[u]`
empties as the elimination proceeds and cubes reach that state sooner.

**Three readings, and the second and third are where the next work is.**

- **`REPORT.md` finding 3 is dead on both halves.** It recorded AMD2's extras as a net loss with
  the hash 72 to 92 percent of the penalty. On cubes `AMD2` is now FASTER than `AMD1`, 5.45 ms
  against 5.69 at 26 a side and 12.33 against 13.03 at 32, so carrying aggressive absorption and
  hash supervariable detection is cheaper than not carrying them. The fill half was reversed on
  2026-08-08 by the entry-4 filing defect; both halves were measurements of defects rather than of
  mechanisms.
- **The families have swapped roles.** The extras were free in 2D and ruinous on cubes; they are
  now free on cubes and cost about 25 percent in 2D, `AMD1` 17.47 against `AMD2` 19.86 and `AMD3`
  21.80 at 400 a side. Whatever remains of the amd gap is a SQUARE-grid effect, which is the
  opposite of where these tables pointed a few hours earlier.
- **And the parked constant-factor proposals need re-pricing rather than reading off the section
  above.** It deprioritized them because they change the shared quotient graph, which would help
  `AMD1` equally, and `AMD1` was the half that was already fine. That argument has been consumed:
  the extras are no longer the story, and what is left is per-element work in walks both branches
  share. Taken per pivot, the excess over the vendored routine is 70 ns in 2D and 149 on cubes
  against 6.28 and 12.72 clique members per pivot, so the residual tracks elements walked rather
  than pivots.

**AMD2's fill moves and AMD3's does not**, which is a tie-break consequence rather than a quality
one: `Amd3` refiles every survivor after the hash in a fixed clique order and comes out canonical,
where `Amd2`'s last write to a bucket is the merge's own refile, so the hash partition reaches the
degree buckets. `Amd2` reads `+1.4` percent of fill at 140x140 and `-3.1` at 26^3, two-sided, and
it now beats the vendored routine at six of the seven square sizes, tying at the smallest, and at
four of the six cubic ones. Every AMD2 and AMD2B fill figure recorded above predates this and should
be read with it. `AMD3`'s column is unchanged and still exact, `make amdorder` matching on all 38
cases.

`docs/DESIGN_DECISIONS.md` (2026-08-09) carries the defect and the account of why five separate
oracles were blind to it; `experiments/ordering/AMD3.md` iterations 21 to 24 are the narrative and
the ledger's entry 8 the record.

## Results

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate, 2026-07-31.** Ordering time in
milliseconds, best of three; nnz(L) in entries.

```
grid          n        MMD    MMD1     AMD    AMD1      MMD nnzL   MMD1 nnzL    AMD nnzL   AMD1 nnzL
 32x32     1024       0.19    1.42    0.20    0.30         11822       11972       11900       12074
 64x64     4096       0.60    3.75    0.44    1.51         63219       71709       67200       67950
100x100   10000       0.76    6.22    0.67    2.67        186835      223806      206332      201856
140x140   19600       1.19   11.52    1.24    5.30        412921      492921      474995      455472
```

Three things to read off it.

**Ours are slower by different factors, roughly 10x for MMD1 and 4x for AMD1**, and the split is
informative rather than noise. MMD1 calls the reachable-set query once per refreshed vertex per
round, where AMD1's whole point is not to: it reads one number per clique instead. The gap tracks
allocation and query count, not the algorithms' relative complexity.

**Superseded, 2026-07-31, by the scratch-buffer pass below.** The figures above are the ones the
first measurement produced and are kept because the reasoning attached to them is what the pass
acted on. Rerun on alpamayo to replace them.

### After the scratch-buffer pass

Three allocations went out of the inner loops, all in the shared quotient graph: the exact refresh
counts the reachable set rather than materializing it (`reachableSize`), the eliminator fills a
member scratch whose capacity survives from pivot to pivot instead of building a set per
elimination, and the absorbed-clique list is taken by swap rather than by copy. Every fill figure
is unchanged, which is what says the pass was behavior-preserving; the prototypes agree with
production as before.

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate, 2026-07-31.**

```
grid          n        MMD    MMD1     AMD    AMD1      MMD nnzL   MMD1 nnzL    AMD nnzL   AMD1 nnzL
 32x32     1024       0.13    0.51    0.12    0.40         11822       11972       11900       12074
 64x64     4096       0.39    1.53    0.42    1.21         63219       71709       67200       67950
100x100   10000       0.79    2.98    0.70    2.27        186835      223806      206332      201856
140x140   19600       1.15    5.22    1.23    4.24        412921      492921      474995      455472
```

At 140x140 MMD1 fell from 11.52 ms to 5.22 and AMD1 from 5.30 to 4.24, so the ratios to the
vendored routines are now about 4.5x and 3.4x, from about 10x and 4x. MMD1 gained most, as
expected: its refresh was where the allocation was, and AMD1's never had that call to begin with.

The fill columns did not move, on either machine, which is the check that the pass changed cost
and nothing else.

### After pass 3 in AMD1

`experiments/ordering` gained a grid mode so the layers' counters could be read at a size worth
reading, and it said the AMD gap was algorithmic: `amd1` obtains `|C[c] - C[p]|` by walking each
touched clique's members, where `amd2` subtracts from a maintained clique degree and walks
incidence lists instead. On a 100x100 grid that is 74281 elements against 272646, a factor of 3.7,
against a measured speed gap of 3.4. So the pass was ported into AMD1, along with the observation
that it makes the C[p] membership stamp dead, one loop per pivot that nothing reads afterwards
(the `amd2` prototype still carries it, inherited and dead there too).

**It bought about 12 percent, not a factor of three.** On alpamayo, 2026-07-31, at 140x140: AMD1
falls from 4.24 ms to 3.72, so 3.4x the vendored AMD down to 3.0x. A Linux sandbox agreed at 3.0 to
3.1x down to 2.6 to 2.7x over three repeats each, which is the same gain read against a different
baseline. MMD1 moved from 5.22 to 5.02, which is noise and should be, the pass not touching that
driver. The permutation and the fill are unchanged everywhere, which is what the pass promised.

```
grid          n        MMD    MMD1     AMD    AMD1      MMD nnzL   MMD1 nnzL    AMD nnzL   AMD1 nnzL
 32x32     1024       0.17    0.65    0.16    0.43         11822       11972       11900       12074
 64x64     4096       0.52    1.72    0.47    1.16         63219       71709       67200       67950
100x100   10000       0.81    3.14    0.72    2.08        186835      223806      206332      201856
140x140   19600       1.22    5.02    1.24    3.72        412921      492921      474995      455472
```

That gap between 3.7x fewer elements in one scan and 1.13x less time is the useful part of the
result. It says the scan was not where the time was, so the remaining cost is spread across the
bound loop's own walks, the eliminator's compactions, and allocation, rather than concentrated
anywhere a single mechanism can reach. **A counter is a count, not a cost**, and the two part
company once the counted thing is no longer the dominant one.

**One thing measured and rejected.** A flat weight array beside the supervariable member lists, so
that the innermost loop's `weight(v)` would be a contiguous load rather than a scattered vector
header per element. Written, measured over five runs, and reverted: 2.77x against 2.76x, no
difference at all, while MMD1's ratio moved by 0.3 in the same runs without calling `weight` even
once, which is the noise floor and puts the AMD1 figure well inside it. The prototypes' rule stands
unchanged, that an array is kept when the quantity stops being derivable.

### The MMD fill gap is the tie-break, and that question is closed

Measured 2026-07-31, in `experiments/ordering`, since the answer belongs beside the layers. Four
filing orders of the same algorithm span 443997 to 513689 at 140x140, a 16 percent spread against a
19 percent gap to the vendored routine. So MMD1's fill is not missing a mechanism and is not a
defect: on a grid nearly every live vertex has the same degree, ties decide almost every pick, and
the choice among equals is very nearly the whole ordering.

Nothing was changed on the strength of it. The order that wins at three sizes loses at the fourth,
which is what an arbitrary choice looks like once it is measured, and adopting it would be adopting
an unproven heuristic. What this does close is the earlier suspicion that porting MMD2 would
recover the fill: it will not, and the number to watch when MMD2 lands is unchanged from this one.

### Where the rest of the gap went, 2026-07-31

Pass 3 buying 12 percent where its counters promised 3.7x said the remaining cost was spread
rather than concentrated. So the next measurement counted allocations instead of elements, by
replacing global `new` in a scratch probe and running one ordering. At 140x140:

```
                       allocations        ms
vendored MMD                    17      2.75
vendored AMD                     8      3.14
MMD1  (before)              110388     12.39
AMD1  (before)              111540      8.82
MMD1  (after)                66236      9.39
AMD1  (after)                66549      5.93
```

The vendored routines allocate a workspace and nothing else. Ours allocated 5.6 times per vertex,
which is the whole difference in kind between a flat pool and a container per list. Three of those
sources were ours rather than the design's:

- **The graph was built and then copied into the quotient graph**, one allocation per vertex for a
  structure the caller had just finished with. It is a sink parameter now, moved in.
- **The eliminator rebuilt each list into a shared scratch and swapped it back**, so every list
  inherited some other vertex's buffer and often had to grow it. Both passes only ever remove, so
  they now compact in place, the write cursor trailing the read one, allocating nothing.
- **Every vertex began with a member list holding itself**, n allocations before anything had
  happened, for a supervariable that is usually a singleton. Members are a chain now, three arrays
  allocated once, which is exactly the condition under which the experiment said a weight array
  would stop being redundant, since a chain no longer gives its size away.

Together, 110388 allocations down to 66236. Every permutation and every fill figure is unchanged
throughout, which is what says the pass touched cost alone.

**alpamayo, 2026-07-31, after all three.**

```
grid          n        MMD    MMD1     AMD    AMD1      MMD nnzL   MMD1 nnzL    AMD nnzL   AMD1 nnzL
 32x32     1024       0.22    0.63    0.21    0.43         11822       11972       11900       12074
 64x64     4096       0.65    1.66    0.61    1.04         63219       71709       67200       67950
100x100   10000       0.92    2.71    0.77    1.63        186835      223806      206332      201856
140x140   19600       1.27    4.11    1.25    2.71        412921      492921      474995      455472
```

At 140x140 that is MMD1 at 3.2x the vendored MMD and AMD1 at 2.2x, against 9.7x and 4.3x when the
benchmark was first run this morning. The whole of that came from allocation and from one
algorithmic pass, with nothing about either ordering changing: the permutations are the same
permutations they were at the start.

### Flattening the adjacency, 2026-07-31

The three lists are not alike, and the difference decides which of them can leave the allocator.
**A[u] only ever shrinks**: it starts as A's off-diagonal column and pruning removes from it, so an
entry's home never moves. That is exactly the condition under which a flat array needs none of the
apparatus the vendored pool carries, no copy on growth, no compaction, no `pfree` and no `ncmpa`,
since all of that exists for lists that grow. It is now one array sized once from the pattern, with
a fixed offset per vertex and a length that falls, and `buildGraph` and the graph-shaped type went
with it: `QuotientGraph` takes A's pattern directly.

I[u] and C[c] stay containers, and honestly so. A vertex joins a new clique every time it is
reached, so its incidence list grows; a clique's members are not known until its pivot is reached.

```
                       allocations    ms (sandbox)    ms (alpamayo, 140x140)
MMD1  (before)               66236            9.39                      4.11
AMD1  (before)               66549            5.93                      2.71
MMD1  (after)                46638            8.81                      3.99
AMD1  (after)                46951            5.55                      2.53
```

Another 30 percent of the allocations. **And only about 5 percent of the time**, measured on
alpamayo at 140x140: MMD1 4.11 ms to 3.99 and AMD1 2.71 to 2.53, so 3.1x and 2.0x the vendored
routines. Permutations and fill unchanged again.

**That gap between 30 percent fewer allocations and 5 percent less time is the result**, and it
says allocation has stopped being what binds. The three earlier passes removed allocations that
were also work, a copy of the whole graph, a buffer regrown per list, a vector per vertex holding
one element; this one removed allocations that were only allocations. So the malloc calls
themselves were never the expensive part, and the reasoning that predicted a large win here, having
been wrong here, should not be trusted to predict one for an arena over I and C either.

What is left of the gap is more likely the indirection and the work itself: a container per list
means a pointer chase per list touched, whoever allocated it, and the vendored routines walk one
array. Measuring that means a profile rather than a counter, which is the next thing to do rather
than another guess.

### The cliques onto an arena, 2026-07-31

The remaining lists are C[c] and I[u], and they are not alike either. **A clique's members are
known exactly when it is formed**, at the moment its pivot is eliminated, and the set only shrinks
afterwards. That is enough for a bump allocator and needs none of the vendored pool's machinery:
take the next block, write the members, and it is that clique's for as long as it lives. No copy on
growth, because nothing grows; no `FLIP` encoding, because nothing has to find object boundaries
during a slide.

What it gives up is reclamation. An absorbed clique leaves a hole, which is exactly what the
vendored pool answers with a compaction pass, and we do not need the answer.

**The bound is the factor itself.** A clique's members are the pivot's column pattern at the moment
it is eliminated, written once, so the arena's total is the sum of those patterns over all pivots,
which is at most nnz(L) for the same run. A bad ordering does not break that: the fill that makes
the factor large is the same fill that makes each clique large, so both sides move together. And
the factor has to be resident anyway, so the ordering's scratch can never be the thing that decides
whether a problem runs. That argument is available to us and not to AMD, whose workspace is
supplied by the caller against a budget it does not control, which is why it must defragment and we
need not.

**Measured, it sits far below the bound**, because mass elimination stores one clique where a
supervariable of w columns produces w columns of L:

```
              n     nnz(A)      arena     nnz(L)   arena/L
MMD1  32x32    1024       4992       5189      11972      0.43
MMD1  64x64    4096      20224      22744      71709      0.32
MMD1 100x100  10000      49600      57843     223806      0.26
MMD1 140x140  19600      97440     115263     492921      0.23
AMD1 140x140  19600      97440     113593     455472      0.25
```

The ratio to nnz(L) falls as the problem grows, and against the *input* the arena is close to flat,
about 1.18 times nnz(A) at every size. So the scratch is the size of A rather than of L, and there
is nothing here to reclaim. (Measured with a temporary counter, since this is instrumentation and
not a feature; nothing in the tree carries it.)

The offsets being indices rather than pointers is what lets the arena grow at all: a reallocation
leaves every offset valid, where in `Iw` it would invalidate everything.

```
                       allocations    ms (sandbox)    ms (alpamayo, 140x140)
MMD1  (before)               46638            8.81                      3.99
AMD1  (before)               46951            5.55                      2.53
MMD1  (after)                31917            8.05                      3.53
AMD1  (after)                32258            5.02                      2.43
```

On alpamayo that is MMD1 at 2.7x the vendored MMD and AMD1 at 1.9x.

**And this time the time moved with the count**, where flattening the adjacency did not: about 10
percent on AMD1 against 5 for the adjacency, on a third fewer allocations either way. The
difference is where the list sits. Walking a clique's members is the innermost loop of the
reachable set, so it was paying a dependent load, header first and data second, on every clique of
every refresh. The adjacency is walked once per vertex per refresh, so it was paying that load far
less often. **Count where the loop is, not how many allocations there are.**

**What remains is about 1.6 allocations per vertex, and it is I[u] alone**: a vertex joins a new
clique every time it is reached, so its incidence list genuinely grows, which is the one case the
vendored pool's apparatus is actually for. It is also the list our two structures walk least, once
per refreshed vertex rather than once per clique member, so on the evidence above it is the least
promising of the three to move. Worth a profile before an arena.

**Answered 2026-08-01, and the premise above was wrong.** The profile was worth taking and said
`I[u]` was the largest single line after all, being an append rather than a walk. But it does not
genuinely grow: paired with `A[u]` it cannot, and the section above records what that made possible.
The vendored pool's apparatus is for the element patterns, which do grow, and not for this.

What is left is structural rather than incidental. The vendored routines hold every list in one
flat integer array while ours hold a vector per list, so what remains is per-list allocation and
pointer chasing, which is the question `experiments/ordering`'s garbage-collection section already
answered in principle: keep the containers and change the allocator, `std::pmr` over a monotonic
buffer, rather than hand-rolling AMD's pool.

**Superseded 2026-08-01**, and not by an allocator. The lists were not the problem the shape of
this paragraph assumed: `A[u]` and `I[u]` never needed to be two, and once they are one there is no
per-list allocation to redirect anywhere. `pmr` is off the list without having been built.

**The fill gap is not symmetric.** MMD1 runs about 19 percent above MMD at 140x140, while AMD1
comes out about 4 percent below AMD. The missing pieces are not equivalent: MMD1 lacks the prepass
and the pairwise merging of `mmdupd`'s q2h path, which cost real fill once a graph is large enough
for them to fire often, where AMD1 lacks mechanisms that trade fill for speed in the other
direction. The two-percent agreement seen on small grids does not survive size, which is worth
knowing before quoting it.

**The fill columns are identical to a Linux run of the same benchmark, digit for digit**, under a
different compiler and standard library. That is what an ordering should be, there being no
floating point anywhere in it, but it is now observed rather than assumed, and it makes fill a
usable regression signal across machines in a way timings never are.
