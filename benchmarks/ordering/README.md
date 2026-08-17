# Ordering Benchmark

> **SUPERSEDED IN PART, 2026-08-15.** Every ordering TIME and every ratio against a vendored
> routine in this document was measured before the encoding work of that date, and understates
> our orderings by roughly 20 to 30 percent. `MMD3` moved from 1.35 to 1.48x genmmd on square
> grids to 1.02 to 1.19x, and to 0.81 at 32 cubed; `MMD2`, `AMD2` and `AMD3` all moved with it.
> **Fill figures are unaffected**, nothing about what is computed having changed: every
> permutation and every nnz(L) is identical. The tables are left as they stand because a dated
> measurement is a record of a run. `docs/DESIGN_DECISIONS.md` (2026-08-15) has the account.

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
make phases2d, phases3d    what the vendored AMD's own time is made of, per size
make digest-record, digest    did any driver's permutation move? see below
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

## The digest, 2026-08-16: a change detector, not a correctness check

```
make digest-record     # BEFORE touching anything
   ... edit ...
make digest            # did any driver's output move?
```

`digest-record` runs all eight of our drivers over 73 small grids, squares from 2 to 55 a side and
cubes from 2 to 20, hashes each permutation and writes 584 lines keyed by driver and grid into
`.digest-baseline`. `digest` runs the same drivers over the same grids and compares. Half a second
either way. It reports WHICH driver moved and at which size, which is what makes it useful for a
change to `QuotientGraph`, where eight drivers read the code being edited.

### What it is not, and this is the whole of how to use it

It says "same as last time". It says nothing about correct. The four cases:

| | previously correct | previously wrong |
|---|---|---|
| **moved** | broke it | possibly fixed it |
| **unchanged** | still correct | still wrong |

**It distinguishes the rows and is silent about the columns.** The columns are anchored by
`make amdorder` and `make mmdorder` in `experiments/ordering`, which compare our elimination order
entry for entry against `AMD_2` and genmmd. Those say correct; this says unmoved. Between two
anchor runs you can make twenty changes and know at each step that nothing drifted, which is
exactly what the encoding folds of 2026-08-15 and 08-16 needed.

**THE FAILURE MODE IS RECORDING AFTER A BREAK**, which certifies the break forever. It is the same
shape as the two dead-code leftovers that survived a counter-based revert check on 2026-08-15: an
instrument measuring the wrong thing, confidently. Record first; re-anchor with `make amdorder` and
`make mmdorder` at the start and end of a stretch of work.

**Two drivers are a partial exception.** `MMD3` reproduces genmmd's permutation exactly and `AMD3`
reproduces `AMD_2`'s, so for those two a baseline recorded from a state that passed the anchor runs
IS a proxy for the vendored answer. The other six have no such anchor and are only ever "unmoved".

### Why small grids, and why the baseline is gitignored

**Small and many, which is the opposite of the timing ladder.** A re-encoding that breaks anything
breaks it on a small graph, usually the smallest where the shape first appears; coverage of shapes
catches it, not size. 2 a side earns its place on its own: a 2x2 grid is a clique, and an ordering
that mishandles a fully connected component fails there and nowhere else.

**The baseline is gitignored on purpose**, and the reason is a failure mode rather than a
preference. A committed baseline catches drift over months, which is real value; but when it fails,
the cheapest response is to regenerate it, and an oracle that can be quietly re-blessed is worse
than none. `make amdorder` and `make mmdorder` cannot be re-blessed, so that is where the durable
check belongs.

**A driver added after recording has no baseline entry**, and the run says so rather than passing
silently on rows it never checked.

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

**A SECOND generated copy answers the timing question instead**, since that objection is about the
raw copy rather than about the idea. `make phases2d` splits `amd_order`'s own time into five phases
and reports what share of it is the region Oblio is comparable against, which is the correction
every ratio in this file needs and has never had. See "The vendored routine's own phases" below.

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

**The AMD3 row is superseded twice below and is kept because it is what produced the finding.**
The hash defect was found by taking the count this table asked for; the section after next gives
1.44x on cubes once it was fixed, "What both families look like now, 2026-08-10" gives 0.98 to
1.33 after the two fusions, and 2026-08-14 gives 1.09 to 1.25. Read the row as the evidence that
sent us looking, not as a current figure. The MMD rows have held.

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

## What each family is actually made of, 2026-08-09

With the hash key fixed, the amd branch was profiled on both families for the first time. Every
earlier profile here was square, because until this day `order_profile.cpp` could not build
anything else. **CPU Counters, alpamayo, 140x140 at 3000 repeats and 26^3 at 600, each
configuration recorded twice**; cycles reproduced within 0.8 percent between passes and the two
cubic vendored runs agreed to seven digits, so these are results rather than readings.

`useful cycles = cycles x useful%` is the work; `cycles / useful cycles` is how efficiently it was
done.

```
                      cycles    useful   useful cycles   stall cycles   work vs AMD   eff vs AMD
2D 140x140
  vendored AMD       15.036 G   57.06%      8.580 G        6.456 G         1.000        1.000
  AMD1               22.408 G   45.36%     10.164 G       12.244 G         1.185        0.795
  AMD2               25.434 G   48.68%     12.381 G       13.053 G         1.443        0.853
3D 26^3
  vendored AMD       10.677 G   25.11%      2.681 G        7.996 G         1.000        1.000
  AMD1               14.999 G   25.07%      3.760 G       11.239 G         1.403        0.999
  AMD2               14.078 G   28.45%      4.005 G       10.073 G         1.494        1.133
```

**The two families fail in completely different ways, and that is the finding.**

**3D is pure work at identical efficiency.** `AMD1` is 25.07 percent useful against the vendored
25.11, which is as close as two measurements get, and its cycle ratio and work ratio agree to three
digits, 1.405 and 1.403. **We execute 40 percent more work and stall exactly as much per unit of
it.** No layout change, no width change and no locality work can touch that number. The only lever
on this family is doing less.

**2D is mostly stalls.** `AMD1` spends 7.37 G extra cycles of which only 1.58 G is extra work, so
**four fifths of the gap is not work at all**, and its Instruction Processing bottleneck is 28.1
percent against the vendored routine's 15.3, which is the back end waiting on operands. Split
multiplicatively, 1.19 of the 1.49 is work and 1.26 is lost efficiency.

**And the two are not one phenomenon seen twice.** Work goes 1.19x in 2D and 1.40x on cubes;
efficiency goes 0.79x in 2D and 1.00x on cubes. Something we do scales worse with clique size than
what `AMD_2` does, and separately something in our 2D memory behavior costs a quarter of the time.

### Why we close on cubes, which is not what it looks like

The tempting reading is that we handle dense structure better. The counters say otherwise:
**both codes collapse**, from 57 and 45 percent useful down to 25, with Discarded at about 40
percent for both. Misspeculation on cubes is a property of the family rather than of either code,
and the vendored routine has further to fall. We close because their efficiency advantage
evaporates, not because ours appears.

### What the extras do, which refuted the hypothesis they were measured to test

`AMD1` lacks aggressive absorption, which `AMD_2` has, so the cubic work gap looked like a missing
mechanism. The prediction was that `AMD2`'s useful cycles would fall toward 1.1x. **They rose, to
1.49x.** The extras add work, as they must.

`AMD2` is faster than `AMD1` on cubes anyway, 5.47 ms against 5.72, and the counters say why: it
**stalls 1.17 G cycles less while doing 0.25 G more work.** Aggressive absorption pays by
shortening the lists every later walk touches, which is a stall saving and not a work saving. In 2D
the same trade loses, `AMD2` stalling slightly more than `AMD1` and doing 22 percent more work, and
the timings agree, 2.05 ms against 1.77.

**One number is worth keeping separately.** `AMD2` on cubes is MORE efficient than the vendored
routine, 28.45 percent against 25.11. That is the first measure of any kind on which one of ours
beats `AMD_2`.

### What each family wants, which is the point of the section

- **2D wants memory work**, and there are two candidates, both already written down and neither
  yet tried. `mEliminated` deleted by giving cliques their own mark space removes one dependent
  byte load per element from the hottest walk in the ordering; `docs/TODO.md` carries it and
  `experiments/ordering/REPORT.md` parked it. And the width question is the other: `REPORT.md`
  measured `std::size_t` counts at 17 to 26 percent against int32, deliberately kept, and 26
  percent is what the efficiency gap here comes to. `AMD_2` runs on 32-bit `Int` throughout and
  streams half the bytes through the same caches. Both are hypotheses; the second is the cheaper
  to test and the harder to act on, since the integer model is a design decision rather than a
  tuning knob.
- **3D wants less work**, and nothing else will move it. The key walk is the visible candidate,
  since `AMD_2` accumulates its hash key inside walks it already makes and we build it in a pass of
  its own. That fusion was measured at zero on 2026-08-08 and reverted, but it was measured when
  the pair loop was ninety percent of the pass and on squares only, which is the family where the
  trade loses; `AMD1B` and `AMD2B` are the same bargain and are ahead by up to 11 percent on the
  cubic ladder while behind by 3 percent in 2D.

### The two fusions, and a null result worth as much as the finding above

`AMD_2` accumulates the weighted clique size inside the loops that build the element, `degme +=
nvi`, and takes its minimum degree inside the loop that restores the degree lists. All four amd
drivers did the first in a pass of their own and three of them did the second, so both were ported
on 2026-08-09: `QuotientGraph::cliqueWeight()` accumulates in the walk that already stamps each
member, `massEliminate` decrements it as members leave, and the minimum folds into the refile loop.
Two passes over `C[p]` per pivot removed, about 13 scattered loads per pivot on cubes and 6 in 2D.

**Measured, it bought nothing.** Same repeat counts, same session, before against after:

```
                    cycles          useful cycles        work vs AMD
3D 26^3
  vendored AMD   10.677 -> 10.686    2.681 -> 2.679     1.000
  AMD3           15.164 -> 15.144    4.338 -> 4.324     1.618 -> 1.614
  AMD2           14.078 -> 14.136    4.005 -> 3.995     1.494 -> 1.491
2D 140x140
  vendored AMD   15.036 -> 15.023    8.580 -> 8.569     1.000
  AMD3           28.107 -> 28.179   13.356 -> 13.390    1.557 -> 1.563
  AMD2           25.434 -> 25.400   12.381 -> 12.361    1.443 -> 1.442
```

Everything is unchanged within half a percent. **The reason is that the loads were already cheap**:
`pivotClique` is a contiguous arena run and `mWeight` is small and hot, so those gathers were
served from L1 and cost near nothing. That is this folder's own lesson arriving from the other
direction. Counting elements walked cannot price them, and here the count was right and the work it
counted was already free.

**The change is kept as measured-neutral rather than as an improvement.** It is a faithful port, it
removes two passes and a stale-value hazard, and it costs nothing; it is not a speed fix and should
not be cited as one. It also carried a correctness dependency: `Amd1`, `Amd2` and `Amd2B`
mass-eliminate inside the eliminator and so needed `cliqueWeight` to follow the trim, which is
`AMD_2`'s own `degme -= nvi`. Without it they read a bound too large per merged vertex. **`make
amdorder` and all 283 assertions passed; only `prototype and production agree` caught it**, `Amd3`
mass-eliminating late and so being legitimately untrimmed.

**And a caution about this benchmark that cost a wrong claim.** The scale tables appeared to show
`AMD3` improving from about 1.41 to 1.30 on the cubic ladder after the fusions, and that was
**drift**. `AMD3` at 26 a side has measured 5.83, 5.95 and 5.65 ms across three runs of the same or
nearly-same code, a 5 percent spread, where the cycle counts for those runs agree to 0.2 percent.
**Ratios in these tables move by several percent between runs at fixed code**, so a change that
does not exceed that band is not visible here at all, and the counters are the instrument for
anything smaller. The harness already picks its repeat count from a timed probe for this reason and
it is still not enough at these sizes.

**So the 1.6x work gap is intact and unexplained.** On cubes it is 2.7 G cycles of real
instructions per 600 orderings, and it is not in any pass that reading the code has found. What
that wants is the per-pass inventory, ours against `AMD_2`, element visits per pass per pivot on
both families, counted rather than reasoned. That is the instrument that found the hash key, and it
is the only one that has answered anything on this question.

**And the parked proposals are re-priced by this rather than merely revived.** `NEXT.md`
deprioritized them because they change the shared quotient graph and would help `AMD1` equally,
where `AMD1` was thought to be fine. `AMD1` is now the larger of the two remaining terms and its 2D
gap is four fifths stalls, which is exactly what those proposals address. The same proposals are
worthless on cubes, where the efficiency gap is zero.

## Two things this benchmark now does to itself, 2026-08-10

**It asks for a performance core.** `pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0)`
at the top of `main`, guarded by `__APPLE__`, in both drivers. A command-line process on Apple
Silicon runs at `QOS_CLASS_DEFAULT`, which prefers a performance core but permits the scheduler to
park the thread on an efficiency one, and that placement is STICKY over long stretches rather than
jittering per iteration. Every row here is a MINIMUM over fifteen to thirty repeats, which filters
per-sample noise completely and filters a whole run placed on the wrong core not at all.

That is the shape of the 4 percent disagreement between identical binaries this file has recorded
since 2026-08-08 and treated as irreducible. **It was not irreducible.** With the call in, a change
that had read as a scattered null across five cubic sizes read as a clean 3 to 9 percent, and the
verdict on it reversed. `taskpolicy -c` at the shell does the same thing and has to be remembered
every run, which is why it is in the source instead.

**And it measures its own floor when it can.** See the section below.

## The noise floor of this benchmark, measured 2026-08-10

**Plus or minus 3 percent, at every size above the smallest.** Between the key fusion landing in
`Amd3` and the `Amd3B` vehicle being removed from the tree, that driver held a VERBATIM COPY of
`Amd3`. Its column in `make scale2d` and `make scale3d` is therefore the same code timed twice in
one run, and the difference between the two is the instrument measuring itself:

```
identical code       AMD3    AMD3B              AMD3    AMD3B
32x32                0.09     0.08  -11.1%      6^3     0.03     0.03    0.0%
64x64                0.38     0.38    0.0%     12^3     0.29     0.29    0.0%
100x100              1.00     1.02   +2.0%     16^3     0.94     0.92   -2.1%
140x140              2.04     2.02   -1.0%     20^3     2.17     2.25   +3.7%
200x200              4.23     4.34   +2.6%     26^3     5.23     5.29   +1.1%
280x280              8.69     8.49   -2.3%     32^3    11.18    11.54   +3.2%
400x400             19.42    19.71   +1.5%
```

The two smallest rows are quantization, 0.02 ms resolution on a 0.03 ms run, and say nothing. The
rest is real drift: same binary, same process, seconds apart.

**What follows for reading any other row here.** A single figure under about 4 percent is not a
result. What rescues a small effect is CONSISTENCY OF SIGN ACROSS SIZES, since drift is not
systematically signed: six consecutive sizes all negative is evidence where any one of them alone
is not. And a figure between 0 and 3 percent is a NULL, whatever its sign, and must be written up
as one.

**Free, and worth arranging deliberately.** An idle variant column costs one benchmark column and
gives every other column an error bar taken under identical conditions, which no amount of
repetition of a single column can supply. `AMD3B` is not a standing column here: it is added with a
candidate change and removed with it, so the arrangement to repeat is to run the benchmark ONCE
with the new vehicle still a verbatim copy, before putting anything into it.

## The per-pass inventory, 2026-08-10, and a correction to our half

**The instrument.** One counter per loop in both codes, on the same grids, counting two kinds of
thing separately. A **walk** is one iteration of a loop over a vertex's adjacency, a vertex's
incidence, or a clique's members. A **sweep** is one iteration of a loop whose body is one member
of the new clique. The 2026-08-08 table had a single shared elements column, which could express
touching more arrays per element and could not express executing a different number of loops.

The vendored half is generated from `private/Amd.cpp` by a script that asserts every anchor, the
arrangement `tools/hook_amd.py` uses and for the reason its header gives. Before any number is read
off it, the counted copy's permutation is compared against the unhooked `amd_order` entry for
entry, its `gc` counter against `AMD_NCMPA`, and its `AMD_LNZ` against the six fill figures this
file already records, which it matches digit for digit.

**Our half as previously recorded was wrong, by up to a factor of two on six of fifteen columns.**
The instrumented `QuotientGraph` carries the same symbols as the original, so it is linked instead
of it and both drivers share it, and the probe read its counters after a control run of the
uninstrumented driver. The passes living in the shared class were counted twice and those in the
driver once.

```
                              per pivot, corrected     as previously recorded
AMD3 element visits, 140x140        112.26                     155.15
AMD3 element visits, 26^3           276.25                     374.83
```

Two consequences worth keeping. The prune's incidence compaction is **not** oversized: corrected it
is 14.74 in 2D and 41.49 on cubes, which is `AMD_2`'s scan 2 to the digit, so the proposal that it
walks a longer list than the vendored routine because it runs before absorption describes nothing
real. And the check that caught it was free: the alignment forces five of our passes to equal a
vendored counterpart exactly, and after the correction all five do at all six sizes. **A count
against our own code cannot see a scale error; a count against a code that must agree with it can.**

**The completed table.**

```
                          ours                 AMD_2             ratio
               walks sweeps    all   walks sweeps    all  walks  sweep    all  cycles
2D 140x140     93.41  56.55 149.96   45.25  25.32  70.57  2.06x  2.23x  2.12x   1.56x
3D 26^3       238.08 114.49 352.57  115.41  51.27 166.68  2.06x  2.23x  2.12x   1.61x
```

**2.12x on both families, to two digits**, which is the first quantity on this question to come out
the same on squares and cubes. Where the excess sits, at 140 a side and at 26 cubed:

```
                                        ours   AMD_2   excess  share
2D 140x140
  I[u]: prune, absorb, bound, key      52.90   14.74    38.16    79%
  A[u]: prune, bound, key              14.70    6.68     8.02    17%
  scan 1, over I[u]                    16.41   14.74     1.67     3%
  construct and the hash passes         9.40    9.07     0.33     1%
3D 26^3
  I[u]: prune, absorb, bound, key     140.71   41.49    99.22    81%
  A[u]: prune, bound, key              35.56   14.76    20.80    17%
  scan 1, over I[u]                    43.66   41.49     2.17     2%
  construct and the hash passes        18.15   17.67     0.48     0%
```

We walk `I[u]` four times per pivot and `AMD_2` walks it twice; we walk `A[u]` three times and it
walks it once. We make nine sweeps over `C[p]` and it makes four, five of ours collapsing into its
scan 2, which per member computes the degree, accumulates the key, compacts the list and tests for
mass elimination in a single visit.

**And the ratio that reframes every null in this file.** 2.12x of visits against 1.56x of cycles
means we execute about **0.74x the work per element visit** that `AMD_2` does. Its loops are fat
because each does four things; ours are thin because each does one. A fusion therefore does not
remove work, it trades a visit for per-visit work, and pays only where that work is redundant
rather than merely distributed.

## What the inventory was worth, and the ranking it suggested was wrong, 2026-08-10

Two candidates came out of the table and went through a scratch `Amd3B` one at a time, the second
only after the first was priced. That driver is added while a candidate is being measured and
removed once it has landed or been recorded, so it is not a standing column in the tables above. Together they had measured 5 percent in 2D and 5 to 12 on cubes, and the
split is the whole reason we know what that was.

```
AMD3 -> AMD3B         B1, a sweep deleted     B2, the key folded into the bound
64x64                        0.0%                        -7.0%
100x100                      0.0%                        -4.6%
140x140                     -0.9%                        -6.5%
200x200                     +3.6%                        -5.1%
280x280                     -1.2%                        -4.4%
400x400                     +0.5%                        -5.4%
12^3                        +5.0%                       -14.3%
16^3                        +2.8%                       -12.8%
20^3                        +1.5%                        -8.1%
26^3                        +2.3%                        -5.0%
32^3                        +1.4%                        -8.0%
```

`nnz(L)` identical to `AMD3` at every size in every run. `AMD3B` carrying B2 came out at 11.20 ms
at 32^3 against `AMD2`'s 11.70 while returning the vendored permutation, which is a firmer
statement than any ratio against `AMD` here, both being ours and both having run in one process.

**B1 deletes a sweep, is provably redundant, and measured NOTHING.** `Amd3` re-sums `weight(u)`
over `C[p]` after mass elimination trims it, recovering a number `massEliminate` already maintains,
so the loop cannot compute anything the class does not hold. Removing it reads 0 percent in 2D and
1 to 3 percent slower on cubes, and the whole of that range sits inside the plus or minus 3 percent
floor measured above. **It is a null.**

The sign was positive at all five cubic sizes in two independent runs, which is weak evidence of a
small real cost and is not more than that: 5 of 5 twice over is suggestive, and a 3 percent effect
against a 3 percent floor is not resolvable here. One mechanism would explain it, that the deleted
sweep walked `pivotClique` immediately before the bound sweep walks the same array and was warming
`C[p]`, but that is a hypothesis with a null attached to it rather than a finding, and CPU Counters
would settle it in ten minutes if it ever mattered.

**What it does establish is enough for the ranking.** Deleting a provably redundant sweep bought
nothing measurable, where fusing two walks bought 4 to 7 percent. So the sweeps axis is not the
better bet, which is what the inventory had suggested; that is a statement about which of two
candidates to try next, and it needs no claim about why B1 failed.

**B2 fuses two walks into two others, had been tried and reverted, and is the only thing that has
moved this gap.** The 2026-08-08 version carried the key in a vector of size n and measured nothing
at 140 a side and minus two percent at 400. This one files each vertex into its bucket where its
key completes and stores nothing extra. **The failure was the array, not the fusion.**

**And only `Amd3` can take the fusion, which is worth knowing before it is proposed elsewhere.**
`Amd2` and `Amd2B` form the bound in one pass and call `buckets.refile` inside it, so their bound
loop's direction is already a tie-break input. Their key pass walks `C[p]` backward against that
forward bound, and head insertion into a degree bucket and into a hash chain want opposite
directions, so one walk cannot serve both without a tail pointer per bucket, which is the size-n
array the change exists to avoid. Measured on a scratch copy: the permutation moves on all ten
grids tried, with fill at `-1.33` percent at 140x140 and `+3.24` at 26^3, two-sided and small. So
it is an ordering change there rather than a schedule change. `Amd3` is free of this because
ledger entry 4 split its bound in two and moved the refile below the hash, for an unrelated reason.
`Amd1` and `Amd1B` have no hash detection and so no key.

**So the ranking the table suggested does not hold.** Deleting a pass was preferred to rescheduling
one, on the argument above and on four fusion attempts having measured zero. The fifth fusion is
the one that worked and the one deletion was negative. What survives is narrower: **a null result
measures one implementation, not the idea it implements**, and the four earlier nulls were read as
the latter for two days.

**What is not yet known**, and it is what a profile would say. B2 removes 18 to 19 percent of the
visits and returns 5 percent in 2D. If the saving is instructions in proportion to the visits, the
ceiling on this direction is low; if cycles fall further than instructions, the win is locality
from touching each list once instead of twice, and the three remaining walks of `I[u]` in the
prune, scan 1 and the bound become the large prize, since `AMD_2` makes two where we make four.

**And what B2 does not touch**, which is larger than B2. `AMD3`'s 2D gap grows with n, 1.54x at 64
a side to 2.07x at 400, and a constant five percent off a growing curve leaves it growing. `MMD3`
over the same quotient graph shows no such trend across the same seven sizes, so whatever grows
lives on the amd branch alone and is not made of passes.

## The first scan folded into the prune, 2026-08-10

**`AMD3` reaches the vendored routine on cubic grids.** The driver walked `I[u]` three times per
pivot and `A[u]` twice; `AMD_2` walks them twice and once. `QuotientGraph::eliminate` now
accumulates |C[c] - C[p]| and the bound's adjacency term on the walk the prune is already making,
so scan 1 goes entirely and the bound's adjacency loop with it.

```
                  AMD3   AMD3B                    AMD3   AMD3B
12^3              0.29    0.26  -10.3%   64x64     0.41    0.41    0.0%
16^3              0.93    0.85   -8.6%   100x100   1.03    0.95   -7.8%
20^3              2.22    1.94  -12.6%   140x140   2.07    1.98   -4.3%
26^3              5.28    4.67  -11.6%   200x200   4.27    4.14   -3.0%
32^3             10.85   10.07   -7.2%   280x280   8.87    8.76   -1.2%
                                         400x400  20.12   19.47   -3.2%
```

`nnz(L)` identical at all thirteen sizes. It is faster than `AMD2` at every cubic size and faster
than `AMD1` too, which no layer carrying the extras has managed before.

**TWO CORRECTIONS TO THAT TABLE, both found the same evening, and both about the instrument.**

**The two columns were not timed the same way.** A scratch variant reached as a free function goes
through `orderTimeFn`, which times the bare ordering call; a standing method goes through
`orderTime`, which times `OrderEngine::compute` and so also builds a Permutation, two assigns of
size n and a loop of size n. Timing identical code down both paths puts that at **0 to 2.4 percent
in the free function's favor**, so every figure above is that much too generous. Re-measured with
both down the same path, the fold is worth **10 to 16 percent on cubes**, 16.1 at 12 a side, 16.2
at 16, 9.8 at 20, 10.3 at 26 and 12.8 at 32. Larger than first recorded, and arrived at properly.

**And the ratio's denominator is the noisiest column in the table.** Over eight runs at 16 cubed,
`AMD` reads 0.74 to 0.86 ms, a 16 percent spread, where `AMD3` reads 0.83 to 0.89, a 7 percent one.
**The vendored routine varies more than we do**, so `AMD3 / AMD` per row is mostly a measurement of
it. This section first said `AMD3` reaches 1.01x, 1.02x and 1.07x at 12, 16 and 20 a side, "which
is parity"; that was one run. What reproduces:

```
                 AMD3 ms        AMD ms        ratio over eight runs
12^3           0.26 - 0.28   0.20 - 0.27         0.98 - 1.33
16^3           0.83 - 0.89   0.74 - 0.86         0.97 - 1.18
26^3           4.65 - 4.89   3.98 - 4.13         1.13 - 1.19
32^3          10.29 - 10.75  8.43 - 8.83         1.18 - 1.27
```

So: **at or near the vendored routine to 16 a side, rising to about 1.2x by 32.** 16 cubed and 26
cubed reproduce to a percent; 12, 20 and 32 wobble by five. **Quote absolute milliseconds with the
vendored range beside them, not a ratio per row.**

**The fusion alone did not do this, and that is the part worth keeping.** The first version carried
the two values crossing from the prune to the bound in two fresh vectors of size n, and measured 3
to 9 percent faster on cubes and **12 percent slower in 2D from 200 a side up**. Both fit in arrays
the driver already had and that are dead at that point, `partial[u]` and `hashNext[u]`, with the
key reduced modulo the bucket count as it accumulates so that it fits an `int32`. With the arrays
gone the 2D penalty went with them.

**So the per-pass inventory above is half an instrument.** It counts VISITS and is silent about how
many size-n streams a change adds, and here the streams were the larger term on one family. This is
the third time the footprint has been the answer: the 2026-08-08 key fusion failed for it and was
recorded as a failure of the fusion, and `Amd1B` is on record as slower at large n after being
faster at small for the same reason. **Price a re-schedule by what it walks AND by what it makes
resident.**

## The two VENDORED routines against each other, 2026-08-10

Every other table in this file measures one of ours against the routine it was ported from. This
one measures genmmd against `AMD_2`, which is a different question and the one that says what we
should be aiming at on each family. Fill is exact and reproduces to the digit; **the times are one
run each and this benchmark moves 7 to 11 percent between runs**, so read the ratios and not the
milliseconds.

```
SQUARE, make scale2d
 side       n   MMD ms   AMD ms  MMD/AMD      MMD nnz    AMD nnz  MMD/AMD   ns/n MMD  ns/n AMD
   32    1024     0.05     0.06    0.83x        11822      11900   0.993x       48.8      58.6
   64    4096     0.21     0.27    0.78x        63219      67200   0.941x       51.3      65.9
  100   10000     0.54     0.62    0.87x       186835     206332   0.906x       54.0      62.0
  140   19600     1.09     1.22    0.89x       412921     474995   0.869x       55.6      62.2
  200   40000     2.27     2.52    0.90x       981766    1081911   0.907x       56.8      63.0
  280   78400     4.46     4.93    0.90x      2137410    2440757   0.876x       56.9      62.9
  400  160000     8.88    10.26    0.87x      4862612    5663298   0.859x       55.5      64.1

CUBIC, make scale3d
 side       n   MMD ms   AMD ms  MMD/AMD      MMD nnz    AMD nnz  MMD/AMD   ns/n MMD  ns/n AMD
    6     216     0.02     0.02    1.00x         3279       3265   1.004x       92.6      92.6
   12    1728     0.47     0.20    2.35x        75674      76038   0.995x      272.0     115.7
   16    4096     1.65     0.76    2.17x       295113     281014   1.050x      402.8     185.5
   20    8000     3.60     1.76    2.05x       840560     842282   0.998x      450.0     220.0
   26   17576     8.76     4.02    2.18x      2869267    2836813   1.011x      498.4     228.7
   32   32768    17.63     8.83    2.00x      7898321    7746501   1.020x      538.0     269.5
```

**The two advantages are made of different things.** In 2D genmmd fills 6 to 14 percent less from
64 a side up and the two are level on time. On cubes `AMD_2` orders twice as fast and the fill is
equal, within two percent at every size and with no direction: MMD is ahead at 12 and 20 a side and
behind at 6, 16, 26 and 32.

**Neither ratio moves with n on either family.** 0.87x time and 0.86x fill in 2D at 400 a side are
what they were at 64; 2.0x time and 1.02x fill on cubes at 32 are what they were at 12. So both are
CONSTANT FACTORS, and whichever routine wins at the smallest size wins at the largest. The choice
between them is a property of the problem's dimensionality rather than its scale, at least out to
these sizes.

**Ordering is near linear in n and fill is not**, which is why the ordering phase shrinks as a
share of a solve:

```
                   time        fill
square grids     n^1.01      n^1.18
cubic grids      n^1.15      n^1.58
```

Fitted over the reliable part of each ladder, n = 10^4 to 1.6*10^5 on squares and 4096 to 32768 on
cubes. Both routines give the same exponents to two digits on both families, which is worth knowing
on its own: **they scale identically and differ by a constant.**

**And the largest number here is not a branch difference at all.** Cubes cost 4 to 8 times more
ordering time per vertex than squares at comparable n, 229 ns against 62 for `AMD_2`. That is
clique size: 12.7 members per pivot on cubes against 6.28 in 2D, which the pass inventory above
measures directly. The family gap dwarfs either branch-against-branch gap.

## What both families look like now, 2026-08-10

With the cubic constant factor largely gone, the two families have the same SHAPE for the first
time, and it is not the shape either had before.

```
cubes, AMD3 / AMD        2D, AMD3 / AMD
12^3   0.98 - 1.33       32x32     1.33
16^3   0.97 - 1.18       64x64     1.39
20^3   1.05 - 1.17       100x100   1.56 - 1.73
26^3   1.13 - 1.19       140x140   1.63
32^3   1.18 - 1.27       200x200   1.63 - 1.64
                         280x280   1.68 - 1.75
                         400x400   1.87 - 1.97
```

**A low constant plus a term that grows with n, on both.** The constants differ, near 1.0 on cubes
and about 1.33 in 2D, and the 2D one is the stall penalty already recorded. What is new is the
growth: cubes rise about 25 percent from 12 a side to 32, where the cubic gap was FLAT at 1.3x that
morning and flat at 3.0x three days before. 2D rises about 45 percent from 32 to 400 with a visible
knee past 280, where the working set roughly doubles.

**That is a different question from the one the pass inventory answered.** The inventory addressed
a constant factor and the constant factor is nearly spent. Nothing measured so far touches a term
that scales, and the knee says memory rather than instructions. `MMD3` over the same quotient graph
shows no growth on either family, which is the control that makes this an amd-branch property
rather than a shared-infrastructure one.

## The ladders as they now stand, and how to read the tables, 2026-08-16

**Square**: `32 50 64 100 128 200 256 400 512 800 1024 1600`, twelve points over a 2560-fold range
in n. **Cubic**: `8 10 16 20 32 40 64 80`, eight points. Each ladder is TWO INTERLEAVED SERIES,
`32 64 128 256 512 1024` against `50 100 200 400 800 1600`, and `8 16 32 64` against `10 20 40 80`,
each quadrupling n and offset from the other.

**Why interleaved.** A power-of-two side aligns every array length, so a trend measured on such a
ladder alone cannot be told apart from an addressing artifact. The second series separates them,
and it earned itself on the first run. Note that a power-of-two SIDE gives a power-of-two n in both
families: `m^3 = 2^3k` when `m = 2^k`, so the cubic ladder aliases exactly as the square one does.

**80 cubed is nnz(A) 3545600 against the 800 square's 3196800**, so the two ladders meet at the top
on input size, which is a more meaningful axis than n. **8 cubed is n = 512** and measures startup
as much as ordering; it anchors the low end of its series rather than being read alone.

**MMD1 and AMD1 are out of the scaling lists**, symmetrically so the tables keep their shape, and
still in `allMethods` so `run2d` and `run3d` show them. `MMD1` costs about 40 seconds per ordering
at 80 cubed, five per row plus a fill; the mmd cubic ladder went from minutes to 46 seconds without
it. `AMD1` is not expensive and left only for symmetry. One line each to restore, and they are what
says how much supervariables and aggressive absorption are worth.

### The rule for reading any ratio column here

**Convert to time per vertex before believing a trend**, and treat the STARRED rows as suspect:
only the vendored AMD zigzags, and it does so because `AMD_1` carves six arrays of exactly n ints
out of one block. At 64 cubed it costs 536 ns per vertex against 314 at 80 cubed, which has half
again as many vertices. The unstarred rows are the honest baseline.

### Slopes, twelve square sizes, `time ~ nnz(A)^alpha`

```
                 all 12    powers of two    other sides
MMD (genmmd)      1.062          1.071          1.056
MMD2              1.052          1.051          1.055
MMD3              1.052          1.047          1.058
MMD3B             1.054          1.054          1.056
MMD1              1.169          1.214          1.138

AMD (vendored)    1.084          1.106          1.080
AMD2              n/a            n/a            n/a
AMD3              1.095          1.107          1.088
AMD3B             1.066          1.081          1.056
```

**The mmd branch has no growth term**: every layer but `MMD1` sits at 1.05 against genmmd's 1.06,
with no series split, and what `MMD3` carries is a constant of 1.04 to 1.20 that erodes with n.
**`MMD1` is the one genuinely different slope in either family**, 1.169, and its ratio to genmmd
climbs from 1.60 at 32 a side to 2.89 at 800.

**On the amd branch `AMD3B` is now BELOW the vendored routine on both series.** It reads 0.93x at
1024 squared and 0.97x at 1600, the first square sizes where anything of ours has beaten it. See
`docs/DESIGN_DECISIONS.md` (2026-08-16, later still) for what did it.

## The wider ladders, and the one column that zigzags, 2026-08-16

**The ladders are two interleaved series each.** Square is `32 50 64 100 128 200 256 400 512 800`,
cubic `10 16 20 32 40 64`. Within each, `32 64 128 256 512` and `16 32 64` quadruple n, and so do
`50 100 200 400 800` and `10 20 40`. Ten square points over a 625-fold range in n where there were
seven over 156-fold.

**The interleaving is the point, not the extra rows.** A power-of-two side aligns every array length
and every grid stride, so a trend measured on such a ladder alone cannot be told apart from an
addressing artifact that grows with n. A second series that quadruples n identically while never
aligning separates them.

### Read these tables per vertex, not as a ratio

```
side       32    50    64   100   128   200   256   400   512   800
AMD      58.6  56.0  63.5  59.0  72.0  64.0  75.2  64.0  96.7  77.1     ns/vertex, zigzags
AMD3     68.4  76.0  78.1  82.0  84.2  87.2  89.3  94.2 108.1 127.2     smooth
genmmd   48.8  52.0  51.3  55.0  55.5  58.0  62.9  57.1  61.5  66.2     smooth
MMD3     58.6  60.0  61.0  60.0  59.8  61.5  61.3  61.6  64.0  72.1     smooth
```

**Only `AMD` zigzags**, costing more per vertex at every power-of-two side than at the larger
50-series side beside it. So the two-series pattern in the RATIO columns is a denominator effect,
and the ratio falling at 32, 64, 128, 256 and 512 means `AMD` got worse, not that we got better.

**A ratio hides which side is moving.** This survived a full differential and seven folds while
being read as a property of our code. If a column here looks like a trend, convert it to time per
vertex before believing it.

### The slopes, which is what the ladder was built to support

Fitted least squares on log-log against nnz(A), all ten square sizes, R2 at or above 0.998
everywhere. `time ~ nnz(A)^alpha`:

```
MMD  (genmmd)   1.039        AMD  (vendored)  1.054   confounded, two series
MMD3            1.018        AMD3             1.077

AMD  powers of two   1.080   AMD3  powers of two   1.071
AMD  other sides     1.049   AMD3  other sides     1.081
```

**Everything is superlinear, and the exponents are small.** Read these as slopes to compare, not as
complexity claims: the published bounds are worst-case and dense, and an exponent near 1 is what
theory expects on a grid. `experiments/ordering/README.md`, "What the literature proves about these
algorithms", has the bounds and why they do not bind here.

**`AMD3` has ONE slope**, 1.071 and 1.081, indistinguishable between the series. **`AMD` has two**,
1.080 aligned and 1.049 unaligned. So at power-of-two sides the vendored routine grows at our
exponent and elsewhere at a lower one: the self-aliasing described below costs it about **0.03 in
the exponent**, which is the same finding stated as a slope rather than as a ratio.

**Our amd excess is 1.081 against 1.049**, about 0.032 in the exponent, or 1.23x over the full
range. Real, small, and the sharpest form the open question has taken.

**And `MMD3`'s slope is BELOW genmmd's**, 1.018 against 1.039. We grow more slowly than the
reference on that branch. The container overhead is a constant paid at every size and it is being
eroded: `MMD3` reads 1.20x genmmd at 32 a side and 1.09x at 800.

**Two caveats on the fits.** The smallest rows are quantised, 0.05 to 0.07 ms printed to two
decimals being plus or minus ten percent; refitting from 100 a side up moves every alpha by at most
0.02 and preserves every ordering. And these are wall-clock, so they include the memory effects. An
instruction-count fit under cachegrind would separate algorithmic growth from growth in cost per
instruction, and cachegrind now runs in the sandbox, so that is cheap.

### The cause, established by intervention

`AMD_1` carves `Pe`, `Nv`, `Head`, `Elen`, `Degree` and `W` out of one block at offsets that are
exact multiples of n. `n = m^2`, so a power-of-two side gives a power-of-two n and the six alias
each other in the cache. Cachegrind, one `amd_order` per run: instructions and data reads per vertex
are FLAT across 400, 512 and 800 to a tenth of a percent, while D1 read misses per vertex read 15.3,
40.3 and 17.0. Padding the six arrays apart by one cache line, which changes addresses and nothing
else, removes 56 percent of the misses at 512 and none at 400, with byte-identical permutations.
See `docs/DESIGN_DECISIONS.md` (2026-08-16, later).

**The honest baseline for growth is therefore `AMD`'s UNALIGNED series**, and the aligned rows
should not be read as our columns improving.

### What the two ends measure

At 32 a side `AMD3` is 1.17x `AMD` per vertex and `MMD3` is 1.20x genmmd. That is the
`std::vector` layer, measured twice independently.

Over a 256-fold range within a series, `AMD3` grows 1.55 to 1.58x and `AMD`'s unaligned series grows
1.38x, so about **1.13x of growth is genuinely ours** and is not yet explained. `MMD3` grows 1.09
and 1.20x against genmmd's 1.26 and 1.27x, so the mmd branch has no such term.

### What the top rows cost

At 64 cubed, n = 262144, nnz(L) reaches 247 million for `AMD1` and 283 million for `MMD1`; the
factor is never materialized, only summed, so peak memory stays modest. `scale-amd-3d` runs in
seconds. **`scale-mmd-3d` does not**: `MMD1` reads 49x genmmd there with 45 percent worse fill and
seconds per ordering, which is most of that target's run time.

## The amd branch's 2D growth, and the one change that removed it, 2026-08-16

**The symptom this file recorded for two weeks.** `AMD3` read 1.25x the vendored routine at 32 a
side and 1.82x at 400, rising monotonically, while `MMD3` over the same quotient graph was flat and
`AMD3`'s own cubic column was nearly flat. Six folds were applied to the amd drivers over one
session; five of them moved the whole 2D column down by a constant and left the slope exactly where
it was.

**The one that removed the slope** was moving a clique's descriptor into the dead pivot's own
`mRun` entry, retiring `mCliquePtr` and `mCliqueSize`. Before it, every clique visit in a walk
probed two separate n-arrays at a dead pivot's id, once per element of every `I[u]`; after it, one
16-byte load on a line the walk already touches.

```
AMD3B against AMD    32    64   100   140   200   280   400
before             1.11  1.19  1.26  1.42  1.39  1.43  1.49     rising
after              1.26  1.15  1.38  1.33  1.38  1.38  1.38     flat
```

**Why 2D and not cubes, and it is not grid shape.** The differential in
`experiments/ordering/README.md` shows a 2D pivot doing about 60 visits and a cubic one 149. A
fixed per-pivot cost is therefore 2.5 times more visible in 2D, at the SAME n. That also explains
why the cubic column never showed the symptom and why every fold has paid more here than there.

**What did not cause it.** Recorded because each was a plausible reading and each was killed by
measurement rather than by argument: it was not extra work, the two codes doing the same visits per
pivot at every size in both families; not the clique arena, whose excess misses FALL with n and
which is smaller than the vendored workspace in 2D; not the hash; not `clear_flag`, which fired
zero times on every case; and not single-level cache locality, an L1 model putting the excess at a
constant 0.09 misses per visit, real but far too small at sizes where everything fits in L2.

### What to run, and what a reader should check

`make scale-amd-2d` and `make scale-amd-3d`. The columns to read against each other are `AMD3` and
`AMD3B` while both exist; they are the same algorithm and should agree within the plus or minus 3
percent floor. The ladder's default sides are geometric enough to show a trend but sparse in the
middle, and the driver takes explicit sides, so a denser ladder needs no code:

```
./order_timing_cpp amd 2d 32 45 64 90 128 181 256 362 400
```

**And the open one.** `AMD1` and `AMD2` still climb in 2D, 1.05 to 1.22x and 1.36 to 1.47x, where
`AMD3` is flat. The descriptor fold is in the shared class and reached them for free, so whatever
remains is driver-side. That is the next thing the differential should be pointed at.

## The vendored routine's own phases, and the denominator every ratio here divides by, 2026-08-15

**`make phases2d` and `make phases3d`.** Every time ratio in this file is against `amd_order`
whole, and two of that routine's phases have no counterpart on our side:

- **`AMD_aat`** forms the pattern of `A + A'`, because AMD takes a one-sided matrix. Ours arrives
  full-symmetric with the diagonal present, which is what a `SparseMatrix` holds and what
  `QuotientGraph` reads directly.
- **`AMD_postorder`** relabels the assembly tree. `ElmForestEngine` postorders the exact supernodal
  tree later with real front and update sizes, so the vendored one is redone and discarded, and
  `Amd.cpp`'s own header says that tree "is not guaranteed to be the precise supernodal elimination
  tree" in any case.

So `AMD3 1.82x` at 400 a side is measured against more work than the comparison is about. **This
file has carried one estimate of the correction since 2026-08-01**, 139 ms and 154 ms of a 3.67 s
run at 140 a side, about 8 percent, from a single profile at a single size, with the conclusion
"the true gap is nearer 1.6x than 1.46x" resting on it. The estimate has never been measured per
size, and that is the part that matters: **a correction that FALLS with n would mean part of the
amd branch's 2D growth, 1.25x at 32 a side to 1.82x at 400, is the baseline shrinking rather than
us.**

### The five phases

```
valid   AMD_valid, amd_preprocess where the input is jumbled, the Len and Pinv vectors
aat     AMD_aat                                                     WE DO NOT DO THIS
build   the S workspace, and AMD_1 filling Iw and Pe from that pattern
core    AMD_2, from entry to the end of its main loop
post    assembly-tree path compression, AMD_postorder, the output permutation
                                                                    WE DO NOT DO THIS
comp    build + core, the comparable region
```

**`build` belongs with `core` rather than with the excluded pair**, which is the one judgment in
the split. It is the vendored routine turning a caller's pattern into its own working structure,
and `orderAmd3` does the same thing in `QuotientGraph`'s constructor, so excluding it would flatter
us by exactly the phase we have most recently been folding arrays out of.

### It is a second generated copy, not a second use of the first

`tools/hook_amd.py` now has two modes. The **raw** copy accumulates supervariable membership and
emits the elimination order, which is what `make amdorder` compares against; it carries a vector
per vertex and a push per member, so it is an oracle and timing it would measure the bookkeeping,
which is why this file has always declined to report its time. The **timed** copy carries five
timestamps and nothing else. Both are generated from whatever `private/Amd.cpp` currently says,
both assert every anchor, and both are gitignored and removed by `make clean`.

**The postorder is TIMED rather than removed, and the reason is worth recording because deleting
the call is the first thing to reach for.** `AMD_postorder` runs inside `AMD_2`, and the block after
it builds the output permutation out of the ranks it writes into `W`. Delete the call and `W` is
unwritten, the returned permutation is meaningless, and both of the checks that make a generated
oracle trustworthy go with it.

### Two checks, and one caveat that cannot be checked away

**The copy must reproduce the shipped routine**, and the row says so and refuses rather than
printing figures if it does not: the permutation entry for entry, `Info[AMD_LNZ]` and
`Info[AMD_NCMPA]`, all against unhooked `amd_order` on the same pattern, once per row before
anything is timed. `Control` is `nullptr` in both, which is what `OrderEngine` passes for
`Ordering::AMD`, so the phases describe the same call the `AMD` column measures.

**And the caveat: a clock read is a compiler barrier.** Each of the five sits at a phase boundary,
between two loops rather than inside one, so there is nothing across it the optimizer should have
been moving; but that is an argument rather than a measurement. The measurement is the `AMDt ms`
column, the hooked copy's own wall time, printed beside the unhooked `AMD ms`. **If the two
disagree by more than this benchmark's plus or minus 3 percent floor, the instrument is distorting
what it measures and the phases are not to be believed.** The `other` column is the second
self-check: the call's wall time less the five phases, which is the clock reads themselves and
should read zero.

### A first reading, and it is NOT a result

**Linux sandbox, GCC, 2026-08-15. Directional only**, on the standing rule that this machine is not
a measurement platform:

```
grid          AMD ms   AMDt ms    valid      aat    build     core     post    other    comp%
32x32          0.305     0.316    0.007    0.009    0.010    0.258    0.031    0.000    85.1%
64x64          1.040     1.044    0.026    0.034    0.042    0.850    0.091    0.000    85.5%
100x100        2.272     2.261    0.064    0.083    0.098    1.834    0.182    0.000    85.5%
```

Both self-checks pass: `AMDt` tracks `AMD` within the floor, and `other` is zero to three decimals.
**And the share reads FLAT**, 85.1, 85.5, 85.5 over a tenfold range in n, which if it holds on
alpamayo means the excluded phases are a constant fraction and **the amd branch's 2D growth is
ours**. It also puts the correction nearer 15 percent than the 8 recorded here in August, which
would move every amd ratio in this file by more than the noise floor.

Three sizes on the wrong machine settle none of that. What the alpamayo run has to answer, in
order: do both self-checks pass; is `comp%` flat across the full 32-to-400 ladder; and is it the
same on cubes, where the phase mix should differ, `AMD_aat` being linear in nnz(A) while `core`
grows with fill.

## MMD3's remaining gap is a CONSTANT, not growth, 2026-08-14

`make scale2d` and `make scale3d` on alpamayo, the first full pair since the 2026-08-10 fusions.
Ratios against the vendored routine of the same lineage:

```
            MMD1    MMD3            AMD1    AMD3
32x32       1.53    1.17            1.13    1.25
64x64       2.29    1.43            1.22    1.60
100x100     2.43    1.39            1.50    1.83
140x140     2.33    1.35            1.44    1.73
200x200     2.52    1.37            1.54    2.01
280x280     2.60    1.39            1.63    1.86
400x400     2.72    1.42            1.66    2.03

6^3         1.49    1.19            1.16    1.25
12^3        3.05    1.02            1.33    1.14
16^3        2.72    1.02            1.20    1.09
20^3        2.94    0.94            1.20    1.11
26^3        3.93    0.88            1.38    1.20
32^3        4.80    0.89            1.37    1.17
```

**MMD3 in 2D is flat at 1.35 to 1.43 from 64 a side to 400**, a forty-fold range in n. The 1.17 at
32 is 0.06 ms against 0.07 and sits under this benchmark's noise floor, so the apparent rise from
it is not a trend. That settles the shape of the question: a constant factor, not a term that
scales, so the profile is looking for something we do per vertex or per bucket operation that
genmmd does more cheaply, and not for an extra pass or a worse complexity. In absolute terms at
400 a side it is 12.87 ms against 9.07, so 3.8 ms spread over 160000 vertices.

**In 3D it crosses below 1 at 20 a side** and holds there, 0.94, 0.88, 0.89. That is consistent
with the clique-size account above: as cliques grow, the per-clique work the two-source split does
faster comes to dominate the per-vertex overhead we carry, and the sign of the difference flips.
The same account then says the 2D residual is whatever we pay per vertex or per bucket operation,
since that is the term a smaller clique cannot amortize.

**MMD2 and MMD3 are the same speed everywhere**, 0.83 against 0.77, 3.24 against 3.26, 13.73
against 12.87. The four tie-break reversals cost nothing measurable, so a profile of either
answers for both, and the question is about the extras rather than about mmd3.

**MMD1 at 2.3 to 2.7x against MMD3's 1.4x says the extras are already paying for themselves in
2D.** The residual is not the extras being slow; it is what remains after they have done their
work.

**AMD3 on cubes is 1.09 to 1.25**, against 1.44x recorded after the hash fix and 3.0x before it,
so the two fusions took another quarter out of it. 2D is 1.60 to 2.03 above 64 a side, still the
larger of the two, and still the growth term that "What both families look like now" identified
and nothing has yet addressed.

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

**A scratch driver left in `src/` is linked into every benchmark whether or not anything calls
it**, because `benchmarks/ordering/Makefile` builds its source list with
`$(wildcard ../../src/*.cpp)`. The top-level `Makefile` and `CMakeLists.txt` have explicit lists,
and removing a name from those changes nothing here: a variant is gone when its FILE is gone. One
more reason to keep a vehicle short-lived, found on 2026-08-10 after `Amd3B` kept appearing in the
compile line long after it had been taken out of both build files.

**The fill columns are identical to a Linux run of the same benchmark, digit for digit**, under a
different compiler and standard library. That is what an ordering should be, there being no
floating point anywhere in it, but it is now observed rather than assumed, and it makes fill a
usable regression signal across machines in a way timings never are.
