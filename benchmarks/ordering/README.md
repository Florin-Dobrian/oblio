# Ordering Benchmark

What each ordering method costs: wall time to produce the permutation, and nnz(L) under it. Six
methods, two lineages: the vendored MMD and AMD against Oblio's own MMD1, MMD2, AMD1 and AMD2.

**A benchmark, not an experiment.** The studies under `experiments/` are frozen: each answers one
question with a measurement and is not maintained afterwards, so nothing may depend on their
contents staying current. This is the opposite. It links `../../src` directly and is expected to
keep compiling as the tree moves, which is why `make` builds it without running it: a benchmark
that silently stopped compiling would be discovered on the day it was wanted.

```
make        build
make run    build and run, grid sides 32, 64, 100, 140
make clean

./order_timing_cpp 200 280      any grid sides
```

## What it measures, and how

Grid Laplacians, since fill on a 2D grid is the case every minimum-degree paper reports and the
one where the methods are known to separate. Ordering time only, not analysis or factorization:
`OrderEngine::compute` and nothing else. Best of three after a warm-up run, because a single cold
reading is a reading rather than a result, and `-O3 -DNDEBUG` for the same reason.

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
