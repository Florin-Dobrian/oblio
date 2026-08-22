# Benchmarks

> **SUPERSEDED IN PART, 2026-08-21.** `MMD1`, `MMD2`, `AMD1` and `AMD2` were retired to
> `retired/` and are out of the build, so their rows here are a record of drivers that no longer
> exist in production. The ladder is intact as prototypes in `experiments/ordering/`; see
> `retired/README.md`.

Timing and profiling work, one self-contained folder per subject.

- **`ordering/`** measures one phase against itself: what each ordering method costs and how much
  it fills.
- **`pipeline/`** measures the phases against each other: what share of a solve the ordering is,
  and after how many factorizations a slower-analyzing ordering pays for itself. It also measures
  how the cost GROWS, over a longer ladder and every factorization Oblio has, which is
  `SCALING.md` there.

The second exists because the first cannot answer whether its own subject matters, and the answer
was not the expected one. Ordering is a substantial share of a one-shot solve rather than a
fraction of a percent, so it is worth optimizing: between a tenth and a fifth of a one-shot solve,
17.8 percent at 140 a side in 2D and 10.3 percent at 26 cubed, falling as the factorization
grows.

**Both of those shares are family dependent**, and so is the comparison between the two branches.
On square grids MMD wins on all three terms, ordering no slower, filling up to 13 percent less and
factoring 5 to 7 percent faster. On cubic grids exactly one of the three reverses and it is the
smallest: AMD orders about twice as fast, fill is within a percent either way, and MMD still
factors 4 to 10 percent faster at the same nnz(L). So AMD wins a one-shot solve on cubes and MMD
wins anything factored more than about four times, where in 2D MMD wins at every count.

The old text here said fill differences among good orderings do not propagate into factorization
time at all. That was a square-grid observation stated generally, and the cubic tables show
something more interesting than either version: **two orderings can fill the same and factor
differently**, by 4 to 10 percent, which is about supernode shape rather than fill.
`pipeline/README.md` has both, and the lesson is the one `experiments/ordering` learned first: a
claim measured on one grid family is a claim about that family.

**A benchmark is not an experiment.** The studies under `experiments/` are frozen: each answers one
question with a measurement and is not maintained afterwards. These run against the current tree
forever, so they link `../../src` directly and are expected to keep compiling as the code moves.
That is why `make` builds them without running them: a benchmark that had silently stopped
compiling would be discovered on the day it was wanted.

Each folder carries its own results with the machine, the date and the compiler beside them. A
timing without those cannot be compared against anything later.

## Two things learned the hard way, 2026-08-08

**Counting and profiling answer different questions.** Counting is right for "are we doing more
work" and blind to two things this folder cares about. It cannot see a loop whose cost is decided
by an EARLY EXIT: a counter that added `adjacencySize + incidenceSize` per pair measured what a
short-circuiting test could cost, and the real iteration count was 1.7x higher on one side, which
is where a 42-percent line was hiding. And it cannot see ALLOCATION at all: the largest single
item found in a day of this was `operator new`, an unreserved arena doubling 18 times per ordering.

**A call tree bounds the search; only the source view ends it.** A driver inlines into one symbol,
so the tree can say only how much sits in its self weight. Open the source annotation. And when the
top line is a CALL rather than a line of work, zoom again: call tree, then the call site, then
`allocate.h`. Stopping at either of the first two finds nothing.

**Fix the harness before chasing anything under 10 percent, and this one was fixed on 2026-08-08.**
`order_timing` took the best of three, which at 140 a side is nine milliseconds of measurement, and
two runs of an IDENTICAL BINARY disagreed by 4 percent. Five changes had landed inside that band
during the amd work before anyone noticed the instrument was the constraint. It now sizes its
repeat count from a timed probe targeting 300 ms per method per row, and two runs agree to a tenth
of a percent at 140 and half a percent at 32. `experiments/ordering/width.cpp` does the same and
also prints the count so a row can be judged.

**What Instruments can and cannot give you from the command line, established 2026-08-09 and worth
knowing before anyone spends an hour on it again.**

`xctrace record --template 'CPU Counters' --launch -- ./binary args` works and produces a trace. It
selects the template's DEFAULT mode, `CPU Bottlenecks`, which reports four derived percentages and
cycles: Useful, Instruction Processing Bottleneck, Instruction Delivery Bottleneck, Discarded
Bottleneck. That mix is genuinely informative: it says whether time is going to back-end stall, to
fetch, or to misspeculation: but it gives no raw event counts.

The GUI offers many other guided modes under Configuration: Guided, Mode: including `L1D Cache
Metrics`, `L1D Miss Sampling` and several instruction-characteristic counts. **Those modes are not
reachable from `xctrace`.** Only `templateCPUBottlenecks` exists as a template resource; the others
are constructed by the GUI and there is no field to pass. Setting the mode in the GUI and recording
from there does work, but the target's working directory and arguments have to be set in the
target's Options sheet or the launch fails with `Path not found`.

One lead that was not followed to the end: the default template ALREADY records eight raw per-core
PMCs in the `kdebug-counters-with-time-sample` table of the trace. They are the bottleneck event
set rather than cache events, but identifying which events those eight columns are would give raw
counts out of traces already taken. The event identities live in the PMU database under the
Instruments app bundle, keyed by a cpufamily hash.

**It costs wall time and that is the trade.** A full `make scale2d` is minutes rather than seconds
now. The alternative is a benchmark that cannot answer the questions being asked of it, which is
what it was.

## How to profile, and what each tool is actually good for

Written down because a day was spent finding this out, and two of the three tools gave answers that
were wrong or too coarse to act on.

**Instruments is the profiler of record.** Its CPU Counters template is the only one of these that
separates doing more work from doing it less efficiently, which is the distinction every question
in this folder has turned on. The others have narrower uses: timing says whether anything changed,
counting says what the algorithm does in machine-independent units, and `sample` gives a function
ranking in ten seconds with no setup. Reach for them in that order, but do not conclude from them.
gprof is on the list only to record that it was tried and misled us; it is not recommended.

The order below is the order to reach for them.

### 1. Time it first

`ordering/order_timing.cpp` is the pattern: best of three after a warm-up, `-O3 -DNDEBUG`, and
report ratios against a reference rather than milliseconds, since a loaded machine moves every
absolute number while leaving ratios intact. A single cold reading is a reading and not a result.

### 2. Count what the algorithm does

Cheap, machine-independent, and it answers "are we doing more work" rather than "are we slower".
Add counters to a scratch copy of the source, never to the source itself; the vendored routines can
be copied and instrumented the same way, which is how we compared against `Amd.cpp`.

**But a count is not a cost, and this is the trap.** Counts and time came apart three times in one
day. Pass 3 of AMD removed 3.7x the elements from a scan and bought 12 percent. Fusing three passes
over the same lists into one removed 30 percent of all element visits and bought nothing at all,
because the second and third visits hit data still in L1 from the first. Removing an allocation
saved about 6 ns; removing a dependent load in an inner loop saved far more. So counts point at
where to look and never at how much there is to gain.

### 3. Sample, for a first split

`sample` ships with macOS and needs no setup. The program must run long enough to sample, which
means a loop: `ordering/order_profile.cpp` exists for exactly that.

```
./order_profile_cpp amd1 140 5000 & sleep 1; sample $! 10 -f /tmp/order.sample; wait
sed -n '/Sort by top of stack/,/Binary Images/p' /tmp/order.sample | head -30
```

The `sleep 1` and `$!` matter: `pgrep` on a program that has not started yet returns nothing, and
`sample` then treats the next argument as a process id and fails confusingly. The `.sample` file has
no default application, so read it with `sed` or `grep` rather than `open`.

"Sort by top of stack" is the self-time ranking, which is the part worth reading. The call tree
above it also carries line numbers, so a hot function can be narrowed to a hot line:

```
grep -n "MyClass::myFunction" /tmp/order.sample | head
```

`sample` is enough to say which function dominates. It is not enough to say why, and on one
occasion it was actively misleading: gprof on Linux attributed 3 percent to allocation where
Instruments on macOS attributed 18, because the two allocators differ and gprof under-attributes
time inside library code it did not instrument. **Never conclude from a profile taken on a different
machine than the one being optimized.**

### 4. Instruments, CPU Counters, which is what actually answers "why"

The tool that ended the argument, and the one to reach for whenever a change is meant to make
something faster rather than merely different.

```
xcrun xctrace record --template 'CPU Counters' --launch -- ./order_profile_cpp amd1 140 3000
open Launch_order_profile_cpp_*.trace
```

Select the whole track, and the summary table gives cycles plus a breakdown:

```
Useful                          the fraction of cycles retiring work
Instruction Processing          stalls waiting on data, cache misses among them
Instruction Delivery            stalls waiting on instructions
Discarded                       work thrown away, branch mispredictions
```

Run it for the code under study **and for whatever it is being compared against**, through the same
binary if possible, so that the comparison is not also measuring two different programs. That is why
`order_profile.cpp` accepts the vendored methods as well as ours.

Then split the gap in two:

```
useful cycles = cycles * useful%      how much work is being done
cycles / useful cycles                how efficiently it is being done
```

This is the decomposition worth having, because the two halves call for opposite work. More useful
cycles means the algorithm does more, and no amount of layout or allocator tuning will help. Worse
efficiency means stalling, and then locality, allocation and data types are the levers.

Worked example, 3000 orderings of a 140x140 grid, alpamayo, 2026-07-31:

```
                cycles           useful    useful cycles   reading
MMD1            43.9 G           42.19%    18.5 G          3.05x cycles, 2.99x work
MMD (vendored)  14.4 G           43.02%     6.2 G          efficiency identical
AMD1            29.5 G           49.02%    14.4 G          1.95x cycles, 1.68x work
AMD (vendored)  15.1 G           56.85%     8.6 G          1.16x of the gap is stalling
```

MMD1 is doing three times the work at the same efficiency, so nothing structural will help it and
the missing mechanisms are the whole story. AMD1 is 86 percent work and 14 percent stalling. Neither
conclusion was reachable from timings, from counts, or from `sample`.

## What was tried and did not pay

Kept because a negative result stops the same idea being tried twice, and every one of these was
plausible enough to be worth an hour.

- **A flat weight array** beside the supervariable members, so the innermost loop's lookup would be
  a contiguous load rather than a scattered vector header. 2.77x against 2.76x over five runs each.
  No difference. (It came back later for a different reason: the members became a chain, so the size
  stopped being free.)
- **Reserving the incidence buffers up front**, so they would be allocated together and in vertex
  order rather than scattered across the run. No change: it moved the allocations from the run to
  construction rather than removing them.
- **Narrowing the per-vertex arrays** from `size_t` to `int32_t`, halving their cache footprint.
  About 2 percent, at the edge of noise, and it cost a departure from the convention plus casts.
- **Fusing three passes** over the same two lists into one, in AMD1. 1.44x against 1.45x. Nothing.
  Fusion removes loop setup and a re-fetch of warm data; it does not remove the per-element work,
  which is what costs.
- **Writing the adjacency through a pointer** instead of `push_back` at construction. Slower, since
  `resize` zero-fills what `reserve` does not.

The constructor turned out to be a fifth of AMD1's time and mostly first touch of freshly allocated
memory, which no loop rewriting reaches. It shrinks only by needing less memory.
