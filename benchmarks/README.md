# Benchmarks

Timing and profiling work, one self-contained folder per subject. `ordering/` is the first.

**A benchmark is not an experiment.** The studies under `experiments/` are frozen: each answers one
question with a measurement and is not maintained afterwards. These run against the current tree
forever, so they link `../../src` directly and are expected to keep compiling as the code moves.
That is why `make` builds them without running them: a benchmark that had silently stopped
compiling would be discovered on the day it was wanted.

Each folder carries its own results with the machine, the date and the compiler beside them. A
timing without those cannot be compared against anything later.

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
