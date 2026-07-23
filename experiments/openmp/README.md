# OpenMP Experiment

How much parallelism does Oblio already have without asking, and how much is left for
OpenMP to take? Reference / teaching only, **not** part of the main Oblio build.

`docs/ARCHITECTURE.md` says node parallelism is present and free through Accelerate
while tree parallelism is absent and would be Oblio's to build, and
`docs/DESIGN_DECISIONS.md` (2026-07-22) recommends OpenMP for the second. Both
statements are reasoned from documentation, and nothing in the tree has ever run a
thread. This experiment puts numbers on the first, because until the free part is
measured the second cannot be priced.

## Terminology

Two names for the same two things, and this file uses both deliberately.

**Data parallelism** is the same operation applied across many elements: one loop, split
by index range. It scales with the core count, since a call can be cut as many ways as
there are cores. OpenMP spells it `parallel for`.

**Task parallelism** is independent units of work running at once. It is capped by how
many units exist, so two independent calls have a ceiling of 2x however many cores are
free. OpenMP spells it `task`.

These are the general names, and the experiment uses them for the mechanism itself. The
sparse-direct literature calls the same two things **node-level** and **tree-level**
parallelism, and MUMPS's three node types are two node-level and one tree-level. That
vocabulary is used here wherever the subject is Oblio, a front, or the forest, since it
is what `ARCHITECTURE`, `TODO`, and the MUMPS papers use.

So: data parallelism is node-level and is what Accelerate already does inside one dense
kernel call. Task parallelism is tree-level and is factoring independent forest branches
at once. The mapping is exact; only the audience differs.

## The question that comes first

The tempting first exercise is to thread something and report a speedup. That would be
premature here, because a vendor BLAS already threads a large call on its own. If
Accelerate is using every core inside one `gemm`, then running two `gemm`s side by side
buys nothing, and a speedup near one would say nothing about OpenMP.

So the primary measurement is not a speedup at all. It is a **crossover**: one product,
timed across a range of orders, with the library's own thread cap on and off. Where the
two runs agree, the BLAS decided the call was too small to be worth threading, the
cores are sitting idle, and only we can reach them. Where they diverge, the BLAS is
already using the machine and OpenMP would be redundant there.

That crossover is a direct design input. A future scheduler has to decide, per region
of the forest, between one front on many threads and many fronts on one thread each,
and the crossing point is where that decision flips. Nothing in the tree knows it
today.

It also bounds the whole idea honestly. If Accelerate turns out to thread everything
down to order 32, then tree parallelism has almost nothing to pick up and the `TODO`
item should be downgraded on the evidence. That would be a good outcome for an
afternoon's work.

## The question that comes second

Below the crossover, the follow-up is worth asking: what do two whole products run side
by side buy over the same two run back to back? Two independent products are two
disjoint subtrees in miniature. There is no dependency to express and no data shared
between them, the easiest case task parallelism ever presents, and that is deliberate.
If the easy case does not scale, the hard one will not either.

Note the shape. Splitting one product across cores would be **node** parallelism, the
part we already have and did not write. Running two whole products at once is **tree**
parallelism, the part we would have to build. The experiment measures the second and
uses the first as context.

The construct is `parallel` / `single` / `task` rather than the simpler `sections`,
even though with exactly two independent units `sections` would read more cleanly. The
task skeleton is the one a postorder walk of the forest would use, and `sections`
generalizes to nothing.

## Two kernels, because one number cannot be read

The same products are computed twice over, by a hand-written triple loop and by BLAS.
The hand loop looks like a slow also-ran and is not one. It is the **control**.

It uses only per-core resources: registers, cache, and whatever SIMD the compiler
finds. Nothing inside it threads. So two of them on two cores should scale close to 2x,
and if they do not, OpenMP never gave us two threads and no other row on the page can
be interpreted. With the control in place, a shortfall on the BLAS row has one
remaining explanation.

That explanation, on Apple Silicon, is the AMX coprocessor. It is shared, roughly one
block per CPU cluster rather than one per core, so two calls that both lean on it queue
for the same hardware. `ARCHITECTURE` predicts this and predicts it from documentation.

The hand loop is deliberately untuned: a plain column-major triple loop, no blocking,
no unrolling, no intrinsics. Tuning it would make it a competitor, and it is not
competing.

## Why the size sweep is the experiment

Order runs from 32 to 512, and the sweep carries the argument. A small front is a leaf:
the BLAS barely threads it and barely reaches for AMX. A large front is a root: one call
may already saturate the machine. That is the leaf-versus-root inversion `ARCHITECTURE`
describes, and a curve is what turns it from a claim into a crossing point.

The number of products per task scales so each order does roughly the same arithmetic.
At order 32 that is thousands of tiny products per task, which is precisely the leaf
swarm; at order 512 it is one product, which is a root front.

## Reading the output

Each run prints, per order, the time for one call and its rate in GF/s, then the two
calls back to back, the two calls in parallel, the ratio, and the skew.

**The GF/s column is the primary result, and it is read across runs rather than down
one.** `make test` invokes the same binary capped and uncapped; comparing that column
between the two gives the data-parallelism factor at each order, and where the factor
falls to one is the crossover.

The `tree` column is the secondary result, read down. It is what OpenMP bought.

The `skew` column is a confound detector, not a result. It is the slower of the two
tasks over the faster, measured inside the parallel region. A value near one means the
pair ran at the same speed and the `tree` column can be trusted. A value well above one
means they did not, and the likely cause is core placement: this is a heterogeneous
machine, and a pair landing on one performance and one efficiency core has its wall
time set by the slower of the two, which is a scheduling artifact rather than
contention. `OMP_PROC_BIND` is ignored on Darwin, so there is no way to pin threads,
which is why this is measured rather than prevented.

Every configuration is warmed up once and timed three times at every order, and the
shortest wall time is kept. An earlier version tapered the trial count with order to
keep the sweep short; the large rows came out incoherent, with a single unwarmed call
absorbing cold-page first touch and reporting one product as slower than two. Cheap
measurements are worse than no measurements, because they look like results.

Both kernels sweep the same orders, which is what lets the two tables be read side by
side. An earlier version ran the BLAS to 1024 and stopped the hand loop at 512, since the
hand loop costs over a hundred milliseconds a call at that size; a row present in one
table and absent from the other turned out to be more confusing than a row missing from
both, so the sweep now stops at 512 for both.

## The BLAS cap

`docs/TODO.md` asks that any parallel forest region cap the BLAS to one thread per
front, since a self-threading BLAS otherwise oversubscribes the cores. Here that is a
measurement rather than a warning. The cap is an environment variable set from outside
the process (`VECLIB_MAXIMUM_THREADS` for Accelerate, `OPENBLAS_NUM_THREADS` for
OpenBLAS), so nothing in the code has to guess at a threading API.

**Worth confirming on the target machine rather than trusting here.** Accelerate's
current documentation also mentions a `BLAS_THREADING` control, and which knob it
actually honors is a question for alpamayo. The binary prints the variables it sees, so
a cap that is not honored shows up as a capped GF/s column identical to the uncapped
one, which is a finding about the knob rather than about the parallelism.

## What the serial build is for

`make test` builds the same source twice, with OpenMP and without, and runs both. The
second build is not a fallback. It is the claim in `DESIGN_DECISIONS` that OpenMP
degrades to valid serial code made checkable.

Two findings from building it, both small and both worth having in writing.

**The pragmas degrade; the runtime API does not.** Without `-fopenmp` a `#pragma omp`
is ignored, exactly as the design entry says. A call to `omp_get_max_threads()` is not:
it fails to *link*. So the property holds only if every runtime call sits behind
`#ifdef _OPENMP`, and the design entry is worth a sentence saying so, since it
currently reads as though the whole feature vanishes cleanly.

**A serial build is not warning-clean under `-Wall`.** Every ignored pragma is reported
through `-Wunknown-pragmas`, four of them in one small file. The tree is otherwise
`-Wall -Wextra` clean, so the serial target passes `-Wno-unknown-pragmas`.

## Build

```
make          build everything
make test     the measurement: three configurations, both binaries
make example  the teaching example, with and without OpenMP
make clean
```

`-O3 -DNDEBUG`, matching the other experiments, and for the same reason: the claim is
about what parallelism costs and buys, and at `-O0` that would measure nothing.

**On macOS, OpenMP needs libomp.** This is the one asymmetry with the BLAS, which is a
system framework needing nothing but `-framework Accelerate`. Apple's Clang ships no
OpenMP runtime, so libomp comes from Homebrew and the Makefile adds its include and
library paths. Set `OMP_PREFIX` if Homebrew is not at `/opt/homebrew`. Every other
platform has OpenMP in the compiler, and there the asymmetry inverts: OpenMP is built
in and the BLAS is the thing to install.

## OpenMP basics

Enough of the model to read the code in this folder, and to know why the pragmas sit
where they do.

### The execution model is SPMD

Single Program, Multiple Data. `#pragma omp parallel` starts a team of N threads and
**every one of them executes the block**. It does not create a coordinator plus workers,
and it does not run the block once. Nothing in the team is distinguished.

That single fact explains the shape of everything else. If all N threads run the block,
then anything that should happen once, rather than N times, needs saying so explicitly.

### The four directives used here

```
parallel     start a team               N threads all run the block
for          divide these iterations    distributes work that already exists
single       one thread runs this       the rest skip to the barrier at its end
task         here is a work item        for whichever thread is free
```

`parallel for` is a combined construct for the common case: `parallel` alone would have
every thread run the whole loop, and the `for` is what splits the iterations between
them. Without the `for`, ten threads each run all n iterations.

### Why `single` is not redundant

The same logic applies to task creation, and this is the part that reads as redundant
until it is measured. Ten threads, one `task` directive, counting the tasks actually
created:

```
without single    10 tasks
with single        1 task
```

So in the two-call form used here, omitting `single` gives twenty tasks rather than two,
because each of the ten threads creates both. The results stay *correct*, since each task
recomputes `c[i] = a[i] + b[i]` from unchanged inputs and ten identical writes are
harmless, so the failure mode is not a wrong answer but ten times the work. That is the
expensive kind of bug: it shows up as a missing speedup with everything passing.

Over a real tree it would be worse than wasteful. Every thread would spawn its own copy
of the subtree walk, and work that is not idempotent, an extend-add accumulating into a
parent, would be applied repeatedly. Wrong answers, nondeterministically.

The division of labor: `parallel` creates the threads, `single` picks one to **create**
the work, `task` describes a unit **any** thread may execute. The nine threads that
skipped the `single` block are not idle. They are waiting at its barrier, and a barrier
is a task scheduling point, so that is exactly where they pick tasks up.

### Data sharing defaults

Worth knowing because the defaults are mostly right and the exception is not obvious.

Variables declared *outside* a `parallel` region are shared. Variables declared *inside*
are private to each thread. For a `task`, a variable that is shared in the enclosing
context stays shared, while one that is private there is captured `firstprivate`, that
is, by value at task creation.

The trap is the shared default combined with a loop variable or a scratch buffer hoisted
out of a loop. Hoisting is good serial practice, since reuse is free and allocation is
not, and it is exactly what makes a loop unparallelizable. Oblio's `gblToLcl` is this
pattern.

### The patterns worth knowing

**Independent loop iterations.** One line, nothing else changes:

```c
// DATA parallelism: split ONE call across cores. The pragma goes inside.
// `for` is what divides the iterations; `parallel` alone would have every
// thread run all n of them.
#pragma omp parallel for
    for (size_t i = 0; i < n; ++i) c[i] = a[i] + b[i];
```

**Independent units of work.** The canonical tasking idiom, and what this folder uses:

```c
// TASK parallelism: run TWO whole calls at once. The callee is not modified at all.
#pragma omp parallel              // start a team; all N threads run this block
    {
#pragma omp single                // one thread creates the tasks; without this,
        {                         // all N would, giving 2N tasks instead of 2
#pragma omp task                  // a work item, for whichever thread is free
            doThing(a);
#pragma omp task
            doThing(b);
        }
    }                             // implicit barrier: both tasks done past here
```

**A tree or DAG.** The recursive form, which is what a postorder walk of an elimination
forest would use. `taskwait` waits for the tasks this task created, and nothing else:

```c
// TASK parallelism over a dependency graph. The recursion expresses the
// dependency, so no `depend` clauses are needed.
void walk(int node) {
    for (int child = firstChild[node]; child != NIL; child = nextSibling[child])
#pragma omp task                  // each subtree is an independent unit
        walk(child);
#pragma omp taskwait              // wait for THIS task's children, nobody else's
    process(node);                // runs only after every child is done
}
```

Called once from inside a `parallel` / `single` block. This is the form that needs no
`depend` clauses, because the recursion already expresses the dependency.

**Accumulating into one variable.** Needs declaring, and changes the answer:

```c
// Not independent: every iteration writes the same `sum`. The reduction clause
// gives each thread a private copy and combines them at the end.
#pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < n; ++i) sum += x[i];
```

For floating point the summation order changes, so the result moves in the last bits.
Fine for a norm, not fine where bit-reproducibility is being asserted, as it is in this
experiment's own checks.

### Variants and small print

`masked`, formerly `master`, also restricts to one thread but requires thread 0 and has
no implicit barrier at the end. For task creation `single` is usually what is wanted,
since that barrier is where the other threads go to collect work.

`nowait` removes the implicit barrier from a `for` or `single` when the following code
does not yet depend on the result.

`task depend(in: x) depend(out: y)` expresses a dependency between sibling tasks without
recursion. It is the alternative to the recursive form above, and the one to reach for
when the work is already a flat loop that cannot easily become a tree walk.

Two further pieces of small print, about building with OpenMP switched off, are covered
above under "What the serial build is for": the runtime API needs an `#ifdef _OPENMP`
guard where the pragmas need nothing, and a serial build is not warning-clean without
`-Wno-unknown-pragmas`. `maxThreads()` in `example_parallel.cpp` demonstrates the guard.

## The example: how it is actually done

`example_parallel.cpp` is the teaching half of this folder, separate from the measurement
and built separately with `make example`. It needs no BLAS, only a compiler, because it
is about the pragmas rather than about this machine.

One function adds two lists elementwise. There are two pairs of lists to get through, and
that is what the two styles divide up differently.

The serial version, and the baseline both parallel forms must reproduce exactly:

```c
void sum(const double* a, const double* b, double* c, size_t n)
{
    for (size_t i = 0; i < n; ++i) c[i] = a[i] + b[i];
}
```

**Data parallelism** puts the pragma inside. Each call is split across the whole team,
and the two calls still happen one after the other:

```c
// The pragma goes INSIDE. One call, divided among N threads.
void sumDataParallel(const double* a, const double* b, double* c, size_t n)
{
#pragma omp parallel for          // divide the iterations across the team
    for (size_t i = 0; i < n; ++i) c[i] = a[i] + b[i];
}

// The call site: two pairs of lists, each pair split across every core,
// the second pair starting only once the first is finished.
sumDataParallel(a1, b1, c1, n);
sumDataParallel(a2, b2, c2, n);
```

**Task parallelism** puts the pragmas outside. The two calls run at once, and the
function they call is the plain `sum` above, with no pragma in it:

```c
// The pragmas go OUTSIDE. Two whole calls at once, callee untouched.
void sumTaskParallel(const double* a1, const double* b1, double* c1,
                        const double* a2, const double* b2, double* c2, size_t n)
{
#pragma omp parallel              // start a team; all N threads run this block
    {
#pragma omp single                // one thread creates the tasks; without this,
        {                         // all N would, giving 2N tasks instead of 2
#pragma omp task                  // a work item, for whichever thread is free
            sum(a1, b1, c1, n);   // plain serial sum: no pragma inside it
#pragma omp task
            sum(a2, b2, c2, n);
        }
    }                             // implicit barrier: both tasks done past here
}
```

**Note where the name sits.** `sumDataParallel` is a parallel version *of* `sum`, and
could be swapped in for it anywhere. `sumTaskParallel` is not: it is a caller, and the
thing it calls is `sum` itself, unmodified. There is no task-parallel version of the
kernel, because task parallelism does not need one. That asymmetry is the distinction
itself. In the real file `sumTaskParallel` also takes the kernel as a parameter, so the
same driver serves both the cheap and expensive flavors; that argument is dropped here
for readability.

That is the entire change in both cases. No threads created by hand, no join, no
restructuring, and both compile back to the serial code without `-fopenmp`. For a loop
whose iterations are independent, the claim that OpenMP is one pragma is true.

The ceilings differ, and this is where the two styles stop being interchangeable. Data
parallelism scales with the core count, since one call can be cut ten ways. Task
parallelism is capped by how many independent units exist: two calls means 2x at best,
however many cores are idle.

In Oblio's terms, the first form is what Accelerate already does inside one dense kernel
call, and the second is factoring two independent forest branches at once.

### Why the same kernel runs twice, cheap and expensive

The example runs both parallelizations on two flavors of the kernel, and the contrast is
the actual lesson.

The cheap flavor is the addition itself: one flop per 24 bytes moved, two loads and a
store. It is memory-bound in the extreme, so both cores wait on the same memory and both
parallelizations buy nearly nothing. They are still correct. There is simply nothing to
win.

The expensive flavor does a long arithmetic chain per element, hundreds of flops for the
same 24 bytes. Same pragmas, same places, and now they pay.

So the pragma is easy and whether it helps is a separate question about what the code is
waiting on. That is the same finding the `gemm` measurement produced from the other
direction, where correct pragmas bought nothing because a shared matrix unit was already
saturated.

### What breaks "just add a pragma"

Three things, worth knowing before assuming a loop is one line away.

**Shared scratch.** Serial code hoists a buffer out of a loop because reuse is free and
allocation is not. That is good serial practice and it is exactly what makes the loop
unparallelizable. Oblio's `gblToLcl` is this pattern. The fix is usually small, moving
the declaration inward or making it per-thread, but it costs the allocation the hoist
was avoiding.

**Early exits.** A `break`, `return`, or `goto` out of the loop body. OpenMP requires a
structured block, so these become a flag tested after the region.

**Reductions.** Accumulating into one variable needs `reduction(+:x)`, and for floating
point the result changes, because the summation order changes. Fine for a norm, not fine
where bit-reproducibility is being asserted, as it is in this experiment's own checks.

A fourth thing is not a correctness matter but decides everything: a dependency graph is
not a parallel loop. `parallel for` handles independent iterations. A tree or DAG needs
the code already expressed as tasks, usually through recursion, which is a restructure
rather than an annotation. That is why `docs/TODO.md` lists four code changes for
`factorMultifrontal` rather than a pragma.

## Correctness

Three checks run before any timing, and the third is the one with teeth.

The two kernels are compared against each other, which is what makes the hand loop
usable as a control rather than merely as a second timing.

Each kernel is then run serially and in parallel and the results are required to be
**bit-identical**. This is a stronger demand than the solver usually makes of itself,
and it is the right one here: each task performs the same arithmetic in the same order
as it would serially, and the two tasks touch disjoint memory, so there is no reduction
and no summation-order freedom to appeal to. Anything short of bit-identical is a data
race rather than rounding. Compare the multifrontal update stack, where reordering the
fold genuinely does move the last bits and the tests bound the residual instead.

## What friend-access already suggests, and how far to trust it

`../friend-access/` times a matvec three ways, and its BLAS row is worth rereading in
this light. It is `dgemv`, so level 2: `O(n^2)` arithmetic on `O(n^2)` data,
bandwidth-bound, with about two flops per element loaded. A front is level 3, which is
why this experiment uses `gemm` instead.

On an M4 that experiment reports the hand loop at 3.18 ms and Accelerate at 0.50 ms for
a 2000 by 2000 matrix. Read as bandwidth rather than flops, 32 MB streamed once, those
are roughly 10 GB/s and 64 GB/s. The second figure is more than one core is likely to
manage, which suggests Accelerate spread a single call across cores, and if so
friend-access has been measuring free data parallelism since before we had a name
for it.

**That is an inference and not a measurement.** Better blocking and prefetch remain
live partial explanations, and so does AMX. A capped run of that same experiment would
settle it in one line, and this one deliberately does not rest on it: the capped
against uncapped comparison here measures the same thing directly, on the kernel that
actually matters.

## Results

Measured on alpamayo (Apple M4, `omp_get_max_threads` 10) on 2026-07-23, Apple Clang
with Homebrew libomp 22.1.8, Accelerate for BLAS, `-O3 -DNDEBUG`. All figures below are
empirical, from the capped run unless stated; the mechanism follows each.

### The three runs

Time in milliseconds for the two products, run back to back and then side by side, with
the resulting speedup. This is `make test` as it stands, one table per configuration.

Serial build, no OpenMP:

```
        gemmHand                    gemmBlas
n     serial  parallel  speedup   serial  parallel  speedup
32    25.46    25.92     0.98      1.95     1.95     1.00
64    18.71    18.92     0.99      1.12     1.12     1.00
128   19.21    18.49     1.04      0.92     0.92     1.00
256   20.39    20.16     1.01      0.74     0.75     0.98
512   32.21    32.08     1.00      1.27     1.28     1.00
```

OpenMP, BLAS uncapped:

```
        gemmHand                    gemmBlas
n     serial  parallel  speedup   serial  parallel  speedup
32    26.15    20.18     1.30      1.92     1.34     1.44
64    18.76     9.81     1.91      1.12     0.98     1.15
128   18.91     9.92     1.91      0.92     0.91     1.01
256   19.60    11.28     1.74      0.74     0.76     0.97
512   32.57    23.18     1.41      1.26     1.26     1.00
```

OpenMP, BLAS capped:

```
        gemmHand                    gemmBlas
n     serial  parallel  speedup   serial  parallel  speedup
32    26.44    20.11     1.31      2.04     1.36     1.50
64    18.86     9.81     1.92      1.12     0.99     1.13
128   18.93     9.86     1.92      0.92     0.90     1.03
256   20.30    11.07     1.83      0.74     0.75     0.98
512   32.31    25.05     1.29      1.26     1.28     0.99
```

### The control validates the harness

Hand loop, capped, speedup column: 1.31, 1.92, 1.92, 1.83, 1.29 across orders 32 to 512.
Two independent tasks on two cores return a genuine 1.9x through the middle of the
range, so OpenMP does what it claims and the BLAS rows can be read as measurements of
the kernel rather than of the scheduler. The serial build returns 1.00 throughout, which
is the degradation check passing.

Two rows at the ends fall short of 1.9x, for different reasons.

Order 512 gives 1.29 capped at `skew` 1.34, against 1.41 uncapped at `skew` 1.01, so the
capped figure is partly core placement and 1.41 is the better estimate. Some real
shortfall remains: each task holds about 6 MB and the loop does no cache blocking, so at
that size the two contend for memory bandwidth.

Order 32 gives about 1.30 in both OpenMP runs with `skew` near 1.00, so both tasks were
slowed equally, by roughly 1.5x each. Equal slowing with no imbalance is not core
placement. The cause is not established here.

### The payoff window closes early

Both kernels from the capped run, identical scheduling, only the kernel differing:

```
n      hand    BLAS
32     1.31    1.50
64     1.92    1.13
128    1.92    1.03
256    1.83    0.98
512    1.29    0.99
```

At order 32 the BLAS scales at least as well as the control, and in this run better. By
order 128 it has stopped entirely, and stays stopped. So tree parallelism over dense
fronts pays only on very small fronts, and the crossing point sits between 64 and 128.

`ARCHITECTURE` predicted that the leaves would pay and the root would not, and that is
confirmed. What it did not predict is how early "leaf" ends: order 64, not order 500.

### Mechanism: unresolved, but the disjunction is tight

The observation to explain is narrow. Two products that are independent, that touch
disjoint memory, and that are handed to two OpenMP tasks do not finish sooner than the
same two run one after the other, above order 64. The hand loop under identical
scheduling finishes in half the time. Same tasks, same runtime, same machine, different
kernel.

Two independent tasks on two idle cores must beat two in sequence. Since they do not,
at least one of the following holds.

**The baseline was never serial.** From OpenMP's point of view the two products run back
to back, but if Accelerate spreads each call across cores then that baseline was already
parallel by other means. There would be no idle cores for OpenMP to find, and no speedup
to win, because the work was already spread before we asked.

**Or the cores are free and the bottleneck is elsewhere.** The threads may be available
but all funnel into a matrix unit that is shared per cluster rather than duplicated per
core, and that one call already keeps fed. Adding cores then adds feeders to a full
queue.

These are not exclusive and both may hold at once.

What the evidence does support is that a shared matrix unit is involved somewhere. Peak
is about 448 GF/s in double precision against the hand loop's 21 on one core. Even
granting a well-tuned NEON kernel several times the hand loop's rate, ten cores of
ordinary SIMD do not reach 448 in double precision. So the BLAS is reaching hardware the
hand loop cannot. That establishes the unit exists, not that Accelerate declines to
thread.

Either way the small-size behavior follows: at order 32 there is headroom, whether
because a small call is not worth threading or because it cannot saturate the shared
unit, and OpenMP finds it. That is why the window exists and why it closes.

Whether the unit is the older AMX path or the M4's SME is not determinable from here,
and does not affect the argument. What matters is that it sits above the core.

### The cap changes nothing measurable

GF/s agrees across all three runs at every order to within noise, and
`VECLIB_MAXIMUM_THREADS=1` makes no difference anywhere.

Two readings remain open and this experiment does not separate them. The variable may be
dead on current macOS, with Accelerate threading regardless. Or Accelerate may not thread
`dgemm` at these sizes, leaving the cap nothing to act on.

**A second knob fails the same way, and that does narrow things.** Under
`OMP_NUM_THREADS=1` the single-call rates are 204.85, 354.43, 423.67, 446.35, 436.07,
439.04 across orders 32 to 1024, matching the unconstrained run. So `OMP_NUM_THREADS`
does not reach Accelerate either. That rules out one mechanism: if the library threaded
*through OpenMP*, capping OpenMP would have capped it. Either it does not thread, or it
threads through something honoring neither variable, and Grand Central Dispatch is the
obvious candidate on this platform.

The practical consequence is worse than "the cap is unnecessary". If Accelerate does
thread through a mechanism we cannot reach, then when Oblio eventually runs fronts as
concurrent tasks, each front's BLAS call could spawn its own work with no lever to stop
it. The `TODO` instruction to cap the BLAS per front would then be unimplementable rather
than merely redundant.

One lever remains untried: Accelerate's own threading control, `BLASSetThreading` with a
single-threaded mode, called from code rather than set in the environment. Its exact
spelling and availability want checking against Apple's documentation rather than
asserting from memory. If it exists and works, it settles both questions at once, since
the GF/s column moving at all would prove the threading was there.

One thing the data does locate is the ceiling. Two concurrent processes, sharing no
OpenMP runtime and no address space, so beyond the reach of both our scheduling and our
cap, cannot exceed one process's aggregate throughput above order 128. Whatever limits
us is in the hardware or the library, not in our code.

The practical conclusion holds under either reading: the `TODO` instruction to cap
Accelerate to one thread per front changes nothing on this hardware. Worth knowing before
that code gets written, and worth re-checking on any machine where the GF/s column does
move with the cap.

### What would settle the mechanism

The obvious test is to serialize our own tasks and ask how many cores the process uses,
leaving only the library's parallelism:

```
for i in 1 2 3; do OMP_NUM_THREADS=1 /usr/bin/time ./test_parallel_cpp > /dev/null; done
```

**This test as written does not work, and the reason is worth stating.** It returns a
user-to-real ratio near 1.06, which looks like a single-threaded process. But the BLAS
sweep is only about 60 ms of a 1.34 s run, since the hand loop dominates by more than an
order of magnitude, so the overall ratio is nearly insensitive to what Accelerate does.
Attributing the excess CPU to the BLAS portion alone gives something closer to two cores
active while the BLAS runs, and the measurement is too coarse to say more.

Making it work requires timing the BLAS sweep in isolation, with the `gemmHand` sweep
removed, so the ratio is entirely about Accelerate. That is a small change and has not
been made.

The second test, which does work, asks whether the ceiling is ours or the machine's by
running two processes that share nothing:

```
OMP_NUM_THREADS=1 ./test_parallel_cpp > run-a.txt & \
OMP_NUM_THREADS=1 ./test_parallel_cpp > run-b.txt & wait
```

Compare `gemmBlas` timings against a single run. If two processes together do no more
work than one, the limit is in hardware or the library rather than in our scheduling.
Read the `two serial` column rather than `one call`, since it is longer and less
sensitive to the two processes drifting out of overlap.

### The skew column caught a false result immediately

Measured while the BLAS sweep still ran to 1024, since removed. Capped, order 1024,
`skew` 2.23: one task took more than twice as long as the other, so that row's 0.82x is
core placement rather than contention and had to be discounted. The uncapped 1024 row
showed skew 1.10 and 0.94x, also discountable. The likely cause is the
efficiency cluster having its own slower AMX block, which would produce exactly this
signature when a pair straddles the two clusters.

### Corrections worth recording

Four, and the last is the general one.

**The tapered trial count.** The first run used no warm-up at the large orders and
produced rows where one product timed slower than two. Read at face value it suggested
contention at every size, and that reading was wrong: the window at order 32 to 64 is
real and the bad data hid it. The tapering is gone and `kTrials` now applies everywhere.

**A single timing outlier that reversed a conclusion.** One `/usr/bin/time` run reported
11.63 real / 14.75 user where the stable figure is 1.34 / 1.42, both numbers inflated
roughly tenfold. A user-to-real ratio of 10.3 looks exactly like ten cores at work, and
on the strength of that one sample the coprocessor explanation was withdrawn and
replaced with "Accelerate threads across all cores and ignores the cap". Three repeats
showed the outlier for what it was. The cause is still unidentified: both real and user
inflated, so idle threads spinning does not fit, and efficiency-core placement or
leftover load from a concurrent test fits better. Not chased further.

**A measurement too coarse for the claim it was used for.** The `/usr/bin/time` ratio of
1.06 was read as proof that Accelerate does not thread `dgemm`, and the mechanism
subsection was written around that. The BLAS sweep is about 5% of that run, so the ratio
could not have shown otherwise whatever Accelerate did. The claim was retracted and the
mechanism is now recorded as unresolved. Stability is not accuracy: three consistent
readings of the wrong quantity are still the wrong quantity.

**The general lesson.** The stated mechanism was reversed twice in one session, each time
on one measurement, in both directions. Both reversals would have been avoided by the
rule the harness already applies to its own rows and the reader did not apply to
`/usr/bin/time`: warm up, repeat three times, and only then believe it. A measurement
taken once is a reading, not a result. And a measurement repeated three times is still
only a reading of whatever it actually measures.

### What this means for Oblio

Tree parallelism over dense fronts is worth roughly 1.5x on fronts of order 64 or
smaller and nothing above that. Whether the `TODO` item is worth building therefore
reduces entirely to the question below, which this experiment cannot answer.

## What this cannot answer

Whether the leaf swarm is a large enough share of Oblio's real runtime to be worth
attacking. Small fronts are a small share of the *flops* and a larger share of the
*time*, because they run at poor efficiency, and the ratio depends on the matrix. That
question needs the actual solver instrumented by front size, which is a separate and
later measurement. It should not be conflated with this one.

## Homebrew maintenance

**Lodging here temporarily.** This is general machine upkeep and has nothing to do with
OpenMP beyond the fact that libomp arrives this way. It is recorded here because this
is where the question came up, and it should move to `CLAUDE.md` under Tooling, or to
`CONTRIBUTING.md`, once there is a natural home for it.

```
brew update            # refresh the formula definitions
brew outdated          # read this before committing to anything
brew upgrade           # formulae
brew upgrade --cask    # casks, separately and after
brew cleanup           # reclaim the old versions
```

Reading `brew outdated` before upgrading is the step worth not skipping. It is the only
chance to see what is about to move while `brew pin` is still an option, and major
version bumps hide in that list without announcing themselves.

One machine at a time, with the tree's `make test` run in between. That keeps a
known-good machine as a fallback, and 153 passing assertions is a sharp enough check to
catch a broken toolchain immediately. The argument is the same one that keeps the 0.9
oracle around: a reference is only useful if it has not moved.

Recent Homebrew versions upgrade casks as part of `brew upgrade`, so the separate
`--cask` pass often finds nothing left to do. It costs nothing to run and is worth
keeping in the sequence for the versions that do not.

Avoid pairing this with `brew autoremove`. A formula installed as a dependency, whose
dependents have since been removed, is what autoremove exists to clean up, and libomp
has been in exactly that position on at least one machine here. Installing it
explicitly marks it as on-request and takes it out of reach. If autoremove is wanted,
`--dry-run` first and read the list.

Afterwards, record what the machines are running:

```
brew list --versions libomp suite-sparse gcc cmake
```

Matching output across machines is the thing worth having later, when a numeric result
disagrees between them and the first question is whether the environments are the same.

Two notes specific to this tree. None of Oblio's build depends on Homebrew: the
compiler is Apple's Clang and the BLAS is Accelerate, both system components. And the
`suite-sparse` formula is unrelated to the vendored AMD 3.3.4 compiled into `src/`,
which nothing links against Homebrew's copy.

## Related

- `docs/ARCHITECTURE.md`, the traversal-choice section, for node against tree
  parallelism, the MUMPS node types, and the AMX inversion this measures.
- `docs/DESIGN_DECISIONS.md` (2026-07-22), for why OpenMP first and what the
  alternatives cost.
- `docs/TODO.md`, "Forest (tree) parallelism on Apple Silicon", the work item this
  experiment is reconnaissance for.
- `../friend-access/`, for the level-2 BLAS handoff and the numbers reread above.
- J.-Y. L'Excellent and W. M. Sid-Lakhdar, "A study of shared-memory parallelism in a
  multifrontal solver", *Parallel Computing*, 2014.
  `https://sciencedirect.com/science/article/pii/S0167819114000246`. The reference for
  the two-rung arrangement this experiment keeps running into: a node has two sources of
  threading, the multithreaded BLAS and the solver's own OpenMP, and the two are measured
  separately rather than lumped. It also proposes a performance model to decide, per
  region of the tree, between tree parallelism and multithreaded libraries, which is the
  same crossover measured here. Its finding that the payoff depends on the ratio of large
  fronts to small ones is the front-size question this experiment cannot answer.
  The longer earlier version is INRIA Research Report RR-8227, 2013, "Introduction of
  shared-memory parallelism in a distributed-memory multifrontal solver".
