# NEXT: 3D grids in the benchmark, and one last bounded attempt at the amd constant factor

A handover note between sessions, rewritten 2026-08-09 when its first priority closed. **It is
meant to be deleted** once the items below are done or abandoned. Everything in it that outlives
the task has already been written somewhere durable, and this file only points at those places:

- `docs/DESIGN_DECISIONS.md`, entry "what the vendored AMD's speed is made of", has the
  measurements, the three blocked routes and why each is blocked, and the parked proposal restated
  properly.
- `docs/DESIGN_DECISIONS.md`, entry "the ordering read freed memory for a week", has what the
  widened acceptance test found and the three method notes that came with it.
- `experiments/ordering/README.md`, section "Aligning a layer against a vendored routine", has the
  method both alignments used, and its `make amdorder` section has what the check now covers.
- `benchmarks/README.md` has what Instruments can and cannot be made to do from the command line,
  and the two rules about counting against profiling that cost a day between them.

If you find yourself about to record something here that a later reader would want, put it in one
of those instead. The previous `NEXT.md` went stale precisely because it accumulated content that
belonged elsewhere.

---

## Closed since the last note

**The alignment check is widened, and the alignment holds.** `make amdorder` ran eleven 2D grids
and now runs 38 cases over four shapes: the seven examples, 2D grids to 140, 3D grids to 24, and
nine random patterns at n = 2000. All match on alpamayo, so `Amd3` reproduces `AMD_2`'s raw
elimination order, member order included, on something far wider than the shape it was built
against.

Getting there found three things, none of them visible before:

- a defect in production `Amd3`, ledger entry 7, the stored clique degree not rewritten after mass
  elimination trimmed the clique;
- a use-after-free in the shared `QuotientGraph` that every ordering had, benign until an allocator
  recycled the block, which is why it surfaced as two machines disagreeing about integer code;
- two harness faults that looked like divergences, a dense threshold turned off by undefined
  behavior and a 3D grid builder emitting unsorted columns.

`experiments/ordering/AMD3.md`, iterations 18 to 20, is the narrative.

**And the mmd branch has an acceptance test at last.** `make mmdorder` compares production `Mmd3`
against genmmd's elimination order on the same four shapes, 38 cases, all matching. Until then the
mmd alignment rested on a scratch probe from 2026-08-07 that died with its session and on the
benchmark's fill column, which `MMD3.md` iteration 6 shows is not sufficient. It needs no hook,
genmmd emitting that order directly, which is why it is forty lines against the amd machinery.
`make aligned` runs both.

---

## Priority

**1. 3D grids in the ordering BENCHMARK**, a separate task from the check above and now the top
item. Every fill conclusion in this work is from square grids, which is one problem family and the
flattering one. The specific claim waiting on it is that our tie-break beats AMD's: `Amd2` fills
6.5 percent BELOW the vendored routine at 140 a side, while `Mmd3` matches genmmd exactly. If that
advantage is a 2D artifact we should know, because fill drives the factorization and the
factorization is where the time is. A 2x ordering costs well under 15 percent of a one-shot solve;
a fill difference costs more than that, permanently, on every factorization.

`graphs.h` now holds a 3D builder, so this is a case list rather than any new code.

**2. The same widening is available to `make test`, cheaply.** Its prototype-against-production
comparison still runs on 2D grids at sides 10 and 20 alone, and `graphs.h` holds the 3D and random
builders. Worth knowing before relying on that check: the prototypes carry no maintained clique
degree, so a defect in production's encoding is invisible to it at any size. See `docs/TODO.md`,
the first of the five ordering questions.

**3. And only then the performance work below, which is a PARK rather than a queue.** One bounded
attempt, worth an hour whenever the amd branch is picked up again. The branch is in a good state to
leave without it: `Amd3` reproduces `AMD_2` at 2.32x, `Amd2` runs at 2.28x with better fill, `Mmd3`
at 1.26x with fill matching genmmd. Nothing here blocks anything.

## The state, in one paragraph

`amd3` is aligned: production `Amd3` returns `AMD_2`'s permutation exactly, up to the postorder it
deliberately does not do, on the seven examples, 2D grids to 140, 3D grids to 24 and nine random
patterns. AMD3 runs at about 2.32x the vendored routine at 140 a side and MMD3 at about 1.26x its
own. The gap is NOT algorithmic, NOT integer width, NOT the number of arrays touched, NOT
instruction count, and NOT any single line: all five measured and rejected. What is left is IPC,
and the profile is stall-shaped.

## The parked attempt

Give cliques their own mark space in `QuotientGraph`, then fold liveness into the vertex marks and
delete `mEliminated`.

```
mMark becomes 2n:  vertices at [v], cliques at [c + n]
                   (the amd driver already does exactly this with cliqueStamp)

then the membership test becomes ONE load:
    mMark[v] <  mTag     live, not yet seen this step
    mMark[v] == mTag     seen this step
    mMark[v] == GONE     eliminated, permanently   (GONE above any tag reachable)

merge(), number() and massEliminate() write GONE where they now set the byte
```

This deletes `mEliminated` and the `mLiveMerges` branch beside it: `(!live || mEliminated[v] == 0)`
goes from all forty sites, and with it a dependent byte load that `Amd.cpp` does not have, because
its liveness rides on `Nv`'s sign.

**Why it is expected to pay.** Forty byte loads in compiled `orderAmd3`, zero in the whole of
`Amd.cpp`. Each is `ldrsw` for an index, `ldrb` at that index, branch on the result: a serial chain
that costs cycles without costing instructions, which is the shape the counters describe.

**Why it might not.** This has not been measured, and five structural hypotheses failed on
2026-08-08 before two succeeded. Treat the mechanism as a reason to try, not as a prediction.

## How to verify it

The oracle is exact and it is the whole safety net:

```
make test                        283 with private/, 269 without
experiments/ordering:
  make test                      63 twin and prototype-against-production checks
  make aligned                   both alignment checks, 38 cases each
```

Every ordering's permutation must be UNCHANGED. This is a re-schedule, not an algorithm change.

**`mmd2` is the case that catches a wrong invariant.** A previous attempt at the same goal, folding
liveness into the weight instead, looked exactly equivalent and produced 201 entries for 200
vertices on a random `mmd2` pattern: an eliminated vertex with weight 1, sitting in a live adjacency
list. `make test` catches it, but check `mmd2` specifically and early.

**And run it under a sanitizer once**, which is cheap and is now known to be worth it:
`-fsanitize=address,undefined` over the six drivers on 2D and 3D grids. That is what found the
arena use-after-free, and a change that moves what a mark array means is exactly the shape that
would introduce another.

Also needed: an overflow guard on `mTag`, as `Amd3`'s `w` array has with `clearFlag`, since `mTag`
only ever increments and `GONE` must stay above it.

## If it does not pay

Stop rather than continue. The alternative account of the remaining gap is data layout, one `Iw`
pool against our `mSource` plus a separate clique arena, with per-vertex arrays as independent heap
allocations. The fix for that is a pooled-storage redesign of `QuotientGraph`, which is a much
bigger change, trades away part of what the shared class buys, and should be a deliberate decision
rather than the next thing tried.

The measurement that would choose between the two is L1D miss counts, ours against the vendored
routine's. `benchmarks/README.md` records how far an attempt to get them got and where it stopped.

## Also open, and unrelated to any of the above

- **Whether LIFO is genuinely better than FIFO**, or genmmd is simply a good ordering. The two
  branches now disagree: aligning MMD improved our fill, aligning AMD costs us 6.5 percent on
  grids, which is what makes it worth answering. `experiments/ordering/REPORT.md` has the question,
  the experiment and the caution that both figures are 2D. Priority 1 above is what it waits on.
- **The clique arena still never reclaims.** The reserve added on 2026-08-09 stops it moving under
  a walk, which was a correctness matter, and it does nothing about the arena growing toward
  nnz(L) where the live cliques fit in nnz(A). `Amd.cpp` compacts in place and counts it in
  `AMD_NCMPA`, which reads 1 for a whole 140x140 run. Section 5.3 of
  `archive/sparse_factorization.md` has the conservation argument that makes that possible, and the
  constructor's own comment already names reclaiming as the real fix.
- **The counter sweep.** Loop counters against sizes should be `std::int32_t` with the bound cast
  once; the size ARRAYS stay `std::size_t` because they are arithmetic operands and narrowing them
  buys a cast at every one. Measured at zero for speed at the hottest loop, so this is a clarity
  change and not a performance one.
