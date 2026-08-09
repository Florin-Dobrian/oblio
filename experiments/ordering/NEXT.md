# NEXT: 3D grids first, and one last bounded attempt at the amd constant factor

A handover note between sessions, written 2026-08-09. **It is meant to be deleted** once the change
below is done or abandoned. Everything in it that outlives the task has already been written
somewhere durable, and this file only points at those places:

- `docs/DESIGN_DECISIONS.md`, entry "what the vendored AMD's speed is made of", has the
  measurements, the three blocked routes and why each is blocked, and the proposal restated
  properly.
- `experiments/ordering/README.md`, section "Aligning a layer against a vendored routine", has the
  method both alignments used.
- `benchmarks/README.md` has what Instruments can and cannot be made to do from the command line,
  and the two rules about counting against profiling that cost a day between them.

If you find yourself about to record something here that a later reader would want, put it in one
of those instead. The previous `NEXT.md` went stale precisely because it accumulated content that
belonged elsewhere.

---

## Priority

**1. Widen the alignment check, and finish it before anything else.** `make raworder` compares
production `Amd3` against the vendored AMD's raw elimination order. As committed it runs the seven
example graphs and 2D grids from 4 a side to 140, all matching. That is one SHAPE at many sizes,
which exercises scale and not mechanism.

An attempt to widen it on 2026-08-09 found divergences immediately, and was NOT committed because
the cause is unestablished:

```
3D grids, 4^3 to 24^3      sizes MATCH, contents DIFFER   (4^3 first differs at position 32)
random n=2000, deg 3, 6    sizes match, contents differ at position 0
random n=2000, deg 12      SIZE MISMATCH, vendored short by about 186 of 2000
```

**The size mismatch says the HOOK is incomplete**, not that `Amd3` is wrong. Every pivot's member
count equals `nvpiv`, instrumented and checked, so the main loop is fully covered and yet 186
vertices are numbered on a path `hook_amd.py` does not see. Find that path first: until it is
found, the random cases mean nothing either way.

**The 3D failures are a different matter and may be real.** Sizes match there, so every vertex is
accounted for and only the ORDER differs. That is either a genuine `Amd3` divergence on a shape
never tested, or the hook mis-ordering members while counting them correctly. The smallest failing
case is a 4x4x4 grid diverging at position 32, small enough to read in full.

A re-schedule is only safe if the oracle is sound. Optimizing against an alignment never checked
outside 2D would be building on sand.

The widening also wants a `graphs.h`: the seven example graphs exist twice already, in
`vendored.cpp` and `production.cpp`, and a third copy was about to go into `raworder.cpp`. One
definition, so a graph added for one driver is available to the others. The drivers answer
different questions; their inputs are not a reason to diverge.

**2. Then 3D grids in the ordering BENCHMARK**, which is a separate task from the check above.
Every fill conclusion in this work is from square grids, which is one problem family and the
flattering one. The specific claim waiting on it is that our tie-break beats AMD's: `Amd2` fills
6.5 percent BELOW the vendored routine at 140 a side, while `Mmd3` matches genmmd exactly. If that
advantage is a 2D artifact we should know, because fill drives the factorization and the
factorization is where the time is. A 2x ordering costs well under 15 percent of a one-shot solve;
a fill difference costs more than that, permanently, on every factorization.

**3. And only then the performance work below, which is a PARK rather than a queue.** One bounded
attempt, worth an hour whenever the amd branch is picked up again. The branch is in a good state to
leave without it: `Amd3` reproduces `AMD_2` at 2.32x, `Amd2` runs at 2.28x with better fill, `Mmd3`
at 1.26x with fill matching genmmd. Nothing here blocks anything.

## The state, in one paragraph

`amd3` is aligned: production `Amd3` returns `AMD_2`'s permutation exactly, up to the postorder it
deliberately does not do, ON THE SEVEN EXAMPLES AND ON 2D GRIDS. Outside those two shapes it is
unverified, which is the priority above. Both alignments are committed. AMD3 runs at about 2.32x
the vendored routine at 140 a side and MMD3 at about 1.26x its own. The gap is NOT algorithmic,
NOT integer
width, NOT the number of arrays touched, NOT instruction count, and NOT any single line: all five
measured and rejected. What is left is IPC, and the profile is stall-shaped.

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
  ./production_cpp amd3 grid N   must equal the probe's raworder at 32, 100, 140
```

Every ordering's permutation must be UNCHANGED. This is a re-schedule, not an algorithm change.

**`mmd2` is the case that catches a wrong invariant.** A previous attempt at the same goal, folding
liveness into the weight instead, looked exactly equivalent and produced 201 entries for 200
vertices on a random `mmd2` pattern: an eliminated vertex with weight 1, sitting in a live adjacency
list. `make test` catches it, but check `mmd2` specifically and early.

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
  the experiment and the caution that both figures are 2D.
- **The counter sweep.** Loop counters against sizes should be `std::int32_t` with the bound cast
  once; the size ARRAYS stay `std::size_t` because they are arithmetic operands and narrowing them
  buys a cast at every one. Measured at zero for speed at the hottest loop, so this is a clarity
  change and not a performance one.
