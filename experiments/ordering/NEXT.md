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

**And the ordering benchmark measures cubic grids.** `make run3d` and `make scale3d`, beside
`run2d` and `scale2d`, both axes now named in every target. It also carries an `AMDraw` column, the
vendored AMD's raw elimination order through the same hook the acceptance test uses, so `AMD3` has
something to sit against that agrees by construction rather than nearly. What the first run found
is in `benchmarks/ordering/README.md` under "Cubic grids, 2026-08-09", and the short version is
that three standing claims were square-grid artifacts: our tie-break beating AMD's, MMD being the
ordering to beat, and the LIFO question having an amd-side answer.

**And the mmd branch has an acceptance test at last.** `make mmdorder` compares production `Mmd3`
against genmmd's elimination order on the same four shapes, 38 cases, all matching. Until then the
mmd alignment rested on a scratch probe from 2026-08-07 that died with its session and on the
benchmark's fill column, which `MMD3.md` iteration 6 shows is not sufficient. It needs no hook,
genmmd emitting that order directly, which is why it is forty lines against the amd machinery.
`make aligned` runs both.

---

## Priority

**1. A matrix that is not a grid.** Cubic grids landed on 2026-08-09 and answered what they were
brought in for: the claim that our tie-break beats AMD's was a square-grid artifact, `Amd2` reading
`-5.5, +2.5, -2.6, +2.5, -4.6` percent on cubes against a monotone `-1.7` to `-6.5` on squares. But
two families of structured grid is two points on one axis, and nothing in this tree has yet been
ordered that came from a real problem. That is the oldest open item on `docs/TODO.md` and it is now
the top of this list. `benchmarks/pipeline` is still square grids alone, so every break-even figure
it carries is one family's.

**2. Count the hash pass's pairs, 2D against 3D. DONE 2026-08-09, and it was a defect.** The count
was taken, and against the vendored routine on the same graphs and for the SAME MERGES we were
testing 19.0 pairs per pivot at 140 a side where it tests 0.333, and 155.3 at 26 cubed where it
tests 0.484. The cause is the key: its incidence half was multiplied by a stride of `n + 1` and
then reduced modulo the same number, which annihilates it exactly, so the bucket was a function of
the adjacency alone. Two lines in nine files. On alpamayo `AMD3` on cubes goes from about 3.0x the
vendored routine to 1.44x and `AMD2` from 2.85x to 1.40x, with both controls unmoved. It is ledger
entry 8, `AMD3.md` iterations 21 to 24 are the narrative, and `docs/DESIGN_DECISIONS.md`
(2026-08-09) carries why five separate oracles were blind to it.

**What that leaves of item 4 below, which is a re-pricing rather than a deletion.** The parked
proposals were deprioritized because they change the SHARED quotient graph and would help `AMD1`
equally, and `AMD1` was the half that was already fine. That argument has been consumed: the extras
are no longer the story, and what is left is per-element work in walks both branches share. Taken
per pivot the excess over the vendored routine is 70 ns in 2D and 149 on cubes, against 6.28 and
12.72 clique members per pivot, so the residual tracks ELEMENTS WALKED rather than pivots, which is
what those proposals are about. The families have also swapped roles: the extras now cost about 25
percent in 2D and nothing on cubes, so what remains of the amd gap is a square-grid effect.

**The original text of this item, kept because the reasoning is what produced the finding.** The
cheapest measurement on this page and the one with the most behind it. `AMD1` is flat across
both problem families at about 1.2 to 1.8x while
`AMD3` goes 2.3x to 3.0x, so the amd branch's whole degradation on cubic grids is in the extras,
and the hash pass is the only part of them whose cost is superlinear in clique size: its pair loop
is the sum of squared bucket sizes over `C[p]`. Count pairs tested per iteration and the bucket
size distribution at comparable n on both families. A flat count per pivot kills the hypothesis; a
growing one says the fix is a better filter rather than a faster loop. One counter beside one the
prototypes already keep. `benchmarks/ordering/README.md`, "What the two families say about where
the amd gap is", carries the tables.

**3. The same widening is available to `make test`, cheaply.** Its prototype-against-production
comparison still runs on 2D grids at sides 10 and 20 alone, and `graphs.h` holds the 3D and random
builders. Worth knowing before relying on that check: the prototypes carry no maintained clique
degree, so a defect in production's encoding is invisible to it at any size. See `docs/TODO.md`,
the first of the five ordering questions.

**4. And only then the performance work below, which is a PARK rather than a queue.** One bounded
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

**READ THIS BEFORE STARTING IT, added 2026-08-09.** The cubic benchmark says this may be aimed at
the wrong half. Split the amd branch into its base and its extras, ratios against the vendored
routine: `AMD1` costs about 1.2 to 1.8x on BOTH families, while `AMD3` goes from 2.3x in 2D to
3.0x in 3D. So the degradation is entirely in aggressive absorption and hash supervariable
detection, and this proposal changes the SHARED QUOTIENT GRAPH, which would help `AMD1` exactly as
much and `AMD1` is the half that is already fine. The same objection applies to the locality
hypothesis below. Neither is refuted, and both may still be worth their 2D value, but a fixed cost
in the shared class cannot explain a gap that appears only with the extras on and only on one
family. The suspect is the hash pass, whose pair loop is quadratic in bucket size while `C[p]`
grows with the family, and the measurement that would settle it is a COUNT rather than a profile:
pairs tested per iteration and the bucket size distribution, 2D against 3D at comparable n. That is
one counter beside one the prototypes already keep. `benchmarks/ordering/README.md`, "What the two
families say about where the amd gap is", has the tables and the argument.

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
