# NEXT: a matrix that is not a grid, and what is left of the amd constant factor

A handover note between sessions, rewritten 2026-08-10 when item 0's instrument was
finished and its first candidate landed. **It is
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

**0. The work gap, DONE as an instrument and OPEN as a question, 2026-08-10.** The per-pass
inventory named here is built, both halves, and it produced the first change in five attempts to
move this gap. Everything durable from it is in `benchmarks/ordering/README.md` under "The per-pass
inventory" and "What the inventory was worth", in `AMD3.md` iteration 25, and in
`docs/DESIGN_DECISIONS.md` (2026-08-10). The short version:

- We make **2.12x** `AMD_2`'s element visits on both families against 1.56x and 1.61x its useful
  cycles, so we execute about 0.74x its work per visit. Nine sweeps over `C[p]` per pivot against
  its four, four walks of `I[u]` against its two.
- Our half as previously recorded here was wrong by up to a factor of two on six of fifteen
  columns, the instrumented shared class having been counted twice. The corrected figures are
  112.26 visits per pivot at 140 a side and 276.25 at 26 cubed, not 155.2 and 374.8. The prune's
  incidence compaction is NOT oversized, contrary to what this file said: it equals `AMD_2`'s scan
  2 exactly.
- **The key fusion landed** and is worth 4.4 to 7.0 percent in 2D across six sizes and 5 to 14 on
  cubes. The 2026-08-08 version of it failed because it stored keys in an array of size n, not
  because of the fusion.
- **Deleting the `degme` re-take is a recorded NULL**, 0 percent in 2D and 1 to 3 percent slower
  on cubes, all of it inside this benchmark's plus or minus 3 percent floor. Not a regression. The
  sign was positive at all five cubic sizes in two runs, which is weak evidence of a small cost
  and no more.
- **THE FLOOR IS PLUS OR MINUS 3 PERCENT**, measured between the fusion landing and the vehicle
  being removed, when `Amd3B` held a verbatim copy of `Amd3` and the two benchmark columns were
  the same code timed twice. A single figure under about 4 percent is not a result; what rescues a
  small effect is consistency of sign across sizes. **Worth arranging again**: whenever the next
  vehicle exists, run the benchmark once with it still a verbatim copy before putting anything in
  it, which costs one column and gives every other column an error bar under identical
  conditions.

**AND THE WALKS ARE DONE, 2026-08-10, later the same day.** The first scan is folded into the
prune, so `I[u]` is walked twice per pivot and `A[u]` once, which is `AMD_2`'s count exactly.
`AMD3` is at 1.01x, 1.02x and 1.07x the vendored routine at 12, 16 and 20 a side and 1.09x at 26,
which is parity on cubic grids, with 2D improving 0 to 8 percent as well. `AMD3.md` iteration 26
and `docs/DESIGN_DECISIONS.md` (2026-08-10, "the algorithm was the smaller half") carry it.

**What is open, in the order we would take it.**

- **A profile, now that cubes are at parity and the shape of the gap has changed.** CPU Counters on
  `amd3` against `AMD` at 26 and 32 cubed. The cubic gap used to be pure work at identical
  efficiency; if that still holds at 1.09x there is little left, and if efficiency has moved the
  ranking needs redoing. At 32 a side the gap is 1.16x against 1.09x at 26, so it now GROWS with
  size on this family too, which it did not before and which is a different question.
- **Price a change by its STREAMS as well as its visits.** The pass inventory counts visits and is
  silent about size-n arrays, and on 2026-08-10 the arrays were the larger term in 2D: the same
  fold measured 12 percent slower there with two extra vectors and 0 to 8 percent faster without
  them. There is no instrument for this and one would be cheap.
- **`Amd2` and `Amd2B`, 2026-08-10: not by the cheap route, and PARKED rather than closed.**
  Checked on a scratch copy: fusing there as `Amd3` does moves the permutation on all ten grids
  tried. They refile inside their single-pass bound, so that loop's direction is already a
  tie-break input and cannot also serve the hash chain, which wants the opposite direction under
  head insertion. Tail insertion would preserve the order and is untried; it needs a `hashTail`
  array of size n, which is the footprint that made the first version of this fusion measure
  nothing, so expect nothing. Not worth a measurement before the profile below. `Amd1` and `Amd1B`
  have no hash and so no key. Only `Amd3` fuses free, because entry 4 moved its refile below the
  hash.
- **The stamp and the mass-elimination sweep are DEMOTED, not queued.** Both were next on the
  reasoning that produced the `degme` deletion, and that reasoning now has a counterexample.
- **`tmp/Amd3I.cpp` is the pre-fusion driver.** The pass inventory reproduces the OLD `Amd3` until
  it is re-derived, which is fine for the before-and-after already recorded and wrong for anything
  further.

**And the 2D scaling divergence, which is larger than any of the above and is not made of passes.**
`AMD3` runs at 1.54x the vendored routine at 64 a side and 2.07x at 400, growing monotonically,
where `MMD3` over the same quotient graph shows no trend at all across those seven sizes. A
constant-factor fix takes five percent off a growing curve and leaves it growing. Two profiles of
`amd3` at 64 and at 400, compared against each other rather than against the vendored routine,
would say what grows.

**The standing caution, unchanged and now with a third instance.** A count locates work; it does
not price it. The 2026-08-09 fusions were 6.7 percent of the visits and moved cycles by 0.25
percent; the `degme` deletion was 4 percent of the visits and moved cycles the wrong way.

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

**2b. Taking `AMD1B` and `AMD2B` out of the public enum: DEFERRED 2026-08-10, deliberately.** The
reasoning below still holds and the work is still worth doing eventually. It is not being done now
because a third B name is in circulation that is not an oracle at all: `Amd3B` is a VEHICLE,
holding one candidate change to `Amd3` at a time. It is SCRATCH and was never committed: it existed
on 2026-08-10 while the two candidates below were being priced and was removed with them. Item 2c
says how the next one should be added so that it touches far less. Deciding the fate of the
other two while that name comes and goes would be churn in the middle of an open investigation.
Revisit when the amd branch is quiet.

**2c. HOW TO ADD THE NEXT VEHICLE.** `Amd3B` was wired in as a full ordering on 2026-08-10, an
enumerator with everything that follows from one, and all of it was reverted a day later. One site
was missed on the way in: three files under `examples/` switch over `Ordering`, so `make examples`
warned. Do it as a FREE FUNCTION instead, which is the pattern `order_timing.cpp`'s `AMDraw` column
already uses on the stated grounds that an enumerator "would put a benchmark's oracle into the
library's public enum and into every switch over it". Two new files, two build entries, one column
in each benchmark driver calling `orderAmd3B` directly, and a local identity check against
`orderAmd3`. No enumerator, so no dispatch, no adapter, no `examples/` arm, no `test_order`
assertions, no `test_pipeline` sweep entry, and no move in `docs/TESTING_SPECIFICATION.md`.

**3. The same widening is available to `make test`, cheaply.** Its prototype-against-production
comparison still runs on 2D grids at sides 10 and 20 alone, and `graphs.h` holds the 3D and random
builders. Worth knowing before relying on that check: the prototypes carry no maintained clique
degree, so a defect in production's encoding is invisible to it at any size. See `docs/TODO.md`,
the first of the five ordering questions.

**4. And only then the performance work below, which is a PARK rather than a queue.** One bounded
attempt, worth an hour whenever the amd branch is picked up again. Nothing here blocks anything.

**Re-priced 2026-08-10, and it is now the WEAKER of the two parked ideas rather than the stronger.**
The proposal below deletes `mEliminated` by giving cliques their own mark space, which removes a
dependent byte load from a hot walk. That is a good thing to want. But its sibling argument, that
deleting work from a sweep is reliably worth something, is what the `degme` deletion tested and
falsified: a provably redundant sweep cost one to three percent on cubes when removed. This
proposal is not that deletion and may well pay, since it shortens a walk rather than removing a
pass. The point is only that it no longer inherits an argument it used to.

## The state, in one paragraph

`amd3` is aligned: production `Amd3` returns `AMD_2`'s permutation exactly, up to the postorder it
deliberately does not do, on the seven examples, 2D grids to 140, 3D grids to 24 and nine random
patterns. After the two fusions of 2026-08-10, the hash key into the bound and the first scan into
the prune, `AMD3` runs at 1.01 to 1.16x the vendored routine on cubic grids from 12 to 32 a side,
which is parity, and 1.44 to 2.00x on square grids from 32 to 400. `MMD3` is at 1.25 to 1.46x its
own with no trend in n. `AMD3` is FASTER than both `AMD1` and `AMD2` at every cubic size while
returning the vendored permutation, which `AMD2` does not. What is left is the 2D penalty, which
grows with n and is stalls rather than work, and a cubic gap that has only now started to grow with
size at all. The gap is NOT algorithmic, NOT integer width, NOT the number of arrays
touched, and NOT any single line. It IS partly the number of passes, which is new: 2.12x the
element visits on both families, of which the fusion removed a fifth. What remains is that plus a
2D-only memory penalty that GROWS with n, and the second of those is now the larger question.

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
