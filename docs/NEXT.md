# NEXT: the amd branch's growth is gone, and where the differential goes next

**REWRITTEN AT THE TOP ON 2026-08-16.** The section immediately below is the current state;
everything under "Priority" and after is older and several of its items are now closed, marked
where they are. The durable account of this session is `docs/DESIGN_DECISIONS.md` (2026-08-16),
`benchmarks/ordering/README.md` and `experiments/ordering/README.md`; this file only points at
them.

## Where the amd branch stands, 2026-08-16

**The 2D growth is gone.** `AMD3` read 1.25x the vendored routine at 32 a side and 1.82x at 400,
rising monotonically; it now reads about 1.13 to 1.43x with no trend from 100 a side up. Seven
folds did it, and ONE of them removed the slope: the clique descriptor moved into the dead pivot's
own `mRun` entry, retiring `mCliquePtr` and `mCliqueSize`. The other six moved the column down by a
constant. Seventeen entity-indexed streams became eleven, against `AMD_2`'s nine.

**Nothing computed changed.** Every permutation from all ten drivers is identical to the
`be589f2` baseline, checked as a digest over 73 grid sizes per driver at each step.

**What the differential established, and it killed three of my hypotheses in a row.** Both codes do
the same visits per pivot, 1.011 in 2D and 1.010 on cubes, in every pass. So the growth was never
work. Not the clique arena either, whose excess misses FALL with n. Not the hash. Not `clear_flag`,
which never fires. And not single-level cache locality, an L1 model putting the excess at a constant
0.09 misses per visit. Read `experiments/ordering/README.md`, the 2026-08-16 section, before
forming a new hypothesis: it also records where the first version of that counter was blind.

### The four things to pick up

**1. `AMD1` AND `AMD2` STILL CLIMB IN 2D**, 1.05 to 1.22x and 1.36 to 1.47x, where `AMD3` is flat.
The descriptor fold is in the shared class and reached them for free, so what remains is
driver-side. This is the open front and the differential is the instrument for it.

**2. FOLD THE B LAYERS IN AND RETIRE THEM.** `AMD1B` and `AMD2B` now satisfy the stop condition in
`Amd1B.h` for the first time: permutation-identical AND faster, 4 to 16 percent on cubes and even
or better in 2D. Their 2D penalty was never the fused schedule, it was `ApproximateScan` crossing
three arrays from the prune; with the tagged `W` the crossing is one. So `Amd1` and `Amd2` take the
fused `TaggedScan` eliminate, the four B files go, and `ApproximateScan` and its overload go with
them, having no users left. `Amd3B` and the `AMD3f` control retire in the same sweep, `Amd3` now
being identical to `Amd3B`.

**3. A DENSER LADDER.** The flatness claim rests on seven sides, geometric but sparse in the middle.
The driver takes explicit sides, so `./order_timing_cpp amd 2d 32 45 64 90 128 181 256 362 400`
needs no code. Worth doing before the claim is quoted anywhere outside this tree.

**4. THE TEN-DRIVER DIGEST WANTS A HOME.** It is a throwaway in `/tmp` that hashes every driver's
permutation over 73 grids and compares against a recorded baseline. It caught nothing this session,
but it is what made a shared-class change safe to attempt at all, and it is strictly stronger than
the `Amd1B == Amd1` pair check it is about to replace.

### Two things worth carrying into any next fold

**A fold can measure nothing and still be why another one works.** The dead-clique test measured
nothing. The tagged `W` measured nothing in `Amd2`. Yet the tagged `W` is what made the fused scan
viable in both B layers. Judging a fold by its own column alone would have discarded the
precondition and kept the failure.

**The footprint trade, now seen three times, is a rule.** A schedule change that saves visits but
adds an array crossing a pass boundary tends to lose at large n. Fold 7 is the general answer to it:
put the crossing in a slot that already exists.

---

A handover note between sessions, rewritten 2026-08-11 after three commits added a real-matrix
benchmark, a scaling benchmark and the three reports drawn from them, and updated 2026-08-10 before
that when three commits closed the constant factor on the amd branch. Appended to on 2026-08-14
with two rounds that touched no algorithm, an idiom sweep and a Makefile consolidation, and again
later that day with three commits on the ordering experiments, all under "Closed since the last
note" and none of them changing a permutation or a fill figure. Appended to again on 2026-08-15
with the MMD3 storage investigation, which is uncommitted work and is item 0 below. **It is meant to be
deleted** once the items below are done or
abandoned.
Everything in it that outlives the task is already somewhere durable and this file only points at
those places:

- `docs/DESIGN_DECISIONS.md`, three 2026-08-10 entries: "a null result measures an implementation",
  "the algorithm was the smaller half", and the earlier "what the vendored AMD's speed is made of".
- `experiments/ordering/AMD3.md`, iterations 25 and 26, the two fusions and how they were found.
- `experiments/ordering/REPORT.md`, "The one gap we can explain", now closed, and the new vendored
  against vendored section beside finding 1.
- `benchmarks/ordering/README.md`, the per-pass inventory, the noise floor, and both scaling
  ladders. `benchmarks/pipeline/README.md`, where the ordering phase gets priced.
- `benchmarks/matrices/README.md`, everything the real matrices said, including the two ordering
  investigations of 2026-08-11 in full. Its `ACCURACY.md` and `PERFORMANCE.md` are the reports
  drawn from it, written for readers outside this tree.
- `benchmarks/pipeline/README.md`, the scaling ladders in full including the two rungs beyond the
  committed six, and why those matrices carry a dominance margin the short ladder's do not. Its
  `SCALING.md` is the third report.

If you find yourself about to record something here that a later reader would want, put it in one
of those instead. The previous `NEXT.md` went stale precisely because it accumulated content that
belonged elsewhere.

---

## Read this first, and the rest only if you are picking up that item

**ITEM 0 IS CLOSED, 2026-08-15, and the answer was not what this file spent a day predicting.**
`MMD3` runs at 1.05 to 1.17x genmmd on square grids and 0.77 to 1.12x on cubes, against 1.35 to
1.48x that morning, and `MMD2` reads the same to within the noise. Nothing about what is computed
changed: every permutation and every nnz(L) is identical, on all 137 acceptance cases across nine
orderings.

**That MMD2 and MMD3 coincide is the useful shape of the result.** They are different orderings
with different mechanisms; what they share is the quotient graph and its clique arena, so both are
now bounded by it rather than by anything of their own. And the trade the library makes is worth
stating in one line: our `std::vector` layer costs a few percent that genmmd does not pay, having
raw arrays in registers for a whole run, and the SECOND ARENA more than covers it. `Mmd3B`, which
is `Mmd3` on genmmd's single nnz(A) storage with every encoding fold present in both, reads 1.12 to
1.23x where `Mmd3` reads 1.05 to 1.17x.

**The cause was the number of ARRAYS a vertex's state lives in**, not placement, not chaining, not
the number of passes, and not the algorithm. genmmd indexes five arrays by a vertex and we indexed
eleven, because each of its arrays answers several questions at once, told apart by sign or by a
reserved value. Four folds and one packing closed it, each of them genmmd's own encoding. The full
account, the pass-by-pass differential that located it, and the three lessons are in
`docs/DESIGN_DECISIONS.md` (2026-08-15). **Read that entry rather than the section this replaced.**

**Two things this file said are now known wrong, and both misdirected the search.** It opened with
"we touch 30 per cent fewer arena entries and take about 50 per cent longer", concluding "so it is
not more work"; counting every loop rather than the source arena alone, we made 1.18x genmmd's
visits in 2D. And it recorded the 2D penalty as stalls rather than work, which was right about the
symptom and wrong about the remedy: our own code was 1.13x its instructions all along.

**And the storage question is answered in the opposite direction from the one it was asked in.**
`Mmd3B` carries `Mmd3`'s algorithm on genmmd's storage, cliques in the dead segment of their own
pivot with no second arena. Every encoding fold is now in both files, so storage is the only
difference left, and ours wins on every axis: 1.02 to 1.19x genmmd in 2D against `Mmd3B`'s 1.15 to
1.38x, 14.22M instructions against 16.61M, 119331 D1 read misses against 123510. **Spending nnz(L)
on a second arena buys speed.** `Mmd3B` therefore STAYS, as the standing equal-encoding comparison
against the vendored storage scheme, which is a change to the stop condition its own header states.

**What is open on the mmd branch, and it is one item.** The TAG SCHEME. genmmd's refresh puts the
element tag ABOVE the per-vertex tags, `mt = tag + md0`, so an element member fails
`marker[nd] < *tag` automatically and needs no second comparison; ours draws both from one counter,
so `elementTag` is below `vertexTag` and every entry in the q2h and qxh paths pays an explicit
`m == elementTag` test. It is a change to what is computed rather than to how it is stored, which
is the category that actually paid, and the refresh is about a third of the run. Untried.

**Everything else on this branch was tried and did not pay.** Four attempts at the container layer,
each failing for its own reason, all in the "2026-08-15, later" design entry: the stamping fold
ported to `Mmd3`, an arena cursor in place of `push_back`, raw bases in place of the accessors, and
q2h indexed rather than looped. Read that entry before reaching for any of them.

**And `Mmd3B` has finished answering its question.** It is now genmmd's data structure essentially
exactly: one array of nnz(A), cliques in their pivot's dead segment, negative links, a value
terminator, no clique length array, no liveness array, the degree list in one link array. Its
obligation from here is to stay ENCODING-IDENTICAL to `Mmd3`, so that the only difference between
them is storage; a fold that lands in `QuotientGraph` lands there too, or the comparison quietly
stops being about storage.

**And the amd branch has had NONE of this.** Its five drivers still carry `degrees`, `outside`,
`cliqueDegree`, `explicitPart`, `hashHead` and `hashNext`, plus a 2n `mark`, and `AMD_2` allocates
none of them: it overlays its hash buckets on the degree heads and keeps the running key in a link.
`AMD3` gained from the shared-class folds alone, reading 1.28 to 1.87x in 2D and 0.97 to 1.26x on
cubes, but nothing driver-side has been touched. That is where the remaining headroom is.

**One item to watch there.** The absorbed-clique stamp was deleted from the shared class, replaced
by `mCliqueSize[c] != 0`. It helps mmd, which the two folds it enabled more than pay for, and it
costs amd, whose aggressive absorption makes that loop heavy and which has no compensating fold
yet. Splitting `mMark` to 2n, vertices at `[v]` and cliques at `[c + n]`, would let the stamp come
back hot and is what both references effectively have, neither of them sharing one stamp array
between the two kinds. Costs n extra int32, and footprint has killed three changes in this tree.

**TWO THINGS ABOUT THE INSTRUMENTS, and they outlast the changes.** A two to three percent movement
in instruction count is INVISIBLE on alpamayo, confirmed in both directions with controls in the
same run, so a counter movement of that size decides nothing by itself. And cachegrind's simulated
cache misses are not comparable across shell invocations, shifting about 17 percent with heap
placement; instructions and data references are exact. Compare builds back to back inside one
invocation, and verify a revert with a DIFF rather than with a counter, which is how two leftovers
survived one that day.

**EVERY TIMING FIGURE IN THIS TREE PREDATES 2026-08-15 AND UNDERSTATES US BY 20 TO 30 PERCENT.**
`benchmarks/ordering/README.md`, `benchmarks/pipeline/README.md`,
`experiments/ordering/README.md`, `experiments/ordering/REPORT.md`, `AMD3.md`, `MMD3.md` and
`docs/TODO.md` all carry ratios measured before it. Fill figures are unaffected, nothing having
moved. Each file carries a dated superseding note rather than being rewritten, since a dated
measurement is a record of a run.

**Where the amd branch stands.** `Amd3` returns `AMD_2`'s raw elimination order exactly on all 38
acceptance cases and now runs close to it: 0.83 to 0.89 ms at 16 cubed over eight runs where the
vendored routine reads 0.74 to 0.86, so they overlap, rising to about 1.2x at 32 a side. In 2D it
is 1.33x at 32 a side and about 1.95x at 400. Three days earlier it was 3.0x on cubes. `Mmd3` is
0.87 to 1.05x genmmd on cubes and 1.26 to 1.55x in 2D, with no trend in n on either.

**THE HASH KEY IS NOW `AMD_2`'S EXACTLY, 2026-08-16.** All four divergences listed anywhere in
this file are closed in `Amd2`, `Amd2B` and `Amd3`: no `+ 1` on a term, the pivot's own clique
excluded, modulus `n` rather than `n + 1`, and a `uint32` accumulator that WRAPS at 2^32 the way
Amd.cpp's `UInt hval` does, one reduction at the end rather than one per term. It changed no
permutation, checked over 730 across ten drivers, and it corrected the grouping: our `pair` and
`stamp` counts now equal `AMD_2`'s digit for digit, where before we were UNDER-grouping on cubes.
Worth about 4 percent. See `docs/DESIGN_DECISIONS.md` (2026-08-16).

**What was done, 2026-08-08 to 08-10**, three commits: the hash key defect (ledger entry 8, worth a
factor of two to three on cubes), the key folded into the bound pass, and the first scan folded into
the prune. **The walk axis is finished**: `Amd3` walks `I[u]` twice per pivot and `A[u]` once,
which is `AMD_2`'s count exactly.

**Item 4 of the previous list is done.** Real matrices are in, `benchmarks/matrices/` holds the
benchmark and two of the three reports, and they changed what the other items are measured
against. What they
said about the orderings, in order of how much it should affect the plan:

- **The fill gap between MMD and AMD collapses on real structure**, to one to three percent from up
  to thirteen on grids. Grids are where the gap lives, because nearly every live vertex has the
  same degree and the tie-break decides almost every pick. **So a fill improvement measured on a
  grid is a claim about grids**, and the ladders in `benchmarks/ordering` should be read that way
  from now on.
- **Ordering is a fifth of a one-shot solve at the median**, and the whole analysis is 45 percent
  against the numeric factorization's 47. That is a stronger case for optimizing the ordering than
  the pipeline benchmark's grid figure of a tenth to a fifth, and it is measured on the input a
  caller actually has.
- **Two ordering pathologies turned up, and one has a documented remedy.** Both are written out in
  `benchmarks/matrices/README.md` under "Why some matrices order slowly".

**What is open, shortest first.**

1. **WE STILL DO NOT COMPUTE WHAT THE VENDORED ROUTINE COMPUTES.** FOUR divergences in the hash
   key, the mechanism that already cost several sessions once, a FIFTH in the mass-elimination
   test, and a SIXTH in the hash's storage, where we allocate 2n int32 that `Amd.cpp` does not.
   None of them is visible to any output we check. Nothing measured, and
   the fix is not obvious, but matching comes before anything clever here. Item 2a-hash, and it
   goes above the profile because the profile measures a routine that is not yet the routine we
   are trying to match.
2. **A profile that compares `amd3` against itself at two sizes**, 140 and 400 a side. It decides
   the next item and takes ten minutes. Item 2d.
2b. **THE PER-ARRAY-AT-SETUP CONTAINER EXPERIMENT IS ANSWERED, 2026-08-16, and not in the form it
   was posed.** It asked whether one allocation carved into arrays would beat separate vectors,
   `AMD_2` carving Pe, Nv, Head, Elen, Degree, W and Iw out of one `S`. That was the leading
   hypothesis for the 2D growth for one session and it was WRONG: `MMD3` has ten separate
   allocations against genmmd's five and is flat, so a doubling of allocations does not by itself
   produce growth. What did produce it was two random probes per clique visit into arrays indexed
   by a dead pivot's id, and folding the arrays away removed them. The container question is
   therefore still open as a question and is no longer a suspect.

3. **A dense-row threshold in `QuotientGraph`.** NEW on 2026-08-11 and the largest ordering-time
   item the real matrices found. A single vertex adjacent to everything makes minimum degree
   quadratic; the vendored AMD sets such rows aside before ordering, above `max(16, 10*sqrt(n))`
   entries, and places them last, and its own source says the cost of not doing so is O(n^2). Oblio
   has no such rule anywhere. On `GHS_indef/bloweybq`, one column of degree 10000 among 9992 of
   degree 5, removing it by hand takes MMD from 70.7 ms to 0.83 and AMD3 from 470 to 1.5, while the
   vendored AMD does not move. **Fill is unaffected**, so this is time only. It belongs in the
   shared quotient graph so all six drivers gain it at once.
4. **The descriptor struct**, `docs/TODO.md` question 3, if that profile says stalls. Item 2d.
5. **Narrow the one-dimensional sizes. THE ORDERING IS COMPLETE, 2026-08-11**; the symbolic and
   numeric phases remain and are the larger half. Item 2e.
6. **Why `Amd3` and the vendored `AMD` disagree on fill on real matrices.** NEW on 2026-08-11, and
   **the instrument for it now exists on the other branch, 2026-08-15**. `make mmdmatrices` runs
   the mmd alignment over `data/*/*.mtx` and reported 243 matched, 0 differed, 3 skipped on its
   first run, so `Mmd3` reproduces genmmd on real structure and not merely on grids. An
   `amdmatrices` is the same driver with the hooked oracle and the dense-row knob, and unlike the
   mmd one it should be EXPECTED to find something, which is what this item is about.
   `HB/bcsstk08`, the 4 percent case named below, is in that run and matches on the mmd side.
   `make amdorder` shows exact agreement on all 38 acceptance cases, which are grids and random
   patterns; on the 107-matrix performance set they differ on a minority, once by 4 percent on
   `HB/bcsstk08`, 29922 against 31153. **Ours fills less where they differ**, so it is not urgent,
   but it is a divergence the acceptance tests cannot see and the widening that would see it is the
   41 pattern files already sitting in `data/`.
7. **The dead-code finding generalizes, and nothing detects it.** NEW on 2026-08-14. Three
   inherited-and-redundant constructs turned up in `mmd2` and `mmd3` in one afternoon: `refile`,
   defined and never called in all three mmd layers; the evicted list and its stamp array, built
   and never read, because the element-walked refresh replaced the need for them; and an inert
   ternary whose branches were identical. None is a compiler warning, none moves a trace, and the
   twins agreed throughout because the dead thing was dead in both. All three arrived the same way,
   a layer copied from the one below and then given a pass that made part of the copy redundant.
   `mda2`, `mdam2`, `mdm2` and `amd1` through `amd4` were built that way and have not been looked
   at, nor has whether their display blocks build text outside the `SHOW_THRESHOLD` guard, which
   was wrong in six of the eight layers checked.

8. **What a `std::vector` actually costs us, measured on purpose rather than inferred.** NEW on
   2026-08-15. Our code carries a container layer genmmd does not: its arrays arrive as parameters
   and live in registers for a whole run, ours are members reached through accessors. On a 100x100
   grid that layer is 2.78M of `Mmd3B`'s 15.89M instructions, and roughly half of it, the
   `operator[]` half, is the element access genmmd also performs and merely has credited to its own
   source. The other half is `push_back` capacity tests, `size()` and growth arithmetic, which
   genmmd does not execute at all.

   **THREE ATTACKS ON IT ALL FAILED, each for its own reason**, and they are in
   `docs/DESIGN_DECISIONS.md` (2026-08-15, later): the stamping fold, an arena cursor in place of
   `push_back`, and raw bases in place of the accessors. The last was also TIMED on alpamayo and
   came out flat while costing 371403 instructions, which is what makes the question worth asking
   properly rather than probing at again.

   **THE EXPERIMENT.** The same small computation written twice, once over `std::vector` members
   reached through accessors and once over raw arrays passed as parameters, and timed against each
   other. It belongs in `experiments/`, being a study rather than a benchmark, and it wants two
   axes rather than one, because a grid conflates them:

   - **Per ACCESS**, in a hot loop: the case the three failed attacks were all aimed at, and where
     the compiler appears to be doing the work already.
   - **Per ARRAY, at setup**, which is the axis nothing has isolated. A pure diagonal is the case
     that shows it: with no elimination work at all, an `Mmd3` ordering is roughly a third
     `QuotientGraph` CONSTRUCTION and a sixth `orderAscending`, and construction allocates and
     initializes about ten size-n arrays where genmmd allocates five plus its 1-based copies. That
     is the array-count finding of the same day, moved from the loops into the constructor, and it
     is invisible on a grid because real work amortizes it.

   **Why it matters and where.** `benchmarks/matrices`, `make ordering`, has the evidence: the five
   pure-diagonal matrices, `Boeing/bcsstm39`, `Cunningham/m3plates`, `HB/bcsstm25`,
   `Oberwolfach/t3dl_e` and `Oberwolfach/t2dal_e`, read 2.03 to 2.28x genmmd after the prepass fold,
   where matrices with real work read 0.40 to 0.86x. Absolute times there are tenths of a
   millisecond, so this is not a performance item; it is the cleanest available measurement of what
   the abstraction costs, which is a thing worth KNOWING before the numeric side is tuned.

   **Two conditions on doing it at all**, both learned the same day. It must be timed on alpamayo,
   since a two to three percent instruction movement is invisible there in either direction and
   three of the four attempts were judged on counters alone. And builds must be compared back to
   back inside one invocation, cachegrind's simulated cache misses shifting about 17 percent with
   heap placement while instruction and data-reference counts are exact.

**Three things not to retry without reading why**, each of which cost a day or a wrong claim:

- The key fusion and the scan fold both FAILED first as versions carrying an array of size n. Price
  a re-schedule by what it walks AND by what it makes resident. The pass inventory counts only the
  first.
- `AMD3 / AMD` per row is a poor measurement: the vendored column moves 16 percent between runs
  where ours moves 7. Quote absolute times with the vendored range beside them.
- A benchmark column reached as a free function is timed differently from one reached through the
  enum, by up to 2.4 percent. Two columns compared must go down the same path.

---

## The MMD3 storage investigation, 2026-08-15: CLOSED

Everything this section held is superseded by `docs/DESIGN_DECISIONS.md` (2026-08-15), which
carries the result, the differential that found it, and the five hypotheses that failed on the way.
Two facts from it are worth having here because they are what a later reader would otherwise
re-derive:

- **The five dead hypotheses stay dead**: construction cost, the liveness array in the clique walks,
  the four-pass refresh preamble, index widths, and clique placement. `Mmd3B` implements genmmd's
  placement in full and the time did not move.
- **The renumbering experiment that pointed at placement was CONFOUNDED**, since shuffling adds a
  large cost to both routines and a large common term compresses any ratio toward 1. If that test
  is ever repeated it needs a control.

### Two rounds that touched no algorithm, 2026-08-13 and 2026-08-14

Neither changed a permutation, a factor or a number in any report; both are recorded in full
elsewhere, and are here only so a later session can tell they finished rather than stalled.

**The idiom sweep.** Four non-conformant sites brought onto the 2026-08-12 initializer rule, plus
two pre-C++11 idioms found while there: a zeroing default constructor in eight classes, and
`std::vector<T>().swap(v)` in `UpdateMatrix::discard`. The multifrontal peak was measured identical
to the byte before and after. `PORTING_LEDGER.md` now carries two sections, C++11 with every site
named and C++17 as a menu with no sweep owed; the 2026-08-13 `DESIGN_DECISIONS` entry has the
reasoning and the two defects the sweep exposed, one of which was in `CODING_RULES.md`'s own
example, which read a moved-from parameter and would have set an nnz to zero silently.

**One behavior change came out of it.** `Vector`'s two-argument constructor took a size and a
vector and checked no agreement between them, so `Vector(10, std::vector<double>(3))` built an
object reporting size 10 over three values. It is now `Vector(std::vector<Val> val)`, size derived,
so the state is unrepresentable rather than merely detectable. Nothing called the old form.

**The Makefile consolidation, and note it went the OTHER WAY from how the previous note framed
it.** That note asked for `benchmarks/ordering` to build on a bare `make`. Instead all nine
Makefiles now PRINT their target list on a bare `make`, with `.DEFAULT_GOAL := help` stated
explicitly rather than falling out of target order, `all`, `clean` and `help` everywhere, and
`test` only where something can fail, which is why the three benchmark directories do not have one.
`CLAUDE.md`'s clone checklist says `make all` in the two benchmark directories, so it no longer
passes vacuously. The convention is a rule in `CODING_RULES.md`. Two working targets turned up that
no help block advertised, `profile` and `example`.

**And the ordering ladder, 2026-08-14.** `degrees[pivot] = 0` moved before the refresh in md4 and
mdm2, so the convention holds in all six layers with a degree cache; mdm2's counter moved after the
loop it counts and its refresh set is named from `C[pivot]` as md4 and md5 name theirs. All four
traces byte-identical to before the move. `experiments/ordering/README.md` gained a section on how
the buckets differ from a bucket priority queue, and the external-degree section now gives the
paper's exact relation and the plain reading: `d_i` is the update size per SUPERNODE, `t_i` the
update size per COLUMN, with the front size as the other dimension.

**And three commits on the ordering experiments, 2026-08-14**, `aa7de25`, `49a673` and `20d6480`.
No permutation and no fill figure moved; the ONE output change is a `prepass:` field added to
mmd2's and mmd3's summary line, counting the vertices the prepass numbered, which until now was
visible only in the trace and so invisible on every grid. `experiments/ordering/README.md` carries
the substance and its mmd2 section roughly tripled; what a later session needs from here is:

- **A real bug, and it predated the day.** `mmd2` and `mmd3` failed at n = 1 in both twins, an
  `IndexError` in Python and an ASan report in C++, and at n = 0 a line later. Pass 1's floor files
  a degree-0 vertex under bucket 1 and the buckets were sized n. Now `n + 1`, matching production's
  `mHead(size + 1)` beside `mNext(size)`, plus an n = 0 early return matching `orderMmd2`'s. The
  prototypes' `Buckets` still sizes one degree-indexed array and three vertex-indexed ones from a
  single argument, which is where the fix belongs and is not where it was made.
- **Display work moved inside the `SHOW_THRESHOLD` guard** in md1 through md5 and the three mmd
  layers. It was being built and thrown away on every grid run; md1's Python was already guarded
  and its C++ was not, so the twins disagreed in a way no trace could show. Roughly a third of a
  C++ `grid 200` run in md5.
- **mmd1 aligned to md5** in four ways: the counter beside `num_clique_entries`, the pivot's unfile
  moved from the pick down beside the merged vertices, that block unfused into md5's two loops, and
  `degree` now `len(neighbors)` inside the display guard as md5 has it. mmd2 and mmd3 have had
  neither of the middle two.
- **`q2h` and `qxh` renamed** `two_source_queue` and `many_source_queue`. The vendored spellings
  are `q2h` and `qxh` exactly, `private/Mmd.cpp` line 118, and they are kept wherever genmmd is
  being cited.
- **`ncsub` was documented as impossible to check and is merely unchecked.** `mmd_order` takes it
  as a local and drops it; one temporary line in the wrapper prints it, as `tools/hook_amd.py`
  already does on the AMD side. Corrected in the README and both C++ twins.

### Earlier, 2026-08-09


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

**The numbering is historical and out of order, deliberately.** Items keep their numbers as they
close so that a reference from a commit message or another document still resolves. Read the "Read
this first" list above for what is actually open; below, an item is live only if its heading does
not say DONE, CLOSED, DEFERRED or PARKED.

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
- **The key fusion landed** and is worth about 4 to 7 percent in 2D across six sizes and 5 to 14 on
  cubes, both approximate for the harness reason below. The 2026-08-08 version of it failed because
  it stored keys in an array of size n, not because of the fusion.
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
Over eight runs `AMD3` reads 0.83 to 0.89 ms at 16 cubed where `AMD` reads 0.74 to 0.86, so the two
overlap there, and it is about 1.2x at 32 a side; 2D improved 0 to 8 percent as well. `AMD3.md`
iteration 26 and `docs/DESIGN_DECISIONS.md` (2026-08-10, "the algorithm was the smaller half")
carry it, with the corrections below.

**What is open, in the order we would take it.**

- **THE GROWTH TERM, which is now the whole question**, and item 2d below has the pass that would
  attack it. Both families are a low constant plus
  something that scales: cubes 0.97x at 16 a side rising about 25 percent to 1.2x at 32, 2D 1.33x
  at 32 a side rising about 45 percent to 1.9x at 400 with a knee past 280. The constant is nearly
  spent, and NOTHING measured so far touches the growth. The knee says memory. Two profiles of
  `amd3` compared AGAINST EACH OTHER rather than against the vendored routine, 140 and 400 a side
  in 2D, or 16 and 32 cubed, would say whether it is instructions or stalls. `MMD3` over the same
  quotient graph shows no growth on either family, which is the control that makes this an
  amd-branch property.
- **QUOTE ABSOLUTE TIMES, NOT RATIOS PER ROW.** The vendored column is the noisiest in the table:
  at 16 cubed it spans 16 percent across runs where `AMD3` spans 7, so `AMD3 / AMD` mostly measures
  `AMD`. And a free-function column is timed by `orderTimeFn`, a standing method by `orderTime`
  which also builds a Permutation, a difference of up to 2.4 percent, so those two are not
  comparable. Both cost a wrong claim on 2026-08-10.
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

**2a-hash. WE DO NOT COMPUTE WHAT `Amd.cpp` COMPUTES. FIVE DIVERGENCES, FOUND 2026-08-12, NONE FIXED.**

This is in the mechanism that cost several sessions in 2026-08-09, where a stride and a modulus
annihilated each other and turned the key into a function of the adjacency alone. That defect was
found and fixed. **What was not done then, and should have been, is a line-by-line comparison of the
whole key computation against the vendored one.** Doing that comparison now, while looking at
something else entirely, turns up four places where we compute a different number from `AMD_2`.

Each of the four preserves the EQUIVALENCE CLASSES, since two indistinguishable vertices have
identical lists and therefore identical keys under any of these schemes. So the merges found are
the same, which is exactly why every check we have is blind to all four: permutations match, the
twins agree, the fill columns are identical. **What differs is the BUCKET ASSIGNMENT**, which
decides the collision rate, and the collision rate is what the 2026-08-09 defect showed can cost a
factor of three.

**Divergence 1: the term. We add `v + 1` and `c + 1`; `Amd.cpp` adds `e` and `j`.**

```
Amd.cpp:   hval += e ;                  hval += j ;
ours:      key = (key + v + 1) % mod;   key += c + 1;
```

Our key is therefore the vendored one plus the NUMBER OF TERMS. The `+ 1` is a leftover: it was
part of the `(c + 1) * (size + 1)` stride removed on 2026-08-09, and only the multiplication was
taken out. Nothing now requires it.

**Divergence 2: we include the pivot; `Amd.cpp` does not include `me`.** Our loop takes every entry
of `I[u]` with the pivot among them. In `Amd.cpp`, `hval` is accumulated in scan 2 over `p1..p2`,
the OLD element list, and `Iw [p1] = me` happens AFTERWARDS, in the "add me to the list for i"
block. The new element is not in the sum.

**And the comment at `src/Amd3.cpp` justifying this says the opposite**: "its scan 2 accumulates
`hval += e` over the element list it is compacting, which by then holds the new element". The
vendored source does not support that reading. **Either the comment is wrong or my reading is; that
has to be settled before anything is changed**, and it is the first thing to do on this item.

**Divergence 3: the modulus. Ours is `n + 1`, `Amd.cpp`'s is `n`.**

```
Amd.cpp:   hval = hval % n ;
ours:      hash = key % (size + 1);
```

The `n + 1` is the other half of the removed stride: it was chosen so the modulus would not divide
the stride, and with no stride left there is nothing for it to be coprime to. It is also the
`n + 1` that overflows `int32_t` at `MAX_IDX` and holds n prisoner, so this divergence and the type
question are the same question.

**Divergence 4: overflow. `Amd.cpp` WRAPS at 2^32 and we never do.** Its `hval` is a `UInt`
accumulated unreduced over both lists, so a long list wraps before the reduction; the header
comment says the type is unsigned so that `%` is defined, not to prevent wrap. We reduce as we
accumulate in the prune and sum wide in the driver, so we compute the exact sum. When the true sum
exceeds 2^32 the two land in different buckets, because `2^32 mod n` is not 0. Reachable only on a
long list at large n, but it is a real difference in the function.

**None of this is measured.** The right first step is a collision counter, pairs tested per pivot
and the bucket size distribution, on 2D and 3D at comparable n, against the vendored routine's
0.33 and 0.48 recorded in the `Amd3.cpp` comment. That number is the one that moved by a factor of
three last time, and it is the only output any of these four can change.

**A FIFTH DIVERGENCE, in a different mechanism, found 2026-08-12 while reading the twins.** The
mass-elimination test in `QuotientGraph::massEliminate` has three conjuncts where `AMD_2` has two:

```
Amd.cpp:   if (Elen [i] == 1 && p3 == pn)
ours:      if (mAdjacencySize[u] == 0 && mIncidenceSize[u] == 1 &&
               mSource[mSourcePtr[u]] == pivot)
```

`AMD_2` asks for one element left and no surviving variables, and never checks WHICH element.
Ours adds that check, and **it is provably redundant**: all three prune variants append the pivot
unconditionally, at `src/QuotientGraph.cpp` lines 395, 460 and 533, each followed immediately by
`mIncidenceSize[u] = write - kept`, so a count of one forces the sole entry to be the pivot. The
`mVendoredListOrder` swap is guarded by `write - kept > 1` and cannot move anything at a count of
one, and the first conjunct is what puts the incidence run at `mSourcePtr[u]`.

**Unlike the four above, this one cannot change any output.** Being redundant it always passes, so
permutations, fill and merges are identical with or without it. It is pure cost: `&&` short
circuits, so it runs once per MERGED vertex rather than per candidate, and costs two dependent
loads into an arena the prune has walked past since writing it. On grids, where mass elimination
fires often, that is a modest number of avoidable misses in a hot region. Removing it is free in
every sense; it is listed here rather than done so that it is decided deliberately.

**AND THE HASH'S STORAGE DIVERGES TOO, not only its arithmetic.** Found 2026-08-12 while reading
`Buckets`. `AMD_2` allocates NO arrays for supervariable detection: `Head [hval]` doubles as the
hash bucket head when the degree list there is empty, distinguished by the `FLIP` encoding, and
`Last [i]` doubles as the stored key. Each amd driver here allocates `hashHead(size + 1)` and
`hashNext(size)`, so **about 2n extra `int32`**, roughly 1.3 MB at 400 a side, touched about twice
per survivor per elimination.

**A lead rather than a plan, and it has a reason it may be declined.** Production already takes
half of `Amd.cpp`'s trick, the running key riding in `hashNext` rather than a third array of its
own. The other half, overlaying the hash heads on the degree heads, is the `FLIP` encoding that
`docs/DESIGN_DECISIONS.md` calls an anti-model. What would settle it is a measurement rather than
an argument: the 2026-08-08 note records two extra arrays of size n costing 12 percent in 2D at
400 a side across the phase boundary, which is the same order of footprint in a different place.

**Two smaller differences in `Buckets`, both examined and both left alone.** `AMD_2` removes a
vertex from its degree list with no guard, needing none because `Nv [i] = -nvi` flags `i` the
moment it enters `Lme` and nothing revisits it; our `unfile` carries an `UNFILED` early return
for the MMD branch, where a batch can evict a vertex a later pivot in the same round then merges
away. And `AMD_2` leaves `Next [i]` and `Last [i]` stale after a removal where we clear them, two
stores per unfile and 224054 unfiles on a 140x140 grid. The clear is not droppable: `file` always
rewrites `mNext[u]`, but MMD's `next(u)` walk would follow a stale link if it ever reached an
unfiled vertex, which is the same batch hazard the guard exists for.

**Why this is item 1 and not item 6.** Item 2d proposes profiling `amd3` against itself to decide
whether the descriptor struct is worth building. That profile measures a routine whose hash is not
the one we are trying to match, and the hash sits in the loop the profile is aimed at. **Matching
first, then measuring, then getting clever, in that order**, which is the rule the 2026-08-09
episode was supposed to have taught and evidently only half taught.

**The general lesson, since it is the second time.** Reinvention concentrates where outputs cannot
see the difference. A hash is the purest case of that in the whole ordering: it is a filter and
never the decision, so ANY key that respects the equivalence classes produces identical output and
differs only in speed. **That makes it the one place where "we match the vendored routine" has to be
checked by reading, because no test will ever check it for us.** The same is true of any future
filter, guard, or early-out.

---

**2a-twins. THE PROTOTYPE SWEEP: TEN LAYERS OF FIFTEEN DONE, 2026-08-12. `mda2`, `mdam2` and
`amd1`-`amd4` REMAIN.** `experiments/ordering` was brought onto the current integer rule and given
four proper classes. Done: `md1`-`md5`, `mdm2`, `mmd1`, `mmd2`, `mmd3`. Nothing in it changes
behavior; every layer's output is identical to its pre-change self except for one word.

**The recipe, in the order it works, so the remaining five go the same way.**

1. **`using Graph = ...` becomes two classes**, `AdjacencyGraph` and `IncidenceGraph`. Identical
   bodies, deliberately not shared: both hold one int32 list per vertex and differ only in what
   the entries MEAN, `A[u]` vertices and `I[u]` clique ids. An alias made them one type, so
   nothing but a variable name said which was which. `AdjacencyGraph` needs a second constructor
   taking `std::initializer_list<std::vector<std::int32_t>>`, for the brace-written examples.
2. **`Cliques` and `Buckets` become classes**, `public:` first and members `m`-prefixed in a
   private section, which is production's layout and Core Guidelines NL.16. Members ordered
   scalars first, since declaration order is initialization order; accessors follow the members
   one for one and in the same order.
3. **`n` leaves `std::size_t`.** `A.size()` returning `std::uint32_t` is what carries it: the
   drivers' `n`, `SHOW_THRESHOLD`, the counters, the `pivots` and `order` loops all follow, and
   the `Storage` function loses its range-`for` because the classes are indexed rather than
   iterable.
4. **`alive` becomes `live`**, in identifiers, printed text and prose, in BOTH twins together. A
   bare `alive` list takes the suffix, `live_vertices`.
5. **Comments last**: the reach-tag note, the eliminator's tag paragraph, the merged-tag note, the
   scratch-buffer note.

**What the verification has to be, because two of these checks are blind to different things.**
Compile `-Wall -Wextra`; diff the C++ against the Python twin in examples AND grid mode; diff the
C++ against its own pre-change output with `alive -> live` applied to the baseline, which proves
the word is the ONLY change; and sanitize. Then diff the shared functions against the layer they
came from, which is the check that caught `md3.py` keeping `candidate`/`best` while its own twin
comparison passed.

**Three traps, all hit at least once.**

- **The mechanical `cliqueTag -> pivotCliqueTag` substitution breaks a trailing comment's column.**
  It happened in `mdm2`, `md4`, `md5` and `mmd1`, every time, and only a cross-layer diff finds it.
- **Signatures with extra parameters escape an exact-match replacement.** `mmd1MinimumDegree` takes
  `std::int32_t delta = 0` and its `Graph` survived, failing the build with `'Graph' does not name
  a type`. The amd layers will have their own.
- **`std::max<std::size_t>` and function return types are not caught by any of the above.**
  `mmd2Degree` and `mmd3Degree` returned `std::size_t` and needed narrowing by hand, and one
  `std::max` needed a cast on `superMembers[u].size()` to form the expression in `uint32_t`.

**What the remaining five will need beyond the recipe.** `mda2` and `mdam2` are md-family and
should port like `mdm2`. **The four amd layers will not**: they carry the approximate bound,
aggressive absorption and the hash, so their `Neighbors` and `Eliminate` will not match the mmd
family and want reading rather than substitution. `amd2` onward also carry `hashHead`, `hashNext`
and `usedKeys`, which is the machinery item 2a-hash is about, so the sweep and that item touch the
same lines.

---

**2b. Taking `AMD1B` and `AMD2B` out of the public enum: DONE, 2026-08-15.** `Ordering` is nine
enumerators, `Natural, MMD, MMD1, MMD2, MMD3, AMD, AMD1, AMD2, AMD3`, and the two B layers are
reached as free functions exactly as item 2c specified and as `AMDraw` and `MMD3B` already were.
`Mmd1` and `Amd1` stay in the enum: they are ladder rungs, but they are also complete, correct
orderings a caller may legitimately choose.

**What moved with it**, since the surface is wider than the enum: `OrderEngine`'s dispatch and its
two private methods; the ordering switch and the sweep list in three files under `examples/`, which
is the site item 2c warns about; `test_pipeline`'s sweep, from nine orderings to seven, with its
written-out expectation corrected; `benchmarks/ordering/order_timing.cpp`, which grew a general
free-function column beside the `mmd3b` one; `benchmarks/pipeline/pipeline_timing.cpp`;
`README.md`'s Structure block and ordering paragraph; and `docs/TESTING_SPECIFICATION.md`, which
moves with the suite by invariant.

**And `test_order` keeps all fourteen B assertions**, calling the free functions instead. That was
the one thing that had to survive: the pair have no prototype, are not in `experiments/ordering`'s
`PORTED` list, and dropping out of `test_pipeline`'s sweep leaves this as the ONLY check on them
anywhere.

**The original entry follows, since it is what specified the shape.**

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

**2d. THE NEXT PERFORMANCE PASS, and what of 2026-08-10 the other layers never got.** The walk
axis is FINISHED on `Amd3`: after the fold it walks `I[u]` twice per pivot and `A[u]` once, which
is `AMD_2`'s count exactly, so there is no pass left to remove that the vendored routine does not
also make. What is left is seven sweeps over `C[p]` against its four, and the growth term. Three
items, in the order the evidence supports.

**(i) One profile first, and it decides whether (ii) is the right tool.** CPU Counters on `amd3` at
140 and at 400 a side, compared AGAINST EACH OTHER rather than against the vendored routine. The
constant factor is nearly spent; what remains grows with n and has a knee past 280 a side, which
says memory. If the growth is stalls, (ii) is aimed at it. If it is instructions, this direction is
finished and the descriptor struct is the wrong tool.

**(ii) The descriptor struct, `docs/TODO.md` question 3.** `mSourcePtr[u]`, `mAdjacencySize[u]` and
`mIncidenceSize[u]` are always read together for the same `u` and live in three arrays, so a
vertex's descriptor touches three cache lines where one struct touches one. Each of the seven
sweeps opens by reading exactly that, so this attacks all seven at once rather than deleting one,
and deleting one is what B1 tried and measured nothing. It is also the only candidate left that
targets STREAMS rather than passes, which is what paid three times on 2026-08-10. It lives in
`QuotientGraph`, so all six drivers move together and no vehicle can isolate it: `make amdorder`
and `make mmdorder` are the guard. The known ceiling is small, about 7 percent of a one-shot solve
for closing `AMD1`'s whole remaining gap, and this is a fraction of that. Narrowing the arrays was
tried once and measured nothing despite cutting simulated misses 17 percent.

**(iii) What 2026-08-10 left on the table for the OTHER layers**, which is the part nobody has
looked at:

- **`Amd1B` and `Amd2B` carry `explicitPart(size)` as an extra array**, which is exactly the
  footprint cost removed from `Amd3` that day. `Amd1B` is on record as "slower at large n after
  being faster at small", and that is the signature of precisely this. **But the fix may not
  transfer.** `Amd3` had somewhere to put the value, `partial[u]`, which exists only because
  ledger entry 4 split its bound in two; `Amd1B` and `Amd2B` form the bound in one pass and have
  no dead size-n array at that moment. `Amd2B` has `hashNext`, dead until filing, and `Amd1B`
  appears to have nothing. **Entry 4's split is now the third thing it has bought by accident**,
  after `Amd3`'s immunity to the hash-bucket order and its ability to take the key fusion at all.
- **The tagged W consolidation never propagated.** `Amd1`, `Amd1B`, `Amd2` and `Amd2B` all carry
  `outside(size)` plus a mark plus a clearing pass, `for (c : touchedCliques) outside[c] = 0`,
  where `Amd3` carries one `w` and invalidates the lot with a single addition. That is one size-n
  array and one pass per pivot, in four drivers, and it is iteration 15's change which landed in
  `Amd3` alone as `AMD3C`. Care needed in `Amd2` and `Amd2B`: their `mark` is 2n and serves the
  hash stamp as well, so only `outside` folds away.
- **The key fusion is NOT available to `Amd2` and `Amd2B`**, checked on 2026-08-10 and recorded
  under item 2b's neighbors: their key pass walks `C[p]` backward against a forward bound that
  also refiles, so head insertion into both structures wants opposite directions. `Amd1` and
  `Amd1B` have no hash and so no key.

**None of (iii) is production work.** `Amd1` and `Amd2` are ladder rungs and the B layers are
oracles; their speed has no consumer. It is listed because it is evidence about the SHARED class:
if folding `explicitPart` away is worth something in `Amd2B`, that is a second data point for the
stream hypothesis at no risk to the default.

**2e. NARROW THE ONE-DIMENSIONAL SIZES. THE ORDERING IS COMPLETE, 2026-08-11.** Eleven steps, each
with the permutation checked against the pre-change tree on 2D grids to 80 a side, 3D to 16 and
three random patterns at n = 2000, all identical, and each clean under `-Wall -Wextra` and under
`-fsanitize=address,undefined`. What landed: the four `QuotientGraph` arrays and the accessors over
them, `Buckets`'s signatures, `degrees` with the scalars that travel with it, the four scan arrays
with the two structs that bind them, then `usedKeys` and its hash locals, `sizeU` and `sizeV`,
`reachableSize`, `absorb`'s `vertexCount`, and one entity loop in `Amd3` brought back to the
`int32_t` form. `docs/DESIGN_DECISIONS.md` (2026-08-11) carries the account and
`docs/CODING_RULES.md` now states the three-way rule.

**Three things came out of it that the list below did not anticipate**, all in the design entry:
disjointness rather than arity is what bounds an accumulation, so nine sites needed no widening; a
narrowing cast goes outside an expression and a widening cast must go on an operand, which removed
the last arithmetic depending on the cap on n; and `bound` in four drivers was narrowed and had to
be reverted, being a fixed sum at its declaration and an accumulator three lines later.

**What remains `std::size_t` in the ordering is correct, and the list is short enough to check
against.** The `colPtr` parameters and the `cp` loop over them; `mSourcePtr` and `mCliquePtr`; the
five accumulators, `bound` in the four accumulating amd drivers, `deg` in `Amd3` and the hash `key`
in three, together with their `std::min<std::size_t>` calls and operand casts; `size()`;
`Buckets`'s constructor; the driver-local `size`; and `numFlagSweeps`, a diagnostic counter bounded
by the tag range rather than by n, which the dimension rule does not decide and which was left
deliberately.

**What is left of this item is everything outside the ordering.** The symbolic and numeric phases
are untouched and are the larger half. **And two casts hold n prisoner**, neither reachable and
both left deliberately:
`Amd3`'s `modulus = static_cast<std::int32_t>(size + 1)` fails at exactly the largest n
`checkIndexRange` admits, and `mark[incidenceU[i] + cliqueStamp]` in the three hash drivers reaches
`2n - 1` in signed `int32_t` and is undefined above `n = 2^30`. The design entry says why they were
not fixed with the rest.

**The original text of the item follows, kept because the target list is what the work was driven
from.** The integer model
becomes: **`std::uint32_t` for one-dimensional sizes, `std::size_t` for two-dimensional ones, and
`std::int32_t` for indices and for anything that carries NIL.** This reverses the 2026-08-08
`DESIGN_DECISIONS` entry, which kept the wide one-dimensional types as a considered trade; that
entry named the count sweep as the cheap half and this is it. **Ordering first**, since it is
self-contained and has the shared class in it.

**The operative test, and it is about the type's range rather than about n.** Cast when the result
can exceed what the narrow type holds:

- a PRODUCT of two 1D quantities always can, `f * u` reaching n^2;
- an ACCUMULATION over an unbounded number of 1D terms always can;
- a sum or difference of a FIXED FEW 1D quantities never can, and `frontSize + updateSize` is not
  even that, a supernode's columns and its update rows being disjoint subsets of the same n
  indices.

### The targets, enumerated, 2026-08-11

Taken from a reading of the ordering sources on that date, so an item is either done or not and
nothing has to be re-derived per session. Four buckets. **The split came out cleaner than expected:
nothing currently `std::int32_t` wants to become `std::uint32_t`, so the narrowing set is exactly
the currently-`std::size_t` set minus the two `Ptr` arrays, minus the two accumulators.**

**Bucket 1, stays `std::size_t`, being two-dimensional.** Two arrays, both in `QuotientGraph`, and
both positions rather than counts:

- `mSourcePtr[u]`, where u's run starts in `mSource`; bounded by nnz(A).
- `mCliquePtr[c]`, where c's block starts in `mCliqueArena`; bounded by the sum of every clique
  ever formed, which grows toward nnz(L).

Both would FIT in 32 bits at any size we can factor. They stay wide because of what they mean,
which is the rule doing its job rather than the rule being slack.

**Bucket 2, stays `std::int32_t`, carrying NIL or a sign.** Larger than it looks, and no work:

- Index payloads: `mSource`, `mCliqueArena`, and the drivers' id lists, `pivots`,
  `touchedCliques`, `deadCliques`, `touched`, `batch`, `elementMembers`, `q2h`, `qxh`,
  `refreshed`, `mMerged`.
- Links: `mSuperNext`, `mSuperLast`, `hashHead`, `hashNext`, and `Buckets`'s `mHead`, `mNext`,
  `mPrev`, which carries two sentinels, `NIL` and `UNFILED = -2`.
- Stamps: `mMark` and `mTag`, the drivers' `mark` and `tag`, `touchedRound`.
- `Amd3`'s `w`, `wflg`, `lemax`, `wbig`, signed by requirement rather than by sentinel; see below.

`mEliminated` and `Mmd2`/`Mmd3`'s `outmatched` are `std::uint8_t` and stay.

**Bucket 3, narrows to `std::uint32_t`.** In `QuotientGraph`, four arrays and one scalar, all
bounded by n: `mAdjacencySize`, `mIncidenceSize`, `mCliqueSize`, `mWeight`, and the per-pivot
`mCliqueWeight`. The accessors over them move too, `adjacencySize`, `incidenceSize`, `cliqueSize`,
`weight`, `cliqueWeight`, `reachableSize` and `reachableWeight`, which is what takes the casts out
of the hot loops.

In the drivers, **twenty-three declarations across the eight files**:

| array | where | what it holds |
|---|---|---|
| `degrees` | all eight | the filed degree, at most n |
| `outside` | `Amd1`, `Amd1B`, `Amd2`, `Amd2B` | per clique, `\|C[c] - C[p]\|` weighted |
| `cliqueDegree` | the five amd drivers other than `Amd3`, plus `Amd3` | per clique, `\|C[c]\|` weighted |
| `explicitPart` | `Amd1B`, `Amd2B` | per vertex, weight summed over the pruned `A[u]` |
| `partial` | `Amd3` | the half-formed bound, ledger entry 4 |
| `usedKeys` | `Amd2`, `Amd2B`, `Amd3` | hash values in [0, n]; a count by the range test, though it does not read as one |

`Amd3` has no `outside`, carrying `Amd.cpp`'s tagged `w` in its place, which is why the amd column
is not uniform. The two scan structs hold references to these same arrays, three in
`ApproximateScan` and two in `TaggedScan`, and move with them rather than being separate targets.

And the scalars that move with the arrays: `numEliminated`, `numLive`, `minDegree`, `degme`,
`numLeft`, `batchLimit`, `dg0`, `weightU`, `weightV`, `sizeU`, `sizeV`, `otherSources`, and the
hoisted lengths `pivotCliqueSize`, `adjacencySize`, `incidenceSize`, `membersSize`, `cliqueSize`.

**Bucket 4, must stay wide.** The two accumulators below, which the operative test catches and a
per-array reading does not.

### Three decisions to take before the first edit

Each is a decision rather than a consequence, and each reaches every file, so taking them in order
is what keeps the change mechanical afterwards.

1. **Does `n` itself narrow?** `QuotientGraph::size()` returns `mAdjacencySize.size()`, a container
   size, so it is `std::size_t` whatever the element type becomes. If the accessors return
   `std::uint32_t` while `size()` does not, every comparison between a count and n is a
   mixed-width one. The likely answer is that `size()` returns `std::uint32_t` with the cast in
   that one place, but it should be chosen rather than fallen into.
2. **What `Buckets` becomes.** Its three arrays are `std::int32_t` links and do not change. What
   changes is the signatures: `file`, `unfile`, `head` and `empty` taking a `std::uint32_t`
   degree, and `refile` taking `std::vector<std::uint32_t>&`. Note `mHead` is sized n + 1 because
   MMD files a vertex under its degree plus one, so the parameter's range reaches n exactly.
3. **Whether the eight subtraction sites get anything beyond the note below.** They are correct by
   invariant today and nothing enforces it; narrowing changes the wrap distance and not the
   hazard.

**The NIL-carrying set needs no work.** `mMark`, `mTag`, `mHead`, `mNext`, `mPrev`, `mSuperNext`,
`mSuperLast`, `hashHead`, `hashNext` and `touchedCliques` are already `std::int32_t`.

**One quantity is deliberately SIGNED and must stay so**, and it is the only one. `Amd3`'s tagged W
array with `wflg`, `wnvi`, `lemax` and `wbig`: `wnvi = wflg - weight(u)` is negative whenever the
tag is still small and the weight is not, which is the first few eliminations, and the arithmetic
only comes right at `w[c] - wflg`. Unsigned wraps there and the bound comes out enormous. `Amd.cpp`
is signed for the same reason. `src/Amd3.cpp` says this at its declaration; do not narrow it and do
not make it unsigned.

**Two accumulators cross and must stay wide**, both in the bound pass. `deg` accumulates
`w[c] - wflg` over every clique in `I[u]`, each term up to n and O(n) of them, so the intermediate
reaches O(n^2); it is at most n only because the two caps are applied afterwards, the bound's whole
purpose being to overcount. And the hash `key` sums `c + 1` over `A[u]` and `I[u]`, same shape. The
`TaggedScan` path already avoids the second by reducing modulo `n + 1` as it accumulates, which is
why it fits an `int32_t` there; the driver's half does not.

**Eight subtraction sites are non-negative by INVARIANT rather than by type**, and none is
enforced: `outside[c] = cliqueDegree[c] - weightU` and its `-=`, `bound = explicitPart + degme -
weight(u)`, `partial[u] + degme - weightU`, `numLeft - weight(u)`, `degrees[u] + degme -
weight(u)`, and `degrees[u] - weightV` in `Amd2`'s hash merge. Each holds because the subtrahend is
part of the minuend. They already wrap silently under `std::size_t`, so `uint32_t` inherits the
hazard rather than creating it; what changes is the wrap distance, about 4.3e9 instead of 1.8e19,
which reads as a large degree rather than an absurd one. Worth knowing when one of them ever fires.

**Expect the model and the deleted casts, not milliseconds.** Narrowing six of these arrays was
tried on 2026-08-01: 17 percent fewer simulated D1 misses and zero on alpamayo. What it does buy is
that the hot loops stop casting: `const std::int32_t adjacencySize = static_cast<std::int32_t>(
mAdjacencySize[u])` and its siblings exist only because the array is wider than the loop. **Eight
of the eleven casts in `src/QuotientGraph.cpp` are that shape**, seven hoisting an adjacency or
incidence length and one narrowing `cliqueDegree[c]` to add it to `wnvi`, so the count is a
before-and-after worth recording rather than an estimate. And it composes with item 2d:
`{ size_t sourcePtr; uint32 adjacencySize; uint32 incidenceSize; }` is 16
bytes, four to a cache line, which makes the descriptor struct the natural next step rather than a
separate argument.

**The rule itself is now stated in its derived form**, in `docs/CODING_RULES.md` and in the
2026-08-12 design entry: a size is always unsigned, one dimensional in 32 bits and two dimensional
in 64; an index is signed 32 bits and only because `NIL` has to share a type with the values it
stands in for. Anything that can exceed n is computed in `std::size_t`, not because 2n overflows,
which it does not, but because its safety would otherwise be inherited from a cap enforced
elsewhere for an unrelated reason. **The rest of the tree is drift against that rule, to be closed
incrementally**: n itself, the symbolic phase, the numeric phase.

**The durable half has to be written before this item closes, not after.** This file is meant to be
deleted, so a list worked through here leaves nothing behind unless two other files move with it:
`docs/CODING_RULES.md`, whose index-type section states a TWO-way split, index against position,
where this is a THREE-way one; and `docs/DESIGN_DECISIONS.md`, since narrowing reverses half of the
2026-08-08 integer-model entry, which kept the wide one-dimensional types as a considered trade and
named the count sweep as the cheap half. Nothing detects drift between the three, there being no
invariant here of the kind that binds the suite to `docs/TESTING_SPECIFICATION.md`, so the
`CODING_RULES` edit belongs in the same step as the last source file rather than as a follow-up.

**3. The same widening is available to `make test`, cheaply.** Its prototype-against-production
comparison still runs on 2D grids at sides 10 and 20 alone, and `graphs.h` holds the 3D and random
builders. Worth knowing before relying on that check: the prototypes carry no maintained clique
degree, so a defect in production's encoding is invisible to it at any size. See `docs/TODO.md`,
the first of the five ordering questions.

**4. And only then the performance work below, which is a PARK rather than a queue.** One bounded
attempt, worth an hour whenever the amd branch is picked up again. Nothing here blocks anything.

**Re-priced 2026-08-10, and it is now the WEAKER of the three candidates rather than the stronger;
item 2d has the ranking.**
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
the prune, `AMD3` runs at roughly 1.0 to 1.25x the vendored routine on cubic grids from 12 to 32 a
side, overlapping it at 16, and 1.33 to about 1.95x on square grids from 32 to 400. `MMD3` is at
0.87 to 1.05x its own on cubes and 1.26 to 1.55x in 2D, with no trend in n on either. `AMD3` is
FASTER than both `AMD1` and `AMD2` at every cubic size while
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
- **`DirectSolver` does not forward `NumFactorDynamic::rank()`.** NEW on 2026-08-11. It forwards
  the perturbation and delay counts, the pivot counts and `inertia`, but not the numerical rank,
  which is the one number that names a numerically singular matrix directly. The accuracy benchmark
  stands in with the inertia's zero count, whose own header warns it is least reliable on singular
  input, which is exactly where it is being used. A small library change.
- **What dynamic pivoting costs, which is now two questions with two kinds of evidence.** NEW on
  2026-08-11 and the largest thing this tree learned about the library rather than about a
  benchmark.

  **The search costs even when it finds nothing**, which the scaling ladder measured cleanly
  because its matrices are built to need no pivoting: 1.48x static in 2D and **6.98x in 3D**, with
  zero columns delayed, growing as nnz(L)^1.31 against static's nnz(L)^1.01. Charged per
  front-column against a growing front. That is a cost anyone paying it unnecessarily should be
  able to avoid, and Oblio already lets them, by asking for Cholesky or static LDL.

  **And the delays themselves cascade on a chain**, which the real matrices showed:
  `Oberwolfach/LFAT5000` delays 3.1 million columns for a 232-fold increase in fill on a matrix
  that is POSITIVE DEFINITE, where no pivot needed delaying at all; `GHS_indef/bloweya` does the
  same past 32 GB and cannot be factored. `GHS_indef/bloweybq` is the control: the same chain
  shape, definite, zero delays.

  Ashcraft, Grimes and Lewis name the mechanism and give four routes out, all four recorded in
  `benchmarks/matrices/ACCURACY.md`. The first is free: relax `setPivotThreshold` from its default
  of 0.1, which they study directly and recommend loosening. **It would test both halves at once**,
  since a looser threshold accepts more pivots and so both delays fewer columns and, if the search
  short-circuits on acceptance, searches less.
- **Nested dissection, which Oblio does not have and which one matrix family badly wants.** NEW on
  2026-08-11. `PARSEC/Si87H76`, an electronic structure Hamiltonian at n = 240369 with 10.7 million
  nonzeros, predicts **5.68 billion entries of fill under MMD3**, which is 42.3 GB of values and
  beyond the machine. AMD3 will not rescue it: the two agree to within a few percent everywhere we
  have measured, so an approximation of a metric will not save a matrix the metric itself handles
  badly. What that family wants is a partitioning ordering rather than a greedy one, which is what
  MUMPS and CHOLMOD reach for on exactly these matrices.

  **Worth knowing that every ordering result in the three reports is minimum degree against minimum
  degree**, two variants of one idea, agreeing to within a few percent nearly everywhere. Nested
  dissection is the first thing that would say whether that agreement is because both are good or
  because both share a blind spot, and `PARSEC` is the sharpest available test of it. METIS is
  source-available and would slot in beside the vendored pair through the existing `Ordering` enum
  and `OrderEngine` dispatch, though it is a directory with its own build rather than two files,
  and its licence wants reading against the PolyForm plan. A simple separator-based ordering
  written here would answer the question less well and cost nothing external.

  The three matrices are fetched and sit in `data/`, named by `benchmarks/matrices/extras.txt`;
  the caps keep them out of every run.

- **Left open by the 2026-08-13 idiom sweep**, none of it blocking and all of it recorded in
  `PORTING_LEDGER.md`'s two idiom sections or the 2026-08-13 `DESIGN_DECISIONS` entry.

  - **Four stored lengths**, `SymFactor::mNumNodeIdx`, `NumFactorStatic::mNumNodeIdx` and
    `mNumVal`, and `UpdateMatrix::mSize`, each a member holding a length its container already
    knows. A DECISION NOT TAKEN rather than an oversight: the arguments run both ways, and two
    were removed during the sweep and put back the same hour because nobody had asked. Taking it
    would move the headers' own reasoning, `ARCHITECTURE.md`'s accessor table and
    `test_pipeline`'s seventeen size assertions along with it.
  - **Four small C++17 items**, each with its site named in the ledger. `[[nodiscard]]` is the
    only one that would catch a bug rather than tidy a line, on the query accessors where a
    discarded result is always a mistake. The others are a `_v` spelling, a `(void)` cast that
    wants `[[maybe_unused]]`, and a `const char*` that wants `std::string_view`.
  - **A `clang-tidy` run with the `modernize-*` checks** over `src/` and `include/`. The ledger
    says nothing is KNOWN to remain, which is weaker than nothing remains: the sweep followed an
    audit of four sites, not a systematic search. One command settles it either way.
  - **Em-dashes in three headers**, `ElmForest.h`, `ElmForestEngine.h` and `Permutation.h`, which
    `WRITING_RULES.md` calls a hard rule. `SparseMatrix.h`'s were cleared on 2026-08-13 because
    that file was being edited anyway.
- **Changing a matrix's VALUES while keeping its structure**, which is the refactorization case and
  what lets one `analyze` serve many `factor` calls. Replacing a matrix wholesale already works and
  needs no new code, assignment from a fresh one going through the constructor so the guards run;
  the values-only path does not exist. `TODO.md` has the entry, with the reasons a `setValues` is
  the right shape and a `setStructure` is not.
