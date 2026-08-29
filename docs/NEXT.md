# NEXT: the mmd and amd branches are being aligned against each other; one regression is open

## DONE 2026-08-28: the ordering code took the naming conventions, and three defects fell out

Two commits, `98af0bf` and `64020aa`, both named "Ordering code alignments". Most of it is naming
and comments; the code changes are small and every one of them was checked against a digest of all
five drivers over grids 1 to 40 plus six edge cases.

### What changed in the code

- `order` is `orderAsMerged`, and `QuotientGraphChained::order` is gone, declared and never defined.
  See `docs/DESIGN_DECISIONS.md` for why not `orderAmd`.
- `mSize` is a stored field in all three classes, and the clique arena's guard is a plain doubling
  on the invariant `capacity() >= mSize`. Same file.
- `orderAscending` lost its temporary: the decode is `-(cursor[u] + 1)`, which is the encode written
  again, `-(x + 1)` being its own inverse. `slot` is `cursor`, `pos` is `k`.
- Both prunes and both `massEliminate`s lost duplicate locals. `pruneMmd` runs on ONE cursor across
  both compactions, taking the boundary from the descriptor it just wrote. The compacted
  `massEliminate` calls `trimClique` where it had an eleven-line copy of that method's body.
- `formReachableSet*` and `reachableSetWeight` are the names in all three classes now, and the
  compacted class's five definitions moved to follow its constructor.
- Locals take the entity prefix and cliques take their role. See `docs/CODING_RULES.md`.

### Three stale comments, and none of them was findable by reading

**A comment claimed an ordering was load-bearing when it was not.** `orderAscending`'s first pass
said the members had to be stamped before the pivot was placed, or the stamp would overwrite the
cursor. It would not: the chain walk starts at `mSuperNext[pivot]`, so the pivot is never stamped,
and the two writes are disjoint. Settled by swapping the two halves and diffing 40 grids, which took
one build. The comment now states the invariant instead, and the pivot goes first because that reads
in the order the reader thinks.

**A comment described a swap that was in the other function.** `pruneMmd` carried the explanation of
the `[c1..ck, pivot] -> [pivot, c2..ck, c1]` rotation; the rotation is `pruneAmd`'s and mmd has
none. Moved to the swap, with a closing line saying mmd does not rotate, since the absence is what a
reader arriving from the other prune will ask about.

**And two comments in the amd drivers claimed the sign restore happened there.** It moved into the
prune in `5dcc57c` and those two were not swept with it; `AmdCompacted.cpp` had two copies, one a
leftover of the other, which is presumably how it survived. Found by checking whether
`massEliminate` sees negative weights, which it does not, on either branch. Deleted.

### Two process notes, both earned the hard way

**`-Wshadow` DOES NOT CATCH A LOCAL SHADOWING A MEMBER FUNCTION.** gcc warns for data members and
for locals, and says nothing when a local named `clique` or `cliqueSize` hides an accessor of that
name. `reachableSetWeight` had four such locals and two of them predate this session. The check
that works is a parse: pull every local out of the function body, pull every member name out of the
class, and intersect. That is how the five aligned functions were cleared.

**A WORD-BOUNDARY REGEX IS NOT A SCOPE, THREE TIMES IN ONE SESSION.** `\badjacencySize\b` renamed
the struct FIELD as well as the local, at nine sites, which did not compile and so surfaced at once.
`\bincidence\b` twice rewrote the English word inside comment prose, eleven lines across two driver
files, which compiles fine and is invisible to every gate. The repair that works is to re-run over
comment text only and then assert that no pure comment line differs from the committed version.
Rename inside a function's own text, not the file's.

### Where this leaves the front

The alignment work is close to done on `reachableSetWeight`, the two form functions, `merge`,
`massEliminate` and the prunes: the flat and compacted bodies of `reachableSetWeight` and
`massEliminate` are now identical modulo the arena name, checked by substitution rather than by eye.
What remains is comments rather than code. The compacted `massEliminate` is missing four comment
blocks flat has, and the chained `merge` has none of the two the others carry and writes
`mMark[v] = GONE` directly where they call `markGone`.

`pruneAmd` still holds `key`, which shadows `Buckets::key`, and it is the last shadowing local in
the file.

## PARTLY DONE 2026-08-28: the flat and compacted drivers could be ONE TEMPLATED BODY per branch

Compaction happens inside `QuotientGraphCompacted`, so a driver never sees it. Diffing the driver
bodies with comments stripped shows how little is left once that is true:

```
                     lines            substantive hunks in the shared body
MmdFlat  / MmdCompacted   149 / 136        NONE
AmdFlat  / AmdCompacted   169 / 162        NONE, since the run order was flipped
```

**THE FOURTH OPTION BELOW WAS DONE ON 2026-08-28 AND STEP 1 IS GONE.** `QuotientGraph` now carries
the run order per branch, `[A, I]` for mmd and `[I, A]` for amd, so `AmdFlat`'s detection is a
single loop and the two amd detection regions are BYTE-IDENTICAL. What is left of this entry is the
templating itself, in two steps rather than three.

The precedent is in this tree already: `NumFactorEngine`'s traversals are
`template<class Val, class Factor>`, and `NumFactorStatic.cpp` and `NumFactorDynamic.cpp` are 12 and
9 lines, explicit instantiations and nothing else.

### The mmd pair needs no design work

Three differences, all at the ends: the signature, `QuotientGraph qg` against
`QuotientGraphCompacted qg`, and the reporting tail, `gMmdCompactions = qg.numCompactions()` against
`*numBornCliqueMembers = qg.numBornCliqueMembers()`. Everything between, including the batch loop,
the refresh, and both the q2h and qxh walks, is BYTE-IDENTICAL. Template on the graph, instantiate
twice.

### The amd pair had two hunks; both are closed

Both were in detection, the stamp walk and the compare walk, two spans in the flat driver against
one in the compacted. Neither was algorithmic: the compacted run was contiguous with the pivot at
index 0 while flat's `[A, I]` put the entry to skip in the middle. Flipping the flat class's amd run
order closed both without moving a walk into a class or inventing an accessor.

### Why this is worth doing, and it is not tidiness

These pairs exist to price layouts against each other. A shared body makes STORAGE THE ONLY
DIFFERENCE MEASURED, where today it is storage plus whatever has drifted, and the compiler enforces
it instead of a periodic reading of the alignment table.

### Three costs to weigh first

The header-only requirement carries over: each driver compiles in its own translation unit with its
graph, so the template has to live in a header or the reason those timings are comparable goes away.

`-Wall` on a template reports only what is instantiated, so a mistake in a branch neither
instantiation reaches goes quiet. Minor with exactly two instantiations, both gated.

And the CHAINED driver does not follow. `QuotientGraphChained` carries runtime flags for list order
where `QuotientGraphCompacted` bakes them in, so a three-way template is a harder problem than
either two-way one. Do the pairs first and leave `MmdChained` on its own body.

### THE ORDER MATTERS, AND "ALIGN FULLY FIRST, TEMPLATE LATER" IS THE WRONG ONE

The natural reading of the above is to make the four drivers identical except for the graph type and
template afterwards. That step has almost no content, and where it has content it needs the
template's design decision anyway.

For the mmd pair the body is ALREADY identical; what is left is the signature and the tail. Aligning
those without the template means adding a counter to a class that has no use for it:

```
QuotientGraph            has numBornCliqueMembers, no numCompactions
QuotientGraphCompacted   has numCompactions, no numBornCliqueMembers
```

`numCompactions()` on the flat class is honest and free, a `constexpr 0`, an arena never compacting.
`numBornCliqueMembers()` on the compacted class is not: it does not track that, so it would mean
adding an accumulator to the elimination path for a probe that class does not want, or returning 0
and lying. `QUOTIENT_GRAPH_USAGE.md` item 5 already records this pair as layout-driven and left
alone. The template absorbs it in one place instead, three lines of `if constexpr` on a trait.

For the amd pair there is no cheap alignment step at all: its two hunks ARE the walk decision.

### THE FOURTH OPTION, FOUND AND DONE 2026-08-28

The two amd hunks are not a flat-against-compacted difference. They are ONE CLASS MISSING A KNOB THE
OTHER HAS.

`QuotientGraphCompacted` carries the run order PER BRANCH, two accessor pairs over one segment:
`adjacencyMmd`/`incidenceMmd` give `[A, I]` and `incidenceAmd`/`adjacencyAmd` give `[I, A]`.
`QuotientGraph` has one pair serving both branches, always `[A, I]`. The four cells:

```
                run       pivot in I    pivot in the RUN     detection
flat  + mmd    [A, I]     back          back                 (unused)
flat  + amd    [I, A]     front         FRONT                one loop      <- flipped
comp  + mmd    [A, I]     back          back                 (unused)
comp  + amd    [I, A]     front         FRONT                one loop
```

**AND `[I, A]` IS THE ACCEPTABLE CELL, NOT THE BEST ONE.** `[A, I]` WITH THE PIVOT APPENDED AT THE
BACK also gives a contiguous span, positions `0 .. runSize - 2`, so detection is one loop that stops
one short instead of one that starts at 1. The length tests are unaffected, both sides counting the
pivot in `incidenceSize` so it cancels. That cell has the single detection span AND NO ROTATION AT
ALL:

```
[A, I], pivot appended     1 store per member per elimination     one detection span
[I, A], pivot rotated in   3 stores                               one detection span
```

So the amd convention costs TWO EXTRA STORES per member of C[p] per elimination, bought with
nothing but the reference's list order; mmd pays none. Small and non-zero, and it is the first case
found where following the oracle costs something measurable rather than merely being arbitrary.

IT IS NOT AVAILABLE WHILE THE ORACLE STANDS, changing the raw order and so failing `make amdorder`
entry for entry. It belongs with the other convention questions, a candidate for a layer with no
oracle behind it, beside the detector cells and the four faces of the list-order convention.

Detection has to stamp `A[u] | ( I[u] - {p} )`, which is the whole run minus one element, and that
is CONTIGUOUS EXACTLY WHEN THE PIVOT SITS AT AN END OF THE RUN. Two of the four cells give one span
and one gives the pivot in the middle. `AmdFlat` is the only bad one, and it is bad because the flat
class serves amd's front-of-I convention on mmd's run order.

**AND THE RUN ORDER APPEARS NOT TO BE PERMUTATION-VISIBLE.** `AmdCompacted` walks `[I, A]` and
`AmdFlat` walks `[A, I]`, and they return the SAME permutation, verified on 38 cases and on 246
matrices with `nnz(L)` agreeing throughout. What is permutation-visible is the pivot's position
within `I[u]` and the direction of the walk, and those are identical in the two.

It held. Five methods changed, all amd-side: the two accessors, `beginEliminationAmd`,
`formReachableSetAmd`, `pruneAmd` and `absorbAggressively`. `massEliminate` needed nothing, its
`mAdjIncSrc[srcPtr] == pivot` test being true under both layouts. `make amdorder` stayed at 38 / 0,
so the flip is permutation-neutral in fact and not only by inference.

TWO THINGS WORTH KNOWING FOR THE NEXT PERSON. `pruneAmd` now builds the run with ONE cursor from the
run's start, incidence then adjacency, and inserts the pivot by THREE MOVES rather than an append;
the spare slot that needs is guaranteed because u is in C[pivot], so either A[u] drops the pivot or
I[u] drops an absorbed clique. And `absorbAggressively` compacts I RIGHTWARD, letting the dropped
entries fall off the front and advancing `srcPtr`, so A[u] never has to slide.

### THE GOAL IS ONE PLACE TO READ HASHING, and the walk move gets there on its own

The reason to do any of this is not tidiness: it is that reasoning about hash detection today means
reading both `AmdFlat` and `AmdCompacted` and holding the difference in your head.

Option 1 above delivers exactly that WITHOUT the template. Move the two walks into the classes and
the driver's detection loop becomes one identical text in both files, with the layout detail in
`stampGenerators` inside each storage class, where it belongs and where nobody thinking about the
algorithm has to look. That is a smaller change than templating and it is the one that pays.

So the sequence is:

```
1  MOVE THE AMD DETECTION WALKS INTO THE CLASSES     DONE DIFFERENTLY, by the run-order flip
2  TEMPLATE THE MMD PAIR                             the rehearsal, no unknowns
3  TEMPLATE THE AMD PAIR                             mechanical once 2 is done
```

Step 1 is closed and was not needed in the form written above: flipping the flat class's amd run
order gave the single detection loop without moving any walk into a class. Two steps remain, and
both are now pure bookkeeping, each pair differing only in the signature, the graph type and the
reporting tail.

### The gate

`make amdorder` and `make mmdorder` at 38 / 0, `make digest` unmoved over 365 digests, and for
steps 2 and 3 both `.cpp` files reduced to explicit instantiations. Nothing here changes an
ordering, so a failure is a mistake rather than a trade.

## OPEN, 2026-08-26: the merge detector is orthogonal to the degree rule, and the swaps are open

**THE TWO CHOICES ARE INDEPENDENT AND WE HAVE THEM WELDED TOGETHER.** A branch makes two decisions,
how it computes a degree and how it finds indistinguishable vertices, and today each branch takes
one pairing only:

```
                  degree rule        merge detector beyond mass elimination
mmd               exact              q2h, free and partial
amd               bound              hash, a hash per vertex and the general case
```

Nothing in either algorithm requires that pairing. Mass elimination is already shared, its test
reading a descriptor both prunes write, so the detector is the only part that differs, and it is a
separate axis from the degree.

**WHAT WE STAY ALIGNED TO IS A PHASE, NOT A CONSTRAINT.** `MmdCorrected` and `AMD_2` are how we
prove the implementations are right while we are still finding defects in them by the week. Once we
are comfortable, the oracles stop being the target and become the baseline a variation is measured
against. Nothing below is blocked; it is queued behind that confidence.

### Variations worth trying, roughly by increasing distance from what we have

```
1  mmd + hash          exact degree, general detection      the most merges any pair here can find
2  amd + q2h           bound, cheap partial detection       the cheapest pair
3  amd + both          bound, q2h first then hash on what it missed
4  mmd + neither       exact degree, mass elimination only  the floor, for attribution
```

(1) and (4) are the two that matter first, because together they say how much of mmd's cost and
quality is the exact degree and how much is the detector. Today those are confounded: every mmd
number we have is exact-degree-plus-q2h and every amd number is bound-plus-hash, so no measurement
we own separates them. (4) is nearly free to build, being a switch that skips a call.

(3) is the one with a real question inside it rather than a measurement: q2h's tag proves a subset
for free, so running it first should leave the hash less to do, but whether the hash then costs
meaningfully less is not obvious and neither is whether the bookkeeping of two detectors in one pass
costs more than it saves.

**AND (1) IS CHEAPER TO BUILD THAN THIS ENTRY IMPLIED, corrected 2026-08-28.** The standing reason
for not attempting mmd-with-hashing was that amd can afford a key only because it rides the bound's
walks, and that mmd would pay a pass of its own. That is wrong. `pruneMmd` already walks `A[u]` and
then `I[u]` per member of `C[p]`, so a key is one add per surviving entry on a walk that exists. The
fusion is in fact SIMPLER there: amd splits its key across the prune and the bound pass because
aggressive absorption sits between them and compacts `I[u]`, and mmd has no aggressive absorption at
all. What still stands in the way is the two link slots, `mPrev` and `mNext` not being free when the
overlay would want them, which means moving mmd's filing after detection and so moving the
permutation.

**AND (2) IS THE ONE NOT TO ASSUME.** q2h is described as free because it reuses a tag genmmd's
two-source walk has already assigned; that tag is a product of that walk's structure. Whether amd's
cliques-first update produces anything equivalent is unknown and is a question for an instrument,
not for reading. Every equivalence claim of this shape that we settled by reading this month was
settled wrong at least once.

### What the axis is really for

More merging means shorter reaches, cheaper degrees and more columns per pivot selection, which is
why both branches spend something on it. It also costs: a hash per vertex is added work, and the
ladder's own timing rule is that folding an array moves the clock while rescheduling work does not.
A detector is neither, so it has to earn its place on the 246 rather than by argument. The gate is
the same one every ladder rung has used: quality on the real set, and not slower.

## DONE 2026-08-24: the branch alignment starts, and the documentation was corrected first

**THE FRONT IS NOW mmd AGAINST amd**, which is a different axis from the flat-against-compacted work
`docs/QUOTIENT_GRAPH_USAGE.md` records. That one asks whether two STORES do the same thing the same
way; this one asks whether the two BRANCHES do, wherever nothing about minimum degree against
approximate minimum degree forces them apart. The test is the same and it is worth stating once: a
difference should be traceable to a vendored routine or to the algorithm, and where it is not, it is
an accident of how the two files were written.

**THREE ITEMS CLOSED, none of which moves a permutation.** `docs/DESIGN_DECISIONS.md` (2026-08-24)
has the account.
- **The minimum-degree seed.** `AmdFlat` and `AmdCompacted` computed it with a `min_element` pass of
  their own, which is neither `AMD_2`'s `mindeg = 0` nor the mmd branch's minimum-at-filing, and
  which answers for buckets that do not exist: the degree-zero and dense rows are numbered or set
  aside rather than filed, so one isolated vertex puts the seed at 0 with bucket 0 empty. Both amd
  drivers now take the minimum in the filing loop, over the vertices actually filed, which is
  genmmd's shape and the mmd drivers'.
- **`Buckets::refile` deleted.** No caller in the built tree since the 2026-08-21 retirement, and it
  was the only member of `Buckets` that named a `degrees` array, which is the amd branch's and which
  the mmd branch does not have. The class is now a list structure and nothing else.
- **A stale sentence in `Buckets`**, saying the heads are sized `n + 1` for slack, three lines above
  a constructor that sizes them `n`. Left behind by the 2026-08-23 change to file at the true
  degree.

**WHAT IS NEXT ON THIS AXIS, and the first item is a fork rather than a step.**

- **mmd HAS NO DENSE-ROW RULE.** `setAside` is called by `AmdFlat` and `AmdCompacted` and by nothing
  else. Item 3 below measures MMD at 70.7 ms against 0.83 on `GHS_indef/bloweybq` for want of it.
  Read that item's figures with its date: they are 2026-08-11, so its `AmdFlat` pair, 470 against
  1.5, describes a driver that got the rule on 2026-08-18 and no longer exists. The mmd pair is
  still current, no mmd driver having the rule.

  **THE MECHANISM IS ALREADY BRANCH-NEUTRAL AND THE SHARED CLASS NEEDS NOTHING**, checked
  2026-08-24. `setAside` does two things and each branch reads a different one: `mWeight[u] = 0`,
  which amd's `reachableSet` and prune skip on `nv > 0`, and `markGone(u)`, which mmd's
  `eliminatedMmd` skips on `mMark[u] == GONE`. Both halves are there because one class serves both,
  and each is inert on the branch that does not read it.

  **THREE THINGS TO SETTLE, none of them the mechanism.**

  - **The oracle, which is the fork.** `MmdCorrected` has no such rule, so the moment an mmd driver
    sets a row aside its permutation leaves the oracle and `make mmdorder` fails on every matrix
    with a hub. Either `MmdCorrected` gets the rule, a second deliberate edit to vendored code after
    the degree-scale repair, or the mmd branch stops being checkable by equality.
  - **The output path differs and would fail SILENTLY.** amd appends `denseRows` to `pivots` and
    calls `order`; mmd calls `orderAscending`, which reserves a supervariable's room with `pos +=
    mWeight[pivot]`, so a set-aside row arriving in the pivot list with weight 0 reserves nothing
    and drops out of the permutation. It has to be appended after the fact, as amd does.
  - **The threshold would be OURS on the mmd side.** `max(16, 10 * sqrt(n))` is `AMD_2`'s and we
    match it because it is the thing being matched. genmmd defines nothing, so the same constant on
    the mmd branch is a choice rather than a port, which is a different kind of claim to defend.
    Worth knowing what it is: a fraction that shrinks with n, 32 per cent of n at 1000 and 1 per
    cent at a million, so it fires well below a full row.

- **`unfile(u); restore(u);` WANTS TO BE ONE METHOD, and it is the cheapest item here.** Every mmd
  driver writes the pair adjacent and unconditional over each member of the new clique, which is
  genmmd's single store `bwd[rn] = 0` at `private/MmdVendored.cpp` line 89 spelled as two calls.
  Ours costs an extra load of `mPrev` and an extra branch per member, since each call guards its own
  entry condition: `unfile` reads the slot to decide the splice and returns early on `OUTMATCHED`,
  then `restore` reads it again to decide whether it says `OUTMATCHED`. genmmd asks once, knowing
  the state at the call site.

  **THE NAME AND THE BODY BOTH ALREADY EXISTED.** `BucketsChained` carried `evict(u)` defined as
  exactly those two lines, and went when that class turned out to be production's `Buckets`
  verbatim; the 2026-08-19 note below records it. Reviving it on `Buckets` decides both states from
  one load and reads as the single operation genmmd treats it as. One method, one call site per mmd
  driver, and no permutation can move.

  **NOT MEASURED, and the size of it should not be oversold**: a load and a branch in a loop over
  clique members, which is below what alpamayo can see on its own. The case for it is that it makes
  our code say what the vendored code says, and the guard in `unfile` stays either way, being there
  for the batch hazard rather than for this.

- **The two prepasses differ in three ways that are not cosmetic**, written up in
  `experiments/ordering/README.md` under "The two prepasses". That section concludes they are not
  the same idea, so it is a place to check the conclusion still holds rather than a place to align.
  One thread to pull: `Buckets::next` is called by the mmd drivers and by no amd driver, because the
  mmd prepass walks a whole bucket and the amd prepass walks nothing. Whether that is forced or is
  an artifact of the two prepasses coming from different codes has not been asked.
- **The rest of the driver-level inventory has not been taken.** What exists is a reading of
  `src/MmdFlat.cpp` against `src/AmdFlat.cpp`, 289 lines against 849, which produced the three
  closed items above and the fork. Nothing systematic, no table, and no equivalent of
  `QUOTIENT_GRAPH_USAGE.md` for this axis. **The one table that does exist** is
  `experiments/ordering/README.md`, "The interface each branch actually uses", which is every
  `Buckets` entry point against the two branches and is where both items above were found.

**THE LIST-ORDER CONVENTION IS ONE FACT WEARING FOUR FACES, and this is the item most likely to
cost someone a day.** Each of these looks like an independent branch difference and they are not:
they are the same tie-break convention, inherited from genmmd and `AMD_2` respectively, expressed
four times.

```
(a) where the new clique lands in I[u]     mmd appends at the BACK and leaves it
                                           amd appends then SWAPS it to the front
(b) which end of I[u] is walked            mmd BACKWARD, mReverseIncidence, amd forward
(c) which half of the run comes first      mmd adjacency then cliques, amd the reverse
(d) how the run is laid out, compacted     mmd A then I, amd I then A
```

**NONE OF THE FOUR IS FREE.** Each decides the order vertices enter `C[pivot]`, hence the order they
are filed into the degree buckets, hence which of several equal-degree vertices is picked. Each
branch reproduces its reference's permutation exactly or the comparison stops being one, so touching
any of these alone moves a permutation and `make mmdorder` or `make amdorder` fails.

**IT IS ALSO HALF A FLAG AND HALF A SUFFIX, which is the inconsistency.** (b) is a member,
`mReverseIncidence`; (c) is the `formReachableSetMmd` / `formReachableSetAmd` suffix. They are
always wanted together, and the compacted class's own comment says so, "two conventions under one
switch because they are one fact", while the flat class splits them across a flag and a name.
Neither arrangement is wrong; having both is.

**(d) IS WHY THE COMPACTED PRUNES ARE STRUCTURALLY DIFFERENT** rather than one loop moved, which the
flat pair is. Amd's layout lets ONE cursor sweep both halves, so its walk is one loop with an
`onCliques` flag; mmd's needs two cursors moving in opposite directions and so two loops. Only mmd's
order lets the walk set `adjacencySize = 0` mid-compaction, having finished the adjacency before it
starts the cliques, where amd is still holding an unread one.

**AND ONE THING HERE IS UNVERIFIED, flagged because reading has lost to instrumenting repeatedly in
this tree.** (a) and (b) appear to COMPOSE, so that both branches meet the new clique FIRST: mmd
puts it last and reads backward, amd puts it first and reads forward. If that is right, the two
conventions are one semantic rule reached by two routes, which changes what an alignment here would
even mean. It follows from reading `pruneMmd` and `pruneAmd` against `formReachableSetMmd` and
`formReachableSetAmd`, and nothing has instrumented it.

**THE MARK ARRAY, AND THE FOUR-LINK CHAIN THAT MAKES ITS WORST TEST LOOK ARBITRARY.** Traced
2026-08-24 by reading, not measured. The condition in question is
`formReachableSetMmd`'s adjacency test, `vWeight > 0 && !(mHasNumbered && mMark[v] == GONE)`, which
is three terms where the amd twin has one.

```
orderAscending reads mWeight[pivot] to reserve a supervariable's room
   -> a prepass-numbered vertex must keep weight 1 or it loses its slot
      -> the weight cannot record that it is dead
         -> mMark must carry GONE, and the walks must ask
            -> mHasNumbered exists to skip that load when nothing was numbered
```

Every link is reasonable and the top is nowhere near the bottom, which is why the condition reads as
an accident. **CUTTING THE TOP LINK DROPS THE BOTTOM THREE**, and the cut is small:
`orderAscending`'s first loop ALREADY walks the chain, `for (u = mSuperNext[pivot]; u != NIL; ...)`,
so it can count the members instead of reading the weight, and the count equals `mWeight[pivot]` by
construction. Then a numbered vertex can be weight zero, `vWeight > 0` catches it like any other
death, and `number()`, `mHasNumbered` and the three-term test all go. The gate is the digest.

**WHAT REMAINS AFTER THAT, and it is structural rather than contingent.** `mMark` and `mTag` stay,
doing DISTINCTNESS in `reachableSetWeight`: an exact degree counts each reached vertex once and u
reaches the same v through two cliques. That is what mmd IS, and amd escapes it by not
deduplicating at all, which is exactly what makes its number a bound rather than a degree. The
drivers' q2h stamps stay too, seven sites per mmd driver through `advanceTag`, `setMark` and `mark`.

**AND THE ELIMINATION WALKS NEED NO STAMP EITHER WAY**, on both branches, which is the part that
makes the above worth doing rather than merely tidy. `formReachableSet*` negates as it goes, so a
vertex already in the reach is negative and `vWeight > 0` rejects it. The negation is not the walk
being clever: the PRUNE requires it, reading `mWeight[v] <= 0` as "v is in C[pivot]" to get the
subtrahend of `A[u] = A[u] - C[p] - {p}`. The walk reads back a sign that is written anyway. So the
array is the price of the DEGREE REFRESH alone, the one call that forms a reach without eliminating
anything, and that call is also where the mmd branch's factor of 21 lives.

**A SMALLER ALTERNATIVE, independent of all of it.** `formReachableSetMmd`'s adjacency test can be
written `mMark[v] != GONE` today. It is equivalent THERE and nowhere else: mmd walks adjacency
FIRST, so nothing is negated yet, `vWeight > 0` is doing liveness alone, and every mmd death carries
GONE. It would not be equivalent in `formReachableSetAmd`, which walks cliques first, nor in either
clique loop, where the sign is doing distinctness. It costs one extra scattered load per adjacency
entry in a walk that runs once per elimination, and it buys one flag in place of three terms, drops
`mHasNumbered` to a single reader, and makes the test identical to `reachableSetWeight`'s adjacency
test. Not taken; the argument against is that it narrows defence in depth, the weight no longer
catching a death that stopped being marked.

**AND TWO COMMENTS ARE NOW STALE ON THIS SUBJECT.** Both class headers still call `mHasNumbered`
"load bearing and not merely an optimization", the short circuit that keeps an empty `mMark` safe in
shared bodies. That was true until 2026-08-24, when the two amd sites carrying the guard were
removed; every reader is now mmd-only and `mMark` is always allocated there, so the safety job is
finished and only the load-skipping remains.

**AND THE DOCUMENTATION WAS CORRECTED BEFORE ANY OF IT**, six files, all uncommitted with the three
above. Two were stale, four were wrong.

- `docs/QUOTIENT_GRAPH_USAGE.md` was entirely in the pre-2026-08-21 driver names, `Mmd3`, `Amd3`,
  `Mmd3C`, `Amd3B`, in all three tables and the ledger.
- `README.md` carried the enum count, the assertion counts, the vendored files under `src/` rather
  than `private/`, `MmdCorrected.cpp` nowhere, the alternative-store drivers described as free
  functions, and four ordering figures predating the 2026-08-23 re-run of `ORDERING.md`.
- **`CMakeLists.txt` DID NOT BUILD**, and this is the one that mattered. Its vendored glob named
  `AmdVendored.cpp` and `MmdVendored.cpp` and never gained `MmdCorrected.cpp`, while
  `OrderEngine.cpp` calls all three entry points under one macro, so every executable failed to link
  on any tree with `private/` present. Reproduced and then fixed and verified both ways. The
  Makefile's list was correct and its comment still said two files.
- **THE ASSERTION COUNTS WERE WRONG IN CLAUDE.md AND IN `docs/TESTING_SPECIFICATION.md`**, and they
  are now MEASURED: 265 with `private/` and 237 without, `test_order` at 73 or 45. `CLAUDE.md` said
  251 and 265, which read as the public build having more assertions than the private one. The two
  halves went stale two days apart rather than being swapped, and the impossible pairing went
  unremarked, which is the third instance of the count drift that file already documents twice.

**ONE GAP FOUND WHILE REWRITING THE TEST SPECIFICATION, not acted on.** `AmdCompacted` is the
default ordering and NO SUITE RUNS A WHOLE SOLVE THROUGH IT: `test_pipeline`'s sweep names five
enumerators by hand and `AmdFlat` is the amd one in it. The two return one permutation, so nothing
about the ordering is unchecked; what is unchecked is the default as a default. Adding it is one
list entry and one constant, and the comment there says an added enumerator is meant to make the
assertion fail until someone has ruled on it, which `MmdCorrected` did not trip because both are
hardcoded.

**AND A SECOND BATCH LANDED AFTER COMMIT `937fb30`, all of it uncommitted with this note.** Every
item below was gated on `make test` 265/265 and 237/237, `make digest` all identical over 365
against `937fb30`, and both alignment checks at 38 of 38; the ones that move a store also had ASan
and UBSan.

- **THE COMMENT TRIM.** 1307 comment lines removed across the three quotient graph headers, and NO
  CODE CHANGED, checked by comparing every non-comment line. What went is the four kinds
  `docs/WRITING_RULES.md` names: measurements, dated history, alternatives tried, and justification
  for the arrangement. What stayed is the invariants a reader can break. 224 of those lines were
  literal duplication between the amd and mmd twins.

```
                        lines            comment          >10 blocks   longest
QuotientGraph.h    1922 -> 1370   1239 ->  687  (50%)      38 -> 11     52 -> 27
QuotientGraphChained  980 -> 734   615 ->  369  (50%)      20 ->  5     48 -> 17
QuotientGraphCompacted 1308 -> 1280 533 ->  505  (39%)     13 -> 12     56 -> 41
```

- **THE SIGN RESTORE MOVED INTO THE PRUNE**, both branches, the prune being the last reader of the
  negated form. `restoreWeight` and `restorePivotWeight` are deleted from both classes with their
  four call sites, and `mLateMassElimination` loses its SECOND job, having decided both "do not
  merge here" and "do not restore here" while being named for the first alone. Measured at nothing,
  and in the unfavourable direction: the restore was a free rider on an existing walk and is now a
  pass of its own. `experiments/ordering/AMD3.md` iteration 28 has it.
- **FIVE UNCALLED MEMBERS DELETED.** `reachableSet` and `reachableSize` on the flat class, called
  from nowhere; and THREE ON THE CHAINED CLASS THAT WERE DECLARED AND NEVER DEFINED, which nothing
  could have caught, an uncalled declaration never having to link.
- **Mmd BEFORE Amd at every suffixed pair**, 8 pairs in the flat class and 10 in the compacted one,
  declarations, definitions and the calls inside the wrappers.
- **`absorb` is `absorbAggressively`**, since the only list ever passed to it is the driver's
  `deadCliques`. Elimination-time absorption owns no function: the kill is in `beginElimination*`
  and the removal rides the prune's compaction.
- **`nv` is `vWeight`**, 31 sites; `reachableWeight` is `reachableSetWeight` with a local
  `totalWeight`; the flat walks are `formReachableSetMmd` and `formReachableSetAmd` taking a
  `reachableSet` parameter. **The naming rule this settles**: a noun-named function IS its return
  value, so a query keeps the noun and an operation takes a verb. `form` rather than `compute`,
  matching `formStaticUpper` and the header's own "formed on demand and never stored".
- **THE LOOP COUNTER SCHEME**, 73 loops across seven files. The counter is the letter of what the
  body assigns, plus `k`: `vk` over an adjacency or a clique's members, `ck` over an incidence,
  `uk` over `C[pivot]`. It replaces `k`, `i`, `ii`, `ri`, `a`, `kk`, `uak` and `uik`. The rule is
  script-checkable and the only five loops that disagree are correct.
- **TWO DEAD `++mTag` BUMPS REMOVED** from the two `formReachableSet` walks, which read `mMark`
  only against the CONSTANT `GONE` and never write it. The tag is a consumable, `GONE` being
  `INT32_MAX`, so this is slightly more than tidiness.
- **AND THE TWO amd SITES THAT MENTIONED `mMark` ARE GONE**, in `formReachableSetAmd` and
  `pruneAmd`. Both were provably dead, `mHasNumbered` being written only by `number()` and called
  only by the mmd drivers. **The justification is that `QuotientGraphCompacted` ALREADY DID THIS**:
  its three amd functions mention neither array, so the flat class was the outlier and this is an
  alignment defect rather than a design question.

**THE EXPERIMENT README GAINED FOUR SECTIONS**, all reference material rather than narrative: the
`Buckets` interface against the two branches, the class surface of flat against compacted,
indistinguishability with the three mechanisms that find it, and what `mWeight` and `mMark` encode
between them as a state table per branch.

## DONE 2026-08-23: the mmd branch files at the true degree, and the oracle changed

**Our three mmd drivers no longer return genmmd's permutation, deliberately.** genmmd files a
vertex under its degree in `mmdint` and under its degree PLUS ONE in `mmdupd`, two scales live in
one bucket array, so a refreshed vertex is penalised by one against a vertex no pivot has reached
and the minimum selected is not always the minimum. `private/MmdCorrected.cpp` is genmmd with that
repaired and is the branch's oracle from here; `private/MmdVendored.cpp` is frozen beside it and
nothing compares against it. `docs/DESIGN_DECISIONS.md` (2026-08-23) has the account, the 4 by 4
evidence and the five changes in the corrected copy.

**MEASURED ON REAL STRUCTURE THE SAME DAY, AND IT IS A WASH.** `ORDERING.md` and `ACCURACY.md` have
both been re-run and rewritten. Against the frozen genmmd over the 246, the corrected rule fills
LESS on 92 matrices and MORE on 83, with a median ratio of exactly 1.0000; the total goes the other
way at +1.8 per cent and three PARSEC-class matrices are 98.9 per cent of that. `PERFORMANCE.md`'s
run says the fill difference costs no factorization time at all, a half per cent on the faster side
across all three traversals.

**WHAT THIS LEAVES FOR THE NEXT SESSION.**

- **`PERFORMANCE.md` has not been redone.** Its run is on disk. Its mmd pair is now two different
  ORDERINGS rather than one computed two ways, which `matrix_performance.cpp`'s header explains;
  the report's prose has not caught up. Florin also wants `MmdChained` and `MmdCompacted` rows
  there, which is a benchmark change and a longer run.
- Nothing else. `experiments/ordering/README.md`'s 246-matrix claim was checked and holds: `make
  mmdorder` makes the same per-matrix comparison against `MmdCorrected` and its 2026-08-23 pass
  flagged nothing on any of the three mmd drivers.

## OPEN, FOUND 2026-08-21: the amd default costs `Oberwolfach/LFAT5000`

**`make accuracy` under `AmdCompacted` is killed by the OOM killer on that matrix, where mmd
factors it.** It is the one open defect from this session and nothing else in the tree shows it.

```
                        predicted fill   actual        ratio     delayed     outcome
MmdCompacted, 08-21          67 463      15 678 721    232.4x    3 128 750   bwd 9.7e-20
MmdCompacted, 08-23          59 972      18 797 474    313.4x    6 250 000   bwd 7.6e-20
AmdCompacted, 08-21          54 978          -            -          -       KILLED, signal 9
```

**Amd predicts LESS fill on this matrix and still dies**, so the cause is the delayed-pivot cascade
rather than the ordering's fill. `LFAT5000` is positive definite, every pivot acceptable, and
dynamic LDL delayed millions of columns anyway under mmd, 3.1 in the 2026-08-21 run and 6.25 in the
2026-08-23 one; under amd it delays past the memory budget. `benchmarks/matrices/ACCURACY.md` has
the mmd account, in "The cost of accuracy".

**AND THE 2026-08-23 ROW IS THE SHARPEST EVIDENCE YET, FOUND BY ACCIDENT.** The mmd branch left
genmmd's degree scale that day for reasons that had nothing to do with this matrix. On `LFAT5000`
the new rule predicts ELEVEN PER CENT LESS fill and the factorization produces TWENTY PER CENT MORE,
with the delay count almost exactly doubling, 1.998x.

**So predicted fill and delay count have now moved in opposite directions three times on this one
matrix, and the third time nothing but the ordering changed WITHIN a branch.** The first two
readings could be argued to be a difference between amd and mmd. This one cannot: same code path,
same numeric phase, same input, a better ordering by the rule's own definition, and twice the
delays. Whatever drives the cascade here is not the quantity the analysis minimises.

**Worth one look before it is reasoned about further:** 6 250 000 is 312.5 delays per column over
n = 19994, and 1.998 times the previous count. Round numbers in a delay count may be a cap, a growth
policy or coincidence, and nobody has checked which.

**The accuracy pass is held at `MmdCompacted` meanwhile**, not at the tree's default, and
`benchmarks/matrices/matrix_accuracy.cpp` says why at its head. `MmdCompacted` returns `MmdFlat`'s
permutation exactly, so every published figure in that report stays comparable to the others.

**THAT PERMUTATION MOVED ON 2026-08-23 AND `ACCURACY.md` HAS BEEN RE-RUN AND REWRITTEN.** The hold
on `MmdCompacted` still makes sense and the three mmd drivers still agree with each other. The
classification held exactly, 106 solved and the same 30/18/56/2, with Cholesky agreeing with the
inertia on every row for the fourth consecutive sample; residuals and delay counts moved, `LFAT5000`
most of all. The amd row above is still from 2026-08-21 and has NOT been repeated since the branch
changed, so that one comparison is between two things of different ages.

**THE CONTROL WAS DONE ON 2026-08-21 AND IT ISOLATED THE ORDERING.** `make accuracy` under
`MmdCompacted`, same 119 files, reproduced the published run exactly: 106 solved, 12 singular, 1
killed, the same 30/18/56/2 classification, and the same eight residual rows to every digit. That
was against the old mmd ordering, and `LFAT5000` factored there at 232.4x with 3 128 750 delays. So
the two runs differ in the ordering and in nothing else, and the matrix was lost to amd rather than
to a changed input set. It also confirms in passing that
`MmdCompacted` and `MmdFlat` are interchangeable here, which is what let the accuracy pass move
onto the compacted store without re-measuring anything.

**WHAT WOULD SETTLE THE MECHANISM.** Whether the cascade is amd-specific or a threshold effect mmd
happens to sit under: run the amd pass with a raised memory ceiling and read the delay count
against mmd's, now 6.25 million. If amd merely delays somewhat more, this is a tuning question and
`docs/TODO.md`'s dynamic-pivoting entries are where it belongs. If it delays an order of magnitude
more, the interaction between an approximate-degree ordering and delayed pivots is a real finding
and nobody here has looked at it. `matrix_accuracy_cpp --max-fill=2e8` is the flag.

**NOT A REASON TO MOVE THE DEFAULT BACK**, on one matrix of 119 and with the mechanism unknown.

## Where the ordering subsystem stands

**CURRENT AS OF 2026-08-21, AND THE TREE IS CLEAN** apart from the accuracy work above.

**THE DRIVERS WERE RENAMED AND THE DEFAULT MOVED.** A driver is now its branch plus the clique
store it runs on: `MmdVendored` `MmdFlat` `MmdChained` `MmdCompacted` and `AmdVendored` `AmdFlat`
`AmdCompacted`, eight enumerators with `Natural`. The three alternative-store drivers joined the
enum, having been free functions. `AmdCompacted` is the default and `MmdCompacted` is the mmd
driver to use; the flat pair stay as the unbounded reference the bounded ones must equal, and
`MmdChained` stays as the measured alternative. See `docs/DESIGN_DECISIONS.md` (2026-08-21), twice.

**THE REASON FOR THE COMPACTED PAIR IS THE BOUND, NOT THE CLOCK.** A fresh run of all 246 put the
compacted store within 1.4 per cent of the arena on both branches, winning 188 of 235 on amd and a
coin flip on mmd. The store is bounded from the input and our arena is not. See
`benchmarks/matrices/ORDERING.md`, which is a snapshot of 2026-08-21.

**THE ORDERING SUBSYSTEM IS OTHERWISE SETTLED.** Eight enumerators, each of ours reproducing its
branch's reference exactly; two quotient graph classes whose driver call sequences are IDENTICAL on
the mmd side and differ in two calls on the amd side, `docs/QUOTIENT_GRAPH_USAGE.md` having the
ledger. What is open is below and, apart from `LFAT5000`, none of it is cleanup.

**THE FIGURES IN THAT SENTENCE WERE THREE AND NINE UNTIL 2026-08-21.** They came from an extraction
that read `qg.<method>` out of the drivers without stripping comments, so a method named in prose
counted as a call. All three mmd differences were comment text and one of the nine amd ones was.
What actually remains on the amd side is two `adjacencyAmd` calls in `AmdFlat`'s hash-detection
block, which is the run order and not an accident; `compactions` against `arenaEntries` is excluded
throughout, being layout. Corrected in `docs/QUOTIENT_GRAPH_USAGE.md` and in the 2026-08-21
`DESIGN_DECISIONS.md` entry.

**WHERE TO START: "The inlining bias" below.** It is short and it invalidates figures quoted
elsewhere in this file and in `benchmarks/matrices/ORDERING.md`. Read it with the caution the
section before it carries: the account of what translation units cost does not predict its own
results, and one movement it cannot explain is on the record.

**WHAT THE LAST TWO COMMITS DID.** `19efe5f`: the three quotient graphs became shared and
header-only, `src/QuotientGraph.cpp` deleted with its bodies folded into its header, so every
ordering driver now compiles its graph into its own translation unit; `MmdCompacted` and
`MmdChained` were both
brought onto shared classes; the real-matrix tables gained compaction counts and `MmdCompacted`;
`benchmarks/matrices/ORDERING.md` was regenerated from a fresh run of all 246. `97f4bc6`: the
vocabulary sweep, element to clique and garbage collection to compaction throughout our own code
and the twins, plus four new sections in `experiments/ordering/README.md` reading the 1996 AMD
paper against what we have built.

`MmdCompacted` was brought level with `AmdCompacted` on storage before any of that, in four steps:
the elbow room
off `nzaat`, `cliqueCountBalances`, the walk in positions and cursors, and the mid-walk collector
with its absorbed capture and contraction report.

What got done and is now behind us:

- `AmdFlat` and `AmdCompacted` return `AMD_2`'s PRE-POSTORDER PERMUTATION ON ALL 246 real matrices,
where
  before four kinds of difference showed on about sixty. One of the four fixes was a CORRECTNESS
  BUG in supervariable detection that had been costing fill silently.
- `MmdFlat` and `MmdChained` return genmmd's permutation on all 246, and needed no fixes at all.
- Peak live clique members is published by all five drivers and CHECKED across each branch's pair,
  in `test_order` and in both benchmark tables. It found two defects in `AmdCompacted`'s collector
that   nothing else in the tree could see.
- The layout question is measured: `AMD_2`'s compacted pool beats our arena by five per cent,
  genmmd's chaining loses to it by 56. **BOTH FIGURES WERE TAKEN IN THE MIXED BUILD ARRANGEMENT;
  see "The inlining bias" below before quoting either.**

## The point of the B and C layers, restated because it was being traded away

A B or C layer exists so that a differential against a vendored code is APPLES TO APPLES: identical
algorithm, identical encoding, identical storage layout, so that whatever still differs is either
the layout being priced or a real improvement to carry back into `AmdFlat` and `MmdFlat`. **If the
layout
is not identical the differential is not a differential**, and the file has no purpose. That is a
validity question, not a cost-benefit one, and it was briefly treated as the latter.

```
                 our arena      genmmd's layout    AMD_2's layout
mmd ladder       MmdFlat           MmdChained              MmdCompacted
amd ladder       AmdFlat           Amd3C   NOT BUILT  AmdCompacted
```

Production keeps a separate append-only clique arena with no constraint and is not required to
match anything. `AmdCompacted` and `MmdCompacted` must match `AMD_2` exactly. `MmdChained` must
match genmmd exactly.

## Where faithfulness actually stands, 2026-08-18

**`AmdCompacted` against `AMD_2`: the storage is closed, and so are both schedule items found by the
audit.** Done in this order, each verified alone, and no permutation moved at any point:

- CLOSED, in-place construction, `AMD_2`'s `if (elenme == 0)`.
- CLOSED, the reclaim, `AMD_2`'s `if (elenme != 0) pfree = p`.
- CLOSED, THE RUN ORDER: incidence first, adjacency behind it, which is `AMD_2`'s and the reverse
  of production's. This is what made the rest expressible, a consumed part now being a PREFIX. It
  also let the new clique go in by `AMD_2`'s three-move rotation instead of by holding a vertex
  back and swapping afterwards, and it retired two flags this layout cannot serve two values of,
  the list order and the reverse incidence walk.
- CLOSED, THE WALK IN POSITIONS off a base hoisted once. The pool never reallocates, so a
  collection invalidates a block's OFFSET and never the array's base; a position survives one and a
  pointer does not. The 277 ms this was feared to cost did not appear: it is the same instruction
  either way.
- CLOSED, WHEN THE COLLECTOR RUNS. Per entry, truncating both lists and carrying the half-built
  clique, so `beginElimination` no longer reserves for a worst-case reach of n.
- CLOSED, THE ELBOW ROOM: `nzaat + nzaat/5 + n` off the OFF-DIAGONAL count, which is `AMD_2`'s
  `iwlen` exactly. It had been computed from `nnz` with the diagonal and ran about 1.2n large.
  Compaction counts now read 1 against 1 from 8 to 800 a side, at the same headroom.
- CLOSED, THE CLIQUE TRIM AFTER DETECTION, `AMD_2`'s `Iw [p++] = i` in RESTORE DEGREE LISTS. The
  hash-absorbed had been staying in the clique for the rest of the run.
- CLOSED, DETECTION'S STAMP AND REJECT. One blind loop over the run from index 1, and a reject on
  the two stored lengths before the candidate's list is touched.

**Verified by varying the headroom, which is stronger than the digest alone.** At the shipped
reserve the collector never runs on grids, so the mid-walk path would go untested. Cut to `nzaat`
with no elbow room at all it fires on nearly every large pivot, 3 to 14 compactions from 3 to 200 a
side, and every permutation is byte identical, clean under ASan and UBSan. That also confirms
`AMD_2`'s own claim that it runs with no elbow room, only slowly.

**AND FOUR MORE THINGS CLOSED SINCE, none of them storage:** the empty-row prepass, the dense-row
rule, the rotation in `absorb`, and the stamp base. Those are what took the amd branch from about
sixty differing matrices to none, and the last was a correctness bug. They have their own section
above; this list is the storage half and stops here.

**`MmdChained` against genmmd: clique storage faithful, vertex list NOT.** The links, the chain
following, the terminator rule and the unchecked first loop all match `mmdelm` line for line. But
genmmd keeps ONE MIXED LIST per vertex, variables and elements together, told apart by the sign of
`fwd`, where we keep two sublists in a run with explicit lengths.

**One thing was closed there today: `MmdChained` no longer compacts C[pivot].** `mmdelm`
mass-eliminates
with `qsize [md] += qsize [rn] ; qsize [rn] = 0` and does NOT rewrite the list; the merged vertex
keeps its place and every later reader skips it on `qsize [nb] != 0`. We had been compacting, which
cost a pass genmmd never pays and saved the skipped entries it pays on every read. Neither is
visible in a permutation, so no check in this tree would have caught it; it was found by reading
`mmdelm` while answering a question about clique lifetimes. Measured after: `MmdChained / MMD`
unchanged within noise, so the two costs are about equal, and the file is now honest rather than
faster.

**AND THE MMD BRANCH IS NOW THE OPEN FRONT.** `make mmdorder` in benchmarks/matrices has never had
the marker the amd table has, because genmmd's own order is what `mmd_order` returns and there is
no postorder to hook around. So the mmd side has no equivalent of the check that found four amd
divergences today. Whether the same kinds of difference exist there is UNKNOWN rather than ruled
out, and giving that table an equality check is the obvious next move.

**`MmdCompacted` has NOT followed `AmdCompacted`.** It is mmd on `AMD_2`'s layout and it uses the
plain list
convention, which the flipped run cannot serve: with the incidence part at the front there is no
way to append the new clique behind it without a free slot or a shift, which is why `AMD_2` inserts
by rotation. That needs a decision before the C cell can be brought level, and until it is, the
`MmdCompacted` column is measuring the old layout.

## The one thing left in `AmdCompacted`, and it is a trade rather than a defect

**Our prune has `AMD_2`'s scan 1 fused into it, and that is what forces the separate `absorb`
visit.** `AMD_2` makes two passes over the new clique: scan 1 walks element lists only and computes
`W[e] = |Le \ Lme|`, then scan 2 walks each list once and does the compaction, the degree, the hash
key, the absorption and the mass-elimination test inline. Because absorption cannot be decided
until all of scan 1 is done, and ours is fused into the prune, we cannot absorb in the walk that
rewrites the list and have to revisit it afterwards.

So the fusion saves a pass and costs a pass, and which way it comes out is a measurement rather
than an argument. It is also the last thing standing between `AmdCompacted` and a pass-for-pass
copy, so
the `AMD / AmdCompacted` column is priced against a schedule that differs in exactly this one place.

Two smaller notes attached to it. The slide in `absorb` is already gone, the incidence part being
compacted rightward so the adjacency never moves, so what remains is the VISIT and not the copy.
And un-fusing would also let the bound live where `AMD_2` puts it, in `Degree[i]` rather than in
`w[u]`, which would retire the second stamp counter.

**Do not start this without deciding it is wanted.** It reverses a deliberate change of 2026-08-10
and it lands in the hottest loop in the file.

## Aligning the amd branch to `AMD_2` on real matrices, 2026-08-18

**Nothing on grids could see any of this.** The digest over 73 grids, both scaling ladders and the
38 acceptance cases in `experiments/ordering` were green before and after every one of the four
changes below. All four were found through `benchmarks/matrices`, `make amdorder`, and three of the
four needed a matrix Florin supplied by hand.

**THE MARKER IS WHAT MADE THEM VISIBLE.** `matrix_ordering.cpp` prints `raw order differs` on a row
where our permutation and the hooked pre-postorder `AMD_2` order disagree. It was added because the
table calls the fill column "both routines'" and nothing was checking that claim on this set.
Sixty-odd rows lit up immediately. They are now zero.

### 1. The empty-row prepass

`AMD_2` numbers every degree-zero vertex where it stands during initialization, ascending, in the
same pass that files everything else. We filed them at degree zero and popped them from the head,
so they came out LIFO: a pure diagonal gave `4 3 2 1 0` against its `0 1 2 3 4`. Identical fill,
reversed order, and twelve `m = 0` matrices in the set differed for this reason alone.

It rides in the filing loop and does NOT call `number()`. That function exists for a vertex
numbered while its neighbors still name it, the degree-ONE case; it sets `mHasNumbered` and puts a
test in every walk for the rest of the run. A degree-ZERO vertex is in nobody's adjacency, so
there is nothing to mark and not filing it is the whole of what has to happen.

### 2. The dense-row rule

`AMD_2` sets aside any row above `max (16, 10 * sqrt (n))`: `Nv [i] = 0`, kept out of every
reachable set, and appended to the output permutation at the end. We had none.

**It is a performance change as much as a fidelity one, and the size of it is new information.** A
hub of degree in the thousands that nobody sets aside sits in every reachable set it touches.
Measured on a 100-squared mesh with eight hubs of degree 3000, our ordering alone: 4.08 ms with the
rule, 80.28 ms without. On the real set it is what the worst rows were: `GHS_indef/bloweybq` 20.4 ms
to 0.43, `bloweybl` 41.0 to 1.24, `QY/case9` 12.7 to 1.24.

`setAside` on both quotient graph classes does it, and a ZERO WEIGHT is the whole mechanism, since
`reachableSet` takes a vertex on `nv > 0` and the prune keeps one on the same test. The threshold is
fixed at the vendored default rather than exposed: `AMD_2` reads `Control [AMD_DENSE]` and has a
control structure to read it from, this driver has none, and inventing one to hold a single constant
while that constant is the thing being matched is the wrong trade.

**The hook had to learn about them too.** `AMD_2` does not eliminate a dense row at all, so it never
reached the raw order and `amd_order_raw` came up short by `ndense` on any matrix with a hub.
`tools/hook_amd.py` now collects them at the branch that sets them aside and appends them after
`AMD_order` returns, which is where `AMD_2`'s own output assembly puts them, `Next [i] = nel++` over
i ascending.

### 3. The rotation in `absorb`

Both codes end I[u] as `[pivot][entries 1..k][entry 0]`. `AMD_2` rotates inside scan 2, where an
aggressively absorbed element has ALREADY been dropped, so entry 0 is the first SURVIVOR. Our prune
rotates first and absorbs afterwards, so when the entry parked at the back is one absorption then
removes, the two lists come out in different orders. That changed which vertex of a hash bucket
absorbed the others, which moved the permutation without moving the fill.

The correction is one read and a conditional `std::rotate` inside a pass `absorb` was already
making. On `GHS_indef/aug2d` it took 6 differing positions of 29008 to 0.

**This is the prune fusion showing its teeth**, and it is worth knowing that it can: the fusion of
scan 1 into the prune is what forces absorption to be a later pass, and the list order is a
consequence. See the section on it below.

### 4. THE STAMP BASE, and this one was a correctness bug

`w` holds two kinds of value: the scan's `w[c] = degree[c] + wflg - nvi`, which reaches
`wflg + lemax`, and detection's stamps. A stamp must be ABOVE every scan value of the same step.
`AMD_2` guarantees it with `wflg += lemax` BETWEEN scan 2 and detection. We raised the base at the
END of the step, so this step's stamps started at `wflg` while this step's scan values ran up to
`wflg + lemax`: the two ranges overlapped exactly, and a clique whose scan value landed on the
current stamp read as marked.

**So the exact test could return a false match and two vertices that are not duplicates could be
merged.** On `Grund/meg4`, n = 5860, vertices 5779 and 5780 were merged at pivot 5080 although their
lists differ in six of sixteen entries. That one false merge moved 109 positions and cost 297
entries of fill, 51809 against `AMD_2`'s 51512.

The fix is one statement moved: `stamp = std::max(stamp, wflg + lemax)` before the detection loop
rather than after it, in both drivers. It has presumably been wrong for as long as detection has
stamped into `w`, which is since the fifth fold, and NOTHING WE RUN COULD HAVE SEEN IT: it needs a
clique degree that lands on exactly the right value and it fired on one matrix in 246.

### How they were found, because the method is the transferable part

Reading the code failed twice on item 3 and once on item 4; each time the argument concluded the
two must already agree. What worked was intervention and logging:

- **Turn a feature off in both and recount.** Aggressive absorption off gave 0 differences on both
  `aug2d` and `meg4`, which said absorption was implicated without saying how.
- **Count the events.** Absorptions matched at 582 on `aug2d` and differed, 504 against 497, on
  `meg4`, which separated a different DECISION from a different EFFECT.
- **Log the event streams and diff them.** Pivots, merges and absorptions separately, since the
  interleaving differs harmlessly. That located `meg4` to one merge at one pivot.
- **Dump the actual state at that point.** Printing both vertices' lists from both codes showed
  them identical and NOT duplicates, which proved the test wrong rather than the lists divergent.

## Two bugs in `AmdCompacted`'s collector, found by a benchmark rather than by a test, 2026-08-19

**Both were introduced by the mid-walk collector**, step 3 of the storage alignment, and both were
invisible to everything the suite runs. `AmdCompacted` alone truncates a list while walking it;
`AmdFlat` has no collector and `MmdCompacted` still reserves before the walk, so neither could have
them.

**Neither moved a permutation.** `make digest` over 73 grids, `make amdorder` over 38 shapes,
`test_order`, the sanitizers and both scaling ladders were green throughout, before and after. They
were found by the `AmdCompacted pC differs` marker added to `benchmarks/matrices`'s amd table, which
compares `AmdCompacted`'s PEAK LIVE CLIQUE MEMBERS against `AmdFlat`'s. That figure is a property of
the
algorithm and not of the layout, so two drivers agreeing on the permutation can still be caught
doing different work, and they were.

**1. The absorbed cliques were read from a list the walk had destroyed.** `beginElimination`
re-read `I[pivot]` after the walk to kill them. `reachableSet` truncates that list as it consumes
it, so after a mid-walk collection it is short or empty and the cliques it named were never marked
dead: each kept a non-zero length, stayed a live block for `garbageCollect` to copy, and was never
subtracted from the live count. `eliminate` now captures the ids into `mAbsorbed` in the pass that
already walks the list to zero `w`, and `beginElimination` kills from that copy.

**And the placement of that kill is forced from both sides**, which is worth knowing because it
looks arbitrary. AFTER the walk, because the walk needs those cliques' member lists. BEFORE the new
clique is born, because the peak is a running maximum and an absorbed clique is never live at the
same instant as the one absorbing it; killing later inflates the peak on every step that absorbs
anything. My first attempt moved it after the birth and `test_order` caught it immediately, which
is the pC check earning itself a second time.

**2. The truncation itself was not counted.** Shortening the clique being consumed, `mRun[c]
.adjacencySize = ln`, is a contraction, and `killClique` afterwards then subtracted only what
remained. One line where the truncation happens.

**What is still missing, and it is the lesson.** `AmdCompacted` has no permanent
`cliqueCountBalances`, the debug recomputation the shared class has. Both bugs were diagnosed by
bolting one on temporarily; with it in place the second would have been a `test_order` failure
rather than a
benchmark marker on a matrix someone happened to send. Adding it is small and is the obvious next
step.

## The peak counter reaches all five drivers, and the layout answer, 2026-08-19

**`MmdChained` publishes `pC` now, and it cost an array.** The flat drivers read a clique's size
from
`mRun[c].adjacencySize`, a descriptor they keep anyway; the chained store has NO clique length at
all, ending a list at a terminator, so `MmdChained` carries `mCliqueLiveMembers`, one uint32 per
vertex.
Birth and contraction are free, the reach count and `merged.size()`; DEATH is the only event that
needs the array, a clique's size then being born minus contracted with nothing else carrying it.

**It tracks the NOTIONAL count, not this file's stored one**, which is the odd part and is
deliberate: `MmdFlat` drops the mass-eliminated from C[pivot] and `MmdChained` does not, `mmdelm`
leaving
them in place. Tracking what the flat drivers hold is what makes the comparison possible at all.
Read as a description of the chained store it would be wrong, and the declaration says so.

**Cost measured, not assumed:** O(1) per elimination and nothing in the member walks, so the ladder
shift was about two per cent with an inconsistent sign, inside a harness whose run-to-run spread on
an unchanged driver is up to seven. The real cost is footprint.

**Both benchmark tables now carry the branch's B sibling**, timed and with its order, fill and `pC`
checked against the driver rather than printed. All markers were clear over the 246 on both
branches.

### And the layout answer, which is the point of having the siblings

| | total | against its vendored | against our arena |
|---|---:|---:|---:|
| `AmdCompacted`, `AMD_2`'s compacted pool | 3.92 s | 1.11 | **0.95** |
| `MmdChained`, genmmd's chained segments | 82.5 s | 1.09 | **1.56** |

**The two say opposite things and both are consistent with the mechanism.** Compaction is a rare
sweep, about once per ordering at a fifth of elbow room, so it amortises to nearly nothing and the
pool BEATS our arena by five per cent, on 224 of 242 matrices. Chaining is a link test on every
read of every clique forever, and no headroom reduces it because chaining exists precisely to need
none, so it COSTS 56 per cent and wins on 60 of 237. In the band holding 96 per cent of the mmd
time, `MmdFlat` reads 0.69 against genmmd and `MmdChained` reads 1.09: a 31 per cent win turned into
a nine per cent loss.

**The grid ladders understate chaining badly**, 1.1 to 1.3 there against 1.56 on real matrices,
because grids have short cliques and few links to follow. benchmarks/matrices ORDERING.md carries
the full reading.

## What was borrowed into production, 2026-08-18

Two of the audit's findings were not `AmdCompacted` defects but shared ones, so they went into
`AmdFlat` and `QuotientGraph` as well. Both verified by digest, `make amdorder` and `make mmdorder`.

- **THE CLIQUE TRIM AFTER DETECTION.** `trimClique` on the shared class, called from `AmdFlat`'s
  restore pass, which now writes survivors back as it walks. One difference from `AmdCompacted`:
there is
  NO CURSOR TO PULL BACK, the arena's length being the vector's own and the clique not necessarily
  its last block, so the trimmed tail is left as a hole. The trim buys visits here, not space.
- **DETECTION'S BLIND STAMP AND LENGTH REJECT.** `sizeU` and `sizeV` are gone. Production gets TWO
  stamp loops where `AmdCompacted` gets one, and that is the layout rather than a choice: with
`A[u]` then   `I[u]` the entry to skip sits in the middle, where `AMD_2`'s order puts it at index 0.

**A third was considered and DECLINED: positions in the walks.** It works in `AmdCompacted` because
that
pool never reallocates, so the base is provably invariant and hoisting it is free. Production's
arena can reallocate, so the compiler cannot hoist `mCliqueArena.data()` across a `push_back` and
positions would add a load per member in the hottest loop. The reserve it would retire is not a
behavioural divergence there either: with no collector, reserving early changes capacity policy and
nothing observable.

## The peak clique counter, 2026-08-18

`numPeakCliqueMembers` and `numLiveCliqueMembers` on the shared class, so all six drivers have them.
Add on birth, subtract on death, take the maximum at BIRTH ALONE since nothing else raises the
total. MEMBERS, not entries: an entry is a slot the arena hands out and never takes back, a member
is a vertex in a live clique at this instant, and `arenaEntries` remains the cumulative figure.

**Why it exists.** It is what a CHUNKED clique store would need, and the flat classes can predict it
before any allocator work is done. On grids the ratio is already known, peak about 0.40x nnz(A)
against cumulative 0.94x in 2D and 1.45x in 3D; what it will say on real matrices is the open
question. It is the payload figure only, with no per-clique header or allocator rounding in it.

**Making it correct meant funnelling four sites, and the fourth was missed twice.** A clique's
length was being zeroed or shortened in `beginElimination`, in `absorb`, in `trimClique`, and in
`massEliminate` dropping the members it merged. Births and deaths now go through `killClique` and
`trimClique`. `merge` deliberately does NOT: `v` there is a live supervariable that never formed a
clique, so its length is `A[v]`'s.

**The check is `cliqueCountBalances`, and it lives in the class for a reason.** Two earlier versions
were wrong. Asserting the live count returns to zero is false: a clique dies when a member becomes a
pivot, and at the close of a run the last cliques have had every member mass eliminated instead, so
a handful of entries legitimately survive. Summing `cliqueSize` over a driver's pivot list is exact
for `AmdFlat` and wrong for `MmdFlat`, which pushes prepass vertices onto that list. Only the class
knows
which vertices ever formed a clique, so it keeps a debug-only owner list and does the recomputation
itself.

## The private types are named for what they are, 2026-08-18

```
AmdCompacted    QuotientGraphA  BucketsA  TaggedScanA                  ->  ...Compacted
MmdCompacted    QuotientGraphC  BucketsC  TaggedScanC                  ->  ...Compacted
MmdChained    QuotientGraphB  BucketsB  TaggedScanB  ApproximateScanB ->  ...Chained
```

All three are in anonymous namespaces, so two files declaring `QuotientGraphCompacted` is legal.
The terminology and the four-class scheme behind these names are in docs/DESIGN_DECISIONS.md
(2026-08-18): a SEGMENT is the logical unit, FLAT and CHUNKED name the allocation, and the eventual
set is `QuotientGraphFlat`, `QuotientGraphChunked`, `QuotientGraphCompacted` and
`QuotientGraphChained`.

**`tmp/` STILL GENERATES THE OLD NAMES.** `make_amd3b.py` and `make_mmd3c.py` rename `QuotientGraph`
to `QuotientGraphA` and `QuotientGraphC`, so regenerating either file undoes this. They also carry
sandbox paths from a previous session. Fix both before running them.

## The shared class is built and both drivers are on it, 2026-08-19: CLOSED

`QuotientGraphCompacted` lives in `include/oblio/QuotientGraphCompacted.h`, header-only, and
serves `AmdCompacted` and `MmdCompacted`. Eight suffixed pairs, `eliminate`
gone, `merge` shared through `markGone`, the mark array on demand. See
docs/DESIGN_DECISIONS.md (2026-08-19) for the inventory and the reasoning; the six-way split this
section used to describe was derived from copies that had drifted and is superseded.

**The bucket question that blocked `QuotientGraphChained` is answered and it needed no decision.**
`BucketsCompacted` and `TaggedScanCompacted` were production's types verbatim and are deleted.
`BucketsChained` is production's `Buckets` too: its `evict(u)` is `unfile(u); restore(u);`, which
`src/MmdFlat.cpp` already writes at lines 114 and 115. So promoting the chained graph needs no
second bucket class and no new method.

## The translation-unit alignment is done, 2026-08-19: CLOSED

All five ordering drivers now compile their quotient graph into their own translation unit, the
three graph classes being header-only. `BucketsChained` went the way of the other two copies, and
two dead `eliminate` declarations went with the chained class: overloads taking the amd branch's
scans, never defined and never called.

**WHAT IT BOUGHT AND WHAT IT DID NOT.** Our own drivers can now be compared with each other, which
they could not before: `MmdCompacted / MmdFlat` is 0.92 to 0.98 and `AmdCompacted / AmdFlat` is 0.82
to 0.84, both
sides built alike. It did NOT make anything comparable to `MMD` or `AMD`, and the three-point gain
predicted there never appeared.

**AND ONE UNEXPLAINED RESULT, worth knowing before trusting any of this.** `MmdChained` did not
change translation unit at all, only where its class's source sits, and it slowed by 4 to 7 per
cent. See docs/DESIGN_DECISIONS.md (2026-08-19).

## The inlining bias, and the numbers that need re-measuring, 2026-08-19

A class in an anonymous namespace inside its driver's file is inlined into the pivot loop; a class
in its own `.cpp` is not, and on alpamayo that is worth about 5 per cent. Every ratio published
before today between a private layer and a production one was taken in that mixed arrangement.
Full account and the three-row table in docs/DESIGN_DECISIONS.md (2026-08-19).

**WHAT THIS MEANS FOR THE NUMBERS ON RECORD:**

- `MmdCompacted / MmdFlat` is 0.92 to 0.95 on large grids, not 0.90. Corrected, both fair
arrangements agree.
- **`AmdCompacted / AmdFlat` at 0.95 on the 246 predates the alignment and has NOT been
re-measured.**
  It is now a both-split build, so a plain re-run of `benchmarks/matrices` gives the fair figure.
  Until then it should not be quoted. `benchmarks/matrices/ORDERING.md` and the header note in
  `src/MmdCompacted.cpp` both still carry mixed-build figures.

  **THAT BULLET WAS OVERTAKEN THE NEXT DAY AND IS WRONG AS IT STANDS, marked 2026-08-21.** It was
  written on 2026-08-19; `19efe5f` landed on 2026-08-20 and regenerated the report from a fresh
  run of all 246 with every driver already split. `ORDERING.md` therefore carries FAIR amd figures,
  and its own header says so. The `AmdCompacted / AmdFlat` row did not exist before that commit:
`af829bb`
  had `AmdCompacted / AMD` only, aggregate 0.95 and median 0.90 over 242 matrices, and the re-run
gives
  aggregate 0.94 and median 0.950 over 235. The two 0.95s are different statistics, which is how
  the figure looked unchanged and so looked un-rerun.

  **WHAT IS ACTUALLY STALE ABOUT IT IS ONE DAY LATER.** `ORDERING.md` was measured at `19efe5f`
  and so predates `97f4bc6` and `1da85c5`, which means it does not carry the `AmdFlat` speedup from
  alignment item 4. A re-run of `make amdorder` is still wanted; the reason is the alignment
  commits, not the build arrangement. And the header note in `src/MmdCompacted.cpp` is not mixed-
  build either: it says outright that its figure was 0.90 until 2026-08-19, that this was not a fair
  comparison, and that building both sides alike gives 0.92 to 0.95. So the bullet was wrong about
  both files it named.
- `MmdChained / MmdFlat` at 1.14 to 1.25 is UNDERSTATED, `MmdChained` still being private. Chaining
loses while   holding the advantage, so that conclusion is safe.
- `AmdFlat / AMD` and `MmdFlat / MMD` overstate our gap by roughly three points, permanently: both
  vendored files are single translation units with everything `static` inside. That is a real cost
  of our arrangement, not an error to correct away, and belongs beside the number rather than
  folded into it.

**INSTRUMENTS:** `tmp/unity_mmd3c.cpp` and `tmp/unity_mmd3.cpp`, two include lines each, built by
excluding the sources they subsume. `tmp/` is gitignored. The compile line is in the session log
and should move into `benchmarks/ordering/README.md`.

**AND A STANDING KNOB:** `-DOBLIO_PAD_ORDERING` staggers the padding of every size-n allocation in
`QuotientGraphCompacted` by whole 16 KiB pages, capacity only, so nothing computed changes. It is
the data-placement intervention from 2026-08-17 made repeatable. Run against the merge, it was a
null, which is how placement was ruled out.

## The ladders moved off the power-of-two sides, 2026-08-18

```
square   33 51 65 101 129 201 257 401 513 801 1025 1601
cubic    9 11 17 21 33 41 65 81
```

199 and 200 are gone; 201 stays as the rung 200 would have become. The 6 cube stays at 6, being the
largest size where `AMD` and `AMDraw` disagree on fill and therefore the only row that observes the
postorder's effect.

**WHAT THIS GIVES UP, and it was deliberate.** The interleaving was originally a power-of-two series
against a non-power-of-two one precisely so that an addressing artifact could be told apart from a
real trend. That diagnostic is gone. What is gained is that no row's denominator carries the
artifact, the vendored routine suffering it worst. If the question comes back, put the old sides
back for a run rather than reading these rows harder.

**AND THE PHASE SPLIT CHANGED WHAT THE AMD LADDER MEANS.** `make phases2d`, from the gitignored
`amd_timed.cpp` that `tools/hook_amd.py` regenerates, shows `AMD_aat` and `AMD_postorder` plus
`valid` at about 13 per cent of the vendored column, none of which we do. Against the comparable
region alone, `AmdFlat` reads 1.63x at 401 a side where the raw column says 1.43x. So the amd ladder
does NOT say containers are free; the denominator was inflated. That removes the apparent
contradiction with the mmd ladder's 15 per cent and takes the pressure off `MmdChained` being
misaligned as the explanation.

## What was done since commit 5ea68cc

**All of this is committed. The bullets below are `c21f447`; everything in the 2026-08-18 sections
above, plus the sections that follow this line, are `db377c9`.**

- **`src/MmdCompacted.cpp` and `include/oblio/MmdCompacted.h` REBUILT.** The old file was mmd on the
production
  arena, a transitional vehicle that had served its purpose; this one is mmd on `AMD_2`'s pool,
  the real matrix cell. Regenerated from current sources by `tmp/make_mmd3c.py`, so it inherits the
  shared class's encoding, then the storage hand-ported from `AmdCompacted`.
- **`src/AmdCompacted.cpp`**: in-place construction, the reclaim, and a `-Wunused-variable` under
`NDEBUG`.
- **`src/MmdCompacted.cpp`**: the same two, plus a real bug found by ASan on a 3 by 3 grid.
`AmdCompacted` negates
  the pivot in `beginElimination` and `MmdCompacted` inside `reachableSet`; the new in-place walk
was a   third path and negated nothing, so the restore flipped a positive weight negative and
  `orderAscending` wrote four billion entries past the permutation. Recorded at the site.
- **`experiments/ordering/README.md`**: the three-layouts section on Florin's one-axis
  classification, no constraint and no bookkeeping, a loose constraint and occasional bookkeeping,
  a tight constraint and constant bookkeeping; the headroom cliff; and the divergence account.

**Verification at every step**: 584 digests over the eight committed drivers identical, `test_order`
87/87, 279 across the suite, both build modes warning-clean, and the digest re-run under ASan and
UBSan with assertions live.

**Figures now known to be stale**, marked in place rather than deleted: the pool timings and the
headroom cliff were measured before the two fixes, so they price a variant and the cliff sits
further left than it should. `AmdCompacted / AmdFlat` is 0.885 as of the last run; the mmd figure
has not been re-measured.

### Other things open

**`AMD_2`'s aliasing is in the cubic ladder too**, and the earlier note that it would be rarer there
is WRONG: `m^3 = 2^3k` when `m = 2^k`, so n is a power of two at every starred cubic size as well.
At 64 cubed the vendored routine costs 536 ns per vertex against 314 at 80 cubed, which has half
again as many vertices. That is the largest instance of the artifact the benchmark has produced and
it has not been confirmed with the padded copy, which is the cheap check.

**The 80-cubed amd row is worth a re-run** before it is read as a trend: `AMD` barely slows from 64
to 80 despite twice the vertices and nearly three times the fill, which reverses a pattern that held
over five sizes.

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
with the MmdFlat storage investigation, which is uncommitted work and is item 0 below. **It is meant to be
deleted** once the items below are done or
abandoned.
Everything in it that outlives the task is already somewhere durable and this file only points at
those places:

- `docs/DESIGN_DECISIONS.md`, three 2026-08-10 entries: "a null result measures an implementation",
  "the algorithm was the smaller half", and the earlier "what the vendored AMD's speed is made of".
- `experiments/ordering/AmdFlat.md`, iterations 25 and 26, the two fusions and how they were found.
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
`MmdFlat` runs at 1.05 to 1.17x genmmd on square grids and 0.77 to 1.12x on cubes, against 1.35 to
1.48x that morning, and `MMD2` reads the same to within the noise. Nothing about what is computed
changed: every permutation and every nnz(L) is identical, on all 137 acceptance cases across nine
orderings.

**That MMD2 and MmdFlat coincide is the useful shape of the result.** They are different orderings
with different mechanisms; what they share is the quotient graph and its clique arena, so both are
now bounded by it rather than by anything of their own. And the trade the library makes is worth
stating in one line: our `std::vector` layer costs a few percent that genmmd does not pay, having
raw arrays in registers for a whole run, and the SECOND ARENA more than covers it. `MmdChained`,
which
is `MmdFlat` on genmmd's single nnz(A) storage with every encoding fold present in both, reads 1.12
to 1.23x where `MmdFlat` reads 1.05 to 1.17x.

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
`MmdChained` carries `MmdFlat`'s algorithm on genmmd's storage, cliques in the dead segment of their
own pivot with no second arena. Every encoding fold is now in both files, so storage is the only
difference left, and ours wins on every axis: 1.02 to 1.19x genmmd in 2D against `MmdChained`'s 1.15
to
1.38x, 14.22M instructions against 16.61M, 119331 D1 read misses against 123510. **Spending nnz(L)
on a second arena buys speed.** `MmdChained` therefore STAYS, as the standing equal-encoding
comparison
against the vendored storage scheme, which is a change to the stop condition its own header states.

**What is open on the mmd branch, and it is one item.** The TAG SCHEME. genmmd's refresh puts the
element tag ABOVE the per-vertex tags, `mt = tag + md0`, so an element member fails
`marker[nd] < *tag` automatically and needs no second comparison; ours draws both from one counter,
so `elementTag` is below `vertexTag` and every entry in the q2h and qxh paths pays an explicit
`m == elementTag` test. It is a change to what is computed rather than to how it is stored, which
is the category that actually paid, and the refresh is about a third of the run. Untried.

**Everything else on this branch was tried and did not pay.** Four attempts at the container layer,
each failing for its own reason, all in the "2026-08-15, later" design entry: the stamping fold
ported to `MmdFlat`, an arena cursor in place of `push_back`, raw bases in place of the accessors,
and q2h indexed rather than looped. Read that entry before reaching for any of them.

**And `MmdChained` has finished answering its question.** It is now genmmd's data structure
essentially exactly: one array of nnz(A), cliques in their pivot's dead segment, negative links, a
value terminator, no clique length array, no liveness array, the degree list in one link array. Its
obligation from here is to stay ENCODING-IDENTICAL to `MmdFlat`, so that the only difference between
them is storage; a fold that lands in `QuotientGraph` lands there too, or the comparison quietly
stops being about storage.

**And the amd branch has had NONE of this.** Its five drivers still carry `degrees`, `outside`,
`cliqueDegree`, `explicitPart`, `hashHead` and `hashNext`, plus a 2n `mark`, and `AMD_2` allocates
none of them: it overlays its hash buckets on the degree heads and keeps the running key in a link.
`AmdFlat` gained from the shared-class folds alone, reading 1.28 to 1.87x in 2D and 0.97 to 1.26x on
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
`experiments/ordering/README.md`, `experiments/ordering/REPORT.md`, `AmdFlat.md`, `MmdFlat.md` and
`docs/TODO.md` all carry ratios measured before it. Fill figures are unaffected, nothing having
moved. Each file carries a dated superseding note rather than being rewritten, since a dated
measurement is a record of a run.

**Where the amd branch stands.** `AmdFlat` returns `AMD_2`'s raw elimination order exactly on all 38
acceptance cases and now runs close to it: 0.83 to 0.89 ms at 16 cubed over eight runs where the
vendored routine reads 0.74 to 0.86, so they overlap, rising to about 1.2x at 32 a side. In 2D it
is 1.33x at 32 a side and about 1.95x at 400. Three days earlier it was 3.0x on cubes. `MmdFlat` is
0.87 to 1.05x genmmd on cubes and 1.26 to 1.55x in 2D, with no trend in n on either.

**THE HASH KEY IS NOW `AMD_2`'S EXACTLY, 2026-08-16.** All four divergences listed anywhere in
this file are closed in `Amd2`, `Amd2B` and `AmdFlat`: no `+ 1` on a term, the pivot's own clique
excluded, modulus `n` rather than `n + 1`, and a `uint32` accumulator that WRAPS at 2^32 the way
Amd.cpp's `UInt hval` does, one reduction at the end rather than one per term. It changed no
permutation, checked over 730 across ten drivers, and it corrected the grouping: our `pair` and
`stamp` counts now equal `AMD_2`'s digit for digit, where before we were UNDER-grouping on cubes.
Worth about 4 percent. See `docs/DESIGN_DECISIONS.md` (2026-08-16).

**What was done, 2026-08-08 to 08-10**, three commits: the hash key defect (ledger entry 8, worth a
factor of two to three on cubes), the key folded into the bound pass, and the first scan folded into
the prune. **The walk axis is finished**: `AmdFlat` walks `I[u]` twice per pivot and `A[u]` once,
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
2b. **THE PER-ARRAY-AT-SETUP CONTAINER EXPERIMENT NOW HAS EVIDENCE AGAINST IT, 2026-08-16.** It
   asked whether one allocation carved into arrays would beat separate vectors, on the model of
   `AMD_2`'s `S`. That is exactly the shape that makes the vendored routine alias against itself at
   power-of-two grid sides: six arrays at offsets that are exact multiples of n, landing in the same
   cache sets whenever n is a power of two, costing it 56 percent of its D1 read misses at 512
   squared. Grid benchmarks are where sizes cooperate to produce that, and grids are most of what
   this tree measures. If the experiment is ever run, it needs padding between the arrays from the
   start, and a ladder that can detect the artifact if the padding is wrong.

2c. **AND THE EARLIER FORM OF THE SAME QUESTION IS ANSWERED, 2026-08-16, not in the form it
   was posed.** It asked whether one allocation carved into arrays would beat separate vectors,
   `AMD_2` carving Pe, Nv, Head, Elen, Degree, W and Iw out of one `S`. That was the leading
   hypothesis for the 2D growth for one session and it was WRONG: `MmdFlat` has ten separate
   allocations against genmmd's five and is flat, so a doubling of allocations does not by itself
   produce growth. What did produce it was two random probes per clique visit into arrays indexed
   by a dead pivot's id, and folding the arrays away removed them. The container question is
   therefore still open as a question and is no longer a suspect.

3. **A dense-row threshold in `QuotientGraph`.** NEW on 2026-08-11 and the largest ordering-time
   item the real matrices found. A single vertex adjacent to everything makes minimum degree
   quadratic; the vendored AMD sets such rows aside before ordering, above `max(16, 10*sqrt(n))`
   entries, and places them last, and its own source says the cost of not doing so is O(n^2). Oblio
   has no such rule anywhere. On `GHS_indef/bloweybq`, one column of degree 10000 among 9992 of
   degree 5, removing it by hand takes MMD from 70.7 ms to 0.83 and AmdFlat from 470 to 1.5, while
the    vendored AMD does not move. **Fill is unaffected**, so this is time only. It belongs in the
   shared quotient graph so all six drivers gain it at once.
4. **The descriptor struct**, `docs/TODO.md` question 3, if that profile says stalls. Item 2d.
5. **Narrow the one-dimensional sizes. THE ORDERING IS COMPLETE, 2026-08-11**; the symbolic and
   numeric phases remain and are the larger half. Item 2e.
6. **Why `AmdFlat` and the vendored `AMD` disagree on fill on real matrices.** NEW on 2026-08-11,
and
   **the instrument for it now exists on the other branch, 2026-08-15**. `make mmdmatrices` runs
   the mmd alignment over `data/*/*.mtx` and reported 243 matched, 0 differed, 3 skipped on its
   first run, so `MmdFlat` reproduces genmmd on real structure and not merely on grids. An
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
   grid that layer is 2.78M of `MmdChained`'s 15.89M instructions, and roughly half of it, the
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
     that shows it: with no elimination work at all, an `MmdFlat` ordering is roughly a third
     `QuotientGraph` CONSTRUCTION and a sixth `orderAscending`, and construction allocates and
     initializes about ten size-n arrays where genmmd allocates five plus its 1-based copies. That
     is the array-count finding of the same day, moved from the loops into the constructor, and it
     is invisible on a grid because real work amortizes it.

   **Why it matters and where.** `benchmarks/matrices`, `make mmdorder`, has the evidence: the five
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
- `AmdFlat / AMD` per row is a poor measurement: the vendored column moves 16 percent between runs
  where ours moves 7. Quote absolute times with the vendored range beside them.
- A benchmark column reached as a free function is timed differently from one reached through the
  enum, by up to 2.4 percent. Two columns compared must go down the same path.

---

## The MmdFlat storage investigation, 2026-08-15: CLOSED

Everything this section held is superseded by `docs/DESIGN_DECISIONS.md` (2026-08-15), which
carries the result, the differential that found it, and the five hypotheses that failed on the way.
Two facts from it are worth having here because they are what a later reader would otherwise
re-derive:

- **The five dead hypotheses stay dead**: construction cost, the liveness array in the clique walks,
  the four-pass refresh preamble, index widths, and clique placement. `MmdChained` implements
genmmd's   placement in full and the time did not move.
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
nine random patterns at n = 2000. All match on alpamayo, so `AmdFlat` reproduces `AMD_2`'s raw
elimination order, member order included, on something far wider than the shape it was built
against.

Getting there found three things, none of them visible before:

- a defect in production `AmdFlat`, ledger entry 7, the stored clique degree not rewritten after
mass   elimination trimmed the clique;
- a use-after-free in the shared `QuotientGraph` that every ordering had, benign until an allocator
  recycled the block, which is why it surfaced as two machines disagreeing about integer code;
- two harness faults that looked like divergences, a dense threshold turned off by undefined
  behavior and a 3D grid builder emitting unsorted columns.

`experiments/ordering/AmdFlat.md`, iterations 18 to 20, is the narrative.

**And the ordering benchmark measures cubic grids.** `make run3d` and `make scale3d`, beside
`run2d` and `scale2d`, both axes now named in every target. It also carries an `AMDraw` column, the
vendored AMD's raw elimination order through the same hook the acceptance test uses, so `AmdFlat`
has
something to sit against that agrees by construction rather than nearly. What the first run found
is in `benchmarks/ordering/README.md` under "Cubic grids, 2026-08-09", and the short version is
that three standing claims were square-grid artifacts: our tie-break beating AMD's, MMD being the
ordering to beat, and the LIFO question having an amd-side answer.

**And the mmd branch has an acceptance test at last.** `make mmdorder` compares production `MmdFlat`
against genmmd's elimination order on the same four shapes, 38 cases, all matching. Until then the
mmd alignment rested on a scratch probe from 2026-08-07 that died with its session and on the
benchmark's fill column, which `MmdFlat.md` iteration 6 shows is not sufficient. It needs no hook,
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
inventory" and "What the inventory was worth", in `AmdFlat.md` iteration 25, and in
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
  being removed, when `AmdCompacted` held a verbatim copy of `AmdFlat` and the two benchmark columns
were
  the same code timed twice. A single figure under about 4 percent is not a result; what rescues a
  small effect is consistency of sign across sizes. **Worth arranging again**: whenever the next
  vehicle exists, run the benchmark once with it still a verbatim copy before putting anything in
  it, which costs one column and gives every other column an error bar under identical
  conditions.

**AND THE WALKS ARE DONE, 2026-08-10, later the same day.** The first scan is folded into the
prune, so `I[u]` is walked twice per pivot and `A[u]` once, which is `AMD_2`'s count exactly.
Over eight runs `AmdFlat` reads 0.83 to 0.89 ms at 16 cubed where `AMD` reads 0.74 to 0.86, so the
two
overlap there, and it is about 1.2x at 32 a side; 2D improved 0 to 8 percent as well. `AmdFlat.md`
iteration 26 and `docs/DESIGN_DECISIONS.md` (2026-08-10, "the algorithm was the smaller half")
carry it, with the corrections below.

**What is open, in the order we would take it.**

- **THE GROWTH TERM, which is now the whole question**, and item 2d below has the pass that would
  attack it. Both families are a low constant plus
  something that scales: cubes 0.97x at 16 a side rising about 25 percent to 1.2x at 32, 2D 1.33x
  at 32 a side rising about 45 percent to 1.9x at 400 with a knee past 280. The constant is nearly
  spent, and NOTHING measured so far touches the growth. The knee says memory. Two profiles of
  `amd3` compared AGAINST EACH OTHER rather than against the vendored routine, 140 and 400 a side
  in 2D, or 16 and 32 cubed, would say whether it is instructions or stalls. `MmdFlat` over the same
  quotient graph shows no growth on either family, which is the control that makes this an
  amd-branch property.
- **QUOTE ABSOLUTE TIMES, NOT RATIOS PER ROW.** The vendored column is the noisiest in the table:
  at 16 cubed it spans 16 percent across runs where `AmdFlat` spans 7, so `AmdFlat / AMD` mostly
measures
  `AMD`. And a free-function column is timed by `orderTimeFn`, a standing method by `orderTime`
  which also builds a Permutation, a difference of up to 2.4 percent, so those two are not
  comparable. Both cost a wrong claim on 2026-08-10.
- **Price a change by its STREAMS as well as its visits.** The pass inventory counts visits and is
  silent about size-n arrays, and on 2026-08-10 the arrays were the larger term in 2D: the same
  fold measured 12 percent slower there with two extra vectors and 0 to 8 percent faster without
  them. There is no instrument for this and one would be cheap.
- **`Amd2` and `Amd2B`, 2026-08-10: not by the cheap route, and PARKED rather than closed.**
  Checked on a scratch copy: fusing there as `AmdFlat` does moves the permutation on all ten grids
  tried. They refile inside their single-pass bound, so that loop's direction is already a
  tie-break input and cannot also serve the hash chain, which wants the opposite direction under
  head insertion. Tail insertion would preserve the order and is untried; it needs a `hashTail`
  array of size n, which is the footprint that made the first version of this fusion measure
  nothing, so expect nothing. Not worth a measurement before the profile below. `Amd1` and `Amd1B`
  have no hash and so no key. Only `AmdFlat` fuses free, because entry 4 moved its refile below the
  hash.
- **The stamp and the mass-elimination sweep are DEMOTED, not queued.** Both were next on the
  reasoning that produced the `degme` deletion, and that reasoning now has a counterexample.
- **`tmp/Amd3I.cpp` is the pre-fusion driver.** The pass inventory reproduces the OLD `AmdFlat`
until
  it is re-derived, which is fine for the before-and-after already recorded and wrong for anything
  further.

**And the 2D scaling divergence, which is larger than any of the above and is not made of passes.**
`AmdFlat` runs at 1.54x the vendored routine at 64 a side and 2.07x at 400, growing monotonically,
where `MmdFlat` over the same quotient graph shows no trend at all across those seven sizes. A
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
the adjacency alone. Two lines in nine files. On alpamayo `AmdFlat` on cubes goes from about 3.0x
the
vendored routine to 1.44x and `AMD2` from 2.85x to 1.40x, with both controls unmoved. It is ledger
entry 8, `AmdFlat.md` iterations 21 to 24 are the narrative, and `docs/DESIGN_DECISIONS.md`
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
`AmdFlat` goes 2.3x to 3.0x, so the amd branch's whole degradation on cubic grids is in the extras,
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

**And the comment at `src/AmdFlat.cpp` justifying this says the opposite**: "its scan 2 accumulates
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
0.33 and 0.48 recorded in the `AmdFlat.cpp` comment. That number is the one that moved by a factor
of three last time, and it is the only output any of these four can change.

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
enumerators, `Natural, MMD, MMD1, MMD2, MmdFlat, AMD, AMD1, AMD2, AmdFlat`, and the two B layers are
reached as free functions exactly as item 2c specified and as `AMDraw` and `MmdChained` already
were. `Mmd1` and `Amd1` stay in the enum: they are ladder rungs, but they are also complete, correct
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
because a third B name is in circulation that is not an oracle at all: `AmdCompacted` is a VEHICLE,
holding one candidate change to `AmdFlat` at a time. It is SCRATCH and was never committed: it
existed
on 2026-08-10 while the two candidates below were being priced and was removed with them. Item 2c
says how the next one should be added so that it touches far less. Deciding the fate of the
other two while that name comes and goes would be churn in the middle of an open investigation.
Revisit when the amd branch is quiet.

**2c. HOW TO ADD THE NEXT VEHICLE.** `AmdCompacted` was wired in as a full ordering on 2026-08-10,
an
enumerator with everything that follows from one, and all of it was reverted a day later. One site
was missed on the way in: three files under `examples/` switch over `Ordering`, so `make examples`
warned. Do it as a FREE FUNCTION instead, which is the pattern `order_timing.cpp`'s `AMDraw` column
already uses on the stated grounds that an enumerator "would put a benchmark's oracle into the
library's public enum and into every switch over it". Two new files, two build entries, one column
in each benchmark driver calling `orderAmd3B` directly, and a local identity check against
`orderAmd3`. No enumerator, so no dispatch, no adapter, no `examples/` arm, no `test_order`
assertions, no `test_pipeline` sweep entry, and no move in `docs/TESTING_SPECIFICATION.md`.

**2d. THE NEXT PERFORMANCE PASS, and what of 2026-08-10 the other layers never got.** The walk
axis is FINISHED on `AmdFlat`: after the fold it walks `I[u]` twice per pivot and `A[u]` once, which
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
  footprint cost removed from `AmdFlat` that day. `Amd1B` is on record as "slower at large n after
  being faster at small", and that is the signature of precisely this. **But the fix may not
  transfer.** `AmdFlat` had somewhere to put the value, `partial[u]`, which exists only because
  ledger entry 4 split its bound in two; `Amd1B` and `Amd2B` form the bound in one pass and have
  no dead size-n array at that moment. `Amd2B` has `hashNext`, dead until filing, and `Amd1B`
  appears to have nothing. **Entry 4's split is now the third thing it has bought by accident**,
  after `AmdFlat`'s immunity to the hash-bucket order and its ability to take the key fusion at all.
- **The tagged W consolidation DID propagate, and this item is stale, 2026-08-17.** It said `Amd1`,
  `Amd1B`, `Amd2` and `Amd2B` all carry `outside(size)` plus a mark plus a clearing pass where
  `AmdFlat` carries one `w`. `Amd1` and `Amd2` both carry the tagged `w` and `wflg` today, and the
two
  B layers it names were retired on 2026-08-16. Finding this out is what made the supervariable
  stamp fold into `Amd2` a direct application rather than a prerequisite: the `w` it needed to
  stamp into was already there.
- **AND THE 2n MARK IS GONE FROM BOTH AMD DRIVERS THAT HAD ONE, 2026-08-17.** `AmdFlat` and `Amd2`
  now stamp supervariable detection into `w`, which is `AMD_2`'s `W [Iw [p]] = wflg`, so neither
  allocates a mark of its own and neither asks the shared class for clique marks. The
  `cliqueMarks` constructor argument and `cliqueBase()` went with them. `Amd1` has no detection at
  all and needed nothing.
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
`reachableSize`, `absorb`'s `vertexCount`, and one entity loop in `AmdFlat` brought back to the
`int32_t` form. `docs/DESIGN_DECISIONS.md` (2026-08-11) carries the account and
`docs/CODING_RULES.md` now states the three-way rule.

**Three things came out of it that the list below did not anticipate**, all in the design entry:
disjointness rather than arity is what bounds an accumulation, so nine sites needed no widening; a
narrowing cast goes outside an expression and a widening cast must go on an operand, which removed
the last arithmetic depending on the cap on n; and `bound` in four drivers was narrowed and had to
be reverted, being a fixed sum at its declaration and an accumulator three lines later.

**What remains `std::size_t` in the ordering is correct, and the list is short enough to check
against.** The `colPtr` parameters and the `cp` loop over them; `mSourcePtr` and `mCliquePtr`; the
five accumulators, `bound` in the four accumulating amd drivers, `deg` in `AmdFlat` and the hash
`key` in three, together with their `std::min<std::size_t>` calls and operand casts; `size()`;
`Buckets`'s constructor; the driver-local `size`; and `numFlagSweeps`, a diagnostic counter bounded
by the tag range rather than by n, which the dimension rule does not decide and which was left
deliberately.

**What is left of this item is everything outside the ordering.** The symbolic and numeric phases
are untouched and are the larger half. **And two casts hold n prisoner**, neither reachable and
both left deliberately:
`AmdFlat`'s `modulus = static_cast<std::int32_t>(size + 1)` fails at exactly the largest n
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
- `AmdFlat`'s `w`, `wflg`, `lemax`, `wbig`, signed by requirement rather than by sentinel; see
below.

`mEliminated` and `Mmd2`/`MmdFlat`'s `outmatched` are `std::uint8_t` and stay.

**Bucket 3, narrows to `std::uint32_t`.** In `QuotientGraph`, four arrays and one scalar, all
bounded by n: `mAdjacencySize`, `mIncidenceSize`, `mCliqueSize`, `mWeight`, and the per-pivot
`mCliqueWeight`. **`mWeight` IS NOW WRONG IN THIS LIST, 2026-08-17.** It goes back to
`std::int32_t` so that `AMD_2`'s `Nv` encoding can ride in its sign, which is the rule deriving
correctly rather than an exception to it: the unsigned case rests on having nothing to stand in for,
and that field now has. See `docs/DESIGN_DECISIONS.md` (2026-08-17) and the four conditions in
`docs/CODING_RULES.md`. The other four are unaffected. The accessors over them move too, `adjacencySize`, `incidenceSize`, `cliqueSize`,
`weight`, `cliqueWeight`, `reachableSize` and `reachableSetWeight`, which is what takes the casts
out of the hot loops.

In the drivers, **twenty-three declarations across the eight files**:

| array | where | what it holds |
|---|---|---|
| `degrees` | all eight | the filed degree, at most n |
| `outside` | `Amd1`, `Amd1B`, `Amd2`, `Amd2B` | per clique, `\|C[c] - C[p]\|` weighted |
| `cliqueDegree` | the five amd drivers other than `AmdFlat`, plus `AmdFlat` | per clique, `\|C[c]\|` weighted |
| `explicitPart` | `Amd1B`, `Amd2B` | per vertex, weight summed over the pruned `A[u]` |
| `partial` | `AmdFlat` | the half-formed bound, ledger entry 4 |
| `usedKeys` | `Amd2`, `Amd2B`, `AmdFlat` | hash values in [0, n]; a count by the range test, though it does not read as one |

`AmdFlat` has no `outside`, carrying `Amd.cpp`'s tagged `w` in its place, which is why the amd
column is not uniform. The two scan structs hold references to these same arrays, three in
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

**One quantity is deliberately SIGNED and must stay so**, and it is the only one. `AmdFlat`'s tagged
W
array with `wflg`, `wnvi`, `lemax` and `wbig`: `wnvi = wflg - weight(u)` is negative whenever the
tag is still small and the weight is not, which is the first few eliminations, and the arithmetic
only comes right at `w[c] - wflg`. Unsigned wraps there and the bound comes out enormous. `Amd.cpp`
is signed for the same reason. `src/AmdFlat.cpp` says this at its declaration; do not narrow it and
do not make it unsigned.

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

`amd3` is aligned: production `AmdFlat` returns `AMD_2`'s permutation exactly, up to the postorder
it
deliberately does not do, on the seven examples, 2D grids to 140, 3D grids to 24 and nine random
patterns. After the two fusions of 2026-08-10, the hash key into the bound and the first scan into
the prune, `AmdFlat` runs at roughly 1.0 to 1.25x the vendored routine on cubic grids from 12 to 32
a
side, overlapping it at 16, and 1.33 to about 1.95x on square grids from 32 to 400. `MmdFlat` is at
0.87 to 1.05x its own on cubes and 1.26 to 1.55x in 2D, with no trend in n on either. `AmdFlat` is
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
routine: `AMD1` costs about 1.2 to 1.8x on BOTH families, while `AmdFlat` goes from 2.3x in 2D to
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

Also needed: an overflow guard on `mTag`, as `AmdFlat`'s `w` array has with `clearFlag`, since
`mTag` only ever increments and `GONE` must stay above it.

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
  nonzeros, predicts **5.68 billion entries of fill under MmdFlat**, which is 42.3 GB of values and
  beyond the machine. AmdFlat will not rescue it: the two agree to within a few percent everywhere
we
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
