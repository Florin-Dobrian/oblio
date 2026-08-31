# AMD3: the alignment, iteration by iteration

A full record of how `amd3` was brought from `amd2` to return `AMD_2`'s permutation EXACTLY, up to
the postorder it deliberately does not do, on 2026-08-08. Written to the same shape as
`experiments/ordering/MMD3.md`, so a reader who has read that one knows where to look: what each step ESTABLISHED,
what it DISCOVERED, what was DECIDED and why, and the exact change made.

Companion documents:

- `experiments/ordering/MMD3.md` is the same record for the mmd branch, done the day before. Read it first; this
  document assumes its method and only notes where the two diverge.
- `experiments/ordering/README.md`, section "Aligning a layer against a vendored routine", is the
  method both alignments used. It was written as a handover note, `NEXT.md`, which was deleted once
  this work finished; the method survived into the README.
- `experiments/ordering/README.md`, section "amd3, and the second alignment ledger", is the durable
  record and holds the authoritative ledger.

---

## Read this first: what the whole thing came to

**Iterations 18 to 20 are a later session, 2026-08-09**, and they are the ones to read if the
question is whether the alignment holds: the acceptance test was widened from one shape to four,
and the widening found two defects and one bug. Two of its three initial failures were the harness
rather than the ordering, which is the caution to carry into any similar exercise.

**The alignment is iterations 0 to 9 and it succeeded.** `amd3` returns `AMD_2`'s permutation
exactly, up to the postorder it deliberately does not do, and production `Amd3` is extracted from
it. Six ledger entries, five of them tie-break conventions and one a real defect that had been
costing `Amd2` and `Amd2B` 3 to 9 percent of fill since they were written.

**Iterations 10 to 17 are a performance chase, and they are the part worth reading for anything
other than amd.** AMD3 went from 2.55x the vendored routine at 140 a side to 2.32x, and from 2.81x
to 2.32x at 400. MMD gained too, because half the changes are in the shared class.

### What actually worked, and where each came from

```
entry 6, the guaranteed match at the head of a short-circuiting walk   4.07x -> 2.55x   SOURCE VIEW
the W array, one tagged array where we had three                       2.60x -> 2.42x   READING Amd.cpp
the hoisted stamp, outer vertex stamped once per bucket not per pair    included above   READING Amd.cpp
the clique arena reserved, 18 reallocations per ordering removed        926 ms -> 111    SOURCE VIEW
mFiled deleted, 224054 byte writes per ordering that nothing read       both branches    COUNTING
the reach written straight into the arena, no scratch, no copy          111 ms -> 0      READING Amd.cpp
```

**Every one came from looking at something: the profile's source view, the vendored source, or a
counter.** Not one came from reasoning about our own code.

### What did not work, and this is the more useful list

```
a std::rotate in the prune, O(|A[u]|) per reached vertex     real extra work, worth nothing
the hash key fused into the bound loop                        a tenth of the visits, worth nothing
Amd3B, per-vertex arrays taken as raw pointers                the 1.12 s the trace blamed, zero
Amd1B and Amd2B, the first scan fused into the eliminator     built earlier, also zero
the integer width of the arrays                               10 to 20 percent, and MMD pays it too
folding mEliminated into the weight                           looked identical, broke mmd2
the read cursor hoisted, and the int32 loop counters          zero, measured
```

**Seven, and five of them were mine.** The pattern in the failures is the finding: everything that
scaled with *how many elements are walked* measured zero, which is what said the cost is what
happens *at* each element and how the loops around them are shaped.

### The structural answer

`AMD_2` never adds an array to record a fact it can encode in the range of one it already loads.

```
W    seen this step (>= wflg), absorbed (== 0), and the value |Le \ Lme|
Nv   live (> 0), taken into Lme this step (< 0), and the weight |nvi|
Pe   the list pointer, and FLIPped, the assembly-tree parent
```

We had `markMmd` plus `outside` plus deadness-by-removal where it has `W`; `mMarkMmd` plus
`mEliminated` plus `mWeight` where it has `Nv`. **Two of the three were portable and are taken. The
third is not**, and iteration 17's read establishes why: `Nv` is negated in the construct loop and
restored in the very last pass, so it is negative across the entire body of an elimination and four
separate readers are written to expect that. It works because `AMD_2` is one function. Ours is a
shared class with six drivers, three of them MMD with different invariants: and the `mmd2`
counterexample in iteration 17 is what that costs when a convention crosses the seam.

**So the remaining gap does NOT have a name, and an earlier draft of this line said it did.** It is
not algorithmic: the counts are equal: same eliminations, same reachable-set elements, same prune
elements, same pairs tested, same fill. It is also not the one-fact-per-array-against-three pattern,
which was the leading explanation for most of a day and was measured false: we touch 1.09x as many
arrays per element as `AMD_2` and take 2.32x as long. `notes/DESIGN_DECISIONS.md` (2026-08-08) has
the table and the reason the wrong conclusion was comfortable. What is left is per-TOUCH cost, which
is a locality hypothesis and is untested.

### If someone picks this up again

The profile is diffuse now. `orderAmd3` is 48 percent of the run with no line above 378 ms, and
everything with a name has been taken. Getting from 2.3x to 1.5x is a third of the total and there
is no item of that size left in it. It would have to come from a driver owning its own storage
instead of sharing `QuotientGraphFlat`: which is a real option, and is the trade `AMD_2` made, and
gives up the shared class, the six-driver ladder and the prototype-against-production check.

**And fix the benchmark first.** `order_timing` uses a fixed repeat count and two runs of the same
binary disagree by 4 percent. Five changes have now landed inside that band. `width.cpp` in
`experiments/ordering` picks its count from a timed probe targeting a fixed wall time per row; the
same twenty lines belong here.

---

**One difference from the mmd account, stated up front because its absence is information.** That
one has two wrong turns of the instructive kind, candidate fixes that scored better and were wrong.
This one has none. Every entry was read out of `AMD_2` before it was applied, and no change was
made and then reverted. What this account has instead is two errors of MINE, both in reasoning
rather than in code, both caught by the author rather than by a measurement, and both recorded in
iteration 8 rather than smoothed away.

---

## Iteration 0: the starting position, and a layer split that did not line up

**Established.** `amd2` has the two mechanisms `AMD_2` has beyond the bound, aggressive absorption
and hash supervariable detection, plus the liveness-aware core they require. It had been built one
pass at a time against the vendored source.

**And it returned a different permutation**, on all seven examples and every grid. That is the same
position `mmd2` was in, and the same problem: minimum degree is a tie-break algorithm, so a fill
comparison measures a difference of mechanism and a difference of arbitrary choice at once, with no
way to separate them.

**Decided: a new layer, and it took the digit 3 rather than 4.** The handover note expected the deliverable
to be called `amd4`, since `amd3` was taken by the layer holding dense rows, `amd_aat`, the
postorder and the Control interface. Instead THAT layer was renamed to `amd4` and the aligned one
took `amd3`, so the digit means the same thing on both branches: 3 is the layer aligned to the
vendored code. `amd4` is temporary, kept to read from, and forks sideways from `amd2` rather than
continuing the chain, which its header says.

The rename was its own step and its own commit, verified by the traces being byte-identical before
and after, since the layer's output carries no layer name. Doing it inside the first alignment
change would have put a mechanical rename and a ledger entry in one diff.

**Discovered, and it had to be settled before any comparison: the layer split does not line up.**
`AMD_2` is one function containing everything, including dense-row removal at its line 1321 and
`AMD_postorder` at 2364. `amd2` has neither, and both change the permutation. Aligning against
`amd_order` as it stands could not succeed, and the failure would have looked like a missing
mechanism.

The handover note named three ways out and recommended the first without concluding it.

**Decided, and it is option 1 taken only in part.** The oracle is `amd_order` on a scratch copy with
`Control[AMD_DENSE]` raised ABOVE `sqrt(n)`, which drives the dense threshold into its
`MIN(n, dense)` clamp so `deg > dense` is unreachable. Use that rather than the negative alpha the
header suggests, which leaves `dense = n - 2` and still removes a vertex adjacent to everything; at
the default neither fires here, the floor being `MAX(16, dense)` against a grid's maximum degree of
four. `Control[AMD_AGGRESSIVE]` stays at its default, aggressive absorption being a mechanism
`amd2` already has.

The postorder was NOT handled by adding one to the layer, which was the other half of option 1 and
would have given it a second job. See iteration 7.

**Also established, and it matters later.** The two knobs are the whole of `Control`. The three
functions `amd_order`, `AMD_1` and `AMD_2` are a decomposition rather than switches: entering lower
skips the input path and skips nothing inside the kernel. And `AMD_1`'s A+A' construction, which
Oblio does not need since its matrix is already full-symmetric, is not merely redundant: it decides
the ORDER of entries within each list in `Iw`, which is the content order everything downstream
walks. That is why the oracle is `amd_order` and not a hand-fed `AMD_2`.

---

## Iteration 1: rebuilding the tooling, and validating it against numbers we did not produce

**The situation.** None of mmd's instrumentation was committed; it lived in `/tmp` and died with
the session, and the method says to rebuild it early.

**What was built.** A scratch copy of `private/Amd.cpp` as `/tmp/probe/Amd_p.h`, never the original,
with six trace sites: the pivot, mass elimination, the hash merge, ordinary absorption, aggressive
absorption and `clear_flag`. `AMD_2` and its helpers are `static` in an anonymous namespace, so
including the source is the way in. A driver carries the seven examples and the grid builder,
transcribed vertex for vertex from `amd2.py` so the two sides cannot be diffing different problems.

**Validated on two numbers the probe could not have been fitted to.** `AMD_LNZ` plus `n`
reproduces the published vendored fill exactly, 206332 at 100 a side and 474995 at 140, the offset
being the diagonal `AMD_LNZ` does not count. And at 100x100 it reports 1 aggressive absorption and
2488 hash merges, which is the pair `experiments/ordering/README.md` records for the vendored
routine and our AMD2 alike.

**One reading corrected while it was fresh.** `clearflag` reads 1 at every size, which is the
priming call at line 1366 where `wflg < 2`, not an overflow. The two in-loop guards fired zero
times, which is the same claim our own `tag sweeps` counter makes.

**Established: the oracle is trustworthy and the grids match ours vertex for vertex.**

---

## Iteration 2: the smallest case, and the first entry

**Method.** Smallest divergence first. The seven examples run from `n = 4` to `n = 12`.

**The baseline, and the instrument choice that earned itself immediately:**

```
             pivot sequence      permutation
graph1           differ             differ
graph2           differ             differ
graph3           differ             differ
graph4           differ             differ
graph5           MATCH              differ        postorder only
graph6           differ             differ
graph7           MATCH              differ        postorder only
```

**Discovered.** On permutations the score is 0 of 7 and two of those are not divergences at all:
`graph5` and `graph7` have identical pivot sequences and differ only by the postorder. Comparing
permutations would have sent us hunting a mechanism on two graphs where nothing was wrong. The
score that means something is 2 of 7, which is exactly where `mmd2` started.

**Also discovered, and recorded rather than chased.** Three of the five real divergences share a
shape: `graph3` matches for nine pivots and we take a tenth, `graph4` for five and we take a sixth,
`graph6` for three and we take a fourth. A matching prefix and one extra pivot. One cause is
likelier than three.

**The case chosen: `graph1`, the 4-cycle, `n = 4`, diverging at the second of two pivots.**

```
              pivot sequence     hash merge      permutation      fill
vendored          3, 2            0 into 2         3 0 1 2          9
our amd2          3, 0            2 into 0         3 0 2 1          9
```

Same pair `{0,2}` merged, opposite survivor, identical fill. A pure tie-break.

**Root cause, from the vendored source.** `AMD_2` at line 1937:

```c
hval = hval % n ;
j = Head [hval] ;
if (j <= EMPTY) { Next [i] = FLIP (j) ;  Head [hval] = FLIP (i) ; }
else            { Next [i] = Last [j] ;  Last [j]    = i ;        }
```

Head insertion in both branches, and `i` is walked forward over `Lme` in scan 2, so the bucket comes
out reversed against `C[p]` and the pair loop takes the head. **The vendored survivor is the member
of `C[p]` seen LAST; ours was the one seen FIRST.** The handover note predicted this one by name.

**Change: walk the bucket backward.** The walk rather than the fill, which is the same order at O(1)
per insertion instead of O(bucket). `used_keys` keeps its forward order deliberately: `AMD_2`
reaches a bucket by rescanning `Lme` forward and taking the first member whose bucket is still full,
so buckets are PROCESSED in `C[p]` order while their CONTENTS are reversed. Two independent orders,
and only the second decides a merge.

**Result: 2 of 7, unchanged.** `graph1` closed and `graph7` opened.

---

## Iteration 3: a regression that was the information, and the deep entry

**The situation.** A correct change, read straight out of the vendored source, that moved the score
not at all.

**What it meant, and this is where the mmd account pays for itself.** `mmd3` stalled at 4 of 7 with
three walks reversed, because the CONTENT order of the lists was wrong upstream. The same shape
here would mean `C[pivot]` is ordered differently, and the new rule would then pick a different
survivor for a reason that has nothing to do with the rule.

**Confirmed from the oracle rather than inferred.** A trace of `Lme` was added to the probe. On
`graph7`, at the pivot where both sides merge the pair `{0,4}`:

```
vendored   LME 2 : 4 0        keeps 0, the last member
ours       C[2] = [0, 4]      keeps 4, the last member
```

The same rule, opposite outcomes, because the lists are ordered differently.

**Root cause.** `AMD_2`'s construct-new-element loop:

```c
for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++) {
    if (knt1 > elenme) { e = me ; pj = p ; ln = slenme ; }   /* the supervariables, LAST */
    else               { e = Iw [p++] ; pj = Pe [e] ; ln = Len [e] ; }  /* the ELEMENTS */
```

The elements first and the supervariables only on the final pass, so `Lme` comes out
cliques-then-explicit. `amd3_neighbors` walked `A[u]` then `I[u]`.

**Change: walk the cliques before the explicit adjacency.** Two loops swapped.

**Result: 4 of 7**, from 2. `graph1`, `graph2`, `graph5` and `graph7`.

**Established, and it is the mmd lesson reproduced exactly.** Entry 2 is the deep one and entry 1
cannot be judged without it: entry 2 fixes the content order of `C[pivot]`, which is the order
entry 1's buckets are built from. A change that scores nothing is not thereby a wrong change, and a
regression can be the measurement that locates the layer beneath.

**Worth noting where the two branches differ here.** genmmd puts variables first and elements last,
which is why the md ladder is laid out that way; AMD is the other way round. This is an amd
convention rather than a correction to the ladder, and `notes/DESIGN_DECISIONS.md` had already
recorded the layout difference from the other side without anyone noticing it mattered to an order.

---

## Iteration 4: one cause, two graphs

**The situation.** `graph3`, `graph4` and `graph6` still diverged, and two of them carried the
matching-prefix-plus-one-extra-pivot shape noticed in iteration 2.

**Smallest first: `graph6`, six vertices.**

```
vendored   PIVOT 5, PIVOT 1, PIVOT 3 with MASS 2, MASS 4, MASS 0     3 pivots
ours       pivots 5, 1, 3 merging 0 and 4 only, then a fourth pivot  4 pivots
```

The vendored routine mass-eliminates all three members of the new clique; we merged two and left
vertex 2 behind.

**Discovered, from our own state after that iteration.** `I[2]` ends as `{c3}` alone, so clique
`c5` WAS aggressively absorbed. Our mass-elimination test simply ran earlier, inside the eliminator,
while `c5` was still in `I[2]` and made `I[u] == {pivot}` false.

**Root cause, and the test itself is not the difference.** `AMD_2`'s test is
`Elen[i] == 1 && p3 == pn`, one element left and no supervariables, which is textually our
`not A[u] and I[u] == [pivot]`. What differs is WHERE it runs: `AMD_2` makes it in scan 2, over an
element list from which aggressive absorption has already dropped every clique lying inside the new
one. Its own comment says why that matters: *with aggressive absorption, `deg == 0` is identical to
the `Elen[i] == 1 && p3 == pn` test*.

**Change: mass elimination moves out of the eliminator and into the driver, after the absorption.**
A restructuring rather than a line, and the eliminator's `C[pivot] = reach(pivot)` identity is
restored as a side effect.

**Result: 6 of 7**, from 4. `graph4` and `graph6` closed together, which the shared shape predicted.

**Discovered while reading `AMD_2` for this, and it corrected one of our claims.** The experiment
README said the shrinking `degme` is vendored behavior, so a survivor handled early sees a larger
value than one handled late, and that our front-loaded mass elimination avoided a loss. `degme` IS
decremented inside scan 2, but it is never read there: the term enters a survivor's degree only in
the later pass that restores the degree lists, `deg = Degree[i] + degme - nvi`, by which point it
is final. Every survivor sees the same number in both codes and there was nothing to avoid.

---

## Iteration 5: the defect

**The situation.** One example left, `graph3`. The two sides agreed for eight pivots and then
parted, vendored taking 8 where we took 6.

**What the traces showed.** Both hash-merge 5 into 8 at the pivot before, so the mechanism fires
identically and the survivor is the same. What differed is the degree 8 is then filed at: vendored
files it at 2, we computed 3.

**Root cause.** `8` absorbs `5`, so its weight goes from 1 to 2. Our bound is
`explicit + degme - weight(u)`, capped by `num_left - weight(u)`, and it was written in the degree
pass which runs BEFORE supervariable detection, so `weight(u)` is the value `u` had before
absorbing `v`. `AMD_2` reads `nvi` in the pass that restores the degree lists, which runs AFTER
detection and after `Nv[i] += Nv[j]`:

```c
deg = Degree [i] + degme - nvi ;
deg = MIN (deg, nleft - nvi) ;
```

so its weight is the post-merge one. **A supervariable is filed one bucket too high per original
vertex absorbed**, and is never picked as early as its size has earned.

**Change: the degree pass stores `min(deg, degrees[u])`, which is `AMD_2`'s
`Degree[i] = MIN(Degree[i], deg)`, and a fourth pass after the hash finishes it with `+ degme`,
`- weight` and the `num_left` cap.** The exact-degree instrumentation moved with it, so it now sees
the same graph the bound does.

**Result: 7 of 7 examples.** `bound below exact` stayed 0, which is the invariant a wrong bound
would break.

**This is a DEFECT, not a convention, and it is mmd's entry 5 in a different array.** Nothing is
incorrect: an over-large bound is still a bound. What is wrong is that the code did not do what its
own comment said, which is the same test mmd's entry 5 failed. The comment argued for the behavior
it got wrong, and argued nearly correctly: an external degree does exclude `u`'s own supervariable,
and `u`'s reachable set really is unchanged. What it missed is that the buckets are keyed on a
degree that has the weight subtracted IN it.

**And `experiments/ordering/MMD3.md` had ruled this out in advance.** Its iteration 5, under "Where it could NOT
be", says the amd layers file at an external degree which does not move when a weight changes, and
cites `amd2`'s own comment in support. Both the claim and the comment cited were making the same
mistake. Had the amd bound been read at the time rather than reasoned about, this would have been
found there. That paragraph is now corrected.

**It was NOT fixed at this point.** The mmd work fixed `Mmd2` and production only after `mmd3` was
fully aligned, and the same order was kept: a divergence found mid-alignment is evidence about the
ledger, not yet a verdict about the layer below. It was parked in `notes/TODO.md` and taken up in
iteration 8.

---

## Iteration 6: the examples were never the test

**The situation.** 7 of 7 examples with four entries in, and grids still wrong. At most twelve
vertices each, the examples cannot exercise a mechanism that needs real structure.

```
side  3  4  5  6  7   agree
side  8                 DIFFER, fill 354 against our 355
side  9                 agree
side 10                 DIFFER, fill 648 against our 649
```

**Smallest first: grid 8, `n = 64`.** 42 pivots agree, then vendored takes 17 and we take 16. Both
are weight-2 supervariables at degree 6.

**Discovered.** They are the SAME supervariable: the pair `{16,17}` is merged on both sides, with
opposite survivors. So `C[pivot]`'s order still differs, more finely than entry 2 reaches. At the
pivot that merges them, and one before it:

```
vendored   LME 10 : 5 12 19 16 17 26        LME 2 : 5 12 19 10 16 17
ours       C[10] = [17, 19, 26, 16, 5, 12]  C[2]  = [16, 10, 17, 5, 12, 19]
```

The same runs in a different order of SOURCES, which is the element list's order rather than any
walk direction.

**Root cause, and reading it mattered.** The obvious reading of "the new element goes to the front"
is an insertion. `AMD_2` does something else, because its two lists live back to back in one run
and inserting at the front of the elements would mean shifting everything right:

```c
Iw [pn] = Iw [p3] ;   /* move first supervariable to end of list */
Iw [p3] = Iw [p1] ;   /* move first element to end of element part of list */
Iw [p1] = me ;        /* add new element, me, to front of list. */
```

It lifts the two entries at the boundaries to the two ends, so the elements come out
`[me, e2, ..., ek, e1]` and the supervariables `[j2, ..., jm, j1]`. **A rotation of each list, not a
shift of both.** We appended.

**Change: both rotations.** Both are load-bearing, because entry 2 made the reachable set walk the
cliques and then the adjacency, so both feed `C[pivot]`'s content order, which decides entry 1's
hash survivor.

**Result: fill exact at every grid size tested, and pivot sequences identical to grid 50.**

**Worth recording for the port.** Production is a better fit for this than the prototype:
`QuotientGraphFlat` already holds both lists in one run behind `mSourcePtr`, which is the layout
that motivated the trick, so AMD's three moves transcribe almost literally.

---

## Iteration 7: upgrading the acceptance test

**The situation.** Aligned on pivot sequences, which was the acceptance test settled in iteration 0
because `AMD_2` always postorders and a permutation comparison would fail on a relabeling.

**The question that changed it, and it came from outside the loop.** Whether the vendored order can
be seen at all without the postorder.

**Discovered: it can, and the check gets strictly stronger.** `AMD_2`'s pivot selection and its
supervariable finalization both sit inside the main loop, and `AMD_postorder` runs after it at line
2364. So the raw elimination order can be reconstructed upstream of the relabeling: track
membership alongside, a hash merge moving `j`'s members to `i` and a mass elimination moving `i`'s
to `me`, and emit each pivot's supervariable as its iteration closes. Concatenating those is what
the vendored routine would return if it stopped at the end of its main loop.

**Result: the acceptance test became the FULL PERMUTATION**, member order within each supervariable
included, and `amd3` passes it on the seven examples and on every square grid tested from 3 a side
to 40. That is a stronger claim than the one the layer was built against, and it is mmd3's own
acceptance test in a form the postorder cannot reach.

**Established, and it belongs in the method rather than in this ledger.** The acceptance test is
worth revisiting once the work is nearly done. It was chosen when the postorder looked immovable,
and the postorder never had to move: what had to move was where we read the oracle.

---

## Iteration 8: fixing the defect, and my two wrong turns

**With the layer aligned, the parked defect was taken up.** Measured first on a scratch copy of
`amd2` in `/tmp`, which touches nothing:

```
grid        AMD (vendored)    AMD1      AMD2 before    AMD2 after
 32x32            11900      12074         12364         11900
100x100          206332     201856        212496        199386
140x140          474995     455472        487111        444191
```

**Fixed at all four sites**, `amd2.py`, `amd2.cpp`, `src/Amd2.cpp` and `src/Amd2B.cpp`. One
subtraction rather than a recomputation, because all three terms of the bound's minimum shift by
the same amount when the weight grows. `Amd1` and `Amd1B` were checked and have no `merge` call at
all. The guard is `test_order`'s `Amd2B == Amd2` identity assertion, which is what says the B pair
moved together.

**It reverses the fill half of REPORT finding 3**, that AMD2's extras are a net loss with the hash
almost all of it. Corrected, AMD2 beats AMD1 at every size by 1 to 3 percent: the extras were not
costing fill, the filing was, and the hash was being charged for it. That is precisely the
attribution alignment was supposed to buy, and the handover note said to treat finding 3 as unresolved
until the permutations matched. The time half is untouched.

**And two wrong turns, both mine, both in reasoning.** Neither reached the code as a change that had
to be reverted, and both are here because the corrected text alone would hide that they happened.

**First: I claimed the postorder buys locality for the looking traversals.** It does not.
Correctness needs a TOPOLOGICAL order and nothing more, which holds for any permutation since
`ElmForestEngine` builds parent links from the permuted matrix and a parent's column is numbered
above its child's. Contiguity is a peak-memory property and only multifrontal has a peak; left- and
right-looking hold one update block at a time. And left-looking gains no locality from it either,
its relay hopping to whichever supernode owns a descendant's next remaining row rather than
climbing the tree.

**Second, and it followed from the first: I raised an open question that does not exist.** Having
noted that `setOptimizeMultifrontal` is off by default, I proposed measuring whether moving to a
non-postordering AMD would raise the multifrontal stack peak. It cannot: `DirectSolver` constructs
the forest engine with `mTraversal == Traversal::Multifrontal`, so choosing multifrontal turns the
option on, `labelDepthFirst` relabels the SUPERNODES into Liu's postorder, and the drivers loop
over those labels. The peak is set there and not by the column permutation. The engine default of
off reaches only a caller wiring the engines by hand.

**The shape both share.** I reasoned about a downstream phase from the ordering's point of view
rather than reading the phase. It is the same failure as the sentence in `experiments/ordering/MMD3.md` that ruled
out the amd branch, one level up: an argument that reaches confidently past where it was checked.

---

## What the alignment then made possible

**With the permutation identical, the remaining difference is implementation and nothing else.**
The counter comparison, which is the amd counterpart of the one in `experiments/ordering/MMD3.md`:

```
grid      pivots           hash merges      aggressive absorptions
          vendored/amd3    vendored/amd3    vendored/amd3
  16        200 / 200         53 / 53            1 / 1
  25        477 / 477        145 / 145           1 / 1
  32        777 / 777        244 / 244           1 / 1
  40       1209 / 1209       388 / 388           1 / 1
```

Identical throughout, and `tag sweeps` is 0 on our side at every size, which is the claim the
overflow guard exists to make checkable.

**Two checks the method names and this work should have used earlier.** `AMD_2` fills `Info[]`
natively, so several comparison points need no instrumentation at all. `AMD_LNZ` was used
throughout; the others were only read at the end:

```
grid 16    AMD_NDIV 1862     nnz(L) 2118  - n 256   = 1862    agrees
grid 32    AMD_NDIV 10876    nnz(L) 11900 - n 1024  = 10876   agrees
```

`AMD_NMULTSUBS_LDL` and `AMD_DMAX` are available on the same terms and were not used. They would
have given a second independent fill check from iteration 1 rather than from the end, and any
future alignment should read `Info[]` before writing a single `fprintf`.

**The strongest single check is against numbers this experiment did not produce.** `amd3`'s nnz(L)
is 206332 at 100 a side and 474995 at 140, which is what `benchmarks/ordering/README.md` records
for the vendored AMD, digit for digit.

**Final position.**

```
                  fill vs vendored AMD, grids
AMD1                 +1.1% to -4.1% in 2D, +11% to +14% in 3D
AMD2, fixed          -1.7% to -6.5% at the larger sizes
AMD3                  0.0% exactly, being its permutation
```

**And one result that goes the other way from mmd, which should not be smoothed over.** `amd3` is
aligned, so its fill IS the vendored routine's, 474995 at 140 a side, where the corrected `Amd2`
reaches 444191. On grids our tie-break now beats AMD's by 6.5 percent. Aligning MMD improved our
fill; aligning AMD costs it.

That is the second data point for the question `REPORT.md` parks, whether LIFO is
genuinely better or genmmd merely good, and it points the opposite way. Grids are one problem
family and the flattering one, so it wants the 3D grids `REPORT.md` has been asking for before it
means anything. Unlike mmd, there is therefore no case yet for `AMD3` becoming a default.

---

## The ledger, for reference

Append only. A row is never edited once closed.

```
#  what diverged                    where in ours       AMD_2                       nature
-  -------------------------------  ------------------  --------------------------  ----------
1  hash bucket walk                 the hash pass       the head push, line 1940    convention
2  reachable set layout,            amd3_neighbors      construct new element,      convention
   cliques before explicit                              the knt1 loop
3  mass elimination ran before      the eliminator,     scan 2, after the           convention
   aggressive absorption            now the driver      aggressive absorb
4  the vertex's own weight          the bound loop,     the fourth pass,            DEFECT
   subtracted before the hash       now a fourth        deg = Degree[i]
   merge that grows it              pass                + degme - nvi
5  the new clique appended to       the eliminator      Iw [p1] = me, and the       convention
   I[u] instead of prepended,                           two moves above it
   with a rotation
```

Six entries, thirteen working iterations, no rejected candidate fix, one real defect that had been
costing production fill since `Amd2` was written, one half-ported mechanism that cost 75 percent of
the layer's run time and nothing else, two claims of our own falsified, and three performance
hypotheses of the author's built or argued and falsified in turn.

**Four of the five are one idiom, and it is mmd's idiom.** AMD pushes at the head where we append,
so the entry seen last is processed first. Entry 1 is that in the hash buckets, entry 5 in a
variable's element list, entry 2 in which of two sources is walked first. Entry 3 is a placement
rather than a direction, but it is the same kind of thing: a step made in a different order.

---

## Iteration 9: production, and one self-inflicted wound

**Three of the six entries land in code the six drivers share**, so this was a port rather than a
copy. Entry 1 is `Amd3`'s own and entry 4 was already fixed everywhere.

**Decided: two flags and one method on `QuotientGraphFlat`, all inert for the other five drivers**, in
the shape `setReverseIncidence` already established for `Mmd3`.

```
setVendoredListOrder(bool)     entries 2 and 5, grouped because they are ONE fact, that AMD's
                               lists run the other way round, and are only ever wanted together
setLateMassElimination(bool)   entry 3, making eliminate() stop at the prune
massEliminate(pivot)           the half it no longer does, for the driver to call after absorb
```

Entry 3 got its own flag rather than being smuggled under a name about ordering, because it is a
different fact: the placement of a step rather than the direction of a list. A flag plus a method
was preferred to a third split of the eliminator, which already carries the
`beginElimination`/`finishElimination` seam built for the B variants.

**The cost, and it is real.** The shared class now carries two alignment flags and two parallel
methods, and a reader has to hold all of them. The alternative was `Amd3` owning its own reachable
set and eliminator, duplicating the two hottest functions in the ordering. Worth knowing before a
third alignment: two flags is fine and four would not be, and the right answer at that point is a
conventions object rather than more booleans.

**Verified inert before `Amd3` existed at all.** With the flags off, every suite passed and every
fill figure in `benchmarks/ordering` was identical digit for digit across all nine orderings. That
check is only available while the new driver does not yet exist, which is why it was taken then.

**And the self-inflicted wound, which grid 20 caught and nothing smaller could.** The rotation
helper for entry 5 was written as a `std::rotate`, giving `[pivot, c1, c2, ..., ck]`. AMD moves
the FIRST element to the end of the element part, so the wanted result is
`[pivot, c2, ..., ck, c1]`, which since the prune has already appended the pivot is a SWAP of the
first and last entries rather than a rotation. Only the two boundary entries move, which is the
whole reason `AMD_2` does it in three assignments instead of shifting.

The seven examples agreed, grid 10 agreed, and grid 20 differed by a single adjacent transposition
in 400 entries. `make test`'s prototype-against-production grid check is what found it, which is
the check added after two defects in `amd2` left all seven examples byte for byte identical. It
earned itself again here.

**Result.** Production `Amd3` reproduces the vendored routine's raw permutation directly, on the
seven examples and on grids at 20, 32 and 40, and agrees with the prototype everywhere `make test`
compares them. Its fill column in `benchmarks/ordering` is the vendored AMD's exactly: 11900 at 32
a side, 206332 at 100, 474995 at 140.

**`AMD3` is not the default, and that is a decision rather than an omission.** `MMD3` became one
because reproducing a reference with decades of use is a better bet on unseen inputs than a
tie-break of our own. The same argument would apply here but the evidence does not: with the
filing defect this alignment uncovered now fixed, `AMD2` fills LESS than the vendored routine on
grids, so `AMD3` fills more. Until that is tested outside 2D grids, both stay and neither is
preferred.

---

## Iteration 10: AMD3 measured slower than AMD2, and it was mine

**The observation, on alpamayo.** `make run` at 140 a side: AMD 1.23 ms, AMD2 3.06, AMD3 4.80.
Five entries that were all conventions about the ORDER of lists should cost nothing, so a 57 percent
penalty over AMD2 is not the alignment. It is an implementation defect.

**Root cause, and it is the same mistake this whole exercise exists to catch.**
`rotateForVendoredOrder` did a `std::rotate` over the adjacency, which is `O(|A[u]|)` per reached
vertex per elimination, on top of the prune that has just walked the same list. A whole extra pass
over the structure.

`AMD_2` pays nothing for it, and its three assignments say why:

```c
Iw [pn] = Iw [p3] ;   /* j1 copied one past the end */
Iw [p3] = Iw [p1] ;   /* e1 into j1's slot, so the supervariable list now STARTS one later */
Iw [p1] = me ;
```

It never moves a list. It writes two boundary entries and lets the list's start shift. I had
translated the EFFECT and not the MECHANISM, which is the rule the method states for candidate fixes
and which applies just as much to how a correct entry is implemented.

**The fix, and the half of it that does not work.** The adjacency reordering folds into the prune's
compaction: hold the first survivor in a register and append it last, `O(1)`. The incidence cannot
be done the same way, and finding out why is worth recording. Writing the pivot FIRST looks
symmetric and corrupts the run: the write cursor starts at `kept` and the read at the original
`adjacencySize`, and those are equal whenever nothing was pruned from `A[u]`, so an extra write
before the reads finish clobbers an unread entry. That is exactly why `AMD_2` makes its three
assignments AFTER both compactions rather than during them. The incidence therefore keeps its
append and takes an `O(1)` swap of the two boundary entries afterwards.

Both halves were tried; the suite caught the bad one immediately, 84 of 86 in `test_pipeline` and
the vendored comparison failing at every grid size.

**Result, on the Linux sandbox, which is not a measurement platform**, so read the ratio and not
the milliseconds. Against the vendored AMD at 140 a side, AMD3 went from clearly worse than AMD2 to
level with it, 1.90x against AMD2's 1.69x where the alpamayo run had 3.9x against 2.5x. It wants
re-measuring on alpamayo before anything is claimed.

**And it leaves the honest question about AMD3's cost still open**, since some gap over AMD2
remains. Two candidates, neither measured: `setVendoredListOrder` splits `reachableSet` into two
paths, so the branch predictor sees a shape AMD2 does not; and the fourth pass of entry 4 walks
`C[p]` once more than AMD2 does, which is real work rather than a convention and is the one entry
that is not free by construction.

---

## Iteration 11: the benchmark could not ask the question

`make scale` existed for the MMD branch alone, over the ladder to 400 a side, because that was the
branch being chased. There was no way to run the same comparison for AMD, which is what made the
observation above a `make run` reading at four sizes rather than a trend.

`order_timing.cpp` now takes `amd` beside `mmd`, and the gap columns are generated from whichever
method list is in force rather than being hardcoded to MMD1, MMD2 and MMD3. `make scale` runs both
branches; `make scale-mmd` and `make scale-amd` run one.

The AMD list is AMD, AMD1, AMD2, AMD3, with the vendored routine first since the gap columns take
the first entry as the baseline. **The B variants are deliberately absent**: they are their
originals' permutations on a different schedule, so their fill column carries no information, and
their time column belongs to the question about the seam rather than to the question about the
branch.

---

## Iteration 12: the profiler could not see AMD3 either, and it said nothing about it

`order_timing.cpp` gained `AMD3` when the driver was extracted; `order_profile.cpp` did not, and
nothing noticed. Its dispatch is a chain of `else if` on the method name with **no final else**, so
an unrecognized name ordered nothing and the program exited zero.

**The failure mode is the point.** The resulting trace is not empty and does not look wrong: it is
six milliseconds of dyld and process startup, with a plausible call tree under `main` and the
ordering simply absent. A second trace taken as the comparison showed the VENDORED AMD, since
`amd` is a name the driver does know. Two traces, both real, neither of the thing being
investigated.

That is the same shape as two defects already in this experiment's record: grid mode's whitelist,
where seven filter keys matched nothing and `nnz(L)` was never printed at all, and the
`make test` `grep` filter narrow enough to skip exactly the lines that had drifted. **An
instrument that silently declines to measure is worse than one that is absent**, because its
output is indistinguishable from a measurement.

`order_profile.cpp` now carries `amd3` and refuses an unknown method with the list of valid ones
and a nonzero exit, rather than falling through it.

---

## Iteration 13: the profile named one line, and it was half a mechanism

**Counting had run out.** Eliminations, reachable-set elements, prune elements, scan 1, the bound
pass, the hash key, the pair count, the bucket probes: all equal between AMD2 and AMD3 within a few
percent. The trace said the whole gap was in the driver's own inlined code, +5.70 s of self weight
against callees that were flat. So the source view, which is the only thing left that can attribute
inside one inlined symbol.

**It named one line, at 6.22 s of a 14.90 s run:**

```cpp
if (mark[incidenceU[i] + static_cast<std::int32_t>(size)] != other)
```

the incidence half of the hash pass's exact test.

**Why every count had missed it.** Both walks of that test carry `&& same` and exit on the first
mismatch, so what they cost depends on WHERE the mismatch is. My counter added
`adjacencySize + incidenceSize` for both vertices, which is what the test COULD cost. Equal
potential work, unequal real work. Counting the iterations actually executed:

```
       pairs tested   merges   adjacency iters   INCIDENCE iters   inc per pair
AMD2         227868     4885                 2            276296           1.21
AMD3         244334     4888                 2            507365           2.08
```

**Exactly one extra iteration per pair.** Which is what a guaranteed match at position zero
predicts, and entry 5 is what put one there: it moved the new clique to the FRONT of every `I[u]`,
and `u` and `v` are both members of `C[pivot]`, so both lists begin with the same entry and it can
never discriminate.

**And `Amd.cpp` does not pay it, because it skips that entry deliberately:**

```c
ok = (Len [j] == ln) && (Elen [j] == eln) ;
/* skip the first element in the list (me) */
for (p = Pe [j] + 1 ; ok && p <= Pe [j] + ln - 1 ; p++)
```

`Pe[j] + 1`, with the comment saying why. **That skip is part of the same mechanism as putting `me`
first.** Entry 5 ported one half of it.

**The fix, and it must be symmetric.** Skip index 0 in both walks. The two feed `sizeV` and `sizeU`
and those are compared at the end, so dropping the entry from one side alone would make every pair
fail on the count. Skipping both leaves the comparison exact, the sets differing by that one shared
element on each side.

```
after entry 6   AMD3   244334 pairs   4888 merges   263032 incidence iters   1.08 per pair
```

Below AMD2's 276296, with the pair count, the merge count and the permutation all unchanged, still
matching the vendored raw order at grids 20, 32 and 40.

**Measured on alpamayo afterwards, and it confirms the diagnosis exactly rather than roughly:**

```
                        AMD2      AMD3 before    AMD3 after
trace total            8.90 s        14.90 s        9.44 s
orderAmdN SELF         4.76 s        10.46 s        4.76 s
```

The driver's self weight is now identical to AMD2's to the millisecond. On the scale ladder against
the vendored AMD, AMD3 went from 4.07x to **2.55x** at 140 a side beside AMD2's 2.42x, and from
4.63x to 2.81x at 400 beside 2.58x. So the whole of the 75 percent penalty was that one line, and
what is left is 6 to 9 percent: AMD3 runs a fourth pass AMD2 does not have and tests 7 percent more
pairs, which between them is about the size of the residual. Not worth chasing without a reason.

**Entry 6 is a fourth nature and the ledger needed the word.** Not a convention, not a defect, not
cosmetic: it changes no ordering, no fill and no permutation, and it changes the COST. Call it
COST.

**And the lesson is about porting rather than about AMD.** Half a mechanism can be correct and
still be wrong. Entry 5 alone gives exactly the right answer and pays for it, and nothing in any
output could ever have shown that: the permutations matched, the fill matched, the counts matched.
Only the profile could, and only its source view.

**Three of my own falsified hypotheses got here**, and the sequence is worth keeping as a caution:
the `std::rotate`, which I fixed on a sandbox reading that was noise and which bought nothing;
`partial`, the size-n array I argued for twice from a precedent that fit its shape exactly and
which the trace never implicated; and the pair count, which I measured and which was equal. Each
was plausible, each had a mechanism, and the one that was right was invisible until a line-level
profile was pointed at it.

---

## Iteration 14: the branch gap, and the third walk

**The question that opened it.** MMD3 is 1.19x its vendored routine and AMD3 is 2.55x its own. The
two references cost the SAME, 1.31 ms against 1.27 ms at 140 a side, so it is not that AMD is
intrinsically harder. Ours are 1.56 and 3.25: our AMD is twice our MMD where the references are
level.

**Two answers, and only the second is actionable.**

The bound requires every member of `C[p]` refreshed at every elimination, because batching defeats
its tightness. 108705 refiles at 140 a side against MMD's 23861 degree updates, 4.5x more. That is
inherent and `AMD_2` pays it too, so it explains the absolute cost and not the ratio.

**The ratio is a difference of PASS STRUCTURE, and it is visible in the shape of the two vendored
routines.** `genmmd` is itself decomposed into `mmdelm` and `mmdupd`, so MMD3 mirroring it costs
little. `AMD_2` is one function that fuses everything into two walks of each list: scan 1 touches
element lists to update `W[e]`, and scan 2 touches both lists once and computes the degree,
compacts the pool AND accumulates the hash key, all in registers. Ours walked the same lists three
times: scan 1 over `I[u]`, the bound loop over both, then the hash key over both again.

**The B variants had already tried fusing a seam and got nothing**, which is worth knowing before
reaching for fusion again: `Amd1B` and `Amd2B` fold the driver's first scan into the eliminator and
are level with their originals. They fused the wrong seam. The duplication was never scan 1; it was
the key.

**The change.** Accumulate the key in the bound loop, which already walks both lists, and let the
reverse walk read it. That is `Amd.cpp`'s `hval += e` and `hval += j`, in the pass that was already
there. REPORT.md had specified this and measured the separate pass at 72 percent of AMD2's overhead
in 2D and 92 in 3D, and had never had it built.

```
AMD3 driver list-element visits, 140x140
before      600500     scan 1 + bound 300250, key 300250
after       541589
```

The reverse walk stays separate, because its DIRECTION is a tie-break: it decides which of two
indistinguishable vertices absorbs the other. Only the accumulation moved.

**Permutation untouched at grids 20, 32 and 40**, as it must be: the key is the same value
computed earlier: and both suites green. Applied to `Amd2` as well, where the 72 percent was
measured and which `Amd3` inherits from.

**The caution REPORT attached to it is now paid and not yet measured.** It needs an array of size
`n`, which is the `explicitPart` footprint trade that made `Amd1B` slower at large `n` after being
faster at small. Whether that eats the saving is the next `make scale-amd-2d` and nothing here
predicts it.

---

## Iteration 15: four fusions that bought nothing, and the two that did

**The question.** MMD3 runs at 1.2x its vendored routine and AMD3 at 2.55x, with the two references
costing the same. Where does our AMD spend what our MMD does not?

**Four things were built and measured, and four bought nothing.** They are recorded because each
was plausible, each had a mechanism, and the pattern in their failure is what led to the answer.

```
the hash key fused into the bound loop     a tenth of the driver's element visits removed   0%
Amd3B, per-vertex arrays as raw pointers   the 1.12 s the trace charged to accessors         0%
                                           built, measured, reverted, and deleted
Amd1B, Amd2B, the first scan fused         built earlier, same seam                          0%
the width of the arrays                    10 to 39%, and identical on the MMD branch      n/a
```

**The pattern is the finding.** Removing a tenth of the element visits changed nothing, so the cost
is not how many elements are walked. Hoisting the array bases changed nothing, so the accessors
were an artifact of where the inliner charged time. And width is paid equally by MMD, which is at
1.2x. Everything that scales with *elements touched* had been eliminated, which left what is
touched *at* each element and how the loops around them are shaped.

**Then two structural differences, both loop shape rather than any line, and both in `Amd.cpp` in
plain sight.**

**The W array carries three facts where we carried them in three places.** Ours: `mark[c]` for
seen-this-step, `outside[c]` for the value, deadness implied by removal from the lists, and a
clearing pass over a touched list at the end of every step. `Amd.cpp`: one tagged array, `0` for
absorbed, below `wflg` for stale, at or above `wflg` for seen with the value in the excess, and a
single `wflg += lemax` to invalidate the lot. Scan 1 becomes one load and one store into one array
where ours was two of each into two.

**The stamp is hoisted out of the pair loop.** `Amd.cpp` stamps the OUTER vertex once, before the
inner loop, and tests every candidate against that one stamp. We stamped the inner vertex once per
PAIR, over its whole list, with no short-circuit: so every pair paid a full list of random writes
even after entry 6 got the comparison down to 1.08 iterations.

```
140x140                  Amd3      Amd3C
stamp writes           639083     290473
compare iterations     263032     264510
ratio                     2.4x       1.1x
```

**And the reason we could not see the hoist is worth more than the hoist.** Our stamping carried
`w != u` and the walk carried `w == v`, exclusions that are pair-dependent and therefore appear to
pin the stamp inside the loop. They are VESTIGIAL: `u` and `v` are both members of `C[pivot]`, and
the prune drops every neighbour lying inside the new clique, so `A[u]` cannot contain `v`. There
was nothing to exclude, and `Amd.cpp` has no such guard for exactly that reason. Two dead
conditions had been holding a loop in the wrong shape.

**Both were landed as `Amd3C` first**, a re-schedule with `Amd3C == Amd3` asserted entry for entry
on the seven test matrices and on grids at 20, 32, 40, 64, 100 and 140. That identity is the whole
oracle a re-schedule has, and here it did real work: it is what says the guards really were
vestigial rather than merely looking so. **Once it had, `Amd3C` was folded into `Amd3` and both it
and `Amd3B` were deleted.** A variant that is exactly its parent and faster should be its parent,
and a variant that measured zero should not be a permanent column in every benchmark table.

**It is the first change in this sequence that moves on the Linux sandbox**, 1.98x against Amd3's
2.10x at 140 and 1.89x against 2.14x at 200, where `g++` at `-O3` had hidden every previous one.
The alpamayo number is the one of record.

---

## Iteration 16: the profile pointed at `operator new`

**Zooming, in decreasing order of weight**, which is the discipline this whole sequence should have
followed from the start. `orderAmd3` self 4.19 s, `eliminate` 1.87 s, and everything else together
1.6 s against a vendored kernel of 3.18 s.

**`orderAmd3` turned out to be diffuse.** Its heaviest own lines are 378, 316, 224, 204, 198 ms.
Entry 6 removed the one spike and what is left is spread thin: no line to fix.

**`eliminate` was not.** One line of its 1.87 s:

```
230    926.00 ms    beginElimination(pivot, inClique, absorbed);
```

and zooming again landed in `allocate.h`, on `__builtin_operator_new`. **The whole of it was
memory allocation**, not computation.

**What allocates.** `mCliqueArena` was never reserved. It is append-only: a new clique goes at the
end and a dead one's block is left where it lies: so it grows to the sum of |C[p]| over the whole
elimination, 108705 entries at 140 a side against nnz(A) of 97440. Unreserved, a vector reaches
that by doubling from nothing: 18 reallocations and 131071 entries copied per ordering, and the
last few blocks are large enough that the allocator serves them from mmap and each faults its pages
in on first touch.

`Amd.cpp` does not have this problem at all. Its `Iw` is one pool sized once and COMPACTED IN PLACE
when it fills, which is what `AMD_NCMPA` counts, and that counter is 1 for a whole 140x140 run.

**The change is one line**, `mCliqueArena.reserve(colPtr.back())`, and nnz(A) is a starting size
rather than a bound: it leaves at most one doubling on a 2D grid and a problem whose fill is many
times its input will still grow, now from a large base and so amortized. Reclaiming the dead blocks
the way `Amd.cpp` does is the real fix and is not this.

**It is in `QuotientGraphFlat`, so it reaches every driver**, and the MMD branch forms cliques the same
way. If it is worth what the trace says, MMD3 should move with AMD3.

**And the method note, which is the point of this iteration.** Four structural hypotheses of mine
were built and measured at zero. Two changes that worked came from reading `Amd.cpp`. This one came
from doing what should have been done first: open the profile, take the heaviest thing, zoom, and
zoom again when the answer is a call rather than a line. It took three zooms to get from a call
tree to `operator new`, and nothing about it could have been reasoned out.

---

## Iteration 17: the arena confirmed, and a substitution that looked identical and was not

**The reserve worked and both branches moved.** AMD3 from 2.37x to 2.31x at 140 and 2.30x to 2.22x
at 280; MMD3 from 1.24x to 1.16x at 140 and 1.22x to 1.01x at 100. `beginElimination`'s arena
`insert` went from 926 ms to 111 ms.

Worth recording that the MMD branch gained here at all. Its time gap had been decomposed and its
largest contributor identified back on 2026-08-07, and this found seven percent in it that has
nothing to do with ordering: the gap was measured with the arena doubling underneath it the whole
time.

**Zooming again, in decreasing order.** `beginElimination` is still the top line of `eliminate` at
868 ms, and inside it `reachableSet` is 554 ms and the arena insert now 111. Inside `reachableSet`,
one line:

```
104   235 ms   if (mMark[v] != mTag && (!live || mEliminated[v] == 0)) {
105     3 ms       mMark[v] = mTag;
106    59 ms       reached.push_back(v);
```

**Two random loads into two arrays.** `Amd.cpp` does it with one: `nvi = Nv[i]` and the sign says
whether the entry is still a vertex, so liveness and value share a load of an array the walk needs
anyway.

**The substitution looked exact and was not.** `merge()` zeroes the weight of precisely the vertex
it folds away, and this test exists for precisely those vertices, so `mEliminated[v] == 0` and
`mWeight[v] != 0` appear to be one fact. I checked the two states where they could part: a
mass-eliminated vertex, which is purged from the only clique that named it, and the pivot, whose
cliques all die in `beginElimination`: and concluded both were unreachable.

Under `mmd2` on a random 200-vertex pattern, vertex 152 comes out **eliminated with weight 1,
sitting in a live ADJACENCY list**, and the substitution emits it twice: 201 entries for 200
vertices. `mmd1`, `mmd3` and all three amd layers were unaffected.

**The unexamined premise was about the adjacency, not the cliques.** The prune removes the pivot
from `A[u]` for every `u` in `C[pivot]`, and that is not every list that can still be walked. I
reasoned carefully about the two cases I thought of and did not ask whether the case list was
complete.

**Reverted, with the counterexample recorded at the site**, because the shape is still right and
someone will try it again: folding liveness into the weight is what `Amd.cpp` does and it is the
hottest line in the ordering. It needs the invariant repaired, not asserted.

**And it is the fifth structural guess of mine to fail**, against two wins that came from reading
`Amd.cpp` and two from opening the profile. The pattern is now unambiguous enough to state as a
rule rather than a tendency, and it is written into the method notes above.

---

## Iteration 18: widening the acceptance test, and two of the three failures were the harness

**The situation, 2026-08-09.** `make amdorder` was committed running eleven 2D grids, which is ONE
SHAPE at many sizes. An attempt to widen it to the seven examples, 3D grids and random patterns
had found divergences and had NOT been committed, because a checker that fails for unknown reasons
is worse than one with narrow coverage. Three failures, and they turned out to have three
different causes, only one of which is ours.

**The random patterns: a SIZE MISMATCH, and the hook was not what was incomplete.** The vendored
side came up short, which says a vertex is numbered on a path the hook does not see. There are
four such paths, not three: empty variables in the initialization, mass elimination in scan 2, the
hash merge, and DENSE-ROW REMOVAL at `Amd.cpp:1387`, which the driver believed it had switched off.

It had not. `dense = alpha * sqrt((double) n)` assigns a double to an `Int`, so `1e30` overflows
and the conversion is undefined; on x86-64 it lands on `INT_MIN` and `MAX (16, dense)` then gives
SIXTEEN, which is dense removal fully on at the strictest setting the code can express. The
measurement that settled it in one line:

```
                        maxdeg    raw    deficit  ndense
alpha 1e30  2D grid 40       4   1600          0       0
            3D grid 12       6   1728          0       0
            random d6       21   1915         85      85
            random d12      39     19       1981    1981

alpha = n   random d12      39   2000          0       0
```

**The deficit EQUALS `ndense`, on every row.** That is the whole of it, and it explains the shape:
invisible on 2D grids at degree 4 and on 3D grids at degree 6, both under the threshold, and
appearing only on patterns dense enough to cross it. With a threshold that cannot overflow there is
no fourth path and the hook is complete.

Worth recording that the repair is not to teach the hook the dense path. If dense removal fires,
the oracle is running a mechanism `amd3` does not have, so the two are ordering different problems
and no member-order comparison means anything. The threshold now comes from `n`, and
`Info[AMD_NDENSE]` is read and fails the case, so a mis-set threshold names itself.

**The 3D grids, first cause: the builder emitted an invalid pattern.** Sizes matched there, so
every vertex was accounted for and only the order differed, first at position 32 of a 4x4x4 grid.
Tracing both sides showed it starts at the FIRST elimination. Pivot 63 is the far corner, and the
two refile its neighbours in opposite orders:

```
ours       REFILE 62, 59, 47
vendored   REFILE 47, 59, 62
```

which is the row order inside column 63. The builder pushed `x-1, x+1, y-1, y+1, z-1, z+1`,
descending for that vertex; our side took the pattern as given and the vendored side sorted it. A
column's rows must be ascending, which is the CSC precondition `SparseMatrix` states and which is
also a tie-break input, since the order within a column decides the content order of `C[pivot]`.
The 2D builder is ascending by construction, which is exactly why 2D could never expose it.

**The lesson is the one this ladder keeps relearning, from the other side.** Every earlier entry
was found by asking which line of the vendored routine we failed to reproduce. Here two of three
failures were not divergences at all, and the time went on the assumption that a red check means a
defect in the thing under test. A check is a program too.

---

## Iteration 19: the third cause was real, and it is ledger entry 7

**With the harness corrected, 3D grids matched from 2 a side through 12, and 16 and 24 still
failed.** A different signature: not order but VALUE. Several vertices filed one degree too high,
late in the run, at line 31097 of 33381.

**The cap was the clue.** Some vertices in the same pass agreed and some were plus one, which is
what `deg = MIN (deg, nleft - nvi)` binding looks like, and a cap binds only when a bound
approaches what is left, which is why it surfaces at 88 percent of the run. But `numLeft` and
`degme` agreed at every iteration, so it was not the cap's operands.

**One vertex, all the way through.** Vertex 982 has fourteen bound computations. Thirteen agree
exactly; the fourteenth has our freshly accumulated degree at 40 against 39, at the point where 982
has weight 6. Splitting that number into its halves named a single clique:

```
ours       clique 1509   dext 40
vendored   element 1509  dext 39
```

Both obtain `dext` by subtraction from a maintained clique degree, so the difference is in that
degree. `Amd.cpp` writes `Degree [me] = degme` at its line 1676 AND AGAIN at 1940, after scan 2 has
run `degme -= nvi` for every mass-eliminated vertex; the second write is what every later step
reads as `|C[me]|`. We wrote it once, with the pre-merge value.

**Half a mechanism, again, and from ledger entry 3.** That entry moved mass elimination out of the
eliminator, which is what makes the second write necessary, and did not carry it. `Amd1` and
`Amd2` mass-eliminate inside the eliminator, so their single write already sees a trimmed clique
and they are unaffected. The prototypes cannot have it at all, obtaining the quantity by walking
`C[c]` rather than maintaining anything.

**So the twin check was structurally incapable of catching this**, which is worth separating from
"the cases were too small". A prototype written to read as the algorithm does not carry the
optimization, so it cannot model a hazard that lives in one. `REPORT.md` parked that divergence as
its fifth lead and could not decide whether it mattered. It does.

One line, and the widened check goes to 38 of 38. Mutation-tested both ways: with the write
removed, `grid3d 16^3` and `24^3` fail and nothing else does; with the threshold put back to
`1e30`, six random cases report the dense removal instead of a mismatch.

---

## Iteration 20: and the acceptance test found a bug that is not the ordering's

**`make amdorder` failed on alpamayo at `grid3d 6^3` and passed on the Linux sandbox.** Same
source, same deterministic input, no floating point anywhere in an ordering, so a difference
between two machines cannot be a legitimate result. Something unspecified was being read.

**AddressSanitizer named it in one run**, which is the right instrument for "two machines disagree
about integer code" and was reached for instead of a sixth hypothesis:

```
heap-use-after-free   QuotientGraphFlat::reachableSet   QuotientGraphFlat.cpp:122
  freed by:  vector<int>::_M_realloc_insert        a push_back that outgrew its reserve
```

`beginElimination` calls `reachableSet(pivot, mCliqueArena)`, so the walk appends to the very
buffer it is reading clique members from, through `mCliqueArena.data() + mCliquePtr[c]`. A
`push_back` past the capacity moves the arena and that pointer dangles for the rest of its clique.
The constructor reserves `nnz(A)` and the arena grows to the sum of `|C[p]|`, 108705 against 97440
at 140 a side, so a reallocation is ordinary rather than exceptional.

**It is every driver, not `amd3`**: mmd1, mmd2, mmd3, amd1, amd2 and amd3 all trip it, on 2D grids
as well as 3D. And it is normally harmless, which is why it survived: a vector growth copies and
then frees, so the stale pointer usually still finds the right values sitting in freed memory. It
stops being harmless when the allocator recycles that block, which is what one machine did and the
other did not.

**Where it came from is iteration 16 of this document.** That change wrote the reach straight into
the arena, for a measured 111 ms, and nothing re-checked the pointers already being held across the
append. `beginElimination` even carries a comment saying the reach pointer is taken after the
append "since that is what can move the arena": the hazard was seen there, for the one pointer it
happened to be about, one line later than the one that mattered.

The repair is to make the arena unable to move rather than to re-fetch per element, which would put
a load in the innermost loop of the whole ordering for a hazard that occurs at most once per
elimination. A reach is at most `n` entries, so guaranteeing room for one before the walk
guarantees it for the walk. `notes/DESIGN_DECISIONS.md` (2026-08-09) carries the entry, since it
belongs to the shared class rather than to the alignment.

**Three method notes, and the third is the one that generalizes.** An acceptance test is worth
widening even when it is passing, because what it cannot reach it cannot check. A disagreement
between two machines on integer code is a sanitizer question and not a reasoning question. And a
performance change that alters WHERE something is stored has to be read for lifetime as well as for
speed, because the profile cannot see a pointer that is about to be invalidated.

---

## Iteration 21: counting the pair loop, because the families disagreed about the extras

**The situation, 2026-08-09, later the same day.** The cubic benchmark had split each branch into
its base and its extras and found `AMD1` flat at 1.2 to 1.8x the vendored routine on both families
while `AMD3` went 2.3x to 3.0x. So the degradation was entirely in aggressive absorption and hash
supervariable detection, and the hash pass is the only part of them whose cost is superlinear in
clique size, its pair loop being the sum of squared bucket sizes over `C[p]`.

**A COUNT and not a profile, which is worth stating because this folder's own advice is to
profile.** The question was whether we do more work, not whether we do it less efficiently, and a
count answers that machine independently and reproduces to the digit. It also priced two things
apart that a timing cannot: how many pairs are enumerated, and what each one costs.

**Measured on scratch copies of `src/Amd2.cpp` and `src/Amd3.cpp`**, entry points renamed, counters
added, and every row checked against the real driver's permutation before it was believed. At
comparable n, `AMD3`:

```
                        pivots   cand/pivot   pairs/pivot   pairs/merge   cmp/pair
2D 140x140, n=19600      14709         6.28          19.0          57.3      1.083
3D 26^3,    n=17576      11863        12.72         155.3         324.0      1.095
```

Clique size doubles and the pair count goes up 8.2x, so it is quadratic in bucket size rather than
proportional to it: `sum |B|^2` grows 5.9x where `sum |B|` grows 1.6x. **The cost per pair did not
move at all**, 1.083 against 1.095, so the exact test still dies on its first iteration and the
growth is entirely in the count. And the pass is now dominated by it, 63 percent of its element
work on cubes against 34 in 2D.

**One thing fell out unasked.** `Amd2` stamps per pair where `Amd3` hoists the stamp, and on cubes
that is 495 stamp writes per pivot against 40.7, a factor of 12, where in 2D it is 3.2. So
iteration 15's hoist is worth four times more on this family, which is a mechanism for the thing
the benchmark had reported and not explained, that `AMD2` is slower than `AMD3` on cubes while the
two are level on squares.

---

## Iteration 22: the control, and it was not the reading either branch expected

**What the count could not say.** `AMD_2` runs the same mechanism with the same kind of key and
does not degrade on cubes. So a growing pair count might be something both codes have, and the two
readings called for opposite work: if its count grows too, the difference is the cost per pair; if
it stays flat, our key spreads worse and the filter is the fix.

**The same counters, in a scratch copy of `private/Amd.cpp`**, generated by a script that asserts
every anchor, the same arrangement as `tools/hook_amd.py` and for the same reason. Ten minutes,
because the technique was already written down.

```
                       pivots   cand/pivot   pairs/pivot   merges   pairs/merge
2D 140x140
  vendored AMD_2        14710         6.28         0.333     4890          1.00
  ours, AMD3            14709         6.28        19.034     4888         57.28
3D 26^3
  vendored AMD_2        11886        12.79         0.484     5690          1.01
  ours, AMD3            11863        12.72       155.335     5688        323.97
```

**Same candidates, same merges, 57 times the pairs on squares and 320 on cubes.** Its buckets are
effectively singletons: essentially every pair it tests is a merge. That is not a cost difference,
it is a different bucket partition, and it sent the next look at the key rather than at the loop.

**Recorded because the technique is not new and the use of it is.** `MMD3.md` entries 5 and 6 came
from two `fprintf` calls in a scratch copy of `Mmd.cpp` and iteration 1 of this document built a
six-site probe of `Amd.cpp`. Every one of those asked which line we fail to reproduce. **None asked
how much more of it we do.**

---

## Iteration 23: the key, read beside `AMD_2`'s

**Ours:**

```cpp
key  = sum over live v in A[u] of (v + 1)
     + sum over c in I[u] of (c + 1) * (size + 1);
hash = key % (size + 1);
```

The stride and the modulus are the same number, so the incidence term is annihilated exactly and
**the bucket is a function of the adjacency alone**. `A[u]` empties as the elimination proceeds and
everything a vertex reaches becomes cliques, so the surviving key carries less and less, and a
cubic grid reaches that state sooner than a square one. `AMD_2` sums element and variable ids into
one running value with no stride and takes it mod n, so both halves reach its bucket.

**Measured without changing anything**, by computing the alternative bucket alongside the real one:

```
                    pairs/pivot   altpairs/pivot   maxBucket   alt maxBucket
2D 140x140               19.034            0.481          20               3
3D 26^3                 155.335            0.592         110               6
```

which is the vendored routine's regime. **The stride's argument is correct about the key and says
nothing about the bucket**, and the invariant the two lines have to hold together is that the
modulus must not divide the stride. That is ledger entry 8.

**And the experiment README stated the whole thing as a virtue.** It explained the key as a number
in base `n + 1` with the adjacency sum as the low digit and the incidence sum as the high one,
three lines below defining the hash as the key modulo `n + 1`. Reducing a base-`n + 1` number
modulo `n + 1` keeps the low digit. Both sentences are true; nobody had put them side by side.

---

## Iteration 24: what the fix is worth, and the one place it moves the permutation

**Two lines, nine files**, the stride dropped everywhere the hash pass lives and, in `Amd3` alone,
`AMD_2`'s guard on the same loop: it enters only for a bucket member with a successor,
`while (i != EMPTY && Next [i] != EMPTY)`, where our hoisted stamp otherwise pays a full list of
random writes for a pair that will never be tested. `Amd2` and `Amd2B` need no counterpart, since
they stamp inside the pair loop.

**alpamayo, `make scale3d` and `make scale2d`, against the figures of earlier the same day:**

```
cubic grids        AMD      AMD1      AMD2      AMD3
26^3   before     4.07      5.76     14.88     12.30
26^3   after      4.06      5.69      5.45      5.83
32^3   before     8.81     13.03     26.4*     25.8*
32^3   after      8.81     13.03     12.33     12.71
```

The vendored routine and `AMD1` are unmoved to within a percent, which is the drift and is what
says the rest is real. `AMD3` on cubes goes from about 3.0x to 1.44x and `AMD2` from 2.85x to
1.40x; in 2D the same change is worth about 1.4x, `AMD3` at 140 going from 2.55x to 1.81x, which is
the ratio the count predicted from 155 pairs per pivot against 19.

**`REPORT.md` finding 3 is dead on both halves.** On cubes `AMD2` is now FASTER than `AMD1`, 5.45
against 5.69 at 26 a side, so carrying the extras is cheaper than not carrying them. The fill half
went on 2026-08-08 with entry 4. Both halves were measurements of defects rather than of
mechanisms.

**And the families have swapped roles**, which is where the next look should start: the extras were
free in 2D and ruinous on cubes, and are now free on cubes and cost about 25 percent in 2D.

**`Amd3` keeps its permutation and `Amd2` does not.** `make amdorder` matches on all 38 cases and
the two orders are identical on 31 shapes checked directly; `Amd2` and `Amd2B` move, with fill
`+1.4` percent at 140x140 and `-3.1` at 26^3, two-sided and small. The reason is structural and
worth keeping: `Amd2` files during the bound pass and the hash merge's refile is the last write, so
the order in which hash buckets are processed reaches the degree buckets, while `Amd3` refiles
every survivor afterwards in `pivotClique` order and comes out canonical whatever the buckets were.
That fourth pass is entry 4's, built for the post-merge weight, and it made `Amd3` immune to this
by accident.

**Verified:** 283/283 in the suite, `make amdorder` and `make mmdorder` 38 of 38 each, every twin
and every `PORTED` layer agreeing, and `hash pairs tested` added to all three amd layers as the
standing witness, reading 94 against 88 merges at side 20.

---

## Iteration 25: the per-pass inventory, and the fifth fusion is the one that worked

**The situation, 2026-08-10.** Entry 8 had left `AMD3` at about 1.6x the vendored routine's useful
cycles on both families, and `NEXT.md` named a per-pass inventory of element visits as the next
instrument. Half of it existed. What was missing was the vendored half, so the table was our
numbers rather than an attribution, and the claim that `AMD_2` walks a vertex's element list twice
where we walk it three times was a reading of the source.

**The vendored half, by the same technique as `tools/hook_amd.py`**, one counter per loop, every
anchor asserted, the counted copy's permutation checked against the unhooked `amd_order` and
`AMD_LNZ` matched against six published fill figures before any number was read off it. Two of our
own errors were caught by guards rather than by reading: `AMD_NDENSE` transcribed as index 3, which
is `AMD_SYMMETRY`, and a zeroed `Control` array, which is not the default `Control` and silently
turns aggressive absorption off.

**And our half was wrong, which the vendored half is what found.** `QuotientGraphI.cpp` carries the
same symbols as `QuotientGraphFlat.cpp`, so it is linked instead of it and both drivers share it,
and the probe read its counters after a control run of the uninstrumented driver. Six of fifteen
columns were doubled and nine were not. The check was free and existed only because the vendored
numbers were beside it: the alignment forces our `reachAdj` to equal `AMD_2`'s construct-adjacency
count, and after the correction it does, 2.61 against 2.61 in 2D and 4.18 against 4.18 on cubes,
with four further passes agreeing to the digit at all six sizes.

```
                              per pivot, corrected     as recorded
AMD3 element visits, 140x140        112.26                155.15
AMD3 element visits, 26^3           276.25                374.83
```

**What the completed inventory says.** Two axes, and both are family independent to two digits,
which nothing on this question had been before:

```
                          ours                 AMD_2             ratio
               walks sweeps    all   walks sweeps    all  walks  sweep    all  cycles
2D 140x140     93.41  56.55 149.96   45.25  25.32  70.57  2.06x  2.23x  2.12x   1.56x
3D 26^3       238.08 114.49 352.57  115.41  51.27 166.68  2.06x  2.23x  2.12x   1.61x
```

We make nine sweeps over `C[p]` per pivot and `AMD_2` makes four; we walk `I[u]` four times, in the
prune, scan 1, the bound and the key, where it walks twice. Five of ours collapse into its scan 2,
which per member computes the degree, accumulates the key, compacts the list and tests for mass
elimination in one visit. And 2.12x of visits against 1.56x of cycles means **we already execute
about 0.74x the work per visit that `AMD_2` does**: its loops are fat because each does four
things, ours are thin because each does one.

**The reading that follows from that, and it was wrong.** Fusion does not remove work, it trades a
visit for per-visit work, so a reschedule pays only where the per-visit work is redundant rather
than merely distributed. Four fusion attempts had measured zero, which looked like the evidence.
Deleting a sweep, by contrast, removes `|C[p]|` scattered loads and a loop and makes nothing else
more expensive. So the sweeps axis was ranked first.

**Two candidates, put through a scratch `Amd3B` one at a time.** That driver is a VEHICLE and not
a layer: it carries exactly one change, is priced against `Amd3` by the benchmark, and goes away
once the change has landed or been recorded. **It was never committed**, so nothing below it left
an enumerator, an assertion count or a benchmark column behind. `NEXT.md` item 2c says how to
build the next one.

- **B1, a sweep deleted.** `Amd3` re-sums `weight(u)` over `C[p]` after mass elimination trims it.
  `massEliminate` already maintains that number, decrementing `mCliqueWeight` per vertex it takes,
  so `cliqueWeight()` returns what the loop recomputes. Provable rather than empirical.
- **B2, the key folded into the bound.** The bound walks `A[u]` and `I[u]`; the key walked both
  again. `Amd.cpp` accumulates `hval += e` and `hval += j` inside its scan 2.

**Carried together they measured 5 percent in 2D across six sizes and 5 to 12 on cubes.** Split,
on alpamayo, `make run2d`, `make run3d`, `make scale2d`, `make scale3d`:

```
AMD3 -> AMD3B         B1 alone        B2 alone
64x64                   0.0%            -7.0%
100x100                 0.0%            -4.6%
140x140                -0.9%            -6.5%
200x200                +3.6%            -5.1%
280x280                -1.2%            -4.4%
400x400                +0.5%            -5.4%
12^3                   +5.0%           -14.3%
16^3                   +2.8%           -12.8%
20^3                   +1.5%            -8.1%
26^3                   +2.3%            -5.0%
32^3                   +1.4%            -8.0%
```

**The fusion carried all of it and the deletion carried none.** `nnz(L)` identical at every
size in every run. `AMD3B` at 32^3 came out at 11.20 ms against `AMD2`'s 11.70 while returning the
vendored permutation, which is a cleaner statement than any ratio here since both are ours and both
ran in one process.

**Why B2 worked this time, having failed on 2026-08-08.** That version carried the key in a vector
of size n and measured nothing at 140 a side and minus two percent at 400, which `REPORT.md` had
already named as the footprint trade, the same stream that made `Amd1B` slower at large n after
being faster at small. This one files each vertex into its bucket at the point its key completes
and stores nothing extra, `hashNext` being size n either way. **The failure was the array, not the
fusion**, and the reason that distinction went unexamined for two days is that four other nulls
made the fusion itself the obvious suspect. It was also measured while entry 8 was live, when the
exact comparison ran 19.0 pairs per pivot against the vendored routine's 0.33, so the pass it
shortens was not the one the profile was standing on.

**B1 is a null, and calling it a regression was the first draft's error.** It reads 0 percent in 2D
and 1 to 3 percent slower on cubes, and the whole of that range is inside this benchmark's floor,
which was measured the same day at plus or minus 3 percent: between the fusion landing and the
vehicle being removed, `Amd3B` held a verbatim copy of `Amd3`, so the two benchmark columns were
the same code timed twice in one run and the difference between them is the instrument measuring
itself. Its sign was positive at all five
cubic sizes in two independent runs, which is weak evidence of a small real cost and not more. One
mechanism would explain it, the deleted sweep having walked `pivotClique` immediately before the
bound sweep walks the same array and so warming `C[p]`, but that is a hypothesis attached to a null
and CPU Counters would settle it if it ever mattered.

What it does establish needs no mechanism: **deleting a provably redundant sweep bought nothing
measurable where fusing two walks bought 4 to 7 percent.** That is enough to reorder the queue, and
it means the stamp and the mass-elimination sweep no longer inherit an argument from their shape.

**What that pair of results does to the ranking.** The argument for the sweeps axis rested on four
fusion nulls looking like a rule. The fifth fusion is the only thing that has ever moved this gap
and the one deletion was negative, so the rule was not one. What generalizes instead is narrower
and more useful: **a null is a measurement of one implementation, not of the idea it implements**,
and the four earlier nulls were read as the latter.

**Two things the split bought that a combined measurement could not.** Carried together the pair
would have been recorded as two improvements when one of them measured nothing, and the remaining
sweeps would have been ranked above the remaining walks on the strength of a number the deletion
contributed nothing to. That is what
`Amd3B` is for, and it is why it holds one candidate at a time.

**And leaving it briefly empty was worth more than the split.** With the fusion landed and the
next candidate not yet chosen, `Amd3B` was a verbatim copy of `Amd3`, so the benchmark ran the same
code in two columns and reported the difference: plus or minus 3 percent, on that machine, in that
run, at every size. That is an error bar under identical conditions, which no amount of repeating
one column supplies, and it cost one column. Worth arranging deliberately whenever a vehicle
exists and is idle. It also retired a claim made earlier the same day, that B1 was a regression.

**Landed in `Amd3` on 2026-08-10, and the vehicle went with it.** It had been wired in as a full
ordering, an enumerator with everything that follows from one, which touched ten files and missed
an eleventh, `make examples` warning that three `Ordering` switches did not handle it. The next
vehicle is a free function called from the two benchmark drivers instead; `NEXT.md` item 2c has
the list.
**`Amd2` and `Amd2B` cannot take it the cheap way, and the reason is ledger entry 4.** They form
the bound in ONE pass and call `buckets.refile` inside it, so their bound loop's direction is
already a tie-break input, deciding which vertex sits at a degree bucket's head. Their key pass
walks `C[p]` backward against that forward bound, and HEAD insertion into both structures wants
opposite directions, so one walk cannot serve both. Measured on a scratch copy: fusing there the
way `Amd3` does changes the permutation on all ten grids tried, with fill moving `-1.33` percent at
140x140 and `+3.24` at 26^3, two-sided and small. So it is an ordering change there, not a schedule
change.

**One route is left and was not tried.** TAIL insertion from a forward walk reproduces the chain
order their backward walk produces now, so the permutation would be preserved and it would be a
genuine schedule change. It costs a `hashTail` array beside `hashHead`, which is the same size-n
footprint that made the 2026-08-08 version of this fusion measure nothing, touched at a similar
rate, so the expectation is that it does not pay. That is an expectation and not a measurement.
Left undone deliberately: `Amd3` is the production default, `Amd2`'s speed has no consumer, and the
profile that would size this whole direction has not been taken yet. Entry 4 split `Amd3`'s bound in two and moved the refile below the hash, for the
post-merge weight, and that is what leaves its bound loop free of tie-break duty. **The second
thing that split has bought by accident**, the first being `Amd3`'s immunity to the hash-bucket
order reaching the degree buckets, recorded in iteration 24.

`Amd1` and `Amd1B` have no hash detection and so no key to fuse. The prototypes deliberately do not
take it either, since an oracle that shares an optimization cannot see a defect in it, which is
entry 8's lesson exactly.

**Verified:** `make amdorder`, `make mmdorder` and `make aligned` 38 of 38 each, the whole test
suite green, and `AMD3B == AMD3` on 32 graphs including twelve random ones.

**What the inventory still cannot say**, and what a profile would. B2 removes 18 to 19 percent of
the visits and bought 5 percent in 2D. Whether the saving is instructions, in proportion to the
visits, or locality from touching each list once instead of twice, decides whether the three
remaining walks of `I[u]` are the next large prize or a much smaller one.

---

## Iteration 26: the last two walks, and three things that were not the algorithm

**The situation, 2026-08-10, after the key fusion.** The corrected inventory had `AMD3` at 1.64x
`AMD_2`'s element visits on cubes and about 1.30x its work, with efficiency at parity on that
family, so the cubic gap was pure work and the walk count was the whole lever. The remaining excess
was concentrated: we walked `I[u]` three times per pivot where `AMD_2` walks it twice, and `A[u]`
twice where it walks it once, and those two rows were 96 percent of what was left.

**The transformation for exactly that already existed and had already failed.**
`QuotientGraphFlat::eliminate(pivot, ApproximateScan&)` folds the driver's first scan into the
prune, which is what `Amd1B` and `Amd2B` are. It measured zero on both, five percent slower on
`Amd1B`. Three things made it worth re-running: it had never been tried on `Amd3`; the reading was
2D, and the direction is now 3D; and it predates entry 8, when a different pass dominated the
profile. That is the same pair of conditions that had made the key fusion look worthless the day
before.

**It could not be reused as it stood, and the detour was the right one.** That overload carries the
pre-iteration-15 encoding, a value array plus a separate seen-this-step mark, where `Amd3` carries
`Amd.cpp`'s tagged W. Adopting it would have bundled a revert of that consolidation into the
measurement, and the record prices W only together with the stamp hoist. So `QuotientGraphFlat`
gained a `TaggedScan` overload and the vehicle differed from `Amd3` in the fold alone.

**Result, alpamayo, and the vendored routine is unmoved throughout:**

```
                  AMD3   AMD3B                    AMD3   AMD3B
12^3              0.29    0.26  -10.3%   64x64     0.41    0.41    0.0%
16^3              0.93    0.85   -8.6%   100x100   1.03    0.95   -7.8%
20^3              2.22    1.94  -12.6%   140x140   2.07    1.98   -4.3%
26^3              5.28    4.67  -11.6%   200x200   4.27    4.14   -3.0%
32^3             10.85   10.07   -7.2%   280x280   8.87    8.76   -1.2%
                                         400x400  20.12   19.47   -3.2%
```

Nine of eleven negative and none positive, `nnz(L)` identical at every size. It is faster than
`AMD2` at every cubic size and faster than `AMD1` too, which no layer carrying the extras has been
before. Three days earlier this branch was at 3.0x on cubes.

**CORRECTED THE SAME EVENING, and the correction is the more useful entry.** This section first
read that `AMD3` reaches 1.01x, 1.02x and 1.07x the vendored routine at 12, 16 and 20 a side,
"which is parity". That was ONE RUN of a quotient whose denominator is the noisiest column in the
table. Over eight runs at 16 cubed, `AMD` reads 0.74 to 0.86 ms, a 16 percent spread, where `AMD3`
reads 0.83 to 0.89, a 7 percent one. **The vendored routine varies more than we do, so a
ratio-per-row is mostly a measurement of it.**

What reproduces, and it is still the best result this branch has had:

```
                 AMD3 ms        AMD ms        AMD3 / AMD over eight runs
12^3           0.26 - 0.28   0.20 - 0.27          0.98 - 1.33
16^3           0.83 - 0.89   0.74 - 0.86          0.97 - 1.18
26^3           4.65 - 4.89   3.98 - 4.13          1.13 - 1.19
32^3          10.29 - 10.75  8.43 - 8.83          1.18 - 1.27
```

So the honest statement is that `AMD3` is **at or near the vendored routine to 16 a side and rises
to about 1.2x by 32**. 16 cubed and 26 cubed reproduce to a percent; 12, 20 and 32 wobble by five.
Quote absolute milliseconds with the vendored range beside them rather than a ratio per row.

---

**Three things were not the algorithm, and together they were most of the day.**

**One: the arrays, which were worth more in 2D than the fold was.** The first version carried the
two values crossing from the prune to the bound in two fresh vectors of size n. It measured 3 to 9
percent faster on cubes and **12 percent slower in 2D from 200 a side up**. Both fit in arrays the
driver already had and that are dead at that point: `partial[u]` is not written until the end of
the bound pass, and `hashNext[u]` holds nothing until the vertex is filed, which happens in that
same pass after the key has been read. The key is reduced modulo the bucket count as it
accumulates, which is what makes it fit an `int32` rather than needing a `size_t` array. With the
arrays gone the 2D penalty went with them. **The visit count predicted the fold and said nothing
about the streams**, and the streams were the larger term on one family.

**Two: the scheduler.** The same code, on the same machine, read as a scattered null before a
`pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0)` went into the benchmark and as a
clean 3 to 9 percent after it. A command-line process on Apple Silicon runs at `QOS_CLASS_DEFAULT`,
which permits the scheduler to park the thread on an efficiency core, and that placement is STICKY
over a whole run rather than jittering per iteration. Every row here is a minimum over fifteen to
thirty repeats, which filters per-sample noise completely and a whole run on the wrong core not at
all. That is the shape of the 4 percent disagreement between identical binaries this folder has
been recording since 2026-08-08.

**Three: the near-miss, which a 2D-only check would have shipped.** A version accumulating the
WHOLE key in the prune, including the incidence half, failed 14 of 32 identity cases. Aggressive
absorption runs between the prune and the bound and compacts `I[u]` in place, so the list the key
must sum over does not exist yet at prune time. **Every square grid passed and every cubic and
random graph failed**, absorption firing far more there. Only the adjacency half can move, and
`I[u]` is still walked twice.

**AND THE MEASUREMENTS ABOVE CARRY A HARNESS BIAS, found after the fact.** `AMD3B` was timed
through `orderTimeFn`, which times the bare ordering call, and `AMD3` through `orderTime`, which
times `OrderEngine::compute` and so also builds a Permutation: two assigns of size n and a loop of
size n. Measured by timing identical code through both paths, that is **0 to 2.4 percent** in the
free function's favor. So every "AMD3B is X percent faster" figure this day, the 7 to 13 percent
above and the key fusion's 4.4 to 7.0 percent in iteration 25, is that much too generous.

Re-measured with both through the free-function path, same run, the fold is worth **10 to 16
percent on cubes**: 16.1 at 12 a side, 16.2 at 16, 9.8 at 20, 10.3 at 26 and 12.8 at 32. Larger
than what was first recorded, not smaller, but arrived at properly. **A comparison between two
columns of a benchmark has to go down the same path**, which is obvious once stated and was not
checked when `orderTimeFn` was written: its probe-and-repeat protocol was copied carefully and
nobody asked whether the thing being repeated was the same thing.

**Verified:** `make amdorder`, `make mmdorder` and `make aligned` 38 of 38 each, `Amd3B == Amd3` on
32 graphs including twelve random ones before the port, and 283/283 with 8 examples after it.

**The vehicle is gone again**, with its build entries and its two benchmark columns. What stays in
the shared class is the `TaggedScan` overload, which `Amd3` now uses.

---

**What is left.** Deleting `amd4`, once nothing further is wanted from it, lifting its postorder
block first if a permutation-level check is ever wanted.

---

## Iteration 27: seven folds, and the one that was not a constant, 2026-08-16

The ladder in this file ends with `Amd3` aligned to `AMD_2` and about 1.8x its speed on a 400
square, rising with n. That rise is now gone. The work was done in `src/Amd3B.cpp`, a private copy
carried for the purpose and since folded into `Amd3`, and the full account is in
`notes/DESIGN_DECISIONS.md` (2026-08-16). What belongs here is what the ladder's own method
contributed and where it misled.

**The method held.** Every fold was landed against a control, `AMD3f`, which is `Amd3` reached down
the same harness path as `Amd3B` so that the comparison is the change and not the seam. That
control also retired a figure this tree had been quoting: `Amd3.h` recorded the harness seam at up
to 2.4 percent, and measured directly it is ZERO across the whole ladder in both families. It was a
real reading on a different driver on a different day.

**Where reasoning failed and counting did not.** Three hypotheses were argued for in writing and
each was killed: that the clique arena's elimination-order placement was the growth, that the hash
was, and that single-level cache locality was. The differential in the README settled it by showing
both codes doing the SAME visits per pivot at every size in both families, which left cost per
visit and nothing else.

**Where the counting failed too, which is the more useful half.** The first differential counted
the passes already under suspicion. It did not count the pass that builds `C[pivot]` at all, nor
`absorb`, nor `clear_flag`, and it counted the hash outer loop on the two sides with different
denominators so that a 0.95 ratio read as agreement by coincidence. A flat table is only evidence
about the counters in it.

**And the ladder's own shape was confirmed once more.** Two folds measured nothing on their own,
the dead-clique test and the tagged `W` in `Amd2`, and the tagged `W` is what then made the fused
scan viable in both B layers. `AMD1B` went from 7 to 9 percent slower than `AMD1` in 2D to even or
better, and `AMD2B` from 3 to 9 percent slower to even, with the cubic advantage roughly doubling
in both. The fused schedule was never the problem; `ApproximateScan` crossing three arrays from the
prune was.

---

## Iteration 28: the restore leaves `AMD_2`'s schedule, 2026-08-24

**The change.** Both branches now restore the negated weights AT THE END OF THE PRUNE. Every earlier
entry here describes the other arrangement, and iteration 5 in particular turns on it, so this is a
departure from the reference and not a correction of our reading of it.

**Where the sign lived before, and it was two places for one fact.**

```
AMD_2    negate in the construct loop, restore under RESTORE DEGREE LISTS at line 2083,
         `nvi = -Nv [i]` then `Nv [i] = nvi`, folded into the pass computing the final degree
AmdFlat  the same, ported: restoreWeight rode the driver's fourth pass and returned the
         magnitude that pass needed one line later
MmdFlat  restore rode massEliminate's walk of C[pivot], the pass that finds what merged
```

**genmmd never negates anything.** Not one negative assignment in the file: it marks membership
with `marker[]` and a monotone tag and reads liveness off `qsize[nb] != 0`. So the sign trick is
`AMD_2`'s alone, and the mmd branch did not inherit it, it was GIVEN it on 2026-08-17 to fold away
an `inClique` stamp array. That is why the restore sat in two unrelated places: one reference's
technique retrofitted onto a branch with no opinion about it, each half then solving the placement
locally against whatever walk it already had.

**What made the move legal is that the prune is the LAST READER on both branches.** Its
`mWeight[v] <= 0` is what consumes the mark. Downstream, absorption never touches a weight, mass
elimination's merge test is structural on `adjacencySize == 0 && incidenceSize == 1`, and every
other reader takes a magnitude through `weight()`.

**What it removed.** `restoreWeight` and `restorePivotWeight` from both class surfaces and their
four call sites; and `mLateMassElimination`'s SECOND job, that flag having decided both "do not
merge here" and "do not restore here" while being named for the first alone.

**MEASURED AT NOTHING, and the direction is worth stating because it is the unfavourable one.** The
restore was a free rider on a walk that was happening anyway and is now a pass of its own over
C[pivot]. That should cost something. Over five grid sizes in both families, every driver of ours
lands inside the run-to-run spread of the two vendored controls, which nothing touched between the
runs: ours 0.963 to 1.034, `MmdCorrected` 0.979 to 1.035, `AmdVendored` 0.981 to 1.029. Fill
identical everywhere, all 365 digests identical.

**AND THE LADDER'S SHAPE HOLDS AGAIN, now from the other side.** Every entry above that moved a
clock moved it by REMOVING AN ARRAY: the eleven-to-five fold, the tagged `W`, the `inClique` stamp
this very sign trick bought. Every entry that only rescheduled work measured zero, and the list is
long enough now to be a rule rather than an observation: the dead-clique test, the tagged `W` on
its own, the fused scans, `eliminateAmd`'s fused walk earlier the same day, and this restore, which
is a FUSION UNDONE and measured zero in that direction too. Fusing and un-fusing are both free
here; folding an array is not.

The reason is not mysterious and it is worth writing once. A fold removes a cache line from every
touch of a vertex, forever. A schedule change moves the same loads and stores between passes over
data that is already warm, so it changes the instruction count by a few per cent and the miss count
by nothing, and a few per cent is under this machine's floor.
