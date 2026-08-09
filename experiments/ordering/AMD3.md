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

We had `mark` plus `outside` plus deadness-by-removal where it has `W`; `mMark` plus `mEliminated`
plus `mWeight` where it has `Nv`. **Two of the three were portable and are taken. The third is
not**, and iteration 17's read establishes why: `Nv` is negated in the construct loop and restored
in the very last pass, so it is negative across the entire body of an elimination and four separate
readers are written to expect that. It works because `AMD_2` is one function. Ours is a shared
class with six drivers, three of them MMD with different invariants: and the `mmd2` counterexample
in iteration 17 is what that costs when a convention crosses the seam.

**So the remaining gap does NOT have a name, and an earlier draft of this line said it did.** It is
not algorithmic: the counts are equal: same eliminations, same reachable-set elements, same prune
elements, same pairs tested, same fill. It is also not the one-fact-per-array-against-three pattern,
which was the leading explanation for most of a day and was measured false: we touch 1.09x as many
arrays per element as `AMD_2` and take 2.32x as long. `docs/DESIGN_DECISIONS.md` (2026-08-08) has
the table and the reason the wrong conclusion was comfortable. What is left is per-TOUCH cost, which
is a locality hypothesis and is untested.

### If someone picks this up again

The profile is diffuse now. `orderAmd3` is 48 percent of the run with no line above 378 ms, and
everything with a name has been taken. Getting from 2.3x to 1.5x is a third of the total and there
is no item of that size left in it. It would have to come from a driver owning its own storage
instead of sharing `QuotientGraph`: which is a real option, and is the trade `AMD_2` made, and
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
convention rather than a correction to the ladder, and `docs/DESIGN_DECISIONS.md` had already
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
ledger, not yet a verdict about the layer below. It was parked in `docs/TODO.md` and taken up in
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
`QuotientGraph` already holds both lists in one run behind `mSourcePtr`, which is the layout that
motivated the trick, so AMD's three moves transcribe almost literally.

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

**Decided: two flags and one method on `QuotientGraph`, all inert for the other five drivers**, in
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
faster at small. Whether that eats the saving is the next `make scale-amd` and nothing here
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

**It is in `QuotientGraph`, so it reaches every driver**, and the MMD branch forms cliques the same
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

**What is left.** Deleting `amd4`, once nothing further is wanted from it, lifting its postorder
block first if a permutation-level check is ever wanted.
