# Ordering Experiment

> **SUPERSEDED IN PART, 2026-08-15.** Every ordering TIME and every ratio against a vendored
> routine in this document predates the encoding work of that date and understates our orderings
> by roughly 20 to 30 percent. **Fill figures are unaffected**, every permutation and every nnz(L)
> being identical. `docs/DESIGN_DECISIONS.md` (2026-08-15) has the account.

`OrderEngine` calls a vendored AMD and a vendored MMD, and until now nothing in the port had
looked inside either. Both are hard to read: MMD is translated Fortran with `goto` labels and
1-based arithmetic, AMD is 1800 lines in one function over eleven parallel `int` arrays. This
experiment rebuilds the family one mechanism at a time, so that what each code *does* can be
separated from how it is *written*. Reference / teaching only, **not** part of the main Oblio
build.

Each layer exists twice, as `<name>.py` and `<name>.cpp`, printing the same trace. The Python is
the one to read and the C++ is the one to trust; `make test` checks they have not drifted.

## The layers

Each adds exactly one mechanism to the one before. The right column cites
`archive/sparse_factorization.md`, where the prose lives.

| layer | adds | sections |
|---|---|---|
| `md1` | naive minimum degree, materializing the fill | 5.1, 5.2 |
| `md2` | the quotient graph: cliques instead of fill | 5.3, 5.4 |
| `md3` | supervariables and mass elimination | 5.5, 5.6 |
| `md4` | maintained degrees, degree updated only where they changed | 5.7, 5.8 |
| `md5` | degree buckets, so the minimum is walked to, not scanned | 5.9, 5.10 |
| `mmd1` | multiple elimination: a batch of pivots per degree update pass | 5.11, 5.12 |
| `mmd2` | the rest of genmmd, one pass at a time | 5.11, 5.12 |
| `mmd3` | mmd2 with genmmd's list order, to match its permutation | 5.11, 5.12 |
| `amd1` | approximate degree: a bound instead of a set union | 5.13, 5.14 |
| `amd2` | the rest of amd_1 and amd_2, one pass at a time | 5.13, 5.14 |
| `amd3` | amd2 with AMD_2's list order, to match its permutation | 5.13, 5.14 |

**`amd4` is PARKED, and it used to be called amd3.** It was renamed on 2026-08-08 so that `amd3`
could be the layer aligned to the vendored AMD, which is what `mmd3` already is on the other
branch: the digit now means the same thing on both, 3 is the aligned layer. `amd4` forks sideways
from `amd2` and does NOT contain `amd3`, its header says so, and it is temporary, kept to read
from and scheduled for deletion. It was built to carry amd_1's input path, the dense-row removal
and the postorder, and it answered what it was built to answer. It is not a comparison point and
no question is open against it.

Parked does not mean abandoned. It stays in `make test`, in `GRID_LAYERS` and in the twin check,
because its Python and C++ twins are the only oracle it has, production having no counterpart to
check it against, and a layer left to drift silently costs more later than one kept green. Its
`matrix1` example stays with it for the same reason.

What parked does mean is that it is not evidence. A result that holds in `amd4` and nowhere else
settles nothing about the four, and a difference between `amd2` and `amd4` is not by itself a
defect worth chasing.

## One word, before anything else

An eliminated pivot's fill-in structure is a **clique** here, never an element. genmmd and AMD both
call it an *element*, and so does most of the literature in that lineage, so the word appears
throughout this file wherever the prose is walking through their code: their `ehead` chain, their
`Elen` array, `mmdupd`'s per-clique loop. Element is their name for our clique, and it stands only
inside those quotations, so that an explanation still matches the source it explains. Everywhere
else, including our own comments and identifiers, it is clique. The rule is in
`docs/WRITING_RULES.md`.

Our own identifiers all follow it as of 2026-08-20, when `element_tag`, `elementTag`,
`element_members`, `elementMembers`, `is_element` and `isElement` were renamed in the mmd2, mmd3
and amd4 twins, alongside the same sweep through the production drivers.

## What the layers show

**`md1` through `md4` return the same ordering.** That is the point of the first four sections:
the heuristic was fixed in `md1` and everything after is implementation, so those layers can be
verified by demanding an identical permutation. From `md5` on the permutation moves, and the
reason is worth separating from the reason the fill moves.

Four things can happen when a layer is added, and all four occur here:

```
                              order        fill        what the change is
md1 -> md2                    same         same        a change of representation
md2 -> md3  (mass elim.)      DIFFERENT    same*       a reordering, free iteration by iteration
md3 -> md4                    same         same        a change of implementation
md4 -> md5  (buckets)         DIFFERENT    same-ish    a change of TIE-BREAK
md5 -> mmd1 (multiple elim.)  different     DIFFERENT  a wager
mmd1 -> mmd2 (genmmd's rest)  different     different  fidelity to the vendored routine
```

`md5`'s change is the odd one, and it is a consequence of the data structure rather than of the
heuristic. A bucket is a linked list pushed and popped at the head, which is the only O(1)
structure for the job, so the winner among equal degrees is whatever was filed last rather than
the lowest index. On the seven examples the order differs every time and `nnz(L)` is identical
every time; on 200 random graphs the order differs on 198, `nnz(L)` agrees on 181, and where it
differs `md5` is better on 12 and worse on 7. Different, not worse: the pivots are still exact
minima and only the choice among equals moves.

So `mmd1` is not the first layer to change the permutation, and neither is `md5`. Mass
elimination already does it at `md3`, on nine of twelve test graphs, with identical `nnz(L)` on
all twelve. What `mmd1` is first to change *systematically* is the **fill**, because it is the
first to choose a pivot on stale information.

`mmd2` is the odd one out: its additions are not there to improve anything but to match what the
vendored routine does, so each pass is judged by fidelity rather than by fill. The first of them,
the prepass, costs a little fill on 10 of 207 graphs and gains none, which is the expected shape.

The asterisk on `md2 -> md3`: the individual merge is provably free, but the equality of the
totals across the whole run is measured rather than proved. See the open question below.

**`mmd1` and `amd1` give up different things.** `mmd1`'s pivots are always true minimum-degree
vertices; only the tie among equals is broken differently, because a batch evicts what it has
touched. `amd1`'s pivots may simply not be minimal, because an overcounted bound can hide the true
minimum. MMD perturbs; AMD can be wrong. Both cost well under a percent of fill, in either
direction, and both are noticeably faster.

## The test graphs

**One definition, in `graphs.h`.** The seven examples were written out twice, in `vendored.cpp`
and in `production.cpp`, and a third copy was about to go into `amdorder.cpp`. Three copies of a
literal is not a style complaint here: the drivers answer different questions about the SAME
graphs and `make test` compares their outputs line for line, so a graph that drifted in one copy
would make two drivers disagree for a reason that is not the code, and that would read exactly
like a defect in an ordering. The header also holds the grid builders, so a shape added for one
driver is available to the others, which is what the acceptance test's widening needed. Each
driver keeps its own conversion to CSC, deliberately: ours takes a full-symmetric pattern with the
diagonal present, and the vendored routines take the off-diagonal pattern alone.

`graph1` is a 4-cycle, the smallest graph that fills at all: eliminating any vertex forces its
two neighbors together, for one fill edge. It is also where md3 merges everything that is left
in a single iteration, 1 taking 2 and 3, which makes it the simplest case for reading the cost of
mass elimination, at the price of the run ending there. `graph2` has six vertices, `graph3` has
twelve and is the first whose ordering is not the identity. `graph4` has eight vertices and
fourteen edges and exists for one reason: **it is the smallest graph we found on which AMD's
bound is ever loose.** The bound overcounts only when a vertex belongs to two cliques that
overlap outside the new one, which needs enough eliminations to have made several cliques and
enough fill for them to intersect.
Checked exhaustively, no connected graph on five or six vertices is ever loose anywhere in its
run, and none in thirty thousand samples on seven. Without `graph4` the `amd1` trace would display
the whole algorithm and never once show it approximating.

`graph5` has five vertices and four edges, two paths joined at 4, and is present in `md1`, `md2`
and `md3` only. It is the smallest graph on which **md3's merge test declines a genuine
supervariable**: at the iteration whose pivot is 0, vertex 4 has nothing explicit left but belongs
to `c1` as well as to the new clique, so `I[4] == {pivot}` fails even though everything 4 reaches
lies
inside that clique. It orders as 2 1 3 0 4 with no merge and no fill, where the exact test would
give 2 1 3 (0 4). It also separates amd2's two extra mechanisms: with aggressive absorption on,
`amd2` takes four iterations and reports `merged = 4, absorbed = c1`; with it off, five iterations
and no merge, exactly like `md3`. The hashing plays no part either way. See the section on mass
elimination.

`graph6` has six vertices and eight edges and is also in `md1`, `md2` and `md3` only. One small
graph carries three things at once. **Its supervariable {0, 4} is a supernode but not a
fundamental one**, since the forest is 2 -> 1 -> 4 and 3 -> 0 -> 4 and 4 already has 1 as a child
when 0 merges into it. It orders as 1 5 (0 4) 2 3, where the exact test would give
1 5 (0 2 3 4). The merge lands at iteration 2 of 5, so the run continues afterwards and the
selection degree, 3 over {2, 3, 4}, differs from the external degree, 2 over {2, 3}, by exactly
the size of what merged. And `super_members` ends with a hole in the middle, slot 4 empty between
two used ones, while no pivot equals its own iteration number.

`graph7` has five vertices and six edges and is the **pairwise case**, also in `md1`, `md2` and
`md3` only. At the iteration whose pivot is 0 and whose clique is {2, 4}, vertices 2 and 4 are
indistinguishable from each other, both reaching the same closed neighborhood, yet neither is
absorbable into the pivot since each still reaches 3 from outside the clique. It orders as
1 0 (2 3 4) under both md3's test and the exact one, since no test framed against the pivot can
see such a pair. Catching them needs a comparison between candidates, which is amd2's hashing.

The three small graphs added for the mass elimination story line up as follows.

| example | what it shows | md3 | exact test |
|---|---|---|---|
| `graph5` | a supervariable the cheap test declines | 2 1 3 0 4 | 2 1 3 (0 4) |
| `graph7` | a pair no pivot-relative test can see | 1 0 (2 3 4) | 1 0 (2 3 4) |
| `graph6` | a supervariable that is not a fundamental supernode | 1 5 (0 4) 2 3 | 1 5 (0 2 3 4) |

## Build

```
make        build every C++ prototype
make test   build them, run each, and check the C++ agrees with its Python twin
make clean
```

`make test` compares traces after stripping spaces, brackets, braces, quotes and commas, because
Python renders a dict as `{0: [1,3]}` and the C++ renders `0: {1,3}`. Every number, name and label
still has to match. The alternative was editing fourteen files to unify the spellings, which
seemed the worse trade.

Each prototype also takes an example number, `python3 md3.py 3` or `./md3_cpp 3`, and runs every
example when given none. Every layer above md5 takes a grid instead, `python3 amd2.py grid 22` or
`./amd2_cpp grid 22`, which prints the closing counters and no trace; `make test` diffs
the twins on grids as well as on the examples, since the examples are too small to fire some of
the mechanisms at all. See the grid-mode section below.

Once a layer has been pulled into the main tree it gains a second check, against the production
driver extracted from it:

```
make production       builds production_cpp, linking ../../src
./production_cpp mmd1        MMD1 on all seven graphs
./production_cpp amd2        AMD2 on all seven graphs
./production_cpp mmd1 3      just the third
./production_cpp amd2 grid 20   one 20x20 grid, the same one the prototypes build
```

`make test` runs it for every layer named in the Makefile's `PORTED` list, `mmd1`, `mmd2`, `amd1`
and `amd2`, and requires the order lines to agree with the prototype's, on the seven examples and
then on the grid sides in `GRID_SIDES`. Both halves are permutation checks and the second exists
because the first cannot see a mechanism that needs real structure to fire: two defects in amd2
left all seven examples byte for byte identical while the ordering was wrong on any grid of 10 a
side or more. `PORTED_FILL` names layers to be checked by nnz(L) instead, and is empty; it was how
amd2 was checked until the postorder moved to amd4, since production skips the postorder and so
the permutations legitimately differed while the fill did not. That weaker check runs Oblio's own
symbolic factorization, which is why the target links the whole library.

The target links `../../src` directly rather than copying it, which is the opposite of what the
vendored target does and deliberately so: a copy is right for code that is not ours to edit and
wrong for code being actively changed at both ends, since noticing when the two come apart is the
whole point. The harness feeds each graph as a full-symmetric CSC with the diagonal present, which
is what a `SparseMatrix` holds, so the production path under check includes `QuotientGraph`
dropping the diagonal rather than only the driver.

The vendored routines have their own target, since they are not layers and have no Python twin:

```
make vendored         builds vendored_cpp
./vendored_cpp        both routines on all seven graphs
./vendored_cpp 3      just the third
```

That target compiles with warnings off, because the two files are not ours to clean up.

The two routines are compiled straight from `private/Mmd.cpp` and `private/Amd.cpp` at the repo
root, which is where they live now: they are not ours, are not published, and this directory no
longer keeps copies of them. `vendored.cpp` only feeds them the same seven graphs and prints their
permutations in our format. When `private/` is absent the target reports that it skipped rather
than failing, so every other layer still builds.

The lowercase names are deliberate, and so is the subfolder. Oblio capitalizes source files, but
macOS formats APFS case-insensitive by default, so `Amd.cpp` and `amd.cpp` are one path: dropping
the vendored copy beside the prototype overwrote the amd prototype silently, and `git status` said
nothing because `core.ignorecase` is on. The subfolder keeps the two kinds of file apart, and the
prefix means no vendored name can ever collide with a prototype.

## What is not implemented

The `mmd1` and `amd1` prototypes are deliberately subsets of the vendored routines, and the plan
for closing the gap is the section that follows. Each file header carries its own list, and
sections 5.11 and 5.13 carry the same lists in prose. In brief:

- **`mmd1`** lacks the prepass that numbers degree-0 and degree-1 vertices before the main loop, and
  `mmdupd`'s `q2h` merging of vertices indistinguishable *from each other* rather than from the
  pivot. It also files degrees at a different offset and never uses bucket 0.
- **`amd1`** is the bound and nothing else. Everything beyond it is `amd2`'s, and `amd2` now
  carries all seven passes: aggressive absorption, hash supervariable detection, the two-pass
  degree update, dense-row handling, `amd_aat` and `amd_preprocess`, the postorder, and
  `amd_valid` with the `Control`/`Info` interface.
- **Neither carries the workspace pool**, `Iw`, `Pe`, `pfree` and the compaction counted by
  `ncmpa`, and that one is a decision rather than a gap. The section on it below sets out why a
  hand-rolled pool is the wrong thing to port.

## The plan for mmd and amd

md1 through md5 were a teaching ladder: each isolated one idea and changed exactly one property,
representation, then order, then implementation twice. What is left of the vendored routines is
not a sequence of ideas but the completion of two, so from here the iterations are bigger. Two
versions each.

**mmd1, the idea. Done.** Multiple elimination. A batch is an independent set in the current
elimination graph, enforced by evicting every reached vertex from its bucket with a stale degree,
so no later pivot in the batch can be a neighbor of an earlier one. `delta` widens the batch to
near-minima. Everything else is md5 unchanged: the quotient graph, mass elimination, the buckets,
the expansion.

`delta` lives here, with the full signed range: negative takes one pivot per iteration, which turns
the batching off and reproduces md5's ordering exactly, verified on 100 random graphs; 0 through
n - 1 widen the window; anything larger saturates. No weight array, because mass elimination
merges only into the pivot and the pivot dies in the same call, so no live vertex ever stands for
more than one original vertex, checked over 200 graphs and 1386 eliminations.

**mmd2, genmmd complete.** Six additions, all of them holes rather than ideas, listed against
the vendored routine that carries each:

1. **The prepass** (`genmmd`, the loop over `head[1]` before the main loop). It numbers every
   vertex in the degree-1 list, marks each `marker[mn] = maxint`, and never updates a neighbor's
degree. `mmdint` maps degree 0 to 1, so isolated and degree-1 vertices are numbered together. Two
   details travel with it: `head[1] = 0` afterwards, and the main loop then starts at `mdeg = 2`.

2. **The q2h path** (`mmdelm` stashes `fwd[rn] = nq + 1`; `mmdupd` splits on it). A reached vertex
   with exactly one explicit neighbor left besides the element goes on the `q2h` list, everything
   else on `qxh`, and the two are walked separately because the q2h case can be answered without
   a full union.

3. **The pairwise merge**, inside the q2h walk: `qsize[en] += qsize[nd]` with `fwd[nd] = -en`.
   This is what makes MMD's supervariables coarser than ours, and it folds a vertex into a LIVE
   one, so from here a candidate can stand for several original vertices.

4. **Outmatched marking**, `bwd[nd] = -maxint` in the same walk. It takes a vertex out of the
   degree lists without merging it; `mmdelm` puts it back with `bwd[rn] = 0` the next time it is
   reached.

5. **The filing convention**, `dg = dg - qsize[en] + 1` floored at 1 (`mmdupd`, n2100), so the
   least bucket is 1 and bucket 0 is never used. `dg` there is a weighted sum over the reached
   set, so it arrives with the weighted counting item 3 brings.

6. **The counters**: `ncsub += mdeg + qsize[mn] - 2`, the subscript statistic, and the early
   termination `if (num + qsize[mn] > neqns) goto n1000`, which stops the main loop as soon as
   the last vertex is accounted for rather than eliminating it.

Items 3 and 5 make a candidate stand for several original vertices, so every degree becomes a
weighted count. It is taken from `len(super_members[v])`, which is O(1), so no weight array is
introduced. Our `super_members` already does what `mmdnum` does, so the expansion needs nothing.

**amd1, the idea. Done.** The approximate degree bound, and nothing else. The same rule mmd1
follows: one layer, one idea. Everything md5 had is kept unchanged, and the only thing that moves
is how a degree is estimated after an elimination.

**amd2, amd_1 and amd_2 complete.** Seven additions, the first two of which are mechanisms that
ride along with the bound's own work and the rest of which are holes:

1. **Aggressive absorption.** |C[c] - C[p]| has just been computed for every clique the new one
   touched, and a zero means C[c] lies inside C[p], so that clique is dead. Ordinary absorption
   kills only what the pivot touched; this kills what any reached vertex touched, at no extra
   cost.

2. **Hash supervariable detection.** Mass elimination merges a vertex into the pivot. Two vertices
   can be indistinguishable from each other without either being absorbable into the pivot. AMD
   hashes (A[u], I[u]), compares within a hash bucket, and merges on an exact match; the hash is
   a filter, never the decision.

3. **The two-pass degree update**, scan 1 obtaining |C[c] - C[p]| by subtraction from a maintained
   clique degree rather than by walking the clique's members.

4. **Dense row and column detection** by the alpha ratio, with those vertices held out of the
   ordering and placed last, at the cost of nnz(L) becoming an upper bound once any fire.

5. **`amd_aat`**, forming the pattern of A + A' with the diagonal dropped, and `amd_preprocess`.

6. **`amd_postorder`**, so the output is a postorder of the assembly tree rather than raw
   elimination order.

7. **`amd_valid`** as an input check, and the `Control`/`Info` interface.

Item 2 makes a candidate stand for several original vertices, exactly as mmd2's items 3 and 5 do,
so every degree there becomes a weighted count taken from `len(super_members[v])`.

### What is deliberately excluded

One piece of the vendored codes is a consequence of packing state into reusable integer arrays
rather than a feature of the ordering, and it is not modeled:

- AMD's workspace with `iwlen`, `pfree` and the `ncmpa` garbage collection. That is a flat pool
  being compacted when it fills, and our member lists grow on their own. It is what the ordering
  engine will need once the lists live in one pool.

**MMD's `maxint` overflow reset used to be the second entry here, and it is now implemented.** The
claim it rested on, that the tag cannot wrap at the sizes these prototypes run, is true and is not
the point: a mark array is a set only while its tag is unique, and a wrapped tag makes a stale
stamp read as a match, which is wrong with no symptom. See the section on the tag guard below.

### Tie-breaks, and what the acceptance test is

The vendored routines break ties differently from every layer here, and not by a rule: a degree
list is a singly linked chain prepended at `head[dg]`, in `mmdint` and again in `mmdupd`, and
the pop takes the front. So the bucket is a stack and the winner is whatever was pushed last.
After construction, which walks `nd` upward, that is the highest-numbered vertex of its degree,
which is why the vendored MMD starts graph1 at vertex 3 where we start at 0. There is no
quality claim behind it: prepending is the cheap end of a linked list, and a tail pointer would
have cost another array of size n.

We now do the same, because the alternative costs a log factor on the dominant operation and the
goal is to match the vendored asymptotics. Buckets are linked lists in both twins from md5 on, so
our tie-break is theirs in kind: whatever was filed last. It is not theirs in detail, since the
insertion history differs, and reproducing that exactly would mean reproducing the order in which
`mmdupd` walks `q2h` and then `qxh`, which is coupled to the features rather than separable from
them. So the permutations still differ, and exact equality stays a possible later mode rather
than the acceptance test.

**What the tie-break is worth, measured 2026-07-31, and it is more than expected.** On small
graphs ties are few and the fill agrees with the vendored routine on 7 of 7 examples and 59 of 60
random graphs. On a grid nearly every live vertex has the same degree, so a tie decides almost
every pick, and the effect compounds. Four filing orders of the same algorithm, against the
vendored MMD, counting nnz(L) through Oblio's own symbolic factorization:

```
                              32x32     64x64   100x100   140x140
vendored MMD                  11822     63219    186835    412921
ours, as it stands            11972     71709    223806    492921
initial buckets descending    12093     67109    194505    443997
degree update order reversed  12570     73184    213784    504177
both                          12074     71487    218989    513689
```

Two things follow, and the second matters more.

**The gap to the vendored routine is a tie-break artifact, not a missing mechanism.** At 140x140
the four orders span 443997 to 513689, a 16 percent spread, which is the size of the whole gap. So
nothing is wrong with what we compute: MMD1 and MMD2 are choosing differently among equally
minimal vertices, and on a regular mesh that choice is nearly the entire ordering.

**And no order among these is right.** Filing the initial buckets in descending vertex order beats
ours at three sizes out of four and loses at the fourth, which is what an arbitrary choice looks
like when it is measured rather than reasoned about. So the code is unchanged: adopting a rule that
wins on three grids would be adopting an unproven heuristic, which is what this experiment
otherwise declines to do. Reproducing the vendored fill exactly would mean reproducing its
insertion history, which the paragraph above explains is coupled to `mmdupd`'s q2h and qxh walks
rather than separable from them.

The test is `vendored.cpp`, which links `private/Mmd.cpp` and `private/Amd.cpp` from the repo root,
and runs both routines on the same seven graphs, printing permutations in our format. mmd2 and amd2
are accepted when every feature above is present and exercised, nnz(L) matches the vendored
routine on the seven graphs and on random ones, and every remaining order difference is traceable
to a tie.

## Aligning a layer against a vendored routine: the method

Both alignments used this and it is not specific to either. It was written down in a handover note
between two sessions, `NEXT.md`, which was deleted once the amd work it handed over was finished;
what follows is the part of it that was method rather than history. The two narratives,
`MMD3.md` and `AMD3.md`, are what happened; this is how.

**The loop.**

```
1. run ours against the vendored routine on the SMALLEST case that still diverges
2. find the FIRST pivot where they differ
3. read both traces at that pivot, instrumenting the VENDORED code if needed
4. root-cause it: which line of the vendored routine does ours fail to reproduce
5. change exactly that one thing
6. re-run everything; record the entry in the ledger
7. repeat
```

Deliberately slow, one entry at a time. **Do not batch fixes.** The mmd work reached 4 of 7 examples
with three changes and could not tell which of the three mattered until a fourth was found and took
it to 7 of 7 in one step. The amd work hit the same shape: entry 1 alone moved the score not at all,
and only entry 2 beneath it showed why it had been right.

### Choose the smallest case, always

For both branches the sequence was the seven examples first, `n` from 4 to 12; then small square
grids once every example matched; then larger grids only to confirm, never to diagnose.

A 25-vertex grid diverging at pivot 19 is readable in full. A 1024-vertex grid diverging at pivot
700 is not, and it is the same defect. Entry 1 of the amd ledger was found by hand on a 4-cycle.

**But the examples are not the test.** With four entries in, the amd work was 7 of 7 on the examples
while grids of side 8 and 10 still diverged. Twelve vertices cannot exercise a mechanism that needs
real structure. Finish on grids.

### Compare the RIGHT object

This mattered more than anything else and it is easy to get wrong. Three comparisons, increasing in
strictness:

```
fill                 nnz(L). Coarsest. Equal fill does NOT mean equal ordering.
pivot sequence       the pivots chosen, in order. This is the ALGORITHM.
permutation          the expanded order, including supervariable members.
```

After entry 5 the mmd work looked stuck: fill exact at every size, permutation still diverging at
pivot 700 of 1024. Comparing PIVOT SEQUENCES showed they were identical, 788 of 788. The algorithm
was already aligned and only the numbering of indistinguishable members differed. Comparing
permutations would have spent the next hour hunting a mechanism that was not missing.

**So: if fill matches but the permutation does not, compare pivot sequences before doing anything
else.**

And the object can be upgraded once the work is nearly done. The amd alignment was settled on pivot
sequences because `AMD_2` postorders and we do not, so permutations could not match. Late on, the
raw elimination order turned out to be reconstructible upstream of the postorder, which made the
test a full permutation including member order. **The acceptance test is worth revisiting: it gets
chosen when an obstacle looks immovable, and the obstacle may not be.**

### Instrument the ORACLE, not just our side

The single most productive technique. Entries 5 and 6 of the mmd ledger and entries 2, 3 and 5 of
the amd ledger were all found by adding `fprintf` to a scratch copy of the vendored source.
Reasoning from our own trace alone produced a plausible and wrong diagnosis more than once.

**Never edit `private/`.** Copy the vendored file to `/tmp`, rename it as a header, and `#include`
it from a small driver. `genmmd` and `AMD_2` are `static` in an anonymous namespace and cannot be
linked from outside; including the source is the way in.

```
cp private/Amd.cpp /tmp/probe/Amd_p.h        # scratch copy, never the original
# add #include <cstdio> at the top, counters as file-scope longs
# then a driver:  #include "Amd_p.h"   plus a main() that builds a grid
g++ -std=c++17 -O3 -w -I/tmp/probe /tmp/probe/probe.cpp -o /tmp/probe/probe
```

Two cautions learned the hard way:

- **Brace counting.** Turning `else if(cond)stmt;}}}}` into `else if(cond){COUNT++;stmt;}}}}` eats a
  closing brace. The vendored code is dense and brace-heavy; count them.
- **`goto` crosses initialization.** `genmmd` is transliterated Fortran and full of `goto`, and
  adding `COUNT++;` on its own line before a labeled statement can make a declaration newly
  crossed and break compilation. Use the comma operator or put the counter after the label.
  (`AMD_2` has no `goto` at all, so this one is mmd's problem only.)

**And read `Info[]` before writing a single `fprintf`.** `AMD_2` fills `AMD_LNZ`, `AMD_NDIV`,
`AMD_NMULTSUBS_LDL` and `AMD_DMAX` natively, which are free comparison points. The amd work used
`AMD_LNZ` from the start and discovered the others at the end, having instrumented by hand what was
already being reported.

### The two alignment checks, and why one needs a hook and the other does not

The pair is `make amdorder` and `make mmdorder`, with `make aligned` running both. They ask the
same question of the two branches, does our layer still compute what the reference routine
computes, so they are **named for the branch and not for the mechanism**. The mechanisms are not
alike, and the difference is the reference routines' rather than ours:

**THE MMD REFERENCE IS `MmdCorrected` SINCE 2026-08-23**, not `MmdVendored`. genmmd files a vertex
under its degree in `mmdint` and under its degree PLUS ONE in `mmdupd`, so the minimum it selects
is not always the minimum; `private/MmdCorrected.cpp` is genmmd with that repaired and is what our
drivers match. The frozen copy stays and nothing compares against it. The amd reference is
unchanged.

```
genmmd emits the order DIRECTLY   mmd_order returns perm, the order it eliminates in, and there
                                  is no postorder anywhere in the routine. The vendored output
                                  vector IS the object to compare. No hook, no generated copy, no
                                  anchors to assert, no Control array.

AMD gets it through a HOOK        amd_order returns a vector AMD_postorder has already relabeled,
                                  so the order AMD_2 would emit at the end of its main loop has to
                                  be reconstructed upstream of the relabeling, by a generated
                                  copy of the vendored source with four insertion points in it.
```

The amd target was called `raworder` until 2026-08-09, for exactly the difference above. Naming one
target for how it works and the other for what it checks made a matched pair read as two different
kinds of thing, so the asymmetry moved into prose, which is where a fact about `AMD_2` belongs.

The rest of this section is the amd half, which is the one with machinery in it. The mmd half is
described under "mmd3, and the alignment ledger" and is a single file with no generator beside it.

### The comparison object for amd is the RAW ORDER, and it has to be built

`AMD_2` does not emit the elimination order. `amd_order` returns `Perm`, which `AMD_1` has already
relabeled by `AMD_postorder`, and that postorder is a heuristic tidy of an approximate assembly tree
that Oblio replaces with Liu's rule on the exact supernodal tree. So the vendored output vector is
NOT what to compare against: our permutation will differ from it by construction and always should.

**What to compare is the order `AMD_2` would emit if it stopped at the end of its main loop**, which
is the pivot sequence together with the member order inside each supervariable. That is the whole
algorithm, and it is what `amd3` and production `Amd3` reproduce exactly.

### Running it

```
make amdorder
```

That generates a hooked copy of the vendored source, builds the checker, and compares production
`Amd3` against it on **four shapes**: the seven examples, 2D grids from 4 a side to 140, 3D grids
from 2 to 24, and nine random patterns at n = 2000. `make clean` removes both the generated source
and the binary; `./amdorder_cpp 20 50` after building runs chosen 2D sizes instead, for bisecting a
failure.

**Four shapes and not one shape at many sizes, which is the whole point of the list.** Widening a
square grid from 4 to 140 exercises scale and never mechanism, so a defect that needs a structure
grids do not produce passes every size of it. That is not hypothetical here: the 2D-only version of
this check was green while production `Amd3` carried the stale clique degree of ledger entry 7, and
a 3D grid at 16 a side finds it. Each shape earns its place differently. The examples are small
enough to read by hand and are where most ledger entries were found. 2D grids give regular
structure at size. 3D grids fill faster, make larger cliques and mass-eliminate far more often.
The random patterns give irregular structure at size, and they are the only family here that
reaches a degree the dense threshold could act on.

**Two settings the driver must get right, and both were wrong when the wider shapes first went in.**

- **The dense threshold has to be derived from `n`.** Dense-row removal is the one mechanism amd3
  does not have, so if it fires the oracle has ordered a different problem and no comparison means
  anything. `Amd.cpp` turns it off through `dense = alpha * sqrt((double) n)`, then `MAX (16,
  dense)` and `MIN (n, dense)`, and its header says to pass "a number larger than sqrt (n)" without
  saying how much larger. It matters: `dense` is an `Int`, so an alpha that drives the product past
  `INT32_MAX` makes the conversion undefined. `1e30` lands on `INT_MIN` on x86-64, `MAX (16,
  INT_MIN)` is 16, and the threshold comes out at SIXTEEN rather than `n`, which is dense removal
  fully on at the strictest setting the code can express; on arm64 the same conversion saturates
  the other way and gives `n`. So the two platforms disagreed about what that line did. Neither
  grid family can show it, both being far below degree 16. `alpha = n` cannot overflow at any size
  we run and clamps to `n` throughout.
- **And it is checked rather than trusted.** `Info[AMD_NDENSE]` reports what was removed, so the
  claim costs one read, and the driver fails the case and says so. Without it a mis-set threshold
  arrives as a size mismatch to be diagnosed, which is how the paragraph above was found; with it
  the instrument names the fault instead of reporting that there is one. An instrument that
  silently declines to measure is worse than one that is absent.

**A third thing the widening turned up, and it is in the graph builders rather than the check.** A
column's row indices must be ascending. It is the CSC precondition `SparseMatrix` states in its own
header, and it is also a TIE-BREAK INPUT, since the order within a column decides the content order
of `C[pivot]` and so which of several equal-degree vertices a bucket hands over first. The 2D
builder is ascending by construction; a 3D builder written the natural way, `x-1, x+1, y-1, y+1,
z-1, z+1`, is not, and the far corner of a 4x4x4 grid comes out `62, 59, 47`. Ours consumed the
pattern as given and the vendored side sorted it, so the two parted company at the FIRST
elimination and the divergence looked like a defect in the ordering. `graphs.h` builds ascending
and says why at the site.

The generated copy is `amd_raw.cpp`, written by `tools/hook_amd.py` from whatever `private/Amd.cpp`
currently says. It is gitignored and removed by `clean`, exactly as the int64 copies the width study
uses are, and for the same reason: a checked-in copy of vendored code carrying our edits would be a
third thing, drifting from the original with nothing to notice. **A test comparing against a stale
oracle is worse than no test**, and this tree has three recorded instances of an instrument quietly
declining to do its job.

`tools/hook_amd.py` asserts every anchor it depends on. If the vendored source moves, generation
FAILS with a message naming which anchor went, rather than producing a copy with the hook in the
wrong place. Fix the anchor there; do not loosen it.

### What the hook does, and the three paths it has to cover

Two file-scope containers:

```c
PB_members[i]   what supervariable i currently stands for, seeded as {i}
PB_raw          the raw elimination order
```

**A vertex is numbered on THREE paths, not one**, and this is the part worth keeping. Building the
generator found it empirically: covering only the obvious path left the order short by three
entries on most grids and six on one, and an irregular deficit like that is the signature of an
unhandled path rather than an off-by-one.

```
EMPTY VARIABLES     numbered in the INITIALIZATION, before the main loop runs. A vertex with no
                    off-diagonal entry never forms an element, so it never reaches the finalize
                    marker.  `Elen [i] = FLIP (1) ; nel++ ;`

MASS ELIMINATION    scan 2 folds a variable whose degree has fallen to the pivot's straight into
                    `me`, with `Nv [i] = 0`, and it never goes through the hash merge.

THE HASH MERGE      supervariable detection, `Nv [i] += Nv [j]`. Note the DIRECTION: j is folded
                    into i, so it is i's list that grows.
```

**And the source list must be CLEARED after each merge.** Without that, a member reached through a
chain of merges is copied more than once, which produces sizes that match with contents that do not.

A fourth insertion point, at `FINALIZE THE NEW ELEMENT`, takes the pivot's members into `PB_raw`
once per elimination. Nothing in it comes from `AMD_postorder`, which runs later in `AMD_1` and only
relabels.

**This is the acceptance test to use from now on**, and it holds entry for entry, member order
included.

**Note what this does NOT claim.** Our AMD3 is not a drop-in match for `amd_order`'s output vector.
Same fill, same elimination order underneath, different labels, because we never compute their
postorder. Anyone comparing output vectors directly will see a difference and it is the intended
one. The fill agrees too, 206332 at 100 a side and 474995 at 140, though NOT by construction: a
postorder of the ELIMINATION tree cannot change fill, and AMD postorders its ASSEMBLY tree, which
its own header says need not be that tree because mass elimination under an approximate degree
merges vertices that were never adjacent. Measured on 2026-08-09, the two agree on every square grid
and on cubic grids from 7 a side up, and differ by one to three entries at 4^3, 5^3 and 6^3. This
paragraph claimed the invariance until then. What the acceptance test compares is the PERMUTATION,
so nothing here reaches it, and `benchmarks/ordering` carries an `AMDraw` column so the two fill
figures can be seen side by side.

**And the acceptance test was upgraded late.** It started as pivot sequences alone, because
`AMD_2` postorders and we do not, so permutations could not be compared. The raw order turned out to
be reconstructible upstream of the postorder, which made the test a full permutation including
member order. Worth revisiting an acceptance test when it was chosen around an obstacle: the
obstacle may not be one.

### Discipline that applied throughout

- **Port, don't rewrite.** Every entry was found by asking which vendored line we failed to
  reproduce, never by reasoning about what would be better.
- **Verify at three levels after every change**: twins agree with each other, prototype agrees with
  production, and the whole thing agrees with the oracle. `make test` does the first two.
- **Never touch `private/`.** Scratch copies in `/tmp`.
- **The ledger is APPEND ONLY.** A row is never edited once closed, so the sequence stays a record
  of what was actually wrong rather than a summary written afterwards.
- **One change at a time**, even when three look obvious.

### Two failures that are worth expecting

**A better score is not a correct port.** The mmd work found a change that improved fill and was
wrong; the question to ask is not "did this help" but "which line of the vendored routine does this
correspond to". Recorded in full in `MMD3.md`, iteration 3.

**The symptom does not identify the cause.** Both branches produced a divergence whose obvious
reading was wrong, and in both cases the actual cause was one layer beneath. `MMD3.md` iteration 5
and `AMD3.md` iteration 3 are the two instances.


## mmd3, and the alignment ledger

**mmd3 adds no mechanism. It exists to return genmmd's permutation, and it is an anchor rather
than a candidate ordering.** `experiments/ordering/MMD3.md` is the narrative counterpart of
this section, iteration by iteration. mmd2 has every mechanism genmmd has and still returned a different
permutation, which meant every comparison between them measured two things at once: a difference
of MECHANISM and a difference of ARBITRARY CHOICE. Minimum degree is a tie-break algorithm, so at
almost every iteration several vertices share the least degree and the winner is whichever the
data structure hands over first. Two codes can agree on every rule and part company on the first
tie, and from there they are ordering different graphs.

**It worked, and the result is stronger than a tidy permutation.** mmd3 now returns genmmd's
permutation EXACTLY, on all seven examples and on every square grid tested from 5 a side to 80,
`n = 6400`. Its fill is genmmd's to the digit at every size on the scale ladder, where mmd2 ran
12 to 25 percent above. Six alignments, four conventions and two real defects, and no mechanism
added to mmd2 at all.

### The ledger

Append only. A row is never edited once closed, so the sequence stays a record of what was wrong
rather than a summary written afterwards. The authoritative copy is in `mmd3.py`'s header and
mirrored in `mmd3.cpp`.

```
#  what diverged                where in ours          genmmd                nature
-  ---------------------------  ---------------------  --------------------  ----------
1  element expansion            mmd3_neighbors, I[u]   mmdelm, the el stack  convention
2  q2h walk                     the refresh            mmdupd, q2h           convention
3  qxh walk                     the refresh            mmdupd, qxh           convention
4  batch element order          the driver             genmmd, ehead         convention
5  merged weight in a           the refresh, q2h       mmdupd, dg -          DEFECT
   supervariable's bucket                              qsize[en] + 1
6  supervariable member order   the final expansion    mmdnum, the scan      cosmetic
```

All six closed 2026-08-07. **The nature column is the one to read first**, because the three
kinds carry entirely different consequences.

**DEFECT** means wrong on its own terms, with no appeal to genmmd needed. Only entry 5 qualifies:
it filed a supervariable one bucket too high per vertex merged into it, so it was never picked as
early as its size had earned, and the code did not do what its own comment, `dg - qsize[en] + 1`,
said. A defect found in mmd3 is a defect wherever the same code sits, so it was **fixed in mmd2
and in production `Mmd2` as well**. It had been costing fill in both since they were written:

```
grid          n      MMD2 fill before   MMD2 fill after
32x32      1024           +1.6%              -0.5%
100x100   10000          +16.2%              +5.9%
200x200   40000          +18.4%              +6.8%
400x400  160000          +25.2%              +8.3%
```

`mmd1` has no q2h path and no live merges, so it cannot have it.

**And the claim that stood here about the amd layers was WRONG, corrected 2026-08-08.** It read
that they file at an external degree, which excludes a vertex's own supervariable and therefore
does not move when its weight changes, so the shape could not arise. The external degree does not
move; the `- weight(u)` term inside the bound does, and it is the term that decides the bucket.
`amd2`, `Amd2` and `Amd2B` carried exactly this defect, subtracting the weight before the hash
merge that grows it, and it was costing 3 to 9 percent of fill on grids. Found by aligning `amd3`
against the vendored routine, where it is ledger entry 4, and fixed the same day;
`docs/DESIGN_DECISIONS.md` (2026-08-08) carries the account. `amd1` and `Amd1B` genuinely cannot
have it, having no live merges at all.

Worth keeping as a lesson about this ledger rather than only as a correction. The sentence was
written from a true premise, that AMD's degree is external, and it did not check whether the
quantity actually filed had the weight in it. A defect ruled out by an argument is ruled out only
as far as the argument reaches.

**CONVENTION** means neither side is wrong. A tie-break has no right answer, and entries 1 to 4
only change which of several equal-degree vertices wins. What they cost is not quality but
COMPARABILITY: while they differed, no measurement against genmmd could separate a difference of
mechanism from a difference of arbitrary choice. That is the whole reason mmd3 exists.

**COSMETIC** means it cannot change the answer at all. Entry 6 reorders the members of a
supervariable, which are indistinguishable by construction, so the fill and the elimination forest
are identical either way.

Worth noting for the next layer: **the defect was the hardest of the six to find**, because it
presented as a tie-break, which is what the other five were.

### Running it

```
make mmdorder
```

Production `MmdFlat` against `MmdCorrected`'s elimination order, on the same four shapes
`make amdorder` uses:
the seven examples, 2D grids from 4 a side to 140, 3D grids from 2 to 24, and nine random patterns
at n = 2000. 38 cases. `make aligned` runs this and the amd one together, which is the one word for
"is either ordering still what its reference computes".

### And the same check on real matrices, 2026-08-15

```
make mmdmatrices
```

`mmdmatrices.cpp` makes the identical assertion on matrices from the SuiteSparse Matrix
Collection, read through `benchmarks/matrices/MatrixMarket.h`. **246 matched, 0 differed, 0
skipped**, over `data/*/*.mtx` as fetched for the accuracy and performance sets.

**THAT RUN WAS AGAINST genmmd AND THE ASSERTION STILL HOLDS AGAINST `MmdCorrected`.** The reference
moved on 2026-08-23 and the 38 generated cases here were rerun and all match. The 246 have not been
rerun through this file, and did not need to be: `benchmarks/matrices`' `make mmdorder` makes the
same per-matrix comparison inside its own run, against `MmdCorrected` since that day, and its
2026-08-23 pass over the same 246 flagged nothing on any of the three mmd drivers. See
`benchmarks/matrices/ORDERING.md`.

The first run read 243 with 3 skipped, the three largest sitting past caps carried over from
`benchmarks/matrices`, which excludes them because a FACTORIZATION of them does not fit. That
reason does not reach an ordering: `PARSEC/Ga41As41H72`, `PARSEC/Si87H76` and `Schenk/nlpkkt80`
order in 10 to 14 seconds each and match like the rest, so the caps were raised to cover the
collection as fetched. They remain as `--max-n` and `--max-nnz` for bounding a run.

**Why it was worth building when the 38 cases were already green.** Every one of those 38 is
generated, and a generated matrix has the structure somebody chose to give it; widening a grid
exercises SCALE and never MECHANISM. That is not a hypothetical objection here, the 2D-only version
of the amd check having been green while production `Amd3` carried a stale clique degree that a 3D
grid at 16 a side finds. Real matrices bring what nothing in `graphs.h` produces, and the run swept
through all of it without a divergence:

- **The dense-row pathology.** `GHS_indef/bloweybq`, one column of degree 10000 among 9992 of
  degree 5, the matrix that takes MMD from 0.83 ms to 70.7.
- **Pure diagonals**, every vertex isolated: `Boeing/bcsstm39` at n = 46772 and nnz = 46772,
  `Cunningham/m3plates`, `JGD_BIBD/bibd_81_2`, five of the `HB/bcsstm*`, `Oberwolfach/t3dl_e`.
- **Graphs rather than meshes**, with degree distributions no grid has: `SNAP/email-Enron`,
  `Arenas/PGPgiantcompo`, `Pajek/Reuters911`, the `Gset` family, the `ML_Graph` nearest-neighbor
  graphs.
- **Near-identical siblings**, which is a test of TIE-BREAKING specifically: `Nemeth/nemeth02`
  through `nemeth09`, eight matrices of nearly the same structure.
- **Sizes two orders past the grid check**, to n = 100000 and nnz = 1.8 million.

**And the fill column settled a standing estimate.** `PARSEC/Si87H76` was on record as predicting
5.68 billion entries under MMD3, extrapolated from grids; measured, it is **5,679,875,732**. The
same run turned up a case nothing had recorded, `FlowIPM22/uni_chimera_i1` at **1,179,373,506**
from n = 100000 and nnz(A) = 1100592, which is 1072x nnz(A) from a matrix an eighth of `Si87H76`'s
size and a cheaper subject for the nested-dissection argument.

**Where the matrices come from, and what crosses the boundary.** They are not in the repository and
they are not named in any rule here: the driver takes them as ARGUMENTS, exactly as
`matrix_accuracy_cpp` does, because a Makefile naming another directory's data files is the case
`docs/WRITING_RULES.md` warns about, nothing binding them and the path dying silently. With `data/`
empty, which is the ordinary state, it says so and exits clean. The one thing that does cross is an
`#include` of the benchmark's reader, which is a compile-time dependency: if that header moves this
stops building, loudly.

**The reader gained one thing for this**, and its banner had already predicted the use: it takes
`pattern` files on request. They carry structure and no values, which is useless for a residual and
exactly an ordering's input. The default is still to refuse them.

**And it is deliberately NOT part of `make aligned`.** That target answers "is either ordering
still what the vendored routine computes" and passes on any machine; a target whose result depends
on what somebody has downloaded into a gitignored directory would make it mean something else.

**The amd side is the one to build next**, and unlike this one it should be expected to find
something: `docs/NEXT.md` item 6 records `Amd3` and the vendored AMD differing on fill on a minority
of the 107-matrix performance set, once by 4 percent on `HB/bcsstk08`, and calls it a divergence the
acceptance tests cannot see. `HB/bcsstk08` is in this run and matches on the mmd side.

**It needs no hook, and the asymmetry with amd is genmmd's rather than ours.** `mmd_order` returns
`perm`, the order genmmd eliminates in, and there is no postorder anywhere in the routine, so the
vendored output vector IS the object to compare. `AMD_2` hides its raw order behind
`AMD_postorder`, which is the entire reason `tools/hook_amd.py` exists. There is no Control array
here either, so the dense threshold that cost a day on the amd side has no counterpart to set
wrong; genmmd's one tunable is the `maxint` ceiling our wrapper supplies for the marker sweep, and
`REPORT.md` records that sweep firing zero times at every size we run.

**Until 2026-08-09 this check did not exist**, and that is the gap worth naming rather than the
target. The alignment was established by a scratch probe that lived in `/tmp` and died with the
session, exactly the arrangement `make amdorder` was built to replace on the other branch. What
ran day to day was `MMD3 nnzL == MMD nnzL` in the benchmark, and iteration 6 of `MMD3.md` is the
proof that this is not sufficient: fill was exact at every size while the permutation still
diverged at pivot 700 of 1024. A weaker check passing is not the stronger one passing.

### How it was actually done, including the wrong turns

The method is a loop and it is deliberately slow: run mmd3 against the vendored routine, find the
FIRST pivot where they differ, read both traces at that pivot, root-cause it, align that one
thing, record it, repeat. What follows is the sequence as it happened, wrong turns included,
because two of them were more instructive than the fixes.

**Starting point.** mmd2 matched the vendored MMD on 2 of the 7 examples, and on grids diverged
at pivot 60 of 1024 at 32 a side. MMD1 diverged at pivot 4 at every size from 3x3 up, which is
expected, since it lacks the mechanisms; mmd2 was the one that should have matched and did not.

**Choosing the case.** The smallest divergence was `graph1`, the 4-cycle, `n = 4`, where mmd2
differs at pivot 2 and MMD1 does not. Four vertices is readable by hand. mmd1 gave `[3, 1, 2, 0]`
and mmd2 `[3, 1, 0, 2]`: mmd1 mass-eliminates 0 into pivot 2, while mmd2's q2h pair merge folds 2
into 0 first, so 0 survives and is eliminated next. genmmd agrees with mmd1, so its merge keeps 2
and ours keeps 0.

**Entry 1 to 4, one defect found four times.** genmmd threads a linked list through an integer
array, pushes at the head, `list[nb] = h; h = nb`, and reads from the head, so the entry seen LAST
is processed FIRST. We hold a vector and append: the same set in the opposite order, at the same
cost. Only the winner among equals changes, and minimum degree is settled by exactly that.

Reversing the q2h walk alone took the examples from 2 of 7 to 3 of 7. That looked like slow
progress, and the next observation is what made it fast.

**The first wrong turn, and what it taught.** With q2h reversed, graph2 still diverged at the same
pivot as before. Instrumenting the q2h list showed why:

```
element 1   members=[0, 3]   reversed walk=[3, 0]   ->  3 survives
element 4   members=[3, 0]   reversed walk=[0, 3]   ->  0 survives
```

The same pair `{0, 3}` stored in OPPOSITE orders in two cliques. So reversing the walk gives a
different survivor depending on how that clique happened to get built, which meant the content
order of the lists was wrong upstream, not just the direction they were read in. Sorting the
element members was tried at that point and reached 5 of 7 examples. **It was rejected**, and the
reason matters: genmmd never sorts. It walks `adjncy` in whatever order `mmdelm`'s compaction left
it, which comes out ascending here only because the input is built from sorted columns and
compaction preserves relative order. Sorting reproduces the EFFECT on these graphs without
reproducing the MECHANISM, and costs a sort per element in a routine whose whole design rests on
not needing one. Two things were learned: a scoring improvement is not evidence of a correct port,
and the real problem was the order of `C[pivot]` itself.

That pointed straight at `mmdelm`, which builds the pivot's new list as explicit neighbors written
forward and then the ELEMENTS drained off a stack, `if (fwd[nb] < 0) { list[nb] = el; el = nb; }`.
Our `neighbors` walked `for c in I[u]` forward. Reversing that fourth walk took the examples from
4 of 7 to **7 of 7** in one step, and it is the deepest of the four: it fixes the order of
`C[pivot]`, hence the content order of every list the other three walk. The other three cannot be
judged without it, which is why they stalled at 4 of 7 on their own.

**The second wrong turn: the symptom is not diagnostic.** With all four reversed, grids still
parted company, at 67 to 78 percent of the pivots, a fraction FLAT in `n` rather than falling. The
smallest case was the 5x5 grid at pivot 19 of 25, where our bucket held `6: [13 11 1]`, three
vertices of equal degree, and we took the head, 13, where genmmd took 1. That reads exactly like
another tie-break, and it was not.

**Entry 5, found by instrumenting genmmd itself rather than reasoning.** Adding two `fprintf`
calls to a scratch copy of `Mmd.cpp`, one at each pivot and one at each merge, gave:

```
PIVOT 9  mdeg 6  qsize 1
PIVOT 5  mdeg 6  qsize 1
  PAIR 21 into 1
PIVOT 1  mdeg 5  qsize 2
```

genmmd files vertex 1 in bucket 5, we had it in bucket 6, and its true degree is the same, since
genmmd files at `dg - qsize + 1` and `5 = 6 - 2 + 1`. A supervariable is filed LOWER by its own
size minus one, so it is picked earlier than a singleton of equal degree. We merged 21 into 1 as
well, so the merge was not missing. What differed was WHEN the weight is subtracted:

```
ours     degree = dg0 - weight(u)          snapshot BEFORE the walk
genmmd   dg kept whole, then dg - qsize[en] + 1 at the END
```

Those agree until the walk MERGES a vertex into `u`, because genmmd's merge does
`qsize[en] += qsize[nd]` in that same walk, so the weight it subtracts is the POST-merge one.
Subtracting first files a supervariable one bucket too high per vertex merged into it, so it is
never picked as early as its size has earned. Moving the subtraction to the end matched 5x5 and
7x7 outright and **took the fill gap to zero at every size on the ladder**.

**Entry 6, and the check that reframed what was left.** Fill was now exact everywhere while the
permutation still diverged, at pivot 700 of 1024 at 32 a side. The right question was whether the
PIVOTS differed or only the printed order, and comparing pivot sequences rather than permutations
answered it: **identical at every size**, 788 of 788 at 32x32. So the algorithm was already
aligned. The remainder was the 6x6 tail:

```
vendored   ... 20, 12, 13, 15, 17, 22
mmd3       ... 20, 13, 12, 17, 22, 15
```

The same six vertices, and genmmd's trace showed them to be one supervariable: `PAIR 13/12/22/15
into 20`, then `PIVOT 20 mdeg 1 qsize 6`. What differed is the order its members are listed when
it expands. `mmdnum` numbers them root first, then by ASCENDING VERTEX INDEX, and gets that from a
single ascending scan, `for nd = 1..neqns, if perm[nd] <= 0`, walking to the root and taking the
next number. No sort: ascending falls out of the scan order. We listed them in merge order.

Entry 6 cannot change the fill and does not pretend to, since supervariable members are
indistinguishable by construction. It was closed anyway, because an exact permutation turns the
comparison into an EQUALITY TEST rather than a judgement, and that is the instrument the next
layer gets aligned with.

### What it bought

From `benchmarks/ordering`'s `make scale2d` on alpamayo, against the vendored MMD:

```
grid          n      MMD2 fill   MMD3 fill      MMD2 time   MMD3 time
32x32      1024        1.6%        0.0%           1.30x       1.29x
64x64      4096       12.0%        0.0%           1.17x       1.11x
100x100   10000       16.2%        0.0%           1.23x       1.13x
140x140   19600       21.6%        0.0%           1.26x       1.24x
200x200   40000       18.4%        0.0%           1.41x       1.39x
280x280   78400       21.2%        0.0%           1.44x       1.43x
400x400  160000       25.2%        0.0%           1.47x       1.49x
```

**The entire fill gap was convention and two off-by-one defects, not missing mechanism.** That was
not knowable before: with mmd2 the 20 percent could have been either, and the only way to tell was
to make the permutations equal and see what was left. What is left is time, 1.1x to 1.5x, and it
can now be profiled against an implementation known to compute the identical answer.

### Method notes worth keeping

- **Instrument the oracle, not just our side.** Entries 5 and 6 were both found by adding two
  `fprintf` lines to a scratch copy of the vendored source. Reasoning from our trace alone had
  produced a plausible and wrong diagnosis twice.
- **A better score is not a correct port.** Sorting beat three of the four reversals and was
  wrong. The question to ask of a candidate fix is which line of the vendored routine it
  corresponds to.
- **Compare the right object.** Comparing permutations said 70 percent aligned; comparing pivot
  sequences said fully aligned. The first framing would have sent us looking for a mechanism that
  was not missing.
- **The symptom does not identify the cause.** Entry 5 presented as a tie-break among three
  equal-degree vertices and was a filing defect. A different pivot at equal apparent degree is
  what BOTH look like.

### One production consequence worth knowing

Three of the four reversed walks are in `Mmd3.cpp`. The fourth, the `I[u]` expansion, lives in
`QuotientGraph::reachableSet`, which all six drivers share, so it is a flag, `setReverseIncidence`,
off by default and turned on only by `Mmd3`, with the branch hoisted and per clique rather than per
member. Entry 6 needed the same treatment for `QuotientGraph::order`, and became a second named
method, `orderAscending`, a counting layout with one ascending pass and no sort. A mode flag and a
parallel method on a shared class are not free to the reader; the alternatives were duplicating
the walks or changing the permutation for every driver.

mmd3 is in `LAYERS`, `GRID_LAYERS` and `PORTED`, so `make test` holds the twins to each other and
production to the prototype. Production carries `Mmd3` and `Ordering::MMD3`, and `make scale2d` in
`benchmarks/ordering` reports it beside MMD1 and MMD2.

## amd3, and the second alignment ledger

**amd3 adds no mechanism. It is amd2 with the vendored routine's list order, and it exists to
return `AMD_2`'s permutation.** This section is the durable record and holds the authoritative
ledger; `experiments/ordering/AMD3.md` is the narrative, iteration by iteration, with what each step established
and both of the author's reasoning errors at full length. `experiments/ordering/MMD3.md` is its counterpart on the
other branch. It is mmd3's counterpart, built the same way for the same reason,
and the digit means the same thing on both branches: 3 is the layer aligned to the vendored code.
The layer that used to carry this name, holding dense rows, `amd_aat`, the postorder and the
Control interface, is now `amd4` and is temporary.

**It worked, and it went further than mmd3's acceptance test asked for.** amd3 returns `AMD_2`'s
permutation exactly up to the postorder: not merely the pivot sequence, but the full expanded
order, member order within each supervariable included, on the seven examples and on every square
grid tested from 3 a side to 40, with pivot sequences checked to 50. Six alignments, four
conventions, one real defect and one that cost nothing but time, and no mechanism added to amd2 at
all.

### What the oracle is, and the one way it differs from mmd's

`genmmd` returns the order it eliminates in, so mmd3 could be held to its permutation directly.
`AMD_2` ends with a postorder of the assembly tree, which Oblio does not want and this layer does
not do, so a permutation comparison would fail on a relabeling rather than on a divergence.

Two things make the comparison exact anyway.

**The pivot selection and the supervariable finalization both sit inside `AMD_2`'s main loop**, and
`AMD_postorder` runs after it at line 2364, so the raw elimination order can be reconstructed
upstream of the relabeling: track membership alongside, a hash merge moving `j`'s members to `i`
and a mass elimination moving `i`'s to `me`, and emit each pivot's supervariable as its iteration
closes. Concatenating those is what the vendored routine would return if it stopped at the end of
its main loop, and it is what amd3 returns.

**And the knob that has to be set.** `Control[AMD_DENSE]` raised above `sqrt(n)`, which drives the
dense threshold into its `MIN(n, dense)` clamp so `deg > dense` is unreachable. Dense-row removal
is the one thing amd2 lacks that `AMD_2` does and that the Control array can switch off. Use that
rather than a negative alpha, which the header calls "off" but which leaves `dense = n - 2`, so a
vertex adjacent to everything is still removed. At the default neither matters here, the floor
being `MAX(16, dense)` against a grid's maximum degree of four.

`Control[AMD_AGGRESSIVE]` stays at its default, aggressive absorption being a mechanism amd2
already has.

### The ledger

Append only. A row is never edited once closed, so the sequence stays a record of what was wrong
rather than a summary written afterwards. The authoritative copy is in `amd3.py`'s header and
mirrored in `amd3.cpp`.

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
6  the exact test walked the new    the hash pass       for (p = Pe[j] + 1 ...),    COST
   clique, which entry 5 had just                       "skip the first element
   made a guaranteed match at                           in the list (me)"
   position zero
7  the stored clique degree not     beginElimination,   Degree [me] = degme,        DEFECT
   rewritten after mass             PRODUCTION Amd3     written TWICE, at its
   elimination trimmed the clique   alone               lines 1676 and 1940
8  the hash key's incidence half    the hash pass,      hval += e and hval += j,    COST
   annihilated by taking the        every amd layer     one sum with no stride,
   modulus to be the stride         and Amd2, Amd2B     then hval % n
```

**Entry 7 is production's alone, and that is the interesting half.** These prototypes obtain
`|C[c] - C[p]|` by walking the members of `C[c]` and counting the live ones outside `C[p]`, so
they recompute it from the truth at every step and nothing can go stale. Production maintains a
clique degree and reaches the same quantity by subtraction, which is amd2's pass 3, and only
production carries that encoding. `AMD_2` writes `Degree [me] = degme` twice, before scan 1 and
again after supervariable detection; the second write is the durable one, because by then scan 2
has run `degme -= nvi` for every vertex mass elimination took, so what a later step reads as
`|C[me]|` is the post-merge size. Production wrote it once, with the pre-merge value, so any
pivot that mass-eliminated left a clique degree permanently too large by the merged weight, and
every later bound taken through that clique inherited it.

It is half a mechanism again, exactly as entry 6 was. Ledger entry 3 moved mass elimination out
of the eliminator and did not carry the second write that the move is the whole reason for.
`Amd1` and `Amd2` cannot have it, mass-eliminating inside the eliminator, so their single write
already sees a trimmed clique.

**And the twin check could not have found it**, which is a limit rather than an oversight. A
prototype written to read as the algorithm does not carry the optimization, so it cannot model a
hazard that lives in one, which leaves the prototype-against-production comparison blind to
precisely the class of defect that optimization introduces. That is the divergence `REPORT.md`
parks as its fifth lead, and this is the first time it has cost anything. What found it was the
acceptance test widened to a shape the defect can move: an inflated bound changes an ordering
only when it changes the head of the minimum bucket, which no 2D grid does at any size to 140 a
side, and a 3D grid at 16 does.

**Entry 4's nature read `convention` in the prototype headers until 2026-08-09**, where this
ledger, `AMD3.md` and `docs/DESIGN_DECISIONS.md` have all said DEFECT since the day it closed.
The column is corrected there and the correction is dated rather than made silently: append-only
protects the record from being rewritten as a tidy summary afterwards, not from being wrong about
itself. Worth knowing that the copy called authoritative had drifted from its mirrors in the one
column this section tells a reader to look at first, and that nothing compares them.

**Entry 6 is a fourth NATURE and the word is needed.** It is neither convention nor defect nor
cosmetic: it changes no ordering, no fill and no permutation, only the COST. Entry 5 put the new
clique at the front of every `I[u]`, which is right; it did not also skip that entry in the exact
test, which `Amd.cpp` does in the same breath at `for (p = Pe[j] + 1 ...)`. So a guaranteed match
sat at the head of a short-circuiting walk and every FAILING pair paid one extra iteration of the
hottest line in the program. **Half a mechanism can be correct and still be wrong**, and nothing in
any output could have shown it: permutations, fill and every count matched.

**Entry 8 is the largest of the eight and the only one that was never a divergence.** Entries 1 to
5 are things we failed to reproduce; entry 8 is a defect of our own that `AMD_2` does not have. Our
key adds the incidence half with a stride, `(c + 1) * (n + 1)`, so that a vertex and a clique of
the same index cannot cancel, and then reduces modulo the same number, which annihilates that half
exactly. The hash was therefore a function of the ADJACENCY ALONE, and `A[u]` empties as the
elimination proceeds, so the key carried less and less and the buckets grew enormous.

**What it cost, against the vendored routine on the same graphs and for the same merges:** 19.0
pairs tested per pivot at 140 a side against its 0.333, and 155.3 at 26 cubed against its 0.484,
with the largest bucket 20 and 110 against essentially two. On alpamayo the fix takes `AMD2` at 26
a side from 14.88 ms to 5.45 and `AMD3` from 12.30 to 5.83, with the vendored routine and `AMD1`
unmoved to within the drift.

**The invariant the two lines hold TOGETHER**, which is how the replacement is stated: the modulus
must not divide the stride. `AMD_2` holds it by having no stride, and accepts the collision on
purpose, the hash being a filter and never the decision.

**It changes no output, which is why nothing here could see it.** Twins collide under any function
of the pattern, so the merges were the vendored routine's throughout. The prototypes carry the
identical key in both twins, so the twin check compared two files wrong the same way and the
prototype-against-production check inherited it. `docs/DESIGN_DECISIONS.md` (2026-08-09) carries
the full account of why five separate oracles were blind, and it is the more general form of what
entry 7 showed.

**The witness is a counter.** All three amd layers now print `hash pairs tested` beside `hash
merges`, for the reason `tag sweeps` and `bound below exact` exist: to make a claim checkable
rather than to measure anything. It should read about one pair per merge, and it reads 94 against
88 on a 20x20 grid.

**Production only: `Amd3` keeps its permutation and `Amd2` does not.** `make amdorder` still
matches on all 38 cases, and `Amd2` and `Amd2B` move, two-sided, `+1.4` percent of fill at 140x140
and `-3.1` at 26^3. The difference is where each driver last writes a bucket position: `Amd2` files
during the bound pass and the hash merge's refile is the last word, so the hash partition reaches
the degree buckets, while `Amd3` refiles every survivor afterwards in `pivotClique` order and comes
out canonical. That fourth pass is entry 4's, made for the post-merge weight, and it made `Amd3`
immune to this by accident.

**Four of the first five are one idiom, and it is mmd's idiom.** AMD pushes at the head where we
append,
so the entry seen last is processed first. Entry 1 is that in the hash buckets; entry 5 is it in
a variable's element list; entry 2 is which of the two sources is walked first, since
`for (knt1 = 1; knt1 <= elenme + 1; knt1++)` takes the elements and only its last pass the
supervariables; entry 3 is a placement rather than a direction, but it is the same kind of thing,
a step made in a different order.

**Entry 4 is the one that is not a convention**, and it is written up below.

### How it was actually done

The method is the loop set out above and it is the one that worked for mmd: run against the
vendored routine on the smallest case that still diverges, find the FIRST differing pivot, read
both traces with the ORACLE instrumented, root-cause it to a named line of `AMD_2`, change exactly
that one thing, record it, repeat.

```
                 examples matching
before               2 of 7        graph5 graph7
after entry 1        2 of 7        graph1 graph5
after entry 2        4 of 7        graph1 graph2 graph5 graph7
after entry 3        6 of 7        all but graph3
after entry 4        7 of 7        every example
after entry 5        7 of 7        and every grid tested
```

nnz(L) is exact against the oracle at every step and every size, so all five moved the tie-break
and nothing else.

**Entry 1 alone looked like a regression, and that was the information.** It closed graph1 and
opened graph7, leaving the score where it started. On graph7 both sides merge the same pair
`{0,4}` at the same pivot and keep opposite survivors, which under the new rule means their
`C[pivot]` orders differ. The oracle settled it: `LME 2 : 4 0` where ours was `[0, 4]`. So entry 1
was right and unjudgeable, exactly as mmd's three reversed walks stalled at 4 of 7 until the
element expansion was fixed beneath them. Entry 2 is that fix here, and it took the score to 4.

**Entry 3 closed two graphs at once, which its shape predicted.** graph4 and graph6 each matched
the oracle for their whole prefix and then took one extra pivot, five then a sixth and three then
a fourth. One cause was likelier than two, and it was the same cause.

**The examples were never the test.** At most twelve vertices each, they cannot exercise a
mechanism that needs real structure. With four entries in they were 7 of 7 while grids of side 8
and 10 still diverged. Grid 8 gave entry 5: 42 pivots agreed, then vendored took 17 and we took
16, both weight-2 supervariables at degree 6, and they were the SAME supervariable, the pair
`{16,17}` merged on both sides with opposite survivors. One pivot earlier, vendored had
`LME 2 : 5 12 19 10 16 17` against our `[16, 10, 17, 5, 12, 19]`: the same runs in a different
order of SOURCES, which is the element list's order rather than any walk direction.

### Entry 5 is not a prepend, and reading it mattered

The obvious reading of "the new clique goes to the front" is an insertion. `AMD_2` does something
else, because its two lists live back to back in one run and inserting at the front of the
elements would mean shifting everything right:

```c
Iw [pn] = Iw [p3] ;   /* move first supervariable to end of list */
Iw [p3] = Iw [p1] ;   /* move first element to end of element part of list */
Iw [p1] = me ;        /* add new element, me, to front of list. */
```

It lifts the two entries sitting at the boundaries to the two ends, so the elements come out
`[me, e2, ..., ek, e1]` and the supervariables `[j2, ..., jm, j1]`. A rotation of each list rather
than a shift of both. Both rotations are load-bearing here, because entry 2 made the reachable set
walk the cliques and then the adjacency, so both feed `C[pivot]`'s content order, and that order
decides entry 1's hash survivor.

Worth noting that production is a better fit for this than the prototype: `QuotientGraph` already
holds both lists in one run behind `mSourcePtr`, which is the layout that motivated the trick.

### Entry 4, the defect, and what it was costing

A hash merge folds `v` into a live `u` and grows `u`'s weight. The bound written moments earlier
has that weight subtracted inside it, in `degme - weight(u)` and in the `num_left - weight(u)`
cap, so absorbing `v` leaves the filed value one bucket too high per original vertex taken and `u`
is never picked as early as its size has earned. `AMD_2` has no such problem: it subtracts `nvi`
in the pass that restores the degree lists, which runs AFTER supervariable detection.

**It is mmd's entry 5 in a different array**, and this README used to say that could not happen
here, because AMD files at an external degree which does not move when a weight changes. The
external degree does not move; the `- nvi` term does, and it is the term that decides the bucket.

**A defect found in one place is a defect wherever the code sits**, so it was fixed in `amd2` and
in production `Amd2` and `Amd2B` as well, where it had been costing fill since they were written.
`amd1` and `Amd1B` cannot have it, having no live merges at all.

**SUPERSEDED A SECOND TIME, 2026-08-09.** The corrected AMD2 column below is itself now out of
date: ledger entry 8 fixed the hash key, which moves `Amd2`'s tie-break and therefore its fill,
to 11900, 199591 and 450190 at 32, 100 and 140 a side. `Amd3` is unaffected, its permutation
being unchanged. See `docs/DESIGN_DECISIONS.md` (2026-08-09).

```
grid        AMD (vendored)    AMD1      AMD2 before    AMD2 after
 32x32            11900      12074         12364         11900
100x100          206332     201856        212496        199386
140x140          474995     455472        487111        444191
```

**That reverses the fill half of REPORT finding 3**, which recorded AMD2's extras as a net loss
with the hash almost all of it. Corrected, AMD2 beats AMD1 at every size by 1 to 3 percent: the
extras were not costing fill, the filing was, and the hash was being charged for it. The time half
is untouched, the hash still being 72 to 92 percent of AMD2's time penalty. That attribution is
exactly what alignment was supposed to buy and what the method said to wait for.

### What it bought, and one result that goes the other way from mmd

The strongest single check is against numbers this experiment did not produce: amd3's nnz(L) comes
out 206332 at 100 a side and 474995 at 140, which is what `benchmarks/ordering/README.md` records
for the vendored AMD, digit for digit.

Beyond that, the same two things mmd3 bought. Any future divergence is a named pivot in a small
grid rather than a fill number somebody has to interpret. And the gap against amd2 is attributable
at last: whatever remains, it is not a difference of arbitrary choice.

**And one result that should not be smoothed over.** amd3 is aligned, so its fill IS the vendored
routine's, 474995 at 140 a side, where the corrected `Amd2` reaches 444191. On grids our tie-break
now beats AMD's by 6.5 percent. Aligning MMD improved our fill; aligning AMD costs it. That is the
second data point for the question `REPORT.md` parks, whether LIFO is genuinely better or
genmmd merely good, and it points the opposite way. Grids are one problem family and the
flattering one, so this wants the 3D grids `REPORT.md` asks for before it means anything.

### Method notes, and the two corrections this work forced

- **Instrument the oracle.** Every entry after the first was found by reading `AMD_2`'s own trace,
  never by reasoning from ours alone. The scratch copy lives in `/tmp` and `private/` is never
  edited; the method section above carries the recipe and the two hazards, brace counting and `goto` crossing an
  initialization.
- **Compare the right object.** With four entries in, the seven examples were exact while grids
  were not, and the permutation comparison was useless throughout because of the postorder. The
  raw-order reconstruction is what made the check an equality test.
- **A regression can be evidence.** Entry 1 lowered nothing and closed nothing on net, and it was
  correct. What it exposed was the layer beneath it.
- **`AMD_2` is the authority on `AMD_2`, including about itself.** This README said the shrinking
  `degme` is vendored behavior, so a survivor handled early sees a larger value than one handled
  late, and that our front-loaded mass elimination avoided a loss. `degme` is decremented in scan 2
  but never read there: the term enters a survivor's degree only in the later pass, by which point
  it is final. Every survivor sees the same number in both codes and there was nothing to avoid.
- **A defect ruled out by an argument is ruled out only as far as the argument reaches.** The
  entry-4 sentence, that the amd branch could not carry mmd's entry 5, was reasoned from a true
  premise and never checked against the quantity actually filed. It stood for a month.

### Where amd3 sits in the build

In `LAYERS` and `GRID_LAYERS`, so `make test` holds its twins to each other on the seven examples
and on grids, and in `PORTED` since production `Amd3` was extracted from it on 2026-08-08, so
production is held to it entry for entry as well.

**The extraction was a port rather than a copy**, three of the six entries landing in code the six
drivers share. `QuotientGraph` gained `setVendoredListOrder` for entries 2 and 5, grouped because
they are one fact and are only ever wanted together, and `setLateMassElimination` with
`massEliminate` for entry 3. Both flags are off for every other driver, and were verified inert
before `Amd3` existed at all: every suite passed and every fill figure across all nine orderings
was identical digit for digit. Entry 1 is `Amd3`'s own and entry 4 was already fixed everywhere.

**One defect the extraction introduced, and only a grid caught it.** Entry 5's rotation was written
as a `std::rotate`, which yields `[pivot, c1, ..., ck]` where AMD moves the FIRST element to the end
and so wants `[pivot, c2, ..., ck, c1]`, a swap of the two boundary entries rather than a shift of
everything. The seven examples agreed, grid 10 agreed, and grid 20 differed by one adjacent
transposition in 400 entries. That is the third time the prototype-against-production grid check has
paid for itself.

**`AMD3` is not the default.** `MMD3` became one because reproducing a decades-old reference beats
a tie-break of our own on unseen inputs; the same argument applies here and the evidence does not,
since the corrected `AMD2` fills less than the vendored routine on grids and `AMD3` therefore fills
more. Both stay.

## Grid mode, and what the missing features are actually worth

Every prototype takes an example number, and since 2026-08-07 every one also takes a grid:

```
./amd1_cpp grid 22        one 22x22 grid Laplacian, counters only
```

The grid is not an eighth example: nothing about it illustrates a mechanism and its trace is
unreadable. It exists so the counters can be read at a size the seven cannot reach.

**Nothing is filtered. What is not printed is not produced.** `SHOW_THRESHOLD`, a constant at the
top of every layer, decides it: above that `n` the run prints nothing from inside itself, no
initial state and no per-iteration trace, and below it everything prints as before. The guard sits
inside the `show` and `show_state` functions rather than at their call sites, so no site can be
missed. The threshold is 32, comfortably above the largest example at n = 12 and far below the
smallest grid at n = 100, so it never fires on anything read by eye. To watch a larger run, raise
it.

The rule it encodes is worth stating on its own: **per-iteration output is O(n) lines of O(n) each
and is for a human reading a small case; end-of-run output is O(1) lines and is what the twin
comparison uses.** So the first is bounded by a threshold and the second always prints.

**And, since 2026-08-03, so that the layers can be CHECKED at that size.** `make test` diffs the
whole output on sides 10 and 20: between the twins for every layer, and against production for
every layer in `PORTED`. Why the check is worth having is the amd2 subsection below: two defects
there left all seven examples byte for byte identical while the ordering was wrong on any grid of
10 a side or more.

**This replaced a filtering sink, on 2026-08-07, and the replacement fixed a defect.** Grid mode
used to run the trace in full and discard it through a whitelist of line prefixes, a `CounterSink`
in Python and a streambuf in the C++, one hardcoded list per file and ten files. Seven of those
keys matched nothing. `"nnz(L)"` was dead in all five layers that had the mode, because the line
begins `n = 100, nnz(L) = ...` and the test was `startswith`, so **grid mode had never once printed
the fill**, which is the number an ordering experiment exists to produce. `"degree computations"`
was dead in all three amd layers, which print `bound computations`. Nothing detected either, since
both twins dropped the same lines and agreed about output neither was producing.

**All thirteen layers have the mode, and the same seven examples.** It used to be five layers, in
the C++ only, so the twin check could compare only on the examples for those and not at all for the
rest. Both sides now spell it the same way, `python3 amd2.py grid 22` against `./amd2_cpp grid 22`,
sharing a `grid_graph` that must match vertex for vertex or the two would be diffing different
problems. The one deliberate exception is `matrix1` in amd4, which is an eighth run there and
nowhere else: it is malformed input, unsymmetric with duplicates and an unsorted column, and it
exists to exercise `amd4_preprocess` and `amd4_aat`. No other layer has an input path for it to
test, and cleaning it up to hand it over would delete the thing under test.

**Why it was added.** `benchmarks/ordering` shows our production MMD1 and AMD1 running about 4.5x
and 3.4x slower than the vendored routines. Two explanations were available and they call for
opposite work: our per-list allocation against their flat array, or the mechanisms these layers do
not yet have. Counters separate them, in units no allocator can move.

**Measured on a 100x100 grid, n = 10000.** The number that matters per branch is the dominant
inner quantity, not the degree update count.

```
                          mmd1      mmd2         amd1      amd2
degree computations      37240     23157        64413     56105
clique-member visits         .         .       272646         .      what amd1 pays per vertex
incidence entries            .         .            .     74281      what amd2 pays instead
nnz(L)                  223806    217102       201856    212496
```

**For AMD the gap is mostly algorithmic, and pass 3 is the whole of it.** `amd1` obtains
`|C[c] - C[p]|` by walking each touched clique's members; `amd2` obtains it by subtraction from a
maintained clique degree, walking incidence lists instead. That is 272646 elements against 74281,
a factor of 3.7, against a measured speed gap of 3.4. So porting pass 3 into AMD1 should close most
of it, and the earlier subsection records that the pass is output-neutral: same orders, same fill,
same looseness. Cost and nothing else.

**For MMD the gap is mostly ours.** `mmd2` runs 1.6 times fewer degree updates, and its `q2h`
shortcut also makes many of the remaining degree updates cheaper, which nothing here counts, so 1.6x
is a lower
bound on the mechanism's worth. Against a 4.5x speed gap that leaves a large remainder, and the
remainder is per-list allocation and pointer chasing.

**And one finding that was not what we expected.** MMD1's fill runs about 19 percent above the
vendored MMD at 100x100, and completing the algorithm barely moves it: `mmd2` gives 217102 against
`mmd1`'s 223806, where the vendored routine gives 186835. So three points of nineteen are the
missing features and sixteen are something else, which can only be the insertion history inside a
degree bucket, since that is the one thing our filing does not reproduce and the tie-break section
above says will differ. Minimum degree being famously sensitive to ties, a systematic 16 percent is
a great deal to attribute to one, and it is worth confirming rather than assuming.

The AMD side has no such problem: `amd1` at 201856 is *better* than the vendored 206332, and
`amd2` at 212496 is worse than both, which is the coarser-supervariables cost the pass 1 and 2
subsection already measured on small graphs.

## The tag guard, and where a sweep is allowed to land

The mark array is a set and the tag names it, so a tag must never repeat. Ours only ever climbs,
so the way it repeats is by overflowing its type, and the failure is the quiet kind: a wrapped tag
makes a stale stamp read as a match, so a vertex is counted as already seen when it was not. No
crash, no assertion, a slightly wrong ordering. Both vendored routines carry a guard against this,
`genmmd` at its line 42 and AMD in `clear_flag`, and the port had dropped it.

Every layer now carries one. The mechanism is three lines at each check site:

```
if tag >= TAG_CEILING:
    mark = [-1] * n          # 2n in amd2 and amd4, which stamp cliques at c + n
    tag = 0
```

**The sweep is unconditional, which neither vendored routine can afford.** `genmmd` sweeps
`if (marker[i] < maxint) marker[i] = 0` and AMD sweeps `if (W[x] != 0) W[x] = 1`, both selective
because both park a permanent state in the same array, `maxint` for a numbered vertex and zero for
an absorbed element, so neither may rest at its own sentinel. Ours carries nothing but tags, so the
sweep is a fill back to exactly what the constructor left and the invariant afterwards is the one
at startup. AMD's `wbig = Int_MAX_VAL - n` headroom is not needed either, since that exists because
`W[e]` can hold `wflg + size` and ours holds only tags.

**`TAG_CEILING` is `2^30 - 1`, and it is a placeholder rather than a derivation.** Nothing in these
files stores anything but a tag, so the true ceiling is the type's own maximum; half the positive
range leaves room against a later layer wanting some of it, and is far enough above the observed
rate that nothing fires. The two directions are not symmetric, which is worth knowing before the
constant is revised: too low costs one `O(n)` fill per sweep, amortized to nothing, while too high
is a correctness cliff with a silent wrong answer past it.

**Where a check may go is the whole of the difficulty, and it is per layer.** A check is legal only
where nothing in `mark` is live, and the eliminators are the obvious hazard: each holds
`pivot_clique_tag` and `absorbed_cliques_tag` live across the prune loop and then the merged set
across the `C[pivot]`
compaction, so the guard goes before the call and never inside. Several other regions hold a stamp
across a span the same way, and each was found by reading rather than assumed:

- `mmd2`'s refresh stamps `clique_tag` once per element and reads it throughout that element's
  q2h walk, where it decides both the pair merge and the outmatched case, with `vertex_tag` fresh
  per vertex nested inside it. Two levels live at once, which is `mmdupd`'s `mt` against its `tag`.
- the amd bound pass stamps `in_clique` and still reads it inside the `outside[c]` loop, with
  `seen_clique` stamped and consumed in between.
- `mdam2_refresh_bounds` stamps once and reads that stamp across all three of its passes.

So the count of sites varies, and the rule that produces it does not: **one check before each
region that advances the tag, placed where nothing is live.**

```
one site      md1, mda2                    only the eliminator advances a tag
two sites     md2 md3 md4 md5 mdm2 mdam2   the eliminator and the degree or bound pass
              mmd1 mmd2 amd1
three sites   amd2 amd4                    those two, plus the hash pair loop
```

**The third site in `amd2` and `amd4` is about a rate rather than about liveness.** Its `other` tag
advances once per PAIR TESTED rather than once per pass, and the pair count is quadratic in the
bucket sizes with no clean bound, so a check before the pass would leave the gap between two checks
unbounded. With the check at the top of each pair, every inter-check advance in every layer is
`O(n)` or better.

**The `tag sweeps` counter is the witness.** It is expected to read 0 at every size we run, so it
is there to make that claim checkable rather than to measure anything, which is the same reason
`bound below exact` sits beside it in the amd layers. Being a closing counter, it prints at every
size, including the grids, which is the only place the number could ever be interesting.

**Verified by forcing it.** With `TAG_CEILING` dropped to 2 or 3, so the sweep fires at nearly
every check, every layer's whole trace is identical to the normal-ceiling run apart from the
counter itself, in both twins, and the twins still agree with each other. On grids the forced
sweeps run from 98 for `mmd2` to 2028 for `amd2` at side 20, and `bound below exact` stays 0
throughout, which is the assertion that would catch a sweep landing inside the bound pass.


### What the ceiling should be, worked out one layer at a time

The code uses `2^30 - 1` everywhere and continues to. This subsection is the investigation
running ahead of it: what the ceiling would be if it were derived rather than chosen, filled in
as each layer is worked through. Nothing here has been applied.

**The quantity to bound is the advance between two consecutive checks, not anything stored.** The
test is `tag > TAG_CEILING` at the START of a region, so a tag that just passes the test still
climbs through that whole region before the next check. If the region advances the tag by at most
`A`, the largest value ever reached is `TAG_CEILING + A`, and the requirement is

```
TAG_CEILING = INT32_MAX - max over regions of A
```

which is AMD's `wbig = Int_MAX_VAL - n` at `Amd.cpp:1349` when the maximum is `n`. Note the
maximum over ALL of a layer's regions, not one of them: a layer whose regions cost `n` and `3`
needs `max(n, 3)`, and the constant is not redundant, since at `n = 2` it is the constant that
binds and `INT32_MAX - n` overflows.

```
              eliminate   degree work   hash pairs   ceiling
md1           n           N/A           N/A          INT32_MAX - n
md2           3           n             N/A          INT32_MAX - max(n, 3)
md3           3 or 4      n             N/A          INT32_MAX - max(n, 4)
md4           3 or 4      n             N/A          INT32_MAX - max(n, 4)
md5           3 or 4      n             N/A          INT32_MAX - max(n, 4)
mdm2          3           n             N/A          INT32_MAX - max(n, 3)
mda2          ?           N/A           N/A          ?
mdam2         ?           ?             N/A          ?
mmd1          ?           ?             N/A          ?
mmd2          ?           ?             N/A          ?
amd1          ?           ?             N/A          ?
amd2          ?           ?             ?            ?
amd4          ?           ?             ?            ?
```

**The n in the middle column is rounded, and from three different exact values.** Worth keeping
straight, since the rounding is what makes the rows look alike:

- **md1, exactly n - 1.** Its eliminate advances once per NEIGHBOR of the pivot, so the bound is
  |A[pivot]|, largest at the first elimination of a complete graph. It appears in the eliminate
  column rather than the middle one because md1 has no degree pass: the picker reads
  `A[u].size()` and spends no tag.
- **md2 and md3, exactly n.** The pivot search calls the neighbors function once per LIVE VERTEX,
  and at the first iteration every one of the n is live.
- **mdm2, exactly n - 1, and ATTAINED.** Maintenance replaces the pivot search with a refresh
  over the members of C[pivot], and a reach excludes the pivot itself, so |C[pivot]| <= n - 1.
  mdm2 has no mass elimination, so a first pivot adjacent to everything gives a refresh set of
  exactly n - 1.
- **md4 and md5, an UPPER BOUND of n - 1 that is NOT attained**, and the entry is n for safety
  rather than because anything reaches it. Their refresh set is C[pivot] copied AFTER the
  eliminator, so the merged vertices are already gone from it, and the two quantities are coupled:
  pushing |C[pivot]| toward n - 1 is exactly what forces its members to merge.

  At the first elimination this is provable. There are no cliques yet, so I[u] is empty for every
  u and the prune leaves I[u] == {pivot}, which satisfies the second merge conjunct for every
  member automatically. If |C[pivot]| = n - 1 then C[pivot] + {pivot} is every vertex, so
  A[u] - C[pivot] - {pivot} is empty for every member and the first conjunct holds too. All n - 1
  merge and the refresh set is EMPTY. A complete graph and a star both do this.

  Later steps are harder and not worth settling. A construction that keeps the refresh set large
  must give every member of C[pivot] a neighbor outside the clique AND have the pivot be a
  minimum-degree vertex at that moment, and those pull against each other. The exact maximum is a
  real combinatorial question whose answer nobody would use: the column records what the ceiling
  must cover, so an upper bound is the correct entry and a loose one costs one unit out of two
  billion.

Rounding all of them to n costs nothing and keeps every row in `Amd.cpp`'s own form,
`Int_MAX_VAL - n`, rather than introducing an n - 1 a reader would have to re-derive to trust. It
does hide the layer's whole point in the mdm2 row: the refresh saves one advance in the WORST case
and a great deal in the average, and a column recording maxima cannot show that.

**Where a column reads "3 or 4" the advance is CONDITIONAL**, and the ceiling takes the larger.
The pair is kept in the table rather than just the maximum because the two cases are different
facts about the layer: md2 and mdm2 advance exactly 3 every elimination, md3 onward advance 3 when
nothing mass-eliminates and 4 when something does.

**The eliminate column is where the layers actually differ**, and it tracks what a tag is ABOUT.
md1's tags are about each neighbor in turn, one per neighbor because the fill is pairwise and each
neighbor's missing edges are its own. From md2 on, the quotient graph gives the fill a name, so
the eliminator stamps sets belonging to the PIVOT: its reach, the members of the clique that reach
becomes, and the ids of the cliques it absorbs. Three, whatever the pivot's degree. md3 adds a
fourth, the merged vertices of mass elimination, and that one is conditional, so its row is a
bound rather than an exact count. Each layer's eliminator carries the same account in its own
comments.

**md2 and mdm2 could each be 2 rather than 3**, and are not. The eliminator's two tags stamp
opposite sides of a clique, `pivotCliqueTag` its members and `absorbedCliquesTag` the absorbed
ids, and the two prune loops query opposite sides, so one value would still name each set
correctly. Two are used so that each test is exact by the tag alone rather than by an invariant
about which lists hold live vertices and which hold eliminated ones. Parked at md3, whose merged
set is on the MEMBER side alongside C[pivot] and so is the first case where two tags would share
one.

The eliminate column is the only universal one, every layer having a check there. The middle
column is one idea under two names, the pivot search in md2 and md3 and the refresh from md4 on,
which is the same degree work moved to the other side of the pivot; `mda2` has no entry because
its bound reads lengths and takes no tag, which is what the recomputing column of the square
means here.

**md1 spends `n` in the eliminator.** `tag += 1` sits inside the loop over the pivot's neighbors,
because each neighbor's adjacency has to be tested against the pivot's set separately, so one
call advances by `|A[pivot]|`, at most `n - 1`. Rounding to `n` is what leaves the margin: the
worst case then lands at `INT32_MAX - 1` rather than exactly on `INT32_MAX`.

**md2 spends 3 in the eliminator and `n` in the search, and the swap is the point of the layer.**
The eliminator advances three times and none of them is in a loop: once inside `md2_neighbors`
for the reachable set, then `pivot_clique_tag`, then `absorbed_cliques_tag`. Every loop after that
only reads `mark` against those two stamps. That is the quotient graph paying off, one set built
and many membership tests against it, where md1 built one set per neighbor. The cost moves to the
pivot search, which calls `md2_neighbors` once per live vertex, so up to `n`.

**md3 spends 3 OR 4 in the eliminator and `n` in the search, and only the fourth is new.** Its
neighbors function and its pivot search are BYTE-IDENTICAL to md2's, in both twins, so the middle
column is the same by construction rather than by coincidence: one `tag += 1` at the top of the
call, one call per live vertex, `n` at the first iteration.

Its eliminator is md2's too, for the first 32 statements. Stripping comments, the two are the same
line for line from the signature through the prune loop, and the first divergence is md3's
`merged_vertices` declaration; the five statements after the mass elimination block are shared
again. So the first three advances are md2's unchanged, in the same order and for the same sets:
the reach inside the neighbors call, then `pivot_clique_tag`, then `absorbed_cliques_tag`.

**The fourth is conditional, and that is what makes the row a bound.** It fires inside
`if merged_vertices`, so an elimination that merges nothing advances 3 and one that merges
advances 4. The table takes the maximum, which is what the ceiling has to cover, but the two cases
are worth naming because md2's and mdm2's 3 is EXACT where this is not. On grids the fourth fires
at most eliminations; on a path graph it never fires and md3's eliminator costs exactly what md2's
does.

It is also the first tag on the MEMBER side alongside `pivot_clique_tag`, which is what ends the
two-into-one saving described above.

**md4 and md5 have the SAME eliminator as md3, verified rather than assumed**: 53 statements each,
identical to md3's after stripping comments, with every difference being a layer-name prefix or a
tag name this sweep has not reached yet. Four `++tag` sites in each of the three files. So both
rows will read 3 or 4 in the eliminate column, and their middle column is the only thing left to
work out: md4 replaces the pivot search with a refresh the way mdm2 does for md2, so `n - 1` is
what to expect there and to check.

**mdm2 spends 3 in the eliminator and `n - 1` in the refresh, and the eliminator is md2's
verbatim.** Its four shared functions are byte-identical to md2's after the name substitution, so
maintenance changes nothing but where the degree work happens: the picker reads cached integers
and spends no tag, and the `tag += 1` per call moves to a refresh over the members of `C[pivot]`.
A reach excludes the pivot, so that is `n - 1`.

**md4 is that same relationship one layer up, and it was checked rather than assumed.** Its
neighbors function and its eliminator are byte-identical to md3's in both twins, the only
difference anywhere being one docstring line saying so. So it inherits md3's 3 or 4, and its
middle column is mdm2's `n - 1`: the picker reads cached degrees and spends no tag, and the
refresh walks `refreshed_vertices`, a copy of `C[pivot]` taken AFTER the eliminator and therefore
already trimmed of the merged vertices. Rounded, its row is md3's.

**The mmd layers have SIX `++tag` sites, not four, and their rows are not worked out.** Recorded
2026-08-12 from reading, not from finishing the argument. Four are the eliminator's, unchanged
from md3. The other two are in the degree refresh, and they are the first case in the family where
TWO TAG LEVELS ARE LIVE AT ONCE:

- `cliqueTag`, stamped once per element of the batch and read all the way through that element's
  q2h walk, where it decides both the pair merge and the outmatched case;
- `vertexTag`, fresh per vertex and nested inside it, so one q2h vertex cannot hide a neighbor
  from the next.

That is mmdupd's `mt` against its `tag`, and it is why a sweep inside an element would erase marks
about to be read: the guard sits before the element loop rather than inside it. The eliminator's
guard sits in the batch loop, which is a second thing to settle, since a batch performs several
eliminations between two visits to it.

**So mmd1, mmd2 and mmd3 need two regions reasoned about rather than one**, and the second has a
nesting the md layers never had. Their rows stay open until that is done. What is already known:
mmd2 and mmd3 also read `buckets.filed`, which no md layer does, because a batch WITHHOLDS
outmatched vertices rather than refiling them.

**md5 spends 3 or 4 in the eliminator and at most n in the refresh, like md4.** Its eliminator and
neighbors function are byte-identical to md3's, and its refresh set is built the same way md4's
is, `const std::vector<std::int32_t> refreshedVertices = C[pivot];` taken after the eliminator.
What md5 changes is how the MINIMUM is found, buckets instead of a linear scan, and the picker
walks buckets rather than calling the neighbors function, so it spends no tag either way.

**So the maintenance idea costs the same wherever it lands.** md2 -> mdm2 and md3 -> md4 are the
same edit and move the same `tag += 1` from one region to the other, which is why both pairs of
rows differ only by an n against an n - 1 that the rounding then erases. The one asymmetry is
which layers ATTAIN their bound: mdm2 does, md4 and md5 do not, and the reason is mass
elimination, the very mechanism that makes their eliminate column a bound rather than an exact
count.

**One consequence to note before this is applied.** The ceiling column differs by layer, so
`TAG_CEILING` stops being one shared constant and becomes thirteen expressions. Today every file
carries the same literal, which is exactly why the pragmatic value is easy to keep and the
derived one will not be a search and replace.

## Four bugs this found, all ours

Worth recording, because every one was invisible to the checks in place at the time.

**Two of them are amd2's and share one cause**, so they are written up together in the amd2
section below, under "Two defects the permutation check found, both inherited from amd1": the
bound's live-vertex cap taken from the wrong counter, which drove the bound below the true degree,
and the fill accounting taking an unweighted count of the new clique. Both were correct lines in
amd1 and became wrong the moment amd2 added a merge into a live vertex.

**amd did not shrink the new clique on mass elimination.** When a vertex is mass-eliminated it
joins the pivot's supervariable, so it stops being outside the new clique and must stop
contributing to `|C[p]|`; `Amd.cpp` does this at `degme -= nvi`. We computed `|C[p]|` once before
the loop. The effect is nearly invisible: identical results on all four test graphs and on every
grid, surfacing only on a five-vertex bowtie where a bound came out one too large.

**`mmd1.cpp` printed display lines its Python twin did not**, left over from the `md5` file it was
derived from. This survived because the verification at the time used a `grep` filter narrow
enough to skip exactly those lines, which is not a test. `make test` exists because of this one:
it compares whole outputs, and it found the drift immediately.

## Two axes, and where the chain runs through them

Everything below md2 sits in a two-dimensional space, and the ladder walks it in a particular
way that is worth stating before the layers rather than discovering across them.

**The first axis is representation, and it is what md1 through md5 climb.** md1 materializes fill,
md2 stops materializing it, md3 merges indistinguishable vertices, md4 stops recomputing degrees
that did not change, md5 stops scanning for the minimum. Each of those changes what is stored or
when it is touched. None of them changes what a degree MEANS.

**The second axis is the degree itself.** From md2 on, a degree is the size of a union of an
explicit adjacency with some cliques, and a union can either be counted or estimated. Counting is
exact and costs a pass over every source. Estimating decomposes the union, adds the sizes of the
pieces, and overcounts where they overlap, at a cost of one number per clique. That estimate is
bound(u), and it is derived in the md2 subsection below, since md2 is where the union it
decomposes first exists.

**The estimate itself never changes.** One decomposition, written once, unaltered from md2 through
amd2. md3 changes what a size means, but identically for degree(u) and for bound(u), so the
relation between them is untouched. md4 adds a cap. md5 does nothing to it at all. So the second
axis is not four variants of an idea; it is one definition and a question of which layer can use
it.

**And the answer to that is md4, not md2.** bound(u) decomposes reach(u) against the new clique, so
it is a statement about the vertices that clique reached and about no others. md2 and md3 repick
from scratch at every iteration, so most of what they update is outside the bound's domain and would
have to be counted exactly anyway. md4 is the layer that narrows the degree update set to exactly
C[p], which is exactly the bound's domain. The bound is inherently incremental, md4 is where the
algorithm becomes incremental, and those two facts meeting is why amd belongs at the top of the
chain rather than beside md2.

So the grid is smaller than it first looks: two exact rungs where the estimate is expressible but
not usable, then a genuine choice from md4 up. Each layer from md2 carries a short subsection
saying where it stands, and amd1 builds the approximate cell once, at the top, as running code.

**And the axis stops being orthogonal at the top, which is why the fork is real.** mmd and amd are
not two knobs that happen to be turned separately. They attack the same cost from opposite ends,
mmd making the degree update rare and amd making it cheap, so their gains overlap rather than add.
They also interfere: the anchored bound is anchored at ONE new clique, and a batch produces several,
so
the anchoring would have to be redone against their union; and the bound stays tight only when
degrees are updated often, which is exactly what batching gives up. That is measured here rather
than assumed, and the measurement is in the delta section. So md1 through md5 is one chain, and
the fork at the top into mmd, which stays exact, and amd, which goes approximate and does not
batch, is a real division rather than a filing convention.

## What the 1996 AMD paper covers: the degree bound in full, the layout in one paragraph

The paper is Amestoy, Davis and Duff, *An Approximate Minimum Degree Ordering Algorithm*, SIAM J.
Matrix Analysis & Applications 17(4), 886 to 905, and there is a copy in the project. It is worth
knowing what it does and does not settle, because several things this ladder treats as `AMD_2`
implementation detail are in fact in the paper, and several things it treats as published are not.

Quotations below use the paper's vocabulary, so *element* is its word for our clique and Lp is the
clique the pivot p forms. See "One word, before anything else" above.

**WHAT THE PAPER IS ABOUT.** The approximate degree bound, its derivation, the proof that
d <= dhat <= dtilde, and the measurements. Sections 4 and 5 develop the bound and Algorithm 2
computes |Le \ Lp|; section 6 is thirty pages of results across four codes. Supervariable
detection by hash is there too, with the hash function written out and the mod (n - 1) term. That
is the paper's subject and it occupies almost all of it.

**WHAT IT SAYS ABOUT STORAGE IS ONE PARAGRAPH**, in section 5 on page 17, immediately after the
algorithm is stated. It is credited to MA27 rather than presented as new: the data structure is
"the quotient graph data structure used in the MA27 minimum degree algorithm". In order, that
paragraph gives:

- the sets A are stored first, followed by a small amount of elbow room;
- when Lp is formed it goes into the elbow room, **or in place of Ap when |Ep| = 0**;
- garbage collection occurs if the elbow room is exhausted;
- during it, the space for Ai and Ei is reduced to exactly |Ai| + |Ei| per supervariable, the space
  for Ae and Ee of every element is fully reclaimed, and so is Le of any absorbed element;
- each collection takes time proportional to the size of the workspace, normally theta(|A|);
- **"In practice, elbow room of size n is sufficient."**

**So two things we had been reading as `AMD_2`'s own choices are the published design.** The
in-place branch, which our drivers spell `elenme == 0` and which takes 62 to 68 per cent of
eliminations on a grid, is in that paragraph. So is the size-n elbow room, which is the `n` term in
`iwlen = nzaat + nzaat/5 + n` and which the code enforces as a floor while treating the `nzaat/5`
as a recommendation.

**AND WHAT IS NOT IN THE PAPER IS EVERYTHING ABOUT HOW THE COLLECTION RUNS.** No FLIP encoding, no
parking of a displaced head, no two-pass sweep, no carrying of a half-built clique. The paragraph
says what the space looks like before and after and says nothing about the mechanism between. All
of that is `AMD_2`, and porting it was reading code rather than reading the paper.

Nor is the run order in the paper: that I[u] comes before A[u], that the new clique goes to the
FRONT of I[u], or the three-move rotation that puts it there. Those decide which of several
equal-degree vertices is picked, so they decide the permutation, and the paper is silent on all
three.

### The complexity bound rests on an assumption about compaction

Page 18 states the bound and then names two conditions it rests on: few or no hash collisions in
supervariable detection, and **a constant number of garbage collections**. The authors say the
assumptions seem to hold in practice and that the asymptotic time would be higher otherwise.

That is a testable claim rather than a theorem, and the `AMD cmp` column in
`benchmarks/matrices/ORDERING.md` is the test, over 246 matrices the authors did not have.

Their own measurement is on page 23, with elbow room of size n: usually one collection, at most two
for AMD and at most three for the other codes, the sentence finishing across the table on page 25
by crediting the difference to aggressive absorption reducing the memory required. They also state
flatly that "garbage collection has little effect on the ordering time obtained".

Ours agrees on the shape and disagrees on the tail: a median of one, 122 of 246 never compacting at
all, but a maximum of ten rather than two. The matrices that do it are interior-point and social
graphs, classes that did not exist in the 1996 test set.

### Chaining gets half a sentence, and no measurement

The one place the paper mentions the other layout is page 25, comparing codes, where MMD is
described as storing the element patterns "in a fragmented manner" and requiring no elbow room.
That is chaining, named accurately and passed over. Nothing follows about what the fragmentation
costs, and nothing could, since the tables that follow compare whole codes and confound the layout
with the algorithm.

**That gap is what `Mmd3B` exists to fill.** Same algorithm, same encodings, same C++, chained
storage against ours, which is the comparison the paper's own tables cannot make. The result is in
`benchmarks/matrices/ORDERING.md`.

## The whole inventory: three families of mechanism, and which branch has what

Minimum degree spends nearly all its time refreshing degrees after a pivot. Everything either
branch adds beyond the naive algorithm attacks that cost, and every one of those additions falls
into one of three families. The families are worth having because they answer "is that all of
them?" structurally rather than by enumeration.

| family | operates on | mechanisms |
|---|---|---|
| detecting indistinguishability | vertex identity | mass elimination, q2h, hash, pre-compression |
| removing redundancy | cliques | natural absorption, aggressive absorption |
| scheduling the work | the refresh itself | multiple elimination, incomplete update, the bound |

Nothing in either code sits outside those three. "Two axes, and where the chain runs through them"
above places the same material on the representation and exactness axes and follows the ladder
through it; this section is the flat inventory.

### Family 1: detecting indistinguishability

Two vertices are INDISTINGUISHABLE when each sees exactly what the other sees plus the other
itself, `reach(u) + {u} == reach(v) + {v}`. Eliminating one right after the other then creates no
fill whatever, so they can be taken as a single pivot, a SUPERVARIABLE. The md3 section below
derives it.

**Four mechanisms find such pairs, and they differ only in AGAINST WHOM and WHEN.**

| mechanism | against whom | when | who has it |
|---|---|---|---|
| mass elimination | the pivot, which is leaving now | at every elimination | everyone from md3 up |
| q2h | another live reached vertex | during the refresh | mmd2, mmd3 |
| hash | another live vertex in `C[p]` | after the prune | amd2, amd3, MA27 |
| pre-compression | another vertex in the ORIGINAL A | once, before elimination | Ashcraft, CMMD |

**Only mass elimination produces pivots directly**, its partner being the vertex leaving in the
same breath. The other three build supervariables that are eliminated together later, so they
change the answer only through the degrees they alter, which is why coarser is not automatically
better.

The two live-vertex routes are the ones that differ by branch, and "Two ways a vertex disappears"
below has that comparison; "Supervariable detection in production" has the mechanism for each.

### Family 2: removing redundancy

Nothing here is about vertex identity. A clique that can never contribute again is pure cost in
every incidence list that names it, and both mechanisms delete such cliques.

| mechanism | what dies | who has it |
|---|---|---|
| natural absorption | every clique in `I[pivot]`, inside `C[pivot]` by construction | everyone |
| aggressive absorption | any clique wholly inside `C[pivot]`, touched or not | amd, MA27 |

### Family 3: scheduling the work

Nothing here merges or deletes anything. These decide how much degree computation is done and when.

| mechanism | what it decides | who has it |
|---|---|---|
| the approximate bound | how expensive ONE refresh is | amd, and it is amd's whole idea |
| multiple elimination | how many pivots go before a refresh | mmd |
| incomplete update | which reached vertices that refresh bothers with | mmd |

**The two branches take opposite answers within this family, and that is the headline difference.**

**MMD REFRESHES LESS OFTEN, AND REFRESHES FEWER.** It keeps the degree EXACT and attacks the cost
from two directions. Both come from Liu, *Modification of the minimum-degree algorithm by
multiple elimination*, ACM TOMS 11 (1985), 141 to 153, which the AMD paper cites as [25]:

- **Multiple elimination.** A whole independent set of minimum-degree pivots goes in one batch, and
  the refresh that follows is amortised over all of them. The Multiple in Multiple Minimum Degree.
- **Incomplete degree update.** Within the refresh, a reached vertex whose reach lies inside
  another reached vertex's cannot be the strict minimum before that one is eliminated, so its
  degree is not computed AT ALL. It is withheld rather than refiled, and put back only when a later
  elimination reaches it.

**AMD MAKES EACH REFRESH CHEAPER INSTEAD.** It eliminates one pivot at a time, refreshes every
vertex of `C[p]`, and replaces the exact degree with a bound computed by decomposition, so a clique
is never opened and never walked. The approximation is the whole of the idea and everything else
follows from it.

The two are mutually exclusive in practice: batching needs exact degrees to know the set is
independent, and the bound is cheap precisely because it does not open the cliques a batch would
have to.

### When each one runs, and why the order is forced

The families say what each mechanism operates on. They say nothing about when it can run, and that
is decided by something else: what has to have been computed already.

**Everything except pre-compression happens inside a single pivot's step.** The order within that
step differs between the branches:

```
amd3, one iteration            mmd3, one iteration
----------------------------   ----------------------------
form C[p], NATURAL ABSORPTION  form C[p], NATURAL ABSORPTION
prune: the bound               prune
  AGGRESSIVE ABSORPTION        MASS ELIMINATION
MASS ELIMINATION               refresh over C[p]
hash detection: PAIRWISE         Q2H PAIRWISE
trim                             OUTMATCHING
```

**Each position is forced by a prerequisite, and naming them is more useful than the timeline.**

- **Natural absorption** needs only `I[pivot]`, which is in hand before anything, so it goes first
  and could not go later: the walk that builds `C[p]` consumes that list.
- **Aggressive absorption** needs `|C[c] - C[p]|` for every clique the new one touched, and that is
  exactly what the bound's walk produces. It cannot run before that walk and there is no reason to
  wait after it.
- **Mass elimination** needs the pruned lists, its test being that nothing explicit is left and the
  only clique is the new one.
- **q2h and hash** need the pruned state of `C[p]` too, and amd's additionally needs the keys the
  bound pass accumulated.
- **Pre-compression** needs only `A`. That is why it can run once before everything, and also why
  it can only ever find pairs that were alike in `A` already.

**And one ordering is a correctness constraint rather than a convenience.** Amd's mass elimination
runs in the driver, after aggressive absorption, and not inside the eliminator where every other
driver here puts it. Absorption is what makes the cheap structural test agree with the true one: a
clique whose members all lie inside `C[p]` contributes nothing to what `u` can reach, yet its
presence in `I[u]` makes the test fail. `Amd.cpp` relies on the same thing, making the test in its
scan 2 over a clique list absorption has already compacted. Asking first declines merges the
vendored routine makes; see `AMD3.md`, ledger entry 3.

### The families are not independent of each other

They are separable as ideas and coupled through one object, the REACH SET.

Family 3's choice decides whether the reach set is ever formed, and that decides what families 1
and 2 can afford. An exact degree materializes it, so q2h merging and outmatching fall out of a
walk that was happening anyway. A bound decomposes it and never forms it, so amd must go looking
for its pairs, which is what the hash is for, and cannot outmatch at all.

**Aggressive absorption is the exception that shows the coupling is not a rule.** It looks like a
consequence of the bound, and it is not: MA27 does it with true degrees. See its subsection below.

### Incomplete update is outmatching, and it is not the batch

The two are easy to conflate because both are Liu's and both defer work, but they are independent.
Multiple elimination decides HOW MANY pivots go before a refresh; incomplete update decides WHICH
of the reached vertices that refresh bothers with. You could batch and still refresh every reached
vertex, and you could outmatch one pivot at a time.

**The eviction during a batch is the first, not the second.** A reached vertex is unfiled so it
cannot be picked again before its degree is known, and it is refiled at the end of the round. That
is bookkeeping for the batch.

**Outmatching is the second, and it is stronger than deferral.** If `reach(v)` is contained in
`reach(u) + {u}`, then `v` cannot be the strict minimum before `u` goes, so its degree is not
computed for that round at all; if it stays outmatched across several rounds it is not computed for
those either. genmmd writes `bwd[nd] = -maxint` and our `Buckets` spells it `outmatch`, with
`restore`, `mmdelm`'s `bwd[rn] = 0`, putting it back.

**And that is what the qxh list is for.** The reached vertices divide in two: `q2h`, cheap enough
to decide exactly and where pairwise merging happens, and `qxh`, where outmatching happens instead.
Both lists exist because of this split, not because of the batch.

### Amd has no incomplete update, and hashing is not the reason

It has none: no amd driver calls `outmatch`, and neither does `mmd1`, the mechanism arriving with
`mmd2` alongside the q2h and qxh split.

**Two reasons, and the second is the same root cause as q2h.** The bound is cheap, so there is
little to be saved by skipping one; and detecting that `v` is outmatched needs `reach(v)` and
`reach(u)` compared, which amd never forms. It has bounds, not reaches. The same wall that stops
amd deciding a q2h merge stops it deciding an outmatch, and for the same reason: a bound never
opens a clique.

So hashing is not what displaces incomplete update. Hashing replaces the q2h half of the split;
nothing replaces the qxh half, because amd computes every vertex of `C[p]` every time and is fast
enough not to mind.

### Natural absorption, which both branches have and neither calls by that name

When `C[pivot]` is formed, every clique in `I[pivot]` is inside it by construction and can never
contribute again. Page 6 of the 1996 paper: the elements adjacent to the pivot are absorbed into
the new element and deleted, and reference to them is replaced by reference to the new one. In our
code it is the `killClique` loop in `beginElimination`, and every driver on both branches runs it.

**Mmd has this and nothing else.** Its cliques die only when a pivot touches them.

### Aggressive absorption, which amd has and mmd does not

Page 17 states it: on top of the natural absorption of the elements in Ep, any element whose
external subset is empty, |Le \ Lp| = 0, is absorbed into p as well, **even if e is not adjacent to
p**. That last clause is the whole of it.

A clique that lies wholly inside the new one is dead whether or not the pivot ever touched it, so
it could never have appeared in `I[pivot]` and natural absorption cannot reach it. The paper's
worked example is a 4-by-4 where element 2 absorbs element 1 although a12 is zero.

**Amd gets the test for nothing, which is not the same as being the only one able to do it.**
`|Le \ Lp|` is exactly what Algorithm 2 computes for every clique the new one touched, so for amd
the test is a comparison against zero on a quantity already in hand. An exact-degree algorithm that
walks a clique's members can also see that all of them are already in the new clique; it just has
to be looking. **MA27 does aggressive absorption with true degrees**, so the mechanism is not tied
to the bound. genmmd simply does not do it.

**The paper's stated purpose is the BOUND, not space.** It "improves the degree bounds by reducing
|E|", and the results confirm it: aggressive absorption tends to give slightly lower fill-in, since
reducing |E| improves the accuracy of the bound. Space is a side effect, mentioned in a parenthesis.
Our own code makes the same point from the other end: it pays twice over, shortening the lists the
bound walks and the lists a later scan walks.

**Its frequency is wildly matrix-dependent and the paper says so plainly**: for many matrices it
rarely occurs, but in some cases up to half of the elements are aggressively absorbed.

**Our grid ladders are at the very bottom of that range.** `amd2` fires aggressive absorption ONCE
across the whole test set, against the hash's 2488. Grids do not produce contained cliques, so
nothing in this directory can see the rule work. That is the same shape as the dense-row finding:
a mechanism invisible to every synthetic ladder here and material on real matrices.

### MA27 is a third point in the design space, and it separates the families

The paper compares against two other codes, and the second is the useful one for seeing which
mechanisms travel together. Page 25: MA27 uses the true degree and the same data structures as AMD,
detects supervariables whenever two variables adjacent to the current pivot have the same structure
in the quotient graph, uses the true degree AS the hash function, and does aggressive absorption.
Neither AMD nor MA27 uses multiple elimination or incomplete update.

| | degree | pairwise merging | aggressive absorption | multiple elim | incomplete update |
|---|---|---|---|---|---|
| MMD | exact | q2h, no hash | no | yes | yes |
| MA27 | exact | hash, keyed on the degree | yes | no | no |
| AMD | bound | hash | yes | no | no |

Reading it by family: all three differ in family 3, MMD alone has both of Liu's techniques and AMD
alone has the bound; MA27 and AMD agree in families 1 and 2 despite disagreeing in 3.

**So the axes are independent, and MA27 is the proof.** Hash detection and aggressive absorption
work with exact degrees; they are not consequences of approximating. And multiple elimination and
incomplete update are not consequences of exactness either, since MA27 has exact degrees and
neither technique. MMD and AMD differ on five things at once, which is what makes a two-code
comparison hard to attribute; MA27 sits between them and holds three of the five fixed.

**MA27's hash is the true degree itself**, which is a nice economy: a quantity it already maintains,
used as the filter, with the exact comparison behind it as in AMD. AMD needs a separate key because
its degree is a bound and two vertices with equal bounds need not be alike.

### What each branch is missing, and whether it could have it

Two questions rather than one, and they have different answers.

#### Amd: neither of Liu's two techniques

Neither of Liu's two techniques is in AMD, and they are not equally out of reach.

**Multiple elimination looks available.** AMD forms `C[p]` in full when it eliminates, so after the
first pivot it knows exactly which vertices were reached; any minimum-bound vertex outside that set
has an unchanged bound and can be eliminated in the same round. Nothing about the bound blocks it.
Whether it would pay is a different question, since AMD's refresh is already cheap and the batch's
value is amortising an expensive one.

**Incomplete update does not look available**, for the reason two subsections above: outmatching
needs a containment between two reach sets and amd has bounds rather than reaches. No filter
rescues it either, since a hash finds candidates for EQUALITY and outmatching is a CONTAINMENT
test, and containment is not decidable from a key.

**IT NEEDS THE REACH SETS, NOT THE EXACT DEGREES, and the difference is worth keeping straight.**
The two travel together because forming the reach set is the expensive part of an exact update and
the degree is a count taken on the way past, so anything in a position to outmatch is already
paying for exactness. But the implication runs one way only:

- **Reach sets without outmatching happens**, and MA27 is the example sitting in the table above:
  true degrees, so the sets are there, and it still does not do incomplete update.
- **Outmatching without reach sets does not**, and that is amd's position. The bound is computed by
  decomposition precisely so that no clique is ever opened, so no reach set is ever formed, so
  there is nothing to compare.

So "incomplete update needs exact degrees" is the wrong statement of it. What it needs is the sets,
and exactness is what you get for free once you have them.

**WHICH MAKES IT SELF-DEFEATING FOR AMD RATHER THAN MERELY UNAVAILABLE.** The reach set is the
expensive object, and not building it is the entire saving: the bound decomposes |reach(u)| into
|A[u] \ Lp| plus a sum of |Le \ Lp| taken by subtraction from a maintained clique degree, so no
clique is opened and no set is materialized. Adding outmatching means building the very thing the
bound exists to avoid, after which the bound has nothing left to buy. That road arrives back at
MMD by a longer route.

**And it is what separates the two techniques cleanly.** Multiple elimination needs the SET OF
REACHED VERTICES, which is `C[p]` and which amd forms anyway. Incomplete update needs THEIR REACH
SETS, which amd forms for nobody. One is portable and the other is not, and the difference is one
word.

That also reads the paper's sentence correctly. "Neither AMD nor MA27 take advantage of multiple
elimination or incomplete update" is one statement covering two different situations: for AMD it is
structural, and for MA27 it is a choice, MA27 having true degrees and therefore the sets and simply
not using them that way.

#### Mmd: aggressive absorption, and it would be a cost play rather than a quality one

**It is available**, and MA27 is the existence proof: true degrees and aggressive absorption at
once. Nothing about exactness forbids it, and any algorithm that walks a clique's members can see
that all of them lie in the new clique.

**But the reason to want it is not amd's reason.** For amd the argument is quality: a shorter
incidence list makes the bound tighter, which is the paper's stated purpose and shows up as
slightly lower fill. That mechanism does not exist for mmd. An exact degree is exact whatever
`I[u]` holds, since the cliques absorption removes have all their members in `C[p]` already and so
change no reach set. **Mmd's degrees would be identical, its pivots would be identical, and its
fill would be identical.** What it would buy is shorter lists to walk and fewer entries to store.

**That makes the prediction sharp, which is the good part.** Less work and a byte-identical
permutation, which `make digest` checks directly, so a wrong implementation announces itself at
once. One caveat: removing cliques from `I[u]` changes the order members are visited in during the
refresh, hence the order vertices enter the buckets, hence which of several equal-degree vertices
is picked. Identical degrees do not quite give an identical permutation.

**Two complications specific to mmd, and neither is fatal.**

- **Outmatching means incomplete information.** The refresh skips outmatched vertices, so a clique
  reachable only through one of them is never walked that round and its containment is never
  tested. Amd walks every vertex of `C[p]` every time and has no such hole.
- **A batch has several new cliques, not one.** Aggressive absorption is defined against `C[p]`;
  under multiple elimination there are several, so the test becomes containment in any of them, or
  in their union, and that is a decision rather than a transcription.

As with hashing, it could not be `Mmd3`, which reproduces genmmd exactly. It would be a new layer
with its own twin.

Both halves of this subsection are reasoning rather than measurement, and are recorded as
directions rather than plans.

### By driver

Grouped by family: detection first, then redundancy, then scheduling.

| | mass elim | pairwise | natural abs | aggressive | bound | batch | incomplete |
|---|---|---|---|---|---|---|---|
| `md3` | yes | no | yes | no | no | no | no |
| `mmd1` | yes | no | yes | no | no | yes | no |
| `mmd2`, `mmd3` | yes | q2h | yes | no | no | yes | qxh |
| `amd1` | yes | no | yes | no | yes | no | no |
| `amd2`, `amd3` | yes | hash | yes | yes | yes | no | no |

Every driver has mass elimination and natural absorption, which is the two families' shared floor;
everything else is a branch or a layer.

`amd1` is the layer that isolates the bound: it has the approximation and none of the mechanisms
that ride along with it, which is why it is the control the two-mechanism cost is measured against.

## Two ways a vertex disappears: mass elimination and pairwise merging

Both branches fold vertices into other vertices, in two distinct operations that are easy to
conflate because both end with one vertex carrying another's weight. They are not the same
operation, they fire at different moments, and the branches implement the second one completely
differently. "Supervariable detection in production" below has the mechanism for each; this section
is the level above it, which is what the two operations ARE and which layer has which.

### The two operations

**MASS ELIMINATION merges into the pivot being eliminated in the same breath.** A vertex whose
reach is exactly the new clique is indistinguishable from the pivot, so it can be eliminated next
at no fill, and both leave the graph together. The test is cheap and local: nothing explicit left
and no clique but the new one, `|A[u]| == 0 && I[u] == {pivot}`.

**PAIRWISE MERGING merges two vertices that are indistinguishable FROM EACH OTHER**, neither of
them the pivot, and the survivor stays LIVE and carries the other's weight onward. That is the
difference that matters downstream: `QuotientGraph::merge` therefore does not call `killClique`,
the absorbed vertex never having formed a clique, and the absorbed vertex is left where it lies at
weight zero rather than purged from every clique that names it.

### Which layer has which

| layer | mass elimination | pairwise merging |
|---|---|---|
| `md3` and up, `mmd1` | yes | no |
| `mmd2`, `mmd3` | yes | yes, by the q2h test |
| `amd1` | yes | no |
| `amd2`, `amd3` | yes | yes, by hash |

**So `mmd1` and `amd1` are the two layers with mass elimination and nothing else**, and pairwise
merging arrives with the second layer on both branches. It arrives by two different routes.

### The routes differ, and q2h is a subset of what the hash can find

`mmd` decides by the q2h test, `|A[v]| + |I[v]| - 1 == 1`, which costs nothing because it is a
by-product of a union it has to compute anyway. `amd` decides by hashing the vertices of `C[p]` and
comparing exactly within a bucket, which costs a key per vertex and a comparison per colliding
pair.

**The population each can reach is not the same.** q2h catches pairs where both vertices have
exactly two sources and both sources coincide; the hash catches any indistinguishable pair in
`C[p]`, whatever is left of either. So the hash finds strictly more, and the reason is not that MMD
is being cleverer but that it is not asking: an exact degree opens every clique, so the pairs fall
out of a walk it was making regardless, while a bound never opens one and has to go looking.

### The paper's third route: pre-compression, and where it comes from

The 1996 paper reports a third arrangement, and it is the SAME HASH FUNCTION used differently
rather than a different idea. Page 18 says Ashcraft applies it as "a preprocessing step on the
entire matrix", dropping the mod (n - 1) term and sorting in O(|V| log |V|) rather than filing into
|V| hash buckets, where the paper's own use is during the ordering and over the variables adjacent
to the current pivot only.

Three axes separate the two:

| | pre-compression | detection during ordering |
|---|---|---|
| when | once, before elimination | repeatedly, at every pivot |
| over what | the whole matrix | `C[p]` only |
| what it finds | identical rows of the ORIGINAL A | any pair, including ones cliques created |

**The motivation is structural engineering**, and it is worth knowing because it explains why the
idea exists at all: those matrices tend to have many rows of identical nonzero pattern, several
degrees of freedom sharing a node, so the initial supervariables are there in A before anything is
eliminated. Ashcraft found that detecting them up front significantly improves MMD's total ordering
time. The authors built CMMD, "compressed" MMD, to measure it.

**Pre-compression does almost nothing for AMD**, and the paper says why in one clause: AMD "finds
these supervariables when their degrees are first updated". Detection during the ordering already
covers the initial case, so compressing first buys nothing, and they report that AMD on compressed
matrices plus the cost of compressing was never faster than plain AMD.

### More detection is not better fill, and that is measured rather than argued

The intuition is that a strictly larger population found means a better ordering. Three
measurements say otherwise.

**The paper's own**: "AMD, MMD, and CMMD find orderings of about the same quality." Three schemes
with very different detection reach, comparable fill.

**Ours, on grids**: `amd2` carries both the hash and aggressive absorption, fires the hash 2488
times across the test set against aggressive absorption's 1, and still fills 7 percent worse and
orders 65 percent slower than `amd1`, which has neither.

**And the reason, which is stated elsewhere in this file and belongs here too:** a supervariable is
a COMMITMENT to eliminate its members consecutively. The commitment is fill-free by construction,
but it removes choices the picker would otherwise have had. Coarser is a smaller search space, not
a better one.

**So the value of stronger detection is in TIME rather than fill**, which is exactly what Ashcraft
reported: fewer, larger supervariables mean fewer pivots and fewer degree updates. Any future
experiment here should be posed as a timing question, with pivot count and `pC` beside it, and
should expect fill to move in both directions across a set.

### If an mmd layer with hashing is ever built

Two constraints, so that it is not attempted inside an existing driver.

**It cannot be `Mmd3`.** That driver reproduces genmmd's permutation exactly and the differential
depends on it; coarser supervariables change pivot choice. It would be a new layer with its own
twin.

**And the port is not the interesting part.** amd's hash key is accumulated inside the bound's
walks, `hval += e` and `hval += j` in scan 2, which is what makes it affordable there. mmd has no
bound walk to fuse it into, so a naive version pays the whole key accumulation as a pass of its
own, and a measurement of that would be measuring the pass rather than the idea.

## The two prepasses: which drivers have one, and why they are not the same idea

Both branches number some vertices before the main loop starts, and it is easy to read that as one
feature appearing twice. It is not. The two prepasses take different vertices, for different
reasons, and only one of them uses `number`.

**WHICH DRIVERS HAVE WHICH.**

| driver | prepass | what it takes |
|---|---|---|
| `Mmd1` | none | files everything |
| `Mmd2`, `Mmd3`, `Mmd3B`, `Mmd3C` | mmd's | the degree-1 bucket, degree 0 with it via the floor |
| `Amd1`, `Amd2` | none | files everything |
| `Amd3`, `Amd3B` | amd's | degree 0, plus dense rows set aside |

**Within each branch the code is identical.** The four mmd drivers run the same eight lines, the
same after comments are stripped; `Amd3` and `Amd3B` run the same filing loop, likewise. So the
question is only ever mmd's against amd's.

### mmd's: the degree-1 bucket, and the floor is what makes it cover degree 0 too

The loop is genmmd's, over `head[1]` before the main loop: take the successor before unfiling
invalidates it, unfile, `number`, push the pivot, count it. It is described in full in "Pass 1: the
prepass" below, at the twin level.

**The floor is the part that is easy to miss.** `mmdint` files a degree-0 vertex under degree 1,
`if (dg == 0) dg = 1`, and our drivers spell it `std::max(qg.adjacencySize(u), 1)`. So isolated and
degree-1 vertices sit in the same bucket and are numbered together, and mmd needs one rule where
amd needs two.

### amd's: degree 0 and dense rows, riding inside the filing loop

`AMD_2` numbers a degree-zero vertex where it stands during initialization, in the same ascending
pass that files everything else, and sets aside any row whose degree exceeds `max(16, 10*sqrt(n))`.
Ours does both in one loop for the same reason: it is one pass either way, and neither vertex is
ever filed.

### Three differences that are not cosmetic

**They take different vertices, and the difference is real rather than a convention.** A degree-1
vertex has a neighbor, so numbering it leaves that neighbor holding a stale degree; a degree-0
vertex has nothing to leave stale. That is why mmd's prepass is a concession bought for speed,
measured below, and amd's costs nothing at all.

**Only mmd's calls `number`, and amd's deliberately does not.** `number` exists for a vertex that
is numbered while still being NAMED by its neighbors: it marks the vertex GONE and sets
`mHasNumbered`, and that flag then puts a test in every walk for the rest of the run. A degree-zero
vertex is in nobody's adjacency, so no walk can reach it and there is nothing to mark. `AMD_2`
writes `W [i] = 0` here for the same non-reason and keeps `Nv [i] = 1`. Not filing it is the whole
of what has to happen.

**Only amd's moves the live count.** That is `nel++` in `AMD_2` and `nleft = n - nel` at the degree
bound: the bound caps on `numLeft - weight(u)`, so leaving a numbered empty row in the count would
make the cap one too large per empty row. Mmd has no such bound and needs no such adjustment.

### Both exist because the permutation differed, not the fill

Neither prepass changes what gets filled. What it changes is which permutation comes out, and in
both cases the reference's answer was the one we were failing to reproduce.

On the mmd side the prepass IS genmmd's, so having it is a condition of the differential existing
at all. On the amd side, without it a degree-zero vertex was filed at degree zero and popped from
the head, so those vertices came out LIFO where `AMD_2` gives them ascending: a pure diagonal gave
`4 3 2 1 0` against `0 1 2 3 4`. Twelve matrices in `benchmarks/matrices` are entirely of that kind
and every one of them differed for that reason alone.

**The dense-row rule is the same shape of finding and a larger one.** Grids have no vertex anywhere
near `max(16, 10*sqrt(n))`, so no digest and no scaling ladder in this directory can see the rule
at all; it took real matrices to expose it, and it accounts for most of the remaining differences
against `AMD_2` on social and power-law graphs.

## The vendored storage scheme, and what it is worth

**WHAT THIS SECTION DOES NOT COVER, 2026-08-16.** It prices our clique arena against GENMMD's
dead-segment scheme, both being candidates for the mmd branch. It says nothing about `AMD_2`'s
storage, which is a different design again: one workspace, compacted and reused, with an element
taking over the slots of the variable that formed it. The amd branch has no equivalent of `Mmd3B`,
so our arena has never been compared against that.

That distinction mattered on 2026-08-16, when the vendored AMD turned out to cost more per vertex at
power-of-two grid sides while genmmd, `MMD3` and `AMD3` did not. Storage was the first suspect and
this section was cited in support of an explanation it cannot support. See
`docs/DESIGN_DECISIONS.md` (2026-08-16, later).


> **ANSWERED 2026-08-15, AND THE ANSWER IS NO.** This section was written while the measurements
> attributed most of our 2D time difference to clique PLACEMENT. That attribution is dead twice
> over. `Mmd3B` implements genmmd's scheme in full and the time did not move; and the renumbering
> experiment the placement claim rested on was CONFOUNDED, since shuffling a grid's numbering adds
> a large cost to both routines and a large common term compresses any ratio toward 1.
>
> What the difference actually was is the number of ARRAYS a vertex's state lives in, five for
> genmmd against eleven for us, each of its arrays answering several questions at once. With every
> one of those encodings now folded into BOTH files, so that storage is the only difference left,
> **our two-arena scheme beats genmmd's single nnz(A) array on every axis**: 1.02 to 1.19x genmmd
> in 2D against `Mmd3B`'s 1.15 to 1.38x, 14.22M instructions against 16.61M, 1.66M data writes
> against 2.14M, 119331 D1 read misses against 123510. Spending nnz(L) on a second arena buys
> speed.
>
> `Mmd3B` therefore stays, as the standing equal-encoding comparison against the vendored storage.
> Everything below about the MECHANISM stands and is worth reading: the chaining, the worked 5x5
> example, the entry counts, the peak-live table and the two failed experiments. Only the
> attribution of time to placement is withdrawn. `docs/DESIGN_DECISIONS.md` (2026-08-15).

Both vendored routines keep their whole quotient graph in ONE array the size of the input pattern,
and neither allocates anything for cliques. We keep two arenas and allocate a second pattern's
worth. That difference was measured in August 2026 and is the largest single thing between our
time and theirs on 2D grids, so it is written up here rather than in a layer section: it belongs
to neither branch and improves both.

**They agree on the arena and differ on two axes**, and the second one is easy to miss:

```
             placement                        density                  cost per read
genmmd   the pivot's own segment, ALWAYS   reuse, exactly nnz(A)   a sign test on every entry
AMD_2    the pivot's own space when the    compaction, 1.2 nnz(A)  none
         pivot is in no element yet;
         APPENDED otherwise
```

Both keep one arena and both stay nnz(A)-scale, so on storage they are close. On PLACEMENT they
are not. `AMD_2` has two branches, `if (elenme == 0)` builds the new element in place at `Pe[me]`
and the else builds it at `pfree`; only the first is vertex-id placement, and it applies to a pivot
that belongs to no element yet, which is an early-run case. Everything after that is appended, in
creation order, exactly as ours is. Its compaction preserves address order, so it restores DENSITY
and never restores placement.

That matters because the measurements below attribute most of the time difference to placement,
not to density. Both schemes beat what we do now; only genmmd's targets the part that was measured
to cost time. The sections describe genmmd's first, then AMD's, then what each would cost to port.

**Nothing about either is Fortran.** The negative entries and the zero terminator are 1985
packaging. Underneath is an allocation argument that would be worth making in any language, and
the sections below separate the two deliberately.

### One array, and why cliques are free

`genmmd` allocates `adjncy` at nnz(A) and nothing else structural. A **segment** is the space
between `xadj[v]` and `xadj[v+1]-1`, one per vertex, fixed at construction, and it is in one of
two states: while `v` is live it holds `A[v]` then `I[v]`; once `v` is eliminated it holds clique
members instead. So a clique is stored in its own pivot's dead segment, and is IDENTIFIED BY THAT
PIVOT: `xadj[p]` and `xadj[p+1]` delimit `C[p]` exactly as they delimited `p`'s own list. There is
no clique arena, no clique offset array, and no clique id space.

The transition happens inside `mmdelm`, which is why the compaction of the pivot's list and the
construction of its clique are the same pass. The write cursor starts at `xadj[md]` and overwrites
the pivot's own entries as it goes; it is safe because the cursor never overtakes the read index,
every entry written having been read first.

**Two conservation arguments make it fit, and they are the whole scheme.**

For a vertex: an elimination that reaches `u` removes at least one entry from `u`'s run and adds
exactly one. Either the pivot was a direct neighbour, so it leaves `A[u]` and its id enters
`I[u]`; or `u` was reached through a clique `c`, and then `p` is a member of `c`, so every member
of `c` is reachable from `p`, so `c` is absorbed and leaves `I[u]`. `|A[u]| + |I[u]|` never grows.
Ours conserves the same way, and this half is not the difference between us.

For a clique: `C[p]` is contained in the pivot's live neighbours together with the members of the
cliques it absorbs, and every one of those contributed a whole segment. Capacity is a sum over a
multiset, content is its union, so capacity is at least content BY CONSTRUCTION. No bound is
checked and no allocation can fail.

### The chain, and how the next segment is chosen

A clique frequently does not fit in its own segment: measured on grids, 65 to 69 per cent do, so
about a third overflow. `mmdelm` handles that by continuing the list in another segment and
linking the two.

The first loop over the pivot's segment sorts its entries. Live neighbours are compacted in place
at the write cursor. Eliminated entries -- the cliques the pivot belongs to, all of which are
about to be absorbed -- are pushed onto a stack, `list[nb] = el; el = nb`. The second loop walks
that stack, and before absorbing each one it writes

```
adjncy[rm] = -el;          rm is the LAST entry of the segment being filled
```

so **the continuation is the next clique to be absorbed**, taken in LIFO order. Nothing is
searched for and no free list is kept: the pivot's own incidence list IS the list of available
segments, because it is exactly the set that is dying. When the cursor reaches that entry it reads
the value back and jumps:

```
while (rl >= rm) { lk = -adjncy[rm]; rl = xadj[lk]; rm = xadj[lk+1]-1; }
```

The link costs one entry of a segment that is being superseded anyway. A reader then has three
cases per entry, and every walk in the file has this shape:

```
positive   a member
negative   a link; jump to xadj[-value] and continue
zero       the clique ends
```

**A real example, from a 5 by 5 grid, 1-based as genmmd is.** Vertex 4 is on the boundary with
three neighbours, so its segment is entries 9 to 11. Its neighbours 3, 5 and 9 have already been
eliminated, so by this point the segment holds the three cliques it belongs to:

```
entry     9   10   11
value     5    3    9        clique ids, not vertices
```

Eliminating 4: there are no live neighbours to compact, so nothing is written yet. Clique 9 is
absorbed first, `-9` goes into entry 11, and vertex 9's segment -- entries 25 to 28, four entries
since 9 is interior -- becomes the continuation. Two live members of `C[9]` are written at entries
9 and 10; the cursor reaches 11, follows the link, and continues at 25. The result:

```
[C4] segment of vertex 4, entries  9..11:   8  10  -9
[C4] segment of vertex 9, entries 25..28:  14   2   0   <- zero terminates
```

`C[4] = {8, 10, 14, 2}`. Note entry 28 holds a stale link written by a later absorption that the
walk never reaches, because the zero at entry 27 always comes first on a correct walk.

**Chains can exceed two segments.** Later in the same run pivot 10 fills its own segment, follows
`-24` into vertex 24's, fills that, and follows `-4` into vertex 4's. That is why the negative test
sits in the inner loop rather than at a boundary check.

**And the chain is almost never followed.** Measured hops as a fraction of entries read: 2.79 per
cent at 100 a side, 1.59 at 200, 0.84 at 400 -- one hop per thirty clique walks at the largest
size, and the rate FALLS with n. A clique overflows only when the reach exceeds the pivot's
original degree; later in a run pivots are supervariables that have absorbed many vertices and
inherit long chains of absorbed segments, so capacity outruns content. So the branch is
predicted-not-taken essentially always. The machinery is not a price paid for the layout; it is a
price that turns out not to be charged.

### AMD's answer: one pool, compacted in place

`AMD_2` reaches the same place without a single link. Its arena is `Iw`, sized

```
slen  = nzaat + nzaat/5 + 7n
iwlen = slen - 6n  =  1.2 nnz(A) + n
```

so one pool with 20 per cent slack. `Pe[i]` is the offset of object `i` and `Len[i]` its length,
for rows and elements alike -- the same identified-by-its-pivot trick, since an element takes over
the supervariable's entry in `Pe`. Note that `Len` is exactly our `mCliqueSize`: AMD keeps a length
array and does NOT use a terminator.

**Placement, and the branch that is easy to miss.** `AMD_2` builds a new element in one of two
places:

```
if (elenme == 0)   pme1 = Pe[me];    in the pivot's OWN space
else               pme1 = pfree;     appended at the tail
```

**THE BRANCH IS A TEST, AND IT IS NOT A FIT TEST.** It asks whether the PIVOT belongs to any
element, and nothing is measured against anything: no length is compared, and the clique's size is
not known when the choice is made. It behaves like a fit test only because of what the condition
implies. With no elements, `reach(me)` is drawn from `A[me]` alone and so cannot be larger than it,
which makes the fit a guarantee rather than a check. **A real fit test would keep more cliques in
place**, since a pivot WITH elements often has a reach that would fit in its own segment and is
sent to the tail regardless. `AMD_2` is deliberately the more conservative of the two, and `Amd3B`
copies the code rather than the behavior it approximates.

The first case is vertex-id placement, and it applies when the pivot belongs to no element yet,
which is an early-run condition: measured 62 to 68 per cent of eliminations on grids. Everything
after is appended in creation order, the same as ours. So AMD does not have genmmd's placement
property, and its compaction does not restore it: the sweep copies live blocks toward the front
preserving their address order, which is creation order.

When `pfree` reaches `iwlen`, it garbage-collects rather than chaining:

```
for (j = 0 ; j < n ; j++)                 mark every live object by FLIPping its first
    if (Pe[j] >= 0) { Pe[j] = Iw[pn] ; Iw[pn] = FLIP(j) ; }

while (psrc <= pend)                      sweep, copying live blocks toward the front
    j = FLIP (Iw[psrc++]) ;
    if (j >= 0) { Iw[pdst] = Pe[j] ; Pe[j] = pdst++ ; ... copy Len[j] entries ... }
```

One O(n + nnz) pass: mark, sweep, rewrite `Pe`. Dead blocks vanish, live ones come out contiguous
and still in ascending address order, and the partially built element is moved last so the
elimination can continue where it left off. `Info[AMD_NCMPA]` counts them, and the figure recorded
elsewhere in this tree is ONE for a whole 140 by 140 run: the 20 per cent slack absorbs nearly
everything, so the sweep is a rare event rather than a per-step cost.

**SO THE ORIGINAL SEGMENTS DO MOVE, and this is the point of the pass rather than a side effect.**
Between collections nothing shifts: a run is rewritten in place by the prune and only ever shrinks.
At a collection every live block slides down, `A[u]` and `I[u]` runs included, and every `Pe` is
rewritten to where its block landed. genmmd never does this to a segment at any time, which is one
of the two real differences between the schemes.

**AND THAT IS THE ONLY WAY THE FREE AREA IS REFILLED.** The space is already free before the sweep
runs, in the gaps left by shortened runs and abandoned blocks; it is unusable because a clique must
be contiguous and the gaps are scattered. The sweep does not create space, it GATHERS space that
already exists into one usable run at the top. Which is why the scheme works with no elbow room at
all: the pool starts full, the first pivot has no elements and so builds in place, and by the time
anything needs the tail there is dead space to gather. Confirmed by intervention rather than by
argument, `src/Amd3B.cpp` sized to exactly `nzaat`: 3 to 14 collections from 3 to 200 a side, every
permutation byte identical to the run with the 20 per cent slack.

**The two schemes trade the same way as a linked list against a vector.** genmmd never copies and
never over-allocates, and pays a branch on every entry read forever. AMD reads a flat block with no
test at all, and pays a full sweep on the rare occasion it fills. The measured hop rate on genmmd's
side, 0.84 to 2.8 per cent, and the compaction count on AMD's, one per run, say both prices are
small; they are small in different places. What they do NOT trade evenly is placement, which only
genmmd keeps, and which is the axis the timings below turn on.

### What it is worth, measured

Same graph, same permutation, `orderMmd3` against `mmd_order`, on square grids.

**We do LESS work and take MORE time.** Arena entries touched over the whole ordering:

```
                vendored elim   refresh    TOTAL      ours QG   refresh    TOTAL    ratio
64x64                  126485    117128    243613       89862     89962    179824   0.74x
140x140                618587    526131   1144718      406662    408269    814931   0.71x
200x200               1271621   1054314   2325935      818496    820571   1639067   0.70x
```

Thirty per cent fewer entries, about fifty per cent more time: roughly 2.1x cost per entry. That
closes the question of whether the gap is work or memory.

**The mechanism is placement.** Their offsets are vertex ids, so the cliques a vertex names --
pivots that were near it in the graph -- land near it in `adjncy` under natural numbering. Ours
are append positions in elimination order, and minimum degree eliminates a vertex's neighbours at
unrelated times. Mean gap between consecutive clique blocks read:

```
                        ours              genmmd         ratio
2D 100          1001 entries (63 lines)   438  (27)      2.3x
2D 200          2073 entries (130 lines)  844  (53)      2.5x
2D 400          4197 entries (262 lines) 1648 (103)      2.5x
```

**The decisive test is renumbering.** If the advantage is placement inheriting the caller's
locality, shuffling the input should destroy it. It does:

```
                        vendored    Mmd3     ratio
2D 200, natural            6.15     9.69     1.57x
2D 200, RANDOM            31.15    36.03     1.16x
```

Stable over repeats, permutations identical throughout. Most of their advantage is the input
numbering, not their code. **1.16x is the honest size of everything else** -- two arenas, index
widths, extra arrays, extra passes -- and 0.41 of the ratio is placement.

**Which is also why 2D is behind and 3D is not.** The placement penalty is per BLOCK; 3D reads far
more members per block, so the same penalty is diluted:

```
              members per clique lookup      locality's share of the ratio
2D 200                  11.67                      0.36   (1.50 -> 1.14)
3D 34                   33.88                      0.05   (0.97 -> 0.92)
```

Two independent quantities moving together in both families, including the crossing where MMD3
overtakes genmmd on cubes.

### What we do instead, and what it costs

`mSource` conserves as `adjncy` does: nnz(A) reserved once, nnz(A) - n used, never grown, offsets
written at construction and never moved. That half is settled.

`mCliqueArena` is the difference. It is append-only in elimination order and NOTHING IS EVER
RECLAIMED -- `mCliqueSize[absorbed] = 0` with the comment "dead, its block left behind". So it
holds every clique ever formed rather than the live ones:

```
            nnz(A)     cumulative        PEAK LIVE        held but dead
2D 200      199200    186730 (0.94x)    78812 (0.396x)      2.4x
2D 400      798400    752920 (0.94x)   317612 (0.398x)      2.4x
3D 26       118976    163671 (1.38x)    48573 (0.408x)      3.4x
3D 32       223232    322820 (1.45x)    92043 (0.412x)      3.5x
```

Live clique storage is 0.40 nnz(A), flat across both families and a sixteenfold range in n. We
hold 2.4 to 3.5 times that. 2D landing just under nnz(A) is a coincidence; 3D crosses because
bigger cliques abandon more per step.

### Two experiments that failed, and why they are worth recording

**One arena, cliques appended to `mSource`.** Removes an allocation and a stream, keeps append
order. Paired interleaved timing, which is the only method that resolves effects this small:
0.867 of base at 140, 0.955 at 200, 0.961 at 280. Real, consistently signed, and far short of the
gap. It does not move the placement, which is what the gap turned out to be about.

**Vertex-ordered arena, blocks at the pivot's own offset, overflow appended.** This was meant to
be the cheap way to buy the placement without chaining. It made things WORSE: 1.022 at 200, 1.044
at 400. The measurement says why -- the mean gap went from 2073 entries to 4980, more than double,
in the direction of the append arena it replaced.

**The lesson is that vertex ORDER is not vertex-order DENSITY.** `adjncy` is dense because a dying
vertex hands over its WHOLE segment, so consecutive ids are consecutive data. An arena of n
segments in which only the fitting two thirds hold anything, each block smaller than its segment,
is correctly ordered and further apart than appending. **Reclaiming is not a second, optional step
after placement; it is the precondition that makes placement worth anything.**

### If this is ported

Three notes for whoever builds it.

**The zero terminator does not survive the move to 0-based.** genmmd can use 0 because vertex ids
start at 1; ours start at 0. We already keep `mCliqueSize`, and with chaining a size and the links
do different jobs: the size says when to stop counting members, the links say where the members
continue. Walk while counted is below size, treating a negative as a jump that does not count. No
terminator needed and no second sentinel introduced.

**The sign bit earns a second use.** Indices are signed `int32_t` in this tree because `NIL` has to
share a type with the values it stands in for. A link is a second reason for the same choice, and
the two do not collide: `NIL` never appears in an arena and a link never appears in a field that
could hold `NIL`.

**It is shared-class work, so it lands on both branches.** `QuotientGraph` serves `Mmd1` through
`Mmd3` and `Amd1` through `Amd3`, and the placement result was measured on mmd only. The amd
branch's own gap is 1.60 to 2.03 in 2D, larger than mmd's, and nothing has tested whether the same
cause is behind it.

**Either vendored scheme beats what we do now, so the choice is which to build first and not
whether to.** Today's arena is appended AND never reclaimed, so it loses to genmmd on placement and
to AMD on density, and to both on total storage. That is the state to leave, and the comparison
between the two vendored answers is empirical: each has to be built and measured, and the second
one is then measured against the first rather than against today.

**On effort, AMD's is the smaller change.** It appends as we already do, keeps blocks contiguous as
we already do, and its `Len` is our `mCliqueSize`. What it adds is one mark-and-sweep pass and a
rewrite of `mCliquePtr`, in one place, called rarely -- `Info[AMD_NCMPA]` records a single
compaction for a whole 140 by 140 run.

**On effect, genmmd's is the one that targets what was measured.** Placement is worth 0.41 of the
time ratio in 2D and vanishes under random renumbering; density is worth the 2.4 to 3.5x of dead
blocks and no measured time. AMD's compaction fixes the second and leaves the first alone, because
its elements are appended in creation order except for the early `elenme == 0` case.

**The cost of chaining is call sites, not novelty.** The scheme is fully described above and
porting it is translation, which is what this tree does. What it touches is the reader: a clique
walk is a flat pointer-and-length loop at about ten sites today, and each becomes a walk that can
resume in another segment, with a sign test in the inner loop. The terminator question does not
arise -- keep `mCliqueSize`, walk while the count is below it, and treat a negative as a jump that
does not count. The sign bit is free because indices are already signed for `NIL`, and the two
cannot collide: `NIL` never appears in an arena and a link never appears where `NIL` could.

## md1: the elimination game


Minimum degree is a greedy heuristic on one observation: eliminating a vertex makes its
neighbors pairwise adjacent, so the cheapest vertex to eliminate now is the one with fewest
neighbors. md1 does exactly that and nothing else. Pick the live vertex of least degree, add
every missing edge among its neighbors, delete it, repeat.

The two things worth taking from md1 are what it stores and what it counts. It stores the
elimination graph explicitly, so the fill edges are materialized as they are created, which
is why md1_eliminate can return them and why the layer can print each one. And it counts
nnz(L) for free: at the moment a vertex is chosen its adjacency is exactly the off-diagonal
pattern of its column of L, so the running degree sum plus n is the answer.

What it does badly is storage. Every fill edge is a real edge in a real set, so the graph
grows as the run proceeds, and on a large problem the fill dominates the original structure
by orders of magnitude. That is the problem md2 exists to solve.

### md1 in set operations

There is only one place set algebra appears at this layer, and it is md1_eliminate. The pick is
`degree(u) = |A[u]|`, a lookup, and there is nothing to say about it. The elimination is four
lines:

```
for u in A[p]:
    fill(u) = A[p] - A[u] - {u}                what was not already there
    A[u]    = ( A[u] | fill(u) ) - {p}
A[p] = {}
```

`fill(u)` is separated out rather than folded into the union because md1 reports it: the
difference is the set of new edges, and printing it is what makes the layer legible. The union
could be written `A[u] = ( A[u] | A[p] ) - {u, p}` in one line, and that is the same thing, but
then the fill would have to be recovered afterwards by a second difference.

**A[p] appears on the right-hand side of a line that assigns to A[u], and u is in A[p].** That is
the one hazard in the block. If A[p] were being modified inside the loop the later iterations
would see a different set from the earlier ones, so the pivot's adjacency is captured once, before
the loop, and cleared once, after it. Every layer from md2 on does the same thing with the
reachable set for the same reason.

The code computes the difference without building a set: stamp A[u] with the mark array, then walk
A[p] and keep whatever is unstamped. Two passes over lists, one stamp and one comparison per
element, rather than a hash per element. The set view is in the docstring, and md5's header states
the convention that every layer follows.

**What the four lines cost is the whole argument for md2.** `fill(u)` is materialized, so A[u]
grows, so the next elimination's difference is taken against a larger set, and the growth
compounds. The set operations themselves never change from here: md2's eliminator does a union and
three differences too. What changes is what they are taken over.

## md2: the quotient graph

The insight is that the fill edges created by one elimination are not independent: they form
a clique, and a clique of d vertices can be stored as a d-element list instead of d(d-1)/2
edges. md2 stores it that way and never materializes fill at all.

The representation splits the neighbor relation in two. A[u] holds the vertices u is still
explicitly adjacent to; I[u] holds the cliques that contain u; C[c] holds the members of clique c. C
carries the same double duty A does, bare C being the set of live cliques and C[c] the vertex set of
one of them, so a clique and its name are the same thing here. The true neighborhood is the union,
computed on demand by md2_neighbors, and George and Liu's reachable-set theorem is the guarantee
that this union equals what md1's adjacency would have been. Same degrees, same pivots, same order,
on every graph.

Two mechanisms keep it from growing. Pruning: when a new clique is formed, any explicit edge
with both ends inside it is deleted, since the clique now implies it. Absorption: the cliques
the pivot itself belonged to are entirely contained in the new one, so they are deleted
outright. The result is that A[u] only ever shrinks and the clique count stays small, so
total storage falls monotonically where md1's rises.

The price is that the degree is no longer a lookup. Every query unions the explicit adjacency
with every clique the vertex belongs to, and the pivot search does that for every live vertex
at every iteration. Making the degree cheap again is what md4 and md5 are for.

### md2 with the exact degree, which is what it does

Every layer from md2 on faces the same choice, since from here a degree is a union rather than a
lookup, and a union can be counted or estimated. Written out in set operations, one iteration of md2
as it stands is:

```
the pick

for every live u:
    reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}
    degree(u) = |reach(u)|
p = the u minimizing degree(u)

the elimination

C[p] = reach(p)                   absorb into C[p]
C    = C - I[p]                   reclaim I[p]
for u in C[p]:
    A[u] = A[u] - C[p] - {p}      prune
    I[u] = ( I[u] - I[p] ) | {p}  absorb into C[p] and reclaim I[p]
```

The union on the second line is the expensive object. It cannot be counted by adding sizes,
because the sources overlap, so it has to be deduplicated, which is the mark pass. The cost per
vertex is |A[u]| + sum |C[c]| over c in I[u].

**Absorption is the first line, not the second.** The union in `reach(p)` pulls in every C[c] with
c in I[p], so each of those cliques is a subset of C[p] the moment C[p] exists. That IS the
absorption, and it has already happened when the second line runs. The second line only reclaims
what has become redundant.

**The last line is the first two written on the I side**, which is why its comment repeats them
verbatim rather than describing it again. The incidence relation is stored twice, once in each
direction: u is in C[c] exactly when c is in I[u]. So `- I[p]` is the same subtraction as
`C = C - I[p]`, on the per-vertex side, and `| {p}` is the same statement as `C[p] = reach(p)`,
read the other way. Four lines, two facts.

That pairing also settles what is optional. The two reclamations go together, `C = C - I[p]` and
the `- I[p]` beside it: drop both and the algorithm is still correct, only wasteful, since an
absorbed clique still sitting in I[u] contributes a subset of C[p] to the union and changes no
degree. Drop only one and `reach` walks a clique that is no longer in C. The `| {p}` half is not
optional at all, since without it the new clique reaches nobody.

**Pruning is the same idea on the other half of the representation.** An explicit edge with both
ends inside C[p] is implied by C[p], so dropping it changes no reachable set. Absorption reclaims
redundant cliques and pruning reclaims redundant edges, and in both cases the line that appears to
do the work is only freeing what the union already made redundant.

**The new clique is written C[p] and not given a name of its own**, which is a small choice with
two consequences worth pointing at. The order of the first two elimination lines is now forced,
and it is the right order: reach(p) reads the cliques in I[p], and the next line deletes them, so
the new clique is built out of the old ones before they go. And `C[p] = reach(p)` turns out to be
the whole of what an elimination IS at this layer. The pivot stops being a vertex with a reachable
set and becomes a clique holding that same set; everything under the loop is bookkeeping, stripping
what the new clique now implies from the explicit edges and from the incidence lists.

Both readings rest on C carrying the same double duty A does: bare C is the set of live cliques,
subscripted C[c] is the vertex set of one of them, exactly as bare A is the matrix and A[u] is one
adjacency. So `C - I[p]` is a difference of two sets of cliques and not an analogy, and the new
clique sits inside the naming scheme the rest of the file already uses rather than beside it.

That identity is exact here and only here. From md3 on, mass elimination trims the merged vertices
out in the same call, so the line becomes `C[p] = reach(p) - merged` and the loop runs over the
trimmed set.

### md2 with an approximate degree, which it can define but not use

The pair above is what md2 computes:

```
reach(u)  = ( A[u] | C[c] for c in I[u] ) - {u}
degree(u) = |reach(u)|
```

and the second line is the expensive one, because the size of a union is not the sum of the sizes
of its parts. Deduplication is what costs. So: decompose reach(u) into pieces, add their sizes
anyway, and accept whatever error the overlaps cause. That is the entire idea, and the rest is
choosing the decomposition well.

**One invariant first**, because it says where an overlap can and cannot be. When clique c is
formed, that is when the vertex c is eliminated, every u in C[c] has
A[u] = A[u] - C[c] (equivalently, A[u] = A[u] - C[p] - {p}, which is how the prune above writes
it, since at that moment c IS the pivot p), and neither set grows afterwards. So

```
A[u] & C[c] = {}     for every c in I[u]
```

The explicit part never overlaps a clique. All overlap is clique against clique.

**The decomposition, anchored at the new clique**, where p is the pivot just eliminated:

```
reach(u)  = ( A[u] - C[p] ) | ( C[p] - {u} ) | ( C[c] - C[p]  for c in I[u] - {p} )
bound(u)  = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]|
```

**Where the first line comes from**, since it is not obvious. The precondition is that u is a
member of C[p], so p is in I[u], and that is what lets C[p] be pulled out as a term of its own.
Four iterations, and only the third does any work.

Start:

```
reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}
```

**1. Split off the new clique.** p is in I[u], so C[p] is one of the C[c]:

```
reach(u) = ( A[u] | C[p] | C[c] for c in I[u] - {p} ) - {u}
```

**2. Distribute the removal of u.** Nothing subtle:

```
reach(u) = ( A[u] - {u} ) | ( C[p] - {u} ) | ( C[c] - {u}  for c in I[u] - {p} )
```

and A[u] - {u} is A[u], since a vertex is never its own explicit neighbor. The middle term is
already in its final form.

**3. Replace C[c] - {u} by C[c] - C[p], which is the whole trick.** Removing more from each old
clique, and the union is unchanged. Both directions:

```
C[c] - C[p]  is a subset of  C[c] - {u}     because u is in C[p], so nothing new is lost
```

and for the other direction, take any x in C[c] - {u}:

```
x in C[p]        ->  x is in C[p] - {u}, which is ALREADY the second term
x not in C[p]    ->  x is in C[c] - C[p], the third term
```

Either way x is still in the union, so nothing is lost. **The whole of C[c] & C[p] can be
subtracted from every old clique, because the middle term already covers it**, and that is the
entire reason the decomposition is anchored at C[p] rather than anywhere else.

**4. Replace A[u] by A[u] - C[p], and this one changes nothing.** By the invariant above,
A[u] & C[p] = {} already, since p is in I[u]. It is written as a subtraction only to make the three
terms VISIBLY disjoint, which is what the next line needs.

**Result.**

```
reach(u) = ( A[u] - C[p] ) | ( C[p] - {u} ) | ( C[c] - C[p]  for c in I[u] - {p} )
```

That is the first line, and it is an identity, exact, for any u that C[p] reached. The second
line is the first with the union replaced by a sum, which is one substitution and the only one, so

```
degree(u) <= bound(u)
```

with equality exactly when the pieces happen to be disjoint. The first two terms always are: the
first by the invariant, the second by construction. So the entire error is the overlap between two
old cliques outside C[p], which is the smallest place the error could have been put, since C[p] is
the largest region they share and it has been factored out.

bound(u) IS the approximate degree, and AMD is what comes of putting it where degree(u) stood: pick
the u minimizing bound(u) rather than the u minimizing degree(u), and change nothing else. The
elimination is untouched, the representation is untouched, and the algorithm is the same algorithm
reading a different number.

Two caps ride along, and they are a different kind of thing, since they bound the same set from
outside rather than measuring a decomposition of it:

```
bound(u) = min( n - k - weight(u),                     nothing exceeds what remains
                degree_old[u] + |C[p] - {u}|,          it can only grow by the new clique
                |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]| )
```

k is the number of original vertices eliminated so far, so n - k is what is still live, and the
weight(u) comes off because reach(u) excludes u and u stands for a whole supervariable. The AMD
paper writes this cap as n - k, which is the md2 reading where every weight is 1; from md3 on the
subtraction is needed and the code has it, as `num_left - len(super_members[u])`.

That k is a count of ORIGINAL vertices while the code increments by `1 + len(merged_vertices)`,
a count of representatives, looks like a units mismatch and is not. An original leaves the live set
exactly once: it is counted at the iteration it is merged, or at the iteration it is chosen as
pivot, and if it later sits inside some other supervariable's member list it is not counted again.
So
`num_eliminated` really does count originals, and the loop condition `num_eliminated < n` compares
like with like.

**Why it is cheap** is one property and worth stating separately from why it is correct:
|C[c] - C[p]| depends on the clique c alone and not on the vertex asking, so it is computed once
per clique per iteration and then read once per vertex. The exact degree recomputes a union per
vertex.

**And here is what md2 cannot do with any of it.** Look again at the first line: it holds for u in
C[p], because it puts C[p] - {u} inside reach(u). For a vertex the iteration never reached, the
decomposition is not merely loose, it is wrong. So bound(u) covers exactly the vertices the new
clique reached, and nothing else.

md2 recomputes every live vertex at every iteration. Most of what it updates lies outside C[p] and
would have to be counted exactly regardless, so the estimate would apply to a shrinking minority of
the work and save nothing worth having. The layer that narrows the degree update set to exactly C[p]
is md4, and C[p] is exactly the bound's domain. bound(u) is inherently incremental; md4 is where the
algorithm becomes incremental; amd sits at the top of the chain because those two facts meet there
and not earlier.

So md2 is where the estimate can be DEFINED, since it is where the union it decomposes first
exists, and md4 is where it can be USED. That gap is the reason this subsection is a derivation and
the md4 one is a decision.

### Fill per elimination, which md1 reports and md2 gives up

md1_eliminate returns its fill edges and the trace prints them. md2_eliminate returns pruned edges
instead, and the two are not the same quantity. The substitution is worth explaining, because it
is the first thing the quotient graph costs that is not storage.

**The fill at u has the same shape as md1's, taken against the reachable set:**

```
fill(u) = ( C[p] - {u} ) - reach_before(u)
```

and reach_before(u) is the union md2 exists to avoid computing. So reporting fill would mean one
full reachable-set query per member of C[p], at the moment of elimination, which is exactly the
operation the layer defers. md1 could report it because its neighborhood was already materialized:
the difference was a byproduct of a set it had lying around.

**Every pair inside C[p] falls into three classes at md2, where md1 had two.** Explicit before,
which is pruned; implicit before through some shared clique; and new, which is fill. md1 has no
third class because it has no third representation, so for md1 "not already an edge" and "new" are
the same statement. At md2 they are not, and separating them is what costs.

**And the third class is expensive for the reason everything else at md2 is expensive.** Counting
the pairs inside C[p] that were already implicit means counting the union of the old cliques
restricted to C[p], which is inclusion-exclusion over overlapping sets. That is the same obstruction
that makes the exact degree expensive and that the bound above sidesteps by overcounting. Fill and
degree are hard at md2 for one reason, not two.

**One detail that looks like an opportunity and is not.** At md4 and md5 the degree update pass
already computes reach(u) for exactly the members of C[p], so the query seems to be sitting there
already.
It is not: the degree update pass runs after the clique is installed, so it yields the new reach,
and fill needs the old one. Getting it would mean a second query per member, before the update,
which doubles
the layer's remaining cost.

**What survives is the aggregate, and it survives intact.** At the moment a vertex is chosen, C[p]
is the off-diagonal pattern of its column of L, so |C[p]| accumulates nnz(L) at every layer for
free, exactly as md1's degree does. Total fill is nnz(L) - nnz(tril A) and is reported by every
layer in the ladder. What md2 gives up is per-iteration, per-vertex fill. The total never depended
on it.

## md3: mass elimination

Two vertices are INDISTINGUISHABLE when they have the same closed neighborhood, that is when

```python
md3_neighbors(A, I, C, u) | {u} == md3_neighbors(A, I, C, v) | {v}
```

so each sees exactly what the other sees plus the other itself. Such vertices are necessarily
adjacent, and eliminating one right after the other creates no fill whatever: the second
one's neighbors were made a clique by the first. Mass elimination is the decision to take
them together as a single pivot, called a SUPERVARIABLE, rather than rediscovering the fact
by search.

The definition is about the REACHABLE sets and not about the containers. The structural test
`A[u] | {u} == A[v] | {v}` together with `I[u] == I[v]` implies it and is far cheaper, but it
is only sufficient. Two stale singleton cliques, one holding u and one holding v, contribute
nothing to either reachable set yet make the I lists differ, and the structural test then
declines a genuine pair. graph5 shows exactly that at its fourth iteration.

md3 detects the case positionally rather than by comparison. Immediately after a clique is
formed it checks each member for `not A[u] and I[u] == {pivot}`, which says the new clique is
everything u can still reach, hence that u and the pivot were indistinguishable before the
iteration. The test is cheap, two container reads per neighbor, and CONSERVATIVE, meaning
sufficient rather than necessary. It misses supervariables, in two independent ways, and both
are a deliberate trade of coverage for cost.

**Miss one, inside the pivot's own neighborhood.** The test fails when A[u] is non-empty but
contained in the clique, or when I[u] holds another clique whose members happen to lie inside
it. Such a u is genuinely indistinguishable from the pivot and is skipped anyway. The exact
test is `md3_neighbors(A, I, C, u) <= C[pivot]`, and it costs a reachability query per
candidate, O(d) set unions per iteration where the cheap test is O(d) constant checks. Measured on
600 random graphs, the exact test merges strictly more on 415 of them and never changes the
fill.

Five vertices suffice to see it. With adjacency 0:{3,4}, 1:{2,4}, 2:{1}, 3:{0}, 4:{0,1}, md3
orders 2 1 3 0 4, five iterations and no merge at all. Iteration 3 is the one to look at: its pivot
is 0 and its clique is {4}. At that moment A[4] is empty and I[4] is {c0, c1}, the second clique
left over from eliminating 1 at iteration 1, with members {4}. So everything 4 can reach lies
inside the new clique and 4 is indistinguishable from the pivot, but `I[4] == {pivot}` is
false and the cheap test declines. The exact test merges it, making the order 2 1 3 (0 4).

**Miss two, between two members of the clique.** Even the exact test compares each candidate
only AGAINST THE PIVOT. Two members can be indistinguishable from each other while neither is
absorbable into the pivot, and no pivot-relative test finds them. An exhaustive pass over the
clique would: all pairs, O(d * d) comparisons at O(d) each, so O(d * d * d) per iteration. That is
what hashing reduces, bucketing by structure and comparing only within a bucket, which is
what amd2 does and what the section on detecting supervariables against each other covers.

Five vertices again, and this one is graph7. With adjacency 0:{1,2,4}, 1:{0,4}, 2:{0,3,4},
3:{2,4}, 4:{0,1,2,3}, md3 orders 1 0 (2 3 4), three iterations, the last of which merges 3 and 4
into 2. Iteration 1 is the one
to look at: its pivot is 0 and its clique is {2, 4}. After pruning, A[2] is {3} and A[4] is
{3}, and both have I equal to {c0}. Neither is absorbable into the pivot, since each still
reaches 3, which is outside the clique. But their closed neighborhoods are equal, 2 reaching
{3, 4} and 4 reaching {2, 3}, so 2 and 4 are indistinguishable from each other already at
this iteration and could be merged here. No test framed against the pivot will ever see it, and
the exact test does not help either: both orders are 1 0 (2 3 4). md3 does group them one
iteration later, when 2 becomes the pivot and both 3 and 4 fall to its own test, which is luck
rather than detection: a graph where 2 and 4 never share a later iteration would keep them apart
for good.

The two misses are independent. Sharpening the pivot test does nothing for the pairwise
population, and a pairwise pass does not subsume the pivot test, since by then the pivot is no
longer live to be compared against.

**A third route: stop the cheap test from being fooled.** Both misses above are about
detection, but graph5's case is not really a detection failure at all. What defeats the cheap
test there is a redundant clique: c1 holds only 4, its single member lies inside the new
clique, and it contributes nothing to what 4 can reach, yet its presence in I[4] makes
`I[4] == {pivot}` false. Ordinary absorption cannot remove it, since absorption takes only the
cliques the PIVOT belonged to and 0 was never in c1. AGGRESSIVE ABSORPTION removes any clique
whose members lie inside the new one, which deletes c1, leaves I[4] == {c0}, and lets the
unchanged cheap test fire.

That is exactly what amd2 does with it. On graph5, amd2 with aggressive absorption on takes four
iterations and reports `merged = 4, absorbed = c1`; with it off it takes five and merges nothing,
behaving like md3. The hash detection reports nothing in either run, so this case is entirely
about absorption and not about comparison. So there are three ways to recover a missed
supervariable, and they are genuinely different: sharpen the test, compare vertices pairwise,
or clean up the structure until the cheap test is no longer fooled.

Three consequences follow, each covered in its own section below. The degree becomes
WEIGHTED, since a neighbor now stands for several original vertices, and it must be EXTERNAL,
excluding the supervariable's own members. The nnz(L) count stops being a running degree sum,
because an iteration is now w consecutive columns rather than one. And the order changes, since a
merged vertex is eliminated immediately where md2 would have reached it later.

The structure this produces is related to a supernode but neither contains nor is contained
in one, and the relation is worth stating carefully because it is easy to get backwards.

A supervariable always has exactly nested patterns, since its members were indistinguishable
when they merged, so it is always a supernode in the general sense of consecutive columns
with the same structure outside the block. What it need NOT be is a FUNDAMENTAL supernode,
which additionally requires each column after the first to have exactly one child in the
elimination tree. Indistinguishability is a property of the current elimination graph and
says nothing about children acquired earlier, on branches the iteration never touched. graph5 is
the counterexample: 0 and 4 are indistinguishable at the fourth iteration, yet the tree is
2 -> 1 -> 4 and 3 -> 0 -> 4, so 4 already has 1 as a child and {0, 4} fails the single-child
test. Amalgamating it anyway costs nothing, one entry stored against one actual nonzero and
zero explicit zeros, so it is a relaxed supernode that a threshold of 0 would accept.

Measured across 800 random graphs, all 1012 supervariables of size greater than one
amalgamate at zero cost, which is what the exact nesting predicts. Whether the converse holds
is a different question and the answer is no: graph5's own {0, 4} is a zero-cost pair that
md3 does not merge at all. So supervariables are a subset of the zero-cost relaxed supernodes,
not the same set.

In the other direction a supernode can be strictly larger than any supervariable, since it
admits columns that were never indistinguishable and merely nest along the tree. Take graph6,
the six-vertex graph with adjacency 0:{2,3,4}, 1:{3}, 2:{0,3,4,5}, 3:{0,1,2,4}, 4:{0,2,3}, 5:{2},
which md3 orders as 1 5 (0 4) 2 3, five iterations with 4 merging into 0 at iteration 2, and no fill
at all. The exact test takes it further, to 1 5 (0 2 3 4), a single supervariable over all four.
Its column patterns nest all the way up, [4,2,3], [2,3], [3], [], so the matrix shows one
dense block over 0, 4, 2, 3, yet the fundamental supernodes are {0,4}, {2} and {3}: nesting
holds everywhere and the single-child condition fails where 5 joins 2 and 1 joins 3. So here
the exact test reaches a grouping the fundamental supernodes do not, which is another way of
saying the two notions cut differently rather than one refining the other.

None of this is something md3 knows. Supernodes are defined against a forest that ordering
does not build: the point of the ordering phase is to reduce fill, and mass elimination is a
batching trick that makes it faster, not a structural analysis. The forest is computed
afterwards, from the finished permutation. Whether a useful forest could be maintained DURING
ordering, at what cost and to what end, is an open question worth a look out of curiosity
rather than a gap in what md3 is trying to do.

## Reading the mass elimination block line by line

Everything md3_eliminate does beyond md2_eliminate is one block, placed after the prune loop
and before the pivot's own containers are cleared:

```python
    merged_vertices = []
    for u in neighbors:
        if not A[u] and len(I[u]) == 1 and I[u][0] == pivot:
            I[u] = []
            eliminated[u] = True
            merged_vertices.append(u)
    if merged_vertices:                 # one compaction pass, not a removal each
        tag += 1
        for u in merged_vertices:
            mark[u] = tag
        C[pivot] = [v for v in C[pivot] if mark[v] != tag]
```

**The candidates come from the snapshot.** `neighbors` was computed before anything was
touched, so it still lists every member of the new clique even though A and I have since been
rewritten by the prune loop. Its order is the order the query produced, which is the same on both
sides, and it matters: the order of merged_vertices becomes the order of super_members[pivot],
which becomes the order of those vertices in the returned permutation.

**The test reads post-prune state.** By this point A[u] has been compacted against the clique
and the pivot dropped, and I[u] has lost the absorbed cliques and gained the pivot. So `not A[u]`
means every explicit neighbor u had was inside the clique, and `I[u] == [pivot]` means the new
clique is its only remaining route out. **So "nothing left" means "everything was inside", and
the prune is what turns one into the other**: a test framed before it would have to compare two
sets, where this one compares two lengths.

**The third conjunct is redundant.** The prune appends the pivot unconditionally, so I[u] always
contains it and `len(I[u]) == 1` already forces `I[u] == [pivot]`. `Amd.cpp` tests only
`Elen [i] == 1` for this reason. Ours is documentation of intent, not a check.

**And what the two live conjuncts give is an EQUALITY, not a containment**, which is worth having
exactly because the stronger statement is what makes the merge free:

```
reach(u) = A[u] union ( C[c] for c in I[u] ) - {u}
         = {} union C[pivot] - {u}
         = C[pivot] - {u}
```

Membership in C[pivot] supplies the adjacency, so u's reach is not merely inside the clique, it
is exactly the rest of it. Read back into the graph as it stood before the elimination, where
reach(pivot) = C[pivot] with u among its members and reach(u) contained the pivot:

```
N[u]     = reach(u) + {u}          = C[pivot] + {pivot}
N[pivot] = reach(pivot) + {pivot}  = C[pivot] + {pivot}
```

The closed neighborhoods are equal, which is the definition of indistinguishable. Eliminating u
next creates no fill the pivot has not already created, and that is the sense in which md2 -> md3
is "a reordering, free iteration by iteration" in the table near the top.

**The last two lines of the first block are the pivot's own.** `I[pivot] = []` and
`eliminated[pivot] = True` appear verbatim a few lines further down, on the pivot itself. That is
the whole content of a merge: u stops being a vertex on exactly the terms the pivot does, in the
same step.

**There is no weight array.** md3 merges only into the pivot, which is eliminated in the same
call, so no live vertex ever stands for more than one original vertex, and the size of a
supervariable is `len(super_members[pivot])` whenever it is wanted. The literature carries a
weight per vertex and later layers will need one, but here it would be a cached length that
nothing on the hot path reads. It earns its place once the members are held as chains over a
flat array rather than as lists, where a size stops being free.

**Clearing I[u] and the second block are two halves of one thing.** An incidence is stored twice,
the clique id in I[u] and the member u in C[c]. Emptying I[u] removes u's side of it; the
compaction removes the cliques' side. They are separate because one is per vertex and the other
per clique, not because of any ordering hazard: C[pivot] is a copy of neighbors, so rewriting it
cannot disturb the iteration over neighbors.

**A[u] is not cleared, and does not need to be.** The test has already established that it is
empty. Nor does anything have to remove u from other vertices' adjacency, because A is
symmetric: if A[u] is empty then no live v has u in A[v].

**Two things the block deliberately does not do.** It does not touch super_members; the driver
does that from merged_vertices, which keeps the eliminator free of expansion bookkeeping. And
at the instant the first loop finishes, C[pivot] still contains the merged vertices, which is
what the second loop fixes. That is also why the driver computes external_degree from
C[pivot] AFTER the call and not before.

**The test is sufficient, not necessary, and the gap has a name.** `len(I[u]) == 1` can fail
while u is genuinely indistinguishable: let a clique c contain u but NOT the pivot, with C[c]
lying inside C[pivot]. Ordinary absorption removes only I[pivot], the cliques containing the
pivot, so c survives, `len(I[u]) == 2`, and the merge is declined although u reaches nothing
outside the clique. That is `graph5` exactly, described in the test graphs section above. The cost
is a little quality in the bound, never correctness, and `Amd.cpp` behaves the same way. `amd2`
recovers this particular case with aggressive absorption, which kills any clique whose members
all lie inside C[p] rather than only those the pivot touched. A DIFFERENT gap needs a different
mechanism: no test framed against the pivot can see two vertices indistinguishable from each
other but not absorbable into it, which is `graph7` and takes the hashing.

**The second block touches C[pivot] only.** `I[u] == [pivot]` guarantees that u belongs to no
other clique, so compacting C[pivot] against the merged set is the whole of it, one pass rather
than a removal each. An earlier version scanned every
clique as a defense against a test that admits a u with more than one clique, which is what the
exact containment test would do; the measurement in the complexity section below made that
defense untenable, since it cost more than the work it accompanied. Should the test ever be
widened, the loop has to be widened with it.

## What mass elimination costs and saves

Two small runs, measured. The unit is elements touched by a neighbor query, `len(A[u])` plus
the sizes of the cliques in `I[u]`, since that union dominates everything else in both layers.
`scan` sums it over every candidate the pivot search evaluates. Neither layer needs anything
beyond the scan: both take the pivot's own degree from the set md*_eliminate returns, which
is computed before any merging and is therefore the same set in both.

graph1, the 4-cycle, where md3 merges at the last iteration:

```
        scan   iterations   candidate touches
md2      21      4            10
md3      16      2             7
```

The merge takes 2 and 3 into the pivot at iteration 1, which ends the run, so md3 never runs md2's
iterations 2 and 3. Nothing downstream is visible here, which is the limitation of this example.

graph6, where the merge lands at iteration 2 of 5 and the run continues:

```
        scan   iterations   candidate touches
md2      62      6            21
md3      53      5            18
```

Iteration by iteration, the two runs are identical through iteration 2: same live sets, same pivots,
same scan costs of 16, 15 and 14. Vertex 4 merges into 0 at iteration 2, and from there md3 scans
[2, 3] then [3] where md2 scans [2, 3, 4], [3, 4] and [4]. One fewer candidate in every
remaining scan, and one fewer iteration at the end, which is 62 down to 53.

Two effects are folded into that number and worth separating. Fewer candidates: a vertex
merged at iteration k is absent from every scan after k, so it stops being touched once per
remaining iteration rather than once. And cheaper candidates: the merged vertex is stripped from
every clique, so the queries that remain walk shorter member sets. On graph6 the surviving
queries at iteration 3 cost 3 each against md2's 4, purely because C[c0] lost a member.

## External degree, new in md3

External degree is the weighted count of the original vertices a supervariable can still
reach, not counting the ones it already stands for. Two exclusions are packed into the word,
and both matter.

The first is the supervariable's own members. A supervariable of size 3 is three mutually
adjacent original vertices. Each of them has two neighbors inside the group, but those are
not a choice the ordering makes: the three are eliminated consecutively and the entries they
create are the w(w-1)/2 triangle that is going to be there whatever else happens. Counting
them would inflate the degree of exactly the vertices mass elimination has just rewarded. In
the code the exclusion is automatic rather than tested for: md3_neighbors drops u itself, and
a merged vertex is stripped from every clique's member set.

The second is the weighting, and here md3 is a special case worth being explicit about. In
general a neighbor is a supervariable too, so it should count for the number of original
vertices it stands for rather than for 1, which is how the literature defines the external
degree. In md3 that never bites: the only vertex that could stand for several is the pivot,
which is eliminated in the same call, and merged vertices are stripped from the cliques, so
every neighbor a scan ever sees stands for exactly one original vertex. md3 therefore has no
weight array at all and its pivot search uses len, exactly as md2's does. Measured by
swapping a weighted key for a plain len over 1000 random graphs, before the array was
removed, neither the order nor nnz(L) changed on any of them.

The quantity appears twice in md3, at different instants, and the two values differ. The
selection degree is evaluated over the neighbors as they stand before the iteration.
external_degree in the nnz(L) block is recomputed from C[pivot] after md3_eliminate returns,
because vertices that were neighbors a moment ago are members now. The difference between the
two is exactly the total weight that merged during the iteration, and confusing them is the
double-counting trap described in the previous section.

The star makes it concrete. At the iteration that eliminates hub 0 with leaf 4 still live, the
selection degree is 1, since 4 is a live neighbor. After the merge the supervariable has size
2 and external_degree is 0, since C[0] no longer contains 4. The supervariable's two columns
then hold one entry and zero entries below their diagonals, and the closed form gives
2 * 0 + 2 * 1 / 2 + 2 = 3, which matches.

The term comes from the AMD literature, where the contrast is with the TRUE degree, which
does count the internal members. The paper's relation is exact, `d_i = t_i - |i| + 1` for a
supervariable `i` of size `|i|`, so `t_i = d_i + |i| - 1`, and the two coincide when `|i| = 1`.

**Both names are worse than the quantities deserve, and the plain reading is worth having beside
them.** "True" suggests correctness and "external" suggests a contrast with an internal degree
that does not exist. What the pair actually encodes is a granularity:

- `d_i`, the external degree, is the reach outside the supervariable, so it is the BORDER of the
  supernode's front, and the update matrix is `d_i` by `d_i`. It is the update size per SUPERNODE.
- `t_i` is that border plus the `|i| - 1` other members sitting below the first column, so it is
  the count below the diagonal in the supernode's FIRST COLUMN. It is the update size per COLUMN.

Both are therefore update sizes at two granularities, and the other dimension is the front size
`|i|`, which is 1 for a plain column. That is the same pair the numeric phase already names,
front and update, so the ordering vocabulary and the factorization vocabulary can describe one
thing rather than two.

We keep "external degree" only because the papers and the vendored routines use it and the amd
layers are read line by line against both; renaming would cut that thread for no gain. But `t_i`
is never computed here, so the qualifier distinguishes nothing INSIDE this codebase, and where
the exact quantity has to be contrasted with AMD's approximation the phrase is EXACT DEGREE and
never "true degree", which in the paper names something else and larger.

Minimum degree with supervariables uses the external one,
and the choice is not cosmetic: it changes which pivot is selected, not merely what number is
printed beside it. Where the weighting first has that effect is amd2, whose hash detection
folds a vertex into a LIVE supervariable, `weight[i] += weight[j]` with i still a candidate.
Through mmd1 every transfer is into the pivot, so the invariant that a live vertex stands for
one original vertex holds all the way, and the weighting stays inert.

## Counting nnz(L), md1 through md3

All three layers report nnz(L) without ever forming L, and all three are evaluating the same
identity at different granularity. The identity is that at the moment a vertex is eliminated
its neighbor set is exactly the off-diagonal pattern of its column of L. Every live neighbor
is eliminated later, so it sits below the pivot in the permuted order; every neighbor already
eliminated was discarded from the set at that earlier iteration, and its entry belongs to that
earlier column. Fill created at this iteration lands among the neighbors, in later columns, and
cannot alter a set that was read before the elimination ran.

md1 and md2 use it directly, one iteration being one column. The pivot degree is the column count,
so a running degree_sum plus n diagonals is nnz(L), with no second pass and nothing else
stored. The two layers compute the neighbor set differently, md1 reading A[pivot] and md2
unioning the explicit adjacency with the cliques, but the sets are equal at every iteration, so
the counts are too.

md3 cannot do that, because one iteration is now w columns rather than one. A supervariable of size
w stands for w original vertices, eliminated consecutively, and their columns differ. Let ext
be the external degree, the count of live neighbors remaining after the merges. The first of
the w columns sees ext outsiders plus the w - 1 members of its own supervariable still to
come, the second sees ext + w - 2, and so on down to the last, which sees ext alone.
Summing:

```
off-diagonals = w * ext + (w-1) + (w-2) + ... + 0 = w * ext + w(w-1)/2
diagonals     = w
```

which is the expression in md3_minimum_degree, term for term. The middle term is fill inside
the supervariable, the dense triangle among its own members, and it is certain to be there
because indistinguishable vertices are mutually adjacent by definition.

Three things follow, and the second is a trap worth naming.

The formula reduces to md1's rule. With w = 1 the triangle term vanishes and ext is the
ordinary degree, giving degree + 1 per iteration, which summed over all iterations is degree_sum +
n. So md1 and md2 are the special case of md3's accounting, not a different scheme.

The degree computed at selection is the wrong number to accumulate. At that moment the
vertices about to merge are still neighbors and still counted, so degree = ext + w - 1, the
pivot itself being the one member it already had. Accumulating from it double counts them,
once as neighbors and once as members of the supervariable. This is why md3 recomputes
external_degree from C[pivot] after md3_eliminate returns, rather than reusing the value it
printed in the iteration title. Getting this wrong reports 43 against the true 37 on graph3.

The accumulation has to be per iteration. Both ext and w vary from iteration to iteration, and
neither is recoverable at the end from a single running total, so md3 evaluates a closed form each
iteration
instead of summing a per-column quantity. What it avoids is materializing the expansion: it
never walks the w columns one at a time, only their sum.

The check on all of this is external rather than internal. On 34 graphs the nnz(L) md3
reports matches a symbolic factorization run independently on the order md3 returns, and
equals what md2 reports on the same graph. That equality is the claim the layer exists to
make, that mass elimination is free and not merely cheap. md1 has a second, cheaper check
built in: it counts fill edges one at a time as it creates them, and prints that total
alongside nnz(L) - nnz(tril A), which arrives by the degree route. The two agree only if the
column-count identity holds and no fill edge is missed or double counted. md2 and md3 have no
such cross-check, since neither materializes fill; computing it would cost a neighbor query
per neighbor plus a quadratic pair scan, which is work the representation exists to avoid.

### md3 with an approximate degree, where only the counting changes

The four lines from md2, unchanged as text:

```
reach(u)  = ( A[u] | C[c] for c in I[u] ) - {u}
degree(u) = |reach(u)|

reach(u)  = ( A[u] - C[p] ) | ( C[p] - {u} ) | ( C[c] - C[p]  for c in I[u] - {p} )
bound(u)  = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]|
```

**What changes is what |S| means.** A vertex now stands for a set of original vertices, so every
size is a weighted count:

```
|S| = sum of weight(v) over v in S,   where weight(v) = |super_members[v]|
```

That substitution is uniform. It lands on degree(u) and on bound(u), in the same places, so
degree(u) <= bound(u) survives it and no term of the decomposition moves.

**Mass elimination touches the anchor, and this is the one place md3 is not merely a substitution.**
C[p] is no longer reach(p): the vertices the pivot absorbs are trimmed out in the same call, so

```
C[p] = reach(p) - merged
```

and both C[p] terms of the bound must use the trimmed set. A merged vertex has joined the pivot's
supervariable, so it is inside the diagonal block rather than outside it, and counting it in
|C[p] - {u}| inflates every bound the iteration produces. `Amd.cpp` does the same thing at
`degme -= nvi`, and our own version of this was a real bug, recorded in the bug section above.

**Mass elimination itself is indifferent.** A[u] == {} and I[u] == {p} is a structural test that
never consults a degree, so it fires the same whichever of the two numbers the layer is computing.

md3 also still repicks from scratch, so it is where md2 is: the estimate is definable and not worth
using. The layer that changes that is md4.

## md4: maintained degrees

Every layer so far rebuilds a neighbor set for every live vertex at every iteration, keeps the
smallest and throws the rest away. md4 keeps the degrees in an array and updates only what
can have changed, so the pivot search reads integers instead of building set unions.

What can have changed is exactly the surviving members of the new clique, and the argument
comes in three parts. PRUNING and clique membership are rewritten only for members of that
clique. ABSORPTION deletes the cliques the pivot belonged to, and every member of such a
clique is reachable from the pivot, hence inside the new one. MERGING removes a vertex, but a
vertex merges only because everything it could see lay inside the new clique, so nobody
outside can tell it is gone. Nothing else in the graph can detect that an elimination
happened.

The layer is three fragments in the driver, and nothing else changes: md4_neighbors,
md4_storage and md4_eliminate are md3's word for word.

**Construction.** One pass, correct because no clique exists yet, so the whole neighborhood
is still explicit.

```python
    degrees = [len(A[u]) for u in range(n)]
    num_degree_computations = n
```

**The pivot search.** md3's key was `len(md3_neighbors(A, I, C, u))`, a full union per
candidate. Here it is an array read, and the scan over range(n) is what md5 removes.

```python
        pivot = min((u for u in range(n) if not eliminated[u]),
                    key=lambda u: degrees[u])
```

**The degree update pass.** After the elimination, and only over the clique's survivors: the merged
vertices were stripped from C[pivot] inside md4_eliminate, and the pivot was never in it. The
zeroing that follows is hygiene rather than logic, since eliminated[u] already keeps dead
vertices out of the search, but it keeps a dead slot holding a neutral value rather than a
stale one. In md5 the same two lines become load-bearing.

```python
        refreshed_vertices = list(C[pivot])
        for u in refreshed_vertices:
            neighbors_u, tag = md4_neighbors(A, I, C, mark, tag, u)
            degrees[u] = len(neighbors_u)
        num_degree_computations += len(refreshed_vertices)
        degrees[pivot] = 0
        for u in merged_vertices:
            degrees[u] = 0
```

Nothing is sorted. C[pivot] is a list in both twins and the degree update pass walks it as it
stands, which is the same sequence on each side. The tag comes back from the query because a degree
computation
advances it and Python has no reference parameters.

Two query sites remain in the whole file, and that is the layer's claim: md4_eliminate's
first line, which becomes the clique, and the degree update pass above. Everything else reads
integers. The run prints the count, and on graph3 it is 34 against the 85 neighbor queries md3
makes,
on graph4 20 against 35, on graph6 11 against 23.

The display follows the same rule. md4_show prints the stored degree and never recomputes it,
so the trace shows what the algorithm believes rather than what is true. An earlier draft
printed both, as a check on the cache, and the two columns were always equal; the check moved
out to a test rather than sitting in the output.

### md4 with an approximate degree, where it becomes usable

The same four lines, still unchanged as text:

```
reach(u)  = ( A[u] | C[c] for c in I[u] ) - {u}
degree(u) = |reach(u)|

reach(u)  = ( A[u] - C[p] ) | ( C[p] - {u} ) | ( C[c] - C[p]  for c in I[u] - {p} )
bound(u)  = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]|
```

**What changes is WHEN they are evaluated, and for whom.** md2 and md3 evaluate the first pair for
every live u at every iteration. md4 evaluates it only for u in C[p], because nothing else can have
changed. Written as the layer sees it:

```
for u in C[p]:  degrees[u] = degree(u)        md4, exact
for u in C[p]:  degrees[u] = bound(u)         md4, approximate
```

and those two loops range over the same set. That is the coincidence the whole placement rests on.
The second line of the decomposition is a statement about members of C[p] and about no others, so
at md2 and md3 most of what gets updated is outside its domain and has to be counted exactly
anyway. At md4 the degree update set IS the domain. bound(u) is inherently incremental, md4 is where
the algorithm becomes incremental, and there is nothing left over.

**A second thing arrives with it.** Of the two caps in the md2 subsection,

```
degree_old[u] + |C[p] - {u}|
```

needs degree_old[u], and a cache is what makes a previous degree exist. It is often the tightest of
the three, since a vertex whose reach barely grew is bounded by what it was plus the new clique. It
also stays valid when the cached value is itself a bound rather than a degree, which is what makes
the chain work inductively: degree_old[u] >= degree(u) at every earlier iteration means it is still
an upper bound now. Nothing has to be exact anywhere.

**And the costs finally differ.** Both loops walk C[p], but degree(u) unites the members of every
clique in I[u], per vertex, while bound(u) reads one number per clique, and that number,
|C[c] - C[p]|, is computed once for the whole iteration. So md4 is the first layer where the two
columns are doing measurably different amounts of work rather than the same work on different
populations.
The decomposition is md2's, the counting is md3's, and what md4 contributes is a place to stand.

## md5: degree buckets

md4 left one O(n) per iteration in place: the scan over every live vertex to find the smallest
cached degree, now over integers rather than unions, but still a full pass. md5 files each
live vertex in a bucket indexed by its degree, so the minimum is found by walking UP from the
last known bound rather than by looking at everything. Both MMD and AMD do this; neither
invented it.

The whole of it is five fragments in the driver plus a three-line helper. md5_neighbors,
md5_storage and md5_eliminate are md4's word for word.

**Fragment 1, construction.**

```python
    buckets = [[] for _ in range(n)]           # buckets[d] holds the live degree-d
    filed = [False] * n                        # whether u is in a bucket at all
    for u in range(n):
        md5_file(buckets, filed, degrees[u], u)
    min_degree = min(degrees) if n else 0
    num_bucket_probes = 0
```

A bucket is a LINKED LIST, not a set. The C++ twin holds `head[d]` with `next` and `prev` over
n, which is what MMD's fwd/bwd and AMD's Next/Last are, and push, pop and splice are all O(1).
An ordered container cannot give that: `std::set` costs O(log n) per file and unfile, and those
happen once per degree change, which is the dominant operation of the whole algorithm. The
Python mirrors the same sequence with a list whose position 0 is the head, so both twins hold
the same order at every iteration and pick the same pivot; it pays O(bucket) for insert and remove,
which is the one place it is asymptotically behind its twin.

n slots is exactly right, indices 0 through n - 1. A live vertex counts only live neighbors,
so its degree is at most n - 1, and the walk stops at the first non-empty bucket, which
exists while anything is live, so no index above n - 1 is ever filed or probed. SuiteSparse
AMD sizes its Head array the same way; MMD's head array looks one longer only because it is
indexed from 1, and it files an isolated vertex under degree 1 rather than 0, so it has no
degree-0 bucket at all. Every vertex is filed here, isolated ones included. min_degree starts
at the true minimum, the tightest legal value for a lower bound; starting at 0 would also be
correct and would cost one extra walk on the first iteration.

**Fragment 2, the pop.**

```python
        while not buckets[min_degree]:         # walk up to the first live bucket
            min_degree += 1
            num_bucket_probes += 1
        num_bucket_probes += 1
        pivot = buckets[min_degree][0]         # the head, whatever was filed last
```

The walk only ever climbs, and min_degree is never reset between iterations, so the work is
amortized across the run rather than paid per iteration. Termination rests on the outer loop
condition: some live vertex exists, it is filed under its own degree, and that degree is at
or above the bound, so a non-empty bucket is found before the array ends.
The pop takes the head, and that is the tie-break: whatever was filed last, not the lowest
index. md1 through md4 scan range(n) and keep the first strict minimum, so ties there go to the
lowest index, and this does not reproduce it. That is the price of an O(1) bucket, and it is why
md5's permutation differs from md4's; see the section on what the layers show for the
measurements.

**Fragment 3, the deletions.**

```python
        md5_unfile(buckets, filed, degrees[pivot], pivot)   # the pivot has left
        degrees[pivot] = 0
        for u in merged_vertices:               # and so have the merged vertices
            md5_unfile(buckets, filed, degrees[u], u)
            degrees[u] = 0
```

This is where md4's harmless zeroing becomes necessary. There a dead vertex's entry was never
read again because eliminated[u] filtered it out of the scan; here there is no scan, and a
dead vertex left in a bucket would be popped as a pivot on a later iteration. The order within
each pair matters: the bucket index is read from degrees[u], so the removal has to come
before the zeroing, or the vertex is erased from buckets[0] and left where it was.

**Fragment 4, the degree update pass.**

```python
        refreshed_vertices = list(C[pivot])
        for u in refreshed_vertices:
            neighbors_u, tag = md5_neighbors(A, I, C, mark, tag, u)
            md5_refile(buckets, filed, degrees, u, len(neighbors_u))
        num_degree_computations += len(refreshed_vertices)
```

```python
def md5_refile(buckets, filed, degrees, u, new_degree):
    md5_unfile(buckets, filed, degrees[u], u)
    degrees[u] = new_degree
    md5_file(buckets, filed, new_degree, u)
```

Same degree update set as md4, and the only change is that the new degree goes through the helper so
the bucket moves with it. The helper exists so the three iterations cannot be written half-way. A
vertex whose degree did not change is removed and reinserted into the same set, which is
harmless. Removal from the middle of a bucket must be O(1), which is why a bucket is a linked
list rather than an ordered container. `filed[u]` makes unfiling idempotent, which matters in
mmd1 rather than here: a vertex evicted early in a batch can be merged away by a later pivot in
the same batch, and unfiling it twice would splice a list it is no longer in. With sets that was
a harmless `discard`; with a linked list it would corrupt the bucket head.

**Fragment 5, lowering the bound.**

```python
        min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])
```

Only the vertices just updated can have moved down, so that is the complete candidate set. The
bound is lowered only here and raised only by the walk in fragment 2, which is what keeps the
total climbing work bounded over the run. Lowering is the safe direction: a bound that is too
low costs a few extra probes, while a bound that is too high skips a non-empty bucket and
picks the wrong pivot. Degrees that rise need no attention at all, since such a vertex is
filed higher and the walk will reach it.

Three invariants hold after every iteration, and they are what the buckets row and the min degree
line in the trace exist to show: every live vertex is in the bucket matching its degree, no
cached degree is stale, and no live vertex has a degree below min_degree. Checked across 400
random graphs, no violations.

The run prints both metrics. Degree computations are unchanged from md4, as they must be,
since this layer touches only how the minimum is found. Bucket probes replace what was an
n-per-iteration scan: 12 on graph3, 5 on graph4, 7 on graph6.

### md5 with an approximate degree, where only the filing changes

The same four lines once more, and now a fifth that is the layer:

```
reach(u)  = ( A[u] | C[c] for c in I[u] ) - {u}
degree(u) = |reach(u)|

reach(u)  = ( A[u] - C[p] ) | ( C[p] - {u} ) | ( C[c] - C[p]  for c in I[u] - {p} )
bound(u)  = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]|

for u in C[p]:  refile(u, <the number>)       and the picker takes the head of the
                                              first non-empty bucket from min_degree
```

**The layer never looks at where the number came from.** refile takes a value and moves u between
two buckets; the walk takes the least filed value. Neither operation can tell degree(u) from
bound(u), so md5 needs no adjustment for the estimate and the estimate needs none for md5.

The three invariants above hold verbatim with "degree" read as "stored value": every live vertex is
in the bucket matching its stored value, no stored value is stale, and no live vertex has a stored
value below min_degree. The one requirement is uniformity. min_degree is a lower bound on the
minimum STORED value, so a mixture, some vertices holding degree(u) and others holding bound(u),
would break it, since the two are not comparable as estimates of the same thing. Choose one and
every argument in this section goes through unchanged.

**What the estimate costs, wherever it is used, is the pivot itself.** Every layer in the exact
chain picks a true minimum-degree vertex, and md5's own tie-break note above is careful that only
the choice among equals moves. bound(u) gives that up: an overcount can hide the true minimum, so
the head of the first non-empty bucket may simply not be minimal. That concession belongs to the
estimate and not to any layer, it becomes payable at md4 where the estimate becomes usable, and
amd1 is where it is actually paid.

## The buckets look like a bucket priority queue and are not one

The structure invites a comparison it does not survive, and the difference is worth stating
because it decides the implementation.

A textbook bucket priority queue needs two operations, insert and extract-min, and extract-min
only ever touches the HEAD of the lowest non-empty bucket. If the pivot were the only departure
each iteration, a stack per bucket would be enough, since the pivot already sits at a head:
`buckets[min_degree][0]`.

**Three things here need removal from the MIDDLE, which is what rules a stack out.**

- The merged vertices, dead but not the pivot, sitting wherever their old degrees put them.
- Every refreshed vertex, which has to leave the bucket for its old degree and enter the one for
  its new degree. That is the whole of `md5_refile`, and it is the dominant operation of the
  algorithm, running once per degree change rather than once per pivot.
- The pivot, which happens to be at a head but is not treated specially.

So a bucket is a DOUBLY linked list and the links are indexed BY VERTEX rather than by slot: a
vertex is its own node, so `unfile` unlinks through `prev` and `next` in O(1) with no search. A
singly linked list would give O(1) at the head and O(bucket) in the middle, which is exactly what
the Python mirror pays and why its own note calls that the one place it is asymptotically behind
the C++.

**`degrees[u]` does double duty, and that is what makes the ordering of two statements matter.**
It is both the value and the index saying which bucket holds `u`, so the two are consistent only
while both are current. Hence

```python
md5_unfile(buckets, filed, degrees[u], u)   # reads degrees[u] to know WHICH bucket
degrees[u] = 0                              # so this must come second
```

Zero first and the unfile searches `buckets[0]`, where `u` is not. In md4 the same assignment is a
lone cache write with nothing reading the old value, which is why md4 can order those two lines
freely and md5 cannot.

### What `min_degree` actually is, and why it is not the minimum

```python
min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])
```

The seeded list is not a trick to avoid `min([])`, though it does that too when the clique
emptied entirely. **The seed is a stand-in for every live vertex the iteration did not touch**,
and it is a lower bound on them rather than their minimum, because nothing re-examines them.
Three ways an iteration can end:

- **The new minimum is at a refreshed vertex.** The expression saw it. Exact.
- **The old bucket still has an occupant** after the pivot and the merged vertices left. The
  minimum is unchanged and the seed returns it. Exact, and for free.
- **The old bucket emptied and nothing refreshed that low.** Now it is a strict lower bound, and
  the next walk-up climbs to the truth.

The second case is the common one on any graph of size, which is why the looseness is rare.
Measured on md5's own examples by comparing `min_degree` against a full scan of the live
vertices: three iterations out of the run are strictly low, by 1, 1 and 2.

The asymmetry is the entire correctness argument. **Too low costs iterations of a `while` loop;
too high would skip a non-empty bucket and return a non-minimal pivot**, silently producing a
worse ordering with nothing to detect it. So the bound is relaxed only by the refreshed vertices,
the ones known to have moved, and never by anything more thorough.

Stated as one loop invariant: `min_degree` is a lower bound on the degree of every live filed
vertex, carried across iterations, tightened upward by the walk at the top and relaxed downward by
the `min` at the bottom. Neither end needs it to be exact. The first bound is the one exception,
`min(degrees)` at construction, which costs one O(n) scan and is exact, since starting at 0 would
only make the first walk climb a distance the degrees were already sitting there to tell us.

The economy is that total walking over the run is bounded by the total distance the bound ever
travels UPWARD, not by iterations times range. And the tension is that every downward correction
is distance that will be re-walked later, which is the reason the relaxation is kept as narrow as
it is.

## mmd1: multiple elimination

Updating degrees is the expensive part, so do it less often: eliminate a whole INDEPENDENT
SET of least-degree vertices before updating any degree. Non-adjacent pivots cannot disturb
each other's degrees, so every pivot in a batch is still a true minimum-degree vertex when it
is taken. That is Liu's M in MMD, and it is the first layer whose ordering differs from md1's
for a reason other than a tie in the same graph state.

**This is the fork, and it is a real one.** md1 through md5 is a single chain because bound(u) is
one definition that no layer alters, and because the layer that makes it usable, md4, does so
without any change to itself, as the notes at the end of md2, md3, md4 and md5 set out. At the top
that independence ends. mmd and amd attack the same cost, the degree update, from opposite ends: mmd
makes it rare, amd makes it cheap, so their gains overlap rather than add. They also interfere. The
anchored bound is anchored at ONE new clique, and a batch produces several, so an amd-flavored mmd
would have to re-anchor against their union and would lose the per-clique reuse that makes the bound
worth having. And the bound stays tight only when degrees are updated often, which is precisely
what batching gives up: the delta measurement below shows the bound's own failure mode in the exact
setting, where a wide batch leaves a whole evicted set invisible for an iteration. So mmd1 and mmd2
stay exact and take the batch; amd1 and amd2 stay one pivot per iteration and take the bound.
Combining them
is a question rather than an iteration, and it is not one this experiment answers.

Six of the seven functions are md5's with the prefix changed: mmd1_neighbors, mmd1_storage,
mmd1_eliminate, mmd1_refile and the two display functions. That is the pattern across the
whole ladder from md2 onward. The elimination itself has not changed since the quotient graph
appeared, and the degree cache and buckets have not changed since md4 and md5. What each layer
varies is the SELECTION POLICY: recompute per candidate, cache, bucket, and now batch. So for
mmd1 the whole layer is the driver.

**Eviction, since the word is used throughout and means one narrow thing.** To EVICT a vertex is
to unfile it from its bucket and nothing else. It stays live, it keeps its stale value in
degrees[], it is simply not visible to the picker until the degree update pass files it back at
its new degree. Every member of C[pivot] is evicted, which is to say every vertex the pivot
reached.

md5 has no eviction, and that is the cleanest way to see what the word is for. There the unfile
and the file are one call, mmd1_refile, three lines apart with the degree computed between them,
so a vertex is never out of the buckets across a decision point. Eviction is what appears when the
file is DEFERRED: mmd1 unfiles during the batch and files at the pass, the two come apart, and the
interval between them needs a name.

It is not extra work. md5 pays the same unfile inside its refile, and the pass files without
unfiling precisely because the eviction already did. What batching saves is the middle, the degree
update, done once for a vertex several pivots reached rather than once per pivot.

One consequence worth stating: a vertex in two pivots' cliques is evicted twice in one iteration,
and the second unfile is a no-op guarded by `filed[u]`. That guard is load-bearing, and it also
covers a vertex evicted and then mass-merged away later in the same iteration.

**The independent set is never searched for.** It falls out of the eviction:

```python
            for u in C[pivot]:                 # EVICT, with a stale degree
                buckets[degrees[u]].discard(u)
                touched.add(u)
```

Every vertex the pivot reached leaves the buckets and stays out until the iteration ends. So
whatever is still filed was not reached by any pivot taken so far, hence is not adjacent to
any of them, and draining a bucket drains an independent set. The eviction is also what makes
the deferred degree update safe: a vertex with a stale degree is not a candidate, because it is not
in a bucket to be found in.

**The batch loop takes at least one pivot, then consults the limit.**

```python
        batch_limit = min_degree + delta
        while True:
            if not buckets[min_degree]:        # this degree is drained
                if min_degree >= batch_limit:
                    break
                min_degree += 1
                num_bucket_probes += 1
                continue
            pivot = buckets[min_degree][0]
            ...
            if delta < 0:                      # one pivot per iteration, as md5 does
                break
```

The shape matters. On entry the outer walk has left buckets[min_degree] non-empty, so the
first iteration always takes a pivot, whatever the limit says. Only after that does delta
decide whether the iteration continues. Written the other way iteration, as a `while min_degree <=
batch_limit` guard, a negative delta makes the loop body unreachable, the batch comes out
empty, nothing is eliminated and the driver spins forever. That was a real bug in the first
draft of this rewrite.

**delta is the whole control, and its sign selects between two behaviors.** delta = 0 keeps
the batch to true minima. delta > 0 admits vertices up to delta above the minimum, which are
not minimal, so that is a concession in quality for still fewer degree updates. delta < 0 takes one
pivot per iteration, which is md5 reached through this code path, and it is the check on the
batching rather than a feature of it: mmd1 at delta = -1 reproduces md5's ordering exactly,
verified on the seven graphs and 150 random ones.

**The degree update pass runs once per iteration, over everything the batch touched.**

```python
        refreshed_vertices = [u for u in touched if not eliminated[u]]
        for u in refreshed_vertices:
            neighbors_u, tag = mmd1_neighbors(A, I, C, mark, tag, u)
            degrees[u] = len(neighbors_u)
            mmd1_file(buckets, filed, degrees[u], u)
        num_degree_computations += len(refreshed_vertices)
        min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])
```

Note the asymmetry with md5's degree update, which called mmd1_refile: here the vertex is already
out of its bucket, evicted during the batch, so the degree update only writes the degree and files
it. The `not eliminated[u]` filter matters because a vertex evicted early in the iteration can be
merged away by a later pivot in the same iteration.

**Choosing delta, measured.** The vendored driver passes 0, which is also what SPARSPAK does and
what Liu's paper treats as the default. On grids, with fill and degree update count both reported:

```
grid 22x22, n=484
   delta   -1: nnz(L) 4773   degree computations 2690   degree updates 2206   iterations 367
   delta    0: nnz(L) 4684   degree computations 1859   degree updates 1375   iterations  36
   delta    1: nnz(L) 4754   degree computations 1733   degree updates 1249   iterations  22
   delta    2: nnz(L) 4706   degree computations 1756   degree updates 1272   iterations  21
   delta    4: nnz(L) 4747   degree computations 1601   degree updates 1117   iterations  14
   delta    n: nnz(L) 5964   degree computations 1514   degree updates 1030   iterations   9
```

**Both middle columns are here because they differ by exactly n, on every row, and n is 484.** The
initial build is one degree computation per vertex and no vertex has been updated yet, so degree
computations is degree updates plus n forever after. Degree updates is therefore the column to
compare, carrying only work the layer can change, and degree computations is the column to quote,
being what the closing line prints and what every figure in `REPORT.md` and
`benchmarks/ordering/README.md` already reports. Neither is wrong and they are not the same
number; a comparison drawn from the wrong one understates every saving by a constant.

Two things in that table are worth more than the recommendation they support.

The first is that delta = 0 beats delta = -1 on BOTH axes, here and on the 10 by 10 and 16 by 16
grids as well. Batching is not trading quality for speed at that setting; it is simply better.
The wager appears only as delta grows, and by delta = n it is decisively bad: 27 per cent more
fill for 19 per cent fewer degree updates.

The second is why it goes bad, which is not stale degrees. Every pivot in a batch has a CORRECT
degree, since anything whose degree could have changed was evicted. What the iteration cannot see is
the evicted set, and those are exactly the vertices whose degrees typically FELL, so they are the
candidates that should be picked next. With delta = 0 the iteration ends as soon as the minimum
bucket drains and they come back at once. With delta = n the iteration keeps climbing through the
degrees, taking vertices of degree 4, 5, 6 while better candidates wait until the end.

**The same comparison across five sizes, which is the other axis.** The table above fixes the grid
and sweeps delta; these fix delta at its two ends and sweep the grid, so the 22x22 rows are the
same two runs seen from the other direction. delta = -1 is md5 reached through this code path, so
the left half is the layer below and the right half is this one.

```
                delta = -1 (md5)          delta = 0
             nnz(L)  iters  degupd     nnz(L)  iters  degupd
10x10           657     78     376        636     19     290    -3.20%
16x16          2195    197    1116       2088     31     751    -4.87%
22x22          4773    367    2206       4684     36    1375    -1.86%
30x30         10436    679    4327      10757     40    2533    +3.08%
40x40         22495   1207    8132      21614     56    4498    -3.92%
```

**The fill is two-sided, and that is the finding rather than the size of it.** Four of the five
grids come out ahead and 30x30 comes out 3 per cent behind, with no trend in n, which is what an
arbitrary choice looks like when it is measured. So delta = 0 is not spending fill to buy degree
updates. The wager is real and it is what delta > 0 spends, as the previous table shows at
delta = n.

The second half of the same runs, counting the work rather than the saving:

```
                 delta = -1 (md5)                    delta = 0
              elims   sum|C|    fill  probes      elims   sum|C|    fill  probes
10x10            78      376     377      89         82      395     356      31
16x16           197     1116    1459     217        201     1108    1352      55
22x22           367     2206    3365     399        373     2198    3276      68
30x30           679     4327    7796     723        687     4358    8117      90
40x40          1207     8132   17775    1273       1215     8010   16894     121
```

**Eliminations rise about half a per cent on every grid, and the mechanism is worth knowing
because it is the two batchings competing.** A vertex the batch takes as a pivot is sometimes one
a later pivot would have absorbed for free by mass elimination, and a merge costs a compaction
where a pivot costs a whole elimination, so where the two overlap the batch takes the more
expensive route. graph1 is the case small enough to read: md5 eliminates 3, refreshes, and then 2
absorbs both 1 and 0 as a supervariable of size 3, where mmd1's batch takes 3 and 1 together, 1
being the opposite corner of the 4-cycle and so not evicted, leaving 2 only 0 to absorb. Two
eliminations against three, `sum |C[p]|` 2 against 4, and identical fill, iterations and degree
updates. So graph1 is not merely a graph where batching buys nothing; it is one where batching
costs and buys nothing. At 1600 vertices the same effect is still there and is half a per cent.

`sum |C[p]|` is a wash, up on two grids and down on three, so those extra pivots are cheap ones.
**Bucket probes are the column that moves most, a factor of ten at 40 a side**, and that is
structural rather than a saving of the same kind: md5 re-enters the walk-up once per pivot, so it
pays a probe per elimination plus the climbs, while the batch enters once per iteration and then
walks inside the drained bucket, so the count tracks iterations rather than pivots.

**delta is total, and the top end is clamped.** Any negative value means one pivot per iteration; 0
through n - 1 widen the window; anything larger saturates at n - 1, since a live vertex's degree
cannot exceed that. That range also fixes its type: delta is signed, it is compared against a
degree, and it stops being meaningful at n - 1, which is itself bounded by the index type. So it
is an index-like quantity by Oblio's rule, a std::int32_t rather than a count.

The clamp is not cosmetic: without it the walk indexes past the last bucket and crashes, which it
did the first time delta = n was tried. The vendored genmmd has the same
latent bug and no clamp, `mdlmt = mdeg + delta` with the walk climbing until it passes mdlmt while
indexing head[mdeg]; it never bites because mmd_order always passes 0.

**What is given up is not what one would guess.** The pivots are exact: every one was a true
minimum when it was taken. What the batch loses is the vertices it evicted, which are
invisible for the rest of the iteration, so the choice is made among the untouched remainder
rather than among all candidates. The batch does not pick a worse vertex, it picks a different
vertex of the same degree. Minimum degree is famously sensitive to tie-breaks, so the fill
moves by a fraction of a percent, in either direction. Batching across connected components is
free, since the components cannot interact at all; batching within one is the wager.

**The metric is iterations against pivots.** The closing line prints degree computations, bucket
probes and iterations, and the ratio of pivots to iterations is the average batch size, which is
what the batching buys. md5 and mmd1 pay the same per degree update; mmd1 pays for fewer of them. On
graph3, md5 makes 36 degree computations over 11 iterations and mmd1 makes 30 over 6 iterations; on
graph4, 20 over 5 against 15 over 3. graph1 is the case where batching buys nothing, 6 either
way, since md5 already finishes it in two iterations.

**Iteration, degree update, degree update pass: the words, said exactly.** Three levels, and they
need three names because we count one and schedule another.

- A **degree update** is ONE vertex's degree recomputed. It is the unit of cost, and the expense
  is real: a union over A[u] and the cliques in I[u].
- A **degree update pass** is the procedure over a set of them. It is the unit of scheduling.
- An **iteration** is the eliminations between two passes.

An iteration is one batch of eliminations followed by one degree update pass, and the two are in
bijection: one iteration, one pass, always. They are not synonyms, they are the two halves of one
cycle, the iteration being the eliminations and the pass the degree work that closes it. That is
why the closing line can print iterations alone; the number of passes is the same number.

The batch is k >= 1 ELIMINATIONS, and an elimination is a pivot together with whatever mass-merges
into it, so k counts pivots and not vertices: in the worked example below a batch of 3 removes 5
vertices. md5 is the case k = 1, and mmd1 at delta < 0 is md5 exactly. The word is deliberately the
same one it has in maximum cardinality matching, where a plain augmenting-path algorithm does one
path per iteration and Hopcroft-Karp does a maximal set of vertex-disjoint ones, for the same
reason: the expensive part is a global recomputation, so one unit of progress per recomputation
wastes it. The correspondence is close, iteration to iteration, elimination to augmenting path,
non-adjacent pivots to vertex-disjoint paths, degree update pass to the BFS layering. It stops at
one place worth remembering. Hopcroft-Karp is exact and batching is a pure speedup with a
sqrt(V) bound on iterations; minimum degree is a heuristic and batching CHANGES the answer, with no
bound on iterations, only measurement.

Two things follow from keeping the levels apart. The bare plural always means vertices: "two degree
updates" is two vertices, never two passes. And iterations and degree updates are separate savings.
Batching removes an iteration by construction; it removes degree updates only when the batch's
pivots reach overlapping neighborhoods, since a vertex touched by three pivots is updated once
rather than three times.

The qualifier is never dropped. UPDATE alone means the supernodal update, as in `UpdateMatrix` and
`SymFactor::updateSize`, and the degree operation is always both words. Which matters immediately,
because of what the next paragraph says.

**What gets updated is the update, not the front, and that is one statement rather than two.**
C[p] excludes p, so it is the off-diagonal column structure. The eliminator then trims it,
`C[pivot] = C[pivot] - merged` in one compaction pass, so by the time the driver sees it the pivot
and everything merged into it are gone. In supernodal terms the pivot plus its merges is the
FRONT, the fully summed block, and what is left in C[p] is the UPDATE, the rows that front will
update. The degree update set is the update and never the front, and it could not be otherwise: a
vertex in the front is being eliminated now and has no future degree to compute, while a vertex in
the update is exactly one whose neighborhood just changed and that is still a candidate.

In md5 that makes the count exact. **One iteration performs |C[p]| degree updates**, with C[p]
read after the trim, and no subtraction: the eviction loop reads the already-trimmed clique, so
the merged vertices were never in the set to begin with. Over an iteration the set is the union of
the batch's updates. Evictions are the other half of that arithmetic and are NOT deduplicated:
they equal the sum of |C[p]| over the eliminations, always, in md5 and in mmd1 alike. Measured on
a 32x32 grid, mmd1 at delta = 0 does 781 eliminations in 49 iterations, evicting 4946 times and
performing 2872 degree updates; at delta < 0 the last two columns are equal, 5053 and 5053. The
gap between evictions and degree updates IS the batching.

The driver still filters that union with `if (!eliminated[u])`, and what that filter guards is
narrow: a vertex evicted by one pivot and then mass-merged away by a later pivot in the same
iteration. Nothing else can eliminate a touched vertex mid-iteration, since eviction unfiles it and
an unfiled vertex cannot be selected. **We have never observed it fire.** Instrumented, the count
is 0 on all seven examples and 0 on grids of 10, 22, 45 and 64 a side. Whether it is reachable at
all is unverified: it would need a later pivot's clique to absorb the earlier one's, leaving I[u] a
singleton, which is not obviously impossible and is not obviously attainable either. The filter is
cheap and correct, so it stays; it is recorded here as an unverified branch rather than as a
saving.

**The two counters are independent, and the closing line prints both.** They are incremented two
lines apart:

```
numDegreeComputations += degreeUpdated.size();   // += m, the vertices in this pass
++numIterations;                                  // += 1, the pass itself
```

So `degree computations: 30, bucket probes: 8, iterations: 6` on graph3 means six passes performing
thirty vertex degree updates between them, twelve of which are the initial build. It does not mean
thirty passes. Neither counter measures the other, which is what lets the two move apart, and the
worked example below turns on exactly that: iterations go from 4 to 3 while degree computations
stay at 10.

**The reported figure includes the initial build.** `numDegreeComputations` starts at n, since
every vertex's degree is computed once before the loop, so the closing line is n plus the sum of
the per-iteration degree updates. Every layer does the same, so comparisons between them are
unaffected, but on a small graph the number looks larger than the iterations suggest. The initial
pass is also the one cheap case, `A[u].size()` rather than the union every later degree update
costs, which is what the whole ladder above md1 exists to avoid.

`degree computations` keeps its name rather than taking the new vocabulary, because it is what
every figure in `REPORT.md` and in `benchmarks/ordering/README.md` already quotes, and because it
would be off by n under the new one. Read it as counting IN UNITS OF one vertex's degree update.

### How the degree update set is accumulated, here and in genmmd

The two codes form the same union and form it in different places. Recording the difference
because it looks like a candidate whenever MMD1's speed comes up, and it is not one.

**Ours accumulates during the batch.** A second stamp array, `touchedIteration`, alongside the
`mark` and `tag` pair the eliminator already uses:

```
for u in C[pivot]:                          # EVICT, with a stale degree
    buckets.unfile(degrees[u], u)
    if touched_iteration[u] != num_iterations:
        touched_iteration[u] = num_iterations
        touched.append(u)                   # a marker, so O(1) per eviction
```

It has to be a second array rather than the shared one. `mark` and `tag` are scratch inside
`mmd1Neighbors` and inside the eliminator, and a batch calls those many times, so `mark` cannot
survive one elimination while `touchedIteration` must survive the whole iteration. Two lifetimes,
two arrays. Using the iteration number as the stamp is the frugal part: it is monotone and already
maintained, so the array never needs clearing and there is no counter to forget to bump. The union
is then built with no set and no sort, `touched` accumulating in first-eviction order, and the
degree update pass walks it as it stands.

**genmmd accumulates during the pass.** It bumps the same shared tag per elimination, `tag++`
immediately before each `mmdelm`, so it has the identical interference and solves it differently:
it never builds a vertex set at all. It chains the new ELEMENTS instead, in the `list` array it
already owns,

```c
mmdelm(mn, ...) ;
num += qsize[mn] ; list[mn] = ehead ; ehead = mn ;   /* chain the new element */
...
mmdupd(ehead, ...) ;                                  /* one pass, given the elements */
```

and `mmdupd` re-derives the vertices by walking each clique's reach, deduplicating as it goes with
`if (marker[nd] < *tag)`. It also gets two logical marks out of the one array by offsetting the
stamp, `mt = *tag + md0`, so a vertex stamped `mt` fails the `< *tag` test and is excluded without
a second array to say so.

So the trade is an n-sized array and a compare per eviction, against no extra storage and a
re-walk in the pass. **Neither is where the time is**, which is the point of recording it.
Instrumented on production MMD1:

```
                  evictions   degree updates   clique-member visits   visits per update
2D  32x32             4,946            2,872                 76,051                 26
2D 140x140          110,403           51,736              3,199,105                 62
3D  12x12x12         15,877           12,245              1,374,615                112
3D  26x26x26        224,975          166,522             65,506,086                393
```

The third column is the same quantity the amd layers print as "clique-member visits an exact
degree would need", so it is comparable with them. The accumulation is one compare and at most one
push_back per eviction, so 225 thousand operations on the 3D case, against 65.5 million
clique-member visits in the degree updates: a ratio of 291 to 1. Adopting genmmd's scheme would
not even remove that, since the pass would have to re-walk the same reaches to find the vertices.
It would trade one array for nothing measurable.

**What the same table does show is where MMD1's time goes.** The last column is the average reach
walked per degree update, and it runs from 26 to 393. `REPORT.md` records MMD1 at 3.8x to 5.7x the
vendored routine on 3D grids while MMD2 comes back to 0.86x to 1.06x, and this is why: genmmd's
`mmdupd` does not recompute each vertex's reach from scratch. It walks per CLIQUE, uses the
two-scale marking above to exclude what the clique already covers, and computes the degrees of the
vertices sharing that clique together. mmd2's `q2h` shortcut is part of the same idea. The gap is
per-clique against per-vertex degree computation, and not how the update set is formed.

### The smallest graph where batching loses, and what it costs

None of the seven examples shows the wager being lost: on all seven the fill is identical either
way, and only the iteration count moves. That is not an accident of the set, it is a fact about
small graphs. Searching every graph up to six vertices exhaustively, 64 of them on four vertices,
1024 on five and 32768 on six, **not one** has a different fill at delta = 0 than at delta = -1. On
seven vertices 1330 of the 2097152 do, and the fewest edges any of those has is eight. So seven
vertices and eight edges is the true minimum, and one of the minimal witnesses is worth looking at
because it is legible:

```
edges: 0-1 0-5 0-6 1-6    2-3 2-4 3-4    2-5

two triangles, {0,1,6} and {2,3,4}, joined by a path of length two through vertex 5
```

Every vertex except 0 and 2 has degree 2, so the minimum bucket holds 1, 3, 4, 5 and 6. Running the
same code twice, changing nothing but delta:

```
delta = -1, which is md5                      delta = 0
  iteration 0: eliminate 6, merge 1             iteration 0: eliminate 6, merge 1
  iteration 1: eliminate 0                                   eliminate 5
  iteration 2: eliminate 5                                   eliminate 4, merge 3
  iteration 3: eliminate 2, merge 3, 4          iteration 1: eliminate 2
                                                iteration 2: eliminate 0
  4 iterations, 10 deg updates, nnz(L) = 15     3 iterations, 10 deg updates, nnz(L) = 16
  fill = 0                                      fill = 1
```

**Vertex 5 is an articulation point**, and that is the whole story. It is the only vertex joining
the two triangles, it sits at the minimum degree along with four others, and it is adjacent to
neither 6 nor 4. So it passes the only test batching applies, non-adjacency, and the batch takes
it. Eliminating 5 while it still has degree 2 creates the clique {0, 2}, and that fill edge is the
entire difference. md5 never faces the choice: it takes 6 and then 0, which drops 5 to degree 1,
and by the time 5 comes up it is free.

So the shape to look for is not two branches of the elimination forest. Taking 6 from one triangle
and 4 from the other is free in the exact sense given above, since those two cannot interact. It is
**two branches plus the vertex that joins them, all sitting at the same minimum degree.** The batch
has no way to tell the joint apart from the leaves.

**And the degree update count does not move.** Both settings report 10, which is the initial build
of 7 plus 3 actual degree updates, and the 3 is the same 3: md5 updates 0, then 5, then 2, then
none, and the batched run updates 0 and 2, then 0, then none. So this graph gives up a fill edge
for a saved iteration and buys nothing at all on the metric the layer exists to improve, which
makes it the worst case rather than the typical one. Set against graph4, in reported figures:

```
                md5: iterations / deg updates    mmd1: iterations / deg updates    nnz(L)
graph4                 5 / 20                           3 / 15                    26, same
this graph             4 / 10                           3 / 10                    15 against 16
```

graph4 is what batching is for. This graph is what it risks. Neither is representative on its own,
and the grids say the expected value is positive: at delta = 0 they come out ahead on fill and on
degree updates at the same time.

One reason this example is so clean is also a reason not to over-read it. The graph is chordal, so
a perfect elimination ordering exists and md5 finds it, fill = 0. A case where the exact order is
attainable is a case where batching has the most to lose. On a graph with no perfect elimination
ordering, which is every interesting one, both orderings are already leaving fill and the gap
between them is the fraction of a per cent described above.

It is quoted here rather than added to the example set. Adding it would mean an eighth graph in all
thirteen layers in both twins, plus `production.cpp` and `vendored.cpp`, since the examples are one
shared set and the production check diffs every order line of a whole run. amd4's `matrix1` is not
a precedent for doing it in one layer only: amd4 can carry its own example because it has no
production counterpart, and mmd1 and mmd2 both have one.

## mmd2: the extras

mmd1 is multiple elimination and nothing else. genmmd is multiple elimination plus six other
things, none of them an idea and all of them load-bearing for matching the vendored routine.
mmd2 adds them one pass at a time, and this section grows a subsection per pass. The checklist
and its vendored references are in the plan section above. The numbering here is the source
file's, so a subsection and the `mmd2.py` header agree; passes 3 and 4 share one subsection
because both are decided inside the same walk, and the heading says so.

**Where the six land in the file, since four of them never leave one function.** Passes 2, 3, 5
and 6 are entirely inside `mmd2_minimum_degree`, which is why that function is 131 of the 158
changed lines. Only two reach the helpers, and they reach different ones. The prepass touches
`mmd2_neighbors`, which gains an `eliminated` parameter and skips numbered vertices in both
halves of its walk, and `mmd2_eliminate`, which drops them when the prune rewrites an adjacency
list; both are genmmd's `marker[mn] = maxint`, which no `marker[nb] < tag` test can pass.
Outmatched marking touches `mmd2_eliminate` alone, one line clearing the flag for everything the
new clique reaches, which is `mmdelm`'s `bwd[rn] = 0`. `mmd2_file`, `mmd2_unfile`, `mmd2_storage`
and the two display functions are byte-identical to mmd1's. `mmd2_degree` is new, and pass 3 is
the reason: once a candidate can stand for several original vertices, a degree has to count
those rather than entries.

**Evictions continue; the evicted LIST goes.** The unfile per member of `C[pivot]` is still there
and still leaves the reached vertices out of the buckets. What mmd1 also did was accumulate who
was unfiled, stamped by iteration, and walk that list at the refresh. mmd2's refresh walks the
new elements instead and re-derives the same vertices from `C[element]`, deduplicating with the
`filed` flag, so the list and its stamp array have no reader and are not carried. genmmd makes
the same trade, chaining its new elements in `list` and building no vertex set at all.

### Pass 1: the prepass

genmmd numbers every vertex in the degree-1 list before the main loop runs, and never updates
a neighbor:

```c
    int nextmd = head[1];
    while (nextmd > 0) { int mn = nextmd; nextmd = invp[mn];
                         marker[mn] = maxint; invp[mn] = -num; num++; }
```

Three things travel with it, and the middle one is the whole of the pass.

**The floor.** `mmdint` files a degree-0 vertex under degree 1, `if (dg == 0) dg = 1`, so
isolated and degree-1 vertices sit in the same bucket and are numbered together. From that
point the bucket a vertex sits in is max(degree, 1) rather than its degree, and MMD compares and
files by that value rather than by the true one. Our degrees[] now holds it too, which is the
same choice made explicit.

**Numbered is not eliminated.** A prepass vertex is ordered and then skipped. No clique is
formed for it, nothing is pruned, and its neighbors keep degrees that still count it. That
staleness is the point rather than an oversight: the prepass buys its speed by not paying for
the updates. In genmmd the skipping is `marker[mn] = maxint`, which no later `marker[nb] < tag`
test can pass; here it is `eliminated[u]`, which the neighbor query checks and which the prune
loop uses to drop the vertex when it compacts. The compaction drop is the same line in
`mmdelm`, where a numbered vertex fails `marker[nb] < tag` and is left out of the rewritten
adjacency.

**Then bucket 1 is empty.** `head[1] = 0` and the main loop starts at `mdeg = 2`. Ours sets
min_degree to 2 for the same reason. Vertices can return to bucket 1 later, since the filing
convention floors there, and min_degree falls back to meet them.

What it costs, measured on 207 graphs including the seven examples: the permutation is always
valid and the reported nnz(L) always matches an independent symbolic factorization, the fill is
unchanged on 197 and worse on 10, and better on none. On the three grids it is unchanged, 636,
2088 and 4684 against mmd1's. So the prepass is a small concession bought for a real saving,
which is the shape the rest of genmmd has too.

Three of the seven examples exercise it: graph3 numbers vertex 11, graph5 numbers 3 and 2, and
graph6 numbers 5 and 1. The other four have no vertex of degree 0 or 1 to begin with, so the
trace shows no prepass line at all.

### Pass 2: the q2h split

`mmdupd` does not walk a flat list of reached vertices. It walks the CLIQUES this iteration
created, `el = list[el]`, and for each one computes `dg0` once, the weighted size of that clique,
before visiting its members. A member is classified by what it has left BESIDES the new clique:
`mmdelm` stashes `fwd[rn] = nq + 1` where `nq` counts the survivors of the compaction, which in
our split representation is `len(A[u]) + len(I[u]) - 1`. `nq == 1` puts the vertex on the `q2h`
list, anything else on `qxh`.

**What the two names mean, and what is being counted.** They are `mmdupd`'s `q2h` and `qxh`,
two chains threaded through `list[]`, `list[nb] = q2h; q2h = nb`, and the number is the vertex's
TOTAL source count: q2 has exactly two, qx has more. Ours are vectors and are named
`two_source_queue` and `many_source_queue`, for the criterion rather than for the abbreviation. Both lists are built from the members of ONE new clique, so every
vertex on either already belongs to it; the split is on what else it has. Two things are easy to
misread here. The count is not of vertices in a bucket and has nothing to do with bucket length:
a vertex sitting in two enormous cliques is q2h, and a vertex in three tiny ones is qxh. And a
source is not only a clique, since `nq` counts explicit adjacency entries too, the `- 1` removing
the new clique from `I[u]`. A vertex in the new clique with one surviving explicit edge is q2h
exactly as one in the new clique and one other clique is.

**Why the split pays.** Everything a q2h vertex reaches is either inside the new clique, already
counted in `dg0`, or comes from its one other source. So its degree is `dg0` plus what that
source contributes, and the union is never built. The qxh case pays for the full union as
before. Measured on grids, the q2h path takes 36 per cent of the degree updates at 10 by 10, 42 at
16 by 16 and 44 at 22 by 22, so it is not a rare case.

**The qxh loop IS mmd1's refresh, and reading it that way is the quickest route in.** It calls
`mmd2_degree` for the true union and files, which is what mmd1 did for every reached vertex
without exception. Two differences, both small and both consequences of the pass rather than
changes to it. It files at `degree + 1` where the q2h site files at `degree - len(super_members[u])
+ 1`, because `mmd2_degree` returns the EXTERNAL weighted degree while the q2h walk starts from
`dg0`, which carries `u`'s own members and has to take them back out; the two expressions produce
the same quantity by different routes, and that asymmetry is where the entry-5 defect hid. And it
opens with two guards mmd1 has no need of, `eliminated[u]` and `outmatched[u]`, because the q2h
walk over the same clique ran a moment earlier and may have merged `u` away or withheld it.

**`refreshed_vertices` is a byproduct here, not the work list.** In mmd1 it was the thing the
refresh walked. In mmd2 the work happens inline as each vertex is classified and filed, and the
list is accumulated only to count the degree updates, to lower `min_degree`, and to print. Nothing
iterates it to do anything.

**Two mark levels, and the bug that taught me why.** My first version marked the other source's
members with the clique's tag, which made a second q2h vertex in the same clique skip them as
already counted and report a degree too small. `mmdupd` avoids this with two levels: members of
the new clique carry `mt`, above every tag used in the iteration, while each vertex gets a fresh
`(*tag)++` for its own walk. Ours does the same with `clique_tag` and `vertex_tag`, testing both.

**The check on this pass, and the check that was wrong.** Every degree the shortcut produces must
equal the full union, and it does, on 307 graphs. My first attempt at that check called
`md2_neighbors` to get the true value, which advances the tag and overwrites the mark array, so
it destroyed the state it was inspecting and reported failures that were its own doing. The
verification has to compute the true neighbors independently, with plain sets and no marks. Worth
remembering for the passes to come: an instrument that shares state with the thing it measures
is not an instrument.

What moves is the filing order, since the degree update is now clique by clique with q2h first, and
filing order decides what a bucket holds. So the permutation changes even though every degree is
identical: on the grids, nnz(L) goes 636 to 633, 2088 to 2101, and 4684 to 4684 against mmd1.
A vertex reached by two pivots in the same iteration is still updated once, skipped on the second
visit by the `filed` flag, which is `if (bwd[en] != 0) goto n2200` there.

**graph1 end to end, which is small enough to hold in one view and still fires both mechanisms.**
The 4-cycle `0-1 1-2 2-3 3-0`. Every degree is 2, so the prepass finds nothing. The batch takes 3
and 1, opposite corners and so non-adjacent, and each forms a clique on `{0, 2}`:

```
0: A = {}  I = {c3, c1}          c3 = {0, 2}
2: A = {}  I = {c3, c1}          c1 = {0, 2}
```

The refresh then walks the two new cliques rather than the reached vertices, and what each
contributes is:

```
clique c3: members [0, 2]   dg0 = 2   q2h = [0, 2]   qxh = []
clique c1: members [0, 2]   dg0 = 2   q2h = []       qxh = []
```

`dg0` is c3's weight, 2, computed once for the whole clique. Classifying its members,
`len(A[u]) + len(I[u]) - 1` is `0 + 2 - 1 = 1` for both 0 and 2, so both have exactly one other
source, c1, and both are q2h. Nothing reaches qxh at all.

Take u = 0. The walk starts at `dg0 = 2`, which is already the whole of c3, and then looks only at
c1. Its members are 0 itself, and 2, which is already marked as belonging to c3 and so already
inside `dg0`. Nothing is added. **The union of c3 and c1 is never formed**, and that is the entire
pass in one line.

That same marked vertex is what triggers the merge of the next pass. 2 lies in both cliques, so it
sees at least what 0 sees, and it is q2h too, so its one other source is c3: identical reach, and
0 absorbs it. So `pair merges: 1`, and 2 leaves without ever being refreshed. c1's own walk then
finds nothing left to do, which is why the whole iteration reports `degree updates: 1`. Vertex 0
files at `dg0 - len(super_members[0]) + 1 = 2 - 2 + 1 = 1`, the weight being 2 AFTER the merge,
which is the subtraction order that was wrong here until 2026-08-07.

**The two grid tables from the mmd1 section, with mmd2 alongside.** Those fix delta at its two
ends and sweep the size; this adds a third column for mmd2 with all six passes in, so the three
groups read as md5, then batching, then the vendored extras. The mmd2 figures elsewhere in this
section are snapshots taken as each pass landed and do not match these, which are the finished
layer. Grids have no vertex of degree 0 or 1, so the prepass numbers nothing here and every
difference in the third group comes from passes 2 to 6.

```
            md5 (delta = -1)         mmd1 (delta = 0)       mmd2, all six passes
           nnz(L)  iters  degupd    nnz(L)  iters  degupd    nnz(L)  iters  degupd
10x10         657     78     376       636     19     290       618     20     192
16x16        2195    197    1116      2088     31     751      2052     31     474
22x22        4773    367    2206      4684     36    1375      4553     39     795
30x30       10436    679    4327     10757     40    2533     10059     44    1398
40x40       22495   1207    8132     21614     56    4498     20536     58    2421
```

```
             md5 (delta = -1)            mmd1 (delta = 0)          mmd2, all six passes
         elims  sum|C|   fill probes   elims  sum|C|   fill probes   elims  sum|C|   fill probes
10x10       78     376    377     89      82     395    356     31      81     359    338     29
16x16      197    1116   1459    217     201    1108   1352     55     204    1034   1316     46
22x22      367    2206   3365    399     373    2198   3276     68     376    2028   3145     61
30x30      679    4327   7796    723     687    4358   8117     90     692    3924   7419     74
40x40     1207    8132  17775   1273    1215    8010  16894    121    1224    7131  15816    104
```

**Degree updates fall again, and by more than batching bought.** mmd1 took 8132 to 4498 at 40 a
side; mmd2 takes 4498 to 2421, another 46 per cent, on top. Two mechanisms do it and the split is
only one: the q2h path answers most updates without a union, and outmatched marking withholds
vertices from both lists so they are not refreshed at all. `sum |C|` falls with them, 8010 to
7131, which is the merges making the cliques smaller.

**The fill improves on all five grids**, 636 to 618, 2088 to 2052, 4684 to 4553, 10757 to 10059,
21614 to 20536. That is worth reading carefully rather than as a win: mmd2's job is fidelity, and
the fill goes where the vendored algorithm puts it. On the 307-graph set the same passes came out
better on 2, same on 289 and worse on 16, so five grids agreeing is a property of grids and not a
general improvement.

**Iterations go slightly UP**, 19 to 20 and 56 to 58, against every other column falling. The
batch is drained by degree, and outmatched vertices are held out of the buckets, so a bucket
empties sooner and the batch ends sooner. Fewer candidates per pass, more passes, and still far
fewer degree updates in total.

**Nothing about this pass is genmmd's.** It lands in mmd2 because that is the layer where the
vendored extras arrive, but the shortcut needs only cliques, weights and a source count, all of
which exist from md5 on. The same is true of the other five, and the closing subsection sets them
out together with what each would cost at md5.

### Passes 3 and 4: the pairwise merge, and outmatched marking

Both live in the same branch of the two-source walk, reached when a member of the one other
source is ALSO a member of the new clique:

```c
    else if(bwd[nd]==0){
        if(fwd[nd]==2){qsize[en]+=qsize[nd];qsize[nd]=0;marker[nd]=maxint;
                       fwd[nd]=-en;bwd[nd]=-maxint;}
        else if(bwd[nd]==0)bwd[nd]=-maxint;}
```

**The merge.** If nd is two-source too, its only other source is that same clique, so en and nd
reach exactly the same vertices: indistinguishable, and en absorbs nd. This is the first merge in the
whole sequence that folds a vertex into a LIVE one. Every earlier merge went into the pivot,
which died in the same call, so a supervariable never survived in the buckets. Now one does, and
a degree has to count original vertices rather than entries. It is also, concretely, what makes
MMD's supervariables coarser than md3's, whose test only ever compares a vertex against the
pivot.

There is still no weight array. The count is `len(super_members[v])`, which is O(1), and the
invariant `weight == len(super_members)` was checked across 307 graphs before the array was
removed again. The rule has not changed since md3: an array is kept when it stops being
derivable, not when the quantity it holds starts varying. MMD keeps qsize because its members are
a chain rather than a list, and that is the condition under which ours will need one too.

**Outmatched, and which of the two the word is about.** If nd is not two-source it has sources
besides these two, so its reach CONTAINS en's and it can never be the minimum before en. The word
is the competitive one, outclassed: en outmatches nd, and the flag lands on nd, the loser. The
AMD paper states the general form, a variable being outmatched when another's adjacency is
contained in its own, and files it under incomplete degree update rather than under pruning,
because the saving is the degree update nd never pays for. The direction reads backwards at
first, since the flag marks the vertex with MORE reach; under minimum degree more is worse.

Note that our test is not the general containment predicate. en and nd already share the new
clique, and nd being many-source implies the containment here, so the split hands it over for
free rather than making us compare two adjacency sets.

MMD withdraws it from the degree lists rather than refiling it, `bwd[nd] = -maxint`. It is not
merged and not eliminated, just held out until something reaches it again, at which point
`mmdelm` restores it with `bwd[rn] = 0`. Ours is an `outmatched` flag, cleared in mmd2_eliminate
for every vertex the new clique reaches. Both are withheld from the two-source and many-source
queues while the flag is set.

**Weights everywhere.** dg0 becomes a weighted sum, the two-source walk starts at
`dg0 - len(super_members[u])` rather than `dg0 - 1`, every count adds `len(super_members[v])`,
the many-source path goes through mmd2_degree, and the nnz(L) accounting sums the same over the
live members of C[pivot].

**graph1 again, where only one of the two fires.** The same run as in the previous subsection.
Vertex 0 is being refreshed from clique c3, `dg0 = 2`, and the walk over its one other source c1
reaches vertex 2, which is marked as belonging to c3 as well. That is the branch. Vertex 2's own
count is `len(A[2]) + len(I[2]) - 1 = 0 + 2 - 1 = 1`, so it is two-source, so this is the equality
case and 0 absorbs it:

```
[probe] MERGE: u=0 absorbs v=2
pair merges: 1, outmatched: 0
```

**Nothing is outmatched on graph1, and the reason is structural rather than luck.** Outmatching
needs a vertex reached by the walk with THREE or more sources, and the 4-cycle after two
eliminations has only two cliques in it, both on `{0, 2}`, so no vertex can have a third. It is
the smallest graph that fills, and it is too small for a containment that is strict. The
`outmatched: 0` in its summary is a fact about the graph, not about the pass being dormant.

**graph3 has one of each, which is the contrast worth seeing.** By iteration 3 the live vertices
are 1, 5, 6 and 8 over three cliques:

```
 1: {6} {c3}          c3: {5, 1, 8}
 5: {}  {c3 c9}       c7: {6, 8}
 6: {1} {c7 c9}       c9: {5, 8, 6}
 8: {}  {c7 c3 c9}
```

Vertex 5 has two sources, c3 and c9; vertex 8 has three, c7, c3 and c9. Refreshing 5 from c9, the
walk over its other source c3 reaches 8, which is in c9 too. So 8 sees at least what 5 sees, but
it is many-source, so the containment is strict:

```
[probe] OUTMATCH: u=5 outmatches v=8; v other sources = 2
[probe] MERGE: u=1 absorbs v=5
[probe] MERGE: u=1 absorbs v=8
pair merges: 2, outmatched: 1
```

The two merges after it are worth noting: 8 is withheld, not discarded, and it comes back as soon
as a clique reaches it, at which point 1 absorbs both. Withholding delays a degree update; it
does not decide the ordering by itself.

**What it does.** On the same 307 graphs: 245 pair merges and 126 outmatched markings, all
permutations valid, all reported nnz(L) correct against an independent symbolic factorization,
and the fill against mmd1 better on 2, same on 289, worse on 16. On grids the mechanisms fire
steadily: 18 merges and 55 outmatched at 10 by 10, 56 and 135 at 16 by 16, 111 and 276 at 22 by
22, with nnz(L) 631, 2105 and 4783 against mmd1's 636, 2088 and 4684.

That last column is worth reading honestly. These two mechanisms do not improve the ordering
here; they make it match the vendored routine. mmd2 is fidelity, and the fill goes where the
vendored algorithm puts it.

### Pass 5: the filing convention

`mmdupd` does not file a vertex under its degree. It files under `dg = dg - qsize[en] + 1`,
floored at 1, where `dg` was the weighted reach INCLUDING en's own members. So the bucket index
is the external degree plus one, and the floor catches a vertex that reaches nothing outside
itself.

`mmdint` meanwhile files at the plain degree, with only the zero case lifted to 1. So MMD runs on
two scales: the initial buckets hold degrees, every refiled bucket holds degree + 1. That is
genuine rather than a misreading, and it tilts the pivot choice slightly against degree updated
vertices, which sit a bucket higher than an untouched vertex of the same reach. From here
degrees[] holds the FILED value, which is what the picker compares and what min_degree tracks;
the nnz(L) accounting is unaffected, since it sums weights over the live members of C[pivot].

**Three places write a bucket index, and only the last two are this pass.**

```
initial       degrees[u] = max(degrees[u], 1)                          pass 1's floor
two-source    degrees[u] = max(degree - len(super_members[u]) + 1, 1)
many-source   degrees[u] = max(degree + 1, 1)
```

The two refresh sites look different and produce the same quantity. `mmd2_degree` already returns
the EXTERNAL weighted degree, so the many-source site adds one and is done; the two-source walk
starts from `dg0`, which carries u's own members, so it has to take them back out first. Same
`dg - qsize[en] + 1` by two routes, and the asymmetry is a consequence of where each path gets
its degree rather than a difference in the convention.

**The two-source site is where the one fill-visible defect lived.** The subtraction has to use
the weight AFTER the walk, because that walk can merge a vertex into u and genmmd does
`qsize[en] += qsize[nd]` inside it. Subtracting first files a supervariable one bucket too high
per merged vertex, so it is never picked as early as its size has earned. That was wrong here
until 2026-08-07 and was found only by aligning mmd3 against genmmd; it is entry 5 in that file's
ledger, and the comment above the line carries the full account.

**Two things this pass does not touch**, worth knowing so they are not searched for. The nnz(L)
accounting reads weights over the live members of C[pivot] and never a bucket index. And
min_degree tracks filed values on both sides, at the initial build and after each refresh, so it
stays consistent without an adjustment for the two scales.

### Pass 6: the counters

Two small things in genmmd's main loop, no mechanism between them: one accumulator and one guard,
neither touching a degree or a permutation.

`ncsub`, accumulated per pivot as `*ncsub += mdeg + qsize[mn] - 2`, is the statistic genmmd
returns alongside the permutation, an estimate of the subscript storage the factor will need.
Ours accumulates the same and prints it with the other counters:

```
ncsub += degree + len(super_members[pivot]) - 2
```

one line in the elimination block, just after the eviction loop.

The early termination, `if((num+qsize[mn])>neqns)goto n1000`, is checked after the pivot is
numbered and before it is eliminated. When the last supervariable is reached there is nothing
left to update, so genmmd skips the elimination entirely and goes to the numbering. Ours is:

```
if num_eliminated_vertices >= n:
    break
```

at the foot of the batch loop, immediately above the `delta < 0` break. The shape differs from
genmmd's and the effect does not: it tests BEFORE eliminating, using `num + qsize[mn]` because
its counter has not yet taken the pivot in, while ours breaks AFTER, its counter already updated.
Worth checking rather than taking on trust, since that is the shape an off-by-one hides in.

**`ncsub` is the one item on the checklist with no oracle behind it, and that is a gap rather
than an impossibility.** It has not been compared against the vendored number, because
`mmd_order` declares it as a local, passes `&ncsub` to `genmmd`, and then drops it. One temporary
line in the wrapper prints it, which is the same move `tools/hook_amd.py` makes on the AMD side,
and the vendored sources are ours to instrument whenever a check needs it. So the honest status
is unchecked, not uncheckable, and it stays on the list until someone runs it.

### Where mmd2 got to, and what is still unchecked

Everything on the six-item list is implemented, in both twins, with the traces agreeing byte for
byte on the seven examples and on forty random graphs. Two things are deliberately absent, both
tag arithmetic that a monotone tag makes unnecessary: MMD's `maxint` overflow reset, and
`md0 = mdeg + delta` feeding the marker window `mt`.

The changes landed almost entirely in the driver, 131 of 158 changed lines. The buckets, the
display and the storage function are untouched, as they have been since md4 and md5. Three
functions outside the driver did change, and each marks a place where an extra was not confined
to policy: mmd2_neighbors skips numbered vertices, which is what the prepass requires;
mmd2_eliminate clears the outmatched flag for everything the new clique reaches, which is
mmdelm's `bwd[rn] = 0` and the first change to the eliminator since md2; and mmd2_degree is new,
because a candidate can now stand for several original vertices.

Against the vendored routine, which is the acceptance test the plan section sets:

```
                       identical order   identical nnz(L)
seven examples               2 of 7            7 of 7
60 random graphs            23 of 60          59 of 60
```

The one disagreement was worth chasing rather than accepting, since a fill difference can hide a
missing feature. It does not here. On that graph, 14 vertices, the two runs agree through the
prepass and the first three pivots, and then diverge on a bucket holding six vertices at degree
4: ours pops the head, 4, and the vendored routine pops 2, which sits second in our list. Same
degrees, different insertion history. Orders differ far more often than fill does for the same
reason, which is what the tie-break section predicts.

For contrast, mmd1 matches the vendored nnz(L) on 54 of the same 60, so the six passes moved the
agreement from 54 to 59 rather than from nothing. That is the honest size of what they buy: MMD's
extra machinery changes the ordering in small ways, and most of the fill was already there.

On grids, mmd2 against mmd1: 624 against 636 at 10 by 10, 2078 against 2088 at 16 by 16, 4681
against 4684 at 22 by 22, with the mechanisms firing steadily, 108 pair merges and 285 outmatched
markings on the largest.

What is still unchecked about mmd2 is listed with everything else, in the section on open items
near the end.

**Where these six could have gone instead, which is not where they went.** They land in mmd2
because it is the layer that starts matching a vendored routine, not because anything in them
needs mmd1 underneath. Taken one at a time against what each requires:

```
prepass              buckets                                available from md5
two-source split     cliques, and dg0 per clique            available from md5
pairwise merge       weights, since a live vertex absorbs   available from md5
outmatched marking   a way to withhold from a bucket        available from md5
filing convention    nothing, it is a convention            available anywhere
counters             nothing                                available anywhere
```

So the whole of mmd2 could sit at md5, and the split in particular, which is the one with a real
saving behind it: 36 to 44 per cent of degree updates answered by an addition rather than a union,
and nothing in that argument mentions batching.

Two things keep it from being a free observation. Each pass changes md5's output, and md5 is the
reference the layers above are read against, so folding six mechanisms into it would move the
ground everything else stands on; the ladder's one-mechanism-per-layer rule would give each its
own layer in any case. And the split is not quite as cheap below mmd1: its classifier comes from
`mmdelm`'s compaction, which counts survivors while it is already rewriting the list, so at md5
the count would be computed on demand as `len(A[u]) + len(I[u]) - 1`. Still cheap, but the pass
whose whole case is that it costs one test would be paying for the test.

We are not moving any of it. The point is that the position of these six in the ladder records
when we started matching genmmd, not when the ideas became available.

## amd1: the approximate degree

md5 leaves one expensive thing standing. After each elimination, every reached vertex needs a new
degree, and getting it means uniting the members of every clique in its list and counting the
result. mmd1 attacks the frequency of that union. amd1 attacks its price.

bound(u) is derived in the md2 subsection above, which is where the union it decomposes first
exists, and it becomes usable at md4, where the degree update set narrows to exactly the vertices it
covers. In brief: reach(u) is decomposed against the new clique, the union is replaced by a sum,
and what comes out is bound(u), with degree(u) <= bound(u). amd1 picks the u minimizing bound(u)
and changes nothing else.

**Why that is cheap, which is the whole point and is easy to miss.** |C[c] - C[p]| depends on the
clique c alone and not on the vertex u. So it is computed once per clique and then read by every
vertex whose incidence list contains c. The exact degree costs, per vertex, a walk over the members
of all its cliques; bound(u) costs, per vertex, one addition per clique. Both are counted in the
trace, and the gap widens with the size of the cliques, which is to say with the fill, which is to
say exactly where it matters:

```
              exact would visit   the bound read
grid 10x10            958              370
grid 16x16           3980             1503
grid 22x22           8426             2955
```

**What is given up, and it is a different kind of loss from mmd1's.** Every layer up to here picks
a true minimum-degree vertex, and md1 through md5 differ only in how they find one while mmd1
differs only in how ties fall. amd1 can pick the wrong vertex outright, because an overcounted
bound can hide the true minimum. It is the first layer whose heuristic changes rather than its
implementation.

That makes the trace do something the other layers do not need: it prints the exact degree beside
the bound at every iteration, and the closing lines count how often the two disagreed. The exact
value is computed with amd1_exact_degree, which is md5's degree update kept alive for no other
purpose. It is
instrumentation, and a real engine would not carry it.

On small graphs the bound is nearly always tight, 12 loose out of 1748 checks over the seven
examples and 200 random graphs, and the ordering comes out identical to md5's on all 207. On grids
it is loose far more often, 798 of 2219 checks at 22 by 22, and the fill still lands within half a
percent of md5, 4762 against 4773. Both facts are the AMD paper's claim in miniature: the bound is
wrong often and wrong by little.

**What amd1 does not have.** Everything else. Aggressive absorption and hash supervariable
detection ride along in the same sweep in the vendored routine and would cost almost nothing here,
which is exactly why they are held back: this file is one idea, as mmd1 was. They are amd2's pass
1.

## amd2: the extras

amd1 is the bound and nothing else. amd_1 and amd_2 are the bound plus seven other things, and
amd2 adds them one pass at a time. The checklist and its vendored references are in the plan
section above.

### Passes 1 and 2: aggressive absorption and hash supervariable detection

The two mechanisms that ride along with the bound. Neither is about the degree, and both are cheap
only because the bound's own work has already been done.

**Aggressive absorption.** |C[c] - C[p]| has just been computed for every clique the new one
touched. A zero means C[c] lies entirely inside C[p], so that clique is dead and can be absorbed at
once. Ordinary absorption, which every layer since md2 has, kills only the cliques the pivot
touched. This kills cliques that any reached vertex touched, and the quantity it tests was needed
anyway.

**Hash supervariable detection.** Mass elimination merges a vertex into the pivot, and md3 through
amd1 can only see pairs where one of them is absorbable into the pivot. Two vertices can be
indistinguishable from each other with neither being absorbable, and no pivot-relative test will
ever find them. AMD hashes (A[u], I[u]), compares only within a hash bucket, and merges on an
exact
match. The hash is a filter and never the decision, so a collision costs a comparison and not a
wrong answer. mmd2 reaches the same population through mmdupd's q2h list, so the goal is shared
and the mechanism is not; the section on detecting supervariables against each other covers both.

**One consequence, and it is the same one mmd2's pairwise merge had.** Hash detection folds a
vertex into a LIVE one, so from here a member of the new clique can stand for several original
vertices. The nnz(L) accounting had used len(C[pivot]) for the external degree, which is correct
exactly while no live vertex stands for more than one original. It is now a weighted sum. Without
that fix the reported fill was too low on 110 of 207 graphs, and every one of them was a graph
where a hash merge fired, which is what made it findable: the wrong number was wrong in a pattern.

The mechanisms fire steadily. Over the seven examples and 200 random graphs, 193 aggressive
absorptions and 198 hash merges. On grids the hash does most of the work, 112 merges against 1
aggressive absorption at 22 by 22, which says that on a regular mesh the cliques rarely nest but
the boundary vertices are often twins.

**SUPERSEDED 2026-08-08: they cost nothing, and the table below was measuring a defect.** `amd2`
filed a supervariable one bucket too high per vertex a hash merge absorbed, because the bound
subtracts the vertex's own weight before the merge that grows it, where `AMD_2` subtracts it in
the pass that restores the degree lists, after supervariable detection. Corrected, amd2 beats
amd1 at every grid size by 1 to 3 percent rather than losing to it: 11900 against 12074 at 32 a
side, 199386 against 201856 at 100, 444191 against 455472 at 140. So the two extras are a net
gain on fill and the coarser-supervariable cost this section describes was the filing, not the
supervariables. Found by aligning `amd3` against the vendored routine, where the same timing is
ledger entry 4; `docs/DESIGN_DECISIONS.md` (2026-08-08) carries it. The table below is kept as
the record of the run that produced it.

**And superseded a second time, 2026-08-09**, by ledger entry 8: the hash key fix moves
`amd2`'s tie-break and so its fill again. The figures here are two corrections behind.

**What these two cost in fill, which is not nothing.** amd2 against amd1 on grids, re-measured
2026-08-03 after the two defects below were fixed, and every figure checked against a symbolic
factorization of the emitted permutation:

```
              amd1     amd2
 10x10         657      659
 16x16        2175     2181
 22x22        4762     4753
 32x32       12074    12364
 64x64       67950    68822
100x100     201856   212496
```

and on the 207 small graphs the same on 206 and worse on 1. Coarser supervariables mean fewer,
larger pivots, so the ordering is cheaper to produce and slightly worse. That is the trade the
vendored routine also makes, and the comparison to make is against the vendored routine rather
than against amd1.

Two things are worth reading off that table beyond the trade itself. The gap is not monotone in
size, amd2 winning at 22 by 22 and losing on either side of it, which is the two-sided noise a
tie-break difference produces rather than a systematic cost. And the last two rows are exactly
production's AMD1 and AMD2 in `benchmarks/ordering/README.md`, digit for digit, which is the
strongest evidence the port is faithful that either file carries: the prototype and the engine
agree on the permutation and on what it fills, at a size the seven examples cannot reach.

Against the vendored amd_order on the seven examples, both layers already match nnz(L) 7 of 7 and
neither matches the permutation on any, which is the tie-break story the mmd2 section tells at
more length. So these two do not move the acceptance test on their own; they move the mechanism
count. The figures with all seven passes in are at the end of this section.

### Two defects the permutation check found, both inherited from amd1

Found 2026-08-03, and they are one finding rather than two. Splitting the old amd2 into amd2 and
amd4 made amd2 comparable to production's `Amd2` by PERMUTATION rather than by fill, and that
stronger check immediately went red on grids of 10 by 10 and larger while staying green on the
seven examples. Both defects are lines that were correct in amd1 and became wrong the moment amd2
added a merge into a LIVE vertex. amd4 already had both fixed, and so did production, so amd2 was
the only file with them.

**The bound's live-vertex cap, taken from the wrong counter.** The first cap is
`num_left - weight(u)`, and amd2 derived `num_left` as `n - num_eliminated`, which is amd1's line.
`num_eliminated` counts what has left the SELECTION, and a hash merge folds `v` into a live `u`, so
`v` stops being selectable while every vertex it stands for is still live inside `u`. The count
therefore overstates what has gone, the cap comes out too tight, and the bound falls below the true
degree, which is the one thing a bound must never do:

```
                bound below exact
grid  10x10          22
grid  16x16          45
grid  25x25          82
```

zero in all three after the fix. So this was never a tie-break difference, and the instrumentation
could not see it: the layer counts how often the bound is LOOSE and never how often it is wrong in
the other direction, which is the same blind spot the pass 5 bug had.

**The vendored routine is the oracle, and it settles this without appeal to amd4.** `AMD_2`
advances `nel` at exactly two sites inside its main loop, `nel += nvpiv` at the pivot and
`nel += nvi` at mass elimination inside scan 2, both of which are original vertices genuinely
leaving. Its supervariable absorption moves the weight and leaves `nel` alone:

```c
Nv [i] += Nv [j] ;
Nv [j] = 0 ;
Elen [j] = EMPTY ;
```

so `nleft = n - nel` there counts live originals. That is production's `numLive` and amd4's, and
it is what amd2 now carries.

**The fill accounting, taking an unweighted count.** `external_degree = len(C[pivot])` is amd1's
line and is correct there, since without live merges every member of the new clique stands for
exactly one original vertex and none of them is dead. With hash detection a member can stand for
several, and a merged one stands for none while still lying in the list. Measured against a
symbolic factorization of the emitted permutation, the reported figure was too low by 12 percent at
10 by 10 and 30 percent at 32 by 32; it is exact at every size now.

This one is worth separating from the first, because the two have opposite reach. The cap is an
ordering defect and changes the permutation. The accounting is instrumentation and changes only
the number printed beside it, so a reader comparing amd2's fill against amd1's was comparing a
wrong number with a right one, and the apparent conclusion, that hash detection buys a large fill
saving that grows with size, was entirely the defect. The corrected table is in the subsection
above and says the opposite, which is what the vendored routine's own trade predicts.

**And the harness is the other half of the finding, again.** Both defects leave the seven examples
byte for byte identical, before and after, so no check in `make test` moved in either direction.
Seven graphs of at most twelve vertices cannot exercise a mechanism that needs enough structure to
fire, and the first defect fires 22 times on a 10 by 10 grid.

**So the harness was changed rather than only the code, in three places.**

- **The instrument was half-blind, and is not now.** The amd layers counted how often the bound was
  LOOSE and never how often it went the other way, which is the one direction that is a defect
  rather than a measurement. All three now print `bound below exact N times, which must be zero`
  beside the looseness count, in both twins. That line reads 22 on a 10 by 10 grid with the cap
  defect in place and 0 without it, so the instrument that missed this now reports it in the
  ordinary trace.
- **`make test` compares on grids as well as on the seven examples**, prototype against production
  for every layer in `PORTED`, at sides 10 and 20.
- **And the Python twins gained the grid mode the C++ ones already had**, so that comparison runs
  between the twins too. Without it the twin check could only ever see the examples, which is a
  gap the mutation test below exhibits directly.

Confirmed to catch what it was built for, by putting the cap defect back into the C++ twin alone
and leaving the Python correct:

```
amd2   py and cpp agree                          the seven examples cannot see it
amd2   DIFFER on grid 10                         the twin grids can
amd2   DIFFER on grid 20
amd2   prototype and production agree            nor can the example check
amd2   DIFFER on grid 10                         the production grids can
amd2   DIFFER on grid 20                         make test exits 2
```

The first line is the gap stated in one line: the Python twin was correct and the C++ was not, and
the seven examples could not tell them apart. Four checks catch it now where none did.

The checks went from 17 to 35 and the run from 0.5 s to 7.8 s with everything already built,
nearly all of the increase being Python at side 20 and nearly all of that the amd layers'
exact-degree instrumentation, which computes per updated vertex the very union the bound exists
to avoid. The Makefile's `GRID_SIDES` comment carries the figures and says what to edit if the
wait ever grates.

### Pass 3: the two-pass degree update

The checklist called this "so every survivor sees the same final degme rather than a shrinking
one", and reading `Amd.cpp` showed that description had it backwards. Worth recording, because the
correction is the interesting part.

**We were already two-pass, in both senses the phrase could mean.** Scan 1 computes every
|C[c] - C[p]| before any vertex's bound is formed, and amd1 already did that.

**The second half of this paragraph was wrong too, corrected 2026-08-08.** It said the shrinking
degme is vendored behavior, so a survivor handled early sees a larger degme than one handled late,
and that our front-loaded mass elimination avoids a loss the vendored file flags in its own
comment. `Amd.cpp` does decrement `degme` inside scan 2, but it never READS it there: the term
enters a survivor's degree only in the later pass that restores the degree lists,
`deg = Degree[i] + degme - nvi`, by which point degme is final. So every survivor sees the same
number in both codes and there was nothing to avoid. The loss the vendored comment flags is a
different one, about mass elimination merging non-adjacent vertices, and it is why AMD's own
nnz(L) is an upper bound.

**What is actually different is how the per-clique quantity is obtained, and it is a cost
difference rather than an output difference.** amd1 computes |C[c] - C[p]| by walking the members
of C[c] and skipping those inside C[p]. Scan 1 never looks at a clique's members at all:

```
|C[c] - C[p]| = |C[c]| - sum of weight(u) over u in C[c] & C[p]
```

The first term comes from a maintained clique degree; the second is accumulated by walking the
INCIDENCE lists of the new clique's members, since c is in I[u] exactly when u is in C[c]. So the
scan pays sum |I[u]| over u in C[p], where amd1 paid sum |C[c]| over the touched cliques. The
vendored code is the same two lines, `we = Degree[e] + wnvi` on first sighting and `we -= nvi`
after.

Measured on grids, with amd1's cost shown for comparison:

```
            member visits amd1 paid   incidence reads scan 1 pays
grid 10x10           1626                        351
grid 16x16           7006                       1427
grid 22x22          14556                       2810
```

**The maintained degree is exact, not an estimate, and that is not obvious.** A stale value would
still leave the bound a bound, since it would only overcount, so the vendored code could afford to
be sloppy here. Ours is not, and three invariants are the reason. A live clique never holds an
eliminated vertex, because eliminating v absorbs every clique in I[v]. Mass elimination only
removes a vertex whose I[u] is exactly {pivot}, so no other clique is touched. And a hash merge
folds v into u where I[u] == I[v], so every clique holding v also holds u, and the weight moves
from one member to another without the clique's weighted size changing. Checked against the
definition with plain sets on 310 graphs including grids, and it holds at every iteration.

The output is unchanged: same orders, same fill, same looseness. Pass 3 buys cost and nothing else,
which is what it should buy.

### Pass 4: dense rows and columns

One dense row touches nearly everything, so it inflates the degree of nearly everything, and the
ordering spends its effort avoiding a vertex that cannot be avoided. AMD takes such rows out before
it starts and puts them last. The threshold is

```
dense = min( n, max( 16, alpha * sqrt(n) ) )       alpha < 0 means n - 2
```

with alpha defaulting to 10, and it is applied to the INITIAL degrees only. A vertex that becomes
dense later is not caught, which is deliberate: this is a property of the matrix, not of the
elimination.

**First, a measurement, because it decides how this can be tested at all.** At the default alpha
the threshold is at least 16, and nothing in our test set comes close: the grids top out at degree
4, the random graphs at fourteen vertices cannot exceed 13, and the seven examples are smaller
still. Across the seven examples, 200 random graphs and three grids, the mechanism fires zero
times and every order and every fill is exactly what pass 3 produced. So alpha is exposed as a
parameter, as `mmd1` exposes delta, and the pass is exercised with stars and with dense random
graphs where it does fire.

**Removing a vertex means zero weight, not deletion, in the vendored code.** `Amd.cpp` sets
`Nv[i] = 0` and `Pe[i] = EMPTY`, so a dense vertex stays in everyone's lists while contributing
nothing to any degree. In our representation a weight of zero and an absence are indistinguishable
to every count the file makes, so we clear its own lists, empty its `super_members`, mark it
eliminated and purge it from the rest. The one visible difference is that our degrees are
recomputed after the purge while `Amd.cpp` files by the degree it computed before, a transient the
first degree update corrects in either case.

**And this is where nnz(L) stops being a count.** A dense row was taken out but it still sits below
every column of L, so `Amd.cpp` uses `r = degme + ndense` for every pivot and then counts the
trailing block of dense rows as completely full. Both are worst cases. We now do the same, so with
`ndense > 0` the reported nnz(L) is an upper bound rather than a count: on two disjoint stars of
thirty vertices it reports 177 against a true 118, since neither centre is actually nonzero in the
other's columns. Checked on 120 random graphs with the threshold forced low, the reported figure
was never below the truth.

That is the first place this experiment gives up an exact number, and it is given up on purpose.
The alternative was to keep the weights and withhold the dense rows from selection only, which
would have preserved the count and kept the label while doing much less than AMD does. Aligning
with the vendored routine matters more here than keeping a number we can check.

### Pass 5: amd_aat and amd_preprocess, the input path

Every layer up to here takes a graph: symmetric, no duplicates, no diagonal. `amd_order` takes a
matrix, which is none of those things, and this pass is the two routines that bridge them.

**It is input conditioning, not a feature for unsymmetric matrices.** Four things happen at once:
the pattern is symmetrized, the diagonal is dropped, duplicates are removed, and unsorted columns
are tolerated. Only the first is about symmetry, and the other three apply to a perfectly symmetric
matrix.

**Dropping the diagonal means dropping it from the GRAPH, not from the matrix.** The matrix keeps
its diagonal, which the factorization needs whether or not an entry is zero. The graph is a
different object: an edge i-j says that eliminating one forces fill between the other and its
neighbors, and a self loop says nothing, since i is not its own neighbor. Left in, it would put u
inside reach(u) and make every degree one too large. So amd4_aat builds A[u] with u excluded, which
is the `i == j` case in the code. It is a property of the conversion, not a modification of
anything.

**For Oblio almost none of this pass is needed, and that is worth stating plainly.** Oblio's rule
is that a matrix is valid on construction: sorted, no duplicates, a diagonal present even where it
is zero, and both triangles stored. Those invariants kill the symmetrization, the deduplication and
the tolerance of unsorted columns outright. What survives is skipping i == j while building
adjacency, one line inside whatever turns a SparseMatrix into a graph, not a routine worth porting.

The asymmetry is deliberate on both sides. AMD is a library that must take whatever a caller hands
it, so `amd_order` can assume nothing and pays for two preprocess passes on every call. Oblio owns
its matrix type and enforces validity once, at construction, which is the right place to pay. Pass
5 is the price of AMD's generality, and Oblio has already bought its way out of it.

The symmetrization is still worth having correct, and it is agnostic to storage. Checked on 200
graphs, a full pattern with a diagonal and a lower-triangle pattern with a diagonal produce the
identical adjacency, both equal to the original graph. Their reported symmetry differs, 1.0 against
0.0, and that is right rather than a defect: symmetry describes the STORED pattern, and one
triangle shares nothing with its transpose. `Amd.cpp` reports it the same way, so a 0.0 there means
half storage rather than trouble.

For genuinely unsymmetric input, ordering A + A' is a heuristic and not a derivation: the LU factors
of A sit inside the Cholesky factor of A + A' under the usual no-cancellation assumption, so it
bounds the fill without predicting it. Oblio is a symmetric solver, so that case is AMD's
generality rather than anything on Oblio's path.

**amd4_preprocess** produces R, the row form of the pattern with duplicates removed, which is the
pattern of A transposed. `Amd.cpp` calls it whenever the input may be unsorted or duplicated, since
R + R' is A + A' and R is clean. The deduplication is the usual stamp, `flag[i] == j` meaning row i
has already appeared in column j.

**amd4_aat** turns the two forms into adjacency lists and reports what the input looked like. Two
things are dropped: the diagonal, a self loop that says nothing about fill, and the distinction
between A(i,j) and A(j,i), since either one forces the same elimination. The symmetry it reports is
the vendored definition, with B the strictly triangular parts:

```
sym = nnz(B & B') / nnz(B),   or 1 when nnz(B) is zero
```

**One deliberate deviation.** `Amd.cpp` computes these counts with a two-pointer scan that walks
the upper and lower triangles together, which saves a pass but needs sorted columns. Ours asks the
row form whether the transposed entry exists, one stamp and one comparison, which is what every
other function in the file does and works on unsorted input. Checked against the definition on 400
random patterns with duplicates: same adjacency, same nz, same nzdiag, same symmetry, same
nz(A+A').

`matrix1` is the one example given as a matrix rather than a graph, six by six and deliberately
awful: unsymmetric, with a diagonal, with duplicate entries and with one column whose rows are out
of order. It reports symmetry 0.222 and is the only example that exercises this path.

### The bug pass 5 found, which had nothing to do with pass 5

Running matrices rather than graphs produced a bound of -1, and tracking it down turned up a defect
that had been latent since pass 2.

`num_eliminated` was doing two jobs: terminating the loop, and standing in for k in the first cap,
`n - k - weight(u)`. Those were the same number until hash detection arrived. A hash merge folds v
into a LIVE u, so v stops being selectable while every vertex it stands for is still live inside u.
The counter increments; the live count does not change. `n - num_eliminated` then understates what
remains, the cap comes out too tight, and **the bound can fall below the true degree**, which is
the one thing it must never do.

Two reasons it stayed hidden. The Python shows a negative bound while the C++ computes the same
quantity in `std::size_t`, where it wraps to something enormous and `min` quietly discards it, so
the twins agreed on the ordering and the trace diverged only where the number itself was printed.
And our own instrumentation counts how often the bound is LOOSE, never how often it is wrong in the
other direction, so 1748 checks said nothing about it.

The fix is to separate the two readings: `num_eliminated` keeps the selection count, and `num_live`
counts live original vertices, reduced only by an elimination and by a dense removal, never by a
hash merge. Re-checked across 210 graphs and 200 matrices: no bound below its exact degree
anywhere.

### Pass 6: amd_postorder

Every layer so far has emitted its pivots in the order it chose them. `amd_order` does not: it
builds the assembly tree and emits a postorder of it. The tree is already there and costs nothing
to record, since a clique's parent is whichever clique absorbed it, which the eliminator already
reports and aggressive absorption already knows.

**Why this is free in quality terms.** Any postorder of the assembly tree gives the same factor. A
node is numbered after all of its descendants either way, so no fill moves and nnz(L) cannot
change. Checked directly: the reported nnz(L) is accumulated during elimination and knows nothing
about the output order, and it still matches a symbolic factorization of the emitted permutation on
all 210 test graphs. The permutation itself differs from raw elimination order on 155 of them, so
the pass is doing work rather than reordering nothing.

**What it buys is locality.** Children finish before their parent starts, so a child's update is
consumed while it is still warm, and a supernode's columns come out contiguous rather than
scattered through the ordering. For a multifrontal or supernodal solver that is the difference
between a stack of pending contribution blocks and a heap of them, which is exactly Oblio's
concern in `ElmForest`.

**Two details from `Amd.cpp`, both about which child goes first.** The child lists are built by
walking the elements downward, from n - 1 to 0, so a list comes out ascending. Then the biggest
child by front size is moved to the END of its list, so the largest subtree is traversed last and
the stack of pending updates stays as small as possible for as long as possible. Front size is
`nvpiv + degme`, the pivots plus what they reach, which we already compute for the nnz accounting.

The traversal is an explicit stack rather than recursion, as the vendored `AMD_post_tree` is, and
for the same reason: the tree can be n deep on a path graph.

One caveat the vendored file states about itself and that applies here too. Mass elimination merges
vertices into a pivot that were not necessarily adjacent to each other, so the result is not an
exact elimination tree postorder. It is a postorder of the assembly tree AMD actually built, which
is what the factorization will follow.

### Pass 7: amd_valid and the Control/Info interface

The last pass, and the one that shows most clearly what `amd_order` has taken on.

**amd_valid is not a yes or no.** It checks that the column pointers start at zero and ascend and
that every row index is in range, and then returns a THIRD answer for unsorted columns or duplicate
entries: `AMD_OK_BUT_JUMBLED`, which is not an error but a request. It is what tells `amd_order` to
run `amd_preprocess` instead of using the pattern directly. So passes 5 and 7 are one mechanism
seen from two sides: the check tolerates what the conditioning pass exists to fix.

**Control is two knobs**, and that is the whole of it: `AMD_DENSE`, the alpha of pass 4, and
`AMD_AGGRESSIVE`, which switches aggressive absorption off. Both are now parameters here, which
makes pass 1 measurable for the first time. Measured, it barely registers: over the seven examples
and 200 random graphs the fill is identical with it on and off, all 207, and on grids it fires once
per run and changes nothing. Cheap, and on this test set worth nothing, which is a useful thing to
know before porting it.

**Info is where the bundling shows.** Fourteen fields in three groups:

- input statistics, `N`, `NZ`, `SYMMETRY`, `NZDIAG`, `NZ_A_PLUS_AT` and `NDENSE`, all of which
  passes 4 and 5 already print,
- workspace telemetry, `MEMORY` and `NCMPA`, the peak workspace and the count of garbage
  collections in the vendored flat pool, which have no meaning here since there is no pool to
  compact and are deliberately omitted,
- and a factorization cost PREDICTION: `LNZ`, `NDIV`, `NMULTSUBS_LDL`, `NMULTSUBS_LU` and `DMAX`.

That third group is the fifth thing `amd_order` bundles, and the largest. It predicts the fill and
the operation counts for both LDL' and LU, plus the largest frontal matrix. Two observations. Our
nnz(L) turns out to have BEEN `AMD_LNZ` all along, computed from the same expression, so part of
Info was implemented in pass 4 without being labeled. And `AMD_NDIV` comes out exactly nnz(L)
minus n on all 207 graphs, which is the internal consistency check the two counts owe each other.

It is also a prediction with two approximations already baked in, an upper bound wherever dense
rows fired and inexact wherever mass elimination merged non-adjacent vertices. A symbolic phase
computes the same quantities exactly, from the permuted matrix, which is the argument of the
section below in miniature.

### Three complexity defects in amd2, and their fixes

Checking whether amd2 matched the vendored asymptotic cost turned up three places where it did not.
All three were in pass 2, the hash detection, and none was visible in any output.

**Sorting to build the hash key.** The key was `(tuple(sorted(A[u])), tuple(sorted(I[u])))`, and
the exact test sorted again, so the detection cost O(d log d) per vertex against `Amd.cpp`'s O(d).
The fix is what the vendored code does: a SUM of the indices, because addition has no order and
neither do the sets, reduced modulo n. The exact test is now stamp-one-side-and-count-the-other,
the same membership test as everywhere else in the file, so it costs one pass and no sort.

**Purging the merged vertex, which was the serious one.** On every hash merge the old code ran a
pass over every live clique and every adjacency list to remove v, which is O(n + total structure)
per merge against `Amd.cpp`'s O(1). The vendored routine never removes a merged variable from
anything: it sets `Nv[v] = 0` and leaves it in place, weighing nothing, which is the same move the
dense prepass makes.

Nothing is lost by leaving it, and the reason is the merge condition itself. The merge required
A[u] - {v} == A[v] - {u}, so every list holding v holds u as well. v is therefore redundant
wherever it appears and never the only way to reach anything, so the walks can simply skip it. That
is one comparison per entry in `amd2_neighbors`, and no asymptotic cost.

The clique degree invariant survives unchanged, which is worth noting since pass 3 depends on it. v
keeps its slot with weight zero while u's weight grows by v's, and since I[u] == I[v] licensed the
merge, every clique holding v holds u, so no clique's weighted size moves.

**One bug found while fixing this**, of the kind that produces silence rather than a wrong answer.
Seeding the hash sum with the vertex's own index makes two distinct vertices unable to collide,
ever, so the detection quietly stopped firing: 0 merges where there had been hundreds. The fill and
the permutations all stayed valid, because a missed merge is a missed opportunity and not an error.
Only the merge counter showed it. A mechanism that can fail silently needs a counter watched, which
is the same lesson the bound-below-exact bug taught two passes earlier.

Every merge is now verified against the definition on 311 graphs, with the reachable sets computed
from plain sets so the check shares no state with the test it audits.

**A third, smaller one in the C++ only.** The hash groups were held in a `std::map` keyed by the
hash value, which costs a log per insertion and a node per group, where the Python used a dict and
paid neither. `Amd.cpp` uses `Head[hval]`, an array indexed by the hash value, which is the obvious
structure once the key has already been reduced modulo n. Both twins now use an array of buckets
allocated once and cleared only where it was used, which removes the log AND the last structural
difference between the two files.

The per-iteration cost accounting that settles whether all this adds up to parity with the vendored
routine is in the complexity section below, since that is where it belongs as a reference.

### Where amd2 lands, with all seven passes in

Against the vendored `amd_order`, on the seven examples: same nnz(L) 7 of 7, same permutation 0 of
7. On 60 random graphs: same nnz(L) on 54, ours lower on 4, ours higher on 2, and again the
permutation never matches.

The fill agreement is the meaningful number and it is high. The permutation disagreement is
expected and says little: both routines now postorder, so a single tie broken differently early on
reorders whole subtrees downstream, and the six passes that change the ordering all change it in
ways that compound. mmd2 showed the same pattern for the same reason, and that section works
through one such tie in detail.

Being lower on 4 and higher on 2 is not a claim to be better. It is the expected spread when two
implementations of an approximate heuristic break ties differently, and on a test set this small
the difference is noise rather than signal.

### What Oblio would take from amd2, and what it would leave

The passes are not equally useful to a solver that has its own symbolic phase, and it is worth
saying which is which while the reasons are fresh rather than rediscovering them at port time.

| pass | what it does | Oblio |
|---|---|---|
| 1 | aggressive absorption | probably yes |
| 2 | hash supervariable detection | probably yes |
| 3 | scan 1 by subtraction from a maintained clique degree | yes |
| 4 | dense row and column removal | maybe |
| 5 | `amd_aat`, `amd_preprocess` | almost nothing |
| 6 | `amd_postorder` | no |
| 7 | `amd_valid`, `Control`/`Info` | mostly no |

**Passes 1, 2 and 3 are ordering, and stay.** All three are internal: they change what the ordering
costs or which pivot it picks, and nothing downstream can do them instead. Pass 3 is unambiguous,
being pure cost for identical output. Passes 1 and 2 change the ordering, so they are worth having
but worth measuring first, which is what the fill figures in the earlier subsections are for.

**Pass 4 is a judgment call.** Removing dense rows is an ordering decision and belongs here, but it
is matrix-dependent and never fires on anything in this test set. Its cost is that nnz(L) becomes
an upper bound, which matters not at all for Oblio, since `SymbolicEngine` computes the exact
structure later and never consults an ordering's estimate.

**Pass 5 is almost entirely dissolved by Oblio's invariants**, as the pass 5 subsection sets out.
What survives is skipping i == j when a matrix becomes a graph.

**Pass 7 splits three ways.** The validity check is a `SparseMatrix` constructor invariant, so it
does not belong to the ordering at all. Of Control, alpha travels with pass 4 and the aggressive
switch travels with pass 1, so neither is a separate decision. Of Info, the input statistics are
free diagnostics worth keeping, the workspace telemetry does not apply, and the cost prediction is
`SymbolicEngine`'s work done early with worse information.

**Pass 6 is redundant, and this is the clearest case.** Oblio computes the elimination tree in
`SymbolicEngine` and postorders it in `ElmForest`, after supernode detection, using real front and
update sizes. Two permutations differing only by a tree postorder yield the same tree up to
relabeling, so Oblio's phase reconstructs the same forest from either input and AMD's postorder is
simply redone. Oblio's is also better informed: it runs on the supernodal tree with real sizes,
where AMD runs on the assembly tree it happened to build and admits, in its own comment, that mass
elimination makes that not an exact elimination tree postorder. Same idea, better information,
later stage.

The one thing worth keeping from pass 6 is not code. AMD moves the biggest child last so the stack
of pending updates stays small, and that is the same stack-minimization argument `ElmForest` makes,
arrived at independently. Useful as corroboration when tuning that traversal, and nothing more.

### Why the ordering does any of this, and why the answer is that it should not

A full solver has a pipeline: order, elimination forest, symbolic factorization, numeric
factorization, solve. Each stage takes the previous stage's artifact and produces its own. Several
of the amd2 passes are the second and third stages done early, inside the first, and the case
against that is not stylistic.

**The merge is lossy, which is the actual defect.** AMD builds the assembly tree internally, in
`Pe`, postorders it, returns a permutation and throws the tree away. The work is done and the
result discarded, with only a side effect kept. A stage boundary that works passes its structure
forward: the ordering hands over a permutation, the forest hands over a tree, the symbolic phase
hands over a pattern. AMD does the second stage's work and then declines to emit the second stage's
output, so the next stage has to do it again from scratch.

**What it does emit is weaker than what the next stage would compute anyway.** Its tree is the
assembly tree, which mass elimination has made approximate by its own admission, and its nnz(L) is
an upper bound rather than a count once dense rows appear. So a downstream stage gets two options
and both are bad: trust an approximation it never asked for, or recompute and waste the work. There
is no reading of the merge on which it helps.

**The confusion is not incidental either, and this experiment produced an instance of it.** The bug
in pass 5 was exactly this failure mode in miniature: `num_eliminated` carried two meanings, loop
termination and the live vertex count, and they diverged silently the moment a third mechanism
arrived. Merged stages breed shared state with more than one reading, and those readings come apart
under change. A pipeline where each stage owns its own artifacts does not have that failure mode
available to it.

**One qualification, and it is about packaging rather than structure.** The layering does exist
inside SuiteSparse. `amd_2` takes and returns the working arrays, `Pe` among them, so a caller that
wants the tree can reach it. It is `amd_order`, the convenience wrapper, that conditions the input,
orders, postorders, estimates the fill and then returns only a permutation. So the criticism lands
on the entry point that most callers use, not on the algorithm's own decomposition.

For Oblio none of the assumptions behind that wrapper hold. The matrix is valid by construction,
the elimination forest is a stage, supernode detection is a stage, and the exact fill is known
downstream. The ordering can stay an ordering: take a graph, return a permutation, and let each
later stage compute its own artifact with the better information it has. That is the reading the
catalogue above follows.

## Why md3 reorders, and what would align it

md1 and md2 are guaranteed to produce the same order on every graph, and the guarantee is
structural rather than empirical. The invariant is that md2's neighbor set equals md1's
adjacency at every iteration: the union of the explicit adjacency with the cliques containing u
is exactly the neighborhood the filled graph would have, which is George and Liu's
reachable-set theorem, and it holds by induction because each elimination adds to md1
precisely the edges md2 records implicitly as the new clique. Equal neighbor sets give equal
degrees, and both drivers select with the same rule, a scan over ascending vertex numbers
keeping the first strict minimum, so the tie-break matches too. Same pivot, same resulting
state, and the induction closes.

md3 breaks that, and not by a tie. On 400 random graphs with 4 to 12 vertices and density
0.3, md2 and md3 disagree on the order for 75 of them, about one in five, while nnz(L) is
identical on all 400. The disagreement is systematic and fill-neutral.

The mechanism is a policy difference. When u merges into pivot p, md3 eliminates u
immediately. In md2 the degree of u after p is gone is d - 1, where d was the pivot's degree,
and nothing says that is the minimum: some vertex elsewhere in the graph can be smaller, so
md2 eliminates several others first and reaches u later. On the eight-vertex example below,
md2 gives [2, 5, 6, 0, 1, 3, 4, 7] and md3 gives [2, 5, 6, 0, 4, 1, 7, 3], both with
nnz(L) = 25. Vertex 4 merges into 0 and jumps from position 6 to position 4. No tie-break can
move it forward in md2, because 4 is not tied for the minimum at that iteration; it is simply not
selected.

```
0: {1, 3, 4, 6, 7}    4: {0, 1, 3, 6, 7}
1: {0, 4, 5, 6}       5: {1, 3, 7}
2: {3}                6: {0, 1, 4}
3: {0, 2, 4, 5, 7}    7: {0, 3, 4, 5}
```

So alignment has to run the other way: give the earlier layer md3's policy rather than give
md3 a tie-break. For md2 that works exactly. The change is one loop after the prune loop,
using md3's own test, no weights and no supervariables, since the vertices are eliminated
rather than merged and the degree stays unweighted:

```python
    merged_vertices = []
    for u in neighbors:
        if not A[u] and len(I[u]) == 1 and I[u][0] == pivot:
            I[u] = []
            eliminated[u] = True
            merged_vertices.append(u)
    if merged_vertices:                 # one compaction pass, not a removal each
        tag += 1
        for u in merged_vertices:
            mark[u] = tag
        C[pivot] = [v for v in C[pivot] if mark[v] != tag]
```

with the driver emitting `order.append(pivot)` then `order += merged_vertices`, and its loop
becoming `while not all(eliminated)` since an iteration can now retire several vertices. Checked
against md3 on the same 400 graphs: zero disagreements. The cost is a constant factor. The
prune loop already touches every neighbor once, the test is O(1) per neighbor, and this is
one more pass over the same set.

For md1 the same alignment is not available cheaply, and the reason is the useful part of
this. In md1's flat graph the natural test is exact: u is indistinguishable from the pivot
when A[u] == neighbors - {u}, meaning everything u still sees lies inside the clique. On the
example above that fires for all four of 1, 3, 4 and 7 at the iteration that eliminates 0, where
md3 merges only 4.

md3's test is conservative. Requiring `not A[u] and I[u] == {pivot}` demands that the
explicit adjacency be empty and that the new clique be the only clique containing u. Vertex 1
also belongs to c5, and although c5's members happen to lie inside the new clique, the test
does not look, so 1 is not merged. The exact test would merge it.

Two things follow. First, md3's order is not "the minimum degree order with mass
elimination"; it is the order produced by one particular cheap sufficient condition, and a
different sufficient condition would give a different order at the same fill. Second,
matching md3 from md1 would mean reproducing that condition, which means reproducing I and C,
which means being md2. The conservatism is representation-dependent, not incidental.

That also names a real gap rather than a defect: the test misses genuine supervariables.
Catching them needs a second mechanism, comparing vertices against each other rather than
against the pivot, which is what the next section is about.

## Detecting supervariables against each other, in mmd2 and amd2

md3's test is positional: it asks whether a neighbor u is indistinguishable from the PIVOT,
in the iteration that just created the clique. Two vertices can be indistinguishable from each
other with neither absorbable into the pivot, and no sharpening of a pivot-relative test will
find them. Both mmd2 and amd2 carry a second mechanism for exactly that population, and the two
mechanisms are different, which is worth recording because the goal is shared and nothing
else about them is.

amd2 hashes. During the pass it is already making over the reached set, it computes a hash of
each survivor's structure, buckets by it, and runs an exact comparison only within a bucket:

```python
        survivors = [u for u in pivot_clique if not eliminated[u]]
        by_hash = {}
        for u in survivors:
            key = (tuple(sorted(A[u])), tuple(sorted(I[u])))
            by_hash.setdefault(hash(key), []).append(u)
```

then, inside a bucket, `(A[i] - {j}) == (A[j] - {i}) and C[i] == C[j]`. Note that this is the
structural test, so amd2 inherits its conservatism: a pair whose reachable sets agree while
their incidence lists differ is missed here too, which is one reason aggressive absorption
travels with the hashing, since removing a contained clique is one way such a difference
disappears. graph5 is the extreme case of that: there the absorption alone recovers the
supervariable and the hashing is not needed at all.

Three things about that shape. The hash is a filter and never the decision, so a collision
costs a comparison and nothing more; without it the pass would be quadratic in the reached
set. The test removes
each vertex from the other's adjacency before comparing, since indistinguishable vertices are
adjacent to each other and would otherwise never match. And the absorption goes the other
way from mass elimination: j is folded into the live supervariable i, which stays a candidate
carrying the combined weight, rather than being eliminated with a pivot.

mmd1 reaches the same vertices through the vendored code's q2h path. mmdelm stashes each
reached vertex's pruned adjacency count as fwd[rn] = nq+1, and mmdupd routes the nq == 1
cases into a separate list where it merges indistinguishable pairs. Same population, entirely
different route, and no hashing.

Neither is implemented in our mmd1, which is why the not-implemented list above names the q2h
path. Our amd2 does implement the hash detection, along with aggressive absorption, and the
amd2 file header calls those out as the two mechanisms beyond md5 that are not about the
degree at all.

Neither mechanism is what its layer is named for, and both are separable from it. mmd1 is
multiple elimination: eliminate a whole independent set of minimum-degree vertices before
updating any degrees, so one expensive degree update pass serves many pivots. amd1 is the
approximate degree: replace the exact size of the union with a bound computable in one pass.
Those are the ideas. Hash detection, the q2h path and aggressive absorption ride along
because both layers already sweep the reached set and the information is at hand. The
independence runs both ways: the approximate degree works with no hashing at all, which is
what amd1 is, and hash detection could be bolted onto md5
with no approximation anywhere.

The practical consequence is a coarseness ordering. md3, md4 and md5 merge only against the
pivot, so their supervariables are the finest. mmd1's are at least as coarse as ours and
sometimes coarser. amd2's are coarser again where the hash finds pairs the elimination never
brings into the same iteration. Coarser is not automatically better: a supervariable is a
commitment to eliminate its members consecutively, and while that commitment is fill-free by
construction, it removes choices the picker would otherwise have had.

## Open question: is the fill guaranteed equal across md1 to md3?

md1 and md2 produce the same order on every graph, so their fill is trivially the same. md3
reorders, as the section above records, yet every measurement so far says its nnz(L) is
identical to md2's. The equality is not something we have proved.

What the measurements cover. 1500 random graphs, 4 to 18 vertices, densities from 0.15 to
0.8: the order differed on 411 of them and nnz(L) on none. Then 62 structured graphs, all
grids from 2 by 2 through 5 by 5, thirty random trees, and chained-clique families built
specifically to force heavy merging: again no difference. With the earlier runs this is on
the order of two thousand graphs and zero disagreements.

What is actually proved is local and weaker. When u is indistinguishable from pivot p, its
neighborhood after p is eliminated is exactly the new clique, so eliminating u next creates no
fill whatever. That is the classical justification for mass elimination, and it establishes
that md3's extra iteration is free at the moment it is taken.

The gap is that this does not settle the global claim. The two orders diverge, md3 taking u
immediately where md2 takes it several iterations later, and from that point the runs face
different graphs and make different subsequent choices. A local exchange argument does not
obviously extend across a run whose pivot sequences have separated, so the equality observed
is stronger than the statement available to justify it.

Two ways to close it. Find the invariant that makes the divergent runs produce identical
filled graphs, which would explain why deferring u never leads md2 into a better or worse
configuration. Or push the search where it is most likely to break: a merge occurring early,
with the deferred vertex separated from its group by several unrelated eliminations before
md2 reaches it. Until one or the other lands, this belongs here as a question and not in the
claims above.

## What the literature proves about these algorithms, 2026-08-16

Looked up rather than recalled, because a claim that "MMD has no published complexity" was made in
this project and is false.

**Heggernes, Eisenstat, Kumfert and Pothen, "The Computational Complexity of the Minimum Degree
Algorithm", NIK 2001 and ICASE 2001-42.** Bounds for quotient-graph implementations under the
`O(n + m)` space constraint every practical code obeys. **`n` is vertices and `m` is EDGES**,
`m = |E|`, so for a `SparseMatrix` holding both triangles with the diagonal, `m = (nnz(A) - n)/2`.
On our 800 square that is 1278400 against nnz(A) of 3196800: `m` is about 2n in 2D and 3n in 3D.

| algorithm | bound | |
|---|---|---|
| MD | `O(n^2 m)` | tight |
| MMD | `O(n^2 m)` | tight |
| AMD | `O(n m)` | tight |

**Where the asymmetry comes from, and it is exactly what our differential measures.** In MD and
MMD the degree update may examine one enode's s-adjacency once for EVERY snode in the reachable
set, since "in the worst case, the same enode e can belong to the e-adjacency of every snode r in
the reachable set". That inner term is `O(nm)`, over n steps. In AMD "each edge in this local graph
is examined at most twice, once from each of its endpoints", so the whole thing is `O(nm)`.

**MMD gets no better bound than MD**, which is worth knowing before optimizing it: multiple
elimination reduces how OFTEN a degree update runs, not what one costs. That is consistent with
what this tree measures, `MMD3` beating `MMD1` by about 2.5x through doing fewer updates on
identical asymptotics.

**All three bounds are tight**, met by one graph family in the paper: a hub, 2k outer and 2k inner
vertices, edges between `x_i` and `y_j` where `|i - j| < k`. MMD needs a clique added to it.

### The three qualifications, which are the ones that matter for us

- "These bounds are for nearly dense graphs. Fortunately, these bounds are not often observed for
  problems that are solved in practice."
- For AMD on BOUNDED DEGREE graphs, Amestoy, Davis and Duff prove a tighter bound. Their own
  statement is that when the nonzeros per row of A are constant and independent of n, AMD takes
  `O(|L|)` time. **That is the bound that applies to our grids**, at 4.9 and 6.9 nonzeros per row,
  and it is the only published bound our benchmark can be read against.
- The paper's closing line is an open problem: "A further development of this work is to identify
  graph classes with provably better MD time complexities."

### The AMD paper's own bound, and its assumptions

`O( sum over k of |L_k*| . |(PAP')_k*| )`, which reduces to `O(|L|)` for bounded row counts. It
holds under two assumptions the paper states explicitly and warns about: **"no (or few)
supervariable hash collisions"** and **"a constant number of garbage collections"**, each collection
costing `Theta(|A|)`. "In practice these assumptions seem to hold, but the asymptotic time would be
higher if they did not." Both are now measured in this tree: `AMD_2` tests 0.33 pairs per pivot at
every size on both grid families, and `Info[AMD_NCMPA]` reads zero on every case run.

One discrepancy between the paper and the shipped code, since we aligned to the code: the paper's
hash is `(sum of j in A_i + sum of e in E_i) mod (n - 1)` plus 1, and `Amd.cpp` computes
`hval % n`.

### Newer theory, for completeness

Cummings, Fahrbach and Fatehpuria, SODA 2021, give an EXACT minimum degree ordering in `O(nm)` with
a matching lower bound, noting that before it no exact algorithm better than `O(n^3)` had been
proven. Fahrbach et al. compute a `(1 + eps)`-approximate greedy minimum degree ordering in
`O(m log^5(n) eps^-2)`, which is a theoretical milestone rather than a practical one.

### What this means for the measurements in benchmarks/ordering

**None of the worst-case bounds bind on grids**, so the exponents fitted there, 1.02 to 1.08 against
nnz(A), are not to be compared with `O(nm)` or `O(n^2 m)`. They are constant-factor behaviour plus
whatever log lives inside `|L|`. What the literature does supply is the STRUCTURAL reason AMD is
cheaper, every edge twice against once per reaching snode, which is the same fact the per-pass
differential shows as identical visit counts between our code and the vendored one.

## Complexity: matching the vendored cost without the vendored style

The goal is the same asymptotic cost as `Mmd.cpp` and `Amd.cpp`, reached in our own style
rather than theirs. Style and complexity are separate questions, and it is possible to write
perfectly modern code that is asymptotically worse, which is what happened here in two places.
Both were found by counting elementary operations rather than by reading.

**The real work was always right.** `md*_neighbors` makes one pass over `A[u]` and one over each
element's member list, which is what `mmdelm` and AMD's inner loop do. No set is unioned twice.
Hashing instead of a mark array is a constant factor. That is the economy symbolic factorization
gets from the elimination tree, in the form available during ordering.

**Two things were asymptotically wrong.**

The driver loop condition, `while not all(eliminated)`, is an O(n) scan per iteration, so the
condition alone is quadratic over the run. On a path of 400 vertices it cost 80800 elementary
iterations against 1596 of real neighbor work. It is now a counter, incremented by
`1 + len(merged_vertices)` per elimination. In md3 and md4 this was only a constant factor,
since their pivot search is already O(n), but in md5 and mmd1 it defeated the buckets outright.

The mass elimination block stripped a merged vertex from every clique:

```python
    for u in merged_vertices:
        for clique_members in C.values():
            clique_members.discard(u)
```

which costs O(number of cliques) per merged vertex. That wide scan was deliberate, described in
this README as defensive against a test that admits a u belonging to more than one clique. The
measurement made the defense untenable: 4247 elementary iterations on a 20 by 20 grid against 28283
of real work, growing faster than the work it accompanied. It is now

```python
    if merged_vertices:                 # one compaction pass, not a removal each
        tag += 1
        for u in merged_vertices:
            mark[u] = tag
        C[pivot] = [v for v in C[pivot] if mark[v] != tag]
```

sound because the merge test requires `I[u] == [pivot]`, so no other clique holds u. A later pass
replaced the per-vertex removal with the single compaction above, since erasing m vertices one at
a time from a clique of size d costs O(m d) where one pass costs O(d + m). Both fixes
are in md3, md4, md5 and mmd1, in both twins. Measured on the 20 by 20 grid: loop 14800 down to
34, clique scan 4247 down to 47, real work 26408. No ordering moved, checked on the seven graphs
and 200 random ones, and mmd1 at delta = -1 still reproduces md5 exactly.

**One place the Python is asymptotically worse than the C++, and it stays.** Both twins hold a
bucket as a list pushed and popped at the head, so the pop is O(1) in each. Filing and unfiling
are not: the C++ splices a doubly linked list in O(1) through `head`, `next` and `prev`, while
the Python does `insert(0, u)` and `remove(u)` on a list, which are O(bucket). Closing that would
mean mirroring the link arrays in Python and giving up the readable list, which is more machinery
than these files should carry. Documented rather than fixed, and noted in the md5 and mmd1
headers.

### Matching the vendored cost: what the containers became

With the two fixes above the prototypes performed the same NUMBER of graph operations as the
vendored routines, but not at the same cost per operation. `std::set` and `std::map` add a factor
of log, which is asymptotic rather than a constant and therefore not something tuning removes,
and a sorted vector trades that for a merge per union. Both are now gone, in both twins.

Where the costs were, and what replaced them:

| operation | before | now | vendored |
|---|---|---|---|
| build a neighbor set | `std::set`, then a merge | one pass per source, marks | the same |
| find a clique's members | `std::map` lookup | index a vector by clique id | the same |
| test v in the new clique | `count(v)`, then search | `mark[v] == tag` | `marker[v] < tag` |
| delete a pruned edge | `erase(v)` each | compaction in place | the same |
| drop absorbed cliques | `set_difference` | stamp and compact | the same |
| add the pivot to I[u] | `lower_bound` then insert | `push_back` | the same |
| file or unfile a bucket | `std::set` insert, erase | splice a list | `fwd`/`bwd`, `Next`/`Last` |

**A, I and the clique member lists are plain vectors, UNSORTED.** Nothing needs the order.
Membership comes from a mark array of size n stamped with a monotone tag, so a query is one
comparison and the array is never cleared, which is `mmdelm`'s `marker[nb] < tag` and Oblio's own
`SymFactorEngine`, where the comment reads "as an index is added to the index set of supernode kk
it is marked with kk, which makes is it already there a single comparison".

**C is indexed by clique id.** An id is a vertex number, so C is a vector over 0..n-1 with a
liveness flag, and the lookup is direct.

**The buckets are head plus next and prev arrays over n**, with a `filed` flag so unfiling is
idempotent, which matters because a vertex evicted early in an mmd1 batch can be merged away by a
later pivot in the same iteration. Filing, unfiling and popping are O(1).

Two things follow from the tie-break rather than the cost. An intrusive list has no order to take
a minimum from, so the pop takes the head, which is whatever was filed last. That is the vendored
convention, it is why md5 and mmd1 order differently from md1 through md4, and it was accepted
deliberately: best complexity wins over a nice ordering, and the ordering is different rather
than worse.

**The Python moved with it.** A, I and C are lists there too, with their own mark array and tag,
because the twins are checked by comparing traces and a trace shows the order the structure
holds. The set algebra that made md2 read as mathematics is gone from the code; that was the
price of the check, and the check is what catches drift.

**But not from the file.** The set view is what the flat containers cost, and it is worth more
than the code was, so it is kept in comments. Every place a set operation is computed by stamping
and compacting now carries the set expression it stands for, as close to the loop as it can sit:
`reach(u) = ( A[u] | C[c] for c in I[u] ) - {u}` over the neighbor walk, the six lines of set
algebra that are the eliminator over the eliminator, `filed = live - reached` over mmd's batch,
`dead = { c : C[c] <= C[p] }` over amd2's aggressive absorption. The reading stays union, difference
and containment while the code stays flat. md5's header states the convention and every layer
from md1 up follows it.

Tags are threaded through return values, since Python has no reference parameters. One sort
remains, at construction, because the input is given as sets in Python and as ascending literals
in C++; it is outside every loop, which is the same place `SymFactorEngine` puts its final sort.

What is left for the ordering code proper is the weight array, and it is a consequence rather than
a choice: once a supervariable's members are a chain rather than a list its size stops being O(1)
to read, which is where `weight` earns its place. The shared pool is not on that list, and used to
be: see the correction at the end of the next subsection, where the part of `Iw` that is algorithm
rather than Fortran is separated from the part that is not.

### The garbage collection, which is not garbage collection

`ncmpa` counts them and the vendored comments call them garbage collections, which is misleading
enough to be worth a section, since this is the piece our prototypes have not built. The conclusion
is that they should not build it either, and the reasoning is worth keeping.

**The mechanism.** Every variable-length list lives in one flat array `Iw`: each variable's
adjacency, each element's pattern, with `Pe[i]` the offset and `Len[i]` the length. A list cannot
grow in place, since the next list starts immediately after it, so a list that grows is COPIED to
the end at `pfree` and its old slot becomes a hole. When `pfree` reaches `iwlen` everything live is
slid down to close the holes. `Amd.cpp` suggests `iwlen = 1.2 * pfree + n`, so the pool starts with
twenty percent slack and a compaction is what happens when that runs out.

Nothing is reachability-analyzed and nothing is freed. It is defragmentation of a bump allocator,
and the name is the confusing part.

**Why it exists, which is not a property of the algorithm.** `amd_2` allocates nothing. The caller
supplies the workspace and the routine promises to fit inside it, which is why `iwlen` is a
parameter rather than an implementation detail. Every piece of the apparatus follows from that one
promise: `pfree`, the copy-on-growth, the compaction, `ncmpa`, the `FLIP` encoding that lets the
sliding pass find object boundaries without a second index. None of it is minimum degree. Remove
the promise and the whole subsystem evaporates.

A library cannot know its caller's memory situation, so AMD accepted the constraint and paid for
it. Oblio allocates its own memory and is not buying anything with it.

**What the pool is actually worth, weighed honestly.** Three benefits get claimed for it and they
do not survive equally.

Contiguity is the weakest, and stating it carelessly overstates it. A vector of vectors is
contiguous WITHIN each list, and that is what the inner loops walk. The indirection is one pointer
per list, not one per element, and at a few hundred lists per iteration it is noise. The pool is not
rescuing element-level locality, because nothing lost it.

Footprint is real and does not matter here. Three pointers per vector is 24 bytes per list before
a single edge is stored, so 24n of pure overhead. At n in the millions that is tens of megabytes.
It is also exactly the kind of cost Oblio spends on purpose elsewhere: the full matrix is kept so
that the elimination forest can walk a row by reading above the diagonal, and fronts and update
blocks are rectangles and squares rather than trapezoids and triangles. Spending memory to keep
loops simple and to hand BLAS the shapes it likes is a deliberate trade, and a solver that makes it
everywhere else has no reason to refuse it in the ordering.

Allocation churn is the one with a genuine TIME cost. Lists grow constantly here, and a growth that
goes through the general allocator is a real cost repeated often. That is the argument worth
taking seriously, and it is not an argument for a hand-rolled pool.

**Because the modern answer is not to write a pool, it is to keep the containers and change the
allocator.** `std::pmr::vector` over a `monotonic_buffer_resource` gives arena allocation, almost
no churn and one contiguous region, while the code still reads `A[u].push_back(v)`. The algorithm
stays visible and the memory strategy becomes a line at construction. That option did not exist when
`Amd.cpp` was written, and the pool is what a library had to do without it.

**The general principle, which this is a case of.** Explicit memory management moves attention from
the algorithm to bookkeeping, and the bookkeeping then needs its own encodings, its own invariants
and its own bugs. Containers exist to stop that happening. Compilers, standards and allocators keep
closing constant-factor gaps, and `pmr` is one of those gaps already closed, so a hand-rolled pool
buys less every year while costing the same to read forever. Human time is the scarce resource, and
a design that spends it on `FLIP` tricks to save memory nobody is short of has the trade backwards.

Explicit management earns its place where the interface demands a fixed budget, which is AMD's
case, or where the footprint is the binding constraint, which is nobody's case here.

**One property worth appreciating anyway.** The pool can only shrink in aggregate over a run,
because absorption destroys more than elimination creates: an eliminated clique's pattern replaces
all the patterns it swallowed. So the compaction reclaims space the algorithm's own progress keeps
making available, and `ncmpa` stays small for exactly that reason. It is a collector for a heap
that mostly empties itself, which is a pleasing thing to have built and a poor reason to build one.

**What this means for the ordering engine.** Do not port `Iw`, `Pe`, `pfree` or `ncmpa`. The one
consequence to keep in view is that `qsize` becomes necessary only if members ever become a chain,
which they since have, so production carries a weight array where these prototypes do not.

**Superseded in part, 2026-08-01.** This section used to end by recommending per-list containers on
an arena, `std::pmr` over a monotonic buffer, if measurement showed allocation mattering. It did
matter, and the answer turned out not to be an allocator. `A[u]` and `I[u]` are the two kinds of
source of `reach(u)`, an elimination destroys one for each it creates, and so their sum never
exceeds `u`'s original degree: the pair fits in one block sized once from the pattern, and there is
nothing left to allocate per list or to put on an arena. Section 5.3 of
`archive/sparse_factorization.md` carries the argument and 5.15 records that both vendored codes
rely on it. The rest of this section stands: the compaction really is an artifact of a fixed
caller-supplied workspace, and it really should not be ported. What was wrong was reading the whole
of `Iw` that way, when only the pool around a vertex's block is archaeology.

The prototypes are unchanged and keep their per-vertex containers, which is right for files whose
job is to read as the algorithm. The engine took the block.

### amd2 against Amd.cpp, iteration by iteration

The table is the reference version of the argument above, for the completed amd2. It counts only
what a production version would run, so `amd2_exact_degree` is left out: it computes the union the
bound exists to avoid, once per updated vertex, and exists only so the trace can print the exact
degree beside the bound.

| iteration | cost | Amd.cpp |
|---|---|---|
| pick from the buckets | O(1) plus the walk, amortized | the same |
| build the new clique | \|A[p]\| + sum \|C[c]\| over c in I[p] | the same |
| prune and absorb | sum over C[p] of \|A[u]\| + \|I[u]\| | the same |
| scan 1 | sum over C[p] of \|I[u]\| | the same |
| aggressive absorption | sum over C[p] of \|I[u]\| | the same |
| the bound | sum over C[p] of \|A[u]\| + \|I[u]\| | the same |
| hash and compare | sum over C[p] of \|A[u]\| + \|I[u]\|, plus the within-bucket pairs | the same |
| refile, clique degree | O(1) | the same |

No sorts remain anywhere except one at construction in the matrix path, outside every loop, which
is where the graph path already sorts its input. No ordered containers remain. Every membership
test is a stamp and a comparison, every list edit is a compaction in place, and every per-clique
quantity is read rather than recomputed.

So the core is at parity as far as this reading goes. What is NOT at parity is everything the
prototypes deliberately do not have: the flat pool and its compaction, which the section on garbage
collection argues against porting, and the per-list allocation our containers pay instead.

Three defects had to be fixed to get here, all in pass 2 and all invisible in the output: sorting
to build the hash key, purging a merged vertex from the whole structure, and holding the hash
groups in an ordered map. The amd2 section works through each.

### Two things the seven examples were hiding

The conversion turned up a latent divergence that had nothing to do with the containers. The
Python's `for u in neighbors` iterated a set, so its order was a hash artifact depending on the
insertion and deletion history, while the C++ iterated ascending. On the seven examples the two
happened to coincide, which is why `make test` never caught it; on random graphs of a dozen
vertices they did not, and the fill and pruned edge lists printed in different orders. It is
moot now that both sides hold lists, but it was fixed twice on the way here, first with `sorted()`
at the loops that feed a printed list and then by the conversion itself.

The lesson is about the test rather than the code: seven small graphs are not enough to keep the
twins honest. The stress harness that found it builds a variant of each `.cpp` whose main runs
forty random graphs, and diffs it against the same graphs through the Python. It caught two more
divergences during the conversion, one per layer, and it is worth running after every change to
either twin.

## Translation choices

The Python is where the thinking happens; the C++ twin is written alongside it and exists to keep
it honest, since two implementations that print the same trace are unlikely to be wrong in the
same way. That was once a near-transliteration, one identifier per identifier, with the C++ using
ordered containers so the two would agree without effort. It is not any more: the C++ is written
for the vendored asymptotics, and the Python mirrors whatever structure decides an answer. The
mapping below is what survives of the correspondence.

The mapping is mechanical wherever it can be. Names translate from `live_vertices` to
`liveVertices`, one identifier to one identifier, and docstrings become the comment block above
the same function. The containers now correspond directly: a list of lists is
`std::vector<std::vector<std::int32_t>>`, the clique store is a vector indexed by clique id with
a liveness flag, a bucket is `head`/`next`/`prev` over n in C++ and a list whose position 0 is
the head in Python, and the mark array is the same array in both.

Six choices are worth recording, since none of them is forced.

**No sets anywhere.** The C++ holds `A`, `I` and the clique members as plain unsorted vectors, and
the Python holds them as plain lists. Membership is a mark array with a monotone tag in both.
That is what the vendored codes do and what `SymFactorEngine` does, and it is the reason the two
twins agree on the order a structure holds: neither is imposing one.

**Cliques are a vector, not a map.** A clique id is the pivot that created it, so the id space is
the vertex space and the lookup is direct. The display walks ids ascending, which is the one
place an order is imposed, and it is imposed identically in both twins.

**Set algebra is gone from both.** `A[u] & neighbors` became a stamp and a compaction pass, in the
Python as much as in the C++, so the two now read alike again. What the Python loses in
expressiveness it gains in saying exactly what the engine does.

**Optional arguments become defaults and a pointer.** `title=None` is `const std::string&
title = ""`, and `eliminated=None` is `const std::vector<bool>* eliminated = nullptr`, so
both defaults behave as Python's do. `std::vector<bool>` is the bit-packed specialization,
which is not obviously right, but the flags are only read one at a time here.

**Multiple return values become tuples.** `neighbors, fill_edges = ...` is a `std::pair` or
`std::tuple` unpacked with a structured binding, which keeps the call sites the same shape.
The returned sets are copies in both languages.

**The driver copies the input.** Both twins take the graph by const reference or value and
build their own working copy, so a graph can be run through several layers in one session.
The Python does `A = [set(adjacency) for adjacency in G]` and the C++ does `Graph A = G`.

Two smaller things that exist only to keep the traces identical. The order prints in Python
list syntax, `order: [0, 1, 2, 3]`, which the C++ formats by hand. And the field widths come
from `std::setw` matching Python's `f"{u:>{width}}"`, so the columns line up the same way in
both.

## Containers

These prototypes exist in two languages already and are likely to grow more, so it is worth
recording what each language actually offers, and where the vocabularies fail to line up. Three
different disciplines all get called "ordered", which is the whole source of the confusion:
sorted by key, insertion order, and no order at all. Lined up that way the correspondence is
clean.

| discipline | Python | C++ | Rust |
|---|---|---|---|
| sorted by key | none | `std::set/map` (red-black tree) | `BTreeSet/BTreeMap` (B-tree) |
| sorted, flat | none | `std::flat_set/flat_map` (C++23) | `Vec` plus `binary_search` |
| sorted, duplicates | none | `std::multiset/multimap` | none |
| hashed, unordered | `set`, `frozenset` | `std::unordered_set/map` | `HashSet`, `HashMap` |
| hashed, insertion order | `dict`, `OrderedDict` | none in stdlib | none in stdlib |

Three cells are where the mismatch bites. Python's `dict` occupies a row the other two leave
empty: a hash map that iterates in insertion order, guaranteed by the language since 3.7 and true
in CPython since 3.6. Neither C++ nor Rust has a standard container there, so the substitutes are
`boost::multi_index`, the `indexmap` crate, or a `Vec` of keys kept alongside the map. Conversely,
C++ and Rust have a sorted row that Python leaves empty; `std::map` and `BTreeMap` iterate in key
order and pay O(log n) lookup for it, where Python's answer is `sorted()` at the point of use,
which is exactly what these prototypes do. And `OrderedDict` has no counterpart anywhere, because
what distinguishes it is not its iteration order but `move_to_end` and order-sensitive equality.

The hashed row hides an asymmetry in how seriously each language means "unordered". C++ leaves the
order unspecified, Rust randomizes it per process by seeding `RandomState`, and Python leaves it
unspecified but in practice stable within a run for a fixed insertion history. Only Rust makes the
non-guarantee impossible to depend on by accident. The Python trap is sharper than it looks: small
non-negative integers hash to themselves, so a set of vertex numbers usually looks sorted, right up
until a resize or a collision modulo the table size. That is an artifact and not a promise, which
is why every display in these files calls `sorted()`. The discipline is load-bearing, since the
traces must match their C++ twins line for line, so anywhere iteration order could leak into a
decision, tie-breaks above all, it has to be pinned rather than left to the container.

The ordering guarantee in `dict` is a side effect of its layout rather than bookkeeping. CPython's
compact dict is two arrays: a sparse array of indices, which is the hash table proper, and a dense
array of hash/key/value entries appended in insertion order. Lookup hashes the key, probes the
sparse array, and follows the index into the dense one; iteration simply walks the dense array
front to back. It also saves memory, because only the sparse array has to be oversized to keep the
load factor low, and its slots are 1, 2, 4 or 8 bytes rather than full entries. Deletion marks the
index slot as a dummy and nulls the dense slot, leaving a hole, and the holes are squeezed out on
resize. One consequence worth knowing: deleting a key and reinserting it moves it to the end, so
the order is insertion order as of the last insertion of that key. Updating the value of a live key
does not move it. Sets never got that treatment; a Python `set` is still one open-addressed table
of hash/key pairs, entries living wherever they probe to, with no dense array to walk and so no
order to promise.

What remains between `dict` and `OrderedDict`, none of it about iteration order, is four things.
Equality: `OrderedDict` compares order-sensitively against another `OrderedDict`, while `dict`
ignores order, and a mixed comparison falls back to order-insensitive. Reordering: `move_to_end`
exists only on `OrderedDict`, and its `popitem` takes an end, so it works as either a stack or a
queue, where `dict.popitem` always takes the last. Reverse: `OrderedDict` has supported
`reversed()` since 3.5 and `dict` since 3.8, which is no longer a difference in practice. Cost:
`OrderedDict` carries a doubly linked list to make those reorderings O(1), so it is larger per
entry and slower to build. Plain `dict` covers most uses now, and `OrderedDict` earns its place
only for `move_to_end`, for eviction from either end, or when order-sensitive equality is the
property being tested.

For these prototypes the choice was settled twice, in opposite directions. The Python layers began
as lists of sets, because the algorithm reads as set algebra that way: `md1` tested `w not in
graph[u]`, `md2` said `A[u] & neighbors` and `I[u] -= absorbed` directly. That lasted until the
goal was stated plainly as a performant ordering engine, at which point the C++ moved to flat
unsorted vectors with a mark array, and the Python followed so the traces would still match. Both
sides now hold lists and test membership with a stamp, which is what `mmd1` and `amd1` need for
their speed and what `SymFactorEngine` already does. What was lost is the notation; what was
gained is that the Python says what the engine does.

## What is still unchecked

Collected in one place, because several of these were noticed at different points and would
otherwise be scattered through the sections that produced them. None is a known defect. Each is a
claim resting on a reading rather than on a measurement, or a mechanism exercised too narrowly to
be trusted at size.

**Scaling on grids, for both branches.** The complexity claims, that mmd2 matches genmmd and that
amd2 matches Amd.cpp, rest on counting operations in both codes and comparing the expressions. The
table in the complexity section is that accounting, not a measurement. The check is a grid family
of growing size with the counters plotted against n, and it is cheap now that both are complete.

**Everything measured is small.** At most 484 vertices, and mostly far less. The q2h share, the
pair merge count, the outmatched count and the looseness of the bound all moved with size on grids,
so the behavior at real problem sizes is extrapolation from a short line.

**delta beyond the batch limit.** The vendored mmdupd also uses `mdeg + delta` in its tag window.
That belongs to machinery we excluded, so the omission is consistent, but it is the one place our
reading of delta is narrower than theirs and it has not been shown to be harmless.

**ncsub against the vendored value.** Not compared, and the reason is that mmd_order takes it as
a local and drops it rather than that it is out of reach: one temporary line in the wrapper
prints it, as `tools/hook_amd.py` does on the AMD side. Ours comes from the same expression, so
today it is a reading rather than a check, and the check is a few minutes' work whenever it is
wanted.

**Aggressive absorption measures at zero benefit.** Identical fill with it on and off across all
207 small graphs, and on grids it fires once per run and changes nothing. That is a real result on
this test set and no result at all about larger ones, where cliques nest more often. It is the
first thing to measure rather than the first thing to port.

**The dense path is only exercised with alpha forced.** At the default the threshold is never
reached by anything here, so stars and dense random graphs with alpha at 1 or -1 are the only
evidence that the pass works. Its interaction with the other passes at a realistic threshold is
untested.

**nnz(L) stops being checkable once dense rows fire.** It becomes an upper bound, deliberately,
following Amd.cpp. On graphs where no dense row appears it is still verified against a symbolic
factorization, which is every graph in the test set at the default alpha, so the oracle covers the
default path and nothing else.

**Timing either amd file measures the wrong thing.** `amd*_exact_degree` computes the union the
bound exists to avoid, once per updated vertex, so it dominates. It is instrumentation and would
not ship, but any timing run needs it removed first.

**The seven examples plus random graphs are not a structured test set.** No mesh from a real
problem, no banded matrix, no matrix with a natural supernodal structure. Two of this session's
findings came from graphs added specifically because the existing set hid something, which is
evidence the set is still thin rather than evidence it is now sufficient.

## Zooming in on md2: four ways to pick a pivot

The ladder from md1 to md5 changes one thing per rung, and by md5 four ideas have accumulated. Two
of them are independent of each other and get entangled by the order they arrive in: **whether the
degree is exact or bounded**, which the ladder introduces at amd1, and **whether it is recomputed
or maintained**, which arrives at md4 alongside mass elimination and supervariables.

They are separable, and separating them at the earliest layer where both are expressible gives four
files rather than two rungs:

```
                        recomputed            maintained
exact degree            md2                   mdm2
bounded                 mda2                  mdam2
```

Each neighbor in that square differs from its partner in exactly one thing, so any difference in
order or in fill is attributable to that thing alone. md2 is the earliest layer where the square
exists at all: md1 has no cliques, so there is nothing to bound and nothing to maintain.

**Where each file spends its work, per iteration.** `n` is the vertex count, `live` the number not
yet eliminated, and `|C[pivot]|` the size of the clique this iteration forms. Every picker walks all
`n`
slots and skips the eliminated, so the SCAN is `O(n)` in all four; what changes is how many of those
slots cost real work and what that work is.

md2, exact and recomputed:

```
picker      md2_neighbors(u) for every live u        O(n) scan, O(live) unions
eliminate   md2_neighbors(pivot)                     1 union, becomes C[pivot]
```

mdm2, exact and maintained:

```
picker      read the cached degrees array            O(n) scan, no set work
eliminate   mdm2_neighbors(pivot)                    1 union, becomes C[pivot]
deg update  mdm2_neighbors(u) for u in C[pivot]      O(|C[pivot]|) unions, after
```

mda2, bounded and recomputed:

```
picker      mda2_bound(u) for every live u           O(n) scan, O(live) additions
                                                     PIVOT-FREE bound
eliminate   mda2_neighbors(pivot)                    1 union, becomes C[pivot]
```

mdam2, bounded and maintained:

```
picker      read the cached bounds array             O(n) scan, no set work
eliminate   mdam2_neighbors(pivot)                   1 union, becomes C[pivot]
deg update  mdam2_refresh_bounds(C[pivot])           O(|C[pivot]|) additions, after
                                                     bound AGAINST C[pivot]
```

Read down the right column and the two axes are the two words that change: `live` becomes
`|C[pivot]|` going right, and `unions` becomes `additions` going down. The `O(n)` scan is in every
box and is what the third axis removes.

**What each box actually costs, per iteration.** The right column above says the shape; this says
the work, with `live` the number of vertices not yet eliminated:

```
md2      sum over live  of ( |A[u]| + sum over c in I[u] of |C[c]| )   a union per vertex
mdm2     sum over C[p]  of ( |A[u]| + sum over c in I[u] of |C[c]| )   the same union, fewer of them
mda2     sum over live  of ( 1 + |I[u]| )                             lengths, one pass
mdam2    3 * sum over C[p] of |I[u]|                                  lengths, three passes
```

Four things to read off that, and the third and fourth are easy to get backwards.

**The exact pair differ only in how many.** md2 and mdm2 compute the identical quantity by the
identical means; maintenance changes the SET of vertices it is computed for and nothing about the
computation. That is the whole of what the first axis does on the exact side.

**mda2 already captures the whole of the second axis.** Per vertex:

```
exact reach(u)   |A[u]| + sum over c in I[u] of |C[c]|      walks every member of every clique
mda2 bound(u)    1      + |I[u]|                            reads lengths
```

The bound's cost is proportional to the NUMBER of cliques in I[u], where the union is proportional
to their MEMBERS. Note the leading `1` rather than `|A[u]|`: the bound needs `len(A[u])`, a size, so
it does not walk the adjacency either. So the bound is nowhere near the price of the union, and it
does not need maintenance to be cheap.

**And mdam2's bound costs MORE per vertex than mda2's, not less.** Per vertex:

```
mda2    len(A[u]) + one length read per clique      ~ |I[u]|,   one pass
mdam2   three passes over I[u], plus outside[]      ~ 3|I[u]|
```

Both sum one number per clique in I[u]; the difference is where that number comes from. mda2 reads
`len(C[c])` straight off the clique, available immediately. mdam2 needs `|C[c] - C[p]|`, which
exists nowhere and has to be MANUFACTURED: one pass to set `outside[c] = len(C[c])`, a second to
subtract the overlap with C[p], and only then can a third read it. So each clique is touched three
times instead of once, and the extra work is not the summing but the building.

The three passes are three phases with real barriers between them: every clique must be initialized
before any member decrements it, and every decrement must be done before any sum is read, since a
clique can be named by a member appearing later in the group.

**And the reason to pay it** is that the manufactured number is smaller, so the bound is tighter,
and the manufacturing is SHARED. `outside[c]` is built once for the whole group even though it is
read once per member: a clique named by ten members of C[p] costs ten decrements in total, not ten
per member.

**So mdam2 wins on the count of vertices and loses on the cost of each.** It pays a constant factor
for a better number, and it pays it on a small group. That is the interaction the square exists to
expose: on the exact side maintenance is purely a reduction, and on the bounded side it is a trade,
slightly worse per vertex in exchange for a tighter bound and far fewer vertices.

**The two bounded boxes do not add the same things.** Same count of additions, different numbers
added. The two bounds differ in their third term and nowhere else:

```
PIVOT-FREE            sum over c in I[u]              of ( |C[c]| - 1 )
AGAINST C[pivot]      sum over c in I[u], c != pivot, of |C[c] - C[pivot]|
```

**The second is strictly the better bound**, and it is one line of set algebra. `u` is one member of
`C[pivot]`, so `|C[c]| - 1` is `|C[c] - {u}|`, and `{u}` is a subset of `C[pivot]`, so removing
`C[pivot]` removes at least as much:

```
|C[c] - C[pivot]|   <=   |C[c] - {u}|   =   |C[c]| - 1
```

**And the gap is exactly the double counting the new clique creates.** Every vertex in
`C[c] & C[pivot]` other than `u` is already counted by the middle term `|C[pivot] - {u}|`; the
pivot-free form counts it a second time and the other form does not. On a graph doing any real
filling that is most of the overcount, since a vertex reached through two cliques is usually reached
through the new one as well.

Neither removes all of it. Both still double count a vertex lying in two of `I[u]`'s cliques
OUTSIDE `C[pivot]`, which is the residual the approximation is named for and what graph4 exhibits.

**Term by term, and this is the other angle on the same pair.** Written as the code computes them:

```
PIVOT-FREE          bound(u) = |A[u]|                              -> len(A[u])
                             + sum over c in I[u] of ( |C[c]| - 1 ) -> sum of len(C[c]) - 1

AGAINST C[pivot]    bound(u) = |A[u] - C[pivot]|                    -> len(A[u])
                             + |C[pivot] - {u}|                     -> len(C[pivot]) - 1
                             + sum over c in I[u], c != pivot,
                                          of |C[c] - C[pivot]|      -> sum of outside[c]
```

**The pivot-free form has two terms, not three**, because there is no designated clique to separate
out: where `u` belongs to the new clique, it is simply one more `c` in the sum. And the first term
is `len(A[u])` in both, for the same reason, that the eliminator has already pruned `C[pivot]` out
of `A[u]`.

**The difference is one array.** `|C[c]| - 1` is a length read straight off `C[c]`, so the
pivot-free form precomputes nothing and stores nothing. `|C[c] - C[pivot]|` is a subtraction, held
in an `outside[]` array indexed by clique id, computed once per iteration and read by every member
of the group whose incidence list names that clique.

So **the tighter bound is cheap but not free**, which is worth stating plainly: it costs one pass
over the group's incidence lists per iteration, `sum over u in C[pivot] of |I[u]|`, to fill
`outside[]`. The pivot-free bound costs nothing at all beyond the additions it shares with the
other. That is a
real price for a better number, and it is small because the number is shared: without the sharing it
would be an intersection per vertex, which is the union cost the whole approximation exists to
avoid.

The pivot-free form exists in this ladder for one reason only, that mda2 recomputes and so has no
group to state the other against.

**The eliminator's line is identical in all four**, and that is not an incidental symmetry. The new
clique IS the reachable set, so forming it needs the members and not a count, and no approximation
can help. One union per pivot is the floor, and every idea in the square is about the other calls.

**The two axes cut different things, and they multiply.** One reduces HOW MANY vertices are worked
on, the other reduces WHAT EACH ONE COSTS:

```
                        recomputed              maintained
exact            n unions per iteration      |C[pivot]| unions per iteration
bounded          n additions            |C[pivot]| additions
```

Neither is a refinement of the other, and each is worth having alone. That is what mdm2 and mda2
show separately and mdam2 shows together.

**The second axis is not a constant factor**, and this is the part most easily lost. The exact side
costs, per vertex,

```
|A[u]| + sum over c in I[u] of |C[c]|          proportional to the MEMBERS of every clique
```

and the bounded side costs

```
|A[u]| + |I[u]|                                proportional to the COUNT of cliques
```

So the ratio between them is the average clique size, which grows with fill. The approximation is
cheapest to skip on a sparse graph with tiny cliques and pays most on exactly the matrices where the
ordering is expensive. Section 5.13 of `archive/sparse_factorization.md` makes the same point from
the reuse side: `|C[c] - C[pivot]|` depends on the clique and not on the vertex, so one number
serves every vertex naming it, where a union cannot be decomposed that way at all.

**A third axis exists and this square does not touch it: the SCAN.** All four files find the minimum
by walking an array of length n and skipping the eliminated, so the scan is `n` in every box, even
in mdm2 and mdam2 where the *work* has fallen to `|C[pivot]|`. Maintained degrees remove the work,
not
the walk. Removing the walk is degree buckets, which is md5, and it is a third independent idea:
file each live vertex under its degree and read the minimum off the front. Extended that way the
square would be a cube, and **amd1 sits at its bounded, maintained and bucketed corner.**

That is a statement about these three axes and nothing more. amd1 is md5 plus the bound, so it also
carries what the md ladder accumulated below md5, supervariables and mass elimination from md3; and
it carries none of the extras, aggressive absorption and hash supervariable detection, which arrive
only at amd2. What the cube does say is that the bound is not the whole of amd1's advantage over
md2: two of its three axes were already in place at md5, and the bound is the third.

**The two axes are not independent in one direction**, which is the finding this square was built to
expose. The bound is

```
reach(u)  = ( A[u] | C[c] for every c in I[u] ) - {u}
degree(u) = |reach(u)|

bound(u)  = |A[u] - C[pivot]|
          + |C[pivot] - {u}|
          + sum over c in I[u], c != pivot, of |C[c] - C[pivot]|
```

and it is stated against `C[pivot]`, which is what ties it to the other axis. A recomputing picker
must produce a number for every live vertex at every iteration, including vertices this elimination
never touched, and for such a vertex `C[pivot]` is not among `I[u]`'s cliques and there is no group
to
state the bound against. So md2's column gets a weaker form,

```
bound(u) = |A[u]| + sum over c in I[u] of ( |C[c]| - 1 )
```

which needs no designated clique. Its terms are exact for two reasons that are md2's doing rather
than the bound's: `u` is in every `C[c]` with `c` in `I[u]`, so `|C[c] - {u}|` is `|C[c]| - 1` with
no test; and `A[u]` is disjoint from all of them, because joining `c` pruned `C[c]` out of `A[u]`
and `A[u]` only ever shrinks after. So its only overcount is one vertex lying in two cliques of
`I[u]`.
It is looser than the bound above, which removes the whole of `C[pivot]` from each other clique
rather than just `{u}`.

**Maintenance buys two things, not one**, and they arrive together because they follow from the same
fact: the degree update set narrows from every live vertex to the members of `C[pivot]`.

1. **Fewer vertices worked on**, `|C[pivot]|` instead of all the live ones. This is the obvious one
   and it is what md4 is introduced for.
2. **The tighter bound becomes expressible.** Every vertex being worked on is now a member of
   `C[pivot]`, so a bound stated against `C[pivot]` applies to all of them.

The second also makes that bound CHEAP rather than merely tighter, which is a third thing hiding
inside it. `|C[c] - C[pivot]|` is one number per clique shared by the whole degree update group,
obtained by walking `C[pivot]`'s members and decrementing a counter per clique in each one's
incidence list.
No clique is opened, and the cost is `sum over u in C[pivot] of |I[u]|` for the entire group. A
recomputing picker has no group, so even if it could state the bound it would have to intersect
cliques per vertex, which is the union cost the bound exists to avoid.

The second goody is easy to miss because md4 and md5 do not use it. They compute exact degrees, so
the tighter bound is available to them and simply never taken; it is taken at amd1, which is md5
with the picker minimizing `bound(u)` instead of `degree(u)` and nothing else changed.

**Which fixes a rule for the whole approximate branch.** A layer in the recomputing column gets the
pivot-free bound because it can have no other. A layer in the maintained column gets the bound
against `C[pivot]`. So `mda2` is the only file in the ladder that uses the pivot-free form, and
`mdam2`, and any later `mda4` or `mda5`, use the tight one. The two bounds belong to the two columns
rather than to two stages of a progression.

**And there is a timing consequence.** In md2 and mda2 the degree work happens in the picker, before
the elimination. In mdm2 and mdam2 it happens in a degree update pass after it, because `C[pivot]`
does not exist until the eliminator has formed it. So the maintained pair have three phases where
the
recomputing pair have two, and a printed number for an untouched vertex is a cached value from
whenever it was last updated rather than something computed for the display.

**Results.** All four files run over the seven examples, each Python twin agreeing with its C++ twin
byte for byte:

```
          order against md2      loose picks, graph3    loose picks, graph4
mdm2      same on all seven      exact                  exact
mda2      differs on 3 and 4     4 of 12                2 of 8
mdam2     differs on 4           0 of 12                1 of 8
```

**`nnz(L)` is identical in all four boxes on all seven graphs.**

**mdm2 reproducing md2 exactly is the check that it is correct**, not a finding: maintaining a
number must not change it. If those two ever diverged, the degree update set would be wrong.

**The two bounded boxes are the finding, and they separate cleanly.** mda2's pivot-free bound goes
loose four times on graph3 and reorders; mdam2's bound against `C[pivot]` is never loose there at
all and reproduces md2's order exactly. **The tighter bound recovers what the looser one loses**,
which is the comparison of the two forms showing up as a result rather than as an argument. graph4
stays loose once even for the tighter bound, and that is the residual the approximation is named
for, two cliques of `I[u]` overlapping outside `C[pivot]`, which is what that graph was added to
exhibit.

**And the fill column says the reordering costs nothing here.** That is the same shape as the
supervariable result in 5.5, permutation differing and fill not, and a stronger version of it, since
mda2 uses the weakest bound in the square with nothing to fall back on.

Seven graphs is not evidence. What the square does establish is that the two axes are separable and
that the bounds differ in a way that is visible on graphs this small.

## Supervariable detection in production: how AMD does it, and how MMD does it

The section "Detecting supervariables against each other, in mmd2 and amd2" above describes the
mechanism as the PROTOTYPES write it, with Python tuples and `sorted`. This one describes what
`src/Amd2.cpp` and the vendored `genmmd` actually do, which is a different implementation of the
same idea and worth having written down beside the code rather than only beside the teaching
version.

Both find the same population: pairs of live vertices that have become indistinguishable from EACH
OTHER without either being indistinguishable from the pivot. Mass elimination, which every layer
from md3 up carries, finds only the pivot-relative case. What follows is two entirely different
routes to the rest.

### AMD: hash, then compare exactly

Three phases, all over the members of the new clique, in `src/Amd2.cpp`.

**What the mechanism means, before how it is done.** The specification is one loop:

```
for every pair (u, v) in C[p]:
    if A[u] - {v} == A[v] - {u} and I[u] == I[v]:
        merge v into u
```

That is `|C[p]|^2 / 2` pairs, and it is what the hash is an implementation of rather than a change
to.

**And it is the same test mass elimination makes, on a different pair.** Inside the eliminator the
question is whether `u` is indistinguishable from the PIVOT; here it is whether `u` is
indistinguishable from another member `v`. Both pairs are drawn from `C[p]`.

What makes the two look so unlike in code is that **the pivot version does not have to be run as a
set comparison.** Its test is

```
A[u] == {} and I[u] == {p}
```

two length checks and one entry. The reason is that the prune has already done the subtraction:
immediately before the test the eliminator has computed

```
A[u] = A[u] - C[p] - {p}
I[u] = ( I[u] - I[p] ) | {p}
```

so "everything u can still reach lies inside C[p]" is literally "A[u] is empty and I[u] names only
the new clique". The symmetric form `A[u] - {p} == A[p] - {u} and I[u] == I[p]` is what it MEANS,
and the prune has already reduced it to a two-field check. No such reduction exists for a pair,
since nothing has subtracted `A[v]` from `A[u]`, so those sets have to be compared directly.

```
              pair       test                                 cost
mass          (p, u)     A[u] == {} and I[u] == {p}           O(1) after the prune
hash          (u, v)     A[u] - {v} == A[v] - {u},            |A| + |I| per pair
                         I[u] == I[v]
```

**And the cheap version is also the conservative one.** graph5 is the case: vertex 4 has nothing
explicit left but belongs to c1 as well as to the new clique, so `I[u] == {p}` fails even though
everything 4 reaches lies inside `C[p]`, and the exact test would have merged it. That is one
reason aggressive absorption travels with the hashing, since removing a contained clique is one way
such a difference disappears.

The hash partitions `C[p]` so that any pair which WOULD pass the test lands in the same bucket,
which makes the pairs in different buckets skippable without changing the outcome. **It is a pruning
of the pair enumeration, not a different question.** If the key spread badly and put everything in
one bucket, the result would be identical and the cost would be the specification's; if it spread
perfectly, the result would still be identical and the cost would be linear.

As implemented:

```
INPUT   C[p], the clique this elimination formed

# pass 1: one key per vertex, |C[p]| keys, no pairs
bucket = empty map from hash value to list
for u in C[p], u live:
    key = sum over v in A[u], v live, of (v + 1)
        + (n + 1) * sum over c in I[u] of (c + 1)
    bucket[key mod (n + 1)].append(u)

# pass 2: pairs, but only inside a bucket
for each nonempty bucket B:
    for each pair (u, v) in B, u before v:
        if u dead or v dead: continue
        if A[u] - {v} == A[v] - {u} and I[u] == I[v]:
            merge v into u
```

Pass 1 costs the sum of `|A[u]| + |I[u]|` over `C[p]`. Pass 2 costs the sum of squared bucket sizes,
where the specification costs `|C[p]|^2`; that ratio is the whole saving. The two liveness guards in
the pair loop are not decoration: `u` may have been merged away by an earlier pair in the same
bucket, and so may `v`.

**The key, at the level of what it must satisfy.**

```
key(u)  = sum over v in A[u] of (v + 1)  +  sum over c in I[u] of (c + 1)
hash(u) = key(u) mod (n + 1)
```

- **A SUM, because the sets are unordered.** The map has to be invariant under permutation of each
  list. A sum is; a concatenation is not. Any symmetric function would do and the sum is the
  cheapest.
- **`+ 1`, because vertex 0 must contribute.** Otherwise `{0, 5}` and `{5}` are indistinguishable.
- **NO STRIDE, and a vertex and a clique of the same index are simply allowed to collide.** That
  costs one exact comparison and nothing else, the hash being a filter and never the decision, and
  it is what `AMD_2` does: `hval += e` and `hval += j` into one running value, then `hval % n`.

**This paragraph used to argue the opposite, and the correction is ledger entry 8.** It read that
the incidence half is multiplied by a stride of `n + 1` so that a vertex and a clique of the same
index cannot cancel, and explained the key as a number in base `n + 1` with the adjacency sum as
the low digit and the incidence sum as the high one. Every sentence of that is true of the KEY.
Three lines above it, the hash is defined as the key modulo `n + 1`, which keeps the low digit and
discards the high one, so the bucket was a function of the ADJACENCY ALONE. Written apart, both
statements read as correct, and the first read as a justification for the thing that broke the
second.

**So the rule the two lines have to satisfy is a joint one, and it is what the text should say
rather than describing each line on its own: the modulus must not divide the stride.** Having no
stride is the cheapest way to hold it. What it was costing, on the same graphs and for the same
merges, was 19.0 pairs tested per pivot at 140 a side against the vendored routine's 0.333 and
155.3 at 26 cubed against its 0.484; the amd section's entry 8 has the account and
`docs/DESIGN_DECISIONS.md` (2026-08-09) has why nothing here could see it.

**And the property that licenses all of it is one-directional:**

```
A[u] = A[v] and I[u] = I[v]   =>   key(u) = key(v)   =>   hash(u) = hash(v)
```

Both implications are trivial, equal multisets having equal sums and equal integers equal residues.
The converse fails freely, since `mod (n + 1)` maps a range of size `O(n^2)` onto `n + 1` values. So
collisions are common by construction, a false positive costs one exact comparison, and a false
negative cannot occur. **That is why the key can be this crude** with no avalanche and no attempt at
uniformity, which a general-purpose hash would need and this one does not.

**1. Build a key and file the vertex.** For each live `u` in `C[p]`:

```
key  = sum over live v in A[u]  of  ( v + 1 )
     + sum over c in I[u]       of  ( c + 1 ) * ( size + 1 )
hash = key % ( size + 1 )
```

then push `u` onto the chain at `hashHead[hash]`, recording the hash in `usedKeys` the first time
it is seen so the next iteration visits only the buckets this iteration touched.

The sum, the `+ 1` and the stride are the key's design, argued above. One further decision belongs
to the implementation alone.

**The chain is filled in REVERSE**, iterating `C[p]` from its last member to its first. A chain
pushed at the head comes out reversed, and the order within a bucket decides which of two
indistinguishable vertices absorbs the other, so filling forward is a tie-break change wearing a
data-structure change's clothes. It moved the permutation on four of the test graphs before the
loop was turned around. Same hazard the degree buckets carry.

**2. Compare, exactly, every pair within a bucket.** The test is

```
A[u] - {v} == A[v] - {u}    and    I[u] == I[v]
```

and it is decided by the mark scheme rather than by sorting, like every other membership test in
this code. Stamp all of `v`'s entries with a fresh tag, counting them into `sizeV`; then walk `u`'s
entries, checking each is stamped and counting them into `sizeU`. The sets are equal exactly when
every one of `u`'s entries was stamped and `sizeU == sizeV`. One pass each, no allocation, and an
early exit on the first mismatch.

The removal of `v` from `A[u]` and of `u` from `A[v]` is not a detail: indistinguishable vertices
are adjacent to each other, so without it no pair would ever match.

`mark` is sized `2 * size` for this pass, with cliques stamped at `c + size`. That is the same
separation the key's stride achieves, done a second time and for the same reason: a vertex and a
clique of the same index must not be confused.

**One elimination uses the mark-and-tag scheme in three different contexts**, which is worth listing
in one place because the array and the counter are shared and the questions are not:

```
eliminator      inClique, absorbed        membership of C[p] and of I[p]
absorption      dead_tag                  membership of the dead-clique set
hash            other                     the exact set comparison
```

**The hash's use is the unusual one, in two ways.** It stamps BOTH vertices and cliques under one
tag, which is why `mark` is `2n` long in this file and nowhere else; the eliminator gets away with a
single range because a clique's id IS a vertex and only one kind is ever live under a given tag.

And it is the only place where the mark answers EQUALITY OF TWO SETS rather than membership in one.
The scheme does not do equality directly, so it is built from two membership passes and a count:
stamp all of `v`, then walk `u` checking each entry is stamped and counting as it goes, then compare
the counts. Membership catches anything in `u` that is not in `v`; the counts catch anything in `v`
that is not in `u`. Neither test alone would do. That is the same construction the union at md2
uses, one stamp of one side and one pass over the other, put to a different question, and it is why
the exact test costs `|A| + |I|` rather than a sort.

**3. Merge.** `qg.merge(u, v)` folds `v` into `u`: `v`'s weight moves to `u`, `v`'s lists are
emptied, and `v` is marked eliminated. From then on `v` is invisible to every walk.

**The two mechanisms differ in their TARGET, and the difference is total rather than incidental.**
Both draw from the same set, `C[p]`, the clique the elimination just formed, whose members are the
vertices the pivot reached. The pivot is not one of them, since `reach` excludes the vertex it is
computed for.

- **Mass elimination always folds into the PIVOT.** Its test is that `u`'s whole remaining reach
  lies inside the new clique, which is a statement about the pivot, so the pivot is the only target
  the test can have. The merged vertex leaves the graph with the pivot in the same iteration.
- **The hash never folds into the pivot.** It compares pairs drawn from `C[p]`, and the pivot is not
  a member of its own clique, so both vertices are live and neither is the pivot. The survivor stays
  in the graph as a candidate, carrying the combined weight.

That is also why the note above about degrees applies only to the hash. After a mass merge there is
nothing to keep, the target being eliminated in the same iteration; after a hash merge the target is
still a live candidate and needs a usable degree, which is the one it already has.

A worked case, graph1 under amd1 and amd2, where the two mechanisms reach the same partition by
different routes:

```
amd1   iteration 0   pivot 3   merged -                    supervariables  {3}  {2} {1} {0}
       iteration 1   pivot 2   merged 1, 0  into 2  MASS   supervariables  {3}  {2,1,0}

amd2   iteration 0   pivot 3   merged 2     into 0  HASH   supervariables  {3}  {0,2}  {1}
       iteration 1   pivot 0   merged 1     into 0  MASS   supervariables  {3}  {0,2,1}
```

Note amd2's first iteration: the target is 0 while the pivot is 3. Two pivots either way, the same
two vertices merged away, the same fill, and the same fundamental supernode partition, one singleton
and
one group of three. amd2 moves one merge earlier and takes it out of the later iteration, with no
gain and no loss on a graph this small. graph3 is where the hash finds a pair mass elimination
cannot,
and there amd2 finishes in ten pivots against amd1's eleven.

**And no degree is recomputed afterwards**, which is worth stating because it is the first thing a
reader checks. `u` keeps the degree the bound pass wrote a few lines earlier, and it is still
correct: an external degree excludes `u`'s own supervariable, the two vertices were adjacent to
each other, and `v` leaves the graph entirely, so `u`'s reachable set is exactly what it was. What
the merge changes is `u`'s WEIGHT, and the buckets are keyed on degree. Every other member of
`C[p]` is unaffected as well, since `v` was in `C[p]` and its weight has moved to `u`, so the
weighted `|C[p]|` is unchanged and with it the middle term of their bounds. The order of the three
passes is therefore bound-then-hash and never the reverse.

**What it costs.** A key is `|A[u]| + |I[u]|` per member of `C[p]`. A comparison is the same again
per pair tested. The bucket loop is quadratic in bucket size, which is what the key's construction
works to avoid, and it is also why AMD2's tag counter has no clean quadratic bound (see item 4 of
the ordering questions in `docs/TODO.md`).

**And it is opportunistic rather than exhaustive.** The pass runs over the reached set at each
iteration, so it finds pairs that are indistinguishable AT THAT MOMENT and both present. A pair that
becomes indistinguishable later and is never in a reached set together again is not found. That is
a different kind of miss from a collision, and AMD accepts it too.

### MMD: no hash, because it is not asking the same question

`genmmd` finds SOME of the same population without any of the above, and the tempting summary, that
an exact degree hands you the answer for free, is wrong. It is worth stating carefully, because the
free-lunch reading is the one a reader arrives at.

**A union walk computes `reach(u)`, never `reach(v)`.** So it cannot on its own conclude that `u`
and `v` are indistinguishable: that needs both sets, and only one of them is being built.

What it CAN conclude is narrower, and the code says exactly when. Our own `Mmd2.cpp`, in the middle
of computing `reach(u)`:

```cpp
if (mark[v] == cliqueTag) {                                   // v is in the new element too
    if (qg.adjacencySize(v) + qg.incidenceSize(v) - 1 == 1)    // and v has no OTHER source
        qg.merge(u, v);                                        // identical reach
    else
        outmatched[v] = 1;                                     // v reaches more, never minimal first
}
```

The walk found `v` through one of `u`'s own sources, so `v` is in the new element AND in that same
other source. The second line then checks that `v` has no source `u` lacks. **Only then is
`reach(v)` fully determined by what the walk has already seen**, and only then can the merge be
made without computing anything more.

And this runs only for a q2h vertex, one with exactly two sources, which is what put `u` on that
list. So the population is: pairs where both vertices have exactly two sources and both sources
coincide. `genmmd` does the same thing by a different route, stashing each reached vertex's pruned
adjacency count as `fwd[rn] = nq + 1` and routing the `nq == 1` cases into its own q2h list.

**So MMD is not avoiding the hash by being cleverer. It is avoiding it by not asking.**

```
MMD    pairs where both have exactly two sources, both shared
       free, a by-product of a union it must compute anyway
       cannot generalize: the walk gives reach(u) and never reach(v)

AMD    any pair in C[p]
       a key per vertex, an exact comparison per colliding pair
       must go looking, because a bound never opens a clique
```

**AMD's mechanism does strictly more work because it answers a strictly harder question.** Its
degree is a bound computed by decomposition, so it never opens a clique and never sees which other
vertices share `u`'s sources; it has to go looking, looking at every pair is quadratic, and the
filter is the hash. That AMD finds strictly more is the capability difference, and the cost of the
hash is what that capability costs.

**Which is the same trade as the degree itself, one level up.** MMD pays a full union and gets
detection thrown in; AMD pays a cheap bound and then pays again for detection. Whether the total
comes out ahead is not decidable by reasoning, and the measurement has been unkind to the
intuition: our AMD2 carries both of AMD's mechanisms, fires the hash 2488 times across the test set
against aggressive absorption's 1, and still fills 7 percent worse and orders 65 percent slower
than our AMD1, which has neither. Every matrix behind those numbers is a grid, which is where a
tie-break decides almost every pick.

Section 5.5 of `archive/sparse_factorization.md` states the same fork in one place, since a reader
meeting supervariables there would otherwise take hashing to be the definition rather than one of
two routes.

## What each branch has to remember, and why an encoding does not transfer, 2026-08-17

The section above compares the two DETECTION MECHANISMS. This one compares what each branch has to
REMEMBER about a vertex while it works, which is a different question and is the one that decides
whether an encoding found on one branch can be carried to the other. It was written while planning
the port of `Amd3B`'s five array folds, three of which live in the shared quotient graph and so
reach every driver whether or not anyone aims them.

Every walk in either branch asks three things of a vertex it meets: **is it dead**, **have I already
seen it this step**, and **what does it weigh**. Both branches ask all three. They differ in where
the answers come from, and the differences are not stylistic.

### The three ways a vertex leaves the graph, and the asymmetry is the whole finding

|  | weight | mark | still named by other lists |
|---|---|---|---|
| `merge`, a live hash merge | zeroed | `GONE` | yes, left exactly where it lies |
| `massEliminate` | zeroed | `GONE` | no, and it was in the pivot's clique alone |
| `number`, the mmd prepass | UNCHANGED, and one for a fresh vertex | `GONE` | yes, everywhere |

**Amd has only the first two.** So on that branch a zero weight and death are the same thing, and
`Amd.cpp` reads exactly that: `nvj = Nv [j] ; if (nvj > 0)` answers all three questions off one
load, positive being the weight, negative meaning already taken into the new clique, zero meaning
absorbed. There is no fourth state to represent because there is no fourth way to die.

**Mmd has the third, and it is deliberate rather than incidental.** A prepass vertex keeps its
weight precisely so that its neighbors' degrees still count it, and it stays in every list that
named it. Under the sign encoding it reads as live and unseen, so a reachable-set walk emits it. The
substitution was tried and reverted on 2026-08-08, where it produced 201 entries for 200 vertices on
a random `mmd2` pattern. So on the mmd branch the weight is a PARTIAL flag and the mark is what
carries liveness.

genmmd holds both at once rather than choosing, and that is the shape of the answer rather than a
way around the problem. Its `marker[v] = maxint` is the general liveness test, and `qsize[nd] != 0`
is used only inside element walks, where a prepass vertex cannot appear because the mark has already
kept it out of every clique. Two encodings covering two regions, where amd needs one covering both.

### The clique side, where the asymmetry runs the other way

Amd stamps clique ids and mmd never does, which is why `mMark` is 2n on one branch and n on the
other, and why the width is a constructor argument rather than always the larger.

The reason is the section above: amd's exact test is set equality over `I[u]` and `I[v]`, two
arbitrary lists, and set equality is done by stamping one side and walking the other, so the stamp
has to reach the clique id space. Mmd asks nothing of that shape. Its merges fall out of a walk it
is already making, testing whether a vertex under the cursor is also in the new element, which the
vertex mark answers. Neither `w[c] == 0` for a dead clique nor `w[c] >= wflg` for one seen this step
is a mark at all; both are the tagged W, which is amd's own and which mmd does not carry.

So the two branches are each frugal where the other is not, and neither pays for what it does not
ask.

### The lifetime, which is where the port gets hard

The negation is not a flag written once and left. It IS the reachable set's insertion operation:
`nv > 0` means "not yet emitted" and `mWeight[v] = -nv` is the emit, so the array is deliberately in
a non-canonical state for the length of an elimination and something has to restore it.

- **Amd negates in `reachableSet`, called ONCE PER PIVOT, and restores in the bound pass it already
  makes.** `Amd.cpp` does the same, `Nv[i] = -nvi` under CONSTRUCT NEW ELEMENT and `Nv[i] = nvi`
  under RESTORE DEGREE LISTS. Neither is a pass of its own.
- **Mmd calls `reachableSize` and `reachableSetWeight` PER VERTEX in the refresh**, and has no bound
  pass to hide a restore in. Each call must therefore leave the array as it found it, which is a
  per-call lifetime rather than a per-pivot one, and it is the design problem `Mmd3C` exists to
  work out.

### Where the sign can serve on the mmd branch, and what it costs there

Two halves, pulling in opposite directions. Neither is measured.

**The good half: a prepass vertex can never appear in a clique member list.** The prepass runs to
completion before the main loop, so no clique exists at the moment `number` is called, and every
clique formed afterwards comes from a reach set that already skipped numbered vertices through the
mark. The exclusion is structural rather than probabilistic. `A[u]` is the opposite case, since
`number` leaves the lists alone: a numbered vertex sits in the adjacency of every neighbor for the
rest of the run.

So the split follows the SOURCE KIND, and `reachableSet` and `reachableSetWeight` already write the
two loops separately:

```
clique walk, C[c] for c in I[u]    sign carries liveness, dedup and value      one array
adjacency walk, A[u]               sign carries dedup and value, and the
                                   mark is still needed for liveness          two arrays
```

Dedup has to be one mechanism across both loops, or a vertex reached first through `A[u]` and then
through a clique is emitted twice, so the negation goes in both and only the liveness test differs.
The clique walk is the half that matters: by the conservation lemma `|A[u]| + |I[u]|` is bounded by
u's original column for the whole run, so the adjacency half shrinks monotonically under the prune,
while each clique in `I[u]` opens a `C[c]` that grows with fill.

**IMPLEMENTED IN `Mmd3C` ON 2026-08-17, and the split above is exactly what landed.** The clique
loops of `reachableSet` read one array, the adjacency loops read two, `beginElimination` loses its
stamping pass, and the prune answers three of its four questions from the sign. `reachableSetWeight`
was DECLINED rather than blocked, for the reason in the next paragraph. Permutations identical to
`Mmd3` on 24 cases and across all nine drivers under `make digest`.

**The bad half: the tag scheme invalidates in O(1) and the sign scheme cannot.** `++mTag` retires
every mark in the array at once, which is exactly why `reachableSetWeight` can count a reach without
materializing it and without cleaning up after itself. A negation has no such trick. Every one must
be undone individually, and these two functions keep no record of what they negated, having been
written deliberately not to build a list.

Amd never meets this. `reachableSet` runs once per pivot and writes its result into the arena, and
the bound pass walks exactly that result and restores as it goes, so the restore rides in a
traversal that already exists. Mmd's refresh calls `reachableSize` and `reachableSetWeight` PER
VERTEX and walks the result again nowhere, so there a restore is a new traversal rather than a free
rider.

**What landed instead, and it is cheaper than either.** `massEliminate` already walks C[pivot],
which is exactly the set `reachableSet` negated, so the restore rides there and costs no pass at
all. That is the mmd analogue of amd's bound pass, and it is why the fold was affordable on the
elimination path. It does NOT extend to `reachableSetWeight`, which is called per refreshed vertex
and whose result nothing walks again; that one is left on the mark.

**SUPERSEDED 2026-08-24 ON WHERE THE RESTORE LIVES, and the two paragraphs above stand as the
account of the day.** Both branches now restore AT THE END OF THE PRUNE, which is the last reader of
the negated form in either: the prune's `mWeight[v] <= 0` is what consumed the mark, absorption
never touches a weight, mass elimination's merge test is structural, and everything downstream takes
a magnitude through `weight()`. So the restore no longer rides in `massEliminate` on one branch and
in the driver's bound pass on the other; it is one place, unconditional, in both classes.

What that removed: `restoreWeight` and `restorePivotWeight` from both class surfaces and their four
call sites in the amd drivers, and `mLateMassElimination`'s SECOND job. That flag had been deciding
both "do not merge here" and "do not restore here", and only the first is what it is named for.

It is still not a saving. The restore rode on an existing walk before and is its own short loop over
C[pivot] now, so this trades a free rider for a pass in exchange for one lifetime instead of two.
The gate was that nothing moves: all 365 digests identical, both alignment checks 38 of 38.

The cheapest form for `reachableSetWeight`, if it is ever wanted, is a member scratch holding the
emitted vertices, cleared per call so its capacity survives, walked afterwards to restore. Per member that gives up one scattered mark load and one
scattered mark store, and takes on one contiguous scratch store and one scattered weight store in
the restore, against a gain of one fewer array touched in the clique walk. The scattered count comes
out about even and lands on one array instead of two, which is the reason to expect anything at all.
Whether it nets out is a measurement, and this is precisely the shape the tree has misjudged three
times: a schedule change that saves visits and adds a pass. See the footprint trade in
`docs/DESIGN_DECISIONS.md` (2026-08-16).

### What this decides about the five folds

| fold | reaches mmd | why |
|---|---|---|
| 1, weight sign as the membership mark | yes, with a caveat | shared class, but `number` breaks the premise |
| 3, `eliminated()` off a zero weight | yes, with the same caveat | same premise, same break |
| 4, mass elimination merges before it compacts | yes | independent of the encodings above |
| 2, the restore rides in the bound pass | no | mmd has no bound pass |
| 5, detection stamps into `w` | no | it retires the clique half of the mark, which mmd never had |

Two of the five are amd-only by construction rather than by difficulty, which is worth knowing
before either is attempted. Folds 1 and 3 are the ones with real work in them, and the work is not
the encoding but finding which regions of the mmd walks a prepass vertex provably cannot reach.

**So the success condition for `Mmd3C` is not that the mark array goes.** It probably cannot go, for
the reason in the table above. It is that the hot walks stop reading it, which is what genmmd
achieves while still carrying `marker`, and which is where the measured cost is: two scattered loads
per clique member where the vendored routines make one.

## Three clique layouts, and one axis: the tax is the price of the constraint, 2026-08-17

A quotient graph has to put C[c] somewhere, and the three codes in this tree answer differently.
With the array encodings finally equal across all six drivers, the difference between a driver and
its B or C sibling is the layout alone, so the columns can be read as the price of a layout.

**THE WHOLE THING IN THREE LINES, and the rest of the section is the detail behind them.** Written
2026-08-18, and each line says where the cliques sit, what is paid, and how much space it takes:

```
genmmd   the original space only, nothing ever shifted        pays in chaining        no extra space
AMD_2    the original space where the pivot allows, else a    pays in compaction      extra space,
         free area; everything shifts at a compaction         but rarely              which is what
                                                                                      makes it rare
ours     a region of its own, nothing ever shifted            pays nothing            unbounded
```

**Placement follows from that, and only `AMD_2` is a mixture.** genmmd is vertex-ANCHORED rather
than uniformly vertex-ordered: a clique starts in its own pivot's segment and, when it does not fit,
continues into the segments of the cliques it absorbed, which sit at other ids. About a third
overflow on grids, so a third are scattered across the ids they came from. `AMD_2` puts 62 to 68 per
cent at the pivot's own id and the rest in elimination order at the cursor. Ours is uniformly
elimination order. **So a locality claim about `AMD_2` has to say which population it is about**,
where the same claim about either of the others does not.

**And the extra space is a knob rather than a constant.** At `AMD_2`'s 20 per cent the collector
runs about once for a whole ordering; at about 1.25 nnz it runs on nearly every pivot. The cliff
between those is measured below, and the space beneath it is what genmmd buys by chaining.

**THE CLASSIFICATION THAT MAKES SENSE OF THEM IS ONE AXIS, NOT THREE SCHEMES.** Storage is
constrained or it is not, and bookkeeping is what a constraint costs:

```
   no constraint      and no bookkeeping           ours
   a loose constraint and occasional bookkeeping   AMD_2
   a tight constraint and constant bookkeeping     genmmd
```

Chaining and compaction are then the same kind of thing at two points on that axis rather than two
different taxes to be compared. That framing is Florin's and it is better than the one this section
first carried, which treated them as rival mechanisms; the measurement below is what settled it.

### The axis is really two, and they trade against each other, 2026-08-18

**MEMORY AND MACHINERY, AND MORE OF ONE BUYS LESS OF THE OTHER.** The one-axis statement above is
the diagonal of a two-axis picture: how much space beyond the pattern are we willing to spend, and
how much bookkeeping are we willing to run given that space. The extremes name themselves.

- **No extra space at all.** The machinery has to be constant, because every clique must be fitted
  into ground already occupied. genmmd chains; `AMD_2` compacts. These are the only two answers
  either code gives, and they are the expensive end by construction.
- **Space is free.** No machinery is worth writing. Append and never look back, which is what ours
  does. There is nothing to reclaim because nothing needs reclaiming.
- **Between them**, and this is where compaction earns its keep: with real but bounded headroom the
  collector is IMPLEMENTED AND RARELY RUN. At `AMD_2`'s twenty per cent it fires about once for a
  whole ordering, and the elbow-room measurement below shows it degrading gracefully to `sum(Len)`
  rather than failing. There is no guarantee it will not be needed, so it must be ready; there is
  every expectation it will not be needed often.

**Chaining does not have a middle.** Its cost is a link test on EVERY READ, paid whether or not the
space is tight, so extra headroom buys it nothing. Compaction's cost is a sweep that extra headroom
makes rarer. So chaining is worth considering only at zero extra memory, and compaction across the
whole range from zero upward. That is a pragmatic narrowing rather than a claim that chaining
cannot work with slack.

**What the constrained end is FOR, given that the free end is simpler.** Prediction. A run under a
memory bound either completes or reports what it would have needed, which says something useful
about whether the factorization is feasible at all. The free-memory version can only fail later and
larger: if the arena does not fit, the factor would not have fitted either, so failing now or
failing then is the same failure with less information.

### An algorithm is never forbidden by a layout, only priced

**A LAYOUT CONSTRAINS IMPLEMENTATIONS, NEVER ALGORITHMS.** Any minimum-degree variant can be run on
any of these layouts; enough indirection will always express it. What a layout can do is make a
particular implementation of a step unavailable, and then the algorithm is reached by a more
expensive route instead. The coupling that remains is a PRICE, not a prohibition, and the pragmatic
question is whether it is worth paying.

The case that made this concrete, 2026-08-18. With the run laid out incidence-first, `AMD_2`'s way,
there is no free slot in which to append a new clique at the BACK of `I[u]`, which is genmmd's
convention. That does not forbid the convention: a shift, a copy aside, or one spare word per vertex
each deliver it. `AMD_2` declines all three and inserts by ROTATION instead,
`Iw[pn] = Iw[p3]; Iw[p3] = Iw[p1]; Iw[p1] = me`, which puts the clique at the FRONT and costs
nothing. So the convention changed to suit the layout, and the tie-breaking changed with it. That is
the shape of every such coupling here: the layout did not decide what could be computed, it decided
what was cheap, and the cheap thing was taken.

**Which is why `Mmd3C` has not followed `Amd3B`** through the run-order flip. It is mmd on
`AMD_2`'s layout and it uses the back-of-list convention, so following would mean either paying for
one of the three routes above or changing the convention and moving the permutation. Neither is
wrong; the choice has simply not been made.

### The three stores, measured on 246 real matrices, 2026-08-19

The layout program's first real answer, from `benchmarks/matrices`. Each sibling is its branch's
driver on the vendored routine's clique store instead of our arena, so the ratio in the last column
is the store and nothing else in the algorithm:

| | store | against our arena |
|---|---|---:|
| `AMD3B` | `AMD_2`'s unified workspace, bounded at a fifth over, compacted | **0.95** |
| `MMD3B` | genmmd's dead segments, chained, bounded at exactly nnz | **1.56** |

**The two disagree, and the mechanism says why.** A compaction is a rare sweep: at a fifth of elbow
room it runs about once for a whole ordering, so its cost amortises to nearly nothing and the pool
beats the arena, on 224 of 242 matrices. A chain is a LINK TEST ON EVERY READ of every clique,
forever, and no amount of headroom reduces it, because chaining exists precisely so that none is
needed. In the band holding 96 per cent of the mmd time, `MMD3` reads 0.69 against genmmd and
`MMD3B` reads 1.09.

**The grid ladders understate chaining badly**, 1.1 to 1.3 there against 1.56 here, because grids
have short cliques and few links to follow. That is worth remembering generally: a scheme whose
cost is per-link is invisible on inputs whose lists are short.

**And read both as directions rather than coefficients.** Each sibling differs from its driver in
more than the store, so the store is the largest of those differences and not the only one.

### genmmd never shortens a clique, and `AMD_2` always does, 2026-08-18

A contraction removes members from a live clique. What the three codes do with the SPACE those
members occupied is a separate question, and they answer it differently:

- **genmmd leaves them where they lie.** `mmdelm`'s mass elimination is
  `qsize [md] += qsize [rn] ; qsize [rn] = 0` and it does NOT rewrite the list; the merged vertex
  keeps its slot and every later reader skips it on `qsize [nb] != 0`. Same for the further
  absorption inside `mmdupd`. So a clique is placed once and never shortened.
- **`AMD_2` compacts.** Its RESTORE DEGREE LISTS pass writes survivors back with `Iw [p++] = i`,
  sets `Len [me] = p - pme1`, and hands the tail back with `pfree = p`.
- **Ours compact**, in every class.

**`Mmd3B` was corrected to match genmmd**, its compaction removed. It cost a pass genmmd never pays
and saved the skipped entries genmmd pays on every later read; measured afterwards, its column
against genmmd is unchanged within noise, so the two costs are about equal. NEITHER IS VISIBLE IN A
PERMUTATION, a skipped member and an absent one reading alike, so no check in this tree would have
caught it. It was found by reading `mmdelm` while answering a question about clique lifetimes.

**What it costs is space, and it is genmmd's space:** a chained clique keeps its dead members for as
long as it lives, so live clique STORAGE there is strictly above what the compacting classes report
for the same ordering. The live MEMBER count is unaffected and is the same for all three, being a
property of the algorithm.

### The life of a clique: birth, contraction, death, 2026-08-18

**THREE EVENTS, AND ONLY THREE.** They are worth naming because the counter that reports peak live
clique members has to see every one of them, and because two of the five sites below were found by
reading rather than by any check the tree runs.

**BIRTH.** One site, `beginElimination`, in every class. The reach of the pivot is placed and the
pivot's own descriptor becomes the clique's. Nothing else creates a clique, which is why the
running maximum is taken here alone: no other event can raise the live total.

**CONTRACTION.** The clique survives and loses members.

| cause | amd | mmd |
|---|---|---|
| mass elimination, the merged leave C[pivot] | yes | yes |
| supervariable detection, the hash-absorbed leave | yes | no |

**DEATH.** The clique ceases to exist.

| cause | amd | mmd |
|---|---|---|
| absorbed into the new clique, being in I[pivot] | yes | yes |
| absorbed aggressively, its external degree reaching zero | yes | no |

So mmd has one of each and amd has one birth, two contractions and two deaths. The first death is
the common case and the reason cliques do not accumulate: every clique the pivot belonged to is
subsumed by the one the pivot forms.

**A CONTRACTION IS NOT A DEATH, and the distinction is not pedantry.** A contracted clique keeps its
identity, its descriptor and its remaining members, and is still read. Conflating the two was what
made the first version of the counter wrong.

**TWO THINGS THAT LOOK LIKE EVENTS AND ARE NOT.** `merge` zeroes a vertex's adjacency length, but
that vertex never formed a clique, so the length is A[v]'s and nothing is subtracted; feeding it to
a clique counter corrupts the total. And in the CHAINED layout the two contractions change the count
without touching storage at all: genmmd leaves the dead members where they lie and skips them on
read, so live members fall while occupied space does not.

### What puts each scheme where it sits, and it is ONE question

**MAY A CLIQUE BE SPLIT?** Everything else follows from the answer, and the three sections above
describe the consequences without ever naming the cause. Written out 2026-08-18.

**genmmd says yes.** A clique may run across several blocks, joined by a negative link at each
block's end, so it can be laid into the dead block of the pivot that formed it and continue into
the blocks of the cliques it absorbed. **Every dead block is therefore usable for every clique**,
which is what buys the exact `nnz` bound. The price is that no walk of a clique is a straight run:
the link test is paid on every read, forever.

**`AMD_2` says no.** A clique is one contiguous block, so a dead block is usable only when the
whole clique fits inside it. Exactly one case qualifies, and it is a property of the PIVOT rather
than of the space: with `I[pivot]` empty the reach is a subset of `A[pivot]`, so the clique is
compacted where `A[pivot]` already sits. Every other clique goes to a free area past the runs.
**So the contiguity requirement is what forces the collector**: with chaining unavailable, dead
space in the middle cannot be reached at all, and the only way to get it back is to move
everything, ORIGINAL SEGMENTS INCLUDED. The elbow room and the sweep are consequences of the
constraint, not a scheme chosen beside it.

**Ours says no, and does not try the dead block either.** Every clique is appended to a region of
its own, so no original space is reused and the region grows toward nnz(L).

Two things follow that are worth having in one place:

- **`AMD_2` cannot chain, so it must collect.** These are not two designs that happened to be
  paired. Given contiguity, a collector is the only remaining answer, which is why the headroom
  cliff below is a property of the constraint rather than a tuning accident.
- **`AMD_2`'s cliques live in TWO POPULATIONS.** Roughly two thirds sit at the pivot's original
  address, which is where the vertices that read them next already are, and one third sits at the
  cursor in elimination order. genmmd's are uniformly at vertex-id placement and ours are uniformly
  in elimination order, so `AMD_2` is the only one of the three whose locality is a mixture, and
  any locality argument about it has to say which population it is about. The share is measured:
  62 to 68 per cent in place on grids, rising with n, because the qualifying condition is an
  early-run one and most pivots are eliminated early.

```
                     extra space      what it pays instead        does that cost grow?
ours, an arena       unbounded        nothing                     n/a
AMD_2, a pool        45% on grids     compacting when it fills    NO, a handful of sweeps
genmmd, segments     none             following chains            YES, with every clique read
```

**Ours** is a separate append-only region in elimination order. Nothing is ever reclaimed, so there
is no bookkeeping of any kind and no bound either.

**genmmd's** puts a clique in the dead pivot's own block and chains onward through the blocks of the
cliques it absorbed, following a negative link at a segment's end. It fits in exactly `nnz`, and the
reason is the conservation argument the prune already relies on: a source is destroyed for each one
created.

**`AMD_2`'s** keeps variable runs and clique blocks in ONE workspace with elbow room, builds each
new clique at a free cursor, leaves absorbed space dead, and compacts when the cursor reaches the
end.

### The asymmetry, which is the point of this section

**Chaining is a PER-ACCESS cost.** Every walk of a clique runs the link test, and a clique is walked
many times over its life, so the total scales with how often cliques are read, which grows with fill.

**Compaction is a ONE-OFF cost AT THE HEADROOM `AMD_2` CHOOSES.** `Mmd3C` compacts 4 times, flat,
from 1,600 vertices to 160,000. Four sweeps of the pool for a whole run, against a number of clique
walks that grows superlinearly, so the cost per pivot tends to zero.

### But that is one point on a curve, and the curve is a CLIFF

The pool's headroom is a knob, so the price of the constraint can be measured directly rather than
argued. `Mmd3C` at 300 squared, headroom varied and nothing else changed:

```
pool / nnz     1.05    1.10    1.20    1.25    1.26   1.27   1.28   1.30   1.45   2.00
compactions   49704   48583   46340   44749     58     29     20     13      4      1
```

**Four orders of magnitude across one percent of headroom.** Below about 1.25 the collector runs on
essentially every pivot and the count scales with n; above about 1.26 it runs a handful of times and
stops scaling. Cubes put the threshold in the same place, 7 compactions at 1.28 and 3 at 1.45, so
the count is not a constant of the algorithm but the CLIFF LOCATION is stable across families.

So the tax is not "more space, proportionally less work". It is: buy enough headroom to hold the
working set and the bookkeeping nearly vanishes; buy slightly less and you pay it per pivot.
`AMD_2`'s `nnz + nnz/5 + n` lands at 1.45 on grids, comfortably past the cliff, and the "4 sweeps,
flat" figure above was taken on the safe side of it without knowing the cliff was there.

### What that says about genmmd, and it is the interesting part

genmmd sits at 1.00, which is BELOW the cliff by a wide margin: a compacting scheme at that headroom
would compact on every pivot and be unusable. So chaining is not an alternative tax at the same
price. **It is what makes the tightest point reachable at all**, buying the last quarter of space
that compaction cannot. The per-access cost is what that space costs.

"genmmd pays chaining, `AMD_2` pays collection, call it even" is therefore the wrong reading twice
over: the costs have different shapes, and the two schemes are not even trying to buy the same
thing.

### What the timings say, and they are consistent with that

Geometric means over the square and cubic ladders, each pair differing in layout alone:

```
MMD3B / MMD3   genmmd's segments COST     1.079 square    1.288 cubic
AMD3B / AMD3   AMD_2's pool EARNS         0.914 square    0.855 at large n
```

The cubic figure for genmmd's scheme being much the worse fits a per-access cost: cubes carry far
more fill and so far more clique reads per pivot. The pool's advantage GROWING with n fits the
removal of one.

**These are wall clock on alpamayo and should be treated as provisional** until they have been seen
across more runs. The vendored columns in the same tables move several percent between invocations
of unchanged code.

**AND THEY PREDATE THE FAITHFULNESS FIXES BELOW.** Both pool figures were measured on a port that
was missing `AMD_2`'s in-place construction and its reclaim, so they are the price of a variant.
`AMD3B / AMD3` has since moved to 0.885; the mmd figure has not been re-measured. The headroom
curve further down has the same problem: it was taken on a scheme that consumed the pool about three
times faster than `AMD_2`'s does, so the cliff it locates sits further left than it should.

### The space figure is 45 percent on grids, not 20, and it is space rather than time

`AMD_2` sizes its workspace `nnz + nnz/5 + n`. On a 5-point grid `nnz(A+I)` is about 4n, so the `+n`
term adds another quarter on top of the fifth and the pool comes out at 1.451, 1.453, 1.452, 1.451,
1.451, 1.451 over sides 40 to 400: flat, and nowhere near the fifth the formula suggests at a
glance. On a denser matrix the `+n` term shrinks against `nnz` and the ratio tends toward 1.2.

Worth keeping distinct from the timing figures above. A layout that costs 45 percent more space and
saves 9 to 15 percent of time is a trade; one that costs 45 percent more space and saves nothing is
not, and only the measurement says which.

**And 1.45 is not a number to tune downward casually.** The cliff sits at about 1.26 on both
families measured, so there is real slack between it and `AMD_2`'s choice, but the penalty for
crossing is four orders of magnitude rather than a few percent. Anything that trims the headroom
wants the compaction count watched, which is what `gMmd3CCompactions` and `gAmd3BCompactions` are
exported for.

### Porting a layout faithfully is harder than it looks, 2026-08-17

The B and C layers exist to price a layout, so a port that drops part of one prices something else.
Three divergences from `AMD_2`'s storage were found in `Amd3B` by reading it against
`private/Amd.cpp` on 2026-08-17, two days after the storage was ported into it. **Finding one was
reason to look for more, and there were two more.**

**The failure was not doing that reading at the time.** The storage was ported by writing something
with the same shape rather than by translating `AMD_2` line by line, which is the one thing
"port, don't rewrite" exists to prevent, and the figures the file produced in between were quoted as
the price of `AMD_2`'s layout when they were the price of something else. An audit against the
vendored source belongs in the port, not after the numbers look interesting.

**What is NOT a divergence, checked the same day:** the vertex list. `AMD_2` splits one run into
elements then variables, `Iw [Pe[i] .. Pe[i]+Elen[i]-1]` and `Iw [Pe[i]+Elen[i] .. Pe[i]+Len[i]-1]`,
which is our run with `adjacencySize` and `incidenceSize` in the opposite order, and
`mVendoredListOrder` is what selects that order. genmmd is the one that differs, keeping variables
and elements MIXED in one list and telling them apart by the sign of `fwd`; that is a `Mmd3B`
question and is open rather than closed.

**IN-PLACE CONSTRUCTION, `AMD_2`'s `if (elenme == 0)`.** A pivot with no elements has a reach that
is a SUBSET of `A[pivot]`, so the clique is compacted where the adjacency already sits and the pool
is not touched. We built at the cursor every time. Measured, **62 to 68 percent of eliminations
qualify on grids, rising with n**, so this was not a corner case: we were putting three times as
much through the pool as `AMD_2` does, and throwing away the locality of a clique landing where its
pivot's adjacency was.

**THE RECLAIM, `AMD_2`'s `if (elenme != 0) pfree = p`.** After supervariable detection removes
newly nonprincipal variables the element is shorter, and the cursor is pulled back. Ours only ever
advanced. It rests on the clique still being the last block in the pool, which is true because
nothing writes to the pool between construction and trimming; that is now an `assert` rather than a
comment, and it is checked on all 73 digest grids in the root build.

Both landed on 2026-08-17. **Compactions fell from 4 per run to 1 on squares and 2 on cubes** in
both files, and `AMD3B / AMD3` moved from 0.929 to 0.885 as a geometric mean, evenly at both ends
of the ladder. A faithful port coming out FASTER is the expected direction here and is the thing
most likely to be accepted without checking; what says faithful is `make amdorder`, the assert, and
the compaction count, not the clock.

### The third divergence, which stays, and what it would cost to close

`AMD_2` tests for room PER ENTRY inside the construction, `if (pfree >= iwlen) garbage_collection`,
and when it collects mid-element it carries the half-built element along: it slides every live block
down, moves the partial element separately, and resets `pme1` and its read pointers. We test ONCE,
before the walk, for the worst case, room for a reach of n.

**The reason is pointers, and it is a trade we took deliberately elsewhere.** `AMD_2` walks in
INDICES, so after a collection it can recompute where it was and resume. Our walks hold
`const std::int32_t*` into the pool, and the collector moves every block, so those pointers would
dangle. The pointers are not an accident: hoisting one out of a loop header was worth 277 ms of an
8.53 s run, which is why the walks are written that way. So the trade is a faster inner loop against
the ability to relocate mid-walk.

**It can only cost us.** Reserving for n when the reach is usually far smaller means we collect
where `AMD_2` would not, never the reverse, so the column is pessimistic rather than flattering. It
does not touch locality: the collector produces the same compacted layout either way. With
compactions now at 1 or 2 per run the window is small.

Three ways to close it, if it is ever worth closing. Rework the walks to index arithmetic, which is
faithful and touches the hottest loop in the file. Keep pointers and re-derive them after a
collection, which needs the collector to report where the partial element moved. Or leave it,
documented, which is where it stands.

### One thing the taxes do not explain, and it is a candidate for the rest

The pool is not merely bounded, it is COMPACT: variable runs and clique blocks share one region,
where ours keeps two and the second only grows. A walk that reads a vertex's run and then the
cliques it names stays inside one region in the pool and crosses between two in ours. That is a
locality effect independent of either tax, and it is the natural candidate for the part of the amd
gap that counting could not find: cachegrind puts `Amd3B` at 10 percent MORE simulated L1 read
misses than `Amd3` at both sizes measured, flat, while the clock has it 9 to 15 percent faster. A
single-level L1 model with no TLB and no prefetcher is the wrong instrument for that difference, and
this is the shape of what it is missing.

**`Amd3C` would test it.** amd on genmmd's segments is the fourth cell, and if the per-access cost
shows up there too, the account above holds on both algorithms rather than one.

## What the sandbox can and cannot answer, 2026-08-16

Established the hard way, by nearly discarding a fold that turned out to be worth 30 percent.

**Cachegrind's counts travel. The sandbox's wall clock does not.** Instructions, data references,
branch counts and simulated misses are machine-independent and were right every time. Wall clock on
the x86 sandbox disagreed with alpamayo in DIRECTION, not merely in magnitude, for any change that
traded one counter against another:

```
first version of the weight-sign fold, +1.4% instructions for -6.3% reads
  x86 sandbox    1.18  1.34  1.27  1.22  1.31    slower at every size
  alpamayo       0.95  0.94  0.97  0.92  0.91    faster at every size
```

x86 charged for the instructions; the M-series paid for the reads. **So a fold that trades counters
must be timed on alpamayo, and a sandbox timing that contradicts a counter improvement is not
evidence against it.** What the sandbox is good for is the counts themselves, and for anything that
improves every counter at once.

**A second use, cheap and often decisive: profile per function and per line.** `cg_annotate` on a
clean harness put 17.3 percent of our instructions in inlined `stl_vector.h` accessors, 65.8 M
against 236 M of algorithm, where the vendored routine runs 337 M of algorithm and no accessors.
That is the container layer measured directly, and it is much larger than the 1.5 percent the
totals suggest, two effects having been cancelling. The same annotation ranked source lines by D1
read misses and put the prune's `incidence[i]` at 12.5 percent of the whole run, which is how the
next fold was chosen.

**One harness warning, since it cost a reading.** Build the test graph without allocating per
vertex. The first version used a `std::vector` per vertex and put 17 percent of the profile in
`malloc`, larger than most of what the profile was being asked about.

## The padded copy, and why intervention beat three inferences, 2026-08-16

The scaling ladder showed the vendored AMD costing more per vertex at power-of-two grid sides while
genmmd, `MMD3` and `AMD3` were all smooth. Three accounts of that were reasoned out and written
down; all three were wrong, and the fourth attempt was not an account at all but a change.

**What failed.** Stride aliasing in a single array indexed by vertex id, refuted by genmmd having
the same storage shape and not zigzagging. Our clique arena buying immunity, withdrawn with it, and
unsupported anyway since nothing compares our arena against `AMD_2`'s storage. The hash modulus
being a power of two, refuted by counting: `AMD_2` tests 0.33 pairs per pivot at every size in both
series.

**What worked, in one run.** Cachegrind on the vendored routine alone, three sizes back to back in
one invocation, showed instructions and data reads per vertex FLAT to a tenth of a percent while D1
read misses jumped 2.6x at the aligned size. That said the cause was addressing rather than work but
not WHICH addresses. So: a copy of `private/Amd.cpp` with sixteen ints of padding inserted between
the six arrays `AMD_1` carves out of `S`, generated by a `sed` and nothing else, changing addresses
and nothing about the algorithm. Same permutations. **56 percent of the read misses at 512 squared
gone, and nothing changed at 400.**

**The technique generalises and is cheap.** A hypothesis about layout can usually be tested by
perturbing the layout rather than by explaining it, and a perturbation that leaves the output
identical is self-checking: if the permutation moves, the probe is wrong and the reading is void.
This one took one command and settled a question three careful arguments had got wrong.

**AND THE PERTURBATION HAS TO BE BIG ENOUGH TO PERTURB, 2026-08-17.** The same technique applied to
`Mmd3C`, whose 200^2 column read 1.28x `MMD3` with 199 and 201 flat, padded by sixteen ints and
moved NOTHING. That was read as evidence against data placement and was evidence of nothing: at 200
a side those arrays are 160,000 bytes, forty pages, and 160,064 rounds to the same forty, so the
allocator very probably returned the same addresses and no intervention took place. Padded by 1024
ints, exactly one page, the ratio fell to 0.99 with the neighbours unchanged and the permutations
byte-identical.

Sixteen ints was the right size for the six arrays `AMD_1` carves out of one workspace and useless
for separate vectors past the page threshold. **A perturbation that does not perturb is
INCONCLUSIVE, not negative**, and the way to tell them apart is to know what granularity the
allocator works at before choosing the amount. Two rounds went on other hypotheses first.

It also refutes this section's own conclusion that our separate allocations buy immunity: they do
not, being page aligned and page rounded once large. See `docs/DESIGN_DECISIONS.md` (2026-08-17).

**And cachegrind now runs in the sandbox**, installed 2026-08-16. Instruction counts and data
references are exact and machine-independent, so this whole class of question no longer needs
alpamayo. Only wall-clock confirmation does. The standing caution about simulated MISS counts
shifting between invocations still applies: compare builds back to back inside one run.

## The differential against AMD_2, and what counting settled, 2026-08-16

The amd branch's 2D ratio against the vendored routine had been RISING with n while `MMD3` over the
same quotient graph was flat. Five folds moved the column down and left the slope alone, so the
cause was not any of them, and reasoning had produced three wrong hypotheses in a row. This is the
instrument that settled it.

### What it is

Two counting copies, both generated by anchor-asserting scripts so a moved source fails the
generation rather than counting the wrong loop: one of `private/Amd.cpp`, one of the driver. The
counters carry the SAME NAMES over the passes that correspond, and everything is normalised per
pivot. The two codes produce the same permutation, so the pivot counts agree exactly and the
normalisation is honest.

Getting the correspondence right is the whole design, and it is not obvious in two places. `scan1`
is a pass of its own in `AMD_2` and is folded into the prune here, so the count is the prune's
incidence entries. `member` is counted per VISIT rather than per pivot, both sides walking `C[p]`
twice, once to prune and once to finish the bound.

### What it found, per pivot, on square grids from 32 to 400 and cubes from 12 to 32

| | |
|---|---|
| visits per pivot | 1.011 to 1.012 in 2D, 1.009 to 1.010 on cubes |
| build, member, scan1, scan2v, restore | 1.00 at every size, both families |
| hash chain, pair, exact, stamp | identical, after the key was aligned |
| degree-bucket search | ours 6 to 100 times CHEAPER |
| `clear_flag` | fired zero times on every case |
| clique arena against the vendored workspace | 0.77 to 0.83 in 2D, 1.01 to 1.13 on cubes |

**So the growth was never work.** Same visits, same passes, fewer bucket-search steps, less storage
in the family where we were slowest.

### What that ruled out, and what it left

It ruled out the clique arena, which had been the standing suspicion: its excess cache misses FALL
with n in both families. It ruled out the hash, which had been the second suspicion. It ruled out
`clear_flag`, which nobody had checked. And a simulated L1 fed from the address streams put the
excess at a CONSTANT 0.09 misses per visit, real but far too small to explain a 1.8x ratio at sizes
where the whole working set fits in L2.

What it left was cost per visit, and the answer turned out to be two random probes per clique
visit, `mCliquePtr[c]` and `mCliqueSize[c]`, which `AMD_2` does not make at all: an element takes
over the `Pe` and `Len` of the variable that formed it. See `docs/DESIGN_DECISIONS.md` (2026-08-16).

### The lesson about the instrument, not the result

**The first version of this counter was blind in the place that mattered.** It counted the passes
already under suspicion and did not count the pass that builds `C[pivot]` at all, nor `absorb`, nor
`clear_flag`, and it counted the hash outer loop on the two sides with DIFFERENT denominators so
that a 0.95 ratio looked like agreement by coincidence. Each omission was found by asking what a
flat table could still be hiding rather than by reading it. A differential is only as good as its
least-considered counter, and the counters worth adding are the ones for passes nobody suspects.

## The special encodings, and the list-status dimension, 2026-08-23

`Buckets` is one class shared by all five drivers, and it holds three arrays of `n` `std::int32_t`
where a naive implementation holds five or six. What was folded away is a `degrees` array, an
`outmatched` flag on the mmd side and a `filed` flag, and what pays for the fold is a set of
reserved values in `mPrev`. This section is the whole encoding in one place.

### The three arrays, and what each code calls them

All three codes keep the same structure, a doubly linked list per degree pushed and read at the
head, in the same three arrays. Only the names differ, and two of them are misleading enough to be
worth a table before anything else in this section is read.

```
                        ours              AMD_2         genmmd
                        Buckets           AmdVendored   MmdVendored

  head per degree       mHead[d]          Head[deg]     head[dg]
  successor             mNext[u]          Next[i]       fwd[nd]
  predecessor           mPrev[u]          Last[i]       bwd[nd]

  indexed by            degree, vertex,   same          same
                        vertex
  size                  n, n, n           n, n, n       n+1 each, 1-based
  filing site           file()            line 1401     line 98
                                          and 2101
  unfiling site         unfile()          line 1500     lines 81 to 84
                                          and 1644
```

**`Last` IS THE PREDECESSOR, NOT THE TAIL**, which is the name most likely to send a reader looking
for something that is not there. `AMD_2`'s filing at line 1401 is our `file` with the names changed,
and its unfiling at 1500 is our `unfile`. `fwd` and `bwd` are forward and backward, the C port's
spelling of Liu's `DFORW` and `DBAKW`.

**AND BOTH VENDORED CODES PARK A LINK IN AN OUTPUT ARRAY, which is where the odd names come from.**
Neither allocates the link arrays at all. `AMD_2` declares `Int Last[]` as `/* the output
permutation */` and uses it as the predecessor link for the whole run, writing the permutation over
it at the end, `Last [k] = i` at line 2428. genmmd does the same thing twice over: `genmmd` calls
`mmdint(neqns, xadj, head, invp, perm, ...)` against a signature reading
`mmdint(..., int head[], int fwd[], int bwd[], ...)`, so **`fwd` is the caller's `invp` and `bwd` is
the caller's `perm`**, and `mmdnum(neqns, perm, invp, qsize)` at line 49 reads them back as the
permutation once the ordering is done. So each of those two arrays is named for its SECOND job and
used for its first.

We do not do this and the cost is visible: `Buckets` owns 3n of its own, the drivers carry a
`pivots` list, and the permutation is built from it by `orderAsMerged` or `orderAscending`. The
vendored codes get the links inside storage the caller supplied anyway. It is an aliasing trade
rather than an encoding one, so it sits outside the rest of this section, and it is recorded because
"why is the predecessor called `Last`" has no answer inside the degree lists and an obvious one
outside them.

### Reading genmmd's names

The C port shortened Liu's Fortran identifiers to two letters and kept nothing that says so. The
expansions below are inferred from the Fortran rather than documented in the file, and they are
worth having because the two-letter forms are unreadable on first contact:

```
nd   NODE       the vertex being visited, and the loop variable over 1..neqns
nb   NABOR      a neighbor of it
en   ENODE      the eliminated node whose clique is being walked
rn   RNODE      the node reached through it
md   MDNODE     the minimum degree node, the pivot
mn   MDNODE     the same, in genmmd's own pivot loop
pv   PVNODE     the previous node in a degree list, our mPrev's value
nx   NXNODE     the next node in one
fn   FNODE      the first node, the old head a filing pushes in front of
dg   DEG        a degree
nq   NQNBRS     the number of qualifying neighbors, and the value fwd carries
el   ELMNT      an element, which is our clique
```

**THE DIMENSION THAT ORGANISES IT IS LIST STATUS**, not the branch. In a list the two branches are
identical and the encoding is closed. Out of every list they diverge, and everything either branch
adds is added there.

```
                                  mmd                          amd
IN A LIST
  mNext[u]   successor, or NIL                 same
  mPrev[u]   [0, n-1]         predecessor      same
             [-(n+1), -2]     head             same
                              d = -mPrev[u] - 2

OUT OF EVERY LIST
  mNext[u]   NIL, left by unfile               NIL, or the HASH CHAIN link
  mPrev[u]   -1               UNFILED          -1  UNFILED
             INT32_MAX        OUTMATCHED       a raw HASH KEY

EITHER WAY
  mHead[d]   a vertex id, or NIL               same
```

### The head form carries the bucket, and that is what deletes the degree array

A head has no predecessor, so its `mPrev` slot is otherwise wasted. It spends it on the one fact the
head alone knows and nothing else records: which bucket the list belongs to. `unfile` needs `d` to
write `mHead[d]`, and takes no degree argument because it reads it back out of the slot.

```
d           0             1        ...   n-2          n-1
mPrev[u]   -2            -3        ...  -n           -(n+1)
```

so the head range is `[-(n+1), -2]`, contiguous, one value per degree, `n` values for `n` degrees.

**THE SHIFT IS BY TWO AND NOT BY ONE, AND DEGREE 0 IS WHY.** An isolated vertex is real and files at
0, so a shift of one would put its head at `-1`, which is `UNFILED`. Shifting by two keeps the head
range clear of it. This is not hypothetical: a scheme with no shift at all was tried, it builds and
passes every 2D case, and it corrupts a bucket list on a 3D grid at 6 a side.

**AND THE BOTTOM IS EXACT.** At the largest admissible `n`, which is `MAX_IDX = 2147483647`, degree
`n - 1 = 2147483646` encodes as `-2147483648`, which is `INT32_MIN` precisely. One more degree would
leave the type. On the other side a predecessor tops at `2147483646`, leaving `INT32_MAX` free for
`OUTMATCHED`. Both ends are reached and neither is exceeded.

That is only true because every branch now files at the TRUE DEGREE. While the mmd refresh filed at
the degree plus one the bucket index reached `n`, the head range was one value too long, and no
arrangement of the two sentinels covered it: putting `OUTMATCHED` at the bottom made it collide
there instead of overflow. See `private/MmdCorrected.cpp` and the 2026-08-23 entry in
`docs/DESIGN_DECISIONS.md`.

### The two sentinels are equality-only, and that is what let them move

Nothing in `Buckets` orders `mPrev` against `UNFILED` or `OUTMATCHED`; every use is `==` or `!=`. So
they sit at the two values a predecessor and a head cannot take rather than at either end of a
scale, and `UNFILED` at `-1` is what removes the bias that used to sit on every predecessor. genmmd
biases by one instead, storing `u + 1`, because its ids are 1-based and zero is therefore free to
mean unfiled; ours are 0-based, so copying that would cost a `+ 1` on every predecessor store and a
`- 1` on every read.

**CONTRAST `GONE`, WHICH IS ORDERED AND CANNOT MOVE.** It lives in `mMark`, not here, and the hot
line is `mMark[v] < mTag`: one load answering two questions, not seen this step and not dead. That
works only if every tag sorts above the initial value and below the dead value, so `NIL` is the
floor, the tags are the middle and `GONE` is the ceiling. A negative `GONE` would sort with `NIL`
and a dead vertex would read as live-and-unseen. `GONE` and `OUTMATCHED` are both `INT32_MAX` today,
in different arrays, and only one of them had a choice about it.

### Each branch adds exactly one thing, out of list, and they are mutually exclusive

**mmd adds `OUTMATCHED`**, genmmd's `bwd[nd] = -maxint`: a vertex withheld from the buckets while
still live and reachable, so it cannot be the minimum before the vertex that outmatched it. It
returns to `UNFILED` when an elimination reaches it, which `restore` does. No amd driver calls
`outmatch`.

**amd adds the hash overlay**, `AMD_2`'s supervariable detection reusing `Last` and `Next`:
`setKey`/`key` on `mPrev` and `setChain`/`chain` on `mNext`. No mmd driver touches them.

**ONLY AMD OVERLOADS `mNext`.** On the mmd side it is a pure link with one meaning; on the amd side
it is a link while filed and a chain link while not.

**AND THE AMD OVERLAY IS NOT SENTINEL-ENCODED, IT IS STATE-SCOPED**, which is the sharpest
difference in this section. A hash key is an arbitrary `int32` with no reserved range, so it can
look like a predecessor, like a head, or like `OUTMATCHED`. What keeps it safe is not its value but
when it is written: only after every member of `C[pivot]` has been unfiled, and read before any of
them is filed again. mmd's meanings are told apart by the value itself and need no such rule.

That asymmetry explains a detail in `unfile`, which returns early on both `UNFILED` and `OUTMATCHED`
and has no guard against a hash key. Those two are out-of-list values it can recognise; a key is one
it cannot. Calling `unfile` on a vertex holding a key would read the key as a predecessor and splice
a list that does not exist.

### The same slot in `AMD_2`, and where the two encodings part

`AMD_2`'s `Last[]` is our `mPrev`: one `Int` per supervariable, carrying the backward link of the
degree lists. Reading the two against each other is the clearest way to see what our encoding buys
and what it costs.

```
                              AMD_2  Last[i]              ours  mPrev[u]

IN A LIST
  predecessor                 a supervariable id          [0, n-1], a vertex id
  at the head of its list     EMPTY (-1)                  [-(n+1), -2], and it
                                                          encodes the bucket:
                                                          d = -mPrev[u] - 2

OUT OF EVERY LIST
  unfiled                     not represented             -1  UNFILED
  withheld                    no such state               INT32_MAX  OUTMATCHED
  holding a hash key          an arbitrary hval           an arbitrary key
```

Two meanings against four, in one slot of the same width.

**A HEAD IS WHERE THEY PART FIRST.** `AMD_2` writes `EMPTY` and finds the bucket somewhere else,
`Head [Degree [i]] = inext` at line 1510 of `private/AmdVendored.cpp`: the degree comes from
`Degree[i]`, which is live. We write `-(d + 2)` and read the bucket back out of the slot itself,
which is why `unfile` takes no degree argument and why the mmd side carries no `degrees` array.

**THE FOLD IS NOT AVAILABLE TO `AMD_2` AND WOULD NOT PAY IF IT WERE.** It keeps `Degree` regardless,
because its degree BOUND reads the value rather than merely filing by it, so removing the array is
not on the table. Our mmd side had no other reader once `unfile` stopped needing one, which is what
made the fold worth doing there and pointless here. The same fold reaches our amd drivers only
because they share the class.

**OUT OF A LIST, `AMD_2` REPRESENTS NOTHING.** Whether `i` is in a degree list is not answered from
`Last` at all; the walks know from context and from `Nv`'s sign. There is no unfiled value and no
withheld value, the second because `AMD_2` has no outmatching at all: aggressive absorption does
that work by another route entirely.

### The hash overlay, which is the one place they agree, and where `AMD_2` goes further

The idea is `AMD_2`'s and ours is copied from it. During supervariable detection a vertex is out of
every degree list, so both slots are free, and both codes park the hash there: the KEY in the
`Last`/`mPrev` slot and the CHAIN LINK in the `Next`/`mNext` slot. `Last [i] = hval` is line 1936;
ours is `setKey(u, hash)` beside `setChain(u, hashHead[hash])`.

**WHERE THE BUCKET HEADS GO IS THE DIFFERENCE.** `AMD_2` refuses to spend an array on them and
overlays `Head` itself, which is already the degree-list head array, telling the two apart by sign:

```
Head [hval] == EMPTY        no hash bucket and no degree list there
Head [hval] <  EMPTY        FLIP(i), a hash head, the degree list being empty
Head [hval] >= 0            a degree-list head j, so the hash head is parked in Last [j]
                            instead, that slot being free while j heads a list
```

Lines 1916 to 1929 write it and 1973 to 1990 read it back, restoring `Last [j] = EMPTY` when the
bucket is emptied. So in the third case `Last` carries a fourth meaning of its own, a borrowed hash
head, and it is legal only because a degree-list head's `Last` is `EMPTY` and therefore spare.

**THAT IS EXACTLY THE SLOT OUR HEAD FORM OCCUPIES**, so the trick does not compose with our
encoding: a head's `mPrev` holds `-(d + 2)` and is not free. `AmdFlat` therefore keeps a separate
`hashHead` array, which is `AMD_2`'s own alternative, present in that file as a commented-out
`Hhead` variant. Two coherent encodings that cannot be combined, and `n` `int32` is the cheaper side
of the trade; `src/AmdFlat.cpp` states this at the allocation.

**AND THE OVERLAY IS STATE-SCOPED IN BOTH CODES, NOT SENTINEL-ENCODED.** A key is an arbitrary value
with no reserved range: it can look like a predecessor, like a head, or like `OUTMATCHED`. Neither
code makes it distinguishable and neither needs to, because it exists only between the unfiling of
`C[pivot]`'s members and their refiling. `AMD_2` clears the borrowed slot back to `EMPTY` when done,
which is the same discipline written as an operation.

**ONE FALSE FRIEND.** `FLIP` appears in both codes and in three unrelated places: `AMD_2` uses it
for the hash head above and for its own pool compaction, and our `QuotientGraphCompacted::FLIPPED`
uses it for pool compaction only. Same `-x - 2` arithmetic, three separate jobs, no shared
invariant. `docs/DESIGN_DECISIONS.md` calls it an anti-model for that reason.

### Why the arrays are `n` and not `n + 1`

The heads are indexed by a DEGREE and the links by a VERTEX, two ranges of the same extent now that
filing is at the true degree: a degree reaches `n - 1`.

The heads were `n + 1` until 2026-08-23 and the extra slot answered two different questions at two
different times. It was needed while the mmd refresh filed at the degree plus one, the index then
reaching `n`. After that stopped it was still needed, for a reason nothing about filing explains:
the mmd PREPASS read `head(1)` unconditionally, before knowing whether any vertex has degree 1, so
at `n == 1` it read a bucket that did not exist. Shrinking the array without bounding that loop
builds, passes every other case, and dies in ASan on the one-vertex tridiagonal in `test_order`. The
loop is bounded now and the array is `n`.

**MEASURED AT BOTH STAGES RATHER THAN ARGUED.** Probing all three of `file`, `head` and `empty`
across 57 grids and 216 random patterns, five drivers each, the largest index ever passed to `mHead`
was `n - 1`, and every exception was `n == 1` in the mmd prepass. The amd drivers never appeared
once, which is what identified the loop as the cause rather than the array.

### The interface each branch actually uses, 2026-08-24

The sections above describe what the slots MEAN. This one is what the drivers CALL, which is the
same division seen from the other side and is the shorter statement of it.

```
                        mmd       amd      slot        what it is

  MUTATORS
    file(d, u)           x         x       both        into bucket d, at the head
    unfile(u)            x         x       both        out, bucket read from mPrev
    outmatch(u)          x         .       mPrev       withheld, still live
    restore(u)           x         .       mPrev       withheld -> unfiled
    setKey(u, k)         .         x       mPrev       hash key, while out of list
    setChain(u, v)       .         x       mNext       hash chain, while out of list

  QUERIES
    head(d)              x         x       mHead       first vertex of bucket d
    empty(d)             x         x       mHead       is bucket d empty
    next(u)              x         .       mNext       successor, for the prepass walk
    filed(u)             x         .       mPrev       is u in any list
    outmatched(u)        x         .       mPrev       is u withheld
    key(u)               .         x       mPrev       the stored hash key
    chain(u)             .         x       mNext       the next vertex in the hash bucket
```

**UNIFORM WITHIN EACH BRANCH.** All three mmd drivers call the same nine entry points and both amd
drivers call the same eight, so no driver differs from its branch sibling. That is what the clique
store not mattering here should look like, and it is worth checking after a change rather than
assumed: `MmdFlat`, `MmdChained` and `MmdCompacted` reach this class identically, and so do
`AmdFlat` and `AmdCompacted`.

**FOUR ARE SHARED AND THEY ARE THE LIST ITSELF**, `file`, `unfile`, `head` and `empty`. Whatever
else a minimum degree ordering does, it puts a vertex in a bucket, takes it out, and asks the lowest
bucket for a candidate.

**AND THE REST DIVIDES ON ONE LINE, THE SAME LINE THE ENCODING DIVIDES ON.** Every branch-specific
entry point reads or writes `mPrev` or `mNext` while the vertex is OUT of every list. mmd spends
that window on withholding and amd spends it on the hash. Neither branch can want both, since a
withheld vertex is not in a hash bucket and a hashed one is not withheld, which is the whole reason
one slot can serve two branches with no arbitration and no flag.

`next` is the one entry point that does not fit that reading. It is a plain successor read on a
FILED vertex, and only mmd makes it, in the prepass, taking the successor before the unfile
invalidates the link it is standing on. So the honest statement is that the shared part is the list,
the divided part is the out-of-list window, and one query sits outside both because a walk of a
whole bucket is an mmd-only shape.

**HOW THE TABLE WAS BUILT, because the obvious way is wrong.** Extracting `buckets.<method>` from
the five sources reports `refile` as an `AmdFlat` call, and `refile` was deleted on 2026-08-24 and
has no caller anywhere. It is a comment at `src/AmdFlat.cpp` naming the method to explain why the
driver does NOT use it. Comments have to be stripped before the match, which is the identical trap
`docs/QUOTIENT_GRAPH_USAGE.md` records against itself: it published call-difference counts of three
and nine that were really zero and two, from prose naming a method.

### The surface each class actually has, 2026-08-24

The table above is one class against two branches. This one is one branch's machinery against two
CLASSES: every member `QuotientGraph` and `QuotientGraphCompacted` declare, side by side. It exists
to be re-read after a change, since the whole claim of the pair is that they differ in STORAGE and
nothing else, and that claim is checkable only by looking.

```
                                     flat   compacted

  size / setup
    size                              x         x
    enableMarks                       x         x

  layout accessors
    adjacencyAmd / adjacencyMmd       x         x
    incidenceAmd / incidenceMmd       x         x
    adjacencySize                     x         x
    incidenceSize                     x         x

  cliques
    clique                            x         x
    cliqueSize                        x         x
    cliqueWeight                      x         x
    trimClique                        x         x
    killClique                        x         x
    bearClique                        .         x
    captureAbsorbed                   .         x

  weights
    weight                            x         x
    setAside                          x         x

  marks
    mark / setMark                    x         x
    advanceTag                        x         x
    markGone                          x         x
    number                            x         x
    eliminatedAmd / eliminatedMmd     x         x

  elimination
    eliminateMmd / eliminateAmd       x         x
    beginEliminationMmd / ...Amd      x         x
    pruneMmd / pruneAmd               x         x
    finishElimination                 x         x
    massEliminate                     x         x

  merging and absorption
    merge                             x         x
    absorbAggressively                x         x

  reachable set
    formReachableSetMmd / ...Amd      x         x
    formReachableSetInPlaceMmd / ..   .         x
    reachableSetWeight                x         x

  output
    orderAscending                    x         x
    orderAsMerged                     x         x

  configuration
    setReverseIncidence               x         x
    setLateMassElimination            x         x

  counters
    numLiveCliqueMembers              x         x
    numPeakCliqueMembers              x         x
    cliqueCountBalances               x         x
    numBornCliqueMembers              x         .
    numCompactions                    .         x
    compact                           .         x
```

**Thirty-five shared, no flat-only, six compacted-only, as of 2026-08-24.** It was thirty-eight,
four and six when this table was first taken the same day; what closed the flat-only column is
below.

**EVERY REMAINING DIFFERENCE IS THE LAYOUT.** `bearClique`, `captureAbsorbed`,
`compact` and `numCompactions` are pool machinery: a store that fills has to be compacted, and a
walk that destroys `I[pivot]` while it runs has to copy the ids out first. `numBornCliqueMembers`
against `numCompactions` is the same fact from the two sides, one store paying for what it holds and
the other for what it reclaims. `formReachableSetInPlaceMmd` and `formReachableSetInPlaceAmd` are
the pool's special case: with `I[pivot]` empty the reach is a SUBSET of `A[pivot]` and is compacted
where it stands, so the pool is not touched at all. The flat class has nothing to spare and no use
for the case.

**AND THE FLAT-ONLY COLUMN CLOSED THE SAME DAY, by deletion.** `reachableSet` and `reachableSize`
were on the flat class and called from NOWHERE: not `src/`, not `tests/`, not `examples/`, not
`benchmarks/`. `reachableSize`'s last caller was `retired/Mmd1.cpp`, which left the build on
2026-08-21, so it had been dead for three days.

**THE CHAINED CLASS WAS WORSE AND NOTHING COULD HAVE CAUGHT IT.** It declared three of them,
`reachableSet` in both forms and `reachableSize`, and DEFINED NONE. An uncalled declaration never
has to link, so a header can carry a promise with no body indefinitely. Five members gone in all.

What they were is worth keeping, because it is not what the first reading said. They are not the
shape `Buckets::refile` had, a member that outlived its caller. They never had one: they are the
readable form of a computation the production path only ever runs fused, which is the next section.
The argument for deleting them is simply that nothing calls them.

### What the reachable set is, and why it is five members and not one

The prototypes have ONE function. `mmd3_neighbors(A, I, C, eliminated, mark, tag, u)` returns the
list, and the file uses it twice: `sum(len(super_members[v]) for v in neighbors)` in the degree, and
as the new clique in `mmd3_eliminate`. Production split that one function five ways, and only two
of the five materialize anything:

```
prototype                       production                   called by
mmd3_neighbors, in the degree   reachableSetWeight(u)           MmdFlat, MmdCompacted, MmdChained
mmd3_neighbors, in eliminate    formReachableSetMmd / ...Amd     beginEliminationMmd / ...Amd,
                                                             appending into the clique's own block
mmd3_neighbors, literally       reachableSet(u)              nobody, DELETED 2026-08-24
                                reachableSize(u)             nobody, DELETED 2026-08-24
```

**NEITHER LIVE FORM BUILDS THE SET.** `reachableSetWeight` sums as it walks and returns a number.
The appending pair writes straight into `mCliqueSrc` at the block the clique will occupy, because
the reach IS the clique. So the prototype's return value has no counterpart in either hot path, and
`reachableSet` is the readable form kept beside the fused ones. Its own comment says as much: a
convenience with no caller inside the elimination, which undoes the sign negation the appending form
deliberately leaves for `massEliminate`, so that a reader reaching for it gets a query rather than a
half-finished elimination.

The open question this raised was why the flat class had the query form and the compacted class did
not, which is not the layout either. It was settled by deleting rather than by matching: adding the
pair to the compacted class would have added a pass that exists nowhere in either branch's real
path, purely so two surfaces agree.

### And the branch split on the reachable set is the whole 21x

```
                        during elimination                for the degree
mmd    formReachableSetMmd, into the clique           reachableSetWeight, per vertex, per refresh
amd    formReachableSetAmd, into the clique           never
```

**BOTH BRANCHES FORM THE REACH ONCE PER PIVOT**, to build the clique. They are identical there.

**ONLY MMD FORMS IT AGAIN.** An exact degree IS `|reach(u)|` weighted, so there is no way to get it
except by walking u's whole neighborhood, once for every vertex whose degree the elimination
changed. An approximate degree decomposes the same quantity into

    |A[u] - C[p]|  +  |C[p] - {u}|  +  sum over c of |C[c] - C[p]|

every term of which the prune has already walked, so the amd branch never forms the set outside
`beginElimination` and calls NOTHING between one elimination and the next.

So the amd branch reads counts off a pass it was making anyway where the mmd branch makes a fresh
pass per reached vertex. That is why `reachableSetWeight` has exactly three callers and all three
are mmd drivers, and it is where the branch's factor of 21 on the benchmark set comes from: not a
better data structure, the same reach formed a different number of times.

`reachableSet` returned the neighborhood of a live vertex BY VALUE, which is what an external user
of `QuotientGraph` would reach for, and that was the one argument for keeping it. It lost to the
plainer one: an uncalled member is read by everyone who opens the file and used by nobody. If a
caller ever wants the query, it is five lines over `formReachableSetMmd` and the restore the prune
now does anyway.

**HOW THE TABLE WAS BUILT, and the first two attempts were both wrong.** A regex over the whole
file collects definitions and call sites along with declarations, and reports 62 members against
52. Restricting it to the class body still misses every member returning a pointer, `adjacencyAmd`
and its three siblings among them. The version above walks the class body and then checks the
known-awkward returns by hand. Anyone regenerating it should assume the naive extraction is wrong
and diff against this list rather than trust it.

## Indistinguishable vertices, and the three ways they are found, 2026-08-24

Two live vertices are INDISTINGUISHABLE when their reachable sets agree, each counting the other:

    reach(u) | {u}  ==  reach(v) | {v}

If that holds it holds forever. Every elimination from here treats them identically, they keep the
same degree, they are eliminated at the same moment, and they produce the same fill. So one is kept
as the PRINCIPAL and the other stops existing as a separate entity: its lists are emptied, its
weight goes to zero, and it is threaded onto the principal's chain through `mSuperNext` and
`mSuperLast`.

**`mWeight[u]` IS HOW MANY ORIGINAL VERTICES u STANDS FOR**, which is what makes it the thing a
degree sums:

    degree(u) = sum of weight(v) over v in reach(u)

A neighbour standing for three originals contributes three, because eliminating u will fill against
all three. Counting entries instead of weights under-counts the fill by exactly the merging. At
output, `orderAsMerged` walks the chain and emits the principal followed by its members, so a
supervariable of weight three occupies three consecutive slots of the permutation.

**THE PAYOFF IS COMPOUND, which is why both branches spend effort finding them.** A merge removes a
vertex from the graph, so every later reach is shorter, every later degree cheaper, and one pivot
selection eliminates a whole group at once.

### Three mechanisms, and they are not alternatives

```
mass elimination     reach(u) == C[pivot] - {u}      structural, free       both branches
q2h merge            two sources, same tag           free, partial          mmd only
hash detection       any two in C[pivot]             a hash per vertex      amd only
```

They are increasing cost for increasing coverage. Mass elimination is the case the elimination
hands you and neither branch would decline it. Beyond that, mmd takes what the tag has already
proved and stops, catching a subset; amd pays a hash and catches the general case. **Neither branch
does both.**

**MASS ELIMINATION IS A MERGE, and the weight arithmetic is the same one.** `massEliminate` appends
u's chain to the pivot's, adds the weights and zeroes u's, which is `merge`'s body written out:

```
massEliminate    mSuperNext[mSuperLast[pivot]] = u;  mSuperLast[pivot] = mSuperLast[u];
                 mWeight[pivot] += mWeight[u];       mWeight[u] = 0;

merge(u, v)      mSuperNext[mSuperLast[u]] = v;      mSuperLast[u] = mSuperLast[v];
                 mWeight[u] += mWeight[v];           mWeight[v] = 0;
```

What differs is only WHERE the merged vertex lands: into the pivot, which is on its way out, rather
than into a live vertex that stays. So the question "does the weight apply during mass elimination"
answers yes, and it has to: without it the pivot would occupy fewer slots in the permutation than
the originals it stands for.

**AND `massEliminate` DOES NOT CALL `merge`, it repeats it.** The two other things it does around
those four lines, `markGone(u)` and zeroing the incidence size, are also `merge`'s. So the body
could be `merge(pivot, u)` plus the one line `mCliqueWeight -= weight(u)`, with the first loop
reduced to detection alone; the detection test reads only u's own segment, so deferring the zeroing
to the second loop cannot change what the first loop decides. Not done, and it wants the digest
rather than an argument.

### What genmmd's q2h branch actually does, and why only two sources

The two-source queue is not merely a cheaper degree update. Inside its walk, when it meets a vertex
`nd` already stamped with this step's tag:

```
if (fwd[nd] == 2)  qsize[en] += qsize[nd]; qsize[nd] = 0; fwd[nd] = -en;   merge
else               bwd[nd] = -maxint;                                      outmatch
```

**A two-source vertex sharing the tag with another two-source vertex is indistinguishable and is
merged.** One that shares the tag but has more sources is only OUTMATCHED: its reach is contained
in the other's, so it can never be the unique minimum and is withheld from selection, but it is not
absorbed and it is not gone.

**Why the two-source restriction is exact rather than cautious.** With exactly two sources, the new
clique plus one other, sharing a tag is enough to conclude the reaches are EQUAL, since there is
nothing else either vertex could reach through. With more sources the tag proves only CONTAINMENT,
which is why the `qxh` path does neither: it recomputes the degree and nothing else.

**Outmatching is not a merge**, and the difference is why it needs state of its own. A merge removes
a vertex from the graph, so zeroing its weight is enough to make every walk skip it. An outmatched
vertex is still live, still in every list that names it, still contributing to its neighbours'
degrees; only its candidacy is withheld. That is `OUTMATCHED` in `mPrev`, and it is mmd's alone.

## What `mWeight` and `mMark` encode between them, 2026-08-24

Every question either branch asks about a vertex is answered from these two arrays. Below is what
each combination means, over the flat class; the compacted and chained ones agree.

**One table per branch, because neither branch has all the states.** mmd never sets a row aside,
having no dense-row rule; amd never numbers, its prepass taking degree 0 alone. And `mMark` is
EMPTY on the amd branch, `enableMarks` being called by the three mmd drivers alone, so it has no
column there at all.

```
MMD                         mWeight[u]        mMark[u]        eliminatedMmd
                                                              (mMark[u] == GONE)

live, principal             > 0, its size     NIL or old tag  false
live, in the reach          < 0, negated      unchanged       false
numbered by the prepass     1, unchanged      GONE            true
merged into another         0                 GONE            true
pivot, eliminated           > 0, unchanged    GONE            true
```

```
AMD                         mWeight[u]                        eliminatedAmd
                                                              (mWeight[u] == 0)

live, principal             > 0, its size                     false
live, in the reach          < 0, negated                      false
set aside, dense row        0                                 true
merged into another         0                                 true
pivot, eliminated           > 0, unchanged                    false, AND IT IS DEAD
```

**ON MMD THE TEST IS EXACT.** All three ways out of the graph write `GONE`, so one comparison
answers for a numbered vertex, a merged one and a pivot alike, with no exceptions to remember.

**ON AMD IT IS EXACT ONLY FOR THE DEATHS AMD CAN PRODUCE.** The last row is the one to hold onto:
a pivot is eliminated and its weight is still POSITIVE, so `eliminatedAmd` returns false about a
dead vertex. The accessor says as much, and it is safe rather than lucky: the amd driver unfiles a
pivot when it chooses it and never revisits it. So the predicate is not "is u eliminated" but "is u
eliminated, given you would never ask about a pivot".

That makes the amd test cheaper AND narrower at once. It is not a weaker way of asking the same
question; it is a different question whose answer happens to agree everywhere amd asks it.

**AND THE REASON THE WEIGHT CANNOT SIMPLY BE ZEROED IN THOSE TWO CASES IS THE PERMUTATION.**
`orderAscending` reads `mWeight[pivot]` to reserve a supervariable's room, and `orderAsMerged` walks
the chain the weight counts. A numbered vertex and a pivot both still need a place in the output, so
their weight has to keep saying how many originals they stand for. That is the whole of it: the
weight is two facts, how many originals and whether alive, and where the two disagree the first
wins.

**`markGone` IS CALLED ON EVERY DEATH ON BOTH BRANCHES**, from `setAside`, `finishElimination` for
the pivot, `massEliminate` and `merge`, with `number` writing `GONE` directly. On amd every one of
those calls is a no-op against an empty vector. So **the mark is a complete death record and the
weight is not**, and mmd is simply the branch that kept the record.

### Where each array is asked, and why the tests differ

**THE NEGATION IS A PHASE AND NOT A STATE**, which is why it appears in both tables as a row and in
neither as a death. A negative weight is a live vertex currently in the reach being built. It exists
between the walk and the end of the prune, which restores it, and nothing outside that window ever
sees one. That is what lets `vWeight > 0` answer two questions in one load: v is live, AND v is not
already in this reach.

```
formReachableSetMmd   adjacency   vWeight > 0 && !(mHasNumbered && mMark[v] == GONE)
                      cliques     vWeight > 0
formReachableSetAmd   adjacency   vWeight > 0
                      cliques     vWeight > 0
reachableSetWeight    adjacency   mMark[v] != GONE
                      cliques     mMark[v] < mTag
```

**ALL FOUR ASK THE SAME TWO QUESTIONS. What differs is the currency each has to pay with.**

```
                        liveness              distinctness
form*, adjacency        vWeight > 0           vWeight > 0     one load, both answers
        plus GONE, because a numbered vertex is live AND positive
form*, cliques          vWeight > 0           vWeight > 0     one load, both answers

reachableSetWeight
        adjacency       mMark[v] != GONE      not needed      A[u] has no duplicates,
                                                              and it runs first
        cliques         mMark[v] < mTag       mMark[v] < mTag one load, both answers
```

**THE NEGATION IS THE DISTINCTNESS MARK IN THE `form*` WALKS.** A vertex already in this reach is
negative, so `vWeight > 0` rejects it, and the SAME comparison rejects merged and set-aside
vertices at zero. Two questions, one load, no tag and no stamp anywhere.

**`reachableSetWeight` CANNOT SPEND THAT CURRENCY, because it is a query and must leave the graph as
it found it.** It emits nothing, so there is no reach to negate and nothing to restore against. It
falls back to the tag, and the tag then absorbs the liveness question too, `GONE` sorting above it.

**AND IT COULD NOT CHOOSE OTHERWISE.** The `form*` walks get their restore for free: the prune walks
C[pivot] at the end anyway and that is where the signs come back. `reachableSetWeight` runs per
REFRESHED VERTEX and nothing walks its result afterwards, so negating would buy it a second pass
over the reach purely to undo itself. The cheapest negating form would need a member scratch holding
the emitted vertices, cleared per call so its capacity survives, walked afterwards to restore, which
is more than the mark costs.

So the tests are not four conventions. They are one pair of questions answered from whichever array
is already being touched, and the single residual asymmetry is `GONE` in the mmd adjacency test,
covering the one state where the weight says live and the vertex is not.

**GONE IS ASKED IN ADJACENCY LISTS AND NEVER IN CLIQUES, and the asymmetry is exact.** A numbered
vertex LINGERS in `A[v]`, since `number()` marks and returns and nobody rewrites its neighbours'
lists. It can NEVER appear in a `C[c]`: the prepass completes before the first elimination, so no
clique exists yet, and every clique since is built from a reach whose adjacency loop already
skipped it. A guard in the clique loops would be dead code.

**`reachableSetWeight` IS THE ODD ONE AND HAS A REASON.** It leads with the mark rather than the
weight because it must touch `mMark[v]` anyway, writing `mMark[v] = mTag` on the next line for
DISTINCTNESS: an exact degree counts each reached vertex once, and u can reach the same v through
two cliques. The load is happening regardless, so `!= GONE` rides on it free, and guarding it with
`mHasNumbered` would add a test to save nothing.

Its clique loop's `mMark[v] < mTag` then answers THREE questions from one load, `GONE` sorting above
every tag: distinct this call, not numbered, not merged. That fold is why `GONE` cannot simply be
retired in favour of a zero weight. Merged vertices linger in other cliques' member lists, so
dropping `GONE` would cost that loop a second test, in the hottest loop in the ordering.

**AND AMD NEEDS NO DISTINCTNESS STAMP AT ALL**, which is the real reason it has no mark array. Its
bound is `|A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]|`, summed term by term over passes the
prune already made. A vertex in two cliques is counted twice and that is not a defect: it is what
makes the result a bound rather than a degree.

## What each file is, and what it adds

The layers below md5 are a line, each rung adding one idea, and the section headings above describe
them in order. From md5 the ladder forks, and from there the names carry two suffix conventions at
once, which is easy to lose track of. This table is the inventory.

```
name     where         base or extras   what it adds
------------------------------------------------------------------------------------------
mmd1     experiment    base             N/A
mmd2     experiment    extras           the q2h path, outmatching
mmd3     experiment    aligned          nothing; genmmd's list order
amd1     experiment    base             N/A
amd2     experiment    extras           aggressive absorption, hash supervariable detection
amd3     experiment    aligned          nothing; AMD_2's list order
amd4     experiment    extras           dense rows by the alpha ratio, amd4_aat and
                                        amd4_preprocess forming A + A', the postorder,
                                        amd4_valid and the Control/Info interface
                                        PARKED and temporary; forks from amd2, not from amd3

Mmd1     production    base             N/A
Mmd2     production    extras           the q2h path, outmatching
Mmd3     production    aligned          nothing; genmmd's list order. The default ordering
Amd1     production    base             N/A
Amd2     production    extras           aggressive absorption, hash supervariable detection
Amd1B    production    base             N/A
Amd2B    production    extras           aggressive absorption, hash supervariable detection
```

**"Base" is base for its BRANCH, not bare.** Every file in the table sits on md5, so all of them
already carry supervariables, mass elimination, maintained degrees and degree buckets. A base file
adds its branch's one idea and nothing else: multiple elimination for the MMD branch, the
approximate degree for the AMD branch.

**A THIRD suffix meaning, and it is the digit 3 on both branches.** `mmd3` and `amd3` add no
mechanism at all: each is the layer below it with the vendored routine's list order, so that the
permutations can be compared as an equality test rather than judged as a fill number. They are
different orderings from their parents, so the digit is honest, but what they add is a set of
tie-break conventions rather than anything the ordering does. `amd4` is the exception that proves
the rule and says so in its own header: it forks sideways from `amd2`, does not contain `amd3`,
and is temporary.

**A trailing digit and a trailing B are different axes.** A digit means a DIFFERENT ORDERING: mmd2
has mechanisms mmd1 lacks, so their permutations and their fill legitimately differ and both are
correct. A B means the SAME ORDERING computed on a different schedule, so Amd1B must return exactly
Amd1's permutation and a difference is a defect in one of them. That is why the B rows repeat their
counterpart's extras rather than adding anything: they carry the same mechanisms and change only
when the work is done, folding the driver's first scan into the eliminator's walk. They have no
prototype, being a re-schedule rather than a layer, and their oracle is the identity check in
`tests/test_order.cpp` rather than a Python twin.

**amd4 has no production counterpart and is not expected to get one.** Its four items are
not ordering ideas: three are input conditioning and output ordering that Oblio does elsewhere or
does not need, and the fourth is a control interface Oblio has no equivalent of. So it is checked
against its own C++ twin and against nothing else, and it is not expected ever to be ported.

That split is what lets **amd2 be checked by PERMUTATION** rather than by fill. Before it, the
experiment's amd2 ended with a postorder that production does not do, so the two could only be
compared on `nnz(L)`. With the postorder moved to amd4, every ported layer is on the strong oracle
and `PORTED_FILL` in the Makefile is empty.

**One thing amd2 takes from amd4 rather than from amd1**, and it is not optional. A hash merge
leaves the merged vertex in place with weight zero rather than removing it from every list, so every
walk has to skip eliminated vertices. amd1 has no such vertices and its core functions take no
`eliminated` argument; amd2's do. That is the prototype's version of what production calls live
merges, and it is why amd2 is amd1 plus two mechanisms AND a liveness-aware core.

**What the extras are.** For MMD: the q2h path finds vertices with exactly two sources and merges
indistinguishable pairs as a by-product of the exact-degree union; outmatching withholds a vertex
that reaches strictly more than another until an elimination puts it back in the running. For AMD:
aggressive absorption kills a clique whose members all lie inside the new one, free because
`outside[c] == 0` is already computed for the bound; hash supervariable detection finds pairs
indistinguishable from each other rather than from the pivot. The section above walks the AMD pair
in detail and contrasts it with the MMD route.

## One example end to end: the 3 by 3 grid, 2026-08-26

The sections above describe the supervariable encodings a piece at a time. This one runs a single
graph all the way through, from the prototypes' containers to the production arrays to the
permutation, so the pieces can be read against each other rather than in turn.

**One example, not one per branch.** The elimination below is `MmdFlat`'s on this grid. Both output
functions are then applied to that same state, including the one the mmd drivers do not call, so
that `orderAsMerged` and `orderAscending` differ in the emission and in nothing else. Running the
amd drivers on this grid would give a different pivot sequence and a different chain, since the
branches select differently, and comparing two states would put a real difference in selection next
to the one point at issue here. The hypothetical is only that a state mmd produced is handed to the
function amd calls. Nothing in either function knows or cares which branch built the arrays.

The graph is a five-point square grid at 3 a side, which `gridGraph` builds and every layer here
already runs:

```
0 1 2
3 4 5
6 7 8
```

The elimination takes pivots 8, 6, 2, 0, 7, 1 and then 5, at which point vertices 3, 4 and 5 have
become indistinguishable and 5 absorbs 4 and then 3.

### The prototypes

`super_members` is a list of lists indexed by root, and `pivots` is the order over supervariables.
`python3 mmd3.py grid 3` prints both at every iteration, and after the last elimination:

```
members: [[0] [1] [2] [] [] [5 4 3] [6] [7] [8]]
pivots:  [8, 6, 2, 0, 7, 1, 5]
```

One list is non-trivial, `[5 4 3]`, holding its root first and then the vertices absorbed into it,
in the order they were absorbed. It is not ascending, and that is the whole reason the two output
functions can be told apart on this example. The absorbed vertices are left as `[]`.

### Production

The same fact is three arrays plus the driver's pivot vector:

```
                       v =   0   1   2   3   4   5   6   7   8

mSuperNext:                 -1  -1  -1  -1   3   4  -1  -1  -1
mSuperLast:                  0   1   2   3   4   3   6   7   8
mWeight:                     1   1   1   0   0   3   1   1   1

pivots, positional and not indexed by v:  [8, 6, 2, 0, 7, 1, 5]
```

The list `[5 4 3]` is four entries across two arrays:

```
mSuperNext[5] = 4          head                          tail
mSuperNext[4] = 3            5 -----> 4 -----> 3 -----> NIL
mSuperNext[3] = NIL          ^                 ^
mSuperLast[5] = 3            |                 |
                             the pivot         mSuperLast[5]
```

No element of the list is stored as an element. Each is stored as the successor of the one before
it, at the index of the one before it, so 4 lives at `mSuperNext[5]` and 3 lives at `mSuperNext[4]`.

**The head lives nowhere in the three arrays, and that is what `pivots` is for.** Nothing marks a
pivot and no member points back to one, so `mSuperNext` on its own is a set of links with no entry
points. Read in index order it gives `3` at entry 4 and `4` at entry 5, which is neither the list
nor its reverse. The pivot vector supplies the heads, and only then does the row resolve into
chains.

Note what `pivots` is NOT indexed by. It is dense and positional, one entry per elimination in
elimination order, where the three arrays are indexed by vertex. Its fifth entry being 7 says the
fifth pivot was vertex 7, and says nothing about vertex 5. Mixing the two readings is the easiest
mistake to make against a dump like this, which is why it is not aligned under the `v =` columns.

Term for term:

```
prototype                   production

super_members[u]            the chain from u through mSuperNext, NIL-terminated,
                            with mSuperLast[u] naming its tail
len(super_members[u])       mWeight[u]
super_members[u] == []      mWeight[u] == 0
pivots                      pivots, the driver's own vector, identical
```

Three things the arrays say that the lists do not have to.

- **The pivot is an index, not a value.** `super_members[5]` is a lookup by root in the prototype
  and the three arrays are the same lookup: 5 being the pivot means entry 5 is where the chain
  begins and where the tail and the weight live. There is no flag saying so and no back-pointer from
  a member, which is why `orderAscending` has to build one for the length of one call and `mmd3.py`
  builds `root_of` in the same shape. Both construct, temporarily, an inverse the data structure
  does not keep. The code says pivot throughout, `root` being the prototypes' word.
- **`mWeight` is `len()` made affordable.** A list knows its own length and a chain does not, and
  the degree computation asks for it once per reached vertex. `mmd3.py` says this where it declines
  to keep a weight array of its own: MMD keeps `qsize` because its members are a chain, not a list.
  It doubles as the liveness test, `mWeight[u] += mWeight[v]; mWeight[v] = 0;` making the size of
  the absorbing list and the emptiness of the absorbed one one pair of stores.
- **`mSuperLast` is the splice made affordable.** The prototype's `super_members[u] +=
  super_members[v]` attaches at a tail Python already holds. Production needs BOTH tails at once,
  `u`'s to attach to and `v`'s to inherit, which is exactly the two lines of `merge`:
  `mSuperNext[mSuperLast[u]] = v; mSuperLast[u] = mSuperLast[v];`. Two stores whatever the sizes,
  order preserved on both sides by construction. This is the one array with no counterpart in the
  prototypes at all, and a member's stale entry never has to be cleared because the value is copied
  out of it and the entry is then never read again. In the dump above `mSuperLast[4]` still reads 4
  and `mSuperLast[3]` still reads 3.

### The final order

Both functions do conceptually the same thing: walk `pivots` in order and splice in each pivot's
members after it. They differ in the member order, and in how each avoids an actual insertion.

**`orderAsMerged`, which the two amd drivers call.** One pass, appending. For each pivot, walk its
chain from the pivot to NIL and push. The chain is already pivot-then-members, so the members arrive
in merge order because that is the order the chain holds:

```
8 | 6 | 2 | 0 | 7 | 1 | 5 4 3
```

```
8 6 2 0 7 1 5 4 3
```

**`orderAscending`, which the three mmd drivers call.** Two passes over one array, `cursor`, which
carries a position at a pivot and an encoded pivot at a member. The first reserves: place each
pivot, advance by `mWeight[pivot]` to leave a gap of exactly the right size behind it, and stamp
each member of that chain with the pivot it belongs to, as `cursor[u] = -(pivot + 1)`.

```
pivot   placed at   room reserved   members stamped
8       0           1
6       1           1
2       2           1
0       3           1
7       4           1
1       5           1
5       6           3               cursor[4] = cursor[3] = -6
```

The second fills: sweep `u = 0 .. n - 1`, decode with `-(cursor[u] + 1)`, and drop each member into
its pivot's advancing cursor. The decode is the encode written again, `-(x + 1)` being its own
inverse, and its sign is the test: a pivot's entry is a position, so the decode comes out negative
and the sweep skips it. Vertex 3 is reached before vertex 4, so 3 takes position 7 and 4 takes
position 8:

```
8 | 6 | 2 | 0 | 7 | 1 | 5 3 4
```

```
8 6 2 0 7 1 5 3 4
```

The ascending order is not sorted for and is not chosen at the chain. It falls out of the sweep's
direction, which is `mmdnum` exactly: genmmd gets the same result from a `1 .. neqns` scan chasing
parent pointers. And the reservation is needed because the function must know where a pivot's
members will end up before it has seen any of them, which `mWeight` answers in one load instead of
a walk.

So the two differ in one transposition on this example, `5 4 3` against `5 3 4`, and in nothing
else. Same pivots, same chains, same fill of 5, same forest, since the members of a supervariable
are indistinguishable by construction. Only the permutation moves.

### Why the branches close differently

This is a question about the oracles, not about the branches. `genmmd` returns its own numbering and
`mmdnum` produces the ascending member order, so the mmd drivers reproduce it. The amd oracle is not
`amd_order`'s output but the raw order the hook reconstructs upstream of `AMD_postorder`, and the
hook emits each pivot's membership in the order it accumulated it, which is merge order. `AMD_2`'s
own output assembly, `Next[i] = Next[e]; Next[e]++` walked over `i` ascending, numbers non-principal
variables ascending and places the element last within its supervariable. So on this one point the
two vendored routines agree with each other and our amd drivers are the odd ones out, because the
amd oracle is taken at the finalize marker, upstream of an output assembly we do not reproduce, and
the hook had to invent a member order there.

### Reproducing it

```
python3 mmd3.py grid 3
./production_cpp mmd3 grid 3
```

`MmdCompacted` and `MmdChained` return what `MmdFlat` returns here, and `AmdCompacted` returns what
`AmdFlat` returns; the store is not what this section is about. To see the amd branch's own
elimination on this grid rather than the emission studied above, run `python3 amd3.py grid 3` and
`./production_cpp amd3 grid 3`, which take pivots 8, 6, 2, 0, 3, 4 and leave the single chain
`4 -> 1 -> 7 -> 5`, mass-eliminated rather than merged in pairs.

**The input convention is load-bearing and this example is small enough to get it wrong by hand.**
`gridGraph` emits each adjacency list in ascending order and `toCsc` writes the diagonal in its
sorted position rather than at the front. A hand-built driver that pushes the four neighbors in
compass order, or the diagonal first, gives a different graph in the only sense that matters here.
mmd is insensitive to it on this grid and amd is not: the wrong list order takes pivot 1 where the
right one takes pivot 3. Same fill, different permutation, and nothing in the output says which
convention produced it. See the note above `gridGraph` on ascending order, which is the 3D builder's
version of the same trap.


## Every state a vertex can be in, and the two questions they answer differently, 2026-08-26

A vertex is in one of six states, and reading the code as though there were one live state and one
dead one gets four of them wrong. They differ along two independent questions, "is it still in the
graph" and "does it have a position in the permutation":

```
                           in graph   has a position       weight     mark (mmd)   mark (amd)

live                       yes        not yet              1 .. n     < mTag       none
live, in the current reach yes        not yet              -(1 .. n)  < mTag       none
pivot, eliminated          no         yes, at that step    KEPT       GONE         none
merged                     no         no, rides later      0          GONE         none
prepass-numbered           no         yes, at that step    1          GONE         none
set aside                  no         no, appended last    0          n/a          none
```

**THE SECOND ROW IS A WINDOW, NOT A STATE THE GRAPH RESTS IN.** A weight is negated by the form walk
that puts the vertex into C[pivot] and restored by that elimination's prune, so it is live only
between those two points and no reader outside them ever meets it. The sign is the reach's
membership mark, which is why the prune's test is `mWeight[v] <= 0` and why the restore cannot move
ahead of the loop that consumes it.

**AND THE amd COLUMN IS EMPTY ALL THE WAY DOWN, WHICH IS NOT AN OVERSIGHT.** That branch never
calls `enableMarks`, so `mMark` is empty, and `markGone` is `if (!mMark.empty()) mMark[u] = GONE;`,
a no-op there. Every `markGone` call in the class is therefore mmd's alone at run time even though
the code is shared. What amd has instead is the zero weight, in every case but one: a pivot keeps
its weight, so `eliminatedAmd` must never be asked about a pivot, and the driver never does.

**A WEIGHT'S RANGE IS WHAT MAKES THAT WORK, AND IT DIFFERS FROM AN INDEX'S.** A weight counts
original vertices, so it runs 1 to n and n is capped at `MAX_IDX`, which is `INT32_MAX`. Zero is
outside the live range and is therefore free to mean dead. An index has no such spare value, which
is why the two encodings are not the same map even though they share an image:

```
index    live range  0 .. 2^31-2      encode -(x + 1)   ->  -(2^31-1) .. -1
weight   live range  1 .. 2^31-1      encode -x         ->  -(2^31-1) .. -1
```

The `+ 1` an index needs is exactly the offset a weight already has. `-2^31` is unreachable in both,
and for weights that is load-bearing rather than tidy: negating `INT32_MIN` is undefined, and it
cannot arise because a weight never reaches `2^31`.

**MERGED IS THE ONLY STATE WITH NO POSITION OF ITS OWN, AND IT NEVER GETS ONE.** It is threaded onto
its principal's chain through `mSuperNext`, and its slot is issued at output time when
`orderAsMerged` or `orderAscending` expands that chain. So it is waiting on the OUTPUT rather than
on an elimination, and the pivot it rides behind may be selected many steps after the merge that
folded it away. It is also the only exit that leaves a vertex named by lists it is no longer in:
nothing purges an absorbed vertex from the adjacency of everyone who knew it, and the zero weight is
what keeps it out of every later reachable set.

**THE TWO ODDITIES IN THE WEIGHT COLUMN ARE WHY THE BRANCHES NEED DIFFERENT PREDICATES.** A pivot
keeps its weight, that weight being the supervariable's size and `orderAscending` needing it to
reserve room. A prepass-numbered vertex keeps weight one deliberately so its neighbors' degrees
still count it. On the amd branch neither case is reachable through the predicate: a pivot is
unfiled when chosen and never revisited, and a degree-zero vertex is in nobody's adjacency, so no
walk can meet either. That is the whole reason `eliminatedAmd` can be `mWeight[u] == 0` while mmd
needs a mark array.

### Where the branches could converge, and it is the same axis as the detector

**SET ASIDE IS THE DENSE-ROW RULE AND IT IS amd's ALONE.** A row whose degree exceeds
`max(16, 10 * sqrt(n))` is not eliminated, not available, kept out of every reachable set by a zero
weight, and appended to the permutation at the end. Nothing about it needs an approximate degree.
mmd has no rule at all, and on a matrix with a hub of degree in the thousands that costs mmd exactly
what it cost amd before the rule went in: the hub sits in every reachable set it touches, which was
20.4 ms against 0.43 on `GHS_indef/bloweybq` and 41.0 against 1.24 on `bloweybl`. Adding it to mmd
is a variation, and one with a measured prize rather than a guess.

**AND THE TWO PREPASSES ARE NOT THE SAME PASS.** Read from the code rather than the comments:

```
mmd    buckets 0 AND 1, drained together, and it calls qg.number on each
amd    degree 0 only, pushed straight onto pivots, no number call, plus the dense-row fork
```

genmmd's `if (dg == 0) dg = 1` lumps degree zero and degree one into one bucket, which is why mmd
takes both; `AMD_2`'s initialization pass tests `deg == 0` alone. The difference is not cosmetic. A
degree-ZERO vertex is in nobody's adjacency, so no walk can reach it and nothing has to be marked;
a degree-ONE vertex is still named by its neighbor, so numbering it requires a mark and that is what
`number` writes. Everything mmd pays for `mMark` on this path follows from taking bucket 1.

**WHICH IS THE TOP OF THE CHAIN RECORDED IN `docs/NEXT.md`**, seen from the other end. That chain
runs `orderAscending` reads `mWeight[pivot]` -> a numbered vertex must keep weight one -> the weight
cannot say dead -> `mMark` must carry GONE -> `mHasNumbered` guards the load. The prepass is where
the numbered vertex comes from, so a variation taking only bucket 0 would remove the case the chain
exists to serve. It would also change the order, bucket 1 being genmmd's, so it belongs with the
other variations rather than with the alignment work.

## Absorb, reclaim and prune, and where each half lives, 2026-08-26

The elimination is three tasks. Prune is one place; absorb and reclaim each have a C side and an I
side, and the two sides are in different functions:

```
absorb    C side   C[p] = reach(p)          formReachableSet*, in beginElimination
          I side   | {p}                    the append after the incidence loop, in prune

reclaim   C side   C = C - I[p]             killClique over I[pivot], in beginElimination
          I side   I[u] - I[p]              the drop in the incidence loop, in prune

prune     one side A[u] = A[u] - C[p] - {p} the adjacency loop, in prune
```

So the incidence loop is the I half of both tasks and decides neither. Both decisions were taken in
`beginElimination`, and the loop records them in u's list. The adjacency loop is the only one of the
three that both decides and performs its task.

Two things this split does not say, and both matter.

**The incidence loop tests DEAD, not IN I[p].** Its condition is `adjacencySize != 0`, so it drops
any clique that has been retired, whoever retired it. Whether that is the same set as `I[p]` is a
separate question from where the work lives.

**AND ON THE amd BRANCH RECLAIM HAS A THIRD SITE.** `absorbAggressively` kills cliques of its own
and compacts `I[u]` itself, outside both functions above. So the two-sided picture is mmd's exactly
and amd's plus one.

## A clique that survives the prune, worked on five vertices, 2026-08-26

The incidence loop keeps some of `I[u]` and drops the rest, and the test is narrower than "is it
still useful". A clique `c` survives in `I[u]` exactly when the PIVOT IS NOT ONE OF ITS MEMBERS,
which is the same as `c` not being in `I[p]`. `reach(p)` unions the cliques `p` belongs to, so those
are the ones subsumed; a clique holding `u` and not `p` was never in that union and nothing covers
it.

Five vertices are enough to see it happen:

```
1 --- 0 --- 2 --- 3            edges: 0-1, 0-2, 2-3, 2-4, 3-4
            |    /|
            +-- 4-+            1 - 0 - 2 is a path, {2,3,4} is a triangle
```

Eliminating in the order 0, 3, 1, 2, 4:

```
eliminate 0:   C[0] = [1, 2]      cliques absorbed: []
    I[1]: []       -> [0]
    I[2]: []       -> [0]
eliminate 3:   C[3] = [2, 4]      cliques absorbed: []
    I[2]: [0]      -> [0, 3]
    I[4]: []       -> [3]
eliminate 1:   C[1] = [2]         cliques absorbed: [0]
    I[2]: [0, 3]   -> [1, 3]
eliminate 2:   C[2] = [4]         cliques absorbed: [1, 3]
    I[4]: [3]      -> [2]
eliminate 4:   C[4] = []          cliques absorbed: [2]
```

Each `eliminate p` line is `beginElimination`: it forms `C[p]` and kills the cliques of `I[p]`, so
"cliques absorbed" is `I[p]`, empty for the first two because neither pivot belonged to a clique
yet. The indented lines under it are the prune, one per `u` in `C[p]`, rewriting `I[u]` as
`( I[u] - I[p] ) | {p}`.

**PIVOT 1 IS THE CASE.** Vertex 2 names two cliques there and they part ways:

- `C[0] = {1, 2}` is ABSORBED, because 1 is a member. `reach(1)` unions it, so it becomes a subset
  of `C[1]` and is killed.
- `C[3] = {2, 4}` SURVIVES, because 1 is not a member. Nothing unioned it, and it still carries 4,
  which `C[1] = {2}` does not cover.

The survivor is carrying real information rather than being left behind by an incomplete test:
`C[3] - C[1] = {4}` is the only way vertex 2 still reaches 4 after this step. That difference is
literally the term the amd bound sums over the survivors.

The matrix, permuted to the elimination order, with F for fill:

```
A                          L + L^T
      0  3  1  2  4              0  3  1  2  4
  0   X  .  X  X  .          0   X  .  X  X  .
  3   .  X  .  X  X          3   .  X  .  X  X
  1   X  .  X  .  .          1   X  .  X  F  .
  2   X  X  .  X  X          2   X  X  F  X  X
  4   .  X  .  X  X          4   .  X  .  X  X
```

One fill entry, from eliminating 0 with its two non-adjacent neighbors 1 and 2.

**WHERE THE TEST STOPS SHORT, and it is the gap between the branches.** Survival is decided by
`p not in C[c]`, not by `C[c]` failing to be contained in `C[p]`. Those differ: a clique can have
every member inside `C[p]` while not containing `p` itself, and it is then redundant but not
absorbed. mmd keeps it. That is exactly what aggressive absorption catches on the amd branch,
testing `outside[c] == 0`, which the bound has already computed.

**PROVENANCE.** The trace above is a direct simulation of the four set-operation lines, not a run of
a production driver: the drivers choose their own pivots and would not pick this order. The
arithmetic is the algorithm's; it has not been through `MmdFlat`.

## Related

- `archive/sparse_factorization.md` section 5, the prose, pseudocode and worked examples.
- `private/Mmd.cpp` and `private/Amd.cpp`, the vendored routines these are read against.
- `src/OrderEngine.cpp`, the glue that calls them.
