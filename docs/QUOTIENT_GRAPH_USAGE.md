# Quotient graph usage: which driver calls what

Two quotient graph classes and four drivers, and the drivers PARTITION across the classes: `Mmd3`
and `Amd3` use `QuotientGraph` and nothing else, `Mmd3C` and `Amd3B` use `QuotientGraphCompacted`
and nothing else. So there are two tables and neither loses anything to the other.

**WHAT THIS DOCUMENT IS FOR.** Not the interfaces, which are in the headers and are authoritative
there. This is the usage: which of each class's entry points each driver actually reaches, so that
what the two branches share and where they part is visible in one place. It exists to support
aligning them further, so a blank cell is a question rather than a fact.

**THE SAME FACTS APPEAR IN THREE TABLES**, one per class and one joining them by concept, and all
three carry an `align` column so progress is visible wherever you happen to be reading. They are
generated from one list and every cell is checked against the four drivers' call sites, so they
cannot drift apart; if one is edited by hand the others must be edited with it.

`QuotientGraphChained` and `Mmd3B` are out of scope. Chaining lost the layout comparison and is kept
as a permanent alternative rather than a candidate, so aligning it buys nothing.

**A blank cell means the driver does not call that entry point.** It does not mean the class lacks
it: every method in the flat table exists in `QuotientGraph`, and every method in the compacted
table exists in `QuotientGraphCompacted`.

**THE `align` COLUMN, in all three tables.**

| value | meaning |
|---|---|
| `aligned` | the two classes already agreed when this started |
| `done` | brought into agreement by this process; `since` carries the date |
| `1` to `5` | pending, and the number is its position in the ledger's order below |
| `layout` | stays apart, the storage genuinely differs |

A number or `layout` on a row means the FLAT and COMPACTED classes differ there. It says nothing
about mmd against amd: a row where one branch calls something and the other does not is `aligned`
as long as both classes agree about it, which is the usual case and the expected one.

## The flat class: `QuotientGraph`

| API | `Mmd3` | `Amd3` | align | since |
|---|:---:|:---:|:---:|:---:|
| **layout accessors** | | | | |
| `adjacencyMmd` / `adjacencyAmd` | x | x | done | 2026-08-21 |
| `incidenceMmd` / `incidenceAmd` | x | x | done | 2026-08-21 |
| `adjacencySize` | x | x | aligned |  |
| `incidenceSize` | x | x | aligned |  |
| **cliques** | | | | |
| `clique` | x | x | aligned |  |
| `cliqueSize` | x | x | aligned |  |
| `cliqueWeight` |  | x | aligned |  |
| `trimClique` |  | x | aligned |  |
| **weights** | | | | |
| `weight` | x | x | aligned |  |
| `restoreWeight` / `restorePivotWeight` |  | x | done | 2026-08-21 |
| `setAside` |  | x | aligned |  |
| **marks** | | | | |
| `enableMarks` | x |  | done | 2026-08-21 |
| `advanceTag` | x |  | aligned |  |
| `mark`, `setMark` | x |  | aligned |  |
| `number` | x |  | aligned |  |
| `eliminatedMmd` / `eliminatedAmd` | x | x | done | 2026-08-21 |
| **elimination** | | | | |
| `eliminate` | x | x | done | 2026-08-21 |
| `massEliminate` |  | x | aligned |  |
| **merging and absorption** | | | | |
| `merge` | x | x | aligned |  |
| `absorb` |  | x | aligned |  |
| **degree** | | | | |
| `reachableWeight` | x |  | aligned |  |
| **output** | | | | |
| `order` |  | x | aligned |  |
| `orderAscending` | x |  | aligned |  |
| **configuration** | | | | |
| `setReverseIncidence` | x |  | aligned |  |
| `setLateMassElimination` |  | x | aligned |  |
| **counters** | | | | |
| `arenaEntries` | x | x | done | 2026-08-21 |
| `numPeakCliqueMembers` | x | x | aligned |  |
| `cliqueCountBalances` | x | x | aligned |  |

Fourteen entry points are called by both, six by `Mmd3` alone and eight by `Amd3` alone.

## The compacted class: `QuotientGraphCompacted`

Where a name is suffixed the two halves share a row, since they are the same entry point split
where the vendored routines disagree.

| API | `Mmd3C` | `Amd3B` | align | since |
|---|:---:|:---:|:---:|:---:|
| **layout accessors** | | | | |
| `adjacencyMmd` / `adjacencyAmd` | x |  | done | 2026-08-21 |
| `incidenceMmd` / `incidenceAmd` | x | x | done | 2026-08-21 |
| `adjacencySize` | x | x | aligned |  |
| `incidenceSize` | x | x | aligned |  |
| **cliques** | | | | |
| `clique` | x | x | aligned |  |
| `cliqueSize` | x | x | aligned |  |
| `cliqueWeight` |  | x | aligned |  |
| `trimClique` |  | x | aligned |  |
| **weights** | | | | |
| `weight` | x | x | aligned |  |
| `restoreWeight` / `restorePivotWeight` |  | x | done | 2026-08-21 |
| `setAside` |  | x | aligned |  |
| **marks** | | | | |
| `enableMarks` | x |  | done | 2026-08-21 |
| `advanceTag` | x |  | aligned |  |
| `mark`, `setMark` | x |  | aligned |  |
| `number` | x |  | aligned |  |
| `eliminatedMmd` / `eliminatedAmd` | x | x | done | 2026-08-21 |
| **elimination** | | | | |
| `eliminate` | x | x | done | 2026-08-21 |
| `massEliminate` |  | x | aligned |  |
| **merging and absorption** | | | | |
| `merge` | x | x | aligned |  |
| `absorb` |  | x | aligned |  |
| **degree** | | | | |
| `reachableWeight` | x |  | aligned |  |
| **output** | | | | |
| `order` |  | x | aligned |  |
| `orderAscending` | x |  | aligned |  |
| **configuration** | | | | |
| `setReverseIncidence` | x |  | aligned |  |
| `setLateMassElimination` |  | x | aligned |  |
| **counters** | | | | |
| `arenaEntries` | x | x | done | 2026-08-21 |
| `compactions` | x | x | layout |  |
| `numPeakCliqueMembers` | x | x | aligned |  |
| `cliqueCountBalances` | x | x | aligned |  |

## The same content joined by concept

The two tables above answer "what does this driver call". This one answers "do the two classes do
this the same way", by putting all four drivers side by side and keying the rows by CONCEPT rather
than by method name. A row whose left pair and right pair have the same shape is aligned; where
they differ, the `align` column carries the item's position in the ledger below, or `done`, or
`layout` for the one pair that stays apart by design.

Both shapes are kept deliberately, to see which is the more useful to work from.

| what | `Mmd3` | `Amd3` | `Mmd3C` | `Amd3B` | align | since | note |
|---|:---:|:---:|:---:|:---:|:---:|:---:|---|
| **layout accessors** | | | | | | | |
| adjacency accessor | x | x | x |  | done | 2026-08-21 | identical in flat, not in compacted |
| incidence accessor | x | x | x | x | done | 2026-08-21 | identical in flat, not in compacted |
| adjacencySize | x | x | x | x | aligned |  |  |
| incidenceSize | x | x | x | x | aligned |  |  |
| **cliques** | | | | | | | |
| clique | x | x | x | x | aligned |  |  |
| cliqueSize | x | x | x | x | aligned |  |  |
| cliqueWeight |  | x |  | x | aligned |  |  |
| trimClique |  | x |  | x | aligned |  |  |
| **weights** | | | | | | | |
| weight | x | x | x | x | aligned |  |  |
| restoring negated weights |  | x |  | x | done | 2026-08-21 | a late run leaves it to the driver |
| setAside |  | x |  | x | aligned |  |  |
| **marks** | | | | | | | |
| mark array on demand | x |  | x |  | done | 2026-08-21 | amd allocates none |
| advanceTag | x |  | x |  | aligned |  |  |
| mark, setMark | x |  | x |  | aligned |  |  |
| number | x |  | x |  | aligned |  |  |
| eliminated | x | x | x | x | done | 2026-08-21 | a zero weight on amd, a tag on mmd |
| **elimination** | | | | | | | |
| performing an elimination | x | x | x | x | done | 2026-08-21 | one call on both |
| mass elimination |  | x |  | x | aligned |  | one method; see restoreWeight |
| **merging and absorption** | | | | | | | |
| merge | x | x | x | x | aligned |  |  |
| absorb |  | x |  | x | aligned |  |  |
| **degree** | | | | | | | |
| reachableWeight | x |  | x |  | aligned |  |  |
| **output** | | | | | | | |
| order |  | x |  | x | aligned |  |  |
| orderAscending | x |  | x |  | aligned |  |  |
| **configuration** | | | | | | | |
| setReverseIncidence | x |  | x |  | aligned |  |  |
| setLateMassElimination |  | x |  | x | aligned |  |  |
| which half of the run comes first |  |  |  |  | done | 2026-08-21 | the flag went |
| **counters** | | | | | | | |
| arenaEntries | x | x | x | x | done | 2026-08-21 | closed 2026-08-21 |
| compactions |  |  | x | x | layout |  | the two stores report different things |
| numPeakCliqueMembers | x | x | x | x | aligned |  |  |
| cliqueCountBalances | x | x | x | x | aligned |  |  |

The concept rows above hide the spellings. In full, where they differ between the classes:

| concept | flat | compacted |
|---|---|---|
| adjacency accessor | `adjacency` | `adjacencyMmd`, `adjacencyAmd` |
| incidence accessor | `incidence` | `incidenceMmd`, `incidenceAmd` |
| eliminated | `eliminated` | `eliminatedMmd`, `eliminatedAmd` |
| performing an elimination | `eliminate` | `beginEliminationMmd`, `beginEliminationAmd`, |
| | | `pruneMmd`, `pruneAmd`, `finishElimination` |
| which half of the run comes first | `setVendoredListOrder` | the suffixed accessors and walks |
| mark array on demand | always allocated | `enableMarks` |

## What the tables say

**Eight entry points are called by all four drivers**: `adjacencySize`, `incidenceSize`, `clique`,
`cliqueSize`, `weight`, `merge`, `numPeakCliqueMembers`, `cliqueCountBalances`. That is the floor
of what an ordering driver over a quotient graph needs whatever branch and whatever layout.

**The branch split is the same in both classes**, which is the reassuring part: marks, `number`,
`reachableWeight`, `orderAscending` and `setReverseIncidence` are mmd's on both sides, and `absorb`,
`trimClique`, `setAside`, `cliqueWeight`, `order` and `setLateMassElimination` are amd's on both.
Moving between layouts changed the storage, not who wants what.

**And the class split is one idea**: the flat pair calls `eliminate`, the compacted pair spells the
three steps, `beginElimination`, prune, `finishElimination`. That wrapper is the only entry point in
either table that exists on one side and not the other for a reason that is not the layout. It went
when the classes merged, because amd's prune takes a `TaggedScan` and mmd's takes nothing.

### Five blanks that are not what they look like

Each of these was checked rather than assumed.

**`Amd3B` does not call `adjacencyAmd`.** SUPERVARIABLE DETECTION IS THE WHOLE OF IT, and this entry
said something else until 2026-08-21: it blamed the prune moving into the class, which is not the
reason and was probably true of an earlier arrangement. Both drivers reach the prune through one
`eliminate` call today and neither walks an adjacency in it. `Amd3`'s only two `adjacencyAmd` calls
in the file are in the hash block.

What the layout decides there is where the entry to skip sits. Under `AMD_2`'s order the new clique
is rotated to index 0, so the run minus its first entry is ONE SPAN and the walk is a single loop
off the run's base, which under that order is the incidence accessor. The flat class puts the new
clique at the front of `I[u]`, which is the MIDDLE of the run, so `Amd3` stamps in two loops and
names both accessors. The same asymmetry gives `Amd3` a second `adjacencyAmd` on the candidate side.

**Detection is the only operation in either driver that is blind to the kind of an entry**, which is
why it is the only one the run order reaches. Everywhere else a clique id is looked up in the store
and a vertex id is weighted, so the two halves are two lists and the order is an offset. Detection
asks only whether two vertices name the same set, so the halves being one span is worth something.
It is sound because the two id populations cannot collide inside one run: a clique's id is the pivot
that formed it, and if that pivot were in `A[u]` then `u` was in its reach and the prune removed it.

**`Mmd3C` does not call `massEliminate`.** Mmd is not late, so `finishElimination` calls it
internally. `Amd3B` calls it from the driver because absorption has to run first; see
`experiments/ordering/README.md`. `Mmd3` and `Amd3` divide the same way, which is why the row is
`aligned`.

**`restoreWeight` and `restorePivotWeight` are compacted-only, and both classes negate.** Membership
in the clique being built is marked by flipping a weight's sign in both. The difference is only who
undoes it: the flat class restores inside `massEliminate` for both branches, while the compacted one
exposes the two accessors so amd's driver can restore in a pass it already makes over `C[pivot]`.

**`eliminated` is ONE method in the flat class and TWO in the compacted one, and the layout has
nothing to do with it.** The flat class answers "is u dead" with `mMark[u] == GONE` for both
branches, allocating the mark array at construction and writing GONE at every retirement site. The
compacted class allocates the array on demand, so the amd branch has none and answers with
`mWeight[u] == 0` instead. That is a footprint decision taken in one class and not the other, worth
n int32 on the amd branch, and it is the clearest candidate in either table for alignment: see the
section below.

**`setVendoredListOrder` is flat-only.** The compacted class has no such flag: the walk order is a
suffixed pair instead. `Mmd3C` carried a dead copy of the flag until 2026-08-19.

**`arenaEntries` was called by three of the four and not by `Amd3B` until 2026-08-21.** That was an
oversight rather than a design difference: `orderAmd3B` had no out-parameter for it where the other
three drivers all had one. Closed; the driver now has the `Impl`-plus-two-overloads shape the others
use.

## Aligning the two classes

**THE TEST TO APPLY IS WHETHER THE LAYOUT CAUSED IT.** A difference between mmd and amd is expected
and is usually a vendored routine's; a difference between the flat class and the compacted one
should be traceable to arena against pool, and where it is not, it is an accident of how the two
were written and is worth removing. The compacted class is the flat one with positions into a
different layout, and it should not differ from it at an algorithmic level.

**Genuinely layout-driven, and to be left alone.** `arenaEntries` against `compactions`: one class
has a store that only grows and the other a workspace that is compacted, so they report different
quantities. Nothing else in either table is of this kind.

### The ledger

A summary of the marked rows above, and deliberately redundant with them. The number in the first
column is the `align` value the tables carry.

| | difference | at `97f4bc6` | now | since |
|---|---|---|---|---|
| 1, `done` | `arenaEntries` unreported by `Amd3B` | apart | **aligned** | 2026-08-21 |
| 2, `done` | `restoreWeight`, `restorePivotWeight` | apart | **aligned** | 2026-08-21 |
| 3, `done` | the `eliminate` wrapper | apart | **aligned** | 2026-08-21 |
| 4, `done` | `eliminated` split, with `enableMarks` | apart | **aligned** | 2026-08-21 |
| 5, `done` | `setVendoredListOrder` against the suffixes | apart | **aligned** | 2026-08-21 |
| `layout` | `arenaEntries` against `compactions` | apart | apart, by design | |

**A closed item keeps its number and its row**, showing `done` in the tables rather than
disappearing, so the numbering never shifts under a reader who has been following it.

**The differences, in the order it makes sense to take them**, smallest and most independent first.
The number is the `align` value the tables carry.

**1. `arenaEntries` unreported by `Amd3B`. CLOSED 2026-08-21.** A missing out-parameter and simply
an oversight: the other three drivers all have an `Impl` taking a pointer plus two public overloads,
and `Amd3B` had a single entry point. It now has the same shape. Both compacted drivers report a
pool of 9088 entries on a 40-square grid, which is the check that the figure means the same thing.

**2. `restoreWeight` and `restorePivotWeight`. CLOSED 2026-08-21.** Both classes negate weights to
mark membership in the clique being built, and the question was who undoes it. The answer is
neither branch nor class: WHOEVER RUNS MASS ELIMINATION restores, the eliminator when it runs
eagerly and the driver when it runs late and has a refile pass to ride the store on. `AMD_2` is the
late case, mass-eliminating with the weights still negative and restoring under RESTORE DEGREE
LISTS, and `Amd3` now matches it as `Amd3B` already did.

The rule is `mLateMassElimination`, which already existed and already meant exactly this, so both
classes now have ONE `massEliminate` and the compacted class LOSES a suffixed pair. Its
`finishElimination` collapsed with it, the two halves having differed only by a `markGone` that is
a no-op when the mark array is absent.

**A first attempt split the flat `massEliminate` by branch instead, and `Amd2` segfaulted.** That
was the useful failure: `Amd1` and `Amd2` are amd drivers that mass-eliminate EAGERLY and need the
restore, so amd-against-mmd was the wrong axis and eager-against-late was the right one.

**3. The `eliminate` wrapper. CLOSED 2026-08-21.** The flat class had one, by overloading on the
scan; the compacted class asked its drivers to spell the three calls. The compacted class now has
the same pair, so all four drivers make one call. The three steps stay public and stay suffixed:
the wrapper is a sequence rather than a replacement, and its value is that the ORDER lives in the
class where nothing in a driver could enforce it.

**4. `eliminated` split in two, and `enableMarks` with it. CLOSED 2026-08-21.** The flat class
answered "is u dead" with `mMark[u] == GONE` for both branches and allocated the array always; it
now does what the compacted class does, a zero weight on the amd branch and the tag on the mmd one,
with the array allocated by `enableMarks`. `Amd3` drops n int32 it never read.

**The two predicates are NOT equivalent, which is why this needed checking rather than renaming.**
They differ on a vertex `number` retired, which is mmd's alone, and on a PIVOT, which is retired
with its weight intact because that weight is the supervariable's and `order` needs it. So the amd
predicate is correct only if the amd driver never asks about a pivot. It does not, a pivot being
unfiled when chosen and never revisited, and that was established by asserting the two agree on
every call `Amd3` makes and running it: the digest's 73 grids at `-O0`, cubes to 33 a side, random
patterns at degree 6, 12 and 40, a star, a diagonal, a dense block that fires the dense-row rule,
and `test_order`. No disagreement anywhere.

**5. `setVendoredListOrder` against the suffixed accessors. CLOSED 2026-08-21.** The flat class
selected the branch with a flag where the compacted one names it. The flag is deleted:
`reachableSet` and `beginElimination` are suffixed pairs there now, and the prune's four flag tests
resolved statically, the flag having been known per `eliminate` overload all along. That also
removed a `heldVertex` in the mmd prune that could never be set.

**It was unblocked by retiring the earlier ladder layers**, which is why it went last rather than
first. Until then the flat class served SIX drivers across THREE list-order conventions, the third
being ours and predating both vendored ones, so a two-way split had nothing to offer `Mmd1`,
`Mmd2`, `Amd1` and `Amd2`. With those four in `retired/` the flat class sees two conventions, the
same two the compacted class does.

**And the accessors were suffixed in all three classes**, `adjacencyAmd` and `adjacencyMmd` with
the incidence twins. Only the compacted class's halves differ, since it reproduces `AMD_2`'s run
order and genmmd's. The flat class has one physical layout, so its four names sit over two
identical bodies; the chained class has one layout AND one driver, `Mmd3B`, so it carries the mmd
half alone and there is nothing there to duplicate. **The names are split so that a driver and its
counterpart read the same**, which is what they are bought for, the flat class having no layout
difference of its own to record.

**Every step must leave every permutation unmoved**, so `make digest` and the two alignment checks
are the gate at each one, and the baseline below is the second gate.

**ALL FIVE ARE CLOSED AS OF 2026-08-21.** What remains between the two classes is `compactions`
against `arenaEntries`, which is layout, and two `adjacencyAmd` calls in `Amd3`'s hash-detection
block that `Amd3B` does not make.

**THE COUNTS THIS PARAGRAPH CARRIED WERE TOO HIGH, corrected 2026-08-21.** It said the mmd pair
differs in three calls of about forty and the amd pair in nine of about thirty-eight. Both figures
came from extracting `qg.<method>` from the sources without stripping comments first, so prose
naming a method counted as a call:

```
                comments in    comments out
Mmd3 vs Mmd3C        3 calls      0 calls
Amd3 vs Amd3B        9 calls      2 calls
```

The three mmd differences are `qg.mark(v)`, `qg.weight(v)` and `qg.eliminatedMmd(v)` written into
comments in `Mmd3.cpp` that `Mmd3C` does not carry, and one of the nine is `qg.cliqueBase()` named
in the comment that records its removal. **Stripped of comments the mmd pair's call sequences are
IDENTICAL apart from `compactions`, and the amd pair's differ only in the two `adjacencyAmd` calls
named above and the order of two adjacent rejects.** So the sentence above the table was right and
the arithmetic beside it was not.

**One residual the tables do not show, and it is the tail.** Three of the four drivers publish
their counters in the same order and `Mmd3C` does not:

```
Mmd3    balances, peak, arena
Amd3    balances, peak, arena
Amd3B   balances, compactions, peak, arena
Mmd3C   balances, arena, compactions, peak
```

`Amd3B` is `Amd3` with one line inserted; `Mmd3C` is neither `Mmd3` with one line inserted nor a
match for `Amd3B`. Nothing observable turns on it, which is exactly what the five closed items had
in common.

### Baseline before alignment, 2026-08-21

Taken at `97f4bc6`, before any alignment step, on alpamayo, all drivers in the same translation unit
as their quotient graph. **Alignment changes nothing a permutation can see, so these ratios should
not move.** A drift outside the harness floor, which is about 3 percent run to run and closer to 5
for a cross-build comparison, says a step changed something it should not have.

| | 801^2 | 1025^2 | 1601^2 | 65^3 | 81^3 |
|---|---:|---:|---:|---:|---:|
| `MmdCompacted / MmdFlat` | 0.938 | 0.927 | 0.941 | 0.959 | 0.911 |
| `AmdCompacted / AmdFlat` | 0.850 | 0.840 | 0.831 | 0.895 | 0.897 |
| `MmdChained / MmdFlat` | 1.216 | 1.194 | 1.207 | 1.357 | 1.542 |

The absolute milliseconds behind them, from the same run:

| | 801^2 | 1025^2 | 1601^2 | 65^3 | 81^3 |
|---|---:|---:|---:|---:|---:|
| `MmdFlat` | 44.81 | 92.69 | 218.19 | 145.92 | 306.56 |
| `MmdCompacted` | 42.05 | 85.94 | 205.27 | 139.96 | 279.17 |
| `MmdChained` | 54.49 | 110.63 | 263.34 | 198.01 | 472.79 |
| `AmdFlat` | 78.39 | 164.56 | 358.09 | 93.15 | 207.41 |
| `AmdCompacted` | 66.67 | 138.16 | 297.43 | 83.36 | 186.00 |

`MmdChained` is here as a control rather than a candidate: chaining is kept as a permanent
alternative and is out of scope for this exercise, so its column should not move at all.

**The real-matrix figures are the other half of this baseline** and are in
`benchmarks/matrices/ORDERING.md`, measured 2026-08-20: `MmdCompacted / MmdFlat` at 0.982 median and
`AmdCompacted / AmdFlat` at 0.950, both over 246 matrices. Those are the ones to quote; the grid
table above is for catching drift between steps, being cheap to rerun.
