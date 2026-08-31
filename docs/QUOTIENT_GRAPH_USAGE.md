# Quotient graph usage: which driver calls what

Two quotient graph classes and four drivers, and the drivers PARTITION across the classes: `MmdFlat`
and `AmdFlat` use `QuotientGraphFlat` and nothing else, `MmdCompacted` and `AmdCompacted` use
`QuotientGraphCompacted` and nothing else. So there are two tables and neither loses anything to the
other.

**A DRIVER IS AN ENGINE INSTANTIATION AS OF 2026-08-30**, `MmdEngine<QuotientGraphFlat>` and its
three siblings, the two bodies each being written once and instantiated per store. The four names
survive as the free functions the `Ordering` enum dispatches to, and the tables below still read the
same way: a row says which of the four reaches that entry point. What changed is that the two
members of a pair now reach it from ONE body, so a difference between them is no longer something a
reader has to find by diffing two files. Where a store genuinely cannot answer, the engine resolves
it by an overload beside the instantiations rather than in the body; `numCompactions` is the only
such case.

**WHAT THIS DOCUMENT IS FOR.** Not the interfaces, which are in the headers and are authoritative
there. This is the usage: which of each class's entry points each driver actually reaches, so that
what the two branches share and where they part is visible in one place. It exists to support
aligning them further, so a blank cell is a question rather than a fact.

**THE SAME FACTS APPEAR IN THREE TABLES**, one per class and one joining them by concept, and all
three carry an `align` column so progress is visible wherever you happen to be reading. They are
generated from one list and every cell is checked against the four drivers' call sites, so they
cannot drift apart; if one is edited by hand the others must be edited with it.

`QuotientGraphChained` and `MmdChained` are out of scope. Chaining lost the layout comparison and is
kept as a permanent alternative rather than a candidate, so aligning it buys nothing.

**A blank cell means the driver does not call that entry point.** It does not mean the class lacks
it: every method in the flat table exists in `QuotientGraphFlat`, and every method in the compacted
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

## The flat class: `QuotientGraphFlat`

| API | `MmdFlat` | `AmdFlat` | align | since |
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
| `advanceTagMmd` | x |  | aligned |  |
| `resetMarkAndTagMmd` | x |  | aligned | 2026-08-30 |
| `numTagResets` | x |  | aligned | 2026-08-30 |
| `markMmd`, `setMarkMmd` | x |  | aligned |  |
| `number` | x |  | aligned |  |
| `eliminatedMmd` / `eliminatedAmd` | x | x | done | 2026-08-21 |
| **elimination** | | | | |
| `eliminateMmd` / `eliminateAmd` | x | x | done | 2026-08-24 |
| `massEliminate` |  | x | aligned |  |
| **merging and absorption** | | | | |
| `merge` | x | x | aligned |  |
| `absorb` |  | x | aligned |  |
| **degree** | | | | |
| `reachableSetWeight` | x |  | aligned |  |
| **output** | | | | |
| `orderAsMerged` |  | x | aligned |  |
| `orderAscending` | x |  | aligned |  |
| **configuration** | | | | |
| `setReverseIncidence` | x |  | aligned |  |
| `setLateMassElimination` |  | x | aligned |  |
| **counters** | | | | |
| `numBornCliqueMembers` | x | x | aligned |  |
| `numPeakCliqueMembers` | x | x | aligned |  |
| `cliqueCountBalances` | x | x | aligned |  |

Counting ROWS rather than names, since a suffixed pair is one row and two entry points: thirteen are
called by both, nine by `MmdFlat` alone and eight by `AmdFlat` alone.

## The compacted class: `QuotientGraphCompacted`

Where a name is suffixed the two halves share a row, since they are the same entry point split
where the vendored routines disagree.

| API | `MmdCompacted` | `AmdCompacted` | align | since |
|---|:---:|:---:|:---:|:---:|
| **layout accessors** | | | | |
| `adjacencyMmd` / `adjacencyAmd` | x |  | done | 2026-08-28 |
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
| `advanceTagMmd` | x |  | aligned |  |
| `resetMarkAndTagMmd` | x |  | aligned | 2026-08-30 |
| `numTagResets` | x |  | aligned | 2026-08-30 |
| `markMmd`, `setMarkMmd` | x |  | aligned |  |
| `number` | x |  | aligned |  |
| `eliminatedMmd` / `eliminatedAmd` | x | x | done | 2026-08-21 |
| **elimination** | | | | |
| `eliminateMmd` / `eliminateAmd` | x | x | done | 2026-08-24 |
| `massEliminate` |  | x | aligned |  |
| **merging and absorption** | | | | |
| `merge` | x | x | aligned |  |
| `absorb` |  | x | aligned |  |
| **degree** | | | | |
| `reachableSetWeight` | x |  | aligned |  |
| **output** | | | | |
| `orderAsMerged` |  | x | aligned |  |
| `orderAscending` | x |  | aligned |  |
| **configuration** | | | | |
| `setReverseIncidence` | x |  | aligned |  |
| `setLateMassElimination` |  | x | aligned |  |
| **counters** | | | | |
| `numCompactions` | x | x | `layout` |  |
| `numBornCliqueMembers` | x | x | aligned | added 2026-08-30 |
| `numPeakCliqueMembers` | x | x | aligned |  |
| `cliqueCountBalances` | x | x | aligned |  |

## The same content joined by concept

The two tables above answer "what does this driver call". This one answers "do the two classes do
this the same way", by putting all four drivers side by side and keying the rows by CONCEPT rather
than by method name. A row whose left pair and right pair have the same shape is aligned; where
they differ, the `align` column carries the item's position in the ledger below, or `done`, or
`layout` for the one pair that stays apart by design.

Both shapes are kept deliberately, to see which is the more useful to work from.

| what | `MmdFlat` | `AmdFlat` | `MmdCompacted` | `AmdCompacted` | align | since | note |
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
| advanceTagMmd | x |  | x |  | aligned |  |  |
| tag guard | x |  | x |  | aligned | 2026-08-30 | the amd guard is the driver's, not the class's |
| tag reset count | x |  | x |  | aligned | 2026-08-30 | published on the same `ElmOrder` field |
| mark, setMarkMmd | x |  | x |  | aligned |  |  |
| number | x |  | x |  | aligned |  |  |
| eliminated | x | x | x | x | done | 2026-08-21 | a zero weight on amd, a tag on mmd |
| **elimination** | | | | | | | |
| performing an elimination | x | x | x | x | done | 2026-08-21 | one call on both |
| mass elimination |  | x |  | x | aligned |  | one method; see restoreWeight |
| **merging and absorption** | | | | | | | |
| merge | x | x | x | x | aligned |  |  |
| absorb |  | x |  | x | aligned |  |  |
| **degree** | | | | | | | |
| reachableSetWeight | x |  | x |  | aligned |  |  |
| **output** | | | | | | | |
| order |  | x |  | x | aligned |  |  |
| orderAscending | x |  | x |  | aligned |  |  |
| **configuration** | | | | | | | |
| setReverseIncidence | x |  | x |  | aligned |  |  |
| setLateMassElimination |  | x |  | x | aligned |  |  |
| which half of the run comes first |  |  |  |  | done | 2026-08-21 | the flag went |
| **counters** | | | | | | | |
| numBornCliqueMembers | x | x | x | x | aligned | 2026-08-30 | compacted keeps a counter |
| numCompactions |  |  | x | x | layout |  | only a bounded store can run out |
| numPeakCliqueMembers | x | x | x | x | aligned |  |  |
| cliqueCountBalances | x | x | x | x | aligned |  |  |

The concept rows above hide the spellings. In full, where they differ between the classes:

| concept | flat | compacted |
|---|---|---|
| adjacency accessor | `adjacency` | `adjacencyMmd`, `adjacencyAmd` |
| incidence accessor | `incidence` | `incidenceMmd`, `incidenceAmd` |
| eliminated | `eliminated` | `eliminatedMmd`, `eliminatedAmd` |
| performing an elimination | `eliminateMmd`, `eliminateAmd` | the same pair, over |
| | | `beginElimination*`, `prune*`, `finishElimination` |
| which half of the run comes first | `setAmdListOrder` | the suffixed accessors and walks |
| mark array on demand | always allocated | `enableMarks` |

## What the tables say

**Eight entry points are called by all four drivers**: `adjacencySize`, `incidenceSize`, `clique`,
`cliqueSize`, `weight`, `merge`, `numPeakCliqueMembers`, `cliqueCountBalances`. That is the floor
of what an ordering driver over a quotient graph needs whatever branch and whatever layout.

**The branch split is the same in both classes**, which is the reassuring part: marks, `number`,
`reachableSetWeight`, `orderAscending` and `setReverseIncidence` are mmd's on both sides, and
`absorb`, `trimClique`, `setAside`, `cliqueWeight`, `orderAsMerged` and `setLateMassElimination` are
amd's on both.
Moving between layouts changed the storage, not who wants what.

**And the class split WAS one idea**: the flat pair called `eliminate` where the compacted pair
spelled the three steps, `beginElimination`, prune, `finishElimination`. That was the only entry
point in either table existing on one side and not the other for a reason that is not the layout,
and it closed on 2026-08-21, both classes now carrying the wrapper. The wrapper is a suffixed PAIR
rather than an overload as of 2026-08-24, so a reader can tell from a call site which branch is
running; the three steps stay public. It went
when the classes merged, because amd's prune takes a `TaggedScan` and mmd's takes nothing.

### Five blanks that are not what they look like

Each of these was checked rather than assumed.

**CLOSED 2026-08-28. `AmdCompacted` does not call `adjacencyAmd`, AND NOW NEITHER DOES `AmdFlat`.**
Supervariable detection was the whole of it. The layout decided where the entry to skip sat: under
the amd run order the new clique is rotated to index 0, so the run minus its first entry is ONE SPAN
and the walk is a single loop off the run's base, which under that order is the incidence accessor.
The flat class used to store `[A, I]` for BOTH branches, which put the new clique in the MIDDLE of
the run and forced `AmdFlat` to stamp in two loops and name both accessors.

THE FLAT CLASS NOW CARRIES THE RUN ORDER PER BRANCH, as the compacted one always did: `[A, I]` for
mmd and `[I, A]` for amd. `AmdFlat`'s detection is a single loop over one span and the two drivers'
detection regions are BYTE-IDENTICAL. There was no reason for the flat class to serve amd's
front-of-`I` convention on mmd's run order; it was inherited rather than chosen.

**Detection is the only operation in either driver that is blind to the kind of an entry**, which is
why it is the only one the run order reaches. Everywhere else a clique id is looked up in the store
and a vertex id is weighted, so the two halves are two lists and the order is an offset. Detection
asks only whether two vertices name the same set, so the halves being one span is worth something.
It is sound because the two id populations cannot collide inside one run: a clique's id is the pivot
that formed it, and if that pivot were in `A[u]` then `u` was in its reach and the prune removed it.

**`MmdCompacted` does not call `massEliminate`.** Mmd is not late, so `finishElimination` calls it
internally. `AmdCompacted` calls it from the driver because absorption has to run first; see
`experiments/ordering/README.md`. `MmdFlat` and `AmdFlat` divide the same way, which is why the row
is `aligned`.

**`restoreWeight` and `restorePivotWeight` are compacted-only, and both classes negate.** Membership
in the clique being built is marked by flipping a weight's sign in both. The difference is only who
undoes it: the flat class restores inside `massEliminate` for both branches, while the compacted one
exposes the two accessors so amd's driver can restore in a pass it already makes over `C[pivot]`.

**`eliminated` is ONE method in the flat class and TWO in the compacted one, and the layout has
nothing to do with it.** The flat class answers "is u dead" with `mMarkMmd[u] == GONE` for both
branches, allocating the mark array at construction and writing GONE at every retirement site. The
compacted class allocates the array on demand, so the amd branch has none and answers with
`mWeight[u] == 0` instead. That is a footprint decision taken in one class and not the other, worth
n int32 on the amd branch, and it is the clearest candidate in either table for alignment: see the
section below.

**`setVendoredListOrder` is flat-only.** The compacted class has no such flag: the walk order is a
suffixed pair instead. `MmdCompacted` carried a dead copy of the flag until 2026-08-19.

**`numBornCliqueMembers` IS NO LONGER FLAT-ONLY, 2026-08-30.** It answers how many members were ever
put into a clique, which the benchmark prints as `cC`. The flat class reads it off its store's
length, which it can because nothing there is reclaimed; the compacted class keeps a counter,
incremented at its single birth site beside the live and peak counters it was already maintaining
there.

**IT WAS ADDED FOR THE CROSS-CHECK RATHER THAN THE FIGURE.** Members born is a property of the
ALGORITHM and not of the layout, exactly as the peak is: a branch's two drivers return one
permutation, so they form the same cliques with the same sizes at the same moments and MUST report
the same `cC`. That check did not exist before. The objection recorded here previously, that a
compacted pool is sized at construction and reused so its class should not pay for the figure, was
answered by looking: the increment is one line at a site that already does the add, O(1) per
elimination and not per member.

An earlier version of this method on the compacted classes reported the POOL's size instead, a
different question with no caller, and both lost it; that is a separate thing from the counter added
now. `numCompactions` remains that class's alone.

## Aligning the two classes

**THE TEST TO APPLY IS WHETHER THE LAYOUT CAUSED IT.** A difference between mmd and amd is expected
and is usually a vendored routine's; a difference between the flat class and the compacted one
should be traceable to arena against pool, and where it is not, it is an accident of how the two
were written and is worth removing. The compacted class is the flat one with positions into a
different layout, and it should not differ from it at an algorithmic level.

**Genuinely layout-driven, and to be left alone.** `numCompactions`: only a store with a bounded
pool can run out and be compacted, so an arena has no such quantity to report. THE FLAT CLASS
BRIEFLY ANSWERED IT WITH A `constexpr 0` and no longer does, that being an accessor for a thing the
class does not have; the difference is absorbed instead by an overload pair local to each engine's
unit, so the shared body asks `numCompactionsOf(qg)` and does not know which store it holds.

`numBornCliqueMembers` was listed here as the other half of this pair until 2026-08-30 and is now
aligned; see above. Nothing else in either table is of this kind.

### The ledger

A summary of the marked rows above, and deliberately redundant with them. The number in the first
column is the `align` value the tables carry.

| | difference | at `97f4bc6` | now | since |
|---|---|---|---|---|
| 1, `layout` | `arenaEntries` unreported by `AmdCompacted` | apart | apart, by design | 2026-08-23 |
| 2, `done` | `restoreWeight`, `restorePivotWeight` | apart | **aligned** | 2026-08-21 |
| 3, `done` | the `eliminate` wrapper | apart | **aligned** | 2026-08-21 |
| 4, `done` | `eliminated` split, with `enableMarks` | apart | **aligned** | 2026-08-21 |
| 5, `done` | `setVendoredListOrder` against the suffixes | apart | **aligned** | 2026-08-21 |
| `layout` | `numBornCliqueMembers` | apart | **aligned** | 2026-08-30 |
| `layout` | `numCompactions` | apart | apart, by design | |

**A closed item keeps its number and its row**, showing `done` in the tables rather than
disappearing, so the numbering never shifts under a reader who has been following it.

**The differences, in the order it makes sense to take them**, smallest and most independent first.
The number is the `align` value the tables carry.

**1. `arenaEntries` unreported by `AmdCompacted`. CLOSED 2026-08-21 AND REOPENED AND REVERSED
2026-08-23, which makes it the one item here that was closed the wrong way.** It was read as a
missing out-parameter: the other three drivers had an `Impl` taking a pointer plus two public
overloads and `AmdCompacted` had a single entry point, so `AmdCompacted` was given the same shape.

**The shape was uniform and the QUANTITY behind it never was.** The flat class returns the clique
store's size, which is what that store cost and which the benchmark prints as `cC`. The compacted
class returned the POOL's size, fixed at construction and derivable from nnz(A) without running
anything. The compacted class's own comment said as much, "Not the same quantity and not meant to
be", which is a comment apologizing for a name.

**And nothing ever called it.** `matrix_ordering.cpp` reads `cC` from the flat pair; the compacted
drivers are reached through the plain overload. So the shape bought here had no consumer on the day
it landed or since.

**The correct reading is that only a store that GROWS has to pay for what it holds.** The compacted
pool is reused, and a chained clique lives in its own pivot's dead segment and costs nothing extra,
so neither has a clique-storage quantity that varies. That is the property the compacted layout was
chosen for rather than an omission. `QuotientGraphCompacted::arenaEntries` is deleted, with the
`Impl` and the two public overloads on both compacted drivers, each of which then collapsed to one
entry point; `numCompactions` is that class's whole storage figure and says whether the fixed pool
sufficed.

**AND THE FLAT METHOD IS NOW `numBornCliqueMembers`, renamed the same day.** `arenaEntries` named
neither a member that still exists nor the quantity it returns: every entry of the store is a clique
member and nothing else is ever written there, so the count is members born, exactly. That makes it
the third of a family, born, live and peak, differing only in WHEN rather than in what is counted.
It reads the store's length rather than a counter, which is a coincidence of nothing ever being
reclaimed and would need a counter the day that changes.

**Which also means the other two classes COULD publish it**, the birth site being in all three, and
it would be identical across a branch's pair by construction, exactly as `numPeakCliqueMembers` is.
That is not ruled out by any of this. What is ruled out is a clique-STORAGE column for a layout
whose clique storage does not vary.

**2. `restoreWeight` and `restorePivotWeight`. CLOSED 2026-08-21.** Both classes negate weights to
mark membership in the clique being built, and the question was who undoes it. The answer is
neither branch nor class: WHOEVER RUNS MASS ELIMINATION restores, the eliminator when it runs
eagerly and the driver when it runs late and has a refile pass to ride the store on. `AMD_2` is the
late case, mass-eliminating with the weights still negative and restoring under RESTORE DEGREE
LISTS, and `AmdFlat` now matches it as `AmdCompacted` already did.

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
answered "is u dead" with `mMarkMmd[u] == GONE` for both branches and allocated the array always; it
now does what the compacted class does, a zero weight on the amd branch and the tag on the mmd one,
with the array allocated by `enableMarks`. `AmdFlat` drops n int32 it never read.

**The two predicates are NOT equivalent, which is why this needed checking rather than renaming.**
They differ on a vertex `number` retired, which is mmd's alone, and on a PIVOT, which is retired
with its weight intact because that weight is the supervariable's and `orderAsMerged` needs it. So
the amd predicate is correct only if the amd driver never asks about a pivot. It does not, a pivot
being unfiled when chosen and never revisited, and that was established by asserting the two agree
on every call `AmdFlat` makes and running it: the digest's 73 grids at `-O0`, cubes to 33 a side,
random patterns at degree 6, 12 and 40, a star, a diagonal, a dense block that fires the dense-row
rule, and `test_order`. No disagreement anywhere.

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
identical bodies; the chained class has one layout AND one driver, `MmdChained`, so it carries the
mmd half alone and there is nothing there to duplicate. **The names are split so that a driver and
its counterpart read the same**, which is what they are bought for, the flat class having no layout
difference of its own to record.

**Every step must leave every permutation unmoved**, so `make digest` and the two alignment checks
are the gate at each one, and the baseline below is the second gate.

**FOUR OF THE FIVE ARE CLOSED AS OF 2026-08-21, AND THE FIFTH WAS REVERSED ON 2026-08-23**, item 1
having been a difference the layout causes rather than an accident; see its entry. What remains
between the two classes is `numCompactions` alone, which is layout: `numBornCliqueMembers` joined it
in that row until 2026-08-30 and is now on both classes. The two
`adjacencyAmd` calls in `AmdFlat`'s hash-detection block are GONE as of 2026-08-28: the flat class
took the amd run order and the two detection regions are now byte-identical.

**THE COUNTS THIS PARAGRAPH CARRIED WERE TOO HIGH, corrected 2026-08-21.** It said the mmd pair
differs in three calls of about forty and the amd pair in nine of about thirty-eight. Both figures
came from extracting `qg.<method>` from the sources without stripping comments first, so prose
naming a method counted as a call:

```
                         comments in    comments out
MmdFlat vs MmdCompacted        3 calls        0 calls
AmdFlat vs AmdCompacted        9 calls        2 calls
```

The three mmd differences are `qg.mark(v)`, `qg.weight(v)` and `qg.eliminatedMmd(v)` written into
comments in `MmdFlat.cpp` that `MmdCompacted` does not carry, and one of the nine is
`qg.cliqueBase()` named in the comment that records its removal. **Stripped of comments the mmd
pair's call sequences are IDENTICAL apart from `numCompactions`, and the amd pair's differed only in
the two `adjacencyAmd` calls named above and the order of two adjacent rejects.** So the sentence
above the table was right and the arithmetic beside it was not. Both of those amd differences are
closed as of 2026-08-28, the rejects having been reordered and the run order flipped.

**The tail residual recorded here is CLOSED as of 2026-08-30.** Three of the four drivers used to
publish their counters in one order and `MmdCompacted` in another; there is now no tail to differ,
the four drivers being two engine instantiations apiece and the publication being four assignments
in one shared body:

```
eo.mOrder                = qg.orderAscending(pivots);   // orderAsMerged on the amd branch
eo.mNumPeakCliqueMembers = qg.numPeakCliqueMembers();
eo.mNumBornCliqueMembers = qg.numBornCliqueMembers();
eo.mNumCompactions       = numCompactionsOf(qg);
```

The counters left the drivers' signatures and three globals on the same day and now ride on
`ElmOrder`, the object an ordering produces.

### The current reading, 2026-08-24, and it supersedes the baseline below

Taken at `937fb30` on alpamayo, `make scale2d` and `make scale3d`, TWO consecutive runs so that each
figure carries a spread rather than a point.

| | 801^2 | 1025^2 | 1601^2 | 65^3 | 81^3 |
|---|---:|---:|---:|---:|---:|
| `MmdCompacted / MmdFlat` | 0.931-0.939 | 0.927-0.935 | 0.914-0.916 | 0.948-0.955 | 0.920-0.938 |
| `AmdCompacted / AmdFlat` | 0.912-0.927 | 0.935-0.960 | 0.955-0.995 | 0.916-0.930 | 0.943-0.957 |
| `MmdChained / MmdFlat` | 1.135-1.141 | 1.221-1.224 | 1.157 | 1.316-1.354 | 1.503-1.566 |

**THE MMD PAIR AND THE CHAINED CONTROL HAVE NOT MOVED** since 2026-08-21, both sitting inside the
harness floor against the table below. **The AMD PAIR HAS**, from 0.83-0.90 to 0.91-0.99, and the
cause is the denominator: `AmdFlat` got faster over the six commits between the two tables, where
`AmdCompacted` did not. **Which commit is not established here**, and this table cannot say: it
spans `1da85c5` through `937fb30`.

**AND THE RUN-TO-RUN SPREAD IS MEASURED RATHER THAN ASSUMED**, which is the more useful half. Two
runs of the same binaries on the same machine, per driver over the five sizes:

```
MmdCorrected   0.957 to 1.000        AmdVendored    0.953 to 1.033
MmdFlat        0.982 to 1.043        AmdFlat        0.996 to 1.006
MmdChained     0.977 to 1.031        AmdCompacted   0.956 to 1.013
MmdCompacted   0.986 to 1.023
```

So the floor is 4 to 8 per cent on most columns, which is wider than the 3 per cent this file has
been quoting, and the vendored columns are among the noisiest, which is the standing caution in
`docs/NEXT.md` about quoting `AmdFlat / AMD` per row. `AmdFlat`'s one per cent is the outlier and
should not be read as the general case.

**A RATIO BETWEEN TWO OF OUR OWN DRIVERS IS THE THING TO WATCH**, since both move together with the
machine and the common term cancels. A ratio against a vendored routine carries that routine's noise
in full.

### Baseline before alignment, 2026-08-21

**HISTORICAL, and no longer a control.** It was taken at `97f4bc6`, six commits back, and the tree
has moved for reasons that are not alignment steps; the current reading is above. Kept because a
dated measurement is a record of a run.

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
