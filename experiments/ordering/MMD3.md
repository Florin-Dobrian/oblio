# MMD3: the alignment, iteration by iteration

A full record of how `mmd3` was brought from `mmd2` to return `genmmd`'s permutation EXACTLY, on
2026-08-07. Written for the next person doing this to another layer, so it includes what each step
ESTABLISHED, what it DISCOVERED, what was DECIDED and why, and the exact changes made. The wrong
turns are here at full length, because two of them taught more than the fixes.

Companion documents:

- `experiments/ordering/README.md`, section "Aligning a layer against a vendored routine", is the
  method stripped of this narrative.
- `experiments/ordering/README.md`, section "mmd3, and the alignment ledger", is the durable record
  and holds the authoritative ledger.
- `experiments/ordering/REPORT.md` holds what the alignment then made measurable.

---

## Iteration 0: the starting position, and why a new layer

**Established.** `mmd2` has every mechanism `genmmd` has. It carries the prepass, the q2h pair
merge, outmatched marking, and the filing convention. It had been built one pass at a time against
the vendored source and was believed complete.

**And it returned a different permutation.** Against the vendored MMD on the seven examples it
matched 2 of 7. On grids it diverged at pivot 60 of 1024 at 32 a side. Fill sat 12 to 25 percent
above genmmd, growing with `n`:

```
grid        n      MMD1 fill   MMD2 fill
32x32     1024        1.3%        1.6%
64x64     4096       13.4%       12.0%
100x100  10000       19.8%       16.2%
140x140  19600       19.4%       21.6%
280x280  78400       23.3%       21.2%
400x400 160000       21.7%       25.2%
```

**The problem this poses, and it is the reason for the whole exercise.** Minimum degree is a
tie-break algorithm: at almost every iteration several vertices share the least degree, and which
is taken is decided by whatever the data structure hands over first. So a fill gap between two
implementations measures TWO things at once, a difference of mechanism and a difference of
arbitrary choice, with no way to separate them. Nobody could say how much of the 20 percent was a
missing feature.

**Decided: build a new layer, `mmd3`, rather than change `mmd2`.** `mmd2` is a rung on a teaching
ladder and its job is to hold genmmd's mechanisms one pass at a time. A layer whose job is to match
genmmd's PERMUTATION is a different job, and the ladder's discipline is one idea per layer. `mmd3`
would add no mechanism at all.

**Also established, and it matters later.** `mmd1` diverged at pivot 4 at every size from 3x3 up.
That is expected and was not chased: `mmd1` lacks the mechanisms by design. The layer that should
have matched and did not was `mmd2`.

---

## Iteration 1: choosing the smallest case, and the first entry

**Method.** Smallest divergence first. The seven examples run from `n = 4` to `n = 12`, far smaller
than any grid.

```
graph    n    mmd1 vs MMD    mmd2 vs MMD
graph1   4    MATCHES        diverges at 2
graph2   6    at 3           at 3
graph3   12   at 1           at 6
graph4   8    at 4           at 4
graph5   5    at 2           MATCHES
graph6   6    at 2           MATCHES
graph7   5    at 2           at 2
```

**Discovered, and it shaped the diagnosis.** `mmd2` is not uniformly closer than `mmd1`. It matches
where mmd1 does not on graph5 and graph6, and on **graph1 mmd1 matches and mmd2 does not**. A
missing mechanism could not do that. A tie-break could.

**The case chosen: `graph1`, the 4-cycle, `n = 4`, diverging at pivot 2.** Four vertices, readable
by hand in full, and the one case where the simpler layer gets it right.

**What the traces showed.** Both eliminate 3 then 1, forming two cliques over `{0, 2}`. Then:

```
mmd1   iteration 1: eliminate 2, merged vertices: 0        order [3, 1, 2, 0]
mmd2   pair merge 2 into 0 during the refresh, then
       iteration 1: eliminate 0, size 2                    order [3, 1, 0, 2]
vendored MMD                                               order [3, 1, 2, 0]
```

mmd1 reaches the answer by MASS elimination, folding 0 into pivot 2. mmd2's q2h pair merge fires
first and folds 2 into 0, so 0 survives and is eliminated next. Both merges are correct: the pair
really is indistinguishable. **They differ in which vertex SURVIVES**, and the survivor is the one
that stays in the bucket and gets picked later.

**Root cause, from the vendored source.** `mmdupd` builds its q2h list as a linked list pushed at
the head, and reads from the head:

```c
if (bwd[nb] == 0) { if (fwd[nb] != 2) { list[nb] = qxh; qxh = nb; }
                    else              { list[nb] = q2h; q2h = nb; } }
...
en = q2h;                              /* the head: the LAST one pushed */
...
qsize[en] += qsize[nd];                /* en survives, nd is absorbed */
```

So the member seen LAST is processed FIRST. Our `q2h` is a `std::vector` built with `append` and
walked forward, so we process the member seen FIRST. On graph1 the element members are `0, 2`:
genmmd takes `en = 2` and absorbs 0; we take `u = 0` and absorb 2.

**Change: walk `q2h` backwards.** `reversed(q2h)` in Python, `rbegin()/rend()` in C++.

**Result: 3 of 7 examples**, up from 2. graph1 now matches. graph7 got WORSE, moving from a
divergence at pivot 2 to a longer disagreement, which is the normal signature of a tie-break
change: it reshuffles which arbitrary choice you make rather than making them right.

**Established: the defect class is a list traversal direction, not a missing mechanism.**

---

## Iteration 2: the same idiom, one level up

**Method.** graph2 still diverged at the same pivot, 3, so something upstream of the merge.

**What the trace showed.** At the divergent iteration our bucket held `3: [0 3 4]`, three vertices
of equal degree, and we took the head, 0. genmmd took 4. The refresh order had been `4, 3, 0`, and
filing pushes at the head, so we end with 0 at the head; genmmd must be filing in the opposite
order.

**Root cause: the driver's element list, the same idiom.** In `genmmd`'s main loop:

```c
mmdelm(mn, ...);
num += qsize[mn];
list[mn] = ehead; ehead = mn;          /* push each pivot at the head */
...
mmdupd(ehead, ...);                    /* then walk from the head */
```

So the LAST pivot taken in a batch is the FIRST element refreshed. We append to `batch` and walk
forward.

**Change: walk `batch` backwards.**

**Result: 4 of 7.** Reversing `qxh` as well, which is the same idiom in `mmdupd` and correct for
consistency, changed nothing on these seven, but was kept: leaving one of three wrong knowingly is
not a defensible state.

**Established: three instances of one defect. Still stuck at 4 of 7.**

---

## Iteration 3: THE FIRST WRONG TURN, and what it revealed

**The situation.** Three walks reversed, 4 of 7, and graph2, graph3 and graph4 all still diverged.

**What was tried, and it looked like progress.** Sorting `clique_members` ascending reached **5 of
7**. Sorting descending gave 2 of 7.

**Why it was REJECTED.** genmmd never sorts. It walks `adjncy` in whatever order `mmdelm`'s
compaction left it, which comes out ascending on these graphs only because the input is built from
sorted columns and compaction preserves relative order. Sorting reproduces the EFFECT on this test
set without reproducing the MECHANISM, and costs an `O(d log d)` sort per element in a routine
whose entire design rests on not needing one.

**The rule this produced, and it is in the README's method section.** The question to ask of any candidate fix is
**which line of the vendored routine does this correspond to?** If the answer is "none, but it
scores better", it is wrong. A better score is not a correct port.

**What the failed attempt revealed, which was the actual value.** Instrumenting the q2h list on
graph2 showed:

```
element 1   members=[0, 3]   reversed walk=[3, 0]   ->  3 survives
element 4   members=[3, 0]   reversed walk=[0, 3]   ->  0 survives
```

**The same pair `{0, 3}` stored in OPPOSITE orders in two different cliques.** So reversing the
walk gives a different survivor depending on how that clique happened to get built. The CONTENT
order of the lists was wrong upstream, not merely the direction they were read in. That is what
sorting was papering over, and it pointed straight at how `C[pivot]` gets built.

---

## Iteration 4: the deep entry, and 7 of 7

**Root cause.** `mmdelm` builds the pivot's new list in two parts:

```c
for (i = is; i <= it; i++) {
    nb = adjncy[i];
    if (marker[nb] < tag) { marker[nb] = tag;
        if (fwd[nb] < 0) { list[nb] = el; el = nb; }   /* an ELEMENT: push on a stack */
        else             { adjncy[rl] = nb; rl++; } } }/* a VARIABLE: write forward */
while (el > 0) { ... }                                 /* then drain the stack */
```

Explicit neighbors written forward, then the ELEMENTS expanded in REVERSE order of encounter,
because `el` is a stack. Our `neighbors` walked `for c in I[u]` forward.

**Change: walk `I[u]` backwards in `mmd3_neighbors`.**

**Result: 7 of 7 examples, in one step, from 4 of 7.**

**Why this one is the deep entry and the other three could not be judged without it.** It fixes the
order of `C[pivot]`, which is the CONTENT order of every list built downstream, including the q2h
and qxh lists that entries 2 and 3 walk. With entries 2, 3 and 4 alone the examples stall at 4 of
7; entry 1 alone would not be right either. They had to be closed together, and that is the one
place where the one-change-at-a-time discipline was correctly relaxed: not to batch fixes, but
because the four are one defect with four instances.

**Established: all four entries are conventions.** genmmd pushes at the head of a linked list
threaded through an integer array and reads from the head, so the entry seen LAST is processed
FIRST. We hold a vector and append. Same set, opposite order, same cost. Only the winner among
equals changes, and minimum degree is settled by exactly that.

**Grids at this point: 67 to 78 percent of pivots matching**, a fraction FLAT in `n` rather than
falling, which says one more mechanism was unaccounted for rather than drift accumulating.

**What it bought already.** Fill against genmmd roughly halved, from about 20 percent to about 10,
with the time gap unchanged. That was the first evidence that the gap was not all mechanism.

---

## Iteration 5: THE SECOND WRONG TURN, and the first real defect

**The situation.** Seven of seven on the examples, grids diverging in the last third.

**Smallest case: the 5x5 grid, `n = 25`, diverging at pivot 19 of 25.** Six pivots from the end,
small enough to read whole.

**What our trace showed, and it read like a tie-break.** At the divergent iteration the bucket held
`6: [13 11 1]`, three vertices of EQUAL degree, and we took the head, 13. genmmd took 1. That is
indistinguishable from the four entries already closed.

**It was not a tie-break, and only the ORACLE could show that.** Two `fprintf` calls added to a
scratch copy of `Mmd.cpp`, one at each pivot and one at each merge:

```
PIVOT 9  mdeg 6  qsize 1
PIVOT 5  mdeg 6  qsize 1
  PAIR 21 into 1
PIVOT 1  mdeg 5  qsize 2
```

**genmmd files vertex 1 in bucket 5. We had it in bucket 6.** Its true degree is the same:
genmmd files at `dg - qsize + 1`, and `5 = 6 - 2 + 1`. A supervariable is filed LOWER by its own
size minus one, so it is picked earlier than a singleton of equal degree.

**And we performed the same merge.** Instrumenting our side printed `OURS PAIR 21 into 1`. So the
merge was not missing. What differed was WHEN the weight is subtracted:

```python
# ours
degree = dg0 - len(super_members[u])      # snapshot BEFORE the walk
...  the walk MERGES 21 into u, so super_members[u] grows to 2  ...
degrees[u] = max(degree + 1, 1)

# genmmd
dg = dg0                                   # kept whole
...  qsize[en] += qsize[nd] during the same walk  ...
dg = dg - qsize[en] + 1                    # subtracted at the END, POST-merge
```

**Root cause.** The two agree until the walk merges a vertex into `u`. When it does, we subtract the
pre-merge weight and genmmd subtracts the post-merge one, so **a supervariable is filed one bucket
too high per vertex merged into it** and is never picked as early as its size has earned.

**Change: keep `dg0` whole and subtract at the end.** One line moved.

**Result.**

```
5x5 and 7x7 grids    MATCH outright
grid divergence      10x10 from pivot 68 to 89, 12x12 from 111 to 133
fill against genmmd  ZERO at every size on the ladder, 32 through 400
```

**This is a DEFECT, not a convention.** The code did not do what its own comment said,
`dg - qsize[en] + 1`. It produces a worse ordering on its own terms, with no appeal to genmmd
needed.

**A defect found in one place is a defect everywhere the code sits.** Grepping for the pattern found
it in exactly three files, and all three were fixed: `mmd2.py`, `mmd2.cpp` and `src/Mmd2.cpp`. It
had been costing fill in production since `Mmd2` was written:

```
grid          MMD2 fill before   after
32x32              +1.6%         -0.5%
100x100           +16.2%         +5.9%
200x200           +18.4%         +6.8%
400x400           +25.2%         +8.3%
```

**MMD2's own gap fell from about 20 percent to about 7, and at 32x32 it now beats genmmd.**

**Where it could NOT be.** `mmd1` has no q2h path and no live merges.

**And this paragraph used to rule the amd branch out too, WRONGLY. Corrected 2026-08-08.** It read
that the amd layers file at an external degree, which excludes a vertex's own supervariable and
therefore does not move when its weight changes, and cited `amd2`'s own comment at its hash merge
as saying the same. The external degree does not move; the `- weight(u)` term inside the amd bound
does, and it is the term that decides the bucket. `amd2`, `Amd2` and `Amd2B` carried exactly this
defect, worth 3 to 9 percent of fill on grids, and the comment cited in support of the claim was
making the same mistake. It was found only when the amd branch got this same treatment, where it is
`experiments/ordering/AMD3.md`'s entry 4, and fixed the same day.

**That makes this the most expensive sentence in the document, and worth reading twice.** It is the
one place mmd's ledger reached beyond its own branch, and it reached with an argument rather than a
check. Had the amd bound been read at the time instead of reasoned about, the defect would have
been found here rather than months later. A defect ruled out by an argument is ruled out only as
far as the argument reaches.

**The lesson, and it is in the README's method section.** The symptom does not identify the cause. A different pivot at
equal apparent degree is what a tie-break AND a filing defect both look like. Only the oracle's
trace separated them, and reasoning from our own trace had produced a plausible and wrong diagnosis
twice by that point.

---

## Iteration 6: comparing the right object, and the last entry

**The situation, which looked like a wall.** Fill was exact at every size and the permutation still
diverged, at pivot 700 of 1024 at 32 a side. A fill of exactly `4862612` against `4862612` at
400x400, from a different permutation, at seven sizes. That is not coincidence.

**The right question, and it is the most transferable thing in this document.** Three comparisons of
increasing strictness:

```
fill               nnz(L). Coarsest. Equal fill does NOT mean equal ordering.
pivot sequence     the pivots chosen, in order. This is the ALGORITHM.
permutation        the expanded order, including supervariable members.
```

**Comparing PIVOT SEQUENCES instead of permutations answered it immediately:**

```
grid   6   pivots: genmmd  31, mmd3  31   IDENTICAL
grid   8   pivots: genmmd  54, mmd3  54   IDENTICAL
grid  10   pivots: genmmd  81, mmd3  81   IDENTICAL
grid  20   pivots: genmmd 314, mmd3 314   IDENTICAL
grid  32   pivots: genmmd 788, mmd3 788   IDENTICAL
```

**The algorithm was already aligned.** 788 of 788 pivots at 32x32. Comparing permutations would
have sent the next hour hunting a mechanism that was not missing.

**What actually differed, from the 6x6 tail:**

```
vendored   ... 20, 12, 13, 15, 17, 22
mmd3       ... 20, 13, 12, 17, 22, 15
```

The same six vertices. The oracle's trace showed them to be ONE supervariable: `PAIR 13/12/22/15
into 20`, then `PIVOT 20 mdeg 1 qsize 6`. What differs is the order its members are listed when it
expands.

**Root cause: `mmdnum`.** It scans vertices in ASCENDING index order, and each merged member takes
the next number after its root:

```c
for (nd = 1; nd <= neqns; nd++)
    if (perm[nd] <= 0) { fa = nd; while (perm[fa] <= 0) fa = -perm[fa];
                         rt = fa; nm = perm[rt] + 1; invp[nd] = -nm; perm[rt] = nm; ... }
```

Root first, then members by ascending vertex index, and **no sort**: ascending falls out of the scan
order. We listed them in merge order.

**Change: one ascending pass, no sort.** Build a root map, then scan `v = 0..n-1` ascending
appending each member to its root's block.

**Result: EXACT permutation match.** Seven of seven examples, and every grid tested from 5 a side to
80, `n = 6400`.

**Decided to close it even though it cannot change the answer.** Supervariable members are
indistinguishable by construction, so fill and the elimination forest are identical either way.
Entry 6 buys no quality at all. It buys an **equality test** instead of a judgement, which is the
instrument the next layer gets aligned with. That was judged worth one function.

---

## Iteration 7: production, and one self-inflicted wound

**Three of the four reversed walks are in `Mmd3.cpp`.** The fourth, the `I[u]` expansion, lives in
`QuotientGraphFlat::reachableSet`, which all six drivers share.

**Decided: a flag, `setReverseIncidence`, off by default and turned on only by `Mmd3`.** The branch
is hoisted like `live` and is per CLIQUE, not per member, so it sits outside the loop that does the
work. The alternatives were duplicating the walk or changing the permutation for every driver.
A mode flag on a shared class is not free to the reader, and that is recorded where it lives.

**Entry 6 needed the same treatment**, since `QuotientGraphFlat::order` is shared too. It became a
second named method, `orderAscending`, rather than a second flag: an output convention reads better
as a named method than as a boolean at a call site.

**And the self-inflicted wound.** `orderAscending` as first written allocated FOUR arrays of size
`n` and made six passes. In the Instruments trace it cost **244 ms of a 4.94 s profile**, where
genmmd's `mmdint` and `mmdnum` TOGETHER cost 116 ms, while doing strictly less work than either.

Rewritten to one scratch array and two passes, with a root's cursor and a member's root marker
sharing that array by SIGN, it is **60 ms**, level with `mmdnum`. That single function was 14
percent of the whole time gap over vendored MMD.

A second attempt on the same file was measured and REJECTED: replacing the constructor's `reserve`
plus `push_back` with a sized write and an index is SLOWER, because `resize` value-initializes the
whole arena and it is then overwritten. The capacity check being removed is cheaper than the
zero-fill being added.

---

## What the alignment then made possible

**With the permutation identical, the remaining difference is implementation and nothing else.**
That is what turned the time gap from an argument into a measurement. Instrumenting genmmd with the
same counters our prototypes print:

```
                 mmd3      genmmd
degree updates   equal      equal      1575, 5515, 12644, 23861 at sides 32, 64, 100, 140
pair merges      equal      equal
outmatched       equal      equal
ncsub            equal      equal
bucket probes      289        193      the only one that differs, O(n) either way
```

No algorithmic difference left. The time gap was then decomposed, and the largest single
contributor found by WIDENING the oracle rather than narrowing ourselves: genmmd retyped from `int`
to `int64_t`, doing byte-for-byte identical work, costs 17 to 26 percent. See
`experiments/ordering/REPORT.md`.

**Final position.**

```
                 fill vs genmmd    time vs genmmd (alpamayo)
MMD1                +1 to +23%       1.6x to 2.8x
MMD2, fixed         -0.5 to +8%      1.06x to 1.47x
MMD3                 0.0% exactly    1.03x to 1.49x
```

MMD3 became the default ordering the same day, on the argument that reproducing a reference with
decades of use behind it is a better bet on unseen inputs than a tie-break of our own tested on
grids. Recorded in `docs/DESIGN_DECISIONS.md`.

---

## The ledger, for reference

Append only. A row is never edited once closed.

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

Six entries, seven working iterations, two wrong turns, one real defect that had been costing
production fill for months, and one performance defect introduced during the work and caught by
profiling the same day.
