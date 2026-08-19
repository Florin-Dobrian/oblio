# Decisions

Durable record of structural choices, newest first. Each entry: date, decision,
why. This is the file to open after a gap to reconstruct the project's shape.

---

## Why this file exists

**In a solver, the bottleneck was never writing the code.** A Cholesky is a day's work. A
supernodal Cholesky is a week's. What takes years is deciding *how the pieces should fit*: where
the factor's values live, whether the engines hold state or the data does, how symmetry is
expressed, what the API promises, which phase owns which fact. Those are the decisions Oblio was
made of, and they were the hard part in the late 1990s.

**They are expensive because they propagate.** Choose the storage and you have chosen the API;
choose the API and you have chosen what every engine can see; and by the time a mistake shows up
it is behind ten thousand lines that assume it. A wrong turn is not a wrong function, it is a
redesign. So each decision carries risk, so you deliberate, so they queue, so the project slows.
That is the real reason a library like this takes a decade rather than a year.

**And the tooling of the time made it worse.** C++98 with shaky compiler support: templates you
could not trust, no move semantics, containers that might or might not be there, and no cheap way
to *try* an abstraction and measure what it cost. The design space was constrained and the
instruments for exploring it were bad. So you chose carefully, once, and lived with it.

**What has changed is not that decisions became easier. It is that exploring one became
disposable.**

A design question can now be turned into a measurement in an hour (`experiments/storage-options`
answered flat-versus-VV with a number and a symbol table). A ported kernel can be tested standalone
before it touches anything (0.9's LDL was reconstructing `A` before the engine had ever seen it). A
clever solution can be built, looked at, called a hack, and rewritten before lunch (the
identity-offsets array in `SymFactorEngine`, which existed for about twenty minutes). A claim can
be made, and superseded the same day when it turns out to have been true only of the case that
prompted it ("symmetry is determined, not chosen", which LDL falsified within hours).

**So the thing to optimize for is not choosing correctly the first time. It is making a wrong
choice cheap to discover and cheap to reverse.**

That reframes what infrastructure is *for*:

- **Per-phase oracles.** Each phase is checkable against something that shares no code with it, so
  a mistake is localized rather than inferred.
- **An end-to-end residual.** `||Ax - b|| / ||b||` in one number, across six phases, tells you the
  pipeline is *consistent*, which no per-phase test can.
- **`experiments/`.** A design argument that can be settled by measurement should be, and the
  measurement should be kept, so nobody re-litigates it from memory.
- **This file.** A decision that was reversed is more instructive than one that was right, and a
  reversal that is not written down will simply be made again. Entries here supersede rather than
  overwrite, and say what was wrong and why it looked right.
- **A willingness to throw work away.** The cost of a wrong turn is now measured in minutes, which
  means the correct response to "that is a hack" is to do it again, not to defend it.

**And the human's role changes rather than shrinks.** Code is fast; *judgment* is not. Which claim
is load-bearing and which is decoration, whether an argument is sound or merely plausible, whether
an abstraction is elegant or merely clever, these do not get cheaper, and they are what decides
whether the fast part was worth doing. Nearly every real correction in this project came from
someone looking at a plausible answer and saying "no, that is not quite it".

The best example is the smallest. Cholesky spent years as an open question here, "in real it is
`CC^T`; in complex, can I have both `CC^T` and `CC^H`?", and it would not close, because it was
the wrong question. `CC^T` in complex does not exist: positive definiteness requires `x* A x` to
be real, which requires Hermitian. Once that is on the table, the design collapses to a sentence:
**Cholesky is `CC^H`, always, and in real that *is* `CC^T`.** No option, no flag, no forbidden
combination to reject. The answer was not hard. Asking the right question was.

## 2026-08-19: a shared class splits where the vendored codes disagree, and nowhere else

**The decision.** `QuotientGraphCompacted` will be one public class serving both `Amd3B` and
`Mmd3C`. Where the two vendored codes disagree, the class carries TWO METHODS with a branch suffix,
`pruneAmd` and `pruneMmd`, rather than one method with a flag. Where they agree, one
implementation. The count of splits is not a constraint; contorting shared code to avoid a suffix
is worse than having one.

**Why suffixes and not flags.** The class had flags, `mVendoredListOrder` and `mReverseIncidence`,
and `Amd3B` removed them in favour of one hardcoded layout. Flags put both vendored conventions in
one body, which is where they are hardest to check against the sources they came from; two bodies
can be read against `AMD_2` and `mmdelm` directly. The runtime branch is the smaller argument, being
constant for a whole run and perfectly predicted.

**What the split is FOR, and it is narrow.** Two conventions, both about ORDER and neither about
storage: where the new clique is written into I[u], `AMD_2` first and genmmd last; and which of a
vertex's two halves `reachableSet` walks first, with the incidence direction. Both decide the order
members enter C[pivot], hence bucket order, hence which of several equal-degree vertices is picked.
Neither is free to move: each driver reproduces its reference's permutation exactly or it is not a
differential.

**And what is NOT the split.** The placement rule is common and stays shared: build the clique in
place when the incidence list is empty, otherwise at the free cursor, compact when the cursor runs
out. So is the collector, including the mid-walk truncation, which was the part that had two bugs.

**One general point worth keeping.** It was tempting to look for a single layout serving both, and
the search was settled by TESTING rather than by argument: four walk variants over `AMD_2`'s
insertion, none reproducing genmmd past a 5-squared grid. The reasoning that preceded it had been
wrong twice in both directions on the same day. Where a claim about equivalence can be run, run it.

**AND A CAUTION ON HOW THIS DECISION WAS REACHED, recorded because it may yet reverse it.** The
split was derived by diffing `Amd3B` against `Mmd3C`, two copies taken into private namespaces and
then evolved apart for months. The SHARED class already serves `Amd3` and `Mmd3` on one layout with
two flags, so two vendored conventions provably coexist in one body. How much of the catalogued
divergence is the references disagreeing and how much is drift between two unreconciled copies is
NOT settled. A private copy is a fine instrument and a poor baseline: what it shows is where one
line of work went, not what the problem requires. See docs/NEXT.md (2026-08-19) for the one
experiment that separates the two.

## 2026-08-19: a derived counter is only as good as the sites it is funnelled through

**The claim being made when `numPeakCliqueMembers` went in** was that funnelling every birth,
contraction and death through `killClique` and `trimClique` made the counter correct by
construction. That was wrong twice within a day, both times in `Amd3B`, and the shape of both
mistakes is the same: a site that changes a clique's length WITHOUT going through either funnel.

- `beginElimination` re-read `I[pivot]` after the walk to find the cliques to kill, and the
  mid-walk collector truncates that list as it consumes it. Cliques lost that way were never killed
  at all.
- The truncation itself shortens the clique being consumed. That is a contraction, and nothing told
  the counter.

**The durable point.** A funnel is a claim about the CALL GRAPH, and it holds only while every
writer of the underlying state goes through it. Neither of these did, and neither was visible by
reading, because both sites look like storage bookkeeping rather than like liveness. The way to
make such a claim checkable is a RECOMPUTATION from independent state, which is what
`QuotientGraph::cliqueCountBalances` does and what `Amd3B` still lacks.

**And on what found them: not a test.** The digest, the vendored acceptance checks, `test_order`
and the sanitizers were all green throughout. What caught them was a cross-driver comparison in a
BENCHMARK, `Amd3B`'s peak against `Amd3`'s, on matrices from a real collection. A quantity that is
a property of the algorithm rather than of the implementation makes two implementations comparable
in a way that output equality does not: the permutations agreed the whole time.

**The general rule this supports**, and it is now twice-earned in two days: an oracle that asks
whether an answer is LEGAL cannot see a defect that changes only how the answer was reached. The
2026-08-18 stamp entry is the other instance, and it cost fill rather than work.

## 2026-08-18: a tag array holding two kinds of value needs the ranges kept apart

**The bug.** `Amd3`'s `w` array carries the scan's `w[c] = degree[c] + wflg - nvi`, which reaches
`wflg + lemax`, AND supervariable detection's stamps. A stamp has to sit above every scan value of
the same step, or a clique whose scan value happens to equal the current stamp reads as marked. We
raised the stamp base at the END of the step, which left this step's stamps starting at `wflg` while
this step's scan values ran up to `wflg + lemax`. The ranges overlapped exactly.

**What it did.** The exact test in supervariable detection could return a FALSE MATCH, merging two
vertices that are not duplicates. On `Grund/meg4` it merged 5779 into 5780 although their lists
differ in six of sixteen entries: 109 positions of the permutation moved and fill went to 51809
against `AMD_2`'s 51512.

**The fix is one statement moved**, `stamp = std::max(stamp, wflg + lemax)` before the detection
loop rather than after it. `AMD_2` does exactly this with `wflg += lemax ; wflg = clear_flag (...)`
between scan 2 and SUPERVARIABLE DETECTION, and the reason is the same one.

**The durable lesson, which is why this is here and not only in NEXT.md.** When one array holds two
populations of value, the invariant separating them is a design obligation and belongs written at
the declaration, not inferred at each use. This tree has folded several arrays into one on purpose
and will fold more; every such fold creates exactly this kind of obligation. The 2026-08-11 entry
narrowing arrays and the fold work generally should be read with this attached.

**And on what could have caught it: nothing we run.** The digest over 73 grids, both scaling
ladders, the 38 vendored acceptance cases and the sanitizers were green throughout. The overlap
needs a clique degree landing on one specific value; it fired on ONE matrix in 246, and only
because `benchmarks/matrices` had just been given a check comparing our permutation against the
vendored one entry by entry. A defect that changes fill without changing validity is invisible to
every oracle that asks whether an answer is legal rather than whether it is the same.

## 2026-08-18: flat and chunked, and the four quotient graph classes

**The words.** A SEGMENT is the logical unit of storage: one vertex's run of `A[u]` and `I[u]`, or
one clique's members. It exists identically under every scheme here. What differs is the
ALLOCATION holding it, and that is what the pair of names is for:

- **FLAT**, one array holds every segment end to end, each segment carved out of it.
- **CHUNKED**, each segment is its own allocation.

Both grow; only chunked can hand memory back. The distinction is the one NumFactor draws today with
`Static` and `Dynamic`, and those words are wrong for it: NumFactor's flat case really is static,
the symbolic factor fixing the size before anything is written, while NEITHER clique store has that
property. Both evolve. So flat and chunked will eventually replace static and dynamic there too.
Not yet, and not in this entry.

**Two words considered and declined.** `segmented` as the opposite of flat, because a segment is
what BOTH schemes hold, so the name would describe the thing they share; "the segments in the flat
space" beside "the segmented space" is a sentence a reader stumbles on. And `blocked`, because a
block is already the dense submatrix in NumFactor and `UpdateBlock` is a class.

**Saying flat or chunked implies C IS SEPARATE from A and I.** When it is not said, the store is
unified and always flat, since a unified store is one array by construction.

### The four classes

```
QuotientGraphFlat        C separate, flat        today's, shared by six drivers
QuotientGraphChunked     C separate, chunked     not built
QuotientGraphCompacted   C unified, compacting   private copies in Amd3B.cpp and Mmd3C.cpp
QuotientGraphChained     C unified, chaining     private copy in Mmd3B.cpp
```

**All four are meant to ship**, as NumFactor's two do, so eventually none is unqualified and
today's `QuotientGraph` becomes `QuotientGraphFlat`. That rename waits until the chunked one is
actually being added, since it touches six drivers and a great deal of prose and is worth doing
once.

**The private copies stop being private and become these classes**, in this order, one at a time:
documentation first, then `QuotientGraphChained`, then `QuotientGraphCompacted`.

**Chained and compacted are separate classes rather than one with a policy flag.** They share a
storage shape, one array sized at nnz, and differ in what happens when a clique will not fit. But
chaining also changes the READ path: every walk of a clique pays a link test, forever, where a
compacted walk is a straight run. A runtime flag would put that branch in the hottest loop in the
ordering.

### What blocks `QuotientGraphCompacted`, and it is a decision rather than work

`Amd3B` and `Mmd3C` both hold a copy called `QuotientGraphA` and THEY HAVE DIVERGED. `Amd3B`'s
carries the run flip, the rotation, positions and the mid-walk collector; `Mmd3C`'s carries none of
it. One class cannot serve both as they stand: the flipped run leaves no free slot for `Mmd3C`'s
back-of-list convention, which is exactly why the flip stopped at `Amd3B`. Three ways out, none
taken yet:

- `Mmd3C` adopts the rotation, and its permutation moves, so its fill column stops being comparable
  to `MMD3`'s and its digest baseline has to be re-recorded.
- The class carries both conventions, costing a shift or a spare word per vertex on the mmd path
  and putting back a branch the flip removed.
- Promote the amd side alone, with `Mmd3C` joining once the convention is decided.

## 2026-08-18: the layout program, two axes and what Oblio ships on each

**The decision.** Oblio's shipping ordering uses the FREE MEMORY layout with no reclamation
machinery at all. The constrained layouts, compaction and chaining, stay private for now and are
built as alternatives to be tightened one at a time.

**The frame this sits in, and it is two axes rather than one.** How much space beyond the pattern
are we willing to spend, and how much bookkeeping are we willing to run given that space. They
trade against each other, and the ends are not symmetric:

- **Free memory, no machinery.** Append every clique and reclaim nothing. Simplest to write,
  fastest per operation, and the one Oblio ships. Its cost is CUMULATIVE rather than peak: nothing
  is reclaimed, so the arena grows by every clique ever built and not by the largest set alive at
  once. Measured at 0.94x nnz(A) in 2D and 1.45x in 3D.
- **Bounded headroom, compaction.** The collector is IMPLEMENTED AND RARELY RUN. At `AMD_2`'s
  twenty per cent it fires about once for a whole ordering, and it degrades to zero elbow room
  rather than failing: `src/Amd3B.cpp` sized to exactly `sum(Len)` completes with 3 to 14 sweeps
  from 3 to 200 a side and the same permutation. There is no guarantee it will not be needed, so it
  has to be ready; there is every expectation it will not be needed often. That asymmetry is what
  makes it worth writing.
- **Zero extra memory, chaining.** Considered ONLY here. Its cost is a link test on every read,
  paid whether the space is tight or not, so headroom buys it nothing, where headroom makes a sweep
  rarer. Both ladders currently rank the chained layout slowest.

**Why a constrained mode is worth having when the free one is simpler.** Prediction. A run under a
memory bound either completes or reports what it would have needed, which is information about
whether the factorization is feasible at all. The free version can only fail later and larger: if
the arena does not fit then the factor would not have fitted either, so failing now or failing then
is the same failure carrying less information.

**Algorithms are never forbidden by layouts, only priced.** Any minimum-degree variant runs on any
of these layouts; enough indirection always expresses it. What a layout removes is a cheap
IMPLEMENTATION of some step, and the algorithm is then reached by a more expensive route. So the
coupling that remains is a price, and the engineering question is whether it is worth paying rather
than whether it is possible. The worked case is in experiments/ordering/README.md: the
incidence-first run leaves no free slot for appending a clique at the back of `I[u]`, which a
shift, a copy aside, or one spare word per vertex would each supply, and which `AMD_2` instead
sidesteps by rotating the clique to the FRONT for nothing, changing the convention and the tie
breaking with it.

**What this does not license.** Building every cell of the cross product. The layouts are three and
the algorithms are two, and the point of the private layers is to price the axes, not to fill a
matrix. Chaining with headroom is the cell explicitly not worth building.

## 2026-08-17: a one dimensional size may go signed, but only to buy an encoding

The 2026-08-11 entry narrowed four `QuotientGraph` arrays to `std::uint32_t` and named `mWeight`
among them, on the reasoning that a size has nothing to stand in for and so spends no sign bit.
`Amd3B` then found that `mWeight` does have something to stand in for, and the rule bends to admit
it rather than being weakened.

### The reasoning, which is the same reasoning

The rule was written to DERIVE rather than assert: unsigned because there is no sentinel, signed
because there is one. `AMD_2`'s `Nv` gives the weight a sentinel and a flag at once, so the premise
changes for that field and the conclusion follows. **Nothing about the rule is being excepted.**
What has to be stated is when the premise is allowed to change, since otherwise "it would be handy
to have a flag here" reopens every field in the tree.

Four conditions, in `docs/CODING_RULES.md`: per FIELD and never per category; the encoding written
out at the declaration; justified by MEASUREMENT rather than by tidiness; and no range lost, checked
rather than assumed.

### What the encoding is, and that it is two things

```
mWeight[v] >  0    live, not yet taken into the clique being built; the weight
mWeight[v] <  0    live, taken into it this step; the weight is -mWeight[v]
mWeight[v] == 0    dead, by a hash merge or by mass elimination
```

**The SIGN and the ZERO are different in kind, and only the sign is general.** The sign is a spare
bit: it carries an orthogonal fact while the magnitude still carries the value, which is unlike
`NIL`, which destroys the value it replaces. The zero is a true sentinel and rests on no live
supervariable weighing less than one. That holds on the amd branch, whose only death sites are
`merge` and `massEliminate`, and NOT on the mmd branch, where `number()` leaves a prepass vertex
live at weight one and in every list that named it. So a tree-wide "signed weights" change is safe
and a tree-wide "zero means dead" change is not, and the two must not travel together on one
justification.

### No range is lost, and this does not contradict the half-range clause

A weight is bounded by n and n is capped at `MAX_IDX = 2^31 - 1`, so the sign bit was never
reachable. Every accumulation of weights is over DISJOINT sets and so also bounded by n, which is
the exception the 2026-08-11 entry already carries. `CODING_RULES.md` says the half range "is not
headroom to spend", and a reader will think that collides with this: it does not, because that
clause is about ARITHMETIC that can exceed n, where the safety would be inherited from a cap
enforced elsewhere. Spending a bit that no value reaches is a different act from relying on one.

### One value is left spare, and it is worth writing down

The negated range is `[-(2^31 - 1), -1]`, so **`INT32_MIN` is unreachable from either direction**
and is free as a fourth state.

**The obvious use for it was proposed here and is REFUTED, same day.** The obstacle to taking this
encoding to the mmd branch is that `number()` must say "numbered" without destroying a weight, and
since `number()` runs once per driver in the prepass, before any merge, every vertex it touches is
at weight one; a `numbered` state at `INT32_MIN` would then be excluded by the same `nv > 0` test
that already excludes dead and taken. That reasoning is sound and the conclusion is still wrong.
**`orderAscending` reads `mWeight[pivot]` for every pivot**, to size each supervariable's run in the
emitted permutation, and a prepass vertex is a pivot. Its weight is load-bearing to the very end of
the ordering and cannot be spent. Found by reading `orderAscending` before implementing, which is
the only reason it cost nothing.

So `Mmd3C` keeps a mark array, and the slot stays free for some state carried by a vertex whose
weight is never read, should one appear.

### What moves with this

`docs/CODING_RULES.md`, whose one dimensional bullet now carries the exception and its four
conditions, and whose index bullet now says `NIL` is the common sentinel rather than the only one,
`Buckets` having carried `UNFILED` and `OUTMATCHED` all along. `docs/NEXT.md` bucket 3 lists
`mWeight` among the arrays narrowed on 2026-08-11 and is now wrong about that one.

---

## 2026-08-17: we self-alias too, and separate allocations do not prevent it

`MMD3C` read 1.28x to 1.47x `MMD3` at exactly 200 a side across five builds, while sitting at 0.95x
to 1.10x at 199 and 201 and at every other size on the ladder. The two compute the same permutation
from the same algorithm, so nothing algorithmic could produce it. Padding every size-n allocation in
`Mmd3C` by ONE PAGE, and changing nothing else, gives byte-identical permutations and:

```
MMD3C / MMD3            199^2    200^2    201^2
as shipped               0.98     1.28     0.99
padded by 1024 int32     0.95     0.99     0.97
```

**So it is data placement**, established by intervention rather than by argument, and the mechanism
is the one found in `AMD_2` on 2026-08-16: size-n arrays landing in the same cache sets at a
particular n.

### What this refutes

The 2026-08-16 entry concluded "our separate allocations are why we never had this." **We do have
it.** A separate allocation large enough is page aligned and rounded up to whole pages, so a set of
same-sized arrays lands at addresses congruent modulo the page size just as surely as a carve does.
The carve makes it easier to reason about and does not cause it. That sentence has been marked in
place; this entry is what it points at.

### The probe that read as a negative result and was not

The first attempt padded by SIXTEEN INTS, one cache line, and moved nothing, and I read that as
evidence against data placement. It was not evidence of anything. At 200 a side these arrays are
160,000 bytes, which rounds to 40 pages; 160,064 rounds to the same 40 pages, so the allocator very
probably returned the same addresses and the intervention never happened. **A perturbation that
does not perturb is inconclusive, not negative**, and the way to tell the difference is to know what
granularity the thing under test works at. One cache line is right for a small array and useless for
one past the page threshold.

Two rounds were spent on other hypotheses because of it: heap history and position in the run, which
were excluded by timing one function under two names at two positions and finding the spike at both.
That exclusion is still sound and is recorded in `benchmarks/ordering/order_timing.cpp`.

### What is not yet known, and matters more than this file

`Mmd3C` is transitional and about to be replaced, at which point every address changes and this
particular point disappears. **The open question is whether the shared class has such points**, since
`QuotientGraph` allocates the same shaped set of size-n vectors and `Mmd3` and `Amd3` would collide
by the same mechanism at whatever n happens to align. The ladder found this one only because 200 is
a rung; a spike at a size nobody runs would be invisible, and would show up as an ordering that is
mysteriously slow on one customer's matrix.

Two follow-ups, cheap, and the second is the important one:

- **Bisect which arrays collide**, by padding one at a time. Whether it is `mWeight` against `mMark`
  or the arena against `mSource` says whether production is exposed at all.
- **Sweep n finely for `MMD3` alone**, every side from say 190 to 210, and look for the same shape.
  One loop, and it is the only thing that would tell us whether this is a property of one
  transitional file or of the library.

**THE SWEEP WAS RUN THE SAME DAY AND CAME BACK CLEAN.** `tools/sweep.cpp` walks consecutive sides
and flags any whose cost per vertex stands more than ten percent above BOTH neighbors. Over 190 to
210 on alpamayo, neither `MMD3` nor `AMD3` has a candidate, and side 200 is unremarkable in both:
mmd reads 62.5 ns/vertex against 60.7 and 62.0, amd reads 85.0 against 85.7 and 87.2, which is
BELOW both. So **the point belongs to `Mmd3C`'s allocations and not to the shared class's**, which
is the answer the entry above wanted. `Mmd3C` differs in the one way that touches addresses: it is
a single translation unit, so its vectors are constructed in a different order and land elsewhere.

Two limits on that. It is one band of 21 sides, so a point at some other n is UNOBSERVED rather
than excluded. And the amd series scatters 12 percent over the range against mmd's 6, close enough
to the threshold that a genuine amd candidate would need a second invocation to separate from noise.

The other follow-up, bisecting which arrays collide, was NOT run and is not worth running: it is
about a file that is being replaced, and replacing it changes every address.

### The lesson, which is a sharpening of the 2026-08-16 one

That entry said prefer intervention to inference, and this one agrees and adds a condition:
**an intervention has to be large enough to intervene.** Both probes here were interventions; only
one of them changed anything, and the failed one looked exactly like a result.

---

## 2026-08-16: the layout matrix, which is where the next few days of work go

Two algorithms, three clique layouts, so six cells. Four exist or are planned; the two that exist
today are both on the diagonal, which is the limitation this entry is about.

```
                 our arena      genmmd's layout    AMD_2's layout
mmd ladder       Mmd3           Mmd3B              Mmd3C   PLANNED
amd ladder       Amd3           Amd3C   PLANNED    Amd3B
```

Ours is the separate append-only arena in elimination order. genmmd's is dead segments chained in
place. `AMD_2`'s is one pooled workspace with a free cursor and a garbage collection. Both vendored
schemes hold the ordering inside `O(n + m)`; ours does not.

**WHAT THE SUFFIXES MEAN, stated 2026-08-17 because the matrix implies a rule that was never
written down and both C cells were about to be built without one.** **B is a driver on its OWN
branch's vendored layout; C is a driver on the OTHER branch's.** So `Mmd3B` is mmd on genmmd's
scheme and `Amd3B` is amd on `AMD_2`'s, both same-lineage and both permanent, which is why they
were built first and why each is the natural differential vehicle for its own branch. `Amd3C` is
amd on genmmd's scheme and `Mmd3C` is mmd on `AMD_2`'s, both cross-lineage, and neither has a
vendored routine to be aligned against.

Two things follow that are easy to get wrong. **The two C files are not built from the files they
are named after**: `Amd3C` copies `Mmd3B`'s storage and `Mmd3C` copies `Amd3B`'s, so each is a port
from the sibling on the same LAYOUT rather than on the same branch. And the eventual pairing for
any shared code runs down the columns, `Mmd3B` with `Amd3C` and `Mmd3C` with `Amd3B`, not along the
suffixes.

**The `Mmd3C` that exists today does not follow the rule and is transitional**, being mmd on the
PRODUCTION arena. It was built to work the amd array folds out on the mmd side without disturbing
the class six drivers run, which it did; it is to be replaced by the real cell. See `src/Mmd3C.cpp`.

### Why the diagonal is not enough

**Each branch is currently aligned only with ITS OWN vendored routine**, which is exactly what a
same-lineage differential needs and is why the two B layers were built. But it means every reading
is a diagonal one, and a diagonal cannot separate the two variables. **With four cells the layout
effect reads DOWN a column and the algorithm effect ACROSS a row.**

**And it asks a question nobody has asked: is genmmd's layout better for amd, or `AMD_2`'s better
for mmd?** The two schemes are different animals, not two spellings of the same idea, and there is
no reason to assume each vendor chose the better one for its own algorithm. Dead segments cost
chaining and never compact; a pooled workspace costs a garbage collection and reuses everything.
Which wins may well depend on how the algorithm walks, and the two algorithms walk differently:
mmd's refresh runs over clique members contiguously, amd's prune alternates between the incidence
list and the clique blocks.

### The reason it matters right now, and not only for coverage

**The shared-class port has to work under every layout.** Three of the five folds from the entry
below land in `QuotientGraph`: the weight sign as the membership mark, `eliminated()` off a zero
weight, and mass elimination merging before it compacts. That class serves every cell of the matrix.
Four cells exercise all three layouts against both algorithms; two exercise each layout once and
each algorithm once, which is the thinnest coverage that can still be called coverage.

### What each of the two planned cells costs, since they are not equal

**`Amd3C` is the cheap one**: `Amd3` on genmmd's layout. The amd side has already proved the
encodings, so it is a storage swap into a driver that is otherwise settled, and it answers whether
genmmd's layout behaves the same under a second algorithm.

**`Mmd3C` carries a real design problem** and it is the blocker for the mmd port rather than a
coverage exercise. Amd calls `reachableSet` ONCE PER PIVOT and restores the negated weights in a
pass it already makes. Mmd calls `reachableSize` and `reachableWeight` PER VERTEX in the refresh, so
the negation needs a per-call lifetime: each call must clean up after itself. The cheap form is to
un-negate while walking the result, those callers already traversing it to count, but it has to be
right at every call site or the state leaks into the next call.

**So the order is a real choice.** `Amd3C` first confirms a layout under a second algorithm before
the hard question is opened; `Mmd3C` first puts the mmd port's blocker on the table before anything
ships. Both are defensible.

### What protects six copies of a quotient graph from drifting

`make digest` hashes every driver's permutation over 73 grids in half a second and names which one
moved, and `tests/test_order.cpp` asserts each B layer against its original. A copy that stops
reproducing its original is caught immediately, which is what makes the maintenance cost of this
matrix tolerable. The alternative shape, one class with the layout behind a compile-time policy, is
recorded in the entry below; three copies that exist beat one abstraction that does not, and this
matrix makes that four.

---

## 2026-08-16: the three storage schemes are all KEPT, and why staying in A's space is a real goal

From a conversation with Alex Pothen, one of the authors of the minimum-degree complexity paper
(`experiments/ordering/README.md`, "What the literature proves about these algorithms"), and it
explains something this tree had been reading as mere frugality.

### The argument

**Given a machine, you know whether A fits. You do not know whether L fits.** nnz(L) depends on the
ordering, which is the thing being computed, so it cannot be bounded before the run. A method that
works within `O(n + m)` therefore has a property no amount of speed substitutes for: **if the input
fits, the answer is reachable.** A method that grows as it goes may simply fail to produce an answer
at all, on a matrix whose input fitted comfortably.

That is why genmmd chains its dead segments and why `AMD_2` carries a garbage collection. Neither is
frugality for its own sake, and neither is an artifact of the era they were written in. They are
machinery bought deliberately, in exchange for a bound that can be stated from the input.

### And that is only HALF of what the B layers are for

The framing above treats them as storage alternatives kept for the trade-off, which is true but is
the outcome rather than the working reason. **Their operational role is to align the storage with
the vendored routine so that a differential is clean.** Comparing our ordering against genmmd or
`AMD_2` is confounded by layout until the two hold their cliques the same way; a B layer removes
that confound, and then whatever still differs is either LAYOUT, whose price these files measure, or
an IMPROVEMENT, which is carried back into our own ladder.

**That is not a hypothetical.** With storage held equal on the amd side, the differential surfaced
five array folds that have nothing to do with layout, which are the subject of the entry below and
were ported to `Amd3` and to `Amd2` on 2026-08-17, and the set is closed; see `src/Amd3B.cpp`'s
header for where each landed and why the fifth could not. The 2n mark shape came out
of `Mmd3B` the same way. So a verdict on storage does not retire either file: storage was never the
only thing they were for, and the next differential will want the same alignment.

**The port target is every applicable layer, not just the 3s.** Three of the five folds live in the
shared class and reach every driver whether or not anyone aims them; the two driver-side ones are a
decision per layer.

### How big our pool actually gets

Measured, because "grows as needed" is not a size and the argument above is easy to read as "ours
stores nnz(L)". It does not.

`tril(A) = n + m` is A in its stored form. `A+I` is `mSource`, which is `2m`: each edge appears in
both endpoints' lists. `nnz(L)` includes the diagonal.

**`C` IN THE TABLES BELOW IS TAKEN OFF THE ELIMINATION FOREST AND OVERCOUNTS.** It sums update parts
over SUPERNODES, where the arena holds one clique per pivot SUPERVARIABLE and each clique's entries
are supervariable representatives rather than rows. Against the arena the ordering actually fills it
reads 0.88 to 0.82 in 2D and 0.81 to 0.60 on cubes.

`benchmarks/matrices` now reports the MEASURED figure instead, through an overload on `orderMmd3`
and `orderAmd3`; the tables here are left as the run that produced them, and should be read as an
upper bound.

```
2D, five point
  grid          n           m      tril(A)          A+I            C         nnz(L)   C/tril(A)  C/nnz(L)
 50^2        2500        4900         7400         9800        12810          35913      1.73     0.357
100^2       10000       19800        29800        39600        54842         206332      1.84     0.266
200^2       40000       79600       119600       159200       227813        1081911      1.90     0.211
400^2      160000      319200       479200       638400       934527        5663298      1.95     0.165
800^2      640000     1278400      1918400      2556800      3783674       27361926      1.97     0.138

3D, seven point
  grid          n           m      tril(A)          A+I            C         nnz(L)   C/tril(A)  C/nnz(L)
 10^3        1000        2700         3700         5400         8668          32190      2.34     0.269
 20^3        8000       22800        30800        45600        90275         842282      2.93     0.107
 40^3       64000      187200       251200       374400       920867       20614676      3.67     0.045
 64^3      262144      774144      1036288      1548288      4447800      184222155      4.29     0.024
 80^3      512000     1516800      2028800      3033600      9172720      535760269      4.52     0.017
```

**C is about 2x tril(A) in 2D, flat, and 2.3 to 4.5x on cubes, rising. Against L it is 0.14 and
0.017.** So the arena tracks the input, not the factor.

### Why the compression is that large, and it is the same idea as symbolic factorization

A clique holds the UPDATE part of a supervariable's column and nothing else; `SymFactor` holds front
and update both. And mass elimination means one clique per pivot SUPERVARIABLE rather than one per
column: without it each clique would be that column's off-diagonal pattern in L.

**The factor is the update-weighted supernode size, not the plain average.** Mean front size on an
800 square is 1.33 columns, but `nnz(L)` sums `f(f+1)/2 + f.u` and is dominated by the few large
supernodes near the root, which have both a large front and a large update; the arena takes one
update from each. The effective factor is 7.2 at 800 squared and 58 at 80 cubed, and it rises with
n, which is why `C/nnz(L)` falls down both ladders while `C/tril(A)` barely moves in 2D.

**AND THERE IS A SECOND COMPRESSION THE FOREST FIGURE CANNOT SEE**, which is why it overestimates. A
clique's ENTRIES are supervariable representatives, one each; `updateSize` counts rows, one per
vertex. So the arena is compressed on both axes, blocks and entries, where the forest sum captures
only the blocks. That is the whole of the 0.82 in 2D and 0.60 on cubes above, and it explains the
direction: the gap widens exactly where supervariables are larger.

### And this is exactly why the guarantee still matters

**Read the tables as "often close to nnz(A)", not "usually".** The compression is bought entirely by
mass elimination, so it is a property of the MATRIX rather than of the method. Where supervariables
barely form, `f` stays near 1, the weighted factor collapses toward 1, and the arena approaches
nnz(L). Grids are a friendly case and both ladders here are grids.

**We cannot tell in advance which matrix we have.** That is precisely the thing the vendored schemes
do not need to know: they are bounded from the input whatever the matrix does. So ours being 2.5x
nnz(A) on a grid is a good number and not a bound, and the two facts sit together rather than one
displacing the other.

### Two working styles, both valid, and the tree should hold both

**Theirs is for DISCOVERY**: a matrix nobody has ordered before, on a machine of known size, where
the question is whether an answer exists at all. Tight and predictable wins, and the extra work
buys the guarantee.

**Ours is for PRODUCTION**: a known shape on a known machine, solved repeatedly, where it has
already been established that everything fits. Then the machinery is paid for a guarantee already
in hand, and the simpler, faster arrangement is the right one.

Neither is a better engineering choice in the abstract. They answer different questions.

### So all three stay, and TWO OF THEM STOP BEING TEMPORARY

| | space | measured cost against ours |
|---|---|---|
| ours: separate arena, grow as needed | not predictable before the run | baseline, and the production default |
| genmmd's dead segments, `Mmd3B` | one pool, predictable from the input | 4 to 10 percent slower on mmd, at every size |
| `AMD_2`'s pool with GC, `Amd3B` | one pool, predictable from the input | reads 2.6 percent fewer, D1 misses 6 to 8 percent more; a wash. Two compactions at every size from 50 to 1600 a side, constant |

Ours stays the production scheme. `Mmd3B` and `Amd3B` are no longer experiments awaiting a verdict:
they are **maintained alternatives**, kept so the choice remains available with its trade-off
measured rather than argued. Their stop conditions are withdrawn; the language in both headers that
said they go when a question closes is obsolete and has been replaced.

Names may change later. What matters is the status.

### The cost of keeping them, stated plainly because it is real

Each carries a **private copy of `QuotientGraph`**, so every shared fold has to be written three
times. The five folds of 2026-08-16 would have been fifteen. That is precisely what the stop-condition
language existed to prevent, and it is being accepted on purpose.

**The alternative shape, recorded so it is not rediscovered as an annoyance:** one class with the
storage behind a compile-time policy, which costs a template parameter instead of two copies. Its
price is that the hot loops become templated and the file harder to read, in a class whose
readability has been load-bearing all along. Three copies that exist beat one abstraction that does
not, so this is not the choice for now; but if the copies start diverging silently rather than
deliberately, that is the signal to revisit.

**What protects against silent divergence** is `make digest`, which hashes every driver's
permutation including both B layers, and the pair assertions in `tests/test_order.cpp`. A B layer
that stops reproducing its original is caught in half a second, which is what makes carrying three
copies tolerable.

---

## 2026-08-16, later still: five folds of one shape, and the amd slope goes below the vendored one

The seven folds below removed arrays. These five remove ARRAY READS FROM CONDITIONALS, which turns
out to be a different and much larger thing. The distinction was Florin's: it was never many arrays
against one, it was consulting a second array in a test where the vendored code answers from a load
it is already making.

### The pattern, in one example

Per adjacency element, `AMD_2` reads one array:

```c
nvj = Nv [j] ;
if (nvj > 0) { deg += nvj ; Iw [pn++] = j ; hval += j ; }
```

`Nv[j]` carries three facts: positive is the weight, negative means already taken into the new
element, zero means absorbed. One load answers "is it dead", "is it inside the new clique" and
"what does it weigh". We read `mMark[v]` for the two tests and `mWeight[v]` for the value, and
`mMark` is an array `AMD_2` DOES NOT HAVE, consulted in the innermost loop of the ordering.

### The five, all in src/Amd3B.cpp

1. **The sign of the weight is the membership mark**, `AMD_2`'s `Nv [i] = -nvi`. Removes `mMark`
   from the reachable-set walk and from the prune.
2. **The restore rides in the bound pass** rather than being a pass of its own. The first version
   negated and restored in two extra traversals per pivot and lost on x86 for that reason alone;
   `AMD_2` pays neither, its restore being inside RESTORE DEGREE LISTS.
3. **`eliminated()` answers from a zero weight**, `AMD_2`'s `Nv [i] == 0`. Removes `mMark` from the
   hash detection loops. It rests on an argument that has to be re-derived per branch: no list can
   retain an eliminated PIVOT, since if w is in A[u] then u was in C[w] and A[u] was pruned when w
   went.
4. **Mass elimination merges before it compacts**, so the compaction reads the zero weight the
   merge had to write anyway instead of stamping a tag of its own.
5. **Supervariable detection stamps into `w`**, which is `AMD_2`'s `W [Iw [p]] = wflg` over the
   whole list, variables and elements alike. This RETIRES mMark outright, 2n int32. The stamps must
   interleave with the tag protocol rather than clobber it: they start above `wflg + lemax` and
   `wflg` ends the step past all of them, so next step they read as alive-and-unseen.

### What they are worth

Cachegrind at 400 a side, exact and machine independent:

```
                 instructions   data reads   D1 read misses
AMD3  / AMD             1.024        1.179            1.016
AMD3B / AMD             1.018        1.081            0.968
```

The read excess falls from 17.9 percent to 8.1, and we now take FEWER D1 misses than the vendored
routine. On alpamayo, `AMD3B / AMD3` runs 0.95 to 0.81 across the square ladder with the advantage
growing in n.

**And the slope goes under.** Fitted as `time ~ nnz(A)^alpha` over twelve square sizes, on the
unaligned series where the vendored routine does not alias against itself:

```
AMD 1.080     AMD3 1.088     AMD3B 1.056
```

So the excess growth that this file spent a day chasing, quoted first at 0.032 and then at 0.016,
is now NEGATIVE. `AMD3B` is 0.93x the vendored routine at 1024 squared and 0.97x at 1600, the first
square sizes where anything of ours has beaten it.

### The storage change, which was the point of the file and is not the answer

`Amd3B` also carries `AMD_2`'s clique storage: one pool with elbow room, a free cursor, absorbed
space left dead, and its garbage collection ported with its own FLIP trick so the compaction is one
linear scan. Two compactions at every size from 50 to 1600, constant, which is what the complexity
bound assumes.

**Measured alone it is a wash**: 2.6 percent fewer data reads, 6 to 8 percent more D1 misses, both
constant across the range. So the arena is not the growth and never was, which is what Florin said
before it was built. It should NOT port; it drags a garbage collector behind it for nothing.

### The instrument finding, which is the one with the widest reach

**The sandbox's wall clock is not alpamayo's for any change that trades one counter against
another.** The first version of fold 1, +1.4 percent instructions for -6.3 percent reads, measured
SLOWER on x86 at every size and 3 to 12 percent FASTER on alpamayo. Cachegrind's counts are
machine-independent and remain trustworthy; the sandbox's timings are not, and a fold judged there
would have been discarded.

**And the container layer is 17 percent of instructions, not the 1.5 percent the totals suggest.**
A per-function profile puts 65.8 M instructions in `stl_vector.h` inlined accessors, against 236 M
of algorithm; the vendored routine runs 337 M of algorithm and no accessors. Two effects were
cancelling: we execute substantially less algorithm and pay it back in the container layer. The
misses are not there, `stl_vector.h` holding 0.3 percent of D1 read misses against `eliminate`'s
42.5.

### Where this leaves the port

Three of the five are in the SHARED class, so mmd sees them: the weight sign, `eliminated()`, and
the mass-elimination reordering. The other two are `Amd3` driver code. The shared three are also
the interesting ones for mmd, because `Mmd3`'s degree refresh is the same walk in mmd's hottest
loop, and mmd's constant against genmmd is 1.04 to 1.20.

**The lifetime is what makes it harder there, not the encoding.** Amd calls `reachableSet` once per
pivot and restores in a pass it already makes; mmd calls `reachableSize` and `reachableWeight` PER
VERTEX in the refresh, so each call must clean up after itself. Two quotient graphs would be the
failure rather than the fallback: the amd folds were cheap precisely because the shared class had
already paid for the encodings.

---

## 2026-08-16, later: the vendored AMD aliases against itself, and our excess growth is smaller than any ratio said

A wider scaling ladder went in after the folds below: ten square sides and six cubic, each built as
TWO INTERLEAVED SERIES rather than one geometric run. Square is `32 64 128 256 512` against
`50 100 200 400 800`, cubic `16 32 64` against `10 20 40`. Each series quadruples n; the two are
offset by about 2.4x and 2x.

**Why interleave.** A power-of-two side makes every array length and every grid stride a power of
two, which is the family that hits cache set conflicts, so a trend measured on such a ladder alone
cannot be told apart from an addressing artifact that grows with n. A second series that quadruples
n identically while never aligning gives the discrimination for free. It paid for itself on the
first run.

### The finding

**`AMD_2` carves six arrays of exactly n ints out of one block, and at a power-of-two side they
alias each other.** `AMD_1` does this:

```c
    s = S ;
    Pe = s ;  s += n ;   Nv = s ;     s += n ;   Head = s ; s += n ;
    Elen = s ; s += n ;  Degree = s ; s += n ;   W = s ;    s += n ;
    Iw = s ;  s += iwlen ;
```

so the six live at offsets that are exact multiples of n. `n = m^2`, so a power-of-two side gives a
power-of-two n, and the six then land in the SAME CACHE SETS at every index while the main loop
touches several of them per vertex.

**Measured with cachegrind, one vendored `amd_order` per run, three sizes back to back:**

```
side       n      I/vertex   Drd/vertex   D1 read misses/vertex
400   160 000       2285.1       603.4          15.32
512*  262 144       2279.0       598.1          40.27
800   640 000       2277.1       594.1          17.04
```

Instructions per vertex are flat to a TENTH OF A PERCENT and data reads are flat too: the routine
does identically the same work at 512 squared as at 400 and 800. Only the misses move, 2.6x at the
aligned size.

**And confirmed by intervention rather than by inference.** Inserting sixteen ints of padding
between the six arrays, changing addresses and nothing else, gives byte-identical permutations and:

```
D1 read misses per vertex     400      512*      800
as shipped                  15.32    40.27    17.04
padded                      15.27    17.56    15.64
change                      -0.4%   -56.4%    -8.2%
```

**56 percent of the read misses at 512 squared are the routine colliding with itself**, and the
padding changes nothing at 400.

### What this corrects, and how the wrong account was reached

An earlier version of this entry said the artifact was stride aliasing in a single array indexed by
vertex id, a 2D grid's vertical neighbour sitting at stride m, and drew from it that our clique
arena's elimination order buys immunity. **Both halves were wrong**, and the first was refuted by
this tree's own data before any measurement: genmmd is also one array indexed by vertex id and does
NOT zigzag, because its arrays are separate allocations rather than a carve. The class of cause was
right and the attribution was not.

Two other accounts died on the way, both quickly and both by measurement. The hash modulus is `n`
and so a power of two at those sides, which should have collided more: counted, and `AMD_2`'s pairs
per pivot are 0.33 at every size in both series. And the phase table localised the artifact to
`core`, `AMD_2`'s main loop, with `valid`, `aat`, `build` and `post` flat to three digits across
both series, which disposed of the workspace-sizing idea since the allocation happens in `build`.

### What the corrected picture is

**Per vertex, and the aligned rows are now known to be an artifact of the baseline:**

```
side       32    50    64   100   128   200   256   400   512   800
AMD      58.6  56.0  63.5  59.0  72.0  64.0  75.2  64.0  96.7  77.1
AMD3     68.4  76.0  78.1  82.0  84.2  87.2  89.3  94.2 108.1 127.2
genmmd   48.8  52.0  51.3  55.0  55.5  58.0  62.9  57.1  61.5  66.2
MMD3     58.6  60.0  61.0  60.0  59.8  61.5  61.3  61.6  64.0  72.1
```

`AMD3`, genmmd and `MMD3` are all smooth; only `AMD` zigzags. So the honest baseline for growth is
`AMD`'s unaligned series.

**The container overhead, measured twice independently.** At 32 a side `AMD3` is 1.17x `AMD` per
vertex and `MMD3` is 1.20x genmmd. That is the `std::vector` layer: raw arrays live in the vendored
routines' registers for a whole run where ours are members behind accessors. The mmd entry below
puts it at "a few percent" in prose; two measurements now put it near 20.

**And the excess growth that is genuinely ours is SMALLER than any ratio column suggested.** Fitted
as `time ~ nnz(A)^alpha` on the ten square sizes, least squares on log-log, R2 at or above 0.998
throughout:

```
MMD  (genmmd)   1.039        AMD  powers of two   1.080     AMD3  powers of two   1.071
MMD3            1.018        AMD  other sides     1.049     AMD3  other sides     1.081
```

**`AMD3` has ONE slope and `AMD` has two.** Ours is 1.071 and 1.081, indistinguishable between the
series; theirs is 1.080 aligned and 1.049 unaligned. So the self-aliasing above costs the vendored
routine about **0.03 in the exponent**, and at power-of-two sides it grows at our rate.

**Our excess over the honest baseline is 1.081 against 1.049**, about 0.032 in the exponent, or
1.23x over the 625-fold range. Not the 1.35x an earlier reading gave, and not the 1.5x the raw ratio
implied.

**And `MMD3`'s slope is BELOW genmmd's**, 1.018 against 1.039: on that branch we grow more slowly
than the reference, and the container constant is being eroded rather than compounded. `MMD3` reads
1.20x genmmd at 32 a side and 1.09x at 800. So the growth term is amd-specific and small, and it is
the only one either branch has.

**These are slopes to compare, not complexity claims.** The published bounds are worst-case and
dense and do not bind on grids; see `experiments/ordering/README.md`, "What the literature proves
about these algorithms". The fits are wall-clock and so include the memory effects; an
instruction-count fit would separate algorithmic growth from growth in cost per instruction.

### What follows for us

**This said "our separate allocations are why we never had this". THAT IS REFUTED, 2026-08-17, and
the correction is below.** It remains a caution about `NEXT.md` item 2b, which asks whether one
allocation carved into arrays would beat separate vectors on the model of `AMD_2`'s `S`: this is
what that costs when the sizes cooperate, and grid benchmarks are exactly where they do. But the
immunity claimed for the separate form does not exist.

### Three lessons, and the middle one is the expensive one

**A ratio hides which side is moving.** This pattern survived a full differential and seven folds
while being read as a property of our code.

**An account that fits the column it was written for can still be refuted by a column beside it.**
The stride account was consistent with everything in the amd table and contradicted by the mmd table
on the same page. Check an explanation against every code in the same run before writing it down.

**Prefer intervention to inference.** Correlation put the cause in the right class and the wrong
place three times. One padded copy settled it in a single run, and cost less than any of the
arguments did.

---

## 2026-08-16: the amd branch, seven folds, and what actually made the 2D curve flat

The two entries below did this to the mmd branch. This is the same work on the amd branch, and it
found something the mmd branch never showed, because the mmd branch never had the symptom.

### The symptom

`AMD3` sat at 1.25x the vendored routine on a 32 square and **1.82x on a 400 square**, rising
monotonically down the ladder, while `MMD3` over the SAME quotient graph was flat at 1.03 to 1.17x
and `AMD3`'s own cubic column was nearly flat at 1.02 to 1.34x. Something on the amd path cost more
per unit of work as n grew, in one problem family only.

### What it was not, and how that was established

A differential was built: counting copies of `Amd3B` and of the vendored `AMD_2`, generated by
anchor-asserting scripts in `tmp/`, counting the same passes under the same names and normalising
per pivot. Both codes produce the same permutation, so the pivot counts agree exactly and the
comparison is honest.

**Every pass matched.** Visits per pivot came out at 1.011 to 1.012 in 2D and 1.009 to 1.010 on
cubes, at every size: the reachable-set build, the prune's adjacency and incidence walks, the
bound, the hash chain, the pair loop, the exact test. `clear_flag` fired ZERO times on every case.
Our degree-bucket search was 6 to 100 times CHEAPER than `AMD_2`'s. Our clique arena was SMALLER
than its workspace in 2D, 0.77 to 0.83 of it, and `AMD_2` compacted on none of these runs.

So it was never work, never the arena, and never the hash. **Same visits, growing time, therefore
cost per visit.** A simulated L1 over the address streams put the excess at a CONSTANT 0.09 misses
per visit in both families, which is a real difference and far too small to explain a 1.8x ratio at
sizes where everything fits in L2, so single-level locality was not it either. Each of those was
a hypothesis of mine that the instrument killed, and two of them I had already argued for in
writing.

### What it was

**`mCliquePtr` and `mCliqueSize`, two n-arrays holding a clique's block position and length.** Every
clique visit in a walk probed both at a dead pivot's id, scattered across the whole id space, once
per element of every I[u]. `AMD_2` allocates neither: an element takes over the `Pe` and `Len` of
the variable that formed it. Folding ours into the dead pivot's own `mRun` entry turned two random
probes into one 16-byte load on a line the walk already touches, and

**the 2D ratio went from 1.11x rising to 1.49x, to FLAT at about 1.38x from 100 a side up.**

That is the finding. It also explains the family split without appealing to grid shape: a 2D pivot
has about 60 visits to amortise its per-pivot array cost over and a cubic one has 149, so the same
fixed cost is 2.5 times more visible in 2D. Same n, different ratio.

### The folds, in order, with what each measured

Measured on alpamayo against `AMD3f`, which is `Amd3` reached down the same harness path as
`Amd3B`, so the comparison is the fold alone. `AMD3f` against `AMD3` measured ZERO throughout,
which retires the "up to 2.4 percent" harness-seam figure in `Amd3.h`: it was a real reading on a
different driver on a different day and it does not reproduce.

| | fold | measured |
|---|---|---|
| 1 | the driver's `mark` into the graph's, `mMark` at 2n | 4.5% |
| 2 | the dead-clique test off the tagged `W` | nothing |
| 3 | the hash key, `AMD_2`'s arithmetic exactly | ~4% |
| 4 | the hash bucket scheme, `AMD_2`'s | the three together, 13% |
| 5 | `cliqueDegree` into `degrees`, which is `AMD_2`'s `Degree` | with 6 |
| 6 | the clique descriptor into the dead pivot's `mRun` | **the flatness** |
| 7 | `partial` into `w[u]`, the tagged `W`'s fourth question | 3.9% |

Seventeen entity-indexed streams down to eleven, against `AMD_2`'s nine.

### Three things worth carrying forward

**A fold can measure nothing and still be the reason another one works.** Fold 2 measured nothing.
The tagged `W` measured nothing when it went into `Amd2` on its own. Yet the tagged `W` is what made
the FUSED SCAN viable in both B layers, turning `AMD1B` from 7 to 9 percent slower than `AMD1` in 2D
into even-or-better, and `AMD2B` from 3 to 9 percent slower into even. Judging a fold only by its
own column would have discarded the precondition and kept the failure.

**The footprint trade is now confirmed three times and should be treated as a rule.** A schedule
change that saves visits but adds an array crossing a pass boundary tends to lose at large n. It
was recorded for the key array on 2026-08-08, for `Amd1B` at large n, and here for
`ApproximateScan`, which crossed `explicitPart`, `outside` and `mark`. With the tagged `W` the
crossing is one array and the same fusion becomes a win. Fold 7 is the general answer: the crossing
rides in a slot that is already there.

**Predictions were wrong in an instructive direction.** Fold 1 was predicted, in writing and before
the run, not to pay, on the grounds that the driver's `mark` is read only inside a hash collision
and is therefore cold. It paid 4.5 percent. The cost was never the number of touches. What that
suggests, and what is not yet confirmed, is that folds should be ranked by how COLD an array is
rather than how often it is read.

### What is left, and it is structural

`mSuperNext` and `mSuperLast` have no counterpart in `AMD_2` at all: it never materializes
supervariable members, which is exactly why `tools/hook_amd.py` has to reconstruct them to get its
elimination order out. Those two chains are what let us emit an expanded permutation directly.

`hashHead` cannot fold. `AMD_2` overlays its hash heads on `Head[]` with FLIP, parking a second head
in `Last[Head[hval]]` where both kinds are live at one index; our `Buckets` carries genmmd's `bwd`
encoding, in which a head's `mPrev` holds `-(degree + 1)` rather than being free. Two coherent
encodings that do not compose, and paying n int32 to keep the two list kinds apart is the cheaper
side of that trade. FLIP remains an anti-model here.

### And the question this leaves open

`AMD1` and `AMD2` still climb in 2D, 1.05 to 1.22x and 1.36 to 1.47x, where `AMD3` is now flat. The
shared descriptor fold reached them for free, so what remains is driver-side and lives in what those
two do that `Amd3` does not. That is where the differential goes next.

---

## 2026-08-15, later: what the container layer costs, and four attempts that did not reduce it

The entry below closed the ordering-speed question by folding arrays away. This is the second half
of that day, and it is mostly negative results, recorded because each one is a thing a later
reader would otherwise try.

### Where the branch finished

`MMD3` runs at 1.05 to 1.17x genmmd on square grids and 0.77 to 1.12x on cubes, and `MMD2` reads
the same to within the noise. That the two coincide is worth noticing rather than passing over:
they are different orderings with different mechanisms, and what they share is the quotient graph
and its clique arena, so the arena is what both are now bounded by.

**The trade the library makes, stated plainly.** Our code carries a `std::vector` layer genmmd
does not: raw `int` arrays arrive as its parameters and live in registers for a whole run, while
ours are members reached through accessors. That costs a few percent and it is not recoverable,
for the reasons below. What pays for it is the SECOND ARENA. genmmd keeps every clique in the dead
segment of the pivot that formed it, one array of nnz(A) and nothing more; we spend a separate
arena growing toward nnz(L). `Mmd3B` exists to price exactly that, and with every encoding fold now
present in both files so that storage is the only difference left, it reads 1.12 to 1.23x genmmd in
2D where `Mmd3` reads 1.05 to 1.17x, on 15.89M instructions against 14.22M. **Spending storage on
the arena more than covers what the containers cost**, which is the opposite of the premise the
storage investigation began from.

### What `Mmd3B` gained, and it was the same kind of thing as the morning

Four alignments, all of them genmmd's own shape, worth 4 to 11 percent on alpamayo with every sign
negative across seven sizes:

- **The stamping pass deleted.** The reach build already writes `mMark[v] = mTag` on every member,
  so `inClique` is that tag and the second walk of C[pivot] goes. `cliqueWeight` went with it,
  having no reader in that file.
- **Mass elimination compacts on `GONE`** rather than stamping a fresh tag over it.
- **`evict` in place of `unfile` then `restore`**, which is `mmdelm`'s single `bwd[rn] = 0`.
- **`mCliqueSize` retired for a value terminator**, `INT32_MIN`, since a link is `-(c+1)` and lies
  in `[-n, -1]`. That one needed the mark split to 2n first, vertices at `[v]` and cliques at
  `[size + c]`, because the dead-clique test had been the size.

`Mmd3B` is now genmmd's data structure essentially exactly, which is what makes it a fair
instrument rather than an approximation.

### Four things that did not work, and why each failed differently

**The stamping fold ported to `Mmd3`: +74000 instructions, +142000 reads.** The shared class cannot
delete the pass outright, `cliqueWeight` having five amd readers, so the accumulation moved into
`reachableSet`'s emit sites. That is the wrong loop: the emit runs over every candidate EXAMINED
and the stamp walk ran only over members EMITTED, which is fewer, since a vertex reached through
two sources is examined twice and emitted once.

**An arena cursor in place of `push_back`: +109000 instructions.** The capacity test can never fire,
`beginElimination` having reserved room for a whole reach, so it looks like free removal. It is not:
the vector's own length can then no longer be the arena's length, and `reserve` does not touch
memory while `resize` value-initializes, so the constructor zeroes nnz(A) entries and every growth
zeroes its new region.

**Raw bases in place of the accessors: +371000 instructions, and FLAT on the clock.** Hoisting
`mMark.data()` and `mWeight.data()` once per refresh round is what genmmd gets for free. GCC was
already keeping them where they belonged and two more live values cost more in the register
allocator than the reloads. Timed on alpamayo with `MMD1` and `MMD2` as controls: three sizes down,
three up, nothing. **The accessors are therefore both the cleaner and the cheaper form**, which is
a rare way for that argument to end.

**And q2h indexed instead of looped: a wash.** genmmd finds the single other source of a q2h vertex
by indexing, since the vertex has two sources and one is the new element; we ran two loops to find
an entry already known to be there. Removing them bought 103000 instructions of an 870000 pass and
nothing on the clock. Kept as the faithful shape, labelled a null. The comment defending the loops
was closer to right than the reasoning that overrode it: they were short, and the 2.5x in that pass
is not there.

### The one optimization that did work, and it came from a matrix with nothing to order

The 246-matrix timing run of the same day, `benchmarks/matrices`, `make mmdorder`, put the shape of
the remaining gap in plain view: where there is elimination work we run at 0.40 to 0.86x genmmd,
and where there is none we run at 2.5 to 2.8x. The worst five rows are matrices with NO OFF-DIAGONAL
ENTRIES, `nnz(A) = n` and `nnz(L) = n`. A higher per-vertex constant and a lower per-unit-of-work
cost, and on a pure diagonal the constant is the whole run.

**The prepass was one term of it.** It collected the degree-1 bucket into a vector and then walked
the vector, the list existing only because unfiling a vertex while walking the bucket destroys the
link the walk stands on. Reading the successor BEFORE the unfile removes the need, which is what
genmmd's numbering loop does. Worth **8.6 percent of a pure-diagonal ordering and 0.2 percent of a
grid**, and it moved those five rows to 2.0 to 2.3x. Applied to all three layers that have a
prepass, `Mmd2`, `Mmd3` and `Mmd3B`, measuring 8.8, 8.6 and 7.8 percent.

**What remains of the constant is CONSTRUCTION**, and this is the part worth carrying forward. With
no elimination work at all, an `Mmd3` ordering is roughly a third `QuotientGraph` construction and a
sixth `orderAscending`. Construction allocates and initializes about ten size-n arrays where genmmd
allocates five plus its 1-based copies. **That is the array-count finding of the same morning,
moved from the loops into the constructor**, and the reason nobody saw it is that a grid amortizes
it against real work while a diagonal has none to amortize against.

It is recorded as an experiment rather than an optimization, `docs/NEXT.md` item 8, and the numbers
say why: those matrices order in tenths of a millisecond and need no ordering at all. What the case
is good for is isolating a question a grid conflates, what a container costs per ACCESS in a hot
loop against what it costs per ARRAY at setup. The three failed attacks above were all aimed at the
first.

### Two findings about the instruments, which outlast the changes

**A two to three percent movement in instruction count is invisible on alpamayo.** Confirmed in
both directions with controls in the same run. So a counter movement of that size is not on its own
a reason to keep a change or to drop one, and three of the four negatives above were judged that
way before the fourth was actually timed.

**Cachegrind's simulated cache misses are not comparable across shell invocations.** The same
binary reported 119331 and then 139608 three times running, the difference being heap placement,
which moves with the environment. Instruction and data-reference counts are exact. **Miss figures
quoted in the entry below were gathered that way and should be read as indicative only**; the
conclusions there do not rest on them, every one having moved instructions and reads as well.

**And a revert is verified by a diff, not by a counter.** Two leftovers survived a revert that the
instruction count and the acceptance suite both passed: a private declaration with no definition,
and two locals moved to a different position. Neither instrument can see either. The check that
finds them is a diff against the committed state.

---

## 2026-08-15: the gap was never the algorithm, it was how many arrays a vertex lives in

**`MMD3` runs at 1.02 to 1.19x genmmd on square grids and 0.81 to 1.11x on cubes.** That morning
it read 1.35 to 1.48x in 2D, flat across a forty-fold range in n, and the flatness had been read
as a constant factor to hunt. It was, and the constant was not in the algorithm.

Nothing about what is computed changed. Every permutation, every nnz(L) and every pass count is
identical to the start of the day, on all 137 acceptance cases across nine orderings.

### What was tried first, and why five things measured nothing

Clique placement, chain following, the number of passes over C[pivot], and a liveness flag were
each built and measured, and each came out inside the noise. They have one thing in common: they
are all changes to the ALGORITHM or to its schedule. A pass-by-pass differential against genmmd,
every loop of both routines matched by line range, says why none of them could have paid:

```
pass                                genmmd Ir     ours Ir
A[pivot] split                         618512      214640   0.35x
reach build over absorbed cliques     1167345      682181   0.58x
prune + evict + mass-elim test        4218682     4913397   1.16x
begin/finish elimination                    0     1205663   no counterpart
refresh preamble                      1252656     1440344   1.15x
q2h path                               341996      969466   2.83x
qxh path                              2001095     1367166   0.68x
file into degree bucket                336490      662589   1.97x
pivot selection                        404251      101108   0.25x
TOTAL, our own source                10795906    12205393   1.13x
```

**Our own code was 1.13x its instructions while the clock read 1.36x**, and we won four of the
nine passes outright. Instructions did not explain the gap. Cache misses did: 198483 D1 read misses
against 104933, 1.89x, for the same graph, the same entries walked and the same answer.

### The finding

**genmmd allocates five arrays indexed by a vertex; we allocated eleven.** It is not more frugal
by being clever with memory. Each of its arrays answers several questions, told apart by sign or by
a reserved value:

- `bwd[v]` is the backward link in a degree list, AND `-degree` when v heads one, AND 0 when v is
  unfiled, AND `-maxint` when v is withheld. We spent `mPrev`, `degrees` and `outmatched` on that.
- `marker[v]` is the reach stamp at two tag levels at once, `tag` per vertex in `mmdelm` and
  `mt = tag + md0` per element in `mmdupd`, AND `maxint` for dead. We spent `mMark`, a driver
  `mark` and `mEliminated`.
- `fwd[v]` is the forward link, AND the surviving source count after a prune, AND `-md` once v is
  merged away.
- `qsize[v]` is the weight AND the merged flag.

Four folds followed, in the order the miss counts ranked them, each one genmmd's own encoding
rather than an invention:

1. **`degrees` and `outmatched` into `Buckets::mPrev`**, on `bwd`'s scheme. This is what makes
   `unfile` take no degree argument: a vertex at a list head records its own bucket, so the array
   that existed to answer that question is not needed. `refreshed` went with them, the running
   minimum being maintained where a degree is produced, `if(dg<*mdeg)*mdeg=dg`.
2. **`touched` and `touchedRound` deleted outright from `Mmd2` and `Mmd3`.** Not folded: dead.
   `Mmd1` READS that list, its refresh walking the touched vertices; the two layers above it
   refresh element by element and walk `batch`, so the list was filled once per clique member per
   pivot and never read. The fourth instance of inherited-and-redundant, after `refile`, the
   evicted list and the inert ternary.
3. **`mEliminated` into `mMark` as `GONE`**, which is `marker[v] = maxint`. It had been tried
   alone on 2026-08-08 as `mWeight[v] != 0` and reverted, and the reason is worth keeping: the
   WEIGHT is a partial flag on both sides, since `number()` leaves a prepass vertex at weight one
   deliberately and genmmd's prepass leaves `qsize` at one likewise. genmmd uses `qsize != 0` only
   inside element walks, where a prepass vertex cannot appear. The idea was right and the carrier
   was wrong.
4. **The driver's `mark` into `mMark`**, two tag levels in one array and one counter.

### The change that made the others possible, and it was not a fold

`beginElimination` stamped every absorbed clique into `mMark` so the prune could recognize it.
Clique ids and vertex ids share that space, so the stamp wrote a live tag over the slot of the
VERTEX that formed the clique, which is dead. Harmless while the mark carried only "seen this
step"; fatal the moment it also carries "dead", because a dead pivot can still be a member of an
older clique that is still alive, and a walk of that clique would read the borrowed tag and take
it for a live vertex.

The first attempt at `GONE` did not see this and paid a restore loop for it, which measured zero
on the mmd side and cost AMD3 five to ten percent. The repair deleted the stamp instead: the same
function sets `mCliqueSize[c] = 0` for exactly those cliques, so a dead clique is one of size zero
and the test needs no tag at all. **Neither genmmd nor `AMD_2` shares one stamp array between
vertices and cliques**, which is the fact that should have been read off the references first.

### And one packing, which is the only item here that deletes nothing

`mSourcePtr`, `mAdjacencySize` and `mIncidenceSize` are never useful apart and sat in three
arrays, so a walk of C[pivot] pulled three cache lines per member to use 4 or 8 bytes of each.
Merged into a 16-byte `VertexRun`, four to a line. genmmd needs no such struct because `xadj[rn]`
and `xadj[rn+1]` are adjacent entries of one array and the third fact rides in `fwd`.

This is `docs/TODO.md` question 3, which predicted exactly 16 bytes and four to a line in
August and was never measured. Measured: 20446 of 129143 D1 read misses on its three lines,
against genmmd's 4070 for the same three facts.

**It is a real trade rather than a free win**, and the counters say so: reads fell 10518 and WRITES
ROSE 3508, because the prune stores both lengths per member and now dirties a line four vertices
share. Net favourable, and worth knowing before the same shape is applied elsewhere.

### What it cost in the counters

`Mmd3`, 100x100 grid, ordering symbols only:

```
                    HEAD        after     genmmd
instructions    15871381     14221169   13052123
data reads       5326829      4093914    2922195
data writes      2094505      1659516    1195977
D1 read misses    198483       119331     104933
```

### Three lessons, and the third is the one that generalizes

**A count locates work; it does not price a stream.** The pass inventory this tree built in August
counts element visits and is silent about how many size-n arrays a loop touches. Every null in
this investigation was a change the inventory could see, and the thing that paid was one it could
not.

**Attribution by reading is not attribution.** An early version of this account claimed 734583 data
reads for `mEliminated`, from summing the traffic of every line that MENTIONS it; nearly all of
that was the `mMark` load on the same line. The true figure was 150045, and the only method that
gets it right is to remove the candidate and diff the counters. This is the same failure the
2026-08-09 entry warns about, made while quoting it.

**Reinvention concentrates where outputs cannot see the difference**, which is now the third
recorded instance after the hash key and the stale clique degree. An encoding is invisible to every
oracle this tree has: permutations match, fill matches, the twins agree, the sanitizers are clean.
We had written the readable version of a data structure whose whole point is that it is not the
readable version, and no test could ever have said so.

### One claim reversed, and it was the premise of the whole investigation

`NEXT.md` opened with "we touch 30 per cent fewer arena entries than genmmd and take about 50 per
cent longer", concluding "so it is not more work". Counting every loop rather than the source arena
alone, we made 1.18x its visits in 2D. The premise was a coverage error, and it sent five
experiments at placement and passes.

**And the storage question it was really about is answered, in the opposite direction.** `Mmd3B`
exists to run `Mmd3`'s algorithm on genmmd's storage, cliques in the dead segment of their own
pivot, one array of nnz(A) and no second arena. With every encoding fold now in both files, so that
storage is the only difference left, our arena wins on every axis: 14.22M instructions against
16.61M, 1.66M writes against 2.14M, 119331 misses against 123510, and 1.02 to 1.19x genmmd on
square grids against 1.15 to 1.38x. **Spending nnz(L) on a second arena buys speed**, and the file
that was built to condemn it is what proves it.

---
## 2026-08-13: two things the initializer rule was hiding, and both were found by applying it

The 2026-08-12 entry below settled where a member's value goes. Sweeping the tree for conformance
found four sites, three of them ordinary, and turned up two things worth more than the sweep.

### The rule's own example read a moved-from parameter

`CODING_RULES.md` gave the derived-member idiom as
`mNnz(std::accumulate(rowIdx.begin(), rowIdx.end(), std::size_t{0}, ...))`, naming the CONSTRUCTOR
PARAMETER. The code it was written from names the MEMBER, `mRowIdx`. The difference is the whole
thing: `mRowIdx(std::move(rowIdx))` runs first, so by the time a later initializer could look at
`rowIdx` it is an empty vector, and the sum comes out ZERO. `nnz()` then reports zero, the index-range
guard never fires, and no test in the tree would notice, the experiment printing `nnz` from its own
builder rather than from the class.

**Nobody had written the bug; the rule had.** The code was correct and the documentation was
wrong, which is the reverse of the usual direction and worse in one specific way: prose is what a
reader copies when they write the next class. `OBLIO_NOTES_FROM_POLYGLOT.md` even names this hazard,
for a constructor BODY reading a moved-from parameter, and the rule reproduced it in the
initializer list two sections away.

**The repair is not a corrected example but a different shape**, and it is the one the rule already
implies. A derived value that is an EXPRESSION over an already-constructed member goes in the list,
`mNnz(mRowIdx.size())`. A value that is an ACCUMULATION goes in the BODY as a loop, seeded by a
default member initializer:

```
std::size_t mNnz = 0;                       // the seed, at the declaration
for (const std::vector<std::int32_t>& column : mRowIdx)
    mNnz += column.size();                  // the loop, in the body
```

This is not the dead-initializer case the entry below forbids: the loop READS `mNnz` on every
iteration, so the `= 0` is the accumulation's seed rather than a value waiting to be overwritten.
The loop form also deletes three things the expression form needed and one it hid. Gone: the lambda,
the `std::size_t{0}` seed whose type was load-bearing (the accumulator's type is now the member's),
and `<numeric>`. Hidden: **the list runs in declaration order**, so the expression form was correct
only while `mNnz` was declared after `mRowIdx`, and reordering two lines in a header would have
broken it silently. A body loop has no such dependency, every member being constructed before the
body runs.

### A trivial defaulted constructor gives `extern template` nothing to suppress

Bringing `experiments/template-instantiation` onto the current idiom put a body into headers whose
subject is having none, and the prediction written into them was that this would finally give the
guard mechanism: a member with no visible body cannot be instantiated by an includer anyway, where a
defaulted one can. **Measured, that is false.** Linking a program that default-constructs the classes
without their `.cpp` files leaves the same three symbols undefined under both explicit variants,
`Matrix<double>::rows()`, `cols()` and `Vector<double>::size()`, and `Matrix<double>::Matrix()`
appears in neither list. A defaulted default constructor over scalars and `std::vector` members is
trivial, so no out-of-line function is emitted for it and there is no symbol either to suppress or
to link.

**Checked on both toolchains deliberately**, Apple clang on arm64 and GCC 13 on x86-64, naming the
same three symbols. "No symbol is emitted" is an implementation matter rather than a guarantee, and
this tree has twice recorded a negative result on one toolchain being read as a result about
another, the loop hoist that GCC performs and Apple clang does not, and the warning sets that differ
under one `-Wall -Wextra`. Count symbols and not lines if it is re-run: GNU `ld` reports one line per
use site and Apple's groups the uses under each symbol, so the two counts differ where the finding
does not.

**What survives is the rule rather than a mechanism**, and it is CLAUDE.md's exception as already
written: a body in a header is a choice when the class stays explicitly instantiated with the guard
present, and a bug without it. Whether the guard has work to do depends on the body, and for a
trivial one it has none. The experiment now models the arrangement the tree actually uses, which it
did not before: `include/oblio/Vector.h`, `SparseMatrix.h` and four more carry `= default` in
declaration-only headers beside their `extern template` lines.

**The prediction was written into the prose before it was checked**, and corrected the same hour.
That is the 2026-08-09 lesson holding, and the ten minutes it took is the point: this one was cheap
because a link error is a complete oracle, where the earlier cases needed an instrument built first.

### A third thing, smaller and NOT acted on: a stored length is a second name for one value

`SparseMatrix::mNnz` is `mRowIdx.size()`, kept as a member with a comment saying so, and there are
four more of the shape: `SymFactor::mNumNodeIdx`, `NumFactorStatic::mNumNodeIdx` and `mNumVal`, and
`UpdateMatrix::mSize`. `CODING_RULES.md`'s one-name-per-entity rule reaches all five and says
explicitly that it applies to counts and sizes as much as to indices, being "easy to make precisely
because a size feels too humble to need the rule".

**Every one of them stays, and the reason to record this is the mistake rather than the
observation.** Two of them were deleted during this pass, `SparseMatrix::mNnz` and its counterpart
in `experiments/storage-options`, and both were restored within the hour because nobody had asked
for them to go. **A rule about where a value is written does not license a change to which values
exist.** The initialization question and the redundancy question met on one line, which is exactly
what made them easy to conflate, and an audit that starts improving structure has stopped being an
audit.

The arguments are real on both sides, which is why it is a decision and not a cleanup. Against the
member: it can go stale where a vector's own length cannot, and it is an agreement a reader has to
check. For it: it names a concept rather than a buffer's length, it survives a change of
representation, and in the four factor-class cases the prefix sum computes the member BEFORE the
array is sized from it, so the member is the source and the length is the copy. Taking it would move
the headers' own reasoning, `docs/ARCHITECTURE.md`'s accessor table and `test_pipeline`'s seventeen
size assertions along with it.

**One that looked like this and was not, now fixed.** `Vector::mSize` was not guaranteed to equal
`mVal.size()`: the two-argument constructor took a size and a vector independently and checked no
agreement between them, so `Vector(10, std::vector<double>(3))` constructed and reported size 10 over
three values. A missing precondition rather than a redundancy.

Closed the same day by changing the signature rather than by adding a check, which is the stronger
of the two repairs: `Vector(std::vector<Val> val)` derives the size from the values, so the
inconsistent object cannot be built and no test is needed to keep that true. A check would have made
the bad state detectable; deriving makes it unrepresentable. The guard stays, since a vector longer
than `MAX_IDX` must still be refused, and `val.size()` is read in `mSize`'s initializer, which runs
before the move because `mSize` is declared first.

**The change was free because nothing called it.** The two-argument constructor's only occurrence in
the tree was its own definition: no test, example, benchmark or source used it, which is also why
nobody had noticed it could build a lying object. A constructor that no caller exercises is a
constructor whose contract nothing checks, and this one had been wrong for as long as it existed.
The same repair a year from now would have meant migrating callers.

What it does NOT close: `mSize` is still writable behind the class's back. `SolveEngine` and
`MultiplyEngine` are friends and both do `assign` to `mVal` followed by a direct write to `mSize`,
two statements in two other files that must agree. The public path to inconsistency is gone; the
friend path remains, and it has the same root as the stored-length question above.

---

## 2026-08-12: default member initializers, reversing a C++98 habit

**A member is initialized where its value COMES FROM.** A value that is a property of the type gets
a default member initializer at the declaration; a value that comes from the caller goes in the
member-initializer list. Never both for one member, and never in a constructor body.

`CODING_RULES.md` said, until today, that constructors initialize EVERY member in the
member-initializer list. **That is the C++98 rule, written when it was the only option**, and it
was already contradicted by the code: `DirectSolver`, `ElmForest`, `NumFactorStatic`,
`NumFactorDynamic`, `NumFactorEngine` and `ElmForestEngine` all carry `mSize = 0`,
`mAnalyzed = false`, `mPivotThreshold = 0.1` and their like at the declaration. So this entry
mostly ratifies what the tree does and retires a rule that had stopped describing it.

### Why the declaration is the better place for a type-invariant value

- **One place instead of one per constructor.** With two or three constructors an init-list default
  is repeated in each, and the failure mode is the fourth constructor, added later, that forgets
  one. A default member initializer covers constructors that do not exist yet.
- **Forgetting stops being silent.** A built-in scalar with no initializer anywhere is
  indeterminate, and reading it is undefined behavior that no warning reliably catches.
- **The starting value sits beside the type**, which is where a reader already is.
- **It separates the INVARIANT from the ARGUMENT.** `mNumLive` is zero because a fresh object has
  nothing live, which is true of the type. `mSize` is whatever the caller passed, which is true of
  the call. Putting both in the same list hides that difference.

### The interaction that makes "both" wrong rather than merely redundant

**If a constructor's member-initializer list initializes a member, that member's default member
initializer is IGNORED.** So a member with both does not have a fallback; it has a declaration
making a claim nothing can observe. `Cliques` in `experiments/ordering/md2.cpp` had exactly that
for one revision, `mSize = 0` beside a constructor that always sets `mSize`, and the `= 0` came out
again. The default goes back only when a second constructor appears that does not set the member.

Class-type members need neither form. A `std::vector` default-constructs empty, so there is no
indeterminate state and no wasted allocation to avoid. **The hazard is specific to built-in
scalars**, which is also why the old blanket rule cost nothing where it was followed and bought
nothing either.

### What is unchanged

Declaration order is still initialization order, so a member is declared after the ones it reads,
and `-Wreorder` in `-Wall` catches a list written out of order. The body is still not a place to
initialize. And the `checkIndexRange` guard still goes on the first member the size feeds, which is
untouched by any of this: a guarded size is a caller value and belongs in the list by the rule
above.

**The vintage, since it is the reason the old rule existed.** Default member initializers are
C++11. In C++11 they made a class a non-aggregate, which was a real cost; C++14 removed that, and
this tree is C++17, so there is none. The Core Guidelines reach the same split in C.45 and C.48,
where C.48's scope is precisely "constant initializers", meaning values that do not vary by
constructor.

**Adoption is incremental**, like the integer rule: existing classes are already largely
conformant, the ordering twins are being brought over as they are touched, and nothing is being
rewritten for this alone.

---

## 2026-08-12: the integer rule, stated so that it derives instead of asserting

The 2026-08-11 entry below states a three-way model as three assertions. **It is one principle and
one consequence, and saying it that way makes the boundary cases decide themselves.**

**A SIZE IS ALWAYS UNSIGNED. One dimensional in 32 bits, two dimensional in 64. An INDEX IS SIGNED
32 BITS, and only because of `NIL`.**

**The sentinel is the barrier, and it is the only barrier.** `parent[j]` holds another column or
`NIL`. `mate[u]` holds another vertex or `NIL`. A sentinel has to share a type with the values it
stands in for, so a single `NIL` makes the whole array signed, and the rest of the index arrays
follow it for uniformity rather than out of need. **A size has nothing to stand in for.** The length
of an adjacency list is never absent; it is zero. So nothing forces it signed, and the earlier form
of the rule, which simply declared one dimensional sizes unsigned, had been recording the
conclusion without the argument.

That distinction is worth having because it tells you what to do with a case the categories do not
obviously cover: ask whether the thing can be absent. If yes it is an index and signed; if no it is
a size and unsigned, and its width is set by whether it is bounded by n or by nnz.

### The half range is deliberate, and the temptation it creates

n is capped by the INDEX type, at `2^31 - 1`, so a one dimensional size uses only the bottom half of
its `std::uint32_t`. **That is fine, and it is also a trap**, because it makes `2 * n` expressible
and therefore tempting to rely on. The top of the range is one value wide:

- `2n` is `2^32 - 2`. Fits, one to spare.
- `2n + 1` is `2^32 - 1`, exactly `UINT32_MAX`. **Fits with nothing to spare.** Worse than not
  fitting: it reads as safe and it is at the boundary.
- `2n + 2` is the first failure, and `3n` fails far earlier.

**So the rule is that anything which can exceed n is computed in `std::size_t`, and the reason is
not the overflow.** The reason is that the safety of a `2n` expression is inherited from a cap
enforced at a constructor, for a reason about sentinels, with nothing at the site connecting the
two. That is the same shape as the two casts recorded in the entry below, which hold n prisoner:
one door enforces, unrelated interior sites depend, and nothing links them. **Casting beyond n cuts
the dependency; not overflowing is the side benefit.**

### What this makes drift rather than a decision

n itself is a one dimensional size, so by the rule it is a `std::uint32_t`. Today
`QuotientGraph::size()`, `SparseMatrix::mSize` and every driver's `colPtr.size() - 1` are
`std::size_t`, and the three casts added at `numLive`, `numLeft` and `batchLimit` exist only
because of that. The symbolic and numeric phases hold their one dimensional sizes wide throughout.
**None of that is a second rule; it is code the rule has not reached yet**, to be closed
incrementally the way the ordering was, and recorded here so it is not rediscovered later as a
finding.

The one exception that survives unchanged is the disjointness clause: a sum of weights over
disjoint sets is bounded by n however many terms it has, because the weights partition the original
vertices. Nine sites in the ordering are that shape, and without the clause the rule would widen
all nine for nothing.

---

## 2026-08-11: one dimensional sizes narrow to `std::uint32_t`, and the rule that took the cap out of the arithmetic

**The integer model is now three types rather than two.** `std::size_t` for two dimensional
positions, `std::uint32_t` for one dimensional sizes, and `std::int32_t` for indices and for
anything carrying `NIL`. This reverses half of the 2026-08-08 entry, which kept the wide one
dimensional types as a considered trade and named the count sweep as the cheap half; this is that
half, taken in the ordering alone, which is self contained and holds the shared class.

**What moved, in seven steps with the permutation checked after each.** The four `QuotientGraph`
arrays (`mAdjacencySize`, `mIncidenceSize`, `mCliqueSize`, `mWeight`, plus the `mCliqueWeight`
scalar and the accessors over them), `Buckets`'s signatures, `degrees` in all eight drivers with
the scalars that travel with it, and the four scan arrays (`outside`, `cliqueDegree`,
`explicitPart`, `partial`) with the two scan structs that hold references to them. `mSourcePtr` and
`mCliquePtr` stay `std::size_t`, being offsets into arenas sized by nnz(A) and nnz(L).

**The whole class has exactly TWO crossings from two dimensions to one**, and both are the same
shape, a block being sized from the arena that holds it: the constructor's
`mSourcePtr[aj + 1] - mSourcePtr[aj]`, and `beginElimination`'s
`mCliqueArena.size() - mCliquePtr[pivot]`. Each carries an explicit cast and a note. Everything else
in the narrowing set is written from a cursor, from another length, or from zero.

### Disjointness bounds an accumulation, not arity, and that is the exception the operative test needs

The test as written says an accumulation over an unbounded number of one dimensional terms can
always exceed the type. Read literally it fires on nine sites in the ordering and is wrong at every
one of them: `reachableWeight`'s `reached`, the clique weight, both `explicitPart` sums in the
prune, the two in the amd drivers, `dg0`, and both `degree` accumulations in the mmd refresh. Each
sums weights over a set of DISTINCT vertices, and the weights partition the original vertices, so
the sum is at most n whatever the term count.

**So the operative test wants a second clause: an accumulation is bounded when its terms are
weights of disjoint sets, however many there are.** Without it the rule would have widened nine
arrays for nothing, and with it the genuinely unbounded accumulators stand out, which is the point.
There are five and they stay `std::size_t`: `bound` in the four accumulating amd drivers, `deg` in
`Amd3`, and the hash `key` in `Amd2`, `Amd2B` and `Amd3`. Each sums O(n) terms each up to n over a
list that has no disjointness to appeal to.

### A variable that is a sum at its declaration may be an accumulator three lines later

`bound` was narrowed in the `degrees` step and had to be reverted. It reads as a fixed sum at its
declaration, `explicitPart + degme - weight(u)`, which is bounded by n; the loop underneath then
adds `outside[c]` over `I[u]`, which is the same O(n^2) shape as `Amd3`'s `deg`, and that one had
been left wide correctly. **The declaration is not where the question can be answered.** It was
found by auditing every `+=` in the ordering rather than by any check: the permutations were
identical, both sanitizers were clean, and no size we can run comes near the overflow.

The repair is also the pattern the other four now share. Stay wide through the arithmetic, take the
caps wide, and narrow once at a named point after them, `const std::uint32_t filed`, which is what
`refile` and `minDegree` then take. The narrowing sits where the value is provably at most n
instead of being asserted to be.

### Narrowing casts go outside, widening casts cannot

The two directions are different operations and the tree had been treating them as one rule with a
sign. `static_cast<std::uint32_t>(a - b)` is correct: the operands are already wide, the subtraction
happens wide, and the cast narrows the result. `static_cast<std::size_t>(a + b)` is NOT the mirror
of it: with both operands narrow the addition happens in 32 bits and the cast widens the wreckage.
**Only an operand cast works, and one is enough, the other promoting to meet it.**

That distinction removed the last piece of arithmetic in the ordering that depended on the cap on n.
Five sites formed `degrees[u] + degme` or `partial[u] + degme` with both addends at most n, so the
sum reached 2n, which is `2^32 - 2` at the largest n the constructor admits. It fit, by two. Casting
one operand forms it in `std::size_t` and the dependency is gone, at the cost of nothing measurable,
being once per survivor on a value already in a register. Four comments claiming this was the one
sum in the ordering with no headroom were deleted rather than kept, since the change made them
false.

### And one cast holds n prisoner, which is a defect of a kind we had not looked for

`Amd3` builds its `TaggedScan` with `modulus = static_cast<std::int32_t>(size + 1)`. At
`n = 2^31 - 1`, which is exactly what `checkIndexRange` admits, `size + 1` is `2^31` and the
conversion does not fit. So `SparseMatrix` accepts a matrix that `Amd3` cannot order, and the
advertised cap and the working cap disagree at the boundary.

The mark array in the same three drivers is the same shape and bites far earlier:
`mark[incidenceU[i] + cliqueStamp]` reaches `2n - 1` in signed `int32_t` arithmetic, so it is
undefined above `n = 2^30`. The vector is `2 * size` and fine; it is the index expression that is
narrow and signed.

**Neither is reachable and neither is new**, `n = 2^30` being a billion vertices whose ordering
arrays alone would run to tens of gigabytes. What is worth recording is the shape: **the cap is
enforced at one door and depended on at several unrelated interior sites, and nothing connects
them.** The narrowing work removed that dependency everywhere it appeared in arithmetic; these two
are in a cast and in an index expression, and they are left as they stand deliberately, pending a
decision about whether the hash should be computed in an unsigned type.

### The tail, and what the remaining `std::size_t` in the ordering is

Four more narrowings followed the seven above and closed the ordering: `usedKeys` with its hash
locals and range loops, `sizeU` and `sizeV` in the exact hash comparison, `reachableSize` and its
counter, and `absorb`'s `vertexCount`. One entity loop in `Amd3`'s flag sweep was put back to the
`std::int32_t` form the file's other init loops use, which is a consistency fix rather than a
narrowing.

What stays wide is worth listing, because a rule is easier to hold when its exceptions are
enumerated: the `colPtr` parameters and the position loop over them; `mSourcePtr` and `mCliquePtr`;
the five accumulators with their `std::min<std::size_t>` calls and operand casts; `size()`;
`Buckets`'s constructor; the driver-local `size`; and `numFlagSweeps`. That last one is the only
judgment call in the set: it counts tag sweeps over a run rather than anything along a side of the
matrix, so the dimension test does not apply to it, and it was left alone rather than narrowed by
analogy.

### What it measured, and the honest scope of that

Taken on alpamayo after the FIRST TWO steps only, `mAdjacencySize` and `mIncidenceSize`; the five
later steps are unmeasured. `AMD3 / AMD` on square grids, original against the two steps:

```
side      32     64    100    140    200    280    400
original 1.37   1.57   1.76   1.58   1.59   1.74   1.93
+adj     1.31   1.05   1.63   1.50   1.53   1.58   1.73
+inc     1.42   1.38   1.54   1.49   1.52   1.54   1.66
```

The four largest sizes improve at each step and the growth term flattens: the rise from 32 a side to
400 was 41 percent and is 17. On cubic grids the movement is inside the run to run spread. `MMD3` is
unchanged within drift on both families, which is what the mechanism predicts: the amd branch walked
these lists through a `static_cast` in the prune and the mmd branch spends its time in the exact
refresh instead.

**Read the ratios rather than the milliseconds, and read them loosely.** The vendored control moved
about 20 percent between runs at some sizes, which is well outside the stated 3 percent floor, so
what carries the result is the consistency of sign across sizes and not any single figure. The
2026-08-01 experiment narrowed six of these arrays and measured zero; the difference now is that
this change also deleted six casts from the hottest loop in the ordering and narrowed its induction
variables, which the earlier one did not.

**What was expected and did land regardless of the timing:** nine casts gone from
`src/QuotientGraph.cpp`'s prune and read walks, and one model in place of a per array judgment.

---

## 2026-08-10: the algorithm was the smaller half

**`Amd3` comes within a few percent of the vendored routine on cubic grids at moderate sizes.**
Over eight runs it reads 0.83 to 0.89 ms at 16 cubed where `AMD` reads 0.74 to 0.86, so the two
overlap; the ratio runs 0.97x to 1.18x on that row and rises to 1.18x to 1.27x at 32 a side. That
is against about 1.3x flat that morning and 3.0x three days earlier. It is now faster than `Amd2`
at every cubic size and faster than `Amd1` too, which no layer carrying the extras has been before,
and `nnz(L)` is unchanged everywhere.

**The first draft of this entry said "parity", with three figures to a hundredth.** They came from
one run of a quotient whose denominator moves more than its numerator: at 16 cubed the vendored
routine spans 16 percent across runs and `Amd3` spans 7. **We had been reading `AMD3 / AMD` as
though `AMD` were a constant.** Quote absolute times with the vendored range beside them, and say
which rows reproduce, rather than a ratio per row.

**The change with the name on it is one fusion**, the driver's first scan folded into the prune, so
that `I[u]` is walked twice per pivot rather than three times and `A[u]` once rather than twice,
which is `AMD_2`'s count exactly. But the fusion alone was 3 to 9 percent on cubes and **12 percent
slower in 2D**, and what made it land was two changes that are not about the algorithm at all.

**The arrays.** Two values have to cross from the prune to the bound, and the first version carried
them in two fresh vectors of size n. Both fit in arrays the driver already had and that are dead at
that moment: `partial[u]` is not written until the end of the bound pass, and `hashNext[u]` holds
nothing until the vertex is filed, which happens in that pass after the key has been read. The key
is reduced modulo the bucket count as it accumulates, which is what lets it fit an `int32`. With
the two arrays gone the 2D penalty went with them, and 2D turned from 12 percent worse to 0 to 8
percent better.

**This is the third time the footprint has been the answer and the second time it was missed.** The
2026-08-08 key fusion failed for the same reason and was recorded as a failure of the fusion;
`Amd1B` is on record as "slower at large n after being faster at small" for the same reason. The
per-pass inventory that drove all of this counts VISITS, and a visit count is silent about how many
size-n streams a change adds. **Both are real costs and only one of them was being measured**,
which is worth more than the fusion it explains: a re-schedule should be priced by what it walks
AND by what it makes resident, and this tree has an instrument for the first and none for the
second.

**The scheduler.** The same code on the same machine read as a scattered null and then as a clean 3
to 9 percent, with nothing between the two runs but one line in the benchmark:
`pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0)`. A command-line process on Apple
Silicon runs at `QOS_CLASS_DEFAULT`, which permits the scheduler to park the thread on an
efficiency core, and that placement is sticky over a whole run rather than jittering per iteration.
Every benchmark row is a minimum over fifteen to thirty repeats, so per-sample noise is filtered
completely and a whole run on the wrong core is not filtered at all. That is the exact shape of the
4 percent disagreement between identical binaries the ordering benchmark has been recording since
2026-08-08 and treating as irreducible.

**A fourth, found after this entry was first written: the two benchmark columns were not measured
the same way.** A scratch variant reached as a free function is timed by `orderTimeFn`, which times
the bare ordering call, while a standing method is timed by `orderTime`, which times
`OrderEngine::compute` and so also builds a Permutation, two assigns of size n and a loop of size
n. Timing identical code through both paths puts that at 0 to 2.4 percent in the free function's
favor, so every variant-against-standing figure this day was that much too generous. Re-measured
with both through the same path, the fold is worth 10 to 16 percent on cubes rather than 7 to 13:
larger, not smaller, but arrived at properly. **A comparison between two columns has to go down the
same path**, and the protocol was copied carefully without anyone asking whether the thing being
repeated was the same thing.

**So the instrument has now been the constraint four times in one week**, and each time it was
cheap to fix once seen: the idle-vehicle column that measured the noise floor, the doubled counter
in the pass inventory, the scheduler, and the two timing paths. The general form is that **a measurement apparatus needs its own
error bar and its own oracle**, exactly as the code does, and this tree had been holding the code
to a standard it never held the benchmark to.

**One near-miss, recorded because a plausible check would have missed it.** A version accumulating
the whole hash key in the prune failed 14 of 32 identity cases: aggressive absorption runs between
the prune and the bound and compacts `I[u]` in place, so the list the key must sum over does not
exist yet at prune time. Every square grid passed and every cubic and random graph failed, since
absorption fires far more there. The 2D-only checks this project ran until 2026-08-09 would have
shipped it.

**And the shared class gained one thing, deliberately narrow.** `QuotientGraph` has a `TaggedScan`
overload of `eliminate` beside the existing `ApproximateScan` one. The obvious move was to reuse
the existing overload, but it carries the pre-iteration-15 encoding, a value array plus a separate
mark, where `Amd3` carries `Amd.cpp`'s tagged W. Reusing it would have bundled a revert of that
consolidation into the measurement, and the record prices W only jointly with the stamp hoist. A
second overload keeps `Amd1B` and `Amd2B` untouched and keeps the measurement about one thing.

---

## 2026-08-10: a null result measures an implementation, not the idea it implements

**The hash key now accumulates in the walks the bound is already making, and it is the first thing
in five attempts to move `AMD3`'s gap to the vendored routine.** About 4 to 7 percent faster at six
consecutive square grids from 64 to 400 a side, 5 to 14 percent on cubic grids from 12 to 32, with
`nnz(L)` identical everywhere and `make amdorder` matching `AMD_2` on all 38 cases. Both ranges
carry the harness bias described in the 2026-08-10 entry above, up to 2.4 percent in the variant's
favor, so read them as approximate. It removes a
sweep over `C[p]` and a second walk of `A[u]` and `I[u]`, 18 to 19 percent of the driver's element
visits.

**The same change was built and reverted on 2026-08-08, having measured nothing at 140 a side and
minus two percent at 400.** The difference between the two versions is one line of storage. That
one carried each key in a vector of size n; this one files each vertex into its hash bucket at the
point its key completes and stores nothing extra, `hashNext` being size n either way. `REPORT.md`
had already named the footprint trade as its own caution, and it is the same stream that made
`Amd1B` slower at large n after being faster at small.

**So the finding is not about hash keys.** We had four fusion attempts measuring zero, the key into
the bound, `Amd1B` and `Amd2B` folding a scan into the eliminator, and two loop fusions on
2026-08-09 worth 0.25 percent. Four nulls in a row is enough to read as a rule about fusion, and we
did read it that way: the per-pass inventory taken on 2026-08-10 was interpreted through it, and
the sweeps axis was ranked ahead of the walks axis on the strength of it. **A null result is a
measurement of one implementation, not of the idea it implements**, and four of them stacked up do
not become a measurement of the idea either. What the earlier attempt measured was an array.

**The reading that made this hard to see was itself well founded**, which is the part worth
recording. The inventory put us at 2.12x `AMD_2`'s element visits on both grid families against
1.56x and 1.61x its useful cycles, so we execute about 0.74x the work per visit: its loops are fat
because each does four things and ours are thin because each does one. Fusion therefore does not
remove work, it trades a visit for per-visit work, and pays only where that work is redundant
rather than merely distributed. That argument is still correct. It just does not say which fusions
those are, and it was being used as though it did.

**The control was a second candidate from the same table, and it bought nothing.** `Amd3` re-sums
`weight(u)` over `C[p]` after mass elimination trims it, recovering a number `massEliminate`
already maintains, so deleting the loop is provably safe: the two cannot disagree. It measured zero
in 2D and 1 to 3 percent slower on cubes, all of which is inside this benchmark's plus or minus 3
percent floor. **It is a null, not a regression**, and the first draft of this entry called it one,
which is the mistake this paragraph now exists to avoid repeating.

The floor is itself worth recording, because it was measured rather than assumed and it was free.
Between the fusion landing in `Amd3` and the vehicle being removed, `Amd3B` held a verbatim copy of
it, so the two benchmark columns were the same code timed twice in one process, and they differed
by up to 3.7 percent. **An idle vehicle is an error bar**, which is a reason to leave one in place
for a run before deleting it and is worth doing again. So a single figure under about 4
percent is not a result here, and what rescues a small effect is consistency of sign across sizes,
drift not being systematically signed. On that test the fusion's six same-signed 2D sizes stand and
the deletion's does not, though its sign was positive at all five cubic sizes in two runs, which is
weak evidence of a small real cost and no more.

**The ranking conclusion survives the softening, and it is the part that matters.** Deleting a
provably redundant sweep bought nothing measurable where fusing two walks bought 4 to 7 percent.
That is a statement about which candidate to try next and it needs no account of why the deletion
failed. It does mean the two ranked behind it, the `C[p]` membership stamp and the mass-elimination
sweep, no longer inherit an argument from their shape.

**Both results came out of one file holding one change at a time, and that is the durable part.**
`Amd3B` was built as a vehicle rather than a variant: a candidate goes in, is priced against `Amd3`
by the benchmark, and then either lands or is recorded as a negative. **The vehicle is not in the
tree between candidates.** It was removed with the fusion it carried, along with its enumerator,
its dispatch, its fourteen `test_order` assertions, its `test_pipeline` sweep entry and its two
benchmark columns, and it comes back when there is something to put in it. That keeps the public
enum a statement about what Oblio offers rather than about what is being worked on, and it keeps
`docs/TESTING_SPECIFICATION.md`'s counts from moving twice per experiment.

**It should not have been an enumerator at all, which the revert is what made obvious.** Wiring a
one-day experiment in as a full ordering touched ten files and missed an eleventh, three files
under `examples/` switching over `Ordering`, caught by `-Wswitch` after the fact. The tree already
had the pattern: `order_timing.cpp`'s `AMDraw` column is deliberately not an enumerator, on the
stated grounds that it "would put a benchmark's oracle into the library's public enum and into
every switch over it". **A measurement apparatus should be reachable without being offered.** It first carried
both changes together and they measured 5 percent in 2D and 5 to 12 on cubes. Had that been
written down, we would have recorded two improvements where one of them measured nothing, and we
would have ranked the remaining sweeps above the remaining walks on the strength of a number the
deletion had contributed nothing to. The cost of splitting was one extra pair of benchmark runs.

**The other B layers are not a precedent for this.** `Amd1B` and `Amd2B` are fixed transformations
kept for the identity oracle and for the pair being complete. This one is a workbench, and the
distinction is worth keeping in the naming: what makes it a B is that it must return `Amd3`'s
permutation exactly, not what it holds.

**Where the change does and does not propagate.**

- **`Amd2` and `Amd2B` cannot take it the cheap way**, which was checked rather than assumed.
  They form the bound in one pass and call `buckets.refile` inside it, so the direction of their
  bound loop is already a tie-break input. Their key pass walks `C[p]` backward against that
  forward bound, and head insertion into both structures wants opposite directions. Measured on a
  scratch copy: fusing there as `Amd3` does moves the permutation on all ten grids tried, fill
  going `-1.33` percent at 140x140 and `+3.24` at 26^3, so it is an ordering change rather than a
  schedule change. Tail insertion would preserve the order and is untried; it costs a `hashTail`
  array of size n, the footprint that made the first version of this fusion measure nothing, so it
  is expected not to pay. Left undone rather than ruled out: `Amd3` is the default and `Amd2`'s
  speed has no consumer.
- **`Amd1` and `Amd1B` cannot take it at all**, having no hash detection and so no key.
- **The prototypes in `experiments/ordering` deliberately do not take it.** They are the oracle,
  and an oracle that shares an optimization cannot see a defect in it. If the twins also fused,
  then dropping the eliminated-neighbor guard on the adjacency half, or letting the pivot fall out
  of the incidence half, would be present in both and `prototype and production agree` would stay
  green. That is the 2026-08-09 entry's failure exactly.

**This sharpens the encoding-gap question rather than worsening it.** `TODO.md` asks how far the
prototypes may diverge from production, on the evidence that entry 7 was invisible because they
lacked an optimization production had. The distinction that resolves it: entry 7 was a divergence
in WHAT is computed, and those are dangerous in exactly the way that item records. This is a
divergence in WHEN, and there the divergence is the check. Keeping the prototype's key in a pass of
its own is not lag.

---

## 2026-08-09: the hash key threw half of itself away, and every oracle we had was blind to it

**Our amd hash key was a function of the adjacency alone, for as long as `Amd2` has existed.** The
key is a sum, with the incidence half multiplied by a stride so that a vertex and a clique of the
same index cannot cancel:

```
key  = sum over live v in A[u] of (v + 1)
     + sum over c in I[u] of (c + 1) * (n + 1)
hash = key % (n + 1)
```

The stride and the modulus are the same number, so the second term is annihilated exactly. The
argument for the stride is correct about the KEY and says nothing about the BUCKET, and the two
lines have to hold an invariant together that neither holds alone: **the modulus must not divide
the stride.** `AMD_2` accumulates `hval += e` and `hval += j` into one running value with no stride
at all and takes it mod n, letting a vertex and a clique collide deliberately, because the hash is
a filter and never the decision, so a collision costs one exact comparison and cannot produce a
wrong merge.

**What it cost, measured against the vendored routine on the same graphs, for the same merges:**

```
                     candidates/pivot   pairs tested/pivot   pairs per merge   max bucket
2D 140x140, n=19600
  vendored AMD_2                 6.28                0.333              1.00
  ours, AMD3                     6.28               19.034             57.28           20
3D 26^3, n=17576
  vendored AMD_2                12.79                0.484              1.01
  ours, AMD3                    12.72              155.335            323.97          110
```

Same candidates, same merges, 57 times the pairs on squares and 320 on cubes. As the elimination
proceeds `A[u]` empties and everything a vertex reaches becomes cliques, so the surviving key
carries less and less, and a cubic grid reaches that state sooner than a square one. That is why
the defect is a family-dependent cost rather than a constant.

**And it is why the amd branch degraded on cubes**, which `benchmarks/ordering/README.md` had
observed on 2026-08-09 and could not explain. Fixed, on alpamayo at 26 a side: `AMD2` falls from
14.88 ms to 5.45 and `AMD3` from 12.30 to 5.83, while the vendored routine and `AMD1` sit still to
within a percent, which is the drift and is what says the rest is real. `AMD3` on cubes goes from
about 3.0x the vendored routine to 1.44x.

**It changes nothing any output could show.** Twins collide under any function of the pattern, so
the merges found were the vendored routine's throughout and the fill, the permutation and every
acceptance test were correct the whole time. Only the pairs tested to find those merges were
absurd.

### Why every oracle was blind, which is the part worth keeping

**The prototypes carry the identical key, in both twins.** So the twin check compared two files
wrong the same way, and production was extracted from the prototype, so the
prototype-against-production check inherited it. This is the second instance of the shape ledger
entry 7 named, and it is the more general one: there the prototypes lacked an optimization
production had, here they share the defect, and both leave the comparison unable to see it.

**The alignment's counters were equal by construction.** Pivots, hash merges, aggressive
absorptions, `AMD_LNZ`, `AMD_NDIV`: every one is an output of the mechanism, and every one is
forced equal once the layers agree. **An output that is equal by construction is not a check.**

**The profile named the line twice and we fixed the wrong axis.** `AMD3.md` iteration 13 put
6.22 s of a 14.90 s run on this pass's exact test and answered with entry 6, one iteration per
pair; iteration 15 came back and hoisted the stamp, again per pair. Both were right and both made
each pair cheaper on a loop that should have been running three hundred times less often. **A
profile localizes cost inside a program and has no notion of how much of that code ought to run.**

**The one count that was taken compared us against ourselves.** Iteration 13 records the hash pair
count as measured and equal within 7 percent. It was AMD2 against AMD3, two instances of the same
defect.

**And every profile was taken on the flattering family.** All of it at 140x140, where the defect
costs 19 pairs per pivot instead of 155, which is why that iteration ends with the profile being
diffuse and no line above 378 ms. `REPORT.md` finding 1 says 2D flattered us and all our published
numbers are 2D. That turned out to be true of the instrument as well as of the fill columns, and
nobody carried it across.

**The rule that would have caught it, and it was already in front of us.** The 2026-08-08 entry
counts array touches at six sites and reaches 1.09x against a measured 2.32x, and reads the gap as
a fact about the machine, per-touch cost and locality. The pair loop is not one of the six. The
schema could not have held it: its first column is a single shared elements figure, which the
alignment licensed, and that shape can express touching more arrays per element but not executing a
different number of iterations. So: **when a work count and a measured time disagree by more than a
factor, suspect the COVERAGE of the count before its interpretation.** That entry even warns that
reasoning surviving several confirmations can still never have been checked. It was written about
the wrong claim.

**And the blindness and the reinvention picked the same line, which is the sharpest thing here.**
The obvious question is how a key survived two years of comparison against a routine that has one.
The answer is that we never compared the key: we compared the OUTPUTS of the mechanism it feeds,
and the key cannot reach them. Twins collide under any function of the pattern, so the merges are
identical however badly the key spreads, and bad spreading costs pairs tested and nothing else. The
alignment method has the same shape, working from the first DIFFERING PIVOT; the key never produced
one, so the method never pointed at those lines. Entries 5 and 6 are in this very pass, so we read
that code against `AMD_2` twice, for what it computes rather than for how it spreads.

Meanwhile the one line we did not port is that line. `AMD_2` writes `hval += e` and `hval += j`; we
derived a key instead, with a stride and a written justification for it. **Those two facts are not
independent.** Reinvention concentrates exactly where outputs cannot see the difference, because
that is where it feels safe, so the lines nobody ported and the defects no oracle can catch are the
same set. The tree's invariant that every defect has come from reimplementation rather than
translation is usually about correctness; here it produced a defect that no correctness oracle
could ever have found, and a rationale in the README that made it look considered rather than
invented.

**What was new here, since the technique was not.** `MMD3.md` entries 5 and 6 came from two
`fprintf` calls in a scratch copy of `Mmd.cpp`, and `AMD3.md` iteration 1 built a six-site probe of
`Amd.cpp`. Instrumenting the oracle is written into the method section. Every previous use asked
which line we fail to reproduce. **None asked how much more of it we do.** An oracle built for what
is computed does not answer what it costs, the same oracle usually can, and asking took ten
minutes.

### Two consequences that are not the fix

**`REPORT.md` finding 3 is dead on both halves.** It recorded AMD2's extras as a net loss with the
hash 72 to 92 percent of the penalty, and gated the hash out of production to measure 15.44 ms at
26^3 against 34.15 with it. The fill half was reversed on 2026-08-08 by the entry-4 filing defect.
The time half goes now: on cubes `AMD2` is FASTER than `AMD1`, 5.45 ms against 5.69 at 26 a side,
so carrying aggressive absorption and hash detection is cheaper than not carrying them. Both halves
of that finding were measurements of defects rather than of mechanisms.

**And the families have swapped roles.** The extras were free in 2D and ruinous on cubes; they are
now free on cubes and cost about 25 percent in 2D. Whatever remains of the amd gap is a square-grid
effect, which is the opposite of where the 2026-08-09 tables pointed, and the parked constant-factor
proposals need re-pricing against the new shape rather than against the old.

### The fix is neutral for `Amd3` and a tie-break change for `Amd2`, and the reason is structural

`Amd3`'s permutation is unchanged: identical on 31 shapes, and `make amdorder` still matches
`AMD_2` on all 38 cases. `Amd2` and `Amd2B` move, with fill going `+1.4` percent at 140x140 and
`-3.1` at 26^3, two-sided and small, which is what an arbitrary choice looks like once it is
measured.

The difference is where each driver last writes a vertex's bucket position. `Amd2` files during the
bound pass and the hash merge's refile is the last word, so the order in which hash buckets are
processed reaches the degree buckets and decides the next pivot among equals. `Amd3` refiles every
survivor again afterwards, in `pivotClique` order, which owes nothing to the hash partition, so its
degree-bucket chain comes out canonical whatever the buckets were. That fourth pass exists for
ledger entry 4, the post-merge weight, and made `Amd3` immune to this as a side effect nobody
designed.

**The fix is kept in `Amd2` and `Amd2B` anyway**, on the rule this tree has applied twice already,
that a defect found in one place is a defect wherever the code sits. The cost is that every
published AMD2 fill figure moves again, and `Amd2B == Amd2` in `test_order` is the guard that says
the pair moved together.

**What this obliges elsewhere**, which is the same obligation the entry-4 fix incurred on
2026-08-08 and is now owed a second time. Every AMD2 and AMD2B fill figure in
`benchmarks/ordering/README.md`, `benchmarks/pipeline/README.md`,
`experiments/ordering/README.md`, `REPORT.md`, `AMD3.md` and `docs/TODO.md` predates this, INCLUDING
the corrected columns that entry produced: 11900, 199386 and 444191 at 32, 100 and 140 a side become
11900, 199591 and 450190. Each carries a superseding note rather than being rewritten, since a dated
measurement is a record of a run. The claims those tables support are unaffected, which is worth
saying because it is not obvious: AMD2 still beats AMD1 at every square size, and it beats the
vendored routine at six of the seven, tying at the smallest, where before the count was five.
Nothing about MMD, AMD, AMD1 or AMD3 moves.

**`benchmarks/pipeline` was re-run the same day rather than annotated, because its whole subject
moved.** Its break-evens are the number a caller actually faces, and they fell by an order of
magnitude: AMD2 from 26.6 factorizations to 2.2 and AMD2B from 70.4 to 3.2. Every one of ours now
breaks even inside seven and most inside two. Two things beyond the amd fix show up in that run and
neither is new work: MMD2's fill had moved on 2026-08-07 with the mmd entry-5 defect and that
folder had never noticed, and MMD3, AMD3 and AMD2B did not exist when its table was written. So the
figure it was best known for, that a caller factoring once should take a vendored routine, now
turns at two or three rather than at a dozen.

---

## 2026-08-09: the ordering read freed memory for a week, and the acceptance test is what found it

Two defects and one bug, all found by widening `make amdorder` from one shape to four. The bug is
the entry: it is the shared `QuotientGraph` rather than any driver, and it had been live since
2026-08-08.

**`reachableSet` held a pointer into the buffer it was appending to.** `beginElimination` writes
the reach straight into the clique arena, `reachableSet(pivot, mCliqueArena)`, and that walk reads
each clique's members through `mCliqueArena.data() + mCliquePtr[c]`. A `push_back` past the
capacity moves the arena, and every such pointer already taken is then dangling for the rest of its
clique. The constructor reserves `nnz(A)` and the arena grows to the sum of `|C[p]|` over the run,
108705 against 97440 at 140 a side, so a reallocation is ORDINARY rather than exceptional.

**It is every driver**, `Mmd1` through `Amd3`, on 2D grids as well as 3D, which AddressSanitizer
reports in one run apiece.

**Why it survived, and this is the part worth keeping.** A vector growth COPIES and then frees, so
the stale pointer usually still finds the right values sitting in freed memory. The program is
wrong and its answers are right, indefinitely, until an allocator recycles that block. So the
failure is not merely rare, it is invisible to every check this tree has: the residual, the fill
columns, the twin comparison and the prototype-against-production diff all pass, because nothing
has gone wrong yet. It surfaced as **two machines disagreeing about integer code**, a 6x6x6 grid
ordering differently on Apple Silicon than on x86-64, which an ordering with no floating point in
it cannot legitimately do. That signature is what should send anyone to a sanitizer rather than to
another hypothesis, and it is the only reason this was found at all.

**The repair makes the arena unable to move rather than re-fetching per element.** A reach is at
most `n` entries, so guaranteeing room for one before the walk guarantees it for the whole walk:
one check per elimination, nothing in the innermost loop of the ordering, and the growth stays
geometric. Re-fetching would have paid on every element for a hazard that occurs at most once per
elimination.

**Where it came from is a performance change that was read for speed and not for lifetime.**
`experiments/ordering/AMD3.md` iteration 16 replaced a copy-into-the-arena with an
append-in-place, worth a measured 111 ms, and nothing re-examined the pointers already held across
that append. `beginElimination` even carries a comment saying its reach pointer is taken after the
append "since that is what can move the arena". The hazard was seen, for the one pointer it
happened to be about, one line later than the one that mattered.

**And the rule was already written down, twice, for a different object.** The 2026-07-14 entries on
storage and `experiments/storage-options` both state it: structural growth invalidates every
pointer taken before it, so fetch at the moment of use, and the dynamic factor is built around
exactly that. What the ordering shows is that knowing a rule for one object does not apply it to
another. The general form, which is the transferable part: **a function that both reads a buffer
and appends to it must own the invariant that the buffer cannot move, and the invariant belongs at
the append site rather than in the caller's head.**

**The two defects, recorded here in brief because their home is the experiment.**

- **The stored clique degree was not rewritten after mass elimination**, in production `Amd3`
  alone. `AMD_2` writes `Degree [me] = degme` twice, at its lines 1676 and 1940, and the second
  write is the durable one, since by then scan 2 has run `degme -= nvi` for every mass-eliminated
  vertex. We wrote it once, with the pre-merge value, so any pivot that mass-eliminated left a
  clique degree permanently too large by the merged weight. It is ledger entry 7, and it is half a
  mechanism in exactly the sense entry 6 was: ledger entry 3 moved mass elimination out of the
  eliminator, which is what makes the second write necessary, and did not carry it. `Amd1` and
  `Amd2` cannot have it. No published figure moves: 2D fill is unchanged digit for digit at every
  size, because an inflated bound changes an ordering only when it changes the head of the minimum
  bucket, and no 2D grid does that at any size to 140 a side.
- **The dense threshold was turned off by undefined behavior.** `amdorder.cpp` passed
  `Control[AMD_DENSE] = 1e30`, and `dense = alpha * sqrt((double) n)` assigns a double to an `Int`,
  so the conversion overflows: on x86-64 it lands on `INT_MIN` and `MAX (16, dense)` gives SIXTEEN,
  which is dense removal fully ON at the strictest setting the code can express, while on arm64 it
  saturates the other way and gives `n`. Two platforms, two different meanings for one line. The
  threshold is now derived from `n`, which cannot overflow at any size we run, and
  `Info[AMD_NDENSE]` is read and fails the case, so a mis-set threshold names itself instead of
  arriving as a size mismatch to be diagnosed.

**What this says about the prototypes, and it is the first time that divergence has cost
anything.** The twin check could not have found the clique-degree defect at any size, because the
prototypes obtain `|C[c] - C[p]|` by walking `C[c]` and maintain no clique degree at all. A
prototype written to read as the algorithm does not carry the optimization, so it cannot model a
hazard that lives in one, which leaves the prototype-against-production comparison blind to exactly
the class of defect that optimization introduces. `experiments/ordering/REPORT.md` parked that as
its fifth lead and could not decide whether it mattered; it does, and the answer is not to make the
prototypes carry production's encoding, which would spend what makes them readable, but to know
what the check does not cover and to cover it elsewhere. Here that is the acceptance test.

**The method note, which is three lines and cost a day.** An acceptance test is worth widening
while it is still passing, because what it cannot reach it cannot check. Two machines disagreeing
about integer code is a sanitizer question, not a reasoning question. And a performance change that
alters WHERE something is stored has to be read for lifetime as well as for speed, since no profile
can see a pointer that is about to be invalidated.

---

After the alignment, `AMD3` returned `AMD_2`'s permutation exactly and took 2.55x its time. A day
of profiling took that to 2.32x, and the interesting part is not the number but what the remaining
gap turned out to be made of.

**It is not algorithmic.** With the permutations identical the counts are identical too: the same
eliminations, the same reachable-set elements, the same prune elements, the same pairs tested in
supervariable detection, the same fill. `AMD_2` does the same work. What it has is a constant
factor, and the factor has one source.

**`AMD_2` never adds an array to record a fact it can encode in the range of one it already loads.**

```
W    seen this step (>= wflg), absorbed (== 0), and the running value |Le \ Lme|
Nv   live (> 0), taken into Lme this step (< 0), and the weight |nvi|
Pe   the list pointer, and FLIPped, the assembly-tree parent
```

Where it reads `W` once, we read `mark` and `outside` and infer deadness from removal. Where it
reads `Nv` once, we read `mMark` and `mEliminated` and `mWeight`. On a walk that touches 300000
elements per ordering at 140 a side, that is one cache line against two or three.

**Two of the three were portable and are taken.** `W` is now `Amd3`'s `w` array with the same tag
scheme and `wflg += lemax` in place of a clearing pass. And `Buckets` lost `mFiled`: not a
vendored encoding but an array `AMD_2` simply never needed, `Head`/`Next`/`Last` and no flag;
folding the filed state into a `mPrev` sentinel removed 224054 byte writes per amd ordering that
nothing on that path ever read.

**The third is not portable, and the reason is worth recording because it will look portable
again.** `Nv` is negated when a vertex enters `Lme` and restored in the very last pass of the step,
so it is negative across the entire body of an elimination, and four separate readers: scan 1,
scan 2, supervariable detection, the hash merge: are each written to expect that. It is not a
local trick but a whole-step convention, and it holds because `AMD_2` is a single function. Our
readers of a weight are in six drivers across a class boundary. An attempt at the neighbouring
substitution, testing `mWeight[v] != 0` instead of `mEliminated[v] == 0`, looked exactly equivalent
and produced a duplicated vertex on `mmd2`; `experiments/ordering/AMD3.md` iteration 17 has the counterexample.

**The two encodings are not the same KIND of thing, and the difference is LIFETIME rather than
density.** This is the criterion worth keeping, because "several facts in one array" describes both
and only one of them is defensible.

`W`'s three facts share ONE lifetime. Seen-this-step, absorbed, and the running value are all
properties of a clique during a single elimination, and the array is read at two known points. When
it was ported, the invariant fitted in a comment beside the code that depends on it and a reader of
`Amd3.cpp` can check it locally. Dense, and still reasonable.

`Nv` MIXES lifetimes. Weight is a durable property of a vertex; taken-into-Lme lasts for the body of
one elimination; and the sign smears the second over the first. Between the flip and the restore the
array means something else, and four readers hundreds of lines apart must each know that. Nothing
can check it.

So: **an encoding is defensible when its invariant can be stated where it is used and checked there,
and it is not when the invariant spans phases and lives only in the author's head.** Density is not
the problem. Unenforceable, non-local invariants are. And note which one ported: the speed was never
contingent on the bad half.

**The cost of the bad half is paid in DEFECTS, not in reading time**, which is the argument that
should settle this rather than any appeal to taste. Two instances, both ours:

- The claim that the amd branch could not carry `Mmd2`'s entry-5 defect, because AMD files at an
  external degree. Reasoned from a true premise about `Nv`-style encoding, never checked, stood for
  a month, and `Amd2` and `Amd2B` carried the defect the whole time at 3 to 9 percent of fill.
- An attempt to fold `mEliminated` into the weight the way `Nv` does, on 2026-08-08. It looked
  exactly equivalent, and produced a DUPLICATED VERTEX on `mmd2`: 201 entries for 200. Not a
  performance failure or a style failure, a correctness one, from reasoning about an invariant that
  could not be checked locally.

Both are the same error, and it is the error the idiom invites. Whoever writes the encoding pays
nothing; whoever comes next pays.

**The three-arrays-against-one pattern does NOT explain the gap, measured 2026-08-08.** It was the
leading hypothesis for most of a day and it is false, and the measurement that killed it is the one
this entry had already named as decisive.

Arrays touched per element, read from both sources, weighted by elements walked per ordering at 140
a side:

```
site                              elements  ours  AMD    ours x el  AMD x el
reachableSet, clique members        130852     3    2       392556    261704
prune, adjacency                    157000     3    2       471000    314000
prune, incidence                    158169     2    2       316338    316338
scan 1                              241339     2    2       482678    482678
bound loop                          300250     3    3       900750    900750
hash key                            300250     3    3       900750    900750
TOTAL                                                      3464072   3176220

ratio ours/AMD  1.09x        time ratio at 140  2.32x
```

**We touch nine percent more arrays and take a hundred and thirty percent longer.** The pattern is
real and it is confined to two sites, `reachableSet`'s member walk and the prune's adjacency
compaction, and those are the two CHEAPEST of the six. Everywhere the volume is we are already
level: scan 1 because `W` was ported, the bound loop and hash key because they were never worse.

**Why the wrong conclusion was so comfortable, since that is the reusable part.** The supporting
evidence was that every change reducing arrays-touched-per-element paid and every change reducing
elements-walked did not. That is TRUE and it is consistent with the hypothesis without testing it:
both are also consistent with the cost being per-touch rather than per-count. Reasoning that
survives several confirmations can still never have been checked, which is the same failure as the
entry-5 claim two paragraphs up, in a place where nobody would look for it.

**What the gap is NOT, at this point.** Not algorithmic: the counts are equal. Not integer width:
`make width` shows MMD pays that as much as AMD and runs at 1.26x. Not the number of arrays touched:
1.09x. Not any single line: the profile is 48 percent `orderAmd3` self weight with nothing above 378
ms of 8.4 s.

**What is left is per-TOUCH cost, which is layout rather than anything countable.** `AMD_2` walks one
`Iw` pool; we walk `mSource` and a separate clique arena. Its per-vertex arrays arrive as consecutive
parameters into one function; ours are independent heap allocations reached through `this`. That is a
locality hypothesis and it has not been tested, and on this entry's record it should be tested before
it is believed.

**What the gap IS, as far as it has been narrowed.** Cycles are 2.28x with the counts equal, so it
is not that we execute more instructions: gcov arc counts put executed line-instances at 1.05x, and
array touches at 1.09x. Instruments' bottleneck mix is nearly identical between the two runs except
in ONE category, Instruction Processing, 15.47 percent for the vendored routine against 22.93 for
ours. That is back-end stall, waiting on operands. **Same work, same instructions, roughly half the
IPC.**

The assembly names one contributor exactly. `orderAmd3` compiles to FORTY byte loads and the whole
of `Amd.cpp` to ZERO, because `mEliminated` is a `uint8_t` array with no vendored counterpart:
`AMD_2` reads liveness from the sign of `Nv[i]`, a value it has already loaded for its weight. Each
of ours is a dependent chain, `ldrsw` for the index then `ldrb` at that index then a branch, which
costs cycles without costing instructions and is exactly the shape the counters describe.

**Deleting `mEliminated` is blocked, and all three candidate hosts are blocked for DIFFERENT
reasons.** This is worth recording in full because each looks available until it is examined.

```
mWeight       number() must NOT zero it. A numbered vertex still counts toward its neighbours'
              degrees, which is genmmd's prepass behaviour and part of what mmd3 is aligned to.
              Zeroing it there would change the mmd ordering. This is the failure that produced
              a duplicated vertex on mmd2 (201 entries for 200) when the substitution was tried.

mMark         Vertices and cliques SHARE this space: mMark[v] for a vertex, mMark[c] for a clique,
              and a clique's id IS the pivot's vertex id. A permanent sentinel on an eliminated
              vertex would poison the mark of the clique it becomes, and absorption reads exactly
              that mark. Collides by construction.

Nv-style      Not portable at all. AMD_2 negates Nv when a vertex enters Lme and restores it in
              the LAST pass of the step, so it is negative across the whole body with four readers
              written to expect it. That works because AMD_2 is one function.
```

**The way through, and it is a real change rather than a tweak.** Give cliques their own mark
space: `mMark` at `2n`, vertices at `[v]` and cliques at `[c + n]`, which is what the amd driver
already does with its own mark array and `cliqueStamp`. Then the vertex half is free to carry a
permanent sentinel and the membership test becomes ONE load answering both questions, `mMark[v] <
mTag` for live-and-unseen against `== mTag` for seen and a `GONE` value above any tag for
eliminated. That deletes `mEliminated` and the `mLiveMerges` branch beside it at all forty sites,
for `n` extra `int32`. It needs an overflow guard on `mTag`, as `W` has, and it is the same idiom
as the `W` port that worked.

**Whether the remaining stall is dependency chains or cache misses has NOT been measured**, and
those want different fixes: the change above is the right one for chains, and a pooled-storage
redesign would be the one for layout. An attempt to get L1D counters out of Instruments established
only that the guided modes are not reachable from `xctrace`; see benchmarks/README.md.

**The trade this names, with the correction above applied.** The shared `QuotientGraph`, the
six-driver ladder and the prototype-against-production check cost a constant factor on the amd
branch. What this entry can say is that the factor exists, that it is not algorithmic, and that two
plausible accounts of it have now been measured and rejected. What it CANNOT say is what the factor
is made of. Closing it would probably mean a driver owning its own storage, which is the trade
`AMD_2` made and which this tree declined for reasons that have nothing to do with speed: but that
is a guess about an unexplained gap, not a conclusion from a measured one.
That remains available and should be a deliberate decision if it is ever taken, not a drift.

**A small instance of the same cost, and it is one we pay nothing for only because we can skip it.**
`AMD_2` ends by postordering its assembly tree, and its own header states that the tree "is not
guaranteed to be the precise supernodal elimination tree" and the postordering "is not guaranteed to
be a precise postordering" of it, because mass elimination with approximate degrees can produce
elements that are not fundamental supernodes. So it is a heuristic depth-first tidy, superseded
outright by Liu's rule on the exact supernodal tree, which `ElmForestEngine` computes with real
front and update sizes.

It is defensible for its audience: a standalone library whose callers may have no symbolic phase,
where a depth-first clustering beats the raw elimination order and nothing downstream would do
better. But it is NOT OPTIONAL, so every caller who does have a symbolic phase pays for a result
they discard. In a decomposed library that would be a flag or a separate entry point. It is a design
cost rather than a bug, and it is the same cost as above wearing different clothes.

**And the defence that a frozen reference implementation pays no maintenance cost is FALSE**, which
matters because it is the argument that would justify copying `AMD_2`'s structure rather than only
its devices.

The premise is that `AMD_2` has been essentially unchanged since the mid-nineties, so the cost of
being unreadable never came due. Two things falsify it.

**Code is re-tuned for architectures, continually.** `W` and `Nv` trade loads for encoding density,
and that trade was made for machines where a cache miss cost far less relative to an instruction
than it does now. This entry's own measurements are a case in point: the binding constraint on
Apple silicon turned out to be back-end stall, not instruction count, which is not the regime those
encodings were tuned for. A codebase that cannot be re-tuned has not avoided the maintenance cost.
It has converted it into a permanent tax that nobody is able to pay down.

**And parallelism is the sharp case.** There are active discussions about multithreading minimum
degree. What blocks that in a monolith is not the data structures, it is the CONVENTIONS: `Nv`
negative for the duration of a step, `W` meaningful only relative to `wflg`, `Pe` flipped to mean a
tree parent. To split work across threads you must know which of those hold over which regions, and
that information exists nowhere but in the head of someone who has understood all 1664 lines at
once. Modular boundaries are precisely what makes such a question answerable, because they are where
you can say what is shared.

So `AMD_2` is correct, fast, and it foreclosed its own future. That it survived thirty years
unchanged is not evidence the trade worked; it is equally consistent with nobody being able to
change it, and the parallelism discussions are what makes those two readings distinguishable. **This
is the reason the shared class is worth its constant factor**, and it is a better-supported claim
than the readability argument alone, which is a matter of taste where this is a matter of what can
be built next.

**Two method notes that cost a day between them and are not about amd at all.**

**Counting and profiling answer different questions, and counting misled twice.** It is the right
instrument for "are we doing more work" and it cannot see a loop whose cost is decided by an early
exit: a counter that added `adjacencySize + incidenceSize` measured what a short-circuiting test
COULD cost, and the real iteration count was 1.7x higher on one side. Nor can it see allocation ,
the largest single item found all day was `operator new` under `beginElimination`, an unreserved
arena doubling 18 times per ordering, and no counter would ever have shown it.

**A call tree bounds the search; only the source view ends it.** Everything in a driver inlines
into one symbol, so the tree could say only that 5.70 s sat in `orderAmd3`'s self weight. The
source view named the line. And when the answer is a CALL rather than a line, it takes a second
zoom: call tree, then line 250, then `allocate.h`. Stopping at either of the first two would have
found nothing.

---

## 2026-08-08: amd3 aligns to the vendored AMD, and the alignment found a defect worth 3 to 9 percent of fill

**`amd3` is the amd counterpart of `mmd3`: `amd2` with the vendored routine's list order, adding
no mechanism.** It returns `AMD_2`'s permutation exactly up to the postorder, on the seven
examples and on every square grid tested from 3 a side to 40, member order within each
supervariable included. The old `amd3`, which carried dense rows, `amd_aat`, the postorder and
the Control interface, is renamed `amd4` and is temporary; the digit now means the same thing on
both branches, 3 is the layer aligned to the vendored code.

Five ledger entries. Four are tie-break conventions, all of the same family: AMD pushes at the
head where we append. The hash bucket is built by head-push while scan 2 walks `Lme` forward, so
its chain is reversed against `C[p]`; the reachable set is laid out cliques-then-explicit, since
`for (knt1 = 1; knt1 <= elenme + 1; knt1++)` takes the elements first and the supervariables on
its last pass; mass elimination runs in scan 2 *after* aggressive absorption, which is what makes
the cheap structural test agree with the true one, as `AMD_2`'s own comment says; and the new
element goes to the front of a variable's list by a rotation, `Iw[p1] = me` with two boundary
entries lifted to the two ends rather than everything shifted.

**The fifth was a defect in `amd2`, `Amd2` and `Amd2B`, and it is the entry that matters.** A hash
merge folds `v` into a live `u` and grows `u`'s weight, and the bound written moments earlier has
that weight subtracted inside it. So a supervariable was filed one bucket too high per original
vertex absorbed. `AMD_2` subtracts `nvi` in the pass that restores the degree lists, which runs
after supervariable detection, so its weight is the post-merge one.

```
grid        AMD (vendored)    AMD1      AMD2 before    AMD2 after
 32x32            11900      12074         12364         11900
100x100          206332     201856        212496        199386
140x140          474995     455472        487111        444191
```

**Three things follow, and the second is the reason this entry exists.**

The comment at all four sites argued for the behavior it got wrong, and argued nearly correctly:
an external degree does exclude `u`'s own supervariable, and `u`'s reachable set really is
unchanged. What it missed is that the buckets are keyed on a degree that has the weight
subtracted *in* it, so the term moves even though the reach does not. A comment that reasons and
concludes wrongly is more dangerous than none, because it answers the question a reader would
otherwise ask.

**It was found by making the code comparable, not by making it fast or by reading it.** This is
the second time: `Mmd2` filed a supervariable one bucket too high for months, worth 13 percent of
fill, and alignment found it in an afternoon. Both defects are the same shape, a weight
subtracted before the merge that grows it, in two different arrays. This file previously recorded
that the amd branch could not carry that shape, because AMD files at an external degree that does
not move when a weight changes; the external degree does not, and the `- nvi` term does.

**And alignment costs fill here where it bought fill on the mmd branch.** `amd3` is aligned, so
its fill is the vendored routine's, 474995 at 140 a side, where the corrected `Amd2` reaches
444191. On grids our tie-break now beats AMD's by 6.5 percent, the opposite direction from mmd.
That is a real number on one problem family, and the family known to flatter us; it is the second
data point for the LIFO-against-FIFO question parked in `experiments/ordering/REPORT.md`, and it
should not be quoted without the 3D grids `REPORT.md` has been asking for.

**What this obliges elsewhere.** Every AMD2 fill figure in `benchmarks/ordering/README.md`,
`benchmarks/pipeline/README.md`, `experiments/ordering/README.md` and `REPORT.md` predates the
fix; each carries a superseding note rather than being rewritten, since a dated measurement is a
record of a run. REPORT finding 3, that AMD2's extras are a net loss with the hash almost all of
it, reverses on fill and stands on time.

---

## 2026-08-08: The integer model, one dimension against two, and what the convenience costs

**Oblio splits integers by DIMENSION, not by role, and that is the whole design.** A quantity that
counts along one side of the matrix is bounded by `n`; a quantity that offsets into an `n x n`
object is bounded by an area. The two have different fates as problems grow, so they get different
treatment:

```
one dimensional   indices and counts   bounded by n, held below 2^31
two dimensional   positions            bounded by an area, held in 64 bits
```

The bound on the first is `2^31` rather than `2^32` because an index may be `NIL`, and giving up
the sign bit costs one power of two. That bound is then enforced at the door, on `n` and on
`nnz(A)`, and never has to be revisited.

**Why the split is at that line and not another.** One-dimensional sizes above `2^32` are not
reachable in any foreseeable machine: even `2^31` indices squared is `2^62` entries, about `5.5e10`
GB. Two-dimensional sizes above `2^32` are ROUTINE. A real factor crosses `2^31` entries at 25.8
GB and `2^32` at 51.5 GB, so both are reached on an ordinary large workstation, never mind a
cluster. So one dimension can be bounded safely and two cannot, and a model that treats them alike
must be wrong about one of them.

**How the established libraries do it, and where it leaves them.** SuiteSparse and MUMPS carry one
source with a width parameter, `using Int = int32_t` in SuiteSparse's case, and the caller selects
32 or 64 at build time. Both dimensions move together. For `A` that works: `nnz(A)` is an INPUT, so
a 32-bit build caps it in the constructor and refuses anything larger, which is a clean contract.
For `L` it does not, because `nnz(L)` is an OUTPUT. A static factorization at least predicts it in
the symbolic phase, so one test before allocation would serve. Under PIVOTING it is not
predictable at all: a delayed column grows its parent's front, so the bound would have to be
re-checked every time a column is delayed, in the most delicate code in the numeric factorization,
and it would fire after real work had been done on an input that was accepted.

**That leaves a fuzzy band, and the band is not exotic.** Between a factor of about 26 GB and about
52 GB, `nnz(L)` has crossed what a 32-bit build can address while `n` is nowhere near any limit.
Precisely the range where a user would reasonably choose the 32-bit build for speed, and precisely
where that choice is wrong. Oblio never enters that band, because two-dimensional quantities are
always 64 bits, so there is no configuration for a caller to get wrong and no failure that arrives
mid-factorization.

**And Oblio goes one step further than the split requires: one-dimensional sizes are `std::size_t`
too.** Only the BOUND is one-dimensional; the TYPE is wide. That is convenience rather than
necessity, and it has a price and a benefit, both now measured.

**The price is width in the hot loops.** In `QuotientGraph` six arrays are `std::size_t` where the
vendored genmmd's equivalents are `int`. Measured two ways in `experiments/ordering/REPORT.md`:
widening genmmd to `int64_t` costs it 17 to 26 percent doing byte-for-byte identical work, and
narrowing our four one-dimensional counts recovers most of it, 1.159x of genmmd against 1.373x at
140x140. Roughly a fifth of the ordering's runtime, paid always.

**The benefit is that the 1D-to-2D crossing cannot be got wrong.** `SymFactor` computes a
supernode's storage as `f * (f + 1) / 2 + f * u`, two one-dimensional counts multiplied into a
two-dimensional size. With `f` as a 32-bit type:

```
front f = 46341    f*(f+1)/2 = 1.07e9   fits
front f = 65536    f*(f+1)/2 = 2.15e9   OVERFLOWS
```

A dense root of 65536 is ordinary for a large 3D problem, and the overflow is silent, in an
arithmetic expression rather than at any boundary anyone would guard, producing a plausible small
number instead of a fault. Narrow the one-dimensional sizes and every such expression needs a cast
that a reader has to recognize as load-bearing. Keeping them wide makes the intermediate already
wide, so the crossing is structurally safe.

**The decision: stay wide, deliberately.** Consistency and a crossing that cannot silently
overflow, against roughly a fifth of the ordering's time. The alternative is available and is
written down rather than forgotten: narrow the one-dimensional COUNTS to `std::uint32_t`, which is
what the type rules are missing a category for, keeping every position at `std::size_t`. That is
most of the measured gain and none of the consistency cost.

**And the reason it is not being done now is not only the effort, though the effort is real.** It
would touch many places, since one-dimensional counts feed arithmetic all through the analysis and
the numeric phases, and every site where a count is combined into a two-dimensional quantity would
need a widening cast. **Each of those casts is a place to introduce the very overflow the wide
types currently make impossible**, and one that is missed is silent. So the work is not only
extensive, it spends effort adding hazards to code that is correct today, in exchange for time.
That trade is worth taking eventually, with the arithmetic audited site by site, and it is not
worth taking in a hurry. Correct code with a measured performance tax beats slightly faster code
with a correctness question in it.

**The general form of this decision, which outlives the particular question.** Every exception to a
uniform rule buys computer time and spends CODER time, and the two are not interchangeable here.
The tax is paid once per run by a machine; the exception is paid again by every reader on every
visit, forever. On a decade-scale project with one author, coder time is the binding resource and
machine time is the abundant one, so a design that spends the scarce resource to save the abundant
one is simply wrong for THIS project, whatever it would be for a team of ten shipping to a
supercomputer centre.

The arithmetic supports it rather than merely excusing it. The width tax is about a fifth of the
ordering, and the ordering is a fraction of a solve: at 140x140 the pipeline benchmark puts MMD3's
ordering at 1.55 ms against a left-looking factorization of 7.01, so a fifth of the ordering is
roughly four percent of an analyze-plus-factor, and less again whenever the factor is reused. That
is what a singularity would buy.

And the project's own history says which failures are expensive. The defect that cost the most was
not slow code: `Mmd2` filed a supervariable one bucket too high for months, worth 13 percent of
fill, MORE than the width tax, and it survived because a subtraction sat somewhere nobody
re-derived. It was found by making the code COMPARABLE to a reference, not by making it fast.

**One asymmetry worth keeping in view: the tax is reversible and the exception is not.** A uniform
`std::size_t` can become a uniform `std::uint32_t` for counts later, mechanically, precisely
BECAUSE it is uniform. A singularity introduced now is load-bearing forever. The convenient choice
keeps the option; the fast one spends it.

**What is NOT claimed here.** That the others are wrong. Their tax is paid only on large problems
and ours is paid always, which is a real advantage for them on small ones. What they buy it with is
a band in the middle where the contract is unclear, and a configuration choice the caller has to
make correctly without being able to predict the quantity it depends on. Oblio pays a known,
uniform, measured tax to have neither.

---

## 2026-08-07: The default ordering moves to MMD3, which reproduces genmmd exactly

**MMD3 is mmd2 with genmmd's list order and one defect fixed, and it returns the vendored
routine's permutation EXACTLY** on all seven examples and every square grid tested from 5 a side
to 80. Six alignments got it there, four tie-break conventions, one cosmetic numbering, and one
real defect; `experiments/ordering`'s README carries the full account and the ledger.

**The default moves in `OrderEngine::mOrdering` and the `DirectSolver` constructor**, MMD2 to
MMD3.

**And the reason is not that it measured best, which is worth stating because the numbers invite
the wrong conclusion.** MMD3's fill gap against genmmd is zero at every size, but that is
guaranteed rather than earned: it IS genmmd's permutation, so identical fill is a tautology. The
comparison that carries information is against the fixed MMD2, and it does not favor MMD3
uniformly:

```
grid          n      MMD2 fill   MMD3 fill
32x32      1024        -0.5%       0.0%
100x100   10000        +5.9%       0.0%
200x200   40000        +6.8%       0.0%
400x400  160000        +8.3%       0.0%
```

MMD2 is very slightly BETTER at 32 a side and a few percent worse above it. So on the evidence we
have, which is grids only, the two are close and MMD3 is ahead on most of the range.

**The argument for the change is about the cases nobody has run.** A default is a bet on unseen
inputs. MMD3 reproduces a reference implementation that has been in use for decades and exercised
on far more than our grids, so its behavior on an unfamiliar matrix is whatever that history has
established. MMD2's tie-break is ours and has been measured on square grid Laplacians and seven
tiny examples. Reproducing the reference is the better bet, and it costs nothing measurable here.

**A second reason, weaker but real.** With MMD3 as the default, any future divergence from genmmd
shows up as a permutation difference in `make test` rather than as a fill number somebody has to
judge. The default is then also the alignment check.

**The defect was fixed in MMD2 as well**, since a defect found in one place is a defect wherever
the code sits. It filed a supervariable one bucket too high per vertex merged into it, by
subtracting the vertex's own weight before a walk that could increase it, where genmmd subtracts
after. MMD2's own gap against genmmd fell from about 20 percent to about 7. `Mmd1` cannot have it,
having no q2h path and no live merges, and the AMD drivers cannot, since AMD files at an external
degree that excludes a vertex's own supervariable.

**What did not change.** `Ordering::MMD2` is still there and still selectable, the examples still
name their ordering explicitly rather than relying on the default, and the vendored `MMD` and
`AMD` remain optional, so the default is still an ordering that is present in every build.

---

## 2026-08-04: The vendored orderings move to private/ and become optional

`src/Amd.cpp` (SuiteSparse AMD 3.3.4) and `src/Mmd.cpp` (Sparspak genmmd, via 0.9) are not our code
and are not ours to publish. They now live in `private/`, which is gitignored, and both builds
detect that directory rather than requiring it.

**Detection, not configuration.** The Makefile takes `$(wildcard private/Amd.cpp private/Mmd.cpp)`
and defines `OBLIO_VENDORED_ORDERINGS` when it is non-empty; CMake does the same with `file(GLOB)`
and reports which mode it configured. There is no flag to pass: a tree that has the directory
behaves as it always has, and a clone behaves without it. The alternative, an explicit option,
would have to be set correctly by whoever builds, and the one person who needs it set is the one
person who would never see the failure.

**One override, `OBLIO_PUBLIC=1`, for building the way everyone else does.** It empties
`VENDOR_SRCS`, so `private/` is ignored even when present. It is spelled the same in the root
Makefile, both benchmark Makefiles and `experiments/ordering`, which is the point: the person who
needs it is the only person who has the directory, and one word in front of any make command in any
of those directories is a habit, where a different target name per directory is a lookup. Each of
the four keeps a `.build-mode` stamp so the two interleave without a clean, rebuilding on a switch
and not on a repeat, and neither discards what the other built.

**What absence costs.** In `OrderEngine` the two `extern` declarations, the two method declarations,
the two definitions and the two `switch` cases sit behind the guard, and the cases return `false`
instead. Every enumerator is still named, so `-Wswitch` still protects the enum against a new
ordering being forgotten. `Ordering::MMD` and `Ordering::AMD` remain in the enum and refuse, which
is what the phase's `bool` is for.

**The default moves to `Ordering::MMD2`**, in `OrderEngine::mOrdering` and the `DirectSolver`
constructor. (Superseded 2026-08-07: the default is now `Ordering::MMD3`. The reasoning below
still applies, since MMD3 is also always present.)
A default that refuses on a fresh clone is the first thing a newcomer meets and gives
no hint why. MMD2 is the closest of ours to genmmd, which is what the default was, so the change is
within the family rather than across it. The examples name `MMD2` explicitly rather than relying on
the new default, so a reader sees a choice being made.

**One copy, not two.** `experiments/ordering/vendored/` held byte-identical duplicates of both
files; those are deleted and `experiments/ordering/Makefile` compiles straight from
`../../private/`, reporting "skipped" rather than failing when it is absent.

**Test counts now depend on the build**, 252 with and 238 without, the difference being fourteen
assertions in `test_order` that check the vendored routines themselves. That is recorded in
`docs/TESTING_SPECIFICATION.md`. Everything else asserts the same thing either way.

**Not settled here: the licenses.** AMD carries a BSD-3-Clause notice, which permits redistribution
provided the copyright, the conditions and the disclaimer are retained; the file has the copyright
line but neither of the other two, so shipping it would need the full text added. MMD carries no
notice at all and its Sparspak provenance is unclear. Neither question is answered by moving the
files, and both remain open.

## 2026-08-01, The B suffix, a seam, and two hypotheses built and falsified

The largest measured item on the ordering list was AMD1's three visits per element against the
vendored routine's one. It was built, it works, it is not faster, and the diagnosis behind it was
wrong. Recording that fully, because a reversed conclusion is worth more than a confirmed one and
because the seam it required is reusable.

**The naming axis, which is new and needs stating once.** A trailing digit means a different
ordering: AMD2 has mechanisms AMD1 lacks, so their permutations and their fill legitimately differ.
A trailing B means the **same ordering computed differently**, so AMD1B must return exactly AMD1's
permutation and any difference is a defect in one of them. Two axes, one enum, and a reader six
months out cannot tell them apart without being told, which is why `OrderEngine.h`, `Amd1B.h` and
`TESTING_SPECIFICATION.md` all say it.

That buys an oracle the digit pairs cannot have. AMD1 against AMD2 can only be compared on fill;
AMD1 against AMD1B must agree entry for entry, on both maps. It is the strongest check in the
ordering suite: 136 graphs and seven assertions, zero mismatches.

**The seam, and why it is not the template it was designed as.** Option C was a
`template <class Scan> eliminate(pivot, Scan&)`, and it collides with the CLAUDE.md invariant that
template definitions live in `.cpp`: the definition must be visible where it is instantiated, so
either it moves to a header or the shared core includes its drivers. Reading the rule's own
rationale settled it rather than leaving it to taste. That entry argues entirely about the `Val`
axis, a closed set of scalar types and the build cost of scattered instantiation, and reaches
nothing about a function templated on a caller-supplied policy type. **The invariant as written is
broader than its own rationale**, which is a genuine gap, and the right response was not to widen
the rule under pressure from a change that did not need it.

So: a non-template overload taking an `ApproximateScan` struct by reference, with `eliminate` split
into private `beginElimination` and `finishElimination` so the half-eliminated state never escapes.
No template, no friendship, one parameter. If a future variant needs a real policy point, the
invariant should be narrowed on its own merits first.

**The result. AMD1B is five percent slower at 140x140** and one to two percent faster at the three
smaller sizes, with identical fill and identical permutations. The visits were never the cost.

**And the second hypothesis failed the same way.** A comparative profile, the first in this
project, put AMD1's gap at 46 percent data stalls, 27 percent work and 17 percent branch
mispredicts. Work being third explains AMD1B. Data stalls being first pointed at footprint, so the
six `std::size_t` arrays in `QuotientGraph` were narrowed to `int32_t` as an experiment: cachegrind
reported D1 misses down 17 percent with instructions flat, and alpamayo reported nothing, with both
vendored controls unmoved. Reverted. **`std::size_t` for a position stands, and now has a
measurement behind it rather than only a rule.**

**What this says about instruments.** Cachegrind was wrong three times in one afternoon in three
directions, and cache simulation on a different machine does not transfer even comparatively, where
instruction counts did. More generally, every one of today's five successful changes came from
reading a trace and fixing the top line, and every failure came from reasoning ahead of one. Two
elaborate hypotheses with good counter evidence were built and falsified in an afternoon; that is
cheap and is the point of the infrastructure, but the ratio is worth remembering before the next
one.

**AMD2B followed, and is the more useful of the pair despite being the easier to write.** The seam
was built, so applying it to AMD2 was driver work: three places it could have differed were checked
rather than assumed, and none does (AMD2's scan 1 is textually AMD1's, absorption runs after the
elimination and touches only incidence lists, and the hash merges run after the bound). It is worth
having because AMD2 carries two mechanisms AMD1 does not, so `AMD2B == AMD2` guards an absorption
pass and a hash comparison that nothing else in the suite checks entry for entry.

**And the pair shows a pattern neither alone would establish.** Both are faster at small n, about 20
percent at 32x32, and both lose that advantage as n grows, ending at plus 4 and plus 1 percent at
140x140. The likely mechanism is that the fusion is not free of memory: it adds `explicitPart`, an
array of size n, 156 KB at n = 19600, which is another stream competing with the dozen already
there. **The fusion trades element visits for footprint, which is the wrong direction for a branch
whose largest gap component is data stalls.** Hypothesis, not finding; testing it would mean a third
variant nobody wants.

**Whether they stay.** The collapse condition was written before either was built: a B variant
replaces its original when permutation-identical and faster. Neither is faster, so neither fires and
the tree carries two implementations of each of two orderings. Kept anyway, as measured negatives in
the shape `experiments/` keeps its rejected ideas, and because the seam and the oracles are worth
more than the duplication costs. Revisit when someone needs those files to be simpler.

---

## 2026-08-01, Fine-grained allocation in the ordering is closed, and all four sites had one cause

Two last sites, both small: AMD2's hash buckets were a vector per bucket over n + 1, constructed and
destroyed per ordering, and `QuotientGraph::eliminate` returned its `merged` list by value. The
buckets are now `hashHead` and `hashNext`, the idiom `Buckets` already carries and `Amd.cpp` uses
for the same job, and `merged` is a member scratch returned by const reference. Allocations per
ordering at 140x140 went from 31915, 29499, 32256 and 46351 to 70, 105, 62 and 56. Measured on
alpamayo at about 4 percent for AMD2 and 1.5 for the two MMDs once a 2 percent drift is removed, and
nothing for AMD1, which had already stopped being allocation-bound.

**The pattern is the entry, not the two sites.** Four things were fixed today and all four had the
same origin: `I[u]` as a vector per vertex, `std::vector<bool>` for the flag arrays, `hashGroup` as a
vector per bucket, `merged` returned by value. **Every one is a data structure that crossed from the
prototypes into production along with the logic**, and the prototypes are right to hold all four that
way, since their job is to read as the algorithm and they order one graph per run. The production
drivers were extracted from them and kept the shapes past the conditions that justified them.

So the lesson is not that we allocate too much. It is that **extracting a driver from a teaching
implementation carries its containers as well as its algorithm**, and the containers are exactly the
part that should not survive the move. Worth one deliberate pass the next time a prototype becomes
production, asking of each container whether the algorithm needs it or the reader did. All four of
these were found one accident at a time instead.

**Two hazards that came out of the same day and generalize.** A container swapped for a chain can be
a tie-break change in disguise: a head-pushed hash chain reverses each bucket, and the bucket order
decides which of two indistinguishable vertices absorbs the other, so the first version moved four
permutations. Filling in reverse restores them, and the only thing that caught it was diffing 76
orderings against the previous tree. And a cache simulation on another machine over-predicted once
and under-predicted once in the same afternoon, so it is worth having to decide whether a change is
worth trying on the target, and worth nothing as an estimate of what it will buy.

**What is not closed.** Coarse allocation remains and is fine, being a handful of whole arrays per
ordering, which is what the forest and the symbolic factorization have always had. Memory is not
closed either: the `std::vector<bool>` result was about initialization rather than allocation, and
the constructor's first touch is still 3 to 4 percent of MMD2 and is not reachable by allocating
less.

---

## 2026-08-01, Two micro-optimizations in the ordering, and both mechanisms were mis-guessed first

Follow-ons to the shared-run entry below, both found by reading a Time Profiler trace rather than by
reasoning, and both worth recording for the same reason: the change was right and the explanation we
brought to it was wrong, twice.

**Hoisting the loop bounds.** `for (i = 0; i < qg.incidenceSize(u); ++i)` re-loads the bound every
iteration, because the loop bodies store through arrays the compiler cannot prove disjoint from the
size vector. The accessor showed as 300 ms of self time in an inlined one-liner. Hoisted, AMD1 fell
6.31 s to 6.09 s.

The refinement is the part worth keeping: **hoist where the loop is long, leave it where the loop is
short or exits early.** Hoisting everywhere made AMD2 2.3 percent worse and MMD2 3.2 percent worse,
in loops that exit on the first mismatch or that run over at most two elements by construction.
Eight sites decline the hoist and say so.

Also: `g++` at `-O3` performs this hoist itself, so a Linux measurement reports zero. Apple Clang
does not. **A negative result on one toolchain is not a result about another**, which the gprof
entry in `benchmarks/README.md` already said about profilers and is now said about optimizers.

**The boolean flags off `std::vector<bool>`.** Three arrays, `QuotientGraph::mEliminated`,
`Buckets::mFiled` and `Mmd2`'s `outmatched`, moved to `std::vector<std::uint8_t>`. Measured on
alpamayo at 140x140, the four branches gained 8 to 12 percent against controls that drifted 3, so
call it 5 to 9.

**The mechanism is construction, not access**, and the obvious story is wrong. A `std::vector<bool>`
of n false is built bit by bit through the proxy machinery; a byte vector of n zeros is a `memset`.
Two agreeing traces put the graph constructor at 710 ms before and about 186 ms after, with the flag
reads moving by less than the noise. Both arrays are constructed once per ordering, so the saving is
per construction and nearly uniform across branches, which is exactly the pattern the timings showed
and which the access story cannot explain.

The spelling is `std::uint8_t` rather than `char`, since the tree pins fixed-width types from
`<cstdint>` and writes no bare fundamental type for storage. There was no precedent to follow: the
only other small-value array in the tree is `NumFactorDynamic::mPivotType`, which is a
`std::vector<std::int32_t>`.

Three `std::vector<bool>` remain, in `ElmForestEngine` and `Permutation`, and stay: none is
constructed per ordering or read per element. The prototypes in `experiments/ordering` keep the
packed form too, and correctly, since they construct one graph per run.

**One number elsewhere is now suspect.** The 2026-07-31 entry credits `mLiveMerges` with about a
quarter of MMD1's time for skipping the `mEliminated[v]` load. That was measured against the packed
array, and the read cost has since measured as noise, so the flag is probably worth much less than
recorded. It still skips a load and a branch and should stay; the number should not be re-quoted
without re-measuring.

**And the two passes bought efficiency where the shared run bought work.** CPU Counters on AMD1 at
three points: 29.50 G cycles at 49.02 percent useful on 2026-07-31, 25.93 G at 43.27 after the
shared run, 22.49 G at 47.05 after these two. The dip in the middle was deleting allocator
bookkeeping, which is high-IPC work, so the remainder was a harder mixture; these two put 3.8 points
back. A timing table cannot tell those apart and this instrument can, which is the whole argument
for it in `benchmarks/README.md`.

The consequence is a work item being repriced. AMD1 against the vendored routine went from 1.95x
cycles and 1.68x useful cycles to 1.49x and 1.23x, so its remaining gap is now roughly half work and
half stalling where it was 86 percent work. The driver restructuring removes two of three visits per
element and so attacks the work half alone: its ceiling is now about a fifth of AMD1 rather than the
two fifths it looked like that morning. Still the largest item on the list, and no longer the whole
answer.

**And the instrument is spent.** Four of the last five decisions landed in the 1 to 4 percent band
against 1.7 to 4 percent of drift between runs. Below roughly 5 percent, `make run` cannot decide
anything, and the Time Profiler line is what settled both changes above.

---

## 2026-08-01, A vertex's two lists share one block, because their sum is conserved

The ordering's last per-list allocation was `I[u]`, one `std::vector` per vertex, and a Time
Profiler trace put its `push_back` in `eliminate` at 1.12 s of `amd1`'s 7.08 s, the largest single
line in the program. The change that was written down for it was an arena, on the reading that
`A[u]` only shrinks and `C[c]` is written once while `I[u]` genuinely grows.

**`I[u]` does not grow, once it is not looked at alone.** `A[u]` and `I[u]` are the two kinds of
*source* that `reach(u)` is a union over, one per explicit neighbor and one per clique, and each
elimination that reaches `u` replaces at least one source with the new clique rather than adding
one. Where `u` was reached through `A[pivot]`, the prune drops `pivot` from `A[u]`; where it was
reached through a clique, that clique is absorbed and leaves `I[u]`. So

```
|A[u]| + |I[u]| <= the off-diagonal entries in u's column of A
```

for the whole run, the pair fits in one block sized once from the pattern, and nothing grows,
moves or is reclaimed. Every later mechanism preserves it: mass elimination, live merging and
numbering all empty a list or destroy a clique, and none creates a source.

**Both vendored routines already do it, differently, and neither says why.** They differ on two
separate choices, which the first version of this entry merged. The LAYOUT: MMD puts variables
first and elements last, as we do, and AMD puts elements first. The BOOKKEEPING: AMD records the
split in `Elen`, we record it in `mAdjacencySize`, and MMD records no split at all, recovering the
kind per entry from `invp[nb] < 0`, the sign of the inverse permutation it must produce anyway.

**Neither choice is forced by its algorithm**, which was the other thing the first version implied.
MMD could keep a boundary and AMD could classify per entry; AMD even carries a liveness flag
already, `Elen[e] < EMPTY`, and declines to use it that way. We took the boundary on the counting,
one read per list against one per entry, and on the precedent that `mLiveMerges` exists precisely
to avoid a per-entry liveness load. **Which is actually faster here is not measured**, and MMD's
loads land in an array the degree lists have already threaded, so they may be warm. Section 5.15
carries the full comparison.

The order within our run is forced rather than chosen: the prune compacts the adjacency and then
the incidence in two passes, and only adjacency-first keeps the write cursor behind the read one.

**The correction this forced is worth more than the change.** Section 5.15 of
`archive/sparse_factorization.md` had filed the whole of `Iw` under archaeology, alongside the
compaction, on the grounds that a flat workspace is what a language without allocation forces. That
is right about the pool and wrong about the run. `Iw` answers two questions at once and only one is
about Fortran: where a vertex's sources live, which the bound settles for any language, and where
the element patterns live, which genuinely grow and are what the pool, `pfree` and `ncmpa` exist
for. The lemma is now 5.3 and the reading of the two codes is 5.15.

**A general hazard, and the reason this sat unclaimed.** A flat workspace forces the question at
the first line of code, since a run that can be outgrown needs a copy and a compaction, so both
authors had to answer it before writing anything. A `std::vector` grows silently, so the question
their encoding put first, ours never raised. That is a real cost of the better encoding rather than
an anecdote about `I[u]`, and it is worth expecting again.

**Measured, alpamayo, 2026-08-01.** Allocations at 140x140 fell from 31915 to 2457 (MMD1) and
32256 to 2760 (AMD1); time fell 6 to 18 percent, MMD2 reaching 1.58x the vendored MMD, which is the
closest any of ours has come and the first time the MMD branch has passed the AMD branch. Every
permutation is bit-identical to before, checked over 76 orderings, and every fill figure is
unchanged on two machines. `benchmarks/ordering/README.md` carries the numbers.

**And the two tables disagree again**, 92 percent fewer allocations for 12 percent less time on
AMD1, which is the third time in this work that allocation counts have not predicted time. What
paid was removing a line from the innermost loop, not removing mallocs, and even that came in under
what the trace attributed to it. AMD's remaining gap is work rather than memory, as the cycle
profile said before this change and says after it, so the driver restructuring is unaffected.

**What is retired without having been built.** `std::pmr` over a monotonic buffer, recommended in
`experiments/ordering`'s garbage-collection section and in TODO, as the answer to per-list
allocation in the ordering. There is no longer a per-list allocation to redirect.

---

## 2026-08-01, MMD2 needed less of the core than expected, and one sizing bug it exposed

MMD2 is MMD1 plus the four things genmmd does around the batch: the prepass, the filing
convention, the element-by-element q2h refresh, and pairwise merging with outmatched marking. A
seventh ordering method, beside MMD1 rather than replacing it.

**The eliminator did not change after all.** The expectation since the first hour of this work was
that `mmd2_eliminate`'s `bwd[rn] = 0`, clearing the outmatched flag for everything a new clique
reaches, would force a second eliminator or a flag through the shared one. It did not: nothing
inside the eliminator reads that flag, and the driver already walks the new clique immediately
afterwards to evict it from the buckets. So the clear lives in the driver, on a flag the driver
owns, and the shared core kept its shape.

What the core did need was smaller: a weighted reachable size, since pairwise merging makes a LIVE
vertex stand for several originals; a `number` operation for the prepass, which marks a vertex
without forming a clique or pruning anything; and the prune skipping vertices numbered that way,
which linger in the adjacency by design.

**A one-line sizing bug the prepass exposed.** `Buckets` indexed its heads by degree and its links
by vertex from one size. Every branch until now filed a vertex under its degree, at most n - 1, so
one array served both. MMD files under degree plus one, which reaches n, and on a one-vertex matrix
that is an out-of-bounds write in the constructor's first filing. Caught by `test_order`'s n = 1
tridiagonal, which is exactly the edge case the specification says to watch, and which the
experiment's seven graphs cannot reach.

**It is the one prediction of the day that paid.** The cycle profile had said MMD1's gap was work
rather than efficiency, 42.2 percent useful cycles against the vendored 43.0, so the mechanisms
were where the time was. MMD2 came in 30 percent faster than MMD1 at 140x140, taking the MMD branch
from about 2.6x the vendored routine to 1.8x. Fill moves either way with the tie-breaks and is not
the point.

---

## 2026-07-31, AMD2 takes two of AMD's mechanisms and not its postorder, and live merging changes the core

AMD1 is the bound alone. AMD2 adds the two mechanisms that make the vendored routine what it is,
both acting on the reached set the bound was just computed over: aggressive absorption, which kills
a clique that scan 1 shows lies wholly inside the new one, and hash supervariable detection, which
finds vertices indistinguishable from each other rather than from the pivot. It is a sixth ordering
method beside the others, not a replacement for AMD1, so any difference between them stays
measurable.

**The postorder is not ported, and that costs the exact harness check.** `Amd.cpp` and the `amd2`
prototype both end by returning a postorder of the assembly tree they build while eliminating.
`ElmForestEngine` builds its own forest from the permuted matrix and orders it, so this would be
work done twice and then overwritten. A postorder relabels an elimination tree without changing it,
so skipping it moves the permutation and not the fill, which means production AMD2 cannot be
diffed against its prototype by permutation. The harness therefore checks it by nnz(L) instead,
through Oblio's own symbolic factorization, and the prototype's `matrix1` example is excluded
because it exists to exercise the input conditioning we deliberately do not have.

**Live merging is the first thing to change the shared core rather than a driver.** Everything
before it merged only into the pivot, which is eliminated in the same breath, so no eliminated
vertex could linger in a live clique. A hash merge folds one live vertex into another and leaves it
where it lies at weight zero, since purging it from every clique that names it would cost a pass
over the structure per merge and buy nothing, the merge test having established that it is
redundant wherever it appears. `QuotientGraph::merge` is that operation, and the reachable set now
has to skip dead vertices, which is what `amd2Neighbors` does and `amd1Neighbors` does not.

**And that filter had to be made conditional, on measurement.** Adding it unconditionally cost MMD1
about a quarter of its time, the exact refresh being its hot loop and the check being a load per
element for a case that branch can never produce. A flag set by the first live merge, hoisted into
a local, short-circuits the condition so the load never happens until a merge has actually
occurred. Same correctness, no cost to the branches that cannot need it.

**What AMD2 buys, measured on grids: not obviously anything yet.** Its fill is about 5 percent worse
than AMD1's (487111 against 455472 at 140x140), which is the coarser-supervariable cost the
experiment already measured on small graphs, and it runs about twice AMD1's time, the hash pass and
the filter both being new work. It refreshes fewer degrees, which is what it was meant to buy, and
on a grid that does not pay for itself. Whether it pays on a problem with real supernodal structure
is exactly the question the thin test set cannot answer, and is now the strongest argument for
widening it.

---

## 2026-07-31, Oblio's own orderings sit over a shared quotient graph, and the drivers return an order

The ordering experiment has run far enough that the prototypes can start becoming production code.
The plan is two parallel orderings beside the vendored pair, `MMD1` and `AMD1`, built from the
`mmd1` and `amd1` layers, kept aligned with those prototypes while both continue to develop, and
deprecating the vendored routines only once they can replace them. This entry records the shape
those first two steps settled on. Both have landed, MMD1 first and AMD1 over the same seam, which
needed three accessors on the quotient graph that an exact degree has no use for (the adjacency,
the incidence list and the supervariable weight) and a refile on the buckets, and nothing else.

**One core, two drivers.** `Cliques`, `Buckets`, `Neighbors`, `Eliminate` and `Refile` are
byte-identical between `mmd1.cpp` and `amd1.cpp`, the duplication being deliberate in an experiment
where each layer must read standalone. That reason does not survive the move, so production has one
`QuotientGraph` unit holding the adjacency, incidence and clique lists, the degree buckets, the
reachable-set query, the eliminator and the CSC-to-adjacency builder, with the drivers holding only
what the two algorithms actually disagree about: a batch loop with eviction against one pivot per
step with the bound.

**The eliminator stays a plain member rather than a policy point.** `mmd2` changes it, clearing the
outmatched flag for everything a new clique reaches, which is `mmdelm`'s `bwd[rn] = 0` and the first
change to that function since md2. A second eliminator beside the first will be the right answer
when it arrives. What is declined now is a parameter or a virtual on the strength of a divergence we
can describe but have not written: one extra function beside its sibling is cheaper than a seam
guessed at in advance.

**The write grant stays with the engine.** `Permutation` declares `friend class OrderEngine`, and
the rule elsewhere in this file is that an object is filled by exactly one engine. So the drivers
return a `std::vector<std::int32_t>` elimination order and `OrderEngine::orderMMD1` writes the two
maps from it, exactly as `orderAMD` already does with the vendored `P`. Making the drivers
first-class fillers would have meant granting friendship to functions that are not engines, and
would have cost the property that a driver can be called without a `Permutation` in hand, which is
what lets the experiment check production against the prototype directly.

**The dispatch became a switch, and it had to.** `OrderEngine::compute` was an if-chain ending in an
unguarded fall-through to MMD, so a new enumerator would have become MMD silently and produced a
valid permutation that nothing in the suite would have questioned. It is now a switch naming every
enumerator with no `default`, which is the house rule and is what makes the next enumerator safe.
The diagonal-stripping rebuild moved down into `orderMMD` with it: that is the vendored routine's
requirement alone, AMD symmetrizing internally and MMD1 skipping `i == j` while building its own
adjacency, and it was in the shared path only because MMD happened to be the fall-through.

**What did not come across.** The printing, the counters, and `nnzL`. The last is not instrumentation
in the way the others are, being the column-count identity computed from quantities the elimination
already holds, but nothing downstream would read it: `SymFactorEngine` computes the exact structure
later and never consults an ordering's estimate, which is the same argument the experiment's README
makes against AMD's `Info` block. It is a few adds to restore if a use ever appears.

**One defect did not come across.** `amd1` allocates and zeroes its per-clique array inside the
pivot loop, which is O(n) per step and O(n * n) over the run in bookkeeping alone, independent of
the graph, and would have swamped the very cost the bound exists to save. It is hoisted in
production and in the prototype alike, the step clearing only what it wrote, with the prototype's
trace unchanged to the character. The Python twin never had it, its own array being a dict over the
cliques the step touched, which is already the right shape.

**What these two are not.** Each carries its idea and not the whole of its vendored counterpart.
MMD1 has no prepass, no `q2h` merging of vertices indistinguishable from each other, no outmatched
marking, and our filing convention rather than MMD's. AMD1 is the bound and nothing else: no
aggressive absorption, no hash detection, no dense-row removal, no postorder. So they are different
orderings rather than reimplementations, and nothing asserts that they agree with anything.

Measured on grid Laplacians through the same pipeline, nnz(L) lands within about two percent of the
vendored routines throughout (at 32 by 32: MMD1 11972 against MMD's 11822, AMD1 12074 against AMD's
11900), which is the two-sided noise the experiment predicts for a heuristic whose ties fall
differently. Ordering time is five to ten times theirs, on per-list allocation rather than on any
asymptotic difference: across a 20-fold size range both grow very close to linearly in n. Those
numbers are the honest starting point rather than a target, and the allocation question is the one
the experiment's garbage-collection section already anticipated, keep the containers and change the
allocator.

---

## 2026-07-31, Project roots are a second map over the tree, laid after the source layout

Oblio is a C++ tree with a Python half in one experiment, and the two are read in different IDEs,
CLion and PyCharm. Placing them raises a question that looks like a layout question and is not one,
and getting it right depends on separating two things that are easy to conflate.

**The constraint is one IDE per project root, and that is the whole of it.** A JetBrains IDE claims
a directory by writing `.idea/` into it, so two IDEs opening the same directory contend over
`misc.xml` and the module type, and one will quietly reconfigure the other. The rule says nothing
about languages. Reading it as "one language per folder" is the mistake that follows from it, and it
is the mistake worth naming, because the fix it suggests is a source reorganization that nothing
actually requires.

**So there are two independent maps over one tree.** The source layout answers to the code: what
belongs beside what, what a reader needs in view at once, what a build treats as a unit. The root
layout answers to the tooling: which directories an IDE is pointed at. Choosing them together, so
that every language sits in its own folder and every folder is one IDE's root, collapses the two
into one map and lets the weaker consideration win. The discipline is to lay the sources out for the
code's reasons and then place roots over whatever layout that produced.

**An IDE project is not a language project, and conflating them is how the second map acquires
weight it should not have.** An IDE project is editor state: which folder is the root, which
interpreter, run configurations, code style. It is `.idea/`, it is generated, it is ignored, and no
build or test ever reads it, so deleting it breaks nothing. A language project is a build and
dependency declaration: `CMakeLists.txt`, `pyproject.toml`, `build.sbt`, `Cargo.toml`, or a plain
`Makefile`. It is authored, it is tracked, and it is what actually compiles or resolves imports, so
deleting it breaks everything. One is a view of the code and the other is part of it.

The two coincide often enough to hide the difference. In the data-structures repo they land on the
same four folders, `python/` holding both a `pyproject.toml` and an IDE root, which is exactly the
case that makes them look like one thing. The ordering experiment separates them in both directions
at once: a language project with no `pyproject.toml`, since the `Makefile` is the thing that builds
and runs and the Python needs nothing to resolve, and an IDE project at the same folder that does
nothing but remember where to open files. `compile_commands.json` is the sharpest case of all,
being neither: generated, ignored, machine-specific, and read only by an indexer, a bridge that
lets an IDE understand a build it does not run.

The consequence is the one that governs here. **PyCharm needs an IDE project and does not need a
language project**, because nothing in the experiment imports anything or depends on anything.
Adding a `pyproject.toml` to satisfy an editor would be answering an IDE question with a build
artifact, which is the same category error as splitting a folder to satisfy an editor, one map
reaching into the other.

**Which kind of language project, then, is a second question, and it turns on whether the ecosystem
is needed rather than on which is nicer.** A `Makefile` is a task runner over a generic dependency
engine: it knows nothing about the language, resolves no dependencies, and finds no toolchain, so
everything it does it does because someone wrote the command out. `CMakeLists.txt`,
`pyproject.toml`, `build.sbt` and `Cargo.toml` are declarations that *other tools read*, which is
the whole of what they buy. So reach for the ecosystem when something only it provides is actually
wanted:

- resolving third-party dependencies
- installing or publishing
- being consumed by someone else's build
- feeding an IDE or a CI that reads that format

A `Makefile` suffices when dependencies are absent or vendored, nothing is installed, one toolchain
is invoked directly, and the work is a few compile-and-link commands plus some run-and-check tasks.

**By that test the experiments are the clean `Makefile` case and should stay there.** No
dependencies, nothing installed, nothing consumed by anyone. And `make test` in the ordering
experiment is mostly not compiling at all: it runs nine C++ binaries, runs nine Python scripts,
strips punctuation from both, and diffs the pairs. That is task running, which is what Make is for,
and an ecosystem build would express it worse rather than better. The same holds for every other
experiment, and it is the reason none of them has anything but a `Makefile`.

**Where Oblio proper sits today is both, with the `Makefile` canonical.** The root carries a
`CMakeLists.txt` that does real work, a `find_package` for BLAS and LAPACK, a
`check_c_source_compiles` probe for the Fortran underscore convention, a static library target and
`ctest`, and the `Makefile` beside it that CLAUDE.md names as the single source of truth. That is a
defensible position for a tree nobody outside consumes, since the `Makefile` is what the
one-unit-at-a-time loop actually runs, but it is two descriptions of one build and the cost of
keeping two in step is real and ongoing.

**Where it has to go is CMake primary, and the trigger is already written down under another
heading.** The moment outside researchers consume the library they need `find_package(Oblio)`,
install targets and an exported config, and none of that is a `Makefile`'s to give. That is the
same event that triggers the PolyForm license and filling out CONTRIBUTING.md, so it is one
milestone rather than three, and the right time to decide whether the `Makefile` stays as a
convenience wrapper over `cmake --build` or goes. Until then the ordering stands: the `Makefile`
leads, CMake follows, and the experiments stay out of it entirely.

**Why the ordering experiment is not split by language.** The tempting move is
`experiments/ordering/cpp/` and `experiments/ordering/python/`, which would give each IDE an
uncontested root. It would be wrong here, and the reason is what the unit of work is. In the
data-structures repo the unit is a language: a session is spent writing the Rust versions, and the
four implementations of a problem never need to see each other, so a folder per language is the
honest layout there. In this experiment the unit is a *layer*, and a layer is a pair. `md3.py` and
`md3.cpp` are one artifact expressed twice, and `make test` exists to diff their traces against each
other. A directory boundary through the middle of that separates precisely the two files whose
agreement is the point.

**The placement rule: each IDE's root is the smallest directory containing all of that language's
code, and no directory is the root of two.** Oblio falls out of it without any adjustment, because
the two languages have different reach. C++ is the product and spans the whole tree, so its root is
the repo root. Python exists in exactly one folder, so its root is that folder. Broad root for the
broad language, narrow root for the narrow one, and the two nest rather than collide.

**The axis that decides between nesting and splitting is reach, not importance.** Nesting is
available exactly when two languages have unequal reach, so that one of them has a proper
subdirectory to itself. Splitting is forced exactly when their reach is equal, both spanning the
same directory, since there is then no narrow root to nest inside the broad one. That this
correlates with a major language and a prototyping one is incidental rather than the rule. If Oblio
grew a second production language spanning the whole tree, splitting would be forced at the root
even though C++ would still be the major language; and conversely a prototyping language confined to
one folder nests whatever its standing.

**The data-structures split is not evidence that tooling should drive layout.** Its four languages
all reach the repo root, so the rule above forces a split there, but it is also the honest source
layout on its own merits, for the unit-of-work reason given above. The two considerations agreed,
which is the comfortable case and probably the common one. What this entry warns against is only the
case where they disagree, and splitting `ordering/` is exactly that: the source layout has an
answer, the tooling would prefer a different one, and the tooling is the weaker claim.

**One escape neither repo takes: a single IDE covering both languages.** CLion has a Python plugin
and IntelliJ Ultimate handles both, which dissolves the constraint outright, since with one IDE
there is no second `.idea/` to contend for and roots may sit wherever the sources put them. The
price is
the native experience in each, which is the reason for running two IDEs in the first place, so this
is a real option rather than the obvious one. It is recorded because it is what turns same-level
coexistence from a conflict to route around into a non-problem, and anyone weighing a future split
should know it is on the table.

**Nesting obliges the outer IDE to exclude the inner folder, and ours already did.** CLion's
`.idea/misc.xml` lists `experiments/ordering` under `excludeRoots`, along with every other
experiment and `reference/`, so the two IDEs are not merely avoiding a shared `.idea/`, they are not
indexing the same files at all. That exclusion was there for its own reasons before this question
came up, which is a small piece of corroboration for the arrangement rather than a coincidence: an
experiment is not part of the main build, and the same fact that keeps it out of CLion's index is
what makes it a clean root for something else.

**Nothing about either root is tracked.** The root `.gitignore` carries `.idea/` with no leading
slash, so it matches at any depth and covers the experiment's copy as well. `vcs.xml` in particular
is derivable, the IDE writes it by walking up to find `.git`, so it records what git already knows.
Worth knowing that this cuts against JetBrains' own convention: the `.idea/.gitignore` they ship,
present in our tree, ignores `workspace.xml` and `shelf/` on the assumption that the rest of
`.idea/` is shared by a team. Our outer ignore sits above it and wins, which leaves that file inert.
For a solo pre-release tree that is the right call, and if part of `.idea/` is ever tracked the
candidate is `codeStyles/`, which holds real formatting decisions and is the JetBrains counterpart
to `.clang-format`, not `vcs.xml`.

**What bounds all of this is contiguity, and naming the bound matters more than the arrangement.**
The narrow root works because Python occupies exactly one folder. What would break it is not Python
spreading to more experiments, which the next paragraph covers and which costs nothing, but Python
appearing in places with no common parent short of the repo root: a matrix generator in `tests/`, a
script in `tools/`, the prototypes here. No narrow directory covers that set, and the choice would
then be several roots or one root colliding with CLion's. It is a real fork rather than a
hypothetical one, since a generator for test matrices is exactly the kind of thing this project
would reach for. When it arrives, the answer is to consolidate the Python rather than to split the
experiment, for the same reason the split is wrong today.

**Python in more than one experiment simply moves the root up.** The rule gives the smallest
directory containing all of it, so Python in every experiment puts the root at `experiments/`, and
Python in three experiments out of ten puts it there too, since no one experiment holds all three.
Three separate roots is a legitimate departure from the rule rather than what the rule gives, and
the question deciding it is the one that decided against splitting `ordering/`: whether a unit of
work spans the three. Prototypes that are standalone and never see each other are honestly three
roots, each matching the self-contained-folder convention experiments already follow. Anything
shared between them, a graph generator or a trace-diffing harness, the sort of thing that tends to
appear by the third prototype, is honestly one root at `experiments/`, and the shared code then has
somewhere to live. Under one root the non-Python experiments can be marked excluded, one click
each, so the middle path costs little.

**A root is cheap to move, and that reversibility is the real payoff of keeping the two maps
separate.** Moving one is two steps: open the new folder, delete the old `.idea/`. The IDE writes a
fresh `.idea/` on open without being asked, so nothing is created by hand. No source file moves, no
import breaks, no Makefile rule changes, no README reference goes stale, and git sees nothing at
all, `.idea/` being generated and ignored. The only care needed is the constraint this entry opened
with, that the new root is not already some other IDE's.

Set that against the alternative it replaces, because the asymmetry is the point. Had `ordering/`
been split into `cpp/` and `python/` to give each IDE an uncontested root, changing our mind later
would mean moving eighteen source files, editing the Makefile's `%_cpp` rule and its test loop,
rewriting the README's file references, and carrying the git history across the renames. **Root
decisions are reversible and layout decisions are not**, so a root should be chosen for what is
convenient now and revisited freely, while a layout should be argued.

What a move does lose is whatever accumulated inside the old `.idea/`: the interpreter setting, run
configurations, the inspection profile, the code style. Here that is close to empty, the interpreter
being the system `python3` and the scripts needing no run configuration, and where it is not, those
are individual XML files that can be copied across rather than rebuilt.

**Where this stands today.** PyCharm's root is `experiments/ordering/`, since that is the only
Python in the tree. If prototyping spreads to other experiments it moves up to `experiments/`, on
the reading above, and that move is the two steps above and nothing else.

---

## 2026-07-31, `python3` over `uv` for the Python twins, because nothing here needs resolving

A companion to the entry above and a narrower question: given that Python now runs in this tree,
which interpreter runs it, in the shell, in the `Makefile`, and in PyCharm. The answer is bare
`python3` in the two places that matter to the repo, and whatever we like in the shell.

**The trigger for uv is a dependency graph or a pinned version, and this has neither.** uv is a
resolver and an environment manager: it earns its place where packages must be fetched and pinned,
where a project is installed or published, or where a Python version has to be held still. The nine
ordering layers import `sys` and `math`, import nothing from each other, install nothing and publish
nothing. An environment is a place to put resolved dependencies, so here uv would create one and it
would be empty. Reach for uv when something needs resolving and skip it when nothing does, which is
a sharper line than prototype against production: a three-package prototype wants uv immediately,
and a shipped tool with no dependencies does not want it at all.

**And the absence is a property of these files rather than a stage they will grow out of.** The
twins exist so that the Python can be read as the specification and the C++ trusted as the
implementation, with `make test` diffing their traces. A dependency would undercut exactly that. The
contrast is the data-structures repo, whose `python/` folder carries `pyproject.toml`, a venv, an
editable install and `ipykernel`, on far less code, because there something genuinely has to be
resolved.

**Three places, and only two of them are the repo's business.**

- **The shell.** Whichever tool is preferred; uv is a fine default there and the repo neither knows
  nor cares.
- **The `Makefile`.** `PYTHON ?= python3`, so `make test` runs under a `PATH` name rather than a
  particular tool: the prerequisite is then only that some Python 3 exists, where `uv run` requires
  uv installed first and then resolves a Python of its own.
- **PyCharm.** Decline the offered uv environment and select the existing `python3`, which is the
  same binary the `Makefile` resolves, so an upgrade moves both together and a `.venv/` this tree
  does not ignore is never created.

**What is deliberately not claimed.** Not that uv is less stable, which it is not, since it pins
exactly where Homebrew's `python3` follows whatever `brew upgrade` last installed. Not that the two
could produce different traces here, which they cannot: measured on 2026-07-31, uv's managed
interpreter and Homebrew's differed first in version (3.14.3 against 3.14.6) and then, after
`uv python upgrade`, only in build (Jul 23 / Clang 22.1.3 against Jun 10 / Clang 21.0.0), and
nothing in stdlib-only scripts with no `hash()` reaching a trace and an arithmetic key modulo n can
see either difference. The choice is about how narrow the prerequisite is and about the editor and
the build agreeing by construction, not about safety. `make test PYTHON="uv run"` is a supported
override and the `Makefile` documents it.

**Changing later is cheap, and nothing is arranged to make it so.** Should a dependency arrive,
three steps:

1. `uv init` and `uv add` in the experiment folder.
2. Point PyCharm's interpreter at the resulting `.venv`.
3. Set `PYTHON ?= uv run`, a one-word edit to a variable that exists for it.

Two chores travel along. Adding `.venv/` and `__pycache__/` to the root `.gitignore`, which this
tree lacked until Python arrived, and tracking `pyproject.toml`, which is an authored language
project rather than generated editor state. Nothing else moves: no file imports another, so no
import path changes, and the README's one `python3 md3.py 3` is a usage example rather than a
dependency. The genuinely irreversible part is not the tooling but the decision itself, since the
twins would stop being stdlib-only and the property that makes them readable as a specification is
what would be spent.

---

## 2026-07-24, The 2x2 pivot is trivial in the triangular sweeps and the whole cost of the diagonal pass

**Splitting `D` into its own solve pass looked like a small tidiness and turned out to be what keeps
the sweeps simple.** A dynamic LDL factor has 1x1 and 2x2 pivots, and the solve applies three
operators to the right-hand side in turn: `L`, then `D^-1`, then `L^H`. The question is where the
2x2 block makes trouble, and the answer is not where the block is, it is where the block is a
*coupled* object.

**In the forward and backward sweeps it is not coupled, because `L` is unit.** A 2x2 pivot does not
change that `L` is unit triangular; it only parks one number, `D`'s off-diagonal, in the storage
slot below the first column's diagonal, the slot that would otherwise hold an `L` entry. So the only
thing the sweep must know is that a column opening a 2x2 has its true `L` entries start one row
lower, because the row between holds a piece of `D`. That is one ternary, `pivotType[lj] != 2 ?
j + 1 : j + 2`, and nothing else. The second column of the block is an ordinary unit-`L` column and
needs no special case at all: a column of `L` is a column of `L` whether or not it sits inside a
block. `forwardDynamic` and `backwardDynamic` differ from their static twins by exactly that
ternary.

**In the diagonal pass the block is genuinely a block.** `D^-1` against a 1x1 is a scalar divide;
against a 2x2 it is a coupled 2x2 solve, and it has to be stable, so `diagonalDynamic` carries a
full 2x2 with *its own* partial pivoting, the `abs(a11) >= abs(a21)` row swap and the small LU
beneath it. Every consequence of the block structure is quarantined here, which is why the
conceptually trivial pass, a diagonal scaling, is the one with the real pivot logic in it, and the
substantial passes, the triangular sweeps, are the ones a 2x2 barely perturbs. The irony is worth
stating plainly: the *shorter* operator carries the *harder* code.

**Two independent layers of pivoting meet in the same 2x2, at different phases and for different
reasons.** The factorization pivoted once, at the threshold test (Section 7 of
`sparse_factorization.md`), to *decide* the block should be 2x2 at all, because the diagonal was too
weak to eliminate on: that is about growth and fill during elimination. The solve pivots again,
inside `diagonalDynamic`, to *invert* that now-frozen block against the right-hand side without
dividing by its smaller entry: that is just stable inversion of a fixed 2x2. Neither layer knows
about the other, and they must not be conflated, since they are answering different questions about
the same four numbers.

**So the decision is: keep `D` as a separate pass, do not fold it into the sweeps.** Cholesky folds
its scalar divide into the sweep and loses nothing, which is the tempting precedent. A 2x2 `D`
cannot be folded the same way without the coupling leaking into the sweep loop and turning that one
ternary into a block solve inside the hot triangular path. Pulling `D` out is what lets `L` and
`L^H` stay unit-triangular sweeps that treat a 2x2 as a one-row offset. The pass split pays for
itself not in the diagonal pass, which it complicates, but in the two sweeps, which it keeps
trivial.

---

## 2026-07-22, Forest parallelism, when we build it: portable threading first, architecture-targeted second

ARCHITECTURE records that Oblio gets node parallelism for free through Accelerate but that tree
parallelism, factoring independent forest branches on separate cores, is Oblio's to build. This entry
is the how, worked out ahead of the need so the reasoning is not lost. It is a recommendation, not a
commitment: nothing here is built. Apple Silicon is the immediate target, but portability is a
first-class constraint here, because Oblio's defining property is that it runs on any CPU with some
BLAS and LAPACK, and forest parallelism must not cost that. So the recommendation leads with what is
portable and treats the architecture-specific options as a deliberate trade, not the default.

**The shape of the work dictates the tool.** Forest parallelism is not a parallel loop. The
elimination tree is a dependency graph, a parent's task waiting on its children's, and the subtrees
are wildly uneven, one near the root enormous and a leaf subtree trivial. Static assignment, thread i
takes subtree i, load-balances badly. What the work needs is dynamic work-stealing over an irregular
task DAG, and that single requirement sorts every option below.

**The real fork is portable versus architecture-targeted, and portability wins by default.** The
options divide cleanly. Portable: OpenMP, and the C++ task-graph libraries (Taskflow, oneTBB), which
run on Linux, macOS, and Windows, x86 and arm64 alike. Architecture-targeted: Apple's GCD, and the
equivalent vendor stacks elsewhere, a different one per platform. The architecture-targeted route buys
tuning, P-core placement on Apple and the analogous thing on an AMD or Intel box, at the cost of a
separate implementation per platform and the loss of the any-CPU property. For Oblio that trade is
usually the wrong way round: the whole point of the solver is that it runs anywhere a BLAS does, so
the parallel layer should too, and the architecture-specifics belong as tuning on top of a portable
base, exactly as the choice of BLAS (Accelerate, MKL, OpenBLAS) is tuning on top of a portable
BLAS/LAPACK interface.

**For shared memory, the model is threads, not processes, which means OpenMP, not MPI.** Both are
portable, so both keep the any-CPU property, but they are not interchangeable here. MPI is a model of
separate processes that communicate by messages; on one multicore machine the ranks do run across the
cores, and a good MPI routes intra-node messages through a shared-memory buffer rather than a network,
but each rank still has its own address space, so data is replicated across ranks and a message is a
copy through a buffer. Threads share one address space, so a subtree handed to another core is a
pointer, not a copy, and nothing is replicated. On genuinely shared memory the thread model is simply
lighter, which is exactly why MUMPS, MPI-first by origin, grew an OpenMP layer for multicore nodes.
Oblio starts where MUMPS had to retrofit: it is already single-process, so a threading runtime is the
natural and lighter fit, and MPI would be the wrong tool for single-node forest parallelism even
though it is portable.

**pthreads: no.** It gives raw threads, mutexes, and condition variables and nothing else, so forest
parallelism on it means hand-building a work-stealing task pool with DAG dependencies, which is the
exact wheel every alternative already provides and a hard one to get right (deque design, termination
detection, load balancing). Reserve pthreads for a hard no-dependency constraint, and even then macOS
offers a better low-level primitive in GCD.

**Recommended, ordered for Oblio: OpenMP first.** It is a compiler feature rather than a library
dependency, so on Linux and Windows there is nothing new to link, and its `task` construct with
`depend` clauses maps onto the tree DAG directly, with the runtime doing the work-stealing. Its
decisive property for Oblio is graceful degradation: with OpenMP switched off the pragmas are ignored,
so a postorder tree walk annotated with tasks compiles and runs as valid serial code everywhere, which
preserves the any-CPU property exactly. It is also the layer MUMPS uses, a proven path. The
portability is the same kind Oblio already accepts for BLAS and LAPACK: a standard interface with a
per-platform backing chosen at link time, Accelerate or MKL or OpenBLAS behind the one and libomp or
libgomp or the MSVC runtime behind the other, with the code compiling against the interface rather
than the backing. The one asymmetry from the BLAS case is how the backing is obtained on Apple:
Accelerate is a macOS system framework that needs no install (`-framework Accelerate` and nothing
else), whereas Apple's Clang ships no OpenMP runtime, so macOS needs libomp from Homebrew or a
Homebrew toolchain. Every other platform has OpenMP in the compiler, and there the asymmetry can
invert, OpenMP built in while the BLAS is the thing to install.

**Then the C++ task-graph libraries, Taskflow or oneTBB, for better ergonomics at the price of a
dependency.** The elimination forest maps onto a task graph one to one, each supernode a task and each
child-to-parent edge a dependency, and the work-stealing scheduler absorbs the uneven subtree sizes,
which is the whole difficulty. Both are portable and build native on arm64 macOS; Taskflow is
header-only, oneTBB more mature, with `task_group` and `flow_graph`. They express the DAG more cleanly
than OpenMP pragmas, but they are a hard dependency to link rather than a compiler feature, and they
do not degrade to nothing when disabled. Reach for them if the pragma ergonomics grate and the
dependency is acceptable.

**Architecture-targeted, only if portability is being deliberately traded for tuning: GCD on Apple.**
Grand Central Dispatch is zero external dependency on macOS, callable from C++, and P-core / E-core
aware through QoS classes, so high-priority work lands on the performance cores, which the portable
options cannot request, since macOS gives no hard thread-to-core affinity, only QoS hints. The cost is
that it is Apple-only and expresses DAG dependencies more manually. It is the right choice only when
the platform is fixed and the P-core placement is worth a per-platform implementation, and the
equivalent statement holds for a vendor-specific stack on any other architecture.

**The caveat that bites regardless of choice, and on every platform: nested parallelism against the
BLAS.** Running N subtrees on N cores while each subtree's fronts call a self-threading BLAS
oversubscribes the machine. The fix is the tree-versus-node decision made concrete: inside a parallel
forest region, cap the BLAS to one thread per front (Accelerate, MKL, and OpenBLAS each expose a
thread cap), and let it thread out only when factoring one large front alone near the root. On Apple
Silicon the AMX-contention inversion makes this clean rather than a compromise: forest parallelism
pays at the leaves, small fronts that run on NEON and barely thread in the BLAS anyway, while the big
root fronts are left to Accelerate and the shared AMX single-stream. Forest across cores at the
bottom, the BLAS at the top, and the two do not fight. The oversubscription caveat is
platform-independent; only the AMX detail is Apple's.

**The GPU is the wrong axis for this, on any platform.** Forest parallelism is many small irregular
tasks, the opposite of what a GPU wants; small fronts on a GPU would drown in launch and transfer
overhead. The GPU is a separate axis, offloading the large dense root fronts, which is node
parallelism; on Apple that overlaps with what AMX already does, and everywhere it carries the
host-to-device transfer cost. A much larger and later lift, orthogonal to forest parallelism, and idle
by design for it.

**Bottom line.** OpenMP first, because it is portable, is a compiler feature rather than a dependency,
and degrades to valid serial code, so it keeps Oblio's any-CPU property intact; the C++ task-graph
runtimes (Taskflow, oneTBB) as the portable, nicer-ergonomics option that costs a dependency; GCD or
another vendor stack only when portability is knowingly traded for platform tuning. Not MPI, which is
portable but process-based and memory-heavy on a single shared-memory node; not pthreads, which means
hand-rolling the work-stealing DAG; not the GPU, which is the wrong shape for the task. And any of
these requires multifrontal, whose subtrees are already independent tasks with a per-subtree stack;
the looking traversals, with their cross-tree update flow, do not decompose this way.

---

## 2026-07-22, The update stack is random-access by supernode, not a true LIFO

The multifrontal update stack was, in 0.9 and in the classical design, an actual stack: each
supernode pushes its contribution block as it finishes, in first-to-next child order, and a parent
pops its children off the top, which is last-to-previous order since the last child pushed sits on
top. One arena, one stack pointer, the peak bounded into contiguous memory, and that contiguity is
the natural on-ramp to an out-of-core factorization.

**Our stack is not that.** It is a `std::vector<UpdateMatrix>` indexed by supernode: `allocate` fills
an arbitrary slot when its supernode is reached, `discard` frees an arbitrary slot once its parent
has assembled it, and the assembly walks children first-to-next. We keep neither the push/pop discipline
nor the reverse pop order. Read strictly, this departs from the original design, and the departure is
worth naming rather than leaving implicit.

**Correctness is untouched, on two independent counts.** First, the postorder numbering means a true
LIFO would in fact work here; the structure is present, we simply do not lean on it. Second, the
assembly is order-free regardless: `assembleUpdateMatrix` accumulates (`+=`, commutative) and
`assembleDelay`
assigns into disjoint positions, so assembling a parent's children in any order lands the same block,
give or take the last bit of the update sum. That last-bit freedom is not hypothetical; it is why the
tier-2 multifrontal pivot counts are bounded rather than pinned, and why multifrontal residuals
differ from left-looking in the low bits. The same summation-order freedom is recorded from the
testing side in TESTING_SPECIFICATION.

**What we give up is the contiguous arena, and only that.** The set of blocks live at any instant is
identical to what a true stack would hold, since a block lives from its supernode's turn until its
parent assembles it, so the peak in bytes is the same. It is merely scattered across per-supernode heap
allocations instead of managed by one stack pointer over one buffer. This is the same decision
recorded from the other side in the UpdateMatrix row of PORTING_LEDGER: we dropped 0.9's abstract
`UpdateStack` and its out-of-core concrete `UpdateStackDynamic` as unneeded in core. The departure
here is precisely that in-core choice, restated in stack-discipline terms.

**The one order that is not free, and must stay first-to-next, is the delay prepend.** The two orders
are decoupled, and that is the trap in any future "reverse everything to make it a real stack"
refactor. The assembly is free because it is commutative. The delay prepend is structural: it writes kk's
index set, which fixes the front's column order, and the pivot sequence reads that order. `gblToLcl`
keeps the values self-consistent for whatever prepend order is used, so reversing it is not a
correctness fault, but it would give kk a different though still valid column ordering, diverge from
0.9's canonical index set, and move the pinned tier-1 pivot counts. So a LIFO refactor would be right
for the assembly and wrong for the delays; the two cannot be swept together.

**When to revisit.** If a contiguous arena or an out-of-core path is ever wanted, reintroduce the
stack discipline for the assembly, with 0.9's `UpdateStack` / `UpdateStackDynamic` split as the reference
for the arena and the out-of-core spill. Leave the delay prepend first-to-next.

---

## 2026-07-19, The per-supernode lookups drop `Ptr`: `nodeIdx(jj)` and `val(jj)`

The naming entry below (07-17) settled the arrays and, in passing, gave `Ptr` a job: an offset is
named for its entity and takes a payload qualifier only when a sibling offset exists, which is why
the numeric factor carries `snodeNodeIdxPtr` and `snodeValPtr` while A collapses to `colPtr`. **The
accessors never followed that rule.** They were `nodeIdxPtr(jj)` and `valPtr(jj)`, and a reader
applying the rule would expect an offset array. They return neither an offset nor an array.

**What they are is lookups**, the same thing the matrix's per-column accessors are: an O(1) address
into storage the object already holds, computed at the moment of use, nothing allocated and nothing
owned. `experiments/storage-options` is where that shape was worked out, and its accessors are
named `rowIdx(j)` / `val(j)` / `colSize(j)`, with no suffix, over exactly the same reasoning. The
factor was the odd one out.

So the lookups are now `nodeIdx(jj)` and `val(jj)`, and the two spellings line up:

| | the lookup | the offset array behind it |
|---|---|---|
| **A** | `rowIdx(j)`, `val(j)` | `colPtr` |
| **L, numeric, static** | `nodeIdx(jj)`, `val(jj)` | `snodeNodeIdxPtr`, `snodeValPtr` |
| **L, numeric, dynamic** | `nodeIdx(jj)`, `val(jj)` | none: each supernode owns its vector |

**The rule, in one line: an offset array keeps `Ptr`; a lookup does not.**

Two things fell out of the rename rather than being designed into it. The static factor already
exposed whole-vector `nodeIdx()` and `val()`, so the renamed forms became overloads of them, which
is precisely `SparseMatrixStatic`'s `rowIdx()` / `rowIdx(j)` pairing: the whole array for a caller
who wants it, one entity's slice for a caller who wants that. And the read-public / write-private
split is untouched, because name lookup gathers every overload before access is checked, so the
`std::as_const` idiom behaves exactly as it did.

The shared read surface is now seven functions, and it is what makes one solve serve both storages:
`size`, `snodeSize`, `factorization`, `frontSize(jj)`, `updateSize(jj)`, `nodeIdx(jj)`, `val(jj)`.
Five are the same field read on both classes and only the last two differ, an offset into a flat
buffer against a per-supernode vector that already holds the pointer. The correspondence with the
matvec is close but not exact: `colSize(j)` has no single counterpart, because a supernode's value
block is a rectangle rather than a run and takes two numbers (`frontSize` for columns,
`frontSize + updateSize` for rows, which is also the leading dimension), and `factorization()` has
no counterpart at all, because a matrix carries no method identity while a factor must say whether
to conjugate.

**Earlier entries in this file keep the old spelling**, since they record what was decided when they
were written: the bulk-versus-direct entry (07-14) argues about `valPtr` as it then stood, and the
rejected extractor it describes materializes arrays it calls `rowIdxPtr[j]` / `valPtr[j]` / `len[j]`.
Those array names are a different thing from the accessor and are correct as written.

---

## 2026-07-17, Naming A and L: `val` stays `val`, and `nodeIdx` unifies row and column

A and L are both column oriented, and the names should show them to be the same idea, with L adding
one thing on top: a compressed structural view, for efficiency. Seeing exactly how far that
compression reaches, and where it stops, fixes every name.

**A is the base case: a column points at two parallel arrays.** `colPtr` brackets, per column, a run
of `rowIdx` (which rows are nonzero) and a run of `val` (their values). One offset serves both,
because the arrays are parallel: entry k of `rowIdx` and entry k of `val` are the same nonzero.

**L keeps the idea and compresses it, and the compression goes two separate ways.** A supernode holds
a set of columns together, and that has different consequences for values and for indices.

- **Values are held together, but not compressed.** Inside a supernode there are still row values,
  but now over more than one column, so `rowVal` would misread it: the run is a dense block spanning
  several columns, not one column's rows. So the value array stays **`val`**. And because we want
  `val` on L, we keep `val` on A as well, rather than split a matching pair of arrays into `rowVal`
  on A and `val` on L for no gain. In isolation A's values could read `rowVal`, and arguably that is
  a hair more descriptive, but consistency across A and L wins, and L's block reading settles it.

- **Indices are compressed.** Values cannot be compressed; indices can, and in L they are. The
  compressed set is then a mix of row and column indices. `rowColIdx` would be perfectly valid, and
  we shorten it: represent both a row and a column by a **node**, giving **`nodeIdx`**. Beyond being
  short, the word does two jobs. It unifies row and column, since in a symmetric structure a node is
  neither and both. And it hints at L's special structure, with nodes coming from the elimination
  forest view and supernodes from node compression.

**The offsets follow from what they bracket, and a qualifier is added only to tell siblings apart.**
A collapses to a single offset (`colPtr`) precisely because its two arrays are parallel: one offset
serves both `rowIdx` and `val`, so it is named for its entity alone, `colPtr`, not `colRowIdxPtr`.
The symbolic factor is in the same position. It holds the index array and nothing else, no values,
so it too needs one offset, and it too is entity-named: `snodePtr`. Moving from symbolic to numeric
is what forces a split. The numeric factor's indices are compressed and its values are not the same
shape (one `nodeIdx` entry is a whole row of a dense block, not a single value), so a single offset
cannot serve both and it carries two. Only now, with two offsets to keep apart, does each earn a
payload qualifier: `snodeNodeIdxPtr` (a supernode's slice of `nodeIdx`) and `snodeValPtr` (a
supernode's slice of `val`). So the rule is one line: an offset is named for its entity, and takes a
payload qualifier only when a sibling offset exists to disambiguate. The inverse map, from a node to
the supernode that owns it, is `nodeToSnode`.

So it is one spine, with the numeric factor the only place a second offset is needed:

| | index offset | indices | value offset | values | inverse map |
|---|---|---|---|---|---|
| **A** | `colPtr` | `rowIdx` | `colPtr` | `val` | |
| **L, symbolic** | `snodePtr` | `nodeIdx` | | | `nodeToSnode` |
| **L, numeric** | `snodeNodeIdxPtr` | `nodeIdx` | `snodeValPtr` | `val` | `nodeToSnode` |

This is settled and now carried by both factors, the symbolic (`SymFactor`) and the numeric, and
the elimination forest shares the same names for what it hands across, `snodeSize` and
`nodeToSnode`, so copying the forest's attributes into the factor is a name-for-name transfer
(`sf.mNodeToSnode = ef.nodeToSnode()`). The engine's own scratch offsets are entity-named for the
same reason, single-purpose so unqualified, and the transpose in `sortIndices` reads as a clean
mirror: `snodePtr -> nodeIdx` (a supernode and the nodes it holds) against `nodePtr -> snodeIdx` (a
node and the supernodes that hold it).

---

## 2026-07-15, An in-header throwing constructor slowed a hot loop in the same translation unit

A concrete, measured finding from the storage-options experiment, worth recording because it is
counterintuitive and it informs where numfact puts its code.

**Symptom.** After the matrix constructors gained the nnz/dimension guard (a `throw`), the static
direct `multiply()` slowed from ~1.03x to ~1.20x of the hand-written baseline on the M4 (Accelerate,
`g++ -O3`). Consistent across many runs, so not noise. The dynamic and baseline rows were unaffected.

**Diagnosis.** The multiply loop's own source had not changed, and the constructor is never called
anywhere near it. What changed is that the constructor, defined *inline in the header*, became a
potentially-throwing body visible in the translation unit that compiles
`multiply<SparseMatrixStatic>` (`MultiplyEngine.cpp`, which includes the matrix header for the
explicit instantiations). The exception path it introduced perturbed the optimizer's treatment of
the hot loop in that unit, even though the loop neither throws nor constructs anything.

**Confirmation, both directions.** Commenting out the `throw` restored 1.03x. Moving the whole
constructor into a `.cpp` (so the throwing body is no longer in the header, hence no longer in
`MultiplyEngine.cpp`'s unit) also restored 1.03x, with the guard intact. Object-file check: with the
constructor in its own `.cpp`, `MultiplyEngine.o` contains no `length_error` or `__cxa_throw`
machinery at all; the exception path is gone from that unit.

**Fix.** Both matrix constructors moved from their headers into `SparseMatrixStatic.cpp` /
`SparseMatrixDynamic.cpp`; the headers declare only. Accessors and the non-throwing mutators
(`setValues`, `setColumn`, which return false rather than throw) stay inline, since only the
constructor carries an exception path. This mirrors the main-code `SparseMatrix`, whose constructor
is already in `SparseMatrix.cpp`, so the experiment now matches the real tree here.

**The lesson, and why it is a third reason.** The usual "small bodies in headers, large bodies in the
`.cpp`" rule (CLAUDE.md, and the explicit-instantiation entry) is argued from inlining and build
time. This adds a distinct, optimization-quality reason: an exception path anywhere in a translation
unit can degrade the codegen of unrelated hot code in that same unit, at a cost that is not small
(17% here). So a throwing body is "heavy" for header purposes even when it is textually short. The
rule for numfact follows directly: keep throwing and heavy bodies out of the translation units that
compile the numeric kernels, since the factorization and solve loops are exactly the hot code this
would silently tax.

---

## 2026-07-15, The nnz cap is an ordering constraint, not a representability one; A is capped, L is not

A follow-on to the index-types entry (2026-07-09), sharpening one thing it left blurred: it speaks
of "the ~2.1-billion index cap" as if a single ceiling covered both indices and nnz. There are two
ceilings, with different origins, equal today only by coincidence.

**Two separate ceilings.** The *dimension/index* ceiling holds every id to `2^31 - 1`, because ids
are `std::int32_t`. But note the dimension `n` is capped one below what index representability alone
would allow, and for a distinct reason (see the next paragraph). The *nnz(A)*
ceiling is not that. nnz is a count, stored as `std::size_t` everywhere (`colPtr`, block offsets,
every position), so nothing internal caps it at `2^31`. What caps it is the ordering handoff: the
vendored AMD/MMD are the `int`-based build, and A's pattern reaches them as `int` arrays, `Ap` (with
`Ap[n] = nnz`) and `Ai`. So `nnz(A)` must fit `int`, `<= 2^31 - 1`, only to be handed to the orderer.
The two ceilings coincide at `INT32_MAX` solely because both are "largest value that fits a signed
32-bit int," one as an id, one as a count in `Ap[n]`. That coincidence is why one constant,
`MAX_IDX = INT32_MAX = 2^31 - 1`, does both jobs. It is a coincidence, not a shared fact.

**Why the dimension cap is `2^31 - 1`, not `2^31`.** Pure index representability would allow
`n = 2^31`: an `n x n` matrix then has largest index `n - 1 = 2^31 - 1`, which fits `int32`. What
removes that last value is the *loop counter*, not the index. Entity loops run
`for (std::int32_t j = 0; j < static_cast<std::int32_t>(size); ++j)`, and at `size = 2^31` the cast
is `INT32_MIN` (negative), so the loop runs zero times, silently; an `int32` counter also cannot reach
`2^31` without overflowing on the final increment. So `n <= 2^31 - 1` is the largest dimension an
`int32` counter can walk, and the dimension cap binds through the counter, not the index. The cost is
one unusable value: the index `INT32_MAX` is representable but never reached (the largest valid index
is `n - 1 = 2^31 - 2`), sacrificed to keep the `int32` entity-loop convention sound. This is why the
guarded constructors reject `size > MAX_IDX` (not `>=`), and it is the same ceiling nnz(A) meets from
the other direction.

**The cap is a forward-iteration limit; a reverse loop could reach `2^31`.** The overflow above is
specific to counting up. A reverse entity loop,
`for (std::int32_t jj = static_cast<std::int32_t>(size) - 1; jj >= 0; --jj)`, has no such problem: `jj`
runs from `size - 1` down to the exit value `-1`, and both ends fit `int32` even at `size = 2^31` (the
top is `2^31 - 1 = INT32_MAX`, the bottom is `-1`). The asymmetry is that each direction, on exit,
spills one step past its range of valid indices, and the two spills land on opposite sides of the
`int32` range: the reverse loop spills to `-1`, which `int32` represents fine, while the forward loop
spills to `size` itself, which at `size = 2^31` is `2^31`, one past `INT32_MAX` and not representable.
So the exit value, `-1` counting down versus `n` counting up, is what decides whether `n = 2^31` is
reachable. A codebase iterating entities purely in reverse could cap at `2^31`; ours iterates them
predominantly forward, so `2^31 - 1` is the honest cap and we keep it.

**The `size_t` counter trick, declined in both directions.** The spill can be absorbed by holding the
counter in `size_t`, which represents the spill value that `int32` cannot, and converting to an `int32`
index inside the body where only valid values occur. It takes a different form each way. Reverse:
`for (std::size_t t = size; t > 0; --t) { const std::int32_t j = static_cast<std::int32_t>(t - 1); ... }`,
with `t` walking `size` down to `1` and `j = t - 1` the index. Forward:
`for (std::size_t t = 0; t < size; ++t) { const std::int32_t j = static_cast<std::int32_t>(t); ... }`,
with `t` walking `0` up to `size` and `j = t` the index. In each, `t` absorbs the spill value (both
`size` and `0` fit `size_t`) and `j` converts only the valid indices `0 .. 2^31 - 1`, never the spill,
since the spill sits at the exit test after the last body. So either form would lift the cap to `2^31`.

We use neither, for opposite reasons. The reverse form is unnecessary: a direct `int32` reverse loop
already reaches `2^31` (its spill is `-1`, which `int32` holds), so `t` buys nothing and only adds a
variable. That is why the reverse countdowns in ElmForestEngine and SolveEngine were rewritten from this
`t` form to the direct `int32` reverse loop. The forward form is unwanted: it is the only way a forward
loop reaches `2^31` (a direct `int32` forward loop cannot), but it would add a second variable and a
cast to every entity loop, permanently and everywhere, to recover a single index at the very top of the
range that no real problem reaches (the memory wall sits far below it). Either way the answer is the
direct `int32` loop, and one unusable index is a cheap price for keeping every loop in its clean
one-variable form, which is the shape the whole convention is built on.

**int32 indices and size_t offsets are a self-consistent CSC pairing.** A fully dense matrix at int32
dimension has `nnz = n^2 <= (2^31)^2 = 2^62`, and `size_t` holds `2^64`, so `colPtr[n] = nnz` is
always representable. The condition is `offset_bits >= 2 * index_bits`, and `2 * 32 = 64` fits
exactly (int64 indices would need 128-bit offsets and would *not* fit, which is a further reason int32
is the comfortable resting point). Consequence: no storage or representability guardrail is needed on
A at all; a dense A is representable, and memory is the only real limit. The dense int32 ceiling is
~48 EiB (real), a formality no machine reaches this century. The only guardrail that exists is the
ordering one.

**Where the guard lives, and why only on the matrix.** The `SparseMatrix` constructor throws
`std::length_error` when `mSize > MAX_IDX || nnz > MAX_IDX`. It belongs to the matrix because A is
what gets ordered. It is a limit borrowed from a third-party interface, held at the one seam between
the `int`-based orderer and our `size_t` world, not a limit of our own storage.

**L is not capped, and the factor must not import the cap.** L, and the dynamic factor during
elimination, is never handed to an orderer, so the nnz ceiling does not apply. Its offsets are
`size_t`, its row ids are int32 (same id space, so L inherits the *dimension* bound but not the *nnz*
bound), and it grows to memory. The factor's constructor carries no nnz cap and throws nothing. The
one discipline: nothing on the numeric or solve path may narrow an L count or offset to `int` out of
symmetry with A. The correct design for the factor is the *absence* of the guard. L cannot approach
the dense ceiling regardless: it is triangular (`nnz(L) <= n^2/2`), and a dense block within it is
bounded by `dimension^2` with per-side dimensions bounded by `n <= 2^31 - 1`, exactly what BLAS's
`int` m/n/k accept, so no block dimension can overflow the `int` BLAS wants. That is foreclosed by the
id type, not by any check we write.

**Liftability.** The nnz(A) cap is cleanly removable: SuiteSparse ships a 64-bit
(`SuiteSparse_long`) AMD build, so lifting it is a link-and-widen of the ordering interface, not an
algorithm change, and the constructor's nnz guard then becomes a clean deletion (indices stay int32;
only the orderer's integer width changes). The reason AMD/MMD store nnz in an `int` at all is not
carelessness but a uniform-width choice: the orderer builds and mutates a quotient graph whose ids and
offsets share one integer type, which is cache-friendly for a pointer-chasing kernel, and it ships two
instantiations (`int` and `_long`) rather than mixing widths. So the "better design" is the `_long`
build, a link choice, not a rewrite. Keeping our cap localized and honestly labeled (A-side,
ordering-driven) is what makes its eventual removal a deletion rather than an archaeology problem.

Recorded now because the storage-options experiment's matrix constructors just gained this guard
(mirroring `SparseMatrix`), and it is the model for numfact: the matrix keeps the cap, the factor
analog drops it.

**The BLAS integer width is a link-time contract, not a compile-time one.** All of the above assumes
BLAS's `int` is 32-bit, which is what lets `static_cast<int>` of a `size_t` dimension be safe: the
per-side block dimensions are bounded by `n <= 2^31 - 1`, and a 32-bit signed `int` holds them. That
assumption is correct on every build we target (LP64: reference BLAS, Accelerate, OpenBLAS default,
MKL's `lp64` interface all use a 32-bit Fortran `INTEGER`, which the `int` prototypes in
`BlasLapack.h` express faithfully). It is worth naming because it can be violated without a compiler
error. An ILP64 BLAS (MKL's `ilp64` interface, OpenBLAS built `INTERFACE64=1`, Accelerate's newer
ILP64 variant) uses a 64-bit `INTEGER`; linking one against our 32-bit `int` prototypes is an ABI
mismatch that both sides compile clean and that surfaces only at run time, as wrong results on large
problems or a crash. The danger runs `size_t -> int`, never `int32_t -> int`: our row ids are
`int32_t` and equal `int` in width on these targets, so that cast loses nothing; the real narrowing
is the `size_t` dimension, and it is safe only while the BLAS `int` is at least 32 bits. If Oblio ever
needs the ILP64 range, the change is localized to `BlasLapack.h` (the raw prototypes and the inline
wrappers swap `int` for a 64-bit integer) because the rest of the code already speaks `size_t` and
narrows only at those wrapper calls; nothing else would move. This pairs with the `SuiteSparse_long`
liftability note above: both the orderer's integer and the BLAS integer are link choices, widened
independently, neither an algorithm change.

---

## 2026-07-14, Bulk versus direct access: direct wins, and bulk was never worth its one advantage

A consumer that reads a CSC-style object can do it two ways. *Direct*: hold the object and ask it
for one column's (or one supernode's) pointer and length at the moment of use, through the storage's
own lookup (`rowIdx` / `val` / `colSize` on the matrix, `valPtr` on the factor). *Bulk*:
extract every pointer up front into three plain arrays, then run a kernel that takes nothing but
those arrays, so the object is out of the picture and its dimensions must be extracted alongside.

Two questions decide between them, and they have different answers.

**Is bulk safe? Only for a consumer that does not mutate the object during its own sweep.** A
read-only sweep (matvec, the solve, both read a frozen structure) may be bulk: nothing grows, so no
extracted pointer goes stale. A mutating sweep must be direct: the numeric factorization grows a
front under itself (a delayed pivot reallocates an ancestor's buffer), so a pointer extracted up
front dangles, silently. This is the hard constraint, and it puts `NumFactorEngine` on direct by
necessity rather than preference.

**Is bulk worth it, where it is safe? No.** Bulk's single advantage is one storage-blind compiled
kernel, verifiable in the symbol table, against direct's one instantiation per storage. For us that
trades a real cost for nothing we need:

- The instantiation it saves is cheap. We have a small fixed set of storages, and the arithmetic
  kernels stay shared regardless (they take raw pointers); only the traversal, which is bookkeeping,
  monomorphizes. One extra copy of the bookkeeping per storage is noise.
- It is slower. Bulk streams three extra arrays that direct never touches, plus an `O(n)` extraction
  pass, and on a memory-bound sweep that traffic shows at full price. Measured, the matvec: bulk runs
  1.05x to 1.07x of hand-written flat, while direct matches hand-written flat.
- It is more API and more hazard. The extractor is method surface a reader must understand, and it
  carries an invalidation warning (a stale extracted pointer) that exists only because something was
  extracted. Direct deletes the method, the warning, and the paragraph explaining the warning, all
  at once.

So the rule is direct, and the storage's lookup is the whole interface: a fact about the layout,
answered by the class that owns it, called by a consumer templated on the storage type. No consumer
carries an extractor of its own, which was the repetition the lookup-versus-view corollary set out to
kill and then, by leaving `columnPointers` on the engine, half-kept. That is corrected in the same
entry: `columnPointers` was a bulk copy of the lookups, not a genuine view, and once the lookups
moved onto the matrix it was removed and the matvec reshaped to direct access.

| engine | access | why |
|---|---|---|
| NumFactorEngine | direct | must be: the factor grows under it, so a bulk pointer would dangle |
| SolveEngine | direct | frozen factor, so bulk would be safe, but direct is faster and smaller |
| MultiplyEngine | direct | frozen matrix, same reason |

One thing stays genuinely open, and only one: the deferred multi-RHS solve. There the kernel is
BLAS-3, `O(n^3)` work on `O(n^2)` data, so an extraction's cost amortizes over heavy block work
rather than showing at full price as it does on the memory-bound single-RHS sweeps. If bulk ever
earns its keep it is there, and even then the arrays would be built once from the same storage
lookups, in one shared helper, never a method per engine. The no-per-consumer-extractor rule holds
regardless of how that measures.

This aligns the matrix and the factor on one pattern, which is the point: whether an object is static
or dynamic because the user changes it (the matrix) or because the algorithm changes it (the factor),
the consumer sees one lookup interface and templates over the storage. See "an object offers what its
storage makes cheap" for the lookup-versus-view rule this completes, and the flat-versus-VV entry for
the storage taxonomy.

## 2026-07-14, An object offers what its storage makes cheap, and nothing else

A third rule in the same family as the two below, and it came out of naming the storage-options
experiment's classes for their *purpose* (`SparseMatrixStatic`, `SparseMatrixDynamic`) rather than
their layout (`Csc`, `Vv`).

**The two classes do not have the same API, and that is the design.**

| | static (flat) | dynamic (vector of vectors) |
|---|---|---|
| `setValues`, same structure, new numbers | **yes**, cheap: nothing moves | **yes**, cheap |
| `setColumn`, one column's structure | **absent by design** | **yes**, cheap: the column owns its buffer |
| restructure | build a new one | `setColumn` |

**`setValues` is the mutation a solver actually performs most often**, a Newton iteration, a time
step: same pattern, new numbers, refactorize. The flat layout is perfectly happy with it, so both
classes have it. This is also why "we need a mutable matrix" is a weaker argument for VV than it
first appears: the common case needs no VV at all.

**Amendment (2026-07-14): `setValues` landed at column granularity.** The experiment settled on
`setValues(std::int32_t j, const std::vector<double>& val)`, one column, with an *identical
signature on both classes*: static overwrites the contiguous run `mVal[colPtr[j]..]`, dynamic
overwrites `mVal[j]` in place, and neither invalidates a pointer. Setting every value is a loop over
columns, the write-twin of reading every column through the accessors. This strengthens the point
rather than changing it: the cheap operation is now not merely present on both classes but the *same
call* on both, and the asymmetry stays exactly where it belongs, on `setColumn` (structural), which
remains dynamic-only.

**`setColumn` is absent from the flat one, and its absence is the point.** Changing a column's
*structure* there means shifting every later column: `O(nnz)`, not `O(column)`. An API that looks
cheap and is secretly linear in the whole matrix is a trap, and the caller who writes it in a loop
will not find out until their program crawls. **Refusing to offer it is not a limitation; it is
telling the truth about the storage.**

And it puts the decision where it belongs. The caller knows whether they are changing one column or
rebuilding, and can pick the object that suits: *want a column replaced for you, use the dynamic
one; want to shift data around a flat buffer, do it yourself, this class will not pretend it is
cheap.*

**This is the API-side argument for having no common base class**, and it is stronger than the
performance one already recorded. A shared interface would force one of two lies:

- `setColumn` on the flat matrix, pretending an `O(nnz)` shift is a column operation, or
- `setColumn` on **neither**, crippling the dynamic one to match its sibling's weakness.

So the asymmetry between `NumFactorStatic` and `NumFactorDynamic` is not a wart to be tidied away
when dynamic LDL lands. It is what the two storages **are**, and the interfaces should say so.

**The rule: an object offers what its storage makes cheap, and nothing else.** Which is the same
instinct as the two rules below, seen from the API rather than from the conventions: *do not name a
thing you cannot honour.*

**A second corollary: a lookup belongs to the storage; a view belongs to the consumer.**

This one took a wrong turn first, which is worth recording because the wrong version *looked*
right.

The storage-options experiment puts `columnPointers` on the **engine**, not on the matrix, and that
is correct. The natural generalization is "access functions live in the algorithm class", and on
that reasoning `blockOf` was first written as a private helper on `NumFactorEngine` and again on
`SolveEngine`. Two copies of one expression, which is exactly the duplication the rule below
forbids.

The generalization was too broad. The two things are not the same kind of thing:

| | `columnPointers` | `valPtr` |
|---|---|---|
| what it does | **materializes** three new arrays | **computes an address** in the existing storage |
| cost | `O(n)`, called once | `O(1)`, called in the hot loop |
| ownership | somebody owns the arrays | nothing is owned |
| shaped by | the traversal that wants it | the storage that holds it |

`columnPointers` is a **view**: a new structure, built in the shape one algorithm wants. A different
algorithm would want a different one, so putting it on the matrix would grow the matrix a method
per consumer, each committing it to a format chosen for somebody else. **A view belongs to the
consumer.**

`valPtr` is a **lookup**: it answers "where does supernode kk's block live", and the answer is a
fact about the layout, not a shape chosen by anyone. The layout (one flat buffer with offsets, or
one vector per supernode) is the factor's own business, and no consumer should have to restate it.
**A lookup belongs to the storage.**

So `valPtr` now lives on `NumFactorStatic` and `NumFactorDynamic`, private, with the engines as
friends. One definition per storage rather than one per engine, and the engines cannot tell the two
factors apart. It still inlines to nothing (no `valPtr` symbol survives in either object file), so
the move is free.

**Amendment (2026-07-14, later the same day): the view was not a view, and the right answer was no
extractor at all.** The rule above is right and `valPtr` is right, but it mis-cast
`columnPointers`. Calling it a view and leaving it on the engine treated it as a genuine
algorithm-shaped structure, which it is not: the three arrays it builds are a bulk copy of the
per-column lookups the matrix already answers (`rowIdx` / `val` / `colSize`, the matrix-side
twin of `valPtr`, added to the storage after this entry was written). Once those lookups exist a
consumer reads a column directly, at the point of use, and the extractor is pure redundancy: an
`O(n)` up-front copy of what the storage already holds, carrying an invalidation hazard that exists
only because something was extracted. Its one property, a single storage-blind compiled kernel, has
no value to us: a small fixed set of storages, so direct costs one instantiation per storage, and
bulk measures a few percent slower on the memory-bound matvec besides. So `columnPointers` was
removed and the matvec reshaped to direct access, three lookups per matrix and one templated
`multiply` that calls them, no extractor. What stands: a lookup belongs to the storage; a genuine
view (a transpose, an algorithm-specific reordering) still belongs to the consumer, but
`columnPointers` was never one. See the bulk-versus-direct entry above.

**And a corollary that dynamic LDL will live or die by: structural mutation invalidates every
pointer previously extracted; value mutation does not.**

```
setValues   does NOT invalidate.  The buffer stays put; only its contents change.
setColumn   DOES invalidate.      The column's buffer is replaced; anything into it dangles.
```

The experiment demonstrates it rather than asserting it (`testInvalidation` extracts the pointers,
mutates, and observes which still point where they did). The rule holds in both storages, and it is
exactly the rule the dynamic factor needs: a delayed pivot grows an ancestor's front, which
reallocates its buffer, which dangles every pointer previously taken into it. **Silently.**

```cpp
eng.blockPointers(f, block);        // extracted once
for (kk) {
    ... factor kk ...
    f.mVal[pp].resize(bigger);      // a delayed column grows ancestor pp
    ...                             // block[pp] is now DANGLING
}
```

Nothing in C++ enforces this. The remedy is to **fetch a supernode's block pointer at the moment of
use rather than up front**, which is one indirection and which `storage-options` measured at
essentially nothing. The alternative, re-extracting after every growth, is more code and more to
forget.

This is the one thing the experiment does *not* rehearse, since its structures change only between
runs of the algorithm, never during one. Worth knowing before writing dynamic LDL rather than
after.

## 2026-07-14, Two rules for conventions: one predicate, and a sparse choice matrix

Two ideas that keep recurring, worth stating once rather than rediscovering.

### One predicate, in one place, applied everywhere

**When several places must agree about a convention, they must not each decide it.** Give them one
named thing to ask, and let it be the only place the decision exists.

**Example 1, and it is the reason for the rule: `hermitian()`.** The factorization and the solve
must both know whether the factor conjugates. 10.12 lets each site decide for itself:

```cpp
// in the factorization
SYRK('L','N', ...)              // decides: no conjugate
GEMM('N','T', ...)              // decides: no conjugate, again

// in the solve, a different file
y[lc] -= y[lr] * val[lij];      // decides: no conjugate, a third time
```

Nowhere is it written down that a Cholesky factor is Hermitian. Each site *inferred* the
convention from context, and each inherited LDL's (symmetric) habit. So Cholesky is wrong in
**all three**, independently, and fixing one would not fix the others. Nothing reveals the
disagreement, because there is nothing to disagree *with*.

Ours has one predicate:

```cpp
static bool hermitian(Factorization f) {
    return f == Factorization::Cholesky
        || f == Factorization::StaticLDLH
        || f == Factorization::DynamicLDLH;
}
```

The factorization asks it; the solve asks it. Add a factorization type and both follow. Get the
predicate wrong and *everything* breaks loudly, which is far better than one site being quietly
wrong.

**Example 2, and it shows the other half of the rule: `exactPatterns()`.**

Symbolic factorization may read one front column per supernode when the supernodes' columns share
a pattern exactly, and must read them all when they do not. Who decides?

Not the caller. A `bool useExactPatterns` parameter on `SymFactorEngine` would put the decision in
the hands of whoever wires the phases together, and that person **cannot know the answer**:
amalgamation at threshold zero merges only free merges, so it *may or may not* store a zero
depending on the matrix. Threshold alone does not determine it. Only the forest engine, which did
the merging, knows what it actually did.

So the forest **records** it (`ElmForest::exactPatterns()`, set false by `compressThreshold` only
when a merge actually pays fill), and symbolic factorization **reads** it. One place decides, at
the moment it can, and everyone downstream asks rather than infers.

**The two examples differ in an instructive way.**

| | `hermitian()` | `exactPatterns()` |
|---|---|---|
| what it is | a pure function of a setting | a fact about what happened at runtime |
| where it lives | computed on demand, anywhere | **recorded** by the code that knows |
| why | the answer is determined by the input | the answer is not determined by any input |

So the discipline is one rule with two shapes: **if the convention is a function of the inputs,
make it a predicate; if it depends on what happened, record it at the moment it happens. Either
way, decide once and let everyone else ask.** The failure mode both avoid is the same: a
convention re-derived at each use, drifting apart, with nothing to compare against.

`Blas<Val>::conjTrans` is the same rule at the smallest scale: one trait, not `'T'` and `'C'` typed
out at each call site.

### The choice matrix is sparse. Enumerate its nonzeros, not the product.

The obvious modeling of "what do we factor" is a product:

```
(Cholesky, LDL)  x  (transpose, conjugate-transpose)   ->  four combinations
```

and then one of the four (**complex symmetric Cholesky**) has to be *forbidden*, because positive
definiteness is meaningless for it. A runtime rejection, guarding a combination the API invited
the caller to ask for.

**We did not do that.** The enum lists only the combinations that exist:

```cpp
enum class Factorization {
    Cholesky,      // always A = CC^H. There is no CC^T variant to name.
    StaticLDLT, StaticLDLH,
    DynamicLDLT, DynamicLDLH
};
```

Cholesky carries no transpose suffix because **it has no choice to offer**: it is always the
conjugate transpose, and over the reals that *is* the plain transpose. LDL carries one because it
genuinely has two forms. The forbidden cell is not rejected; **it does not exist**.

The cost is one letter of naming convention. The gain is that a whole class of runtime error is
unrepresentable, which is the same principle as `BlasLapack`'s operation-named wrappers: **make the
wrong thing unwriteable, not merely refused.**

Worth naming the general shape, since it will recur: a design space is rarely a full product. Model
the combinations that mean something and leave the rest unnameable. In a sparse-matrix library the
metaphor is right there: **the choice matrix is sparse, so store its nonzeros.**

## 2026-07-13, The solve, and the first test that checks the pipeline rather than a phase

**Every test before this one checks a phase against an oracle.** The forest against a
recomputation, the symbolic factor against a dense pattern, the numeric factor against a dense
Cholesky or by reconstruction. Each says *this phase computed what it should*. **None says the
phases compose.**

```
|| A x - b ||  /  || b ||
```

says that, in one number, through ordering, elimination forest, symbolic factorization, numeric
factorization, triangular solve and sparse matvec. Six phases. It is where two phases disagreeing
about a convention would show, an ordering, a conjugate, an index base, and nothing else we have
would catch that. It comes out at 3e-16 for every factorization, both traversals, both scalar
types.

The right-hand side is manufactured from a known solution (`b := A x`), so the test needs no
reference solver. And it checks the **residual**, not the distance to the manufactured `x`: those
differ by the conditioning of `A`, and the residual is the honest thing to require of a direct
solver.

**`MultiplyEngine` exists for this**, and only for this. It is fifteen lines and it is what makes
the check possible.

**One right-hand side, and the scalar solve that follows from it.** 0.9 has two vector classes and
the split is principled: with a *single* right-hand side there is no level-3 BLAS to be had, so its
`SingleVector` solve is scalar and works directly on the vector through indirect indexing. With
*many*, a supernode's rows become a dense block and the solve becomes TRSM and GEMM, which is worth
the gather and scatter that packing demands. We take the first. The multi-column path is real and
worth adding; it is a performance path, not a correctness one.

**And the conjugate in the backward pass is where 10.12 repeats its bug.** Its backward solve is

```cpp
y[col] -= y[row] * val[...];        // no conjugate: this applies L^-T, not L^-H
```

which is right for its complex-symmetric LDL and **wrong for its Cholesky**, exactly as its `SYRK`
is. Ours conjugates when the factorization does, using the *same* `hermitian()` predicate the
factorization uses, so it is one rule stated once rather than two rules that must be kept in step.
That is the pattern to hold on to: **when two places must agree about a convention, give them one
predicate, not two copies of a decision.**

## 2026-07-13, Static LDL: three kernels BLAS does not have, and why the traversals did not change

**LAPACK has no unpivoted LDL.** `?sytrf` is Bunch-Kaufman, which pivots, and pivoting is exactly
what a *static* factorization refuses to do. So the kernel is ours. 0.9 wrote it; we port it.

**Three kernels, all from 0.9, all absent from BLAS:**

- **`ldl`** (0.9's `OBLIO_POTRF2`). Recursive: split the block in half, factor the leading part,
  solve for the off-diagonal, form its upper counterpart, take the Schur complement, recurse.
- **`formUpper`** (`OBLIO_COMPUTE_U`). `U := D L^T` (or `D L^H`), into the upper triangle.
- **`gemmLower`** (`OBLIO_GEMM`). `A -= L U`, filling only the lower triangle, because the product
  is symmetric. BLAS has nothing for this: `syrk` computes `A A^T`, and there is no "`A B` where
  the product is known symmetric".

**A nice factoring in 0.9 worth naming.** `OBLIO_POTRF1` and `OBLIO_POTRF2` are the *same
algorithm*, differing only at `n == 1`: the first fails on a non-positive pivot, the second
replaces a tiny one and counts it. So 0.9's `NO_LAPACK` "Cholesky" is really an LDL that refuses
to proceed unless `D > 0`, which is legitimate (`LDL^T` with positive `D` *is* Cholesky, with
`L_chol = L_ldl sqrt(D)`; same factorization, different storage). We need only the second, since
LAPACK gives us a real Cholesky.

**The storage, which the statistics predicted.** In an LDL block:

```
the diagonal        holds D          where L's implicit 1s would be
the lower triangle  holds L          unit lower triangular
the upper triangle  holds U = D L^T  which Cholesky leaves as explicit zeros
```

So an LDL block uses the whole rectangle. `U` is not redundant: the recursion needs `D L^T` in two
places (to solve for the next `L`, and to form the Schur complement), so it is computed once and
kept. We reasoned to exactly this layout from 0.9's `numberOfEntries` versus
`numberOfAllocatedEntries` counts, before seeing the code.

**Perturbation is not a refinement; it is the only recourse.** A static factorization cannot pivot,
so a pivot too small to divide by has no remedy but replacement. We then factor a matrix slightly
different from the one we were given, which is a real thing to have done, so the count is reported
on the factor (`NumFactorStatic::numPerturbations`) rather than hidden. It belongs on the factor,
not the engine: it is a property of what was computed, not of the thing that computed it.

**And the update is not a rank-k operation, which is the deep difference from Cholesky.**

```
Cholesky:  T -= L21 L21^H          one HERK. Exactly what a rank-k routine computes.
LDL:       T -= L21 D L21^H        no BLAS routine at all: the D in the middle rules it out.
```

So LDL forms `U := D L21^H` into a **scratch block**, then multiplies. That scratch is the whole
price of the `D`, and it is why Cholesky's update is one call and LDL's is three.

**The traversals did not change by a line, and that was the design working.** `factorSupernode` and
`updateSupernode` dispatch on the factorization type internally, so left-looking and right-looking
decide *when* to factor and *when* to update, never *how*. Adding LDL touched those two functions
and one `switch`. That separation is worth protecting when multifrontal arrives: a traversal is a
schedule, not an arithmetic.

**Complex Hermitian LDL (`StaticLDLH`) is an extension, not a port.** 0.9's complex LDL is
symmetric only, which is why its complex LDL correctly uses `SYRK` and `'T'`. We support both, and
the `T`/`H` distinction runs through all three kernels as a single `bool hermitian`, which for
`double` is a no-op in every branch. That is the honest expression of "real symmetric and real
Hermitian are the same case".

**Checked by reconstruction, which is a better oracle than a second implementation.** Cholesky is
compared against an independently written dense Cholesky. LDL is checked by multiplying the factor
back out: `L D L^H == P A P^T`. That needs no second implementation, validates `D`, `L`, the
storage layout and the supernodal assembly in one statement, and works unchanged across all three
symmetries, which a dense oracle would not.

## 2026-07-13, Symmetry is part of the factorization, not a setting beside it

*(Supersedes the "symmetry is determined, not chosen" half of the entry below, written earlier
today. That claim was true for Cholesky and false the moment LDL arrived. Recording the
correction rather than editing the original, because the reason it was wrong is worth keeping.)*

**What changed.** The earlier entry argued that `Val` and the factorization type fix the symmetry
between them, so symmetry never needs to be named:

```
real    + anything   ->  symmetric
complex + Cholesky   ->  Hermitian
complex + LDL        ->  symmetric
```

The third line is where it breaks. It is what **0.9** does, not what LDL *is*. A complex matrix
may be symmetric (`A = A^T`, so `D` comes out complex) or Hermitian (`A = A^H`, so `D` comes out
real), and both are legitimate, useful factorizations. `LDL^T` and `LDL^H` are different
computations, and nothing in `Val` or in "static LDL" says which one is wanted. **Symmetry becomes
a genuine choice.**

**The choice goes into the factorization type, not beside it.**

```cpp
enum class Factorization {
    Cholesky,      // A = CC^H  (CC^T for real). Positive definite, hence Hermitian.
    StaticLDLT,    // A = LDL^T          complex: D complex
    StaticLDLH,    // A = LDL^H          complex: D real
    DynamicLDLT,
    DynamicLDLH
};
```

**Why here and not as a `Symmetry` flag beside it**, which was the obvious alternative: the same
principle that shaped `BlasLapack`, **make the wrong thing unwriteable.** A separate flag would let
a caller ask for `Cholesky` + `Symmetric` + complex, which we would then reject at runtime. In the
enum above there is no such value. The combination is not *forbidden*; it does not *exist*. That is
strictly stronger, and it costs one letter in a name.

A flag would also be **inert for real**, where symmetric and Hermitian coincide. An API parameter
that is sometimes meaningless is a small lie, and worth avoiding for the price of an enumerator.

The cost, stated honestly: the enum grows, and `T`/`H` is a suffix convention a reader must learn.
Both seem cheap against a runtime rejection that exists only to guard a combination we could
simply not have named.

**And complex Hermitian LDL is an extension, not a port.** 0.9 does only the symmetric one, which
is why its complex LDL correctly uses `SYRK` and `'T'`. That its *Cholesky* also uses them is the
bug we have already recorded; the two are related, since Cholesky is the one place 0.9 needed the
Hermitian convention and did not reach for it.

## 2026-07-13, The factorization space, and a BLAS layer that names operations rather than routines

Entering the numeric phase, three questions had to be settled together: which combinations of
scalar type, symmetry and factorization we support, how the code selects among them, and what
shape the BLAS wrapper takes. They turn out to be one question.

**The space, and it collapses.** Three axes, `Val` in {real, complex}, symmetry in
{symmetric, Hermitian}, factorization in {`CC^T`, `LDL^T`}:

| `Val` | symmetry | `CC^T` | `LDL^T` |
|---|---|---|---|
| real | symmetric **is** Hermitian | yes, SPD | yes, indefinite |
| complex | Hermitian | yes, HPD (`zpotrf`) | yes; 0.9 does not do it, we will |
| complex | symmetric | **forbidden** | yes, the standard case |

**For real the symmetry axis does not exist.** Conjugation is the identity, so `A^H = A^T` and
the two conditions are the same condition. That row is one row wearing two hats.

**Complex symmetric Cholesky is forbidden, and not merely unimplemented.** Positive definiteness
means `x* A x > 0`, which requires that quantity to be *real*, which happens for all `x` exactly
when `A` is Hermitian. For a complex *symmetric* `A`, `x^T A x` is a complex number and the
inequality does not even typecheck. Concretely: `A = [[0, 1], [1, 0]]` is symmetric and
nonsingular, and `a11 = 0` kills the very first square root. Cholesky has no pivoting to recover
with, because not needing pivoting is the entire point of Cholesky, and that guarantee comes from
positive definiteness. LAPACK has no complex-symmetric Cholesky for precisely this reason. So the
API rejects it, as a hard error rather than a to-do.

The stable factorization for complex symmetric matrices is `LDL^T` with 2x2 pivots, which is
already on the plan. Complex symmetric `CC^T` would be an unstable duplicate of a thing we are
building anyway.

**Complex Hermitian `LDL^T` (`A = LDL^H`, `D` real) is a gap in 0.9, not a gap in the
mathematics.** It is a perfectly good factorization, indefinite Hermitian, and we intend to
support it. 0.9 simply never wrote it. Worth distinguishing sharply from the cell above: one is a
mathematical impossibility and must be rejected forever, the other is work not yet done. They
must not report the same error.

**So symmetry is determined, not chosen, and there is no third parameter.** Given `Val` and the
factorization type, everything else follows:

```
real    + anything   ->  symmetric  ->  'T',  syrk
complex + CC^T       ->  Hermitian  ->  'C',  herk
complex + LDL^T      ->  symmetric  ->  'T',  syrk
```

`Val` is the template parameter and `FactorType` is an engine setting (10.12 already has the
enum: `eCC`, `eStaticLDL`, `eDynamicLDL`). Between them the transpose character and the rank-k
routine are fixed. No symmetry flag, no extra argument.

**The BLAS layer names operations, not routines. This is the part 0.9 gets wrong.** A wrapper is
needed regardless: BLAS is a Fortran interface, everything by pointer, and without a layer every
call site carries `&uplo, &trans, &n, &k, &alpha, ...` plus a branch to choose `d` or `z`. 0.9
has such a layer, overloaded inline wrappers on `Real*` and `Complex*`, and that much is right.

What is wrong is *what it wraps*. 0.9 wraps BLAS **routine by routine**: `SYRK`, `GEMM`, `TRSM`,
`POTRF`, leaving the caller to choose `'T'` versus `'C'` and `SYRK` versus `HERK`. The convention
therefore leaks into the engine, and the engine gets it wrong: **0.9's complex Cholesky calls
`SYRK` and `TRSM('T')` and `GEMM('N','T')`, all of which are the complex *symmetric* pattern,
while `POTRF` maps to `zpotrf_`, which is Hermitian.** There is no `HERK` anywhere in 0.9's
`BlasLapack`. For Hermitian `A = CC^H` the update must be `L21 L21^H` and the off-diagonal solve
must be against `L11^H`; using `'T'` is correct only when `L11` is real. Almost certainly never
exercised on a genuinely complex Hermitian matrix.

Our layer therefore exposes an operation whose *meaning* is fixed and lets the type pick the
routine:

```cpp
// "A times A-conjugate-transpose", whatever that means for this Val.
void herk(char uplo, char trans, ...);   // double -> dsyrk_ ;  complex -> zherk_

template<class Val> struct Blas;
template<> struct Blas<double>               { static constexpr char conjTrans = 'T'; };
template<> struct Blas<std::complex<double>> { static constexpr char conjTrans = 'C'; };
```

so the Cholesky kernel is one piece of code, correct for both:

```cpp
potrf('L', f, block, ld);
trsm ('R', 'L', Blas<Val>::conjTrans, 'N', u, f, ...);
herk ('L', 'N', ...);
gemm ('N', Blas<Val>::conjTrans, ...);
```

No branch, no `if constexpr`, and **0.9's bug becomes unwriteable**: the engine never names
`SYRK` or `HERK`, so it cannot pick the wrong one. `syrk` and a literal `'T'` remain available in
the header for `LDL^T`, where plain transpose is what is wanted, and there the *algorithm* asks
for them explicitly rather than inheriting them by accident.

**Storage: this is what the storage-options experiment was for.** Yesterday's study
(`experiments/storage-options/`) established that one compiled algorithm serves both a flat CSC
storage and a vector-of-vectors, through nothing but a pointer array, and that the abstraction
costs about nothing (1.07x flat, 1.10x packed VV, one `multiply` symbol in the binary). It was
run against exactly this moment. `NumFactor` uses both:

- **static** (`CC^T`, static `LDL^T`): **flat**. Symbolic has already sized every block, nothing
  grows, so one buffer with per-supernode offsets. 0.9 does the same (`FactorsStatic` allocates
  one array and points into it).
- **dynamic** (dynamic `LDL^T`): **VV**. Delayed pivoting grows a front at runtime by an amount
  symbolic never predicted, and the growth is local. 0.9 does the same (`FactorsDynamic`
  allocates one array per supernode).

**Static and dynamic, not flat and VV**, in the naming. The layout is a *consequence* of
mutability, not the thing being chosen, and 0.9 names its classes the same way. Flat-versus-VV
describes the bytes; static-versus-dynamic describes why.

**Status.** Settled: the table above, symmetry determined rather than chosen, operation-named
BLAS, static-flat and dynamic-VV. The objects are `NumFactor` and `NumFactorEngine`, the engine
taking a factorization type (`CC^T` first, then `LDL^T` static and dynamic) and a traversal
(left-looking, right-looking, multifrontal). Cholesky first, left- and right-looking first.

## 2026-07-12, Amalgamation: a second compression, and the third bug in 10.12

**Two compressions, and they are not the same algorithm with a knob.** Fundamental compression
contracts *paths*: it merges a supernode with its child when that child is an **only** child
and shares its pattern. Amalgamation contracts *stars*: it merges a supernode with any number
of its children, and will pay for the privilege in explicitly stored zeros, up to a budget.
They are orthogonal settings in the engine (`Supernodes::Fundamental` and an optional
threshold), as they are in both references, and both may run, fundamental first.

**Absent is not zero.** The threshold is a `std::optional<std::size_t>`, because "do not
amalgamate" and "amalgamate but pay nothing" are different instructions. At threshold zero it
still merges, and it merges strictly more than fundamental compression does.

**Which is the interesting fact: fundamental supernodes are not maximal.** The zero-fill
condition and the fundamental pattern condition are the *same* test. The only difference is
the only-child requirement, which exists to make supernodes paths, and hence unique, not
because merging would cost anything. Drop it and free merges appear. The smallest case is a
three-column star: fundamental gives 3 supernodes, amalgamation at threshold zero gives 2, at
no cost. On the grid of the notes it takes 7 supernodes to 6, column 3 joining the separator
for nothing.

**Uniqueness is what is being traded away, and it is worth stating plainly.** Where two
children could each merge for free, only one can (absorbing the first widens the front, which
prices out the second), so the algorithm must **break ties**, and a tie-break is a convention
rather than a theorem. 0.9's rule is: least fill, then largest front, then first in the child
list. That last clause is arbitrariness made deterministic, which is exactly what a canonical
algorithm never needs. Fundamental supernodes are unique because they refuse to make the
choice; amalgamation is not, because it makes it.

**The third 10.12 bug, and the most instructive.** 0.9 updates a parent's front size after it
absorbs children:

```
frontSizeArray[kk] += frontIncrement;
```

10.12 has that line, **commented out**, with the explanatory comment above it left intact:

```
// update the front size of supernode s2.
//numFrntIdxsArr[s2] += frntInc;
```

It matters. Parents are processed in increasing label order, and a child's label is below its
parent's, so by the time a parent is reached its children **have already been parents
themselves** and may have absorbed children of their own. Pricing a merge needs the child's
*current* front size. 10.12 reads the stale one, understating both the fill and the resulting
block, silently. And the cause is visible in its own source: the array is declared `const`, so
the update would not compile, and someone commented out the statement rather than fix the
constness. We restore it.

That is now three bugs found in 10.12 (the sibling-link copy in `SymbolicEngine`, the wrong
front-size operand in the fundamental merge test, and this), and two of them left a corpse in
place: a commented-out line whose comment still promises the behavior. The pattern is worth
naming, since it will recur: **10.12 transcribes 0.9's prose faithfully and its code
approximately.** Read the comments for intent; verify the code against 0.9.

**The unsigned subtraction is safe by theorem, and says so.** The merge cost contains
`|indexSet(K)| - |update(J)|`, which on unsigned types would wrap if it went negative. It
cannot: by the containment theorem `update(J)` is a **subset** of `K`'s index set, so this is
a set-difference size. We name it (`zerosPerCol`) and cite the theorem, rather than leaving a
subtraction that happens not to wrap. Unlike the fundamental merge test, it cannot be
rearranged into an addition, because the running fill total genuinely needs the difference.

**Testing what survives tie-breaking.** Since the partition is not canonical, the tests assert
only what is invariant under the tie-break: at threshold zero the stored-zero count is exactly
zero; the factor's true nonzeros are all still present; the supernode count never rises, and a
larger budget never raises it; the links, height and topological labeling stay valid. The
specific partition is deliberately not asserted, except on the star and the grid, where it is
forced.

## 2026-07-12, Friendship is a write grant, and reading needs no friend

**Every argument an engine takes falls into one of three cases, and the third is the only one
that is subtle.**

**Written: friend, and pass the object.** An engine is declared `friend` by exactly the object
it fills, and by no other. `ElmForestEngine` writes `ElmForest`, `SymFactorEngine` writes
`SymFactor`, `OrderEngine` writes `Permutation`. That is the complete list, and friendship exists
for no other purpose. Having granted it, pass the object rather than its fields: the engine can
reach them regardless, so enumerating them in the signature restricts nothing and only makes it
long. This is the subject of the entry above.

**Read: pass the object, through its public API, with no friendship.** We considered granting
the engines read friendship on `SparseMatrix` and `Permutation`, on the theory that it would
compact the signatures further. It compacts nothing. Every read already goes through an
accessor returning a `const&`, so there is no copy to avoid and no access to gain: friendship
would only let us write `A.mColPtr` where we write `A.colPtr()`, at the price of new friends on
objects whose headers say plainly that they have no writer. `SymFactorEngine` reads a dozen
fields of `ElmForest` and is not its friend, which is exactly right.

**Read, and only part of the object is needed: take that part.** `SparseMatrix` is the one
object an engine does not simply take whole, and the reason is not access but need. Ordering,
the elimination forest and the symbolic factorization are graph algorithms: they read a
sparsity pattern and never touch a value. So they take one. `SparseMatrix` offers two overloads,
one taking the matrix and one taking `colPtr` and `rowIdx`; the second is the implementation,
the first a one-line adapter over it. See the entry below.

**Friendship was never the constraint, and it is worth being precise about that.** The
field-taking signatures are not there because access is blocked, it is not. They are there
because a structural algorithm has no business asking for values. No amount of friendship would
change that, and granting it would not have shortened a single signature.

## 2026-07-12, Engine helpers take whole objects, or exactly the part they need

**Passing an object's pieces to a function that is already its friend restricts nothing.**
The engine can reach every field regardless; unpacking the arrays into the signature does not
narrow that access, it only spreads it across nine parameters. We had been writing helpers
that took the forest apart, on the theory that a signature listing `const&` inputs and `&`
outputs documents the data flow. It does not, and the cost of pretending otherwise was real.

**The signature could not carry the contract that actually matters.** `compressFundamental`
leaves `mFirstChild`, `mLastChild`, the sibling links and the roots **stale**: they still
describe the nodal forest and are wrong for the compressed one until rebuilt. A `const&`
parameter says "not written", which a reader takes as "still valid", which is exactly
backwards. The read set, the write set and the stale set always lived in a comment. We were
paying nine parameters for documentation we never received.

**Meanwhile the parameters were manufacturing a bug class.** `finalizeLinks` took five
`std::vector<std::int32_t>&` arguments in a row. Transposing `firstChild` and `lastChild`
compiles silently and builds a subtly wrong forest. Same type, same reference-ness, no help
from the compiler. `finalizeLinks(f)` cannot express that error at all.

**So: helpers take the object.** Across `ElmForestEngine` this took the private helpers from
38 parameters to 11 (`finalizeLinks` 9 to 1, `compressFundamental` 8 to 1, `computeHeight` 5
to 1). The gain is not brevity. It is that `compute()` now shows the shape of the algorithm,
parent links, child links, sizes, compress, child links again, height, instead of the
plumbing, and a structurally interesting fact like `finalizeLinks` running twice is legible
at a glance rather than buried in argument lists.

**The read/write/stale sets move to the comment, where they can be said properly.** That is
not a loss. A signature can express "I do not write this"; only prose can express "I leave
this stale, and you must rebuild it before anything else reads it".

**Three cases, and the third dissolves rather than trading off.**

- **Objects we write through friendship** (`ElmForest`, `SymFactor`): pass whole. This is what
  the entry is about, and it is also 0.9's shape, whose `compress_()` is a member taking no
  arguments at all. Our array-passing came from following 10.12.
- **Objects we only read, and take whole** (`Permutation`): pass whole. No friendship needed,
  the public accessors suffice, and it is free. This collapses `oldToNew, newToOld` into one
  parameter with no downside whatever.
- **Objects we only read, and of which we need only a part** (`SparseMatrix`): take the part,
  and adapt at the boundary. See below.

**Pass only the structural part of a matrix when only structural work is done.** Ordering, the
elimination forest and the symbolic factorization are graph algorithms. They read a sparsity
pattern; they never touch a value. So the implementation takes a pattern, and the overload
taking a matrix is a one-line adapter over it:

```
template<class Val>
bool compute(const SparseMatrix<Val>& A, const Permutation& p, ElmForest& f) const
{ return compute(A.colPtr(), A.rowIdx(), p, f); }         // adapter, one line

bool compute(const std::vector<std::size_t>&  colPtr,
             const std::vector<std::int32_t>& rowIdx,
             const Permutation& p, ElmForest& f) const;   // the engine
```

**This is about honest dependencies, not about C++.** A function's parameters should say what
it actually consumes. A structural algorithm that demands a matrix is lying about what it needs,
and it forces a caller holding only a graph to fabricate numbers to satisfy a signature that
will ignore them. Both overloads are public for that reason: the lower one is not an internal
shortcut, it is the **graph interface**. The whole structural pipeline,
`OrderEngine -> ElmForestEngine -> SymFactorEngine`, runs on a bare graph today: no `SparseMatrix`,
no scalar type, no numbers. The rule would be right in any language, and we should not restate
it as a fact about templates.

**The template mechanics are a consequence, not the motive.** In C++ the structural overload
happens to be non-templated, so it is compiled once rather than once per scalar type, and the
adapter is a forwarding line that inlines away. Pleasant, and worth nothing on its own. We first
wrote this entry as though `Val` were the reason, which inverted cause and effect and made an
honest interface look like a workaround for a language wart. It is not: the pattern-taking
overload would earn its place even if every type in the codebase were concrete.

**Adapt once, at the public boundary.** The adapter belongs on the entry point and nowhere
else. We first tried it on the private helpers too, and those adapters immediately became dead
code: once `compute` has unpacked the matrix, every helper below it already holds the pattern,
and a second layer has no callers. All three engines now have exactly one adapter each, on
`compute`.

This is not a new idea in the code, only a newly named one: `compute` already pulled `colPtr`
and `rowIdx` out of `A` and handed them to the helpers. We had simply not seen that the same
move, applied at the entry point, makes the pattern a public capability rather than an internal
convenience.

**`SparsityPattern` is a packaging question, not a design one.** `SparseMatrix` is two things
under one name: a pattern, and values indexed by it. Since the interface already passes the
pattern, introducing a type to name it would only replace two array parameters with one. The
interface is already right; whether the pattern travels as a named type or as two arrays is a
matter of convenience. Recorded as an improvement available, not a debt owed.

**Entry points are named `compute`.** One verb across every engine, since an engine's job is
to derive a fact from its inputs. `OrderEngine::order` was renamed to match. 0.9 calls them
all `run`, which is uniform but says nothing; 10.12 uses `ComputeElmForest` and
`ComputeSymbolic`, right about the verb but appending a noun the class already carries
(`ElmForestEngine::ComputeElmForest` stutters). `ElmForestEngine::compute(A, p, f)` says what
is computed three times over: in the class, in the arguments, and in the output.

## 2026-07-11, Choice objects are constructible; derived-fact objects are engine-filled

**The pipeline holds two kinds of object, and they deserve opposite rules.** A permutation
is a choice. Given a matrix there are n! valid permutations, and which one is wanted is
policy: AMD, MMD, nested dissection, the problem's own numbering, an ordering read back from
a file, an ordering composed from two others. An elimination forest is not a choice. Given a
matrix and a permutation there is exactly one correct forest, one correct symbolic
factorization, one correct numeric factor. They are derived facts, and any value other than
the derived one is simply wrong.

**So: choice objects are freely constructible, derived-fact objects are engine-filled only.**
`SparseMatrix` and `Permutation` are inputs in the sense that matters, they encode decisions
the caller is entitled to make, so they take a full public interface for building and setting
them. `ElmForest`, `SymFactor` and the numeric factor to come are written only by their engine,
through friend access, and expose read-only accessors. A caller-supplied elimination forest
is meaningless at best and silently wrong at worst, so nothing is gained by permitting one,
while a whole class of bugs is prevented by forbidding it.

**Input versus output is the wrong test, and would have misled us here.** A permutation is
`OrderEngine`'s output and `ElmForestEngine`'s input, so that framing gives no answer. The
question that does give an answer is whether the object is determined by what precedes it.
Determinacy, not position in the pipeline.

**A consequence worth stating: an ordering engine is optional.** `OrderEngine` is a
convenience for computing a good permutation; it is not the sole authority on what a
permutation may be, which is why `setOldToNew`/`setNewToOld` are public API rather than a
testing backdoor. This is not an innovation, 0.9 had `set`, `get`, `read`, `write`, `compose`
and the grid orderings, all public. We had merely under-ported the class, and the missing
setter first showed up as an inability to test the permutation maps against a known answer.
The friend grant to `OrderEngine` survives, but only as an optimization (it can skip
re-checking a bijection it just built), not as the mechanism by which a permutation gets
filled.

## 2026-07-11, Symbolic factorization: 10.12's design, 0.9's behavior

**The source-of-truth rule got its first real test here, and both halves of it fired.** The
union recurrence at the heart of symbolic factorization is identical in 0.9 and 10.12: the
same two-part union (the sparsity patterns of a supernode's own front columns, then the
update indices of its children), the same two skip tests, the same marker array. There was
no behavioral choice to make. Every difference between the two references was a matter of
shape, and in one case a matter of one of them being wrong.

**We took 10.12's naming, and it is a real improvement.** `s1`/`s2` for a child and parent
supernode, `lc`/`lr` for a local (factor-order) column and row, `ac`/`ar` for their
counterparts in the original matrix. 0.9 uses `jj`/`kk`/`i`/`di_`, which is hard to read and
harder to review. The 10.12 vocabulary also lines up with the `lc1`/`lc2` names already in
our `computeParent`, so the port reads consistently across engines.

**Front size is computed by counting the map, not by filling ones.** 10.12's
`rComputeNumIdxs` fills the front sizes with 1 unconditionally. We count how many columns
map to each supernode. Both run only on the nodal forest, where the map is the identity, so
both give all ones and the choice is stylistic: counting derives the value from the map
rather than asserting it. We originally justified this as generality, on the grounds that it
would stay correct after compression. That was wrong, and worth recording as such:
compression derives the merged front sizes itself, and the rest of the function is
nodal-only anyway (its update-size walk indexes by column, which coincides with supernode
only while supernodes are trivial), so it is never re-run afterwards. The map-count buys no
capability. This is also not a departure from 0.9, whose symbolic factorization recomputes
front sizes from the map with exactly this loop; we merely compute them one stage earlier,
on the forest, where they belong as an attribute.

**SymFactor stores its index sets flat, as 0.9 does, not as a vector of vectors.** This
follows directly from the flat-vs-VV decision below: the symbolic factor is written once
into a structure whose size the forest already knows (`frontSize + updateSize` per
supernode), so a flat buffer with per-supernode offsets is the right shape and stays
comparable to 0.9 buffer-for-buffer. 10.12 uses a vector of vectors here and we do not
follow it. One modernization within the flat layout: 0.9's `pointerToIndex` holds absolute
pointers into the index array, while we hold `std::size_t` offsets. Same layout, but the
offsets bracket each supernode's block as `[s]`/`[s + 1]` and match the convention already
used by `SparseMatrix`.

**SymFactor copies three links, not five.** The forest is doubly linked (parent, first and
last child, next and previous sibling), because `computeHeight` walks it backward. Symbolic
factorization only ever walks a child list forward, so it needs parent, first child and next
sibling, and that is exactly what both 0.9 and 10.12 store in their symbolic object. We
follow. The other two are cheap to add later if numeric factorization wants a backward walk.

**10.12's symbolic engine has a bug, and finding it is the point.** It builds its sibling
vector with `std::copy(elmForest.mLstChldVec.begin(), ..., nxtSblgVec.begin())`, copying the
last-child links into the next-sibling links. The factorization itself still comes out right,
because the recurrence reads the sibling links straight off the forest rather than off the
copy, but the sibling links stored in the resulting symbolic object are corrupt, and those
are what numeric factorization would later traverse. This is the concrete reason the rule is
"favor 10.12's design, verify against 0.9's behavior" and not "follow 10.12". We took 0.9's
correct three-link copy.

**The generality for compression is written but currently dead.** The recurrence unions the
patterns of every front column of a supernode, not just the first, because threshold-based
compression groups columns whose patterns are merely similar rather than identical. While
supernodes are trivial every supernode has exactly one front column, so that loop runs once
and the generality is untested. It is deliberately kept (both references keep it) and will
get its first real exercise when fundamental-supernode compression lands.

**Verification is against an independent oracle, not against ourselves.** The tests compare
every supernode's index set to a dense simulation of Cholesky fill (eliminate a column, make
its subdiagonal rows pairwise adjacent), which shares no code path with the forest or the
symbolic factorization. Agreement is therefore evidence rather than tautology. It runs under
both the natural ordering and AMD, so the permutation maps are exercised too, and it
incidentally reconfirms the forest's update sizes on the same matrices.

## 2026-07-11, Flat vs vector-of-vectors storage, and the descriptor view that spans both

*Partly measured (2026-07-12, `experiments/storage-options/`): the abstraction cost and the
layout cost are now numbers, not guesses, and one of this entry's original claims is superseded
by them. What remains provisional is the dynamic-numfactor storage, still recorded ahead of having
0.9's numfactor code, and the flatten-or-not question the measurement reopened. See Status.*

**Not an API question, a data-structure question.** Flat storage (one contiguous buffer
plus a per-column/per-supernode offset array, CSC style) versus vector-of-vectors (VV, one
inner `std::vector` per column/supernode) is not about the public interface. A clean API
could wrap either; that part is cosmetic. What matters is the friend algorithms: they reach
into the internals directly for speed (the friend-access decision), so the storage layout
is part of each algorithm's contract. Change the layout and the friends adapt; they cannot
go through a high-level per-element accessor in a hot loop.

**The dividing line is mutability.** VV earns its keep only when the structure is mutable
and edits must stay local: inserting into one column without reflowing every downstream
offset in a flat buffer. Flat wins everywhere the structure is write-once with sizes known
up front: one allocation, contiguous streaming, cheap offsets to hand to BLAS. The cost of
flat is that each column's size must be known before filling; the cost of VV is scattered
heap allocations (a cache miss at every column boundary) plus per-vector overhead.

**Symfactor and static numfactor are flat.** Symbolic never grows: the elimination forest
supplies `frontSize + updateSize` per supernode, so the flat buffer is sized exactly and
each supernode's block is written once with a local cursor. This is what 0.9 does:
`pointerToIndex` from the known sizes, then the union loop fills. A mutable symfactor would be
overkill; nobody patches a symbolic factor in place, one recomputes it after the matrix
structure changes. Static numeric factorization is likewise write-once into a known
structure. Both stay directly comparable to 0.9 buffer-for-buffer.

**Dynamic numfactor is why VV matters.** Dynamic LDL pivoting delays an unstable pivot and
passes it up to an ancestor, so that ancestor's front acquires columns symbolic never
predicted; its index set and value block grow at runtime by an amount unknowable until the
numerics run. That is genuine mutability, arriving through the numerics rather than through
matrix edits, and it is local, independent growth (one ancestor grows, its siblings do
not), which is precisely what VV does cheaply. Flat for dynamic needs contortions:
preallocate every front to its worst case (mostly wasted memory), or reflow the buffer on
each delay (O(nnz) per delay). Flat for dynamic is overkill; VV is the natural structure.

**A hybrid can keep the stored object flat.** The dynamic structure is unpredictable only
*during* elimination; once numfactor finishes, every delay has resolved and the factor is a
fixed structure with known counts, so it is flattenable, we just cannot size the flat
buffer until the numerics have run. The option on the table: dynamic (growable) working
storage inside the numfactor engine, then one flatten pass into the persistent flat factor.
The persistent factor is then uniformly flat whether it came from static or dynamic
numfactor, so the solve phase and the verification against 0.9 are identical either way, and
dynamism becomes an implementation detail of the engine's scratch. The copy is O(nnz(L)),
one pass, noise against numfactor's O(flops), and we factor once and solve many. Two honest
caveats: the flat factor is sized from the numfactor *result*, not from symbolic (delayed
pivoting adds fill, so symbolic stays the static lower bound and the dynamic flat factor is
a separate, post-hoc-sized object); and the flatten canonicalizes, it is copy-and-order
into the fixed front-then-update sorted layout the solve expects, not a pure memcpy. An
incremental flatten, supernode by supernode, releasing each dynamic front as it is copied,
avoids a peak that holds both representations at once.

**The descriptor view unifies the code across both layouts.** A VV is, in essence, an array
of (base pointer, length) pairs, one per column/supernode; flat presents the identical
shape by handing out pointers into its one buffer plus the counts. Write each friend against
that view (`{std::int32_t* ptr; std::size_t len}` for an index run, `{std::int32_t* rowIdx;
Val* val; std::size_t len}` for a matrix column) and the algorithm does not know which
layout produced the base pointer. 0.9 is already in this shape: `pointerToIndex[kk]` is an
absolute pointer into the flat index buffer and the consumer walks `(base, len)`. This is
the friend-access decision generalized from "member access" to "a layout-agnostic view over
members," and it applies to any CSC-style object, the matrix included.

**The view must not be a per-element virtual call. It need not be a template either.**
*(Written before the measurement; the second half of this paragraph is what changed.)* Our
worry was that unifying the layouts would reintroduce the non-inlinable per-element cost the
friend-access experiment clocked at ~6x, so we planned a compile-time abstraction: one algorithm
source, `template<class Storage>`, instantiated per layout. Two monomorphizations, not one binary
branching on layout inside the loop.

The measurement says even that is more machinery than the problem needs. Hoist the descriptors
out of the loop entirely, into three plain arrays (`rowIdxPtr[j]`, `valPtr[j]`, `len[j]`), and
the algorithm becomes an ordinary non-template function taking nothing but pointers. **One**
compiled kernel, no instantiation per layout, no virtual call, and it matches hand-written CSC.
The per-element indirection we feared never arises, because the indirection happens once per
column, not once per entry. So the rule survives in its negative half (never a virtual accessor
in the hot loop) and is superseded in its positive half (a template is sufficient but not
necessary).

The write side unifies too: because sizes are known up front for the flat cases, both flat and a
pre-sized VV inner accept the same `base[cursor++] = i` fill.

**One VV-only hazard: growth invalidates descriptors.** When a VV inner vector reallocates
on growth its base pointer moves, so any cached `{ptr, len}` goes stale. For flat and
symfactor this never happens (read-only after fill), so the view is stable. For dynamic the
rule is grow, then fetch the descriptor, then stream; never hold a descriptor across a
growth event. A `reserve()` per front to a symbolic-size-plus-margin reduces reallocations
but cannot guarantee zero, so the fetch-after-growth discipline still governs. Delays land
at assembly and the BLAS streaming comes after, so the ordering is natural.

**Measured, in `experiments/storage-options/`, and the result is stronger than predicted.**
Two sparse matrix classes (CSC and VV) holding identical content, one `MultiplyEngine`, and a
sparse matvec. Each class fills the same three arrays, `rowIdxPtr[j]`, `valPtr[j]`, `len[j]`.
CSC points into its single buffer; VV reads each inner vector's `data()`. The multiply takes
nothing but those arrays, and on an M4:

```
hand-written CSC (baseline)     1.362 ms
multiply(), CSC pointers        1.454 ms   1.07x
multiply(), VV pointers         1.499 ms   1.10x
multiply(), VV scattered        8.723 ms   6.41x
```

Bit-identical results. Two findings, and the first supersedes what this entry originally
claimed.

**The abstraction is free, and it needs no template.** We had assumed the unifier would have to
be a compile-time template, monomorphized per layout, to avoid a virtual call's cost. It does
not. There is exactly **one compiled multiply**, verified in the symbol table, and its signature
names neither matrix class:

```
T MultiplyEngine::multiply(unsigned long, int const* const*, double const* const*,
                           unsigned long const*, double const*, double*) const
```

It cannot tell CSC from VV because by the time it runs there is nothing left to tell apart. The
storage question and the algorithm question are **separable**: the layout decides where the
pointers come from, and nothing else. Two layouts therefore cost us zero duplicated kernels, not
"two monomorphizations" as this entry first supposed.

**The interface costs the same on both layouts**, 1.07x and 1.10x, a three percent spread. What
costs is *locality*: the two VV rows are the same class with the same content, differing only in
the order their inner vectors were allocated. **A flat buffer guarantees consecutive columns are
adjacent; a vector of vectors only ever borrows that from the allocator.**

**Two caveats on the 6.41x, both important.** It is *constructed*, shuffled allocation plus
interleaved spacers, engineered to remove the allocator's help entirely, so read 1.10x as VV's
structural cost and 6.41x as an upper bound. And it is hardware-dependent: the same code measures
8.87x on a machine with smaller caches. More importantly, **this kernel is the harshest possible
setting for a cache miss**: a sparse matvec does about two flops per element loaded, so a miss
shows at full price with nothing to hide behind.

**Which layout for which object, with reasons rather than symmetry.** Having both layouts is now
free, so the temptation is to offer both everywhere. Most of that would be storage nobody uses.

- **`A`: CSC.** It is built once and read forever, so VV buys it nothing but scatter. A VV `A`
  is defensible only for *structural* mutation, and that case is weaker than it looks. Value-only
  mutation (same pattern, new numbers, refactorize) simply overwrites `mVal` in place and needs
  no VV at all. Structural mutation invalidates the ordering, the forest, the symbolic
  factorization and the factor's size, so it forces the whole analysis phase again, against which
  rebuilding `A` in CSC is one `O(nnz)` pass and therefore noise. Which is also why one *batches*
  structural changes: each forces re-analysis anyway, so nobody applies them one at a time. And
  incremental assembly, the one place VV genuinely helps, has a better answer already: triplets,
  then convert once, with no per-column allocation at all.
- **`SymFactor`: CSC.** Write-once into a size the forest already knows
  (`frontSize + updateSize` per supernode). The textbook case for flat. VV buys it nothing.
- **Static numfactor: CSC**, for the same reason: write-once into a structure symbolic has sized.
- **Dynamic numfactor: VV.** Delayed pivoting grows a front at runtime by an amount symbolic never
  predicted, and the growth is local. This is the one place the algorithm genuinely mutates.

**Flatten-or-not is now a real fork, not a foregone conclusion.** This entry originally treated
the flatten-to-flat hybrid as settled. At 1.10x for a packed VV, "do not flatten; let the solve
read the VV factor" is defensible, and it saves an `O(nnz(L))` copy plus the peak memory of
holding both representations at once. Against that: the solve *streams* the factor, which is
exactly the cache-hostile shape the experiment measures, so it is the phase most likely to want
flat. Both options are now open, and the experiment is what opened them.

**The measurement we owe, and it is load-bearing.** The matvec is a toy for this purpose: about
two flops per element loaded, so a cache miss shows at full price. Numeric factorization is the
opposite, a dense front handed to BLAS level 3, `O(n^3)` arithmetic on `O(n^2)` data, where the
same miss amortizes over far more work. **We have predicted that VV-during-elimination is
therefore affordable, and we have not measured it.** That prediction is what the entire dynamic
design rests on, and it deserves its own study against a realistic front, not a sparse column
scatter. Until then, treat "VV is affordable in numfactor" as a hypothesis, not a result.

(Note also that the ~6x virtual-call figure this entry cites was carried over from the
friend-access experiment and is *not* tested by storage-options, which has no virtual path. The
question turned out not to arise, since the unifier needed no polymorphism at all.)

**Status.** Settled: CSC for `A`, for symfactor, and for static numfactor. Settled: VV for dynamic
numfactor's working storage. Open: whether the persistent factor is flattened, and whether a
VV-all-the-way numfactor is viable, both of which want the BLAS-3 measurement above.

**Confirmed against 0.9 (2026-07-13).** This entry predicted the split and flagged it as
unverified, pending sight of 0.9's numeric code. That code has now arrived, and 0.9 does exactly
what we reasoned it must:

- **`FactorsStatic`** sums `numberOfAllocatedEntries` over every supernode, allocates **one
  buffer**, and points `pointerToEntry[jj]` into it. **Flat.**
- **`FactorsDynamic`** runs `for (jj) pointerToEntryArray[jj] = new Real[jjEntrySize];`. **One
  allocation per supernode. A vector of vectors in all but name.**

So the split is a *port*, not a modernization: 0.9 already separates the two storages by
mutability, for the same reason we would. The prediction and the oracle agree, and the entry's
provisional status on this point is closed. What remains open is the *hybrid* question, whether
the dynamic factor is flattened afterwards, which 0.9 answers only implicitly (it keeps a
`FactorsDynamic` object, and `SolveEngine` reads it through the same abstract base), and the
BLAS-3 measurement that would tell us whether flattening is worth it.

**Amendment (2026-07-14): the taxonomy stands; the descriptor view in the title did not.** The
storage split above (CSC for the matrix, the symbolic factor, and the static numeric factor; a
vector of vectors for the dynamic factor) is the settled, 0.9-confirmed core and is unchanged. What
did not survive is the mechanism this entry pairs with it, the "descriptor view that spans both": a
bulk extractor that materializes `rowIdxPtr[j]` / `valPtr[j]` / `len[j]` arrays and feeds one
storage-blind compiled multiply. That was superseded by **direct access**. Each storage now exposes
per-column accessors (`rowIdx(j)` / `val(j)` / `colSize(j)`, one pointer or size per call), a
consumer templates over the storage and calls them at the moment of use, and there is no extractor
at all, so nothing is materialized, owned, or left to go stale. Direct is also faster (it streams no
extra arrays) and it is what the numeric factorization must use regardless, since a growing dynamic
factor would dangle any pointer extracted up front. So the array-valued names in this entry describe
the retired bulk form; the current interface is the three per-column accessors. See the
bulk-versus-direct entry above, and the storage-options README.

## 2026-07-09, Index types: `std::int32_t` IDs (NIL = -1), `std::size_t` offsets

Two kinds of integer, two types:

- **IDs**, a value that *names* a vertex/row/column/supernode, and may need a "none"
  marker, are **`std::int32_t`**, with sentinel **`NIL = -1`** (`constexpr std::int32_t`).
  E.g. `SparseMatrix::rowIdx`, the permutation maps, and the forest's parent / child /
  sibling / supernode-map arrays.
- **Offsets / counts / sizes**, row-pointers, `nnz`, dimensions, anything that indexes
  into or measures, are **`std::size_t`**. E.g. `SparseMatrix::colPtr`, `size()`.

**Why signed int32 for IDs.** The forcing function was the sentinel. The forest needs a
"no parent"/"no child" marker; on an unsigned `std::size_t` that can only be the max
value, spelled `static_cast<std::size_t>(-1)`, defined behavior but an ugly wraparound
smell, and easy to misuse in arithmetic. A signed `std::int32_t` gives a clean, obvious
`-1`. This also matches (a) the graph code, which already uses `int32_t` vertex IDs with
`static const int NIL = -1;` so companion arrays like `mate` hold `-1` naturally, and (b)
the vendored AMD/MMD, which are `int`-based, so `rowIdx` as `int32_t` largely removes the
`size_t -> int` conversion at that boundary. The cost is a ~2.1-billion index cap (int32 vs
size_t's range). Accepted deliberately: cleaner and more agile, and Oblio isn't targeting
matrices past 2^31 structural indices for now. If that changes, widen the ID type to
`std::int64_t` in one place.

**Why `size_t` stays for positions.** A position is never negative and can legitimately exceed
2^31 even when the index *count* does not: a row-pointer is a position into an `nnz`-length
array, and `nnz` outgrows `n`. So positions keep the full unsigned range and never carry a
sentinel. This mirrors the graph's `GeneralGraph` exactly: `idx` (the row-pointer) is `size_t`,
`adj` (the neighbor indices) is `int32_t`.

(CODING_RULES now calls these **indices** and **positions** rather than IDs and offsets, and
bans the phrase "index into a vector": `index` names a matrix or graph entity, and reusing it
for array access collides with that meaning exactly where the distinction matters. Same
distinction, sharper words.)

**No signed/unsigned friction, because loop counters stay `std::size_t`.** g2_csr
demonstrates the discipline: counters that enumerate positions are `size_t` (so they
compare against `.size()` cleanly), and `int32_t` appears only as *stored values* that
may be `NIL`, never as a loop variable. `size_t` is a safe superset of the non-negative
`int32_t` range, so viewing an index as `size_t` in a bounded, non-negative loop loses
nothing. The two types reconcile with explicit casts at exactly two crossings:
`static_cast<std::int32_t>(counter)` when storing a counter into an ID array (the
narrowing where the 2^31 cap lives), and `static_cast<std::size_t>(id)` when an ID
subscripts an array (widening; guard against `NIL` first if it could be a sentinel). Casts
are few and mark the ID <-> offset boundary.

Spelling: `std::int32_t` from `<cstdint>`, `std::size_t` from `<cstddef>`, std-qualified,
C++ headers, same rule as every other stdlib type. Applied first to `SparseMatrix`
(`rowIdx` -> `int32_t`, `colPtr` -> `size_t`); `Permutation` (maps -> `int32_t`) and `ElmForest`
(parent/etc. -> `int32_t` with a shared `NIL`) follow. The graph code uses bare `int32_t`;
it's the code to bring in line later, like the matching codebase was for `size_t`.

**Later note (2026-07-15):** the "~2.1-billion index cap" above is really two separate ceilings that
happen to meet at `INT32_MAX`, a dimension/index cap (representability, intrinsic to `int32`) and an
nnz(A) cap (the `int`-based ordering interface, external and A-only). nnz(L) is not capped. See the
ordering-constraint entry for the distinction.

## 2026-07-09, Friend grants write access; reads are public

The engine <-> data access rule:

- **`friend` = write access.** An engine befriends exactly the data class(es) it
  *produces/mutates*. There is no public mutation API for those internals, only the
  producing engine reaches the private members. E.g. `OrderEngine` writes
  `Permutation`; `ElmForestEngine` writes `ElmForest`. Friendship is declared by the
  data class (`friend class FooEngine;`), not by the engine.
- **Reads are public, everywhere, including hot paths.** Engines, tests, and users
  all read via the public const API. `friend` is *not* needed for reading.
- **This corrects an earlier overstatement** (the friend/BLAS entry below implied
  friend was needed for hot-path reads). The `experiments/friend-access/` study
  measured a *per-element, cross-translation-unit accessor call*, which can't inline,
  so it blocks vectorization, against direct friend member access; friend won ~6x.
  But that gap comes from calling the accessor *per element*, not from using a public
  accessor at all. Two public-accessor patterns match friend's performance exactly:
  (a) hand BLAS the block via a single `.data()` call, BLAS then owns the O(n^3) loop,
  so the one call is free; (b) for a hand-written element loop, bind the returned
  container once (`const auto& v = A.val();`) and loop over *that*, one non-inlined
  call total, then a vectorizable loop over contiguous memory. So the hot-path
  discipline is **"bind-once / pass the pointer," not "be a friend."**
- **Consequence for `SparseMatrix`:** A is input and nothing writes it (its
  construction path is TBD), so it needs *no* friends. Its earlier
  `friend class OrderEngine` is removed; `OrderEngine` reads A via
  `colPtr()`/`rowIdx()`/`size()` and remains a friend only of `Permutation` (which
  it writes). `ElmForestEngine` already followed this (friend of `ElmForest` only).
- **Numeric data classes** (`Factors`, ...) still befriend their producing engines for
  *writes*; their hot-path reads also go through public accessors with the
  bind-once/pointer discipline.
- **Exposure stance (pragmatic, not purist):** exposing internals read-only is fine;
  we don't design curated "won't-break" read APIs up front. A representation change
  already forces editing the friends (writers); public read exposure adds only the
  *tests* to that blast radius, cheap, and tests *should* feel a rep change. Curate a
  narrower public API later only if a structure's representation proves unstable. For
  the canonical structures we have (CSC, etree `parent[]`), the representation is the
  settled standard form, so exposing it by reference is low-risk.

## 2026-07-09, Sparse matrix storage: flat CSC, stored FULLY (both triangles)

`Matrix` (input A) stores its structure and values as flat **compressed sparse column
(CSC)**: three contiguous `std::vector`s, `colPtr` (size+1), `rowIdx` (nnz), `val`
(nnz), row indices sorted ascending per column. A symmetric matrix is stored **fully
(both triangles)**, each column holding its complete neighbor list plus the diagonal.

Storage layout (flat CSC vs vector-of-vectors) and triangle (full vs lower) are two
separate decisions:

**Layout, flat CSC**, across the three generations:
- **0.9** stored CSC via four manually `new`/`delete`d `Array*` pointers, the manual
  memory management the port removes.
- **10.12** modernized to `std::vector` but chose **vector-of-vectors** (one inner
  vector per column), RAII but the wrong layout: columns scatter across the heap.
- **PoC / port** use flat CSC (`mColPtr`/`mRowIdx`/`mVal`), vectors *and* contiguous,
  satisfying the "contiguous storage to BLAS via `.data()`" invariant.

**Triangle, full, not lower.** Both 0.9 and 10.12 store A **fully**: 0.9's
`getNumberOfNonzeroEntries() = size + 2*numOffDiagonals` (each off-diagonal in both
triangles) and its "storing A within the structure of A+A^T"; 10.12 has `SymmetrizeStrc`
and its etree reads full per-column neighbor lists. The **PoC diverged**, storing the
lower triangle only ("stored as lower triangle in CSC"). That divergence is what forced
every structural consumer to expand lower -> full first (the MMD path in `OrderEngine`, the
etree in `ElmForestEngine`), and the etree bug where lower-triangle input silently
produced an empty tree until expansion was added. **The port matches the oracle: A is
stored fully.** Consequences: structural phases (ordering, elimination forest, symbolic)
read each column's neighbors directly with no expansion (the etree's diagonal
self-skips via `lc1 < lc2`; MMD just strips the diagonal; AMD ignores it); it's the
faithful port (lower-triangle was a rewrite of the data structure); and it's the natural
substrate for a future **unsymmetric extension**, factor the symmetrized structure
A+A^T while carrying asymmetric values. Cost: ~2x off-diagonal storage for A, and the
numeric phase carries a redundant triangle. Accepted, matching 0.9/10.12; A is the input
and is far smaller than the factors, where the real memory lives.

Open for the port: the PoC exposed this as a public `struct` with a `fromCOO` builder
and a weak `isValid`. The modern `SparseMatrix` keeps the flat-CSC layout but is a
`class` with `friend` engines and a structural interface; 0.9 is the oracle for the
COO -> CSC assembly details (zero-diagonal insertion, duplicate merging, symmetrization),
10.12 shows which operations the solver actually calls.

## 2026-07-09, Two layers of modernization: rules prevent, clang-tidy catches

The coding rules and `.clang-tidy` are complementary layers catching different
failures at different times, not redundant work:
- **Coding rules** (CODING_RULES + CLAUDE invariants), preventive and broad. They
  shape code as it's written and cover what no tool can judge: port-verbatim
  discipline, friend -> BLAS, when to split a header, `.cpp`, `mFoo`, `std::size_t`.
- **`.clang-tidy` `modernize-*`**, a mechanical safety net, narrow but certain. It
  can't reason about intent, but within scope it catches the idiom slips that get
  through (a stray `NULL`, a `typedef`) and can auto-fix them.

So a `NULL` written by mistake is exactly what the rules *say* and the tool
*guarantees*, belt and braces, each doing what the other can't.

`modernize-*` is aligned with this project's purpose at the concept level: porting 0.9
(late-90s C++) forward *is* turning `typedef` -> `using`, `NULL` -> `nullptr`, raw loops ->
range-for, and so on across old code, exactly what that checkset does. So it can do
part of the mechanical modernization *for* you on each ported file (`--fix`), leaving
your attention on the algorithmic faithfulness no tool can verify. Hence the per-unit
workflow (in CLAUDE.md Process): port faithfully -> `clang-tidy --fix` -> verify vs 0.9.

## 2026-07-08, `experiments/` convention (runnable design studies)

`experiments/<name>/` holds self-contained, runnable studies that establish or
validate a coding standard before it is applied in the main tree. Each is its own
folder with its own sources, `Makefile`, and `README.md`; builds standalone
(`make test`); and is reference/teaching material, **not** part of the main Oblio
build. Executables carry the `_cpp` suffix and are gitignored.

Distinct from its two neighbors:
- `archive/`, frozen history (superseded PoC devlog, 0.9-analysis notes, old
  harnesses). Not maintained, not built.
- `examples/`, usage samples showing how to *call* the library (`examples/basic.cpp`, since
  renamed `examples/example_basic.cpp`).

An experiment answers a design question with code you can run and measure, then feeds
a decision here. Current studies: `template-instantiation/` (how to instantiate the
`Val` template, implicit vs plain/guarded explicit) and `friend-access/` (public API
vs `friend`-direct access, with timing). Experiments use the already-settled standards
(guarded explicit, `.cpp`, `mFoo`, `Oblio` namespace), so they double as worked
references for those standards.

## 2026-07-08, Numeric hot path: `friend` access, then BLAS (carried from 0.9)

The 0.9 design for numeric work, which the port preserves: **an engine reaches the
data's contiguous block via `friend`, then hands that raw block to BLAS** wherever
BLAS applies (gemv/gemm/syrk/trsm/potrf, via Accelerate on macOS). Not a new choice,
this two-step (`friend` -> BLAS) is how 0.9 does dense numerics.

This *is* supernodal numeric factorization: supernodes are dense blocks embedded in
the sparse structure, and the numeric phase is a long sequence of dense BLAS calls on
them, `syrk`/`gemm` for Schur-complement updates, `potrf`/`getrf` to factor the pivot
block, `trsm` for the off-diagonal solve, repeated per supernode across the whole
elimination. `FactorEngine` reaches each supernode's contiguous storage via `friend`
and passes the pointer straight to BLAS, no copy, thousands of times per factorization.
So `friend` isn't an optimization detail; it's the access mechanism the entire numeric
phase is built on. The `experiments/friend-access/` mat-vec is the single-block toy of
this pattern.

Important distinction, **the supernode blocks live in the factors, not in A.** `A`
(the input `Matrix`) is never handed to BLAS block-by-block; it's *read*, its structure
by the ordering and symbolic phases, its values once when they're scattered into the
factor. The dense blocks that BLAS operates on are created during factorization and
stored in `Factors`. So the `friend` -> BLAS hot path is specifically a `Factors` /
`FactorEngine` story. `A`'s storage (CSC) is chosen for a different reason: cheap,
cache-friendly *sequential column traversal* by the structural phases, plus being the
standard interchange format (what AMD/MMD expect). Both `A` and the factors favor flat,
contiguous storage over 10.12's vector-of-vectors, but for `A` the reason is streaming
structural reads, and for the factors it's the contiguous block BLAS needs.

Data classes (`Matrix`, `Vector`, `Factors`, `Symbolic`) expose a public,
bounds-checked API (`operator()`, `operator[]`) for reads, by all callers, on hot
paths too (see the "friend = write access" entry above: reads are public; the hot-path
discipline is bind-once / pass `.data()` to BLAS, not friendship). Engines befriend a
data class only to *write* it. 0.9 grants `FactorEngine` friend access into `Matrix`,
`Factors`, `Symbolic`; in the port that friendship is retained only where the engine
writes (the factor storage), and reads go through the public API.

Why (performance): the public-operator path is one non-inlined, cross-translation-unit
call per element (data-class body in its `.cpp`, loop in the engine's), which blocks
vectorization. Direct `friend` access fetches the raw block pointer once and walks
contiguous memory, which vectorizes, measured ~6x on Apple Silicon (M4/AppleClang),
~3x on x86/g++, over the API path. But the raw block should then go to **BLAS**, not a
hand loop: on the M4, Accelerate's `dgemv` ran the 2000x2000 mat-vec ~6x faster than
the vectorized hand loop (and ~36x over the API path), because it breaks past a single
core's bandwidth ceiling (multithreading + prefetch/blocking) that a hand loop can't.
The advantage only grows for the compute-bound O(n^3) kernels (`dgemm`/`syrk`/`trsm`/
`dpotrf`) the factorization leans on. So the hand loop is a fallback/baseline; where
BLAS applies, use it.

Why `friend`, not public getters: `friend` is *tighter* encapsulation, not looser, it
grants access to exactly the named engine classes, where a public `data()`/getter
exposes internals to the whole program. It honestly encodes that a data class and its
engines are one subsystem split for organization, not two modules talking through a
narrow API. Deliberate pragmatic choice over OO-purist accessors, and it's faster.

Consequence for porting: `friend` couples a data class to its engines, so the natural
port/verify unit is the *cluster* (e.g. `Factors` + `FactorEngine`), not the data class
in isolation, the friend boundary sets porting granularity.

Measured note: the gap is structural (non-inlined calls / vectorization), not assertion
overhead, toggling bounds-check asserts barely moves it. See `experiments/friend-access/`.

## 2026-07-08, Source extension: `.cpp` (headers stay `.h`)

Switch Oblio source files from `.cc` to `.cpp`; headers remain `.h`. All
extensions (`.cc`, `.cpp`, `.cxx`, `.C`) are identical to the compiler, so this is
convention, not correctness.

Why `.cpp`: cross-project consistency. The matching codebase uses `.cpp` (alongside
`.rs`, `.py`), one extension per language across the ecosystem, so Oblio being
`.cc` made it the odd one out. `.cpp` is also the more common choice in the wider
world (`.cc` is mainly the Google-style corner). Done now because the tree is still
scaffold + PoC with no ported units yet, so the rename is at its cheapest.

Blast radius (all mechanical, filename-level, no semantics): `git mv *.cc -> *.cpp`
across `src/`, `tests/`, `examples/`, `archive/`, `experiments/`; the CMake source and
executable lists; the manual build glob in CLAUDE.md / README (`src/*.cpp`); `.clang-format`;
the example `Makefile`s and the build-command comments in the example files. Headers and
object files are untouched, `#include`s point at `.h`, and `Foo.o` derives from the
source basename regardless of extension, so no `.o` reference changes. This makes
the rename strictly safer than the `exp` -> `ext` one, which touched `#include`s.

## 2026-07-08, Explicit instantiation over header-only templates (rationale)

Decision (already active in CLAUDE.md; this entry records *why*, which otherwise
lives only in the `experiments/template-instantiation/` example comments): Val-dependent classes keep a single
`Val` template, but their definitions live in `.cpp` files with explicit
instantiation for the supported scalar types, and headers carry declarations plus
`extern template`. Fuller treatments: `archive/oblio_modernization_notes.md` section "Why
explicit instantiation still works" (the Val-surface table) and
`archive/oblio-new-devlog.md` Session 3 (adoption + the link-failure proof). The
`_tpl` / `_ext` example files are the compact head-to-head (see naming below).

Mental model: a template is a recipe, not code, it generates code once a type is
plugged in. Header-only templates plug in *late and everywhere*: every translation
unit that includes the header re-runs the recipe for each type it uses, and the
linker discards the duplicates (N files x 2 types -> the same bodies compiled 2N
times). That was the real cost of 0.9's header-heavy templating, not the templates,
but the repeated late instantiation. Explicit instantiation is "one template,
applied early, once per type": the recipe runs exactly twice, in one `.cpp`, at
library-build time, and every other file links the existing result instead of
re-running the recipe. (Two mechanisms achieve that "instead of re-running",
declaration-only headers and/or `extern template`, see History below; they are
not the same feature and have different dates.) Generality is preserved (adding `float` is one
early application), but instantiation collapses from scattered-2N to centralized-2.

Key framing: build cost that scales with the number of scalar types is *incidental*
to the C++ compilation model, not *inherent* to supporting real and complex.
Nothing about "a matrix can hold real or complex" requires recompiling matrix code
in every includer, that only happens because header-only defers instantiation to
include time. Explicit instantiation removes the accident and keeps the capability.

The tradeoff, and why it's nearly free here: explicit instantiation gives up
instantiating *arbitrary* types at use sites, a consumer can't spin up
`Matrix<long double>` unless that line exists in the `.cpp`. For a maximally-generic
header library (Eigen) that's a real loss. For Oblio it isn't, because the scalar
world is closed and tiny: a type only makes sense if a dense BLAS/LAPACK kernel
exists for it, which bounds the space to BLAS's four, `float`, `double`,
`complex<float>`, `complex<double>` (s/d/c/z). We use two, might add the others.
Enumerating even all four is a handful of lines, far below the build cost of keeping
them implicit. So the one thing explicit instantiation costs is something a
closed-world numerical solver doesn't want anyway.

Bonus, and it matters for a verification-focused port: the `template class Foo<...>;`
lines *are* the list of supported types, in one place. Support becomes a declared,
reviewable fact rather than an emergent property of whatever anyone happened to
instantiate, and adding a type is a deliberate act that forces confronting whether
kernels and tests exist for it, exactly the gap the appendix flagged when complex
was "a new code path with zero test coverage."

Three mechanisms (not two, this frame makes the history click). Three distinct
tools, three natures, three dates:
- **Forcing**, `template class Foo<double>;`, emit all of `Foo<double>` in this
  TU. C++98.
- **Suppressing**, `extern template class Foo<double>;`, do *not* implicitly
  instantiate here; link it from elsewhere. C++11 (GCC extension earlier).
- **Definition-hiding**, not a keyword but a code-organization move: put member
  bodies in a `.cpp`, leave declarations in the header. A TU that can't *see* a body
  can't implicitly instantiate it. Works in every era.

Definition-hiding is the hinge, and it pairs with forcing. The three configurations:

| Case | Available | Approach | Build cost |
|---|---|---|---|
| 1 | neither forcing nor suppressing | Inclusion model: all definitions in headers, every TU re-instantiates what it uses, linker merges duplicates. No way to move bodies out and still get symbols. | High (~2N), unavoidable |
| 2 | forcing only (C++98) | Definition-hiding + forcing: bodies to `.cpp`, declaration-only headers, `template class Foo<double>;` in the `.cpp`. Other TUs see declarations only -> can't implicitly instantiate -> link the forced symbols. | Low |
| 3 | forcing + suppressing (C++11) | Case 2 still works and stays the choice; *additionally* you may keep bodies in headers and use `extern template` to suppress re-instantiation, needed only when definitions must stay header-visible. | Low |

Key insight: the big jump is **1 -> 2, not 2 -> 3**. Forcing is what unlocks the whole
technique (hide a definition, still guarantee the symbol). Suppressing is the
incremental step that only adds a second route for the case where you insist on
header-visible definitions. If you move bodies to the `.cpp` (Oblio does), you never
need it.

Precondition for all of it: an **enumerable** type set. For genuinely arbitrary
types you can't force or suppress anything (you don't know the list), so you're in
Case 1 regardless of language version. Forcing/suppressing are tools for closed type
sets, which Oblio's is (BLAS s/d/c/z), the same fact that makes the tradeoff above
nearly free.

Dates and 0.9: forcing is C++98, suppressing is C++11, so when 0.9 was written (late
90s) only Case 1 was portably reachable, header-only, inclusion model. That was the
*correct* choice for the era, and not because suppressing was missing (Case 2 doesn't
need it) but because template separate-compilation was a portability minefield then
(the `export` saga, inconsistent two-phase lookup, compilers disagreeing on
inclusion-vs-separation). Header-only-everything was the safe default. The modern
refactor applies matured portability; it does not correct a 0.9 error.

Where Oblio sits: current `ext` code is **Case 3**, bodies in `.cpp` (declaration-only
headers) *plus* `extern template`. But because the headers are already
declaration-only, the build win is really Case 2's (definition-hiding + forcing); the
`extern template` lines suppress nothing here (no visible header body to instantiate),
so they are documentation, not mechanism, a header annotation of intent, latent
unless a body is later (wrongly) added to a header. See the naming note below for the
full implicit / plain explicit / guarded explicit framing. So Oblio's pattern was achievable in C++98; C++11 was
not strictly required.

Template-instantiation example, naming (one algorithm, dense mat-vec, built
three ways). File names are `<Class><Variant>` with no separator, e.g.
`MatrixImplicit.h`, `MatrixPlainExplicit.{h,cpp}`, `MatrixGuardedExplicit.{h,cpp}`
(same for `Vector`, `MultiplyEngine`). The three variants:
- `Implicit`, body in the header; instantiated implicitly per TU. Stands in for
  what 0.9 effectively was.
- `PlainExplicit`, bodies in the `.cpp`, header signatures only (explicit
  instantiation, forcing only).
- `GuardedExplicit`, plain explicit + `extern template` in the header. The
  pattern used in the real tree.

Conceptual framing (two axes):

- **Axis 1, where the body lives (implicit vs explicit).** *implicit* = body in
  the header (`.h`), instantiated implicitly per translation unit -> `Implicit`.
  *explicit* = body outside the header, forced in the `.cpp`; header carries
  signatures only -> `PlainExplicit` and `GuardedExplicit`.
- **Axis 2, applies only within explicit; guarded vs plain.** *plain explicit* =
  `.cpp` bodies, declaration-only header, nothing more -> `PlainExplicit`.
  *guarded explicit* = same, plus `extern template` in the header -> `GuardedExplicit`.

  The three named layers:
  - **implicit**, body in `.h` (`Implicit`)
  - **plain explicit**, bodies in `.cpp`, header signatures only (`PlainExplicit`)
  - **guarded explicit**, same as plain explicit, plus the guard (`GuardedExplicit`)

  Important: the guard is NOT a sub-kind of "more correct" explicit, plain explicit
  is a complete, valid design. `extern template` only ever acts on *visible* header
  bodies (the implicit failure mode); in a declaration-only design there is nothing
  for it to suppress, so it is pure documentation, a header annotation reading
  "instantiated elsewhere," aimed back at the implicit branch it guards against. It
  gains mechanical effect only if someone reintroduces a header body (which the
  "definitions live in `.cpp`" invariant in CLAUDE.md forbids). So "guarded" (a
  reminder/guard), not "suppressed" or "enforced", in this design it suppresses
  nothing and enforces nothing; the invariant does the enforcing.

All three are built and tested together via the example's `Makefile` (`make test`)
against one shared source, `test_multiply.cpp`, they must produce identical
results, and the plain explicit and guarded explicit variants share the same link-failure behavior when
their `.cpp` files are omitted (empirical confirmation that with declaration-only
headers `extern template` suppresses nothing, it is documentation, not mechanism).
Selector macros: `OBLIO_TI_IMPLICIT` / `OBLIO_TI_PLAIN_EXPLICIT` /
`OBLIO_TI_GUARDED_EXPLICIT`.

Naming history: suffixes were once `_tpl`/`_exp`/`_ext`, where `_exp` inaccurately
labeled the extern-template variant. Renamed in two steps, first `_exp` -> `_ext` to
free `_exp` for the genuine forcing-only variant, then all three to the conceptual
`Implicit`/`PlainExplicit`/`GuardedExplicit` once the two-axis framing settled.

## 2026-07-08, Matrix naming: explicit `SparseMatrix` / `DenseMatrix`

Rename the sparse matrix type from `Matrix` to **`SparseMatrix`**; keep
**`DenseMatrix`**. Both are plain concrete types with **no shared base class**.

Why: the old `Matrix` (implicitly sparse) + `DenseMatrix` (explicitly marked) is a
half-committed convention, it privileges sparse as the unmarked default but marks
dense as the exception, so a reader must already know "unmarked = sparse" to parse
it. Going fully explicit removes that. It also makes matrix naming consistent with
the graph naming used elsewhere (`GeneralGraph` / `BipartiteGraph`, no bare
`Graph`). Note the typed-library mainstream (Eigen, scipy) actually defaults the
*other* way (unmarked = dense, `SparseMatrix` marked), so `Matrix` = sparse
actively misleads anyone arriving from there. The sparse-first precedent
(CSparse/SuiteSparse `cs`) is the one respectable counter-argument, but explicit
wins on cross-codebase consistency and reader-independence.

No speculative base `Matrix` interface. If something ever needs to consume sparse
and dense polymorphically, add the interface then, on top of the two concrete
types, not before a caller forces it.

Scope: this is a rename (wrapping/API), not an algorithm change, it's on the
port-and-modernize track, not the rewrite track. Do it as one deliberate
mechanical pass **before** porting proper, since it's cross-cutting (every `friend`
decl, FactorEngine, SolveEngine, tests name `Matrix`) and only gets more expensive
as code solidifies around the name. **Oracle mapping: 0.9 `Matrix` <-> modern
`SparseMatrix`**, record this so output comparisons against 0.9 stay unambiguous.

## 2026-07-08, Minimal abstraction; containers are the structure that matters

Design stance for the port: concrete types, minimal OO ceremony, `std::vector` as
the spine. Don't build base classes, single-implementation interfaces, or
inheritance ahead of an actual polymorphic caller (YAGNI). A direct solver is a
pipeline of concrete transforms, not a class hierarchy.

One guard against a tempting over-correction: "AI makes code cheap to generate, so
structure matters less" is only half true. Generation got cheap; **verification did
not**, and this project is the proof (every PoC bug was cheap to write, expensive
to trust, and fixed only by slow comparison against the 0.9 oracle). Structure's
real job was never to save typing; it's to keep code readable and checkable and to
localize where a bug can hide. That matters *more* when code is machine-generated,
because the bottleneck shifts onto review. So: drop OO flavor that only served
human writing time; keep the structure that serves verification (clear module
boundaries, one concern per unit, the friend seams that let a supernode block be
diffed against 0.9).

This is why "proper containers everywhere" (`std::vector` over `Array`) is not a
style preference but the load-bearing invariant: a vector carries its own size (no
drifting length variables), self-frees (no leak/double-free surface), bounds-checks
under sanitizers, and hands clean pointers to BLAS via `.data()`. Each removes a
bug class you'd otherwise verify the absence of by hand. For this codebase, the
container discipline *is* the architecture that matters.

## 2026-07-07, Align with standard project files; adopt clang tooling

The doc set maps to established conventions rather than being bespoke:
CODING_RULES.md ~ a style/conventions guide, DESIGN_DECISIONS.md ~ a lightweight ADR log
(single-file variant of the Nygard/ADR pattern), CLAUDE.md ~ the agent-instructions
file (AGENTS.md is the emerging cross-tool equivalent, if we ever run more than
Claude Code, keep content in AGENTS.md and make CLAUDE.md a one-line `@AGENTS.md`).
PORTING_LEDGER.md stays bespoke, it's specific to a port, no standard analog.

Renamed CPP_RULES.md -> **CODING_RULES.md** and made it language-general (Rust is a
likely future scope), with per-language sections.

Adopted **`.clang-format`** and **`.clang-tidy`** so mechanical rules (`nullptr`,
`using`, `enum class`, `= delete`, narrowing) are tool-enforced instead of written
as prose an agent may skip. `.clang-tidy`'s `modernize-*` family does part of the
port mechanically. The prose rules doc now holds only judgment calls tools can't
check.

Stubbed ahead of need (deliberate, not yet load-bearing): CONTRIBUTING.md and
CHANGELOG.md ("Keep a Changelog" format). CONTRIBUTING fills out on going public
(and becomes the canonical build/test source then); CHANGELOG gets its first real
entry at the first tagged release, which requires settling the tree's version
identity (Oblio 11 vs a fresh 1.0), tied to the working-tree question below. Still
outstanding: LICENSE (bundled AMD is BSD-3-clause, so a project license decision is
eventually needed).

Tension worth remembering: canonical ADR is many small numbered files, which suits
humans browsing history but fights the "always in context" goal, a growing pile
can't all stay loaded. The single-file log + distilled always-on layer is a
deliberate adaptation of the standard, not an oversight.

## 2026-07-07, What stays always in context vs on demand

Goal: coding rules and active design constraints present in every session, without
letting context grow unbounded.

Mechanism: only CLAUDE.md and its `@`-imports load every session; everything else
is read on demand. `@import` costs the same context as inlining, it just keeps
the file separate and editable.

Split:
- **Always-on** (in CLAUDE.md or imported): invariants (inline), conventions
  (`@docs/CODING_RULES.md`), and a distilled *Active design constraints* summary.
- **On demand**: the full DECISIONS log (this file, grows over time, so importing
  it would erode context), PORTING_LEDGER (porting-specific; read after a gap),
  README, archive history.

Rule of thumb: the always-on set is the distilled essence; narrative, dates,
history, and open questions stay here. Keep CLAUDE.md well under ~200 lines, past
that, adherence drops. If the always-on set grows, distill harder rather than
importing the growing logs.

## 2026-07-07, Invariants live in CLAUDE.md, not CODING_RULES.md

The C++ **invariants** (port-verbatim, `std::vector` default, no signed -> unsigned
index slips, `.data()` to BLAS) are written directly in CLAUDE.md. Only the
**conventions** (style preferences) stay in CODING_RULES.md.

Why: Claude Code auto-loads CLAUDE.md every session (directory walk up to repo
root, concatenated) but does **not** auto-load other files, a pointer from
CLAUDE.md to CODING_RULES.md does not pull its contents in. So a rule meant to be
active every session must physically live in CLAUDE.md, or it isn't loaded until
Claude happens to read the file it's in. Invariants are always-on; conventions
are fine to read on demand.

Do not "consolidate" the invariants back into CODING_RULES.md, that silently
disables them.

## 2026-07-07, Documentation structure

Established this doc set: CLAUDE.md (operating contract + index),
PORTING_LEDGER.md (per-unit status), CODING_RULES.md (conventions a linter can't
enforce), DESIGN_DECISIONS.md (this log), plus `.clang-format` / `.clang-tidy` for
tool-enforced mechanical style. The existing md files, devlog, modernization
notes, appendix, README, are kept as history/rationale for now; decide later
whether to fence or retire them.

Why: after a multi-month gap the project needs to be reconstructable from three
things (history, layout, decisions) without re-deriving them. CLAUDE.md holds
the operational index (it's what Claude Code loads each session); the reasoning
lives here so CLAUDE.md stays lean. The two live in different files on purpose and
should not duplicate each other.

(This entry documents a list that includes DESIGN_DECISIONS.md itself. The self-reference
is intentional, not an oversight.)

## 2026-07-07, OPEN: which tree is the working copy

Unresolved; resolve before the first port step and record the outcome here.

Candidates:
- **(a) Fresh port from 0.9** into a new tree, treating the PoC the way 10.12 is
  treated, learn from, don't build on. "Replacing Array gradually" only applies
  to a tree that still has `Array` (0.9 or 10.12), which points here.
- **(b) Continue the PoC tree**, where Array -> vector is already done, in which
  case the remaining work is different (coverage, finishing, cleanup), and this
  is not really an Array migration.

The stated plan ("port carefully from 0.9, one function at a time, replace Array
gradually") reads as (a). Confirm.

## 2026-03-07 (PoC), Choices carried from the proof of concept

Recorded for continuity. The PoC was exploratory, so revisit each on its merits
rather than inheriting it unquestioned.

- **One `Val` template** instead of 0.9's separate `*Real.h` / `*Complex.h` file
  pairs. Cholesky treated as Hermitian input; LDL^T as complex-symmetric.
- **`std::vector` storage** instead of the hand-rolled `Array`.
- **Explicit template instantiation** for `double` and `std::complex<double>`,
  headers declare, `.cpp` files define and instantiate. Faster builds; `float` /
  `long double` remain one line away.
- **Namespaced `include/oblio/` headers**, declarations only.
- **Flat `src/`**, all sources, including `Mmd.cpp` / `Amd.cpp`, directly in `src/`.
