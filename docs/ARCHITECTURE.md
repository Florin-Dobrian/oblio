# Architecture

How the pieces fit, as they are today. This is the *living* description: what the units are, what
they promise each other, and which of those promises are load bearing. It is meant to be edited
whenever the code changes, which is what separates it from its neighbors.

- **README** is the outside view: what Oblio does and how to build it.
- **CLAUDE.md** holds the invariants, the rules that must not be broken while working.
- **DESIGN_DECISIONS** is history: dated entries recording what was chosen and why, superseded
  rather than rewritten. It explains how we got here; this file explains where here is.
- **PORTING_LEDGER** tracks unit by unit what has been ported and checked.
- **experiments/** are frozen studies. Each answers one question with a measurement and is not
  maintained afterwards, so nothing here should depend on their contents staying current.

Sections are added as they earn their place. Each describes one structure or contract that a reader
would otherwise have to reconstruct by grepping.

## The storage-agnostic factor contract

The numeric factor comes in two storages. `NumFactorStatic` holds flat index and value buffers with
an offset per supernode; `NumFactorDynamic` holds one index vector and one value vector per
supernode, so that a delayed pivot can expand a front without moving its neighbors. Dynamic LDL
requires the second; everything else is cheaper in the first.

**There is deliberately no common base class**, so the thing that lets one algorithm serve both is
not a type. It is a set of accessors that both classes spell the same way, with the same signatures,
differing only in body. The contract is structural rather than nominal, which is exactly why it
needs writing down: no header declares it, and a reader would otherwise have to infer it from what
compiles.

Five of the shared accessors are the same field read on both classes. Only the two *lookups* differ,
and they differ in the way the storages differ: an offset into a flat buffer against a vector that
already holds the pointer.

| accessor | returns | static | dynamic |
|---|---|---|---|
| `size()` | `std::size_t`, the matrix order | `mSize` | `mSize` |
| `snodeSize()` | `std::size_t`, how many supernodes | `mSnodeSize` | `mSnodeSize` |
| `factorization()` | `Factorization`, which method produced this | `mFactorization` | `mFactorization` |
| `frontSize(jj)` | `std::size_t`, supernode jj's own columns | `mFrontSize[jj]` | `mFrontSize[jj]` |
| `updateSize(jj)` | `std::size_t`, its update rows | `mUpdateSize[jj]` | `mUpdateSize[jj]` |
| `nodeToSnode()` | `const std::vector<std::int32_t>&`, which supernode owns a node | `mNodeToSnode` | `mNodeToSnode` |
| `nodeIdx(jj)` | `const std::int32_t*`, where jj's node indices start | `mNodeIdx.data() + mSnodeNodeIdxPtr[jj]` | `mNodeIdx[jj].data()` |
| `val(jj)` | `const Val*`, where jj's value block starts | `mVal.data() + mSnodeValPtr[jj]` | `mVal[jj].data()` |

### Who reads what

The contract is not one set. Each consumer uses the part it needs, and the differences are
informative: the solve needs to know which factorization produced the factor, because that decides
whether it conjugates, while the traversals never ask, because the engine already knows what it is
computing. Conversely the traversals need the node-to-supernode map to find an ancestor, and the
solve never does. The column is the static traversals (`factorStaticLeftLooking` and its
right-looking twin), which are the pair templated on the factor; dynamic pivoting reaches further,
for the expansion and contraction verbs below.

| | `SolveEngine` | the static traversals in `NumFactorEngine` |
|---|---|---|
| `size`, `snodeSize`, `frontSize(jj)`, `updateSize(jj)`, `nodeIdx(jj)`, `val(jj)` | read | read |
| `factorization()` | read | not used |
| `nodeToSnode()` | not used | read |
| non-const `nodeIdx(jj)`, `val(jj)`, `numPerturbations()` | not used | write, through friendship |

**Reading is public; writing is not.** The const overloads are public, so a consumer that only
reads the factor needs no friendship, which is why `SolveEngine` is a plain reader. The non-const
overloads that hand out a writable value block are private, reached only by `NumFactorEngine` as a
friend. Name lookup gathers every overload before access is checked, so a non-const object outside
the friend context selects the private overload and fails to compile; `std::as_const(nf).val(jj)`
is the idiom for reading through such an object.

**Inside the friend context it is unnecessary, and the engine does not use it.** `NumFactorEngine`
may select either overload, and the writable one converts to a `const` local anyway, so the cast
would change nothing a reader can observe. What says "this block is only read here" is the type of
the local, `const Val*` in left-looking, where `jj` is an already-factored descendant, against
plain `Val*` in right-looking, where both blocks are written. Sprinkling `as_const` over some of
those sites and not others makes the asymmetry itself the thing a reader has to explain.

**The dynamic factor's expansion and contraction verbs are private and one-sided.** `expandNodeIdx`,
`resetVal`, `expandVal`, `swap` and `contractVal` exist only on `NumFactorDynamic`, because only its storage
makes them cheap. That asymmetry is the same rule `experiments/storage-options` arrived at for the
matrix pair, where `setColumn` is absent from the flat class by design: an object offers what its
storage makes cheap, and nothing else.

### The life of a delayed column

Two flows run over the elimination forest during dynamic LDL, and they have different reach. An
**update** goes from a supernode to any ancestor holding one of its update rows, which may be far up
the tree and may skip levels; this is the flow the traversals disagree about, left-looking pulling
and right-looking pushing. A **delay** goes from a child to its parent, one edge, never further.

It helps to see the asymmetries between left- and right-looking in two layers. The first layer is
already present in the *static* drivers, where there are no delays at all: left-looking pulls updates
from descendants and factors last, right-looking factors first and pushes updates to ancestors, and
that difference forces right-looking's A-assembly prepass (a push lands on a front not yet reached, so
A must be there first). These asymmetries are essential to the two orderings and visible in the
simplest code. The second layer appears only in the *dynamic* drivers, and it is the delay flow. The
delay flow is driven by the child-to-parent relationship, which is not the relationship the update
flow follows, so laying it on top of the update flow adds asymmetries that have nothing to do with
left-versus-right: they come from fitting a tree-shaped migration into a loop whose shape was chosen
for updates. Much of the dynamic drivers' apparent divergence is this second layer, and separating it
from the first is what keeps the two drivers legible.

Only one flow is contested, so the delay machinery is identical in both drivers, which fell out
rather than being designed. A delay can still travel further than one edge, but only by being
delayed again: a column passed up to the parent becomes an ordinary front column there and may fail
to pivot a second time. At a root there is no next hop, and that is an error rather than a numeric
failure.

From this the storage behavior follows, and it is worth stating because it bounds the allocation:

| | when | how often |
|---|---|---|
| `expandNodeIdx` | the supernode is reached, if its children delayed | at most once, never undone |
| `resetVal` / `expandVal` | the same moment | at most once |
| `contractVal` | the *parent* is reached, if this supernode delayed | at most once, always after |

Once and not more, because delays arrive only from children and a supernode is reached once. "At
most" rather than "exactly" for two reasons that are different in kind: structurally, a leaf never
expands and a root never contracts; numerically, whether a supernode delays at all depends on the
*values*, not the pattern. The same matrix under a different ordering delays a different set of
columns, which is the whole difference between static and dynamic pivoting.

The order cannot interleave, since a parent's index is above its child's, so a block is at its
largest between its own expansion and its parent's collection. That interval is the storage high-water
mark.

**Growth counts and contract counts need not match**, and where they diverge says something about the
forest: several children may delay into one parent, giving one expansion and several contracts. A chain
forest, which a banded matrix under Natural ordering produces, makes them equal; a branching one
does not. Measured on a saddle-point matrix under AMD: 15 supernodes grew, 24 shrank, identically in
both traversals.

**The index set never contracts.** `contractVal` truncates the value block's columns and keeps every
row, because a delayed column remains a row of L: its entries under this supernode's pivots are
genuine. That is the height invariant, `frontSize + delaySize + updateSize`, constant
throughout. It also makes the three regions a partition, so every row belongs to exactly one flow,
which is what the traversals rely on when they advance past `frontSize + delaySize`
rather than past the front alone. Advancing by the front alone would push the delayed rows to the
parent as an update *as well as* handing them over as delays, and the parent would double-count
them silently.

**Where a delayed column gets its values.** `factorStaticSupernode` pivots every front column, so
nothing is left behind. `factorDynamicSupernode` pivots only the columns it accepts, its *post-factor*
front, but eliminating those pivots applies the Schur update across the whole trailing block, and the
delay region sits in that block. So a delayed column is brought fully current during the factor call;
it is simply not used as a pivot there. The update kernel never touches it: `updateDynamicSupernode`
reads only past `frontSize + delaySize`, so it carries the update area to ancestors and leaves the
delay region alone. A delayed column is therefore factored-into where it sits, migrated intact to the
parent by `assembleDelay` and `contractVal`, and pivoted there. It is never short-changed and never
double-updated: the factor step updates it once, the update kernel never does, and the parent pivots
it once. This falls out rather than being arranged, since the pivot loop's trailing update does not
care whether a trailing column is destined for delay or for update, it updates the whole tail either
way.

**This set has expanded.** The solve reads `delaySize(jj)` and `pivotType()` from the
dynamic factor, the first because a delayed column leaves its row behind so the leading dimension is
`frontSize + delaySize + updateSize`, the second because a 2x2 pivot puts D's
off-diagonal where L's first sub-diagonal entry would be. Both are dynamic-only, which is why
`SolveEngine`'s three passes are paired: `forwardStatic` and `forwardDynamic`, and so on, the static
three templated on the factor and the dynamic three naming `NumFactorDynamic` outright. Same split,
same reason, as the traversals in `NumFactorEngine`.

### The shape of an update: the (jj, kk) edge

An update from `jj` to `kk` is a rank-`f` outer product, and its two other dimensions, the width `w`
(`jjKkWidth` in the code) and height `h` (`jjKkHeight`), are properties of the **(jj, kk) edge**, not
of either supernode alone. Only `f`, jj's post-factor front, is jj's on its own: it is the number of
pivots eliminated at `jj`, and so the rank of every update `jj` makes.

`w` is how many of jj's rows land in `kk`. Those rows are nodes owned by `kk`, and a node owned by
`kk` is a front column of `kk`, so the same `w` nodes are at once a horizontal band of jj's rows and a
subset of kk's columns. That identity is why an update is possible at all: one index is a row in jj's
block and a column in kk's, the structure being symmetric. Change `jj` or change `kk` and `w` moves,
which is what makes it a two-supernode quantity rather than either one's.

`h` is that band plus everything below it in jj: all of jj's rows from the band's start to the end of
jj's index set. So `w` measures a subset of kk's columns and `h` a subset of jj's rows, each the set
it is most naturally a subset of, and `w <= h` always, with equality exactly when nothing of jj
reaches past `kk`.

Both are contiguous, which is what lets the kernel take one start position and two lengths instead of
an index list. jj's index set is sorted and the supernodes partition it in increasing order, so kk's
band is a single run and `h` is a suffix: the `w` band on top, jj's higher-reaching rows beneath. The
update block is `h` by `w`. Its top `w`-by-`w` square is a block against its own conjugate transpose,
symmetric, filled by `gemmLower`; the `(h - w)`-by-`w` rectangle below is jj's higher-reaching rows
against kk's rows, not symmetric, a plain `gemm`.

So the unit of work is the edge: subsets of jj's rows updating subsets of kk's columns, `f` summed
over. The kernel forms only the `f`-by-`w` slice `D L21_kk^H` for this one pair, which is why it runs
once per (jj, kk) the traversal visits rather than once per `jj`, and why the same function serves both
traversals. The edge is symmetric to which end one stands at: `kk` held fixed while its descendants
arrive (left-looking), or `jj` held fixed while its ancestors are visited (right-looking).

### The life of an update

The section above followed a delayed column, which travels one edge. An update travels differently:
from a supernode to *any* ancestor holding one of its update rows, possibly several of them, and
possibly skipping levels. How that gets scheduled is the one place the two traversals genuinely
disagree, and it is worth setting down because the mechanism is not obvious from either driver.

**Two independent flows share these drivers.** The update flow moves a supernode's update area to its
ancestors, driven by the update (descendant-to-ancestor) relationship. The delay flow moves a
supernode's unpivotable columns to its parent, driven by the child-to-parent relationship. They are
not aligned: an ancestor that receives an update is usually not the parent, and the parent that
receives a delayed column is usually not an update target. The two lives in these two sections are
that pair of flows. A useful consequence is that the update itself is exactly what the static twin
does, in static and dynamic alike: the update area out of `jj` updates `kk`, and the delay machinery
does not touch that path. In the dynamic drivers this shows up as the update being *contract-agnostic*.
Left-looking updates from `jj` after `jj` has been contracted, right-looking before; the result is
identical, because the update reads only `jj`'s update area, which is disjoint from both delay regions
(`jj`'s own delayed columns, which the update walk starts past, and `kk`'s inbound delays, which
arrive by `assembleDelay` instead) and which never changes size. The delay flow runs alongside the
update flow without perturbing it.

**The structural fact underneath everything here** is the elimination tree's absorption property:
every update row of `jj` also appears in `parent(jj)`'s index set. Measured across grids and random
matrices, nodal and fundamental, roughly 6000 update rows with no exception. So if `ii` has a row in
some distant ancestor `kk`, then `parent(ii)` has that row too, and so does every supernode above it
on the way to `kk`. The relay therefore always has a next hop and never strands an update.

**Note the direction, because the convenient-looking converse is false.** Absorption says `ii`'s
update rows are a *subset* of its parent's, and a proper subset is normal. It does not say `ii`
updates everything its parent updates, so "`ii` updates `jj`, `jj` updates `kk`, therefore `ii`
updates `kk`" does not hold. A four-node chain is enough to break it, nodal, etree `0->1->2->3`,
no branching anywhere:

```
  snode 0  nodeIdx=[0 1 3]   updates: 1 3      <- skips 2
  snode 1  nodeIdx=[1 2 3]   updates: 2 3
  snode 2  nodeIdx=[2 3]     updates: 3
```

`0` updates `1` and `1` updates `2`, yet `0` does not update `2`: it has no entry in row 2, and no
fill creates one, since that would need a path from 0 to 2 through nodes numbered below 0. Skipping
comes from sparsity within a column, not from branching in the forest.

**Right-looking delivers data.** Standing at `jj`, it walks `jj`'s own index set, finds each ancestor
in turn, forms the update and assembles it immediately. Having done the work it stores nothing, so
the full ancestor walk costs it only the walk.

**Left-looking relays intent.** Standing at `kk`, it needs to know which descendants will update
`kk` before `kk` is reached, and that is a record it must carry. It keeps one queue per supernode,
`descendantUpdateQueue[kk]`, and a supernode sits on exactly one queue at a time, the next ancestor
it must update, hopping to the following one each time it delivers an update. `nextUpdateSp[jj]` is
the bookmark saying how far through `jj`'s index set it has got.

This is **bookmarking**, and it is a technique in its own right, not a local trick. Whenever several
consumers walk the same sorted sequence at different times, each keeps a position that only ever
advances, so resuming is O(1) instead of a rescan from the start. The same shape appears in a k-way
merge (one bookmark per input list) and in a two-pointer sweep over two sorted sets; here there is
one bookmark per supernode, over that supernode's own sorted index set. Sorted order is what makes
the single saved position sufficient: nothing already passed can become relevant again.

**And here is the crux of the relay: where a descendant goes next does not depend on where it is.**

```cpp
    nextUpdateSp[jj] = from + width;                                    // jj's own bookmark
    if (nextUpdateSp[jj] < jjNumNodeIdx)
        descendantUpdateQueue[nf.nodeToSnode(jjNodeIdx[nextUpdateSp[jj]])].push_back(jj);
        //          ^ jj's index set          ^ jj's bookmark
```

Nothing in that expression mentions the supernode currently being factored, or `parent(jj)`, or the
forest at all. It reads `jj`'s own index set at `jj`'s own bookmark and asks which supernode owns
that row. The ancestor being updated is merely the *occasion*, the reason `jj` is awake and
being processed, not an input to the decision.

Two things follow immediately, which are otherwise separate facts to be checked.

The tree is never walked, so there are no levels to skip and skipping is not an optimization. The
relay moves along an index set, and the supernodes it names are whichever ones happen to own those
rows. On the 8x8 grid the 111 hops jump over 92 etree levels, up to 9 at once, and none of that is
deliberate, it simply falls out of where the rows are.

And disjointness is forced rather than merely observed: `jj` has exactly one bookmark, therefore at
most one next destination, therefore at most one queue membership. The single `nextUpdateSp[jj]` is what
makes the whole family of queues encodable as flat arrays.

It is also the sharpest way to state the difference from right-looking. Both read the same index
set. Right-looking reads all of it at once and acts immediately; left-looking reads one run at a
time, and the bookmark is the only thing it carries between visits.

**The relay does not plod up the tree.** Each hop reads the supernode owning `jj`'s next *remaining*
row, so it jumps straight to the next ancestor `jj` genuinely updates and skips any it has no rows
in. On the 8x8 grid Laplacian, 111 hops skip 92 etree levels between them, as many as 9 in one hop.
Each `(descendant, ancestor)` pair is popped exactly once: 164 pairs, 164 pops, no repeats.

Two things about the relay are easy to get wrong.

**It is not avoiding redundant updates.** `ii` still has to update `kk` itself; its contribution is
not subsumed by `jj`'s, and both arrive at `kk`. The point is not that `ii` need only tell `jj`. It
is that `ii` will *reach* `kk` eventually, and can be told where to go next when it gets there. The
absorption property is what guarantees the relay never drops anyone: the ancestor `ii` is re-queued
onto always has more work upward too, so there is always a next hop until the work is done.

**And the eager alternative is not rejected on computation, it costs no more operations at all.**
Enumerating all of `ii`'s ancestors when `ii` is factored is one walk of its index set,
`O(|Idx(ii)|)`, and right-looking does exactly that walk. The enqueue count is identical either way:
lazily, `ii` is queued on `jj`, then `ii` and `jj` are queued on `kk`; eagerly, `ii` is queued on
`jj` and `kk`, then `jj` on `kk`. Three enqueues both times. Measured on the 8x8 grid: 164 pairs,
164 pops, **no pair popped twice**.

What eager enqueueing costs is memory and structure. Every membership would be live at once,
`O(pairs)` against `O(N)` in the number of supernodes, 164 against at most 54 on that grid, and the
gap widens with the matrix. It would also destroy the disjointness of the queues, and disjointness
is what allows them to be flat arrays rather than lists.

So the relay does not save work. **It saves peak storage**, and preserves disjointness as a
consequence.

**The relay does not plod up the tree, either**, which is the crux above, seen from the outside. A
natural misreading is that `ii` is re-queued at every supernode between itself and `kk`. It is not.

Those gaps come from the structure of the index sets, not from branching in the forest. A chain
`0-1-2-3-4-5` with one extra edge `0-5` has a path etree with no branching anywhere, and supernode 0
still updates only 1 and 5, skipping three levels. The tempting implication, `ii` updates `jj` and
`jj` updates `kk`, therefore `ii` updates `kk`, is false, and false in that example: 0 updates 1,
1 updates 2, 0 does not update 2. Absorption runs the other way, and that direction is the one the
relay depends on: if `ii` updates `kk`, every ancestor of `ii` below `kk` updates `kk` too, so there
is always a next hop.

So: **left-looking relays intent, right-looking delivers data.** Right-looking can afford the full
ancestor walk because it acts immediately and stores nothing; left-looking only records that a supernode still
has someone to update, so it needs the relay to keep that record down to one slot per supernode.

The same push-versus-pull direction has a second, structural consequence: **right-looking assembles A
into every front in a prepass, before its main loop; left-looking assembles A lazily, at the top of each
supernode's own iteration.** Because right-looking pushes, a supernode's factorization writes updates into
ancestors the loop has not yet reached, so those ancestors' fronts must already hold their A values when the
push lands. `assembleFromA` assigns rather than accumulates, so it cannot run after a push has landed without
overwriting it; every front is therefore filled up front. Left-looking pulls, so nothing writes into a front
until its own turn, and A can wait until then. This is the one asymmetry between the two dynamic drivers that
is essential rather than incidental: every other difference has been aligned away, but the prepass belongs to
right-looking and the lazy assembly to left-looking, and neither can adopt the other's shape without breaking.

There is a clean way to hold the two loops in the head, and then a small improvement the code takes on top
of it. **The clear description: both drivers expand a front by preserving what is in it.** Per supernode,
each loop widens the front to take its children's delayed columns and keeps whatever the front already holds,
shifting it to make room; that is `expandVal`, and it is the honest mirror-image form. Read this way the two
drivers differ only in the direction of update flow, nothing else: right-looking assembles A up front and
carries it through the expand, left-looking assembles A after the expand and carries an empty front through
it, but both *expand-by-preserving*.

**What the code actually does, with a small improvement: left-looking discards instead of preserving.** When
left-looking expands a front, nothing has been written into it yet, its A is assembled after the expand, its
descendants have not pulled through. So the front to preserve is empty, and preserving an empty front is
wasted work: `expandVal` would allocate the wider block and copy nothing worth copying, then A would be
assembled at an offset past the delayed columns. Left-looking instead calls `resetVal`, which re-sizes the
block to the widened shape and zeroes it in one assign, and then assembles A directly at its final offset. It
is the same result as expand-by-preserving an empty front, reached more cheaply. Right-looking cannot take
this shortcut: by the time it expands, the front already holds A and every update pushed in by an
already-factored descendant, so it must preserve. This is the only place the two verbs diverge, and whether to
give the divergence up for a single uniform loop skeleton is left open in `docs/TODO.md`.

In pseudocode, the symmetric pair first, both expanding by preserving. Right-looking:

```
for J in supernodes:
    assemble A into J
for J in supernodes:
    expand J from its children
    assemble children's delays, contract children
    factor J
    for K that J updates:
        update J -> K
```

Left-looking, the honest mirror, also expanding by preserving:

```
for K in supernodes:
    assemble A into K
    expand K from its children
    assemble children's delays, contract children
    for J that updates K:
        update J -> K
    factor K
```

The two read as reflections: right-looking factors then pushes up, left-looking pulls up then factors, and the
prepass is the price right-looking pays for pushing into fronts it has not reached. What the code actually runs
is this left-looking loop with the one improvement: the front A would be assembled into, then preserved through
the expand, is instead expanded empty and filled afterward, so `expandVal` becomes `resetVal` and the A
assembly moves to after the expand:

```
for K in supernodes:
    reset K at its expanded shape
    assemble A into K
    assemble children's delays, contract children
    for J that updates K:
        update J -> K
    factor K
```

One asymmetry is deliberately *not* visible in any of these blocks, and its invisibility is the point.
At the `update J -> K` step, J has already been contracted in left-looking but not yet in right-looking:
left-looking contracts a child when it processes that child's parent, which for the update target is the
same iteration and comes first; right-looking contracts J only later, when it reaches J's own parent,
which is after J has pushed. So the update runs from a contracted J in one driver and an uncontracted J
in the other. The pseudocode does not distinguish the two, and neither does the code, because the update
does not need to: it reads only J's update area, which is disjoint from the delayed columns that
contraction trims (see *The life of an update*). This is where the two flows come apart. Contraction is
on the delay flow's schedule (child to parent), the update is on its own (descendant to ancestor), and
the update operation is written so it never has to know which side of J's contraction it is on.

One consequence, recorded because it looks like a bug when first noticed. The order in which
descendants arrive on a queue is by *previous ancestor*, not by their own index, so a queue can be out of index
order, commonly, on branching forests. Along a single root-path it is always sorted; every
inversion observed was between different branches. Order does not affect correctness, since the
updates are summed, but it does affect the last bits of that sum.

#### What the two traversals actually cost

**Both do work proportional to the same count: the number of (descendant, ancestor) pairs**, which
is the right unit here rather than supernodes. On an 8x8 grid Laplacian under AMD, 54 supernodes,
that count is 164, supernodes update about 3 ancestors each, up to 5. Neither traversal avoids any of
those pairs; they spend the per-pair cost differently.

| | per pair | shape |
|---|---|---|
| left-looking | one list node allocated and freed, plus a `nextUpdateSp` write | O(1), but an allocation and a pointer chase |
| right-looking | global-to-local map set and cleared | O(\|Idx(kk)\|), but a contiguous sweep over an array already warm |

So it is O(1)-with-an-allocation against O(\|Idx\|)-with-good-locality, and which wins depends on how
large the index sets are. Fat supernodes should favour left-looking, since \|Idx(kk)\| grows while the
list node does not; thin ones may favour right-looking, since a malloc can cost more than sweeping
twenty integers. **That is a prediction and nothing has been measured**; the work item lives in
docs/TODO.md.

**What the relay does not buy is operations.** Enqueueing every ancestor eagerly at `ii`'s factor
time would perform the same number of enqueues and dequeues, one per pair either way. For `ii`
updating `jj` and `kk`, with `jj` updating `kk`: lazily, `ii` is enqueued for `jj`, then `ii` and
`jj` are enqueued for `kk`. Eagerly, `ii` is enqueued for `jj` and `kk`, then `jj` for `kk`. Three
enqueues both times. Nor is enumerating the ancestry expensive, it is one walk of `ii`'s index set,
which right-looking performs without difficulty.

**What it buys is peak storage**, and that is the whole of it: `O(N)` live memberships against
`O(pairs)`, where `N` is the number of supernodes, 54 against 164 on that grid, and the gap widens
with the matrix. The disjointness below follows from it.

#### The queues do not need to be lists

Disjointness, a supernode is on at most one queue at a time, is not just a description of the
schedule. It is the licence for a data structure: one successor slot per supernode covers every
queue simultaneously, so the whole family collapses to flat arrays.

```cpp
    std::vector<std::int32_t> queueHead(snodeSize, NIL);   // first descendant queued for kk
    std::vector<std::int32_t> queueTail(snodeSize, NIL);   // last, so pushes stay FIFO
    std::vector<std::int32_t> queueNext(snodeSize, NIL);   // jj's successor on whatever queue holds it
```

Push becomes three or four integer writes, pop becomes two, and the loop tests
`queueHead[kk] != NIL`. It is the same idea as the `firstChild`/`nextSibling` pair already used for
the forest, applied to a queue instead of a tree. The general pattern: **a family of disjoint lists
over a fixed universe of n elements encodes as `head[numLists]` plus `next[n]`, with no allocation
at all.** If an element could sit on two lists at once the encoding breaks.

Keep the tail array if this is ever done. Dropping it gives LIFO with two arrays instead of three,
but that reverses the order updates are summed into an ancestor, which perturbs the last bits and
would break `test_pipeline`'s assertion that the two traversals agree bit for bit.

**And be precise about what it would buy.** The complexity does not change: push and pop are O(1)
either way. What changes is the constant, and the reason the constant is worth caring about is that
one side of it is a call into the allocator. Related but distinct, the *number of allocations* does
drop asymptotically, from one per pair to three for the whole factorization. Likely larger than
either: the traversal then walks contiguous `int32_t` arrays rather than chasing heap nodes, which
is the part neither an operation count nor an allocation count expresses.

The port uses `std::vector<std::list<std::int32_t>>` today. 0.9 used an `ArraySingleLinkedList`,
whose source did not come across with the rest of the reference, so whether it pooled its nodes is
not known here.

### The parallel with the matrix, and where it stops

The same shape appears one level down, in `SparseMatrixStatic` and `SparseMatrixDynamic` in
`experiments/storage-options`, which is where the pattern was measured before it was adopted here.
The correspondence is close enough that the two kernels read alike:

```cpp
const std::int32_t* nodeIdx = nf.nodeIdx(jj);   // a supernode's value block, in the factor
const Val*          val     = nf.val(jj);

const std::int32_t* rowIdx  = A.rowIdx(j);      // a column, in the matrix
const double*       val     = A.val(j);
```

Two places it stops. The matrix's `colSize(j)` has no single counterpart, because a column is a run
and a supernode's value block is a rectangle: its extent takes two numbers, `frontSize(jj)` for the
columns and `frontSize(jj) + updateSize(jj)` for the rows, which is also the leading dimension.
That pair is what lets the numeric kernels take `(Val* block, rows, cols, ld)` and never see either
factor class. And `factorization()` has no counterpart at all, because a matrix carries no method
identity while a factor must say whether the solve conjugates.

### Reading the names

An offset array keeps the `Ptr` suffix (`colPtr`, `snodeNodeIdxPtr`, `snodeValPtr`); a lookup does
not (`rowIdx(j)`, `val(j)`, `nodeIdx(jj)`, `val(jj)`). Where a class holds both the whole array and
the per-entity slice, they are overloads of one name: `nodeIdx()` returns the flat buffer and
`nodeIdx(jj)` the address of one supernode's slice within it. See the 2026-07-19 and 2026-07-17
entries in DESIGN_DECISIONS for how those names were settled.
