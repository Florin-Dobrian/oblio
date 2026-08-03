# Ordering Experiment

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
| `md4` | maintained degrees, refreshed only where they changed | 5.7, 5.8 |
| `md5` | degree buckets, so the minimum is walked to, not scanned | 5.9, 5.10 |
| `mmd1` | multiple elimination: a batch of pivots per refresh | 5.11, 5.12 |
| `mmd2` | the rest of genmmd, one pass at a time | 5.11, 5.12 |
| `amd1` | approximate degree: a bound instead of a set union | 5.13, 5.14 |
| `amd2` | the rest of amd_1 and amd_2, one pass at a time | 5.13, 5.14 |

## What the layers show

**`md1` through `md4` return the same ordering.** That is the point of the first four sections:
the heuristic was fixed in `md1` and everything after is implementation, so those layers can be
verified by demanding an identical permutation. From `md5` on the permutation moves, and the
reason is worth separating from the reason the fill moves.

Four things can happen when a layer is added, and all four occur here:

```
                              order        fill        what the change is
md1 -> md2                    same         same        a change of representation
md2 -> md3  (mass elim.)      DIFFERENT    same*       a reordering, free step by step
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

`graph1` is a 4-cycle, the smallest graph that fills at all: eliminating any vertex forces its
two neighbors together, for one fill edge. It is also where md3 merges everything that is left
in a single step, 1 taking 2 and 3, which makes it the simplest case for reading the cost of
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
supervariable**: at the step whose pivot is 0, vertex 4 has nothing explicit left but belongs to
`c1` as well as to the new clique, so `I[4] == {pivot}` fails even though everything 4 reaches lies
inside that clique. It orders as 2 1 3 0 4 with no merge and no fill, where the exact test would
give 2 1 3 (0 4). It also separates amd2's two extra mechanisms: with aggressive absorption on,
`amd2` takes four steps and reports `merged = 4, absorbed = c1`; with it off, five steps and no
merge, exactly like `md3`. The hashing plays no part either way. See the section on mass
elimination.

`graph6` has six vertices and eight edges and is also in `md1`, `md2` and `md3` only. One small
graph carries three things at once. **Its supervariable {0, 4} is a supernode but not a
fundamental one**, since the forest is 2 -> 1 -> 4 and 3 -> 0 -> 4 and 4 already has 1 as a child
when 0 merges into it. It orders as 1 5 (0 4) 2 3, where the exact test would give
1 5 (0 2 3 4). The merge lands at step 2 of 5, so the run continues afterwards and the
selection degree, 3 over {2, 3, 4}, differs from the external degree, 2 over {2, 3}, by exactly
the size of what merged. And `super_members` ends with a hole in the middle, slot 4 empty between
two used ones, while no pivot equals its own step number.

`graph7` has five vertices and six edges and is the **pairwise case**, also in `md1`, `md2` and
`md3` only. At the step whose pivot is 0 and whose clique is {2, 4}, vertices 2 and 4 are
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
example when given none.

Once a layer has been pulled into the main tree it gains a second check, against the production
driver extracted from it:

```
make production       builds production_cpp, linking ../../src
./production_cpp mmd1        MMD1 on all seven graphs
./production_cpp amd2        AMD2, reporting nnz(L) rather than the order
./production_cpp mmd1 3      just the third
```

`make test` runs it for every layer named in the Makefile's `PORTED` list, `mmd1` and `amd1`, and
requires the order lines to agree with the prototype's. `PORTED_FILL` names the layers checked by
nnz(L) instead, `amd2` today: production skips the postorder that prototype ends with, since
`ElmForestEngine` does that work itself, so the two permutations legitimately differ while the fill
does not. That check runs Oblio's own symbolic factorization, which is why the target links the
whole library. The prototype's `matrix1` is left out of it, being the example that exists to
exercise `amd2Aat` and `amd2Preprocess`, which production has no counterpart to. It links `../../src/QuotientGraph.cpp` and
`../../src/Mmd1.cpp` directly rather than copying them, which is the opposite of what the vendored
target does and deliberately so: a copy is right for code that is not ours to edit and wrong for
code being actively changed at both ends, since noticing when the two come apart is the whole point.
The harness feeds each graph as a full-symmetric CSC with the diagonal present, which is what a
`SparseMatrix` holds, so the production path under check includes `buildGraph` dropping the
diagonal rather than only the driver.

The vendored routines have their own target, since they are not layers and have no Python twin:

```
make vendored         builds vendored_cpp
./vendored_cpp        both routines on all seven graphs
./vendored_cpp 3      just the third
```

That target compiles with warnings off, because the two files are not ours to clean up.

`vendored/vendored_mmd.cpp` and `vendored/vendored_amd.cpp` are copies of `src/Mmd.cpp` and
`src/Amd.cpp`, kept here so the comparison is self-contained, and never edited. `vendored.cpp`
only feeds them the same seven graphs and prints their permutations in our format.

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
not a sequence of ideas but the completion of two, so from here the steps are bigger. Two
versions each.

**mmd1, the idea. Done.** Multiple elimination. A batch is an independent set in the current
elimination graph, enforced by evicting every reached vertex from its bucket with a stale degree,
so no later pivot in the batch can be a neighbor of an earlier one. `delta` widens the batch to
near-minima. Everything else is md5 unchanged: the quotient graph, mass elimination, the buckets,
the expansion.

`delta` lives here, with the full signed range: negative takes one pivot per round, which turns
the batching off and reproduces md5's ordering exactly, verified on 100 random graphs; 0 through
n - 1 widen the window; anything larger saturates. No weight array, because mass elimination
merges only into the pivot and the pivot dies in the same call, so no live vertex ever stands for
more than one original vertex, checked over 200 graphs and 1386 eliminations.

**mmd2, genmmd complete.** Six additions, all of them holes rather than ideas, listed against
the vendored routine that carries each:

1. **The prepass** (`genmmd`, the loop over `head[1]` before the main loop). It numbers every
   vertex in the degree-1 list, marks each `marker[mn] = maxint`, and never refreshes a neighbor.
   `mmdint` maps degree 0 to 1, so isolated and degree-1 vertices are numbered together. Two
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

Two pieces of the vendored codes are consequences of packing state into reusable integer arrays
rather than features of the ordering, and neither is modeled:

- MMD's `maxint` overflow reset, `if (tag >= maxint) { tag = 1; for each i with marker[i] < maxint,
  marker[i] = 0; }`. We use the same mark-and-tag idiom, so the reset would be needed only if the
  tag could wrap, which it cannot at the sizes these prototypes run.
- AMD's workspace with `iwlen`, `pfree` and the `ncmpa` garbage collection. That is a flat pool
  being compacted when it fills, and our member lists grow on their own.

Both are what the ordering engine will need once the lists live in one pool; the mark arrays
themselves are already here, in `mark` and `touched_round`.

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
refresh order reversed        12570     73184    213784    504177
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

The test is `vendored.cpp`, which links `vendored/vendored_mmd.cpp` and
`vendored/vendored_amd.cpp`, copies of `src/Mmd.cpp` and `src/Amd.cpp` that are never edited, and
runs both routines on the same seven graphs, printing permutations in our format. mmd2 and amd2
are accepted when every feature above is present and exercised, nnz(L) matches the vendored
routine on the seven graphs and on random ones, and every remaining order difference is traceable
to a tie.

## Grid mode, and what the missing features are actually worth

Every prototype takes an example number; `mmd1`, `mmd2`, `amd1` and `amd2` also take a grid:

```
./amd1_cpp grid 22        one 22x22 grid Laplacian, counters only
```

The trace is discarded as it is written rather than afterwards, a filtering streambuf keeping only
the closing counter lines, because at n = 10000 a full trace runs to gigabytes and a process
holding it is killed. The grid is not an eighth example: nothing about it illustrates a mechanism
and its trace is unreadable. It exists so the counters can be read at a size the seven cannot
reach.

**Why it was added.** `benchmarks/ordering` shows our production MMD1 and AMD1 running about 4.5x
and 3.4x slower than the vendored routines. Two explanations were available and they call for
opposite work: our per-list allocation against their flat array, or the mechanisms these layers do
not yet have. Counters separate them, in units no allocator can move.

**Measured on a 100x100 grid, n = 10000.** The number that matters per branch is the dominant
inner quantity, not the refresh count.

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

**For MMD the gap is mostly ours.** `mmd2` refreshes 1.6 times less often, and its `q2h` shortcut
also makes many of the remaining refreshes cheaper, which nothing here counts, so 1.6x is a lower
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

## Two bugs this found, both ours

Worth recording, because both were invisible to the checks in place at the time.

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
from scratch at every step, so most of what they refresh is outside the bound's domain and would
have to be counted exactly anyway. md4 is the layer that narrows the refresh set to exactly C[p],
which is exactly the bound's domain. The bound is inherently incremental, md4 is where the
algorithm becomes incremental, and those two facts meeting is why amd belongs at the top of the
chain rather than beside md2.

So the grid is smaller than it first looks: two exact rungs where the estimate is expressible but
not usable, then a genuine choice from md4 up. Each layer from md2 carries a short subsection
saying where it stands, and amd1 builds the approximate cell once, at the top, as running code.

**And the axis stops being orthogonal at the top, which is why the fork is real.** mmd and amd are
not two knobs that happen to be turned separately. They attack the same cost from opposite ends,
mmd making the refresh rare and amd making it cheap, so their gains overlap rather than add. They
also interfere: the anchored bound is anchored at ONE new clique, and a batch produces several, so
the anchoring would have to be redone against their union; and the bound stays tight only when
degrees are refreshed often, which is exactly what batching gives up. That is measured here rather
than assumed, and the measurement is in the delta section. So md1 through md5 is one chain, and
the fork at the top into mmd, which stays exact, and amd, which goes approximate and does not
batch, is a real division rather than a filing convention.

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
at every step. Making the degree cheap again is what md4 and md5 are for.

### md2 with the exact degree, which is what it does

Every layer from md2 on faces the same choice, since from here a degree is a union rather than a
lookup, and a union can be counted or estimated. Written out in set operations, one step of md2
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
formed, every u in C[c] has A[u] = A[u] - C[c], and neither set grows afterwards. So

```
A[u] & C[c] = {}     for every c in I[u]
```

The explicit part never overlaps a clique. All overlap is clique against clique.

**The decomposition, anchored at the new clique**, where p is the pivot just eliminated:

```
reach(u)  = ( A[u] - C[p] ) | ( C[p] - {u} ) | ( C[c] - C[p]  for c in I[u] - {p} )
bound(u)  = |A[u] - C[p]| + |C[p] - {u}| + sum |C[c] - C[p]|
```

The first line is an identity, exact, for any u that C[p] reached. The second line is the first
with the union replaced by a sum, which is one substitution and the only one, so

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
exactly once: it is counted at the step it is merged, or at the step it is chosen as pivot, and if
it later sits inside some other supervariable's member list it is not counted again. So
`num_eliminated` really does count originals, and the loop condition `num_eliminated < n` compares
like with like.

**Why it is cheap** is one property and worth stating separately from why it is correct:
|C[c] - C[p]| depends on the clique c alone and not on the vertex asking, so it is computed once
per clique per step and then read once per vertex. The exact degree recomputes a union per vertex.

**And here is what md2 cannot do with any of it.** Look again at the first line: it holds for u in
C[p], because it puts C[p] - {u} inside reach(u). For a vertex the step never reached, the
decomposition is not merely loose, it is wrong. So bound(u) covers exactly the vertices the new
clique reached, and nothing else.

md2 recomputes every live vertex at every step. Most of what it refreshes lies outside C[p] and
would have to be counted exactly regardless, so the estimate would apply to a shrinking minority of
the work and save nothing worth having. The layer that narrows the refresh set to exactly C[p] is
md4, and C[p] is exactly the bound's domain. bound(u) is inherently incremental; md4 is where the
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

**One detail that looks like an opportunity and is not.** At md4 and md5 the refresh already
computes reach(u) for exactly the members of C[p], so the query seems to be sitting there already.
It is not: the refresh runs after the clique is installed, so it yields the new reach, and fill
needs the old one. Getting it would mean a second query per member, before the update, which doubles
the layer's remaining cost.

**What survives is the aggregate, and it survives intact.** At the moment a vertex is chosen, C[p]
is the off-diagonal pattern of its column of L, so |C[p]| accumulates nnz(L) at every layer for
free, exactly as md1's degree does. Total fill is nnz(L) - nnz(tril A) and is reported by every
layer in the ladder. What md2 gives up is per-step, per-vertex fill. The total never depended on it.

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
declines a genuine pair. graph5 shows exactly that at its fourth step.

md3 detects the case positionally rather than by comparison. Immediately after a clique is
formed it checks each member for `not A[u] and I[u] == {pivot}`, which says the new clique is
everything u can still reach, hence that u and the pivot were indistinguishable before the
step. The test is cheap, two container reads per neighbor, and CONSERVATIVE, meaning
sufficient rather than necessary. It misses supervariables, in two independent ways, and both
are a deliberate trade of coverage for cost.

**Miss one, inside the pivot's own neighborhood.** The test fails when A[u] is non-empty but
contained in the clique, or when I[u] holds another clique whose members happen to lie inside
it. Such a u is genuinely indistinguishable from the pivot and is skipped anyway. The exact
test is `md3_neighbors(A, I, C, u) <= C[pivot]`, and it costs a reachability query per
candidate, O(d) set unions per step where the cheap test is O(d) constant checks. Measured on
600 random graphs, the exact test merges strictly more on 415 of them and never changes the
fill.

Five vertices suffice to see it. With adjacency 0:{3,4}, 1:{2,4}, 2:{1}, 3:{0}, 4:{0,1}, md3
orders 2 1 3 0 4, five steps and no merge at all. Step 3 is the one to look at: its pivot is
0 and its clique is {4}. At that moment A[4] is empty and I[4] is {c0, c1}, the second clique
left over from eliminating 1 at step 1, with members {4}. So everything 4 can reach lies
inside the new clique and 4 is indistinguishable from the pivot, but `I[4] == {pivot}` is
false and the cheap test declines. The exact test merges it, making the order 2 1 3 (0 4).

**Miss two, between two members of the clique.** Even the exact test compares each candidate
only AGAINST THE PIVOT. Two members can be indistinguishable from each other while neither is
absorbable into the pivot, and no pivot-relative test finds them. An exhaustive pass over the
clique would: all pairs, O(d * d) comparisons at O(d) each, so O(d * d * d) per step. That is
what hashing reduces, bucketing by structure and comparing only within a bucket, which is
what amd2 does and what the section on detecting supervariables against each other covers.

Five vertices again, and this one is graph7. With adjacency 0:{1,2,4}, 1:{0,4}, 2:{0,3,4},
3:{2,4}, 4:{0,1,2,3}, md3 orders 1 0 (2 3 4), three steps, the last of which merges 3 and 4
into 2. Step 1 is the one
to look at: its pivot is 0 and its clique is {2, 4}. After pruning, A[2] is {3} and A[4] is
{3}, and both have I equal to {c0}. Neither is absorbable into the pivot, since each still
reaches 3, which is outside the clique. But their closed neighborhoods are equal, 2 reaching
{3, 4} and 4 reaching {2, 3}, so 2 and 4 are indistinguishable from each other already at
this step and could be merged here. No test framed against the pivot will ever see it, and
the exact test does not help either: both orders are 1 0 (2 3 4). md3 does group them one
step later, when 2 becomes the pivot and both 3 and 4 fall to its own test, which is luck
rather than detection: a graph where 2 and 4 never share a later step would keep them apart
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
steps and reports `merged = 4, absorbed = c1`; with it off it takes five and merges nothing,
behaving like md3. The hash detection reports nothing in either run, so this case is entirely
about absorption and not about comparison. So there are three ways to recover a missed
supervariable, and they are genuinely different: sharpen the test, compare vertices pairwise,
or clean up the structure until the cheap test is no longer fooled.

Three consequences follow, each covered in its own section below. The degree becomes
WEIGHTED, since a neighbor now stands for several original vertices, and it must be EXTERNAL,
excluding the supervariable's own members. The nnz(L) count stops being a running degree sum,
because a step is now w consecutive columns rather than one. And the order changes, since a
merged vertex is eliminated immediately where md2 would have reached it later.

The structure this produces is related to a supernode but neither contains nor is contained
in one, and the relation is worth stating carefully because it is easy to get backwards.

A supervariable always has exactly nested patterns, since its members were indistinguishable
when they merged, so it is always a supernode in the general sense of consecutive columns
with the same structure outside the block. What it need NOT be is a FUNDAMENTAL supernode,
which additionally requires each column after the first to have exactly one child in the
elimination tree. Indistinguishability is a property of the current elimination graph and
says nothing about children acquired earlier, on branches the step never touched. graph5 is
the counterexample: 0 and 4 are indistinguishable at the fourth step, yet the tree is
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
which md3 orders as 1 5 (0 4) 2 3, five steps with 4 merging into 0 at step 2, and no fill at
all. The exact test takes it further, to 1 5 (0 2 3 4), a single supervariable over all four.
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
clique is its only remaining route out. Together they say u's reachable set is contained in the
clique, which is the fill-free condition.

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

graph1, the 4-cycle, where md3 merges at the last step:

```
        scan   steps   candidate touches
md2      21      4            10
md3      16      2             7
```

The merge takes 2 and 3 into the pivot at step 1, which ends the run, so md3 never runs md2's
steps 2 and 3. Nothing downstream is visible here, which is the limitation of this example.

graph6, where the merge lands at step 2 of 5 and the run continues:

```
        scan   steps   candidate touches
md2      62      6            21
md3      53      5            18
```

Step by step, the two runs are identical through step 2: same live sets, same pivots, same
scan costs of 16, 15 and 14. Vertex 4 merges into 0 at step 2, and from there md3 scans
[2, 3] then [3] where md2 scans [2, 3, 4], [3, 4] and [4]. One fewer candidate in every
remaining scan, and one fewer step at the end, which is 62 down to 53.

Two effects are folded into that number and worth separating. Fewer candidates: a vertex
merged at step k is absent from every scan after k, so it stops being touched once per
remaining step rather than once. And cheaper candidates: the merged vertex is stripped from
every clique, so the queries that remain walk shorter member sets. On graph6 the surviving
queries at step 3 cost 3 each against md2's 4, purely because C[c0] lost a member.

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
selection degree is evaluated over the neighbors as they stand before the step.
external_degree in the nnz(L) block is recomputed from C[pivot] after md3_eliminate returns,
because vertices that were neighbors a moment ago are members now. The difference between the
two is exactly the total weight that merged during the step, and confusing them is the
double-counting trap described in the previous section.

The star makes it concrete. At the step that eliminates hub 0 with leaf 4 still live, the
selection degree is 1, since 4 is a live neighbor. After the merge the supervariable has size
2 and external_degree is 0, since C[0] no longer contains 4. The supervariable's two columns
then hold one entry and zero entries below their diagonals, and the closed form gives
2 * 0 + 2 * 1 / 2 + 2 = 3, which matches.

The term comes from the AMD literature, where the contrast is with the TRUE degree, which
does count the internal members. Minimum degree with supervariables uses the external one,
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
eliminated was discarded from the set at that earlier step, and its entry belongs to that
earlier column. Fill created at this step lands among the neighbors, in later columns, and
cannot alter a set that was read before the elimination ran.

md1 and md2 use it directly, one step being one column. The pivot degree is the column count,
so a running degree_sum plus n diagonals is nnz(L), with no second pass and nothing else
stored. The two layers compute the neighbor set differently, md1 reading A[pivot] and md2
unioning the explicit adjacency with the cliques, but the sets are equal at every step, so
the counts are too.

md3 cannot do that, because one step is now w columns rather than one. A supervariable of size
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
ordinary degree, giving degree + 1 per step, which summed over all steps is degree_sum + n.
So md1 and md2 are the special case of md3's accounting, not a different scheme.

The degree computed at selection is the wrong number to accumulate. At that moment the
vertices about to merge are still neighbors and still counted, so degree = ext + w - 1, the
pivot itself being the one member it already had. Accumulating from it double counts them,
once as neighbors and once as members of the supervariable. This is why md3 recomputes
external_degree from C[pivot] after md3_eliminate returns, rather than reusing the value it
printed in the step title. Getting this wrong reports 43 against the true 37 on graph3.

The accumulation has to be per step. Both ext and w vary from step to step, and neither is
recoverable at the end from a single running total, so md3 evaluates a closed form each step
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
|C[p] - {u}| inflates every bound the step produces. `Amd.cpp` does the same thing at
`degme -= nvi`, and our own version of this was a real bug, recorded in the bug section above.

**Mass elimination itself is indifferent.** A[u] == {} and I[u] == {p} is a structural test that
never consults a degree, so it fires the same whichever of the two numbers the layer is computing.

md3 also still repicks from scratch, so it is where md2 is: the estimate is definable and not worth
using. The layer that changes that is md4.

## md4: maintained degrees

Every layer so far rebuilds a neighbor set for every live vertex at every step, keeps the
smallest and throws the rest away. md4 keeps the degrees in an array and refreshes only what
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

**The refresh.** After the elimination, and only over the clique's survivors: the merged
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

Nothing is sorted. C[pivot] is a list in both twins and the refresh walks it as it stands, which
is the same sequence on each side. The tag comes back from the query because a degree computation
advances it and Python has no reference parameters.

Two query sites remain in the whole file, and that is the layer's claim: md4_eliminate's
first line, which becomes the clique, and the refresh above. Everything else reads integers.
The run prints the count, and on graph3 it is 34 against the 85 neighbor queries md3 makes,
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
every live u at every step. md4 evaluates it only for u in C[p], because nothing else can have
changed. Written as the layer sees it:

```
for u in C[p]:  degrees[u] = degree(u)        md4, exact
for u in C[p]:  degrees[u] = bound(u)         md4, approximate
```

and those two loops range over the same set. That is the coincidence the whole placement rests on.
The second line of the decomposition is a statement about members of C[p] and about no others, so
at md2 and md3 most of what gets refreshed is outside its domain and has to be counted exactly
anyway. At md4 the refresh set IS the domain. bound(u) is inherently incremental, md4 is where the
algorithm becomes incremental, and there is nothing left over.

**A second thing arrives with it.** Of the two caps in the md2 subsection,

```
degree_old[u] + |C[p] - {u}|
```

needs degree_old[u], and a cache is what makes a previous degree exist. It is often the tightest of
the three, since a vertex whose reach barely grew is bounded by what it was plus the new clique. It
also stays valid when the cached value is itself a bound rather than a degree, which is what makes
the chain work inductively: degree_old[u] >= degree(u) at every earlier step means it is still an
upper bound now. Nothing has to be exact anywhere.

**And the costs finally differ.** Both loops walk C[p], but degree(u) unites the members of every
clique in I[u], per vertex, while bound(u) reads one number per clique, and that number,
|C[c] - C[p]|, is computed once for the whole step. So md4 is the first layer where the two columns
are doing measurably different amounts of work rather than the same work on different populations.
The decomposition is md2's, the counting is md3's, and what md4 contributes is a place to stand.

## md5: degree buckets

md4 left one O(n) per step in place: the scan over every live vertex to find the smallest
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
the same order at every step and pick the same pivot; it pays O(bucket) for insert and remove,
which is the one place it is asymptotically behind its twin.

n slots is exactly right, indices 0 through n - 1. A live vertex counts only live neighbors,
so its degree is at most n - 1, and the walk stops at the first non-empty bucket, which
exists while anything is live, so no index above n - 1 is ever filed or probed. SuiteSparse
AMD sizes its Head array the same way; MMD's head array looks one longer only because it is
indexed from 1, and it files an isolated vertex under degree 1 rather than 0, so it has no
degree-0 bucket at all. Every vertex is filed here, isolated ones included. min_degree starts
at the true minimum, the tightest legal value for a lower bound; starting at 0 would also be
correct and would cost one extra walk on the first step.

**Fragment 2, the pop.**

```python
        while not buckets[min_degree]:         # walk up to the first live bucket
            min_degree += 1
            num_bucket_probes += 1
        num_bucket_probes += 1
        pivot = buckets[min_degree][0]         # the head, whatever was filed last
```

The walk only ever climbs, and min_degree is never reset between steps, so the work is
amortized across the run rather than paid per step. Termination rests on the outer loop
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
dead vertex left in a bucket would be popped as a pivot on a later step. The order within
each pair matters: the bucket index is read from degrees[u], so the removal has to come
before the zeroing, or the vertex is erased from buckets[0] and left where it was.

**Fragment 4, the refresh.**

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

Same refresh set as md4, and the only change is that the new degree goes through the helper so
the bucket moves with it. The helper exists so the three steps cannot be written half-way. A
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

Only the refreshed vertices can have moved down, so that is the complete candidate set. The
bound is lowered only here and raised only by the walk in fragment 2, which is what keeps the
total climbing work bounded over the run. Lowering is the safe direction: a bound that is too
low costs a few extra probes, while a bound that is too high skips a non-empty bucket and
picks the wrong pivot. Degrees that rise need no attention at all, since such a vertex is
filed higher and the walk will reach it.

Three invariants hold after every step, and they are what the buckets row and the min degree
line in the trace exist to show: every live vertex is in the bucket matching its degree, no
cached degree is stale, and no live vertex has a degree below min_degree. Checked across 400
random graphs, no violations.

The run prints both metrics. Degree computations are unchanged from md4, as they must be,
since this layer touches only how the minimum is found. Bucket probes replace what was an
n-per-step scan: 12 on graph3, 5 on graph4, 7 on graph6.

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

## mmd1: multiple elimination

Refreshing degrees is the expensive step, so do it less often: eliminate a whole INDEPENDENT
SET of least-degree vertices before refreshing anything. Non-adjacent pivots cannot disturb
each other's degrees, so every pivot in a batch is still a true minimum-degree vertex when it
is taken. That is Liu's M in MMD, and it is the first layer whose ordering differs from md1's
for a reason other than a tie in the same graph state.

**This is the fork, and it is a real one.** md1 through md5 is a single chain because bound(u) is
one definition that no layer alters, and because the layer that makes it usable, md4, does so
without any change to itself, as the notes at the end of md2, md3, md4 and md5 set out. At the top
that independence ends. mmd and amd attack the same cost, the refresh, from opposite ends: mmd
makes it rare, amd makes it cheap, so their gains overlap rather than add. They also interfere. The
anchored bound is anchored at ONE new clique, and a batch produces several, so an amd-flavored mmd
would have to re-anchor against their union and would lose the per-clique reuse that makes the bound
worth having. And the bound stays tight only when degrees are refreshed often, which is precisely
what batching gives up: the delta measurement below shows the bound's own failure mode in the exact
setting, where a wide batch leaves a whole evicted set invisible for a round. So mmd1 and mmd2 stay
exact and take the batch; amd1 and amd2 stay one pivot per step and take the bound. Combining them
is a question rather than a step, and it is not one this experiment answers.

Six of the seven functions are md5's with the prefix changed: mmd1_neighbors, mmd1_storage,
mmd1_eliminate, mmd1_refile and the two display functions. That is the pattern across the
whole ladder from md2 onward. The elimination itself has not changed since the quotient graph
appeared, and the degree cache and buckets have not changed since md4 and md5. What each layer
varies is the SELECTION POLICY: recompute per candidate, cache, bucket, and now batch. So for
mmd1 the whole layer is the driver.

**The independent set is never searched for.** It falls out of the eviction:

```python
            for u in C[pivot]:                 # EVICT, with a stale degree
                buckets[degrees[u]].discard(u)
                touched.add(u)
```

Every vertex the pivot reached leaves the buckets and stays out until the round ends. So
whatever is still filed was not reached by any pivot taken so far, hence is not adjacent to
any of them, and draining a bucket drains an independent set. The eviction is also what makes
the deferred refresh safe: a vertex with a stale degree is not a candidate, because it is not
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
            if delta < 0:                      # one pivot per round, as md5 does
                break
```

The shape matters. On entry the outer walk has left buckets[min_degree] non-empty, so the
first iteration always takes a pivot, whatever the limit says. Only after that does delta
decide whether the round continues. Written the other way round, as a `while min_degree <=
batch_limit` guard, a negative delta makes the loop body unreachable, the batch comes out
empty, nothing is eliminated and the driver spins forever. That was a real bug in the first
draft of this rewrite.

**delta is the whole control, and its sign selects between two behaviors.** delta = 0 keeps
the batch to true minima. delta > 0 admits vertices up to delta above the minimum, which are
not minimal, so that is a concession in quality for still fewer refreshes. delta < 0 takes one
pivot per round, which is md5 reached through this code path, and it is the check on the
batching rather than a feature of it: mmd1 at delta = -1 reproduces md5's ordering exactly,
verified on the seven graphs and 150 random ones.

**The refresh is one pass per round, over everything the batch touched.**

```python
        refreshed_vertices = [u for u in touched if not eliminated[u]]
        for u in refreshed_vertices:
            neighbors_u, tag = mmd1_neighbors(A, I, C, mark, tag, u)
            degrees[u] = len(neighbors_u)
            mmd1_file(buckets, filed, degrees[u], u)
        num_degree_computations += len(refreshed_vertices)
        min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])
```

Note the asymmetry with md5's refresh, which called mmd1_refile: here the vertex is already
out of its bucket, evicted during the batch, so the refresh only writes the degree and files
it. The `not eliminated[u]` filter matters because a vertex evicted early in the round can be
merged away by a later pivot in the same round.

**Choosing delta, measured.** The vendored driver passes 0, which is also what SPARSPAK does and
what Liu's paper treats as the default. On grids, with fill and refresh count both reported:

```
grid 22x22, n=484
   delta   -1: nnz(L) 4773   degree computations 2690   rounds 367
   delta    0: nnz(L) 4684   degree computations 1859   rounds  36
   delta    1: nnz(L) 4754   degree computations 1733   rounds  22
   delta    2: nnz(L) 4706   degree computations 1756   rounds  21
   delta    4: nnz(L) 4747   degree computations 1601   rounds  14
   delta    n: nnz(L) 5964   degree computations 1514   rounds   9
```

Two things in that table are worth more than the recommendation they support.

The first is that delta = 0 beats delta = -1 on BOTH axes, here and on the 10 by 10 and 16 by 16
grids as well. Batching is not trading quality for speed at that setting; it is simply better.
The wager appears only as delta grows, and by delta = n it is decisively bad: 27 per cent more
fill for 19 per cent fewer refreshes.

The second is why it goes bad, which is not stale degrees. Every pivot in a batch has a CORRECT
degree, since anything whose degree could have changed was evicted. What the round cannot see is
the evicted set, and those are exactly the vertices whose degrees typically FELL, so they are the
candidates that should be picked next. With delta = 0 the round ends as soon as the minimum
bucket drains and they come back at once. With delta = n the round keeps climbing through the
degrees, taking vertices of degree 4, 5, 6 while better candidates wait until the end.

**delta is total, and the top end is clamped.** Any negative value means one pivot per round; 0
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
invisible for the rest of the round, so the choice is made among the untouched remainder
rather than among all candidates. The batch does not pick a worse vertex, it picks a different
vertex of the same degree. Minimum degree is famously sensitive to tie-breaks, so the fill
moves by a fraction of a percent, in either direction. Batching across connected components is
free, since the components cannot interact at all; batching within one is the wager.

**The metric is rounds against pivots.** The closing line prints degree computations, bucket
probes and rounds, and the ratio of pivots to rounds is the average batch size, which is what
the batching buys. md5 and mmd1 pay the same per refresh; mmd1 pays for fewer of them. On
graph3, md5 makes 34 degree computations over 10 steps and mmd1 makes 26 over 5 rounds; on
graph4, 20 over 5 against 15 over 3. graph1 is the case where batching buys nothing, 6 either
way, since md5 already finishes it in two steps.

## mmd2: the extras

mmd1 is multiple elimination and nothing else. genmmd is multiple elimination plus six other
things, none of them an idea and all of them load-bearing for matching the vendored routine.
mmd2 adds them one pass at a time, and this section grows a subsection per pass. The checklist
and its vendored references are in the plan section above.

### Pass 1: the prepass

genmmd numbers every vertex in the degree-1 list before the main loop runs, and never refreshes
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

`mmdupd` does not walk a flat list of reached vertices. It walks the ELEMENTS this round created,
`el = list[el]`, and for each one computes `dg0` once, the size of that element, before visiting
its members. A member is classified by what it has left BESIDES the new element: `mmdelm` stashes
`fwd[rn] = nq + 1` where `nq` counts the survivors of the compaction, which in our split
representation is `len(A[u]) + len(I[u]) - 1`. `nq == 1` puts the vertex on the `q2h` list,
anything else on `qxh`.

**Why the split pays.** Everything a q2h vertex reaches is either inside the element, already
counted in `dg0`, or comes from its one other source. So its degree is `dg0` plus what that
source contributes, and the union is never built. The qxh case pays for the full union as
before. Measured on grids, the q2h path takes 36 per cent of the refreshes at 10 by 10, 42 at
16 by 16 and 44 at 22 by 22, so it is not a rare case.

**Two mark levels, and the bug that taught me why.** My first version marked the other source's
members with the element's tag, which made a second q2h vertex in the same element skip them as
already counted and report a degree too small. `mmdupd` avoids this with two levels: element
members carry `mt`, above every tag used in the round, while each vertex gets a fresh `(*tag)++`
for its own walk. Ours does the same with `element_tag` and `vertex_tag`, testing both.

**The check on this pass, and the check that was wrong.** Every degree the shortcut produces must
equal the full union, and it does, on 307 graphs. My first attempt at that check called
`md2_neighbors` to get the true value, which advances the tag and overwrites the mark array, so
it destroyed the state it was inspecting and reported failures that were its own doing. The
verification has to compute the true neighbors independently, with plain sets and no marks. Worth
remembering for the passes to come: an instrument that shares state with the thing it measures
is not an instrument.

What moves is the filing order, since the refresh is now element by element with q2h first, and
filing order decides what a bucket holds. So the permutation changes even though every degree is
identical: on the grids, nnz(L) goes 636 to 633, 2088 to 2101, and 4684 to 4684 against mmd1.
A vertex reached by two pivots in the same round is still refreshed once, skipped on the second
visit by the `filed` flag, which is `if (bwd[en] != 0) goto n2200` there.

### Pass 3: the pairwise merge, and outmatched marking

Both live in the same branch of the q2h walk, reached when a member of the one other source is
ALSO a member of the new element:

```c
    else if(bwd[nd]==0){
        if(fwd[nd]==2){qsize[en]+=qsize[nd];qsize[nd]=0;marker[nd]=maxint;
                       fwd[nd]=-en;bwd[nd]=-maxint;}
        else if(bwd[nd]==0)bwd[nd]=-maxint;}
```

**The merge.** If nd is q2h too, its only other source is that same element, so en and nd reach
exactly the same vertices: indistinguishable, and en absorbs nd. This is the first merge in the
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

**Outmatched.** If nd is not q2h it has sources besides these two, so its reach contains en's and
it can never be the minimum before en. MMD withdraws it from the degree lists rather than
refiling it, `bwd[nd] = -maxint`. It is not merged and not eliminated, just held out until
something reaches it again, at which point `mmdelm` restores it with `bwd[rn] = 0`. Ours is an
`outmatched` flag, cleared in mmd2_eliminate for every vertex the new clique reaches. Both are
withheld from the q2h and qxh lists while the flag is set.

**Weights everywhere.** dg0 becomes a weighted sum, the q2h walk starts at
`dg0 - len(super_members[u])` rather than `dg0 - 1`, every count adds `len(super_members[v])`,
the qxh path goes through mmd2_degree, and the nnz(L) accounting sums the same over the live
members of C[pivot].

**What it does.** On the same 307 graphs: 245 pair merges and 126 outmatched markings, all
permutations valid, all reported nnz(L) correct against an independent symbolic factorization,
and the fill against mmd1 better on 2, same on 289, worse on 16. On grids the mechanisms fire
steadily: 18 merges and 55 outmatched at 10 by 10, 56 and 135 at 16 by 16, 111 and 276 at 22 by
22, with nnz(L) 631, 2105 and 4783 against mmd1's 636, 2088 and 4684.

That last column is worth reading honestly. These two mechanisms do not improve the ordering
here; they make it match the vendored routine. mmd2 is fidelity, and the fill goes where the
vendored algorithm puts it.

### Pass 4: the filing convention

`mmdupd` does not file a vertex under its degree. It files under `dg = dg - qsize[en] + 1`,
floored at 1, where `dg` was the weighted reach INCLUDING en's own members. So the bucket index
is the external degree plus one, and the floor catches a vertex that reaches nothing outside
itself.

`mmdint` meanwhile files at the plain degree, with only the zero case lifted to 1. So MMD runs on
two scales: the initial buckets hold degrees, every refiled bucket holds degree + 1. That is
genuine rather than a misreading, and it tilts the pivot choice slightly against refreshed
vertices, which sit a bucket higher than an untouched vertex of the same reach. From here
degrees[] holds the FILED value, which is what the picker compares and what min_degree tracks;
the nnz(L) accounting is unaffected, since it sums weights over the live members of C[pivot].

### Pass 5: the counters

Two small things in genmmd's main loop, and one of them cannot be checked.

`ncsub`, accumulated per pivot as `mdeg + qsize[mn] - 2`, is the statistic genmmd returns
alongside the permutation, an estimate of the subscript storage the factor will need. Ours prints
it with the other counters. It cannot be compared against the vendored number: `mmd_order`
computes it and drops it, and `genmmd` is static, so the only way to see it would be to edit a
file we do not edit.

The early termination, `if((num+qsize[mn])>neqns)goto n1000`, is checked after the pivot is
numbered and before it is eliminated. When the last supervariable is reached there is nothing
left to update, so genmmd skips the elimination entirely and goes to the numbering. Ours is the
same test on num_eliminated.

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

## amd1: the approximate degree

md5 leaves one expensive thing standing. After each elimination, every reached vertex needs a new
degree, and getting it means uniting the members of every clique in its list and counting the
result. mmd1 attacks the frequency of that union. amd1 attacks its price.

bound(u) is derived in the md2 subsection above, which is where the union it decomposes first
exists, and it becomes usable at md4, where the refresh set narrows to exactly the vertices it
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
the bound at every step, and the closing lines count how often the two disagreed. The exact value
is computed with amd1_exact_degree, which is md5's refresh kept alive for no other purpose. It is
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

**What these two cost in fill, which is not nothing.** amd2 against amd1: 681 against 657 at 10 by
10, 2283 against 2175 at 16 by 16, 4898 against 4762 at 22 by 22, and on the 207 small graphs the
same on 206 and worse on 1. Coarser supervariables mean fewer, larger pivots, so the ordering is
cheaper to produce and slightly worse. That is the trade the vendored routine also makes, and the
comparison to make is against the vendored routine rather than against amd1.

Against the vendored amd_order on the seven examples, both layers already match nnz(L) 7 of 7 and
neither matches the permutation on any, which is the tie-break story the mmd2 section tells at
more length. So these two do not move the acceptance test on their own; they move the mechanism
count. The figures with all seven passes in are at the end of this section.

### Pass 3: the two-pass degree update

The checklist called this "so every survivor sees the same final degme rather than a shrinking
one", and reading `Amd.cpp` showed that description had it backwards. Worth recording, because the
correction is the interesting part.

**We were already two-pass, in both senses the phrase could mean.** Scan 1 computes every
|C[c] - C[p]| before any vertex's bound is formed, and amd1 already did that. And the shrinking
degme is the VENDORED behavior, not ours: `Amd.cpp` detects mass elimination inside scan 2 and does
`degme -= nvi` as it goes, so a survivor handled early sees a larger degme than one handled late.
Our eliminator front-loads mass elimination, so C[p] arrives already trimmed and every survivor
sees the same number. The vendored file flags its own version as a loss in the comment above that
line, one that costs analysis quality rather than ordering quality.

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
definition with plain sets on 310 graphs including grids, and it holds at every step.

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
first refresh corrects in either case.

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
inside reach(u) and make every degree one too large. So amd2_aat builds A[u] with u excluded, which
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

**amd2_preprocess** produces R, the row form of the pattern with duplicates removed, which is the
pattern of A transposed. `Amd.cpp` calls it whenever the input may be unsorted or duplicated, since
R + R' is A + A' and R is clean. The deduplication is the usual stamp, `flag[i] == j` meaning row i
has already appeared in column j.

**amd2_aat** turns the two forms into adjacency lists and reports what the input looked like. Two
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
Info was implemented in pass 4 without being labelled. And `AMD_NDIV` comes out exactly nnz(L)
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

The per-step cost accounting that settles whether all this adds up to parity with the vendored
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
adjacency at every step: the union of the explicit adjacency with the cliques containing u
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
move it forward in md2, because 4 is not tied for the minimum at that step; it is simply not
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
becoming `while not all(eliminated)` since a step can now retire several vertices. Checked
against md3 on the same 400 graphs: zero disagreements. The cost is a constant factor. The
prune loop already touches every neighbor once, the test is O(1) per neighbor, and this is
one more pass over the same set.

For md1 the same alignment is not available cheaply, and the reason is the useful part of
this. In md1's flat graph the natural test is exact: u is indistinguishable from the pivot
when A[u] == neighbors - {u}, meaning everything u still sees lies inside the clique. On the
example above that fires for all four of 1, 3, 4 and 7 at the step that eliminates 0, where
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
in the step that just created the clique. Two vertices can be indistinguishable from each
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
refreshing any degrees, so one expensive update pass serves many pivots. amd1 is the
approximate degree: replace the exact size of the union with a bound computable in one pass.
Those are the ideas. Hash detection, the q2h path and aggressive absorption ride along
because both layers already sweep the reached set and the information is at hand. The
independence runs both ways: the approximate degree works with no hashing at all, which is
what amd1 is, and hash detection could be bolted onto md5
with no approximation anywhere.

The practical consequence is a coarseness ordering. md3, md4 and md5 merge only against the
pivot, so their supervariables are the finest. mmd1's are at least as coarse as ours and
sometimes coarser. amd2's are coarser again where the hash finds pairs the elimination never
brings into the same step. Coarser is not automatically better: a supervariable is a
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
that md3's extra step is free at the moment it is taken.

The gap is that this does not settle the global claim. The two orders diverge, md3 taking u
immediately where md2 takes it several steps later, and from that point the runs face
different graphs and make different subsequent choices. A local exchange argument does not
obviously extend across a run whose pivot sequences have separated, so the equality observed
is stronger than the statement available to justify it.

Two ways to close it. Find the invariant that makes the divergent runs produce identical
filled graphs, which would explain why deferring u never leads md2 into a better or worse
configuration. Or push the search where it is most likely to break: a merge occurring early,
with the deferred vertex separated from its group by several unrelated eliminations before
md2 reaches it. Until one or the other lands, this belongs here as a question and not in the
claims above.

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

The driver loop condition, `while not all(eliminated)`, is an O(n) scan per step, so the
condition alone is quadratic over the run. On a path of 400 vertices it cost 80800 elementary
steps against 1596 of real neighbor work. It is now a counter, incremented by
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
measurement made the defense untenable: 4247 elementary steps on a 20 by 20 grid against 28283
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
later pivot in the same round. Filing, unfiling and popping are O(1).

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
per list, not one per element, and at a few hundred lists per step it is noise. The pool is not
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

### amd2 against Amd.cpp, step by step

The table is the reference version of the argument above, for the completed amd2. It counts only
what a production version would run, so `amd2_exact_degree` is left out: it computes the union the
bound exists to avoid, once per refreshed vertex, and exists only so the trace can print the exact
degree beside the bound.

| step | cost | Amd.cpp |
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

The mapping is mechanical wherever it can be. Names translate from `alive_vertices` to
`aliveVertices`, one identifier to one identifier, and docstrings become the comment block above
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

**ncsub against the vendored value.** It cannot be compared through the public interface, since
mmd_order drops it and genmmd is static. Ours comes from the same expression, which is a reading
rather than a check.

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
bound exists to avoid, once per refreshed vertex, so it dominates. It is instrumentation and would
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

**Where each file spends its work, per step.** `n` is the vertex count, `live` the number not yet
eliminated, and `|C[pivot]|` the size of the clique this step forms. Every picker walks all `n`
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
refresh     mdm2_neighbors(u) for u in C[pivot]      O(|C[pivot]|) unions, after
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
refresh     mdam2_bound(u) for u in C[pivot]         O(|C[pivot]|) additions, after
                                                     bound AGAINST C[pivot]
```

Read down the right column and the two axes are the two words that change: `live` becomes
`|C[pivot]|` going right, and `unions` becomes `additions` going down. The `O(n)` scan is in every
box and is what the third axis removes.

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
in an `outside[]` array indexed by clique id, computed once per step and read by every member of
the group whose incidence list names that clique.

So **the tighter bound is cheap but not free**, which is worth stating plainly: it costs one pass
over the group's incidence lists per step, `sum over u in C[pivot] of |I[u]|`, to fill `outside[]`.
The pivot-free bound costs nothing at all beyond the additions it shares with the other. That is a
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
exact            n unions per step      |C[pivot]| unions per step
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
by walking an array of length n and skipping the eliminated, so the scan is `n` in every box, even in
mdm2 and mdam2 where the *work* has fallen to `|C[pivot]|`. Maintained degrees remove the work, not
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
must produce a number for every live vertex at every step, including vertices this elimination never
touched, and for such a vertex `C[pivot]` is not among `I[u]`'s cliques and there is no group to
state the bound against. So md2's column gets a weaker form,

```
bound(u) = |A[u]| + sum over c in I[u] of ( |C[c]| - 1 )
```

which needs no designated clique. Its terms are exact for two reasons that are md2's doing rather
than the bound's: `u` is in every `C[c]` with `c` in `I[u]`, so `|C[c] - {u}|` is `|C[c]| - 1` with
no test; and `A[u]` is disjoint from all of them, because joining `c` pruned `C[c]` out of `A[u]` and
`A[u]` only ever shrinks after. So its only overcount is one vertex lying in two cliques of `I[u]`.
It is looser than the bound above, which removes the whole of `C[pivot]` from each other clique
rather than just `{u}`.

**Maintenance buys two things, not one**, and they arrive together because they follow from the same
fact: the refresh set narrows from every live vertex to the members of `C[pivot]`.

1. **Fewer vertices worked on**, `|C[pivot]|` instead of all the live ones. This is the obvious one
   and it is what md4 is introduced for.
2. **The tighter bound becomes expressible.** Every vertex being worked on is now a member of
   `C[pivot]`, so a bound stated against `C[pivot]` applies to all of them.

The second also makes that bound CHEAP rather than merely tighter, which is a third thing hiding
inside it. `|C[c] - C[pivot]|` is one number per clique shared by the whole refresh group, obtained
by walking `C[pivot]`'s members and decrementing a counter per clique in each one's incidence list.
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
the elimination. In mdm2 and mdam2 it happens in a refresh phase after it, because `C[pivot]` does
not exist until the eliminator has formed it. So the maintained pair have three phases where the
recomputing pair have two, and a printed number for an untouched vertex is a cached value from
whenever it was last refreshed rather than something computed for the display.

**Results so far.** mda2 against md2 across the seven examples: the order differs on graph3 and
graph4 and `nnz(L)` differs on none. On graph4 one pivot was chosen at a bound of 6 against a true
degree of 3. That is the same shape as the supervariable result in 5.5, permutation differing and
fill not, and a stronger version of it, since mda2 uses the weakest bound in the square with nothing
to fall back on. Seven graphs is not evidence, and mdm2 and mdam2 are not written yet.

## Related

- `archive/sparse_factorization.md` section 5, the prose, pseudocode and worked examples.
- `src/Mmd.cpp` and `src/Amd.cpp`, the vendored routines these are read against.
- `src/OrderEngine.cpp`, the glue that calls them.
