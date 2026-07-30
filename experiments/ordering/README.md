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
| `amd` | approximate degree: a bound instead of a set union | 5.13, 5.14 |

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

The asterisk on `md2 -> md3`: the individual merge is provably free, but the equality of the
totals across the whole run is measured rather than proved. See the open question below.

**`mmd1` and `amd` give up different things.** `mmd1`'s pivots are always true minimum-degree
vertices; only the tie among equals is broken differently, because a batch evicts what it has
touched. `amd`'s pivots may simply not be minimal, because an overcounted bound can hide the true
minimum. MMD perturbs; AMD can be wrong. Both cost well under a percent of fill, in either
direction, and both are noticeably faster.

## The test graphs

`graph1` is a 4-cycle, the smallest graph that fills at all: eliminating any vertex forces its
two neighbors together, for one fill edge. It is also where md3 merges everything that is left
in a single step, 1 taking 2 and 3, which makes it the simplest case for reading the cost of
mass elimination, at the price of the run ending there. `graph2` has six vertices, `graph3` has
twelve and is the first whose ordering is not the identity. `graph4` has eight vertices and
fourteen edges and exists for one reason: **it is the smallest graph we found on which AMD's
bound is ever loose.** The bound overcounts only when a vertex belongs to two elements that
overlap outside the new one, which needs enough eliminations to have made several elements and
enough fill for them to intersect.
Checked exhaustively, no connected graph on five or six vertices is ever loose anywhere in its
run, and none in thirty thousand samples on seven. Without `graph4` the `amd` trace would display
the whole algorithm and never once show it approximating.

`graph5` has five vertices and four edges, two paths joined at 4, and is present in `md1`, `md2`
and `md3` only. It is the smallest graph on which **md3's merge test declines a genuine
supervariable**: at the step whose pivot is 0, vertex 4 has nothing explicit left but belongs to
`c1` as well as to the new clique, so `I[4] == {pivot}` fails even though everything 4 reaches lies
inside that clique. It orders as 2 1 3 0 4 with no merge and no fill, where the exact test would
give 2 1 3 (0 4). It also separates amd's two extra mechanisms: with aggressive absorption on,
`amd` takes four steps and reports `merged = 4, absorbed = c1`; with it off, five steps and no
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
see such a pair. Catching them needs a comparison between candidates, which is amd's hashing.

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
the vendored copy beside the prototype overwrote `amd.cpp` silently, and `git status` said nothing
because `core.ignorecase` is on. The subfolder keeps the two kinds of file apart, and the prefix
means no vendored name can ever collide with a prototype.

## What is not implemented

The `mmd1` and `amd` prototypes are deliberately subsets of the vendored routines, and the plan
for closing the gap is the section that follows. Each file header carries its own list, and
sections 5.11 and 5.13 carry the same lists in prose. In brief:

- **`mmd1`** lacks the prepass that numbers degree-0 and degree-1 vertices before the main loop, and
  `mmdupd`'s `q2h` merging of vertices indistinguishable *from each other* rather than from the
  pivot. It also files degrees at a different offset and never uses bucket 0.
- **`amd`** performs its degree update in one pass where `Amd.cpp` uses two, which changes the
  ordering on four of eleven test graphs. It also lacks dense-row handling and the postorder, both
  of which change the output, along with `amd_aat`, `amd_valid`, the `Control`/`Info` interface and
  the workspace compression, none of which do.

## The plan for mmd and amd

md1 through md5 were a teaching ladder: each isolated one idea and changed exactly one property,
representation, then order, then implementation twice. What is left of the vendored routines is
not a sequence of ideas but the completion of two, so from here the steps are bigger. Two
versions each.

**mmd1, the idea.** Multiple elimination. A batch is an independent set in the current
elimination graph, enforced by evicting every reached vertex from its bucket with a stale
degree, so no later pivot in the batch can be a neighbor of an earlier one. `delta` widens the
batch to near-minima, which is a real concession rather than a free one. Everything else is
md5 unchanged: the quotient graph, mass elimination, the buckets, the expansion.

`delta < 0` belongs here, taking a single pivot per round, which is what `genmmd` does when the
tolerance is negative. It is one line at the bottom of the batch loop, and it is worth having
because it turns the batching off: mmd1 at `delta < 0` should reproduce md5's ordering exactly,
which is a check on the batching rather than a feature of it. The current code loops forever on
a negative delta, so the line has to be written either way.

No weight array in mmd1, for the reason md3 through md5 have none: mass elimination merges only
into the pivot, the pivot dies in the same call, and a supervariable's size is
`len(super_members[pivot])`. Checked rather than assumed, over 200 graphs and 1386
eliminations with batching on: no live vertex ever stands for more than one original vertex.

**mmd2, genmmd complete.** Four additions, all of them holes rather than ideas:

- the PREPASS over `head[1]`, which by `mmdint`'s mapping of degree 0 to degree 1 numbers
  isolated and degree-1 vertices together, before the main loop and without refreshing their
  neighbors;
- `mmdelm`'s `fwd[rn] = nq + 1` stash and `mmdupd`'s split into `q2h` and `qxh`, with the
  pairwise merge inside the `q2h` walk, which is what makes MMD's supervariables coarser than
  ours;
- OUTMATCHED marking, `bwd[nd] = -maxint`, which takes a vertex out of the degree lists without
  merging it;
- the filing convention, `dg - qsize[en] + 1` floored at 1, hence no bucket 0, together with the
  `ncsub` subscript statistic.

The weight array returns in mmd2, and this is where it earns its place: the `q2h` merge folds a
vertex into a LIVE one, so from then on a candidate can stand for several original vertices and
every degree that reaches it has to count them.

**amd1, the ideas.** The approximate degree bound, computed in one pass as now, plus the two
mechanisms that ride along in the same sweep, aggressive absorption and hash supervariable
detection.

**amd2, amd_1 and amd_2 complete.** The two-pass degree update, so every survivor sees the same
final degme rather than a shrinking one, which is the one difference that already shows in our
output; dense row and column detection by the alpha ratio, with those vertices held out of the
ordering and placed last; `amd_aat` forming the pattern of A + A' with the diagonal dropped,
and `amd_preprocess`; `amd_postorder`, so the output is a postorder of the assembly tree rather
than raw elimination order; `amd_valid` as an input check; and the `Control`/`Info` interface.

### What is deliberately excluded

Two pieces of the vendored codes are consequences of packing state into reusable integer arrays
rather than features of the ordering, and neither is modeled:

- MMD's tag and marker machinery with its `maxint` overflow reset. The marks exist because the
  same integer arrays are reused across eliminations; we use sets.
- AMD's workspace with `iwlen`, `pfree` and the `ncmpa` garbage collection. That is the flat
  pool being compacted when it fills.

Both are exactly the kind of thing the C++ rewrite will need and the prototypes do not.

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

The test is `vendored.cpp`, which links `vendored/vendored_mmd.cpp` and
`vendored/vendored_amd.cpp`, copies of `src/Mmd.cpp` and `src/Amd.cpp` that are never edited, and
runs both routines on the same seven graphs, printing permutations in our format. mmd2 and amd2
are accepted when every feature above is present and exercised, nnz(L) matches the vendored
routine on the seven graphs and on random ones, and every remaining order difference is traceable
to a tie.

## Two bugs this found, both ours

Worth recording, because both were invisible to the checks in place at the time.

**`amd` did not shrink the new element on mass elimination.** When a vertex is mass-eliminated it
joins the pivot's supervariable, so it stops being outside the new element and must stop
contributing to `|L|`; `Amd.cpp` does this at `degme -= nvi`. We computed `|L|` once before the
loop. The effect is nearly invisible: identical results on all four test graphs and on every grid,
surfacing only on a five-vertex bowtie where a bound came out one too large.

**`mmd1.cpp` printed display lines its Python twin did not**, left over from the `md5` file it was
derived from. This survived because the verification at the time used a `grep` filter narrow
enough to skip exactly those lines, which is not a test. `make test` exists because of this one:
it compares whole outputs, and it found the drift immediately.

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

## md2: the quotient graph

The insight is that the fill edges created by one elimination are not independent: they form
a clique, and a clique of d vertices can be stored as a d-element list instead of d(d-1)/2
edges. md2 stores it that way and never materializes fill at all.

The representation splits the neighbor relation in two. A[u] holds the vertices u is still
explicitly adjacent to; I[u] holds the ids of the cliques that contain u; C[c] holds the
members of clique c. The true neighborhood is the union, computed on demand by md2_neighbors,
and George and Liu's reachable-set theorem is the guarantee that this union equals what md1's
adjacency would have been. Same degrees, same pivots, same order, on every graph.

Two mechanisms keep it from growing. Pruning: when a new clique is formed, any explicit edge
with both ends inside it is deleted, since the clique now implies it. Absorption: the cliques
the pivot itself belonged to are entirely contained in the new one, so they are deleted
outright. The result is that A[u] only ever shrinks and the clique count stays small, so
total storage falls monotonically where md1's rises.

The price is that the degree is no longer a lookup. Every query unions the explicit adjacency
with every clique the vertex belongs to, and the pivot search does that for every live vertex
at every step. Making the degree cheap again is what md4 and md5 are for.

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
what amd does and what the section on detecting supervariables against each other covers.

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

That is exactly what amd does with it. On graph5, amd with aggressive absorption on takes four
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
printed beside it. Where the weighting first has that effect is amd, whose hash detection
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

## mmd1: multiple elimination

Refreshing degrees is the expensive step, so do it less often: eliminate a whole INDEPENDENT
SET of least-degree vertices before refreshing anything. Non-adjacent pivots cannot disturb
each other's degrees, so every pivot in a batch is still a true minimum-degree vertex when it
is taken. That is Liu's M in MMD, and it is the first layer whose ordering differs from md1's
for a reason other than a tie in the same graph state.

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

## Detecting supervariables against each other, in mmd1 and amd

md3's test is positional: it asks whether a neighbor u is indistinguishable from the PIVOT,
in the step that just created the clique. Two vertices can be indistinguishable from each
other with neither absorbable into the pivot, and no sharpening of a pivot-relative test will
find them. Both mmd1 and amd carry a second mechanism for exactly that population, and the two
mechanisms are different, which is worth recording because the goal is shared and nothing
else about them is.

amd hashes. During the pass it is already making over the reached set, it computes a hash of
each survivor's structure, buckets by it, and runs an exact comparison only within a bucket:

```python
        survivors = [i for i in sorted(L) if not eliminated[i]]
        by_hash = {}
        for i in survivors:
            h = (hash(frozenset(A[i])), hash(frozenset(C[i])))
            by_hash.setdefault(h, []).append(i)
```

then, inside a bucket, `(A[i] - {j}) == (A[j] - {i}) and C[i] == C[j]`. Note that this is the
structural test, so amd inherits its conservatism: a pair whose reachable sets agree while
their element lists differ is missed here too, which is one reason aggressive absorption
travels with the hashing, since removing a contained element is one way such a difference
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
path. Our amd does implement the hash detection, along with aggressive absorption, and the
amd file header calls those out as the two mechanisms beyond md5 that are not about the
degree at all.

Neither mechanism is what its layer is named for, and both are separable from it. mmd1 is
multiple elimination: eliminate a whole independent set of minimum-degree vertices before
refreshing any degrees, so one expensive update pass serves many pivots. amd is the
approximate degree: replace the exact size of the union with a bound computable in one pass.
Those are the ideas. Hash detection, the q2h path and aggressive absorption ride along
because both layers already sweep the reached set and the information is at hand. The
independence runs both ways: the approximate degree works with no hashing at all, which is
what our amd would be with that block deleted, and hash detection could be bolted onto md5
with no approximation anywhere.

The practical consequence is a coarseness ordering. md3, md4 and md5 merge only against the
pivot, so their supervariables are the finest. mmd1's are at least as coarse as ours and
sometimes coarser. amd's are coarser again where the hash finds pairs the elimination never
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
holds. The set algebra that made md2 read as mathematics is gone; that was the price of the
check, and the check is what catches drift. Tags are threaded through return values, since
Python has no reference parameters. One sort remains, at construction, because the input is
given as sets in Python and as ascending literals in C++; it is outside every loop, which is the
same place `SymFactorEngine` puts its final sort.

What is left for the ordering code proper is the shared pool with its garbage collection, and
the weight array. Both are consequences rather than choices: member lists in a pool can outgrow
their slots when a clique grows, which is where AMD's `iwlen`, `pfree` and `ncmpa` come from, and
once a supervariable's members are a chain rather than a list its size stops being O(1) to read,
which is where `weight` earns its place.

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
sides now hold lists and test membership with a stamp, which is what `mmd1` and `amd` need for
their speed and what `SymFactorEngine` already does. What was lost is the notation; what was
gained is that the Python says what the engine does.

## Related

- `archive/sparse_factorization.md` section 5, the prose, pseudocode and worked examples.
- `src/Mmd.cpp` and `src/Amd.cpp`, the vendored routines these are read against.
- `src/OrderEngine.cpp`, the glue that calls them.
