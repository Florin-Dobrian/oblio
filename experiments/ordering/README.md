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
| `mmd` | multiple elimination: a batch of pivots per refresh | 5.11, 5.12 |
| `amd` | approximate degree: a bound instead of a set union | 5.13, 5.14 |

## What the layers show

**`md1` through `md5` return the same ordering.** Every one of them. That is the point of the
first five sections: the heuristic was fixed in `md1` and everything after is implementation, so
the layers can be verified by demanding an identical permutation.

One refinement, because "same ordering" is not quite true across the whole run. Three things can
happen when a layer is added, and all three occur here:

```
                              order        fill        what the change is
md1 -> md2                    same         same        a change of representation
md2 -> md3  (mass elim.)      DIFFERENT    same*       a reordering, free step by step
md3 -> md4 -> md5             same         same        a change of implementation
md5 -> mmd  (multiple elim.)  different     DIFFERENT  a wager
```

So `mmd` is not the first layer to change the permutation. Mass elimination already does, on nine
of twelve test graphs, with identical `nnz(L)` on all twelve. What `mmd` is first to change is the
**fill**.

The asterisk on `md2 -> md3`: the individual merge is provably free, but the equality of the
totals across the whole run is measured rather than proved. See the open question below.

**`mmd` and `amd` give up different things.** `mmd`'s pivots are always true minimum-degree
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

## What is not implemented

The `mmd` and `amd` prototypes are deliberately subsets of the vendored routines. Each file header
carries its own list, and sections 5.11 and 5.13 carry the same lists in prose. In brief:

- **`mmd`** lacks the prepass that numbers degree-0 and degree-1 vertices before the main loop, and
  `mmdupd`'s `q2h` merging of vertices indistinguishable *from each other* rather than from the
  pivot. It also files degrees at a different offset and never uses bucket 0.
- **`amd`** performs its degree update in one pass where `Amd.cpp` uses two, which changes the
  ordering on four of eleven test graphs. It also lacks dense-row handling and the postorder, both
  of which change the output, along with `amd_aat`, `amd_valid`, the `Control`/`Info` interface and
  the workspace compression, none of which do.

## Two bugs this found, both ours

Worth recording, because both were invisible to the checks in place at the time.

**`amd` did not shrink the new element on mass elimination.** When a vertex is mass-eliminated it
joins the pivot's supervariable, so it stops being outside the new element and must stop
contributing to `|L|`; `Amd.cpp` does this at `degme -= nvi`. We computed `|L|` once before the
loop. The effect is nearly invisible: identical results on all four test graphs and on every grid,
surfacing only on a five-vertex bowtie where a bound came out one too large.

**`mmd.cpp` printed display lines its Python twin did not**, left over from the `md5` file it was
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
    for u in sorted(neighbors):
        if not A[u] and I[u] == {pivot}:
            I[u].clear()
            eliminated[u] = True
            merged_vertices.append(u)
    for u in merged_vertices:
        for clique_members in C.values():
            clique_members.discard(u)
```

**The candidates come from the snapshot.** `neighbors` was computed before anything was
touched, so it still lists every member of the new clique even though A and I have since been
rewritten by the prune loop. The `sorted` is not cosmetic: the order of merged_vertices
becomes the order of super_members[pivot], which becomes the order of those vertices in the
returned permutation. Set iteration order would make the output depend on hashing.

**The test reads post-prune state.** By this point A[u] has had `redundant` subtracted and
the pivot discarded, and I[u] has had the absorbed cliques removed and the pivot added. So
`not A[u]` means every explicit neighbor u had was inside the clique, and `I[u] == {pivot}`
means the new clique is its only remaining route out. Together they say u's reachable set is
contained in the clique, which is the fill-free condition.

**There is no weight array.** md3 merges only into the pivot, which is eliminated in the same
call, so no live vertex ever stands for more than one original vertex, and the size of a
supervariable is `len(super_members[pivot])` whenever it is wanted. The literature carries a
weight per vertex and later layers will need one, but here it would be a cached length that
nothing on the hot path reads. It earns its place once the members are held as chains over a
flat array rather than as lists, where a size stops being free.

**Clearing I[u] and the second loop are two halves of one thing.** An incidence is stored
twice, the clique id in I[u] and the member u in C[c]. `I[u].clear()` removes u's side of it;
the second loop removes the cliques' side. They are separate loops because one is per vertex
and the other per clique, not because of any ordering hazard: C[pivot] is a copy of
neighbors, so mutating it cannot disturb the iteration over neighbors.

**A[u] is not cleared, and does not need to be.** The test has already established that it is
empty. Nor does anything have to remove u from other vertices' adjacency, because A is
symmetric: if A[u] is empty then no live v has u in A[v].

**Two things the block deliberately does not do.** It does not touch super_members; the driver
does that from merged_vertices, which keeps the eliminator free of expansion bookkeeping. And
at the instant the first loop finishes, C[pivot] still contains the merged vertices, which is
what the second loop fixes. That is also why the driver computes external_degree from
C[pivot] AFTER the call and not before.

**One conservatism in the second loop.** It scans every clique, `for clique_members in
C.values()`, where `I[u] == {pivot}` already guarantees that u belongs to no clique but
C[pivot], so `C[pivot].discard(u)` would do. The wide scan is defensive against a test that
admits a u with more than one clique, which is exactly what the exact containment test would
do.

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
Through mmd every transfer is into the pivot, so the invariant that a live vertex stands for
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
        refreshed_vertices = sorted(C[pivot])
        for u in refreshed_vertices:
            degrees[u] = len(md4_neighbors(A, I, C, u))
        num_degree_computations += len(refreshed_vertices)
        degrees[pivot] = 0
        for u in merged_vertices:
            degrees[u] = 0
```

The `sorted` is for the trace, not the algorithm: the refresh is order-independent, but the
printed list has to match the C++ twin, whose std::set is sorted by construction.

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
    buckets = [set() for _ in range(n)]        # buckets[d] holds the live degree-d
    for u in range(n):
        buckets[degrees[u]].add(u)
    min_degree = min(degrees) if n else 0
    num_bucket_probes = 0
```

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
        pivot = min(buckets[min_degree])       # lowest index, as the scan did
```

The walk only ever climbs, and min_degree is never reset between steps, so the work is
amortized across the run rather than paid per step. Termination rests on the outer loop
condition: some live vertex exists, it is filed under its own degree, and that degree is at
or above the bound, so a non-empty bucket is found before the array ends.
`min(buckets[min_degree])` is the tie-break, and it is why the orders match md1 through md4:
those scan range(n) keeping the first strict minimum, so ties go to the lowest index, and
this reproduces it exactly.

**Fragment 3, the deletions.**

```python
        buckets[degrees[pivot]].discard(pivot)  # the pivot has left the graph
        degrees[pivot] = 0
        for u in merged_vertices:               # and so have the merged vertices
            buckets[degrees[u]].discard(u)
            degrees[u] = 0
```

This is where md4's harmless zeroing becomes necessary. There a dead vertex's entry was never
read again because eliminated[u] filtered it out of the scan; here there is no scan, and a
dead vertex left in a bucket would be popped as a pivot on a later step. The order within
each pair matters: the bucket index is read from degrees[u], so the removal has to come
before the zeroing, or the vertex is erased from buckets[0] and left where it was.

**Fragment 4, the refresh.**

```python
        refreshed_vertices = sorted(C[pivot])
        for u in refreshed_vertices:
            md5_refile(buckets, degrees, u, len(md5_neighbors(A, I, C, u)))
        num_degree_computations += len(refreshed_vertices)
```

```python
def md5_refile(buckets, degrees, u, new_degree):
    buckets[degrees[u]].discard(u)
    degrees[u] = new_degree
    buckets[new_degree].add(u)
```

Same refresh set as md4, and the only change is that the new degree goes through the helper so
the bucket moves with it. The helper exists so the three steps cannot be written half-way. A
vertex whose degree did not change is removed and reinserted into the same set, which is
harmless. Removal from the middle of a bucket must be O(1), which is why a bucket is a set
here; the vendored codes use doubly linked lists, MMD's fwd/bwd and AMD's Next/Last, because
Fortran and C of that era had no such container.

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
    for u in sorted(neighbors):
        if not A[u] and I[u] == {pivot}:
            I[u].clear()
            eliminated[u] = True
            merged_vertices.append(u)
    for u in merged_vertices:
        for clique_members in C.values():
            clique_members.discard(u)
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

## Detecting supervariables against each other, in mmd and amd

md3's test is positional: it asks whether a neighbor u is indistinguishable from the PIVOT,
in the step that just created the clique. Two vertices can be indistinguishable from each
other with neither absorbable into the pivot, and no sharpening of a pivot-relative test will
find them. Both mmd and amd carry a second mechanism for exactly that population, and the two
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

mmd reaches the same vertices through the vendored code's q2h path. mmdelm stashes each
reached vertex's pruned adjacency count as fwd[rn] = nq+1, and mmdupd routes the nq == 1
cases into a separate list where it merges indistinguishable pairs. Same population, entirely
different route, and no hashing.

Neither is implemented in our mmd, which is why the not-implemented list above names the q2h
path. Our amd does implement the hash detection, along with aggressive absorption, and the
amd file header calls those out as the two mechanisms beyond md5 that are not about the
degree at all.

Neither mechanism is what its layer is named for, and both are separable from it. mmd is
multiple elimination: eliminate a whole independent set of minimum-degree vertices before
refreshing any degrees, so one expensive update pass serves many pivots. amd is the
approximate degree: replace the exact size of the union with a bound computable in one pass.
Those are the ideas. Hash detection, the q2h path and aggressive absorption ride along
because both layers already sweep the reached set and the information is at hand. The
independence runs both ways: the approximate degree works with no hashing at all, which is
what our amd would be with that block deleted, and hash detection could be bolted onto md5
with no approximation anywhere.

The practical consequence is a coarseness ordering. md3, md4 and md5 merge only against the
pivot, so their supervariables are the finest. mmd's are at least as coarse as ours and
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

## Translation choices

The Python is where the thinking happens; the C++ twin is generated alongside it and exists
to keep the Python honest, since two implementations that print the same trace are unlikely
to be wrong in the same way. That makes the twin a correctness instrument and not yet a
performance one. Every choice below was made to keep the two texts recognizably the same
program, and several of them are known to be the wrong choice for speed. Revisiting them is
a separate activity, to be done when the ordering code moves out of this directory, and the
likely conclusion is that most of the sets become plain vectors.

The mapping is mechanical wherever it can be. Names translate from `alive_vertices` to
`aliveVertices`, one identifier to one identifier, and docstrings become the comment block
above the same function. The container mapping is: a list of sets is `std::vector<std::set<int>>`,
a dict from clique id to member set is `std::map<int, std::set<int>>`, a list of pairs is
`std::vector<std::pair<int, int>>`, and a list of lists is `std::vector<std::vector<int>>`.

Six choices are worth recording, since none of them is forced.

**Sets stay sets.** `std::set` is a red-black tree, so a membership test costs O(log n) against
Python's O(1), and a traversal chases pointers between separately allocated nodes. It is used
anyway because it iterates in key order, which is what makes the two traces match without the
C++ having to sort anything, and because `A[u]`, `I[u]` and the clique members are written in
set algebra in the Python. The performance answer later is almost certainly a sorted
`std::vector<int>` with `std::set_intersection` and friends, or an index array over a flat
pool with a scratch mark array, which is what the real minimum degree codes do.

**Cliques are a `std::map`, not an `unordered_map`.** The reason is again iteration order:
`for c in sorted(C)` in the display becomes a plain range-for, because the map is already in
key order. An `unordered_map` would be faster to probe and would then need an explicit sort
at every display.

**Set algebra is spelled out.** `A[u] & neighbors` has no operator, so the C++ builds
`redundant` in a loop and then erases from it. This is the one place the twin is visibly
longer than the Python rather than the same shape.

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

For these prototypes the choice is settled and worth stating. The Python layers use a list of sets,
because the algorithm is written in set algebra and reads that way: `md1` tests `w not in graph[u]`
inside a doubly nested loop over a neighborhood, which is O(1) hashed and O(deg) on a list, and
`md2` says `A[i] & neighbors` and `C[i] -= absorbed` directly. On lists those become membership
scans and stop resembling the mathematics they illustrate. None of this carries to the C++ side,
where the same objects are index arrays over a flat pool: `mmd` and `amd` get their speed from a
scratch mark array and in-place compression, not from hashing, and that is what the real ordering
code will look like.

## Related

- `archive/sparse_factorization.md` section 5, the prose, pseudocode and worked examples.
- `src/Mmd.cpp` and `src/Amd.cpp`, the vendored routines these are read against.
- `src/OrderEngine.cpp`, the glue that calls them.
