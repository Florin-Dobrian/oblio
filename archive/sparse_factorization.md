# Sparse Factorization

Notes on the algorithms behind Oblio's numerical phase. The goal is to keep the
mathematics and the traversal schedules in one place, so the code (left-looking,
right-looking, multifrontal) can be read against a common reference.

Conventions: matrices are `n x n`, indexed `1 .. n` (math style). `A` is the input
matrix, `L` the computed factor. Column indices are `j` and `k`; row indices are `i`.
Throughout, **`k` is the column currently being produced**, and **`j` is the other
column**, to the *left* (`j < k`, an earlier column being gathered) in the
left-looking view, or to the *right* (`j > k`, a trailing column being updated) in the
right-looking view. Sums written `sum_{j<k}` run over `j = 1 .. k-1`.


## 1. Cholesky

### 1.1 The factorization

For a symmetric positive-definite (SPD) matrix `A`, Cholesky computes a lower-triangular
`L` such that

```
A = L L^T
```

`L` is unique once the diagonal is taken positive. Because `A` is symmetric, only its
lower (or upper) triangle is needed, the two triangles carry the same information.

### 1.2 Scalar recurrences

Multiplying `L L^T` and matching entries gives, for the current column `k` and rows
`i > k`:

```
L[k][k] = sqrt( A[k][k] - sum_{j<k} L[k][j]^2 )

L[i][k] = ( A[i][k] - sum_{j<k} L[i][j] * L[k][j] ) / L[k][k]
```

Each entry of `L` depends only on columns `j < k`, everything to its left. Cholesky
therefore sweeps left to right; the two pseudocodes below are just two *schedules* for
performing the same arithmetic, differing in *when* the cross-column contributions are
applied.

The SPD condition surfaces as the square-root arguments staying positive: `A[k][k]`
minus the accumulated sum must be `> 0` at every pivot (each is a Schur-complement
pivot; see 1.6).

### 1.3 Left-looking Cholesky (dot product / "gather")

Column `k` is left untouched until its turn, then it *pulls in* all contributions from
the columns to its left in one pass. Each entry of `L` is written exactly once.
Indexing: `k` is the column being produced, `j` ranges over the earlier columns
(`j < k`) looked at to the left, and `i` ranges over rows.

```
cholesky_left(A) -> L:
    L = zero(n, n)
    for k = 1 .. n:
        for i = k .. n:
            L[i][k] = A[i][k]
    for k = 1 .. n: # current column k
        for j = 1 .. k-1: # earlier column j (to the left)
            for i = k .. n: # row
                L[i][k] = L[i][k] - L[i][j] * L[k][j]
        L[k][k] = sqrt(L[k][k])
        for i = k+1 .. n:
            L[i][k] = L[i][k] / L[k][k]
```

Accumulated over the middle loop, the total update subtracted from entry `L[i][k]` is

```
sum_{j<k} L[i][j] * L[k][j]
```

which is the **dot product of rows `i` and `k`** restricted to the earlier columns
`j < k`. Read this as a *collection* of earlier work: each term `L[i][j] * L[k][j]` is
the contribution that eliminating pivot `j` made to entry `(i,k)`. A term is nonzero
only when *both* `L[i][j]` and `L[k][j]` are nonzero, i.e. when rows `i` and `k` both
have structure in column `j`. So the left-looking update **gathers, at the moment
column `k` is finalized, all the contributions that earlier pivots produced for
`(i,k)`.**

Organized as the loop nest above (`k` then `j` then `i`), each earlier column `j`
contributes to the *whole* of column `k` in one inner sweep, a scaled column
subtraction `col_k <- col_k - L[k][j] * col_j` (an axpy). This is the "column
modification" form, and it is the shape the supernodal left-looking traversal uses:
each contributing column (or supernode) updates the target column (or supernode) as a
vector/block operation, rather than one entry at a time.

### 1.4 Right-looking Cholesky (outer product / "scatter")

The moment column `j` is finished, it *pushes* its contribution into every trailing
column at once. A trailing entry is written many times, once per relevant earlier
pivot, accumulating updates as pivots complete. Indexing: `j` is the column just
finalized, `k` ranges over the later columns (`k > j`) to the right, and `i` ranges
over rows. As in the left-looking version, `A` is only read: the factor is built in `L`.

```
cholesky_right(A) -> L:
    L = zero(n, n)
    for j = 1 .. n:
        for i = j .. n:
            L[i][j] = A[i][j]
    for j = 1 .. n: # current column j
        L[j][j] = sqrt(L[j][j])
        for i = j+1 .. n:
            L[i][j] = L[i][j] / L[j][j]
        for k = j+1 .. n: # later column k (to the right)
            for i = k .. n: # row
                L[i][k] = L[i][k] - L[i][j] * L[k][j]
```

Both algorithms have the same shape: seed `L` from the lower triangle of `A` up front,
then run a single column loop, `A` is read-only, the factor is built in `L`. They
differ only in what that loop does: left-looking **gathers then finalizes** (pull the
earlier columns' contributions into column `k`, then finalize it), right-looking
**finalizes then scatters** (finalize column `j`, then push its contribution into the
later columns). Gather versus scatter is the whole distinction; the seeding is
deliberately identical so the two read as mirror images.

The double loop applies the **rank-1 outer product** of column `j` with itself to the
later columns of `L` (lower triangle only, the factor stays symmetric in structure, so
the upper half is redundant; the `i == k` case updates the diagonal
`L[k][k] -= L[k][j]^2`; this is the shape BLAS `dsyrk` computes):

```
L_trailing <- L_trailing - L[*][j] * L[*][j]^T
```

Entry `(i,k)` of that outer product is `L[i][j] * L[k][j]`, nonzero only when column
`j` has nonzeros in *both* rows `i` and `k`. So finishing column `j` **scatters a
contribution into every later `(i,k)` whose rows both touch column `j`**, updating the
lower triangle only. This is the block operation LAPACK/BLAS accelerate as a symmetric
rank-`k` update (`dsyrk`, with `dgemm` for the off-diagonal supernode blocks) in the
supernodal version.

### 1.5 Left = right = the same contributions, on different schedules

Both algorithms compute the identical `L`. They differ only in *when* the
cross-column contributions are applied:

- **Left-looking gathers (pulls):** many reads from earlier columns, one write per
  entry, at the time the entry's column is finalized.
- **Right-looking scatters (pushes):** one "source" event per pivot, many writes into
  trailing entries, as each pivot completes.

Under the fixed convention `j < k`, the left-looking dot-product term and the
right-looking outer-product entry are literally the **same expression** into the same
target `(i,k)`:

```
    contribution of column j to entry (i,k):   L[i][j] * L[k][j]

    left-looking  (current column k): SUMMED   at read time, over all j < k
    right-looking (current column j): DEPOSITED at write time, into (i,k)
```

Left-looking, finalizing column `k`, gathers `sum_{j<k} L[i][j] L[k][j]` into `(i,k)`.
Right-looking, finalizing column `j`, deposits `L[i][j] L[k][j]` into `(i,k)` for each
later `k`. Same term, same target, one sums at read time, the other deposits at write
time. Pull-when-needed versus push-when-ready.

### 1.6 Connection to vertex elimination and fill

View the nonzero structure of `A` as a graph `G`: one vertex per row/column, and an
edge joining the two indices of each off-diagonal nonzero. Factorization eliminates
vertices in order `1, 2, .., n`.

Across this section the three indices keep fixed roles, matching the `j < k` convention
of the algorithms: **`j` is the vertex being eliminated** (the pivot), and **`i` and `k`
are two later neighbors of `j`** that its elimination connects, so the fill it creates
lands at entry `(i,k)`.

**Eliminating a vertex is a rank-1 Schur-complement update.** Partition `A` around
pivot `j` (the vertex being eliminated):

```
A = [ a_jj   a_j^T ]
    [ a_j    B     ]
```

where `a_j` is the pivot column (nonzero exactly at `j`'s neighbors). Eliminating `j`
replaces the trailing block by its Schur complement:

```
B' = B - (a_j a_j^T) / a_jj
```

The correction `a_j a_j^T / a_jj` is precisely the **rank-1 outer product** of 1.4
(with the finished column playing the role of `a_j`).

So **right-looking factorization is the elimination game, step for step**: its outer
loop `for j = 1 .. n` eliminates vertex `j` at iteration `j`, and the scatter that
follows is exactly the Schur update `- a_j a_j^T / a_jj` pushed into the later columns.
Left-looking computes the same result but *defers* each pivot's effect until the
affected column is reached, the gather at column `k` replays, all at once, the scatters
that pivots `j < k` would have made. Same eliminations, same fill; scatter does them
eagerly, gather does them lazily.

**Fill = the outer product connecting a pivot's neighbors.** The rule in one line: when
`j` is eliminated, if column `j` has nonzeros in both rows `i` and `k` (both `> j`), then
`L[i][k]` fills, equivalently, in the graph `j` is connected to both `i` and `k`, and
eliminating it adds a fill edge between `i` and `k`. This is the `a_j a_j^T` outer
product read entrywise: entry `(i,k)` is `a_j[i] * a_j[k]`, nonzero iff `a_j[i]` and
`a_j[k]` are both nonzero, iff `i` and `k` are *both neighbors of `j`*. If `B[i][k]` was
zero (i and k were not connected), it becomes nonzero: a **fill** edge. Do this for every
pair of `j`'s still-uneliminated neighbors and they all become mutually adjacent, a
clique, and each new clique edge is a fill entry.

**The dot-product term is that same fill, observed later.** In the left-looking view,
the term `L[i][j] * L[k][j]` is nonzero iff both `i` and `k` descend from pivot `j`'s
elimination, i.e. iff `j` was a shared neighbor. So the three statements below are the
same condition in three languages:

1. **Graph:** some earlier vertex `j` is adjacent to both `i` and `k` (eliminating `j`
   connects them).
2. **Outer product / Schur (right-looking):** some pivot column `a_j` is nonzero at both
   `i` and `k`, so `a_j a_j^T` deposits a nonzero at `(i,k)`.
3. **Dot product (left-looking):** some `j < k` has both `L[i][j]` and `L[k][j]`
   nonzero, so that term of `sum_{j<k} L[i][j] L[k][j]` is nonzero.

All three say **"`i` and `k` share an earlier neighbor `j`."** The right-looking outer
product *creates* the fill by scattering; the left-looking dot product *observes* it by
gathering; the graph calls it a shared adjacency.

**Structural, not numerical.** Symbolic factorization predicts fill from *structure*
alone: it declares `(i,k)` filled whenever the patterns of rows `i` and `k` overlap in
some column `j < k`, without evaluating any products. It assumes **no numerical
cancellation**, if the structure admits a nonzero contribution, the position is
treated as nonzero even if the values would happen to cancel. This is why the fill
pattern can be computed ahead of the values, as a graph computation rather than an
arithmetic one.

**Fill depends on the ordering.** The clique formed at each elimination depends on which
vertex is eliminated when, so different permutations of `A` yield different amounts of
fill. Minimizing fill is the job of the ordering phase (AMD / MMD): reorder `A` so the
induced cliques stay small.

### 1.7 Worked examples

Two small examples, both traced by the elimination game of 1.6. Legend: `X` = original
nonzero of `A`, `F` = fill (new nonzero in `L`), `.` = structural zero. Graphs show the
undirected structure of `A` (an edge per off-diagonal nonzero); the diagonal is
implicit. Vertices and matrix indices run `1 .. n`.

#### Example 1, Arrowhead (ordering determines fill)

Vertex 1 is a hub adjacent to all others; there are no other off-diagonal entries. The
graph is a star centered at 1.

```
graph G(A):                 A (natural order 1..6):

     2   3   4                    1  2  3  4  5  6
      \  |  /               1  [  X  X  X  X  X  X ]
       \ | /                2  [  X  X  .  .  .  . ]
        (1)                 3  [  X  .  X  .  .  . ]
       / | \               4  [  X  .  .  X  .  . ]
      5  6                  5  [  X  .  .  .  X  . ]
                            6  [  X  .  .  .  .  X ]
```

Eliminate 1 first. Its later-neighbors are {2,3,4,5,6}, and the rule connects *every
pair* of them, turning the five leaves into a complete graph `K5`. The trailing 5x5
block fills in entirely:

```
after eliminating 1:        L (natural order):

    2 --- 3                       1  2  3  4  5  6
    | \ / |                 1  [  X                ]
    |  X  |   (K5: all      2  [  X  X             ]
    | / \ |    10 pairs     3  [  X  F  X          ]
    6 --- 4    connected)   4  [  X  F  F  X       ]
       5                    5  [  X  F  F  F  X    ]
                            6  [  X  F  F  F  F  X ]
```

That is **10 fill entries**, the maximum possible for this pattern. Every off-diagonal
among rows/cols 2-6 is fill.

Now reorder so the hub is eliminated **last** (elimination order 2,3,4,5,6,1; the hub
becomes the highest-numbered vertex). Each leaf now has only the hub as a later
neighbor, so no elimination ever has two later-neighbors to connect, **no clique
forms, and there is zero fill**:

```
reordered (hub last):       L (reordered; vertex 6 = old hub):

    1                             1  2  3  4  5  6
     \                      1  [  X                ]
    2 \                     2  [  .  X             ]
     \ \                    3  [  .  .  X          ]
    3--(6)*  * = old hub    4  [  .  .  .  X       ]
     / /                    5  [  .  .  .  .  X    ]
    4/                      6  [  X  X  X  X  X  X ]
    5                            (no F anywhere)
```

Same matrix, same graph, only relabeled: **10 fills vs 0 fills.** This is why the
ordering phase (AMD / MMD) exists, for the arrowhead, eliminating the hub last is
optimal, and a good ordering finds that automatically.

#### Example 2, Path with long edges (fill cascades)

A path `1-2-3-4-5-6-7` plus three "long" edges `1-5`, `2-6`, `3-7`. Here fill is not a
single event: each elimination is enlarged by edges the *previous* eliminations
created, so fill propagates down the path.

```
graph G(A):                          A (natural order):

  long edges:                              1  2  3  4  5  6  7
     1 --------- 5                   1  [  X  X  .  .  X  .  . ]
         2 --------- 6               2  [  X  X  X  .  .  X  . ]
             3 --------- 7           3  [  .  X  X  X  .  .  X ]
                                     4  [  .  .  X  X  X  .  . ]
  path:                              5  [  X  .  .  X  X  X  . ]
     1 - 2 - 3 - 4 - 5 - 6 - 7       6  [  .  X  .  .  X  X  X ]
                                     7  [  .  .  X  .  .  X  X ]
```

Eliminating in natural order, step by step (later-neighbors -> new fill):

```
eliminate 1: later {2,5}       -> fill (2,5)
eliminate 2: later {3,5,6}     -> fill (3,5), (3,6)      [5 is present via fill (2,5)]
eliminate 3: later {4,5,6,7}   -> fill (4,6), (4,7), (5,7)  [5,6 present via earlier fill]
eliminate 4: later {5,6,7}     -> (no new fill)          [clique now USES fill edge (5,7)]
eliminate 5: later {6,7}       -> (no new fill)
eliminate 6: later {7}         -> (no new fill)
eliminate 7: later {}          -> done
```

The cascade is explicit: the fill edge `(2,5)` created when eliminating 1 makes 5 a
neighbor of 2, which enlarges 2's clique and produces `(3,5)`; that in turn feeds 3's
clique, and by the time vertex 4 is eliminated the clique is operating on the fill edge
`(5,7)` itself. Fill begets fill. The result is **6 fill entries**:

```
L (natural order):

     1  2  3  4  5  6  7
 1 [ X                    ]
 2 [ X  X                 ]
 3 [ .  X  X              ]
 4 [ .  .  X  X           ]
 5 [ X  F  F  X  X        ]
 6 [ .  X  F  F  X  X     ]
 7 [ .  .  X  F  F  X  X  ]

fill = (2,5), (3,5), (3,6), (4,6), (4,7), (5,7)
```

The elimination tree for this ordering is the single path
`1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7` (each `parent[k] = k+1`), which is exactly the chain
along which the fill propagates: column `k`'s update indices flow to its parent `k+1`,
carrying the accumulated fill forward, the etree union rule of the symbolic phase.

### 1.8 The fill-path theorem

Section 1.6 described fill mechanically, eliminating a vertex cliques its neighbors.
The **fill-path theorem** (Rose-Tarjan-Lueker 1976) turns that into an exact,
value-free test for *which* entries fill, stated purely on the graph of `A`. It is the
tool the later proofs (2.4, 2.5) rest on, so it is worth establishing once here.

Model elimination as a sequence of graphs. Let `G = G_0` be the graph of `A`, and form
`G_k` from `G_{k-1}` by making the neighborhood of vertex `k` a clique and then deleting
`k`; so `G_k` lives on `{k+1, ..., n}`. Cliquing `k`'s neighbors is precisely the fill
`k` creates (1.6). One observation links this sequence to `L`: for `i > j`,

```
L[i][j] != 0   iff   (i,j) is an edge of G_{j-1}
```

the graph just before `j` is eliminated, because `1..j-1` are already gone, `j`'s
neighbors there are exactly its higher-numbered neighbors, i.e. the below-diagonal
nonzeros of column `j` of `L`.

**Lemma (reachability).** For `u, v > k`, `(u,v)` is an edge of `G_k` if and only if `G`
has a path `u - ... - v` whose interior vertices are all `<= k`.

**Proof (induction on `k`).** Base `k = 0`: `G_0 = G` and "interior vertices `<= 0`"
means no interior vertices, i.e. a direct edge, so the claim is `(u,v) in E(G)` iff the
edge exists. Step: assume the claim for `k-1`. `G_k` is `G_{k-1}` with `N(k)` cliqued and
`k` deleted. For `u, v > k`, `(u,v)` is an edge of `G_k` iff either `(u,v)` was already an
edge of `G_{k-1}`, or both `(u,k)` and `(v,k)` were edges of `G_{k-1}` (cliquing joins two
neighbors of `k`). In the first case the hypothesis gives a `u - v` path with interior
`<= k-1 <= k`. In the second, the hypothesis gives a `u - k` path and a `k - v` path each
with interior `<= k-1`; concatenating them yields a `u - v` path whose interior is those
vertices (`<= k-1`) together with the join vertex `k`, so `<= k`. Conversely, take a
`u - v` path in `G` with interior `<= k`. If it avoids `k`, its interior is `<= k-1`, so
`(u,v)` is an edge of `G_{k-1}`, and it survives into `G_k` (elimination only adds edges
among survivors, never removes one). If it uses `k`, split the path at `k` into `u - k`
and `k - v`; each half has interior `<= k-1`, so by hypothesis `(u,k)` and `(v,k)` are
edges of `G_{k-1}`, and cliquing puts `(u,v)` into `G_k`. QED

**Theorem (fill path).** For `i > j`, `L[i][j] != 0` if and only if `G(A)` has a path
`i - ... - j` whose interior vertices are all `< j`.

**Proof.** Specialize the lemma to `k = j - 1`: for `i, j > j-1` (both `>= j`, and `i > j`
so `i, j` are distinct survivors of `G_{j-1}`), `(i,j)` is an edge of `G_{j-1}` iff `G`
has a path `i - ... - j` with interior `<= j-1`, i.e. `< j`. Combined with the observation
`L[i][j] != 0` iff `(i,j)` is an edge of `G_{j-1}`, this is the claim. QED

The reachability lemma is the same "reach through lower-numbered vertices" relation that
reappears in 2.5 as the edges of the update DAG, the fill-path theorem, that
characterization, and the elimination-graph sequence are three views of one fact: **an
edge of the filled graph is a path in `A` that dips only through already-eliminated
vertices.**


## 2. The elimination forest

The elimination forest (often called elimination tree, but in the most general case of a
structurally reducible A matrix, it is a forest) is the combinatorial backbone of sparse
factorization: an `O(n)` structure, one parent pointer per column, that encodes the entire
dependency structure of `L`. Which columns update which, where fill lands, and what can run
in parallel, all read off it. Every later phase (symbolic factorization, supernode
detection, scheduling) is a walk over this forest.

### 2.1 Definition

Vertices, rows, and columns coincide (the matrix is square with symmetric structure), so
"vertex `j`", "row `j`", and "column `j`" name the same object. A column's
below-diagonal pattern can be written

```
Struct(j) = { i : i > j and L[i][j] != 0 }
```

basically the rows where column `j` of the factor has nonzeros beneath the diagonal. The
**parent** of `j` is the first of them:

```
parent(j) = min Struct(j)          (NIL if Struct(j) is empty)
```

That single array `parent[1..n]` *is* the forest. A column with no below-diagonal
nonzero is a **root**. Why *forest* and not *tree*: if `A` is structurally reducible,
its graph has several connected components, each component contributes its own tree, so
the general object is a forest; an irreducible matrix gives a single tree. Read off `L`,
`parent(j)` is simply the row of the topmost off-diagonal nonzero in column `j`, and the
height of the forest is the length of the longest child-to-parent chain.

### 2.2 Computed from A, used to build L

If we already had `L`, the forest would be trivial, one scan per column for its first
below-diagonal nonzero. But that is backwards in practice: computing `L` is the goal, and
we want the forest *first*, precisely so we can build `L`. Reading `parent(j)` off `L`
needs the very thing it is meant to help produce.

The escape is that `parent(j)` never depends on the *values* in `L`, only on where its
nonzeros are, and that structure is recoverable from `A`'s structure alone, without ever
forming `L`. So the dependency runs one way, not in a circle:

```
A (structure)  --climb-->  forest  --union up the tree-->  Struct(L)  -->  numeric L
```

The forest is the pivot of the whole architecture: derivable from `A` cheaply
(near-linear), and once we have it, `L`'s pattern is a tree walk (symbolic
factorization) and the numeric values are scheduled along it. The deeper statement is

```
Struct(L) is a function of Struct(A) alone
```

independent of the numeric values (assuming no lucky cancellation, the "structural, not
numerical" principle of 1.6). The forest is the compact `O(n)` witness of that function.
This is what lets the ordering phase score fill and the symbolic phase allocate storage
before a single floating-point operation.

### 2.3 The containment theorem

One theorem justifies both computing the forest and using it. It is the clique fact of
1.6, restated for columns.

**Theorem (column containment).** If `Struct(j)` is nonempty, let `k = parent(j)` (so
`k > j`). Then

```
Struct(j) \ {k}   is a subset of   Struct(k)
```

equivalently `Struct(j)` is a subset of `{k} union Struct(k)`: a column's pattern is
contained in its parent, together with the parent's pattern.

**Proof.** The elements of `Struct(j)` are `j`'s higher-numbered neighbors, and
eliminating `j` makes them pairwise adjacent (the clique of 1.6). Since `k` is the
smallest of them, every other `i` in `Struct(j)` has `i > k` and becomes a neighbor of
`k`, so `L[i][k] != 0`, i.e. `i` is in `Struct(k)`. QED

Two corollaries, the two ways the theorem is used:

**Corollary 1 (nonzeros are ancestors).** `Struct(j)` is a subset of `j`'s proper
ancestors, every below-diagonal nonzero of column `j` lies on the path from
`parent(j)` to the root.

**Proof (induction up the tree).** The theorem is a single hop: `Struct(j) \ {k}` is a
subset of `Struct(k)` for `k = parent(j)`. Chain the hops. Process columns from the
roots downward, so that `k = parent(j)` is handled before `j`; this is well-founded
because `parent(j) > j` strictly, so the recursion always moves to a higher-numbered
column and terminates at the roots. Base case: if `Struct(j)` is empty, `j` is a root
and the claim holds vacuously. Inductive step: assume the claim for `k = parent(j)`.
Every `i` in `Struct(j)` with `i != k` lies in `Struct(k)` by the theorem, hence is a
proper ancestor of `k` by the hypothesis; since `k` is an ancestor of `j` and ancestry
is transitive, `i` is a proper ancestor of `j`. The remaining element `i = k` is
`parent(j)`, a proper ancestor of `j` by definition. QED

So a column's pattern can only ever point at ancestors, which is why the forest can be
*computed by climbing* (2.4): to place a column's connections we walk upward, never
sideways or down.

**Corollary 2 (the symbolic union recurrence).** Read the other direction, a column
absorbs its children's patterns:

```
Struct(j) = { i > j : A[i][j] != 0 }
            union, over each child c of j, ( Struct(c) \ {j} )
```

column `j`'s structure is `A`'s column `j` unioned with each child's structure minus that
child's own entry. This is the recurrence symbolic factorization runs, the
"union-with-marks" loop, in one bottom-up pass over the forest.

**Proof.** Write `Acol(j) = { i > j : A[i][j] != 0 }` for the original below-diagonal
entries of column `j`, and let `RHS = Acol(j) union ( over children c of j: Struct(c) \
{j} )` be the whole right-hand side of the recurrence. The claim is `Struct(j) = RHS`;
prove the two inclusions `RHS is a subset of Struct(j)` and `Struct(j) is a subset of
RHS`.

*`RHS is a subset of Struct(j)`.* Each original nonzero is a nonzero of `L`, so
`Acol(j)` is a subset of `Struct(j)`. For a child `c` of `j`, `parent(c) = j`, so the
containment theorem gives `Struct(c) \ {j}` a subset of `Struct(j)`. Every set making up
`RHS` therefore lies in `Struct(j)`, so their union does too.

*`Struct(j) is a subset of RHS`.* Take `i` in `Struct(j)`, so `i > j` and
`L[i][j] != 0`. If `A[i][j] != 0` then `i` is in `Acol(j)`, done. Otherwise `(i,j)` is
fill: by the elimination mechanism of 1.6 a fill edge `(i,j)` is created by eliminating
some pivot `c < j` whose later-neighbors include both `i` and `j`, that is, some `c < j`
with both `i` and `j` in `Struct(c)`. Since `j` is in `Struct(c)`, Corollary 1 makes `j`
a proper ancestor of `c`, so the tree path from `c` upward reaches `j`; let `c'` be the
child of `j` on that path (`parent(c') = j`, with `c, ..., c'` all `< j` and strictly
increasing up the path). Carry `i` up from `c` to `c'` with the theorem: at each vertex
`x` on the path `parent(x) <= j < i`, so `i != parent(x)`, hence `i` is in
`Struct(x) \ {parent(x)}` and therefore in `Struct(parent(x))`; iterating lands `i` in
`Struct(c')`. As `i > j` we have `i != j`, so `i` is in `Struct(c') \ {j}`, one of the
sets in `RHS`. QED

Same theorem: Corollary 1 climbs (justifies the algorithm below), Corollary 2 accumulates
(drives symbolic factorization).

### 2.4 Computing the elimination forest from A

Sweep columns left to right; for each earlier neighbor, climb to the root of its subtree
and attach. The forest is *defined* from `L`, but this reads only `A` and the partial
forest, the justification (that the sweep still yields the `L`-defined forest) is the
incremental argument below. Conventions: 1-indexed, `NIL` = root.

**Plain climb**, `parent[]` is the forest; each climb walks it.

```
computeForest(A) -> parent:
    for k = 1 .. n:
        parent[k] = NIL
    for k = 1 .. n: # current column k
        for each j < k with A[k][j] != 0: # earlier neighbor j
            r = j
            while parent[r] != NIL: # climb j's subtree to its root
                r = parent[r]
            if r != k:
                parent[r] = k # attach that root under k
```

Correct and transparent, but the climb is `O(height)` per neighbor; a long chain (a
path-shaped forest) makes it `O(n * nnz)` worst case, re-walking the same lengthening
chain at every step.

**Access orientation.** The plain-climb pseudocode is written with `A[k][j]`, `j < k`,
lower-triangle flavored, but it assumes nothing about how `A` is stored; it is purely
algorithmic. The inner test `{ j < k : A[k][j] != 0 }` is the lower-triangle structure of
**row `k`**, so the forest is a *row-oriented* sweep. That makes it the exception among
the phases: symbolic factorization (3.1) and numeric factorization march over columns
(`Acol(j) = { i > j : A[i][j] != 0 }` and its numeric analogue), so the practical storage
format for a sparse solver is compressed-sparse-column (CSC). The forest is the one pass
that wants rows.

Full storage reconciles this with a single subscript. If `A` is stored fully (both
triangles), the part of each column above the diagonal is, by symmetry, a row below the
diagonal: row `k`'s earlier neighbors `{ j < k : A[k][j] != 0 }` are exactly column `k`'s
above-diagonal entries `{ j < k : A[j][k] != 0 }`. So under full CSC the only change the
pseudocode needs, for efficient access, is to read `A[j][k]` (the above-diagonal of stored
column `k`) in place of `A[k][j]` (row `k`), everything else is identical. We write
`A[k][j]` for clarity of exposition and use `A[j][k]` for an efficient implementation.
This is Oblio's approach, full storage of A.

*Lower-triangle-only.* The forest can still be built from lower-only storage at the same
`O(nnz * alpha(n))` cost, but not by the subscript swap above. Path compression is correct
only if edges are processed in increasing order of their *higher* endpoint (so that when a
subtree is attached under `k`, every smaller vertex is already resolved). Lower CSC stores
each edge `(i, j)`, `i > j`, in column `j`, grouped by its *lower* endpoint, so a plain
column sweep presents higher endpoints out of order and gives the wrong tree. Regrouping
edges by higher endpoint first (a transpose) fixes it:

```
for k = 1 .. n: earlier[k] = empty                 # bucket per higher endpoint
for j = 1 .. n:                                     # sweep lower-triangle columns
    for each i > j in stored column j:              # edge (i,j), lower part
        append j to earlier[i]                      # file it under higher endpoint i
for k = 1 .. n:                                     # now standard climb
    for each j in earlier[k]:
        climb from j, attach root under k
```

Same result, same complexity, but two passes and a scatter into `earlier[i]` whose write
target jumps with the higher endpoint, poorer locality than a single streaming column
walk, so a potential cache cost on large problems. Reading the upper triangle of each
column is in fact the standard route (Liu 1986; CSparse's `cs_etree` in Davis 2006 reads
column `k`'s entries with row index below `k`); code that stores only the lower triangle
transposes to that orientation first, which is the regrouping above. Oblio takes the
full-storage route instead: it is the clean, pragmatic choice, one column-oriented access
reused across the etree, symbolic, and numeric phases, no transpose or bucketing, and the
price, storing both triangles, is one Oblio accepts. The lower-only variant is recorded
here as an alternative, not the path taken.

**Why sweeping A gives the L-defined forest.** The forest is *defined* from `L`
(`parent(j) = min Struct(j)`), yet the algorithm reads only `A` and the partial forest.
Two arguments justify this; the first is Liu's and is shorter but leans on an external
theorem, the second is self-contained within this document.

*Proof 1 (path characterization; Liu 1990, on the fill-path theorem of 1.8).* The
fill-path theorem (1.8) states: for `i > j`, `L[i][j] != 0` iff there is a path
`i - ... - j` in `G(A)` whose interior vertices are all numbered `< j`. Substitute this
into the definition:

```
parent(j) = min { i > j : G(A) has a path i - ... - j with all interior vertices < j }
```

The right-hand side mentions only the edges of `A` and the vertex numbering, no `L`, no
fill. So the forest is manifestly a function of `Struct(A)` alone, and the sweep computes
it directly. In this light the algorithm is disjoint-set union: adding `k` merges the
subtrees (components) of `k`'s lower-numbered `A`-neighbors and makes `k` their root,
which is exactly the climb-and-attach, and is why the `ancestor[]` shortcut below is a
union-find with the `alpha(n)` bound. This is clean, but it *cites* the fill-path theorem
rather than proving it. QED

*Proof 2 (incremental; self-contained, using Corollaries 1-2).* Here the fill-path
theorem is not assumed; only this document's own results are used. Let `A_k = A(1:k, 1:k)`
be the leading `k x k` submatrix and `F_k` its elimination forest (defined from *its own*
factor); the full forest is `F_n`. Write `P` for the array the algorithm produces.

**Lemma A (Cholesky factors nest).** The top-left `(k-1) x (k-1)` block of the factor of
`A_k` equals the factor of `A_{k-1}`.

**Proof.** Partition `A_k = [[A_{k-1}, b], [b^T, a]]`. Its Cholesky factor is
`[[L_{k-1}, 0], [c^T, d]]` with `L_{k-1} c = b` and `d = sqrt(a - c^T c)`, forced by the
block product, and its top-left block is `L_{k-1}`, the factor of `A_{k-1}`. QED

So passing from `A_{k-1}` to `A_k` adds vertex `k` (the largest index), and for every
`j < k` the only possible new nonzero in column `j` is at row `k`:

```
Struct_{A_k}(j) = Struct_{A_{k-1}}(j), possibly with k added,   for j < k.
```

Two consequences for parents. (1) If `j` was **not** a root of `F_{k-1}`, its
`min Struct(j)` is already `< k`, so adding the larger `k` cannot change it: existing
parent pointers never move. (2) A **root** `r` of `F_{k-1}` either stays a root or gets
`parent_{F_k}(r) = k` (the only larger index available). Growing the matrix by one vertex
can only hang former roots under `k`, the forest is built once, never rewritten. This is
exactly why `L` is not needed.

**Lemma B (which roots attach).** Let `r` be a root of `F_{k-1}`. Then
`parent_{F_k}(r) = k` **iff** `r`'s subtree in `F_{k-1}` contains some `j` with
`A[k][j] != 0`.

**Proof.** (if) Let `j` in `r`'s subtree have `A[k][j] != 0`. Then `k` is in
`Struct_{A_k}(j)`, so by Corollary 1 (in `A_k`) `k` is a proper ancestor of `j`. Climbing
`j`'s ancestor chain follows `F_{k-1}` while below `k` (consequence 1), reaching `r`;
since `k` must be an ancestor of `j`, the chain continues past `r`, forcing
`parent_{F_k}(r) = k`. (only if) If `parent_{F_k}(r) = k` then `k` is in `Struct_{A_k}(r)`;
apply Corollary 2 in `A_k` repeatedly, `k` in `Struct(x)` means either `A[k][x] != 0`
(an `A`-neighbor in `r`'s subtree) or `k` in `Struct(c) \ {x}` for a child `c`, and
descending to `c` stays in `r`'s subtree. The descent terminates at some `d` in `r`'s
subtree with `A[k][d] != 0`. QED

So the children of `k` in `F_k` are exactly the roots of `F_{k-1}` whose subtrees touch
`k` in `A`. The algorithm realizes precisely this: assuming `P[1:k-1] = F_{k-1}` before
step `k`, it climbs from each `A`-neighbor `j < k` to the root of `j`'s subtree and sets
`P[r] = k`. Every attached root is Lemma-B-certified (some neighbor lies in its subtree),
every such root is reached, and duplicates climb to the same root and attach once (the
`r != k` guard); all other pointers stay put (consequence 1). Hence `P[1:k] = F_k`. The
base case `k = 1` is the empty forest, so by induction `P = F_n`. The sweep needs only
`k`'s `A`-neighbors and the partial forest, fill is never consulted, because the
re-rooting rule depends on `A`'s edges and the tree so far, not on `L`. QED

Both proofs are now self-contained within this document; they differ in which prior
result they lean on. Proof 1 is short because it stands on the fill-path theorem of 1.8
(itself one induction), which makes the `A`-only dependence immediate. Proof 2 is longer
but needs only the containment corollaries of 2.3. They agree, and either certifies that
the plain-climb sweep, and Liu's path-compressed version below, computes the true
elimination forest from `A`.

**Liu's algorithm (path compression)**, keep the true forest in `parent[]`, but climb
along a second array `ancestor[]`, a union-find shortcut flattened toward `k` each pass.

```
computeForest(A) -> parent:
    for k = 1 .. n:
        parent[k]   = NIL
        ancestor[k] = NIL
    for k = 1 .. n: # current column k
        for each j < k with A[k][j] != 0: # earlier neighbor j
            r = j
            while ancestor[r] != NIL and ancestor[r] != k:
                t = ancestor[r]
                ancestor[r] = k # compress: point r straight at k
                r = t
            if ancestor[r] == NIL:
                ancestor[r] = k
                parent[r] = k # true forest edge
```

The two arrays do different jobs: `parent[]` is the forest we keep, `ancestor[]` is
scratch we discard. Path compression drops the cost to near-linear,
`O(nnz * alpha(n))`. Both algorithms produce the identical `parent[]`, the plain form
is the one to reason about for correctness, Liu's is what production code runs.

**The two end states, and four ways to test them.** Every climb finishes in one of two
cases. Either `j` reaches a fresh subtree root, which then takes `k` as its parent
(**attach**), or `j` already lies under `k` from an earlier neighbor this same pass, so
nothing is added (**skip**). Detecting the skip can live either in the loop's exit or in
the attach test, which yields four equivalent forms, two climbing `parent[]` (plain) and
two climbing `ancestor[]` (compressed). Writing `X` for the climbed array:

```
form                        loop continues while         attach when
--------------------------  ---------------------------  ------------------
plain,    single (classic)  parent[r] != NIL             r != k
plain,    double            parent[r] != NIL and != k    parent[r] == NIL
compress, single            ancestor[r] != NIL           r != k
compress, double (classic)  ancestor[r] != NIL and != k  ancestor[r] == NIL
```

All four are correct because `X[k]` is NIL throughout iteration `k`: a column is linked
only by a later column, so its own slot is empty during its turn, and the walk may run
onto `k` without corrupting it. Trace the two end states:

- **Attach** is identical in every form. The climb halts at a node with `X[r] == NIL` that
  is not `k`, and writes `parent[r] = k`.
- **Skip** is where they differ. The single-condition loop climbs all the way to the
  subtree root, which on a skip *is* `k`, so it detects the skip as `r == k` and the attach
  test `if r != k` rejects it. The double-condition loop stops one node short, at the vertex
  already pointing at `k` (`X[r] == k`), so the skip becomes the loop's second exit and the
  attach test `if X[r] == NIL`, false here, rejects it. In the double-condition form the
  attach test is really reading back which exit fired.

The tests cannot be swapped: `if r != k` after a double-condition loop would attach at a
non-root and overwrite a real parent, and `if X[r] == NIL` after a single-condition loop is
never true at `k`, so it would never skip. Loop form and attach test are one choice.

The classic forms sit on the diagonal, plain with the single condition (the plain-climb
pseudocode above) and compressed with the double (Liu, CSparse, Oblio). We follow that
convention, but it may be historical more than principled. Each classic is the natural spelling
of its own algorithm: plain reads as "climb to the root, then test whether it is `k`",
while compressed is a union-find `find` with path compression, whose idiomatic early-out is
"stop on reaching a node already linked to `k`" (`ancestor[r] == k`). There is one genuine
but negligible asymmetry: the extra step a single-condition loop takes onto `k` is a
redundant *write* in the compressed case (`ancestor[r] = k` again) and only a *read* in the
plain case, so the double condition has marginally more to recover for compressed. That
aside, the case for stopping one node short applies to both climbs equally and does not
explain why only the compressed classic uses it. For exposition the symmetric pairing would
be clearer, since the conventional split carries no lesson; the content is the 2x2 above
and the loop-condition/attach-test coupling. Oblio ports the classic compressed form.

### 2.5 The forest as transitive reduction of the update DAG

A clean characterization ties the forest to the factorization dependencies of 1.6. Direct
the "who-updates-whom" relation: an edge `j -> k` whenever eliminating `j` sends a
contribution to `k`, i.e. for every `k` in `Struct(j)`. Every edge points forward in the
elimination order (`j < k`), so this is a **DAG**, the filled graph, directed by
elimination order. Its edges are the *reach* relation: `j -> k` iff `k` is reachable from
`j` through lower-numbered vertices (the fill paths).

The forest keeps only the minimal out-edge of each node, `j -> parent(j)`; every other
edge is transitively implied:

```
j -> k in the DAG   <=>   k is an ancestor of j in the forest
```

and each dropped edge `j -> k` is recovered by climbing from `parent(j)` up to `k`,
exactly the containment theorem. That is the definition of **transitive reduction**: the
unique minimal edge set with the same reachability. So

```
elimination forest = transitive reduction of the filled (update) DAG
                   = Hasse diagram of the elimination-reach order
```

The reduction is unique because the filled graph is a DAG. Two directions close the
circle: the plain-climb algorithm *computes* the reduction (climbing past
transitively-implied edges to the minimal one); symbolic factorization *inverts* it
(unioning up the tree to reconstruct the full reach, `Struct(L)`). Reduction to store,
closure to use. The DAG must be the *filled* graph, not `A`'s graph, the reach relation
runs over paths through eliminated vertices, so fill edges are essential.

### 2.6 Worked example: nested dissection of a grid

A 3x3 grid, small enough to trace, rich enough to branch, and a real sparse pattern.
Nine nodes with grid adjacency; the ordering is the whole story, so take the good one,
**nested dissection**: cut the grid with its middle column as a separator, number the
left block first (1,2,3), the right block second (4,5,6), the separator last (7,8,9).

```
   left   sep   right
    1  --  7  --  4
    |      |      |
    2  --  8  --  5
    |      |      |
    3  --  9  --  6
```

Structure of `A` (`X` = nonzero; full symmetric):

```
     1 2 3 4 5 6 7 8 9
  1  X X . . . . X . .
  2  X X X . . . . X .
  3  . X X . . . . . X
  4  . . . X X . X . .
  5  . . . X X X . X .
  6  . . . . X X . . X
  7  X . . X . . X X .
  8  . X . . X . X X X
  9  . . X . . X . X X
```

The two subdomains `{1,2,3}` and `{4,5,6}` are the diagonal blocks; they have **no
coupling** (the `4:6 x 1:3` block is empty), interacting only through the separator band
`{7,8,9}`.

Compute the forest by the climb of 2.4, the multi-hop climbs are fill reconstructed from
`A` without forming `L`:

```
k=1: no earlier neighbors
k=2: earlier neighbor j=1                        -> parent[1]=2
k=3: earlier neighbor j=2                        -> parent[2]=3
k=4: no earlier neighbors
k=5: earlier neighbor j=4                        -> parent[4]=5
k=6: earlier neighbor j=5                        -> parent[5]=6
k=7: earlier neighbor j=1, climb 1->2->3, attach -> parent[3]=7 fill
k=7: earlier neighbor j=4, climb 4->5->6, attach -> parent[6]=7 fill
k=8: earlier neighbor j=2, climb 2->3->7, attach -> parent[7]=8
k=8: earlier neighbor j=5, climb 5->6->7->8, reaches k=8, skip
k=8: earlier neighbor j=7, climb 7->8, reaches k=8, skip
k=9: earlier neighbor j=3, climb 3->7->8, attach -> parent[8]=9
k=9: earlier neighbor j=6, climb 6->7->8->9, reaches k=9, skip
k=9: earlier neighbor j=8, climb 8->9, reaches k=9, skip
```

The trace has one line per off-diagonal nonzero of a single triangle, each edge seen once
at its higher endpoint `k`, plus one line per column with no earlier neighbor. Those
columns are exactly the forest leaves, so the length is (edges) + (leaves): here 12 + 2 =
14.

The `skip` lines are the guard `if r != k` of 2.4 firing. Once one earlier neighbor of
`k` has attached its subtree under `k`, the remaining earlier neighbors of `k` already sit
in that subtree, so their climb runs straight into `k` with nothing to attach. Skips
appear only at the separators `7,8,9`, where several chains meet; each chain column has a
single earlier neighbor and never skips.

The same run under path compression, the `ancestor` shortcut of 2.4: each climb redirects
the `ancestor` of every node it passes straight to `k`, so a later climb reaches the root
in fewer hops. The `parent[]` result is identical; only the walk lengths change.

```
k=1: no earlier neighbors
k=2: earlier neighbor j=1                        -> parent[1]=2, ancestor[1]=2
k=3: earlier neighbor j=2                        -> parent[2]=3, ancestor[2]=3
k=4: no earlier neighbors
k=5: earlier neighbor j=4                        -> parent[4]=5, ancestor[4]=5
k=6: earlier neighbor j=5                        -> parent[5]=6, ancestor[5]=6
k=7: earlier neighbor j=1, climb 1->2->3, attach -> parent[3]=7 fill, ancestor[1,2,3]=7
k=7: earlier neighbor j=4, climb 4->5->6, attach -> parent[6]=7 fill, ancestor[4,5,6]=7
k=8: earlier neighbor j=2, climb 2->7, attach    -> parent[7]=8, ancestor[2,7]=8
k=8: earlier neighbor j=5, climb 5->7, skip (already at k), ancestor[5]=8
k=8: earlier neighbor j=7, skip (already at k)
k=9: earlier neighbor j=3, climb 3->7->8, attach -> parent[8]=9, ancestor[3,7,8]=9
k=9: earlier neighbor j=6, climb 6->7, skip (already at k), ancestor[6]=9
k=9: earlier neighbor j=8, skip (already at k)
```

The shortcut shows at `k=8`: the walk from `2` reaches root `7` in one hop (`2->7`), not
`2->3->7`, because `k=7` had already pointed `ancestor[2]` at `7`.

The same fill, seen as the elimination game of 1.6 run in label order: eliminating `k`
cliques its higher-numbered neighbors (the still-uneliminated ones), and any clique edge
not already present is fill.

```
eliminate | higher neighbors | new fill edges  | total fill
----------+-------------------+-----------------+-----------
    1     | {2,7}             | (2,7)           |   1
    2     | {3,7,8}           | (3,7) (3,8)     |   3
    3     | {7,8,9}           | (7,9)           |   4
    4     | {5,7}             | (5,7)           |   5
    5     | {6,7,8}           | (6,7) (6,8)     |   7
    6     | {7,8,9}           | -               |   7
    7     | {8,9}             | -               |   7
    8     | {9}               | -               |   7
    9     | {}                | -               |   7
```

The higher-neighbor set already reflects prior fill, e.g. at vertex 3 the neighbors are
`{7,8,9}`, where `7` and `8` arrived as fill from eliminating 2. Three things this makes
visible:

- **Fill is created only at vertices 1-5** (the subdomain columns); eliminating 6-9 adds
  nothing. Every fill lands in the separator rows 7-9, in the "new fill edges" column all
  endpoints are `>= 7` except the smaller endpoints 2,3,5,6, which are the subdomain
  vertices reaching up into the separator.
- **The two subtrees are independent and merely interleave by label.** Fills from
  `{1,2,3}` (rows 1-3: `(2,7),(3,7),(3,8),(7,9)`) and from `{4,5,6}` (rows 4-5:
  `(5,7),(6,7),(6,8)`) never interact, no fill edge joins a `{1,2,3}` vertex to a
  `{4,5,6}` vertex. That is the sibling independence, seen dynamically.
- **`(7,9)` is second-order fill.** Eliminating 3 creates `(7,9)`, but 3 only became
  adjacent to 7 through the earlier fill `(3,7)` (from eliminating 2). So it is fill
  begetting fill, confined to the separator.

giving `parent = [2, 3, 7, 5, 6, 7, 8, 9, NIL]`, a forest that branches:

```
              9            separator top (root)
              |
              8
              |
              7            the two blocks MERGE here
            /   \
           3     6
           |     |
           2     5
           |     |
           1     4
        (left)  (right)
```

Structure of the factor `L` (`F` = fill), 7 fills:

```
     1 2 3 4 5 6 7 8 9
  1  X
  2  X X
  3  . X X
  4  . . . X
  5  . . . X X
  6  . . . . X X
  7  X F F X F F X
  8  . X F . X F X X
  9  . . X . . X F X X
```

Everything reads off this one figure:

- **It branches -> parallelism.** The subtrees `{1,2,3}` and `{4,5,6}` are disjoint,
  neither an ancestor of the other, so by sibling independence they factor completely in
  **parallel**, meeting only at the separator `7`. In the matrix this is the empty
  `4:6 x 1:3` block: the two halves never write into each other. That branch is the
  divide-and-conquer of nested dissection made visible.
- **Chains -> supernodes.** Each straight run (`1->2->3`, `4->5->6`, `7->8->9`) is a chain
  of single-child nodes whose patterns nest, the raw material for fundamental supernodes.
- **Climbs -> fill, localized.** All 7 fills land in the separator rows 7-9; the subdomain
  blocks stay clean. The good ordering shows up as fill *quarantined* into the separator
  corner, the entire point of nested dissection.
- **First nonzero = parent.** Each column's topmost below-diagonal entry is its parent:
  col 3 -> row 7 (fill), col 6 -> row 7 (fill), col 7 -> row 8, reproducing `parent[]`
  straight from the matrix.

Finally, mark each forest edge as original or fill. A forest edge is `j -> parent(j)`; it
is **original** if `A[parent(j)][j] != 0`, otherwise **fill**:

```
forest edge | A[parent][j] | kind
------------+--------------+----------
  1 -> 2    | A[2][1] = X  | original
  2 -> 3    | A[3][2] = X  | original
  3 -> 7    | A[7][3] = .  | fill
  4 -> 5    | A[5][4] = X  | original
  5 -> 6    | A[6][5] = X  | original
  6 -> 7    | A[7][6] = .  | fill
  7 -> 8    | A[8][7] = X  | original
  8 -> 9    | A[9][8] = X  | original
  9 -> root |     -        | (root)

original: 1-2, 2-3, 4-5, 5-6, 7-8, 8-9   fill: 3-7, 6-7
```

- **The two fill forest edges `3-7` and `6-7` are exactly the merge edges** where the left
  and right chains roll up into the separator, the two multi-hop climbs at `k=7` in the
  trace. The within-subdomain chain edges (`1-2, 2-3, 4-5, 5-6`) and the separator chain
  edges (`7-8, 8-9`) are all original. So the forest's original edges follow `A`'s actual
  chains, and its fill edges are precisely where independent subtrees join.
- **A forest edge is fill only when a column's first below-diagonal nonzero is itself
  fill**, here only columns 3 and 6, the subtree roots whose smallest separator-neighbor
  (`7`) is reached through fill rather than an original edge. For every other column the
  minimum of its pattern is an original neighbor.
- **Most fill is not on the tree.** Of the 7 fill entries in `L`, only 2 are forest edges
  (`3-7`, `6-7`); the other 5, `(2,7), (3,8), (5,7), (6,8), (7,9)`, are non-tree fill,
  higher entries in a column, above its parent. This is the transitive-reduction fact of
  2.5 made concrete: the forest keeps one edge per column (the minimal one), so a column's
  other fills are transitively implied by climbing. The 7 fill entries collapse to 2 fill
  tree-edges.

That last note is a nice payoff, it makes "forest = transitive reduction" (2.5)
concrete: 7 fill entries collapse to 2 fill tree-edges.

Contrast: the *row-major* numbering of the same grid collapses the forest to a single
**path** `1->2->...->9` (`parent = [2,3,4,5,6,7,8,9,NIL]`), height 8, no parallelism,
fill spread across the whole band. Nested dissection turns that path into a balanced tree
of height 6 with two independent halves. Same matrix, same theory, better tree, the
ordering is what the forest sees.

### 2.7 A lighter ordering of the same grid

The 2.6 ordering numbered each subdomain column top to bottom (`1,2,3` down the left),
making each a chain. A better choice numbers the two *ends* of every column before its
*middle*, so the middle vertex separates the two ends, a nested dissection applied once
more, inside each line. Relabel the same grid this way (`3`, `6`, `9` are the per-column
sub-separators):

```
    1 -- 7 -- 4
    |    |    |
    3 -- 9 -- 6
    |    |    |
    2 -- 8 -- 5
```

`A` under the new labels:

```
     1 2 3 4 5 6 7 8 9
  1  X . X . . . X . .
  2  . X X . . . . X .
  3  X X X . . . . . X
  4  . . . X . X X . .
  5  . . . . X X . X .
  6  . . . X X X . . X
  7  X . . X . . X . X
  8  . X . . X . . X X
  9  . . X . . X X X X
```

Compute the forest by the climb of 2.4:

```
k=1: no earlier neighbors
k=2: no earlier neighbors
k=3: earlier neighbor j=1                        -> parent[1]=3
k=3: earlier neighbor j=2                        -> parent[2]=3
k=4: no earlier neighbors
k=5: no earlier neighbors
k=6: earlier neighbor j=4                        -> parent[4]=6
k=6: earlier neighbor j=5                        -> parent[5]=6
k=7: earlier neighbor j=1, climb 1->3, attach    -> parent[3]=7 fill
k=7: earlier neighbor j=4, climb 4->6, attach    -> parent[6]=7 fill
k=8: earlier neighbor j=2, climb 2->3->7, attach -> parent[7]=8 fill
k=8: earlier neighbor j=5, climb 5->6->7->8, reaches k=8, skip
k=9: earlier neighbor j=3, climb 3->7->8, attach -> parent[8]=9
k=9: earlier neighbor j=6, climb 6->7->8->9, reaches k=9, skip
k=9: earlier neighbor j=7, climb 7->8->9, reaches k=9, skip
k=9: earlier neighbor j=8, climb 8->9, reaches k=9, skip
```

The same under path compression:

```
k=1: no earlier neighbors
k=2: no earlier neighbors
k=3: earlier neighbor j=1                        -> parent[1]=3, ancestor[1]=3
k=3: earlier neighbor j=2                        -> parent[2]=3, ancestor[2]=3
k=4: no earlier neighbors
k=5: no earlier neighbors
k=6: earlier neighbor j=4                        -> parent[4]=6, ancestor[4]=6
k=6: earlier neighbor j=5                        -> parent[5]=6, ancestor[5]=6
k=7: earlier neighbor j=1, climb 1->3, attach    -> parent[3]=7 fill, ancestor[1,3]=7
k=7: earlier neighbor j=4, climb 4->6, attach    -> parent[6]=7 fill, ancestor[4,6]=7
k=8: earlier neighbor j=2, climb 2->3->7, attach -> parent[7]=8 fill, ancestor[2,3,7]=8
k=8: earlier neighbor j=5, climb 5->6->7, skip (already at k), ancestor[5,6]=8
k=9: earlier neighbor j=3, climb 3->8, attach    -> parent[8]=9, ancestor[3,8]=9
k=9: earlier neighbor j=6, climb 6->8, skip (already at k), ancestor[6]=9
k=9: earlier neighbor j=7, climb 7->8, skip (already at k), ancestor[7]=9
k=9: earlier neighbor j=8, skip (already at k)
```

Here the shortcut shows at `k=9`: `3->8` directly rather than `3->7->8`, since `k=8`
redirected `ancestor[3]` to `8`.

The elimination game (label order), each subdomain column's two ends fill into its
sub-separator, which then fills into the main separator:

```
eliminate | higher neighbors | new fill edges | total fill
----------+-------------------+----------------+-----------
    1     | {3,7}             | (3,7)          |   1
    2     | {3,8}             | (3,8)          |   2
    3     | {7,8,9}           | (7,8)          |   3
    4     | {6,7}             | (6,7)          |   4
    5     | {6,8}             | (6,8)          |   5
    6     | {7,8,9}           | -              |   5
    7     | {8,9}             | -              |   5
    8     | {9}               | -              |   5
    9     | {}                | -              |   5
```

The forest `parent = [3, 3, 7, 6, 6, 7, 8, 9, NIL]`:

```
          9
          |
          8
          |
          7
        /   \
       3     6
      / \   / \
     1   2 4   5
```

`L` structure (`F` = fill), **5 fills**:

```
     1 2 3 4 5 6 7 8 9
  1  X
  2  . X
  3  X X X
  4  . . . X
  5  . . . . X
  6  . . . X X X
  7  X . F X . F X
  8  . X F . X F F X
  9  . . X . . X X X X
```

Forest edges, original vs fill:

```
forest edge | A[parent][j] | kind
------------+--------------+----------
  1 -> 3    | A[3][1] = X  | original
  2 -> 3    | A[3][2] = X  | original
  3 -> 7    | A[7][3] = .  | fill
  4 -> 6    | A[6][4] = X  | original
  5 -> 6    | A[6][5] = X  | original
  6 -> 7    | A[7][6] = .  | fill
  7 -> 8    | A[8][7] = .  | fill
  8 -> 9    | A[9][8] = X  | original

original: 1-3, 2-3, 4-6, 5-6, 8-9   fill: 3-7, 6-7, 7-8
```

Against 2.6, the same matrix under this ordering gives **5 fills instead of 7**, and a
shorter, bushier forest (height 5 rather than 6). The reason is exactly the
sub-separators: 2.6's left block was the chain `1 -> 2 -> 3`, whereas here `1` and `2` are
independent leaves under `3`, so the column that was a path becomes a two-leaf fork, one
fewer coupling to fill, per subdomain. The sub-separator shows in the elimination box:
eliminating the two ends `1, 2` creates `(3,7), (3,8)`, and their meeting point `3` then
creates `(7,8)`, confining each subdomain's interaction to a single vertex. Numbering the
separators last, recursively, is the whole idea of nested dissection, this subsection is
one more level of it applied to 2.6.

One subtlety is worth drawing out: 2.7 has *fewer* fill entries (5 vs 7) yet *more* fill
forest edges (3 vs 2). These count different things, total nonzeros created in `L`,
versus columns whose parent edge (first below-diagonal nonzero) is fill. Split the fill
into on-tree (a column's minimal entry, so a forest edge) and off-tree (higher entries,
above the parent) and they reconcile:

```
            on-tree fill    off-tree fill    total
2.6              2                5            7
2.7              3                2            5
```

2.6's chains produce mostly off-tree fill, each chain column fills into several separator
rows, most of them above its parent, redundant climbs high in the band. 2.7's
sub-separators cut those columns to near-single below-diagonal entries, concentrating fill
at the parent. So a tighter ordering does two things at once: it lowers total fill, and it
pushes the remaining fill toward the forest's transitive reduction (2.5), a larger share
of it becomes minimal (on-tree), leaving fewer transitively-implied off-tree entries. Fill
forest edges going up while total fill goes down is a signature of that.

### 2.8 Connection to the traversals

The two schedules of 1.3-1.4 are tree motions. Left-looking column `k` gathers from its
**descendants**, exactly the columns whose scatters reach `k` (the subtree below it).
Right-looking column `j` pushes to its **ancestors**, `parent(j)` first. "Descendants
feed a column; the column feeds its ancestors" is the forest restatement of
gather-versus-scatter. And because siblings are independent, any topological order of the
forest is a valid elimination order, the freedom the ordering phase exploits.

From here two roads leave the forest, both walks already visible above: **symbolic
factorization** (Corollary 2's union up the tree, producing `Struct(L)`) and
**supernodes** (the nesting chains, grouped into dense blocks). Both are consumers of
`parent[]`.


## 3. Symbolic factorization

Symbolic factorization is Corollary 2 turned into a loop. Given the forest, it builds
`Struct(L)`, the nonzero pattern of every column of `L`, in one bottom-up pass, unioning
each column's original entries with its children's patterns. No arithmetic, no explicit
fill formation; the output is the storage map the numeric phase allocates and fills.

### 3.1 The algorithm

Corollary 2 is the recurrence:

```
Struct(j) = { i > j : A[i][j] != 0 }
            union, over each child c of j, ( Struct(c) \ {j} )
```

Turn it into code. Process columns in increasing order `1 .. n`; because `parent(c) > c`,
every child is finished before its parent, so the natural column order is already a valid
bottom-up (topological) order of the forest, no separate tree traversal is needed. Each
child's pattern is enumerated through the forest's child links.

```
symbolicFactor(A, forest) -> Struct:
    for j = 1 .. n: # increasing order = children before parents
        Struct(j) = { i > j : A[i][j] != 0 } # start from column j of A
        for each child c of j: # absorb each child's pattern
            Struct(j) = Struct(j) union ( Struct(c) \ {j} )
```

**Cost, union-with-marks.** Done naively the set unions would be expensive, but a mark
array (a timestamp per row) makes each union a linear scan with `O(1)` duplicate
rejection: to build column `j`, append a row and stamp `mark[i] = j`; a row already
stamped `j` is a duplicate and is skipped.

```
symbolicFactor(A, forest) -> Struct: # with marks
    for r = 1 .. n:
        mark[r] = 0
    for j = 1 .. n:
        for each i > j with A[i][j] != 0: # column j of A
            append i to Struct(j); mark[i] = j
        for each child c of j:
            for each i in Struct(c) with i != j: # scan child pattern
                if mark[i] != j:
                    append i to Struct(j); mark[i] = j
```

The whole phase is `O(nnz(L))`, optimal, proportional to the output it produces. The key
accounting: each column `c` is a child of exactly one parent, so its pattern `Struct(c)`
is scanned exactly once (when its parent is built). Summed over all columns that is
`nnz(L)`, plus one pass over `A`. So the total is `O(nnz(A) + nnz(L)) = O(nnz(L))`.

What it touches: `A`'s structure (each column once) and the forest. It never forms fill
explicitly, never does arithmetic, never looks at a numeric value. This is the claim of
2.2, `Struct(L)` is a function of `Struct(A)` alone, made operational: `Struct(A)` and
the forest go in, `Struct(L)` comes out.

### 3.2 Worked example: the grid

Take the grid of 2.6, forest `parent = [2, 3, 7, 5, 6, 7, 8, 9, NIL]`. First read each
column's original below-diagonal entries `Acol(j) = { i > j : A[i][j] != 0 }` off `A`, and
the children off `parent[]`:

```
col j:    1     2     3     4     5     6     7     8     9
Acol(j): {2,7} {3,8} {9}   {5,7} {6,8} {9}   {8}   {9}   {}
children: -     1     2     -     4     5     3,6   7     8
```

Now run the recurrence bottom-up. Each line is `Struct(j) = Acol(j)  union  (child terms
Struct(c) \ {j})`; the child term drops `j` itself, since a child's edge to its parent is
not part of the parent's below-diagonal pattern:

```
j=1  children -      Struct(1) = {2,7}                    = {2,7}
j=2  children 1      Struct(2) = {3,8} + {7}              = {3,7,8}
j=3  children 2      Struct(3) = {9}   + {7,8}            = {7,8,9}
j=4  children -      Struct(4) = {5,7}                    = {5,7}
j=5  children 4      Struct(5) = {6,8} + {7}              = {6,7,8}
j=6  children 5      Struct(6) = {9}   + {7,8}            = {7,8,9}
j=7  children 3,6    Struct(7) = {8}   + {8,9} + {8,9}    = {8,9}
j=8  children 7      Struct(8) = {9}   + {9}              = {9}
j=9  children 8      Struct(9) = {}    + {}               = {}
```

(`+` denotes set union; the child terms are `Struct(c) \ {j}`.) These column patterns are
exactly the `X`/`F` positions of `L` in 2.6, column 2 gains row 7 as fill from its child
1, column 3 inherits `{7,8}` from column 2, and at column 7 the two subtrees (children 3
and 6) merge, each contributing `{8,9}`. The recurrence reconstructs the entire fill
pattern from `A` and the forest, with no arithmetic.

Two features of the forest show through the trace. Fill enters exactly at the multi-hop
climbs of 2.6: a column's pattern grows only where a child hands up rows that were not in
`A`'s column. And the two independent subtrees `{1,2,3}` and `{4,5,6}` are built without
reference to each other, column 7 is the first to see both, which is the sibling
independence of 2.8: the symbolic pass over disjoint subtrees is itself independent work.

### 3.3 Output and storage

The product is `Struct(1..n)`, stored **flat** as compressed-sparse-column (CSC-of-L): one
row-index array holding the columns' patterns back to back, and a column-pointer array
marking where each column begins. That flat pattern is the allocation the numeric phase
fills, it fixes where every nonzero of `L` lives before a single floating-point operation
runs, which is the whole purpose of computing it ahead of time.

One structural fact carries into the next section. By the containment theorem a child's
pattern is nearly its parent's (`Struct(c) \ {parent(c)}` is a subset of
`Struct(parent(c))`), so along a chain of single-child columns the patterns *nest*, each
is the next with one fewer row. Storing and processing such a chain column by column
repeats almost identical patterns. Grouping the chain into one unit with a shared pattern,
a **supernode**, removes that redundancy, and is the subject of Section 4.


## 4. Supernodes

Section 3 ended on a redundancy. Along a chain of single-child columns the patterns nest,
each is the next with one fewer row, so storing and processing such a chain column by
column repeats nearly the same pattern over and over. A **supernode** is that chain,
grouped into one unit with a single shared pattern. This section says exactly which
columns may be grouped, how to find the grouping in one pass, and what it buys.

### 4.1 Definition

A **fundamental supernode** is a maximal set of columns `{j, j+1's successor, ..., k}`
forming a path in the elimination forest, such that

1. every column on the path except the bottom one has **exactly one child**, and
2. the columns have **identical sparsity patterns**, in the sense that the pattern of each
   column is the pattern of the next with that next column's index removed.

Condition 2 is the point: the columns share one pattern, so the supernode needs one index
list rather than one per column. Condition 1 is what makes them *contiguous* in the forest:
if a column on the path had two children, the two subtrees would both feed it, and the
column below it on the path would not carry the whole story.

Fundamental supernodes are **unique** for a given forest. There is no choice to make and
no heuristic to tune, unlike the threshold-based merging that trades explicit zeros for
larger blocks. That uniqueness is why they are the natural default.

Write `snode(j)` for the supernode containing column `j`. The supernode's columns are its
**front indices**; the rows below them, common to all of them, are its **update indices**.

### 4.2 The counting test

Condition 2 looks like it needs a set comparison per candidate pair. It does not. Let `j`
be the only child of `k`, so `j < k` and `parent(j) = k`. If they belong to one supernode,
then by the containment theorem of 2.3 the rows below `j` are exactly `k` followed by the
rows below `k`:

```
Struct(j) = {k} union Struct(k)
```

Take sizes of both sides. The union is disjoint (`k` is not in `Struct(k)`), so

```
|Struct(j)| = 1 + |Struct(k)|
```

and the test becomes arithmetic on numbers we already have from symbolic factorization, or
from the column counts computed with the forest:

```
k joins j's supernode  <=>  j is k's only child  and  |Struct(k)| = |Struct(j)| - 1
```

No pattern is ever compared. This is the whole trick: containment guarantees `Struct(j)`
is *at least* `{k} union Struct(k)`, so equality of the sizes forces equality of the sets.

**The supernodal form.** Written for a forest whose vertices are already supernodes rather
than single columns, the same identity reads

```
|update(j)| = |front(k)| + |update(k)|
```

because the rows below supernode `j` are the columns of `k` followed by the rows below `k`.
For a nodal forest `|front(k)| = 1` and this collapses to the test above. The general form
is worth keeping: it is what lets compression run on a forest that has already been
compressed. Note it is the **parent's** front size that appears, not the child's.

### 4.3 Compression

One pass in increasing column order suffices. Because `parent(j) > j`, a column's child is
always numbered before it, so when column `k` is reached its child's supernode is already
decided.

```
compress(forest, Struct) -> snode: # fundamental supernodes
    s = 0
    for k = 1 .. n: # increasing order = children before parents
        j = firstChild[k]
        if j != NIL and j == lastChild[k] # condition 1: j is k's only child
               and |Struct(k)| == |Struct(j)| - 1: # condition 2: one shared pattern
            snode[k] = snode[j] # k continues j's supernode
        else:
            s = s + 1
            snode[k] = s # k starts a new supernode
```

Both conditions of 4.1 are in that test. `j == lastChild[k]` says the child list of `k` has
a single entry, which is why the forest keeps a last-child link and not only a first-child
link. The size equality is the pattern test of 4.2.

That assigns every column to a supernode and counts them. The forest must then be rebuilt
over supernodes rather than columns. A second pass, in *decreasing* order, does it:

```
rebuild(forest, snode) -> supernodal forest:
    for S = 1 .. s:
        front(S) = 0
        done[S] = false
    for k = n .. 1: # decreasing order, so a supernode's topmost column comes first
        S = snode[k]
        front(S) = front(S) + 1 # every column of S adds to its front size
        if not done[S]: # k is the topmost column of S
            if parent(k) != NIL:
                supParent(S) = snode[parent(k)]
            else:
                supParent(S) = NIL # S is a root
            update(S) = |Struct(k)| # the rows below k are the rows below S
            done[S] = true
```

The asymmetry between the two lines inside the `if` is worth reading twice. Front sizes
**accumulate** over every column of the supernode, so that line sits outside the guard.
The parent link and the update size are taken **once**, from the first column seen, which
in decreasing order is the supernode's **topmost** column. That is the right column for
both: its parent is the one link that leaves the supernode, and its rows below are exactly
the supernode's update indices, the ones common to all of its columns.

Finally the child and sibling links are rebuilt from the new parent links exactly as in
2.4, and the height recomputed. Cost is `O(n)` beyond what we already have, so compression
is free next to the factorization it accelerates.

**Columns of a supernode need not be consecutive.** The forest is in topological order, not
necessarily postorder, so a supernode can own scattered column indices. Consumers must
gather its front indices through `snode[]` rather than assume a contiguous range. (A
postordering of the forest would make them consecutive, and is worth doing for other
reasons, but the algorithms above do not depend on it.)

### 4.4 Worked example: the grid

Take the grid of 2.6 again, with `parent = [2, 3, 7, 5, 6, 7, 8, 9, NIL]` and the patterns
computed in 3.2. Tabulate what the test needs, the only child (if any) and the pattern
sizes:

```
col k:        1     2     3     4     5     6     7     8     9
Struct(k):   {2,7} {3,7,8} {7,8,9} {5,7} {6,7,8} {7,8,9} {8,9} {9} {}
|Struct(k)|:  2     3     3     2     3     3     2     1     0
children:     -     1     2     -     4     5     3,6   7     8
only child:   -     1     2     -     4     5     -     7     8
```

Now run the test on each column with an only child:

```
k=2  j=1   |Struct(2)|=3, |Struct(1)|-1=1   3 != 1   no merge
k=3  j=2   |Struct(3)|=3, |Struct(2)|-1=2   3 != 2   no merge
k=5  j=4   |Struct(5)|=3, |Struct(4)|-1=1   3 != 1   no merge
k=6  j=5   |Struct(6)|=3, |Struct(5)|-1=2   3 != 2   no merge
k=7   -    two children (3 and 6)                    no merge
k=8  j=7   |Struct(8)|=1, |Struct(7)|-1=1   1 == 1   merge, 8 joins 7
k=9  j=8   |Struct(9)|=0, |Struct(8)|-1=0   0 == 0   merge, 9 joins 8
```

Nine columns collapse to **seven supernodes**: six trivial ones, and `{7, 8, 9}`.

```
snode:     1   2   3   4   5   6   7
columns:  {1} {2} {3} {4} {5} {6} {7,8,9}
front:     1   1   1   1   1   1   3
update:    2   3   3   2   3   3   0
```

The one non-trivial supernode is exactly the **separator**. That is not a coincidence, it
is nested dissection working as designed: the separator is eliminated last, by then it is
completely filled in, and a completely filled-in block is precisely a set of columns with
one shared pattern. Its index set is `{7,8,9}` with no update rows, a dense 3x3 trailing
block.

The left and right blocks do not compress at all, and the trace shows why. Column 2's
pattern is `{3,7,8}` while its child column 1 has `{2,7}`: dropping the 2 leaves `{7}`,
not `{7,8}`. The patterns differ, so the columns cannot share an index list. Column 2
*acquires* row 8 from `A` that column 1 never had. Only once a column's pattern has stopped
growing do the columns start agreeing, and in this ordering that happens only inside the
separator.

### 4.5 What changes downstream

Symbolic factorization runs unchanged, at supernode granularity. The recurrence of 3.1 is
the same union, with two adjustments. It reads the patterns of *every* front column of the
supernode rather than a single column, and when absorbing a child it drops the child's
front indices rather than a single index `j`:

```
symbolicFactor(A, supernodal forest) -> Idx:
    for K = 1 .. numSupernodes: # increasing = children before parents
        Idx(K) = union, over front columns k of K, { i >= k : A[i][k] != 0 }
        for each child J of K: # absorb each child's update indices
            Idx(K) = Idx(K) union ( Idx(J) \ front(J) )
```

For a nodal forest each supernode has one front column and `front(J) = {j}`, and this is
exactly the algorithm of 3.1. The same code covers both regimes, which is a good reason to
write the supernodal form even when the forest happens to be nodal.

The payoff is not in the symbolic phase, which was already optimal, but in the numeric one.
A supernode's columns share a pattern, so its nonzeros form a **dense rectangular block**,
and the numeric factorization can update it with dense matrix kernels (BLAS level 3)
instead of a scatter over individual columns. That is the entire motivation: supernodes
convert sparse column operations into dense block operations, which run at a large multiple
of the speed on real hardware. The larger the supernodes, the better the ratio, which is
also why threshold-based merging exists, it accepts some explicitly stored zeros to buy
bigger blocks. Fundamental supernodes are the free case: they introduce no zeros at all.


## References

The material above is standard sparse-matrix theory; the grouping below points to the
primary sources for each, roughly section by section. Citation details are given from
knowledge and are worth a spot-check against the originals.

**Cholesky; left- and right-looking (Section 1).**

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins
  University Press, 2013. Dense Cholesky, the left-/right-looking schedules, and the
  block (`dsyrk` / `dgemm`) kernels of 1.4.

**Fill and vertex elimination (Sections 1.6, 2).**

- D. J. Rose, "A graph-theoretic study of the numerical solution of sparse positive
  definite systems of linear equations", in *Graph Theory and Computing*, R. C. Read,
  ed., Academic Press, 1972, pp. 183-217. The elimination graph, chordal fill, and the
  clique view of a pivot's neighbors.
- D. J. Rose, R. E. Tarjan, and G. S. Lueker, "Algorithmic aspects of vertex elimination
  on graphs", *SIAM J. Comput.* 5(2):266-283, 1976. The fill-path theorem (proved in 1.8
  via the elimination-graph reachability lemma) and the graph theory of elimination.

**The elimination tree / forest (Section 2).**

- R. Schreiber, "A new implementation of sparse Gaussian elimination", *ACM Trans. Math.
  Software* 8(3):256-276, 1982. Formalizes the elimination tree.
- J. W. H. Liu, "A compact row storage scheme for Cholesky factors using elimination
  trees", *ACM Trans. Math. Software* 12(2):127-148, 1986. The path-compressed etree
  computation (Liu's algorithm in 2.4) and compact factor storage.
- J. W. H. Liu, "The role of elimination trees in sparse factorization", *SIAM J. Matrix
  Anal. Appl.* 11(1):134-172, 1990. The survey: the containment property, the symbolic
  union recurrence, and subtree independence / parallelism.
- R. E. Tarjan, "Efficiency of a good but not linear set union algorithm", *J. ACM*
  22(2):215-225, 1975. The union-find bound behind the `alpha(n)` in Liu's algorithm.

**Symbolic factorization, ordering, sparse direct methods (Section 3; general).**

- A. George and J. W. H. Liu, *Computer Solution of Large Sparse Positive Definite
  Systems*, Prentice-Hall, 1981. Symbolic factorization, ordering, and the classic
  treatment of the whole pipeline.
- A. George, "Nested dissection of a regular finite element mesh", *SIAM J. Numer. Anal.*
  10(2):345-363, 1973. Nested dissection, the ordering of the grid example in 2.6.
- T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006. A modern,
  code-level reference for the entire order / symbolic / numeric pipeline. Its CSparse
  `cs_etree` computes the tree from the upper triangle of each column, the standard route
  used by Oblio (2.4); single-triangle storage in the other orientation is transposed
  first.

**On the framing.** Two presentational choices here are expository, not lifted from a
single source: casting the forest as the *transitive reduction of the update DAG* (2.5)
is a folklore-standard restatement rather than any one paper's language, and structuring
the 2.4 correctness proof around the "Cholesky factors nest" lemma (Lemma A) is this
document's packaging of the incremental argument. Both are consistent with the sources
above; the specific wording is ours.
