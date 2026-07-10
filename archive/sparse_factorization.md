# Sparse Factorization

Notes on the algorithms behind Oblio's numerical phase. The goal is to keep the
mathematics and the traversal schedules in one place, so the code (left-looking,
right-looking, multifrontal) can be read against a common reference.

Conventions: matrices are `n x n`, indexed `1 .. n` (math style). `A` is the input
matrix, `L` the computed factor. Column indices are `j` and `k`; row indices are `i`.
Throughout, **`k` is the column currently being produced**, and **`j` is the other
column** — to the *left* (`j < k`, an earlier column being gathered) in the
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
lower (or upper) triangle is needed — the two triangles carry the same information.

### 1.2 Scalar recurrences

Multiplying `L L^T` and matching entries gives, for the current column `k` and rows
`i > k`:

```
L[k][k] = sqrt( A[k][k] - sum_{j<k} L[k][j]^2 )

L[i][k] = ( A[i][k] - sum_{j<k} L[i][j] * L[k][j] ) / L[k][k]
```

Each entry of `L` depends only on columns `j < k` — everything to its left. Cholesky
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
only when *both* `L[i][j]` and `L[k][j]` are nonzero — i.e. when rows `i` and `k` both
have structure in column `j`. So the left-looking update **gathers, at the moment
column `k` is finalized, all the contributions that earlier pivots produced for
`(i,k)`.**

Organized as the loop nest above (`k` then `j` then `i`), each earlier column `j`
contributes to the *whole* of column `k` in one inner sweep — a scaled column
subtraction `col_k <- col_k - L[k][j] * col_j` (an axpy). This is the "column
modification" form, and it is the shape the supernodal left-looking traversal uses:
each contributing column (or supernode) updates the target column (or supernode) as a
vector/block operation, rather than one entry at a time.

### 1.4 Right-looking Cholesky (outer product / "scatter")

The moment column `j` is finished, it *pushes* its contribution into every trailing
column at once. A trailing entry is written many times — once per relevant earlier
pivot — accumulating updates as pivots complete. Indexing: `j` is the column just
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
then run a single column loop — `A` is read-only, the factor is built in `L`. They
differ only in what that loop does: left-looking **gathers then finalizes** (pull the
earlier columns' contributions into column `k`, then finalize it), right-looking
**finalizes then scatters** (finalize column `j`, then push its contribution into the
later columns). Gather versus scatter is the whole distinction; the seeding is
deliberately identical so the two read as mirror images.

The double loop applies the **rank-1 outer product** of column `j` with itself to the
later columns of `L` (lower triangle only — the factor stays symmetric in structure, so
the upper half is redundant; the `i == k` case updates the diagonal
`L[k][k] -= L[k][j]^2`; this is the shape BLAS `dsyrk` computes):

```
L_trailing <- L_trailing - L[*][j] * L[*][j]^T
```

Entry `(i,k)` of that outer product is `L[i][j] * L[k][j]`, nonzero only when column
`j` has nonzeros in *both* rows `i` and `k`. So finishing column `j` **scatters a
contribution into every later `(i,k)` whose rows both touch column `j`** — updating the
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
later `k`. Same term, same target — one sums at read time, the other deposits at write
time. Pull-when-needed versus push-when-ready.

### 1.6 Connection to vertex elimination and fill

View the nonzero structure of `A` as a graph `G`: one vertex per row/column, and an
edge joining the two indices of each off-diagonal nonzero. Factorization eliminates
vertices in order `1, 2, .., n`.

Across this section the three indices keep fixed roles, matching the `j < k` convention
of the algorithms: **`j` is the vertex being eliminated** (the pivot), and **`i` and `k`
are two later neighbours of `j`** that its elimination connects — so the fill it creates
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
affected column is reached — the gather at column `k` replays, all at once, the scatters
that pivots `j < k` would have made. Same eliminations, same fill; scatter does them
eagerly, gather does them lazily.

**Fill = the outer product connecting a pivot's neighbors.** Entry `(i,k)` of
`a_j a_j^T` is nonzero iff `a_j[i]` and `a_j[k]` are both nonzero — iff `i` and `k` are
*both neighbors of `j`*. If `B[i][k]` was zero (i and k were not connected), it becomes
nonzero: a **fill** edge. So eliminating `j` makes all of `j`'s still-uneliminated
neighbors mutually adjacent — a clique — and each new clique edge is a fill entry.

**The dot-product term is that same fill, observed later.** In the left-looking view,
the term `L[i][j] * L[k][j]` is nonzero iff both `i` and `k` descend from pivot `j`'s
elimination — i.e. iff `j` was a shared neighbor. So the three statements below are the
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
cancellation** — if the structure admits a nonzero contribution, the position is
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

#### Example 1 — Arrowhead (ordering determines fill)

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
pair* of them — turning the five leaves into a complete graph `K5`. The trailing 5x5
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

That is **10 fill entries** — the maximum possible for this pattern. Every off-diagonal
among rows/cols 2-6 is fill.

Now reorder so the hub is eliminated **last** (elimination order 2,3,4,5,6,1; the hub
becomes the highest-numbered vertex). Each leaf now has only the hub as a later
neighbor, so no elimination ever has two later-neighbors to connect — **no clique
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
ordering phase (AMD / MMD) exists — for the arrowhead, eliminating the hub last is
optimal, and a good ordering finds that automatically.

#### Example 2 — Path with long edges (fill cascades)

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
carrying the accumulated fill forward — the etree union rule of the symbolic phase.
