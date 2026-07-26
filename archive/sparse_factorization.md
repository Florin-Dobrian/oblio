# Sparse Factorization

Notes on the algorithms behind Oblio: the ordering that decides the fill, the elimination
forest and symbolic factorization that predict it, the supernodes that make it fast, and the
numeric traversals that compute it. The goal is to keep the mathematics and the schedules in one
place, so the code (left-looking, right-looking, multifrontal) can be read against a common
reference.

The sections run in dependency order rather than pipeline order: Cholesky and fill first, since
everything is stated in their terms; then the forest, symbolic factorization, and supernodes; and
ordering last (Section 5), because it is the one phase whose *output* the others take as given.

Conventions: matrices are `n x n`, indexed `1 .. n` (math style). `A` is the input
matrix, `L` the computed factor. Column indices are `j` and `k`; row indices are `i`.
Throughout, **`k` is the column currently being produced**, and **`j` is the other
column**, to the *left* (`j < k`, an earlier column being gathered) in the
left-looking view, or to the *right* (`j > k`, a trailing column being updated) in the
right-looking view. Sums written `sum_{j<k}` run over `j = 1 .. k-1`.


## 0. The direct solution of Ax = b

Solving `A x = b` directly is three operations: forming a product `A x`, factoring `A`
into triangular pieces, and solving triangular systems. All three come in the same two
shapes, a **gather** that pulls finished values in and writes each output once, and a
**scatter** that pushes each finished value out into the rest. The factorization
schedules of Section 1 (left- and right-looking) are these same two shapes on a harder
operation, so it helps to meet them first on the easy ones. Conventions are the front
matter's: `A` is `n x n`, indexed `1 .. n`; `i` is a row, `j` and `k` are columns.

### 0.1 Matrix-vector product

`y = A x`. Two loop orders, differing only in which index is outer.

**Row-oriented (dot product / "gather").** Each `y[i]` is the dot product of row `i` of
`A` with `x`, computed and written once.

```
matvec_row(A, x) -> y:
    for i = 1 .. n:          # row i
        y[i] = 0
        for j = 1 .. n:      # column j
            y[i] = y[i] + A[i][j] * x[j]
```

**Column-oriented (axpy / "scatter").** Each column `j`, scaled by `x[j]`, is added into
the whole of `y`; `y` starts at zero and accumulates across columns.

```
matvec_col(A, x) -> y:
    for i = 1 .. n:
        y[i] = 0
    for j = 1 .. n:          # column j
        for i = 1 .. n:      # row i
            y[i] = y[i] + A[i][j] * x[j]
```

Row order reads `x` in full for each output and writes `y[i]` once (gather); column
order reads `x[j]` once and updates every `y[i]` (scatter). Same arithmetic, transposed
schedule, the distinction that runs through everything below.

### 0.2 LU factorization

For a general (nonsymmetric) `A`, LU computes a unit lower-triangular `L` (`L[i][i] = 1`)
and an upper-triangular `U` with

```
A = L U
```

Matching entries of `L U` against `A`, with `k` the current column, gives the column-`k`
entries of `U` (on and above the diagonal) and of `L` (below it):

```
U[i][k] = A[i][k] - sum_{j<i} L[i][j] * U[j][k]                 for i = 1 .. k
L[i][k] = ( A[i][k] - sum_{j<k} L[i][j] * U[j][k] ) / U[k][k]   for i = k+1 .. n
```

`U[k][k]` is the pivot and must be nonzero (pivoting, which reorders rows to keep it away
from zero, is omitted here). Cholesky is the symmetric special case: when `A` is SPD the
two triangles carry the same information and only one is computed, `U = L^T` up to the
diagonal (Section 1).

As with Cholesky the two schedules compute the same `L` and `U` and differ only in *when*
the cross-column contributions are applied. Both overwrite a working copy `LU` of `A`,
leaving `U` in its upper triangle and `L` (minus the unit diagonal) in its strict lower
triangle.

**Left-looking (gather).** Column `k` waits, then pulls in every earlier column's
contribution in one pass before its `L` part is normalized.

```
lu_left(A) -> (L, U):
    LU = A
    for k = 1 .. n:              # current column k
        for j = 1 .. k-1:        # earlier column j
            for i = j+1 .. n:    # LU[i][k] -= L[i][j] * U[j][k]
                LU[i][k] = LU[i][k] - LU[i][j] * LU[j][k]
        for i = k+1 .. n:        # normalize the L part by the pivot
            LU[i][k] = LU[i][k] / LU[k][k]
```

**Right-looking (scatter).** The moment column `k` is finished, its multipliers push a
rank-1 update into the entire trailing submatrix.

```
lu_right(A) -> (L, U):
    LU = A
    for k = 1 .. n:              # pivot column k
        for i = k+1 .. n:        # multipliers: column k of L
            LU[i][k] = LU[i][k] / LU[k][k]
        for i = k+1 .. n:        # rank-1 update of the trailing block
            for j = k+1 .. n:
                LU[i][j] = LU[i][j] - LU[i][k] * LU[k][j]
```

Left-looking gathers `sum_{j<k} L[i][j] U[j][k]` into column `k` when it is reached;
right-looking scatters `LU[i][k] * LU[k][j]` into every trailing `(i,j)` as soon as column
`k` is done. That rank-1 update is the outer product of column `k`'s multipliers with row
`k` of `U`, the unsymmetric analog of the Cholesky rank-1 update (a general `dgemm` in
place of `dsyrk`).

### 0.3 Triangular solves for LU factorization

With `A = L U`, solving `A x = b` is two triangular solves: forward `L y = b`, then back
`U x = y`. Each comes in the same two loop orders as the matvec, a row-oriented dot
product (gather) and a column-oriented axpy (scatter).

**Forward substitution, `L y = b`** (`L` unit lower-triangular, so no diagonal division).

Row-oriented (gather), `y[i]` finished by subtracting the dot product of the already-known
`y[1 .. i-1]`:

```
forward_row(L, b) -> y:
    for i = 1 .. n:
        y[i] = b[i]
        for j = 1 .. i-1:
            y[i] = y[i] - L[i][j] * y[j]
```

Column-oriented (scatter), each `y[j]` pushed into the later equations as soon as it is known:

```
forward_col(L, b) -> y:
    for i = 1 .. n:
        y[i] = b[i]
    for j = 1 .. n:
        for i = j+1 .. n:
            y[i] = y[i] - L[i][j] * y[j]
```

**Back substitution, `U x = y`** (`U` upper-triangular; divide by the diagonal pivot).

Row-oriented (gather):

```
back_row(U, y) -> x:
    for i = n .. 1:
        x[i] = y[i]
        for j = i+1 .. n:
            x[i] = x[i] - U[i][j] * x[j]
        x[i] = x[i] / U[i][i]
```

Column-oriented (scatter):

```
back_col(U, y) -> x:
    for i = 1 .. n:
        x[i] = y[i]
    for j = n .. 1:
        x[j] = x[j] / U[j][j]
        for i = 1 .. j-1:
            x[i] = x[i] - U[i][j] * x[j]
```

Forward and back are mirror images across the diagonal: row order gathers the solved
unknowns into the current one, column order scatters each unknown, once solved, into the
rest. The gather/scatter pair that opened with the matvec closes the direct solve.

### 0.4 LDL^T factorization

For a symmetric `A`, LDL^T computes a unit lower-triangular `L` (`L[i][i] = 1`) and a
diagonal `D` with

```
A = L D L^T
```

It is the LU factorization of 0.2 specialized to a symmetric `A`, with the pivots pulled
out of the upper factor: `U = D L^T`, so `A = L U = L (D L^T) = L D L^T`, and only `L` and
`D` are stored. Equivalently it is Cholesky without the square root, the pivots collected
in `D` rather than taken under a root, so the factorization stays in the arithmetic of `A`
and extends to the symmetric indefinite case (with pivoting, omitted here; the Cholesky
relation `C = L D^{1/2}` is in 0.6). Only the lower triangle is computed; with `k` the
current column,

```
D[k]    = A[k][k] - sum_{j<k} L[k][j]^2 * D[j]
L[i][k] = ( A[i][k] - sum_{j<k} L[i][j] * D[j] * L[k][j] ) / D[k]   for i = k+1 .. n
```

Both schedules leave `D` on the diagonal of the working matrix `LD` and `L` (minus the
unit diagonal) in its strict lower triangle.

**Left-looking (gather).** Column `k` pulls in every earlier column's contribution, each
scaled by that column's pivot `D[j] = LD[j][j]`, before its `L` part is normalized.

```
ldlt_left(A) -> (L, D):
    LD[i][j] = A[i][j] for i >= j, else 0
    for k = 1 .. n:              # current column k
        for j = 1 .. k-1:        # earlier column j
            for i = k .. n:      # LD[i][k] -= L[i][j] * D[j] * L[k][j]
                LD[i][k] = LD[i][k] - LD[i][j] * LD[j][j] * LD[k][j]
        for i = k+1 .. n:        # normalize the L part by the pivot D[k] = LD[k][k]
            LD[i][k] = LD[i][k] / LD[k][k]
```

**Right-looking (scatter).** Once column `k`'s pivot and multipliers are set, they scatter
a `D[k]`-scaled rank-1 update into the trailing lower triangle.

```
ldlt_right(A) -> (L, D):
    LD[i][j] = A[i][j] for i >= j, else 0
    for k = 1 .. n:              # pivot column k
        for i = k+1 .. n:        # multipliers: column k of L
            LD[i][k] = LD[i][k] / LD[k][k]
        for j = k+1 .. n:        # rank-1 update of the trailing lower triangle
            for i = j .. n:
                LD[i][j] = LD[i][j] - LD[i][k] * LD[k][k] * LD[j][k]
```

Left-looking gathers `sum_{j<k} L[i][j] D[j] L[k][j]` into column `k`; right-looking
scatters `LD[i][k] * D[k] * LD[j][k]` into the trailing `(i,j)`. The update is the symmetric
rank-1 form of the LU outer product, touching only the lower triangle.

### 0.5 Triangular solves with LDL^T factorization

With `A = L D L^T`, solving `A x = b` is three steps: forward `L z = b`, a diagonal solve
`D w = z`, then back `L^T x = w`. The two triangular steps use `L` (unit lower) and its
transpose (unit upper), each in the two loop orders of 0.3; the diagonal step between them
is a single loop.

**Step 1, forward `L z = b`** (unit lower, no diagonal division).

```
forward_row(L, b) -> z:
    for i = 1 .. n:
        z[i] = b[i]
        for j = 1 .. i-1:
            z[i] = z[i] - L[i][j] * z[j]
```

```
forward_col(L, b) -> z:
    for i = 1 .. n:
        z[i] = b[i]
    for j = 1 .. n:
        for i = j+1 .. n:
            z[i] = z[i] - L[i][j] * z[j]
```

**Step 2, diagonal `D w = z`** (a single loop).

```
diag(D, z) -> w:
    for k = 1 .. n:
        w[k] = z[k] / D[k]
```

**Step 3, back `L^T x = w`** (unit upper, `L^T[i][j] = L[j][i]`, no diagonal division).

```
backT_row(L, w) -> x:
    for i = n .. 1:
        x[i] = w[i]
        for j = i+1 .. n:
            x[i] = x[i] - L[j][i] * x[j]
```

```
backT_col(L, w) -> x:
    for i = 1 .. n:
        x[i] = w[i]
    for j = n .. 1:
        for i = 1 .. j-1:
            x[i] = x[i] - L[j][i] * x[j]
```

The forward is the `L` solve of 0.3; the back reads `L`'s columns as the rows of `L^T`,
the same entries accessed the other way. The diagonal solve between them is the one step
with a single loop.

### 0.6 Cholesky factorization

For a symmetric positive definite `A`, Cholesky computes a lower-triangular `C` with

```
A = C C^T
```

It is the LDL^T factorization of 0.4 with the diagonal absorbed into the factor:
`C = L D^{1/2}`, that is `C[i][k] = L[i][k] * sqrt(D[k])` and `C[k][k] = sqrt(D[k])`, so
`A = C C^T = (L D^{1/2})(L D^{1/2})^T = L D L^T`. The square root needs `D > 0`, which is
why Cholesky is the SPD case and LDL^T the general symmetric one. Unlike `L`, the factor
`C` is not unit lower-triangular; its diagonal carries the pivots as `sqrt(D[k])`. With
`k` the current column,

```
C[k][k] = sqrt( A[k][k] - sum_{j<k} C[k][j]^2 )
C[i][k] = ( A[i][k] - sum_{j<k} C[i][j] * C[k][j] ) / C[k][k]   for i = k+1 .. n
```

**Left-looking (gather).**

```
chol_left(A) -> C:
    C[i][j] = A[i][j] for i >= j, else 0
    for k = 1 .. n:              # current column k
        for j = 1 .. k-1:        # earlier column j
            for i = k .. n:      # C[i][k] -= C[i][j] * C[k][j]
                C[i][k] = C[i][k] - C[i][j] * C[k][j]
        C[k][k] = sqrt(C[k][k])
        for i = k+1 .. n:        # normalize by the pivot C[k][k]
            C[i][k] = C[i][k] / C[k][k]
```

**Right-looking (scatter).**

```
chol_right(A) -> C:
    C[i][j] = A[i][j] for i >= j, else 0
    for k = 1 .. n:              # pivot column k
        C[k][k] = sqrt(C[k][k])
        for i = k+1 .. n:        # column k of C
            C[i][k] = C[i][k] / C[k][k]
        for j = k+1 .. n:        # rank-1 update of the trailing lower triangle
            for i = j .. n:
                C[i][j] = C[i][j] - C[i][k] * C[j][k]
```

Left-looking gathers `sum_{j<k} C[i][j] C[k][j]` into column `k`; right-looking scatters
`C[i][k] * C[j][k]` into the trailing `(i,j)`, the same schedules as LU and LDL^T with a
unit weight in place of `U[k][j]` or `D[k]`. Section 1 develops this same factorization in
`L L^T` notation.

### 0.7 Triangular solves with Cholesky factorization

With `A = C C^T`, solving `A x = b` is two steps: forward `C y = b`, then back
`C^T x = y`. There is no separate diagonal step, the `sqrt(D)` that was the LDL^T diagonal
solve of 0.5 is folded into `C`'s diagonal, so each triangular solve divides by
`C[i][i] = sqrt(D[i])` (the unit-diagonal `L` solves of 0.5 did not divide).

**Step 1, forward `C y = b`** (lower triangular, divide by the diagonal).

```
forward_row(C, b) -> y:
    for i = 1 .. n:
        y[i] = b[i]
        for j = 1 .. i-1:
            y[i] = y[i] - C[i][j] * y[j]
        y[i] = y[i] / C[i][i]
```

```
forward_col(C, b) -> y:
    for i = 1 .. n:
        y[i] = b[i]
    for j = 1 .. n:
        y[j] = y[j] / C[j][j]
        for i = j+1 .. n:
            y[i] = y[i] - C[i][j] * y[j]
```

**Step 2, back `C^T x = y`** (upper triangular, `C^T[i][j] = C[j][i]`, divide by the
diagonal).

```
backT_row(C, y) -> x:
    for i = n .. 1:
        x[i] = y[i]
        for j = i+1 .. n:
            x[i] = x[i] - C[j][i] * x[j]
        x[i] = x[i] / C[i][i]
```

```
backT_col(C, y) -> x:
    for i = 1 .. n:
        x[i] = y[i]
    for j = n .. 1:
        x[j] = x[j] / C[j][j]
        for i = 1 .. j-1:
            x[i] = x[i] - C[j][i] * x[j]
```

As with LDL^T the back solve reads `C`'s columns as the rows of `C^T`. The two Cholesky
solves are the LDL^T three-step solve with its diagonal split evenly between them: each
divides by `sqrt(D[i])`, and the two together account for the single `D[i]` division.

### 0.8 Block Cholesky factorization

The scalar recurrences of 0.6 are the `1 x 1` pivot case of a block factorization. Partition
symmetric positive-definite `A` into a leading pivot block and a trailing block:

```
A = [ A11  A12 ]
    [ A21  A22 ]
```

with `A11` square and, by symmetry, `A21 = A12^H`. Cholesky produces a block lower-triangular
`C` with `A = C C^H`:

```
C = [ C11   0  ]
    [ C21  C22 ]
```

Multiplying `C C^H` back out,

```
C C^H = [ C11        0  ] [ C11^H  C21^H ]
        [ C21       C22 ] [   0    C22^H ]

      = [ C11 C11^H              C11 C21^H         ]
        [ C21 C11^H       C21 C21^H + C22 C22^H    ]
```

and matching blocks against `A` gives three equations,

```
A11 = C11 C11^H
A21 = C21 C11^H
A22 = C21 C21^H + C22 C22^H
```

which read top to bottom as the factorization:

```
C11 = chol(A11)                  # LAPACK potrf (dpotrf / zpotrf)
C21 = A21 C11^-H                 # BLAS trsm (dtrsm / ztrsm)
C22 = chol( A22 - C21 C21^H )    # Schur update: BLAS herk (zherk, dsyrk in real) + gemm;
                                 #   then LAPACK potrf on the result
```

Factor the pivot block, solve the off-diagonal block against `C11^H`, then form the Schur
complement `A22 - C21 C21^H` and factor it in turn. The scalar `chol_left`/`chol_right` of 0.6
is this with `A11` a single entry: `C11 = sqrt(A11)`, `C21 = A21 / C11`, and the Schur
complement the rank-1 trailing update. The block form is the same recurrence with the pivot a
block and the reciprocal a triangular solve, which is what lets a supernode factor its whole
front at once through BLAS rather than one column at a time.

### 0.9 Block LDL^T factorization

The complex-symmetric LDL of 0.4 has the same block form, with the pivots kept in a
block-diagonal `D` instead of absorbed into the factor, and with `^T` throughout because this
variant does not conjugate. Partition symmetric `A` as before, `A21 = A12^T`, and factor
`A = L D L^T` with `L` unit block-lower-triangular and `D` block-diagonal:

```
L = [ L11   0  ]      D = [ D11   0  ]
    [ L21  L22 ]          [  0   D22 ]
```

Multiplying `L D L^T` out,

```
L D L^T = [ L11 D11 L11^T                  L11 D11 L21^T           ]
          [ L21 D11 L11^T          L21 D11 L21^T + L22 D22 L22^T    ]
```

and matching against `A`,

```
A11 = L11 D11 L11^T
A21 = L21 D11 L11^T
A22 = L21 D11 L21^T + L22 D22 L22^T
```

which read as:

```
(L11, D11) = ldlt(A11)                     # Oblio's own ldl; no library equivalent (see note)
L21       = A21 L11^-T D11^-1              # BLAS trsm against L11^T, then Oblio scales by D11^-1
(L22, D22) = ldlt( A22 - L21 D11 L21^T )    # Schur update L21 D11 L21^T: Oblio's own formUpper
                                          #   + gemmLower + gemm; then Oblio's own ldl
```

The Oblio-own routines here are not stand-ins for a library call we have not gotten to; the
standard libraries do not provide them. LAPACK's `sytrf`/`hetrf` factor a symmetric block, but
as Bunch-Kaufman *with* pivoting, which the static path refuses (its pivots are fixed by the
symbolic structure) and the dynamic path handles under its own delayed-pivot scheme, so neither
can call `sytrf`. The Schur update `A22 - L21 D11 L21^T` has no BLAS routine at all: the `D` in
the middle rules out a rank-k call (`syrk`/`herk` compute `A A^H`, not `A D A^H`), and there is
no primitive for "`A B` whose product is known symmetric", which is why the update forms
`U := D L21^T` in a scratch (`formUpper`) and then multiplies (`gemmLower`, `gemm`). Cholesky
above, by contrast, is BLAS and LAPACK throughout. These are permanent, not a to-do.

The off-diagonal solve splits in two, a unit-triangular solve against `L11^T` and a division
by the pivot block `D11`, which for a `1 x 1` `D11` is a scalar division and for a `2 x 2` pivot
a `2 x 2` solve. The Schur complement carries the pivot explicitly, `A22 - L21 D11 L21^T`, where
Cholesky folded it into `C21 C21^H`. This is the factorization the static and dynamic LDL
kernels compute, the `2 x 2` `D11` being exactly the delayed-pivot case the dynamic path handles.

### 0.10 Block LDL^H factorization

The Hermitian LDL is the same block factorization with `^H` in place of `^T`, so it conjugates
where LDL^T does not. Partition Hermitian `A`, `A21 = A12^H`, and factor `A = L D L^H` with `L`
unit block-lower-triangular and `D` block-diagonal (its `1 x 1` pivots real, its `2 x 2` pivots
Hermitian):

```
L = [ L11   0  ]      D = [ D11   0  ]
    [ L21  L22 ]          [  0   D22 ]
```

Multiplying `L D L^H` out,

```
L D L^H = [ L11 D11 L11^H                  L11 D11 L21^H           ]
          [ L21 D11 L11^H          L21 D11 L21^H + L22 D22 L22^H    ]
```

and matching against `A`,

```
A11 = L11 D11 L11^H
A21 = L21 D11 L11^H
A22 = L21 D11 L21^H + L22 D22 L22^H
```

which read as:

```
(L11, D11) = ldlh(A11)                     # Oblio's own ldl; no BLAS/LAPACK equivalent (see 0.9)
L21       = A21 L11^-H D11^-1              # BLAS trsm against L11^H, then Oblio scales by D11^-1
(L22, D22) = ldlh( A22 - L21 D11 L21^H )    # Schur update L21 D11 L21^H: Oblio's own formUpper
                                          #   + gemmLower + gemm; then Oblio's own ldl
```

The routine choices match 0.9 exactly, `^H` for `^T`: the pivot block and the D-weighted Schur
update are Oblio's own for the same reasons, with no library equivalent, while the triangular
solve is BLAS `trsm`.

The three sections are one factorization seen through three conjugation choices. Cholesky
carries the pivot under a square root and needs `D > 0`; LDL^T and LDL^H keep the pivot in `D`
and so reach the indefinite case, differing only in whether the transpose conjugates. In the
code this is the single `Val` kernel with `hermitian(factorization())` selecting `^H` against
`^T`, the reason the three share their block structure and split only at the conjugation.

### 0.11 The supernode as a local 2x2 block

The block algebra above describes a whole matrix. A supernode is the same algebra applied to a
local `2x2` block, and this is the view the numeric kernels are written in. The `2x2` indices
are the supernode's own, not the global matrix's: the leading block is the supernode's front,
the block below it is the rows its columns reach into ancestors, and the trailing block is the
update it sends upward. Name the two inputs `F11` and `F21`, where `F` is the front the kernel
receives, already carrying whatever was gathered into it (in right-looking, earlier updates were
scattered in before this supernode is reached; in left-looking, `A` is assembled just in time),
not raw `A`. The two outputs of factoring are `C11`, `C21` (Cholesky) or `L11`, `D11`, `L21`
(LDL), and the single update the supernode emits is `U22`. Where `U22` goes, and that its
ancestors subtract it, is a separate matter (Section 1's traversals); here it is just what comes
out.

Two functions split the work along the block. `factorSupernode` turns `F11`, `F21` into the
factor; `updateSupernode` turns the factor into `U22`. For the three factorizations the kernel
supports:

**Cholesky** (`A = C C^H`):

```
factorSupernode:   C11 = chol(F11)      # LAPACK potrf
                   C21 = F21 C11^-H     # BLAS trsm
updateSupernode:   U22 = C21 C21^H      # BLAS herk (+ gemm), the update this supernode emits
```

**LDL^T** (`A = L D L^T`, complex-symmetric, no conjugation):

```
factorSupernode:   (L11, D11) = ldlt(F11)        # Oblio's own ldl
                   L21       = F21 L11^-T D11^-1  # BLAS trsm against L11^T, then scale by D11^-1
updateSupernode:   U22       = L21 D11 L21^T      # Oblio's own formUpper (+ gemmLower + gemm)
```

**LDL^H** (`A = L D L^H`, Hermitian, conjugating):

```
factorSupernode:   (L11, D11) = ldlh(F11)        # Oblio's own ldl
                   L21       = F21 L11^-H D11^-1  # BLAS trsm against L11^H, then scale by D11^-1
updateSupernode:   U22       = L21 D11 L21^H      # Oblio's own formUpper (+ gemmLower + gemm)
```

The shape is one across all three: factor the pivot block, solve the off-diagonal against the
pivot's (conjugate) transpose, emit `U22` as the outer product of the off-diagonal block. They
differ only in the two independent choices the block algebra already drew out. Whether a middle
`D` is present: Cholesky folds the pivot under the root, so `U22` is a plain rank-k product and
BLAS `herk` computes it; LDL keeps `D11` in the middle, which no BLAS routine accommodates, so
the update forms `D11 L21^H` in a scratch with `formUpper` first. And whether the transpose
conjugates: `^H` for Cholesky and LDL^H, `^T` for LDL^T, which in the code is the
`hermitian(factorization())` flag driving `maybeConjugate`. One factor per supernode; one `U22`
emitted, then routed and subtracted by whichever ancestors it reaches.

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
tool the later proofs (2.4, 2.6) rest on, so it is worth establishing once here.

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
reappears in 2.6 as the edges of the update DAG, the fill-path theorem, that
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
prove the two inclusions `RHS is a subset of Struct(j)` and `Struct(j) is a subset of RHS`.

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

`r` is the node the climb has reached. The name is for *root*, which is what it holds when
the loop ends here, and the plain climb is where that reading is literally true. So it names
a postcondition, not an invariant: mid-climb `r` is just the current node, and in the
compressed form below it can finish on a node that is not a root at all. That is the right
trade, and deliberately so. A variable named for the loop's goal (`result`, `best`, `min`
are all false until the last iteration) keeps the purpose in view, while one named for its
interior, `current` say, would be accurate at every step and tell the reader nothing. Where
the postcondition fails, in the compressed form, the failure is a fact about the algorithm
and not a defect in the name: it is exactly the point that path compression halts the climb
short of the root, and a name that quietly throws that into relief is doing its job.

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
scratch we discard. `t` is a *temporary*, and it exists only because compression is
destructive: `ancestor[r]` must be saved before it is overwritten with `k`, or the loop
would lose the very hop it is about to follow. The plain climb needs no such variable
(`r = parent[r]` suffices), so `t` is an artifact of compression rather than of climbing.
Path compression drops the cost to near-linear, `O(nnz * alpha(n))`. Both algorithms
produce the identical `parent[]`, the plain form is the one to reason about for
correctness, Liu's is what production code runs.

**The two end states, and four ways to test them.** Every climb finishes in one of two
cases. Either `j` reaches a fresh subtree root, which then takes `k` as its parent
(**attach**), or `j` already lies under `k` from an earlier neighbor this same pass, so
nothing is added (**skip**). Detecting the skip can live either in the loop's exit or in
the attach test, which yields four forms, two climbing `parent[]` (plain) and two climbing
`ancestor[]` (compressed). Writing `X` for the climbed array:

```
form                        loop continues while         attach when
--------------------------  ---------------------------  ------------------
plain,    single (classic)  parent[r] != NIL             r != k
plain,    double            parent[r] != NIL and != k    parent[r] == NIL
compress, single            ancestor[r] != NIL           r != k
compress, double (classic)  ancestor[r] != NIL and != k  ancestor[r] == NIL
```

The pairing is symmetric, and one sentence explains it. **The double-condition loop detects
the skip in its own exit test, so the attach test merely reads that detection back
(`X[r] == NIL`). The single-condition loop makes no such detection, so the attach test must
derive it (`r != k`).** Loop form fixes the attach test, and the rule is the same whichever
array is climbed. Trace the two end states:

- **Attach** is identical in every form. The climb halts at a node with `X[r] == NIL` that
  is not `k`, and writes `parent[r] = k`.
- **Skip** is where they differ. The single-condition loop climbs all the way to the
  subtree root, which on a skip *is* `k`, so it detects the skip as `r == k` and the attach
  test `if r != k` rejects it. This is why `X[k]` must stay NIL through iteration `k`: a
  column is linked only by a later column, so its slot is empty during its own turn and the
  walk may safely run onto it. (The double-condition loops never reach `k` at all, stopping
  one node short at the vertex already pointing at `k`, so for them that fact is not
  needed.) The double-condition loop turns the skip into its second exit, and the attach
  test `if X[r] == NIL`, false there, rejects it.

**Correctness, and one asymmetry between the arrays.** The single-condition forms need
`r != k`, since `X[r] == NIL` is the loop's own exit condition and would therefore always
be true. The double-condition forms accept `X[r] == NIL`. What about `r != k` after a
double-condition loop? It is correct on `parent[]` and **wrong** on `ancestor[]`, and the
reason is the difference between the two arrays:

- `parent[]` is a tree, so every edge joins **adjacent** levels. `parent[r] == k` therefore
  says `r` is a *child* of `k`, which during pass `k` can only be a root attached this same
  pass. "Already under `k`" and "is a root" coincide, so `r != k` still identifies an attach
  target. It writes `parent[r] = k` over a slot that already holds `k`: redundant, but
  harmless.
- `ancestor[]` is a shortcut structure, and compression rewrites whole paths to point
  straight at `k`, so an edge can span **arbitrarily many** levels. `ancestor[r] == k` then
  says only "`r` lies somewhere under `k`", which is true of interior nodes as well as
  roots. The two notions come apart, `r != k` can no longer tell them apart, and attaching
  at an interior node overwrites a genuine forest edge. Only `ancestor[r] == NIL` still
  identifies a root, because that is the one property compression never forges.

So three of the four cells admit `r != k`; compressed with a double condition is the sole
exception, and it is the form production code runs.

**What optimality means here.** Correct is not the same as clean, and neither is quite the
same as fast. Each single-condition form does one redundant thing on a skip, the price of
keeping the loop condition simple: a redundant *read* in the plain case (stepping onto `k`
to find `parent[k] == NIL`) and a redundant *write* in the compressed case (`ancestor[r] = k`
onto a slot already holding `k`). The double-condition forms do nothing redundant, but their
loop condition is itself twice the test, evaluated on every iteration, to save one iteration
at the end. Which is faster is therefore not obvious, and it hardly matters: this pass is
`O(nnz * alpha(n))` and vanishes beside the factorization it prepares. Trading a redundant
operation for a simpler loop is an ordinary thing to do in code. For *exposition*, though,
the choice is clear: the form with no redundancy is the one that explains itself, because
every step in it has a reason.

Read that way, the four cells are:

- **Plain, single loop, `r != k`** *(classic; pedagogical).* Climb `parent[]` to the very
  top, then ask whether the top is `k`. It is the mental model written out: a subtree has
  one root, find it, and either it is already ours or we adopt it. The redundancy is a
  single extra read, which costs essentially nothing. This is the form to reason about and
  to prove things with, and it is why we introduced it first. It is not the form to run.
- **Plain, double loop, `parent[r] == NIL`** *(the plain optimum, rarely written).* Stops
  one node short, on `parent[r] == k`, and reads that stop back as the skip. No redundant
  read, no redundant write, the only cell in the plain row that does no wasted work. It is
  cleaner than the plain classic and is seldom seen, because plain climbing is the
  pedagogical form and nobody optimizes the version they do not intend to run. Here
  "historical" is the honest explanation.
- **Compressed, single loop, `r != k`** *(the trap).* Legal, and it looks like the natural
  companion to the plain classic: same attach test, same shape. But the extra step onto `k`
  is no longer a free read, it is `ancestor[r] = k` written over a slot that already says
  `k`. Compression turns the plain form's harmless overshoot into a redundant store, on
  every skip, in the inner loop of the hot path. The cell that most tempts one to preserve
  the symmetry is the one where the symmetry stops paying.
- **Compressed, double loop, `ancestor[r] == NIL`** *(classic; principled, and forced).*
  The only correct compressed form with an early exit, for the reason above: `ancestor[r] ==
  NIL` is the one predicate compression cannot forge. The double condition detects the skip
  and the NIL test reads that detection back, so nothing is read or written twice. This is
  what Liu, CSparse, and Oblio run, and the reason is structural rather than conventional.

The classics sit on the diagonal, and the diagonal turns out to hold two cells chosen for
two different reasons wearing the same name. The compressed classic is **principled**: it is
forced, being the only legal compressed form with the early exit, and it happens also to be
the one with no redundancy, hence the cleanest to explain. The plain classic is
**pedagogical**: it is not the plain row's optimum, merely its clearest statement, and since
plain climbing is never what production code runs, nobody has had reason to prefer the
tighter spelling. Oblio ports the classic compressed form, which is both the fast choice and
the clean one.

### 2.5 Column counts of L, without computing L

The forest gives the parent links. The next thing we want is the **size** of each column of
`L`, and we want it **without computing the column itself**.

Be precise about the diagonal, since it is one entry and every count here is off by one if we
are careless. `Struct(j)` is *strictly* subdiagonal (2.1): `Struct(j) = { i > j : L[i][j] != 0 }`.
Column `j` of `L` therefore holds `{j} union Struct(j)`, so its size is `|Struct(j)| + 1`. Both
numbers are wanted, and they differ by one: symbolic factorization stores the column, diagonal
included, while a supernode's *update size* (Section 4) counts only what lies below its columns.

The reason is allocation. If the size of a structure can be computed before the structure
itself, it is nearly always worth doing: allocate exactly once, fill, and never grow. Guessing
and growing costs reallocation and copying, and over-allocating costs memory that a sparse
solver does not have to spare. This is not a fact about any particular layout. One flat buffer
needs the total length up front; a vector of vectors wants each column reserved to its exact
size; either way the counts come first. And they are cheap: one extra pass and `O(n)` of
scratch, against a structure whose size can run to many times `nnz(A)`.

**The pruned row subtree.** Fix a row `k`. Which columns `j < k` have `L[k][j] != 0`? By the
containment theorem (2.3), if `L[k][j] != 0` then `k` is in `Struct(j)`, so `k` is a proper
ancestor of `j`, and moreover every column on the tree path from `j` up to `k` also has a
nonzero in row `k`. Corollary 2 gives the converse direction: row `k`'s nonzeros are generated
by the `A`-neighbors of `k` below it, and propagate upward along tree paths.

So the columns holding a nonzero in row `k` are exactly

```
T(k) = the union of the tree paths from j up to k, over all j < k with A[k][j] != 0
```

a subtree of the elimination forest rooted at `k`. It is called the **pruned row subtree** of
`k`: pruned, because it is not the whole subtree below `k`, only the part reachable upward
from `k`'s own `A`-neighbors.

Two readings of the same set, and both matter:

- **By row:** `T(k)` is the sparsity pattern of row `k` of `L`.
- **By column:** column `j` gains one nonzero for every `k > j` whose pruned row subtree
  contains `j`. So `|Struct(j)|` is the number of subtrees `T(k)` that `j` belongs to.

The second reading is the algorithm. Walk every `T(k)`, and each time a column is visited,
credit it with one more nonzero.

**The algorithm.**

```
columnCounts(A, parent) -> colSize:
    for j = 1 .. n:
        colSize[j] = 1 # the diagonal, which every column has
        mark[j] = NIL

    for k = 1 .. n: # k roots the current pruned row subtree
        mark[k] = k # k is in its own subtree, and this halts every climb below
        for each j < k with A[k][j] != 0: # an A-neighbor of k, below it
            r = j
            while mark[r] != k: # climb until we re-enter the part already marked
                colSize[r] = colSize[r] + 1 # column r gains row k
                mark[r] = k
                r = parent[r]
```

The marker array does two jobs at once, and it is worth separating them. It **prevents double
counting**: two different neighbors of `k` may climb into a shared upper path, and the second
climb must not credit those columns twice. And it **terminates the climb**: since `mark[k]`
was set to `k` before the neighbor loop, a climb that reaches `k` stops there, with no test
for `r == k` needed. One array, one test, both problems.

**This is a companion to 2.4, and it inherits 2.4's access pattern.** The two algorithms have
the same shape: sweep `k` upward, and for each earlier `A`-neighbor `j < k` of `k`, climb the
forest from `j`. In 2.4 the climb attaches a subtree root under `k`; here it credits every
column it passes. Both read the same thing, the *earlier neighbors of row `k`*, and so both are
**row-oriented**, which makes them the exception among the phases (symbolic and numeric
factorization march over columns). The access-orientation argument of 2.4 therefore applies
here unchanged: the pseudocode is written `A[k][j]` for exposition, but with `A` stored fully,
row `k`'s earlier neighbors are exactly column `k`'s **above-diagonal** entries, so an
implementation reads `A[j][k]` and one column-oriented access serves every phase. Full storage
of `A` pays for itself twice here, once in 2.4 and once in 2.5.

**The marker is doing what a union does.** The counting algorithm is a lightweight symbolic
factorization: the same skeleton of *visit an index, mark it so it is handled exactly once,
act on it*, with the action reduced from **store** to **count**. Symbolic factorization unions
index sets, and its marker is what makes the union idempotent, an index already in the set is
skipped. Here the marker makes the *count* idempotent, an index already counted is skipped.
Detecting that two visits refer to the same index is the whole content of a union; the only
difference is that we discard the index's identity and keep only the tally.

That is exactly what makes the two-pass scheme work. Run the walk once to learn the sizes,
allocate, then run it again to store the identities in the space just allocated. The second
pass is the *row-oriented* symbolic factorization, and it is the natural partner of a count
pass, since both traverse the pruned row subtrees. (The symbolic factorization of Section 3
takes the other route: it fixes a column and gathers from that column's children, rather than
fixing a row and climbing. Same object, same marker trick, opposite traversal. It is the one
we implement, and it needs the counts of this section just the same, to size what it fills.)

**The marker is the same array in both, with its two roles exchanged.** This is the sharpest
way to see the duality, and it is worth writing out:

```
here (row-oriented):     mark[column] = row       "column r is already counted for row k"
3.1 (column-oriented):   mark[row]    = column    "row i is already in Struct(j)"
```

One algorithm fixes a row and marks the columns it reaches; the other fixes a column and marks
the rows it reaches. Same array, same size `n`, same `O(1)` duplicate rejection, and the two
are mirror images because they traverse the same object from opposite ends.

**And in both, the marker is stamped rather than set, which is what keeps the cost linear.**
The natural thing to reach for is a boolean array, cleared before each iteration. That would
cost `O(n)` per iteration and `O(n^2)` overall, swamping the `O(nnz(L))` the algorithm is
supposed to run in. Stamping with the *identifier of the current iteration* (the row `k` here,
the column `j` in 3.1) removes the clearing entirely: a stale stamp from an earlier iteration
is simply not equal to the current one, so it reads as unmarked. The array is initialized once,
to a value no iteration uses, and never touched again. That is not a micro-optimization, it is
the reason the mark array is used at all rather than a set.

**The diagonal is a choice.** Initializing `colSize[j] = 1` counts the diagonal, giving
`|Struct(j)|` in full. Initializing to `0` instead gives the count of *subdiagonal* entries,
which is what a supernode's **update size** is (Section 4). The two differ by exactly one and
either is a line's edit; which to compute is a matter of what the consumer wants, not of the
algorithm.

**Cost.** Every step of the inner climb either credits a column (a nonzero of `L` that is
counted exactly once, thanks to the marker) or terminates. So the total work is proportional
to the number of nonzeros of `L`, plus one terminating step per edge of `A`:

```
O(nnz(L) + nnz(A))  =  O(nnz(L))
```

That is optimal in the sense of being linear in what is counted, but note what it implies:
counting `L` costs as much as touching `L` would. The saving is in **memory**, not time. We
learn the size of a structure we never store, using `O(n)` scratch.

There are asymptotically better algorithms. Gilbert, Ng and Peyton (1994) compute the same
counts in nearly linear time in `nnz(A)`, not `nnz(L)`, using the skeleton matrix and least
common ancestors, and never touching a fill entry. For a solver this matters when `L` is much
denser than `A`, which is the usual case. We do not use it here; the point of this section is
the pruned row subtree, and the simple walk makes it visible.

**Worked example: the grid.** Take the grid of 2.7 again, with
`parent = [2, 3, 7, 5, 6, 7, 8, 9, NIL]`. Trace two of the pruned row subtrees.

`k = 7`. Its `A`-neighbors below it are `1` and `4` (grid edges `1-7` and `4-7`). Mark
`mark[7] = 7`, then:

```
j = 1:  r = 1 -> credit column 1, mark it 7
        r = parent(1) = 2 -> credit column 2, mark it 7
        r = parent(2) = 3 -> credit column 3, mark it 7
        r = parent(3) = 7 -> mark[7] == 7, stop
j = 4:  r = 4 -> credit column 4, mark it 7
        r = parent(4) = 5 -> credit column 5, mark it 7
        r = parent(5) = 6 -> credit column 6, mark it 7
        r = parent(6) = 7 -> mark[7] == 7, stop
```

So `T(7) = {1,2,3,4,5,6,7}`: both chains, which is exactly row 7 of `L`, and it says that
every one of the six subdomain columns fills into the separator's first column.

`k = 8`. Its `A`-neighbors below it are `2`, `5` and `7`. Mark `mark[8] = 8`:

```
j = 2:  credit 2, 3, 7 (climbing 2 -> 3 -> 7 -> 8, stopping at 8)
j = 5:  credit 5, 6      (climbing 5 -> 6 -> 7, and 7 is already marked 8: stop)
j = 7:  r = 7 is already marked 8: stop immediately, credit nothing
```

Both jobs of the marker are on display. Column `7` is reached three separate times and
credited **once**. And the third neighbor, `7` itself, terminates before doing anything, since
its climb was already absorbed into the first.

Running every `k` gives:

```
column:      1    2    3    4    5    6    7    8    9
|Struct|:    3    4    4    3    4    4    3    2    1     total 28
subdiagonal: 2    3    3    2    3    3    2    1    0     total 19
```

The 28 is the length of the flat buffer symbolic factorization will need. Note it counts the
9 diagonal entries plus 19 subdiagonal ones, and the fill is already in there: the original
matrix has only 12 edges, so 12 subdiagonal nonzeros, and the other 7 are fill (matching the
count in 2.7).

### 2.6 The forest as transitive reduction of the update DAG

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

### 2.7 Worked example: nested dissection of a grid

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
columns are exactly the forest leaves, so the length is (edges) + (leaves): here 12 + 2 = 14.

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
  2.6 made concrete: the forest keeps one edge per column (the minimal one), so a column's
  other fills are transitively implied by climbing. The 7 fill entries collapse to 2 fill
  tree-edges.

That last note is a nice payoff, it makes "forest = transitive reduction" (2.6)
concrete: 7 fill entries collapse to 2 fill tree-edges.

Contrast: the *row-major* numbering of the same grid collapses the forest to a single
**path** `1->2->...->9` (`parent = [2,3,4,5,6,7,8,9,NIL]`), height 8, no parallelism,
fill spread across the whole band. Nested dissection turns that path into a balanced tree
of height 6 with two independent halves. Same matrix, same theory, better tree, the
ordering is what the forest sees.

### 2.8 A lighter ordering of the same grid

The 2.7 ordering numbered each subdomain column top to bottom (`1,2,3` down the left),
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

Against 2.7, the same matrix under this ordering gives **5 fills instead of 7**, and a
shorter, bushier forest (height 5 rather than 6). The reason is exactly the
sub-separators: 2.7's left block was the chain `1 -> 2 -> 3`, whereas here `1` and `2` are
independent leaves under `3`, so the column that was a path becomes a two-leaf fork, one
fewer coupling to fill, per subdomain. The sub-separator shows in the elimination box:
eliminating the two ends `1, 2` creates `(3,7), (3,8)`, and their meeting point `3` then
creates `(7,8)`, confining each subdomain's interaction to a single vertex. Numbering the
separators last, recursively, is the whole idea of nested dissection, this subsection is
one more level of it applied to 2.7.

One subtlety is worth drawing out: 2.8 has *fewer* fill entries (5 vs 7) yet *more* fill
forest edges (3 vs 2). These count different things, total nonzeros created in `L`,
versus columns whose parent edge (first below-diagonal nonzero) is fill. Split the fill
into on-tree (a column's minimal entry, so a forest edge) and off-tree (higher entries,
above the parent) and they reconcile:

```
            on-tree fill    off-tree fill    total
2.7              2                5            7
2.8              3                2            5
```

2.7's chains produce mostly off-tree fill, each chain column fills into several separator
rows, most of them above its parent, redundant climbs high in the band. 2.8's
sub-separators cut those columns to near-single below-diagonal entries, concentrating fill
at the parent. So a tighter ordering does two things at once: it lowers total fill, and it
pushes the remaining fill toward the forest's transitive reduction (2.6), a larger share
of it becomes minimal (on-tree), leaving fewer transitively-implied off-tree entries. Fill
forest edges going up while total fill goes down is a signature of that.

### 2.9 Connection to the traversals

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

Two things about that array, both established in 2.5 and both easy to pass over. The stamp is
the *column being built*, not a boolean, which is what removes the clearing pass: a stale stamp
from an earlier column is simply not equal to `j`, so it reads as unmarked, and the array is
initialized once rather than `n` times. Without that the phase would be `O(n^2)`, not
`O(nnz(L))`. And the array is the exact mirror of the one in 2.5: here `mark[row] = column`,
there `mark[column] = row`, because that algorithm fixes a row and marks the columns it reaches
while this one fixes a column and marks the rows. Same trick, same array, roles exchanged.

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

Take the grid of 2.7, forest `parent = [2, 3, 7, 5, 6, 7, 8, 9, NIL]`. First read each
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
exactly the `X`/`F` positions of `L` in 2.7, column 2 gains row 7 as fill from its child
1, column 3 inherits `{7,8}` from column 2, and at column 7 the two subtrees (children 3
and 6) merge, each contributing `{8,9}`. The recurrence reconstructs the entire fill
pattern from `A` and the forest, with no arithmetic.

Two features of the forest show through the trace. Fill enters exactly at the multi-hop
climbs of 2.7: a column's pattern grows only where a child hands up rows that were not in
`A`'s column. And the two independent subtrees `{1,2,3}` and `{4,5,6}` are built without
reference to each other, column 7 is the first to see both, which is the sibling
independence of 2.9: the symbolic pass over disjoint subtrees is itself independent work.

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

Fundamental supernodes are **unique** for a given forest. Conditions 1 and 2 are properties,
not choices: there is nothing to decide and no heuristic to tune. That is why they are the
natural default, and it is worth being precise about what it costs.

They are **not maximal**. Condition 1 is stricter than the sparsity requires: a parent with
two children is refused outright, even when one of those children shares its pattern exactly
and could be absorbed at no cost whatever. The only-child condition is there to make the
supernodes *paths*, hence unique, not because merging would introduce fill. Section 4.5 drops
it, and pays for the extra merges in uniqueness.

Write `supernode(j)` for the one containing column `j`. Its columns are the **front indices**;
the rows below them, common to all of them, are the **update indices**.

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
than single columns, the same identity reads (uppercase now, since these are supernodes):

```
|update(J)| = |front(K)| + |update(K)|
```

because the rows below supernode `J` are the columns of `K` followed by the rows below `K`.
For a nodal forest `|front(K)| = 1` and this collapses to the test above. The general form
is worth keeping: it is what lets compression run on a forest that has already been
compressed. Note it is the **parent's** front size that appears, not the child's.

### 4.3 Compression

One pass in increasing column order suffices. Because `parent(j) > j`, a column's child is
always numbered before it, so when column `k` is reached its child's supernode is already
decided.

```
compress(forest, Struct) -> supernode: # fundamental supernodes
    S = 0 # supernodes so far; the first one created is numbered 1
    for k = 1 .. n: # increasing order = children before parents
        j = firstChild[k]
        if j != NIL # k has a child at all
               and j == lastChild[k] # condition 1: and only the one, j
               and |Struct(k)| == |Struct(j)| - 1: # condition 2: sharing one pattern
            supernode[k] = supernode[j] # k continues j's supernode else:
            S = S + 1
            supernode[k] = S # k starts a new supernode
```

Both conditions of 4.1 are in that test, and the guard in front of them is not decoration.
`firstChild[k] == lastChild[k]` says "the child list of `k` has a single entry", which is
why the forest keeps a last-child link and not only a first-child link. But it is *also*
true when both are NIL, that is, at every leaf, so without the `j != NIL` guard every leaf
would report exactly one child and the pattern test would go on to ask for `Struct(NIL)`.
The guard is what makes the idiom mean what it reads as. The size equality is then the
pattern test of 4.2.

That assigns every column to a supernode and counts them. The forest must then be rebuilt
over supernodes rather than columns. A second pass, in *decreasing* order, does it:

```
rebuild(forest, supernode) -> supernodal forest:
    for S = 1 .. numSupernodes:
        front(S) = 0
        done[S] = false
    for k = n .. 1: # decreasing order, so a supernode's topmost column comes first
        S = supernode[k]
        front(S) = front(S) + 1 # every column of S adds to its front size
        if not done[S]: # k is the topmost column of S
            if parent(k) != NIL:
                supParent(S) = supernode[parent(k)] else:
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
gather its front indices through `supernode[]` rather than assume a contiguous range. (A
postordering of the forest would make them consecutive, and is worth doing for other
reasons, but the algorithms above do not depend on it.)

### 4.4 Worked example: the grid

Take the grid of 2.7 again, with `parent = [2, 3, 7, 5, 6, 7, 8, 9, NIL]` and the patterns
computed in 3.2. Tabulate what the test needs, the only child (if any) and the pattern sizes:

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
supernode:   1    2    3    4    5    6    7
columns:    {1}  {2}  {3}  {4}  {5}  {6}  {7,8,9}
front:       1    1    1    1    1    1    3
update:      2    3    3    2    3    3    0
```

The one non-trivial supernode is exactly the **separator**. That is not a coincidence, it
is nested dissection working as designed: the separator is eliminated last, by then it is
completely filled in, and a completely filled-in block is precisely a set of columns with
one shared pattern. Its index set is `{7,8,9}` with no update rows, a dense 3x3 trailing block.

The left and right blocks do not compress at all, and the trace shows why. Column 2's
pattern is `{3,7,8}` while its child column 1 has `{2,7}`: dropping the 2 leaves `{7}`,
not `{7,8}`. The patterns differ, so the columns cannot share an index list. Column 2
*acquires* row 8 from `A` that column 1 never had. Only once a column's pattern has stopped
growing do the columns start agreeing, and in this ordering that happens only inside the
separator.

### 4.5 Amalgamation: buying bigger blocks with explicit zeros

Fundamental supernodes are free but timid. This section relaxes them in two directions at
once, and a single parameter controls both.

**The fill of a merge.** Let `K` be a supernode and `J` one of its children. Merging `J` into
`K` makes `J`'s columns part of `K`'s front, so each of them must now store the whole of `K`'s
index set below it, where before it stored only `J`'s own update rows. By the containment
theorem (2.3), `update(J)` is a **subset** of `K`'s index set, so the new entries per column of
`J` are exactly

```
indexSet(K) \ update(J)      of which there are   |indexSet(K)| - |update(J)|
```

and every one of them is a **zero**, stored explicitly because the block must be rectangular.
Multiply by the number of columns `J` contributes:

```
fill(J -> K)  =  ( |indexSet(K)| - |update(J)| ) * |front(J)|
```

Being a set-difference size, it is never negative: containment guarantees it.

**Zero fill is the pattern test.** The expression vanishes exactly when
`update(J) = indexSet(K)`, which is condition 2 of 4.1, the identical-pattern condition. So a
zero-fill merge and a fundamental merge ask the *same* question about patterns. The only
difference is condition 1: fundamental compression additionally demands that `J` be `K`'s
**only** child. Drop that demand and we merge strictly more, at no cost.

**And that changes the answer.** Take the smallest case, a three-column star, `A` with edges
`1-3` and `2-3`:

```
Struct(1) = {1,3}    Struct(2) = {2,3}    Struct(3) = {3}
parent(1) = 3        parent(2) = 3        parent(3) = NIL
```

Column 3 has **two** children, so fundamental compression merges nothing: three supernodes.
Now compute the fills. `indexSet(3) = {3}`, of size 1, and `update(1) = update(2) = {3}`, also
of size 1, so

```
fill(1 -> 3) = (1 - 1) * 1 = 0
fill(2 -> 3) = (1 - 1) * 1 = 0
```

Both children merge for free. But not both *together*. Absorb column 1 first, and supernode
`{1,3}` now has a front of two columns, so its index set has grown to `{1,3}`. Recompute:

```
fill(2 -> {1,3}) = (2 - 1) * 1 = 1
```

Column 2 would have to store a zero at row 1. So one child merges and the other does not:
**two supernodes, no fill**, where fundamental compression gave three.

**Uniqueness is the price.** Which child? Both were equally free, and the choice mattered:
taking column 2 first yields `{2,3}` and `{1}`, a different partition with the same count and
the same (zero) fill. There is no canonical answer, so the algorithm must **break ties**, and
any tie-breaking rule is a convention rather than a theorem. The usual one is: least fill
first; ties to the largest front size; ties again to the first child in the list. That last
rule is arbitrariness made deterministic, which is exactly what a canonical algorithm never
needs. Fundamental supernodes are unique because they *refuse* to make this choice.

**What threshold zero leaves behind.** The zero-fill expression `|indexSet(K)| - |update(J)|`
is always `>= 0` by containment, and equals zero exactly at an identical-pattern pair, the free
merge. A fundamental forest can still hold such pairs: the star's root has two children, each at
zero, that fundamental compression could not take because neither is an only child.
Threshold-zero amalgamation exists to consume exactly these. It absorbs one free child, which
widens the parent's front and lifts every remaining sibling off zero (in the star, the second
child jumps from zero to one), then looks again, and stops only when no child sits at zero any
longer. So its post-condition is clean: **after threshold-zero amalgamation, no parent-child
pair has zero fill.** Equivalently, no two supernodes still in a parent-child relation share the
identical-pattern property, because if they did the merge would have been free and already
taken. Every surviving pair is strictly positive, which is why raising the threshold above zero
is what it takes to merge anything more.

**Why one pass suffices.** That the equality cases vanish can be seen directly, without appeal
to the loop's mechanics. Run one pass over the children at threshold zero: the merge condition
fires exactly on the zero-fill pairs, so every equality present is caught and merged in that
pass. Could a later pass find a fresh equality? It would have to have either existed during the
first pass, in which case it was not missed and was already taken, or been created by one of the
first pass's merges. But a merge only ever widens the parent's front, and a wider front raises
every remaining child's fill or leaves it unchanged; it can never lower a positive fill back to
zero. So no new equality can be manufactured. Neither source is possible, and the equality cases
are therefore gone after a single pass. The same monotonicity, that a merge never decreases
another child's fill, is what lets a child ruled out of budget be struck off for good.

**Spending fill on purpose.** So far the threshold has been zero. Raising it lets a merge
through even when it costs zeros, and the reason is Section 4.6: a supernode's nonzeros form
a dense block, and dense blocks run on BLAS level 3. A wider block is a better block, and past
some point the arithmetic on a few explicit zeros is cheaper than the overhead of processing
two narrow blocks instead of one. The threshold is the budget: **how many explicitly stored
zeros are we willing to buy a bigger front with**. It is a genuine tuning parameter, with no
right answer independent of the machine.

**The algorithm.** Greedy, one parent at a time.

```
amalgamate(forest, threshold) -> supernode:
    for K = 1 .. numSupernodes:
        candidate[K] = true

    for K = 1 .. numSupernodes: # increasing: children before parents
        fillInc = 0 # zeros already bought for K
        frontInc = 0 # columns already absorbed into K
        kSize = |front(K)| + |update(K)| # K's index set, before growth

        loop:
            best = NIL
            for each child J of K:
                if not candidate[J]:
                    continue
                zerosPerCol = kSize + frontInc - |update(J)| # >= 0, by containment
                cost = fillInc + zerosPerCol * |front(J)|
                if cost > threshold:
                    candidate[J] = false # over budget now, and only ever more so
                    continue
                if best == NIL or cost < bestCost or
                   (cost == bestCost and |front(J)| > |front(best)|):
                    best = J
                    bestCost = cost

            if best == NIL:
                break # nothing more fits the budget

            supernode[best] = K # absorb it
            fillInc = bestCost
            frontInc = frontInc + |front(best)|
            candidate[best] = false

        front(K) = front(K) + frontInc # K has grown; later merges must see this
```

Three details are load-bearing.

**`frontInc` must be inside the cost.** Every child absorbed widens `K`'s front, which raises
the price of the next one. That is what stopped column 2 in the star, and it is why the loop
recomputes the cost on every pass rather than ranking the children once.

**`candidate[J] = false` is permanent.** A child that fails the test can never pass it later,
since `K` only grows. Striking it off keeps the greedy loop from rescanning it, and turns what
looks quadratic into something close to linear in practice.

**Updating `front(K)` at the end is not bookkeeping.** The parents are processed in increasing
order, so by the time `K` is reached its children have *already been parents themselves* and
may have absorbed children of their own. The cost of merging `J` into `K` depends on `J`'s
**current** front size, not the one it started with. Forgetting this line underestimates both
the fill and the resulting block, and the error is silent.

**A consequence for storage.** 4.3 noted that a supernode's columns need not be consecutive,
and that a postordering would make them so. That escape closes here. A fundamental supernode
is a *path*, and a postorder numbers a path consecutively. An amalgamated supernode is a
parent plus a **subset** of its children, and a postorder places the parent after *all* of its
children, including the ones that were not merged. The star gives the smallest instance:
supernode `{1,3}` skips column 2, and no renumbering fixes it, because 2 must precede 3 and
cannot precede 1. So consumers must gather a supernode's front indices through `supernode[]`
and never assume a contiguous range. Both the symbolic factorization of 4.6 and the numeric
one depend on this.

**What does survive, and why the sort is safe.** Contiguity is lost, but a weaker property is
not, and it is the one that actually matters: within each supernode, every front index is
smaller than every update index. A supernode's update indices are the rows below its whole
front, which is to say the ancestors of its top column. Under a postorder the top column is the
parent, the **largest** of the front columns, and an ancestor of the parent is numbered above
it. So every update index exceeds the parent, which in turn is at least every front column:
front-max is below update-min, always. Amalgamation does not touch this. It scatters the front,
a non-merged sibling can sit numerically between two merged columns, but that gap is filled by
some *other* supernode's column, never by this supernode's own update rows, which remain
strictly above the parent. The star shows it: supernode `{1,3}` has a gap at 2, but 2 is a
different supernode's column, not an update row of `{1,3}`.

This is why the indices can be sorted front and update **together**, in one pass, without
interleaving the two groups. It is worth being clear about what is required and what is chosen.
The factor is retrieved entirely through indirection, an entry found by its row and column
index, so **the indices do not need to be sorted at all**; the one real requirement is that the
front indices be **separated** from the update indices, front before update. Sorting is a
convenience laid on top of that requirement, not the requirement itself. We choose to sort
anyway: it is linear time (the transpose-based bucket sort of 4.6, cheap next to the
factorization), and a factor whose indices are in order is easier to read and to check by eye.
The separation comes for free from front-max below update-min; the ordering within each group
is the part we pay a small, optional price for. (One later wrinkle, noted here so it is not a
surprise: dynamic pivoting can reorder the front indices of a supernode, though never the update
indices; the numeric factorization tolerates this precisely because order was never required,
only separation.)

**Where it fits.** Amalgamation runs on any forest, nodal or supernodal, since a nodal forest
is just one whose supernodes are all trivial. In practice one runs it *after* fundamental
compression: fundamental does the free, canonical work cheaply, and amalgamation then does the
tie-broken and the paid work on a smaller forest. Running it directly on a nodal forest at
threshold zero reaches the same fill (none) but by a longer route.

**The same widening happens again later, for a different reason.** Dynamic LDL cannot always
pivot a column where the symbolic factorization placed it: a candidate pivot too small relative
to its column is **delayed** (7.6), leaving its supernode's front and joining its parent's. That
migration is mechanically what a merge is. The column leaves `J`, becomes one of `K`'s front
columns, and must now span `indexSet(K)` below it where before it spanned only what `J`'s index
set gave it. The count of new positions has the shape of the merge formula, and they arrive as
**zeros**, since the widened block is allocated zeroed and the migration copies in only the
values the column already had.

So both are the same operation seen twice: a child's columns moving into a parent's front, and
paying for it in explicitly stored zeros. What differs is everything around it.

Amalgamation is **symbolic and chosen**. It happens before any value exists, under a threshold we
set, and its zeros are *known* zeros of `L`: merging supernodes does not change `L`'s pattern, so
those entries are zero at the end of the factorization as surely as at the start. That is the
whole bargain, storage and flops spent on zeros we can name, bought back in block size.

Delay is **numeric and forced**. It happens while the arithmetic runs, no threshold makes it
optional, and its zeros are not known zeros: the column is now eliminated at `K`, so its pattern
has genuinely changed, and the symbolic factorization did not predict the result. Some of those
positions may fill during `K`'s own elimination. Nothing is bought; the cost is what stability
costs.

Two smaller asymmetries follow the same line. A merge takes **all** of `J`'s columns and removes
`J` from the forest; a delay takes only the columns that could not be pivoted and leaves `J` in
place, narrower than it was. And a merge is decided once, for a whole forest, before the numbers
arrive; a delay is decided per column, per factorization, and the same matrix with different
values delays differently.

### 4.6 What changes downstream

Symbolic factorization runs unchanged, at supernode granularity. It is the same union as 3.1,
and the three regimes differ in exactly two places: how much of `A` the supernode reads, and how
much of a child it drops. Here are all three, written out. Repetition is cheap and the
differences are the point.

Throughout, `Idx(K)` is a supernode's **index set**, `front(K) union update(K)`: its own columns
followed by the rows below them. For a nodal supernode `{k}` that is `{k} union Struct(k)`, the
column with its diagonal. (The doc does not permute, so `j` and `k` are the same indices in `A`
and in `L`.)

**Nodal.** One column per supernode. This is 3.1, restated in the index-set vocabulary.

```
symbolicFactorNodal(A, forest) -> Idx:
    for k = 1 .. n: # increasing = children before parents
        Idx(k) = { i >= k : A[i][k] != 0 } # A's column k: the diagonal, and below it
        for each child j of k:
            Idx(k) = Idx(k) union ( Idx(j) \ {j} ) # drop j's own diagonal; the rest carries up
```

**Fundamental supernodes.** The columns of a supernode share a pattern in `L`, so for any two
front columns `k' < k''` of `K`,

```
Struct(k'') is contained in Struct(k')
```

Reading the **first** (lowest) front column of `A` is therefore enough. But that inference is
less obvious than it looks, and it is worth doing carefully, because two things one might assume
are false.

**Sharing a pattern in `L` does not mean sharing one in `A`.** Take the smallest case that shows
it, three columns with `A`-edges `1-2` and `1-3`:

```
        1  --  2
        |
        3
```

Structure of `A` (`X` = nonzero; full symmetric):

```
     1 2 3
  1  X X X
  2  X X .
  3  X . X
```

Eliminating column 1 makes 2 and 3 adjacent, so `L` gains one fill entry at `(3,2)`, marked `F`:

```
     1 2 3
  1  X . .
  2  X X .
  3  X F X
```

Read the columns off both:

```
A columns:  A(1) = {1,2,3}    A(2) = {1,2}     A(3) = {1,3}
L columns:  L(1) = {1,2,3}    L(2) = {2,3}     L(3) = {3}
```

All three merge into **one** fundamental supernode: the `L` patterns nest, each is the next with
one more row on top. The `A` patterns do not resemble each other at all. And look at column 2: in
`A` it has **nothing below its diagonal**, while in `L` it holds `{3}`. That entry is fill, and
it arrives from column 2's child (column 1), not from `A` at all.

**And `A`'s column is generally a strict subset of `L`'s.** Column `k'` has fill of its own, from
its own children. So `A(k')` does not contain `L(k')`, and cannot be expected to. (In the picture
above column 1 happens to have none, being a leaf, but in general it will.)

Which raises the real question. The shortcut reads only `A`'s column `k'`. What if `A`'s column
`k''` holds an index that `A`'s column `k'` does not?

It can, and the shortcut still finds it. Let `i` be such an index: `A[i][k''] != 0`, `i >= k''`,
and `i` not in `A(k')`. Then `i` is in `L(k'')`, and by the containment above `i` is in `L(k')`.
Now `L(k')` decomposes exactly two ways:

```
L(k') = { rows of A's column k', at or below k' }  union  { update indices of k's children }
```

`i` is not in the first set, so it is in the second, **and the shortcut reads that set too**. The
step that closes it: a fundamental supernode is a *path* in which every column above the bottom
has exactly one child, and that child is the previous column of the same supernode. So the only
edges reaching in from outside land on `k'`, and the supernode's children are exactly `k'`'s
children.

Hence

```
Idx(K) = L(k') = (A's column k', at or below k')  union  (the children's update indices)
```

which is precisely what the algorithm computes. The shortcut is **not** "the later `A` columns
are redundant copies of the first", which is false. It is "the later `A` columns can only
contribute indices that `L(k')` already holds, and every index of `L(k')` arrives through one of
the two doors the union already opens".

The front indices need no seeding either, for the same reason: `L(k')` already contains every
later column of `K`, since those are rows below `k'`. The block diagonal *emerges* from the union.

```
symbolicFactorFundamental(A, forest) -> Idx:
    for K = 1 .. numSupernodes:
        k' = the lowest column of K # its first front index
        Idx(K) = { i >= k' : A[i][k'] != 0 } # A's column k' only
        for each child J of K:
            Idx(K) = Idx(K) union ( Idx(J) \ front(J) ) # drop the child's columns, not one index
```

Read the first line carefully, because the obvious gloss on it, "one column of `A` suffices", is
**false as an argument**, even though it is true as a conclusion. Getting from one to the other
is the whole content of this subsection, and it is worth the space.

**`A`'s columns of a supernode are not nested.** The three-column picture above shows it:
`A(2) = {1,2}` and `A(3) = {1,3}`, and neither sits inside the other. Only `L`'s columns nest.
So the shortcut cannot be "the later columns of `A` are redundant copies of the first". They are
not copies of anything.

**And a later column of `A` really can hold something the first does not.** The smallest case
needs four columns, with `A`-edges `1-2`, `1-4`, `2-3`, `3-4`, and crucially **no** edge `2-4`:

```
A (X = nonzero):          L (F = fill):
     1 2 3 4                   1 2 3 4
  1  X X . X              1  X . . .
  2  X X X .              2  X X . .
  3  . X X X              3  . X X .
  4  X . X X              4  X F X X
```

The supernodes are `{1}` and `{2,3,4}`, so `K = {2,3,4}` and `k' = 2`. Look at what `A` offers
each front column, below its own diagonal:

```
A(2) | >= 2 = {2,3}      <- all the shortcut reads
A(3) | >= 3 = {3,4}      <- holds 4, and A(2) does not
A(4) | >= 4 = {4}
```

`A`'s column 3 has an original entry at row 4. The shortcut never reads column 3. So where does
the 4 come from?

**From the child.** `K`'s only child is the supernode `{1}`, whose update indices are
`L(1) \ {1} = {2,4}`. So the union computes

```
Idx(K) = A(2)|>=2  union  (child's update indices)
       = {2,3}     union  {2,4}
       = {2,3,4}                                    which is L(2), and correct
```

The 4 arrives, but **not as `A`'s entry at (4,3)**. It arrives as *fill*, through the graph path
`2 - 1 - 4`: vertex 1 lies below both, so eliminating it makes 2 and 4 adjacent. That is the `F`
at `(4,2)` in the picture. The original entry and the fill entry are **two different reasons for
the same index**, and the shortcut collects it by the second reason while ignoring the first.

**The general argument, now that the example has made it concrete.** Take any `i` in `A(k'')`
with `i >= k''`. Then `i` is in `L(k'')`, since `A`'s pattern is contained in `L`'s. And
`L(k'')` is contained in `L(k')`, since they are front columns of one supernode. So

```
i is in L(k')
```

**guaranteed, by the shared pattern alone.** And `L(k')` has a complete recipe that never mentions
`A(k'')`:

```
L(k') = { i >= k' : A[i][k'] != 0 }  union  (update indices of k's etree children)
```

So `i` must arrive through one of those two doors. Which one varies, and the example shows both:
the 3 comes straight from `A(2)`, the 4 comes from the child. The shortcut does not need to know
which.

That is the honest form of "one column of `A` suffices". Not *"the other columns of `A` say
nothing new"*, which is false, but *"whatever the other columns of `A` could say about `L` is
already forced into `L(k')` by the shared pattern, and `L(k')` is computed anyway"*. The columns
of `A` are bound together by the pattern their `L` columns share.

**And it must be `k'`, not any other front column.** Two things pick it out, and both are needed.

*It sees the widest window.* `k'` reads rows `>= k'`, which includes every row a later front
column could see, plus the rows between `k'` and `k''` that a later column structurally cannot
reach. So no other column would do.

*It is where the children attach.* A fundamental supernode is a path in which every column above
the bottom has exactly one child, and that child is the previous column of the same supernode. So
`k'` is the only front column whose elimination-tree children are the **supernode's** children.
Its recurrence therefore has both its inputs to hand; a later column's recurrence would refer to
the update indices of `k'' - 1`, a column *inside* the supernode whose index set was never
computed separately.

The first says no other column would do. The second says this one does. Try `k''` on the
three-column picture, a leaf supernode with no children at all, and the failure is immediate:

```
from k' = 1:   {1,2,3} union {}  =  {1,2,3}       correct
from k'' = 2:  {2}     union {}  =  {2}           wrong: 1 and 3 are lost
```

Finally, the front indices need no seeding, for the same containment reason: `L(k')` already
contains every later column of `K`, since those are rows below `k'`. The block diagonal *emerges*
from the union rather than being put into it.

**Amalgamated supernodes.** Amalgamation (4.5) merges columns whose patterns are only *nearly*
identical, paying explicit zeros for the difference. So `Struct(k'')` is no longer contained in
`Struct(k')`, the shortcut above fails, and **every** front column of `A` must be read.

```
symbolicFactorAmalgamated(A, forest) -> Idx:
    for K = 1 .. numSupernodes:
        Idx(K) = union, over every front column k of K, { i >= k : A[i][k] != 0 }
        for each child J of K:
            Idx(K) = Idx(K) union ( Idx(J) \ front(J) )
```

**The three are one algorithm.** The last is the general form and is correct in all three
regimes: on a fundamental forest the extra front columns contribute nothing (the union simply
finds every index already present), and on a nodal forest there is only one front column and
`front(J) = {j}`, which is the first listing exactly. So an implementation can write the third
and get all three, which is what both references do, and what we do.

The saving forgone is real but modest, and it is worth being exact about what it is *not*. Both
forms need a pass over the column-to-supernode map, and neither can avoid it: the map runs the
wrong way, and a supernode's columns need not be contiguous (the forest is topological, not
postordered), so its lowest column cannot be found by inspection. That pass is the price of
admission in either regime.

What differs is what the pass must *produce*. The fundamental form needs only the lowest column
of each supernode (`numSupernodes` entries); the general form needs all front columns gathered
contiguously (`n` entries), and pays an inner loop over them in the union. So the specialization
is one of **size**, not of passes. It is available and not taken.

**What actually changed, in a table.** The nodal case is not a special case of the supernodal
one; it is the uncoarsened one.

```
                     nodal                          supernodal
own contribution     A's column k, rows >= k        A's columns of K, rows >= each, unioned
front space          {k}, one diagonal entry        the block diagonal, |front(K)| entries
children give        each child's update indices    each child's update indices
```

Three things are worth reading off that table.

**The children's half does not change at all.** A child gives its update indices and only those;
its front indices die with it, because they are the child's own columns and lie strictly above
the parent. Nodal and supernodal differ only in *how many* indices are dropped, one, or
`|front(J)|` of them.

**What generalizes is the front space**, from a single diagonal entry to a block diagonal. That
is the whole content of the supernodal step, and it is why a supernode's front indices have to
be found at all: nodally there is nothing to find, since the front of column `k` is `k`.

**And `A` contributes to both spaces, not just the front.** The union seeds a supernode with
`{ i >= k : A[i][k] != 0 }`, the diagonal *and* the subdiagonal. Fill comes from the children;
the original nonzeros come from `A`, and most of them lie below the diagonal. It is tempting to
read the parent's own contribution as "the block diagonal" and the children's as "everything
else", but that is not what happens: `A` reaches into the update space too.

The payoff is not in the symbolic phase, which was already optimal, but in the numeric one.
A supernode's columns share a pattern, so its nonzeros form a **dense rectangular block**,
and the numeric factorization can update it with dense matrix kernels (BLAS level 3)
instead of a scatter over individual columns. That is the entire motivation: supernodes
convert sparse column operations into dense block operations, which run at a large multiple
of the speed on real hardware. The larger the supernodes, the better the ratio, which is
also why threshold-based merging exists, it accepts some explicitly stored zeros to buy
bigger blocks. Fundamental supernodes are the free case: they introduce no zeros at all.


## 5. Ordering

Everything above takes the permutation as given. This section is about where it comes from, and
it is the last black box in the pipeline: `OrderEngine` calls a vendored AMD and a vendored MMD,
and nothing in the port has yet looked inside either.

The section is here for two reasons. The obvious one: the ordering decides the fill, and
therefore the forest, the supernodes, the work, and the memory. Everything downstream is
downstream *of this*. The less obvious one: the classical implementations are hard to read, and
it is worth separating what is hard because the *algorithm* is intricate from what is hard
because the *encoding* is fifty years old. They are not the same, and only one of them is worth
preserving.

**A note on where the pseudocode comes from, because it is not all the same kind.** Sections 5.1
through 5.10 build up minimum degree in layers, and their pseudocode is written *from* the
prototypes in `experiments/ordering`, which are the specification: what the pseudocode says, the
code does, and the worked examples are that code's own output. Sections 5.11 to 5.13 go the other
way. Their pseudocode is a reading of the vendored `Mmd.cpp` and `Amd.cpp`, which are the
specification there, and it describes those routines in full. Our prototypes for the two named
algorithms are deliberately incomplete, and 5.11 records exactly which parts are missing.

The arrow therefore points one way in the first ten subsections and the other way in the last two,
which is the arrangement we want rather than an oversight: it lets the layered sections be exact
about code we control, and lets the MMD and AMD sections be complete about code we do not, with
the gap between the vendored routine and our prototype written down instead of hidden.

One consequence is visible in the notation. The layered sections say *clique* and write `C_i`,
following 5.3, while 5.11 to 5.13 say *element* and write `E_i`, following their sources. The
translation is one word and one letter, given in 5.3, and it is kept deliberately so that a reader
who opens `Mmd.cpp` after 5.11 meets the vocabulary the section just used.

### 5.1 Minimum degree, and why the obvious version does not work

The greedy idea (Tinney and Walker, 1967) is one line: **eliminate the vertex of least degree
next.** A vertex of degree `d` creates a clique on its `d` neighbors, so `d` bounds the fill it
causes. Picking the smallest causes the least immediate fill.

It is a heuristic, not an optimum. Minimizing fill exactly is NP-hard, and minimum degree is
greedy on a local quantity. It is nevertheless very good in practice, and it is what AMD and MMD
both are.

Written out literally, on a graph that stores its edges, the whole algorithm is this:

```
minimum_degree(G) -> elimination order:

    while G has a live vertex:
        p = a live vertex of least degree             # the greedy pick
        for each pair (u, w) of distinct neighbors of p:
            if u and w are not adjacent:
                add the edge (u, w)                   # FILL
        remove p and its incident edges from G
        append p to the order
```

That is the vertex elimination game, and it is the definition the rest of Section 5 optimizes
without changing what it computes.

It is worth counting what a single elimination costs the graph, because the count is the reason
the heuristic is *minimum* degree rather than anything else. Eliminating `p` of degree `d` deletes
exactly `d` edges, the ones joining `p` to its neighbors, and adds at most `d(d-1)/2`. The two
cross at `d = 3`. A degree-2 pivot deletes two edges and adds at most one, so eliminating it
strictly shrinks the graph; a degree-3 pivot is break-even at worst; only from `d = 4` up does
elimination generally grow the edge count. Choosing the least-degree vertex is therefore choosing
a pivot likely to shrink the graph rather than grow it, which is a stronger statement than the
usual one about `d` merely bounding the fill.

The `at most` in that count carries the rest of the story. The pivot only pays for the pairs of
neighbors that were *not already adjacent*, so an elimination whose neighborhood is already a
clique creates no fill at all and is pure deletion, free. Minimum degree cannot see this coming:
it ranks vertices by `d` alone, never by how much of the clique already exists, so it will
sometimes pass over a high-degree vertex that would have cost nothing. That blindness is a large
part of why the heuristic is a heuristic.

The trouble is entirely in the implementation, and there are two problems.

**The graph grows, on balance.** Not monotonically: by the count above a degree-2 pivot shrinks it
and a degree-3 pivot is break-even at worst. But the pivots that matter are the ones with room to
fill, and summed over an elimination the additions win, so the ordering code ends up storing
inside itself the very fill the ordering exists to reduce.

Note what the ordering does and does not get for that price. It gets the total fill nearly free, and
by a route worth knowing: the degree of a pivot at the moment it is eliminated is exactly the
column count of `L` for that column, since the pivot's reachable set is precisely the set of rows
below the diagonal that fill in. So

```
nnz(L) = sum over pivots of (degree at elimination)  +  n
fill   = nnz(L) - nnz(tril A)
```

and any minimum-degree code can report `nnz(L)` by accumulating a number it already computes at
every step. This is the same quantity 2.5 obtains from the elimination forest; the forest is not
needed for the counts when one is performing the elimination anyway, it is needed to get them from
`A` and a permutation *without* performing it.

What the ordering cannot do is consult that total while choosing. The number only exists once the
whole sequence is fixed, and each greedy step has no view of the future, so the choice is made on
a local proxy and the total arrives too late to inform it.

It is worth being precise about what the ordering *does* hold, because it is more than the count.
The reachable set of a pivot at the moment it is eliminated is not merely the size of a column of
`L`, it is that column's nonzero pattern. So an ordering code is already performing a symbolic
factorization, one column at a time, and discarding each column as it goes. Symbolic factorization
(Section 3) exists as a separate phase not because the ordering lacks the information but because
producing it there would be the most expensive possible way to get it:

- **Allocation.** A symbolic factorization wants to size its output exactly, once, and fill a flat
  array. An ordering discovers each column's size only on reaching it, so emitting `L` from inside
  it means growable per-column storage. The forest breaks the deadlock by computing every column
  count *before* any pattern is formed, in nearly `O(nnz(A))` rather than `O(nnz(L))` (2.5), so the
  array can be sized up front and never grown.
- **Waste.** To find each pivot the ordering computes a reachable set for every live vertex and
  keeps one. The other `n - 1` per step are built and thrown away, so the great majority of its set
  work is not the answer.
- **Granularity.** Once supervariables are in play (5.5) the patterns are over principal variables,
  not individual ones, and would need expanding before they were `L`'s pattern.

The ordering, in short, can do symbolic factorization and does do it, in the wrong order and with
structures built for a different question. Separating the phases is about allocating once and
computing each thing where it is cheap, not about who knows what.

The way out is not to fight the arithmetic but to notice what is being stored. The edges an
elimination adds are not arbitrary: they complete a *clique* on the pivot's neighborhood, and a
clique is fully described by its vertex list. Every edge inside it is implicit. So the `d(d-1)/2`
edges need never be written down; keeping the list of `d` vertices says the same thing.

That observation cuts twice, and the second cut is the one that makes the technique bounded rather
than merely cheaper. If membership in the clique implies adjacency, then the edges *already
present* between two of its members are redundant too, and can be deleted from the explicit
adjacency. An elimination therefore adds nothing and removes something. This is why the
representation of 5.3 never needs more room than the original graph did.

**And the degrees must be recomputed.** After eliminating `p`, every neighbor's degree changes,
and computing a degree means taking a set union over the vertex's neighborhood. Done naively this
is the dominant cost, by a wide margin. Note that the pseudocode above hides two separate costs
under one line: finding the least-degree vertex, and knowing the degrees in the first place. A
linear scan for the minimum is `O(n)` per step and is the easier of the two to fix, by filing
vertices in buckets by degree so the pick is a walk up from the last known minimum (this is what
`mmdint` builds, and the buckets must be doubly linked, since a vertex whose degree changes has to
be pulled from the middle of one). But bucketing only helps if the degrees are *maintained*. As
long as each degree is recomputed on demand, the pick is doing a set union per live vertex and
dominates the elimination itself. Maintaining degrees incrementally is the real problem, and the
two codes answer it in opposite ways: AMD makes each update cheap (5.13), MMD makes updates rare
(5.11).

Everything intricate in a real minimum-degree code is one of those two problems being solved.

### 5.2 A worked example: the naive algorithm

Three graphs, used again in 5.4 so the two representations can be compared step for step.

```
   graph1                graph2                graph3
   a 4-cycle             uneven degrees        twelve vertices

   0---1                     0                 a path 0-1-2-...-11 plus eight
   |   |                    / \                extra edges: 0-3 0-8 1-6 1-8
   3---2                   1   2                2-5 5-9 6-10 7-8
                           |   |
   edges:                  3---4               degrees run 1 to 4, and the
   0-1 1-2 2-3 3-0          \ /                elimination order is not the
                             5                 identity
```

**graph1.** All four degrees are 2, so the tie goes to the lowest index.

```
start                                          storage 8
  graph: {0: [1,3], 1: [0,2], 2: [1,3], 3: [0,2]}
step 0: eliminate 0 (degree 2)    fill: 1-3              storage 6
  graph: {1: [2,3], 2: [1,3], 3: [1,2]}
step 1: eliminate 1 (degree 2)    fill: none             storage 2
  graph: {2: [3], 3: [2]}
step 2: eliminate 2 (degree 1)    fill: none             storage 0
  graph: {3: []}
step 3: eliminate 3 (degree 0)    fill: none             storage 0
  graph: {}
order = [0, 1, 2, 3]
fill = 1,  nnz(L) = 9 against nnz(tril A) = 8
```

One fill edge, `1-3`, the diagonal of the square. Watch vertex 1: it enters with neighbors
`[0, 2]` and leaves step 0 with `[2, 3]`. It lost the eliminated pivot and gained a vertex it was
never adjacent to. That gain is the fill, written into the graph. After it, `{1,2,3}` is a
triangle, so nothing further can fill and every later step reports `none`.

**graph2.** Six vertices, two fill edges.

```
start                                          storage 14
  graph: {0: [1,2], 1: [0,3], 2: [0,4], 3: [1,4,5], 4: [2,3,5], 5: [3,4]}
step 0: eliminate 0 (degree 2)    fill: 1-2              storage 12
  graph: {1: [2,3], 2: [1,4], 3: [1,4,5], 4: [2,3,5], 5: [3,4]}
step 1: eliminate 1 (degree 2)    fill: 2-3              storage 10
  graph: {2: [3,4], 3: [2,4,5], 4: [2,3,5], 5: [3,4]}
step 2: eliminate 2 (degree 2)    fill: none             storage 6
  graph: {3: [4,5], 4: [3,5], 5: [3,4]}
step 3: eliminate 3 (degree 2)    fill: none             storage 2
  graph: {4: [5], 5: [4]}
step 4: eliminate 4 (degree 1)    fill: none             storage 0
  graph: {5: []}
step 5: eliminate 5 (degree 0)    fill: none             storage 0
  graph: {}
order = [0, 1, 2, 3, 4, 5]
fill = 2,  nnz(L) = 15 against nnz(tril A) = 13
```

Note the degrees: 2, 2, 2, 2, 1, 0. Vertex 3 began at degree 3 and is eliminated at degree 2,
because two of its neighbors went first and the fill it received did not make up the difference.
The degree a vertex is eliminated at is not the degree it started with, which is the whole reason
degrees have to be recomputed as the algorithm runs.

**graph3.** Twelve vertices, eighteen edges, and the first example whose order is not the
identity.

```
start                                          storage 36
  graph: {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,4], 4: [3,5],
          5: [2,4,6,9], 6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9],
          9: [5,8,10], 10: [6,9,11], 11: [10]}
step 0: eliminate 11 (degree 1)   fill: none             storage 34
  graph: {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,4], 4: [3,5],
          5: [2,4,6,9], 6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9],
          9: [5,8,10], 10: [6,9]}
step 1: eliminate 4 (degree 2)    fill: 3-5              storage 32
  graph: {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,5], 5: [2,3,6,9],
          6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9], 9: [5,8,10], 10: [6,9]}
step 2: eliminate 7 (degree 2)    fill: 6-8              storage 30
  graph: {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,5], 5: [2,3,6,9],
          6: [1,5,8,10], 8: [0,1,6,9], 9: [5,8,10], 10: [6,9]}
step 3: eliminate 10 (degree 2)   fill: 6-9              storage 28
  graph: {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,5], 5: [2,3,6,9],
          6: [1,5,8,9], 8: [0,1,6,9], 9: [5,6,8]}
step 4: eliminate 0 (degree 3)    fill: 1-3, 3-8         storage 26
  graph: {1: [2,3,6,8], 2: [1,3,5], 3: [1,2,5,8], 5: [2,3,6,9],
          6: [1,5,8,9], 8: [1,3,6,9], 9: [5,6,8]}
step 5: eliminate 2 (degree 3)    fill: 1-5              storage 22
  graph: {1: [3,5,6,8], 3: [1,5,8], 5: [1,3,6,9], 6: [1,5,8,9],
          8: [1,3,6,9], 9: [5,6,8]}
step 6: eliminate 3 (degree 3)    fill: 5-8              storage 18
  graph: {1: [5,6,8], 5: [1,6,8,9], 6: [1,5,8,9], 8: [1,5,6,9],
          9: [5,6,8]}
step 7: eliminate 1 (degree 3)    fill: none             storage 12
  graph: {5: [6,8,9], 6: [5,8,9], 8: [5,6,9], 9: [5,6,8]}
step 8: eliminate 5 (degree 3)    fill: none             storage 6
  graph: {6: [8,9], 8: [6,9], 9: [6,8]}
step 9: eliminate 6 (degree 2)    fill: none             storage 2
  graph: {8: [9], 9: [8]}
step 10: eliminate 8 (degree 1)   fill: none             storage 0
  graph: {9: []}
step 11: eliminate 9 (degree 0)   fill: none             storage 0
  graph: {}
order = [11, 4, 7, 10, 0, 2, 3, 1, 5, 6, 8, 9]
fill = 7,  nnz(L) = 37 against nnz(tril A) = 30
```

Minimum degree begins at vertex 11, the only degree-1 vertex, and works inward; the four
degree-4 vertices (1, 5, 6, 8) are all deferred to the end, which is exactly the behavior the
heuristic is for. The degrees at elimination run 1, 2, 2, 2, 3, 3, 3, 3, 3, 2, 1, 0, never
exceeding 3 even though four vertices started at 4.

Two things are visible here that the smaller graphs cannot show. Fill arrives in bursts rather
than singly: step 4 creates two edges at once, `1-3` and `3-8`, because vertex 0 had three
neighbors and two of the three pairs were not yet adjacent. And the graph becomes dense at the
end without any further fill being created, since by step 7 the surviving `{5,6,8,9}` is already
a complete graph, so the last five eliminations report `none`.

The storage line is the clearest statement of the problem 5.3 solves. It falls by 2 per step at
first, each elimination removing a low-degree vertex and adding little, then the decrease slows
through the middle as fill offsets the removals, and only collapses once the graph is nearly
gone. On a graph this small the additions never overtake the deletions, so the peak is the
starting value of 36. On a larger and denser problem they do overtake and the peak sits well
above the start: that is the growth 5.1 warned about, and the reason a real code cannot store the
graph this way.


**graph4.** Eight vertices, fourteen edges, and denser than the others: four fill edges
appear in the first two steps and none after.

```
start                                                storage 28
  graph: {0: [2,3,4,7], 1: [3,4,6,7], 2: [0,3,5], 3: [0,1,2,6,7],
          4: [0,1,5], 5: [2,4,6], 6: [1,3,5], 7: [0,1,3]}
step 0: eliminate 2 (degree 3)    fill: 0-5, 3-5     storage 26
  graph: {0: [3,4,5,7], 1: [3,4,6,7], 3: [0,1,5,6,7], 4: [0,1,5],
          5: [0,3,4,6], 6: [1,3,5], 7: [0,1,3]}
step 1: eliminate 4 (degree 3)    fill: 0-1, 1-5     storage 24
  graph: {0: [1,3,5,7], 1: [0,3,5,6,7], 3: [0,1,5,6,7], 5: [0,1,3,6],
          6: [1,3,5], 7: [0,1,3]}
step 2: eliminate 6 (degree 3)    fill: none         storage 18
  graph: {0: [1,3,5,7], 1: [0,3,5,7], 3: [0,1,5,7], 5: [0,1,3],
          7: [0,1,3]}
step 3: eliminate 5 (degree 3)    fill: none         storage 12
  graph: {0: [1,3,7], 1: [0,3,7], 3: [0,1,7], 7: [0,1,3]}
step 4: eliminate 0 (degree 3)    fill: none         storage 6
  graph: {1: [3,7], 3: [1,7], 7: [1,3]}
step 5: eliminate 1 (degree 2)    fill: none         storage 2
  graph: {3: [7], 7: [3]}
step 6: eliminate 3 (degree 1)    fill: none         storage 0
  graph: {7: []}
step 7: eliminate 7 (degree 0)    fill: none         storage 0
order = [2, 4, 6, 5, 0, 1, 3, 7]
fill = 4
```

### 5.3 The quotient graph

Do not add the clique's edges. **Store the clique as one object.**

An eliminated vertex becomes an **element**, and its pattern `L_p` is the set of vertices it made
mutually adjacent. A live vertex `i` is then adjacent to a mixture:

```
A_i     the original neighbors of i that have not been eliminated  (variables)
E_i     the elements i belongs to                                   (cliques)
```

and its true neighborhood is

```
adj(i) = A_i  union  ( union of L_e over e in E_i )   minus  i itself
```

which is never formed. The clique that eliminating `p` created is represented by the single fact
that its members all list `p` in their `E`.

**Elements are absorbed.** When `p` is eliminated, every element adjacent to `p` is consumed into
`L_p` and dies: whatever those cliques joined, the new one joins too. So the number of elements
stays bounded, and the representation does not grow without limit. This is the whole reason the
quotient graph is affordable.

**A note on the name.** The literature calls these objects *elements*, and writes `E_i` for the
list above. The word is inherited from finite-element assembly, where the cliques really were
mesh elements, and it has stuck ever since. They are cliques, and the rest of this section says
so where it helps: `A_i` for the explicit neighbors, `C_i` for the cliques containing `i`, and
`L_c` for the members of clique `c`. When reading AMD or MMD, read `element` for `clique` and
`E_i` for `C_i`; nothing else changes.

Setting the naive version of 5.1 and 5.2 beside this one, the loop is untouched and only the two
operations underneath it differ:

```
minimum_degree_quotient(G) -> elimination order:

    for each vertex i:
        A_i = the neighbors of i in G          # explicit, as given
        C_i = {}                                # no cliques yet

    while a live vertex remains:

        # DEGREE, on demand: unite the explicit neighbors with the members of
        # every clique i belongs to. Never stored, formed only when asked.
        neighbors(i) = A_i  union  ( union of L_c over c in C_i )   minus i

        p = a live vertex minimizing |neighbors(p)|
        L = neighbors(p)

        # ELIMINATE. The clique on L is recorded once, not written out as edges.
        for each c in C_p:
            delete clique c                     # ABSORPTION: subsumed by the new one
        L_p = L;  p is now a clique

        for each i in L:
            A_i = A_i \ L                        # PRUNING: the new clique implies these
            A_i = A_i \ {p}
            C_i = ( C_i \ absorbed ) union {p}

        append p to the order
```

Three things distinguish it from 5.1. The clique on `L` is never expanded into edges, so no fill
is stored. The absorbed cliques are deleted, which bounds how many can accumulate. And `A_i` is
pruned of everything the new clique now implies, which is what makes each explicit list shrink
monotonically. The result is identical to 5.1's, pivot for pivot.

### 5.4 A worked example: the quotient graph

Take the smallest graph that fills, a 4-cycle:

```
   0---1        edges: 0-1 1-2 2-3 3-0
   |   |
   3---2
```

At the start there are no cliques and every neighbor is explicit, so the quotient graph is
just the graph:

```
A       = {0: [1,3], 1: [0,2], 2: [1,3], 3: [0,2]}
C       = {0: [],    1: [],    2: [],    3: []}
cliques = {}
```

All four degrees are 2, so the tie goes to the lowest index and `p = 0`. Its neighborhood needs
no clique expansion, since `C_0` is empty:

```
neighbors(0) = A_0                     = {1, 3}
```

That set is stored once as clique `c0` and vertices 1 and 3 record their membership:

```
A       = {0: [], 1: [2],    2: [1,3], 3: [2]}
C       = {0: [], 1: [c0],   2: [],    3: [c0]}
cliques = {c0: [1,3]}
neighbors = [1,3], absorbed = none, pruned = none
```

Compare with 5.1, which at this point would have added the edge `1-3` to the graph. Here that
edge exists only as the fact that 1 and 3 both list `c0`. Note also that vertex 2, which was not
reached, is untouched.

Next, `p = 1`. Now the clique does work: `A_1` alone would say the degree is 1, but expanding
`c0` supplies the implicit neighbor 3.

```
neighbors(1) = A_1 union L_c0 minus 1  = {2} union {1,3} minus 1  =  {2, 3}
```

Eliminating 1 does three things at once. The new clique `c1 = {2,3}` is created; `c0` is absorbed,
because 1 was one of its members; and the *real* edge `2-3` is pruned, because `c1` now implies
it. That last step is why both explicit lists go empty:

```
A       = {0: [], 1: [], 2: [],     3: []}
C       = {0: [], 1: [], 2: [c1],   3: [c1]}
cliques = {c1: [2,3]}
neighbors = [2,3], absorbed = c0, pruned = 2-3
```

From here the graph is carried entirely by cliques. For `p = 2` the explicit side contributes
nothing at all:

```
neighbors(2) = {} union L_c1 minus 2   = {2,3} minus 2  =  {3}

A       = {0: [], 1: [], 2: [], 3: []}
C       = {0: [], 1: [], 2: [], 3: [c2]}
cliques = {c2: [3]}
neighbors = [3], absorbed = c1, pruned = none
```

And the last vertex has nothing left to join:

```
neighbors(3) = {} union L_c2 minus 3   = {3} minus 3  =  {}

A       = {0: [], 1: [], 2: [], 3: []}
C       = {0: [], 1: [], 2: [], 3: []}
cliques = {c3: []}
neighbors = [], absorbed = c2, pruned = none
```

Three things are worth reading off the trace. The `neighbors` line gives 2, 2, 1, 0, which by
5.1's identity are the column counts of `L`; they sum with the diagonal to `nnz(L) = 9` against
`nnz(tril A) = 8`, so the fill is the single edge `1-3`, exactly what the naive version stored.
The explicit lists hold 8 entries at the start, 4 after the first elimination and 0 after the
second, never rising. And `cliques` never holds more than one entry, because each elimination
absorbs its predecessor.

**A is the same object as the naive graph.** This is worth stating plainly, because the two
sections can look like different data structures and they are not. `A` starts as the adjacency
of 5.2, entry for entry, and both versions end with it empty. The difference is entirely in
what happens in between. The naive version both adds and deletes, so its graph can swell before
it drains. The quotient version only ever deletes, by three separate routes: no fill is written,
the pivot's own edges go with it, and existing edges are pruned once a clique implies them.

The third route is the one with no counterpart in 5.2, and step 1 above isolates it. In the
naive trace, eliminating 1 leaves vertex 2 holding `[3]`. In the quotient trace the same step
leaves `A_2` empty, because `c1 = {2,3}` now carries that edge. Nothing was lost; the information
moved from `A` to `C`. That migration, not merely the absence of fill, is why the explicit side
drains faster here than there.

**A clique owns its pivot implicitly.** `c1` is named 1, stores `{2,3}`, and stands for the
clique on `{1,2,3}`. The name carries the pivot and the list carries the rest, so the pivot is
never stored twice. This is not a bookkeeping trick but the `+ n` of 5.1's identity: the stored
list `L_c` is the below-diagonal part of column `c` of `L`, and the implicit pivot is that
column's diagonal entry. It is also what lets cliques and vertices share one name space, since a
clique is always named for the vertex whose elimination created it.

**The last clique of a component is empty.** A clique is created empty exactly when its pivot had
no live neighbors left, and an empty clique appears in nobody's `C`, so nothing can ever absorb
it. Each connected component therefore terminates in one empty clique, named for the last pivot
eliminated in that component. In the trace above that is `c3`.

Two notes on reading that, one about labels and one about structure. Cliques are named by *old*
index, because the new ordering does not exist yet: it is the thing being computed, and the only
names available while running are A's. So the terminal clique of a connected graph is `c` of
`order[n-1]` in old labels, which is the same object as position `n-1` in the new ordering. On
`graph3` the survivor is `c9`, and 9 is the twelfth vertex eliminated, so in the permuted
numbering it is exactly `n-1`. The two statements agree; they are one fact in two coordinate
systems, and the code can only use the first.

The structural note is the disconnected case, which is not a labeling artifact. Two triangles
plus an isolated vertex end with three empty cliques, sitting at new positions 0, 3 and 6: the
last vertex of each component. Only the final one is at `n-1`. One terminal clique per component
is the general rule, and the single clique of a connected graph is its special case.

**Reducible matrices, and what the empty cliques are counting.** A symmetric `A` whose graph is
disconnected is structurally *reducible*: some permutation makes it block diagonal, one block per
component. That is exactly the case 2.1 raises when it explains why this document says
elimination **forest** and not tree, since each component contributes its own tree and only an
irreducible `A` gives a single one.

The two observations are the same observation. A clique ends up empty precisely when its pivot
has no below-diagonal nonzero, and `Struct(j)` empty is 2.1's definition of a **root**. So the
terminal empty cliques and the roots of the elimination forest are the same objects, counted by
the same rule: one per connected component. On `graph3` there is a single root, vertex 9, and the
forest is a tree. On two triangles plus an isolated vertex there are three roots, 2, 5 and 6, one
per block.

Nothing in the algorithm needs to detect this. Components never interact, because no fill can
cross between them and no clique can span them, so a reducible `A` is factorized as independent
subproblems without the code being told they are independent. The forest simply comes out with
more than one root.

The independence is stronger than it may look, and 5.11 turns it into a concrete result. Since no
fill crosses between components, the total fill is the sum of the components' fills, and each
component's fill depends only on the relative order *within* it. The order in which we interleave
between components is therefore free: any interleaving gives the same fill. That is what makes one
particular batching strategy provably lossless where the general one is not.

Two more graphs, in the same order as 5.2 so the traces can be set side by side.

**graph2.** Same order, same degrees and same fill as the naive run, which is the point.

```
start                                          storage 14
  A       = {0: [1,2], 1: [0,3], 2: [0,4], 3: [1,4,5], 4: [2,3,5],
             5: [3,4]}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
  cliques = {}
  neighbors = none, absorbed = none, pruned = none
step 0: eliminate 0 (degree 2)                 storage 12
  A       = {0: [], 1: [3], 2: [4], 3: [1,4,5], 4: [2,3,5], 5: [3,4]}
  C       = {0: [], 1: [c0], 2: [c0], 3: [], 4: [], 5: []}
  cliques = {c0: [1,2]}
  neighbors = [1,2], absorbed = none, pruned = none
step 1: eliminate 1 (degree 2)                 storage 10
  A       = {0: [], 1: [], 2: [4], 3: [4,5], 4: [2,3,5], 5: [3,4]}
  C       = {0: [], 1: [], 2: [c1], 3: [c1], 4: [], 5: []}
  cliques = {c1: [2,3]}
  neighbors = [2,3], absorbed = c0, pruned = none
step 2: eliminate 2 (degree 2)                 storage 6
  A       = {0: [], 1: [], 2: [], 3: [5], 4: [5], 5: [3,4]}
  C       = {0: [], 1: [], 2: [], 3: [c2], 4: [c2], 5: []}
  cliques = {c2: [3,4]}
  neighbors = [3,4], absorbed = c1, pruned = 3-4
step 3: eliminate 3 (degree 2)                 storage 2
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [c3], 5: [c3]}
  cliques = {c3: [4,5]}
  neighbors = [4,5], absorbed = c2, pruned = 4-5
step 4: eliminate 4 (degree 1)                 storage 1
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [c4]}
  cliques = {c4: [5]}
  neighbors = [5], absorbed = c3, pruned = none
step 5: eliminate 5 (degree 0)                 storage 0
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
  cliques = {c5: []}
  neighbors = [], absorbed = c4, pruned = none
order = [0, 1, 2, 3, 4, 5]
fill = 2,  nnz(L) = 15 against nnz(tril A) = 13
```

Pruning fires twice, at steps 2 and 3, each time deleting a genuine original edge that a new
clique had made redundant: `3-4` and then `4-5`. Storage runs 14, 12, 10, 6, 2, 1, 0 against the
naive 14, 12, 10, 6, 2, 0, 0, which is nearly identical, because every clique here has exactly
two members and a two-member clique costs precisely what the edge it replaces cost. The quotient
graph's advantage is quadratic in clique size, so a graph this sparse cannot show it.

Follow one vertex to see the migration. Vertex 4 starts with `A_4 = [2,3,5]`, three explicit
neighbors and no cliques. By step 2 it holds `A_4 = [5]` and `C_4 = [c2]`: two of its neighbors
have moved into a clique and only one is still spelled out. By step 3 its explicit list is empty
and everything it knows is carried by `c3`. No adjacency was lost at any point.

**graph3.** Twelve vertices, and the first graph whose cliques exceed two members.

```
start                                          storage 36
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,4], 4: [3,5],
             5: [2,4,6,9], 6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9],
             9: [5,8,10], 10: [6,9,11], 11: [10]}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  cliques = {}
  neighbors = none, absorbed = none, pruned = none
step 0: eliminate 11 (degree 1)                storage 35
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,4], 4: [3,5],
             5: [2,4,6,9], 6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9],
             9: [5,8,10], 10: [6,9], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [c11], 11: []}
  cliques = {c11: [10]}
  neighbors = [10], absorbed = none, pruned = none
step 1: eliminate 4 (degree 2)                 storage 33
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2], 4: [],
             5: [2,6,9], 6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9],
             9: [5,8,10], 10: [6,9], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [c4], 4: [], 5: [c4], 6: [], 7: [],
             8: [], 9: [], 10: [c11], 11: []}
  cliques = {c4: [3,5], c11: [10]}
  neighbors = [3,5], absorbed = none, pruned = none
step 2: eliminate 7 (degree 2)                 storage 31
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2], 4: [],
             5: [2,6,9], 6: [1,5,10], 7: [], 8: [0,1,9], 9: [5,8,10],
             10: [6,9], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [c4], 4: [], 5: [c4], 6: [c7], 7: [],
             8: [c7], 9: [], 10: [c11], 11: []}
  cliques = {c4: [3,5], c7: [6,8], c11: [10]}
  neighbors = [6,8], absorbed = none, pruned = none
step 3: eliminate 10 (degree 2)                storage 28
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2], 4: [],
             5: [2,6,9], 6: [1,5], 7: [], 8: [0,1,9], 9: [5,8], 10: [],
             11: []}
  C       = {0: [], 1: [], 2: [], 3: [c4], 4: [], 5: [c4], 6: [c7,c10],
             7: [], 8: [c7], 9: [c10], 10: [], 11: []}
  cliques = {c4: [3,5], c7: [6,8], c10: [6,9]}
  neighbors = [6,9], absorbed = c11, pruned = none
step 4: eliminate 0 (degree 3)                 storage 23
  A       = {0: [], 1: [2,6], 2: [1,3,5], 3: [2], 4: [], 5: [2,6,9],
             6: [1,5], 7: [], 8: [9], 9: [5,8], 10: [], 11: []}
  C       = {0: [], 1: [c0], 2: [], 3: [c0,c4], 4: [], 5: [c4],
             6: [c7,c10], 7: [], 8: [c0,c7], 9: [c10], 10: [], 11: []}
  cliques = {c0: [1,3,8], c4: [3,5], c7: [6,8], c10: [6,9]}
  neighbors = [1,3,8], absorbed = none, pruned = 1-8
step 5: eliminate 2 (degree 3)                 storage 20
  A       = {0: [], 1: [6], 2: [], 3: [], 4: [], 5: [6,9], 6: [1,5],
             7: [], 8: [9], 9: [5,8], 10: [], 11: []}
  C       = {0: [], 1: [c0,c2], 2: [], 3: [c0,c2,c4], 4: [], 5: [c2,c4],
             6: [c7,c10], 7: [], 8: [c0,c7], 9: [c10], 10: [], 11: []}
  cliques = {c0: [1,3,8], c2: [1,3,5], c4: [3,5], c7: [6,8], c10: [6,9]}
  neighbors = [1,3,5], absorbed = none, pruned = none
step 6: eliminate 3 (degree 3)                 storage 15
  A       = {0: [], 1: [6], 2: [], 3: [], 4: [], 5: [6,9], 6: [1,5],
             7: [], 8: [9], 9: [5,8], 10: [], 11: []}
  C       = {0: [], 1: [c3], 2: [], 3: [], 4: [], 5: [c3], 6: [c7,c10],
             7: [], 8: [c3,c7], 9: [c10], 10: [], 11: []}
  cliques = {c3: [1,5,8], c7: [6,8], c10: [6,9]}
  neighbors = [1,5,8], absorbed = c0, c2, c4, pruned = none
step 7: eliminate 1 (degree 3)                 storage 11
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [9], 6: [], 7: [],
             8: [9], 9: [5,8], 10: [], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [c1], 6: [c1,c7,c10],
             7: [], 8: [c1,c7], 9: [c10], 10: [], 11: []}
  cliques = {c1: [5,6,8], c7: [6,8], c10: [6,9]}
  neighbors = [5,6,8], absorbed = c3, pruned = 5-6
step 8: eliminate 5 (degree 3)                 storage 7
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [c5,c7,c10],
             7: [], 8: [c5,c7], 9: [c5,c10], 10: [], 11: []}
  cliques = {c5: [6,8,9], c7: [6,8], c10: [6,9]}
  neighbors = [6,8,9], absorbed = c1, pruned = 8-9
step 9: eliminate 6 (degree 2)                 storage 2
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [c6], 9: [c6], 10: [], 11: []}
  cliques = {c6: [8,9]}
  neighbors = [8,9], absorbed = c5, c7, c10, pruned = none
step 10: eliminate 8 (degree 1)                storage 1
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [c8], 10: [], 11: []}
  cliques = {c8: [9]}
  neighbors = [9], absorbed = c6, pruned = none
step 11: eliminate 9 (degree 0)                storage 0
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  cliques = {c9: []}
  neighbors = [], absorbed = c8, pruned = none
order = [11, 4, 7, 10, 0, 2, 3, 1, 5, 6, 8, 9]
fill = 7,  nnz(L) = 37 against nnz(tril A) = 30
```

The same order and the same fill as 5.2, as they must be. What differs is the shape of the
storage curve, and it is worth setting the two side by side:

```
naive     36, 34, 32, 30, 28, 26, 22, 18, 12, 6, 2, 0, 0      peak 36
quotient  36, 35, 33, 31, 28, 23, 20, 15, 11, 7, 2, 1, 0
```

They cross twice. The quotient graph is *behind* over the first four steps, since each
elimination costs it a clique entry where the naive version simply drops edges and stores no
replacement. From step 5 it pulls ahead as pruning and absorption compound, reaching 15 against
18 and 11 against 12, then falls behind again at the very end for the trivial reason that it
still holds one empty clique when the naive graph holds nothing at all.

That crossover is characteristic and worth not overreading in either direction. The quotient
graph is not uniformly smaller on small problems; on this one it wins by three entries at its
best. What it is, is *bounded*: `A` never grows, so the total cannot exceed what the original
graph cost, whatever the fill turns out to be. The naive version has no such guarantee, and on a
problem where the fill is the dominant term the difference stops being three entries and becomes
the difference between fitting in memory and not.

**Why these examples cannot show the advantage.** On all three graphs above the naive storage
falls monotonically; it never once exceeds its starting value, so the bound the quotient graph
offers is never tested. That is not because the graphs are small. It is because they are the
wrong *shape*, and the reason is the count from 5.1: a pivot of degree `d` deletes `2d` endpoints
and adds at most `d(d-1)`, so storage can only grow when `d >= 4`. Minimum degree exists
precisely to avoid such pivots, and it will keep avoiding them as long as the graph offers
anything cheaper. A 2D grid always does, since its corners sit at degree 2 and its edges at 3, so
even an 8x8 grid never exceeds its start. `graph3`, with a degree-1 vertex and several degree-2s,
is comfortably in the same easy regime. Enlarging it would not help.

Adding a dimension does. A 3D grid has interior degree 6 and no cheap vertices to fall back on,
which forces minimum degree into the growth regime:

```
3D grid, 5x5x5:  n = 125, 300 edges, fill = 941

                 start   peak
naive             600    1132     1.89x, rising to a peak near the middle
quotient          600     600     1.00x, monotone down: 576, 544, 508, ...
```

The naive graph nearly doubles before it drains; the quotient graph never rises at all, so its
starting size is also its ceiling. That gap is the whole point of 5.3, and it is invisible on any
example small and shapely enough to trace by hand.


**graph4.** The densest of the four, and the clearest case for the quotient graph: storage
falls from 28 to 0 without ever rising, where the naive run of 5.2 had to add fill edges to a
representation that only grows.

```
start                                          storage 28
  A       = {0: [2,3,4,7], 1: [3,4,6,7], 2: [0,3,5], 3: [0,1,2,6,7],
             4: [0,1,5], 5: [2,4,6], 6: [1,3,5], 7: [0,1,3]}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  cliques = {}
  neighbors = none, absorbed = none, pruned = none
step 0: eliminate 2 (degree 3)                 storage 23
  A       = {0: [4,7], 1: [3,4,6,7], 2: [], 3: [1,6,7], 4: [0,1,5],
             5: [4,6], 6: [1,3,5], 7: [0,1,3]}
  C       = {0: [c2], 1: [], 2: [], 3: [c2], 4: [], 5: [c2], 6: [], 7: []}
  cliques = {c2: [0,3,5]}
  neighbors = {0,3,5}, absorbed = none, pruned = 0-3
step 1: eliminate 4 (degree 3)                 storage 20
  A       = {0: [7], 1: [3,6,7], 2: [], 3: [1,6,7], 4: [], 5: [6],
             6: [1,3,5], 7: [0,1,3]}
  C       = {0: [c2,c4], 1: [c4], 2: [], 3: [c2], 4: [], 5: [c2,c4],
             6: [], 7: []}
  cliques = {c2: [0,3,5], c4: [0,1,5]}
  neighbors = {0,1,5}, absorbed = none, pruned = none
step 2: eliminate 6 (degree 3)                 storage 15
  A       = {0: [7], 1: [7], 2: [], 3: [7], 4: [], 5: [], 6: [],
             7: [0,1,3]}
  C       = {0: [c2,c4], 1: [c4,c6], 2: [], 3: [c2,c6], 4: [],
             5: [c2,c4,c6], 6: [], 7: []}
  cliques = {c2: [0,3,5], c4: [0,1,5], c6: [1,3,5]}
  neighbors = {1,3,5}, absorbed = none, pruned = 1-3
step 3: eliminate 5 (degree 3)                 storage 9
  A       = {0: [7], 1: [7], 2: [], 3: [7], 4: [], 5: [], 6: [],
             7: [0,1,3]}
  C       = {0: [c5], 1: [c5], 2: [], 3: [c5], 4: [], 5: [], 6: [], 7: []}
  cliques = {c5: [0,1,3]}
  neighbors = {0,1,3}, absorbed = c2, c4, c6, pruned = none
step 4: eliminate 0 (degree 3)                 storage 3
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  C       = {0: [], 1: [c0], 2: [], 3: [c0], 4: [], 5: [], 6: [], 7: [c0]}
  cliques = {c0: [1,3,7]}
  neighbors = {1,3,7}, absorbed = c5, pruned = 1-7, 3-7
step 5: eliminate 1 (degree 2)                 storage 2
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  C       = {0: [], 1: [], 2: [], 3: [c1], 4: [], 5: [], 6: [], 7: [c1]}
  cliques = {c1: [3,7]}
  neighbors = {3,7}, absorbed = c0, pruned = none
step 6: eliminate 3 (degree 1)                 storage 1
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [c3]}
  cliques = {c3: [7]}
  neighbors = {7}, absorbed = c1, pruned = none
step 7: eliminate 7 (degree 0)                 storage 0
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  cliques = {c7: []}
  neighbors = none, absorbed = c3, pruned = none
order = [2, 4, 6, 5, 0, 1, 3, 7]
```

### 5.5 Supervariables

Two vertices with **identical neighborhoods** will have equal degree at every step, will be
chosen at the same moment, and will fill in identically. They can be merged and eliminated
together.

They are detected by **hashing** the pattern (`A_i` and `E_i` together) and then comparing
candidates within a bucket, since a hash collision is not a match. A merged vertex is a
**supervariable**, and it carries a weight: how many original vertices it stands for.

This is not a marginal optimization. Real matrices are full of indistinguishable vertices,
especially after a few eliminations have merged their neighborhoods, and supervariables are
often the difference between a usable code and an unusable one.

**Mass elimination** is the same idea one step further. After forming the clique `L`, if some
vertex `i` in `L` has nothing left outside `L`, it is adjacent only to the new clique, and can be
eliminated immediately at no cost. Its degree, computed below, comes out zero and says so.

Both mechanisms need the same two additions to 5.3, and nothing else. Each supervariable carries
a weight, and the degree becomes a weighted count:

```
minimum_degree_supervariables(G) -> elimination order:

    for each vertex i:
        A_i = the neighbors of i in G
        C_i = {}
        w_i = 1                                 # original vertices represented

    while a live supervariable remains:

        # DEGREE, now WEIGHTED: a neighbor of weight 3 contributes 3, because it
        # stands for three original vertices and each will fill in.
        degree(i) = sum of w_j over j in ( A_i union ( union of L_c over c in C_i ) ), j != i

        p = a live supervariable minimizing degree(p)
        L = the reachable set of p

        for each c in C_p:  delete clique c      # absorption, as in 5.3
        L_p = L
        for each i in L:
            A_i = A_i \ L  \ {p}                  # pruning, as in 5.3
            C_i = ( C_i \ absorbed ) union {p}

        # MASS ELIMINATION. Everything i can still see lies inside the new
        # clique, so eliminating it next would create no fill. Merge it now.
        for each i in L:
            if A_i is empty and C_i = {p}:
                w_p = w_p + w_i
                delete i from the graph and from every clique
                record i as a member of p

        append p to the pivot sequence

    expand: each pivot p contributes its w_p members, consecutively
```

The merge test is worth reading twice, because its cheapness is the whole point. `A_i` empty and
`C_i = {p}` says that after the pruning just performed, `i` has no explicit neighbor left and no
other clique to reach through, so its entire neighborhood is inside `L`. In the vendored MMD this
is a single integer comparison, `if (nq <= 0)`, on the count left after compacting `i`'s list.
Pruning (5.3) is what makes it that cheap: without it the test would be a containment check
against `L` rather than a test for emptiness.

**Where this changes the permutation, but not the fill.** Every layer up to here was a change of
representation, and 5.4 verified that the orderings came out identical. This one is different.
Merged vertices are numbered consecutively whether or not the degree would have chosen them in
that order, so the permutation genuinely can differ from what 5.3 produces, and the final
expansion pass exists for exactly that reason.

The constraint turns out to cost nothing. Across fourteen test graphs (grids in two and three
dimensions, and random graphs from thirty to seventy vertices) the permutation differed from 5.3's
on nine of them and `nnz(L)` was identical on all fourteen. The merge is not merely sound, it
appears to be free: a merged vertex creates no fill when eliminated next, so promoting it ahead of
whatever the degree would have picked cannot cost anything. Measured rather than proved, and on
prototypes rather than production code, but the result is clean enough to state.

### 5.6 A worked example: supervariables

The same three graphs again, so all four representations can be compared on identical input.

**graph1.** Four vertices, but only two pivots.

```
start                                          storage 8
  A       = {0: [1,3], 1: [0,2], 2: [1,3], 3: [0,2]}
  C       = {0: [], 1: [], 2: [], 3: []}
  cliques = {}
  weights = {0: 1, 1: 1, 2: 1, 3: 1}
  neighbors = none, absorbed = none, pruned = none, merged = none
step 0: eliminate 0 (degree 2, weight 1 -> 1)  storage 6
  A       = {0: [], 1: [2], 2: [1,3], 3: [2]}
  C       = {0: [], 1: [c0], 2: [], 3: [c0]}
  cliques = {c0: [1,3]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1}
  neighbors = [1,3], absorbed = none, pruned = none,
  merged = none
step 1: eliminate 1 (degree 2, weight 1 -> 3)  storage 0
  A       = {0: [], 1: [], 2: [], 3: []}
  C       = {0: [], 1: [], 2: [], 3: []}
  cliques = {c1: []}
  weights = {0: 1, 1: 3, 2: 0, 3: 0}
  neighbors = [2,3], absorbed = c0, pruned = 2-3,
  merged = 2, 3
supervariable pivots = [0, 1]   (2 of 4)
order = [0, 1, 2, 3]
fill = 1,  nnz(L) = 9 against nnz(tril A) = 8
```

Step 1 is the whole idea in one step. Eliminating 1 forms the clique `{2,3}`, and the pruning
that follows empties both `A_2` and `A_3` while leaving each with `C = {c1}` and nothing else.
Both therefore satisfy the merge test, both are absorbed, and vertex 1 becomes a supervariable of
weight 3 standing for `{1,2,3}`. The clique `c1` is left empty, since all of its members have
joined the supervariable that owns it.

Two pivots instead of four, and the same ordering and the same `nnz(L) = 9` as 5.4. The saving is
not in the answer but in the number of times the degree scan has to run.

**graph2.** Four pivots instead of six.

```
start                                          storage 14
  A       = {0: [1,2], 1: [0,3], 2: [0,4], 3: [1,4,5], 4: [2,3,5],
             5: [3,4]}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
  cliques = {}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1}
  neighbors = none, absorbed = none, pruned = none, merged = none
step 0: eliminate 0 (degree 2, weight 1 -> 1)  storage 12
  A       = {0: [], 1: [3], 2: [4], 3: [1,4,5], 4: [2,3,5], 5: [3,4]}
  C       = {0: [], 1: [c0], 2: [c0], 3: [], 4: [], 5: []}
  cliques = {c0: [1,2]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1}
  neighbors = [1,2], absorbed = none, pruned = none,
  merged = none
step 1: eliminate 1 (degree 2, weight 1 -> 1)  storage 10
  A       = {0: [], 1: [], 2: [4], 3: [4,5], 4: [2,3,5], 5: [3,4]}
  C       = {0: [], 1: [], 2: [c1], 3: [c1], 4: [], 5: []}
  cliques = {c1: [2,3]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1}
  neighbors = [2,3], absorbed = c0, pruned = none,
  merged = none
step 2: eliminate 2 (degree 2, weight 1 -> 1)  storage 6
  A       = {0: [], 1: [], 2: [], 3: [5], 4: [5], 5: [3,4]}
  C       = {0: [], 1: [], 2: [], 3: [c2], 4: [c2], 5: []}
  cliques = {c2: [3,4]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1}
  neighbors = [3,4], absorbed = c1, pruned = 3-4,
  merged = none
step 3: eliminate 3 (degree 2, weight 1 -> 3)  storage 0
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
  cliques = {c3: []}
  weights = {0: 1, 1: 1, 2: 1, 3: 3, 4: 0, 5: 0}
  neighbors = [4,5], absorbed = c2, pruned = 4-5,
  merged = 4, 5
supervariable pivots = [0, 1, 2, 3]   (4 of 6)
order = [0, 1, 2, 3, 4, 5]
fill = 2,  nnz(L) = 15 against nnz(tril A) = 13
```

**graph3.** Ten pivots instead of twelve.

```
start                                          storage 36
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,4], 4: [3,5],
             5: [2,4,6,9], 6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9],
             9: [5,8,10], 10: [6,9,11], 11: [10]}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  cliques = {}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = none, absorbed = none, pruned = none, merged = none
step 0: eliminate 11 (degree 1, weight 1 -> 1) storage 35
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2,4], 4: [3,5],
             5: [2,4,6,9], 6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9],
             9: [5,8,10], 10: [6,9], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [c11], 11: []}
  cliques = {c11: [10]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [10], absorbed = none, pruned = none,
  merged = none
step 1: eliminate 4 (degree 2, weight 1 -> 1)  storage 33
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2], 4: [],
             5: [2,6,9], 6: [1,5,7,10], 7: [6,8], 8: [0,1,7,9],
             9: [5,8,10], 10: [6,9], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [c4], 4: [], 5: [c4], 6: [], 7: [],
             8: [], 9: [], 10: [c11], 11: []}
  cliques = {c4: [3,5], c11: [10]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [3,5], absorbed = none, pruned = none,
  merged = none
step 2: eliminate 7 (degree 2, weight 1 -> 1)  storage 31
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2], 4: [],
             5: [2,6,9], 6: [1,5,10], 7: [], 8: [0,1,9], 9: [5,8,10],
             10: [6,9], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [c4], 4: [], 5: [c4], 6: [c7], 7: [],
             8: [c7], 9: [], 10: [c11], 11: []}
  cliques = {c4: [3,5], c7: [6,8], c11: [10]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [6,8], absorbed = none, pruned = none,
  merged = none
step 3: eliminate 10 (degree 2, weight 1 -> 1) storage 28
  A       = {0: [1,3,8], 1: [0,2,6,8], 2: [1,3,5], 3: [0,2], 4: [],
             5: [2,6,9], 6: [1,5], 7: [], 8: [0,1,9], 9: [5,8], 10: [],
             11: []}
  C       = {0: [], 1: [], 2: [], 3: [c4], 4: [], 5: [c4], 6: [c7,c10],
             7: [], 8: [c7], 9: [c10], 10: [], 11: []}
  cliques = {c4: [3,5], c7: [6,8], c10: [6,9]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [6,9], absorbed = c11, pruned = none,
  merged = none
step 4: eliminate 0 (degree 3, weight 1 -> 1)  storage 23
  A       = {0: [], 1: [2,6], 2: [1,3,5], 3: [2], 4: [], 5: [2,6,9],
             6: [1,5], 7: [], 8: [9], 9: [5,8], 10: [], 11: []}
  C       = {0: [], 1: [c0], 2: [], 3: [c0,c4], 4: [], 5: [c4],
             6: [c7,c10], 7: [], 8: [c0,c7], 9: [c10], 10: [], 11: []}
  cliques = {c0: [1,3,8], c4: [3,5], c7: [6,8], c10: [6,9]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [1,3,8], absorbed = none, pruned = 1-8,
  merged = none
step 5: eliminate 2 (degree 3, weight 1 -> 1)  storage 20
  A       = {0: [], 1: [6], 2: [], 3: [], 4: [], 5: [6,9], 6: [1,5],
             7: [], 8: [9], 9: [5,8], 10: [], 11: []}
  C       = {0: [], 1: [c0,c2], 2: [], 3: [c0,c2,c4], 4: [], 5: [c2,c4],
             6: [c7,c10], 7: [], 8: [c0,c7], 9: [c10], 10: [], 11: []}
  cliques = {c0: [1,3,8], c2: [1,3,5], c4: [3,5], c7: [6,8], c10: [6,9]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [1,3,5], absorbed = none, pruned = none,
  merged = none
step 6: eliminate 3 (degree 3, weight 1 -> 1)  storage 15
  A       = {0: [], 1: [6], 2: [], 3: [], 4: [], 5: [6,9], 6: [1,5],
             7: [], 8: [9], 9: [5,8], 10: [], 11: []}
  C       = {0: [], 1: [c3], 2: [], 3: [], 4: [], 5: [c3], 6: [c7,c10],
             7: [], 8: [c3,c7], 9: [c10], 10: [], 11: []}
  cliques = {c3: [1,5,8], c7: [6,8], c10: [6,9]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [1,5,8], absorbed = c0, c2, c4, pruned = none,
  merged = none
step 7: eliminate 1 (degree 3, weight 1 -> 1)  storage 11
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [9], 6: [], 7: [],
             8: [9], 9: [5,8], 10: [], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [c1], 6: [c1,c7,c10],
             7: [], 8: [c1,c7], 9: [c10], 10: [], 11: []}
  cliques = {c1: [5,6,8], c7: [6,8], c10: [6,9]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [5,6,8], absorbed = c3, pruned = 5-6,
  merged = none
step 8: eliminate 5 (degree 3, weight 1 -> 1)  storage 7
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [c5,c7,c10],
             7: [], 8: [c5,c7], 9: [c5,c10], 10: [], 11: []}
  cliques = {c5: [6,8,9], c7: [6,8], c10: [6,9]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1,
             10: 1, 11: 1}
  neighbors = [6,8,9], absorbed = c1, pruned = 8-9,
  merged = none
step 9: eliminate 6 (degree 2, weight 1 -> 3)  storage 0
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [],
             8: [], 9: [], 10: [], 11: []}
  cliques = {c6: []}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 3, 7: 1, 8: 0, 9: 0,
             10: 1, 11: 1}
  neighbors = [8,9], absorbed = c5, c7, c10, pruned = none,
  merged = 8, 9
supervariable pivots = [11, 4, 7, 10, 0, 2, 3, 1, 5, 6]   (10 of 12)
order = [11, 4, 7, 10, 0, 2, 3, 1, 5, 6, 8, 9]
fill = 7,  nnz(L) = 37 against nnz(tril A) = 30
```

All three produce the orderings and the fill counts of 5.2 and 5.4 exactly. What differs is the
pivot count, and on these graphs the saving is modest: 2 of 4, 4 of 6, 10 of 12. The mechanism
rewards density, so it is worth measuring where the fill actually is:

```
                     n     pivots without    pivots with     saving
2D grid 8x8         64            64              54          15%
2D grid 12x12      144           144             118          18%
3D grid 4x4x4       64            64              46          28%
3D grid 5x5x5      125           125              93          25%
```

Each avoided pivot is an avoided pass over every live vertex, so a quarter fewer pivots is close
to a quarter off the dominant cost, for a test that is one integer comparison per clique member.
This is the sense in which 5.5's opening claim is meant: not a marginal optimization.

**Where the three layers stand.** It is worth naming the progression, because it mirrors one the
document has already made elsewhere. 5.1 is minimum degree written out literally, storing the fill
it creates. 5.3 keeps the graph from growing but still eliminates one vertex at a time, which
makes it *nodal*: one pivot, one column, one pass of the picker. 5.5 eliminates a group at a time
and is *supernodal* in the same sense Section 4 uses the word, and the analogy is more than
verbal, since the group in both cases is the same object, a set of columns whose patterns nest.

The two do not cash out in the same currency, though, and the difference is worth keeping
straight. Supernodes in the numeric phase buy arithmetic: dense blocks, BLAS-3 instead of BLAS-2,
cache and vectorization. Supervariables in the ordering phase buy *pivots*, and a saved pivot is a
saved pass over every live vertex. Same idea, handling a group where we handled one, applied to
different bottlenecks.

**And why this is not the whole symbolic factorization.** An ordering code simulates the
elimination and therefore already computes what symbolic factorization computes: 5.1 established
that the reachable set of a pivot is the pattern of its column of `L`, and 5.5 has now grouped
those columns into what are, at least on the examples of 5.6, the fundamental supernodes. The
whole analysis could in principle be read off a single tangled pass.

It is not, and the reason is a dependency chain rather than a preference for tidiness. Allocating
`L` exactly needs `nnz(L)`; the column counts that give `nnz(L)` cheaply need the elimination
forest (2.5); the forest needs a permutation (2.4); and the permutation is what the ordering is
computing. Nothing in that chain can be reordered. So the ordering has the patterns but has them
in the wrong sequence and before any of the sizes are known, and keeping them would mean growable
per-column storage, which forfeits the single exact allocation the staged version exists to
provide. Order, then forest, then symbolic factorization, is the order the dependencies impose.

**What each layer actually bought.** It is worth totalling this honestly, because the natural
assumption, that each layer helps a bit with everything, is wrong. Measuring 5.4's version against
5.6's on grids:

```
                  n |  pivots       degree scans   |  peak storage  |  scans   waste
                    | 5.4    5.6    5.4      5.6   |  5.4    5.6    |  needed
2D grid 12x12   144 | 144    118    10584    10125 |  528    528    |    652   15.5x
3D grid 5x5x5   125 | 125     93     8000     7419 |  600    600    |    748    9.9x
```

Three readings, none of them the obvious one.

Supervariables cut the pivot count by 15 to 28 percent, which is what 5.6 already claimed.

They do **not** reduce peak storage at all: the peaks are identical. That follows from two facts
already established rather than being a surprise. Quotient-graph storage never rises (5.3), so the
peak is always the starting value; and merges fire only once the surviving graph has become dense,
which is late. A mechanism that acts at the end cannot lower a peak that occurred at the
beginning.

And they barely help the degree cost either: 3 to 9 percent, far short of the pivot saving. The
reason is the same lateness. The pivots supervariables remove are the last ones, when the live set
has already collapsed to a handful, so each removed pivot removes a cheap scan. Thirty-two fewer
pivots out of 125 sounds substantial until one notices which thirty-two.

So the ledger reads: 5.3 solved the growth problem properly, since `A` cannot grow and the bound
holds. 5.5 is about pivot count, and its deeper value is structural rather than arithmetic, since
the supervariables it finds are the supernodes the later phases would have had to find anyway.
Neither touches the degree recomputation. The last column above is what remains: the picker
recomputes a reachable set for every live vertex at every step, while only the vertices the pivot
actually reached can have changed, so it does roughly ten to fifteen times the necessary work and
the ratio grows with `n`. That is the whole of the second problem of 5.1, still untouched after
four layers, and it is what 5.11 and 5.13 exist to solve from opposite directions.

One qualification on the 3 to 9 percent. Our version detects only mass elimination, vertices
indistinguishable from the pivot. Real AMD also hashes patterns (5.5) and so catches merges that
happen earlier in the run, not only at the collapse. The figure is therefore a lower bound on what
supervariables are worth in a complete implementation, though the shape of the conclusion, that
this is a pivot-count mechanism and not a degree-cost one, does not change.

**Nothing here improves the ordering.** It is worth saying outright, because the sequence of
increasingly sophisticated versions invites the opposite assumption. Every layer from 5.1 to 5.6
computes the *same heuristic*, and the measurements bear that out: 5.3 reproduces 5.1's
permutation exactly, and 5.5 produces a different permutation with identical fill on every graph
tried. What the layers buy is time and space, never quality.

The two vendored codes then go a step further and spend a little quality to buy a lot of speed.
Measuring 5.11's batching against one pivot at a time on the same base, the fill moves by about
half a percent on average and by a few percent on individual graphs, in *either* direction: it
came out better on two grids and worse on a third. The full table and the mechanism are in 5.11.
AMD's approximate degree is the same kind of trade; its authors report fill within a few percent
of exact (5.13).

The point to carry forward is the sign of the thing rather than its size. Every layer up to 5.10
left the ordering untouched and could be verified by demanding an identical permutation. From 5.11
onward that test no longer applies, and the only way to judge a change is to measure the fill and
accept a distribution rather than an answer.

So the whole of Section 5 is one heuristic implemented five ways. The differences are in cost, and
the fastest versions give back a fraction of a percent of quality for what the previous table
measures in factors of ten. All of the figures in this section come from the prototypes in
`experiments/ordering`, on small graphs, and are meant to show the shape of each effect rather
than to stand as benchmarks; the production implementations are where quality and speed get
measured properly.



**graph4.** Five pivots instead of eight, the largest saving of the four graphs, which fits:
it is the densest, so its cliques swallow the most vertices.

```
start                                          storage 28
  A       = {0: [2,3,4,7], 1: [3,4,6,7], 2: [0,3,5], 3: [0,1,2,6,7],
             4: [0,1,5], 5: [2,4,6], 6: [1,3,5], 7: [0,1,3]}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  cliques = {}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1}
  neighbors = none, absorbed = none, pruned = none, merged = none
step 0: eliminate 2 (degree 3, weight 1 -> 1)  storage 23
  A       = {0: [4,7], 1: [3,4,6,7], 2: [], 3: [1,6,7], 4: [0,1,5],
             5: [4,6], 6: [1,3,5], 7: [0,1,3]}
  C       = {0: [c2], 1: [], 2: [], 3: [c2], 4: [], 5: [c2], 6: [], 7: []}
  cliques = {c2: [0,3,5]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1}
  neighbors = {0,3,5}, absorbed = none, pruned = 0-3, merged = none
step 1: eliminate 4 (degree 3, weight 1 -> 1)  storage 20
  A       = {0: [7], 1: [3,6,7], 2: [], 3: [1,6,7], 4: [], 5: [6],
             6: [1,3,5], 7: [0,1,3]}
  C       = {0: [c2,c4], 1: [c4], 2: [], 3: [c2], 4: [], 5: [c2,c4],
             6: [], 7: []}
  cliques = {c2: [0,3,5], c4: [0,1,5]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1}
  neighbors = {0,1,5}, absorbed = none, pruned = none, merged = none
step 2: eliminate 6 (degree 3, weight 1 -> 1)  storage 15
  A       = {0: [7], 1: [7], 2: [], 3: [7], 4: [], 5: [], 6: [],
             7: [0,1,3]}
  C       = {0: [c2,c4], 1: [c4,c6], 2: [], 3: [c2,c6], 4: [],
             5: [c2,c4,c6], 6: [], 7: []}
  cliques = {c2: [0,3,5], c4: [0,1,5], c6: [1,3,5]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1}
  neighbors = {1,3,5}, absorbed = none, pruned = 1-3, merged = none
step 3: eliminate 5 (degree 3, weight 1 -> 1)  storage 9
  A       = {0: [7], 1: [7], 2: [], 3: [7], 4: [], 5: [], 6: [],
             7: [0,1,3]}
  C       = {0: [c5], 1: [c5], 2: [], 3: [c5], 4: [], 5: [], 6: [], 7: []}
  cliques = {c5: [0,1,3]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1}
  neighbors = {0,1,3}, absorbed = c2, c4, c6, pruned = none, merged = none
step 4: eliminate 0 (degree 3, weight 1 -> 4)  storage 0
  A       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  C       = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
  cliques = {c0: []}
  weights = {0: 4, 1: 0, 2: 1, 3: 0, 4: 1, 5: 1, 6: 1, 7: 0}
  neighbors = {1,3,7}, absorbed = c5, pruned = 1-7, 3-7, merged = 1, 3, 7
order = [2, 4, 6, 5, 0, 1, 3, 7]
fill = 4,  5 pivots instead of 8
```

### 5.7 Maintained degrees

Everything from 5.3 to 5.6 attacked the first problem of 5.1, the graph that grows. This section
begins on the second, and it begins with the observation that costs nothing to make.

The measurement closing 5.6 is that the picker recomputes a reachable set for every live vertex at
every step, and keeps one. The waste is not subtle once stated: **eliminating a pivot can only
change the degrees of the vertices it reached.** Everything else has the same `A`, the same
cliques and the same weights it had a moment ago, so its degree is still whatever it was. There is
nothing to recompute.

So keep a `degree` array and refresh only the reached set. The picker then scans cached integers
instead of building set unions:

```
minimum_degree_maintained(G) -> elimination order:

    initialize A, C, w as in 5.5
    for each vertex i:
        degree[i] = |A_i|                       # once, at the start

    while a live supervariable remains:

        p = a live supervariable minimizing degree[p]      # reads the cache
        L = the reachable set of p

        eliminate p exactly as in 5.5: absorb, prune, merge

        # REFRESH, and only here.
        for each i in L that survived the merges:
            degree[i] = sum of w_j over the reachable set of i

    expand supervariables as in 5.5
```

**Why refreshing `L` is enough**, which is the only thing in this layer that needs an argument.
There are three ways an elimination can alter what a vertex sees, and all three stay inside `L`:

- **Pruning and clique membership.** `A_i` and `C_i` are modified only for `i` in `L`. No other
  vertex's lists are touched.
- **Absorption.** A clique `c` is deleted when the pivot belonged to it, which changes `C_i` for
  every member of `c`. But if the pivot was a member of `c`, every other member of `c` is
  reachable from the pivot, so all of them are in `L` already.
- **Merging.** A merged vertex disappears entirely, which would corrupt the cached degree of any
  neighbor still expecting to see it. But a vertex merges precisely when everything it can see
  lies inside `L`, so it has no neighbors outside `L` to corrupt.

Nothing outside `L` can tell that an elimination happened. That is a stronger statement than it
looks, and it is worth checking rather than assuming when implementing, because the failure mode
is quiet: a stale degree does not crash and does not produce an invalid permutation, it produces a
valid one that is slightly worse, and nothing announces it.

This is the first half of what both vendored codes do before their ideas diverge. The second half
is degree buckets, which remove the remaining scan; that is 5.9.

### 5.8 A worked example: maintained degrees

The state is the same as 5.6's with one array added, so `graph1` is shown in full and the other
two are reduced to the two lines that are new. Nothing else about them has changed.

**graph1.**

```
start                                          degree computations 4
  A       = {0: [1,3], 1: [0,2], 2: [1,3], 3: [0,2]}
  C       = {0: [], 1: [], 2: [], 3: []}
  cliques = {}
  weights = {0: 1, 1: 1, 2: 1, 3: 1}
  degrees = {0: 2, 1: 2, 2: 2, 3: 2}
  refreshed = all (initial)
step 0: eliminate 0 (degree 2, weight 1 -> 1)  total 6
  A       = {0: [], 1: [2], 2: [1,3], 3: [2]}
  C       = {0: [], 1: [c0], 2: [], 3: [c0]}
  cliques = {c0: [1,3]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1}
  degrees = {0: -, 1: 2, 2: 2, 3: 2}
  refreshed = 1, 3
step 1: eliminate 1 (degree 2, weight 1 -> 3)  total 6
  A       = {0: [], 1: [], 2: [], 3: []}
  C       = {0: [], 1: [], 2: [], 3: []}
  cliques = {c1: []}
  weights = {0: 1, 1: 3, 2: 0, 3: 0}
  degrees = {0: -, 1: -, 2: -, 3: -}
  refreshed = none
order = [0, 1, 2, 3]
fill = 1,  nnz(L) = 9 against nnz(tril A) = 8
degree computations: 6
```

Read the `degrees` line across step 0. Vertices 1 and 3 are refreshed because the pivot reached
them; vertex 2 is not, and keeps the degree it was assigned before any elimination happened. Its
cached 2 is still correct, which is the claim of 5.7 in its smallest instance.

**graph2.**

```
start                                          degree computations 6
  degrees = {0: 2, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}
  refreshed = all (initial)
step 0: eliminate 0 (degree 2, weight 1 -> 1)  total 8
  degrees = {0: -, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}
  refreshed = 1, 2
step 1: eliminate 1 (degree 2, weight 1 -> 1)  total 10
  degrees = {0: -, 1: -, 2: 2, 3: 3, 4: 3, 5: 2}
  refreshed = 2, 3
step 2: eliminate 2 (degree 2, weight 1 -> 1)  total 12
  degrees = {0: -, 1: -, 2: -, 3: 2, 4: 2, 5: 2}
  refreshed = 3, 4
step 3: eliminate 3 (degree 2, weight 1 -> 3)  total 12
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -}
  refreshed = none
order = [0, 1, 2, 3, 4, 5]
fill = 2,  nnz(L) = 15 against nnz(tril A) = 13
degree computations: 12
```

**graph3.**

```
start                                          degree computations 12
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 3, 11: 1}
  refreshed = all (initial)
step 0: eliminate 11 (degree 1, weight 1 -> 1) total 13
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  refreshed = 10
step 1: eliminate 4 (degree 2, weight 1 -> 1)  total 15
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  refreshed = 3, 5
step 2: eliminate 7 (degree 2, weight 1 -> 1)  total 17
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: 2, 11: -}
  refreshed = 6, 8
step 3: eliminate 10 (degree 2, weight 1 -> 1) total 19
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  refreshed = 6, 9
step 4: eliminate 0 (degree 3, weight 1 -> 1)  total 22
  degrees = {0: -, 1: 4, 2: 3, 3: 4, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  refreshed = 1, 3, 8
step 5: eliminate 2 (degree 3, weight 1 -> 1)  total 25
  degrees = {0: -, 1: 4, 2: -, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  refreshed = 1, 3, 5
step 6: eliminate 3 (degree 3, weight 1 -> 1)  total 28
  degrees = {0: -, 1: 3, 2: -, 3: -, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  refreshed = 1, 5, 8
step 7: eliminate 1 (degree 3, weight 1 -> 1)  total 31
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: 3, 6: 3, 7: -, 8: 3, 9: 3,
             10: -, 11: -}
  refreshed = 5, 6, 8
step 8: eliminate 5 (degree 3, weight 1 -> 1)  total 34
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: 2, 7: -, 8: 2, 9: 2,
             10: -, 11: -}
  refreshed = 6, 8, 9
step 9: eliminate 6 (degree 2, weight 1 -> 3)  total 34
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -, 8: -, 9: -,
             10: -, 11: -}
  refreshed = none
order = [11, 4, 7, 10, 0, 2, 3, 1, 5, 6, 8, 9]
fill = 7,  nnz(L) = 37 against nnz(tril A) = 30
degree computations: 34
```

Ten steps, thirty-four degree computations, where 5.6's version would perform one per live vertex
per step. The gap widens quickly with size:

```
                     n   degree computations       ratio
                          5.6          5.7
2D grid 8x8         64      2019         302        6.7x
2D grid 12x12      144     10007         770       13.0x
2D grid 16x16      256     31095        1418       21.9x
3D grid 4x4x4       64      1901         362        5.3x
3D grid 5x5x5      125      7326         841        8.7x
3D grid 6x6x6      216     21481        1721       12.5x
```

The orderings are identical on every one, as they must be: this layer changes when a degree is
computed, never what it is. The ratio grows with `n` because the two are in different complexity
classes, `O(n^2)` set unions against `O(nnz(L))`, so the factor of twenty on a 256-vertex grid is
the beginning of the trend rather than the end of it.

What remains is the scan itself. The picker still walks every live vertex to find the smallest
cached degree, which is `O(n)` per step, though now over integers rather than set unions. That is
cheap enough to be invisible on these examples and expensive enough to matter at scale. Filing
vertices in buckets indexed by degree removes it, and that is 5.9.


**graph4.**

```
start                                          degree computations 8
  degrees = {0: 4, 1: 4, 2: 3, 3: 5, 4: 3, 5: 3, 6: 3, 7: 3}
  refreshed = all (initial)
step 0: eliminate 2 (degree 3, weight 1 -> 1)  total 11
  degrees = {0: 4, 1: 4, 2: -, 3: 5, 4: 3, 5: 4, 6: 3, 7: 3}
  refreshed = 0, 3, 5
step 1: eliminate 4 (degree 3, weight 1 -> 1)  total 14
  degrees = {0: 4, 1: 5, 2: -, 3: 5, 4: -, 5: 4, 6: 3, 7: 3}
  refreshed = 0, 1, 5
step 2: eliminate 6 (degree 3, weight 1 -> 1)  total 17
  degrees = {0: 4, 1: 4, 2: -, 3: 4, 4: -, 5: 3, 6: -, 7: 3}
  refreshed = 1, 3, 5
step 3: eliminate 5 (degree 3, weight 1 -> 1)  total 20
  degrees = {0: 3, 1: 3, 2: -, 3: 3, 4: -, 5: -, 6: -, 7: 3}
  refreshed = 0, 1, 3
step 4: eliminate 0 (degree 3, weight 1 -> 4)  total 20
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -}
  refreshed = none
order = [2, 4, 6, 5, 0, 1, 3, 7]
degree computations: 20
```

### 5.9 Degree buckets

5.7 removed the recomputation. What it left is the scan: the picker still walks every live vertex
to read its cached degree and keep the smallest, `O(n)` per pivot over integers rather than set
unions. Cheap per step and the only remaining `O(n)` per pivot.

File each supervariable in a bucket indexed by its degree, and the minimum can be found by walking
**up** from the last known minimum instead of looking at everything:

```
minimum_degree_bucketed(G) -> elimination order:

    initialize A, C, w, degree as in 5.7
    for each vertex i:
        bucket[degree[i]] = bucket[degree[i]] union {i}
    mdeg = min over i of degree[i]              # a LOWER BOUND, not the true minimum

    while a live supervariable remains:

        while bucket[mdeg] is empty:  mdeg = mdeg + 1     # the walk
        p = a member of bucket[mdeg]
        remove p from bucket[mdeg]

        eliminate p as in 5.5; merged vertices leave their buckets too

        for each i in L that survived:
            d = the refreshed degree of i
            move i from bucket[degree[i]] to bucket[d]     # REFILE
            degree[i] = d
            if d < mdeg:  mdeg = d                         # the bound may only fall
```

Two things carry the weight. `mdeg` is a **lower bound** on the current minimum, not the minimum
itself: it is allowed to lag, and the walk corrects it. What it must never do is overshoot, since
buckets below `mdeg` are never examined and a vertex sitting in one would never be chosen. That is
why the refile lowers `mdeg` and nothing ever raises it except the walk.

And a vertex whose degree changes must be pulled out of the **middle** of its old bucket, which is
why MMD's degree lists carry both `fwd` and `bwd` links and AMD's carry `Next` and `Last`. The
doubly linked list is not decoration; singly linked would make the refile `O(bucket size)` and
give the whole structure back.

**This is a priority queue question.** Naming it that way makes the sequence legible. The picker
needs repeated find-min over keys that change constantly, which is exactly what a priority queue
is for, and the last three sections have been walking up the standard implementations:

```
5.5   recompute every key on every query          no queue at all
5.7   maintained keys, linear search              unsorted array:  O(1) update, O(n) find-min
5.9   keys bucketed by value                      bucket queue:    O(1) update, O(1) amortized
```

The conspicuous absence is the **binary heap**, `O(log n)` for both operations, which would be the
obvious middle step and is deliberately skipped. Degrees move in *both* directions here: pruning
and merging lower them, fill raises them, and a single refresh can do either. Heaps handle
decrease-key tolerably and increase-key badly. A bucket queue does not care about direction, since
a refile is a removal and an insertion whichever way the key moved. The heap would be a worse
middle, not a smaller one.

**When bucket queues are actually fast**, which is the part worth being careful about. The `O(1)`
amortized claim needs the keys to be small bounded integers, which degrees are, and it needs the
minimum to be roughly monotone, so that the total distance walked over a whole run is bounded
rather than paid afresh at every step. That second condition is the fragile one, and this
algorithm does not strictly satisfy it: a vertex's degree genuinely falls when its neighbors are
eliminated, so `mdeg` can drop below where it had reached. Measured on our prototype it drops only
three to seven times in an entire run, and the total walk comes to roughly 0.1 to 0.3 buckets per
pivot, so the structure behaves as though it were monotone. It is not, and no clean bound follows;
the evidence is empirical.

**An aside, on three ways to bucket a changing key.** The same structure appears in matching
algorithms under the same pressure, and comparing them isolates the design choice. Gabow's
cardinality matching buckets edges by the dual value at which they are *projected* to become
tight, `L[Delta + p]`, computed from duals that then change; the entries can go stale, so the
extraction loop carries a guard that discards any edge whose duals no longer make it tight. That
is **lazy invalidation**: insert eagerly, verify on extraction, never remove. Micali and
Vazirani's algorithm buckets bridges by tenacity, and when a bridge's tenacity is not yet
computable it refuses to guess, parking the edge in a per-vertex `hanging` list until its level is
known and only then filing it at the correct key. That is **deferred insertion**. Our version does
neither: it moves the entry the moment the key changes, so a stale entry never exists. That is
**eager refiling**.

The choice is forced by monotonicity rather than taste. Both matching codes advance their key in
one direction only, a `Delta` that is only ever incremented and a phase index that only counts up,
so a stale entry is encountered at most once, at a level that was going to be visited anyway.
Minimum degree has a falling key, and lazy entries would accumulate in buckets already passed,
below `mdeg`, where nothing will ever look at them again. Eager refiling is what a non-monotone
key requires, and the doubly linked lists are what make it affordable.

### 5.10 A worked example: degree buckets

The state is 5.8's with the buckets added, so `graph1` is shown in full and the other two are
reduced to the lines that are new.

**graph1.**

```
start
  A       = {0: [1,3], 1: [0,2], 2: [1,3], 3: [0,2]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1}
  degrees = {0: 2, 1: 2, 2: 2, 3: 2}
  buckets = {2: [0,1,2,3]}   mdeg = 2
step 0: eliminate 0 (degree 2, weight 1 -> 1)  walked 0
  A       = {0: [], 1: [2], 2: [1,3], 3: [2]}
  weights = {0: 1, 1: 1, 2: 1, 3: 1}
  degrees = {0: -, 1: 2, 2: 2, 3: 2}
  buckets = {2: [1,2,3]}   mdeg = 2
  refreshed = 1, 3
step 1: eliminate 1 (degree 2, weight 1 -> 3)  walked 0
  A       = {0: [], 1: [], 2: [], 3: []}
  weights = {0: 1, 1: 3, 2: 0, 3: 0}
  degrees = {0: -, 1: -, 2: -, 3: -}
  buckets = {}   mdeg = 2
  refreshed = none
order = [0, 1, 2, 3]
fill = 1,  nnz(L) = 9 against nnz(tril A) = 8
bucket probes: 2 for 2 pivots
```

Every step walks zero buckets: the pivot is found in the first slot examined. That is the common
case and the reason the structure works, since eliminating a vertex tends to lower its neighbors'
degrees and pull `mdeg` straight back down to where the search will start next time.

**graph2.**

```
start
  degrees = {0: 2, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}
  buckets = {2: [0,1,2,5], 3: [3,4]}   mdeg = 2
step 0: eliminate 0 (degree 2, weight 1 -> 1)  walked 0
  degrees = {0: -, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}
  buckets = {2: [1,2,5], 3: [3,4]}   mdeg = 2
  refreshed = 1, 2
step 1: eliminate 1 (degree 2, weight 1 -> 1)  walked 0
  degrees = {0: -, 1: -, 2: 2, 3: 3, 4: 3, 5: 2}
  buckets = {2: [2,5], 3: [3,4]}   mdeg = 2
  refreshed = 2, 3
step 2: eliminate 2 (degree 2, weight 1 -> 1)  walked 0
  degrees = {0: -, 1: -, 2: -, 3: 2, 4: 2, 5: 2}
  buckets = {2: [3,4,5]}   mdeg = 2
  refreshed = 3, 4
step 3: eliminate 3 (degree 2, weight 1 -> 3)  walked 0
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -}
  buckets = {}   mdeg = 2
  refreshed = none
order = [0, 1, 2, 3, 4, 5]
fill = 2,  nnz(L) = 15 against nnz(tril A) = 13
bucket probes: 4 for 4 pivots
```

**graph3.**

```
start
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 3, 11: 1}
  buckets = {1: [11], 2: [4,7], 3: [0,2,3,9,10], 4: [1,5,6,8]}   mdeg = 1
step 0: eliminate 11 (degree 1, weight 1 -> 1) walked 0
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  buckets = {2: [4,7,10], 3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 1
  refreshed = 10
step 1: eliminate 4 (degree 2, weight 1 -> 1)  walked 1
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  buckets = {2: [7,10], 3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 2
  refreshed = 3, 5
step 2: eliminate 7 (degree 2, weight 1 -> 1)  walked 0
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: 2, 11: -}
  buckets = {2: [10], 3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 2
  refreshed = 6, 8
step 3: eliminate 10 (degree 2, weight 1 -> 1) walked 0
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 2
  refreshed = 6, 9
step 4: eliminate 0 (degree 3, weight 1 -> 1)  walked 1
  degrees = {0: -, 1: 4, 2: 3, 3: 4, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [2,9], 4: [1,3,5,6,8]}   mdeg = 3
  refreshed = 1, 3, 8
step 5: eliminate 2 (degree 3, weight 1 -> 1)  walked 0
  degrees = {0: -, 1: 4, 2: -, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [3,9], 4: [1,5,6,8]}   mdeg = 3
  refreshed = 1, 3, 5
step 6: eliminate 3 (degree 3, weight 1 -> 1)  walked 0
  degrees = {0: -, 1: 3, 2: -, 3: -, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [1,9], 4: [5,6,8]}   mdeg = 3
  refreshed = 1, 5, 8
step 7: eliminate 1 (degree 3, weight 1 -> 1)  walked 0
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: 3, 6: 3, 7: -, 8: 3, 9: 3,
             10: -, 11: -}
  buckets = {3: [5,6,8,9]}   mdeg = 3
  refreshed = 5, 6, 8
step 8: eliminate 5 (degree 3, weight 1 -> 1)  walked 0
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: 2, 7: -, 8: 2, 9: 2,
             10: -, 11: -}
  buckets = {2: [6,8,9]}   mdeg = 2
  refreshed = 6, 8, 9
step 9: eliminate 6 (degree 2, weight 1 -> 3)  walked 0
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -, 8: -, 9: -,
             10: -, 11: -}
  buckets = {}   mdeg = 2
  refreshed = none
order = [11, 4, 7, 10, 0, 2, 3, 1, 5, 6, 8, 9]
fill = 7,  nnz(L) = 37 against nnz(tril A) = 30
bucket probes: 12 for 10 pivots
```

Twelve probes for ten pivots, so the walk contributed two steps in the whole run. At scale:

```
                     n   scanned (5.7)   probes (5.9)     ratio
2D grid 8x8         64          2019             63       32x
2D grid 12x12      144         10007            136       74x
2D grid 16x16      256         31095            227      137x
3D grid 4x4x4       64          1901             61       31x
3D grid 5x5x5      125          7326            123       60x
3D grid 6x6x6      216         21481            198      108x
```

The probe count is essentially the pivot count, 227 against 202 pivots on the largest grid, which
says the walk is close to free in aggregate. The orderings are identical to 5.8's on every graph
tried, including random ones, which is the thing to check when replacing a linear scan with a
structure: the tie-break has to survive. Taking the lowest-numbered member of `bucket[mdeg]`
reproduces what `min` over the live set was already doing.

**And still nothing has improved the ordering.** 5.7 changed when a degree is computed, 5.9 changed
how the smallest is found, and neither changed what any degree *is*. The sequence has now spent
four sections on speed and none on quality, which is the whole point: the heuristic was fixed in
5.1 and everything since has been implementation.

That statement needs one refinement, though, because "unchanged ordering" is not quite what holds
across all of it, and the exception is instructive. Three things can happen when a layer is added,
and all three occur in this section:

```
                              order        fill        what the change is
5.1 -> 5.3                    same         same        a change of representation
5.3 -> 5.5  (mass elim.)      DIFFERENT    same        a provably free reordering
5.5 -> 5.7 -> 5.9             same         same        a change of implementation
5.9 -> 5.11 (multiple elim.)  different    DIFFERENT   a wager
```

So 5.11 is not the first layer to change the permutation. Mass elimination already does, and
substantially: on twelve test graphs it produced a different ordering from 5.3's on nine of them,
with identical `nnz(L)` on all twelve. What 5.11 is the first to change is the *fill*.

The middle row deserves its adjective. The reordering is not merely harmless in practice, it is
free by construction: a merged vertex creates no fill when eliminated next, so promoting it ahead
of whatever the degree would otherwise have chosen cannot cost anything. Nor is the promotion
forced. Measured at the moment of each merge, the merged vertex would have been the global
minimum-degree choice only about seventy percent of the time; in the other thirty the ordering is
genuinely different from what plain minimum degree would produce, and the fill is identical anyway.

It follows that making 5.5 agree with 5.3 is easy and pointless. One would recognize
indistinguishability, use it to skip the degree recomputation, and still leave the vertex in its
bucket as a separate pivot. The orderings would match, and the entire pivot saving of 5.6 would be
gone, because that saving *is* the promotion. The differing order is the design working, not a
defect in it.

**Two kinds of free reordering, and one that is not.** Collecting them, since the distinction is
what the rest of the section turns on:

- **Consecutive numbering of indistinguishable vertices** (5.5). Free: the promoted vertex adds no
  fill.
- **Interleaving between connected components** (5.4, and 5.11). Free: total fill is the sum of the
  components' fills and each depends only on the order within it.
- **Reordering ties inside one component** (5.11's batch). *Not* free, because a tie-break is
  precisely where this heuristic's quality lives.

Which sets up what remains. 5.11 and 5.13 both start from exactly this base, quotient graph,
supervariables, maintained degrees, buckets, and both go further by giving something up. MMD gives
up its freedom to reconsider a vertex the current batch has touched; AMD gives up the exactness of
the degree. Each buys speed with a small, two-sided change in the fill, and the fact that they must
give something up is the clearest sign that the free wins have all been taken.


**graph4.**

```
start
  degrees = {0: 4, 1: 4, 2: 3, 3: 5, 4: 3, 5: 3, 6: 3, 7: 3}
  buckets = {3: [2,4,5,6,7], 4: [0,1], 5: [3]}   mdeg = 3
step 0: eliminate 2 (degree 3, weight 1 -> 1)  walked 0
  degrees = {0: 4, 1: 4, 2: -, 3: 5, 4: 3, 5: 4, 6: 3, 7: 3}
  buckets = {3: [4,6,7], 4: [0,1,5], 5: [3]}   mdeg = 3
  refreshed = 0, 3, 5
step 1: eliminate 4 (degree 3, weight 1 -> 1)  walked 0
  degrees = {0: 4, 1: 5, 2: -, 3: 5, 4: -, 5: 4, 6: 3, 7: 3}
  buckets = {3: [6,7], 4: [0,5], 5: [1,3]}   mdeg = 3
  refreshed = 0, 1, 5
step 2: eliminate 6 (degree 3, weight 1 -> 1)  walked 0
  degrees = {0: 4, 1: 4, 2: -, 3: 4, 4: -, 5: 3, 6: -, 7: 3}
  buckets = {3: [5,7], 4: [0,1,3]}   mdeg = 3
  refreshed = 1, 3, 5
step 3: eliminate 5 (degree 3, weight 1 -> 1)  walked 0
  degrees = {0: 3, 1: 3, 2: -, 3: 3, 4: -, 5: -, 6: -, 7: 3}
  buckets = {3: [0,1,3,7]}   mdeg = 3
  refreshed = 0, 1, 3
step 4: eliminate 0 (degree 3, weight 1 -> 4)  walked 0
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -}
  buckets = {}   mdeg = 3
  refreshed = none
order = [2, 4, 6, 5, 0, 1, 3, 7]
bucket probes: 5 for 5 pivots
```

### 5.11 Multiple elimination, which is the idea in MMD

Sections 5.3 to 5.6 disposed of the first of the two problems of 5.1, the graph that grows: the
quotient graph never forms a clique, and supervariables shrink what remains. The second problem,
the expense of the degree, is still open, and it is exactly here that the two codes we vendor part
company.

**What both codes do first, and neither section is about.** Before either of the two mechanisms,
AMD and MMD do the same two things, and both are already in place by the end of 5.10: the
maintained degree of 5.7, refreshing only the vertices the pivot reached, and the degree buckets
of 5.9, turning the search for the minimum into a short walk. Both are common ground rather than
either paper's contribution, which is worth saying plainly, because it means the fork below is
narrower than it first appears. What the two codes disagree about is the refreshes that remain.

There are two ways to make an expensive step affordable: do it less often, or make each one
cheaper. MMD does it less often. It keeps the degree *exact*, the honest set union of 5.1, and pays
for that by recomputing it rarely. AMD makes each one cheaper (5.13). Multiple elimination is the
mechanism behind "less often", and it is the idea the M in MMD names.

Ordinary minimum degree eliminates one vertex, refreshes every degree that elimination disturbed,
then picks again. The refresh is the costly step, and running it after every single pivot is the
waste. Liu's observation (1985) is that several least-degree vertices can usually be eliminated
*before any refresh*, so long as they do not interfere: if two minimum-degree vertices are
non-adjacent and share no vertex about to be eliminated with them, then eliminating one leaves the
other's degree untouched, so both are still current and both can go. A maximal independent set of
minimum-degree vertices is eliminated as a batch, and only then are degrees recomputed, once, for
everything the whole batch reached.

We never search for that independent set. It falls out of the bookkeeping. Eliminating a vertex
evicts every vertex it reaches from the degree buckets, to be re-filed only after the batch, so
whatever is still sitting in the minimum-degree bucket is by construction non-adjacent to everyone
already eliminated in this batch. Draining the bucket drains an independent set. A tolerance
`delta` may widen the batch to vertices within `delta` of the minimum, trading a little ordering
quality for still fewer refreshes; we pass `delta = 0`, which keeps the batch to true minima.

**Multiple elimination is not mass elimination**, and the names invite the confusion, so it is
worth separating them before going further. Both group several eliminations together; they group
opposite things for opposite reasons.

```
                    what is grouped        relation to the pivot   why it is safe
mass (5.5)          vertices of L_p        adjacent, inside it     they are one object
multiple (here)     other pivots           non-adjacent, disjoint  they cannot interact
```

Mass elimination merges vertices indistinguishable *from the pivot*. They lie inside its clique,
so they are as adjacent to it as vertices can be, and merging them is not a wager but the
recognition that there was only ever one object there. It happens *within* one elimination step
and yields one supervariable, one column block of the factor.

Multiple elimination eliminates several *pivots* that cannot see each other. They are non-adjacent
by construction, so eliminating one leaves the others' degrees untouched and all of them may go
before anything is recomputed. It happens *across* several elimination steps and yields several
separate supervariables.

Two consequences follow. Mass elimination changes *what* is eliminated: the graph shrinks, one
supervariable stands in for many vertices. Multiple elimination changes nothing about what is
eliminated, only *when the degrees are refreshed*; strip the batching out of MMD and it would
choose the same pivots, refreshing more often. And mass elimination is exact, which 5.6 verified
by measurement, whereas multiple elimination is a genuine wager: the second pivot in a batch is
chosen on a degree computed before the first was eliminated.

**Batching, and two ways to justify it.** The technique here is general enough to name. When a
loop alternates cheap local work with an expensive global refresh, the refresh can be amortized by
performing many local steps before paying for one refresh, provided something guarantees the
steps do not invalidate each other. Hopcroft-Karp is the familiar instance: rather than one
augmenting path per iteration it augments along a maximal set of *vertex-disjoint* shortest paths
per phase, and the disjointness is precisely what makes the batch safe to apply without
recomputing in between. Multiple elimination is the same shape, with non-adjacency in the role of
vertex-disjointness.

Section 5.5 batches too, and it is worth seeing that the justification is the opposite one. Mass
elimination merges vertices that are *indistinguishable from the pivot*: they sit inside its
clique and are as dependent on it as vertices can be. They may be eliminated together because they
are effectively one object. Multiple elimination batches vertices that cannot see each other at
all. Independence and identity are the two ways a group of eliminations can be safe, and this
section and 5.5 take one each.

Two limits on the analogy. Hopcroft-Karp's batching buys an asymptotic improvement, from `O(V E)`
to `O(sqrt(V) E)`, because the number of phases is provably bounded; multiple elimination is a
practical constant-factor win with no comparable bound. And Hopcroft-Karp must *search* for its
disjoint set, with a layering pass and then a set of disjoint augmentations, where MMD gets its
independent set for free from the bucket eviction described above. The mechanism transfers; the
guarantee does not.

```
mmd(A) -> elimination order:

    for each vertex i:                            # mmdint
        A_i = the neighbors of i in A
        E_i = {}                                  # no elements yet
        w_i = 1                                   # supervariable weight
        degree[i] = |A_i|                         # exact degree
        file i in the bucket for degree[i]

    # PREPASS: degree 0 and degree 1 create no fill, so number them now and     # genmmd
    # leave their neighbors' degrees stale. A stale degree costs ordering
    # quality, never correctness: an eliminated vertex is skipped on sight.
    for each i with degree[i] <= 1:
        number i next;  eliminate i
    if no vertex remains:  goto FINISH

    while vertices remain:                          # genmmd, outer loop

        mdeg = the least non-empty bucket

        # ELIMINATE A MAXIMAL INDEPENDENT SET at mdeg as one batch, refreshing
        # no degree until the batch is complete. This is multiple elimination.
        batch = {}
        while the mdeg bucket still holds a vertex p:
            remove p from its bucket

            # FORM p's ELEMENT: its reachable set in the quotient graph.        # mmdelm
            L = A_p
            for e in E_p:
                L = L union L_e                     # ELEMENT ABSORPTION: e dies into p
            L = L \ {p}
            L_p = L;  p is now an element

            for i in L:
                A_i = A_i \ (eliminated and absorbed vertices)
                E_i = (E_i with dead elements dropped) union {p}
                evict i from its degree bucket      # <- this is what keeps the batch independent
                if A_i \ L is empty:                # i sees nothing outside the new clique
                    w_p = w_p + w_i                 # MASS ELIMINATION: merge i into p
                    absorb i

            reserve the next w_p numbers for p;  add p to batch

        # REFRESH once, for the whole batch: the exact external degree of        # mmdupd
        # every vertex the batch reached, then re-file it.
        for i reached by batch:
            degree[i] = | A_i union ( union of L_e over e in E_i ) | \ {i}
            re-file i in the bucket for degree[i]
            # the code splits this by whether i sees exactly two elements or many
            # (q2h vs qxh); the two-element pass also catches indistinguishable
            # vertices the A_i \ L test above missed. Collapsed here.

FINISH:                                              # mmdnum
    resolve each supervariable: every merged vertex takes one number inside its
    representative's reserved block, giving a consecutive ordering
    return the numbering
```

**Where the independence comes from, in the code.** The pseudocode says "evict `i` from its degree
bucket", and it is worth seeing that this single line is the whole of the batching discipline.
In `mmdelm`, every reached vertex is spliced out of its degree list:

```c
int pv = bwd[rn];
if (pv != 0 && pv != (-maxint)) {
    int nx = fwd[rn];
    if (nx > 0) bwd[nx] = pv;
    if (pv > 0) fwd[pv] = nx;
    if (pv < 0) head[-pv] = nx;      /* rn leaves its bucket */ }
```

That is the doubly linked removal of 5.9, and here is what it is for. Because every vertex the
pivot reached has just left the buckets, anything *still* sitting in the `mdeg` bucket when the
loop comes back for another pivot was not reached, hence is non-adjacent to every pivot already
taken this round. Draining the bucket drains an independent set, and no code anywhere computes
one. The batch is a side effect of the eviction.

The driver makes the two-level structure plain, with the inner loop taking pivots and the outer
one refreshing:

```c
while (1) {
    while (head[mdeg] <= 0) mdeg++;
    int mdlmt = mdeg + delta, ehead = 0; n500: ...
    mmdelm(mn, ...);                 /* eliminate; refresh nothing */
    list[mn] = ehead; ehead = mn;    /* remember it for the batch refresh */
    if (delta >= 0) goto n500;       /* another pivot, degrees still stale */ n900:
    mmdupd(ehead, ...);              /* ONE refresh for the whole batch */ }
```

**What the batch buys and what it costs.** Measured against the same code with the batch limited
to a single pivot, so that only the batching differs:

```
                 pivots | refreshes  batched  saved | batch size |  fill change
2D grid 24x24       447 |      3417     2427    28% |       10.0 |    -0.70%
2D grid 16x16       202 |      1418     1094    22% |        5.8 |    -0.96%
2D grid 12x12       118 |       770      597    22% |        4.5 |    +0.68%
3D grid 5x5x5        93 |       841      558    33% |        6.6 |     0.00%
3D grid 6x6x6       156 |      1721     1458    15% |        4.4 |    +3.70%
random graphs                              15 to 23% |   3.1-4.3 |  0.00 to +0.37%
overall                                                          |    +0.40%
```

Batches run three to ten pivots and remove 15 to 33 percent of the refreshes. The saving is
precisely the overlap between reached sets: a vertex reached by two pivots of one batch is
refreshed once rather than twice, and non-adjacent pivots may still share neighbors.

**That column counts refreshes, not work, and the two are not the same claim.** A refresh is one
vertex having its degree recomputed; the work inside it is the walk over the members of that
vertex's elements. Batching reduces how *often* a vertex is refreshed while increasing what each
refresh *costs*, because more new elements exist by the time the batch ends. Counted in element
members walked rather than in refreshes, the same runs give a saving of roughly six percent, not
fifteen to thirty-three. Both numbers are honest and they answer different questions; the second
is closer to time on a machine. We report the first because it isolates the mechanism, and record
the second here so the first is not read as more than it is.

**What the batch actually gives up**, which is worth pinning down because the obvious answer is
wrong. It is tempting to say that the later pivots of a batch are chosen on stale degrees, and
that is precisely what does *not* happen: the pivots are non-adjacent, so none of them disturbs
another's degree, and every one of them is a genuine minimum-degree vertex when it is taken. The
heuristic is obeyed throughout.

What changes is *which* minimum-degree vertex. Every vertex the batch has already reached was
evicted from the buckets and will not be re-filed until the round ends, so the choice is made
among the untouched remainder rather than among all the candidates. Tracing the two versions
against each other makes this visible: on a 16x16 grid the pivot sequences diverge at step 11
while the sequence of pivot *degrees* stays identical until step 114. For a hundred steps the two
runs take different vertices of the same degree.

So the batch does not make a worse choice, it makes a different one, and minimum degree is
famously sensitive to exactly that. A tie-break is arbitrary by construction, so perturbing it
moves the fill in either direction: two of the grids above came out better than one-at-a-time and
the 6x6x6 came out worse. The mean is slightly bad and the spread exceeds the mean. The honest
description is not "a small penalty" but "the ordering becomes arbitrary in a different way", and
that framing predicts the sign pattern where "penalty" does not.

**Where batching is genuinely free.** The distinction is not between adjacent and non-adjacent
pivots but between pivots in the same connected component and pivots in different ones. The fill
of a reducible matrix is the sum of its components' fills (5.4), and each component's fill depends
only on the relative order *within* that component. Interleaving between components therefore
reorders nothing that matters. Restricting a batch to at most one pivot per component gives an
ordering identical to one-at-a-time, which we checked on six disconnected graphs; draining the
whole bucket, which can take several pivots from one component, differed on four of the six.

That is the clean boundary between the two kinds of batching. Across components it is free, in the
same way mass elimination is free: no wager, just a reordering that provably cannot matter. Within
a component it is a wager, small and two-sided. MMD does not distinguish the cases, because
tracking components would cost more than it saves on a matrix that is usually irreducible, but
knowing which side of the line a given batch sits on explains the whole of its behavior.

It also separates two questions that sound like one. Identical *fill* is reachable: restrict the
batch to one pivot per component and 5.10's guarantee is restored, at the price of maintaining
component labels. Identical *ordering* is not, because within a component the batch would have to
be a single pivot, which is multiple elimination deleted rather than constrained. So the mechanism
can be made lossless but not order-preserving, and since the fill is what matters, that is the
useful half to be able to recover.

These are prototype measurements on small graphs; the production numbers are the ones to trust.

The same four mechanisms as AMD, minus one and plus one. MMD keeps the quotient graph with element
absorption (form `L`, absorb `E_p`), keeps supervariables (the `w` merges), and keeps mass
elimination (the empty `A_i \ L` test). What it drops is AMD's approximate degree: the refresh
computes the exact set union, the very quantity 5.13 will bound rather than compute. What it adds is
multiple elimination, the batching that makes an exact refresh affordable by making it rare. Its
supervariable detection is structural, done during the refresh, where AMD's is a hash (5.13).

So the two codes are one algorithm forking on the second problem of 5.1. Both refuse to form a
clique and both shrink the graph; they differ only on the cost of the degree. AMD makes each update
cheap and runs one per pivot; MMD keeps each update exact and runs one per batch. Approximate
degree and multiple elimination are the same answer, do the expensive thing less, reached from
opposite ends.

The vendored code is five routines, and this table doubles as a reading key, since the pseudocode
above is annotated with the same names:

| routine | role |
|---|---|
| `mmdint` | bucket every vertex by degree; the buckets are doubly linked so a vertex can be pulled from the middle |
| `genmmd` | the driver: the prepass, then the outer loop that finds `mdeg`, drains its bucket as a batch, and calls the refresh |
| `mmdelm` | eliminate one vertex: form its reachable set `L`, splice the quotient graph, absorb dead elements, merge indistinguishable vertices |
| `mmdupd` | the batch refresh: exact external degree of every reached vertex, then re-file (the `q2h`/`qxh` split is the two-element versus many-element case) |
| `mmdnum` | recover the permutation, expanding each `w`-merged group into consecutive numbers |

Two coarse-grainings in the pseudocode above are deliberate: the `q2h`/`qxh` split inside `mmdupd`
(the two-element versus many-element reachable cases, collapsed into one refresh loop) and the
exact number-reservation bookkeeping between the batch and `mmdnum`. Neither changes the shape.

**What the prototype implements, and what it does not.** The pseudocode above describes the
vendored routine. Our `experiments/ordering/mmd` is a subset of it, and the gap is recorded here
rather than only in the file headers, so that it can be closed deliberately.

Present: the quotient graph with element absorption and pruning, supervariables by mass
elimination, maintained degrees refreshed only where they can have changed, degree buckets with a
lagging `mdeg`, and multiple elimination with a `delta` parameter defaulting to zero. That is
enough to reproduce the behavior this section describes, which is why the measurements above come
from it.

Absent, in two categories. Two are **algorithmic**, and closing them would change the ordering:

- **The prepass.** `genmmd` numbers every degree-0 and degree-1 vertex before entering the main
  loop, leaving their neighbors' degrees stale. Such a vertex creates no fill, so eliminating it
  early is safe; the staleness it leaves behind is a small ordering concession for a real saving.
  We enter the main loop immediately.
- **Pair merging in `mmdupd`.** `mmdelm` stashes each reached vertex's pruned adjacency count as
  `fwd[rn] = nq + 1`, and the refresh routes the `nq == 1` cases into the `q2h` list, where it
  merges vertices indistinguishable *from each other*. Our merge test fires only for vertices
  indistinguishable *from the pivot*, so MMD's supervariables are at least as coarse as ours and
  sometimes strictly coarser. This is the larger of the two gaps: it costs both supervariables we
  never find and the pivots we would have saved by finding them.

Two are **conventions**, and they matter only when comparing against the vendored code, though
they are not quite cosmetic:

- **MMD never uses bucket 0.** `mmdint` maps degree 0 to 1, and `mmdupd` floors with
  `if (dg < 1) dg = 1`. Its least bucket is 1 where ours is 0.
- **MMD files at `dg - qsize[en] + 1`**, not at the plain external degree. The offset is a
  monotone shift, so it cannot change which vertex is minimal; but the flooring above collapses
  degrees 0 and 1 into one bucket, and since this section has established that MMD's whole quality
  difference is a tie-break effect, a convention that merges two buckets can break a tie
  differently. Small, and not zero.

The `tag` and `marker` machinery with its overflow reset, and the `ncsub` statistic, have no
counterpart in the prototype and need none: the first is a way to do set membership without sets,
and the second is a diagnostic.

### 5.12 A worked example: multiple elimination

The state is 5.10's, reduced to the degrees and buckets, with the round structure made visible.
Two fields are new and carry the whole idea: `evicted` lists the vertices a pivot pulled out of
the buckets, and `batch` lists what the round managed to take before the buckets at `mdeg` ran dry.

**graph1.** The 4-cycle is the smallest graph that shows batching at all, and it shows it
immediately.

```
start
  degrees = {0: 2, 1: 2, 2: 2, 3: 2}
  buckets = {2: [0,1,2,3]}   mdeg = 2
round 0: mdeg = 2
  eliminate 0 (degree 2, weight 1 -> 1)   merged = none   evicted = 1, 3
  eliminate 2 (degree 2, weight 1 -> 1)   merged = none   evicted = 1, 3
  batch of 2: [0, 2]    refreshed = 1, 3
  degrees = {0: -, 1: 1, 2: -, 3: 1}
  buckets = {1: [1,3]}   mdeg = 1
round 1: mdeg = 1
  eliminate 1 (degree 1, weight 1 -> 2)   merged = 3   evicted = none
  batch of 1: [1]    refreshed = none
  degrees = {0: -, 1: -, 2: -, 3: -}
  buckets = {}   mdeg = 1
order = [0, 2, 1, 3]
rounds = 2, pivots = 3, average batch = 1.5
fill = 1,  nnz(L) = 9 against nnz(tril A) = 8
degree computations: 6
```

Round 0 takes two pivots. Vertices 0 and 2 are the diagonal of the 4-cycle, so they are
non-adjacent, and that is exactly why both survive to be taken: eliminating 0 evicts its neighbors
1 and 3, and 2 is not among them. Nothing tested whether 0 and 2 were adjacent. The bucket simply
still held 2 after 0 was done with it.

The `evicted` lines also show where the saving comes from. Both pivots evict the same pair, 1 and
3, and the refresh runs once over that pair rather than twice. One round, two pivots, one refresh.

**graph2.**

```
start
  degrees = {0: 2, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}
  buckets = {2: [0,1,2,5], 3: [3,4]}   mdeg = 2
round 0: mdeg = 2
  eliminate 0 (degree 2, weight 1 -> 1)   merged = none   evicted = 1, 2
  eliminate 5 (degree 2, weight 1 -> 1)   merged = none   evicted = 3, 4
  batch of 2: [0, 5]    refreshed = 1, 2, 3, 4
  degrees = {0: -, 1: 2, 2: 2, 3: 2, 4: 2, 5: -}
  buckets = {2: [1,2,3,4]}   mdeg = 2
round 1: mdeg = 2
  eliminate 1 (degree 2, weight 1 -> 1)   merged = none   evicted = 2, 3
  eliminate 4 (degree 2, weight 1 -> 1)   merged = none   evicted = 2, 3
  batch of 2: [1, 4]    refreshed = 2, 3
  degrees = {0: -, 1: -, 2: 1, 3: 1, 4: -, 5: -}
  buckets = {1: [2,3]}   mdeg = 1
round 2: mdeg = 1
  eliminate 2 (degree 1, weight 1 -> 2)   merged = 3   evicted = none
  batch of 1: [2]    refreshed = none
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -}
  buckets = {}   mdeg = 1
order = [0, 5, 1, 4, 2, 3]
rounds = 3, pivots = 5, average batch = 1.7
fill = 2,  nnz(L) = 15 against nnz(tril A) = 13
degree computations: 12
```

**graph3.**

```
start
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 3, 11: 1}
  buckets = {1: [11], 2: [4,7], 3: [0,2,3,9,10], 4: [1,5,6,8]}   mdeg = 1
round 0: mdeg = 1
  eliminate 11 (degree 1, weight 1 -> 1)   merged = none   evicted = 10
  batch of 1: [11]    refreshed = 10
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  buckets = {2: [4,7,10], 3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 1
round 1: mdeg = 2
  eliminate 4 (degree 2, weight 1 -> 1)   merged = none   evicted = 3, 5
  eliminate 7 (degree 2, weight 1 -> 1)   merged = none   evicted = 6, 8
  eliminate 10 (degree 2, weight 1 -> 1)   merged = none   evicted = 6, 9
  batch of 3: [4, 7, 10]    refreshed = 3, 5, 6, 8, 9
  degrees = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 2
round 2: mdeg = 3
  eliminate 0 (degree 3, weight 1 -> 1)   merged = none   evicted = 1, 3, 8
  eliminate 2 (degree 3, weight 1 -> 1)   merged = none   evicted = 1, 3, 5
  eliminate 9 (degree 3, weight 1 -> 1)   merged = none   evicted = 5, 6, 8
  batch of 3: [0, 2, 9]    refreshed = 1, 3, 5, 6, 8
  degrees = {0: -, 1: 4, 2: -, 3: 3, 4: -, 5: 4, 6: 3, 7: -, 8: 4, 9: -,
             10: -, 11: -}
  buckets = {3: [3,6], 4: [1,5,8]}   mdeg = 3
round 3: mdeg = 3
  eliminate 3 (degree 3, weight 1 -> 1)   merged = none   evicted = 1, 5, 8
  eliminate 6 (degree 3, weight 1 -> 1)   merged = none   evicted = 1, 5, 8
  batch of 2: [3, 6]    refreshed = 1, 5, 8
  degrees = {0: -, 1: 2, 2: -, 3: -, 4: -, 5: 2, 6: -, 7: -, 8: 2, 9: -,
             10: -, 11: -}
  buckets = {2: [1,5,8]}   mdeg = 2
round 4: mdeg = 2
  eliminate 1 (degree 2, weight 1 -> 3)   merged = 5, 8   evicted = none
  batch of 1: [1]    refreshed = none
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -, 8: -, 9: -,
             10: -, 11: -}
  buckets = {}   mdeg = 2
order = [11, 4, 7, 10, 0, 2, 9, 3, 6, 1, 5, 8]
rounds = 5, pivots = 10, average batch = 2.0
fill = 7,  nnz(L) = 37 against nnz(tril A) = 30
degree computations: 26
```

Ten pivots in five rounds, and 26 degree computations against the 34 that 5.8 needed for the same
graph. Round 1 is the clearest instance: pivots 4, 7 and 10 are mutually non-adjacent, and vertex
6 is evicted twice, once by 7 and once by 10, then refreshed once.

**And here the ordering finally moves.** Every version from 5.1 through 5.10 returns

```
[11, 4, 7, 10, 0, 2, 3, 1, 5, 6, 8, 9]
```

on `graph3`, and this one returns

```
[11, 4, 7, 10, 0, 2, 9, 3, 6, 1, 5, 8]
```

They agree for six positions and then diverge, where the earlier versions take 3 and this one
takes 9. Both are vertices of degree 3, which is the point: the batch did not choose a worse
vertex, it chose a different vertex of the same degree, because 3 had been evicted by an earlier
pivot of the round and was not available to be chosen. The fill is 7 either way.

That single line is the taxonomy at the end of 5.10 made concrete. Mass elimination changed the
ordering without changing the fill; this changes the ordering and puts the fill in play. On
`graph3` the wager happens to cost nothing, and on the graphs of 5.11's table it costs a fraction
of a percent, in either direction.


**graph4.** Three rounds for six pivots, and the ordering diverges from every earlier layer at
position three.

```
start
  degrees = {0: 4, 1: 4, 2: 3, 3: 5, 4: 3, 5: 3, 6: 3, 7: 3}
  buckets = {3: [2,4,5,6,7], 4: [0,1], 5: [3]}   mdeg = 3
round 0: mdeg = 3
  eliminate 2 (degree 3, weight 1 -> 1)   merged = none   evicted = 0, 3, 5
  eliminate 4 (degree 3, weight 1 -> 1)   merged = none   evicted = 0, 1, 5
  eliminate 6 (degree 3, weight 1 -> 1)   merged = none   evicted = 1, 3, 5
  eliminate 7 (degree 3, weight 1 -> 1)   merged = none   evicted = 0, 1, 3
  batch of 4: [2, 4, 6, 7]    refreshed = 0, 1, 3, 5
  degrees = {0: 3, 1: 3, 2: -, 3: 3, 4: -, 5: 3, 6: -, 7: -}
  buckets = {3: [0,1,3,5]}   mdeg = 3
round 1: mdeg = 3
  eliminate 0 (degree 3, weight 1 -> 1)   merged = none   evicted = 1, 3, 5
  batch of 1: [0]    refreshed = 1, 3, 5
  degrees = {0: -, 1: 2, 2: -, 3: 2, 4: -, 5: 2, 6: -, 7: -}
  buckets = {2: [1,3,5]}   mdeg = 2
round 2: mdeg = 2
  eliminate 1 (degree 2, weight 1 -> 3)   merged = 3, 5   evicted = none
  batch of 1: [1]    refreshed = none
  degrees = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -}
  buckets = {}   mdeg = 2
order = [2, 4, 6, 7, 0, 1, 3, 5]
rounds = 3, pivots = 6
```

### 5.13 Approximate degree, which is the idea in AMD

MMD made the exact degree affordable by computing it rarely (5.11). AMD takes the other route: it
computes a degree at every single pivot, and makes that affordable by refusing to compute the
degree exactly at all.

**AMD does not batch**, which is worth saying outright because the code invites the opposite
reading. Its main loop is `while (nel < n)` with one pivot chosen per iteration and the degrees
updated before the next, and multiple elimination appears nowhere. What does appear is
`nel += nvpiv`, which advances the count by the pivot's *weight*, so a single pivot numbers many
columns at once. That is mass elimination, not multiple elimination: many columns from one pivot,
grouped by identity rather than by independence, exactly the pair 5.11 separates. The two codes
share no trick beyond the common ground; each takes one branch of the fork and neither takes the
other's.

The exact external degree of `i` is

```
degree[i] = | A_i  union  ( union of L_e over e in E_i ) |    minus i
```

a set union per vertex per step. That is the bottleneck of true minimum degree, and it is what
AMD refuses to pay.

**Bound it instead.** Sum the elements' contributions rather than uniting them:

```
degree[i] <= min( n - k,                                  # nothing can exceed what remains
                  degree_old[i] + |L \ i|,                # it can only grow by the new element
                  |A_i \ L| + |L \ i| + sum |L_e \ L| )   # over e in E_i, e != the new element
```

The third line is the approximation. It **overcounts**, because two elements may overlap outside
`L` and that overlap is counted twice. So it is an *upper bound* on the degree, not the degree.

**And the payoff is in one word: reuse.** The quantity `|L_e \ L|` depends only on the element
`e`, not on the vertex `i`. It is computed **once per element** and read by every vertex that
touches `e`. The exact degree would need a union *per vertex*.

The reason the exact degree cannot be shared the same way is worth stating, because it is the
whole of the matter: **a union is not decomposable and a sum is.** There is no way to split
`|A_i union L_e union L_f|` into a part that depends on `e` and a part that depends on `f`, since
the answer turns on how much `L_e` and `L_f` overlap, and the overlap that matters is different
for every vertex. So the union has to be recomputed for each `i`, and that is the per-vertex cost
being paid. Adding `|L_e \ L|` and `|L_f \ L|` instead asks a question about each element alone,
which is why one answer serves every vertex that touches it. The cost per vertex falls from a
walk over its elements' members to one addition per element. That is the entire speed argument,
and it says where the gain lands: the ratio between the two grows with element size, which is to
say with fill, which is to say precisely where the ordering is expensive.

**The other two terms are the brake**, and they are easy to read past as bookkeeping. Both are
exact and both are cheap: `n - k` counts what is left to eliminate, and `degree_old[i] + |L \ i|`
says a degree can only have grown by the new element. Neither needs a union. Taking the minimum
of all three lets the loose term be used where it is good and discarded where it is not, so the
overcounting cannot accumulate over a run. That clamp is load-bearing rather than defensive: the
combination result below turns on it, since a bound that is re-tightened at every pivot stays
useful and one that is not degrades quickly.

The algorithm whole:

```
amd(A) -> elimination order:

    for each vertex i:
        A_i    = the neighbors of i in A
        E_i    = {}                        # no elements yet
        w_i    = 1                         # supervariable weight
        degree[i] = |A_i|                  # degree

    while vertices remain:

        # PICK
        p = the live vertex of least degree

        # FORM THE ELEMENT: absorb every element p touches
        L = A_p
        for e in E_p:
            L = L union L_e
            kill e                         # element absorption
        L = L \ {p}
        L_p = L                            # p is now an element
        eliminate p

        # APPROXIMATE THE NEW DEGREES
        for each element e reachable from L:
            we = |L_e \ L|                 # ONCE per element, not per vertex
            if we == 0:                    # AGGRESSIVE ABSORPTION
                kill e                     # e lies wholly inside L, so it is dead

        for i in L:                        # ONE pass here, for legibility. The
            degree[i] = min( remaining,    # vendored code splits this in two;
                                           # see the note below the algorithm.
                             degree[i] + |L| - w_i,
                             |A_i \ L| + |L| + sum of we over e in E_i )

            if degree[i] == 0:             # MASS ELIMINATION
                w_p = w_p + w_i            # i joins the pivot's supervariable
                |L| = |L| - w_i            # ... and so LEAVES the new element,
                                           # which every later i in this same
                                           # pass must see. Easy to omit.
                eliminate i with p
                                           # (the degree[i] == 0 test is
                                           # equivalent to the structural one
                                           # only BECAUSE of the aggressive
                                           # absorption above)
            else:
                h_i = hash(A_i, E_i)       # for the merge below

        # ABSORB INDISTINGUISHABLE SUPERVARIABLES
        for i, j in one hash bucket, with identical patterns:
            w_i = w_i + w_j
            kill j                         # j eliminates whenever i does

    expand the supervariables back into individual vertices
```

Five mechanisms, and each is a solution to one of the two problems of 5.1:

| mechanism | solves |
|---|---|
| quotient graph, element absorption | the graph growing |
| approximate degree | the degrees being expensive |
| supervariables, found by hash | both, by shrinking the graph |
| mass elimination | both, likewise |
| aggressive absorption | the graph growing, more of it |

Aggressive absorption is the one addition beyond what 5.1 to 5.10 built, and it is nearly free.
Ordinary absorption kills the elements the *pivot* touched. Aggressive absorption also kills any
element whose members all lie inside the new one, which is detected by `|L_e \ L| == 0`, a
quantity the degree bound has just computed anyway. The vendored code makes it optional
(`Control[AMD_AGGRESSIVE]`), and the pseudocode above assumes it is on: the `degree[i] == 0` test
for mass elimination is equivalent to the structural test only when it is.

**The two codes give up different kinds of thing**, and the distinction matters more than the
sizes involved. MMD's pivots are exact. A batch never chooses a vertex that is not of true minimum
degree, so its ordering moves only because the tie among equally-minimal vertices is broken
differently (5.11). AMD's pivots are not exact. An overcounted bound can make the true minimum
look worse than it is, so the wrong vertex can be chosen outright. MMD perturbs; AMD can be wrong.

That the second is not obviously worse in practice is the striking result of the AMD paper, and it
survives our own measurement. On nine graphs, comparing the prototype against exact minimum degree
on the same base:

```
                     exact work   bound work   ratio   bound loose   fill against exact
2D grid 16x16              3800         1493    2.5x    464 / 1081        -1.84%
3D grid 6x6x6              7206         2602    2.8x    816 / 1321        +0.47%
2D grid 12x12              1823          731    2.5x    209 / 571         +1.17%
random graphs                                2.4-3.0x                  0.00 to +0.26%
overall                                                                   -0.15%
```

The middle column is the mechanism: element members walked once per *element* against once per
*vertex*, and the ratio grows with element size, which is to say with fill, which is to say
precisely where the cost is. The bound is genuinely loose, on a third to a half of the vertices it
touches, so this is a real approximation and not a formality. And the fill came out slightly
*better* than exact on these graphs, which is the same two-sided noise 5.11 found: choosing
greedily on a bound is a different heuristic, not a degraded one, and its errors are not
systematically in the bad direction.

**One line in that loop is easy to get wrong**, and it is worth flagging because we got it wrong.
When `i` is mass-eliminated it joins the pivot's supervariable, so it is no longer *outside* the
new element and must stop contributing to `|L|`. The vendored code does this at `degme -= nvi`;
our prototype originally did not, having computed `|L|` once before the loop began. The effect is
easy to miss: it changes nothing on any of the four graphs of 5.14 and nothing on the grids, and
it surfaced only on a five-vertex graph where a bound came out one too large. It matters because
the loop is a single pass, so a mass elimination at one `i` is seen by every later `i` in the same
pass, and a bound is only worth having if it bounds the right quantity.

**Can the two mechanisms be combined?** They look orthogonal: one decides how often to refresh,
the other how to compute a refresh, and nothing obviously forbids doing both. The reason they are
not orthogonal is in the formula above. Every term of the bound subtracts `L`, the element just
formed, and that subtraction is what makes it tight rather than trivial. It presumes there is
exactly *one* new element. Batch the pivots and several form at once, a vertex can belong to more
than one of them, and the bound then needs their union, which is the object the approximation
exists to avoid computing.

Running one implementation in all four configurations, so that only the two mechanisms vary:

```
                              work      fill
exact  + one pivot   (5.9)    1.00x    +0.00%
exact  + batch       (5.11)   0.94x    +0.38%
bound  + one pivot   (5.13)   0.14x    -0.98%
bound  + batch                0.16x   +14.67%
```

Combining is worse on both axes: fifteen percent more fill than exact minimum degree, and *more*
work than the bound alone rather than less.

One mechanism explains both failures. The bound stays tight because of its `min` against
`degree_old[i] + |L \ i|`, where `degree_old[i]` is itself a bound from the last time `i` was
touched. With one pivot per refresh that chain is re-tightened at every step, since every reached
vertex is recomputed at once. Batch, and a vertex may be touched by three pivots and refreshed
once, so bounds accumulate on bounds instead of being reset. Looser bounds then choose worse
pivots, worse pivots make more fill, more fill makes larger elements, and larger elements are why
the work rises too. So multiple elimination works by *delaying* refreshes and the approximate
degree depends on refreshes being *frequent*, which is a genuine conflict and not a coincidence of
implementation.

**This result is empirical, and worth labeling as such.** It is one implementation of one
reasonable generalization of the bound, on seven small graphs. The compounding argument explains
the numbers but does not prove they must come out that way, and a scheme that re-tightened bounds
inside a batch might do better, though it would have to compute per-vertex unions over the batch's
elements, which is the expense being avoided. What the measurement does establish is that AMD's
choice to refresh at every pivot is load-bearing rather than incidental, which is the thing worth
knowing when reading the code.

**What the prototype implements, and what it does not.** The pseudocode above describes the
vendored routine. Our `experiments/ordering/amd` is a subset of it, and the gap is recorded here
rather than only in the file header, so that it can be closed deliberately. The same accounting
as 5.11's, for the same reason.

Present: the quotient graph with element absorption, the three-way bound above, aggressive
absorption, mass elimination, hash supervariable detection, and degree buckets. That is enough to
reproduce the behavior this section describes, which is why the measurements above come from it.

Absent, in two categories. Three are **algorithmic**, and closing them would change the output:

- **The update runs in one pass here and two there.** This is a difference in structure rather
  than a missing feature, and it is the important one. `Amd.cpp` computes, in a first loop, a
  bound that deliberately excludes the new element, performing the mass eliminations as it goes
  and shrinking `|L|` at each one; a *second* loop then adds the **final** `|L|` to every
  survivor. We fold both into a single pass, so a vertex handled early sees a larger `|L|` than
  one handled late, where the vendored code gives them all the same value. Since `|L|` only
  shrinks, our early vertices get looser bounds. Measured on eleven graphs, the two forms give
  different orderings on four, all of them grids, with fill moving a few tenths of a percent in
  either direction. One pass is easier to read, which is why the pseudocode above shows it; two
  is faithful.
- **Dense-row handling** (`amd_preprocess`, with the `AMD_DENSE` control) removes rows above a
  density threshold from the ordering entirely and places them last. It changes the result on any
  matrix that has such rows, and it is the first thing to add if the prototype is ever run on real
  problems.
- **The postorder** (`amd_postorder`, `amd_post_tree`) reorders the output within the assembly
  tree. It leaves the fill alone but changes the permutation, and so changes everything downstream
  that reads one, our elimination forest included.

The rest are **packaging**, and closing them would change nothing:

- **`amd_aat`** forms `A + A'`, so that unsymmetric input can be ordered on its symmetric pattern.
  Our prototypes take a symmetric graph directly.
- **`amd_valid`** checks the input, and **`Control` and `Info`** carry parameters in and
  statistics out, `aggressive` among them, where we hardcode the choice.
- **Workspace compression** reclaims the index array when it fills, a memory strategy for a
  flat-array representation with no counterpart in a prototype built on sets.

### 5.14 A worked example: approximate degree

Two lines carry this section. `bounds` is what the algorithm believes and acts on; `exact` is what
the degree really is, computed alongside purely so the two can be compared. In a production
implementation the second line would not exist, which is the whole point of the first.

**graph1.**

```
start: no elements yet, so nothing to approximate over
  bounds  = {0: 2, 1: 2, 2: 2, 3: 2}   (exact)
  exact   = {0: 2, 1: 2, 2: 2, 3: 2}
  buckets = {2: [0,1,2,3]}   mdeg = 2
eliminate 0 (bound 2, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: -, 1: 2, 2: 2, 3: 2}
  exact   = {0: -, 1: 2, 2: 2, 3: 2}
  buckets = {2: [1,2,3]}   mdeg = 2
eliminate 1 (bound 2, weight -> 3)   merged = 2, 3   absorbed = none
  bounds  = {0: -, 1: -, 2: -, 3: -}
  exact   = {0: -, 1: -, 2: -, 3: -}
  buckets = {}   mdeg = 2
order = [0, 1, 2, 3]
fill = 1,  bound was loose 0 of 2
element members an exact degree walks: 0;  element reads used: 0
```

**graph2.**

```
start: no elements yet, so nothing to approximate over
  bounds  = {0: 2, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}   (exact)
  exact   = {0: 2, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}
  buckets = {2: [0,1,2,5], 3: [3,4]}   mdeg = 2
eliminate 0 (bound 2, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: -, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}
  exact   = {0: -, 1: 2, 2: 2, 3: 3, 4: 3, 5: 2}
  buckets = {2: [1,2,5], 3: [3,4]}   mdeg = 2
eliminate 1 (bound 2, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: -, 1: -, 2: 2, 3: 3, 4: 3, 5: 2}
  exact   = {0: -, 1: -, 2: 2, 3: 3, 4: 3, 5: 2}
  buckets = {2: [2,5], 3: [3,4]}   mdeg = 2
eliminate 2 (bound 2, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: -, 1: -, 2: -, 3: 2, 4: 2, 5: 2}
  exact   = {0: -, 1: -, 2: -, 3: 2, 4: 2, 5: 2}
  buckets = {2: [3,4,5]}   mdeg = 2
eliminate 3 (bound 2, weight -> 3)   merged = 4, 5   absorbed = none
  bounds  = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -}
  exact   = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -}
  buckets = {}   mdeg = 2
order = [0, 1, 2, 3, 4, 5]
fill = 2,  bound was loose 0 of 6
element members an exact degree walks: 0;  element reads used: 0
```

**graph3.**

```
start: no elements yet, so nothing to approximate over
  bounds  = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 3, 11: 1}   (exact)
  exact   = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 3, 11: 1}
  buckets = {1: [11], 2: [4,7], 3: [0,2,3,9,10], 4: [1,5,6,8]}   mdeg = 1
eliminate 11 (bound 1, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  exact   = {0: 3, 1: 4, 2: 3, 3: 3, 4: 2, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  buckets = {2: [4,7,10], 3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 1
eliminate 4 (bound 2, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  exact   = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: 2, 8: 4, 9: 3,
             10: 2, 11: -}
  buckets = {2: [7,10], 3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 2
eliminate 7 (bound 2, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: 2, 11: -}
  exact   = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: 2, 11: -}
  buckets = {2: [10], 3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 2
eliminate 10 (bound 2, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  exact   = {0: 3, 1: 4, 2: 3, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [0,2,3,9], 4: [1,5,6,8]}   mdeg = 2
eliminate 0 (bound 3, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: -, 1: 4, 2: 3, 3: 4, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  exact   = {0: -, 1: 4, 2: 3, 3: 4, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [2,9], 4: [1,3,5,6,8]}   mdeg = 3
eliminate 2 (bound 3, weight -> 1)   merged = none   absorbed = c4
  bounds  = {0: -, 1: 4, 2: -, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  exact   = {0: -, 1: 4, 2: -, 3: 3, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [3,9], 4: [1,5,6,8]}   mdeg = 3
eliminate 3 (bound 3, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: -, 1: 3, 2: -, 3: -, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  exact   = {0: -, 1: 3, 2: -, 3: -, 4: -, 5: 4, 6: 4, 7: -, 8: 4, 9: 3,
             10: -, 11: -}
  buckets = {3: [1,9], 4: [5,6,8]}   mdeg = 3
eliminate 1 (bound 3, weight -> 1)   merged = none   absorbed = c7
  bounds  = {0: -, 1: -, 2: -, 3: -, 4: -, 5: 3, 6: 3, 7: -, 8: 3, 9: 3,
             10: -, 11: -}
  exact   = {0: -, 1: -, 2: -, 3: -, 4: -, 5: 3, 6: 3, 7: -, 8: 3, 9: 3,
             10: -, 11: -}
  buckets = {3: [5,6,8,9]}   mdeg = 3
eliminate 5 (bound 3, weight -> 4)   merged = 6, 8, 9   absorbed = c10
  bounds  = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -, 8: -, 9: -,
             10: -, 11: -}
  exact   = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -, 8: -, 9: -,
             10: -, 11: -}
  buckets = {}   mdeg = 3
order = [11, 4, 7, 10, 0, 2, 3, 1, 5, 6, 8, 9]
fill = 7,  bound was loose 0 of 19
element members an exact degree walks: 19;  element reads used: 7
```

**And on all three, the two lines never differ.** The bound is not merely close, it is exact at
every step, and the last line of each trace says so: loose 0 of 2, 0 of 6, 0 of 19. A reader who
stopped here would have watched the whole of AMD's machinery, the elements, the absorptions, the
mass eliminations, and never once seen it approximate anything.

That is not an accident of these particular graphs, and the reason is the useful part. The sum
overcounts only when a vertex belongs to **two elements that overlap outside the new one**. That
needs enough eliminations to have created several live elements, and enough fill for them to
intersect. Small sparse graphs do not get there: a vertex usually has at most one element, and a
sum with one term is a union with one term. Checking exhaustively, **no connected graph on five or
six vertices has a loose bound anywhere in its run**, and a sample of thirty thousand on seven
vertices found none either.

So the fourth graph exists for this section. Eight vertices, fourteen edges, the smallest we
found on which the bound is ever wrong.

**graph4.**

```
start: no elements yet, so nothing to approximate over
  bounds  = {0: 4, 1: 4, 2: 3, 3: 5, 4: 3, 5: 3, 6: 3, 7: 3}   (exact)
  exact   = {0: 4, 1: 4, 2: 3, 3: 5, 4: 3, 5: 3, 6: 3, 7: 3}
  buckets = {3: [2,4,5,6,7], 4: [0,1], 5: [3]}   mdeg = 3
eliminate 2 (bound 3, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: 4, 1: 4, 2: -, 3: 5, 4: 3, 5: 4, 6: 3, 7: 3}
  exact   = {0: 4, 1: 4, 2: -, 3: 5, 4: 3, 5: 4, 6: 3, 7: 3}
  buckets = {3: [4,6,7], 4: [0,1,5], 5: [3]}   mdeg = 3
eliminate 4 (bound 3, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: 4, 1: 5, 2: -, 3: 5, 4: -, 5: 4, 6: 3, 7: 3}
  exact   = {0: 4, 1: 5, 2: -, 3: 5, 4: -, 5: 4, 6: 3, 7: 3}
  buckets = {3: [6,7], 4: [0,5], 5: [1,3]}   mdeg = 3
eliminate 6 (bound 3, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: 4, 1: 4, 2: -, 3: 4, 4: -, 5: 4, 6: -, 7: 3}
  exact   = {0: 4, 1: 4, 2: -, 3: 4, 4: -, 5: 3, 6: -, 7: 3}
  buckets = {3: [7], 4: [0,1,3,5]}   mdeg = 3
eliminate 7 (bound 3, weight -> 1)   merged = none   absorbed = none
  bounds  = {0: 3, 1: 3, 2: -, 3: 3, 4: -, 5: 4, 6: -, 7: -}
  exact   = {0: 3, 1: 3, 2: -, 3: 3, 4: -, 5: 3, 6: -, 7: -}
  buckets = {3: [0,1,3], 4: [5]}   mdeg = 3
eliminate 0 (bound 3, weight -> 4)   merged = 1, 3, 5   absorbed = c6
  bounds  = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -}
  exact   = {0: -, 1: -, 2: -, 3: -, 4: -, 5: -, 6: -, 7: -}
  buckets = {}   mdeg = 3
order = [2, 4, 6, 7, 0, 1, 3, 5]
fill = 4,  bound was loose 1 of 12
element members an exact degree walks: 21;  element reads used: 12
```

The divergence is at the third elimination. After `6` goes, vertex `5` is credited with a bound of
4 while its true degree is 3:

```
  bounds  = {0: 4, 1: 4, 2: -, 3: 4, 4: -, 5: 4, 6: -, 7: 3}
  exact   = {0: 4, 1: 4, 2: -, 3: 4, 4: -, 5: 3, 6: -, 7: 3} ^
```

Vertex `5` sits in two elements by then, and they share a vertex outside the element just formed.
The union counts it once and the sum counts it twice, which is exactly the overcount the third
line of the bound admits to. One vertex, one step, one too many.

Note what does *not* go wrong. The bound is never too small, which is what makes it safe to use:
an upper bound can make a vertex look worse than it is and cost us a good pivot, but it can never
make one look better than it is and tempt us into a bad one. And on this graph the error changes
nothing that matters, since `amd` returns the same fill as every earlier layer.

**Which is the section in miniature.** The approximation is invisible on small graphs, appears
only once fill has accumulated, and even then usually costs nothing. Its looseness rate tracks
density directly: 0 percent here on `graph1` to `graph3`, about 8 percent on `graph4`, 37 percent
on a 2D grid, 62 percent on a 3D one. And the *saving* tracks the same quantity, since both grow
with element size. The approximation bites hardest exactly where it saves the most, which is why
it is worth making.

### 5.15 What is intricate, and what is merely old

We vendor two codes, and both are long for the same two reasons, an intricate algorithm and a
fifty-year-old encoding, in different proportions. The split matters because the two causes have
opposite implications for a rewrite: one is worth preserving and the other only worth deleting.
Take AMD first, the larger.

The vendored SuiteSparse AMD is about 1200 lines of code (plus 990 of comment) in essentially one
function. It is worth being exact about where that goes.

**Inherent, and no style will remove it:**

- the quotient graph, with elements, absorption, and the variable/element split
- the approximate degree bound, and the `|L_e \ L|` computation that makes it cheap
- supervariable hashing, and the pattern comparison that confirms a match
- mass elimination

**Archaeological, and would evaporate in a modern encoding:**

- **One flat `int` array** (`Iw`) holding every adjacency list, every element pattern, and the
  free space, interleaved, indexed by offsets in a second array. This is a Fortran constraint: no
  dynamic allocation, no structs, one workspace handed in by the caller.
- **A garbage collector**, several hundred lines, to compact that array when absorbed lists have
  left it full of holes. This exists *only because of the flat array*. Give each vertex its own
  vector and there is nothing to collect.
- **Parallel arrays instead of a record**: `Pe`, `Len`, `Nv`, `Elen`, `Degree`, `W`, `Head`,
  `Next`, `Last`, all indexed by vertex, all conceptually fields of one struct. A COMMON block,
  transliterated.
- **Sentinel-encoded state**: `Elen[i] < 0` means "i is an element", `Nv[i] < 0` means "absorbed",
  and so on. Negative numbers as tags, because there was nowhere else to put a flag.
- **One 1800-line function**, because subroutine calls were expensive and everything shared the
  workspace anyway.

So the code is long *both* because the mechanisms are intricate *and* because the encoding is
fifty years old. Only the first is worth preserving. A legible version keeps every mechanism above
and loses every item in the second list, and the honest expectation is that it lands at a few
hundred lines and can be read.

MMD sits at a different point on both axes, and the contrast is instructive. Its inherent core is
if anything lighter: it keeps the quotient graph and mass elimination but drops the two subtlest
mechanisms above, the approximate degree bound and supervariable hashing, because it computes the
exact degree and detects indistinguishable vertices structurally during the batch update rather
than by hash. What it adds in their place, multiple elimination, is a batching loop, not a data
structure. The encoding, though, is more concentrated, not less: some 160 lines of algorithm with
none of AMD's near-thousand lines of comment, the Sparspak transliteration entire, numeric goto
labels, the pointer decrement that fakes one-based indexing, sentinel-tagged state, and variables
hoisted to a function's top only to dodge a goto crossing their initialization. The one
archaeological burden it escapes is the garbage collector: MMD compacts each adjacency list in
place as it eliminates, so it needs no separate compaction pass, and that absence is most of why it
is a couple hundred lines rather than two thousand. Between a simpler algorithm to keep and a
shorter, if denser, text to replace, MMD is the more approachable of the two to rewrite first.

The risk, which is real: a rewrite must produce **the same permutation**, or the fill changes and
every downstream measurement drifts. That equivalence is the thing that makes it a piece of work
rather than a stylistic pass, and it is also the thing that makes it testable: run both, compare
the permutations, and the vendored code is its own oracle.

## 6. Multifrontal factorization

Section 1.4 gave two schedules for the same arithmetic, left-looking and right-looking, and
Section 4.6 lifted both to supernodes. This section is the third, and it is different in kind. Left
and right-looking disagree about *when* an update is applied; multifrontal disagrees about *where
it lives in between*. That one change is what carried the method out of core in the 1970s and what
makes it the natural shape for parallelism now, and it is also what makes it cost more memory than
either of the other two. All three are in the same file: our `factorStaticMultifrontal` and
`factorDynamicMultifrontal` sit beside the left- and right-looking drivers and reuse their kernels
unchanged.

The section is arranged around that trade. 6.1 and 6.2 say what the method is and follow one
example through. 6.3 and 6.4 say why it was invented and why it survived its original purpose.
6.5 through 6.8 are about the price, which is the extra storage the method needs, and about the
one classical technique for lowering it.

### 6.1 The frontal matrix, and an update that travels one edge

Supernode `K` owns a set of columns, its **front**, and its index set is those columns followed by
the rows below them, `update(K)` in the vocabulary of 2.3. Gather the whole of that into one dense
block and it is the **frontal matrix** of `K`. Factor its front columns and two things come out:
the finished columns of `C` or `L`, which go to the factor and are never touched again, and a
Schur complement over `update(K)` alone, square and symmetric, which is not part of the factor and
must reach `K`'s ancestors somehow. That square is the **update matrix** of `K`, and it is the
object the whole method is built around. We write `U(K)` for its size in entries, `|update(K)|`
squared.

**A note on the frontal matrix's shape, because implementations differ and the difference is
visible in the memory.** The classical presentation makes it the full `|indexSet(K)|` square, with
the update matrix its trailing submatrix, so the whole square is what goes on the stack and comes
back off. Ours is the `|indexSet(K)|`-by-`|front(K)|` rectangle instead: the columns being factored
and the rows below them, which is precisely the block the factor already needs to store, with the
update matrix allocated separately. Nothing algorithmic turns on the choice, but the scratch does.
A stack of whole squares carries `|indexSet|` squared per supernode; a stack of update matrices
carries `|update|` squared, and the front columns are not duplicated because they live in the
factor. Everything counted in 6.5 is the second arrangement.

The three traversals differ in what happens to it.

Left-looking never forms it. When `K` is reached, `K` asks each descendant `J` that reaches into it
for the part of `J`'s update that lands in `K`, one rectangle per `(J, K)` pair, and each rectangle
is built, assembled and destroyed on the spot. Right-looking is the mirror: `J`, once factored,
walks its ancestors and pushes the corresponding rectangle into each. Both hold one rectangle at a
time, and neither keeps anything between supernodes.

Multifrontal forms the whole square once, and moves it exactly one edge:

```
J finishes  ->  its update matrix waits
                parent(J) is reached  ->  the update matrix is assembled into the parent's
                                          frontal matrix, then freed
```

Nothing is sent to a grandparent. What `J` owed an ancestor higher up is now inside the parent's
frontal matrix, and it travels on inside the parent's own update matrix when the parent finishes.
An update reaches a distant ancestor by riding up the tree one level at a time, in someone else's
luggage.

**The containment theorem is what permits this.** By 2.3, `update(J)` is a subset of the index set
of `parent(J)`. Every row and column of `J`'s update matrix therefore has a place in the parent's
frontal matrix, so the assembly is total: nothing is left over, and no entry has to be held back
for a further ancestor. Without containment the method would not close.

**The assembly is a scatter, not a copy.** `J`'s update matrix is indexed by `J`'s update rows,
which are a subset of the parent's index set but not a contiguous one, so each entry is routed by
where its global index falls in the parent. An entry whose row and column both land in the parent's
front columns goes into the part of the frontal matrix that is about to be factored. An entry that
lands in the parent's update rows goes into the parent's own update matrix instead, to travel on.
Both destinations are one dense block each, and the routing is a lookup through the global-to-local
map. The literature calls this operation *extend-add*; we call it **assembly**, and our
`assembleUpdateMatrix` is the routine that performs it.

### 6.2 A worked example

Take the forest

```
            K
           / \
          P   J
              |
              M
```

with these sizes, in columns:

```
snode   front   update      indexSet
  M       2       3         its 2 columns, then 3 rows shared with J
  J       2       2         its 2 columns, then 2 rows shared with K
  P       1       2         its 1 column,  then 2 rows shared with K
  K       3       0         its 3 columns, root, nothing below
```

Multifrontal walks the supernodes in an order that reaches every child before its parent, so
`M, J, P, K` will do. Following the update matrices:

```
M factored     frontal 5x2, factor 2 columns, leaves a 3x3 update matrix   U(M) = 9
                 live: M
J reached      assemble M's 3x3 into J, routed by index: some entries into the
                 columns J is about to factor, the rest into J's own update
                 matrix; free M
                 J factored, leaves a 2x2 update matrix                    U(J) = 4
                 live: J
P reached      nothing to assemble, P is a leaf
                 P factored, leaves a 2x2 update matrix                    U(P) = 4
                 live: J, P
K reached      assemble J's 2x2 and P's 2x2 into K, both entirely into K's own
                 columns since K is a root; free both
                 K factored, no update rows, nothing left                  U(K) = 0
                 live: nothing
```

Two things to notice, and both matter later.

`M`'s contribution reaches `K` without ever being addressed to `K`. It is assembled into `J`, and
whatever part of it belongs to `K` is carried inside `J`'s update matrix when `J` finishes. `M` and
`K` are never in contact.

And the peak is not at the end. At the moment `K` is reached, `J`'s and `P`'s update matrices are
both live and `K`'s frontal matrix is being built. Earlier, while `P` was being factored, `J`'s
update matrix was already waiting. That waiting is the extra storage the method costs, and 6.5
counts it.

### 6.3 Why it was invented: a file has only one cheap end

The method is older than the memory it now runs in. The frontal method (Irons, 1970) came out of
finite element work, where the matrix did not fit: assemble one element front at a time, factor
what can be factored, write it out, and keep only the active front in core. Duff and Reid (1983)
generalized one front to many, organized by the elimination tree, and the multifrontal method is
that generalization.

What makes the organization work out of core is a property of the lifetimes. A supernode's update
matrix is created when the supernode is factored and destroyed when its parent is reached. If the
supernodes are numbered so that each subtree occupies a contiguous run, then those lifetimes
**nest**: any two are either disjoint, or one contains the other. A set of nested intervals is
exactly what a stack holds, so the update matrices can live on one, pushed as they are created and
popped as they are consumed. That is why the literature calls the storage the **update stack**, and
why the peak is called the stack peak.

And a stack is what a file can be. A file is a linear sequence of bytes with one cheap end: it can
be appended to and truncated there, and neither anywhere else. Push is an append and pop is a
truncate, so a stack spills to a file without further machinery. Every alternative is worse.
Freeing a block in the middle means either compacting, which rewrites the tail on every free, or
tracking holes, which is a free-space allocator on disk written by hand, or giving each update
matrix its own file, which trades the problem for thousands of files with the open, close and
metadata cost of each, for blocks that live a few steps.

So the discipline was not a stylistic preference. Out of core, it is the only arrangement that
works, and the requirement comes from the medium rather than from anything about the
factorization.

**None of that applies in core.** Main memory is randomly addressable, so a block can be freed
wherever it sits and the peak is unchanged: what is live at any instant is fixed by the traversal
and the lifetimes, not by how the bytes were obtained. Our implementation keeps one update matrix
per supernode in a plain vector, allocated when the supernode is reached and freed when its parent
consumes it, and it holds exactly the same set at every moment a stack would. The word *stack*
survives in the name of the quantity, which is standard, and we keep it in prose for that reason;
the code does not claim it.

### 6.4 Why it survived: the shape matches parallelism

The out-of-core motivation is largely historical. What kept multifrontal central is that the same
organization is the one parallelism wants, and for reasons that have nothing to do with disks.

**The heavy arithmetic is concentrated.** Each supernode's numeric work happens inside one dense,
contiguous frontal matrix: a dense factorization of the front columns, a triangular solve, and a
symmetric update, each a single large BLAS-3 call at near peak throughput. The irregular part, the
scatter that assembles a child into its parent, is data movement and sits *outside* the
arithmetic. Left- and right-looking invert this. They chop the same update into many smaller
index-mapped rectangles, each fronted by a gather, so the same flops arrive in smaller kernels with
worse locality. Multifrontal localizes the arithmetic and pushes the sparsity out to the assembly;
the looking traversals smear the sparsity through the arithmetic.

**And the tree is a task graph.** Two kinds of parallelism live in it, and they hand off. Near the
leaves the forest is wide, so independent subtrees run concurrently, one task each: *tree
parallelism*, many small fronts at once. Near the root the tree has narrowed to a few supernodes
and there is little concurrency left, but the fronts are now large, so the parallelism moves
*inside* one front, the big dense kernels threaded or handed to an accelerator: *node
parallelism*. The bottom of the tree parallelizes across fronts, the top within a front. MUMPS is
the reference implementation of exactly this split, sorting fronts into three node types: type 1
handled by a single process, which is tree parallelism; type 2, a large interior front partitioned
by rows across a master and slaves, which is node parallelism in 1D; and type 3, the root,
distributed 2D over a process grid.

**The costs are real and worth stating plainly.** Multifrontal does slightly *more* arithmetic than
the looking traversals, and the surplus is all in assembly. The multiplicative work that forms the
factor is identical, but a contribution travelling from `M` to `K` through `J` is added twice,
once into `J` and once into `K`, where left-looking would add it to each ancestor once. A
contribution once formed is data the tree moves rather than recomputes, so no multiplication is
repeated; it is additions that are.

The larger cost is memory, and it is intrinsic. Concurrency raises the peak by itself: `N` branches
in flight hold `N` live sets where serial execution holds one. That is the price of the wall-clock
time, and no arrangement of storage avoids it. What softens it is that the multiplier applies where
blocks are cheap. Tree parallelism is available near the wide leaves, where update matrices are
small; near the root, where they are enormous, there is almost no tree parallelism left to exploit.
The expensive blocks stay singular.

### 6.5 The extra storage, and the formula that measures it

Left- and right-looking hold one rectangle at a time and fold it straight into the factor. They
never keep a standing Schur complement. Multifrontal keeps a set of update matrices alive, scratch
that is not part of the factor, from each supernode's turn until its parent's. That difference is
the whole memory cost of the method, and it is worth measuring rather than estimating.

Let `I` and `J` range over the children of `K`, with `I < J` meaning `I` is assembled before `J`,
and write `storage(X)` for `U(X)`, the entries one update matrix occupies. The storage a
supernode's whole subtree needs is

```
maxStorage(K) = max( max_J [ sum_{I<J} storage(I) + maxStorage(J) ],   term 1, while J runs
                     sum_J storage(J) + storage(K) )                   term 2, at K's assembly
```

**Term 1 is the worst moment while one of `K`'s subtrees is still running.** When `J` is running,
two things are resident: whatever `J`'s own subtree needs at its worst instant, which is
`maxStorage(J)` recursively, and the update matrices of its earlier siblings, which finished and
are waiting for `K`. It is a maximum over the children, and the sum `sum_{I<J}` is the only place
the assembly order enters the formula at all.

**Term 2 is the single moment at `K` itself**, when every child has finished, none has been freed,
and `K`'s own update matrix is allocated on top of them. It is a plain total, so no order can
change it. A leaf has no children, so both reduce to `storage(leaf)`, which is where the recursion
bottoms out.

The recursion has optimal substructure, which matters for 6.6: `maxStorage(J)` is a single number
summarizing everything below `J`, and `K` needs nothing else from inside `J`. In particular the
internal arrangement of `J`'s own subtree does not change how `J` and its siblings should be
ordered.

**How large is it in practice?** Measured on our implementation, as a percentage of the factor
itself, with the child order left as the ordering produced it:

```
                       |L|          update stack peak
grid2D 80x80         133261         42735    ( 32.1%)
grid3D 16x16x16      340513        357398    (105.0%)
band n=2000 bw=20     41980          1200    (  2.9%)
```

On banded problems it is negligible. On two-dimensional grids it is a third of the factor. On
three-dimensional grids it *exceeds* the factor, so peak working memory is roughly double what
`|L|` alone suggests, and the half that is pure scratch is the half a technique might reduce.

**Dynamic pivoting does not change any of this**, which is not obvious and is worth stating. An
update matrix is `|update(K)|` squared, and `update(K)` is fixed by the symbolic factorization.
Delaying a column moves it between the front and the delayed columns; it never touches the update
size. So the stack peak is identical for a dynamically pivoted factorization and a static one on
the same forest, and it can be predicted before any value is known.

### 6.6 Choosing the assembly order

Term 2 is fixed. Term 1 is not: it depends on the order in which `K`'s children are assembled,
through `sum_{I<J}`. A child that finishes early leaves its update matrix occupying space for as
long as its later siblings take to run, so the question is which children should pay that wait.

A priori this ranges over all orderings of the children. An exchange argument collapses it to a
sort. Take two adjacent children `X` and `Y`, and let `S` be the storage of everything assembled
before them. Their two positions contribute

```
P1  (X first) = max( maxStorage(X),  storage(X) + maxStorage(Y) )
P2  (Y first) = max( maxStorage(Y),  storage(Y) + maxStorage(X) )
```

with `S` dropped, since it appears in both, and everything after the pair untouched, since the
running sum grows by `storage(X) + storage(Y)` either way. Now suppose

```
maxStorage(X) - storage(X)  >=  maxStorage(Y) - storage(Y)
```

and write `T = storage(Y) + maxStorage(X)`, the second argument of `P2`, so that `P2 >= T`. Then
`maxStorage(X) <= T`, because `storage(Y) >= 0`, and `storage(X) + maxStorage(Y) <= T`, which is
the supposition rearranged. Both arguments of `P1` are at most `T`, so `P1 <= T <= P2`, and putting
`X` first is never worse.

Any out-of-order adjacent pair can therefore be swapped without loss, so the key induces a total
order and sorting by it is optimal with no search. This is Liu's rule (1986): **assemble the
children in decreasing order of `maxStorage(J) - storage(J)`.** The intuition reads off the key. A
child whose subtree needs a great deal and leaves little behind should run while little is waiting;
a child that needs little but leaves a bulky update matrix should run last, so that its leftover
waits as briefly as possible.

The cost is `O(m log m)` for a supernode with `m` children, and since the child counts sum to the
supernode count, `O(snodeSize log m)` for the forest, with `m` small in practice. It runs once,
during the symbolic phase, and it is structural, so it is reused by every numeric factorization
sharing the pattern. The `maxStorage` values it needs come free in the same upward pass, because
supernodes are visited in an order that reaches every child before its parent.

**The formula measures an order; it does not choose one.** That distinction is worth keeping,
because the two are separate operations in any implementation: a sort that picks the order and an
evaluation that computes `maxStorage(K)` under the order picked. Ours are
`sortForOptimalMultifrontal` and the loop at its end, run in that sequence.

### 6.7 A worked example: what the ordering buys

Take a root `K` with two children:

```
snode   maxStorage   storage      key = maxStorage - storage
  A        210          10                  200
  B        150         100                   50
  K          0           0                    0
```

`A` heads a deep subtree that needs 210 at its worst but leaves only 10 behind. `B` needs less, 150,
but leaves 100. Evaluating the formula both ways:

```
                     term 1                                      term 2      maxStorage(K)
A first    max( 210, 10 + 150 )  = 210            10 + 100 + 0 = 110              210
B first    max( 150, 100 + 210 ) = 310            100 + 10 + 0 = 110              310
```

Term 2 is identical, as it must be. The whole difference is term 1, and it is 100 entries, a third
of the worse arrangement. Liu's key ranks `A` at 200 and `B` at 50, so decreasing key puts `A`
first, which is the better column. The mechanism is visible in the arithmetic: putting `A` first
means only `A`'s small leftover of 10 is waiting while `B` runs, whereas putting `B` first leaves
`B`'s bulky 100 waiting through the whole of `A`'s expensive subtree.

**A parent all of whose children are leaves has no decision to make.** A leaf's `maxStorage` is its
own `storage`, so every key is zero, every order ties, and the arithmetic agrees: term 1 collapses
to `sum_J storage(J)`, which is order-invariant, and term 2 is that plus `storage(K)`. The ordering
does work only where children are themselves internal nodes with depth beneath them. This is why
banded matrices, whose forests are essentially paths, gain nothing from it while grids gain a third.

### 6.8 The labels have to carry the order

There is a second requirement, easy to miss, and without it the sort of 6.6 is inert.

A supernode numbering is constrained to be **topological**: `parent(K) > K`, which holds
automatically because a parent is the first index below the diagonal in a column of the factor.
That is all correctness needs, and every traversal is happy with it.

It is not all the storage analysis needs. The formula of 6.5 describes a tree traversal, in which a
subtree runs to completion before its siblings begin. An implementation that loops over supernodes
in index order realizes that traversal only if each subtree occupies a **contiguous** run of
indices, which topological order does not imply. The smallest counterexample is four nodes:

```
parent(0) = 2,  parent(1) = 3,  parent(2) = 3,  parent(3) = NIL

        3
       / \
      1   2
          |
          0
```

Every parent exceeds its child, so the numbering is topological and any traversal by index is
legitimate. But the subtree rooted at `2` is `{0, 2}`, and index `1` lies between them while
belonging to the other branch. Walking indices visits `0`, then leaves the branch for `1`, then
returns to `2`. The lifetimes tell the same story: `0`'s update matrix is live over `[0, 2]` and
`1`'s over `[1, 3]`, and those two intervals **overlap without nesting**. A block is held alive
across a branch it has nothing to do with, which is both more storage than the formula predicts and
the exact configuration no stack can hold.

The repair is a **postorder relabeling**: renumber the supernodes by a depth-first walk of the
forest, so that each subtree becomes a contiguous run ending at its own root. On the example this
exchanges the labels of `0` and `1`, after which the subtree at `2` is `{1, 2}` and every interval
nests. Note what does *not* move: the factor, the permutation, the columns. Only which label refers
to which supernode changes, so the relabeling is a pure permutation of the numbering.

**The relabeling does two jobs, and they are separable.** It makes the peak match the formula at
all, which matters even when no sort has been run, because a non-contiguous numbering costs storage
on its own. And it makes the sort's choice reachable, because a depth-first walk that follows the
sorted child links produces labels in that order, whereas the links alone are never read by a
driver that loops over indices.

Measured on our implementation, with the two steps run as a pair and compared against the same
forest without them:

```
                        without         with        saving
grid2D 80x80             42735         26425         38%
grid3D 16x16x16         357398        255433         28%
band n=2000 bw=20         1200          1200          0%
```

The three-dimensional case falls from above `|L|` to three quarters of it. Bands do not move,
having no sibling choice to make. And the sort alone, without the relabeling, moves nothing at all
while appearing to work: the links carry the chosen order and the traversal never reads them. This
is why the two belong together, and why our `ElmForestEngine` runs them as one operation,
`orderForMultifrontal`, rather than exposing either half.

## 7. Pivoting

Every section so far has assumed the pivot sequence is settled before any value is read. Cholesky
earns that assumption, static LDL asserts it, and dynamic LDL gives it up. This section is about
what the third one does instead: why a symmetric factorization cannot pivot the way LU does, what
the dense theory offers, what a sparse solver has to give back, and what our
`factorDynamicNonRootSupernode` actually implements, which is not quite what the literature it
is usually named after describes.

The one idea under all of it is that symmetry is a stability constraint, and everything in this
section is the response to it. Symmetry is a bargain: `A = L D L^H` factors and stores one triangle,
so the arithmetic and the storage are half of LU's. What it takes in exchange is the free parameter
LU leans on for stability. LU permutes rows independently of columns, `P A = L U`, so it can pull
the largest entry in a column straight to the pivot, and partial pivoting is cheap and never fails
because a row is free to move on its own. Symmetry ties the two sides together, `P A P^H`, so a row
swap is forced to be the same column swap, and only a diagonal entry can become a `1 x 1` pivot
(7.1); the off-diagonal one might have wanted is unreachable by any symmetric permutation. So LU
spends more arithmetic to buy easy stability, and LDL^H spends more pivoting machinery to buy back
the stability that symmetry took while keeping the arithmetic saving. The `2 x 2` blocks, the
Bunch-Parlett and Bunch-Kaufman search, the thresholds, and the delays are all that price. Sparse
tightens the same screw once more (7.4): the symmetric search cannot even run in full, only over the
fully summed part of the front, which is why delay had to be added on top. Oblio sits at the strict
end, symmetric and sparse at once, and the subsections below are the full accounting of what those
two constraints cost and how each is paid.

### 7.1 The diagonal is the only place to look

LU with partial pivoting swaps rows and leaves columns alone, computing `P A = L U`. Nothing
constrains the two sides of the matrix to agree, because a general `A` has no symmetry to lose.

A symmetric factorization does have one. `A = L D L^H` says the two triangles carry the same
information, which is why only one is stored and why the arithmetic is half of LU's. Keeping that
means applying the permutation on both sides, `P A P^H`, and a two-sided permutation maps diagonal
entries to diagonal positions and off-diagonal entries to off-diagonal positions. **A `1 x 1` pivot
therefore has to be a diagonal entry.** There is no symmetric permutation that brings an
off-diagonal entry to the pivot position, so the entire search happens along the diagonal.

For a positive definite `A` that costs nothing at all. Every diagonal entry of every Schur
complement is positive, by 1.6, and the largest entry in magnitude of a positive definite matrix
lies on its diagonal, since `|a_jk| < sqrt(a_jj * a_kk)`. The entry we would have wanted to pivot on
is already where it is needed. Cholesky therefore pivots for stability not at all, which is why
Sections 1 through 6 could pass over the subject in silence. The permuting those sections do care
about, the ordering of Section 5, is chosen for fill before any value is seen and is a different
question entirely.

Drop definiteness and the guarantee goes with it. Consider

```
A = [ 0  1 ]
    [ 1  0 ]
```

which is nonsingular, perfectly conditioned, and has no nonzero diagonal entry anywhere. No
symmetric permutation creates one, since such a permutation only reorders the diagonal it already
has. There is no sequence of `1 x 1` symmetric pivots that factors this matrix, and the obstruction
is structural rather than a matter of accuracy: it is not that the pivots would be small, it is that
they do not exist.

The repair is to let the pivot block itself be `2 x 2`, and to take the whole of `A` above as one
such block. Then `D` is block diagonal with `1 x 1` and `2 x 2` blocks and `L` is unit
block-lower-triangular, which is exactly the shape 0.9 and 0.10 set up, and the reason those
sections carry a `2 x 2` `D11` case rather than a scalar division. The inertia comes out correctly
along with it, by Sylvester's law: a `2 x 2` block with negative determinant contributes one
positive and one negative eigenvalue, and counting blocks is how such a factorization reports the
inertia of `A`. That is what Bunch and Kaufman's title is about.

**A catalog of the `2 x 2` cases.** Every question this section asks can be put to a symmetric
`2 x 2`, which is small enough to answer completely. Seven are worth keeping in mind, and the names
below are used only here.

```
identity                  negative identity         signed diagonal
[ 1  0 ]                  [ -1   0 ]                [ 1   0 ]
[ 0  1 ]                  [  0  -1 ]                [ 0  -1 ]

zero diagonal             zero diagonal, flipped    small diagonal, 0 < eps << 1
[ 0  1 ]                  [  0  -1 ]                [ eps  1 ]
[ 1  0 ]                  [ -1   0 ]                [  1   0 ]

rank one
[ 1  1 ]
[ 1  1 ]
```

Their spectra, with the inertia written `(positive, negative, zero)`:

```
                        eigenvalues       eigenvectors           eigenspaces   inertia
identity                +1, +1            every nonzero vector   one, dim 2    (2, 0, 0)
negative identity       -1, -1            every nonzero vector   one, dim 2    (0, 2, 0)
signed diagonal         +1, -1            (1,0) and (0,1)        two lines     (1, 1, 0)
zero diagonal           +1, -1            (1,1) and (1,-1)       two lines     (1, 1, 0)
zero diagonal, flipped  +1, -1            (1,-1) and (1,1)       two lines     (1, 1, 0)
small diagonal          see below         (lambda, 1) for each   two lines     (1, 1, 0)
rank one                2, 0              (1,1) and (1,-1)       two lines     (1, 0, 1)
```

The small-diagonal case has `trace = eps` and `det = -1`, so its eigenvalues are
`(eps +- sqrt(eps^2 + 4)) / 2`, which for small `eps` are `1 + eps/2` and `-1 + eps/2` to first
order. Its eigenvectors are `(lambda, 1)` for each `lambda`, so they sit near `(1,1)` and `(-1,1)`.
The two definite cases are the only ones with a repeated eigenvalue, and a repeated eigenvalue on a
symmetric matrix means a genuine plane of eigenvectors rather than a defective one: every nonzero
vector qualifies, and no basis is preferred.

What each case does to a pivot search is the part that matters here:

```
                        a 1x1 pivot on the diagonal?    what the case shows
identity                yes, either column              the easy case, accepted at any threshold
negative identity       yes, either column              LDL ignores the sign, Cholesky cannot
signed diagonal         yes, either column              indefinite does not imply a 2x2 is needed
zero diagonal           none exists                     why D carries 2x2 blocks at all
zero diagonal, flipped  none exists                     the sign of an off-diagonal reaches no test
small diagonal          yes, and it is a trap           what the threshold is for (7.5)
rank one                none, here or anywhere          a real singularity, caught at a root (7.6)
```

Three things follow, and they are the reason the catalog is worth having.

**Definiteness decides, not conditioning.** The first five all have condition number 1, as
well conditioned as a matrix can be, and their behavior runs the whole range from trivial to
impossible. Conditioning bounds how much an error already made can be amplified; it says nothing
about whether a pivot exists.

**Indefinite is necessary but nowhere near sufficient.** The signed diagonal and the zero diagonal
have the same inertia and opposite fates: one pivots twice with no search, the other cannot pivot at
all. The implication runs one way only. If `a_kk = 0` and some `a_ik` in its column is nonzero, the
principal `2 x 2` on `k` and `i` has determinant `-|a_ik|^2 < 0`, so it has one eigenvalue of each
sign, and interlacing forces the whole matrix to be indefinite. A zero on the diagonal beside a
nonzero off it therefore forces indefiniteness; indefiniteness forces nothing about the diagonal.

**The flipped pair is one matrix.** With the signature matrix `S = diag(1, -1)`, which is symmetric,
orthogonal and its own inverse, `S X S` carries the zero diagonal to its flipped twin. That single
identity is at once a similarity `S X S^-1`, which preserves eigenvalues, and a congruence
`S X S^T`, which preserves inertia by Sylvester's law, so both invariants are pinned by the same
equation and only the eigenvectors trade places. Our kernel reaches the same conclusion by a shorter
route: every quantity the acceptance tests read is a magnitude or a determinant, and the
off-diagonal enters the determinant squared, so its sign cannot reach any comparison.

That last observation generalizes into the classification the `2 x 2` acceptance test of 7.5 is
built on. For a symmetric `2 x 2` the trace is the sum of the eigenvalues and the determinant their
product, so two numbers settle everything:

```
det < 0                 indefinite,          inertia (1, 1, 0)
det > 0, trace > 0      positive definite,           (2, 0, 0)
det > 0, trace < 0      negative definite,           (0, 2, 0)
det = 0                 singular
```

which is why `readPivotBlock2x2` computes a determinant and the test compares against it. It also
settles what a `2 x 2` pivot contributes to the inertia: a block selected because its diagonal is
weak has both `|a|` and `|c|` below `|b|`, hence `|ac| < b^2` and `det < 0`, hence exactly one
eigenvalue of each sign. Counting the blocks is therefore counting the inertia.

**Why counting them is legitimate, and why it costs nothing.** The inertia of a symmetric `A` is
the triple `(n+, n-, n0)`, its counts of positive, negative and zero eigenvalues with multiplicity,
so the three sum to `n`. It carries a name because of **Sylvester's law of inertia**: the triple is
invariant under congruence, so that for every nonsingular `S`,

```
inertia( S A S^H )  =  inertia( A )
```

The eigenvalues themselves move freely under congruence; only how many lie on each side of zero is
pinned. And the law is sharp in both directions, since two real symmetric matrices are congruent
exactly when their inertias agree, which makes inertia a complete invariant of congruence rather
than merely a preserved one.

An `LDL^H` factorization is itself a congruence, with `S = L^-1`, so `inertia(A) = inertia(D)`, and
`D` is block diagonal with blocks small enough to classify by the rule just given:

```
1 x 1, positive                     (1, 0, 0)
1 x 1, negative                     (0, 1, 0)
1 x 1, zero                         (0, 0, 1)
2 x 2, det < 0                      (1, 1, 0)
2 x 2, det > 0, trace > 0 or < 0    (2, 0, 0) or (0, 2, 0)
```

Summing over the blocks gives the inertia of `A` exactly, with no arithmetic beyond the
factorization we wanted anyway and without going anywhere near an eigenvalue solver. The last line
does not arise in practice, for the reason given just above: a `2 x 2` block is selected only when
the diagonal is too weak to pivot on, and that is what drives its determinant negative.

This is why Bunch and Kaufman's title leads with the calculation of inertia rather than with the
factorization. For the problems that motivated it the inertia was the answer wanted and `LDL^H` was
the means, and the same holds wherever a count of negative eigenvalues is the question: confirming
that a stationary point is a minimum rather than a saddle, checking that a KKT system has the
inertia an interior point step requires, or counting a pencil's eigenvalues below a shift by
factoring `A - sigma B` and reading `n-` off the result.

`n0` is the weak entry of the three, and its weakness is numerical rather than mathematical. Exactly
zero is not a decidable property of a computed pivot, so a rank-deficient matrix does not announce
itself as a clean zero block. It appears instead as a column that cannot be accepted at any
threshold, which is the root failure of 7.6.

**We do not compute it.** Nothing in Oblio reports an inertia today, though both ingredients are
already in hand: `NumFactorDynamic::pivotType()` is public and says which columns open a `2 x 2`,
and the diagonal itself carries the signs. A count over the supernodes would be a short read-only
pass over a finished factor. It is recorded here as an observation about what the factorization
already knows, not as a plan.

**One boundary case, a single character away and entirely outside this framework.** Flip only the
lower off-diagonal of the zero-diagonal matrix and the result is skew-symmetric rather than
symmetric:

```
[ 0  -1 ]      eigenvalues +i and -i, no real eigenvectors, no inertia
[ 1   0 ]
```

Real symmetric matrices always have real eigenvalues and an orthogonal eigenbasis; the moment the
two off-diagonal entries disagree in sign rather than agreeing, that guarantee is gone. Nothing in
this section applies, and Oblio rejects such a matrix at input validation as structurally
unsymmetric.

### 7.2 Bunch-Parlett and Bunch-Kaufman

Two questions remain: which diagonal entry, or which pair, and how much work is spent deciding.

Bunch and Parlett (1971) answered the first well and the second badly. Their strategy searches the
entire trailing submatrix at every step, comparing the largest diagonal entry against the largest
off-diagonal entry and taking a `1 x 1` or a `2 x 2` accordingly. The growth bound is of the same
character as complete pivoting for LU, polynomial in `n` rather than exponential. The search costs
`O((n-k)^2)` comparisons at step `k` and so `O(n^3)` over the run, which is the order of the
arithmetic itself. For a factorization that costs `O(n^3/3)` multiplications, a search of that order
is not something anyone can afford.

Bunch and Kaufman (1977) answered the second. Their strategy examines **at most two columns** per
step, so `O(n)` comparisons at each and `O(n^2)` over the run, negligible beside the arithmetic. The
price is a weaker guarantee, and 7.3 is about which one.

The constant that makes it work is

```
alpha = (1 + sqrt(17)) / 8 = 0.6403882...
```

and it is not arbitrary. It is chosen so that the worst-case element growth of two consecutive
`1 x 1` steps equals that of one `2 x 2` step, which is what makes the two branches comparably safe
and what fixes the bound quoted at the end of this section.

With `k` the column being produced, the trailing submatrix `A(k:n, k:n)`, and `r` the candidate row:

```
bunch_kaufman(A, k) -> pivot size and swap:

    omega1 = max |A(i,k)|  over i = k+1 .. n          # largest subdiagonal entry of column k
    r      = the i attaining it

    if omega1 == 0:                                   # column already eliminated below the diagonal
        return 1x1 at k, no swap                      #   A(k,k) may still be zero: singular

    if |A(k,k)| >= alpha * omega1:                    # case 1: the diagonal stands on its own
        return 1x1 at k, no swap

    omegar = max |A(i,r)|  over i = k .. n, i != r    # largest off-diagonal entry of column r

    if |A(k,k)| * omegar >= alpha * omega1 * omega1:  # case 2: column r is no better, so keeping
        return 1x1 at k, no swap                      #   A(k,k) is bounded after all

    if |A(r,r)| >= alpha * omegar:                    # case 3: A(r,r) is a good pivot, bring it up
        return 1x1 at k, swap k <-> r

    return 2x2 at k, k+1, swap k+1 <-> r              # case 4: no acceptable 1x1 anywhere here
```

Every swap is symmetric, rows and columns together, so the factorization stays symmetric throughout.
Case 4 does not move column `k` at all: it stays as the first column of the block, and `r` comes
down to `k+1`.

Case 2 is the one that reads as arbitrary and is not. It says that `A(k,k)` failed the direct test
against its own column, but that column `r` is itself so large that switching to it would buy
nothing, so the growth from keeping `A(k,k)` is bounded anyway. Without it the algorithm would take
a `2 x 2` block every time the diagonal happened to be modest, which is both slower and, since a
`2 x 2` block consumes two columns, structurally different.

The elimination is the block algebra of 0.9 with the pivot block at one of its two sizes. For a
`1 x 1` pivot `d = A(k,k)`:

```
D(k)   = d
L(i,k) = A(i,k) / d                             for i = k+1 .. n
A(i,j) = A(i,j) - L(i,k) * d * L(j,k)           for i >= j > k
```

For a `2 x 2` pivot, write the block and its inverse:

```
E = [ a  b ]     det = a*c - b*b      E^-1 = (1/det) [  c  -b ]
    [ b  c ]                                         [ -b   a ]
```

and then

```
L(i, k:k+1) = A(i, k:k+1) * E^-1                for i = k+2 .. n
A(i,j)      = A(i,j) - L(i, k:k+1) * E * L(j, k:k+1)^H
```

The trailing update is a symmetric rank-2 update rather than rank-1, and it carries `E` in the
middle. That is the same `A22 - L21 D11 L21^H` shape 0.9 described, with the same consequence: no
BLAS routine computes it, which is why `formUpper` and `gemmLower` exist.

The growth bound that follows is roughly `(1 + 1/alpha)^(n-1)`, about `2.57^(n-1)`, against
`2^(n-1)` for partial pivoting on a general matrix. Worse in the bound, comparable in practice, and
the point is that a bound exists at all while symmetry is preserved.

### 7.3 What Bunch-Kaufman does not bound

Element growth in the trailing submatrix is bounded. **The entries of `L` are not.** They can be
arbitrarily large, and this is the known weakness of the method.

It matters because the backward error of a symmetric indefinite factorization is governed not by
`|A|` but by `|L| |D| |L^H|`, which a large multiplier can make much larger. A factorization can
therefore be the exact factorization of a matrix far from the one we handed it, while every
individual step looked acceptable. Growth in the trailing submatrix is bounded and the factor is
still poor.

Bounded Bunch-Kaufman, usually called **rook pivoting**, repairs it by continuing the search:
alternate along a row and then a column until an entry is found that is maximal in both, and accept
only then. The resulting multipliers are bounded, and in practice the extra searching is slight,
though in the worst case it grows to the order of the arithmetic and the cost advantage of 7.2 is
given back.

Ashcraft, Grimes and Lewis (1998) is the paper that made this concrete for solvers of the kind we
are building. Their point for our purposes is that a sparse factorization has *less* freedom to
recover than a dense one, because, as 7.4 is about to say, it cannot even search the whole trailing
submatrix, so a strategy that relies on a good candidate always being somewhere below is relying on
a set it is not allowed to look at.

### 7.4 What a sparse factorization takes away

Everything in 7.2 assumes the whole trailing submatrix is available to search. In a supernodal or
multifrontal factorization it is not, and the reason is the one 6.1 gave from the other side.

When supernode `K` is reached, only its front is **fully summed**. The rows below the front,
`update(K)`, are still owed contributions from other branches of the forest, which is exactly what
the update matrices of 6.1 are carrying around. Their values are partial sums. Pivoting on one would
mean dividing by a number that has not finished being computed, and the result would be wrong rather
than merely inaccurate.

So the candidate set is the front, not the trailing submatrix. In our non-root kernel this is
visible as one
loop bound: the scan that chooses a `2 x 2` partner stops at the end of the front, while the scan
that measures how large the pivot has to be runs the full column height, update rows included.
**The test uses everything available, the choice uses only what is fully summed.** That asymmetry is
the whole of the sparse adaptation, and it is worth carrying as a principle rather than as a detail
of one loop.

Named at the level of Bunch-Kaufman's own machinery, this is the cleanest single reason pure BK
cannot run at a supernode with update rows. BK cases 3 and 4 both turn on `r`, the row carrying the
largest off-diagonal in column `k`: case 3 promotes `A(r,r)` to the pivot by swapping `r` up, case 4
pairs `k` with `r` as a `2 x 2`. Both *use `r` as a pivot column*, and BK is free to pick `r`
anywhere in the trailing submatrix. But at a supernode with update rows, `r` can land on an update
row, which is not fully summed, so using it as a pivot column computes with an unfinished value: the
result is wrong, not merely inaccurate. The supernode is therefore not entitled to `r` even though
BK's search would select it. This is why the partner scan must stop at the front, and why what runs
is not BK but BK *confined to the front*. A confined search can come up empty, which the full search
never does. (It is also why the roots are the exception of 7.7: with no update rows, `r` is always a
front column, "confined to the front" and "the whole block" coincide, and full BK becomes available
again.)

It also breaks Bunch-Kaufman's case ladder. Case 4 in 7.2 is unconditional: once cases 1 through 3
have failed the `2 x 2` is taken, no test is applied to it, and the choice of `alpha` is precisely
what makes taking it provably safe. The ladder terminates because there is always somewhere to go.
Restrict the partner to the front and case 4 can fail, since the only acceptable partner may be an
update row, and the algorithm then runs out of cases with no pivot chosen. A strategy built on a
fixed constant and an unconditional fallback has nothing to offer at that point.

Two things are needed instead: a test that is allowed to refuse a `2 x 2`, and somewhere for a
refused column to go.

### 7.5 Threshold pivoting

The test came from Duff and Reid (1983), in the same paper as the multifrontal method, and it is
what MA27, its successor MA57, and MUMPS all use. Instead of a fixed constant chosen to make a
fallback safe, it takes a parameter `u`, the **threshold**, and asks of each candidate whether
accepting it would keep the factor bounded.

A `1 x 1` pivot `a_kk` is accepted when

```
|a_kk| >= u * max_{i != k} |a_ik|
```

and a `2 x 2` block `E` on rows and columns `k` and `r` is accepted when

```
| E^-1 | * [ max_{i != k,r} |a_ik| ]   <=   [ 1/u ]
           [ max_{i != k,r} |a_ir| ]        [ 1/u ]
```

componentwise, where `|E^-1|` is the matrix of magnitudes. Writing `E` and `E^-1` as in 7.2 and
calling the two maxima `max1` and `max2`, that pair of conditions expands into a pair of scalar
inequalities on `|det|`, and taking whichever binds gives the form our kernel evaluates directly:

```
maxmax = max( |c| * max1 + |b| * max2 ,  |a| * max2 + |b| * max1 )
accept when |det| >= u * maxmax
```

**What this buys is exactly what Bunch-Kaufman does not give.** The two conditions bound the entries
of `L` by `1/u`, by construction, because they are nothing other than the statement that the
multipliers `A21 E^-1` are bounded. Bunch-Kaufman bounds growth in the trailing submatrix and leaves
`L` free; threshold pivoting bounds `L` and lets growth follow. Each bounds the quantity its setting
cares about, and neither is a weaker version of the other.

The parameter is a genuine knob rather than a constant. Raising `u` accepts fewer pivots, which
costs work and fill and buys accuracy; lowering it does the reverse. At `u = 0` everything is
accepted and there is no pivoting left, and `0.5` is the practical ceiling, above which so little is
acceptable that the factorization struggles to progress. MA57 and MUMPS both default to `0.01`.

**One vocabulary trap, and it is worth naming, because it is the source of the confusion this
section exists to clear up.** The Duff literature calls its threshold `alpha` about as often as `u`,
and writes the `1 x 1` condition as `|a_kk| >= alpha * max |a_ik|`. Bunch-Kaufman's case 1 is
`|A(k,k)| >= alpha * omega1`. The same letter, in the same position, in the same inequality, naming
two entirely different objects: a tunable parameter around `0.01` in one, the fixed constant
`0.6403882...` in the other. Read either with the other in mind and threshold pivoting looks like
Bunch-Kaufman with most of its cases missing, which is how our own comments came to say that it was.

### 7.6 Delayed columns

A test that can refuse needs somewhere for the refusal to lead. It leads up the forest.

When no acceptable pivot is found among `K`'s remaining front columns, those columns are
**delayed**: `K`'s front contracts by that count, the columns are handed to `parent(K)`, and they
are tried again there. Containment is what permits this, the same property 6.1 needed: a delayed
column's rows lie inside the parent's index set, so the parent is able to hold it. And the parent is
a better place to ask, because more of the matrix has been summed there and the fully summed block
is larger, so the candidate set the column now faces is bigger than the one that rejected it. A
column may be delayed repeatedly, riding up one level at a time in the way a contribution does.

At a root there is no next hop. A column that cannot be pivoted there is a real failure, and it is
how the factorization detects that `A` is numerically singular.

Three consequences run through the rest of the implementation.

**The front grows.** A parent absorbs its children's delayed columns, so its front is larger than
the symbolic factorization predicted. This is why dynamic pivoting requires dynamic storage, one
vector per supernode, where a static factorization is content with flat buffers whose sizes were
known in advance. The two senses of *dynamic* meet here, and this is the connection between them.
4.5 sets the same widening beside amalgamation, which reaches it by a different road and pays for
it in the same explicitly stored zeros.

**Fill exceeds the prediction.** The structure of `L` is no longer settled before the values are
seen, so `nnz(L)`, the delay count and the counts of `1 x 1` and `2 x 2` pivots are properties of
the numbers rather than of the pattern. Nothing about them can be asserted from the symbolic phase.

**The update size is untouched.** Delaying moves a column between the front and the delayed set and
never changes `update(K)`, which is fixed by the symbolic factorization. This is the invariant that
the height `frontSize + delaySize + updateSize` is preserved, and it is what let 6.5 say that the
multifrontal stack peak is identical for a dynamically pivoted factorization and a static one on the
same forest.

The threshold governs how much of all this happens, which is the same knob of 7.5 seen from the
other end: a stricter threshold delays more columns, which grows fronts, which adds fill and work.

### 7.7 What Oblio does

**Rewritten 2026-07-26.** What follows describes the code as it now stands. Where later parts of
Section 7 still discuss "pass 1" and "pass 2" of a single `factorDynamicSupernode`, they are
describing the state before that date; the changes are listed at the end of this subsection.

Two kernels, `factorDynamicRootSupernode` and `factorDynamicNonRootSupernode`, chosen by the caller
on `parent[jj] == NIL`. That predicate is the same fact as "the supernode has no update rows", read
off the forest rather than off the storage, and it is the same split 0.9 made inside one function.
What it distinguishes is whether there is anywhere to delay a column to at all, and the two kernels
are now different algorithms rather than two settings of one.

**Non-root, threshold pivoting.** An ancestor exists, so nothing is forced. It holds a list of
candidate front columns and rotates through it: a column that fails is moved to the back and the
next is tried, with no elimination in between, so a later candidate sees exactly the values an
earlier one saw, and a round in which nothing is accepted ends the sweep and delays the remainder.
The `2 x 2` partner is chosen by a scan bounded at the front, per 7.4, and the block is accepted on
the threshold test of 7.5, which is AGL Figure 3.3's fourth branch and nothing else. The loop shape
is AGL Figure 3.4.

**Root, bounded Bunch-Kaufman.** `update(K)` is empty, so the front is dense and no ancestor is
waiting. Rather than force an acceptance, this kernel runs the chase of AGL Figure 2.4: take the
next unfactored column, follow its largest off-diagonal to the column holding it, and repeat until
either a diagonal passes its own `1 x 1` test or the maximum becomes mutual, at which point the
`2 x 2` is taken. There is no queue and no failing branch. The threshold here is the dense
`alpha = (1 + sqrt(17)) / 8`, not the tunable `u`, since a root front is dense and stays dense under
any symmetric permutation, so a stricter threshold costs comparisons and never fill.

```
                    root                          non-root
algorithm           AGL Figure 2.4                AGL Figures 3.4 and 3.3
forced 1x1          no                            no
2x2 partner         the chase's mutual maximum    front columns only
2x2 acceptance      max1 == max2, both diagonals  |det| >= u * maxmax
                      already failed their tests
on failure          cannot arise                  delay to the parent
threshold           (1 + sqrt(17)) / 8, fixed     u, tunable
```

**What changed on 2026-07-26, against what the rest of Section 7 describes.** The root kernel used
to be a weakened copy of the non-root one: it forced the last remaining candidate as a `1 x 1`
whatever it looked like, and accepted a `2 x 2` on `max1 == max2` alone, without reading the block's
determinant. Both are gone. The non-root kernel used to accept a `2 x 2` on that same
symmetric-maximum condition as an extra disjunct beside the threshold test; that clause is gone
too, having admitted blocks whose factor entries are unbounded (see the end of 7.8).

**This is not Bunch-Kaufman, and our own comments have said for some time that it is.** The `1 x 1`
condition is Duff's rather than Bunch-Kaufman's case 1, differing in that `u` is tunable where
`alpha` is fixed. The `2 x 2` condition is Duff's determinant test, and Bunch-Kaufman has no `2 x 2`
test at all, its case 4 being unconditional by construction. Cases 2 and 3 of 7.2 have no
counterpart here. The two strategies are certainly related, both choosing between `1 x 1` and
`2 x 2` blocks by comparing a diagonal against off-diagonal maxima, and they descend from the same
pair of papers, but they are different algorithms with different guarantees, and the difference is
the one 7.5 named: Bunch-Kaufman bounds growth, this bounds `L`.

The one place a Bunch-Kaufman comparison is genuinely instructive is case 3, which tests the
*partner's* diagonal, `|A(r,r)| >= alpha * omegar`, and promotes `r` by an immediate swap. We have
no such test. But the rotation reaches the same candidates anyway: `r` is a front column, so it
comes up as its own candidate later in the same round and is tested by the `1 x 1` condition, which
measures the same quantity, over values nothing has changed in the meantime. What differs is the
order in which pivots are accepted, not which columns are eligible to be.

**The threshold defaults to `0.1`**, in `NumFactorEngine` and in `DirectSolver` alike, an order of
magnitude stricter than the `0.01` of MA57 and MUMPS. Oblio therefore accepts fewer pivots and
delays more columns than those codes would on the same matrix, buying accuracy with fill and work.
That is a defensible setting, and it is the value the reference implementation shipped, but nothing
in the test suite exercises any other value, so whether the extra strictness earns its cost on our
matrices is unmeasured. A sweep reporting delays, `2 x 2` counts, `nnz(L)` and residual against `u`
would answer it, and it needs no change to the kernel.

### 7.8 Zero pivots, and an open question about singularity

A `1 x 1` pivot can be exactly zero. This section is about when that happens, what it means, and
the one place where what it means is not yet settled.

**In dense Bunch-Kaufman the story is complete.** Read the four cases of 7.2 against the possibility
that the accepted pivot is `0`. Case 1 accepts on `|A(k,k)| >= alpha * omega1` past the
`omega1 == 0` guard, so the bar is strictly positive and the pivot is nonzero. Case 2 needs
`|A(k,k)| * omegar >= alpha * omega1^2`, a strictly positive right side, which a zero `A(k,k)`
cannot meet. Case 3 accepts `A(r,r)` against `alpha * omegar` with `omegar >= omega1 > 0`, again a
positive bar. Case 4 forms a `2 x 2` whose determinant is `A(k,k) A(r,r) - omega1^2`, which is
`-omega1^2 < 0` even when both diagonals vanish, so the block is nonsingular, not a zero pivot.
Every
proper case guards against a zero by comparing a magnitude against a strictly positive threshold.

The only way through is the short-circuit: `omega1 == 0`, the whole subdiagonal of column `k` zero,
which by symmetry means row `k` is zero too, and if `A(k,k)` is also zero the entire row and column
are zero. That is the door, and it leads to a clean theorem. A zero subdiagonal and a zero diagonal
make a zero column; a matrix with a zero column is singular; so **a zero `1 x 1` pivot in
Bunch-Kaufman implies the matrix is singular.** The contrapositive is the useful direction: a
factorization that completes with every pivot nonzero certifies nonsingularity, which is why
LAPACK's
`sytrf` reports the first zero pivot in `INFO` rather than treating it as a numerical nuisance. And
it propagates, by the dynamic reading of Sylvester's law from 7.1: eliminating a nonsingular pivot
block is a congruence, so a nonsingular `A` has nonsingular trailing submatrices throughout, no
trailing row ever becomes zero, and the door never opens at any step.

The converse holds more tightly than one might expect, and getting it right matters for what
follows. One might guess that singular does not force a zero `1 x 1`, that a singular `A` could
instead surface as a `2 x 2` with `det = A(k,k) A(r,r) - omega1^2` happening to vanish. In dense
Bunch-Kaufman that cannot happen: **a case-4 `2 x 2` is always nonsingular, in fact strictly
indefinite.** Case 4 is reached only after case 1 fails (`|A(k,k)| < alpha * omega1`) and case 3
fails (`|A(r,r)| < alpha * omegar`), so both diagonals are small against the off-diagonal, and the
standard analysis sharpens this to `|A(k,k) A(r,r)| < omega1^2`, which makes
`det = A(k,k) A(r,r) - omega1^2` strictly negative. A negative determinant is one positive and one
negative eigenvalue, so the block is nonsingular by construction, and the constant `alpha` is
exactly what buys that: BK's `2 x 2` needs no determinant test because case 4 cannot form a singular
one. So in Bunch-Kaufman the rank deficiency of a singular `A` has nowhere to hide but a zero
`1 x 1`, reached through the `omega1 == 0` short-circuit and only there. The theorem is therefore an
iff after all: **a zero `1 x 1` pivot in Bunch-Kaufman if and only if `A` is singular**, with the
`2 x 2` structurally excluded as a carrier. (This is the dense luxury 7.5 contrasts against: a
threshold code has to *test* for `2 x 2` safety, where BK gets it free from `alpha`.)

**In Oblio the same question is open on both pivot sizes, because neither forcing path has the
shape Bunch-Kaufman's argument needs.** Take the `1 x 1` first. Pass 2's acceptance tests are
magnitude-against-threshold, `|d11| >= u * max1` and `|det| >= u * maxmax`, so by the Bunch-Kaufman
argument above neither can accept a zero: whatever the non-root kernel admits is nonzero. That half
transfers
cleanly. What does not transfer is the forcing path. Bunch-Kaufman's zero door is an empty column,
`omega1 == 0`, and the singularity theorem rests on that predicate: empty column plus zero diagonal
equals zero row of `A`. Pass 1's forcing path (7.7) is different. It forces on **threshold
exhaustion**: the candidate rotation cycles the whole front, nothing clears `u`, and the last
remaining column is accepted whatever its diagonal, because a root has no update rows and so nowhere
to delay. "No candidate cleared the threshold" is not the statement "a row is identically zero", and
the linear-algebra step that turned Bunch-Kaufman's door into a singularity certificate does not
obviously apply to this one.

The `2 x 2` is where Oblio departs from Bunch-Kaufman, and not in its favor. A zero column never
reaches a `2 x 2`: the `max1 == 0` guard in both passes routes an empty column to a `1 x 1` before
any partner is considered, so a formed block always has nonzero off-diagonals and its `a11` (the
larger first-column entry, after the row swap) is nonzero. So far so safe. But the root kernel used
to accept a
`2 x 2` on `max1 == max2`, the off-diagonal magnitudes alone, **without reading the block's
determinant** (7.7). This was precisely the guarantee BK case 4 has and that root acceptance
lacked: BK tests the
*diagonals* against the off-diagonal, which forces `det < 0`, while `max1 == max2` compares the two
columns' off-diagonal maxima and says nothing about the diagonals, hence nothing about the
determinant. So a singular block with equal off-diagonals, `[[d, d], [d, d]]` with `d` nonzero,
passed that test, was formed, and reached the diagonal solve where `u22 = det / a11 = 0` and the
block solve divided by zero. What BK structurally excludes, that acceptance could construct. The
non-root kernel does not
have the problem, because its `2 x 2` acceptance is the determinant test `|det| >= u * maxmax`,
which bounds `det` away from zero exactly as BK's `alpha` does, only by an explicit test rather than
a fixed constant. So the three routes to `2 x 2` safety line up: BK case 4 forces it with a
constant, Oblio's non-root kernel forces it with a threshold, and Oblio's old root pass did
neither. That root pass was the one
that runs at a dense root, which is exactly where singularity concentrates.

So the honest state is a conjecture with one part proved and two parts open, one per pivot size:

- **Proved.** Pass 2 cannot accept a zero pivot of either size, by the same
  threshold-against-positive-bar argument as Bunch-Kaufman cases 1 through 4.
- **Conjectured, `1 x 1`.** A forced pass-1 `1 x 1` acceptance carries an exactly zero diagonal only
  when the fully summed root block is singular, hence only when `A` is singular. Plausible, because
  a root block that rejects every candidate ordering at every threshold comparison behaves like a
  rank-deficient block, but the implication has not been shown, and "behaving like" is precisely the
  phrase the port discipline says to convert into a proof or a measurement before it is asserted.
  One piece of evidence stronger than "plausible", though short of a proof: 0.9 already *counts*
  this event as a rank loss. Its `factorDynamicLDL_` keeps a `rank_` initialized to `n` and
  decrements it at exactly this site, a forced or isolated `1 x 1` accepted with a zero diagonal, so
  the original
  author treated "forced zero `1 x 1`" and "rank dropped by one" as the same thing. That is a design
  belief embedded in code, not a theorem, but it is the same claim this bullet conjectures, held by
  the person who wrote the algorithm. (The port had dropped this counter and now restores it; see
  the ledger.)
- **Closed, `2 x 2`, on 2026-07-26.** This asked whether a root `2 x 2` accepted on `max1 == max2`
  could have a zero determinant, the `[[d, d], [d, d]]` shape, and divide by zero in the solve. The
  question no longer arises. The root kernel is now the chase of AGL Figure 2.4, which reaches its
  `2 x 2` only after both diagonals have failed their own `1 x 1` tests, so with the shared
  off-diagonal magnitude `m` it holds that `|det| > m^2 (1 - alpha^2)`, bounded away from zero by
  construction. The `[[d, d], [d, d]]` block fails the first diagonal test and is never offered.

The `1 x 1` conjecture above is likewise narrower than it was. Since 2026-07-26 there is no forced
acceptance at a root at all, so the only path to a zero `1 x 1` is the isolated-column branch, a
column whose off-diagonals are all zero, which is a zero row of the permuted matrix and therefore
singular outright. That is the same predicate dense Bunch-Kaufman forces on, so the gap the
conjecture was reaching across has closed and the reading is no longer conjectural for roots. It
remains open for a delayed column that arrives at a root and is accepted there by the ordinary
tests.

Two ways to settle what remains, cheapest first. **Measure it.** Build singular indefinite
matrices with a known null vector, factor them dynamically, and check whether an isolated-column
acceptance with a zero diagonal occurs, and occurs only there. Then the half that
actually matters for the claim: factor a batch of nonsingular indefinite matrices, including ones
constructed to delay heavily, and confirm that no zero pivot ever appears. The second batch is the
empirical stand-in for the theorem, and it is the natural shape of an `experiments/pivoting` study
or a targeted case in the numeric-factorization suite with matrices built to force delay to a root.
**Or prove it.** Show that at a root, with no update rows, a column whose off-diagonals all vanish
and whose diagonal vanishes occurs exactly when the root block is rank deficient. That is now a
small theorem about the isolated-column branch rather than about forced acceptance, the latter
having been removed.

One caveat runs under all of the above, and it is the same one that ends 7.5. This is exact
arithmetic. In floating point an exact zero almost never appears; a nearly singular `A` produces a
tiny pivot rather than a zero one, which divides to a huge multiplier rather than to `inf`. The
clean
iff is a statement about the exact algorithm, and the practical failure mode is smallness, a matter
of conditioning, not the exact zero this section is about. The zero pivot is the mathematically
sharp
case; the threshold of 7.5 exists precisely because the interesting real-arithmetic cases sit just
to
one side of it.

The consequence in the solve is recorded where it happens, in `diagonalDynamic`: a zero `1 x 1`
skips
its divide rather than dividing by zero, leaving that component of the solution to the value the
forward sweep produced. Whether that is a valid free choice (the system consistent, the null
component genuinely free) or a silently accepted contradiction (the system inconsistent) is not
distinguished there, for the reason that the honest test is the post-solve residual and not a
per-pivot comparison against a value the sweep has already modified. That is 0.9's abandoned
`InconsistentSystem` check, and the reasoning for leaving it abandoned is the same reasoning given
here.


### 7.9 Complete and partial pivoting, and the axes the words name

"Partial" and "complete" are worth pinning down, because they mean one thing in the general LU
world and the symmetric-indefinite world borrows the words and shifts them, and because a third
property Oblio actually turns on, the threshold, is often smeared into the same phrase while being
a different axis entirely.

**The general (LU, nonsymmetric) meaning is about how much of the matrix is searched** for the
pivot at step `k`.

- *Partial pivoting* searches one column, the current column below the diagonal, takes the
  largest-magnitude entry, and swaps its row up. `O(n)` comparisons per step, `O(n^2)` over the
  run. It is the standard, what `getrf` does, and it permutes rows only: `P A = L U`.
- *Complete pivoting* searches the whole trailing submatrix, all `(n-k)^2` entries, takes the
  global largest, and swaps both its row and its column into place. `O(n^2)` per step, `O(n^3)`
  total. It permutes on both sides, `P A Q = L U`, gives a stronger growth bound, and is rarely
  worth the cost.
- *Rook pivoting* sits between, searching alternately along a row and a column until an entry is
  maximal in both. Usually cheap, `O(n^3)` in the worst case.

So the LU axis is search extent: one column, a row-and-column zigzag, or everything.

**The symmetric-indefinite meaning keeps the axis but renames its endpoints after people.**
Symmetry forces the two-sided permutation `P A P^H` to preserve it (7.1), so "search a column and
swap only its row" is not even available, since a one-sided swap breaks symmetry. The Bunch family
fills the two ends:

- *Bunch-Parlett* (7.2) is the complete method: it searches the whole trailing submatrix each step
  for its `1 x 1`-or-`2 x 2` choice, at `O(n^3)` search cost. It is the symmetric analogue of
  complete pivoting.
- *Bunch-Kaufman* (7.2) is the partial method: at most two columns examined per step, `O(n)`
  search, `O(n^2)` total. It is called partial by direct analogy to LU partial pivoting, bounded
  cheap search rather than exhaustive.
- *Bounded Bunch-Kaufman* is the symmetric rook, the in-between (7.3).

Same axis, search extent; the endpoints are now complete = Bunch-Parlett, partial = Bunch-Kaufman.

**The axis Oblio actually turns on is a different one, and this is the trap.** Oblio pivots by
*threshold* (Duff-Reid, 7.5), and threshold is orthogonal to partial-versus-complete. Partial and
complete answer *how far is the search*; the threshold answers *what must a candidate pass to be
accepted*. The two are independent: one can pair a threshold acceptance test with a partial search
(MA57, MUMPS, and Oblio) or in principle with a complete one. So "threshold-partial pivoting", the
phrase 0.9 uses, is two choices at once: a partial *search* (bounded, Bunch-Kaufman-shaped) and a
*threshold acceptance* (`|pivot| >= u * max`) in place of Bunch-Kaufman's fixed-`alpha` test.

Oblio narrows the search axis one notch further than "partial", because a sparse supernode cannot
run even a full partial search: it may look only at its fully summed columns (7.4), delaying the
rest to a parent (7.6). So the search is partial *and* front-restricted, which is why the case
ladder can run out of candidates and why delay exists at all.

Laid against each other, the axes are:

```
axis               Oblio's choice              the alternatives
search extent      partial, front-restricted   complete (Bunch-Parlett); rook (bounded BK)
acceptance test    threshold u                 fixed constant (Bunch-Kaufman alpha)
permutation        two-sided (P A P^H)         one-sided is unavailable under symmetry
on failure         delay to the parent         dense methods cannot delay
```

The one-line name is **threshold-partial pivoting**, but the phrase is worth unpacking exactly
because "threshold" and "partial" answer different questions, and conflating them hides where a
design choice actually lives. The full-Bunch-Kaufman-at-the-roots proposal (recorded in TODO) is a
change to the *acceptance* axis alone, swapping the threshold `u` for the fixed `alpha` at the
fronts where delay is impossible, while leaving the *search* axis partial. That proposal is only
stateable once the two axes are held apart: it is not "more complete pivoting", it is a different
acceptance test at the same search extent.

### 7.10 The LU pseudocode, and why complete pivoting is almost never used

7.9 named the search-extent axis in the abstract. This subsection makes it concrete in the general
(nonsymmetric) LU setting, where it originates, then asks the question that decides which end of
the axis anyone actually uses: is the quality difference worth the cost. The answer, partial almost
always, complete essentially never for solving, is the same judgment the symmetric world reaches in
7.2 when it standardizes on Bunch-Kaufman over Bunch-Parlett, and it is worth seeing why in the
simpler LU case first.

Both fragments below are dense LU, 1-indexed, `A` overwritten in place by `L` (unit, below the
diagonal) and `U` (on and above). `k` is the current step.

**Partial pivoting.** Search the current column below the diagonal, swap one row.

```
for k = 1 .. n-1:
    p = argmax over i in k..n of |A(i,k)|      # largest magnitude in column k, rows k..n

    if A(p,k) == 0:
        continue                               # column already zero: singular here, skip

    swap rows k and p of A                     # row interchange only
    record p in the pivot vector

    for i = k+1 .. n:                          # eliminate below the pivot
        A(i,k) = A(i,k) / A(k,k)               # multiplier, stored in L
        for j = k+1 .. n:
            A(i,j) = A(i,j) - A(i,k) * A(k,j)
```

Search cost `O(n-k)` at step `k`, so `O(n^2)` over the run. Produces `P A = L U`, one permutation,
applied to rows.

**Complete pivoting.** Search the entire trailing submatrix, swap one row and one column.

```
for k = 1 .. n-1:
    (p, q) = argmax over i,j in k..n of |A(i,j)|   # largest anywhere in the trailing block

    if A(p,q) == 0:
        break                                  # whole trailing block zero: rank deficient, done

    swap rows    k and p of A                  # row interchange
    swap columns k and q of A                  # column interchange
    record p in the row-pivot vector
    record q in the column-pivot vector

    for i = k+1 .. n:
        A(i,k) = A(i,k) / A(k,k)
        for j = k+1 .. n:
            A(i,j) = A(i,j) - A(i,k) * A(k,j)
```

Search cost `O((n-k)^2)` at step `k`, so `O(n^3)` over the run, the same order as the elimination
itself. Produces `P A Q = L U`, two permutations, one on each side.

**The elimination body is byte-identical between them.** The divide and the rank-1 update are the
same; only the search and the swap differ.

```
                partial                     complete
search          one column, k..n            whole block, (k..n) x (k..n)
                O(n-k)                       O((n-k)^2)
swap            rows only                    rows and columns
factorization   P A = L U                    P A Q = L U
growth bound    2^(n-1) worst case          polynomial in n (much tighter)
```

The singular handling differs in a small telling way. Partial pivoting meets a zero *column* and
can `continue`, that column is finished, move on. Complete pivoting meets a zero *entire trailing
block* and can `break`, everything left is zero and the factorization is done at rank `k-1`.
Complete pivoting therefore reveals rank as a side effect: the step at which the running maximum
reaches zero is the numerical rank. Hold onto that, it is the one job complete pivoting keeps.

**The quality gap in theory is enormous.** Element growth, how large trailing entries get relative
to the original matrix, is what bounds the backward error. Partial pivoting's worst case is
`2^(n-1)`, exponential. Complete pivoting's is Wilkinson's bound, a slowly growing sub-polynomial
function that is small for any realistic `n`. Read only the bounds and one would always pick
complete.

**But the theory almost never bites.** The `2^(n-1)` partial bound is achievable only by
specifically adversarial matrices (the classic is a `-1/+1` triangular-plus-spike construction); on
matrices that actually arise, growth is small, empirically `O(n^{2/3})` or better, often a small
constant. Partial pivoting essentially never fails in practice, its instability is real but, in the
standard phrasing, spectacularly rare. So for the overwhelming majority of real problems the
practical quality difference between partial and complete is nil: both are backward stable, and the
complete search is `O(n^3)` spent defending against a failure that does not occur.

**And the cost always bites.** Complete pivoting's whole-block search is the order of the whole
factorization, so it roughly doubles the flop count, and worse, it destroys the memory-access and
blocking structure that makes dense linear algebra fast. Partial pivoting searches one column,
which is cheap, cache-friendly, and foldable into blocked BLAS-3 panel factorizations (what
LAPACK's `getrf` does). Complete pivoting's search cannot be blocked the same way, so the
real-hardware slowdown is worse than the flop count alone suggests. This is why LAPACK ships no
complete-pivoting `LU` at all: partial pivoting is the default and effectively the only choice for
general dense solve.

**Where the complete-pivoting instinct survives is rank, not stability.** Bringing the largest
remaining entry to the pivot orders the pivots by decreasing magnitude, so the step where the pivot
falls to near-zero is the numerical rank, and partial pivoting cannot do this reliably. When the
question is rank or the null space, that behavior is wanted, but even then complete-pivoting *LU*
is rarely the tool. QR with column pivoting (`geqp3`) is the workhorse rank-revealing
factorization: it permutes columns only, cheaper than full complete pivoting, and orders the
diagonal by magnitude, so it reveals rank at close to ordinary-QR cost. The SVD is the definitive
answer, more expensive still. So complete-pivoting LU as a *solver* is a museum piece, its
stability edge defends against nothing that happens; complete pivoting as an *idea*, largest
element to the pivot to reveal rank, lives on in QR-with-column-pivoting and the SVD.

**The tie back to Oblio.** The symmetric-indefinite field made the identical call for the identical
reason. Bunch-Parlett (7.2) is the complete-pivoting analogue, the whole-block search;
Bunch-Kaufman (7.2) is the partial analogue; and the field standardized on Bunch-Kaufman because
the complete search costs `O(n^3)` and defends against growth that, with the `2 x 2` blocks
handling the genuinely dangerous cases, essentially never materializes. Oblio inherits that
judgment. And the one place the complete-pivoting instinct resurfaces, revealing rank, is exactly
the `rank_` counter and the singularity questions of 7.8: that is Oblio caring about rank at the
roots, the one spot where "is this actually deficient" is the real question rather than "is this
stable". The same partial-beats-complete-except-for-rank story plays out one level down, in the
symmetric setting, for the same reasons.

### 7.11 Bunch-Parlett, and the growth bounds of the symmetric methods

7.10 gave the LU pseudocode and its partial-versus-complete growth story. This subsection does the
same for the symmetric pair, Bunch-Kaufman (partial) against Bunch-Parlett (complete), with the one
extra concern that the symmetric setting has and the LU setting does not. Bunch-Kaufman's own
pseudocode is in 7.2; the missing piece is Bunch-Parlett, given first.

**Bunch-Parlett** is the complete method: it searches the whole trailing submatrix each step for its
`1 x 1`-or-`2 x 2` choice, the symmetric analogue of complete-pivoting LU. As in 7.2, `alpha =
(1 + sqrt(17)) / 8`. Symmetric, 1-indexed, factoring `P A P^H = L D L^H`, working on the trailing
block `A(k:n, k:n)`.

```
bunch_parlett_step(A, k):
    # mu0 = largest magnitude on the diagonal; mu1 = largest magnitude off it
    mu0 = max over i in k..n of |A(i,i)|,     at diagonal position p
    mu1 = max over i > j in k..n of |A(i,j)|, at position (r, s)   # r > s

    if mu0 == 0 and mu1 == 0:
        return done                       # trailing block is entirely zero: rank exhausted

    if mu0 >= alpha * mu1:
        # a diagonal entry dominates: take it as a 1x1 pivot
        swap rows/cols k and p            # symmetric swap, both sides
        return 1x1 at k

    else:
        # no diagonal entry is large enough: take the symmetric 2x2 spanning r and s
        swap rows/cols k   and s          # bring the pair to (k, k+1)
        swap rows/cols k+1 and r
        return 2x2 at k, k+1
```

The elimination that follows each choice is the same block algebra as 7.2, a `1 x 1` scalar step or
a `2 x 2` block step; only the *selection* above differs from Bunch-Kaufman. And the difference is
exactly the search extent of 7.9: Bunch-Kaufman finds its candidate by scanning at most two columns,
while Bunch-Parlett scans the whole trailing block for `mu0` and `mu1`. That whole-block scan is the
`O((n-k)^2)` per step, `O(n^3)` total, that makes it complete and makes it expensive.

**The growth bounds line up with the LU case, partial exponential and complete polynomial.**

```
                          growth, worst case            LU analogue
LU partial pivoting       2^(n-1)          exponential  --- BK
LU complete pivoting      polynomial in n               --- BP

Bunch-Kaufman  (partial)  ~ 2.57^(n-1)     exponential
Bunch-Parlett  (complete) polynomial in n
```

Bunch-Kaufman's `(1 + 1/alpha)^(n-1)` is about `2.57^(n-1)`, from 7.2. It is *worse* than LU
partial's `2^(n-1)`, and the reason is the `2 x 2` blocks: a `2 x 2` step can amplify slightly more
than a `1 x 1`, and `alpha` is chosen to equalize the two, landing the base at 2.57 rather than 2.
Bunch-Parlett's bound is Wilkinson-style, the same slowly growing polynomial character as
complete-pivoting LU. So the headline is identical to 7.10: partial (BK) is exponential worst-case,
complete (BP) is polynomial, and BP is the safe-but-expensive one that is essentially never used for
the same reason complete-pivoting LU is not.

**The twist unique to the symmetric case: growth is not the whole story.** In LU, bounding growth
bounds the backward error, because partial pivoting makes every multiplier in `L` at most 1 in
magnitude automatically (each multiplier is a column entry divided by the largest entry in that
column). Growth is the only failure mode. In the symmetric case growth and the multipliers come
apart, which is the 7.3 point: Bunch-Kaufman bounds *growth* but leaves the entries of `L`
*unbounded*, and a large `L` entry inflates the backward error `|L| |D| |L^H|` even when growth is
controlled. So Bunch-Kaufman has two separate worst-case concerns where LU partial has one, and the
growth bound alone does not capture its stability.

```
              growth, worst      L entries bounded?     verdict
LU partial    2^(n-1)            yes (automatic, <= 1)  good in practice
LU complete   polynomial         yes                    safe, expensive
Bunch-Kaufman ~ 2.57^(n-1)       NO                     good in practice, but L can blow up
Bunch-Parlett polynomial         yes                    safe, expensive
rook (bBK)    close to BK        yes                    the practical fix for BK's L
```

That extra column, whether the multipliers are bounded, is why **rook pivoting (bounded
Bunch-Kaufman) exists in the symmetric world as a distinct third option with no clean LU analogue.**
Rook keeps Bunch-Kaufman's cheap partial search but adds enough extra searching, alternating along a
row then a column until an entry is maximal in both, to bound `L`, plugging exactly the hole the
growth bound leaves open. In LU no "rook" was ever needed because partial pivoting bounds the
multipliers for free; in the symmetric case it is needed because Bunch-Kaufman does not.

**In practice the verdict is the same as LU partial.** Despite the uglier `2.57^(n-1)` bound and the
unbounded `L`, Bunch-Kaufman essentially never blows up on real matrices: the exponential growth
needs an adversarial construction, and the unbounded-`L` case, while more real than LU's growth
problem, is rare enough that Bunch-Kaufman with `2 x 2` blocks is the standard and Bunch-Parlett is
the museum piece. Threshold pivoting (7.5) is the further practical answer a sparse solver reaches
for: rather than Bunch-Kaufman's fixed `alpha`, it tunes `u` to bound `L` *directly*, which is what
a sparse solver actually wants, since a large `L` entry is a fill problem there and not only an
accuracy one.

This is what fixes the ladder of the root proposals (the two TODO entries). If a root front is
numerically nasty, what has gone wrong is almost certainly Bunch-Kaufman's *unbounded `L`*, not its
growth. Rook fixes exactly that at roughly `O(k^2)` search cost, while Bunch-Parlett spends `O(k^3)`
to bound a growth that was probably fine. The extra column of the table is precisely why rook, not
Bunch-Parlett, is the natural next rung above plain Bunch-Kaufman at a hard root.

### 7.12 The reference codes: MA27, MA57, MUMPS, and where Oblio sits

The pivoting of Section 7 is not Oblio's invention; it is the Duff-Reid family, and three codes in
that family are the ones the literature and this document keep naming. They are worth
distinguishing, because they are the same algorithm at three generations and on three classes of
hardware, and Oblio sits at a definite point among them.

All three are multifrontal symmetric-indefinite solvers doing `L D L^H` with `1 x 1` and `2 x 2`
threshold pivoting and delayed columns, exactly the algorithm the rest of Section 7 describes. What
separates them is generation and target machine, not method.

**MA27** (Duff and Reid, early 1980s, Fortran 77) is the original. It is the first multifrontal
solver, and it established the whole approach Oblio descends from: the analyze phase picks a pivot
order from the sparsity pattern alone, and factorization adjusts it numerically through threshold
pivoting and delayed columns. It used minimum degree ordering. Its phases are vectorized but not
parallel, and it predates the widespread use of the Level-2 and Level-3 BLAS, so its dense front
work is hand-written rather than handed to `GEMM`.

**MA57** (Duff, early 2000s, Fortran 77 with a Fortran 90 interface layer) is MA27 modernized for a
single core. Same algorithm, restructured so the dense front factorizations go through Level-3
BLAS, which is the single biggest speed difference from MA27. It moved to approximate minimum
degree (AMD), the ordering Oblio also vendors, and added a large feature set MA27 lacked: multiple
right-hand sides, iterative refinement and error analysis, matrix modification, partial solves,
restart, static pivoting, and the ability to identify and change pivots. It is still a serial code,
still actively maintained. Its threshold default is `u = 0.01`, the value 7.5 and 7.7 compare
Oblio's `0.1` against.

**MUMPS** (Amestoy et al., late 1990s onward, Fortran 90 with C, MPI for parallelism) is the same
method taken to distributed-memory machines. It is designed primarily for distributed computation
under MPI, with shared-memory parallelism through multithreaded BLAS and OpenMP and out-of-core
support, and it uses BLAS and ScaLAPACK for the dense kernels. Two things follow from the
distribution that matter for Section 7. First, MUMPS is more general than the other two: it handles
unsymmetric matrices as well, not only symmetric. Second, distribution is what forces its pivoting
to differ, because a pivot candidate may live on another process and so be unreachable by the same
threshold search a serial code runs; this is exactly why the Duff-Pralet work on pivoting for the
sequential *and parallel* case (in the references) exists, and it is the parallel analogue of the
front restriction of 7.4, the candidate set narrowed one step further, now by process boundary
rather than by fully-summed-ness.

```
code    era      language                parallelism        role
MA27    1980s    Fortran 77              vectorized, serial  the original multifrontal
MA57    2000s    Fortran 77 (+ F90 API)  serial, BLAS-3      MA27 modernized: speed + features
MUMPS   1990s+   Fortran 90, C           distributed (MPI)   the method for parallel machines
```

**Where Oblio sits.** Algorithmically Oblio is MA27/MA57: serial, symmetric, threshold-partial
pivoting with `1 x 1` and `2 x 2` blocks and delayed columns, AMD ordering, dense fronts handed to
the Level-3 BLAS. The pivot test it runs is the later componentwise member of the Duff-Reid fork
specifically (7.13), not the early normwise version and not either trap of 3.4, its search is the
MA27-refined standard ordering, and its blocks are `1 x 1` and `2 x 2` only, not AGL's larger blocks.
It is not MUMPS, it does not distribute, and it is symmetric-only. It is written
in C++17 rather than any Fortran, which buys the storage-agnostic templating and the read-public
write-private accessor discipline that the Fortran codes express differently or not at all. The
threshold default of `0.1` is stricter than MA57's `0.01` (7.7), which is a deliberate inheritance
from the reference implementation rather than a match to the HSL default. The forest-parallelism
work noted elsewhere is the first step in the direction MA57's shared-memory successors (HSL_MA77
out-of-core, HSL_MA97 shared-memory and bit-compatible) and MUMPS took, but it is a step not yet
taken. So Oblio is a modern-C++ serial member of the MA27/MA57 branch of the family, with the
parallel and distributed branches as future territory rather than current ground.

**What the resemblance claims, and what is actually checked.** The family placement is solid, and the
AGL reading made it specific: "Oblio implements Duff-Reid" is right, and it is the middle rung of the
Duff-Reid fork, the later componentwise test, with the MA27-refined search and delayed columns, not
Duff-Reid in general. But the "subject to verification" that phrase invites hides two different
verifications, and only one of them is done.

The read-level check, does the code realize the later Duff-Reid test, is essentially complete:
`factorDynamicNonRootSupernode`'s `2 x 2` test is term-for-term AGL Figure 3.3, which AGL themselves
identify as the later Duff-Reid test, and its search is Figure 3.4, which they call a refinement of
MA27's. That
is the match 7.13 establishes by reading. The run-level check, does Oblio's *output* agree with an
actual Duff-Reid code such as MA27 or MA57, is not done, and it is confounded, so it would not be a
clean check even if run: before any two outputs could agree, three things would have to be matched
first, the threshold (Oblio defaults to `0.1`, MA57 and MUMPS to `0.01`), the ordering (Oblio and MA57
use AMD, MA27 uses minimum degree, and a different order means different fronts and a different
everything downstream), and the root handling. Oblio *is* verified numerically, but against its own 0.9
oracle at a residual near `3e-16`, not against the HSL codes.

Three departures keep Oblio from being identical even to MA27 and MA57, and naming them keeps the
resemblance honest. The root handling was the one real algorithmic departure, until 2026-07-26: its
forced-last
`1 x 1` and its `max1 == max2` `2 x 2` accepted without a determinant test are weaker than the bounding
test Duff-Reid and AGL would apply, and they are the subject of 7.8 and 7.13's pass-1 block. The `0.1`
threshold default against `0.01` is the same algorithm at a different operating point, so it delays more
columns and fills more on the same matrix. And the Cramer-rule `L21` formation in `factor2x2` is a
place where Duff-Reid prescribe a stable block solve, so MA27 and MA57 presumably form the multipliers
accordingly, and Oblio's explicit inverse may be doing less carefully what they do carefully (7.13's
pitfall block).

The three reference codes do not sit at equal distance. For MA27 and MA57 the resemblance is tight: the
same core method, the same later Duff-Reid test, `1 x 1` and `2 x 2` blocks only, delayed columns, and
AMD with MA57. For MUMPS it is looser: MUMPS shares the method, but its pivoting genuinely diverges,
because distribution restricts the candidate set further still, a pivot may live on another process,
the parallel analogue of the front restriction of 7.4 noted above. So Oblio is the serial branch, and
grouping MUMPS with MA27 and MA57 overstates the resemblance. The accurate one-line summary is that
Oblio is a serial Duff-Reid-family solver of the MA27/MA57 class, running the later componentwise test
with the MA27-refined search and delayed columns, differing in the root handling, the `0.1` default,
and the Cramer `L21` formation; and that "almost the same" is a claim about the algorithm family,
established by reading and by the 0.9 oracle, not an output-equivalence claim against the HSL codes,
which is a separate and confounded experiment.

### 7.13 Ashcraft, Grimes, and Lewis, and the source of Oblio's pivoting

7.7 recorded what the dynamic kernels do and left their provenance open, placing it only in
"the Duff-Reid family" without a specific source. The source is the paper of Cleve Ashcraft, Roger
Grimes, and John Lewis (1998), and the match is close enough to settle the question outright: one of
the two tests our kernel runs is identical to theirs, term for term. This subsection records the
paper and lays Oblio beside it, so that the mindset the code was written in sits on the page next to
the code.

**What the paper is.** It has a dense half and a sparse half. The dense half presents two close
cousins of Bunch-Kaufman, bounded Bunch-Kaufman (the rook pivoting of 7.3) and a fast Bunch-Parlett,
both of which fix the unbounded-`L` defect of 7.3 by taking as the `2 x 2` off-diagonal an entry
that is a local maximum in its column. The sparse half is the one that reaches us. Its conclusion is
that Bunch-Kaufman cannot be rescued in the sparse case, because the local-maximum off-diagonal a
bounded `2 x 2` wants may lie in the part of the front not yet assembled and so unavailable as a
pivot (the `A22` of 7.4), and that bounding `L` instead leads to one particular version of the
Duff-Reid algorithm. They then extend Duff-Reid in two directions: a more effective search for pivot
blocks, and blocks larger than `2 x 2`. So the paper both recommends the Duff-Reid bounding-`L` test
and supplies its stability proof from the bounding-`L` side. This is why 7.5 crediting Duff and Reid
and this subsection crediting Ashcraft, Grimes, and Lewis are both correct: theirs is the paper where
the test is written in the form our kernel uses and shown to bound `L`, and they themselves note
(their section 3.4) that this same test is identical to the later Duff-Reid `|D^-1|` magnitude test.
The three coincide.

**Same test, opposite motivation, and why Liu is the foil.** The coincidence just noted is not an
accident of two groups landing on the same inequality; it is one test reached from two directions, and
the paper is explicit about which direction is theirs. Duff and Reid motivated the `2 x 2` test by
bounding growth in the reduced matrix, the trailing Schur complement as elimination proceeds: their
stated requirement is `|delta_ij| <= max|a_uv| / u`, a bound on the modification each elimination makes
to a reduced-matrix entry. AGL bound a different object, the entries of `L` itself, the multipliers
`A21 E^-1`, keeping `|L| <= 1/u`. The one `2 x 2` inequality secures both, which is why the two land on
the same test, but each side's target is the other's byproduct. AGL put it as means and end: for Duff
and Reid, bounding `L` was the *means* to the *end* of bounding reduced-matrix growth; for AGL, bounding
`L` is the end and the growth bound follows, because a bounded `L` together with a backward stable
`2 x 2` solve is exactly what proves backward stability of the whole factorization. Same inequality, the
goal and the byproduct swapped. Their own summary is that their
acceptance test is the same as the later Duff-Reid one but their motivation differs, and that the
bounding-`L` perspective "provides a proof of stability for the Duff-Reid algorithm." So the
contribution of Section 3 is not a new test. It is a new justification for an existing one, offered as
the stability proof Duff-Reid's growth-motivated derivation never made explicit, a theory result in an
algorithm's clothing.

The Liu thread is the other half, and it is the writeup of a regression the authors shipped. Their
Section 3 opening tells the story plainly: they had a working Duff-Reid multifrontal code, then, drawn
by Liu's better reported sparsity and by the popularity of dense Bunch-Kaufman, they replaced the
Duff-Reid pivot selection with Liu's sparse Bunch-Kaufman variant, and it was the failure of that
variant that sent them back to reconsider Bunch-Kaufman in general. So Section 3 is not a neutral
survey. Liu is the foil because Liu is the thing they tried in place of Duff-Reid and it broke, failing
to bound `L`; Duff-Reid is the thing they returned to, now with a proof in hand for why it was right.
The "different direction" is thus both genuine theory, the `L`-bound derivation, and the path a burned
implementer actually walked.

**The strategies the paper weighs, and what each bounds.** The argument runs across a whole family of
pivot rules, dense and sparse, and the single axis that sorts them is whether they bound `L`, since
that is the property Bunch-Kaufman lacks and the one AGL are chasing. Every rule below also bounds
growth in the reduced matrix; `L` is the discriminator, so it is the only bound column.

```
+----------------------------+---------+--------+------------------------------------------+--------------------+
| strategy                   | scope   | bounds | AGL's role / verdict                     | where              |
|                            |         | L?     |                                          |                    |
+============================+=========+========+==========================================+====================+
| Bunch-Kaufman              | dense   | no     | the standard; unbounded L can degrade    | 2.1, 2.2           |
|                            |         |        | accuracy, and LAPACK's blocked form can  |                    |
|                            |         |        | be unstable                              |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| Bunch-Parlett              | dense   | yes    | complete pivoting, the accuracy          | 2.3                |
|                            |         |        | standard; O(n^3) search, too costly      |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| bounded BK (rook)          | dense   | yes    | AGL's new dense algorithm 1; replaces    | 2.4                |
|                            |         |        | blocked LAPACK BK                        |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| fast Bunch-Parlett         | dense   | yes    | AGL's new dense algorithm 2; replaces    | 2.5                |
|                            |         |        | LINPACK-like unblocked                   |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| Liu sparse threshold BK    | sparse  | no     | the one AGL adopted and that failed;     | 3.1, 3.2, Fig 3.1  |
|                            |         |        | case 2's L entry can be gq/g1 times 1/u  |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| Duff-Reid early, normwise  | sparse  | yes    | bounds L in the infinity norm; AGL call  | 3.4 (DR 1983)      |
|                            |         |        | it unnecessarily severe                  |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| Duff-Reid later,           | sparse  | yes    | the recommended test; = AGL Fig 3.3, re- | 3.3, 3.4, Fig 3.3  |
| componentwise              |         |        | derived from the L side; Oblio's         |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| explicit-BK componentwise  | sparse  | no     | trap: one refinement step past DR-later, | 3.4                |
|                            |         |        | collapses to explicit BK                 |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| Liu two-parameter DR       | sparse  | weak   | trap: only a 1/(2u^2) bound on L; AGL    | 3.4 (Liu)          |
| variant                    |         |        | call it a mistake                        |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
| larger blocks, s <= 5      | sparse  | yes    | AGL's extension of the bounding-L test   | 3.6                |
|                            |         |        | to s-by-s blocks                         |                    |
+----------------------------+---------+--------+------------------------------------------+--------------------+
```

Two groupings are worth reading off the table. Among the sparse rules there is one Liu strategy that is
central, the sparse threshold Bunch-Kaufman of Figure 3.1, which does not bound `L` and is the one whose
failure drove the paper; a second Liu entry, the two-parameter Duff-Reid variant, appears only as a
named mistake. And there is a Duff-Reid fork of two, or three if the trap is counted: the early normwise
test, the later componentwise test, and, one refinement step away from the later one, the
explicit-Bunch-Kaufman componentwise test that looks almost identical but bounds only growth. The three
Duff-Reid-family tests differ by a norm or by a single step, and they give three different guarantees,
which is why AGL spend Section 3.4 pulling them apart rather than treating "the Duff-Reid test" as one
thing. The rule Oblio and AGL both use is the middle one, the later componentwise test, explicitly not
the early normwise one and explicitly not either trap. AGL list three pitfalls in 3.4: the two above are
selection traps, and the third is Cramer's rule for the `2 x 2` solve, a solve method rather than a
selection rule, which is the subject of the pitfall block below.

Two things the table leaves out on purpose. The search orderings, Figure 3.4's standard MA27 baseline
and Figure 3.6's exhaustive order, are not here because they bound nothing: they order the candidates,
they do not decide acceptance, and they sit on the driver axis of the split below, not the predicate
axis. And the dense rows are context, the ladder the sparse argument climbs down from: bounded
Bunch-Kaufman and fast Bunch-Parlett are the dense algorithms that bound `L` by taking a
mutual-local-max off-diagonal, the move 7.4 shows is unavailable in the sparse front, which is what
forces the sparse side onto the threshold test instead. The two rows worth marking for later are the
Duff-Reid early normwise and the explicit-BK componentwise: they are precisely the neighboring tests
the 7.14 catalog names as the designed A/B against Oblio's current test, and AGL have already told us
which way each should come out, the normwise stricter, the explicit-BK trap unbounded in `L`.

**The identity.** Their Figure 3.3, the explicit bounding sparse pivot strategy, accepts a `2 x 2`
block on columns `1` and `q` when

```
max{ |aqq| * g1 + |aq1| * gq ,  |a11| * gq + |aq1| * g1 }  <=  |a11 * aqq - aq1^2| / ahat
```

where `g1` and `gq` are the two columns' largest off-diagonal magnitudes, `aq1` is the coupling, and
`ahat` is the threshold. Their symbols map onto 7.5's block `E = [[a, b], [b, c]]` and onto the
names the kernel actually uses as follows:

```
AGL 1998         7.5 / kernel     meaning
ahat             u                the threshold
g1               max1             column 1's largest off-diagonal magnitude
gq               max2             column q's largest off-diagonal magnitude
aq1              b = m            the coupling, the front-restricted partner entry
a11              a = d11          column 1's diagonal, the current candidate
aqq              c = d22          column q's diagonal, the chosen partner
a11 aqq - aq1^2  det              the block determinant
```

Under that identification their inequality is

```
max( |c| * max1 + |b| * max2 ,  |a| * max2 + |b| * max1 )  <=  |det| / u
```

which is 7.5's `maxmax`, cleared of the division, character for character. It is not a restatement of
their test in our terms; it is their test. `readPivotBlock2x2` fills `a`, `b`, `c`, `det`, and the
kernel evaluates exactly this.

**Figure 3.3 is the test, Figure 3.4 is the driver.** The paper gives two figures that are easy to
run together, and they are not one figure at two levels of detail; they are two different components
on two different axes, and they compose. Figure 3.3 (in their section 3.3) is an acceptance predicate,
answering *when do we accept this column or pair*. Figure 3.4 (in their section 3.5, where it is the
standard-ordering baseline that the exhaustive Figure 3.6 then improves on) is a search driver,
answering *which column next, and what on failure*. This is the same split 7.14 calls the seam,
predicate against driver, seen here in the paper's own two figures.

Figure 3.3 takes one column, its `a11` and largest off-diagonal magnitude `g1`, together with its
largest partner in `A11`, `aqq` with magnitude `gq` and coupling `aq1`, and returns one of four
verdicts: a `1 x 1` on `a11`, a `1 x 1` on `aqq`, the `2 x 2` block, or none. It states the `2 x 2`
test concretely, as the L-bounding inequality above. It says nothing about search order: its last
line, "no pivot found, repeat search using next column," hands off to an outer loop it does not
describe.

Figure 3.4 is that outer loop made concrete: the queue, the rotation (pop the front, and on failure
add to the rear), and the two-level termination (inner, until a pivot is found or every column has
been tested since the last success; outer, until the queue empties or nothing is found). It calls an
acceptance test at each column but does not fix which one. Its `2 x 2` branch is written abstractly,
"`Dj,qj` is acceptable as a `2 x 2` pivot," an open slot, and in their section 3.5 text they say
plainly that the slot takes "either of the Duff-Reid tests or Liu's condition." It is deliberately
test-agnostic.

So there is one narrow sense in which 3.4 implements 3.3: the outer loop 3.3 only gestures at is
exactly what 3.4 supplies. But 3.4 does not reproduce 3.3's per-column test, and there are two real
differences.

The first is that 3.4's `2 x 2` test is a slot and 3.3 is one filling of it. Drop 3.3's inequality
into the slot and the two compose; but the slot could equally hold the Duff-Reid-early test or Liu's,
which is why the whole Tier-0 catalog of 7.14 fits behind this one branch.

The second is that 3.3 tests the partner's `1 x 1` inline and 3.4 does not. The per-column bodies,
side by side:

```
+--------------------------------------------+--------------------------------------------+
| Figure 3.3, the predicate (column +        | Figure 3.4, the driver's per-column body   |
| partner)                                   |                                            |
+============================================+============================================+
| g1 == 0  ->  nothing                       | gj == 0  ->  nothing                       |
+--------------------------------------------+--------------------------------------------+
| |a11| >= u*g1  ->  1x1 on a11 (current)    | |ajj| >= u*gj  ->  1x1 on ajj (current)    |
+--------------------------------------------+--------------------------------------------+
| |aqq| >= u*gq  ->  1x1 on aqq (partner)    | no partner-1x1 branch; aqq is reached only |
|                                            | when it rotates to the front later         |
+--------------------------------------------+--------------------------------------------+
| the L-bounding inequality  ->  2x2 block   | "Dj,qj is acceptable as a 2x2 pivot"  ->   |
|                                            | 2x2 block, with the test left unspecified  |
+--------------------------------------------+--------------------------------------------+
| else  ->  repeat with the next column      | else  ->  add j to the rear of the queue   |
|                                            | (delay)                                    |
+--------------------------------------------+--------------------------------------------+
```

Figure 3.3 checks three candidates for the column, the current `1 x 1`, the partner `1 x 1`, and the
`2 x 2`; Figure 3.4 checks two, the current `1 x 1` and the `2 x 2`, then delays. The partner's
`1 x 1` chance in 3.4 arrives later, when `qj` itself rotates to the front of the queue. The same
columns are eligible under both; what differs is the order and the timing in which pivots are
accepted, not which columns can be. Figure 3.4 also wraps this body in the queue and the two-level
loop, machinery Figure 3.3 has none of.
**Of the two figures, Oblio takes 3.4's shape and 3.3's test.** `factorDynamicNonRootSupernode` is
3.4's search machinery, the queue, the rotation, and the delay, with 3.3's concrete `2 x 2`
inequality
dropped into 3.4's abstract slot, and it follows 3.4's leaner per-column body: one current-column
`1 x 1`, with the partner `1 x 1` reached by rotation rather than tested inline, not 3.3's double
`1 x 1`. The symmetric-maximum fast-accept it adds on top is the bounded-Bunch-Kaufman idea of 2.4
and 2.5. This is why the comparison below cites Figure 3.3 for the `2 x 2`-accept row and Figure 3.4
for the queue, rotation, and delay rows: each row names whichever figure owns that piece. And it is
why the partner `1 x 1`, the one place 3.3 and 3.4 differ in eligibility timing, is a separate and
smaller knob that sits with the driver, not with the test.

**A shift, not a refinement.** It is worth being exact about the *kind* of gap between the two figures,
because the natural expectation is the other kind. A textbook move from a mathematical algorithm to an
implementation, an `A[i][j]` matvec down to the `colPtr` / `rowIdx` / `val` loop, is a refinement: it
is output-equivalent on every input, and there is a structure-preserving map between the two, one
operation to one operation, with only the addressing changed. No decision differs; `A[i][j]` and
`val[k]` name the same number, and the two versions make the same choices in the same order. The move
from 3.3 to 3.4 is not that. It fails both tests, and the failure is the point.

The cleanest witness is the partner-`1 x 1` difference already noted, followed to its consequence.
Take a focus column whose own `1 x 1` fails, and suppose both its partner's `1 x 1` and the pair's
`2 x 2` would pass. 3.3, in branch order, reaches the partner's `1 x 1` before the `2 x 2`, so it
eliminates the partner as a `1 x 1` and leaves the focus. 3.4 has no partner-`1 x 1` branch, so it
reaches the `2 x 2` and eliminates the pair. Same matrix, same values, a `1 x 1` against a `2 x 2`:
different pivots, and so a different `D`, a different reduced matrix, and different fill. Two
algorithms that select different pivots on one input are not a refinement pair.

Underneath that witness are three shifts, none of them an addressing change.

**The unit of work is redefined.** 3.3's step is "given a focus column and its best partner, choose
the best pivot among three candidates: focus-`1 x 1`, partner-`1 x 1`, the pair." 3.4's step is "visit
one column; accept it, or accept its pair, else defer it." These are different notions of what a step
is, not one step seen at two resolutions.

**A prioritized move is dropped.** 3.3 can promote the partner as a standalone `1 x 1`, and ranks that
above the `2 x 2`. 3.4 can promote no column other than the one it is visiting, except as the mate of
a `2 x 2`. This is exactly why one cannot plug 3.3's body into 3.4's loop unchanged: 3.4's bookkeeping
is `remove j`, and on a `2 x 2` also `remove qj`, and it has no case for "accept a third thing, the
partner alone, while visiting j." The shapes do not mate, and that missing case is algorithmic
content, not a data-structure detail.

**A capability is added.** 3.4 can delay: a column leaves the front, the front shrinks, and the column
is postponed to the parent (7.6). 3.3's only non-accept outcome is "try the next column." Delay lives
only at 3.4's level, and there is nothing in 3.3 for it to refine from.

The rotation is what papers over the dropped move, and it does so only data-dependently. Within one
round no elimination happens between visits, so if the partner is the next acceptable candidate it is
promoted on the same values and the two coincide; but if some other column accepts first, its
elimination has updated the partner's diagonal, and the partner's later `1 x 1` is now a test on
numbers 3.3 never saw, or that other column simply becomes the pivot instead. Data-dependent
equivalence is the opposite of a refinement, which is equivalent by construction.

So the refinement the matvec model expects does exist in this code, just between a different pair. It
is Oblio's kernel against 3.4 with 3.3's `2 x 2` test dropped into the slot: same loop, same choices,
same output, with only the addressing changed, `gblToLcl`, `val[at(i_, k_)]`, and the flat or dynamic
storage underneath (the seam of 7.14). That pair is the `A[i][j]`-to-CSR relationship, and
it holds cleanly. The gap we sensed is that `(3.3, 3.4)` is not that relationship; `(kernel, 3.4 with
3.3's test)` is.

**On the authors' intent.** The two figures sit in two different sections doing two different jobs, and
that placement is the best evidence of why the move is a shift. Figure 3.3 is in their section 3.3,
"Pivoting to bound L," where its job is to establish the acceptance test the L-bound justifies. Figure
3.4 is in their section 3.5, "An exhaustive pivoting strategy," where its job is to be the standard
search baseline, their MA27 refinement, that the exhaustive Figure 3.6 then improves on. So the move
from 3.3 to 3.4 is not a refinement within one algorithm; it is a change of subject, from the
acceptance criterion to the search ordering, and 3.4's abstract test slot ("either of the Duff-Reid
tests or Liu's condition") is the sign that the search was meant to be independent of the test. One
mechanical reading of the single real omission, the partner-`1 x 1`, is that 3.4's queue reaches it by
another route: every column is tested for its own `1 x 1` when it is visited as a focus, so the dropped
branch removes only the option to promote a partner *early*, which the authors evidently did not judge
worth a special case. Whether that is the whole of their intent, or whether the early promotion was
given up for a reason the paper does not state, is the thing still worth zooming in on.

**Oblio beside the paper.**

```
+------------------------------------+--------------------------------------------+----------------------------+
| Oblio, non-root kernel             | Ashcraft, Grimes, Lewis (1998)             | relationship               |
+====================================+============================================+============================+
| the 2x2 test |det| >= u * maxmax   | Figure 3.3, the explicit-bounding test     | term-for-term identical    |
+------------------------------------+--------------------------------------------+----------------------------+
| candidate list: pop the front,     | Figure 3.4, "standard ordering" (their     | same                       |
| test the current column's 1x1 then | refinement of MA27's search)               |                            |
| one 2x2 with its largest front     |                                            |                            |
| partner, else push to the rear     |                                            |                            |
+------------------------------------+--------------------------------------------+----------------------------+
| outer/inner loop: rotate until a   | Figure 3.4: "until pivot found or all      | same control structure     |
| pivot is found or every column has | columns tested since last successful       |                            |
| been tested since the last         | pivot", "until queue empty or no pivot"    |                            |
| success, then stop                 |                                            |                            |
+------------------------------------+--------------------------------------------+----------------------------+
| leftover columns delayed,          | postponing columns from A11 into the       | same, and for the same     |
| frontSize -= delaySize             | enlarged A22                               | reason (A22 unavailable,   |
|                                    |                                            | not fully assembled)       |
+------------------------------------+--------------------------------------------+----------------------------+
| partner scan confined to [j,       | the "potential second pivot column must be | same, their 7.4 in one     |
| frontSize)                         | taken from the first block column",        | sentence                   |
|                                    | entries in A22 unavailable                 |                            |
+------------------------------------+--------------------------------------------+----------------------------+
| mutual-max fast-accept m == max1   | bounded BK / fast BP: a 2x2 whose off-     | same idea                  |
| && m == max2 (both passes)         | diagonal is a local maximum in both        |                            |
|                                    | columns bounds L for free (2.4, 2.5, and   |                            |
|                                    | the opening of 3.3)                        |                            |
+------------------------------------+--------------------------------------------+----------------------------+
```

The `1 x 1` test is theirs too: our `|d11| >= u * max1` is their `|ajj| >= ahat * gj` of Figure 3.4,
so both of the two acceptance conditions the kernel evaluates are identical to theirs, not only the
harder one.

The candidate list of 7.7 is their queue, popped from the front with failures sent to the rear; the
barren-pass stop is their "all columns in the queue tested since the last successful pivot"; the
front-only partner scan of 7.4 is their rule that the second pivot column must come from the fully
assembled block; and the delay of 7.6 is their postponement of a column into an enlarged `A22`. The
symmetric-maximum disjunct 7.7 named, accepting a `2 x 2` on `m == max1 == max2` without the
determinant test, is the bounded-Bunch-Kaufman idea from their dense half carried into the sparse
one: when the coupling is the largest off-diagonal in both columns, the block bounds `L` for free and
needs no threshold.

**What Oblio did not take, and it is their two headline extensions.** Their more effective search is
Figure 3.6: a bias toward `2 x 2` pivots that reuses each column maximum and tests all pairwise
`2 x 2` blocks among the leading `m2x2` columns (they recommend `m2x2 = 5`) before accepting any
`1 x 1`. We use the simpler standard ordering of Figure 3.4 instead, one `2 x 2` per column against
its largest front partner. Their second extension is blocks larger than `2 x 2` (their section 3.6),
each solved through a QR factorization of the block; we do `1 x 1` and `2 x 2` only. So Oblio took
the core, the bounding-`L` test and the MA27-refined search, and neither extension. If a future pass
wants more `2 x 2` pivots for kernel speed, Figure 3.6 is the recipe, and it changes only the search
order, not the acceptance test.

**Where the old root pass departed from them.** Their framework assumes a refused column can be
postponed into
`A22`. At a root there is no `A22`, and they prescribe no forced acceptance there: an unfactorable
final front is singularity, the reading 7.8 reaches from the other side. Pass 1's `max1 == max2`-only
`2 x 2` and its forced last `1 x 1` are our own root handling, and they are strictly weaker than the
bounding test of Figure 3.3, which bounds `L` at a root as anywhere else. This is the same gap the
fixed-`alpha` root proposal and the 7.8 conjecture circle. One option those notes do not yet name is
to run Figure 3.3's test at the root as well, since it is valid there, and force only the genuinely
last column that no test can accept.

**The three pitfalls, and where Oblio stands.** AGL name exactly three pitfalls in their section 3.4,
and frame them with one theme: all three are traps that focusing only on bounding reduced-matrix growth
leads to, which is to say traps you avoid by bounding `L` instead (the means-and-end point above). Two
concern the acceptance test, one the block solve.

1. **The block solve: Cramer's rule or explicit inversion.** Solving the `2 x 2` systems by explicit
inverse rather than by a normwise backward stable scheme. This is a *solve-method* trap, not a test
one, and it is the one 7.15 works out in full. AGL's fix is a stable solve everywhere `E` is inverted,
their own choice being Gaussian elimination with complete pivoting.

2. **The over-refined test: an explicit Bunch-Kaufman in disguise.** Refining the acceptance test one
step past Duff-Reid to a componentwise magnitude bound yields a test that is nothing more than an
explicit version of the Bunch-Kaufman test: it bounds growth but gives no useful bound on `L`, so it
does not suffice for backward stability. A *selection-test* trap that looks like a sharpening and
silently drops the `L` bound. It is the explicit-BK componentwise row of the catalog above.

3. **The relaxed test: the Bunch-Parlett square.** Taking the `2 x 2` growth bound to be the square of
the `1 x 1` bound, a weaker condition that admits `L` entries as large as `1/(2u^2)`, around `1e6` for
the small `u` a sparse code uses. Liu's two-parameter Duff-Reid variant uses exactly this, and AGL call
it in hindsight a mistake. The other *selection-test* trap, this one by relaxing rather than
over-refining, the Liu two-parameter row of the catalog above.

**Where Oblio stands: clear of 2 and 3, in 1 only, and only on the factor side.** Oblio runs the later
Duff-Reid, AGL Figure 3.3 test, which is exactly the `L`-bounding test that pitfalls 2 and 3 trade
away, so it is clear of both by construction: it is not the explicit-Bunch-Kaufman bound of 2, and not
Liu's weaker bound of 3. It walks into pitfall 1 alone, and there only in the factor: `factor2x2`
forms the multipliers by Cramer, while the solve's `diagonalDynamic` already uses a pivoted LU. So the
planned repair, LU everywhere with partial or complete pivoting, kept as selectable variants, closes
the one pitfall Oblio is in and leaves the acceptance test where it already correctly sits.

**One pitfall we walk into.** They note (their section 3.3, restated in 3.4) that the test bounds
`|L|` but not the condition number of the `2 x 2` block, so the `2 x 2` systems must be solved by a
normwise backward stable scheme, their own code by Gaussian elimination with complete pivoting, and
they name Cramer's rule and explicit inversion as the trap.
`factor2x2` forms the multiplier columns as `(t1 * d22 - t2 * d21) / det`, which is exactly that
explicit inverse. The `diagonalDynamic` solve (DESIGN_DECISIONS, 2026-07-24) carries a `2 x 2` with
its own partial pivoting, so the solve side is already off this path; it is the factor's formation of
`L` that remains on it. The acceptance test has bounded `|L21|` by `1/u` before we arrive here, so it
is likely harmless in practice, but it is precisely the case they predict can lose accuracy, a
`2 x 2` block that is ill-conditioned yet passed the magnitude test, and it is worth a measurement
rather than an assumption.

**The block solve, where symmetry is deliberately dropped.** The pitfall above is the factor's half of
applying `D^-1`; the solve's half is its constructive counterpart, and it is where the symmetric
discipline of the whole factorization is set aside on purpose. Once dynamic pivoting has chosen a
`2 x 2` block, that block is a tiny dense object, and `diagonalDynamic` solves the `D z = y` step
against it by an explicit LU with partial pivoting, the `abs(a11) >= abs(a21)` row swap and the small
LU beneath it. It does not exploit the block's symmetry, and that is the point.

The non-symmetric solve is not a liberty taken because it is cheap; for the blocks this factorization
produces it is required. A diagonal entry becomes a `2 x 2` pivot precisely because it was too small or
too zero to be a stable `1 x 1`, so the `2 x 2` blocks dynamic pivoting emits are exactly the ones with
a problematic diagonal. The archetype is the exchange block `[[0, 1], [1, 0]]`: perfectly symmetric,
and a symmetric solve would pivot on its zero diagonal and divide by zero. Partial pivoting swaps the
nonzero off-diagonal into the pivot position and the divide is safe. Solving such a block
"symmetrically" would reintroduce the division-by-zero that going to a `2 x 2` was meant to avoid, so
the asymmetry is the safety, not a shortcut. (This is the zero-*diagonal* case, which partial pivoting
handles cleanly; it is a separate matter from the zero-*determinant* block of 7.8, which the old
root pass could
construct and which no local solve method can rescue.)

Why this is fully contained, and leaves the overall symmetry untouched, is worth stating exactly. The
outer factorization `P A P^H = L D L^H` is fixed before the solve runs. The block step computes
`z_k = D_k^-1 y_k`, and the mathematical value of `D_k^-1 y_k` does not depend on how `D_k` is
factored; only its accuracy does. So symmetry is a property of the outer factorization, which nothing
here touches, and the block solve is an inner dense problem where the method is a pure accuracy choice.
This is the two independent layers of pivoting the design notes name (DESIGN_DECISIONS, 2026-07-24):
the factor pivots symmetrically to *choose* the block, about growth and fill during elimination, and
the solve pivots non-symmetrically to *invert* the frozen block, about the stable inversion of four
fixed numbers, and neither layer knows the other exists. Both of the reasons symmetry is worth keeping
at scale are absent at the block level: there is no sparsity in a `2 x 2` to preserve, and the
arithmetic saving from symmetry on so small a block is nil, so the only property left to serve is
stability, which is what favors dropping symmetry.

This containment is also what makes larger pivots possible, which is where the `2 x 2` story
generalizes. AGL's `s <= 5` blocks (their 3.6) work because each block is a self-contained dense solve,
so the block may grow to any size the search can afford and the outer sparse-symmetric structure stays
indifferent to it. They factor those larger blocks with a QR of the block and solve the `2 x 2` by
Gaussian elimination with complete pivoting; both are non-symmetric stable dense factorizations, the
same move. And it is their stated design principle: the acceptance test bounds `|L|` but not the
condition number of the block (their 3.3), so the block must be solved by a normwise backward stable
scheme, never Cramer's rule. Oblio's partial-pivoted LU is the `2 x 2` instance of that prescription.
Oblio is `1 x 1` and `2 x 2` only today; `k x k` is AGL's extension and a future direction, but the
containment described here is exactly why it would be a clean one.

The two halves of `D^-1` are worth holding together. The same block inverse is applied in two phases
and two ways: the solve applies it stably, by the partial-pivoted LU above, and the factor applies it
by explicit inverse, `L21 = A21 D^-1` formed as `(t1*d22 - t2*d21)/det`, the Cramer's-rule pitfall of
the block just above. Both are the same four numbers inverted; only the solve half currently inverts
them stably. So "symmetry deliberately dropped, and done stably" describes the solve; the factor's
formation of `L21` drops symmetry too, but by the method AGL warn against.

**Provenance.** Unlike 7.5 through 7.7, which were written from the code with the source uncertain,
this subsection is written from the paper itself, read in full. It closes the question 7.7 left open,
and 7.7's "this is not Bunch-Kaufman" stands unchanged: what is added here is the positive
identification of what the pivoting *is*, not a correction of what it is not.

### 7.14 Multiple pivoting strategies, and the seam that makes them cheap

Pivoting is the crux of dynamic LDL, and it is not one strategy but a family: 7.2, 7.3, 7.5, and 7.13
each named a different rule for the same decision, and we want to be able to experiment with all of
them. Oblio can hold more than one acceptance test, each a different strategy, but the ease
of that is real only for part of the family, and it is worth setting down which part, because that is
what tells us where the cheap experiments are and where the expensive ones hide.

**The seam.** The whole strategy lives at one point: given the current front state and the cursor `j`,
decide the next accepted pivot, a column for a `1 x 1` or a pair for a `2 x 2`, or report none-found.
Everything downstream of that decision is shared and does not vary between strategies: `factor1x1`
and `factor2x2`, the swap, the rank-1 and rank-2 trailing updates, `readPivotBlock2x2`, the
candidate-list rotation, the delay tail, and the height invariant `frontSize + delaySize + updateSize`
of 7.6. A strategy is therefore a decision policy bolted into an otherwise fixed driver. That is what
makes "more than one kernel" true, but only for the strategies that vary the predicate. Some vary the
control flow instead, and those cost more. Four tiers, cheapest first.

**Tier 0, swap the predicate, same driver.** Only the accept tests change, the `|d11| >= u * max1`
and the `|det| >= u * maxmax` at the accept site. Every threshold-family variant is here, and these
are a few lines each. The catalog below is the Tier-0 menu.

**Tier 1, same elimination, different search order.** The candidate enumeration changes but the accept
test and the elimination do not. AGL's exhaustive `2 x 2`-biased order (Figure 3.6, and 7.13) is the
case: a window of the leading `m2x2` columns, each column maximum reused, all pairwise `2 x 2` blocks
tried before any `1 x 1`.

**Tier 2, different block size.** A new elimination path and a new numeric solve. AGL's larger blocks
(their 3.6), `s <= 5`, a principal-minor descent to find an acceptable block, a QR of the block to
apply it.

**Tier 3, different shape entirely.** A posteriori pivoting: factor optimistically on the symbolic
pivots, check the computed factor afterward, delay only what failed. This is Duff, Hogg and Lopez
(2020, in the references), the one strategy we name that we have not yet opened.

The Tier-0 menu, which is where the variety is and where each entry guarantees something different:

```
+----------------------------+------------------------+----------------------+------------------------------------+
| acceptance test            | what it guarantees     | source               | Oblio relation / note              |
+============================+========================+======================+====================================+
| bounding-L threshold,      | bounds |L| by 1/u      | Duff-Reid late, AGL  | current Oblio (7.5, 7.7)           |
| |det| >= u * maxmax        |                        | Fig 3.3              |                                    |
+----------------------------+------------------------+----------------------+------------------------------------+
| Duff-Reid early, normwise  | bounds L, infinity     | AGL section 3.4      | AGL call it unnecessarily severe:  |
|                            | norm                   |                      | rejects blocks the magnitude test  |
|                            |                        |                      | accepts; a clean A/B against the   |
|                            |                        |                      | current test                       |
+----------------------------+------------------------+----------------------+------------------------------------+
| Liu sparse threshold       | does not bound L       | AGL Fig 3.1          | the ahat/tau four-case sparse BK   |
| Bunch-Kaufman              |                        |                      | our comments wrongly claimed we    |
|                            |                        |                      | ran; extra case 2, tau tied to     |
|                            |                        |                      | ahat by a cubic                    |
+----------------------------+------------------------+----------------------+------------------------------------+
| explicit-BK componentwise  | bounds growth, not L   | AGL 3.4, pitfall #2  | a named trap; worth having to show |
|                            |                        |                      | the trap on our own matrices       |
+----------------------------+------------------------+----------------------+------------------------------------+
| bounded BK / rook          | bounds L via a mutual  | AGL 2.4 (7.3)        | our m == max1 == max2 disjunct is  |
|                            | local-max off-diagonal |                      | this, restricted to the front      |
+----------------------------+------------------------+----------------------+------------------------------------+
| fast Bunch-Parlett         | bounds L via a mutual  | AGL 2.5              | dense cousin of rook; same front   |
|                            | local-max off-diagonal |                      | difficulty in the sparse case      |
+----------------------------+------------------------+----------------------+------------------------------------+
| Bunch-Parlett complete     | the strongest bound    | Bunch-Parlett 1971   | O(n^3) whole-submatrix search; an  |
|                            |                        | (7.2)                | accuracy oracle, not a production  |
|                            |                        |                      | path                               |
+----------------------------+------------------------+----------------------+------------------------------------+
```

Three axes are orthogonal to the choice of test and combine with any of them. **Root handling:** the
current forced-last plus `max1 == max2` (7.7), the Figure 3.3 bounding test applied at the root as
7.13 proposes, or a fixed-`alpha` Bunch-Kaufman at the root as the fixed-`alpha` proposal has it; this
is the 7.8 territory. **The threshold value** `u` itself, the sweep 7.5 and 7.7 both call for,
reporting delays, `2 x 2` counts, `nnz(L)`, and residual against `u`. **The `2 x 2` numeric solve,**
Cramer's rule in `factor2x2` against AGL's Gaussian elimination with complete pivoting against the
partial-pivoted solve already in `diagonalDynamic` (7.13); this last is a correctness axis rather than
a strategy, but it rides along with any `2 x 2` policy.

Two cautions, so the seam is not oversold. The two-pass split on `updateSize` (7.7) is itself part of
some strategies and absent from others, so it may have to move behind the seam rather than stay in the
driver. And the accept tests currently read fields the driver computed, `max1`, `max2`, `m`, and the
block from `readPivotBlock2x2`, so the seam is really "the driver hands the policy a small view of the
front, and the policy returns a decision"; that view is worth defining once, so that every Tier-0
variant shares it.

**The cheapest first move that teaches the most** is a Tier-0 pair: keep the current bounding-L test
as the baseline and stand up the Duff-Reid-early normwise test beside it. Both bound L and differ only
in the norm, so any divergence in delays, `2 x 2` counts, `nnz(L)`, and residual is a pure measurement
of componentwise against normwise with nothing else moving, and it forces the seam to be cut cleanly
on an easy case before the harder tiers lean on it. The method is the one the ordering family used:
prototype twins in an `experiments/pivoting` folder against the compiled 0.9 oracle, with whole-output
comparison, before anything enters the main kernel.

**This subsection is a plan, not a record of built code.** Nothing here is implemented yet; it is the
map of what could be, written while the AGL reading of 7.13 is fresh, so the experiments can be picked
up in any order later.

### 7.15 Cramer's rule for the 2x2, and why it is the pitfall

The pitfall block of 7.13 and the block-solve block both name Cramer's rule as the factor's unstable
way of handling a `2 x 2` block, against the LU the solve uses. This subsection writes both out on one
example, because seeing the arithmetic is what shows where the accuracy goes.

Cramer's rule solves a linear system by writing each unknown as a ratio of determinants. For the
`2 x 2` case, take

```
E z = y,   E = [ a  b ]   y = [ y1 ]   z = [ z1 ]
               [ c  d ]       [ y2 ]       [ z2 ]
```

with `det(E) = a*d - b*c`. Cramer's rule gives

```
     y1*d - b*y2               a*y2 - c*y1
z1 = -----------,       z2 = -----------
      a*d - b*c                a*d - b*c
```

Each numerator is the determinant of `E` with the matching column overwritten by `y`: `z1`'s numerator
is `E` with its first column replaced by `y`, and `z2`'s is `E` with its second.

**The two solutions, on one example.** Take `E = [[2, 1], [1, 3]]` and `y = [5, 10]`. Both methods
return `z = [1, 3]`; they differ only in what they divide by.

Cramer, dividing by `det`:

```
det = 2*3 - 1*1 = 5
z1 = (5*3 - 1*10)/5 = 5/5  = 1
z2 = (2*10 - 1*5)/5 = 15/5 = 3
```

LU with partial pivoting, dividing by matrix entries:

```
Column 1: |2| >= |1|, so 2 is already the pivot, no row swap.

Eliminate row 2:   l21 = 1/2 = 0.5          <- divide by a11 = 2
   row2 - 0.5*row1:  a22 -> 3 - 0.5*1 = 2.5
                     y2  -> 10 - 0.5*5 = 7.5

(the a22 update is the factorization, the y2 update is the forward solve; they are
 interleaved here, done in one pass over row 2)

Now  [ 2   1 ] [z1]   [ 5  ]
     [ 0  2.5] [z2] = [7.5]

Back-substitute:
   z2 = 7.5 / 2.5 = 3               <- divide by u22 = 2.5
   z1 = (5 - 1*3) / 2 = 2/2 = 1     <- divide by a11 = 2

btw, above we have the following LU factorization:

E = L U:   [ 2  1 ]   [  1    0 ] [ 2   1  ]
           [ 1  3 ] = [ 0.5   1 ] [ 0  2.5 ]
```

Cramer makes one division, by `det`. LU divides only by matrix entries, `2` and `2.5`, and partial
pivoting arranges that the entry divided by is the larger of the column.

**Cost: LU against the determinant way, and why it does not decide.** The `2 x 2` here is the small
case of a `k x k` block (AGL's `k <= 5`, 7.13's block-solve block), and at block scale the two methods
cost the same order, so cost is not what chooses between them. It helps to separate the two kinds of
cost, arithmetic and search, because only one of them can differ in order and it is not the one the
determinant way is blamed for.

*Arithmetic.* Both are `O(k^3)` for the block, the same leading order. Solving `k x k` is about
`k^3/3` either way: LU factors then substitutes, and the determinant way done sanely is that same
elimination expressed as an explicit inverse, `O(k^3)` to build `E^-1` and `O(k^2)` per right-hand
side to apply it. The difference is a constant factor, not an order. The one trap is Cramer taken
*literally*, a fresh determinant per unknown by cofactor expansion, which is factorial and which
nobody means: the real determinant way is elimination with a shared denominator, exactly the `2 x 2`
code's one `det` and adjugate over it, and that is the same order as LU.

*Search.* This is the only place a genuine order difference can enter, and it is LU's cost, not the
determinant way's. Pivot selection adds comparisons: partial pivoting scans a column, `O(k)` per step
and `O(k^2)` over the block; complete pivoting scans the whole trailing block, `O(k^2)` per step and
`O(k^3)` over the block, the order of the arithmetic itself. So partial-pivot search is a lower order
than the flops and disappears into the constant, while complete-pivot search is the same order as the
flops, can roughly double the constant, and worse on real hardware breaks the blocking. The
determinant way has no search at all, which is exactly why it is cheaper to write and why it was
tempting.

At these block sizes even the constant factor is beside the point. For `k = 2`, or AGL's `k <= 5`,
`k^3` is a rounding error against the sparse factorization around it, so neither the flop constant nor
the partial-pivot search term is measurable in the total. What separates the two methods is therefore
not cost but stability: same order, same ballpark constant, and Cramer exposed on an ill-conditioned
block where pivoted LU is not. The honest summary is that there is no order difference, a
constant-factor arithmetic difference, and a search term that is sub-leading for partial pivoting and
leading for complete, none of it large at `k <= 5`, which leaves stability as the only axis that
actually decides.

Two refinements bear on the sparse setting specifically. In the *solve* the block inverse is applied
once per pivot block, to the right-hand side; in the *factor* it is applied to every trailing row, as
`A21 E^-1`, so the apply cost there is `O(k^2)` times the number of update rows rather than `O(k^2)`
once, though still sub-leading against the rank-`k` trailing update it feeds. And complete pivoting
inside the block would buy the block a better growth bound, but the block is `k <= 5` and its
acceptance test already bounds what matters, so it defends against almost nothing, the same
museum-piece verdict as 7.10 at `k = 5` instead of `k = n`. Partial pivoting in the block is the right
level of care; complete would be the same over-insurance, one level down.

Why this is the pitfall. The `det` sits in every denominator, and it is computed as `a*d - b*c`, a
subtraction of two products. When `a*d` and `b*c` are both large and nearly equal, that subtraction
loses relative accuracy to cancellation, and nothing in the formula arrests it, unlike a pivoted LU,
which reorders to divide by the larger entry. AGL's specific objection (their 3.3) is that the
acceptance test bounds `|L|` but not the condition number of the block, so `E` is allowed to be
ill-conditioned, and for an ill-conditioned block a normwise backward stable solve, their Gaussian
elimination with complete pivoting, keeps the computed answer solving a nearby system, while Cramer's
rule and explicit inversion carry no such guarantee. Concretely, a near-singular block
`E = [[1, 1], [1, 1.0000001]]` has `det = 1e-7` formed by subtracting two numbers equal to seven or
eight digits, so the denominator arrives with most of its significant figures already gone, and
everything divided by it inherits that.

**Which path each phase takes.** The solve is on the LU path: `diagonalDynamic` factors each `2 x 2` by
the LU with partial pivoting above. The factor is on the Cramer path: when `factor2x2` eliminates a
trailing row against the pivot, it forms that row's two multipliers by the explicit-inverse formula
`(t1*d22 - t2*d21)/det` and `(t2*d11 - t1*d12)/det`, with `t1, t2` the row's two entries in the pivot
columns, dividing by `det`. That is Cramer row by row, so the factor is the half still on the pitfall
while the solve is off it. One note on symbols, since it is easy to trip on: the `l21 = 0.5` in the LU
above is the `(2,1)` entry inside the LU of the `2 x 2` block itself, a single number; the multipliers
just named are entries of the *outer* factor `L` below the pivot, a different object at a different
level, despite the shared digits.

One distinction is worth keeping, because it separates this from the block solve of 7.13. The Cramer
danger is a small `det`, not a zero diagonal. The exchange block `[[0, 1], [1, 0]]` has `det = -1`, so
Cramer handles it cleanly, `z = [y2, y1]`; it is a *symmetric* solve that breaks there, dividing by the
zero diagonal. The two failure modes are different: pivoted LU in the solve covers both, while the
factor's Cramer form covers neither hard case by design, which is why it is the one still sitting on
the pitfall.

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
- J. R. Gilbert, E. G. Ng and B. W. Peyton, "An efficient algorithm to compute row and
  column counts for sparse Cholesky factorization", *SIAM J. Matrix Anal. Appl.*
  15(4):1075-1091, 1994. Column counts in nearly linear time in `nnz(A)`, using the skeleton
  matrix and least common ancestors, rather than the `O(nnz(L))` pruned-row-subtree walk of 2.5.

**Symbolic factorization, ordering, sparse direct methods (Section 3; general).**

- A. George and J. W. H. Liu, *Computer Solution of Large Sparse Positive Definite
  Systems*, Prentice-Hall, 1981. Symbolic factorization, ordering, and the classic
  treatment of the whole pipeline.
- A. George, "Nested dissection of a regular finite element mesh", *SIAM J. Numer. Anal.*
  10(2):345-363, 1973. Nested dissection, the ordering of the grid example in 2.7.
- T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006. A modern,
  code-level reference for the entire order / symbolic / numeric pipeline. Its CSparse
  `cs_etree` computes the tree from the upper triangle of each column, the standard route
  used by Oblio (2.4); single-triangle storage in the other orientation is transposed first.

**Ordering (Section 5).** The two orderings we vendor, MMD and AMD, come from this lineage; the
primary sources for each, in the order Section 5 builds them.

- W. F. Tinney and J. W. Walker, "Direct solutions of sparse network equations by optimally
  ordered triangular factorization", *Proc. IEEE* 55(11):1801-1809, 1967. The origin of minimum
  degree: eliminate the vertex of least degree next, greedily.
- A. George and J. W. H. Liu, "The evolution of the minimum degree ordering algorithm", *SIAM
  Review* 31(1):1-19, 1989. **The one to read first.** It explains why a naive minimum degree is
  unusable and how every trick in a real implementation (the quotient graph, element absorption,
  supervariables, mass elimination, incomplete degree update) exists to make it affordable. Read
  this before touching AMD's source: the code is long because these mechanisms are, not because it
  is badly written.
- J. W. H. Liu, "Modification of the minimum-degree algorithm by multiple elimination", *ACM
  Trans. Math. Software* 11(2):141-153, 1985. Multiple minimum degree, which is what our vendored
  MMD is.
- P. R. Amestoy, T. A. Davis, and I. S. Duff, "An approximate minimum degree ordering algorithm",
  *SIAM J. Matrix Anal. Appl.* 17(4):886-905, 1996. AMD. The idea is to bound the degree instead
  of computing it, which turns the expensive step into a cheap one and, surprisingly, does not
  hurt the ordering.
- P. R. Amestoy, T. A. Davis, and I. S. Duff, "Algorithm 837: AMD, an approximate minimum degree
  ordering algorithm", *ACM Trans. Math. Software* 30(3):381-388, 2004. The software paper for the
  SuiteSparse code we vendor (`src/Amd.cpp`); the 1996 entry above is the algorithm, this is the
  implementation.
- T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006, chapter 7. Davis's own
  compact AMD (`cs_amd`), several times shorter than the SuiteSparse implementation. **The
  shortness is in the surrounding machinery, not the algorithm**: statistics, control parameters,
  input validation, debug dumps, and the unsymmetric preprocessing all fall away, while the
  quotient graph with garbage collection and supervariable hashing remains, because it is
  irreducible.

**Supernodes and amalgamation (Section 4).**

- C. Ashcraft and R. G. Grimes, "The influence of relaxed supernode partitions on the multifrontal
  method", *ACM Trans. Math. Software* 15(4):291-309, 1989. Relaxed (amalgamated) supernodes: the
  trade of explicitly stored zeros for larger dense blocks, which is what 4.5 is about.
- E. G. Ng and B. W. Peyton, "Block sparse Cholesky algorithms on advanced uniprocessor
  computers", *SIAM J. Sci. Comput.* 14(5):1034-1056, 1993. The supernodal Cholesky the numeric
  phase implements, and the case for why level-3 BLAS on dense blocks is the whole point.

**The multifrontal method (mentioned in 2.9, not yet implemented here).**

- I. S. Duff and J. K. Reid, "The multifrontal solution of indefinite sparse symmetric linear
  systems", *ACM Trans. Math. Software* 9(3):302-325, 1983. The method, and its update stack.

**Dense kernels.**

- J. J. Dongarra, J. Du Croz, S. Hammarling, and I. Duff, "A set of level 3 basic linear algebra
  subprograms", *ACM Trans. Math. Software* 16(1):1-17, 1990. Level-3 BLAS: `gemm`, `syrk`,
  `herk`, `trsm`. The reason supernodes are worth forming at all.
- E. Anderson et al., *LAPACK Users' Guide*, 3rd ed., SIAM, 1999. `potrf` (Cholesky) and `sytrf`
  (Bunch-Kaufman). Note there is **no unpivoted LDL in LAPACK**, which is why a static LDL kernel
  has to be written by hand.

**Multifrontal (Section 6).**

- B. M. Irons, "A frontal solution program for finite element analysis", *Internat. J. Numer.
  Methods Engrg.* 2(1):5-32, 1970. The frontal method: one active front in core, the rest
  streamed, which is the out-of-core origin of 6.3.
- I. S. Duff and J. K. Reid, "The multifrontal solution of indefinite sparse symmetric linear
  equations", *ACM Trans. Math. Software* 9(3):302-325, 1983. The generalization from one front
  to many over the elimination tree.
- J. W. H. Liu, "The multifrontal method for sparse matrix solution: theory and practice",
  *SIAM Review* 34(1):82-109, 1992. The standard survey, and the source for the frontal and
  update matrix vocabulary of 6.1.
- J. W. H. Liu, "On the storage requirement in the out-of-core multifrontal method for sparse
  factorization", *ACM Trans. Math. Software* 12(3):249-264, 1986. The storage formula of 6.5 and
  the child-ordering rule of 6.6, including the exchange argument.
- P. R. Amestoy, I. S. Duff, J.-Y. L'Excellent, and J. Koster, "A fully asynchronous multifrontal
  solver using distributed dynamic scheduling", *SIAM J. Matrix Anal. Appl.* 23(1):15-41, 2001.
  MUMPS, and the type 1, 2 and 3 node classification of 6.4.
- J.-Y. L'Excellent and W. M. Sid-Lakhdar, "A study of shared-memory parallelism in a multifrontal
  solver", *Parallel Computing* 40(3-4):34-46, 2014. The tree-against-node parallelism tradeoff
  of 6.4, and where the boundary between them is placed.

**Pivoting (Section 7).**

- J. R. Bunch and B. N. Parlett, "Direct methods for solving symmetric indefinite systems of
  linear equations", *SIAM J. Numer. Anal.* 8(4):639-655, 1971. Complete pivoting for the
  symmetric indefinite case: the strong bound, and the `O(n^3)` search that made a cheaper
  strategy necessary.
- J. R. Bunch and L. Kaufman, "Some stable methods for calculating inertia and solving symmetric
  linear systems", *Math. Comp.* 31(137):163-179, 1977. The `1 x 1` / `2 x 2` selection of 7.2,
  its four cases, and the constant `(1 + sqrt(17)) / 8`.
- I. S. Duff and J. K. Reid, "The multifrontal solution of indefinite sparse symmetric linear
  equations", *ACM Trans. Math. Software* 9(3):302-325, 1983. Listed above for the multifrontal
  method; it is also where the threshold conditions of 7.5 and the delayed columns of 7.6 come
  from, the two arriving together and for the same reason.
- I. S. Duff, "MA57: a code for the solution of sparse symmetric definite and indefinite systems",
  *ACM Trans. Math. Software* 30(2):118-144, 2004. The successor to MA27, and where the threshold
  conditions are stated in the form 7.5 uses, with the default `u = 0.01`.
- C. Ashcraft, R. G. Grimes, and J. G. Lewis, "Accurate symmetric indefinite linear equation
  solvers", *SIAM J. Matrix Anal. Appl.* 20(2):513-561, 1998. Why Bunch-Kaufman's unbounded `L`
  matters (7.3), bounded Bunch-Kaufman (rook pivoting), and delayed pivoting in a sparse
  factorization, which is why dynamic LDL grows a front at runtime. It is also the source of Oblio's
  own pivoting (7.13): their Figure 3.3 is the `2 x 2` acceptance test the non-root kernel
  evaluates, and their Figure 3.4 is its search procedure.
- I. S. Duff and S. Pralet, "Strategies for scaling and pivoting for sparse symmetric indefinite
  problems", *SIAM J. Matrix Anal. Appl.* 27(2):313-340, 2005. The pivoting strategies MUMPS
  uses, and the scaling that precedes them.
- I. S. Duff and S. Pralet, "Towards stable mixed pivoting strategies for the sequential and
  parallel solution of sparse symmetric indefinite systems", *SIAM J. Matrix Anal. Appl.*
  29(3):1007-1024, 2007. What changes when the candidate set is restricted further still by
  distribution across processes, which is 7.4's restriction taken a step past ours.
- I. S. Duff, J. Hogg, and F. Lopez, "A new sparse LDL^T solver using a posteriori threshold
  pivoting", *SIAM J. Sci. Comput.* 42(2):C23-C42, 2020. A modern alternative to all of the
  above: factor optimistically and check afterwards, rather than testing each candidate before
  accepting it.
- N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, 2002, chapter 11.
  The growth bounds quoted in 7.2 and the error analysis behind 7.3. Also, for the *static* path,
  the perturbation of a small pivot in a factorization that cannot pivot and so has no other
  recourse.

**On the framing.** Two presentational choices here are expository, not lifted from a
single source: casting the forest as the *transitive reduction of the update DAG* (2.6)
is a folklore-standard restatement rather than any one paper's language, and structuring
the 2.4 correctness proof around the "Cholesky factors nest" lemma (Lemma A) is this
document's packaging of the incremental argument. Both are consistent with the sources
above; the specific wording is ours.
