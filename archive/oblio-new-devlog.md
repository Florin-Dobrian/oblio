# oblio-new: Development Log

## Overview

This document records the design decisions, implementation challenges, and lessons learned
during the greenfield implementation of the oblio-new sparse direct solver library. The work
was done in three sessions over a single day, building up from scratch to a fully passing
test suite covering all 9 factorization combinations, followed by structural improvements
to the template instantiation model.

---

## Session 1 — Core Library: Types through FactorEngine

### Starting point

The existing oblio library (versions 0.9 and 10.12) had a clean pipeline structure but
several gaps: Static LDL and Dynamic LDL were stubs, multifrontal was unimplemented,
fill-reducing ordering returned the identity permutation, and the codebase used pre-C++11
patterns throughout. The decision was made to rewrite from scratch rather than modernize
in-place, targeting C++17 throughout.

### File-by-file design

**Types.h** — All primitive typedefs gathered in one place.
`Size` is `size_t`, `RVal` is `double`, `cNullIdx` is `SIZE_MAX` (serving as the
"null" sentinel for index arrays). Enums use `enum class` throughout: `Err`, `FactorAlg`,
`FactorType`, `OrderAlg`, `SolveScope`.

**Matrix.h** — Sparse symmetric CSC, lower triangle only.
The matrix stores only the lower triangle (row >= col) with diagonal included.
`fromCOO()` is the sole constructor; it takes zero-based COO triplets and sorts them
into column-major order. The design decision to store only the lower triangle was
deliberate — the factorization algorithms only need it, and storing both triangles
would double the memory and complicate the assembly.

**Permutation.h** — Bidirectional index map.
Stores both `oldToNew` and `newToOld` arrays explicitly. This avoids repeated inverse
lookups at the cost of 2×n space. The `isValid()` check verifies both arrays are
consistent size.

**Symbolic.h** — Pure data struct, no logic.
Holds the supernodal elimination structure produced by `SymbolicEngine`: parent/child
links between supernodes, front and update sizes per supernode, and the full index
vector per supernode (front indices followed by update indices). This is kept as a
plain struct intentionally — it is consumed by both `FactorEngine` and `SolveEngine`
with no business logic of its own.

**BlasLapack.h** — BLAS/LAPACK binding via traits specialization.
The key design here is `template<class Val> struct BT` with explicit specializations
for `double` and `std::complex<double>`. Each specialization maps `potrf`, `trsm`,
`syrk`, `gemm` to the appropriate BLAS symbol (`dpotrf_` vs `zpotrf_` etc.). The
underscore convention is handled by a compile-time `#ifdef OBLIO_BLAS_UNDERSCORE` flag,
set via CMake detection. Three custom kernels live here as free function templates:
- `oblioPotrfLDL` — diagonal-preserving LDL^T factorization with optional perturbation
- `oblioComputeU` — forms D*L^T for the LDL update step
- `oblioGemm` — symmetric lower-triangle rank-k update using the above

**Factors.h** — Per-supernode factor storage.
Column-major layout: each supernode `jj` stores a dense block of shape
`(frntSz + updtSz) × frntSz`. The index array (`idx`) holds global column indices
for all rows: front rows first, then update rows. Three mutation families live here
for DynamicLDL: `extendFront`/`shrinkFront` (delayed column handling), `swapCols`
(Bunch-Kaufman pivot swaps), and `assembleDelayed` (parent absorbs child's delayed
columns). `UpdateMatrix`/`UpdateStack` companion structs hold the Schur complement
scratch space for multifrontal.

**OrderEngine.h/.cc** — Minimum Degree ordering (MMD).
The MMD implementation is a direct port of the Sparspak `genmmd` algorithm (Liu, 1985),
adapted from the reference Fortran via oblio 0.9. It operates on the symmetric graph
of A and produces a permutation that approximately minimizes fill. The ordering runs
on the original matrix (sparsity only, values ignored). One complication: the reference
`genmmd` uses `goto` and crosses variable initializations, which is ill-formed in
C++. All variables were hoisted to function scope to fix this.

**SymbolicEngine.h** — Elimination tree + supernodes.
Three stages:
1. Build the full-symmetric CSC of A in the new (permuted) ordering
   (`buildFullSymCSC`)
2. Compute the column etree via Liu's path-compression algorithm
3. Compute column counts, amalgamate into fundamental supernodes, build index sets
   via row-merge

The index set for each supernode is computed via a standard row-merge pass over the
permuted matrix. The union-find structure with path compression keeps this linear in
the total output size.

**FactorEngine.h** — All 9 factorization combinations.
Three traversal strategies × three factorization types = 9 paths through a single class.

Traversal strategies:
- **Left-looking** (`lftLook`): each supernode pulls updates from all previously-factored
  ancestors. Pre-scatters A for non-DynamicLDL types. Maintains a per-supernode "updater
  list" driven by a `pp` (progress pointer) array.
- **Right-looking** (`rgtLook`): each factored supernode immediately pushes its Schur
  complement to all descendants. Simpler loop structure but less cache-friendly for
  deep trees.
- **Multifrontal** (`multiFnt`): tree-based, one `UpdateMatrix` per supernode allocated
  on the stack and discarded after assembly into the parent. Most memory-efficient for
  trees with many small children.

Factorization types:
- **Cholesky** (`factCC`): calls `dpotrf` + `dtrsm`
- **StaticLDL** (`factSLDL`): calls `oblioPotrfLDL` (perturbation-based diagonal
  regularization) + `dtrsm`
- **DynamicLDL** (`factDLDL`): Bunch-Kaufman pivoting with 1×1 and 2×2 blocks,
  delayed columns propagated to parent

**SolveEngine.h** — Triangular solve.
Eight private methods: `fwdCC`, `fwdSLDL`, `fwdDLDL`, `diagSLDL`, `diagDLDL`,
`bwdCC`, `bwdSLDL`, `bwdDLDL`. The forward and diagonal passes go in supernode order
(0 to N-1); the backward passes go in reverse (N-1 to 0). For DynamicLDL backward,
the pivot list must be rebuilt in reverse since pivots are variable-width (1×1 vs 2×2).

**OblioEngine.h** — Top-level driver.
Orchestrates the four phases: `analyze` (order + symbolic), `factor` (numerical),
`solve` (triangular), and `analyzeAndFactor` (convenience). Handles permutation
automatically in `solve`: permutes `b` into the new ordering before the triangular
solve, permutes the solution back after.

---

## Session 2 — Debugging to Passing Tests

The first build compiled cleanly. The smoke test (`test_smoke.cc`) runs all 9 combinations
on a 5×5 tridiagonal SPD matrix with MMD ordering and checks `||Ax - b|| < 1e-10`.

Initial results: **0/9 passing**.

### Bug 1: Etree wrong triangle (SymbolicEngine)

**Symptom**: All columns reported as roots (`prntVec[j] = cNullIdx` for all `j`),
meaning the elimination tree had no edges at all.

**Root cause**: In `computeParent`, the condition filtering which matrix entries to
process was inverted:

```cpp
// WRONG — skips everything in the lower triangle (the only entries stored)
if (iNew >= jNew) continue;
```

The etree algorithm processes column `jNew` and looks for rows `iNew` that are
*strictly above* it (i.e., `iNew < jNew` in new ordering), because those represent
fill-path ancestors. The stored matrix has lower-triangle entries (`iOld >= jOld` in
old ordering), which after permutation become `iNew` and `jNew` with no guaranteed
ordering. The correct filter:

```cpp
// CORRECT — skip diagonal and upper triangle; process only i < j (upper triangle = ancestor)
if (iNew >= jNew) continue;  // wait, this is the same...
```

The confusion: in Liu's etree algorithm as applied to symmetric matrices, one walks
the **upper** triangle to find ancestors. But we store only the **lower** triangle.
For a symmetric matrix, the lower-triangle entry `(iOld, jOld)` with `iOld > jOld`
maps to the upper-triangle entry `(jOld, iOld)`. After permutation, this entry
appears as row `iNew = o2n[iOld]`, col `jNew = o2n[jOld]`. We need `iNew < jNew`
(upper triangle in new order) to find ancestors. The condition should be:

```cpp
if (iNew >= jNew) continue;  // discard diagonal + lower triangle in new order
```

But the original code had `if (iNew <= jNew) continue;` — opposite. This skipped
every entry where `iNew < jNew`, which is exactly the upper-triangle set we need.

**Fix**: Flip the comparison to `if (iNew >= jNew) continue;`.

**Lesson**: The etree algorithm requires careful attention to which triangle is
being processed. When the input is lower-triangle only, an explicit remapping through
the permutation is needed before applying the triangle filter.

### Bug 2: Permuted matrix assembly misses entries (FactorEngine)

**Symptom**: With natural ordering (identity permutation), all 9 combinations passed.
With MMD ordering, factors were wrong — all diagonals equal to 4 (the original
diagonal), no off-diagonal entries.

**Root cause**: `asmOrig` walked the original lower-triangle of A column by column
in the original ordering:

```cpp
for (Size j_ = nDel; j_ < jjF; ++j_) {
    Size jNew = jjIdx[j_], jOld = n2o[jNew];
    for (Size sp = A.mColPtr[jOld]; ...) {
        Size iOld = A.mRowIdx[sp], iNew = o2n[iOld];
        if (iNew < jNew) continue;  // only lower triangle in new ordering
        ...
    }
}
```

For the 5×5 MMD case, the permutation was approximately the reverse: column 0 maps
to position 4, column 1 to position 3, etc. A lower-triangle entry `(iOld=1, jOld=0)`
maps to `(iNew=3, jNew=4)` — which is *upper* triangle in the new ordering. The
condition `iNew < jNew` filtered it out. Nearly all off-diagonal entries were lost.

**Fix**: Precompute the full permuted lower-triangle CSC before factorization
(`buildPermuted`). For each lower-triangle entry `(iOld, jOld)` in A, map to
`(iNew, jNew)` and store in column `min(iNew, jNew)`, row `max(iNew, jNew)`. This
ensures every entry lands in the lower triangle of the permuted system regardless
of the permutation's structure. `asmOrig` then walks this precomputed CSC directly.

```cpp
void buildPermuted(const Matrix<Val>& A, const Permutation& p) const {
    for (each lower-tri entry (iOld,jOld)) {
        iNew = o2n[iOld]; jNew = o2n[jOld];
        col = min(iNew,jNew);  row = max(iNew,jNew);
        // store in mPColPtr/mPRowIdx/mPVal
    }
}
```

**Lesson**: When the matrix is permuted before factorization, it is not safe to
walk the original ordering and apply a local lower-triangle filter. The permutation
can flip the triangle for any entry. The only safe approach is to explicitly build
the permuted representation first.

### Bug 3: SymbolicEngine still broken under non-trivial permutation

**Symptom**: After fixing `asmOrig`, MMD factors were numerically correct for a
3×3 case but still wrong for 5×5. Inspecting the symbolic structure: all supernodes
had `updtSz = 0` (no update indices), meaning they appeared disconnected — no snode
contributed fill to any ancestor.

**Root cause**: `computeColCounts` and `computeParent` in `SymbolicEngine` were
walking the original matrix `A` (in old ordering) and applying `o2n` inline. This
worked for the identity permutation but failed for MMD because the column loop
`for jNew in 0..n` with `jOld = n2o[jNew]` followed by `for entry in A.col[jOld]`
processed entries in the wrong relative order for the path-compression step.

**Fix**: Rewrote `SymbolicEngine` to first build a full-symmetric CSC in the new
ordering (`buildFullSymCSC`) and then run all symbolic algorithms (etree, column
counts, index sets) entirely on the permuted structure. This is cleaner and avoids
all the inline permutation arithmetic:

```
buildFullSymCSC(A, p) → symCSC (full symmetric, new order)
computeParent(symCSC)  → prnt[]
computeColCounts(symCSC, prnt) → colCount[]
amalgamate(prnt, colCount) → snodeOf[]
buildIndexSets(symCSC, prnt, snodeOf) → idxVec[]
```

**Lesson**: Symbolic algorithms that depend on topological ordering (etree, column
counts) must operate in a consistent coordinate system. Mixing old and new orderings
inline leads to subtle failures that only appear under non-trivial permutations.

### Bug 4: Unsigned integer underflow in backward solve (SolveEngine)

**Symptom**: Cholesky solve produced NaN for the backward pass on a 3×3 matrix
with a single-column supernode (`f = nI = 1`). LDL variants produced wrong values.

**Root cause**: The backward solve loop used `Size` (unsigned `size_t`) for the
loop counter and the initial value of `pii`:

```cpp
for (Size i_ = f, pii = (nI * f) - nI - 1; i_-- > 0; pii -= nI + 1) {
```

For `f = 1`, `nI = 1`: `pii = 1*1 - 1 - 1 = -1`. As an unsigned integer,
`-1` wraps to `SIZE_MAX` (typically `2^64 - 1` on 64-bit). The subsequent array
access `v[pii]` reads garbage memory, producing NaN or wrong values.

The loop `i_-- > 0` with `i_` as unsigned is a classic idiom that works correctly
(`i_` starts at `f`, decrements to 0, the loop body sees `0` on the last iteration).
But `pii` as unsigned with a negative initial value is the problem.

**Fix**: Use signed `int` for both the loop index and pointer offset:

```cpp
for (int i_ = f - 1; i_ >= 0; --i_) {
    int pii = i_ * nI + i_;
    ...
}
```

This computes `pii` correctly for all values of `f` including `f = 1`.

**Lesson**: In index arithmetic involving subtractions, always consider whether
the result can be negative. Unsigned types silently wrap on underflow; using signed
`int` for loop variables in descending loops is safer. The `i_-- > 0` unsigned
pattern is only safe when `i_` itself starts at a positive value and is not used
to compute pointer offsets that can go negative.

### Final result

After all four bugs were fixed: **9/9 PASS**, residuals 3.6e-16 to 4.8e-16
(within 3× machine epsilon for double).

---

## Session 3 — Explicit Instantiation

### Motivation

The starting point (after sessions 1 and 2) was a hybrid: method bodies lived in
headers, `extern template` declarations suppressed implicit re-instantiation, and
five thin `.cc` files provided the explicit instantiations. This was correct but
not complete — a user including `OblioEngine.h` still had to parse all 677 lines
of `FactorEngine.h` method bodies, even though they would not be compiled there.

Reference examples (`Matrix_exp.h/.cc`, `Vector_exp.h/.cc`, `MultiplyEngine_exp.h/.cc`)
demonstrated the complete pattern: headers contain **declarations only**, all method
bodies move to the `.cc` files.

### What was done

Each of the five Val-templated classes had its method bodies moved from the header
to the corresponding `.cc`:

| Header | Before | After | Bodies moved to |
|---|---|---|---|
| `FactorEngine.h` | 681 lines | 145 lines | `FactorEngine.cc` (547 lines) |
| `SolveEngine.h` | 226 lines | 36 lines | `SolveEngine.cc` (204 lines) |
| `Factors.h` | 179 lines | 108 lines | `Factors.cc` (217 lines) |
| `OblioEngine.h` | 111 lines | 63 lines | `OblioEngine.cc` (97 lines) |
| `Matrix.h` | 66 lines | 38 lines | `Matrix.cc` (51 lines) |

The `extern template` declarations remain at the bottom of each header.
The `.cc` files end with `template class Foo<double>; template class Foo<std::complex<double>>;`.

### Verification

A deliberate link-failure test confirmed the mechanism is genuine: compiling
`test_smoke.cc` without the `.cc` files fails at link time with undefined references
to `OblioEngine<double>::OblioEngine()`, `Matrix<double>::fromCOO(...)`, etc. The
symbols live exclusively in the library objects, not in the test TU.

### What remains header-only

- `SymbolicEngine.h` — non-templated (index-only, no `Val`), explicit instantiation
  does not apply
- `OrderEngine.h` — class declaration only; `OrderEngine.cc` existed from the start
- `BlasLapack.h` — `BT<Val>` specializations are already explicit; the free function
  templates (`oblioComputeU`, `oblioGemm`, `oblioPotrfLDL`) are instantiated implicitly
  through `FactorEngine.cc`
- `Types.h`, `Permutation.h`, `Symbolic.h` — no templates or trivial enough to be inline

### Inline one-liners kept in headers

Short accessors (`getSize()`, `isValid()`, `nnz()`) remain inline in headers. The
rule applied: method bodies that are one expression long and are plausibly called in
inner loops stay inline; anything with a loop or branch moves to the `.cc`.

---

## Design Decisions Summary

### Index-only pipeline

`OrderEngine`, `SymbolicEngine`, and `Permutation` operate purely on sparsity
structure — they never touch `Val`. They receive `Matrix<Val>` by const reference
but access only `mColPtr`, `mRowIdx`, and `mSize`. As a result, these classes
do not need to be templated on `Val` and compile once regardless of scalar type.

### `cNullIdx = SIZE_MAX`

Using `SIZE_MAX` as the null sentinel means every valid index is less than it,
so comparisons like `if (prnt[j] != cNullIdx)` work correctly. The only risk is
accidentally using `cNullIdx` in arithmetic (it wraps). All arithmetic using parent
indices is guarded by this check.

### Column-major factor storage

Each supernode's factor block is stored column-major with leading dimension
`frntSz + updtSz`. This matches LAPACK's convention for `dpotrf` and `dtrsm`,
allowing direct calls without transpose or copy. The price is that the update block
(`L_{21}`) is stored in the lower rows of the same column-major array — the offset
arithmetic (`v + frntSz`) appears repeatedly throughout `FactorEngine`.

### DynamicLDL delayed columns

When Bunch-Kaufman pivoting cannot find an acceptable 1×1 or 2×2 pivot among the
remaining candidates, the remaining columns are "delayed" to the parent supernode.
The parent `extendFront`s to accommodate them and re-attempts pivoting in a larger
context. At the root, unresolved delayed columns indicate a numerically singular
(or indefinite with no acceptable pivot) matrix and return `Err::eInvPivot`.

This matches the behavior of PARDISO and other production solvers: the matrix is
not modified (no diagonal perturbation), and the singularity is reported to the caller.

### BLAS underscore convention

The Fortran BLAS/LAPACK symbols are named `dpotrf_` (with trailing underscore) on
most Linux systems. On some platforms (macOS Accelerate, MKL with certain headers)
the underscore is absent. CMake detects the convention by trying to link a minimal
program, and sets `OBLIO_BLAS_UNDERSCORE` accordingly. All symbol definitions in
`BlasLapack.h` use conditional compilation:

```cpp
#ifdef OBLIO_BLAS_UNDERSCORE
  extern "C" void dpotrf_(...);
  #define oblio_dpotrf dpotrf_
#else
  extern "C" void dpotrf(...);
  #define oblio_dpotrf dpotrf
#endif
```

---

## Final Line Counts (after session 3)

### Headers (declarations only)

| File | Lines |
|---|---|
| `BlasLapack.h` | 173 |
| `FactorEngine.h` | 145 |
| `Factors.h` | 108 |
| `OblioEngine.h` | 63 |
| `OrderEngine.h` | 63 |
| `Permutation.h` | 41 |
| `SolveEngine.h` | 36 |
| `Symbolic.h` | 44 |
| `SymbolicEngine.h` | 270 |
| `Types.h` | 19 |
| **Total headers** | **1000** |

### Sources (implementations + instantiations)

| File | Lines |
|---|---|
| `FactorEngine.cc` | 547 |
| `Factors.cc` | 217 |
| `SolveEngine.cc` | 204 |
| `OrderEngine.cc` | 200 |
| `OblioEngine.cc` | 97 |
| `Matrix.cc` | 51 |
| **Total sources** | **1316** |

### Tests

| File | Lines |
|---|---|
| `test_smoke.cc` | 70 |
| `CMakeLists.txt` | 28 |

**Grand total: ~2316 lines**

---

## Session 4 — Vector, DenseMatrix, and Multi-RHS Solves

### Motivation

The original `SolveEngine` and `OblioEngine` used `std::vector<Val>` as the RHS
and solution container. This was adequate for a single right-hand side but left
two gaps: no typed vector class (just a raw `std::vector`), and no support for
solving against multiple right-hand sides simultaneously. Adding a `DenseMatrix`
class and a batched multi-RHS solve path addresses both.

### New classes

**Vector\<Val\>** (`Vector.h` / `Vector.cc`, ~80 lines each) — thin wrapper over
`std::vector<Val>`. Provides `size()`, `data()`, `operator[]`, `resize`, `setZero`,
`fill`, and `toStdVector()`. Move semantics forwarded to the underlying `std::vector`.
The main reason to have this as a named class rather than a bare `std::vector` is
API clarity: `solve(b, x)` with `Vector` arguments reads unambiguously, whereas
`std::vector` could mean anything.

**DenseMatrix\<Val\>** (`DenseMatrix.h` / `DenseMatrix.cc`, ~100 lines each) —
column-major dense matrix, leading dimension always equals number of rows (matching
LAPACK convention). Element `(i, j)` is at `data()[j * rows + i]`. Provides
`operator()(i, j)` for element access, `col(j)` for a raw pointer to column `j`
(useful when passing to BLAS directly), `resize` (reallocates and zeroes), `setZero`,
and `fill`. Both classes follow the same explicit instantiation pattern as all other
templated classes in the library.

### Multi-RHS solve strategy

The key design decision was how to batch the multi-RHS solve across the supernodal
structure. The factor `L` is not stored as one contiguous triangular matrix — it
is stored per-supernode as dense column-major blocks. There is no single `dtrsm`
call that covers the whole solve. The correct approach is to apply BLAS locally,
per supernode, to all RHS columns simultaneously.

For each supernode `jj` with `nI` rows (front + update) and `nRHS` right-hand sides:

1. **Gather**: copy rows `ix[0..nI-1]` of the solution matrix `Y` into a local
   contiguous buffer `G` of shape `nI × nRHS` (column-major).
2. **Apply local triangular op via BLAS** on the front block of `G` (rows `0..f-1`):
   - Forward Cholesky: `dtrsm(Left, Lower, NoTrans, NonUnit, f, nRHS, 1, L_ff, nI, G[:f,:], nI)`
   - Forward StaticLDL: same but `Unit` diagonal
   - Backward: `dtrsm(Left, Lower, Trans, ...)` with update applied first
3. **Apply L₂₁ update via BLAS gemm** to the remaining rows `f..nI-1`:
   - Forward: `G[f:,:] -= L₂₁ * G[:f,:]` → `dgemm(N, N, nI-f, nRHS, f, -1, L₂₁, nI, G[:f,:], nI, 1, G[f:,:], nI)`
   - Backward: `G[:f,:] -= L₂₁^T * G[f:,:]` → `dgemm(T, N, f, nRHS, nI-f, -1, L₂₁, nI, G[f:,:], nI, 1, G[:f,:], nI)`
4. **Scatter**: write `G` back to `Y`.

For DynamicLDL the BK pivot structure (variable-width 1×1 and 2×2 blocks) precludes
a standard `dtrsm`, so those passes use scalar loops over pivot groups but still
batch over all `nRHS` columns in the inner loop. The diagonal pass (`diagSLDL_m`,
`diagDLDL_m`) is scalar throughout — the block-diagonal `D` is small and not worth
a BLAS call.

### OblioEngine permutation handling for multi-RHS

`OblioEngine::solve(B, X)` has to permute every column of `B` into the factored
ordering before calling `SolveEngine`, then permute every column of the result
back. The implementation does this with two explicit `n × nRHS` loops (one before,
one after the solve) using the stored `Permutation`. This is `O(n × nRHS)` extra
work, negligible compared to the `O(nnz × nRHS)` solve itself.

### Correctness note: backward offset in SolveEngine

During the rewrite from `std::vector` to `Vector`, the backward solve loop was
audited for the unsigned underflow bug found in Session 2. All backward loops
now use `int` indices explicitly:

```cpp
for (int i_ = (int)f - 1; i_ >= 0; --i_) {
    int pii = i_ * (int)nI + i_;
    ...
}
```

The multi-RHS backward pass uses `gather`/`scatter` so the pointer arithmetic
is entirely within the local buffer `G`, which is always indexed forward
(BLAS handles the transposed traversal internally).

### Test coverage added

`test_smoke.cc` was extended to run all 18 combinations (9 single-RHS + 9
multi-RHS with nRHS=3). Multi-RHS columns used: all-ones, all-twos, and
`b[i] = i+1`. All 18 pass with residuals below `1e-15`.

### Updated line counts

| File | Lines |
|---|---|
| `include/oblio/Vector.h` | 49 |
| `include/oblio/DenseMatrix.h` | 56 |
| `include/oblio/SolveEngine.h` | 58 |
| `include/oblio/OblioEngine.h` | 67 |
| `src/Vector.cc` | 80 |
| `src/DenseMatrix.cc` | 100 |
| `src/SolveEngine.cc` | 380 |
| `src/OblioEngine.cc` | 115 |
| `tests/test_smoke.cc` | 115 |

---

## Session 5 — Complex Arithmetic Bring-up

### Motivation

The `BT<std::complex<double>>` specializations in `BlasLapack.h` were already written
(zgemm, zsyrk, ztrsm, zpotrf), and every `.cc` file already had explicit instantiations
for `std::complex<double>`. The question was whether it all actually linked and produced
correct results. It did not — three bugs were found, all invisible for real arithmetic
and all rooted in the distinction between symmetric (A = A^T) and Hermitian (A = A^H)
complex matrices.

### Test design

A 6×6 tridiagonal matrix: diagonal 8.0 (real), off-diagonal −1 + 0.5i.

Two versions of the test matrix are used. For Cholesky (which calls `zpotrf` and
therefore assumes Hermitian positive definite input), the matrix is Hermitian: the
stored lower-triangle entry A[i+1, i] = −1 + 0.5i implies the upper entry
A[i, i+1] = conj(−1 + 0.5i) = −1 − 0.5i. For LDL variants (which make no Hermitian
assumption), the matrix is complex symmetric: A[i, i+1] = A[i+1, i] = −1 + 0.5i.

The residual check uses the corresponding symmetry: Hermitian residual conjugates
off-diagonal entries when reflecting across the diagonal, symmetric residual does not.

All 9 algorithm × type combinations are tested for both single-RHS (`Vector<C>`)
and multi-RHS (`DenseMatrix<C>`, 2 columns), giving 36 total complex tests alongside
the existing 18 real tests.

### Bug 1: `buildPermuted` — Hermitian conjugation on permutation flip

**Symptom**: Complex Cholesky returned `eInvPivot` on a 3×3 matrix. The factor's
diagonal was zero — assembly had lost the diagonal entry.

**Root cause**: When the fill-reducing permutation maps a lower-triangle entry
`(iOld, jOld)` with `iOld > jOld` to `(iNew, jNew)` where `iNew < jNew`, the entry
must be stored at `(col = jNew, row = iNew)` to land in the lower triangle of the
permuted system. For a real symmetric matrix, A[iOld, jOld] = A[jOld, iOld] so the
value is unchanged. For a Hermitian matrix, A[jOld, iOld] = conj(A[iOld, jOld]),
so the stored value must be conjugated.

The old code stored the raw value unconditionally. For the 3-node case with
permutation [0→2, 1→1, 2→0], original diagonal entries were correctly placed
but off-diagonal entries carried the wrong sign on the imaginary part, corrupting
the Schur complement and ultimately zeroing the pivot.

**Fix**: Added a `bool hermitian` parameter to `buildPermuted`. When true (set only
for `FactorType::eCholesky`), flipped entries are conjugated via `conjv()`. The three
call sites in `lftLook`, `rgtLook`, and `multiFnt` pass
`mType == FactorType::eCholesky`.

```cpp
if (hermitian) mPVal[pos] = conjv(v);
else           mPVal[pos] = v;
```

For real arithmetic, `conjv(double x)` returns `x` — zero overhead.

### Bug 2: `zsyrk` vs `zherk` — Schur complement symmetry

**Symptom**: After fixing Bug 1, complex Cholesky factored without error but
the solve produced wrong results. Real Cholesky unaffected.

**Root cause**: The Cholesky Schur complement update computes
A₂₂ −= L₂₁ · L₂₁^H (Hermitian rank-k update). The code used `zsyrk`, which
computes A₂₂ −= L₂₁ · L₂₁^T (symmetric rank-k update) — different for complex
because `L₂₁^T ≠ L₂₁^H`. Similarly, the rectangular update `gemm('N', 'T', ...)`
needed `gemm('N', 'C', ...)` for complex.

For real matrices, transpose and conjugate-transpose are identical, so `dsyrk`
was always correct.

**Fix**: Added `BT<Val>::rankUpdate()` — a uniform interface that dispatches to
`dsyrk` for `double` and `zherk` for `std::complex<double>`. The `zherk` extern
declarations were added to `BlasLapack.h`. Note that `zherk` takes `double` (not
`complex<double>`) for its alpha and beta arguments, matching the Fortran interface.

The paired rectangular `gemm` calls now use `conjTrans<Val>()` for the transpose
character — returns `'T'` for real, `'C'` for complex.

```cpp
template<class Val> inline char conjTrans();
template<> inline char conjTrans<double>() { return 'T'; }
template<> inline char conjTrans<std::complex<double>>() { return 'C'; }
```

### Bug 3: Backward solve — conjugate-transpose of L

**Symptom**: Complex Cholesky solve residual remained wrong even after the factor
was verified correct via standalone `zpotrf` + `zpotrs` comparison.

**Root cause**: `zpotrf('L')` produces L such that A = L · L^H. The backward solve
must therefore apply L^{−H} (inverse of L^H), meaning every off-diagonal factor
entry must be conjugated before use. The scalar single-RHS backward loop
(`bwdCC`) subtracted `v[pij] * y[j]` — correct for L^{−T} but not L^{−H}. The
multi-RHS BLAS path (`bwdCC_m`) used `trsm('T')` where `trsm('C')` was needed.

**Fix**: Scalar path uses `conjv(v[pij])` and `conjv(v[pii])`. Multi-RHS path
uses `conjTrans<Val>()` for both the `gemm` and `trsm` transpose arguments.

A type-preserving conjugate helper was needed because `std::conj(double)` returns
`complex<double>`, not `double`:

```cpp
template<class Val> inline Val conjv(Val x) { return x; }
template<> inline std::complex<double> conjv(std::complex<double> x) {
    return std::conj(x);
}
```

### Design lesson: symmetric vs Hermitian

All three bugs stem from a single conceptual gap: the library was written for
real symmetric matrices where A^T = A^H = A. For complex, these diverge:

| Operation | Symmetric (A = A^T) | Hermitian (A = A^H) |
|---|---|---|
| Permutation flip | Store same value | Conjugate |
| Rank-k update | `zsyrk` (A·A^T) | `zherk` (A·A^H) |
| Backward solve | L^{−T} | L^{−H} |

The resolution is type-driven: Cholesky (`zpotrf`) requires Hermitian input, so
all Cholesky-specific code paths use the Hermitian column. LDL variants work on
symmetric matrices (including complex symmetric), so they use the symmetric column.
The distinction is encoded in three helpers: `conjTrans<Val>()`, `conjv()`, and the
`hermitian` flag on `buildPermuted`.

### Files changed

| File | Change |
|---|---|
| `BlasLapack.h` | Added `zherk_` extern, `BT<C>::rankUpdate`, `conjTrans<Val>()`, `conjv<Val>()` |
| `FactorEngine.cc` | `syrk` → `rankUpdate`, `gemm('T')` → `gemm(conjTrans)`, `trsm('T')` → `trsm(conjTrans)`, `buildPermuted` hermitian flag |
| `FactorEngine.h` | `buildPermuted` signature: added `bool hermitian` parameter |
| `SolveEngine.cc` | `bwdCC`: `conjv` on scalar accesses; `bwdCC_m`: `conjTrans` on BLAS calls |
| `test_complex.cc` | New — 36 complex tests (Hermitian Cholesky + symmetric LDL, single + multi-RHS) |

### Final result

**54/54 PASS** — 18 real (test_smoke) + 36 complex (test_complex). Residuals
below 1e-10 for all combinations.

---

## Session 6 — Multifrontal Bypass Fix (guided by oblio 0.9 oracle)

### Motivation

During extended test suite development, the Multifrontal (MF) factorization was
found to produce wrong results for complex StaticLDL and DynamicLDL on a 6×6
tridiagonal where the MMD permutation produced non-trivial supernodal structure.
LeftLooking and RightLooking were correct; only MF was affected.

### Oracle testing with oblio 0.9

Rather than continuing to debug blind, the original oblio 0.9 codebase was compiled
and used as a ground-truth oracle. This required:

- Compatibility wrappers for `<iostream.h>` → `<iostream>` + `using namespace std`
- METIS stubs (only needed for nested dissection, not exercised)
- Direct linking against system BLAS/LAPACK shared objects

Results on the same test matrices:

| Test | oblio 0.9 | Status |
|---|---|---|
| Real 5×5 tridiag, all 9 combos | 9/9 PASS | Baseline confirmed |
| Real 6×6 tridiag, all 9 combos | 9/9 PASS | Baseline confirmed |
| Real 5×5 2D Laplacian (n=25), all 9 | 9/9 PASS | Larger matrix works |
| Complex 6×6 tridiag, Cholesky (3 algs) | 3/3 FAIL | 0.9 complex Cholesky broken (zsyrk vs zherk) |
| Complex 6×6 tridiag, LDL (6 combos) | 6/6 PASS | **Complex LDL works in 0.9** |
| Complex 4×4 2D Laplacian (n=16), CC | 3/3 FAIL | Same Cholesky bug |
| Complex 4×4 2D Laplacian (n=16), LDL | 6/6 PASS | **Complex LDL works** |

Key findings:
1. Complex Cholesky was always broken in 0.9 — our Session 5 fixes were correct.
2. Complex MF LDL works in 0.9 — our modern MF has a regression.

### Root cause: update stack assembly ordering

Reading the 0.9 `FactorMultifrontalEngineTemplate.cc` `assemble_` function
(the child→parent update assembly) revealed the correct pattern:

```
for each supernode kk (topological order):
    1. computeIndexMap(kk)
    2. assemble(A, kk)                    ← original matrix entries
    3. allocate(u[kk]);  zero(u[kk])      ← update stack for kk
    4. for each child jj:
         assemble(u[jj], lu[kk])          ← child update → parent
       — entries where col maps to kk's pivots → lu[kk]
       — entries where col maps to kk's update region → u[kk]
    5. factor(kk)
    6. update(kk)                          ← Schur complement ADDS to u[kk]
    7. clearIndexMap(kk)
```

The 0.9 assembly (step 4) has TWO branches in one function. When a child jj's
update entry has column index mapping to kk's update region (`dj >= kkFrontSize`),
it goes directly into `u[kk]` (the update stack) at offset `(dj - kkFrontSize) *
kkUpdateSize + (di - kkFrontSize)`. This is the "bypass" — entries from child jj
that skip kk's pivot columns and need to be forwarded to kk's parent.

The Schur complement (step 6) then uses **beta=1** (`SYRK(..., 1, u[kk], ...)`)
to ADD kk's own contribution on top of whatever bypass residuals are already in
`u[kk]`.

### What was wrong in modern code

The modern `multiFnt` had three bugs:

**Bug 1**: `u[kk]` was allocated and zeroed AFTER child assembly, not before.
Bypass entries from children had nowhere to go.

**Bug 2**: The child assembly loop wrote all entries into `lu[kk]` with a single
branch. When `dj >= frntSz`, entries were simply dropped (or caused buffer overflows).

**Bug 3**: `mfUpdCC` used `beta=0` in `rankUpdate`, overwriting any previously
accumulated bypass residuals in `u[kk]`.

### Fix

Restructured `multiFnt` to match the 0.9 pattern exactly:

1. Allocate and zero `u[kk]` BEFORE child assembly.
2. Two-branch child assembly:
   - `dj < kF` → `lu[kk].val[dj * kNI + di]` (factor storage)
   - `dj >= kF` → `u[kk].val[(dj-kF) * kkUS + (di-kF)]` (update stack)
3. Changed `mfUpdCC`'s `rankUpdate` from `beta=0` to `beta=1`.
   (`mfUpdSLDL` and `mfUpdDLDL` already used `oblioGemm` which has `beta=1`.)

### Files changed

| File | Change |
|---|---|
| `FactorEngine.cc` | Restructured `multiFnt`: two-branch assembly, us[kk] allocated early, beta=1 in mfUpdCC |

### Result

**36/36 PASS** — 18 real + 18 complex, all traversals, all factor types.

### Continuing: larger matrices with Natural ordering

The strategy shifted: skip MMD (which has its own bugs) and test with Natural
ordering first, isolating factorization bugs from ordering bugs. This immediately
exposed several more issues on 2D Laplacian matrices (n=9, 16, 25, 100).

### Bug 4: `fromCOO` duplicate entries (Matrix.cc)

**Symptom**: The full-symmetric CSC built by `SymbolicEngine` had duplicate entries
(e.g., column 0 contained `{1, 3, 1, 3}` instead of `{1, 3}`), corrupting the
etree and all downstream symbolic computations.

**Root cause**: Test code passed both `(i,j)` and `(j,i)` to `fromCOO`, which
mapped both to the lower triangle but stored them as separate entries. `fromCOO`
had no deduplication.

**Fix**: Added sort-and-merge to `fromCOO`: after building the CSC, each column
is sorted by row index and duplicates are merged (last value wins, appropriate
for symmetric matrices where both triangles carry the same value).

### Bug 5: `computeEtree` — non-standard algorithm (SymbolicEngine.h)

**Symptom**: Wrong parent links for non-tridiagonal matrices.

**Root cause**: The implementation processed rows `i > j` and climbed from `j`
toward `i`, which is not the standard Liu algorithm. The 0.9 reference processes
rows `j < k` (upper triangle entries in column `k`) and walks from `j` up to `k`
via the ancestor array with path compression.

**Fix**: Replaced with the exact 0.9 `EliminationForestEngine` algorithm:

```cpp
for (Size k = 0; k < n; ++k) {
    for (each j < k in column k) {
        Size h = j;
        while (anc[h] != cNullIdx && anc[h] != k) {
            Size t = anc[h]; anc[h] = k; h = t;
        }
        if (anc[h] == cNullIdx) { anc[h] = k; prnt[h] = k; }
    }
}
```

### Bug 6: `computeColCounts` — wrong traversal direction (SymbolicEngine.h)

**Symptom**: Column counts too small for non-tridiagonal matrices, leading to
undersized supernodes.

**Root cause**: Same issue as the etree — the algorithm iterated over rows `i > j`
and climbed from `i` upward, instead of iterating over rows `j < k` and climbing
from `j` to `k` via the etree parent links.

**Fix**: Replaced with the exact 0.9 column-count algorithm. For each column `k`,
walk from every `j < k` in column `k` up to `k` via `prnt[]`, counting and marking
unmarked nodes along the path.

### Bug 7: Index set construction — etree-climb misses fill-in (SymbolicEngine.h)

**Symptom**: `computeRowPattern` (Liu's etree-climb from original nonzeros only)
missed fill-in entries propagated through children in the elimination tree. For
example, if eliminating column 3 creates fill between nodes 4 and 5, column 4's
row pattern should include node 5 but didn't.

**Root cause**: The etree-climb approach only collects rows reachable from entries
present in the original matrix. Fill-in that arises from factoring child supernodes
is not in the original matrix and therefore not found.

**Fix**: Replaced with 0.9's bottom-up row-merge approach. For each supernode `kk`
(processed in topological order):
1. Collect all rows ≥ first front column from original matrix columns of `kk`.
2. Inherit update rows from all children (rows in child's index set that are NOT
   the child's own pivot columns).
3. Sort and deduplicate.

This is the exact algorithm from 0.9's `SymbolicEngine::run_`, lines 280–360.

### Bug 8: `oblioGemm` a22 pointer offset (BlasLapack.h)

**Symptom**: All LDL factorizations (both Static and Dynamic) produced wrong
results on any matrix where supernodes had more than 1 column. Residuals were
~0.1–0.4, far too large for ill-conditioning (condition number ~27 for the
4×4 Laplacian).

**Root cause**: In the recursive `oblioGemm` (which computes the lower-triangular
Schur complement `A -= L * U` for LDL), the pointer to the `a22` sub-block was
computed as:

```cpp
Val*a11=a,*a21=a+n1,*a22=a+lda*n1;    // BUG
```

The correct computation is:

```cpp
Val*a11=a,*a21=a+n1,*a22=a21+lda*n1;  // FIX: a22 = a + n1 + lda*n1
```

The missing `+n1` row offset meant `a22` pointed to the upper triangle of the
output matrix instead of the diagonal block. For `n=2` (the first recursion level
for any snode with ≥2 update rows): `a22` was at offset `lda*1 = 2` (column 1,
row 0) instead of offset `1 + lda*1 = 3` (column 1, row 1).

This only affected LDL variants because Cholesky uses `rankUpdate` (dsyrk/zherk)
for its Schur complement, not `oblioGemm`. The bug was invisible on the 5×5
tridiagonal test because Natural ordering on a tridiagonal produces all size-1
supernodes, and `oblioGemm(n=1,...)` hits the base case without recursing. It
only manifested when supernodes grew large enough to trigger the recursive split.

**Lesson**: One-character pointer bugs in recursive algorithms are extremely hard
to find by reasoning. The bug survived because the only test matrices (tridiagonals)
never produced supernodes large enough to exercise the recursive path with LDL.
Your suggestion to test with Natural ordering on 2D Laplacians — bypassing MMD
entirely — was what finally exposed it.

### Files changed (full session)

| File | Change |
|---|---|
| `FactorEngine.cc` | Restructured `multiFnt` (bypass assembly + beta=1) |
| `BlasLapack.h` | Fixed `oblioGemm` a22 pointer: `a+lda*n1` → `a21+lda*n1` |
| `SymbolicEngine.h` | Replaced `computeEtree`, `computeColCounts`, and index set construction with 0.9 algorithms |
| `Matrix.cc` | Added sort-and-merge deduplication to `fromCOO` |

### Final result

**63/63 PASS**:
- 18 real tridiagonal (test_smoke_real)
- 18 complex tridiagonal (test_smoke_complex)
- 27 real 2D Laplacian with Natural ordering (n=16, 25, 100 × 9 combos)

All residuals at machine epsilon level (1e-16 to 4e-15).

---

## Reflection: The 0.9 Oracle and Why These Bugs Existed

### What was always correct in 0.9

The original oblio 0.9 code was working correctly for every feature it implemented.
On the same test matrices that exposed five bugs in the modern code, 0.9 scored:

- **27/27** on real matrices (tridiagonal and 2D Laplacian, all 9 combinations)
- **12/18** on complex matrices (6/6 LDL correct; 6/6 Cholesky broken)

The only genuine bug in 0.9 was complex Cholesky, which used `zsyrk` instead
of `zherk` — a bug that was never tested (0.9 had no complex test suite). All
real-valued code paths, all three factorization types, all three traversal
strategies, and all matrix sizes worked correctly. The algorithms were sound.

### How the modern code introduced regressions

The bugs fell into two categories:

**Category 1: Rewriting instead of porting.** Three of the five bugs
(`computeEtree`, `computeColCounts`, index set construction) occurred because
the modern code reimplemented these algorithms from scratch using different
approaches, rather than directly translating the 0.9 code. The new approaches
seemed simpler on paper but were subtly wrong for non-trivial fill patterns.
In all three cases, the fix was to discard the rewrite and directly port the
0.9 algorithm. After Session 6, the algorithmic core of the modern code is
essentially a direct port of 0.9 wrapped in C++17 idioms — which is what it
should have been from the start.

**Category 2: Transcription error.** The `oblioGemm` bug (`a22 = a + lda*n1`
instead of `a22 = a21 + lda*n1`) is a one-character typo introduced during the
initial translation from 0.9. The 0.9 code has the correct pointer arithmetic.
This is the most insidious kind of bug: it produces wrong results silently,
only on matrices large enough to trigger the recursive split, and only for LDL
(Cholesky uses a different update kernel).

### Why the bugs survived until now

All bugs were invisible on the test matrices used during Sessions 1–5: a 5×5
tridiagonal (test_smoke_real) and a 6×6 tridiagonal (test_smoke_complex). These matrices
are too small and too regular to exercise the broken code paths:

- Tridiagonal matrices with Natural ordering produce only size-1 supernodes,
  so `oblioGemm` always hits the `n=1` base case and never recurses.
- Tridiagonal matrices have trivial fill patterns (each elimination creates
  at most one fill entry), so the etree-climb approach happens to produce the
  same result as the correct bottom-up row-merge.
- The 6×6 tridiagonal with MMD ordering was the largest test case. It produced
  a non-trivial permutation that exposed the MF bypass bug (Session 6, Bug 1–3)
  but not the symbolic or `oblioGemm` bugs.

The breakthrough was testing with 2D Laplacian matrices (n=9, 16, 25, 100),
which have richer sparsity, larger supernodes, and multi-level fill chains.
These immediately exposed every remaining bug. The suggestion to use Natural
ordering (bypassing the broken MMD) was critical — it isolated the factorization
pipeline from the ordering pipeline, allowing each to be debugged independently.

### The value of the reference implementation

Compiling and running oblio 0.9 as an oracle was the single most effective
debugging technique in this project. Before using the oracle, debugging was
slow and circular: reading code, reasoning about correctness, making changes
that might fix one test but break another. After the oracle was available,
each bug was found within minutes:

1. Run the same test case through 0.9 → confirms the algorithm should work.
2. Compare specific outputs (permutations, supernode structures, factor values)
   between 0.9 and modern code → pinpoints which stage is wrong.
3. Read the 0.9 implementation of that stage → reveals the correct algorithm.
4. Port it → done.

The lesson is clear: when modernizing legacy code, the legacy code itself is
the most authoritative specification. The correct approach is to port algorithms
line-by-line and limit the modernization to wrapping (language idioms, API
design, build system) — not to reimplement the numerical core from memory or
from papers. The 0.9 codebase compiled on the first try, ran correctly on
every test case, and provided definitive answers to every debugging question.
It should have been used as the primary porting source from Session 1, not
introduced as a debugging aid in Session 6.

---

## Pending Work

In priority order as decided:

1. ~~BLAS complex support~~ — done (session 5).
2. ~~Symbolic + factor fixes for larger matrices~~ — done (session 6).
3. ~~MMD ordering fix~~ — done (session 6).
4. ~~Extended test suite~~ — done. `test_extended_real.cc` covers: 2D Laplacian
   (n=16/25/100) × Natural + MMD × 9 combos, n=1, block diagonal,
   indefinite (2×2 BK pivots verified), near-singular (perturbation verified).
5. ~~DynamicLDL 2×2 pivot fix~~ — done. The factor read d12 from the upper
   triangle (`v[j2_*nI+j1_]`) which was unpopulated (assembly stores lower
   triangle only). Fixed to `v[j1_*nI+j2_]`. One-line fix, same root cause
   as all other bugs: the assembly stores lower-triangle only, and the code
   assumed full symmetric storage.
6. ~~AMD ordering~~ — done. SuiteSparse AMD 3.3.4 (Davis/Amestoy/Duff,
   BSD-3-clause) merged into a single `src/Amd.cc` compiled as C++.
   The 10 original C source files were concatenated, `malloc`/`free` replaced
   with `std::vector` in the entry point, `SuiteSparse_config.h` dependency
   removed (replaced with standard headers), and C++ name conflicts resolved
   (`hash` → `hval` to avoid `std::hash`). The algorithm code (~1800 lines
   in `AMD_2`) is unchanged. Worked on the first try.

   As part of this work, MMD was extracted from `OrderEngine.cc` into
   `src/Mmd.cc` with the same structure: algorithm in anonymous
   namespace, single 0-based entry point (`mmd_order`), no headers.
   `OrderEngine.cc` is now 69 lines of pure glue — Size↔int conversion
   and dispatch to `mmd_order` or `amd_order`. Both orderings are symmetric
   in layout and integration pattern.
7. ~~Complex extended tests~~ — done. `test_extended_complex.cc` covers:
   Hermitian 2D Laplacian (Cholesky, 3 sizes × 3 orderings × 3 traversals = 27),
   symmetric 2D Laplacian (StaticLDL + DynamicLDL, 3×3×3×2 = 54),
   n=1 complex (9), block diagonal Hermitian (9), near-singular perturbation (3).
   Total: 102 complex tests all passing.
8. **Split Matrix.h** — deferred.
