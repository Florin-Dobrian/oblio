# Porting Ledger

Tracks the migration of each 0.9 unit into the modern tree. After a context gap,
read this first, it turns "refresh my context" into a two-minute scan.

## How to use

- One row per unit (a class or a small function-group).
- Status, in increasing order of confidence:
  - `not started`
  - `ported`, translated and compiles, but nothing checks that it is right.
  - `checked`, verified against an independent oracle: an implementation that
    shares no code path with the unit, so agreement is evidence rather than
    tautology. For the structural units this is a dense simulation of Cholesky
    fill, which checks the mathematical property directly.
  - `verified`, output compared against 0.9 itself on real input, buffer for
    buffer.
- `checked` and `verified` are different kinds of evidence, not merely different
  amounts. An oracle says the unit computes the right thing; a 0.9 comparison says
  it computes the same thing 0.9 does. A unit can be `checked` and still differ
  from 0.9 in a way that matters downstream (index order within a block, say), so
  `verified` remains the bar before numeric factorization consumes a unit.
- Notes: index-signedness fixes, departures from 0.9 or 10.12, open questions,
  anything a future session needs.

## Pipeline

`OrderEngine -> ElmForestEngine -> SymFactorEngine -> NumFactorEngine -> SolveEngine`,
orchestrated by `DirectSolver`.

Naming note: the modern tree renames as it ports. 0.9's `Matrix` is `SparseMatrix`
(a `DenseMatrix` will follow), and 0.9's `Symbolic` / `SymbolicEngine` pair is
`SymFactor` / `SymFactorEngine`, with `NumFactor` / `NumFactorEngine` to follow for the
numeric phase. The table gives the modern name, with the 0.9 name in the notes where they
differ.

## Units

| Unit | Uses Val? | Status | Notes |
|---|---|---|---|
| Types | no | ported | sentinels (`NIL`), index typedefs |
| Permutation | no | checked | index-only. Ported: `set` (as `setOldToNew`/`setNewToOld`, replacing 0.9's direction flag) and `compose`. Not yet ported from 0.9: `get`, `read`/`write` (persist an ordering), `initialize2dGrid`/`initialize3dGrid` (structured orderings, useful as test inputs with hand-computable fill) |
| SparseMatrix | yes | checked | 0.9 `Matrix`; flat CSC, stored fully (both triangles). Build and structural symmetry tested |
| OrderEngine | no | checked | AMD (SuiteSparse) and MMD (Sparspak) vendored verbatim; output checked for validity as a permutation, not against 0.9's ordering |
| ElmForest | no | checked | data; supernodal shape throughout, nodal being the uncoarsened case rather than a special one |
| ElmForestEngine | no | checked | parent links, child/sibling links, roots, height, column sizes, fundamental compression, threshold amalgamation. Links and height recomputed independently; sizes and supernodes against the dense oracle, natural and AMD ordered. Amalgamation is greedy and not canonical, so only its tie-break-invariant properties are asserted |
| SymFactor | no | checked | 0.9 `Symbolic`; flat index sets with per-supernode offsets |
| SymFactorEngine | no | checked | 0.9 `SymbolicEngine`; index sets against the dense oracle, natural and AMD ordered. 10.12's design, 0.9's behavior (see DESIGN_DECISIONS) |
| BlasLapack | yes | checked | wraps potrf/trsm/herk/syrk/gemm, overloaded on the scalar type. Named by *operation*, not routine: `herk` means A times A-conjugate-transpose, so it is `dsyrk_` for real and `zherk_` for complex, and the engine cannot pick the wrong one (0.9 does; see DESIGN_DECISIONS). One trait, `Blas<Val>::conjTrans`. Also carries the three kernels BLAS lacks, ported from 0.9: `ldl` (unpivoted LDL, 0.9's `OBLIO_POTRF2`), `formUpper` (`U = D L^T`, `OBLIO_COMPUTE_U`), `gemmLower` (`A -= L U` with the product known symmetric, `OBLIO_GEMM`). Verified on hand-computed factors and by reconstruction, 1x1 to 23x23 |
| UpdateBlock | yes | checked | 0.9 `Temporary`; one supernode's update to one ancestor, dense column-major plus its row indices. Not the multifrontal update matrix, which is a different object |
| NumFactorStatic | yes | checked | 0.9 `FactorsStatic`; SymFactor's structure copied, plus one flat value buffer with per-supernode offsets. Blocks are dense column-major rectangles (the upper front triangle is allocated and zero, so BLAS can take the whole block) |
| NumFactorDynamic | yes | checked | 0.9 `FactorsDynamic`. One index vector and one value vector per supernode, so a front can grow under delayed pivoting. Written by every factorization: the static ones run into it unchanged and produce a factor identical to the flat one, bit for bit; dynamic LDL writes it through the growth verbs (`extendIndex`, `swap`, `shrinkEntry`). No base class shared with the static one: `experiments/storage-options` showed a pointer array does the job a base would |
| NumFactorEngine | yes | checked | **Static factorizations functionally complete**: `Cholesky`, `StaticLDLT`, `StaticLDLH`, each left- and right-looking, real and complex. Cholesky checked against an independent dense Cholesky (4e-16); LDL by reconstruction, `L D L^H == P A P^T` (2e-15), through AMD ordering and supernodes. `StaticLDLH` (complex Hermitian LDL) is an **extension**: 0.9's complex LDL is symmetric only. Gaps, both in Owed: the LDL **perturbation branch has never fired**, and a complex input is **not validated as Hermitian**. `DynamicLDLT`/`DynamicLDLH` and `Traversal::Multifrontal` not started |
| Vector | yes | checked | 0.9 `SingleVector`; one column. 0.9 also has `MultipleVector`, whose solve uses TRSM/GEMM with a gather and scatter; with one right-hand side there is no level-3 BLAS to be had, so the scalar solve is right and the multi-column path is a later performance addition |
| MultiplyEngine | yes | checked | `y = A x` and `r = A x - b`. Exists for the residual: it is what turns per-phase oracles into an end-to-end check |
| DenseMatrix | yes | not needed so far | a supernode's block is a raw pointer plus (rows, cols, ld), handed straight to BLAS. See Owed |
| SolveEngine | yes | checked | forward, diagonal (LDL only), backward, for Cholesky and both static LDLs, real and complex. Scalar, one right-hand side, as 0.9's `SingleVector` path is. **The backward pass conjugates when the factorization does**, which 10.12 omits: its backward solve applies `L^-T` where a Hermitian factor needs `L^-H`, correct for its complex-symmetric LDL and wrong for its Cholesky. Verified by residual, `\|Ax - b\| / \|b\|` at 3e-16 through the whole pipeline |
| OblioEngine | yes | ported as `DirectSolver` | top-level driver: owns the permutation, forest, symbolic and numeric factors, and exposes analyze / factor / solve. Renamed because it is not an engine in this tree's sense (engines are stateless and produce one object) and because `Oblio` names the package, not the method: `DirectSolver` pairs with an `IterativeSolver` if one ever arrives |

Units from 0.9 deliberately not carried over: `Utility` (`ResizeVector` and
friends, obviated by `std::vector`) and `Functional` (pre-C++11 comparators,
obviated by lambdas).

## Owed

### Numeric factorization

- **The input is not checked for Hermitian symmetry, and this is a correctness hole.** For
  complex Cholesky the matrix must be Hermitian (`A[i][j] == conj(A[j][i])`). Nothing enforces
  it. `zpotrf` reads only the lower triangle and *assumes* the upper is its conjugate, and
  `assembleFromA` takes the same triangle, so a complex **symmetric** matrix would be factored as
  though it were the Hermitian matrix agreeing with its lower triangle, and would **succeed**,
  returning a plausible wrong answer. Silent, and complex-only.

  The fix, agreed: compute two flags on `SparseMatrix` at construction, in one pass, and have the
  engine require the first for Cholesky.

  ```
  isHermitian()   A == A^H    what Cholesky needs
  isSymmetric()   A == A^T    what complex LDL will need
  ```

  For `double` the two predicates coincide and are both true together; for complex they are
  genuinely different and a matrix may be one, the other, or neither. Compute both now, since it
  is the same pass and one extra comparison per entry, and LDL then needs no change to
  `SparseMatrix`. 10.12 carries only `mIsNmrclySym`, which sufficed for it because it never
  finished complex Cholesky. 0.9 checks only that the diagonal is real, and only in its
  `NO_LAPACK` path, so its LAPACK build does not check at all.

  Deferred deliberately: correct input gives correct output today, so this guards against misuse
  rather than fixing a defect. It must land before anyone but us runs the solver.

- **The LDL perturbation branch has never executed.** A static factorization cannot pivot, so a
  pivot smaller than the threshold is *replaced* and counted (`ldl`'s `n == 1` case;
  `NumFactorStatic::numPerturbations` reports it). That branch is the only part of static LDL that
  changes *what is computed* rather than how, and every test matrix so far is diagonally dominant,
  so no pivot has ever been small enough to trigger it. It wants a matrix with a deliberately tiny
  pivot, asserting both that the count is nonzero and that the reconstruction then differs from
  `A` by about the perturbation, which is the honest statement of what perturbing means: we
  factored a slightly different matrix, and said so.

- **Statistics are not computed.** `setSymFactor` derives the value-block sizes it needs, but the
  aggregate counts 0.9 keeps (`numberOfEntries`, `numberOfAllocatedEntries`, the multifrontal
  stack high-water marks, the factor and solve flop estimates) are unported. 0.9 computes them in
  `EliminationForest::computeStatistics_`; 10.12's `rComputeStats` is commented out. The stack
  counts become necessary for multifrontal; the rest are reporting.

- **No `DenseMatrix`.** The ledger lists one as a unit. It has not been needed: a supernode's
  block is a raw pointer plus (rows, columns, leading dimension), handed straight to BLAS, which
  is what `experiments/storage-options` argued for and what keeps the kernels blind to the
  storage. Revisit only if something wants a dense matrix as an object.

### Structural

- **Nothing is `verified` yet.** Every structural unit is `checked` against the
  dense oracle, which is real evidence, but none has been compared against 0.9
  buffer for buffer. That comparison is owed before numeric factorization starts
  consuming these structures, because a difference the oracle cannot see (index
  order within a supernode's block, say) is exactly the kind that surfaces as a
  numeric bug much later.
- **`sortForOptimalMultifrontal`** (0.9) / `rOptimizeForMultifrontal` (10.12) is not ported.
  Both references call it after compression to reorder children for the multifrontal stack;
  10.12's call site is commented out. Not needed until the numeric factorization exists.
