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
| NumFactorDynamic | yes | checked | 0.9 `FactorsDynamic`. One index vector and one value vector per supernode, so a front can expand under delayed pivoting. Written by every factorization: the static ones run into it unchanged and produce a factor identical to the flat one, bit for bit; dynamic LDL writes it through the expansion and contraction verbs (`expandNodeIdx`, `resetVal`, `expandVal`, `swap`,
`contractVal`). No base class shared with the static one: `experiments/storage-options` showed a pointer array does the job a base would |
| NumFactorEngine | yes | checked | **Static factorizations functionally complete**: `Cholesky`, `StaticLDLT`, `StaticLDLH`, each left- and right-looking, real and complex. Cholesky checked against an independent dense Cholesky (4e-16); LDL by reconstruction, `L D L^H == P A P^T` (2e-15), through AMD ordering and supernodes. `StaticLDLH` (complex Hermitian LDL) is an **extension**: 0.9's complex LDL is symmetric only. Gaps, both in Owed: the LDL **perturbation branch has never fired**, and a complex input is **not validated as Hermitian**. **Dynamic LDL works for real input in both traversals**, `DynamicLDLT` and `DynamicLDLH` alike (the
same computation over the reals), with delayed columns crossing the forest and 2x2 pivots in the
solve, verified by residual in `test_pipeline` and by right-looking agreeing with left-looking bit
for bit. 0.9's two dynamic kernels are byte-identical between its left- and right-looking engines,
so only the driver differs: right-looking expands a front with `expandVal`, which preserves, where
left-looking uses `resetVal`, which discards. Complex dynamic LDL and `Traversal::Multifrontal` not
started. **Complex `DynamicLDLT` needed no kernel change at all**: 0.9's complex
`factorDynamicLDL_` differs from its real one in six lines, all declaring the pivot magnitudes real
rather than scalar, and this port declared them `double` from the start, so it was already the
complex form; `updateDynamicLDL_` is byte-identical between 0.9's two engines. Only the dispatch
guard had to widen. **Complex `DynamicLDLH` is done too, and it is an extension rather than a port**: 0.9's complex LDL
is symmetric only, so nothing here was transcribed and its oracles are the residual and
reconstruction of `L D L^H`. It needed the conjugate in `readPivotBlock2x2`, conjugated `L` where the
`D L^H` rows are formed, `forceReal` on the diagonal, and one fix in `swap` described below |
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

- **The input is not checked for Hermitian symmetry, and this is a correctness hole.** Moved to
  docs/TODO.md, under "Validate the input matrix", where it joins the two other unchecked
  preconditions found since (a structurally present diagonal, and the absence of duplicate
  entries). All three want one validation pass, so they are one job rather than three, and it is
  not a porting job: neither reference validates its input either. The agreed fix, two flags on
  `SparseMatrix` computed in one construction pass, is recorded there.

- **The LDL perturbation branch executes, but nothing asserts it.** A static factorization cannot
  pivot, so a pivot smaller than the threshold is *replaced* and counted (`ldl`'s `n == 1` case;
  `NumFactorStatic::numPerturbations` reports it). That branch is the only part of static LDL that
  changes *what is computed* rather than how.

  This entry used to say the branch had never run, on the grounds that every test matrix was
  diagonally dominant. **That stopped being true on 2026-07-19** and the trigger turned out to be a
  family already in the suite: the banded matrices with zeroed diagonals that `test_pipeline` uses
  for tier 1. Measured on `bandIndefinite(40, 3, zf, 7)` under Natural ordering with `StaticLDLT`:

  | zero fraction | perturbations | residual |
  |---|---|---|
  | 0.00 | 0 | 1.9e-16 |
  | 0.10 | 1 | 1.1e-03 |
  | 0.30 | 1 | 1.4e-03 |
  | 0.50 | 1 | 3.6e-04 |

  Only one perturbation even at half the diagonal zeroed, for the same reason most zero diagonals
  never delay: they fill in from the Schur complement before they are reached. One column arrives
  still tiny.

  It was found twice by accident, both times as a "failing" assertion that was the branch working
  correctly on input a static factorization cannot handle. The test it still wants is unchanged in
  shape, assert the count is nonzero *and* that the reconstruction differs from `A` by about the
  perturbation, which is the honest statement of what perturbing means: we factored a slightly
  different matrix and said so. What has changed is that the matrix no longer has to be invented.

- **The two pivot bodies in `factorDynamicSupernode` were duplicated, and the debt is now paid.**
  0.9 splits `factorDynamicLDL_` on whether the supernode has update rows and writes both bodies
  out; the port followed it, deliberately, because those lines were 0.9's and merging them before
  they had ever run would have put their first execution in a shape the oracle never had.

  The trigger was the end-to-end residual, and once it fired the merge was done: the eliminations
  are now `applyPivot1x1` and `applyPivot2x2`, called from both selection loops, and the 2x2 block
  is read in one place, `readPivotBlock2x2`. 230 code lines became 204 across three functions, and
  the refactor was covered by 147 assertions throughout.

  Two things fell out of it that the duplication had hidden. **Pass 1 never reads the 2x2 block at
  all**, it accepts on `max1 == max2`, on the magnitudes alone, where pass 2 tests the determinant;
  the compiler said so with an unused-variable warning the moment the shared body stopped reading it
  for both. And `readPivotBlock2x2` turns out to be the single place the *symmetry* of D is decided,
  which is exactly what complex `LDL^H` needs to change.

  The selection loops stay separate, and that is a decision rather than unfinished work. See
  docs/TODO.md.

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
