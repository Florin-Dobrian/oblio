# Porting Ledger

Tracks the migration of each 0.9 unit into the modern tree. After a context gap,
read this first, it turns "refresh my context" into a two-minute scan.

Two axes, tracked separately because they move separately. **What** a unit computes is the port,
and the tables below carry it unit by unit. **How** it is written is the idiom, sweeps that run
across the whole tree at once, listed in the two sections before the tables. A unit can be
`verified` against 0.9 and still hold a 1998 idiom, and the reverse.

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

## Idiom modernization: C++11

The tree compiles `-std=c++17`, and the jump being made is from **pre-C++11**. 0.9 predates C++11
entirely, so porting a unit means carrying its algorithm over unchanged, which is the invariant,
while its idioms come forward.

C++11 is the large jump and the one to finish. Each entry below gives what was there, what replaced
it, and every site it was applied at, so a later reader can tell whether a file has been through
this and does not have to re-derive it by reading.

**Swept 2026-08-13**, in commit `387689d`:

- **A default constructor written only to zero scalars** -> a default member initializer at each
  declaration, plus `= default` on the constructor. Eight classes:
  `experiments/friend-access/{Matrix,Vector}.h` and
  `experiments/template-instantiation/{Matrix,Vector}{Implicit,PlainExplicit,GuardedExplicit}.h`.
  The six bodies in the matching `.cpp` files went with them.
- **A default member initializer made dead by an initializer-list entry** -> delete one of the two.
  Keep the declaration's only where a constructor exists that does not set the member, which is why
  `SparseMatrix` and `UpdateMatrix` correctly carry both. One site:
  `include/oblio/UpdateBlock.h`, `mHeight` and `mWidth`, unreachable behind the only constructor.
- **`std::vector<T>().swap(v)` to force deallocation** -> assignment from a fresh object,
  `*this = T()`. `clear()` does not free and `shrink_to_fit` is only a request, which is why the
  1998 idiom existed. One site: `src/UpdateMatrix.cpp`, `discard()`. Multifrontal peak measured
  identical to the byte before and after.
- **Accumulating in the initializer list to avoid a body** -> a loop in the body, seeded by a
  default member initializer. One site: `experiments/storage-options/SparseMatrixDynamic.cpp`, the
  constructor's `mNnz`.
- **Default construct, then assign field by field** -> aggregate initialization at the construction
  site. One site: `experiments/storage-options/test_multiply.cpp`, `Columns` in `build`.

**Swept earlier**, during the port of each unit rather than as a pass:

- **`NULL`, `typedef`, iterator loops, `0` for null** -> `nullptr`, `using`, range-for. Applied by
  `clang-tidy`'s `modernize-*` checks as each unit came over.
- **Raw owning arrays, manual `new` / `delete`** -> `std::vector`, the default container. This one
  is a CLAUDE.md invariant rather than a convention.

The rules that came out of these live in `CODING_RULES.md`, under C++; the reasoning and the
two defects the sweep exposed are in the 2026-08-13 entry of `notes/DESIGN_DECISIONS.md`.

Nothing here is known to remain, which is weaker than nothing remains: the sweep followed an audit
of four sites plus two idioms found while there, not a systematic search. A `clang-tidy` run with
the `modernize-*` checks over `src/` and `include/` would settle it either way and has not been
done.

## Idiom modernization: C++17

Available to use, and no sweep is owed. C++17 adds far less to this codebase than C++11 did, and the
parts that mattered are already adopted where they earn their place:

- **`if constexpr`**, `src/SolveEngine.cpp`, where a branch is impossible rather than merely
  untaken: the dynamic passes do not exist for a static factor.
- **`std::optional`**, `include/oblio/DirectSolver.h` and `ElmForestEngine.h` with their `.cpp`
  files, for amalgamation and threshold, so absent is a state rather than a sentinel value.
- **`std::is_same_v`**, the variable-template spelling, `src/SolveEngine.cpp`.
- **Pointer plus length** in place of `std::span`, which is C++20; one house convention, tree-wide.

What is left is small, scattered, and worth doing when the surrounding line is being touched anyway
rather than as a pass of its own. Each is listed with its site, so none of them has to be searched
for twice:

- **`[[nodiscard]]`**, used nowhere. Real value on the query accessors, `SparseMatrix::nnz`,
  `DirectSolver::nnz` and `relativeResidual`, `Vector::size`, where a discarded result is always a
  mistake. The largest of the four, and the only one that would catch a bug rather than tidy a
  line.
- **`std::is_same<Val, double>::value`**, `src/DirectSolver.cpp`, in `factor`. The one place still
  using the C++11 spelling where the rest of the tree says `_v`. Consistency only, and not an
  `if constexpr` candidate: the condition mixes a compile-time term with a runtime one.
- **`(void)ldd;`**, `src/BlasLapack.cpp`, the C++98 spelling of `[[maybe_unused]]`.
- **`std::string_view`** would suit `checkIndexRange`'s `const char* what`,
  `include/oblio/Types.h` and `src/Types.cpp`. The `const char*` parameters in
  `include/oblio/BlasLapack.h` must stay: that is the Fortran ABI, not a choice.

Structured bindings have no site here, nothing returning a pair or a tuple. Nested namespace
definitions buy nothing with a single namespace.

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
| UpdateMatrix | yes | checked | the multifrontal contribution block: one supernode's *entire* Schur complement, `u` by `u` (u = update size), symmetric, dense column-major plus its node indices. The per-supernode block of the multifrontal update stack, held in a flat `std::vector<UpdateMatrix>` sized once to the supernode count; `allocate`/`discard` are engine-only, `discard` frees rather than clears so the stack peak stays bounded. Now exercised by multifrontal for every factorization, static and dynamic (residuals at machine precision end to end). Replaces 0.9's abstract `UpdateStack` / concrete `UpdateStackDynamic` (an out-of-core split not needed in core) |
| NumFactorStatic | yes | checked | 0.9 `FactorsStatic`; SymFactor's structure copied, plus one flat value buffer with per-supernode offsets. Blocks are dense column-major rectangles (the upper front triangle is allocated and zero, so BLAS can take the whole block) |
| NumFactorDynamic | yes | checked | 0.9 `FactorsDynamic`. One index vector and one value vector per supernode, so a front can expand under delayed pivoting. Written by every factorization: the static ones run into it unchanged and produce a factor identical to the flat one, bit for bit; dynamic LDL writes it through the expansion and contraction verbs (`expandNodeIdx`, `resetVal`, `expandVal`, `swap`,
`contractVal`). No base class shared with the static one: `experiments/storage-options` showed a pointer array does the job a base would |
| NumFactorEngine | yes | checked | **Static factorizations functionally complete**: `Cholesky`, `StaticLDLT`, `StaticLDLH`, each left- and right-looking, real and complex. Cholesky checked against an independent dense Cholesky (4e-16); LDL by reconstruction, `L D L^H == P A P^T` (2e-15), through AMD ordering and supernodes. `StaticLDLH` (complex Hermitian LDL) is an **extension**: 0.9's complex LDL is symmetric only. Gaps, both in Owed: the LDL **perturbation branch has never fired**, and a complex input is **not validated as Hermitian**. **Dynamic LDL works for real input in both traversals**, `DynamicLDLT` and `DynamicLDLH` alike (the
same computation over the reals), with delayed columns crossing the forest and 2x2 pivots in the
solve, verified by residual in `test_pipeline` and by right-looking agreeing with left-looking bit
for bit. 0.9's two dynamic kernels are byte-identical between its left- and right-looking engines,
so only the driver differs: right-looking expands a front with `expandVal`, which preserves, where
left-looking uses `resetVal`, which discards. **`Traversal::Multifrontal`
is now complete: multifrontal works for every factorization, static and dynamic** (the third
traversal, a postorder pass carrying each supernode's contribution block up a
`std::vector<UpdateMatrix>` stack to its parent). Cholesky is verified against the dense-Cholesky
oracle directly (real and complex, natural and AMD) and the static LDLs by reconstruction across all
three traversals, all also end to end by residual. It reuses `factorStaticSupernode` and the new
`assembleUpdateMatrix` assembly, and forms the contribution block with `updateStaticUpdateMatrix`
(one `herk` for Cholesky, `formStaticUpper` + `gemmLower` for static LDL; the herk resolves to zherk
for complex, fixing 0.9's complex-Cholesky syrk bug for free). **Dynamic multifrontal, where delayed
columns meet the stack, is done.** Its driver is the left-looking dynamic skeleton (expand a front by
its children's delayed columns with `expandNodeIdx` + `resetVal`, assemble A past them, factor, which
may delay again) with the stack in place of the pull queue: assembling a child does both halves of the
assembly, its delayed columns become front columns (`assembleDelay`) and its contribution block
adds in (`assembleUpdateMatrix`), then it is contracted (`contractVal`) and freed. The factor reuses
the two dynamic pivot kernels unchanged, since the block layout at factor time matches
left/right-looking,
and the contribution block is formed by the new `updateDynamicUpdateMatrix` (`formDynamicUpper` +
`gemmLower`, block-diagonal D). Verified by residual at all three tiers, real and complex, symmetric
and Hermitian; tier 1 matches left-looking's delay and pivot counts exactly, heavier tiers match to
tolerance. **Complex `DynamicLDLT` needed no kernel change at all**: 0.9's complex
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

- **0.9's `rank_` was dropped in the port, and is now restored (2026-07-24).** 0.9's
  `factorDynamicLDL_` keeps a rank counter: `rank_ = a.getSize()` at the start, then `rank_--` at
  each `1 x 1` pivot accepted with an exactly zero diagonal, at three sites (the forced `1 x 1`
  when the candidate list empties, and the isolated-column accept `max1 == 0` in each of the two
  passes). The end value is the numerical rank, `n` minus the count of zero pivots. Real and
  complex 0.9 are identical here, three drops each. The port had transcribed the pass split but
  not the counter, so a singular matrix factored without any record that a rank was lost. Restored
  faithfully: `NumFactorDynamic::rank()` (read and write overloads, mirroring `numPerturbations`),
  initialized to the factor size in `initNumFactor`, decremented at the same three zero-guarded
  sites; the templated kernel covers both scalar types from one set of edits. This was found while
  checking whether the pass-1 zero-pivot comment matched 0.9, and it matters beyond bookkeeping:
  0.9 counting rank at exactly the forced-`1 x 1`-with-zero-diagonal site is the original author
  already treating that event as a rank loss, which is direct evidence for the pass-1 singularity
  conjecture in Section 7.8 of `sparse_factorization.md` and connects to the detection half of the
  full-BK-at-the-roots item in TODO. Not yet asserted in a test; that belongs with the pass-1
  verification work.

- **The input is not checked for Hermitian symmetry, and this is a correctness hole.** Moved to
  notes/TODO.md, under "Validate the input matrix", where it joins the two other unchecked
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

- **The two pivot bodies in the dynamic kernel were duplicated, and the debt is now paid.**
  0.9 splits `factorDynamicLDL_` on whether the supernode has update rows and writes both bodies
  out; the port followed it, deliberately, because those lines were 0.9's and merging them before
  they had ever run would have put their first execution in a shape the oracle never had.

  The trigger was the end-to-end residual, and once it fired the merge was done: the eliminations
  are now `factor1x1` and `factor2x2`, called from both selection loops, and the 2x2 block
  is read in one place, `readPivotBlock2x2`. 230 code lines became 204 across three functions, and
  the refactor was covered by 147 assertions throughout.

  Two things fell out of it that the duplication had hidden. **Pass 1 never reads the 2x2 block at
  all**, it accepts on `max1 == max2`, on the magnitudes alone, where pass 2 tests the determinant;
  the compiler said so with an unused-variable warning the moment the shared body stopped reading it
  for both. And `readPivotBlock2x2` turns out to be the single place the *symmetry* of D is decided,
  which is exactly what complex `LDL^H` needs to change.

  The selection loops stay separate, and that is a decision rather than unfinished work. See
  notes/TODO.md.

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
- **`sortForOptimalMultifrontal`** (0.9) / `rOptimizeForMultifrontal` (10.12) **is ported**, as
  `ElmForestEngine::sortForOptimalMultifrontal`, gated by `setOptimizeMultifrontal` and off by
  default as in both references. It runs after both compressions, exactly where 0.9 calls it, and
  reorders each supernode's children by decreasing `maxStorage(c) - updateSize(c)^2`, which is
  Liu's rule. 0.9 reaches that order by selecting the largest key onto the front of a list and then
  popping that list onto the front of the child list, two reversals that cancel; the port sorts
  directly, stably, to match 0.9's tie-breaking. Verified against a symbolic simulation: with the
  option on, the peak the tree model predicts drops to the computed optimum on every grid tested,
  16 to 38 percent below the unsorted order.
- **`labelDepthFirst`** (0.9 `labelDepthFirst_`) **is ported**, called immediately after the sort
  under the same option, as 0.9 calls it. It relabels the supernodes into a postorder, so every
  subtree holds a contiguous run of labels. This is the half that makes the sort pay: the drivers
  loop in increasing label order, so a contribution block is live from its own supernode to its
  parent *in the numbering*, which the child links alone do not affect. The traversal is 0.9's, a
  stack seeded with the roots pushed first to last and extended with each supernode's children the
  same way, so each set pops last to first, while the label counts down from `snodeSize - 1`; the
  two reversals leave the children of any supernode in increasing label order, preserving what the
  sort chose. It rewrites the map and permutes every per-supernode array, and it is contained
  entirely within the forest, since `SymFactorEngine` reads the forest through its accessors and
  derives everything fresh.
  With the pair on, the peak the multifrontal driver actually reaches falls to the computed
  optimum: 1845 to 1218 on grid2D 20x20, 42735 to 26425 on grid2D 80x80, 357398 to 255433 on
  grid3D 16x16x16, and unchanged on banded matrices, which have no sibling choice to make.
