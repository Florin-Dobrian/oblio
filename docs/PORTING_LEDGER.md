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

`OrderEngine -> ElmForestEngine -> SymFactEngine -> FactorEngine -> SolveEngine`,
orchestrated by `OblioEngine`.

Naming note: the modern tree renames as it ports. 0.9's `Matrix` is `SparseMatrix`
(a `DenseMatrix` will follow), and 0.9's `Symbolic` / `SymbolicEngine` pair is
`SymFact` / `SymFactEngine`. The table gives the modern name, with the 0.9 name in
the notes where they differ.

## Units

| Unit | Uses Val? | Status | Notes |
|---|---|---|---|
| Types | no | ported | sentinels (`NIL`), index typedefs |
| Permutation | no | checked | index-only. Ported: `set` (as `setOldToNew`/`setNewToOld`, replacing 0.9's direction flag) and `compose`. Not yet ported from 0.9: `get`, `read`/`write` (persist an ordering), `initialize2dGrid`/`initialize3dGrid` (structured orderings, useful as test inputs with hand-computable fill) |
| SparseMatrix | yes | checked | 0.9 `Matrix`; flat CSC, stored fully (both triangles). Build and structural symmetry tested |
| OrderEngine | no | checked | AMD (SuiteSparse) and MMD (Sparspak) vendored verbatim; output checked for validity as a permutation, not against 0.9's ordering |
| ElmForest | no | checked | data; supernodal shape, trivial supernodes for now |
| ElmForestEngine | no | checked | parent links, child/sibling links, roots, height, column sizes, fundamental compression, threshold amalgamation. Links and height recomputed independently; sizes and supernodes against the dense oracle, natural and AMD ordered. Amalgamation is greedy and not canonical, so only its tie-break-invariant properties are asserted |
| SymFact | no | checked | 0.9 `Symbolic`; flat index sets with per-supernode offsets |
| SymFactEngine | no | checked | 0.9 `SymbolicEngine`; index sets against the dense oracle, natural and AMD ordered. 10.12's design, 0.9's behavior (see DESIGN_DECISIONS) |
| Vector | yes | not started | |
| MultiplyEngine | yes | not started | |
| BlasLapack | yes | not started | traits + underscore handling |
| DenseMatrix | yes | not started | fronts and update blocks |
| FactorEngine | yes | not started | most complex; friend access into SparseMatrix/Factors/SymFact. Dynamic LDL pivoting is what forces the storage question (see DESIGN_DECISIONS) |
| SolveEngine | yes | not started | watch backward-solve index signedness |
| OblioEngine | yes | not started | top-level driver |

Units from 0.9 deliberately not carried over: `Utility` (`ResizeVector` and
friends, obviated by `std::vector`) and `Functional` (pre-C++11 comparators,
obviated by lambdas).

## Owed

- **Nothing is `verified` yet.** Every structural unit is `checked` against the
  dense oracle, which is real evidence, but none has been compared against 0.9
  buffer for buffer. That comparison is owed before numeric factorization starts
  consuming these structures, because a difference the oracle cannot see (index
  order within a supernode's block, say) is exactly the kind that surfaces as a
  numeric bug much later.
- **`sortForOptimalMultifrontal`** (0.9) / `rOptimizeForMultifrontal` (10.12) is not ported.
  Both references call it after compression to reorder children for the multifrontal stack;
  10.12's call site is commented out. Not needed until the numeric factorization exists.
