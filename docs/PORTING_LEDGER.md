# Porting Ledger

Tracks the migration of each 0.9 unit into the modern tree. After a context gap,
read this first — it turns "refresh my context" into a two-minute scan.

## How to use

- One row per unit (a class or a small function-group).
- Status: `not started` → `ported` (translated, compiles) → `verified` (output
  checked against 0.9 on real input).
- Notes: Array→vector boundary, index-signedness fixes, open questions — anything
  a future session needs.
- Advance to `verified` only after comparing output against 0.9.

## Needs reconciliation

The inventory below is inferred from the three-month-old notes, **not** from a
fresh read of the 0.9 tree. Reconcile against the actual 0.9 source before
relying on it — names and groupings may differ.

## Pipeline

`OrderEngine → ElmForestEngine → SymbolicEngine → FactorEngine → SolveEngine`,
orchestrated by `OblioEngine`.

## Units

| Unit | Uses Val? | Status | Notes |
|---|---|---|---|
| Types | no | not started | typedefs, enums, sentinels |
| Utility | no | not started | helpers; `ResizeVector` etc. — removal candidates |
| Functional | no | not started | pre-C++11 comparators; likely deletable |
| Permutation | no | not started | index-only |
| Matrix | yes | not started | largest unit; storage + build + I/O |
| Vector | yes | not started | |
| MultiplyEngine | yes | not started | |
| OrderEngine | no | not started | index-only |
| ElmForestEngine | no | not started | elimination forest; index-only |
| SymbolicEngine | no | not started | etree, col counts, index sets; index-only |
| Symbolic | no | not started | data struct |
| BlasLapack | yes | not started | traits + underscore handling |
| FactorEngine | yes | not started | most complex; friend access into Matrix/Factors/Symbolic |
| SolveEngine | yes | not started | watch backward-solve index signedness |
| OblioEngine | yes | not started | top-level driver |

`Array → std::vector` is the cross-cutting migration, not a unit — it lands unit
by unit as each row is ported, and the boundary between converted and
un-converted code creeps up the call graph. Track that boundary in the Notes
column as it moves.

## Suggested start

Non-numeric, leaf-first: Types, Utility, Functional, Permutation. These carry no
`Val` and the fewest dependencies, so the Array→vector boundary starts small and
grows outward from stable ground.
