# CLAUDE.md

Operating contract for the oblio refactoring effort. Auto-loaded by Claude Code
each session — keep it lean. Rationale lives in docs/DESIGN_DECISIONS.md, not here.

## What this is

oblio is a sparse direct solver (supernodal Cholesky / LDL^T). This effort is a
careful refactor of **oblio 0.9** into modern C++17: migrate the hand-rolled
`Array` container to `std::vector`, modernize idioms, one function at a time.
Porting and modernizing are the same act — a good algorithm re-expressed in
modern style. They are not two passes.

## Source of truth

- **oblio 0.9** — complete, correct, fully functional. The reference and the
  oracle. Every ported unit is verified against 0.9's output before it counts as
  done.
- **oblio 10.12** — a partial, unfinished refactor of 0.9. Reference only; do
  not build on it.
- **PoC tree** (built in one day, 2026-03-07) — a proof of concept. Learn from it;
  do not assume its code is correct.

## Process (every session)

1. **One unit at a time.** One function or small function-group per step. Verify,
   record in the ledger, then advance.
2. **Verify against 0.9.** Same inputs through both codebases; compare outputs.
   Don't advance a unit that hasn't been checked.

## Invariants (breaking one is a bug)

These live here, not behind a link, because only this file loads every session.
Conventions (style preferences) are imported below from docs/CODING_RULES.md.

- **Port, don't rewrite.** Carry 0.9's algorithm over unchanged. Changing *what*
  is computed is a rewrite — a separate, opt-in track, not part of a port. Every
  PoC bug came from reimplementing a 0.9 algorithm instead of porting it.
- **`std::vector` is the default container.** Exceptions, each deliberate:
  fixed-size small blocks (e.g. 2×2 pivots) → `std::array`; non-owning view over a
  column/block → pointer + length (C++17, no `std::span`), one house convention;
  a `std::vector<std::vector<T>>` is a valid port target but **do not flatten** it
  into one buffer with offsets mid-port — that's a layout change (rewrite track).
- **No signed→unsigned index slips.** An index or offset that was `int` and can go
  negative must not become `size_t` arithmetic. Unsigned underflow on a descending
  loop is a bug this project has already hit once. Danger spots: descending loops,
  pointer offsets.
- **Contiguous storage to BLAS via `.data()`.** Never reintroduce the
  `(size==0) ? NULL : &v[0]` guard; `.data()` is well-defined for empty vectors.

## Active design constraints

Decisions that shape code written now. Full rationale and history in docs/DESIGN_DECISIONS.md;
these are the always-on summary.

- **One `Val` template** for scalar type — `double` and `std::complex<double>`.
  Cholesky assumes Hermitian input; LDL^T variants assume complex-symmetric. No
  separate `*Real.h` / `*Complex.h` file pairs.
- **Explicit instantiation.** Headers declare; `.cpp` files define and instantiate
  for `double` and `std::complex<double>`. Adding a scalar type is one line per `.cpp`.
- **Namespaced headers** under `include/oblio/`, declarations only.
- **Flat `src/`** — all sources directly in `src/`, no per-category subdirectories.

## Coding rules (imported)

Loaded every session via import so conventions stay consistent. Edit the file, not
this line. (Path is relative to this file; CODING_RULES.md now lives in docs/, so
the import is `@docs/CODING_RULES.md` — a wrong path imports nothing, silently.)

@docs/CODING_RULES.md

## Build

macOS (alpamayo, Apple Silicon; Accelerate provides BLAS/LAPACK):

```
g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -I include \
    tests/<test>.cpp src/*.cpp -framework Accelerate -o <test>
```

Linux: replace `-framework Accelerate` with `-lblas -llapack -lm`.

## Tooling

- `.clang-format` — formatting. Run routinely; safe to apply whole-file.
- `.clang-tidy` — `modernize-*` and narrowing checks. Run as a *catch* per unit
  (report, no `--fix`) or as a deliberate sweep — not blind `--fix` mid-port, since
  it rewrites whole files and blurs one-unit-at-a-time diffs.
- Both are starting points, not yet calibrated to the tree; see the headers in each
  config before mass-applying.

## Docs

- **CLAUDE.md** (this file) — operating contract + doc index.
- **docs/PORTING_LEDGER.md** — per-unit porting status. Read first after a context gap.
- **docs/CODING_RULES.md** — conventions a linter can't enforce. Imported above, so
  always loaded. Language-general.
- **docs/DESIGN_DECISIONS.md** — full rationale, history, dates, open questions. Read on
  demand; the code-shaping subset is summarized above under Active design
  constraints.
- **`.clang-format` / `.clang-tidy`** — mechanical style + modernization,
  tool-enforced (not context). See Tooling above.
- **CONTRIBUTING.md / CHANGELOG.md** — stubs, ahead of need. Fill CONTRIBUTING on
  going public; start CHANGELOG at the first tagged release.
- **README.md** — public-facing overview.
- **archive/** — history: `oblio-new-devlog.md`, `oblio_modernization_notes.md`,
  `oblio_modernization_appendix.md`. These describe 0.9/10.12 and the PoC;
  rationale, not current guidance.
