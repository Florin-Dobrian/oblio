# Decisions

Durable record of structural choices, newest first. Each entry: date, decision,
why. This is the file to open after a gap to reconstruct the project's shape.

## 2026-07-08 — Source extension: `.cpp` (headers stay `.h`)

Switch Oblio source files from `.cc` to `.cpp`; headers remain `.h`. All
extensions (`.cc`, `.cpp`, `.cxx`, `.C`) are identical to the compiler, so this is
convention, not correctness.

Why `.cpp`: cross-project consistency. The matching codebase uses `.cpp` (alongside
`.rs`, `.py`) — one extension per language across the ecosystem — so Oblio being
`.cc` made it the odd one out. `.cpp` is also the more common choice in the wider
world (`.cc` is mainly the Google-style corner). Done now because the tree is still
scaffold + PoC with no ported units yet, so the rename is at its cheapest.

Blast radius (all mechanical, filename-level, no semantics): `git mv *.cc → *.cpp`
across `src/`, `tests/`, `examples/`, `archive/`; the CMake source and executable
lists; the manual build glob in CLAUDE.md / README (`src/*.cpp`); `.clang-format`;
`archive/Makefile` and the build-command comments in the example files. Headers and
object files are untouched — `#include`s point at `.h`, and `Foo.o` derives from the
source basename regardless of extension, so no `.o` reference changes. This makes
the rename strictly safer than the `exp`→`ext` one, which touched `#include`s.

## 2026-07-08 — Explicit instantiation over header-only templates (rationale)

Decision (already active in CLAUDE.md; this entry records *why*, which otherwise
lives only in the archive template-instantiation example comments): Val-dependent classes keep a single
`Val` template, but their definitions live in `.cpp` files with explicit
instantiation for the supported scalar types, and headers carry declarations plus
`extern template`. Fuller treatments: `archive/oblio_modernization_notes.md` §"Why
explicit instantiation still works" (the Val-surface table) and
`archive/oblio-new-devlog.md` Session 3 (adoption + the link-failure proof). The
`_tpl` / `_ext` example files are the compact head-to-head (see naming below).

Mental model: a template is a recipe, not code — it generates code once a type is
plugged in. Header-only templates plug in *late and everywhere*: every translation
unit that includes the header re-runs the recipe for each type it uses, and the
linker discards the duplicates (N files × 2 types → the same bodies compiled 2N
times). That was the real cost of 0.9's header-heavy templating — not the templates,
but the repeated late instantiation. Explicit instantiation is "one template,
applied early, once per type": the recipe runs exactly twice, in one `.cpp`, at
library-build time, and every other file links the existing result instead of
re-running the recipe. (Two mechanisms achieve that "instead of re-running" —
declaration-only headers and/or `extern template` — see History below; they are
not the same feature and have different dates.) Generality is preserved (adding `float` is one
early application), but instantiation collapses from scattered-2N to centralized-2.

Key framing: build cost that scales with the number of scalar types is *incidental*
to the C++ compilation model, not *inherent* to supporting real and complex.
Nothing about "a matrix can hold real or complex" requires recompiling matrix code
in every includer — that only happens because header-only defers instantiation to
include time. Explicit instantiation removes the accident and keeps the capability.

The tradeoff, and why it's nearly free here: explicit instantiation gives up
instantiating *arbitrary* types at use sites — a consumer can't spin up
`Matrix<long double>` unless that line exists in the `.cpp`. For a maximally-generic
header library (Eigen) that's a real loss. For Oblio it isn't, because the scalar
world is closed and tiny: a type only makes sense if a dense BLAS/LAPACK kernel
exists for it, which bounds the space to BLAS's four — `float`, `double`,
`complex<float>`, `complex<double>` (s/d/c/z). We use two, might add the others.
Enumerating even all four is a handful of lines, far below the build cost of keeping
them implicit. So the one thing explicit instantiation costs is something a
closed-world numerical solver doesn't want anyway.

Bonus, and it matters for a verification-focused port: the `template class Foo<...>;`
lines *are* the list of supported types, in one place. Support becomes a declared,
reviewable fact rather than an emergent property of whatever anyone happened to
instantiate, and adding a type is a deliberate act that forces confronting whether
kernels and tests exist for it — exactly the gap the appendix flagged when complex
was "a new code path with zero test coverage."

Three mechanisms (not two — this frame makes the history click). Three distinct
tools, three natures, three dates:
- **Forcing** — `template class Foo<double>;` — emit all of `Foo<double>` in this
  TU. C++98.
- **Suppressing** — `extern template class Foo<double>;` — do *not* implicitly
  instantiate here; link it from elsewhere. C++11 (GCC extension earlier).
- **Definition-hiding** — not a keyword but a code-organization move: put member
  bodies in a `.cpp`, leave declarations in the header. A TU that can't *see* a body
  can't implicitly instantiate it. Works in every era.

Definition-hiding is the hinge, and it pairs with forcing. The three configurations:

| Case | Available | Approach | Build cost |
|---|---|---|---|
| 1 | neither forcing nor suppressing | Inclusion model: all definitions in headers, every TU re-instantiates what it uses, linker merges duplicates. No way to move bodies out and still get symbols. | High (~2N), unavoidable |
| 2 | forcing only (C++98) | Definition-hiding + forcing: bodies to `.cpp`, declaration-only headers, `template class Foo<double>;` in the `.cpp`. Other TUs see declarations only → can't implicitly instantiate → link the forced symbols. | Low |
| 3 | forcing + suppressing (C++11) | Case 2 still works and stays the choice; *additionally* you may keep bodies in headers and use `extern template` to suppress re-instantiation — needed only when definitions must stay header-visible. | Low |

Key insight: the big jump is **1 → 2, not 2 → 3**. Forcing is what unlocks the whole
technique (hide a definition, still guarantee the symbol). Suppressing is the
incremental step that only adds a second route for the case where you insist on
header-visible definitions. If you move bodies to the `.cpp` (Oblio does), you never
need it.

Precondition for all of it: an **enumerable** type set. For genuinely arbitrary
types you can't force or suppress anything (you don't know the list), so you're in
Case 1 regardless of language version. Forcing/suppressing are tools for closed type
sets — which Oblio's is (BLAS s/d/c/z), the same fact that makes the tradeoff above
nearly free.

Dates and 0.9: forcing is C++98, suppressing is C++11, so when 0.9 was written (late
90s) only Case 1 was portably reachable — header-only, inclusion model. That was the
*correct* choice for the era, and not because suppressing was missing (Case 2 doesn't
need it) but because template separate-compilation was a portability minefield then
(the `export` saga, inconsistent two-phase lookup, compilers disagreeing on
inclusion-vs-separation). Header-only-everything was the safe default. The modern
refactor applies matured portability; it does not correct a 0.9 error.

Where Oblio sits: current `ext` code is **Case 3** — bodies in `.cpp` (declaration-only
headers) *plus* `extern template`. But because the headers are already
declaration-only, the build win is really Case 2's (definition-hiding + forcing); the
`extern template` lines suppress nothing here (no visible header body to instantiate),
so they are documentation, not mechanism — a header annotation of intent, latent
unless a body is later (wrongly) added to a header. See the naming note below for the
full plain-vs-guarded framing. So Oblio's pattern was achievable in C++98; C++11 was
not strictly required.

Template-instantiation example — naming (one algorithm — dense mat-vec — built
three ways). File names are `<Class><Variant>` with no separator — e.g.
`MatrixImplicit.h`, `MatrixPlainExplicit.{h,cpp}`, `MatrixGuardedExplicit.{h,cpp}`
(same for `Vector`, `MultiplyEngine`). The three variants:
- `Implicit` — body in the header; instantiated implicitly per TU. Stands in for
  what 0.9 effectively was.
- `PlainExplicit` — bodies in the `.cpp`, header signatures only (explicit
  instantiation, forcing only).
- `GuardedExplicit` — plain explicit + `extern template` in the header. The
  pattern used in the real tree.

Conceptual framing (two axes):

- **Axis 1 — where the body lives (implicit vs explicit).** *implicit* = body in
  the header (`.h`), instantiated implicitly per translation unit → `Implicit`.
  *explicit* = body outside the header, forced in the `.cpp`; header carries
  signatures only → `PlainExplicit` and `GuardedExplicit`.
- **Axis 2 — applies only within explicit; guarded vs plain.** *plain explicit* =
  `.cpp` bodies, declaration-only header, nothing more → `PlainExplicit`.
  *guarded explicit* = same, plus `extern template` in the header → `GuardedExplicit`.

  The three named layers:
  - **implicit** — body in `.h` (`Implicit`)
  - **plain explicit** — bodies in `.cpp`, header signatures only (`PlainExplicit`)
  - **guarded explicit** — same as plain explicit, plus the guard (`GuardedExplicit`)

  Important: the guard is NOT a sub-kind of "more correct" explicit — plain explicit
  is a complete, valid design. `extern template` only ever acts on *visible* header
  bodies (the implicit failure mode); in a declaration-only design there is nothing
  for it to suppress, so it is pure documentation — a header annotation reading
  "instantiated elsewhere," aimed back at the implicit branch it guards against. It
  gains mechanical effect only if someone reintroduces a header body (which the
  "definitions live in `.cpp`" invariant in CLAUDE.md forbids). So "guarded" (a
  reminder/guard), not "suppressed" or "enforced" — in this design it suppresses
  nothing and enforces nothing; the invariant does the enforcing.

All three are built and tested together via the example's `Makefile` (`make test`)
against one shared source, `test_multiply.cpp` — they must produce identical
results, and the plain/guarded variants share the same link-failure behaviour when
their `.cpp` files are omitted (empirical confirmation that with declaration-only
headers `extern template` suppresses nothing — it is documentation, not mechanism).
Selector macros: `OBLIO_TI_IMPLICIT` / `OBLIO_TI_PLAIN_EXPLICIT` /
`OBLIO_TI_GUARDED_EXPLICIT`.

Naming history: suffixes were once `_tpl`/`_exp`/`_ext`, where `_exp` inaccurately
labeled the extern-template variant. Renamed in two steps — first `_exp`→`_ext` to
free `_exp` for the genuine forcing-only variant, then all three to the conceptual
`Implicit`/`PlainExplicit`/`GuardedExplicit` once the two-axis framing settled.

## 2026-07-08 — Matrix naming: explicit `SparseMatrix` / `DenseMatrix`

Rename the sparse matrix type from `Matrix` to **`SparseMatrix`**; keep
**`DenseMatrix`**. Both are plain concrete types with **no shared base class**.

Why: the old `Matrix` (implicitly sparse) + `DenseMatrix` (explicitly marked) is a
half-committed convention — it privileges sparse as the unmarked default but marks
dense as the exception, so a reader must already know "unmarked = sparse" to parse
it. Going fully explicit removes that. It also makes matrix naming consistent with
the graph naming used elsewhere (`GeneralGraph` / `BipartiteGraph`, no bare
`Graph`). Note the typed-library mainstream (Eigen, scipy) actually defaults the
*other* way — unmarked = dense, `SparseMatrix` marked — so `Matrix` = sparse
actively misleads anyone arriving from there. The sparse-first precedent
(CSparse/SuiteSparse `cs`) is the one respectable counter-argument, but explicit
wins on cross-codebase consistency and reader-independence.

No speculative base `Matrix` interface. If something ever needs to consume sparse
and dense polymorphically, add the interface then, on top of the two concrete
types — not before a caller forces it.

Scope: this is a rename (wrapping/API), not an algorithm change — it's on the
port-and-modernize track, not the rewrite track. Do it as one deliberate
mechanical pass **before** porting proper, since it's cross-cutting (every `friend`
decl, FactorEngine, SolveEngine, tests name `Matrix`) and only gets more expensive
as code solidifies around the name. **Oracle mapping: 0.9 `Matrix` ↔ modern
`SparseMatrix`** — record this so output comparisons against 0.9 stay unambiguous.

## 2026-07-08 — Minimal abstraction; containers are the structure that matters

Design stance for the port: concrete types, minimal OO ceremony, `std::vector` as
the spine. Don't build base classes, single-implementation interfaces, or
inheritance ahead of an actual polymorphic caller (YAGNI). A direct solver is a
pipeline of concrete transforms, not a class hierarchy.

One guard against a tempting over-correction: "AI makes code cheap to generate, so
structure matters less" is only half true. Generation got cheap; **verification did
not** — and this project is the proof (every PoC bug was cheap to write, expensive
to trust, and fixed only by slow comparison against the 0.9 oracle). Structure's
real job was never to save typing; it's to keep code readable and checkable and to
localize where a bug can hide. That matters *more* when code is machine-generated,
because the bottleneck shifts onto review. So: drop OO flavor that only served
human writing time; keep the structure that serves verification (clear module
boundaries, one concern per unit, the friend seams that let a supernode block be
diffed against 0.9).

This is why "proper containers everywhere" (`std::vector` over `Array`) is not a
style preference but the load-bearing invariant: a vector carries its own size (no
drifting length variables), self-frees (no leak/double-free surface), bounds-checks
under sanitizers, and hands clean pointers to BLAS via `.data()`. Each removes a
bug class you'd otherwise verify the absence of by hand. For this codebase, the
container discipline *is* the architecture that matters.

## 2026-07-07 — Align with standard project files; adopt clang tooling

The doc set maps to established conventions rather than being bespoke:
CODING_RULES.md ≈ a style/conventions guide, DESIGN_DECISIONS.md ≈ a lightweight ADR log
(single-file variant of the Nygard/ADR pattern), CLAUDE.md ≈ the agent-instructions
file (AGENTS.md is the emerging cross-tool equivalent — if we ever run more than
Claude Code, keep content in AGENTS.md and make CLAUDE.md a one-line `@AGENTS.md`).
PORTING_LEDGER.md stays bespoke — it's specific to a port, no standard analog.

Renamed CPP_RULES.md → **CODING_RULES.md** and made it language-general (Rust is a
likely future scope), with per-language sections.

Adopted **`.clang-format`** and **`.clang-tidy`** so mechanical rules (`nullptr`,
`using`, `enum class`, `= delete`, narrowing) are tool-enforced instead of written
as prose an agent may skip. `.clang-tidy`'s `modernize-*` family does part of the
port mechanically. The prose rules doc now holds only judgment calls tools can't
check.

Stubbed ahead of need (deliberate, not yet load-bearing): CONTRIBUTING.md and
CHANGELOG.md ("Keep a Changelog" format). CONTRIBUTING fills out on going public
(and becomes the canonical build/test source then); CHANGELOG gets its first real
entry at the first tagged release — which requires settling the tree's version
identity (Oblio 11 vs a fresh 1.0), tied to the working-tree question below. Still
outstanding: LICENSE (bundled AMD is BSD-3-clause, so a project license decision is
eventually needed).

Tension worth remembering: canonical ADR is many small numbered files, which suits
humans browsing history but fights the "always in context" goal — a growing pile
can't all stay loaded. The single-file log + distilled always-on layer is a
deliberate adaptation of the standard, not an oversight.

## 2026-07-07 — What stays always in context vs on demand

Goal: coding rules and active design constraints present in every session, without
letting context grow unbounded.

Mechanism: only CLAUDE.md and its `@`-imports load every session; everything else
is read on demand. `@import` costs the same context as inlining — it just keeps
the file separate and editable.

Split:
- **Always-on** (in CLAUDE.md or imported): invariants (inline), conventions
  (`@docs/CODING_RULES.md`), and a distilled *Active design constraints* summary.
- **On demand**: the full DECISIONS log (this file — grows over time, so importing
  it would erode context), PORTING_LEDGER (porting-specific; read after a gap),
  README, archive history.

Rule of thumb: the always-on set is the distilled essence; narrative, dates,
history, and open questions stay here. Keep CLAUDE.md well under ~200 lines — past
that, adherence drops. If the always-on set grows, distill harder rather than
importing the growing logs.

## 2026-07-07 — Invariants live in CLAUDE.md, not CODING_RULES.md

The C++ **invariants** (port-verbatim, `std::vector` default, no signed→unsigned
index slips, `.data()` to BLAS) are written directly in CLAUDE.md. Only the
**conventions** (style preferences) stay in CODING_RULES.md.

Why: Claude Code auto-loads CLAUDE.md every session (directory walk up to repo
root, concatenated) but does **not** auto-load other files — a pointer from
CLAUDE.md to CODING_RULES.md does not pull its contents in. So a rule meant to be
active every session must physically live in CLAUDE.md, or it isn't loaded until
Claude happens to read the file it's in. Invariants are always-on; conventions
are fine to read on demand.

Do not "consolidate" the invariants back into CODING_RULES.md — that silently
disables them.

## 2026-07-07 — Documentation structure

Established this doc set: CLAUDE.md (operating contract + index),
PORTING_LEDGER.md (per-unit status), CODING_RULES.md (conventions a linter can't
enforce), DESIGN_DECISIONS.md (this log), plus `.clang-format` / `.clang-tidy` for
tool-enforced mechanical style. The existing md files — devlog, modernization
notes, appendix, README — are kept as history/rationale for now; decide later
whether to fence or retire them.

Why: after a multi-month gap the project needs to be reconstructable from three
things — history, layout, decisions — without re-deriving them. CLAUDE.md holds
the operational index (it's what Claude Code loads each session); the reasoning
lives here so CLAUDE.md stays lean. The two live in different files on purpose and
should not duplicate each other.

(This entry documents a list that includes DESIGN_DECISIONS.md itself. The self-reference
is intentional, not an oversight.)

## 2026-07-07 — OPEN: which tree is the working copy

Unresolved; resolve before the first port step and record the outcome here.

Candidates:
- **(a) Fresh port from 0.9** into a new tree, treating the PoC the way 10.12 is
  treated — learn from, don't build on. "Replacing Array gradually" only applies
  to a tree that still has `Array` (0.9 or 10.12), which points here.
- **(b) Continue the PoC tree**, where Array→vector is already done — in which
  case the remaining work is different (coverage, finishing, cleanup), and this
  is not really an Array migration.

The stated plan ("port carefully from 0.9, one function at a time, replace Array
gradually") reads as (a). Confirm.

## 2026-03-07 (PoC) — Choices carried from the proof of concept

Recorded for continuity. The PoC was exploratory, so revisit each on its merits
rather than inheriting it unquestioned.

- **One `Val` template** instead of 0.9's separate `*Real.h` / `*Complex.h` file
  pairs. Cholesky treated as Hermitian input; LDL^T as complex-symmetric.
- **`std::vector` storage** instead of the hand-rolled `Array`.
- **Explicit template instantiation** for `double` and `std::complex<double>` —
  headers declare, `.cpp` files define and instantiate. Faster builds; `float` /
  `long double` remain one line away.
- **Namespaced `include/oblio/` headers**, declarations only.
- **Flat `src/`** — all sources, including `Mmd.cpp` / `Amd.cpp`, directly in `src/`.
