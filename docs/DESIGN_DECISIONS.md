# Decisions

Durable record of structural choices, newest first. Each entry: date, decision,
why. This is the file to open after a gap to reconstruct the project's shape.

## 2026-07-09 — Index types: `std::int32_t` IDs (NIL = -1), `std::size_t` offsets

Two kinds of integer, two types:

- **IDs** — a value that *names* a vertex/row/column/supernode, and may need a "none"
  marker — are **`std::int32_t`**, with sentinel **`NIL = -1`** (`constexpr std::int32_t`).
  E.g. `SparseMatrix::rowIdx`, the permutation maps, and the forest's parent / child /
  sibling / supernode-map arrays.
- **Offsets / counts / sizes** — row-pointers, `nnz`, dimensions, anything that indexes
  into or measures — are **`std::size_t`**. E.g. `SparseMatrix::colPtr`, `numCols()`.

**Why signed int32 for IDs.** The forcing function was the sentinel. The forest needs a
"no parent"/"no child" marker; on an unsigned `std::size_t` that can only be the max
value, spelled `static_cast<std::size_t>(-1)` — defined behaviour but an ugly wraparound
smell, and easy to misuse in arithmetic. A signed `std::int32_t` gives a clean, obvious
`-1`. This also matches (a) the graph code, which already uses `int32_t` vertex IDs with
`static const int NIL = -1;` so companion arrays like `mate` hold `-1` naturally, and (b)
the vendored AMD/MMD, which are `int`-based — so `rowIdx` as `int32_t` largely removes the
`size_t→int` conversion at that boundary. The cost is a ~2.1-billion index cap (int32 vs
size_t's range). Accepted deliberately: cleaner and more agile, and Oblio isn't targeting
matrices past 2^31 structural indices for now. If that changes, widen the ID type to
`std::int64_t` in one place.

**Why `size_t` stays for offsets.** Offsets are never negative and can legitimately
exceed 2^31 even when the ID *count* doesn't (a row-pointer indexes into `nnz`-length
arrays). So they keep the full unsigned range and never carry a sentinel. This mirrors the
graph's `GeneralGraph` exactly: `idx` (row-pointer) is `size_t`, `adj` (neighbour IDs) is
`int32_t`.

**No signed/unsigned friction, because loop counters stay `std::size_t`.** g2_csr
demonstrates the discipline: counters that enumerate positions are `size_t` (so they
compare against `.size()` cleanly), and `int32_t` appears only as *stored values* that
may be `NIL` — never as a loop variable. `size_t` is a safe superset of the non-negative
`int32_t` range, so viewing an index as `size_t` in a bounded, non-negative loop loses
nothing. The two types reconcile with explicit casts at exactly two crossings:
`static_cast<std::int32_t>(counter)` when storing a counter into an ID array (the
narrowing where the 2^31 cap lives), and `static_cast<std::size_t>(id)` when an ID
subscripts an array (widening; guard against `NIL` first if it could be a sentinel). Casts
are few and mark the ID↔offset boundary.

Spelling: `std::int32_t` from `<cstdint>`, `std::size_t` from `<cstddef>` — std-qualified,
C++ headers, same rule as every other stdlib type. Applied first to `SparseMatrix`
(`rowIdx`→`int32_t`, `colPtr`→`size_t`); `Permutation` (maps→`int32_t`) and `ElmForest`
(parent/etc.→`int32_t` with a shared `NIL`) follow. The graph code uses bare `int32_t`;
it's the code to bring in line later, like the matching codebase was for `size_t`.

## 2026-07-09 — Friend grants write access; reads are public

The engine↔data access rule:

- **`friend` = write access.** An engine befriends exactly the data class(es) it
  *produces/mutates*. There is no public mutation API for those internals — only the
  producing engine reaches the private members. E.g. `OrderEngine` writes
  `Permutation`; `ElmForestEngine` writes `ElmForest`. Friendship is declared by the
  data class (`friend class FooEngine;`), not by the engine.
- **Reads are public — everywhere, including hot paths.** Engines, tests, and users
  all read via the public const API. `friend` is *not* needed for reading.
- **This corrects an earlier overstatement** (the friend/BLAS entry below implied
  friend was needed for hot-path reads). The `experiments/friend-access/` study
  measured a *per-element, cross-translation-unit accessor call* — which can't inline,
  so it blocks vectorization — against direct friend member access; friend won ~6×.
  But that gap comes from calling the accessor *per element*, not from using a public
  accessor at all. Two public-accessor patterns match friend's performance exactly:
  (a) hand BLAS the block via a single `.data()` call — BLAS then owns the O(n³) loop,
  so the one call is free; (b) for a hand-written element loop, bind the returned
  container once (`const auto& v = A.val();`) and loop over *that* — one non-inlined
  call total, then a vectorizable loop over contiguous memory. So the hot-path
  discipline is **"bind-once / pass the pointer," not "be a friend."**
- **Consequence for `SparseMatrix`:** A is input and nothing writes it (its
  construction path is TBD), so it needs *no* friends. Its earlier
  `friend class OrderEngine` is removed; `OrderEngine` reads A via
  `colPtr()`/`rowIdx()`/`numCols()` and remains a friend only of `Permutation` (which
  it writes). `ElmForestEngine` already followed this (friend of `ElmForest` only).
- **Numerical data classes** (`Factors`, …) still befriend their producing engines for
  *writes*; their hot-path reads also go through public accessors with the
  bind-once/pointer discipline.
- **Exposure stance (pragmatic, not purist):** exposing internals read-only is fine;
  we don't design curated "won't-break" read APIs up front. A representation change
  already forces editing the friends (writers); public read exposure adds only the
  *tests* to that blast radius — cheap, and tests *should* feel a rep change. Curate a
  narrower public API later only if a structure's representation proves unstable. For
  the canonical structures we have (CSC, etree `parent[]`), the representation is the
  settled standard form, so exposing it by reference is low-risk.

## 2026-07-09 — Sparse matrix storage: flat CSC, stored FULLY (both triangles)

`Matrix` (input A) stores its structure and values as flat **compressed sparse column
(CSC)**: three contiguous `std::vector`s — `colPtr` (size+1), `rowIdx` (nnz), `val`
(nnz) — row indices sorted ascending per column. A symmetric matrix is stored **fully
(both triangles)**, each column holding its complete neighbour list plus the diagonal.

Storage layout (flat CSC vs vector-of-vectors) and triangle (full vs lower) are two
separate decisions:

**Layout — flat CSC**, across the three generations:
- **0.9** stored CSC via four manually `new`/`delete`d `Array*` pointers — the manual
  memory management the port removes.
- **10.12** modernized to `std::vector` but chose **vector-of-vectors** (one inner
  vector per column) — RAII but the wrong layout: columns scatter across the heap.
- **PoC / port** use flat CSC (`mColPtr`/`mRowIdx`/`mVal`) — vectors *and* contiguous,
  satisfying the "contiguous storage to BLAS via `.data()`" invariant.

**Triangle — full, not lower.** Both 0.9 and 10.12 store A **fully**: 0.9's
`getNumberOfNonzeroEntries() = size + 2*numOffDiagonals` (each off-diagonal in both
triangles) and its "storing A within the structure of A+Aᵀ"; 10.12 has `SymmetrizeStrc`
and its etree reads full per-column neighbour lists. The **PoC diverged**, storing the
lower triangle only ("stored as lower triangle in CSC"). That divergence is what forced
every structural consumer to expand lower→full first (the MMD path in `OrderEngine`, the
etree in `ElmForestEngine`) — and the etree bug where lower-triangle input silently
produced an empty tree until expansion was added. **The port matches the oracle: A is
stored fully.** Consequences: structural phases (ordering, elimination forest, symbolic)
read each column's neighbours directly with no expansion (the etree's diagonal
self-skips via `lc1 < lc2`; MMD just strips the diagonal; AMD ignores it); it's the
faithful port (lower-triangle was a rewrite of the data structure); and it's the natural
substrate for a future **unsymmetric extension** — factor the symmetrized structure
A+Aᵀ while carrying asymmetric values. Cost: ~2× off-diagonal storage for A, and the
numeric phase carries a redundant triangle. Accepted, matching 0.9/10.12; A is the input
and is far smaller than the factors, where the real memory lives.

Open for the port: the PoC exposed this as a public `struct` with a `fromCOO` builder
and a weak `isValid`. The modern `SparseMatrix` keeps the flat-CSC layout but is a
`class` with `friend` engines and a structural interface; 0.9 is the oracle for the
COO→CSC assembly details (zero-diagonal insertion, duplicate merging, symmetrization),
10.12 shows which operations the solver actually calls.

## 2026-07-09 — Two layers of modernization: rules prevent, clang-tidy catches

The coding rules and `.clang-tidy` are complementary layers catching different
failures at different times — not redundant work:
- **Coding rules** (CODING_RULES + CLAUDE invariants) — preventive and broad. They
  shape code as it's written and cover what no tool can judge: port-verbatim
  discipline, friend→BLAS, when to split a header, `.cpp`, `mFoo`, `std::size_t`.
- **`.clang-tidy` `modernize-*`** — a mechanical safety net, narrow but certain. It
  can't reason about intent, but within scope it catches the idiom slips that get
  through (a stray `NULL`, a `typedef`) and can auto-fix them.

So a `NULL` written by mistake is exactly what the rules *say* and the tool
*guarantees* — belt and braces, each doing what the other can't.

`modernize-*` is aligned with this project's purpose at the concept level: porting 0.9
(late-90s C++) forward *is* turning `typedef`→`using`, `NULL`→`nullptr`, raw loops→
range-for, and so on across old code — exactly what that checkset does. So it can do
part of the mechanical modernization *for* you on each ported file (`--fix`), leaving
your attention on the algorithmic faithfulness no tool can verify. Hence the per-unit
workflow (in CLAUDE.md Process): port faithfully → `clang-tidy --fix` → verify vs 0.9.

## 2026-07-08 — `experiments/` convention (runnable design studies)

`experiments/<name>/` holds self-contained, runnable studies that establish or
validate a coding standard before it is applied in the main tree. Each is its own
folder with its own sources, `Makefile`, and `README.md`; builds standalone
(`make test`); and is reference/teaching material — **not** part of the main Oblio
build. Executables carry the `_cpp` suffix and are gitignored.

Distinct from its two neighbors:
- `archive/` — frozen history (superseded PoC devlog, 0.9-analysis notes, old
  harnesses). Not maintained, not built.
- `examples/` — usage samples showing how to *call* the library (`examples/basic.cpp`).

An experiment answers a design question with code you can run and measure, then feeds
a decision here. Current studies: `template-instantiation/` (how to instantiate the
`Val` template — implicit vs plain/guarded explicit) and `friend-access/` (public API
vs `friend`-direct access, with timing). Experiments use the already-settled standards
(guarded explicit, `.cpp`, `mFoo`, `Oblio` namespace), so they double as worked
references for those standards.

## 2026-07-08 — Numerical hot path: `friend` access, then BLAS (carried from 0.9)

The 0.9 design for numerical work, which the port preserves: **an engine reaches the
data's contiguous block via `friend`, then hands that raw block to BLAS** wherever
BLAS applies (gemv/gemm/syrk/trsm/potrf, via Accelerate on macOS). Not a new choice —
this two-step (`friend` → BLAS) is how 0.9 does dense numerics.

This *is* supernodal numerical factorization: supernodes are dense blocks embedded in
the sparse structure, and the numerical phase is a long sequence of dense BLAS calls on
them — `syrk`/`gemm` for Schur-complement updates, `potrf`/`getrf` to factor the pivot
block, `trsm` for the off-diagonal solve — repeated per supernode across the whole
elimination. `FactorEngine` reaches each supernode's contiguous storage via `friend`
and passes the pointer straight to BLAS, no copy, thousands of times per factorization.
So `friend` isn't an optimization detail; it's the access mechanism the entire numerical
phase is built on. The `experiments/friend-access/` mat-vec is the single-block toy of
this pattern.

Important distinction — **the supernode blocks live in the factors, not in A.** `A`
(the input `Matrix`) is never handed to BLAS block-by-block; it's *read* — its structure
by the ordering and symbolic phases, its values once when they're scattered into the
factor. The dense blocks that BLAS operates on are created during factorization and
stored in `Factors`. So the `friend`→BLAS hot path is specifically a `Factors` /
`FactorEngine` story. `A`'s storage (CSC) is chosen for a different reason: cheap,
cache-friendly *sequential column traversal* by the structural phases, plus being the
standard interchange format (what AMD/MMD expect). Both `A` and the factors favor flat,
contiguous storage over 10.12's vector-of-vectors — but for `A` the reason is streaming
structural reads, and for the factors it's the contiguous block BLAS needs.

Data classes (`Matrix`, `Vector`, `Factors`, `Symbolic`) expose a public,
bounds-checked API (`operator()`, `operator[]`) for reads — by all callers, on hot
paths too (see the "friend = write access" entry above: reads are public; the hot-path
discipline is bind-once / pass `.data()` to BLAS, not friendship). Engines befriend a
data class only to *write* it. 0.9 grants `FactorEngine` friend access into `Matrix`,
`Factors`, `Symbolic`; in the port that friendship is retained only where the engine
writes (the factor storage), and reads go through the public API.

Why (performance): the public-operator path is one non-inlined, cross-translation-unit
call per element (data-class body in its `.cpp`, loop in the engine's), which blocks
vectorization. Direct `friend` access fetches the raw block pointer once and walks
contiguous memory, which vectorizes — measured ~6× on Apple Silicon (M4/AppleClang),
~3× on x86/g++, over the API path. But the raw block should then go to **BLAS**, not a
hand loop: on the M4, Accelerate's `dgemv` ran the 2000×2000 mat-vec ~6× faster than
the vectorized hand loop (and ~36× over the API path), because it breaks past a single
core's bandwidth ceiling (multithreading + prefetch/blocking) that a hand loop can't.
The advantage only grows for the compute-bound O(n³) kernels (`dgemm`/`syrk`/`trsm`/
`dpotrf`) the factorization leans on. So the hand loop is a fallback/baseline; where
BLAS applies, use it.

Why `friend`, not public getters: `friend` is *tighter* encapsulation, not looser — it
grants access to exactly the named engine classes, where a public `data()`/getter
exposes internals to the whole program. It honestly encodes that a data class and its
engines are one subsystem split for organization, not two modules talking through a
narrow API. Deliberate pragmatic choice over OO-purist accessors — and it's faster.

Consequence for porting: `friend` couples a data class to its engines, so the natural
port/verify unit is the *cluster* (e.g. `Factors` + `FactorEngine`), not the data class
in isolation — the friend boundary sets porting granularity.

Measured note: the gap is structural (non-inlined calls / vectorization), not assertion
overhead — toggling bounds-check asserts barely moves it. See `experiments/friend-access/`.

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
across `src/`, `tests/`, `examples/`, `archive/`, `experiments/`; the CMake source and
executable lists; the manual build glob in CLAUDE.md / README (`src/*.cpp`); `.clang-format`;
the example `Makefile`s and the build-command comments in the example files. Headers and
object files are untouched — `#include`s point at `.h`, and `Foo.o` derives from the
source basename regardless of extension, so no `.o` reference changes. This makes
the rename strictly safer than the `exp`→`ext` one, which touched `#include`s.

## 2026-07-08 — Explicit instantiation over header-only templates (rationale)

Decision (already active in CLAUDE.md; this entry records *why*, which otherwise
lives only in the `experiments/template-instantiation/` example comments): Val-dependent classes keep a single
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
full implicit / plain explicit / guarded explicit framing. So Oblio's pattern was achievable in C++98; C++11 was
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
results, and the plain explicit and guarded explicit variants share the same link-failure behaviour when
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
