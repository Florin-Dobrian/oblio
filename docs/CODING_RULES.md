# Coding Rules

Conventions not enforced by tooling. Mechanical style — formatting and idiom
modernization — is enforced by `.clang-format` and `.clang-tidy`; this file covers
only what a linter can't check. Language-general; add per-language sections as the
project grows (e.g. a Rust section later).

**Invariants** (must not break) live in CLAUDE.md, always loaded. This file is the
softer layer: conventions for consistency, not correctness.

## Project conventions (all languages)

- Symmetric naming for paired files (`*_real` / `*_complex`, not asymmetric pairs).
- Flat source layout — no per-category subdirectories.
- Consistent capitalization across filenames.
- Compiled executables carry a language suffix indicating what built them: `_cpp`
  for C++ (`_rs` for Rust, etc.). E.g. `test_smoke_real_cpp`, `test_multiply_implicit_cpp`.
  Applies to the executable's output name, not source files. Keep `.gitignore`, the
  CMake target/output names, and any Makefile targets in sync when adding one.
- **Optimization level: `-O3`.** Manual builds and Makefiles use `-O3` (the highest
  standard level). Use it consistently across every build file. (CMake derives its
  own level from the build type — Release is `-O3`; see the CMake note when that tree
  is set up.)
- **Package name is `Oblio`** (capital O) — in prose, documentation, and the C++
  namespace (`namespace Oblio`, `Oblio::`). Lowercase `oblio` is reserved for path
  and artifact identifiers only: the `include/oblio/` directory and `#include
  "oblio/…"` paths, `liboblio`, the CMake target, and existing lowercase filenames
  (`oblio-new-devlog.md`, etc.). Macros are `OBLIO_`. Never write bare `oblio` as
  the package name in running text.

## C++

- **Modern spellings — pin one per historical variation** (check this list before
  reintroducing an old form):
  - source files: **`.cpp`**, not `.cc` (headers `.h`)
  - size/index type: **`std::size_t`** for offsets/counts, **`std::int32_t`** for IDs,
    never bare `size_t`/`int32_t` — see the std-types and index-type rules below
  - aliases: **`using`**, not `typedef`
  - null pointer: **`nullptr`**, not `NULL`
  - enums: **`enum class`**, not bare `enum`
  - deleted members: **`= delete`**, not private-undeclared

  The last four are auto-enforced by `.clang-tidy` (`modernize-*`); `.cpp` is
  convention; `std::size_t` needs the header detail below. This list is the single
  place to settle "which spelling" so it doesn't resurface per file.

  Note: `.clang-tidy` does **not** read this document — the tool works from its own
  `Checks:` config (running `modernize-*` AST rewrites), and the two are independent
  records of the same intent. Nothing detects a disagreement between them, so keep
  this list and `.clang-tidy` in sync by hand. "Auto-enforced" also assumes something
  actually runs `clang-tidy` (editor, pre-commit, or CI) against a compilation
  database; until that's wired up, the checks fire only when run manually.

- **Member naming: `mFoo` prefix** — house style, consistent with 0.9 (`mRows`,
  `mCols`, `mVals`). Arguments and locals are bare by default: they can't collide
  with each other, and `mFoo` already separates them from members. If marking an
  argument genuinely aids clarity, use `aFoo` — *not* `foo_` (reads as "member" in
  Google/LLVM style) and *not* leading `_foo` (brushes the reserved-identifier
  rules). *Provisional:* may switch to the more common trailing-underscore `foo_`
  for members; if so, this rule and every existing `mFoo` change together (a
  mechanical rename), and `aFoo` for arguments still stands.
- **Standard-library types: use the `std::`-qualified name, from the C++ header that
  declares it** — not the bare name and not the C-style `<*.h>` header. Include only
  the header for the type you actually use (one header per need):
  - `<cstddef>` → `std::size_t`, `std::ptrdiff_t`, `std::nullptr_t`, `std::byte`
  - `<cstdint>` → fixed-width ints: `std::int32_t`, `std::uint64_t`, `std::int8_t`, …

  e.g. `std::size_t` comes from `<cstddef>` (which is all Permutation needs). The
  bare forms (`size_t`) rely on a global-namespace leak that isn't guaranteed by the
  C++ headers. This is the correct spelling even where a matching codebase uses the
  bare form — that codebase is the one to fix, not this one.
- **Index types — IDs vs offsets** (the graph/matrix convention; see the design
  decision for the full rationale):
  - **IDs** — values that name a vertex/row/column/supernode and may carry a "none"
    sentinel → **`std::int32_t`**. E.g. `SparseMatrix::rowIdx`, permutation maps,
    `ElmForest` parent/child/sibling/supernode-map arrays. Sentinel is
    **`NIL = -1`** (a `constexpr std::int32_t`), never `static_cast<std::size_t>(-1)`.
    Cost: ~2.1 billion index cap — accepted for a clean signed sentinel, matching the
    graph code and the vendored `int`-based AMD/MMD.
  - **Offsets / counts / sizes** — row-pointers, `nnz`, dimensions, anything that
    indexes-into or measures → **`std::size_t`**. Never negative, may exceed 2^31.
    E.g. `SparseMatrix::colPtr`, `numCols()`, `nnz()`.
  - **Loop counters are `std::size_t`.** A counter enumerates real positions (never
    negative, never the sentinel), so the unsigned view is safe and avoids
    signed/unsigned comparison against `.size()`. Do **not** loop with an `int32_t`
    counter compared to a container size.
  - **Cast explicitly at the two crossings, and only there:**
    `static_cast<std::int32_t>(counter)` when storing a counter into an ID array
    (narrowing — this is where the 2^31 cap lives); `static_cast<std::size_t>(id)`
    when an ID subscripts an array (widening — guard against `NIL` first if the ID
    could be a sentinel). These casts mark the ID↔offset boundary; they should be few.
- Beyond the spellings above, `.clang-tidy` (`modernize-*`) and `.clang-format` also
  handle idiom cleanups (e.g. `.data()` over raw-pointer extraction) and formatting.
  Rely on the tools; don't hand-police these.
- The judgment calls tools can't make are invariants, in CLAUDE.md: port-verbatim,
  container choice, index signedness, don't-flatten-mid-port.
