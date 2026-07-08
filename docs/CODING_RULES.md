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
- **Package name is `Oblio`** (capital O) — in prose, documentation, and the C++
  namespace (`namespace Oblio`, `Oblio::`). Lowercase `oblio` is reserved for path
  and artifact identifiers only: the `include/oblio/` directory and `#include
  "oblio/…"` paths, `liboblio`, the CMake target, and existing lowercase filenames
  (`oblio-new-devlog.md`, etc.). Macros are `OBLIO_`. Never write bare `oblio` as
  the package name in running text.

## C++

- **Member naming: `mFoo` prefix** — house style, consistent with 0.9 (`mRows`,
  `mCols`, `mVals`). Arguments and locals are bare by default: they can't collide
  with each other, and `mFoo` already separates them from members. If marking an
  argument genuinely aids clarity, use `aFoo` — *not* `foo_` (reads as "member" in
  Google/LLVM style) and *not* leading `_foo` (brushes the reserved-identifier
  rules). *Provisional:* may switch to the more common trailing-underscore `foo_`
  for members; if so, this rule and every existing `mFoo` change together (a
  mechanical rename), and `aFoo` for arguments still stands.
- Mechanical modernization (`nullptr`, `using` over `typedef`, `enum class`,
  `= delete`, `.data()`) is handled by `.clang-tidy` (`modernize-*`) and
  `.clang-format`. Rely on the tools; don't hand-police these.
- The judgment calls tools can't make are invariants, in CLAUDE.md: port-verbatim,
  container choice, index signedness, don't-flatten-mid-port.
