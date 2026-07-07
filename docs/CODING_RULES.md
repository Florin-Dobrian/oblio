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

## C++

- Mechanical modernization (`nullptr`, `using` over `typedef`, `enum class`,
  `= delete`, `.data()`) is handled by `.clang-tidy` (`modernize-*`) and
  `.clang-format`. Rely on the tools; don't hand-police these.
- The judgment calls tools can't make are invariants, in CLAUDE.md: port-verbatim,
  container choice, index signedness, don't-flatten-mid-port.
