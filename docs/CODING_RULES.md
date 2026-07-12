# Coding Rules

Conventions not enforced by tooling. Mechanical style, formatting and idiom
modernization, is enforced by `.clang-format` and `.clang-tidy`; this file covers
only what a linter can't check. Language-general; add per-language sections as the
project grows (e.g. a Rust section later).

**Invariants** (must not break) live in CLAUDE.md, always loaded. This file is the
softer layer: conventions for consistency, not correctness.

## Project conventions (all languages)

- Symmetric naming for paired files (`*_real` / `*_complex`, not asymmetric pairs).
- Flat source layout, no per-category subdirectories.
- Consistent capitalization across filenames.
- Compiled executables carry a language suffix indicating what built them: `_cpp`
  for C++ (`_rs` for Rust, etc.). E.g. `test_smoke_real_cpp`, `test_multiply_implicit_cpp`.
  Applies to the executable's output name, not source files. Keep `.gitignore`, the
  CMake target/output names, and any Makefile targets in sync when adding one.
- **Optimization level: `-O3`.** Manual builds and Makefiles use `-O3` (the highest
  standard level). Use it consistently across every build file. (CMake derives its
  own level from the build type, Release is `-O3`; see the CMake note when that tree
  is set up.)
- **Package name is `Oblio`** (capital O), in prose, documentation, and the C++
  namespace (`namespace Oblio`, `Oblio::`). Lowercase `oblio` is reserved for path
  and artifact identifiers only: the `include/oblio/` directory and `#include
  "oblio/…"` paths, `liboblio`, the CMake target, and existing lowercase filenames
  (`oblio-new-devlog.md`, etc.). Macros are `OBLIO_`. Never write bare `oblio` as
  the package name in running text.

- **Definition order follows declaration order.** The `.cpp` defines functions in the same
  order the header declares them. The header is the table of contents; the source is the
  book, and they should agree. A reader who has found something in the header knows where
  to look for it, and a diff that reorders one without the other is a review hazard for no
  gain. This is a convention, not a correctness rule, which is exactly why nothing enforces
  it and it drifts unless stated.
- **Keep an overload set together, in both files.** Declare the members of an overload set
  adjacently in the header and define them adjacently in the `.cpp`, and do not bunch all
  the adapters in one place and all the implementations in another. An overload pair is one
  idea in two spellings; splitting it makes the reader hold the connection in their head
  instead of seeing it. The pattern this most applies to:

  - **The overload taking less is the implementation; the one taking more is an adapter.**
    Where a function needs only part of an object, take the part in the body and the object
    at the boundary, as a one-line forward. E.g. `OrderEngine::compute(const SparseMatrix<Val>&,
    Permutation&)` forwards to `compute(colPtr, rowIdx, Permutation&)`. The implementation
    then states its real dependencies, and a caller who has only that part can reach it. See
    the structural-part rule below for why this matters beyond signature length.

    **Adapt once, at the public boundary.** One adapter per engine, on the entry point, and
    nothing below it takes more than it needs. Do not repeat the pattern on each private
    helper: once the entry point has unpacked the object, the helpers already have the part
    they need, and a second layer of adapters would have no callers. All three engines
    (`OrderEngine`, `ElmForestEngine`, `SymFactEngine`) have exactly one adapter, on `compute`.

- **Prefer passing the whole object over passing its pieces**, where doing so does not drag
  in a template parameter the callee has no use for. Passing an object's fields to a function
  that is already its `friend` restricts nothing, it only lengthens the signature, and a long
  run of same-typed reference parameters (five `std::vector<std::int32_t>&` in a row) is a
  transposition bug the compiler cannot catch. What a signature cannot say anyway (which
  fields the callee leaves *stale*) belongs in the comment regardless. See DESIGN_DECISIONS.

- **Access rule: friendship is a write grant.** Three cases, and they cover every argument an
  engine takes.

  - **Written -> friend, pass the object.** An engine is declared `friend` by exactly the
    object it fills, and by no other. `ElmForestEngine` writes `ElmForest`, `SymFactEngine`
    writes `SymFact`, `OrderEngine` writes `Permutation`. Pass the object, not its fields:
    the engine can reach them anyway, so listing them only lengthens the signature. What a
    signature cannot say (which fields the callee leaves *stale*) goes in the comment.
  - **Read -> pass the object, public API, no friendship.** Reading needs no friend. The
    accessors return `const&`, so there is no copy to avoid and no access to gain; friendship
    would spend encapsulation for nothing. `SymFactEngine` reads a dozen fields of `ElmForest`
    through its accessors and is not its friend.
  - **Read, and only the structure is needed -> take only the structure.** See below.

- **Pass only the structural part of a matrix when only structural work is done.** Ordering,
  the elimination forest and the symbolic factorization are graph algorithms: they read a
  sparsity pattern and never touch a value. So they take one, and `SparseMatrix` offers two
  overloads, one taking the matrix and one taking `colPtr` and `rowIdx`, the second holding
  the implementation and the first a one-line adapter over it.

  This is about honest dependencies, not about C++. A function's parameters should say what it
  actually consumes. A structural algorithm that demands a matrix is lying about what it needs,
  and it forces a caller holding only a graph to fabricate numbers to satisfy a signature that
  will ignore them. The rule would be right in any language.

  The mechanism, in C++, is that the structural overload is not templated, which is also why
  it is compiled once rather than once per scalar type. That is a consequence and not the
  motive, and it should not be mistaken for one.

## C++

- **Index naming: `j` and `k` are columns, `jj` and `kk` are their supernodes.** 0.9's
  convention, and worth keeping. Two rules, and the second is the useful one:

  - **`j` is the lower column, `k` the higher**, whichever the loop happens to scan. An
    ascending loop scans `k` and looks down at its child `j` (`j = firstChild[k]`); a
    descending loop scans `j` and looks up at its parent `k` (`k = parent[j]`). Either way
    `j < k`, because a parent's column is numbered above its child's.
  - **Doubling a letter applies the column-to-supernode map.** `jj` is the supernode of `j`,
    `kk` the supernode of `k`. One rule, applied at both levels.

  The payoff is that a statement reads as the fact it derives. In the compression rebuild,
  `parent[jj] = kk` says: the supernode of `j`'s parent column is the parent of `j`'s
  supernode. Both levels are visible, and the relation between them needs no lookup. The rule
  holds even where no column is in scope: in `finalizeLinks` and `computeHeight`, which see
  only supernodes, they are still `jj` and `kk`. A doubled letter means supernode, full stop.

  Where both orderings are in play, prefix with `l` for the factor (permuted) ordering and `a`
  for the original ordering of `A`: `lj`, `lk`, `aj`, `ak`. This keeps 10.12's ordering
  distinction and 0.9's column roles at once, where each reference had only one of them.

  The same prefix marks a **position** in the flat storage: `ap` is an offset into `A`'s
  arrays, `lp` an offset into the factor's. This is 0.9's `p`, the name `p` being taken by the
  `Permutation`. Positions measure rather than name, so they are `std::size_t`, cannot be
  `NIL`, and may exceed 2^31; the ID rules above do not apply to them.

  Row indices are `i`. Scratch is `r` for the node a climb has reached and `t` for a general
  temporary.

  See the index-type rules below for why `j`, `k`, `jj`, `kk` are `std::int32_t` and are used
  to subscript without a cast.

- **Modern spellings, pin one per historical variation** (check this list before
  reintroducing an old form):
  - source files: **`.cpp`**, not `.cc` (headers `.h`)
  - size/index type: **`std::size_t`** for offsets/counts, **`std::int32_t`** for IDs,
    never bare `size_t`/`int32_t`, see the std-types and index-type rules below
  - aliases: **`using`**, not `typedef`
  - null pointer: **`nullptr`**, not `NULL`
  - enums: **`enum class`**, not bare `enum`
  - deleted members: **`= delete`**, not private-undeclared

  The last four are auto-enforced by `.clang-tidy` (`modernize-*`); `.cpp` is
  convention; `std::size_t` needs the header detail below. This list is the single
  place to settle "which spelling" so it doesn't resurface per file.

  Note: `.clang-tidy` does **not** read this document, the tool works from its own
  `Checks:` config (running `modernize-*` AST rewrites), and the two are independent
  records of the same intent. Nothing detects a disagreement between them, so keep
  this list and `.clang-tidy` in sync by hand. "Auto-enforced" also assumes something
  actually runs `clang-tidy` (editor, pre-commit, or CI) against a compilation
  database; until that's wired up, the checks fire only when run manually.

- **Member naming: `mFoo` prefix**, house style, consistent with 0.9 (`mRows`,
  `mCols`, `mVals`). Arguments and locals are bare by default: they can't collide
  with each other, and `mFoo` already separates them from members. If marking an
  argument genuinely aids clarity, use `aFoo`, *not* `foo_` (reads as "member" in
  Google/LLVM style) and *not* leading `_foo` (brushes the reserved-identifier
  rules). *Provisional:* may switch to the more common trailing-underscore `foo_`
  for members; if so, this rule and every existing `mFoo` change together (a
  mechanical rename), and `aFoo` for arguments still stands.
- **Standard-library types: use the `std::`-qualified name, from the C++ header that
  declares it**, not the bare name and not the C-style `<*.h>` header. Include only
  the header for the type you actually use (one header per need):
  - `<cstddef>` → `std::size_t`, `std::ptrdiff_t`, `std::nullptr_t`, `std::byte`
  - `<cstdint>` → fixed-width ints: `std::int32_t`, `std::uint64_t`, `std::int8_t`, …

  e.g. `std::size_t` comes from `<cstddef>` (which is all Permutation needs). The
  bare forms (`size_t`) rely on a global-namespace leak that isn't guaranteed by the
  C++ headers. This is the correct spelling even where a matching codebase uses the
  bare form, that codebase is the one to fix, not this one.
- **Index types, IDs vs offsets** (the graph/matrix convention; see the design
  decision for the full rationale):
  - **IDs**, values that name a vertex/row/column/supernode and may carry a "none"
    sentinel → **`std::int32_t`**. E.g. `SparseMatrix::rowIdx`, permutation maps,
    `ElmForest` parent/child/sibling/supernode-map arrays. Sentinel is
    **`NIL = -1`** (a `constexpr std::int32_t`), never `static_cast<std::size_t>(-1)`.
    Cost: ~2.1 billion index cap, accepted for a clean signed sentinel, matching the
    graph code and the vendored `int`-based AMD/MMD.
  - **Offsets / counts / sizes**, row-pointers, `nnz`, dimensions, anything that
    indexes-into or measures → **`std::size_t`**. Never negative, may exceed 2^31.
    E.g. `SparseMatrix::colPtr`, `size()`, `nnz()`.
  - **Loop counters are `std::size_t`.** A counter enumerates real positions (never
    negative, never the sentinel), so the unsigned view is safe and avoids
    signed/unsigned comparison against `.size()`. Do **not** loop with an `int32_t`
    counter compared to a container size.
  - **Cast only where information can be lost, which is one direction, not two:**
    `static_cast<std::int32_t>(x)` when a `std::size_t` becomes an ID (**narrowing**; this is
    where the 2^31 cap lives, and the compiler warns without it). The other direction needs
    nothing: subscripting a container with an `std::int32_t` is legal, is silent under
    `-Wall -Wextra`, and is correct for every value that is not `NIL`. Write `parent[jj]`, not
    `parent[static_cast<std::size_t>(jj)]`.
  - **A widening cast is not a NIL guard, and must not be mistaken for one.**
    `static_cast<std::size_t>(NIL)` yields `SIZE_MAX` exactly as silently as the implicit
    conversion would: the cast performs the mistake, it does not prevent it. What prevents it
    is a branch, `if (kk != NIL)`, and that branch is required whether or not a cast is
    written. Adding casts "for safety" buys nothing and costs a great deal of noise.
  - **One name per entity.** `j`, `k`, `jj`, `kk` are `std::int32_t`, because they name things
    and may carry `NIL`. Do not introduce a second variable holding the same value as a
    `std::size_t`: a reader then has to check whether the two are really the same thing. Where
    an ascending loop counts entities, the counter is that entity's type (`std::int32_t`) and
    the *bound* is cast once in the loop header, rather than the counter being cast at every
    use. Where a loop counts *positions* rather than entities (an offset into `colPtr`, or a
    descending count-down), it is a `std::size_t` and names nothing.

- Beyond the spellings above, `.clang-tidy` (`modernize-*`) and `.clang-format` also
  handle idiom cleanups (e.g. `.data()` over raw-pointer extraction) and formatting.
  Rely on the tools; don't hand-police these.
- The judgment calls tools can't make are invariants, in CLAUDE.md: port-verbatim,
  container choice, index signedness, don't-flatten-mid-port.
