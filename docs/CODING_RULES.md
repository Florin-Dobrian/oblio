# Coding Rules

Conventions not enforced by tooling. Mechanical style, formatting and idiom
modernization, is enforced by `.clang-format` and `.clang-tidy`; this file covers
only what a linter can't check. Language-general; add per-language sections as the
project expands (e.g. a Rust section later).

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
  "oblio/..."` paths, `liboblio`, the CMake target, and existing lowercase filenames
  (`oblio-new-devlog.md`, etc.). Macros are `OBLIO_`. Never write bare `oblio` as
  the package name in running text.

- **Terminology lives in WRITING_RULES.** Its Terminology section is the one part of that file
  that governs code as well as prose, because names and comments must use the same words the
  documentation does, or a reader ends up carrying a translation table. It settles
  **numeric** factorization (not "numerical", which is reserved for the mathematics: numerical
  stability, numerical pivoting), **front** as the noun against **frontal** as the adjective,
  and when to reach for matrix words, graph words, or the neutral `index` and `snode`.

- **Data structures take the noun; functions may take the adjective.** Fields and accessors are
  strict: a supernode has a **front** (its columns), so the field is `mFrontSize` and the
  accessor `frontSize()`. Function names can afford to be descriptive, so `gatherFrontalIndices`
  reads as prose in a declaration. Both are correct, and the split is the ordinary one: the
  front is a thing, a frontal index is a kind of index. Prose documentation follows the fields,
  since that is what a reader cross-references.

  Function names also spell words out where a field would abbreviate: `sortIndices`, not
  `sortIdx`. Loop-local variables go the other way and abbreviate freely (`sfp`, `sp`, `lk`),
  because there they are read a hundred times and written once.

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
    (`OrderEngine`, `ElmForestEngine`, `SymFactorEngine`) have exactly one adapter, on `compute`.

- **Prefer passing the whole object over passing its pieces**, where doing so does not drag
  in a template parameter the callee has no use for. Passing an object's fields to a function
  that is already its `friend` restricts nothing, it only lengthens the signature, and a long
  run of same-typed reference parameters (five `std::vector<std::int32_t>&` in a row) is a
  transposition bug the compiler cannot catch. What a signature cannot say anyway (which
  fields the callee leaves *stale*) belongs in the comment regardless. See DESIGN_DECISIONS.

- **Access rule: friendship is a write grant.** Three cases, and they cover every argument an
  engine takes.

  - **Written -> friend, pass the object.** An engine is declared `friend` by exactly the
    object it fills, and by no other. `ElmForestEngine` writes `ElmForest`, `SymFactorEngine`
    writes `SymFactor`, `OrderEngine` writes `Permutation`. Pass the object, not its fields:
    the engine can reach them anyway, so listing them only lengthens the signature. What a
    signature cannot say (which fields the callee leaves *stale*) goes in the comment.
  - **Read -> pass the object, public API, no friendship.** Reading needs no friend. The
    accessors return `const&`, so there is no copy to avoid and no access to gain; friendship
    would spend encapsulation for nothing. `SymFactorEngine` reads a dozen fields of `ElmForest`
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

- **A switch over one of our enums names every enumerator, and never uses `default`.** Under
  `-Wall` the compiler warns on a switch that misses one, and that warning is the point: adding a
  sixth `Factorization` should make the build list every place that has to make a decision about it.
  A `default` throws that away. The new enumerator compiles, falls into the catch-all, and returns
  something plausible from a site nobody revisited.

  This holds even where a case is unreachable. `NumFactorEngine::compute`'s dynamic overload
  dispatches the dynamic factorizations before its static switch, so those two enumerators cannot
  reach it, but they are still listed and still return false, because the compiler cannot know they
  are unreachable and neither can the next reader. An unreachable named case costs one line; a
  `default` costs the warning.

  Corollary for the comment beside a refusal: **say whether it is a gap or a decision.** "Not
  implemented" invites someone to implement it. Dynamic pivoting into `NumFactorStatic` returns
  false permanently, because a delayed column expands a front and a flat buffer cannot expand, and the
  comment says so; multifrontal returns false because nobody has ported it yet, and the comment says
  that instead. Same statement, opposite instruction to the reader.

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

  **Two kinds of thing, and they must not be confused: an index and a position.**

  An **index** identifies a matrix or graph entity: a row, a column, a node, a supernode. These
  are the same kind of thing (rows, columns and vertices share one index space, the structure
  being symmetric; supernodes have their own space), and any of them may need to say "none", so
  an index is a **`std::int32_t`** and its sentinel is `NIL`.

  A **position** is an offset into a vector. It identifies nothing; it locates something inside
  one particular flat array, and means nothing outside it. A position is a **`std::size_t`**: it
  measures, so it is never negative, never `NIL`, and free to exceed 2^31.

  **An index becomes a position when we use it to select from a vector.** That is the one place
  the two meet, and it is where care is needed: turning `NIL` into a position yields `SIZE_MAX`
  and reads a mile past the end. The remedy is a guard (`if (kk != NIL)`), never a cast; see the
  index-type rules below, where this is spelled out.

  So the storage rule falls out. **A vector holding indices is a `std::vector<std::int32_t>`, and
  a vector holding positions is a `std::vector<std::size_t>`. Both are accessed by position.**
  `rowIdx` holds indices; `colPtr` holds positions; both are subscripted by a position.

  **Do not say "index into a vector".** It is ordinary English, and that is exactly the problem:
  `index` is reserved here for a matrix or graph entity, and reaching for it to mean array access
  collides with that meaning at the moment the distinction matters most. Say **position**.
  (*Offset* would also do, but it suggests a delta rather than a location, so position is
  better.)

  **And `colPtr` is not a pointer.** Neither are `snodePtr` or `snodeFrontPtr`. They
  hold positions, `std::size_t` offsets into a flat array. The `Ptr` is inherited vocabulary,
  universal in sparse
  matrix codes (0.9, 10.12, CSparse, SuiteSparse all say `colptr` or `Ap`), and we keep it for
  that reason alone. It is a name, not a claim about the type.

  **Indices are global.** There are exactly six, and there will only ever be six, because there
  are two orderings and three roles:

  ```
  aj  ak  ai      a column, a column and a row, in A's ordering
  lj  lk  li      the same three, in the factor's ordering
  ```

  Rows are `i`, columns are `j` and `k`, with `j` the lower and `k` the higher. They do turn into
  rows where the two meet (`A[k][j]` names an entry, so `k` indexes a row there), and that is not
  a problem: rows and columns share one index space, the structure being symmetric, so an index
  is an index and the letter says what it is *for*, not what it *is*. Every structure in the
  solver, `A`, the forest, the symbolic factor, the numeric factor, indexes into that one space.

  **Strictly, the prefix is the ordering, not the matrix.** "`a` for `A`, `l` for `L`" is the
  mnemonic, and a good one, but what the prefix records is *which numbering the index is in*:
  `A`'s original one, or the permuted one the factor is computed in. There is exactly one of
  each, and they are connected by exactly one thing, the `Permutation`: `oldToNew[aj] == lj` and
  `newToOld[lk] == ak`. A line that converts between them should read as the conversion it is.

  The distinction matters because the factor is not one matrix. Cholesky computes `C` with
  `A = CC^T`; static and dynamic LDL compute `L` with `A = LDL^T`, unit lower triangular, not
  absorbing the square root of `D`. Numerically these are different factors, and we write `C`
  for the first precisely to keep `L` for the second (WRITING_RULES). But they
  share one structure and one index space, which is why the elimination forest and the symbolic
  factorization are computed once and serve all three. So `lk` means the same thing regardless of
  which factorization is running, and it should: an index has no business knowing.

  (10.12 writes `lc`/`lr` for column and row and drops the role distinction; 0.9 writes bare `j`
  and `k` and drops the ordering distinction. We keep both, at a cost of one character.)

  **A position is named for the pointer array it walks: its initials, and nothing else.** That is
  the whole rule, and it has no exceptions.

  | position | walks | one entry per | the data is |
  |---|---|---|---|
  | `cp` | `colPtr` / `rowIdx` | column | row indices |
  | `sp` | `snodePtr` / `nodeIdx` | supernode | node indices |
  | `np` | `nodePtr` / `snodeIdx` | node | supernode indices |
  | `sfp` | `snodeFrontPtr` / `frontNodeIdx` | supernode | front (column) indices |

  Read across the table and the CSC convention falls out: **the letter names what the array has
  one entry of, and the data is the other thing.** `colPtr` has one entry per column and holds
  rows. `snodePtr` has one per supernode and holds nodes. `nodePtr` has one per node
  and holds supernodes.

  `sp` and `np` are exact mirrors, and worth holding onto as the clearest case: `sp` is a *snode*
  pointer, so its data is nodes; `np` is a *node* pointer, so its data is supernodes. That is why
  the two counting sorts in `sortIndices` are the same code with the roles exchanged, and why the
  names stay honest across both passes without needing a second set.

  `sfp` is not an exception but the rule applied to a longer name: `snodeFrontPtr` -> `sfp`.
  It is a `sp` restricted to the front, and the letters say so.

  The letters name the *array*, not the matrix. There is only one `colPtr` in the solver, so `cp`
  is unambiguous. Should two objects ever both carry one, qualify then and not before: `acp`,
  `lcp`.

  0.9 calls all of these `p`; we cannot, since `p` is the `Permutation`, and in any case the
  prefix is more useful than the bare letter. Unlike indices there is no fixed set of positions:
  each flat structure brings its own, and `lp` will arrive with the numeric factor.

  **A position is not an index.** It cannot be compared with one, it means nothing outside its
  own array, and it is a `std::size_t` because it measures rather than names: never negative,
  never `NIL`, free to exceed 2^31, and outside the index-type rules below. The two meet only through a
  dereference, `lk = frontNodeIdx[sfp]`, and code should make that step visible rather than nest
  it.

  Scratch is `r` for the node a climb has reached and `t` for a general temporary.

  See the index-type rules below for why `j`, `k`, `jj`, `kk` are `std::int32_t` and are used
  to subscript without a cast.

- **Modern spellings, pin one per historical variation** (check this list before
  reintroducing an old form):
  - source files: **`.cpp`**, not `.cc` (headers `.h`)
  - index types: **`std::int32_t`** for indices (which name something and may be `NIL`),
    **`std::size_t`** for positions (which measure), never bare `size_t`/`int32_t`; see the
    std-types and index-type rules below
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
  - `<cstddef>` -> `std::size_t`, `std::ptrdiff_t`, `std::nullptr_t`, `std::byte`
  - `<cstdint>` -> fixed-width ints: `std::int32_t`, `std::uint64_t`, `std::int8_t`, ...

  e.g. `std::size_t` comes from `<cstddef>` (which is all Permutation needs). The
  bare forms (`size_t`) rely on a global-namespace leak that isn't guaranteed by the
  C++ headers. This is the correct spelling even where a matching codebase uses the
  bare form, that codebase is the one to fix, not this one.
- **Index types.** The same index/position split as in the naming rules above, now as types.
  (The design decision calls these IDs and offsets, which is the older vocabulary for the same
  thing.)

  - **An index is a `std::int32_t`.** It names a vertex, row, column or supernode, and may
    carry the sentinel **`NIL = -1`** (a `constexpr std::int32_t`, never
    `static_cast<std::size_t>(-1)`). E.g. `SparseMatrix::rowIdx`, the permutation maps, the
    forest's parent/child/sibling/map arrays. Cost: a ~2.1 billion cap, accepted for a clean
    signed sentinel, and matching the graph code and the vendored `int`-based AMD/MMD.
  - **A position is a `std::size_t`.** It measures rather than names: an offset into a buffer,
    a count, a dimension. Never negative, never `NIL`, free to exceed 2^31. E.g.
    `SparseMatrix::colPtr`, `size()`, `nnz()`.
  - **A loop counter takes the type of what it counts, and direction does not change the type.**
    A loop over *entities* (columns, supernodes) counts indices, so the counter is `std::int32_t`
    and the bound is cast once in the condition. Ascending:
    `for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk)`. Descending is the
    same rule with the cast on the start:
    `for (std::int32_t jj = static_cast<std::int32_t>(snodeSize) - 1; jj >= 0; --jj)`; the loop exits
    at `-1`, which `int32_t` represents, so a descending entity loop is a *direct* `int32_t`
    reverse, not a `std::size_t` counter with `j = t - 1` inside (that trick is declined; see the
    forward/reverse asymmetry in DESIGN_DECISIONS). The cast is loop-invariant and costs nothing. A
    loop over *positions* (an offset into `colPtr`) counts nothing that has a name, so it is
    `std::size_t`: `for (std::size_t ap = colPtr[ak]; ap < colPtr[ak + 1]; ++ap)`. A descending
    *position* loop is still `std::size_t`, but must not test `p >= 0` (always true, it never
    terminates); use `while (p-- > 0)` or an explicit guard. What is forbidden is an `int32_t`
    counter compared against a container size *without* the cast, which is where the signed/unsigned
    warning comes from.

  Three rules about casting, and they matter more than they look, since we got each of them
  wrong at least once.

  - **Cast in one direction only, the one that loses information.**
    `static_cast<std::int32_t>(x)`, when a `std::size_t` becomes an index. That is where the
    2^31 cap lives, and the compiler warns without it. The other direction needs nothing:
    subscripting a container with an `std::int32_t` is legal, is silent under `-Wall -Wextra`,
    and is correct for every value that is not `NIL`. Write `parent[jj]`, never
    `parent[static_cast<std::size_t>(jj)]`.
  - **A widening cast is not a NIL guard, and must not be mistaken for one.**
    `static_cast<std::size_t>(NIL)` yields `SIZE_MAX` exactly as silently as the implicit
    conversion would: **the cast performs the mistake, it does not prevent it.** What prevents
    it is a branch, `if (kk != NIL)`, and that branch is required whether or not a cast is
    written. Casts added "for safety" buy nothing and cost a great deal of noise.
  - **One name per entity, and this includes sizes.** Do not introduce a second variable holding
    the same value in the other type: a reader then has to check whether the two are really the
    same thing. This applies to **counts and sizes** as much as to indices, and
    `const std::int32_t n = static_cast<std::int32_t>(size);` is the same mistake wearing
    different clothes. It is easy to make precisely because a size feels too humble to need the
    rule.

- **Constructors initialize every member in the member-initializer list, in declaration order.**
  The body is left for work an initializer cannot do: a loop that fills, a call made for its side
  effect (`setIdentity()`). A member set by assignment in the body is default-constructed first and
  then overwritten, which for a `std::vector` is a wasted allocation the list would have avoided.
  Where the initializer is not a plain copy, use an expression, not a body assignment and not a
  helper with internal linkage: `mNnz(std::accumulate(rowIdx.begin(), rowIdx.end(), std::size_t{0},
  ...))`, and note the `std::size_t{0}` seed is load-bearing, it fixes the accumulator type so the
  sum does not overflow an `int`. Two consequences: the list runs in *declaration* order whatever
  order it is written in, so declare a member after the ones it reads; and a guard on the first
  member runs before any later member allocates, which is how check-before-allocate and init-list
  purity coexist (see the next rule).

- **A size-taking constructor guards its size against the index range before it allocates.** A
  dimension or count that will size a `std::vector` of indices, or be stored as one, must fit
  `std::int32_t`. It passes through `checkIndexRange(size, "<what> size")` (declared in `Types.h`,
  defined in `Types.cpp`), which returns the size unchanged or throws `std::length_error` when it
  exceeds `MAX_IDX` (`= 2^31 - 1`). Put the call in the member-initializer list on the *first*
  member the size feeds, so it runs before any allocation:
  `Vector(std::size_t size) : mSize(checkIndexRange(size, "Vector size")), mVal(size, Val(0)) {}`. Where the
  constructor instead moves in already-built vectors (`SparseMatrix`), nothing is allocated from the
  size, so the guard is a plain body check. The throw lives in `Types.cpp` and never inline in a
  header: a throwing body in a widely-included header degrades codegen of unrelated hot loops in the
  same translation unit (the in-header-throw finding in DESIGN_DECISIONS). A structure sized from
  user input or from A carries the guard (`SparseMatrix`, `Vector`, `Permutation`); a structure
  sized from already-guarded internal data does not re-check (the numeric factor, sized from the
  symbolic factor, throws nothing). See DESIGN_DECISIONS for where the cap comes from and why it is
  `2^31 - 1`, not `2^31`.

- Beyond the spellings above, `.clang-tidy` (`modernize-*`) and `.clang-format` also
  handle idiom cleanups (e.g. `.data()` over raw-pointer extraction) and formatting.
  Rely on the tools; don't hand-police these.
- The judgment calls tools can't make are invariants, in CLAUDE.md: port-verbatim,
  container choice, index signedness, don't-flatten-mid-port.
