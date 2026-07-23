# CLAUDE.md

Operating contract for the Oblio refactoring effort. Auto-loaded by Claude Code
each session, keep it lean. Rationale lives in docs/DESIGN_DECISIONS.md, not here.

## What this is

Oblio is a sparse direct solver (supernodal Cholesky / LDL^T). This effort is a
careful refactor of **Oblio 0.9** into modern C++17: migrate the hand-rolled
`Array` container to `std::vector`, modernize idioms, one function at a time.
Porting and modernizing are the same act, a good algorithm re-expressed in
modern style. They are not two passes.

## Source of truth

Two references, two roles: favor 10.12's *design*, verify against 0.9's *behavior*.

- **Oblio 0.9** is the correctness oracle: complete, correct, fully functional. It
  has everything. Every ported unit is verified against 0.9's output before it counts
  as done, and when 0.9 and 10.12 disagree on *what is computed*, 0.9 wins.
- **Oblio 10.12** is a partial, unfinished refactor of 0.9 headed toward a better
  design, and it is the design direction to favor: prefer 10.12's structure,
  decomposition, and interfaces wherever it reaches. It is unverified, so never assume
  its code is correct (it can carry bugs 0.9 does not). Take *shape* from 10.12 and
  *behavior* from 0.9, and fall back to 0.9 wherever 10.12 is silent or incomplete.
- **PoC tree** (built in one day, 2026-03-07) was a proof of concept. It informed the
  effort but was never trusted, and nothing was ported from it. It has since been removed
  from the tree.

**Where things live.** The 0.9 and 10.12 sources live in **`reference/`**, which is
gitignored and so absent from REPO_MAP; files are suffixed `-0.9` and `-10.12`. All work
happens at the top level (`src/`, `include/oblio/`, `tests/`). Nothing else in the tree is
either a reference or a target.

## Process (every session)

1. **One unit at a time.** One function or small function-group per step. Verify,
   record in the ledger, then advance.
2. **Port, then modernize, then verify.** (a) Carry 0.9's algorithm over faithfully,
   human judgment, no tool checks this. (b) Run `clang-tidy --fix` to sweep the
   mechanical idioms (`modernize-*`: `NULL`→`nullptr`, `typedef`→`using`, …). (c)
   Verify output against 0.9, same inputs through both codebases; don't advance a
   unit that hasn't been checked. The machine does the mechanical modernization; the
   human does the correctness verification.

## Invariants (breaking one is a bug)

These live here, not behind a link, because only this file loads every session.
Conventions (style preferences) are imported below from docs/CODING_RULES.md
(code) and docs/WRITING_RULES.md (prose).

- **Port, don't rewrite.** Carry 0.9's algorithm over unchanged. Changing *what*
  is computed is a rewrite, a separate opt-in track, not part of a port. Every PoC
  bug came from reimplementing a 0.9 algorithm instead of porting it. Following
  10.12's *design* is not a rewrite: shape and decomposition may track 10.12, but
  *what* is computed still matches 0.9.
- **The test suite and docs/TESTING_SPECIFICATION.md move together.** Adding, removing or
  changing what a test asserts is a change to both, in the same step. The specification is the
  description of the suite, so a test it does not describe is a test nobody can find, and a
  description with no test behind it is a claim we are not making good on. This is an invariant
  rather than a convention because nothing detects the drift: both files keep compiling, the suite
  keeps passing, and the record quietly stops being true.
- **`std::vector` is the default container.** Exceptions, each deliberate:
  fixed-size small blocks (e.g. 2×2 pivots) → `std::array`; non-owning view over a
  column/block → pointer + length (C++17, no `std::span`), one house convention;
  a `std::vector<std::vector<T>>` is a valid port target but **do not flatten** it
  into one buffer with offsets mid-port, that's a layout change (rewrite track).
- **No signed→unsigned index slips.** An index or offset that was `int` and can go
  negative must not become `size_t` arithmetic. Unsigned underflow on a descending
  loop is a bug this project has already hit once. Danger spots: descending loops,
  pointer offsets.
- **Contiguous storage to BLAS via `.data()`.** Never reintroduce the
  `(size==0) ? NULL : &v[0]` guard; `.data()` is well-defined for empty vectors.
- **Template definitions live in `.cpp`, never in headers.** Value-templated
  classes are declaration-only in the header; member bodies and the explicit
  instantiations (`template class Foo<double>;`) go in the `.cpp`. Moving a body
  into a header silently reintroduces per-translation-unit instantiation, the
  build cost the whole design exists to avoid, with no error to catch it. The
  `extern template` lines in the headers are the in-place *reminder* of this rule
  (they'd suppress the regression if it happened), but they don't enforce it, the
  rule does. Rationale: the explicit-instantiation entry in DESIGN_DECISIONS.
  *Deliberate exception:* a trivial accessor (e.g. `size()`) may be inlined in the
  header for performance, but only if the class stays explicitly instantiated with
  `extern template` present, which is exactly what stops that inline from being
  emitted per-TU. An accidental inline without that guard is the bug; a chosen
  inline under the guard is fine.

## Active design constraints

Decisions that shape code written now. Full rationale and history in docs/DESIGN_DECISIONS.md;
these are the always-on summary.

- **One `Val` template** for scalar type, `double` and `std::complex<double>`.
  Cholesky assumes Hermitian input; LDL^T variants assume complex-symmetric. No
  separate `*Real.h` / `*Complex.h` file pairs.
- **Explicit instantiation.** Headers declare; `.cpp` files define and instantiate
  for `double` and `std::complex<double>`. Adding a scalar type is one line per `.cpp`.
- **Namespaced headers** under `include/oblio/`, declarations only.
- **Flat `src/`**, all sources directly in `src/`, no per-category subdirectories.

## Coding rules (imported)

Loaded every session via import so conventions stay consistent. Edit the file, not
this line. (Path is relative to this file; CODING_RULES.md now lives in docs/, so
the import is `@docs/CODING_RULES.md`, a wrong path imports nothing, silently.)

@docs/CODING_RULES.md

## Writing rules (imported)

Loaded every session via import so prose conventions stay consistent across the
documentation. Edit the file, not this line. (Path is relative to this file;
WRITING_RULES.md lives in docs/, so the import is `@docs/WRITING_RULES.md`; a wrong
path imports nothing, silently.)

@docs/WRITING_RULES.md

## Build

macOS (alpamayo, Apple Silicon; Accelerate provides BLAS/LAPACK):

```
g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -I include \
    tests/<test>.cpp src/*.cpp -framework Accelerate -o <test>
```

Linux: replace `-framework Accelerate` with `-lblas -llapack -lm`.

## Tooling

- `.clang-format`, formatting. Run routinely; safe to apply whole-file.
- `.clang-tidy`, `modernize-*` and narrowing checks. Run as a *catch* per unit
  (report, no `--fix`) or as a deliberate sweep, not blind `--fix` mid-port, since
  it rewrites whole files and blurs one-unit-at-a-time diffs.
- Both are starting points, not yet calibrated to the tree; see the headers in each
  config before mass-applying.

**Where `clang-tidy` comes from.** Apple ships Clang but not `clang-tidy`, so it comes from
Homebrew's `llvm`. Only that one binary is exposed, by symlink:

```
brew install llvm
ln -s /opt/homebrew/opt/llvm/bin/clang-tidy /opt/homebrew/bin/clang-tidy
```

**Do not put `/opt/homebrew/opt/llvm/bin` on PATH instead.** That directory also holds `clang`
and `clang++`, which would shadow Apple's and silently change the compiler the whole tree builds
with, including which BLAS headers and which OpenMP spelling apply. The symlink takes the one
tool we lack and leaves the compiler alone. Homebrew's `clang++` remains reachable by full path
when a second toolchain is wanted deliberately.

### Editor navigation (CLion, terminal build kept)

The build stays in the terminal; CLion is a read-and-navigate surface only, driven by a
`compile_commands.json` so its clangd index sees the exact flags the Makefile uses. Written up one
step at a time as each is confirmed on alpamayo, so every step here is one that actually ran.

1. **Generate the compile database.** `brew install bear`, then from the repo root run the normal
   build prefixed with `bear --`, from clean so every translation unit is captured. On macOS the
   Makefile detects Darwin via `uname` and sets `BLAS_LIBS = -framework Accelerate` itself, so the
   command needs no BLAS argument: `make clean && bear -- make test`. This writes
   `compile_commands.json` (gitignored, machine-specific: it holds absolute paths) and changes
   nothing about the Makefile. Re-run only when the compile commands change, a new source file or a
   changed flag, not on every edit. Confirmed on alpamayo (macOS Tahoe, Apple Silicon): bear
   intercepted cleanly with no SIP trouble, and the database captured all 24 translation units (17
   in `src/`, 7 tests) each carrying `-std=c++17 -Iinclude -DOBLIO_BLAS_UNDERSCORE`, which is what
   clangd needs to index the code the way the build sees it. `Amd.cpp` and `Mmd.cpp` additionally
   carry `-w`; that silences the compiler but not clangd, so expect analysis noise on those two
   vendored/ported files regardless.

2. **Ignore the generated artifacts.** `compile_commands.json` holds absolute paths and the local
   toolchain, so it is per-machine and must not be committed; clangd also writes a `.cache/`
   index directory. Both are in `.gitignore` (`.idea/`, CLion's project folder, was already there).
   Confirmed on alpamayo: after the ignore rule was in place, `git status` did not list
   `compile_commands.json` among untracked files, so the database is present on disk but invisible
   to git.

3. **Open it in CLion as a compilation-database project.** File -> Open, and select the
   `compile_commands.json` file itself (not the folder); CLion offers "Open as Project", take it. It
   imports the database as the source of truth for flags and does not configure or run a build of
   its own, so the terminal stays the only place anything compiles. Confirmed on alpamayo: it
   reported "Compilation database project successfully imported", indexed without error, and
   cmd-click navigation resolves through the templates (jump-to-definition and find-usages both
   work). `reference/` needed no manual exclusion: it is gitignored, and CLion skips it on that
   basis, so cmd-click on a symbol with a `-0.9`/`-10.12` twin still jumps straight to the port
   with no pick-list. If a future setup does surface duplicate symbols, the fix is right-click the
   directory -> Mark Directory as -> Excluded, but it was not needed here.

**Why the build stays in the terminal.** A compilation-database project is deliberately
reduced-capability: the database records *what was compiled*, not *how to build*, so CLion can
compile a single file on demand (replaying one entry) but has no working full-build button. It knows
nothing of link steps, the `make test` orchestration, or the target structure, because that lives in
the Makefile, which CLion does not read. This is the wanted state, not a limitation: the terminal is
the only thing that compiles, and the two never fight over who owns the build. Even where CLion
builds *could* be wired in (registering `make` as an external tool, or migrating to CMake), the
terminal is the better choice on the merits here. The Makefile is the single source of truth and does
real work a second build description would have to duplicate and keep in sync: the `uname` platform
branch, the `-w` on the vendored orderers, the guarded explicit-instantiation object layout. And the
one-unit-at-a-time loop is a terminal loop, reading `PASS`/`FAIL` lines and comparing residuals
against the 0.9 oracle; the hand-rolled test prints are not a framework CLion's runner would parse
into a green/red tree, so routing the build through the IDE would cost the tidy output and gain
nothing. The clean model: CLion is a reading instrument pointed at the code, the terminal is the
workshop, and the compile database is the one-way bridge that lets the instrument understand what the
workshop produces without needing to run it. Wiring builds into CLion earns its keep only if the work
moves into CLion (its debugger, refactoring, inline run buttons); for navigation-only it is a moving
part with no payoff. The lightest middle option, if tab-switching ever grates, is binding `make test`
to a key as a CLion External Tool, which gives terminal-identical output in a panel without CLion
owning anything.

## Docs

- **CLAUDE.md** (this file), operating contract + doc index.
- **docs/PORTING_LEDGER.md**, per-unit porting status. Read first after a context gap.
- **docs/TODO.md**, work we intend to do that is not a porting question. The ledger asks whether a
  0.9 unit is carried over, which is a question with an end; TODO holds what outlives it, including
  things neither reference does. Straddling items live in TODO and the ledger points at them.
- **docs/TESTING_SPECIFICATION.md**, what the suite covers, stated independently of the test
  sources. Tests are written from it, not recovered from it, so an unsupported combination is a
  stated expectation rather than a missing test. Currently a catalog of the suite as it stands.
- **docs/CODING_RULES.md**, conventions a linter can't enforce. Imported above, so
  always loaded. Language-general.
- **docs/WRITING_RULES.md**, prose and documentation conventions (no em-dashes,
  ASCII only, American spelling, minimal formatting). Imported above, so always
  loaded. The prose counterpart to CODING_RULES.md.
- **docs/DESIGN_DECISIONS.md**, full rationale, history, dates, open questions. Read on
  demand; the code-shaping subset is summarized above under Active design
  constraints. When adding an entry, date it with **today's actual date read from
  context** (never copied from the entry above); the git commit date is the
  authoritative record.
- **`.clang-format` / `.clang-tidy`**, mechanical style + modernization,
  tool-enforced (not context). See Tooling above.
- **CONTRIBUTING.md / CHANGELOG.md**, stubs, ahead of need. Fill CONTRIBUTING on
  going public; start CHANGELOG at the first tagged release.
- **README.md**, public-facing overview.
- **archive/sparse_factorization.md**, the algorithm notes: elimination forest,
  symbolic factorization, supernodes and amalgamation, with worked examples. **Current
  guidance, despite the folder.** Actively maintained, and the code cross-references it
  by section number (`SymFactorEngine.h` cites 4.6). Its home is under review.
- **archive/**, otherwise history: `oblio-new-devlog.md`, `oblio_modernization_notes.md`,
  `oblio_modernization_appendix.md`, plus `template_comparison.jsx` and the 0.9-era
  `test09*.cc`. These describe 0.9/10.12 and the PoC; rationale, not current guidance.
- **experiments/**, runnable design studies, one self-contained folder each
  (`template-instantiation/`, `friend-access/`, `storage-options/`). Each validates a
  coding standard; build standalone with `make test`. Not part of the main build. See
  DESIGN_DECISIONS.
