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
   mechanical idioms (`modernize-*`: `NULL` -> `nullptr`, `typedef` -> `using`, ...). (c)
   Verify output against 0.9, same inputs through both codebases; don't advance a
   unit that hasn't been checked. The machine does the mechanical modernization; the
   human does the correctness verification.
3. **Ask before touching a file.** Nothing is created, modified or deleted without an explicit
   go-ahead for that specific change. Reading code, running a measurement, answering a question
   and proposing an edit are all fine unprompted; making the edit is not. Thinking aloud about
   how a vendored routine works, or asking where a check belongs, is a question and not a
   request to change our code, and the difference is often invisible in the wording, so the
   default when it is unclear is to answer and stop. An unasked change costs more than it saves:
   it has to be reviewed, it buries the change that was actually wanted, and a wrong one leaves
   the tree worse than doing nothing.
4. **Present only the files touched in that step, and present all of them.** The output
   directory is cumulative, so earlier uncommitted work sits there and gets presented again by
   accident, which buries the files that actually changed and makes the review ambiguous. If
   older work is still uncommitted, say so in one line rather than presenting the file a second
   time. The other half matters as much: a file that is not presented cannot be collected, so
   presenting a representative sample of a step is useless. If a previous answer failed to
   present a file it touched, present exactly the missing files and say what they are
   correcting, without re-presenting files that were already presented correctly.
5. **Never discard uncommitted work to undo a mistake.** `git checkout` on a path, `git reset
   --hard`, `git clean` and `git stash drop` throw away everything uncommitted in what they
   touch, not just the change being undone, and a path argument that looks narrow may not be:
   `git checkout experiments/ordering/` takes the whole directory, README included. Undo by
   editing the file back, or by copying from the staged copy in the output directory, both of
   which are surgical and reversible. This is not a hypothetical: it cost a full session's
   documentation work once and a second file twenty minutes later, and the only reason the
   second was cheap is that the staged copy was current. Which is the other half of the rule:
   stage after every meaningful edit, not at the end, so a copy always exists to restore from.

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
  fixed-size small blocks (e.g. 2x2 pivots) -> `std::array`; non-owning view over a
  column/block -> pointer + length (C++17, no `std::span`), one house convention;
  a `std::vector<std::vector<T>>` is a valid port target but **do not flatten** it
  into one buffer with offsets mid-port, that's a layout change (rewrite track).
- **No signed -> unsigned index slips.** An index or offset that was `int` and can go
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

`g++` there is Apple Clang, not GCC. On macOS `/usr/bin/g++` and `/usr/bin/clang++` are one file
under two names, dispatching to the active developer directory; `g++ --version` says `Apple clang`.
The README's macOS prerequisites explain it. What matters for this work is that the two compilers
warn about different things under the same `-Wall -Wextra`, so **a warning-clean build here is not
evidence of one on Linux, or the reverse**. Apple Clang has `-Wbitwise-instead-of-logical` and GCC
has no such option, which cost a build break on 2026-08-04; the optimizers differ too, and a GCC
`-O3` loop-hoist once made a real speedup measure as zero. `CXX` uses `?=`, so `make CXX=clang++`
names it explicitly.

Note that `src/*.cpp` no longer picks up the vendored orderings: they live in `private/`, which is
gitignored. A by-hand command line like the one above therefore builds without them, and
`Ordering::MMD` and `Ordering::AMD` refuse. Add `private/*.cpp` and `-DOBLIO_VENDORED_ORDERINGS` to
include them, or use the Makefile, which detects the directory itself.

### Building the way everyone else does

`private/` holds the two vendored orderings and is not published, so this machine builds something
nobody else can. To build as they do, put `OBLIO_PUBLIC=1` in front of any make command:

```
make test                          # this machine: 283 assertions
OBLIO_PUBLIC=1 make test           # everyone else:  269
```

`make help` in any of those directories prints its target list and this note, so the reminder is a
command rather than a file to find.

**One word, every Makefile that links the library**, so there is one thing to remember rather than
one per directory:

```
OBLIO_PUBLIC=1 make test                    # repo root
OBLIO_PUBLIC=1 make test                    # experiments/ordering
OBLIO_PUBLIC=1 make all                     # benchmarks/ordering
OBLIO_PUBLIC=1 make all                     # benchmarks/pipeline
```

`export OBLIO_PUBLIC=1` once and a whole shell session builds that way. On any other machine the
variable changes nothing, `private/` being absent already. It is harmless on every target: `clean`
removes the same things either way, so there is nothing to remember about when it applies.

The two modes interleave with no clean in between: each directory keeps a `.build-mode` stamp
recording which it was built in, and a switch rebuilds while a repeat does not. Nothing is deleted
either way.

### Before a release: build a real clone

`OBLIO_PUBLIC=1` cannot catch two things, because `private/` is still on disk either way: a
published file referring to a path inside it, and a file that is needed but was never committed.
Only a clone proves absence, `git clone` carrying committed history and nothing else.

```
git clone . /tmp/oblio-clone
cd /tmp/oblio-clone
make test
```

Expect `269/269 assertions across 8 suites, 8 examples run`, against 283 in the working tree. Then

```
cd -
rm -rf /tmp/oblio-clone
```

Run 2026-08-04 on alpamayo, after making the vendored orderings private: 252 in the tree, 238 in
the clone, as expected. The counts have since grown to 283 and 269, with `MMD3`, `AMD3` and the amd
alignment work; the figures above are the current ones.

Worth extending to the other directories when something there has changed, all of which should
work in the clone with no `private/` present:

```
cd /tmp/oblio-clone/experiments/ordering && make test    # every layer agrees, 35 comparisons
cd /tmp/oblio-clone/benchmarks/ordering  && make all     # builds; MMD and AMD rows refuse
cd /tmp/oblio-clone/benchmarks/pipeline  && make all     # builds both drivers
cd /tmp/oblio-clone && cmake -S . -B bld && cmake --build bld && (cd bld && ctest)
```

The CMake configure line prints which mode it chose, so `private/ absent, MMD and AMD will refuse`
is itself a check.

Once before a release rather than per push: `make test` and `OBLIO_PUBLIC=1 make test` cover the
day to day.

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
builds *could* be wired in, by registering `make` as an external tool or migrating to CMake, the
terminal is the better choice on the merits here:

- The Makefile is the single source of truth and does real work a second build description would
  have to duplicate and keep in sync: the `uname` platform branch, the `-w` on the vendored
  orderers, the guarded explicit-instantiation object layout.
- The one-unit-at-a-time loop is a terminal loop, reading `PASS`/`FAIL` lines and comparing
  residuals against the 0.9 oracle. The hand-rolled test prints are not a framework CLion's runner
  would parse into a green/red tree, so routing the build through the IDE would cost the tidy
  output and gain nothing.

The clean model: CLion is a reading instrument pointed at the code, the terminal is the workshop,
and the compile database is the one-way bridge that lets the instrument understand what the workshop
produces without needing to run it. Wiring builds into CLion earns its keep only if the work moves
into CLion (its debugger, refactoring, inline run buttons); for navigation-only it is a moving part
with no payoff. The lightest middle option, if tab-switching ever grates, is binding `make test` to
a key as a CLion External Tool, which gives terminal-identical output in a panel without CLion
owning anything.

### Python prototypes in PyCharm (`experiments/ordering`)

The ordering experiment carries a Python twin per layer, read in PyCharm while the C++ is read in
CLion. One rule governs the pair: **one IDE per project root, and never two at the same directory.**
CLion's root is the repo root and PyCharm's is `experiments/ordering/`, so the two `.idea/`
directories never meet, and CLion already lists that folder under `excludeRoots`. Do not open the
repo root in PyCharm; it is CLion's. The reasoning, and what to do if Python spreads to other
experiments, is the 2026-07-31 entry in docs/DESIGN_DECISIONS.md.

Setup, confirmed on alpamayo:

1. File -> Open, select `experiments/ordering/`, trust the project. PyCharm writes
   `experiments/ordering/.idea/` and changes nothing else in the folder, and `git status` lists no
   new file, since the root `.gitignore`'s `.idea/` has no leading slash and so matches at any
   depth.
2. It then asks for an interpreter, offering to set up a uv environment or a custom one. **Take the
   custom one and select the existing `python3`**, the same one `which python3` resolves to in the
   shell that runs `make`.

No venv and no `pyproject.toml`, for three reasons. The layers are standalone stdlib-only scripts
with no cross-imports, so there are no dependencies for an environment to hold. The `Makefile` runs
them with a bare `python3`, so a venv would put `make test` and the editor on different
interpreters, which would surface as the twins disagreeing for a reason that is not the code. And
it would drop a `.venv/` that this tree's `.gitignore` does not cover.

Moving the root later is two steps, open the new folder and delete the old `.idea/`, and no source
file, Makefile rule or README reference is affected by it.

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
- **experiments/ordering/**, the minimum-degree ladder, and the largest of these by some way. Its
  own README is the durable record and holds both alignment ledgers. Beside it, `MMD3.md` and
  `AMD3.md` record the two alignments ITERATION BY ITERATION: what each step established,
  discovered and decided, the wrong turns at full length, and the defects each found in the layer
  below. Narrative, read once. They live here rather than under docs/ because they are the record
  of experiment work, not guidance for a reader of the library. The README's method section is where the
  alignment procedure itself is written down, `REPORT.md` the measurements.
