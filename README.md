# Oblio

A sparse direct solver library for symmetric matrices, real and complex. Supernodal Cholesky and
`LDL^T` / `LDL^H` factorization, definite and indefinite, with fill-reducing ordering and three
traversals over one pipeline.

## Features

- **Five factorizations**: Cholesky (`A = CC^H`); static `LDL^T` and `LDL^H`, whose pivots are
  fixed by the symbolic structure and which perturb a pivot too small to divide by; and dynamic
  `LDL^T` and `LDL^H`, which choose pivots as they go, taking 1x1 and 2x2 blocks and delaying a
  column they cannot use. The `T`/`H` pair coincides over the reals and differs over the complex
  field, where one is symmetric and the other Hermitian
- **Three traversal algorithms**: Left-looking, Right-looking, Multifrontal
- **Two fill-reducing orderings**, both minimum degree: MMD, using the exact degree, and AMD,
  using an approximate degree bound. `Ordering::MMD` and `Ordering::AMD`
- **One right-hand side per solve**, a `Vector<Val>`, with the factorization reused across as
  many as wanted. Many right-hand sides at once, which is where the solve would become a level-3
  BLAS operation, is the one thing on the roadmap rather than in the library; see Status
- **Scalar types**: `double`, `std::complex<double>` (via explicit instantiation).
  Cholesky requires Hermitian input for complex; LDL^T variants handle complex
  symmetric (non-Hermitian) matrices.
- **C++17**, header-declaration / `.cpp`-definition pattern for fast builds

## Quick Start

```cpp
#include "oblio/DirectSolver.h"
using namespace Oblio;

// A symmetric matrix in CSC, both triangles stored.
const SparseMatrix<double> A(n, colPtr, rowIdx, val);

Vector<double> b(n), x(n);
for (std::size_t i = 0; i < n; ++i) b[i] = 1.0;

// Ordering, factorization, traversal: the pipeline order. Each also has a setter.
DirectSolver<double> solver(Ordering::MMD2, Factorization::Cholesky, Traversal::LeftLooking);

// The three phases have different lifetimes: analyze depends only on the pattern,
// factor on the values, solve on the right-hand side.
solver.analyze(A);
solver.factor(A);
solver.solve(b, x);

// A second right-hand side reuses the factorization.
solver.solve(b2, x2);

printf("residual %.3e\n", solver.relativeResidual(A, b, x));
```

Multiple right-hand sides (a dense `B`) are not wired yet; see Status.

## Examples

Eight standalone programs in `examples/`, built by `make examples` and run for exit status by
`make test`. Each is heavily commented and answers one question, in roughly this order:

- **`example_basic`**, what one solve looks like start to finish, through the facade.
- **`example_matrix`**, how to get *your* matrix in: CSC written out by hand, converted from a
  dense array, and assembled from coordinate triplets, with the conversion that sorts, merges
  duplicates and keeps the diagonal.
- **`example_analysis`**, what the ordering buys, in fill, supernode count and forest height,
  across every ordering method on an arrow and two grids. Nothing is factored.
- **`example_indefinite`**, what happens when the matrix is not positive definite: one pattern
  under three value sets, and which factorizations refuse, which perturb and which pivot.
- **`example_reuse`**, what can be kept between solves and what invalidates it: one analysis over
  three matrices, one factor over three right-hand sides, and what each setter discards.
- **`example_amalgamation`**, what merging supernodes costs and buys: fill against factorization
  time against solve time, and the ratio of solves to factorizations where the trade turns.
- **`example_pipeline_real`** and **`example_pipeline_complex`**, what the facade is doing
  underneath, with the engines wired by hand over every ordering, factorization and traversal.

**If you are deciding whether you want a direct solver at all**, read `example_analysis`: it shows
an ordering turning a dense factor into a sparse one on nine vertices. **If you are deciding
whether you want this one rather than a Cholesky**, read `example_indefinite`: one matrix, three
value sets, and the three factorizations refusing, perturbing and pivoting in turn.

## Prerequisites

The whole dependency list is a C++17 compiler, GNU make, and a BLAS and LAPACK. Nothing is vendored
that has to be fetched, nothing is downloaded at build time, and there is no package manager to
satisfy: `git clone` and `make` is the entire path. The one platform difference worth knowing before
starting is that **macOS already has the BLAS**, in Accelerate, and Linux does not.

### macOS

Apple Silicon or Intel.

- **Xcode Command Line Tools**: `xcode-select --install`. Provides `clang++`, `make`, and the
  Accelerate framework, which is the BLAS and LAPACK. There is nothing further to install to build
  and run the tests.
- **CMake**, only for the CMake build: `brew install cmake`. Homebrew itself, if it is not already
  present, installs with
  `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`.

Accelerate is selected automatically: the Makefile links `-framework Accelerate` when `uname` says
Darwin, and CMake finds it through `find_package(BLAS)`.

### Linux

Verified on Ubuntu 24.04 with GCC 13.3 and GNU Make 4.3, both builds from a clean container.

- **Compiler and make**: `sudo apt install build-essential`
- **BLAS and LAPACK**: `sudo apt install libblas-dev liblapack-dev`. The `-dev` packages are the
  ones that matter: a system carrying only the runtime `libblas3` and `liblapack3` links nothing,
  and the failure is `/usr/bin/ld: cannot find -llapack` rather than anything about a missing
  package.
- **CMake**, only for the CMake build: `sudo apt install cmake`

The Makefile links `-llapack -lblas` on anything that is not Darwin. Any BLAS implementation with
the standard Fortran symbols will do, so OpenBLAS or MKL can be substituted by overriding
`BLAS_LIBS` on the make command line.

On distributions other than Debian and Ubuntu the package names differ, `lapack-devel` and
`blas-devel` on Fedora and RHEL, and those have not been tested here.

### Optional, for development rather than for building

None of these is needed to compile the library, run the tests or run the examples.

- **`bear`**, to generate the `compile_commands.json` an editor indexes from: `brew install bear`
  or `sudo apt install bear`. See Editors below.
- **`clang-tidy`**, which Apple does not ship: `brew install llvm` and a symlink to that one
  binary. The reason not to put its directory on `PATH` is in `CLAUDE.md` under Tooling.
- **Instruments**, the profiler of record on Apple Silicon, which comes with Xcode proper rather
  than with the Command Line Tools. Used by `benchmarks/`; see `benchmarks/README.md`.
- **Python 3**, for the ordering prototypes in `experiments/ordering/` only.

## Build

Two builds cover the same sources, and they are not redundant:

- **Makefile**, canonical. What the development loop runs, and when the two disagree it is the one
  to believe.
- **CMake**, the release path. Probes for its dependencies rather than assuming them, and is what a
  consumed library is built and installed through.

Both compile the same `src/`, `tests/` and `examples/`, and both discover what they build by
wildcard, so a new `tests/*.cpp` file needs no edit in either.

### Makefile

```bash
make            # build everything (tests and examples)
make test       # build and run the test suites, then run the examples for exit status
make tests      # build the test binaries only
make examples   # build the example programs
make clean
```

Direct and predictable: one compiler invocation per unit, no configure step, no generated
directory, and nothing to install beyond a compiler and a BLAS. It selects the BLAS by platform,
branching on `uname` to `-framework Accelerate` on macOS and `-llapack -lblas` elsewhere, and sets
`-DOBLIO_BLAS_UNDERSCORE` unconditionally. That is an assumption rather than a discovery, correct
on every platform targeted so far and wrong on a build whose Fortran symbols carry no trailing
underscore.

It builds in place: `.o` files land beside their sources in `src/`, and test binaries in the repo
root with a `_cpp` suffix, all of them gitignored. There is no library, only objects linked
directly into each executable. Nothing to remove but the files themselves, which is what
`make clean` does.

What it buys is the development loop. `make test` prints one `PASS`/`FAIL` line per assertion and
a count, which is what the port is read against, followed by one line per example: each is run and
its exit status checked, with its output discarded so it cannot drown the suites. That catches an
example that crashes or has stopped compiling, and nothing about whether its numbers are
right, and the Makefile's own header lists the
inner-loop targets beyond the ones above: building one suite, or compiling the library alone as
the fastest check that a unit still builds. The experiments under `experiments/` each carry their
own Makefile for the same reasons, and none of them is part of this build.

### CMake

```bash
cmake -B build && cmake --build build && ctest --test-dir build
```

Configures, builds a static `liboblio.a`, and registers each suite with `ctest`. Where the Makefile
assumes, CMake **detects**: `find_package(BLAS)` and `find_package(LAPACK)` locate the libraries,
and a `check_c_source_compiles` probe calling `dpotrf_` decides whether to define
`OBLIO_BLAS_UNDERSCORE`. On macOS that reports Accelerate and the underscore convention; on a
machine with split reference libraries it finds those instead, with no edit.

Two structural differences follow from that, beyond detection. It builds **out of source**, with
everything under `build/`, so the source tree stays clean and `rm -rf build` is a complete reset,
where the Makefile leaves `.o` files among the sources. And it produces a real **library target**,
`liboblio.a`, which each suite links against, where the Makefile links the objects into every
executable. Neither matters much for eight test binaries; both matter for anything installed.

It is also more portable than the Makefile, in the sense that matters for a consumer rather than
for a developer here. The Makefile defaults to `g++` and branches on `uname`, both overridable but
both written out by hand; CMake picks up whatever compiler and generator the host offers, including
Ninja, MSVC and Xcode, none of which this project has needed.

The cost is a configure step, a generated `build/` tree, and a second description of the same
build to keep honest. What it buys beyond detection is the future: install targets and an exported
config, so that a consumer can write `find_package(Oblio)`, are CMake's to give and not a
Makefile's. That is why CMake is expected to become the primary build when the library is consumed
from outside; see the 2026-07-31 entry in `docs/DESIGN_DECISIONS.md`.

### Compiling one unit by hand

Neither build is required. All sources are flat in `src/`, so `src/*.cpp` catches everything:

```bash
# macOS (Accelerate provides BLAS/LAPACK):
g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude \
    tests/test_solve.cpp src/*.cpp -framework Accelerate -o test_solve

# Linux (system BLAS/LAPACK):
g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude \
    tests/test_solve.cpp src/*.cpp -lblas -llapack -lm -o test_solve
```

### Which to use

- **Makefile, while working.** Faster to invoke, carries the single-unit targets, and is what the
  porting loop is built around.
- **CMake, to check the build still configures from nothing**, on a machine whose BLAS the
  Makefile has not hardcoded, and before anything is released.

The division is not a capability gap. CMake builds single pieces perfectly well,
`cmake --build build --target test_order` and `ctest -R test_order`. It is that Make is a better
task runner and a thinner layer, which is what a tight edit-compile-check loop wants: no configure
step in the way, the exact compiler line visible rather than a progress percentage, and arbitrary
shell where a target is really a task rather than an artifact. `make objs`, compiling the library
alone as a fast syntax check, is that shape, and so is the trace-diffing `make test` in
`experiments/ordering`. CMake is the better build *description*, which is what a released library
wants. Neither subsumes the other, which is why both are kept here and why a project without a
porting loop would reasonably keep only one.

Running both occasionally is worth the seconds, since they share no code and a difference between
them is information.

### Editors

Neither build needs an IDE, and the IDE in use does not build. CLion can attach to this tree two
ways, and the choice is made by what is opened, not by a setting:

- **Compilation-database project** (current). `bear -- make test` writes a `compile_commands.json`
  recording what the Makefile actually compiled, flags and all, and CLion indexes from that.
  Navigation and find-usages resolve through the templates. There is no build button, because the
  database records what was compiled rather than how to build it.
- **CMake project.** Opening `CMakeLists.txt` instead lets CLion drive CMake itself, configuring
  into its own `cmake-build-*/` directory, with working build and run buttons, per-target run
  configurations and a `ctest` tree.

The database arrangement is the current one, and the terminal stays the only thing that compiles;
the database is a one-way bridge letting the editor understand a build it never runs. It is
machine-specific and gitignored, and wants regenerating only when the compile commands change, a
new source file or a changed flag, not on every edit. The CMake arrangement is declined for now
because it would be a third configuration of the same sources alongside the two above, because the
Makefile is canonical and does work a second description would have to duplicate, and because the
suites print `PASS`/`FAIL` lines rather than a format an IDE runner would parse into a tree. It
becomes the better choice if the work moves into the IDE, for its debugger and refactoring.

PyCharm is opened separately on `experiments/ordering/`, for the Python prototypes there. The
rule governing both, one IDE per project root and never two at the same directory, is in
`docs/DESIGN_DECISIONS.md` (2026-07-31); the setup steps for each are in `CLAUDE.md` under Tooling.

## Structure

```
include/oblio/      , public headers (declarations only)
  Types.h           , enums (Factorization, Traversal), factorization predicates, typedefs
  SparseMatrix.h    , sparse symmetric matrix (CSC)
  Vector.h          , dense vector (one right-hand side)
  Permutation.h     , bidirectional index map (oldToNew / newToOld)
  OrderEngine.h     , fill-reducing ordering, MMD and AMD behind one enum
  QuotientGraph.h   , the representation Oblio's own orderings run on, and its degree buckets
  Mmd1.h  Mmd2.h    , Oblio's own MMD orderings
  Amd1.h  Amd2.h    , Oblio's own AMD orderings
  Amd1B.h Amd2B.h   , the same two orderings on a fused eliminator schedule (see below)
  ElmForest.h       , elimination forest and supernodes (data)
  ElmForestEngine.h , builds the elimination forest
  SymFactor.h       , symbolic factor: supernodal index structure (data)
  SymFactorEngine.h , computes the symbolic factorization
  NumFactorStatic.h , numeric factor, flat per-supernode storage (data)
  NumFactorDynamic.h, numeric factor, vector-of-vectors storage (data)
  NumFactorEngine.h , computes the numeric factorization (Cholesky, static and dynamic LDL)
  UpdateBlock.h     , one supernode's update to one ancestor, for the looking traversals
  UpdateMatrix.h    , one supernode's whole contribution block, for multifrontal
  BlasLapack.h      , operation-named BLAS/LAPACK wrappers and custom kernels
  MultiplyEngine.h  , sparse matrix-vector product and residual
  SolveEngine.h     , triangular solves and right-hand-side permutation
  DirectSolver.h    , the whole pipeline behind one object (analyze / factor / solve)
src/                , method bodies + explicit instantiations (flat layout)
                      One .cpp per header, plus the vendored orderings, which have none:
  Amd.cpp           , AMD ordering (SuiteSparse 3.3.4, Davis/Amestoy/Duff, BSD-3-clause)
  Mmd.cpp           , MMD ordering (Sparspak/Liu, via Oblio 0.9)
tests/              , test suites (252 assertions; see docs/TESTING_SPECIFICATION.md)
  smoke.cpp                    5,  quick end-to-end sanity
  test_order.cpp              77,  the eight non-trivial orderings, and each B pair against its
                                   original entry for entry
  test_permutation.cpp        11,  permutation maps
  test_forest.cpp             29,  elimination forest and supernodes
  test_symfactor.cpp          29,  symbolic factorization
  test_numfactor.cpp          18,  numeric factorization
  test_solve.cpp              14,  the solve step, residual at machine precision
  test_pipeline.cpp           69,  whole-pipeline combinations, by residual
examples/           , eight usage examples, named example_* as the tests are named test_*
                      (described under Examples above, not repeated here)
  example_basic.cpp            , example_matrix.cpp
  example_analysis.cpp         , example_indefinite.cpp
  example_reuse.cpp            , example_amalgamation.cpp
  example_pipeline_real.cpp    , example_pipeline_complex.cpp
benchmarks/         , timing against the current tree, and expected to keep compiling as it moves
  ordering/         , one phase against itself: what each ordering costs, in time and in fill
  pipeline/         , the phases against each other: what share of a solve the ordering is, and
                      after how many factorizations a slower-analyzing ordering pays for itself
experiments/        , frozen design studies, each answering one question with a measurement
  ordering/         , the minimum-degree family rebuilt one mechanism at a time, in C++ and Python
  storage-options/  , flat against vector-of-vectors, and the accessor that spans both
  template-instantiation/, three ways to instantiate a Val-templated class
  friend-access/    , the access pattern the numeric kernels use
  openmp/           , how much parallelism Accelerate already supplies, and what is left
```

The ordering enum also carries Natural, the identity, and our own minimum-degree implementations
under the names MMD1, MMD2, AMD1, AMD2, AMD1B and AMD2B. Those are work in progress; see
`experiments/ordering/`.

## History

This codebase is a C++17 modernization of Oblio 0.9, a sparse direct solver
developed between 1998 and 2005. The algorithmic core, symbolic factorization,
numeric factorization kernels, custom BLAS routines, solve engines, is
ported directly from the 0.9 source, not reimplemented. What is new is the
wrapping:

- `enum class` instead of bare enums
- a single `Val` template parameter instead of separate `*Real.h` / `*Complex.h` file pairs
- `std::vector` storage instead of hand-rolled arrays
- explicit template instantiation for fast builds
- a namespaced `include/oblio/` header layout

Early development attempted to rewrite some algorithms from first principles
(elimination tree, column counts, supernodal index sets). Every rewrite
introduced subtle bugs that only appeared on non-trivial matrices. The bugs
were found and fixed by compiling the original 0.9 code and using it as a
ground-truth oracle, running the same test matrices through both codebases
and comparing outputs entry-by-entry. In each case, the fix was to replace
the rewritten algorithm with a direct port of the 0.9 code. The lesson was
clear: the 0.9 algorithms were correct and well-tested over years of use;
reimplementing them from scratch added risk with no benefit.

The one genuine bug in 0.9, complex Cholesky using `zsyrk` (symmetric rank-k
update) instead of `zherk` (Hermitian rank-k update), was identified and
fixed in the modern code, along with the related conjugate-transpose issues
in the factor and solve paths.

## Pipeline

```
OrderEngine -> ElmForestEngine -> SymFactorEngine -> NumFactorEngine -> SolveEngine
```

`MultiplyEngine` supplies the sparse matvec and residual. Today the phases are wired by hand
(see `examples/example_pipeline_real.cpp` and its complex counterpart), or driven together by
`DirectSolver<Val>`, which owns the intermediates and exposes the analyze / factor / solve phases
(see `examples/example_basic.cpp`).

## Status

Done:

- [x] MMD and AMD ordering, vendored (AMD from SuiteSparse 3.3.4; MMD via Oblio 0.9)
- [x] Supernodal symbolic factorization (elimination forest + symbolic factor, ported from 0.9)
- [x] Cholesky and static LDL, both LDL^T and LDL^H, left-looking and right-looking
- [x] Single-RHS triangular solve (`Vector`)
- [x] Complex arithmetic: Hermitian Cholesky, complex-symmetric LDL^T, complex-Hermitian LDL^H
- [x] Namespaced headers (`include/oblio/`), explicit instantiation throughout
- [x] Validated against Oblio 0.9 as oracle; end-to-end residual at machine precision
- [x] `DirectSolver<Val>`, the top-level analyze / factor / solve driver
- [x] 252 assertions across 8 suites
- [x] Dynamic LDL, threshold 1x1 / 2x2 pivots: all three traversals, delayed columns and all, at
      machine precision. Non-root supernodes follow Ashcraft, Grimes and Lewis (1998) Figure 3.4
      with the Figure 3.3 acceptance test; roots, which cannot delay, use bounded Bunch-Kaufman
- [x] Multifrontal traversal, for every factorization and both factor storages
- [x] Complex LDL^H with dynamic pivoting, an extension rather than a port, since 0.9's complex
      LDL is symmetric only: no reference to check against, so the oracles are the residual and
      reconstruction of L D L^H on dense fronts

Not yet:

- [ ] Multi-RHS solve (dense right-hand sides)
