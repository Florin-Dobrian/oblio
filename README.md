# Oblio

A sparse direct solver library for symmetric matrices. Implements supernodal
Cholesky and LDL^T factorization with multiple fill-reducing orderings and
traversal strategies.

## Features

- **Three factorization types**: Cholesky, Static LDL^T (with diagonal perturbation),
  Dynamic LDL^T (threshold pivoting with 1x1 and 2x2 pivots)
- **Three traversal algorithms**: Left-looking, Right-looking, Multifrontal
- **Ordering**, nine methods: Natural (identity), the vendored MMD (Multiple Minimum Degree) and
  AMD (Approximate Minimum Degree, Davis/Amestoy/Duff), and Oblio's own, built over a shared
  quotient graph: MMD1 and AMD1, the base algorithms; MMD2 and AMD2, which add the mechanisms
  their vendored counterparts carry; and AMD1B and AMD2B, which are AMD1 and AMD2 computed on a
  different schedule and return the same permutations. Ours are not drop-in replacements for the
  vendored pair and order differently
- **Single and multiple RHS**: `Vector<Val>` for one RHS, `DenseMatrix<Val>` for
  batched solves using per-supernode BLAS (`dtrsm` + `dgemm`)
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
DirectSolver<double> solver(OrderMethod::AMD, Factorization::Cholesky, Traversal::LeftLooking);

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
make test       # build and run the test suites
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
a count, which is what the port is read against, and the Makefile's own header lists the
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
  OrderEngine.h     , fill-reducing ordering, the nine methods below behind one enum
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
tests/              , test suites (241 assertions; see docs/TESTING_SPECIFICATION.md)
  smoke.cpp                    5,  quick end-to-end sanity
  test_order.cpp              77,  the eight non-trivial orderings, and each B pair against its
                                   original entry for entry
  test_permutation.cpp        11,  permutation maps
  test_forest.cpp             29,  elimination forest and supernodes
  test_symfactor.cpp          29,  symbolic factorization
  test_numfactor.cpp          18,  numeric factorization
  test_solve.cpp              14,  the solve step, residual at machine precision
  test_pipeline.cpp           58,  whole-pipeline combinations, by residual
examples/           , usage examples
  pipeline.cpp      , the pipeline by hand, every factorization / traversal / ordering
  basic.cpp         , the same solve through the DirectSolver facade
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

**On the ordering names.** A trailing digit means a different ordering: MMD2 has mechanisms MMD1
lacks, so their permutations and their fill legitimately differ and both are correct. A trailing B
means the same ordering computed on a different schedule, so AMD1B must return exactly AMD1's
permutation and a difference is a defect in one of them. The two axes share one enum, which is why
they are spelled out here and in `include/oblio/OrderEngine.h`.

## History

This codebase is a C++17 modernization of Oblio 0.9, a sparse direct solver
written circa 2003-2005. The algorithmic core, symbolic factorization,
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
(see `examples/pipeline.cpp`), or driven together by `DirectSolver<Val>`, which owns the
intermediates and exposes the analyze / factor / solve phases (see `examples/basic.cpp`).

## Status

Done:

- [x] MMD and AMD ordering, vendored (AMD from SuiteSparse 3.3.4; MMD via Oblio 0.9)
- [x] Oblio's own minimum-degree orderings over a shared quotient graph: MMD1, MMD2, AMD1,
      AMD2, and the AMD1B and AMD2B schedule variants
- [x] Supernodal symbolic factorization (elimination forest + symbolic factor, ported from 0.9)
- [x] Cholesky and static LDL, both LDL^T and LDL^H, left-looking and right-looking
- [x] Single-RHS triangular solve (`Vector`)
- [x] Complex arithmetic: Hermitian Cholesky, complex-symmetric LDL^T, complex-Hermitian LDL^H
- [x] Namespaced headers (`include/oblio/`), explicit instantiation throughout
- [x] Validated against Oblio 0.9 as oracle; end-to-end residual at machine precision
- [x] `DirectSolver<Val>`, the top-level analyze / factor / solve driver
- [x] 241 assertions across 8 suites
- [x] Dynamic LDL, threshold 1x1 / 2x2 pivots: all three traversals, delayed columns and all, at
      machine precision. Non-root supernodes follow Ashcraft, Grimes and Lewis (1998) Figure 3.4
      with the Figure 3.3 acceptance test; roots, which cannot delay, use bounded Bunch-Kaufman
- [x] Multifrontal traversal, for every factorization and both factor storages
- [x] Complex LDL^H with dynamic pivoting, an extension rather than a port, since 0.9's complex
      LDL is symmetric only: no reference to check against, so the oracles are the residual and
      reconstruction of L D L^H on dense fronts

Not yet:

- [ ] Multi-RHS solve (dense right-hand sides)
