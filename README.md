# Oblio

A sparse direct solver library for symmetric matrices. Implements supernodal
Cholesky and LDL^T factorization with multiple fill-reducing orderings and
traversal strategies.

## Features

- **Three factorization types**: Cholesky, Static LDL^T (with diagonal perturbation),
  Dynamic LDL^T (Bunch-Kaufman pivoting with 1×1 and 2×2 pivots)
- **Three traversal algorithms**: Left-looking, Right-looking, Multifrontal
- **Ordering**: Natural (identity), MMD (Multiple Minimum Degree), AMD (Approximate Minimum Degree, Davis/Amestoy/Duff)
- **Single and multiple RHS**: `Vector<Val>` for one RHS, `DenseMatrix<Val>` for
  batched solves using per-supernode BLAS (`dtrsm` + `dgemm`)
- **Scalar types**: `double`, `std::complex<double>` (via explicit instantiation).
  Cholesky requires Hermitian input for complex; LDL^T variants handle complex
  symmetric (non-Hermitian) matrices.
- **C++17**, header-declaration / `.cpp`-definition pattern for fast builds

## Quick Start

```cpp
#include "oblio/OblioEngine.h"
#include "oblio/Vector.h"
#include "oblio/DenseMatrix.h"
using namespace Oblio;

// Build matrix from COO (lower triangle, 0-based).
auto A = Matrix<double>::fromCOO(n, rows, cols, vals);

// Configure solver.
OblioEngine<double> eng;
eng.setOrderAlg(OrderAlg::eMMD);
eng.setFactorAlg(FactorAlg::eMultifrontal);
eng.setFactorType(FactorType::eCholesky);
eng.analyzeAndFactor(A);

// Single RHS.
Vector<double> b(n, 1.0), x;
eng.solve(b, x);

// Multiple RHS, B is n×nRHS column-major, X is filled on output.
DenseMatrix<double> B(n, nRHS), X;
eng.solve(B, X);
```

## Build

The Makefile compiles the units in `src/` once and links each test or example against them.

```bash
make            # build everything (tests and examples)
make test       # build and run the test suites
make tests      # build the test binaries only
make examples   # build the example programs
make clean
```

To compile a single unit by hand (all sources are flat in `src/`, so `src/*.cpp` catches
everything):

```bash
# macOS (Accelerate provides BLAS/LAPACK):
g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude \
    tests/test_solve.cpp src/*.cpp -framework Accelerate -o test_solve

# Linux (system BLAS/LAPACK):
g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude \
    tests/test_solve.cpp src/*.cpp -lblas -llapack -lm -o test_solve
```

CMake is also supported: `cmake -B build && cmake --build build && ctest --test-dir build`.

## Structure

```
include/oblio/      , public headers (declarations only)
  Types.h           , enums (Factorization, Traversal), factorization predicates, typedefs
  SparseMatrix.h    , sparse symmetric matrix (CSC)
  Vector.h          , dense vector (one right-hand side)
  Permutation.h     , bidirectional index map (oldToNew / newToOld)
  OrderEngine.h     , fill-reducing ordering (Natural, MMD, AMD)
  ElmForest.h       , elimination forest and supernodes (data)
  ElmForestEngine.h , builds the elimination forest
  SymFactor.h       , symbolic factor: supernodal index structure (data)
  SymFactorEngine.h , computes the symbolic factorization
  NumFactorStatic.h , numeric factor, flat per-supernode storage (data)
  NumFactorDynamic.h, numeric factor, vector-of-vectors storage (data; stub)
  NumFactorEngine.h , computes the numeric factorization (Cholesky, static LDL)
  UpdateBlock.h     , temporary update block used during numeric factorization
  BlasLapack.h      , operation-named BLAS/LAPACK wrappers and custom kernels
  MultiplyEngine.h  , sparse matrix-vector product and residual
  SolveEngine.h     , triangular solves and right-hand-side permutation
src/                , method bodies + explicit instantiations (flat layout)
  Amd.cpp           , AMD ordering (SuiteSparse 3.3.4, Davis/Amestoy/Duff, BSD-3-clause)
  Mmd.cpp           , MMD ordering (Sparspak/Liu, via Oblio 0.9)
tests/              , test suites (111 tests)
  smoke.cpp                    5,  quick end-to-end sanity
  test_order.cpp              21,  ordering (Natural, MMD, AMD)
  test_permutation.cpp        11,  permutation maps
  test_forest.cpp             23,  elimination forest and supernodes
  test_symfactor.cpp          29,  symbolic factorization
  test_numfactor.cpp          14,  numeric factorization
  test_solve.cpp               8,  end-to-end solve, residual at machine precision
examples/           , usage examples
  pipeline.cpp      , the pipeline by hand, every factorization / traversal / ordering
  basic.cpp         , sketch of the planned OblioEngine facade (not yet compiling)
```

## History

This codebase is a C++17 modernization of Oblio 0.9, a sparse direct solver
written circa 2003–2005. The algorithmic core, symbolic factorization,
numeric factorization kernels, custom BLAS routines, solve engines, is
ported directly from the 0.9 source, not reimplemented. What is new is the
wrapping: `enum class` instead of bare enums, a single `Val` template parameter
instead of separate `*Real.h` / `*Complex.h` file pairs, `std::vector` storage
instead of hand-rolled arrays, explicit template instantiation for fast builds,
and a namespaced `include/oblio/` header layout.

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
(see `examples/pipeline.cpp`); `OblioEngine<Val>` will orchestrate them behind a single header,
as sketched in Quick Start.

## Status

Done:

- [x] MMD and AMD ordering (AMD from SuiteSparse 3.3.4; MMD via Oblio 0.9)
- [x] Supernodal symbolic factorization (elimination forest + symbolic factor, ported from 0.9)
- [x] Cholesky and static LDL, both LDL^T and LDL^H, left-looking and right-looking
- [x] Single-RHS triangular solve (`Vector`)
- [x] Complex arithmetic: Hermitian Cholesky, complex-symmetric LDL^T, complex-Hermitian LDL^H
- [x] Namespaced headers (`include/oblio/`), explicit instantiation throughout
- [x] Validated against Oblio 0.9 as oracle; end-to-end residual at machine precision
- [x] 111 tests across 7 phase suites

Not yet:

- [ ] Dynamic LDL, Bunch-Kaufman 1x1 / 2x2 pivots (`NumFactorDynamic` is a stub)
- [ ] Multifrontal traversal
- [ ] Multi-RHS solve (dense right-hand sides)
- [ ] `OblioEngine<Val>` top-level driver (the Quick Start API)
