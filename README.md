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

```bash
# macOS (Accelerate provides BLAS/LAPACK):
g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -I include \
    tests/test_smoke_real.cpp src/*.cpp -framework Accelerate -o test_smoke_real

# Linux (system BLAS/LAPACK):
g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -I include \
    tests/test_smoke_real.cpp src/*.cpp -lblas -llapack -lm -o test_smoke_real
```

All source files live flat in `src/`, so `src/*.cpp` catches everything.
Replace the test file to build each of the four test suites.

## Structure

```
include/oblio/      , public headers (declarations only)
  Types.h           , primitive typedefs, enums
  Matrix.h          , sparse symmetric CSC (lower triangle)
  Vector.h          , dense vector
  DenseMatrix.h     , dense column-major matrix
  Permutation.h     , bidirectional index map
  Symbolic.h        , supernodal elimination structure (data)
  BlasLapack.h      , BLAS/LAPACK trait wrappers + custom kernels
  OrderEngine.h     , fill-reducing ordering (MMD, Natural)
  SymbolicEngine.h  , etree, column counts, supernode amalgamation
  FactorEngine.h    , 9 factorization combinations
  SolveEngine.h     , single and multi-RHS triangular solves
  OblioEngine.h     , top-level driver (the only header users need)
src/                , method bodies + explicit instantiations
  Mmd.cpp            , MMD ordering (Liu/Sparspak via Oblio 0.9)
  Amd.cpp            , AMD ordering (SuiteSparse 3.3.4, Davis/Amestoy/Duff, BSD-3-clause)
tests/              , test suite
  test_smoke_real.cpp           18 tests, real tridiagonal, quick sanity
  test_smoke_complex.cpp        18 tests, complex tridiagonal, quick sanity
  test_extended_real.cpp       123 tests, Laplacians, edge cases, all orderings
  test_extended_complex.cpp    102 tests, complex Laplacians, all orderings
examples/           , usage examples
```

## History

This codebase is a C++17 modernization of Oblio 0.9, a sparse direct solver
written circa 2003–2005. The algorithmic core, symbolic factorization,
numerical factorization kernels, custom BLAS routines, solve engines, is
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
OrderEngine → SymbolicEngine → FactorEngine → SolveEngine
```

All orchestrated by `OblioEngine<Val>`, the only header users need to include.

## Status

- [x] MMD ordering (fixed `xadj` conversion for general graphs)
- [x] AMD ordering (SuiteSparse AMD 3.3.4, merged into single C++ file)
- [x] Supernodal symbolic factorization (bottom-up row-merge, ported from 0.9)
- [x] Cholesky, Static LDL, Dynamic LDL (Bunch-Kaufman 1×1 and 2×2 pivots)
- [x] Left-looking, Right-looking, Multifrontal
- [x] Single-RHS triangular solve (`Vector`)
- [x] Multi-RHS triangular solve (`DenseMatrix`, batched BLAS per supernode)
- [x] Namespaced headers (`include/oblio/`), explicit instantiation throughout
- [x] Complex arithmetic (`std::complex<double>`, Hermitian Cholesky + symmetric LDL)
- [x] Extended test suite: 2D Laplacian, n=1, block diagonal, indefinite, perturbation
- [x] Validated against Oblio 0.9 reference implementation
- [x] Complex extended tests (Hermitian Cholesky + symmetric LDL on Laplacians, all orderings)
