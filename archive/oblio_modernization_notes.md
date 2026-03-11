# Oblio Modernization — Observations & Recommendations

## Overall Structure

Oblio is a well-architected, header-heavy library. Almost all substantive logic lives in the headers as templated implementations (`Matrix`, `FactorEngine`, `ElmForestEngine`, `OrderEngine`, `SolveEngine` are all header-only), with only `ElmForestEngine.cc`, `Symbolic.cc`, `SymbolicEngine.cc`, and `Timer.cc` as separate translation units. This is a reasonable design for a templated library, but it means modernization touches are mostly in `.h` files.

The class breakdown maps cleanly to the sparse direct solver pipeline:

```
OrderEngine → ElmForestEngine → SymbolicEngine → FactorEngine → SolveEngine
```

...all orchestrated by `OblioEngine`. That layering is clean and should be preserved as-is during modernization.

---

## C++ Modernization

### Issues identical to matchbox's starting point

- `typedef size_t Size` and `typedef double RVal` — can become `using`
- `NULL` used extensively for pointer guards instead of `nullptr` — the `(size == 0) ? NULL : &vec[0]` pattern appears dozens of times
- `bind2nd` in `Matrix.h` (deprecated in C++11, removed in C++17)
- Manual copy-constructor/assignment-operator deletion via private declarations — should become `= delete`
- `ResizeVector` helper exists instead of using `vec.assign()` or direct construction — can be eliminated
- `std::fill` on raw pointer ranges where `std::fill(vec.begin(), vec.end(), ...)` would be cleaner
- Error codes returned as `Err` enum rather than exceptions or `std::expected` — this is a style choice, but it is at least consistent throughout

### New issues specific to Oblio

- **The `(size == 0) ? NULL : &vec[0]` raw-pointer extraction pattern** is pervasive in `FactorEngine` and `ElmForestEngine`. Modern C++ replaces this entirely with `.data()`, which handles the empty-vector case correctly and eliminates the ternary.

- **BLAS/LAPACK bindings** in `BlasLapack.h` use `#define UNDERSCORE` and C-style `extern "C"` with `int*` parameters. The `int` vs. `size_t` mismatch (BLAS uses `int`, Oblio uses `size_t`) should be made explicit with proper casts, and the `#define` guard is fragile.

- **`Functional.h`** contains hand-rolled `ItmLess`, `ItmGreater`, `ValLess`, `ValGreater` comparator functors — these predate C++11 lambdas and can all be replaced.

- **`friend class`** is used broadly (`FactorEngine`, `SolveEngine`, `ElmForestEngine` are all friends of data classes). This is fine for a library of this style, but it is worth noting as a design coupling point.

---

## Incomplete Features (Explicit TODOs)

This is more notable than in matchbox. Several unimplemented branches are stubbed with `TODO` and currently `assert(false)` or return `eErrInvType`:

| Feature | Location | Status |
|---|---|---|
| Static LDL | `FactorType` enum, `FactorEngine`, `OblioEngine` | Stub only |
| Dynamic LDL | Same + `SolveEngine` | Stub only — the indefinite case |
| Multifrontal factorization | `FactorAlgType::eMultifrontal` | Defined, unimplemented |
| Fill-reducing ordering | `OrderEngine::Order` | Returns identity permutation only |

The ordering gap in particular is significant — `OrderEngine::Order` currently just returns the identity permutation. A real fill-reducing ordering (AMD, nested dissection, etc.) is entirely absent and would be essential for practical performance on any non-trivial matrix.

These represent genuine feature gaps, not just style issues. **Modernization should clarify whether these are in scope before touching the stubs** — it would be wasteful to carefully modernize LDL stubs if the intent is to implement them, since you'd be modernizing twice. Conversely, if multifrontal is out of scope, the dead enum values and commented-out code should be cleaned up.

---

## Structural Recommendations

### 1. Start with `Types.h` and `Utility.h`
Same as matchbox: convert `typedef` to `using`, introduce `nullptr`, and remove `ResizeVector` in favor of `.data()` and direct construction.

### 2. `BlasLapack.h` needs its own treatment
The `#define UNDERSCORE` guard is fragile. A cleaner approach is to use CBLAS (which has a stable C interface without the underscore issue) or at minimum move to a `namespace`-scoped inline wrapper and make the `int`/`size_t` truncation explicit via `static_cast`.

### 3. `Matrix.h` is the largest single file (1187 lines)
It mixes storage, validation, printing, norm computation, and I/O. A good candidate for splitting or at minimum organizing into clearly demarcated internal sections.

### 4. The `(size == 0) ? NULL : &vec[0]` pattern
This should be the first mechanical pass across all headers — replace with `.data()`. It is noisy, hides intent, and `.data()` is well-defined for empty vectors since C++11. The payoff per line changed is high.

### 5. Decide on incomplete features first
See the table above. The decision on Static LDL, Dynamic LDL, and multifrontal should be made before starting, as it determines how much of `FactorEngine` and `OblioEngine` is live code vs. dead scaffolding.

### 6. `Functional.h`
Minor cleanup — all four hand-rolled comparator functors can be replaced with lambdas at their call sites or removed in favor of standard library comparators.

### 7. Preserve and run the test suite first
The tests in `tests/` (`testMatrix.cc`, `testOblioEngine.cc`, `testElmForestEngine.cc`, etc.) are a valuable correctness baseline. Worth confirming they compile and pass before touching anything, so regressions are immediately visible.

---

## Key Difference from Matchbox

In matchbox, the algorithms were self-contained enough that the modernized output became standalone `.cpp` files. Oblio is a more tightly coupled library with interdependent class hierarchies — modernizing it will likely stay in header/library form rather than becoming standalone files. The approach will feel more like a systematic header-by-header refactor than the algorithm-by-algorithm rewrite done for matchbox.

---

## Template Strategy

### The requirement

Complex support (`std::complex<double>`) is a hard requirement, and the possibility of other scalar types (`float`, `long double`) should remain open. Dropping the `Val` template parameter entirely is not an option.

### Why explicit instantiation still works

The scalar type set is small and the template surface area is narrower than it first appears. Looking at what actually uses `Val` across the codebase, it reduces to a short list:

| Class / Engine | Needs `Val`? |
|---|---|
| `Matrix<Val>` | Yes — stores numerical values |
| `Factors<Val>` | Yes — stores factor values |
| `Temporary<Val>` | Yes — update matrix scratch space |
| `Vector<Val>` | Yes — numerical vector |
| `FactorEngine` | Yes — `rFactorLftLookingCC`, `rFactorRgtLookingCC`, `rFactorCC`, `rUpdateCC` |
| `SolveEngine` | Yes — `rSolveCC`, `rSolveStaticLDL` |
| `MultiplyEngine` | Yes |
| `ElmForest` | No — index structure only |
| `Symbolic` | No — index structure only |
| `Permutation` | No — index structure only |
| `ElmForestEngine` | No — only touches sparsity pattern |
| `SymbolicEngine` | No — only touches sparsity pattern |
| `OrderEngine` | No — only touches sparsity pattern |

The index-only pipeline (`ElmForestEngine`, `SymbolicEngine`, `OrderEngine`) is only nominally templated today because it receives a `Matrix<Val>` argument, even though it never touches the values. Separating the structural interface from the value interface would remove these from the template machinery entirely — they would compile once, not once per scalar type.

### Recommended approach

**Explicit instantiation on the value-dependent classes only, combined with structural separation of the index pipeline.**

The template definitions stay in headers (preserving full generality for any future scalar type), but explicit instantiations are provided for `double` and `std::complex<double>` in dedicated `.cc` files. This gives fast builds for the common cases while keeping the door open. Adding `float` later means one additional line per `.cc` file — a minimal maintenance burden.

The explicit instantiation list is short and stable, covering roughly five `.cc` files:

```cpp
// e.g. at the bottom of FactorEngine.cc
template class FactorEngine<double>;
template class FactorEngine<std::complex<double>>;
```

### BLAS traits layer

Completing the BLAS bindings for complex is a prerequisite for this approach and should be done as an early milestone. The right structure is a traits specialization rather than the current `#define UNDERSCORE` flat header:

```cpp
template<class Val> struct BlasTraits;

template<> struct BlasTraits<double> {
    static void gemm(...) { dgemm_(...); }
    static void potrf(...) { dpotrf_(...); }
    // ...
};

template<> struct BlasTraits<std::complex<double>> {
    static void gemm(...) { zgemm_(...); }
    static void potrf(...) { zpotrf_(...); }
    // ...
};
```

All scalar-type dispatch lives in one place. Adding `float` or `long double` later means adding one more `BlasTraits` specialization rather than hunting through the codebase. The commented-out complex BLAS stubs already present in `BlasLapack.h` confirm this work was always planned — this is the right moment to complete it.
