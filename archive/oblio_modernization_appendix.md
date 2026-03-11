# Oblio Modernization — Appendix
## Detailed Notes on Each Recommendation

---

## A1. `Types.h` and `Utility.h` — Mechanical Modernization

### `typedef` → `using`

Current:
```cpp
typedef size_t Size;
typedef double LSize;
typedef double RVal;
```

Becomes:
```cpp
using Size  = size_t;
using LSize = double;
using RVal  = double;
```

Purely mechanical, no behavioural change. `using` is preferred in modern C++ because it is more readable, works with templates, and is consistent with the rest of the type alias syntax introduced in C++11.

### On the names `LSize` and `RVal`

- **`LSize`** — Large Size, a floating point type for counts or sizes that could overflow a `size_t`-based integer. Things like flop counts, fill estimates, or memory estimates for the factorization — numbers that can legitimately exceed 2^64 in a large sparse problem when computed as products. Using `double` gives approximate but overflow-safe arithmetic for those statistics.

- **`RVal`** — Real Value, the scalar type for matrix and factor entries. Distinguishes from the template parameter `Val` used in `Matrix<Val>` and `Factors<Val>`, which can be complex. `RVal` specifically means the real-valued scalar — used for norms, tolerances, pivot thresholds, and backward error computation, all of which must be real even when the matrix entries are complex. The distinction becomes essential once complex support is fully implemented: `std::abs` of a `std::complex<double>` returns `double`, so norm computations, pivot comparisons, and error bounds naturally live in `RVal` regardless of what `Val` is.

### `NULL` → `nullptr`

`NULL` is a macro that expands to an integer zero in most implementations, which means it can silently participate in integer arithmetic or overload resolution in ways `nullptr` cannot. `nullptr` has its own type (`std::nullptr_t`) and only converts to pointer types. In practice, most of the `NULL` occurrences dissolve entirely once `ResizeVector` is removed (see below), but anywhere a raw pointer is explicitly null-initialised, `NULL` becomes `nullptr`.

### `ResizeVector` → gone

`ResizeVector` currently does this:

```cpp
template<class T>
inline void ResizeVector(std::vector<T>* vec, Size size) {
    std::vector<T> tmpVec;
    tmpVec.reserve(size);
    tmpVec.resize(size);
    vec->swap(tmpVec);
}
```

The intent is to guarantee a fresh allocation with no excess capacity. That concern is real but the implementation is roundabout. The modern equivalent is:

```cpp
vec->assign(size, T{});  // resize and zero-initialise in place
// or, to genuinely release excess capacity:
*vec = std::vector<T>(size);
```

In practice, most `ResizeVector` calls in oblio are initialising a freshly constructed vector, so the capacity concern rarely applies. For those cases the call collapses into direct construction at the declaration site:

```cpp
// Before
std::vector<Size> prntVec;
ResizeVector<Size>(&prntVec, size);

// After
std::vector<Size> prntVec(size);
```

Which eliminates the follow-up null-guard pattern entirely:

```cpp
// Before — two lines after ResizeVector
Size* prntArr = (size == 0) ? NULL : &prntVec[0];

// After — one line, always safe
Size* prntArr = prntVec.data();
```

### The `cInfSize` double semicolon and infinity constants

A minor fix to catch while in `Types.h`:

```cpp
const Size cInfSize = std::numeric_limits<Size>::max();;  // double semicolon
```

Also, `cPosInfRVal` is currently computed as division by zero. The explicit modern form:
```cpp
const RVal cPosInfRVal = std::numeric_limits<RVal>::infinity();
const RVal cNegInfRVal = -std::numeric_limits<RVal>::infinity();
```

---

## A2. `BlasLapack.h` — Two Separate Problems

### Problem 1: the `#define UNDERSCORE` guard

BLAS and LAPACK are Fortran libraries. When a Fortran compiler builds them, it typically mangles symbol names by appending an underscore — so the Fortran routine `DGEMM` becomes the C-visible symbol `dgemm_`. This is platform and build-system dependent. Some builds (notably Windows, or certain MKL configurations) produce `dgemm` without the underscore, or even `DGEMM` in uppercase. The current code handles this with:

```cpp
#define UNDERSCORE

#ifdef UNDERSCORE
  extern "C" { void dgemm_(...); }
#else
  extern "C" { void dgemm(...); }
#endif
```

This is fragile in two ways. First, `#define UNDERSCORE` is hardcoded in the header itself — whoever builds oblio on a system where the convention differs has to edit the header. Second, a bare `#define UNDERSCORE` in a header is a landmine: any other code that uses `UNDERSCORE` as an identifier will silently break.

The standard fix is to move this to the build system. The compiler flag `-DOBLIO_BLAS_UNDERSCORE` gets passed at compile time, and the header just consumes it:

```cpp
#ifdef OBLIO_BLAS_UNDERSCORE
  extern "C" { void dgemm_(...); }
#else
  extern "C" { void dgemm(...); }
#endif
```

The `OBLIO_` prefix on the macro avoids the collision risk. Alternatively, switching to CBLAS eliminates the problem entirely — CBLAS has a stable, standardised C header (`cblas.h`) with no underscore ambiguity.

### Problem 2: `int` vs `size_t`

BLAS and LAPACK were written in Fortran 77, where all integers are 32-bit. Their C interfaces take `int*` for all dimension arguments. Oblio uses `size_t` (64-bit on any modern platform) for all its sizes. The current wrappers silently truncate, and at call sites in `FactorEngine.h`, oblio's `Size` values are passed directly as `int` arguments. On a 64-bit system, if any dimension exceeds `INT_MAX` (~2.1 billion), this silently wraps around.

The fix has two parts. First, a checked cast utility:

```cpp
inline int ToBlasInt(size_t n) {
    if (n > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("dimension exceeds BLAS int range");
    }
    return static_cast<int>(n);
}
```

Second, use it in the wrappers so they accept `size_t` natively:

```cpp
inline void GEMM(char transa, char transb,
                 size_t m, size_t n, size_t k, ...) {
    int im = ToBlasInt(m), in = ToBlasInt(n), ik = ToBlasInt(k);
    dgemm_(&transa, &transb, &im, &in, &ik, ...);
}
```

Now the wrappers match oblio's internal types, and the truncation is explicit, checked, and will throw rather than silently corrupt.

### The CBLAS alternative

If willing to depend on CBLAS (which ships with OpenBLAS, MKL, Accelerate, and most other BLAS distributions), the underscore problem and much of the boilerplate disappear. The complex equivalents (`cblas_zgemm` etc.) are in the same header, directly unblocking complex support work.

### The `BlasTraits` layer on top

Whichever raw binding approach is chosen, a `BlasTraits<Val>` struct sits on top and is what `FactorEngine` actually calls. The decision between raw BLAS and CBLAS is then contained entirely within `BlasLapack.h` — nothing above it changes either way.

```cpp
template<class Val> struct BlasTraits;

template<> struct BlasTraits<double> {
    static void gemm(...) { dgemm_(...); }
    static void potrf(...) { dpotrf_(...); }
};

template<> struct BlasTraits<std::complex<double>> {
    static void gemm(...) { zgemm_(...); }
    static void potrf(...) { zpotrf_(...); }
};
```

---

## A3. `Matrix.h` — Structure of a 1187-line File

### What is actually in the file

Six distinct responsibilities are mixed together:

**Storage and internal representation** — the private data members (`mRowIdxVecVec`, `mRowValVecVec`, `mDiagIdxVec`) and the column-major CSC-like layout they implement. This is the core.

**Construction and input validation** — `Set()`, `SetCol()`, `SetDiagVals()`, and the constructor that delegates to them. Contains the `bind2nd` usage and does bounds checking on row/column indices before delegating to the private helpers.

**Internal column-building helpers** — `rFormCols()`, `rSortCols()`, `rRemoveDupVals()`, `rCheckNoDupVals()`, `rFormDiag()`. These are the heaviest part of the file by line count. `rSortCols()` alone is a full two-pass bucket sort. Together this is 500+ lines of purely algorithmic code that has nothing to do with the public interface.

**Read accessors** — `GetVals()`, `GetCol()`, `GetDiagVals()`, `GetNormOne()`, `GetNormInf()`, `GetTranspose()`. These expose the stored data in various forms.

**Symmetry handling** — `SymmetrizeStrc()`, which merges the matrix with its transpose in a merge-sort style pass, and the `mIsStrclySym`/`mIsNmrclySym` flags.

**Printing** — `Print()`, `mPrecision`, `mFullPrint`. The `Print()` body is currently almost empty — it prints dimension counts but not values even in full-print mode, which looks like unfinished work.

### Recommended split

With explicit instantiation, the file splits naturally into a `.h` (declaration only, ~60–70 lines) and a `Matrix.cc` organised into clearly separated blocks:

```cpp
// ── Construction and validation ───────────────────────────────────────────
// Set(), SetCol(), SetDiagVals(). Contains the bind2nd → any_of replacement.

// ── Internal column builders ──────────────────────────────────────────────
// rFormCols(), rSortCols(), rRemoveDupVals(), rCheckNoDupVals(), rFormDiag().
// Pure algorithms, no public surface. The bucket sort in rSortCols is the
// most complex piece in the file.

// ── Read accessors ────────────────────────────────────────────────────────
// GetVals(), GetCol(), GetDiagVals(), GetNormOne(), GetNormInf(),
// GetTranspose(). All const.

// ── Symmetry ──────────────────────────────────────────────────────────────
// SymmetrizeStrc(). Self-contained merge pass.

// ── Printing ──────────────────────────────────────────────────────────────
// Print(). Currently minimal — good place to complete the full-print path.

// ── Explicit instantiations ───────────────────────────────────────────────
template class Matrix<double>;
template class Matrix<std::complex<double>>;
```

### The `bind2nd` replacement

Appears three times in `Set()` and `SetCol()`. Replaced with `std::any_of` and a lambda, which expresses the intent more directly:

```cpp
// Before
if (std::find_if(rowIdxVec.begin(), rowIdxVec.end(),
                 std::bind2nd(std::greater_equal<Size>(), numRows))
    != rowIdxVec.end()) {
    return eErrInvRowIdx;
}

// After
if (std::any_of(rowIdxVec.begin(), rowIdxVec.end(),
                [numRows](Size r){ return r >= numRows; })) {
    return eErrInvRowIdx;
}
```

### The `rcVal == 0.0` comparisons

For `double` this is technically correct — checking for an explicitly stored zero. For `std::complex<double>`, the comparison to `0.0` relies on an implicit conversion. The explicit form that works for both:

```cpp
if ((rcVal == Val{0}) && (r != c)) {
```

This is a check for the zero element of the scalar type, not a floating-point near-zero test.

---

## A4. The `(size == 0) ? NULL : &vec[0]` Pattern

### The basic case

```cpp
// Before
const Size* prntArr = (numSnodes == 0) ? NULL : &prntVec[0];

// After
const Size* prntArr = prntVec.data();
```

`std::vector::data()` has been guaranteed since C++11 to return a valid non-null pointer even for empty vectors — it just must not be dereferenced. That is exactly the semantics the ternary was trying to enforce manually.

### The form that disappears entirely

The most common pattern in `FactorEngine.h` and `ElmForestEngine.h`:

```cpp
std::vector<Size> prntVec;
ResizeVector<Size>(&prntVec, numSnodes);
Size* prntArr = (numSnodes == 0) ? NULL : &prntVec[0];
if (prntArr != NULL) {
    std::fill(&prntArr[0], &prntArr[numSnodes], cNullItm);
}
```

After applying both the `ResizeVector` and `NULL` fixes together:

```cpp
std::vector<Size> prntVec(numSnodes, cNullItm);
```

Four lines become one.

### The nested vector form

```cpp
// Before
std::vector<std::vector<Size>>* rowIdxVecArr =
    (numCols == 0) ? NULL : &(*rowIdxVecVec)[0];

// After
std::vector<Size>* rowIdxVecArr = rowIdxVecVec->data();
```

### The inner loop form

```cpp
// Before
for (Size c = 0; c < numCols; ++c) {
    Size cNumVals = rowIdxVecArr[c].size();
    const Size* cRowIdxArr = (cNumVals == 0) ? NULL : &rowIdxVecArr[c][0];
    const Val*  cRowValArr = (cNumVals == 0) ? NULL : &rowValVecArr[c][0];
    ...
}

// After
for (Size c = 0; c < numCols; ++c) {
    const Size* cRowIdxArr = rowIdxVecArr[c].data();
    const Val*  cRowValArr = rowValVecArr[c].data();
    ...
}
```

### The form where the raw pointer can go away entirely

Anywhere a raw pointer is introduced solely to be passed to `std::fill`, `std::copy`, or `std::find_if`, the iterator-based form of those algorithms removes the need for the pointer:

```cpp
// Before
Size* tmpPrntArr = (size == 0) ? NULL : &tmpPrntVec[0];
if (tmpPrntArr != NULL) {
    std::fill(&tmpPrntArr[0], &tmpPrntArr[size], cNullItm);
}

// After — raw pointer eliminated completely
std::fill(tmpPrntVec.begin(), tmpPrntVec.end(), cNullItm);
```

### BLAS call sites

```cpp
// Before
const RVal* sValArr = (sNumVals == 0) ? NULL : &sValVec[0];
GEMM('N', 'N', m, n, k, alpha, sValArr, lda, ...);

// After
GEMM('N', 'N', m, n, k, alpha, sValVec.data(), lda, ...);
```

Once the BLAS wrapper layer is updated to accept `size_t`, these call sites can be further simplified by passing the vector directly.

### Summary of forms and payoff

| Pattern | Replacement | Lines saved |
|---|---|---|
| `(n == 0) ? NULL : &vec[0]` | `vec.data()` | 0, but cleaner |
| `ResizeVector` + null guard + `fill` | `std::vector<T>(n, val)` | 3 → 1 |
| Raw pointer used only for `std::fill` | `vec.begin(), vec.end()` | pointer eliminated |
| Raw pointer used only for `std::copy` | `vec.begin(), vec.end()` | pointer eliminated |
| Nested vector null guard | `outerVec.data()` | same |
| BLAS call sites | `vec.data()` | 0, but cleaner |

Overall effect across the codebase: approximately 150–200 lines removed with no behavioural change anywhere.

---

## A5. Incomplete Features — What Is Actually There

### Static LDL — half done

`rSolveStaticLDL` in `SolveEngine.h` is **fully implemented** — all three passes of the LDL^T solve are present (L^{-1}, D^{-1}, L^{-T} sweeps over the supernodal factor). The solve side exists and works.

What is missing is the factorization. `SetFactorType` rejects `eStaticLDL` with `eErrInvType`, and `Factor()` hits `assert(false)`. The factorization for static LDL^T is structurally similar to Cholesky — instead of A = LL^T, compute A = LDL^T where D is block diagonal. For positive definite matrices with static pivoting, D ends up all 1×1 and the factorization degenerates to Cholesky. The interesting case is indefinite matrices where 2×2 pivots are needed — which is exactly what `mNum1x1Pivots` and `mNum2x2Pivots` in `FactorEngine` are tracking, confirming this was planned from the start.

**Decision required:** if implementing, write `rFactorStaticLDL` and wire into switch statements. If not, remove `rSolveStaticLDL` and `eStaticLDL` from the enum.

### Dynamic LDL — genuinely absent

No `rSolveDynamicLDL`, no `rFactorDynamicLDL`, no data structures. Only the `TODO` comments and the `eDynamicLDL` enum value exist. Dynamic LDL^T (Bunch-Kaufman) differs from static in that pivot selection happens at factorization time — the algorithm inspects numerical values and may swap rows and columns for stability. This requires a dynamic permutation that is not known at symbolic analysis time, making the symbolic factorization and supernodal structure interact with pivoting in a non-trivial way. Significantly more complex than static LDL^T, and what makes an indefinite solver practically useful for saddle-point problems and augmented systems.

**Decision required:** if out of scope, remove `eDynamicLDL` from `Types.h` entirely along with all TODO stubs.

### Multifrontal — partially prepared

`eMultifrontal` exists in the enum, `SetFactorAlgType` rejects it, and `Factor()` hits `assert(false)`. However, commented-out code in `ElmForestEngine.h` is directly relevant:

```cpp
//  if (mOptimizeForMultifrontal == true) {
//    bool success = rOptimizeForMultifrontal(...);
```

`ElmForestEngine` has a private member `mOptimizeForMultifrontal` set to `true` in the constructor but never used because the call is commented out. The supernodal structure preparation for multifrontal was at least partially worked on.

The distinction: left-looking and right-looking (both implemented) differ in traversal order. Multifrontal assembles a dense frontal matrix per supernode from both original values and child contributions, factors it with LAPACK, then stores only the update part. More memory-intensive but typically faster on modern hardware due to better cache utilization.

If multifrontal is in scope, the `rOptimizeForMultifrontal` path in `ElmForestEngine` should be uncommented and examined first.

### The ordering gap

`OrderEngine::Order` currently returns the identity permutation. Without a fill-reducing ordering, the factorization of any non-trivial sparse matrix produces catastrophically dense factors. The entire symbolic analysis pipeline runs fine with the identity permutation — it just produces a useless result in practice. AMD (Approximate Minimum Degree) is the natural addition: well-understood, compact reference implementations, and directly targets the minimum fill criterion that the elimination forest structure is already built around.

### Summary table

| Feature | Current state | Work to implement | Work to remove |
|---|---|---|---|
| Static LDL factorization | Solve done, factor missing | Write `rFactorStaticLDL`, wire into switches | Remove `rSolveStaticLDL`, `eStaticLDL` |
| Dynamic LDL | Entirely absent | Major — dynamic pivoting + symbolic interaction | Remove `eDynamicLDL` from enum |
| Multifrontal | Partially prepared in ElmForestEngine | Moderate — uncomment + implement `rFactorMultifrontal` | Remove `eMultifrontal`, clean up ElmForestEngine |
| Fill-reducing ordering | Identity permutation only | Moderate — implement AMD | Remove `OrderEngine::Order` body |

---

## A6. `Functional.h` — Dead Code

### What is in the file

Five functor classes: `Increment`, `ItmLess`, `ItmGreater`, `ValLess`, `ValGreater`. All pre-C++11 workarounds for the absence of lambdas.

### What is actually used

After searching every source file, test file, and header across the entire codebase: `ItmLess` appears exactly **once**, on line 621 of `Matrix.h` inside `SetCol()`:

```cpp
sort(rowVec.begin(), rowVec.end(), ItmLess<std::pair<Size, Val>>());
```

`ItmGreater`, `ValLess`, `ValGreater`, and `Increment` have **zero call sites**.

### The single live usage

`ItmLess` sorts a vector of `std::pair<Size, Val>` by the first element (the row index). This is exactly what `std::sort` with the default comparator does for pairs — `std::pair::operator<` already compares `first` before `second`. So `ItmLess` is replaceable with nothing at all:

```cpp
// Before
sort(rowVec.begin(), rowVec.end(), ItmLess<std::pair<Size, Val>>());

// After — default comparator on pair does the same thing
std::sort(rowVec.begin(), rowVec.end());
```

Or, to make intent explicit:

```cpp
std::sort(rowVec.begin(), rowVec.end(),
          [](const auto& a, const auto& b){ return a.first < b.first; });
```

Either way, `Functional.h` is entirely eliminated: one call site gets a trivial change, the file is deleted, and the `#include "Functional.h"` in `Matrix.h` goes with it.

---

## A7. Test Suite — Baseline State

### What each test covers

**`testVector.cc`** (38 lines, 1 test) — minimal smoke test. Constructs an empty `Vector<RVal>`, checks `GetInitErr`, `GetSize`, `GetPrecision`, `GetFullPrint`, and `Clear`.

**`testMatrix.cc`** (1625 lines) — the most thorough file. Covers trivial/empty construction, horizontal and vertical matrices, `GetCol`, `GetDiagVals`, `GetNormOne`, `GetNormInf`, `GetTranspose`, `SetCol`, `SetDiagVals`, `SymmetrizeStrc`, and the full COO input pipeline including zero-value removal, duplicate removal, and diagonal insertion.

**`testMultiplyEngine.cc`** (126 lines, 2 tests) — tests `MultiplyEngine::Multiply` on horizontal and vertical matrices against hand-computed results.

**`testElmForestEngine.cc`** (460 lines) — exercises the full pipeline short of `OblioEngine`: constructs symmetric matrices, runs `ElmForestEngine::ComputeElmForest`, `SymbolicEngine::ComputeSymbolic`, `FactorEngine::Factor`, and `SolveEngine::Solve`. Checks specific supernodal structure properties (supernode counts, tree counts, heights, front/update sizes, index sets) as well as backward error on the solve.

**`testOblioEngine.cc`** (224 lines, 1 active test) — exercises `OblioEngine::Solve` on a 9×9 numerically symmetric positive definite matrix, running both `eLftLooking` and `eRgtLooking`, asserting backward error below `1e-15` for both. A commented-out `Test2()` call in `main()` references a function that does not exist — planned but never written.

### Verified baseline: all five pass cleanly

All tests compile and pass under `-std=c++17 -Wall -Wextra` with GCC 13. The complete warning output across all five tests reduces to exactly three distinct warnings:

- **`bind2nd` deprecated** — all instances from the same three calls in `Matrix.h:Set()` and `Matrix.h:SetCol()`. Gone after the `Matrix.h` modernization pass.
- **`unused parameter 'argc'`** — in every `main(int argc, char** argv)`. Minor.
- **`unused parameter 'sIdxArr'`** in `FactorEngine::rAssembleUpdateVals` — a genuine unused parameter; the function body uses `tIdxArr` instead. Fix when touching `FactorEngine`.

### Compile commands

Three tests need only the headers:
```bash
g++ -std=c++17 -Wall -Wextra -I./include tests/testVector.cc -o testVector
g++ -std=c++17 -Wall -Wextra -I./include tests/testMatrix.cc -o testMatrix
g++ -std=c++17 -Wall -Wextra -I./include tests/testMultiplyEngine.cc -o testMultiplyEngine
```

Two need the `.cc` source files and BLAS/LAPACK:
```bash
LIBS="/usr/lib/x86_64-linux-gnu/blas/libblas.so.3 /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3"

g++ -std=c++17 -Wall -Wextra -I./include \
    tests/testElmForestEngine.cc \
    src/ElmForestEngine.cc src/SymbolicEngine.cc src/Symbolic.cc \
    $LIBS -o testElmForestEngine

g++ -std=c++17 -Wall -Wextra -I./include \
    tests/testOblioEngine.cc \
    src/ElmForestEngine.cc src/SymbolicEngine.cc src/Symbolic.cc \
    $LIBS -o testOblioEngine
```

### Gaps in coverage

- `complex<double>` is completely untested. Every test uses `RVal` (`double`). The explicit instantiation work for `complex<double>` adds a new code path with zero test coverage. At minimum, `testMatrix.cc` and `testOblioEngine.cc` should each get a parallel `complex<double>` test case before complex support work begins.
- `eRgtLooking` is tested in `testOblioEngine` but not in `testElmForestEngine`.
- The identity-permutation limitation in `OrderEngine` means backward error tests are testing an unordered factorization. Once AMD is added, existing tests will still pass (factorization of a positive definite matrix is correct regardless of ordering) but the fill and performance characteristics will change.

---

## A8. Comparison with Matchbox — Why This Feels Different

### What made matchbox amenable to standalone files

The matching algorithms in matchbox — Hopcroft-Karp, Micali-Vazirani, Gabow, Edmonds — share a common shape: they take a graph as input, run a self-contained algorithm, and return a matching. The data structures they operate on are either passed in or constructed locally. The algorithms do not call each other at runtime; they are alternatives, not collaborators. That shape made it natural to extract each algorithm into a `.cpp` file with a clean function signature at the top.

### Why oblio is structurally different

Oblio is a pipeline, not a menu. The five engines run in a fixed sequence and each one consumes the output of the previous:

```
Matrix → OrderEngine → ElmForestEngine → SymbolicEngine → FactorEngine → SolveEngine
```

Each stage produces an intermediate object with internal structure tightly matched to the algorithm that produced it. `Symbolic` stores the supernodal index sets, front sizes, and update sizes in exactly the layout that `FactorEngine` expects to traverse. `Factors` stores packed column-major supernode blocks in exactly the layout that `SolveEngine` expects to index into.

This coupling runs deeper than data passing. `FactorEngine` and `SolveEngine` access private members of `Matrix`, `Factors`, and `Symbolic` directly via `friend class` declarations — a deliberate performance choice. The alternative (accessor methods that copy data into caller-owned buffers) would add allocation and indirection into the innermost loops of the factorization. The `friend` access lets `FactorEngine` write directly into supernodal value blocks without going through a public interface.

### What "header-by-header refactor" means in practice

In matchbox, the unit of work was an algorithm. In oblio, it is a concern that spans a small cluster of related headers. Changes tend to ripple. Changing the storage layout of `Matrix` immediately requires checking `FactorEngine.h`, `ElmForestEngine.h`, `SymbolicEngine.h`, and anywhere else that accesses `Matrix` internals via `friend`. Changing `Factors<Val>` member names requires simultaneous updates to `SolveEngine` and `OblioEngine`.

The practical approach is to distinguish two categories of change:

**Mechanical modernization** — `typedef` → `using`, `NULL` → `nullptr`, `ResizeVector` → direct construction, `bind2nd` → lambda, `.data()` replacements. These are local and safe. Each can be done one header at a time with tests run between each change. They do not change any interfaces or data layouts.

**Structural changes** — explicit instantiation split (moving implementation to `.cc` files), splitting `Matrix.h`, reorganising the BLAS layer. These touch interfaces and require coordinated changes across multiple files. They should be done one class at a time with the full test suite run between.

The ordering matters. Mechanical passes first, structural changes after. Doing them simultaneously creates a large change surface where a failure is harder to localise.

### The explicit instantiation split is where it gets asymmetric

In matchbox, adding explicit instantiation was straightforward because each algorithm file was already self-contained. In oblio, most of the value-dependent code lives in `FactorEngine.h` — a 664-line header that is also the most algorithmically complex file in the codebase. Moving its implementation to `FactorEngine.cc` is the right thing to do, but because `FactorEngine` has `friend` access into `Matrix`, `Factors`, and `Symbolic`, the `.cc` file needs to include all of those headers. That is not a problem — it is just the correct dependency graph made explicit rather than implicit.

The result is not standalone algorithm files the way matchbox produced — `FactorEngine.cc` is still deeply coupled to the rest of the library. But it is a real improvement: the implementation is no longer recompiled in every translation unit that includes `FactorEngine.h`, the header is much smaller and faster to parse, and the coupling is now explicit in the `#include` list of the `.cc` file.

### The one place oblio could produce something matchbox-like

`OrderEngine` is the exception. A proper AMD implementation is algorithmically self-contained — it takes a sparsity pattern and returns a permutation. It does not need `friend` access into anything. It does not call BLAS. If implemented as a standalone function with its own `.cc` file, it would look structurally very similar to the matchbox algorithm files. AMD is also index-only and does not actually need `Val` at all — it is the one component of the value-dependent pipeline that would compile exactly once regardless of how many scalar types are supported.
