# Friend-Access Example (mat-vec)

Establishes the access-pattern standard for Oblio's compute kernels: **public
operators as the convenient API, `friend` + direct storage access as the
performance path.** Reference / teaching only, **not** part of the main Oblio build.

## The standard it demonstrates

- Data classes (`Matrix`, `Vector`) expose a public, bounds-checked API
  (`operator()`, `operator[]`).
- Compute engines (`MultiplyEngine`) are declared `friend` and reach the contiguous
  storage (`mVals`) directly for hot loops.
- Both coexist by design: the API for readable / non-hot-path use, `friend`-direct
  for performance. `friend` is *more* encapsulated than public getters, it grants
  access to the named engine only, not the whole program.

All files are guarded-explicit style (declaration-only headers; bodies + explicit
instantiation in `.cpp`), the settled instantiation pattern. No BLAS yet, `friend`
access is precisely what enables a later handoff of the raw block to dgemv/dgemm.

## Three methods, one result

- `multiplyByApi`, element access through the public operators. Each element is a
  non-inlined, cross-translation-unit call, and the compiler cannot vectorize across
  the calls.
- `multiplyDirectly`, `friend` access: fetch `A.mVals.data()` / `x.mVals.data()`
  once, then walk raw contiguous memory. No per-element calls; the inner loop
  vectorizes. Hand-written fast path.
- `multiplyWithBlas`, same `friend` access to obtain the raw block, then hand it to
  BLAS gemv (`dgemv_`/`zgemv_`). This is what `friend` ultimately enables and the
  real solver's fast path. (The matrix here is row-major, so the call uses gemv with
  `TRANS='T'` on the column-major `A^T` view; the real solver stores dense blocks
  column-major and skips the transpose.)

## Build & run

```
make test     # builds ./test_multiply_cpp, runs correctness + timing
make clean
```

## What it shows

- All three methods agree with each other and with hand-computed results (real +
  complex, the complex path uses gemv `TRANS='T'`, i.e. no conjugation).
- On a 2000×2000 dense matrix, `multiplyDirectly` runs several times faster than
  `multiplyByApi`, measured ~6× on Apple Silicon (M4 / AppleClang, release), ~3× on
  x86 / g++. It varies with how well the machine vectorizes the direct path: the API
  path stays call-bound (~17 ms) on both, while the direct path is much faster on
  stronger SIMD (~3 ms on the M4), so the ratio grows on better hardware.
- Built with `-DNDEBUG` (release, asserts off), so the measured gap is **structural**,
  not assertion overhead: the API path's `operator()` lives in a different translation
  unit, so it's a non-inlined call per element that blocks vectorization, while the
  direct path walks contiguous memory and vectorizes. (Toggling asserts barely moves
  the gap, the cost is the calls, not the checks. Worth measuring, not assuming.)
- `multiplyWithBlas` vs `multiplyDirectly` depends heavily on the BLAS. With
  reference (Netlib) BLAS the two tie (reference gemv is an unoptimized loop). With a
  tuned BLAS the difference is large: on Apple Silicon (M4, Accelerate) `multiplyWithBlas`
  ran **~0.5 ms vs ~3.2 ms** for the hand loop, **~6× faster**, and ~36× over the API
  path. Mat-vec is bandwidth-bound *per core*, but Accelerate breaks past a single
  core's ceiling (multithreading across cores, plus better prefetch/blocking), which a
  single-threaded hand loop can't. So even for mat-vec, handing the raw block to BLAS
  wins decisively on real hardware.
- Takeaway for the standard: the hand-written direct loop is a good teaching baseline
  and fallback, but wherever a BLAS call exists, the `friend`-obtained raw block should
  go to BLAS, not a hand loop. The advantage only grows for the compute-bound O(n³)
  kernels (`dgemm`/`syrk`/`trsm`) the factorization actually leans on.

## Not a verdict against the API

`multiplyByApi` is the right tool for readable, non-performance-critical code. The
standard is: use `friend` direct access (and later BLAS) in the hot paths, keep the
public API everywhere else.

## Related

`../../docs/DESIGN_DECISIONS.md`, the friend-access decision and the guarded-explicit
instantiation rationale.
