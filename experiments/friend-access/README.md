# Friend-Access Example (mat-vec)

Establishes the access-pattern standard for Oblio's compute kernels: **public
operators as the convenient API, `friend` + direct storage access as the
performance path.** Reference / teaching only — **not** part of the main Oblio build.

## The standard it demonstrates

- Data classes (`Matrix`, `Vector`) expose a public, bounds-checked API
  (`operator()`, `operator[]`).
- Compute engines (`MultiplyEngine`) are declared `friend` and reach the contiguous
  storage (`mVals`) directly for hot loops.
- Both coexist by design: the API for readable / non-hot-path use, `friend`-direct
  for performance. `friend` is *more* encapsulated than public getters — it grants
  access to the named engine only, not the whole program.

All files are guarded-explicit style (declaration-only headers; bodies + explicit
instantiation in `.cpp`), the settled instantiation pattern. No BLAS yet — `friend`
access is precisely what enables a later handoff of the raw block to dgemv/dgemm.

## Two methods, one result

- `multiplyByApi` — element access through the public operators. Each element is a
  non-inlined, cross-translation-unit call carrying a bounds-check assert, and the
  compiler cannot vectorize across the calls.
- `multiplyDirectly` — `friend` access: fetch `A.mVals.data()` / `x.mVals.data()`
  once, then walk raw contiguous memory. No per-element calls, no asserts, and the
  inner loop vectorizes.

## Build & run

```
make test     # builds ./test_multiply_cpp, runs correctness + timing
make clean
```

## What it shows

- Both methods agree with each other and with hand-computed results (real + complex).
- On a 2000×2000 dense matrix, `multiplyDirectly` runs roughly 3× faster than
  `multiplyByApi` (indicative — varies by machine and build).
- Two effects drive the gap: (a) the API path's per-element bounds-check asserts,
  which a release build (`-DNDEBUG`) removes; and (b) non-inlined cross-TU operator
  calls that block vectorization, present regardless of asserts. The direct path
  avoids both.

## Not a verdict against the API

`multiplyByApi` is the right tool for readable, non-performance-critical code. The
standard is: use `friend` direct access (and later BLAS) in the hot paths, keep the
public API everywhere else.

## Related

`../../docs/DESIGN_DECISIONS.md` — the friend-access decision and the guarded-explicit
instantiation rationale.
