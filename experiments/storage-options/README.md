# Storage-Options Experiment

Can **one multiply** serve both a *static* sparse matrix (fixed structure, stored flat) and a
*dynamic* one (mutable structure, stored as a vector of vectors)? Reference / teaching only, **not**
part of the main Oblio build.

**A note on the names.** The classes are `SparseMatrixStatic` and `SparseMatrixDynamic`, not `Csc`
and `Vv`. The layouts *are* CSC and vector-of-vectors, but the layout is a **consequence** of
whether the structure can change, not the thing being chosen. Static and dynamic name the reason;
flat and VV name the bytes. It is also the vocabulary the solver itself uses (`NumFactorStatic`,
`NumFactorDynamic`), and one vocabulary is easier to hold than two.

Scope: the timed kernel does not mutate, both matrices are built once and read, so the four
timed rows measure what the *layout* costs. The mutation the names promise is exercised separately,
in `testMutation` and `testInvalidation`, which is where `setValues` and `setColumn` earn their
place.

## Why this exists

Oblio needs both layouts, and the reason is mutability. See the flat-vs-VV entry in
DESIGN_DECISIONS; in brief:

- **Static** (flat, CSC) is right wherever the structure is written once into a size known in
  advance: the matrix, the symbolic factorization, the static numeric factor. One allocation,
  contiguous streaming, cheap raw blocks to hand to BLAS.
- **Dynamic** (a vector of vectors) is right for **dynamic LDL**, where a delayed pivot passes
  columns up to an ancestor and that ancestor's front grows at runtime by an amount symbolic never
  predicted. The growth is local, one front grows while its siblings do not, which is what VV does
  cheaply and a flat buffer does not.

The engines reach into storage directly for speed, so the layout is part of each algorithm's
contract. Two layouts would therefore mean two copies of every kernel, unless the algorithm can be
written so that it does not know the layout at all.

## The idea

Each class exposes the same three per-column accessors, named for what they return, not for the
array behind them:

```
rowIdx(j)    where column j's row indices start   (const std::int32_t*)
val(j)       where column j's values start         (const double*)
colSize(j)   how many entries column j has         (std::size_t)
```

In the static matrix a real pointer is one step from `colPtr`: `mRowIdx.data() + mColPtr[j]`. In the
dynamic one it is already there: `mRowIdx[j].data()`. **The same type, from a different place.**
`colPtr` never appears in a signature; it is a static-internal detail.

These are *lookups*: O(1) addresses into storage the matrix already holds, nothing allocated,
nothing owned. A lookup is a fact about the layout, so it lives on the storage, exactly as
`blockPtr` lives on the numeric factor.

The multiply then reads a column directly, at the moment of use:

```cpp
template <class Matrix>
void multiply(const Matrix& A, const double* x, double* y) const;   // y = A*x
```

It is a template over the matrix. It calls `rowIdx(j)` / `val(j)` / `colSize(j)` and names no
member, buffer, or layout, so **one source serves both storages** and the compiler specializes it
per storage. There is deliberately **no extractor**: nothing materializes a pointer-array view, so
nothing is owned, nothing goes stale, and no consumer carries a `columnPointers` of its own to
restate. That per-consumer extractor was the repetition this design avoids, and its only advantage,
a single storage-blind compiled kernel, buys nothing here while costing speed, API, and an
invalidation hazard (see the bulk-vs-direct entry in DESIGN_DECISIONS).

This is the same access shape the numeric engine must use on the factor anyway: a dynamic factor
grows under it, so any pointer extracted up front would dangle. The matvec matches that shape rather
than inventing a second one.

## What is measured

The kernel is `y = A*x`, column by column, scattering into `y`. It is a *pure* multiply (BLAS's
`beta = 0`): it overwrites `y`, zeroing it once up front, because column-outer touches every `y[i]`
across several columns and cannot write each in one shot. That access pattern is exactly
right-looking factorization and every assembly step in the multifrontal method: stream a column's
(index, value) run and scatter into a dense target.

Four rows. The three `multiply()` rows are one templated source, monomorphized per storage. The
first is the honest baseline, **hand-written flat**: it reads the raw CSC arrays through the static
class's public `colPtr()` / `rowIdx()` / `val()` (no friendship, reading is public) and walks
`colPtr` directly, calling no per-column accessor, as one would if CSC were the only storage. If the
templated version matches it, reaching a column through the accessor costs nothing.

## Results

Apple M4, 200000 columns, ~3.2M nonzeros (representative run):

```
hand-written flat (baseline)    1.340 ms
multiply(), static (direct)     1.383 ms   1.03x
multiply(), dynamic, packed     1.438 ms   1.07x
multiply(), dynamic, scatter    9.945 ms   7.42x
```

All results bit-identical.

**One source does serve both**, monomorphized per storage, with the storage-specific part confined
to each matrix's own accessors. Adding a storage is one `template ...` line.

**The accessor costs almost nothing, and costs the same on both layouts.** Static 1.03x, packed
dynamic 1.07x: a few percent between two entirely different storage formats. The static row matches
the hand-written flat baseline within noise, so reaching a column through the accessor is close to
free, and free on *both* layouts rather than only on the one it was written for. On a slower machine
the same code came out at or below the baseline; a fast machine makes an abstraction look more
expensive, because the kernel it wraps got cheaper.

**What costs is the layout, and only when the columns are scattered.** The two VV rows are the *same
class holding the same content*, differing only in the order their inner vectors were allocated.
Built in column order the allocator hands back near-adjacent blocks and VV costs about 1.07x, the
extra pointer hop and nothing else. Built in shuffled order it costs several times that, a cache
miss at every column boundary.

That is the real lesson, and it is about the *layout*, not the accessor: **a flat buffer guarantees
consecutive columns are adjacent in memory; a vector of vectors only ever borrows that from the
allocator.**

## Two caveats on the scattered number

Both matter, and without them the number would be misread.

**The scattered case is constructed, not observed.** To defeat locality we shuffle the allocation
order *and* interleave a spacer vector between every column, an adversarial worst case engineered to
remove the allocator's help entirely. A real program does not usually try that hard to be slow. Read
the packed row as VV's **structural** cost and the scattered row as an **upper bound** on the
locality penalty, with reality between. It is also strongly hardware-dependent: the same code has
measured from roughly 6x to nearly 9x depending on cache size and prefetcher, so it is a ceiling
that moves with the machine, not a constant to quote.

**And this kernel is the harshest possible setting for a cache miss.** A sparse matvec does about
two flops per element loaded: pure streaming, so a miss has nothing to hide behind and shows at full
price. A numeric factorization front is the opposite, a dense block handed to BLAS level 3, doing
`O(n^3)` arithmetic on `O(n^2)` data, so the same miss amortizes over far more work. VV's penalty
inside dynamic factorization will therefore be much smaller than the scattered figure here, though
not zero.

Which is why the experiment is worth having: it isolates the locality cost by holding everything
else fixed (same class, same content, same source, same accessors) and puts a *ceiling* on it rather
than a guess.

## The asymmetry is the design

The two classes do **not** offer the same mutation API, and that is deliberate.

| | static | dynamic |
|---|---|---|
| `setValues(j, val)`, one column, same structure, new numbers | **yes**, cheap: overwrites in place | **yes**, cheap: overwrites in place |
| `setColumn(j, rowIdx, val)`, one column's structure | **absent by design** | **yes**, cheap: the column owns its buffer |
| restructure | build a new one | `setColumn` |

`setValues` is at **column granularity**, and its signature is identical on both classes:
`setValues(std::int32_t j, const std::vector<double>& val)`. Static overwrites the contiguous run
`mVal[colPtr[j]..]`; dynamic overwrites `mVal[j]`'s contents; both in place, so neither invalidates
a pointer. Setting every value is a loop over columns, exactly as reading every column is. Value
mutation is what the solver does most (a Newton step, a time step, refactorize), and both layouts do
it cheaply, so both offer it the same way.

**`setColumn` is absent from the static one, and its absence is the point.** Changing a column's
*structure* in a flat layout means shifting every later column: `O(nnz)`, not `O(column)`. An API
that *looks* cheap and is secretly linear in the whole matrix is a trap. Refusing to offer it is not
a limitation; it is telling the truth about the storage, and it puts the decision where it belongs:
the caller knows whether they are changing one column or rebuilding, and picks the object that suits.

**This also explains why there is no common base class.** A shared interface would force one of two
lies: `setColumn` on the flat matrix, pretending an `O(nnz)` shift is a column operation, or
`setColumn` on neither, crippling the dynamic one to match its sibling's weakness. So the asymmetry
is not a wart to tidy away. It is what the two storages *are*.

**The rule, in one line: an object offers what its storage makes cheap, and nothing else.**

## Pointer validity, which is the rule that will bite

Direct access reads a column at the moment of use, so a read holds no pointer across a mutation and
nothing goes stale *within* the multiply. But the underlying rule still governs anything that does
hold a pointer, and `testInvalidation` demonstrates it rather than asserting it:

```
setValues   pointer UNCHANGED   (buffer reused; contents overwritten in place)
setColumn   pointer MOVED       (buffer replaced; anything held from before now dangles)
```

**So the rule is exactly: structural mutation invalidates; value mutation does not.** It holds in
both storages.

**It is the rule the solver's dynamic factor will live by.** Delayed pivoting grows a front, which
reallocates its buffer, which would dangle any pointer taken into it up front. Direct access is the
remedy already built in: fetch a supernode's block pointer *at the moment of use* (one indirection,
and this experiment says that costs nothing) rather than extracting up front. The numeric engine is
direct for exactly this reason, and the matvec here is direct in the same shape.

**The one thing this experiment does not rehearse** is growth *during* the sweep: its structures
grow only between runs, not mid-loop. Worth knowing before writing dynamic LDL, not after.

## What this settles for Oblio

The storage question and the algorithm question are **separable**. The layout decides where a
column's pointer comes from, and nothing else. So:

- We can have both storages without duplicating a kernel. One templated source calls each storage's
  accessors and monomorphizes per storage; the arithmetic stays written against pointers, which is
  what BLAS wants anyway.
- The accessor is not a tax. Reaching a column through it costs nothing against hand-written CSC.
- The locality cost is real but bounded, and it depends on the kernel. It supports the hybrid in
  DESIGN_DECISIONS (use VV as the numeric engine's growable scratch, because the *algorithm* needs
  it, then flatten once into the persistent factor) but does not by itself prove it. The solve
  streams the factor, the cache-hostile shape this experiment measures, so flattening is likely
  worth it there. The factorization hands dense blocks to BLAS, where the penalty amortizes, so VV
  during elimination is likely affordable. Both claims want their own measurement against a real
  front, not this one.

## Build

```
make test
```

`-O3 -DNDEBUG`. The optimization level is not incidental: the claim is that reaching a column
through an accessor costs nothing against hand-written CSC, and at `-O0` that would measure nothing.
