# Storage-Options Experiment

Can **one multiply function** serve both a *static* sparse matrix (fixed structure, stored flat)
and a *dynamic* one (mutable structure, stored as a vector of vectors)? Reference / teaching only, **not** part of the main Oblio build.

**A note on the names.** The classes are `SparseMatrixStatic` and `SparseMatrixDynamic`, not
`Csc` and `Vv`. The layouts *are* CSC and vector-of-vectors, but the layout is a **consequence** of
whether the structure can change, not the thing being chosen. Static and dynamic name the reason;
flat and VV name the bytes. It is also the vocabulary the solver itself uses (`NumFactorStatic`,
`NumFactorDynamic`), and one vocabulary is easier to hold than two.

Honesty about the gap: **nothing in this experiment mutates.** Both matrices are built once and
read. What is measured is what the *layout* costs, not what the mutation buys. A `setColumn` on the
dynamic one would make the name earn itself; it is not written yet.

## Why this exists

Oblio needs both layouts, and the reason is mutability. See the flat-vs-VV entry in
DESIGN_DECISIONS; in brief:

- **Static** (flat, CSC) is right wherever the structure is written once into a size known in advance:
  the matrix, the symbolic factorization, the static numeric factor. One allocation,
  contiguous streaming, cheap raw blocks to hand to BLAS.
- **Dynamic** (a vector of vectors) is right for **dynamic LDL**, where a delayed pivot passes columns up
  to an ancestor and that ancestor's front grows at runtime by an amount symbolic never
  predicted. The growth is local, one front grows while its siblings do not, which is what VV
  does cheaply and a flat buffer does not.

The engines reach into storage directly for speed (the friend-access decision), so the layout
is part of each algorithm's contract. Two layouts would therefore mean two copies of every
kernel, unless the algorithm can be written so that it does not know the layout at all.

## The idea

In the static matrix, `colPtr[j]` is an **index**. But a real pointer is one step away:
`&mRowIdx[mColPtr[j]]` is where column `j`'s row indices begin, and it is a plain
`const std::int32_t*`.

In the dynamic one, the pointer is already there: `mRowIdx[j].data()`. **The same type, from a different
place.**

So both classes can fill the same three arrays:

```
rowIdxPtr[j]   where column j's row indices start
valPtr[j]      where column j's values start
len[j]         how many entries column j has
```

and the multiply takes nothing but those:

```cpp
void multiply(std::size_t size,
              const std::int32_t* const* rowIdxPtr,
              const double* const*       valPtr,
              const std::size_t*         len,
              const double* x, double* y) const;
```

No matrix. No template parameter. No virtual call. **One function, compiled once.** It cannot
tell CSC from VV, because by the time it runs there is nothing left to tell apart.

The symbol table says so plainly. Compiling `MultiplyEngine.cpp` yields exactly one `multiply`
symbol, and its signature mentions neither matrix class:

```
T MultiplyEngine::multiply(unsigned long, int const* const*, double const* const*,
                           unsigned long const*, double const*, double*) const
```

The two storage classes appear only in the two `columnPointers` overloads, which are the *only*
code that knows a layout exists.

## What is measured

The kernel is `y = A*x`, column by column, scattering into `y`. That is the access pattern of
right-looking factorization and of every assembly step in the multifrontal method: stream a
column's (index, value) run and scatter into a dense target.

Four rows. The last three are the same compiled function, called with pointers from different
places. The first is the honest baseline: **hand-written flat**, which builds no pointer arrays
at all and walks `colPtr` directly, as one would if CSC were the only storage.

## Results

Apple M4, 200000 columns, ~3.2M nonzeros:

```
hand-written flat (baseline)    1.362 ms
multiply(), static pointers     1.454 ms   1.07x
multiply(), dynamic, packed     1.499 ms   1.10x
multiply(), dynamic, scattered  8.723 ms   6.41x
```

All results bit-identical.

**One algorithm does serve both.** That is the answer to the question, and it needs no
templating and no polymorphism.

**The interface costs almost nothing, and costs the same on both layouts.** static 1.07x,
packed dynamic 1.10x: a three percent spread between two entirely different storage formats, through
one compiled function. The generality is not a tax paid for flexibility; it is close to free,
and it is free on *both* layouts rather than only on the one it was written for.

The residual 7% is worth naming honestly rather than rounding away. The baseline reloads
`colPtr[j+1]` each column while the general version already holds `len[j]`, so the two do very
nearly the same work; what the general version adds is a second array to stream. At this speed
there is little latency left to hide it behind. Measured on a slower machine the same code came
out at 0.95x, *faster* than the baseline. Both numbers are true: a fast machine makes an
abstraction look more expensive, because the kernel it is wrapped around got cheaper.

**What costs is the layout, and only when the columns are scattered.** The two VV rows are the
*same class holding the same content*, differing only in the order their inner vectors were
allocated. Built in column order, the allocator hands back blocks that are essentially
back-to-back, and VV costs 1.10x, the extra pointer hop and nothing else. Built in shuffled
order, it costs **6.41x**, and the gap is a cache miss at every column boundary.

That is the real lesson, and it is a lesson about the *layout*, not about the interface: **a
flat buffer guarantees that consecutive columns are adjacent in memory; a vector of vectors
only ever borrows that property from the allocator.**

## Two caveats on the scattered number

Both matter, and without them the number would be misread.

**The scattered case is constructed, not observed.** To defeat locality we shuffle the
allocation order *and* interleave a spacer vector between every column. That is an adversarial
worst case, engineered to remove the allocator's help entirely. A real program does not usually
try that hard to be slow. Read 1.10x as VV's **structural** cost, 6.41x as an **upper bound** on
the locality penalty, and reality as somewhere between.

It is also hardware-dependent, and strongly so: the same code measured 8.87x on a machine with
smaller caches and a weaker prefetcher. The M4 recovers a good deal of the loss. So the
scattered figure is a ceiling that moves with the machine, not a constant to quote.

**And this kernel is the harshest possible setting for a cache miss.** A sparse matvec does
about two flops per element loaded: it is pure streaming, so a miss has nothing to hide behind
and shows up at full price. A numeric factorization front is the opposite: a dense block handed
to BLAS level 3, doing `O(n^3)` arithmetic on `O(n^2)` data, so the same miss amortizes over far
more work. VV's penalty inside dynamic factorization will therefore be much smaller than 6.41x,
though it will not be zero.

Which is exactly why the experiment is worth having: it isolates the locality cost by holding
everything else fixed (same class, same content, same interface, same compiled function), and
it puts a *ceiling* on it rather than a guess.

## The asymmetry is the design

The two classes do **not** offer the same API, and that is deliberate.

| | static | dynamic |
|---|---|---|
| `setValues`, same structure, new numbers | **yes**, cheap: nothing moves | **yes**, cheap: nothing moves |
| `setColumn`, one column's structure | **absent by design** | **yes**, cheap: the column owns its buffer |
| restructure | build a new one | `setColumn` |

**`setValues` is the mutation the solver actually does most often**, a Newton iteration, a time
step, the same pattern with new numbers, refactorize, and the flat layout is perfectly happy with
it. Both classes have it.

**`setColumn` is absent from the static one, and its absence is the point.** Changing a column's
*structure* in a flat layout means shifting every later column: `O(nnz)`, not `O(column)`. An API
that *looks* cheap and is secretly linear in the whole matrix is a trap, and the caller who writes
it in a loop will not find out until their program crawls.

**Refusing to offer it is not a limitation. It is telling the truth about the storage.** And it
puts the decision where it belongs: the caller knows whether they are changing one column or
rebuilding, and can pick the object that suits. Want a column swapped for you? Use the dynamic one.
Want to shift data around a flat buffer? Do it yourself; this class will not pretend it is cheap.

**This also explains why there is no common base class.** A shared interface would force one of two
lies:

- `setColumn` on the flat matrix, pretending an `O(nnz)` shift is a column operation, or
- `setColumn` on **neither**, crippling the dynamic one to match its sibling's weakness.

So the asymmetry is not a wart to be tidied away. It is what the two storages *are*.

**The rule, in one line: an object offers what its storage makes cheap, and nothing else.**

## Pointer validity, which is the rule that will bite

The engine works through extracted pointers. So: **when do they stop being valid?**

```
setValues   does NOT invalidate.  The buffer stays put; only its contents change.
setColumn   DOES invalidate.      The column's buffer is replaced; anything into it dangles.
```

The experiment demonstrates this rather than asserting it (`testInvalidation`): extract the
pointers, mutate, and observe which pointers still point where they did.

```
setValues   pointer UNCHANGED   (buffer reused; contents overwritten in place)
setColumn   pointer MOVED       (buffer replaced; anything held from before now dangles)
```

**So the rule is exactly: structural mutation invalidates; value mutation does not.** It holds in
both storages, and it is not a quirk of this experiment.

**It is the rule the solver's dynamic factor will live by.** Delayed pivoting grows a front, which
reallocates its buffer, which dangles every pointer previously taken into it:

```cpp
eng.blockPointers(f, block);        // extracted once
for (kk) {
    ... factor kk ...
    f.mVal[pp].resize(bigger);      // a delayed column grows ancestor pp
    ...                             // block[pp] is now DANGLING, and silently so
}
```

Nothing in C++ enforces this, and nothing can. The remedy is one of two: fetch a supernode's block
pointer *at the moment of use* rather than up front (one indirection, and this experiment says that
costs nothing), or re-extract after any growth. The first is simpler and is what the numeric engine
should do.

**And that is the one thing this experiment does not rehearse**, because its structures do not grow
*during* the algorithm, only between runs of it. Worth knowing before writing dynamic LDL, not
after.

## What this settles for Oblio

The storage question and the algorithm question are **separable**. The layout decides where the
pointers come from, and nothing else. So:

- We can have both storages without duplicating a single kernel. The kernels stay written
  against pointers, which is what BLAS wants anyway.
- The abstraction is not a tax. Reaching a column through a pointer array costs nothing against
  hand-written CSC.
- The locality cost is real but bounded, and it depends on the kernel. It supports the hybrid
  in DESIGN_DECISIONS (use VV as the numeric engine's growable scratch, because the *algorithm*
  needs it, then flatten once into the persistent factor) but it does not by itself prove it.
  The solve phase streams the factor, which is the cache-hostile shape this experiment
  measures, so flattening is likely worth it there. The factorization phase hands dense blocks
  to BLAS, where the penalty amortizes, so VV during elimination is likely affordable. Both
  those claims want their own measurement against a real front, not this one.

## Build

```
make test
```

`-O3 -DNDEBUG`. The optimization level is not incidental: the claim is that reaching a column
through a pointer array costs nothing against hand-written CSC, and at `-O0` that would measure
nothing.
