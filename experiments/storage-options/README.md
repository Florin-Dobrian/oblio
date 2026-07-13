# Storage-Options Experiment

Can **one multiply function** serve both a CSC sparse matrix and a vector-of-vectors sparse
matrix? Reference / teaching only, **not** part of the main Oblio build.

## Why this exists

Oblio needs both layouts, and the reason is mutability. See the flat-vs-VV entry in
DESIGN_DECISIONS; in brief:

- **CSC (flat)** is right wherever the structure is written once into a size known in advance:
  the matrix, the symbolic factorization, the static numeric factor. One allocation,
  contiguous streaming, cheap raw blocks to hand to BLAS.
- **Vector of vectors** is right for **dynamic LDL**, where a delayed pivot passes columns up
  to an ancestor and that ancestor's front grows at runtime by an amount symbolic never
  predicted. The growth is local, one front grows while its siblings do not, which is what VV
  does cheaply and a flat buffer does not.

The engines reach into storage directly for speed (the friend-access decision), so the layout
is part of each algorithm's contract. Two layouts would therefore mean two copies of every
kernel, unless the algorithm can be written so that it does not know the layout at all.

## The idea

In CSC, `colPtr[j]` is an **index**. But a real pointer is one step away:
`&mRowIdx[mColPtr[j]]` is where column `j`'s row indices begin, and it is a plain
`const std::int32_t*`.

In VV, the pointer is already there: `mRowIdx[j].data()`. **The same type, from a different
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
places. The first is the honest baseline: **hand-written CSC**, which builds no pointer arrays
at all and walks `colPtr` directly, as one would if CSC were the only storage.

## Results

Apple M4, 200000 columns, ~3.2M nonzeros:

```
hand-written CSC (baseline)     1.362 ms
multiply(), CSC pointers        1.454 ms   1.07x
multiply(), VV pointers         1.499 ms   1.10x
multiply(), VV scattered        8.723 ms   6.41x
```

All results bit-identical.

**One algorithm does serve both.** That is the answer to the question, and it needs no
templating and no polymorphism.

**The interface costs almost nothing, and costs the same on both layouts.** CSC pointers 1.07x,
packed VV 1.10x: a three percent spread between two entirely different storage formats, through
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
