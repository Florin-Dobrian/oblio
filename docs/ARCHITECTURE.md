# Architecture

How the pieces fit, as they are today. This is the *living* description: what the units are, what
they promise each other, and which of those promises are load bearing. It is meant to be edited
whenever the code changes, which is what separates it from its neighbors.

- **README** is the outside view: what Oblio does and how to build it.
- **CLAUDE.md** holds the invariants, the rules that must not be broken while working.
- **DESIGN_DECISIONS** is history: dated entries recording what was chosen and why, superseded
  rather than rewritten. It explains how we got here; this file explains where here is.
- **PORTING_LEDGER** tracks unit by unit what has been ported and checked.
- **experiments/** are frozen studies. Each answers one question with a measurement and is not
  maintained afterwards, so nothing here should depend on their contents staying current.

Sections are added as they earn their place. Each describes one structure or contract that a reader
would otherwise have to reconstruct by grepping.

## The storage-agnostic factor contract

The numeric factor comes in two storages. `NumFactorStatic` holds flat index and value buffers with
an offset per supernode; `NumFactorDynamic` holds one index vector and one value vector per
supernode, so that a delayed pivot can grow a front without moving its neighbors. Dynamic LDL
requires the second; everything else is cheaper in the first.

**There is deliberately no common base class**, so the thing that lets one algorithm serve both is
not a type. It is a set of accessors that both classes spell the same way, with the same signatures,
differing only in body. The contract is structural rather than nominal, which is exactly why it
needs writing down: no header declares it, and a reader would otherwise have to infer it from what
compiles.

Five of the shared accessors are the same field read on both classes. Only the two *lookups* differ,
and they differ in the way the storages differ: an offset into a flat buffer against a vector that
already holds the pointer.

| accessor | returns | static | dynamic |
|---|---|---|---|
| `size()` | `std::size_t`, the matrix order | `mSize` | `mSize` |
| `snodeSize()` | `std::size_t`, how many supernodes | `mSnodeSize` | `mSnodeSize` |
| `factorization()` | `Factorization`, which method produced this | `mFactorization` | `mFactorization` |
| `frontSize(jj)` | `std::size_t`, supernode jj's own columns | `mFrontSize[jj]` | `mFrontSize[jj]` |
| `updateSize(jj)` | `std::size_t`, its update rows | `mUpdateSize[jj]` | `mUpdateSize[jj]` |
| `nodeToSnode()` | `const std::vector<std::int32_t>&`, which supernode owns a node | `mNodeToSnode` | `mNodeToSnode` |
| `nodeIdx(jj)` | `const std::int32_t*`, where jj's node indices start | `mNodeIdx.data() + mSnodeNodeIdxPtr[jj]` | `mNodeIdx[jj].data()` |
| `val(jj)` | `const Val*`, where jj's value block starts | `mVal.data() + mSnodeValPtr[jj]` | `mVal[jj].data()` |

### Who reads what

The contract is not one set. Each consumer uses the part it needs, and the differences are
informative: the solve needs to know which factorization produced the factor, because that decides
whether it conjugates, while the traversals never ask, because the engine already knows what it is
computing. Conversely the traversals need the node-to-supernode map to find an ancestor, and the
solve never does. The column is the static traversals (`factorStaticLeftLooking` and its
right-looking twin), which are the pair templated on the factor; dynamic pivoting reaches further,
for the growth verbs below.

| | `SolveEngine` | the static traversals in `NumFactorEngine` |
|---|---|---|
| `size`, `snodeSize`, `frontSize(jj)`, `updateSize(jj)`, `nodeIdx(jj)`, `val(jj)` | read | read |
| `factorization()` | read | not used |
| `nodeToSnode()` | not used | read |
| non-const `nodeIdx(jj)`, `val(jj)`, `numPerturbations()` | not used | write, through friendship |

**Reading is public; writing is not.** The const overloads are public, so a consumer that only
reads the factor needs no friendship, which is why `SolveEngine` is a plain reader. The non-const
overloads that hand out a writable value block are private, reached only by `NumFactorEngine` as a
friend. Name lookup gathers every overload before access is checked, so a non-const object outside
the friend context selects the private overload and fails to compile; `std::as_const(nf).val(jj)`
is the idiom for reading through such an object.

**The dynamic factor's growth verbs are private and one-sided.** `extendIndex`, `swap` and
`shrinkEntry` exist only on `NumFactorDynamic`, because only its storage makes them cheap. That
asymmetry is the same rule `experiments/storage-options` arrived at for the matrix pair, where
`setColumn` is absent from the flat class by design: an object offers what its storage makes cheap,
and nothing else.

**This set will grow.** When dynamic LDL delays columns across a forest, the solve's diagonal pass
needs `pivotType` to apply 2x2 block solves, and the table above gains a row.

### The parallel with the matrix, and where it stops

The same shape appears one level down, in `SparseMatrixStatic` and `SparseMatrixDynamic` in
`experiments/storage-options`, which is where the pattern was measured before it was adopted here.
The correspondence is close enough that the two kernels read alike:

```cpp
const std::int32_t* nodeIdx = nf.nodeIdx(jj);   // a supernode's value block, in the factor
const Val*          val     = nf.val(jj);

const std::int32_t* rowIdx  = A.rowIdx(j);      // a column, in the matrix
const double*       val     = A.val(j);
```

Two places it stops. The matrix's `colSize(j)` has no single counterpart, because a column is a run
and a supernode's value block is a rectangle: its extent takes two numbers, `frontSize(jj)` for the
columns and `frontSize(jj) + updateSize(jj)` for the rows, which is also the leading dimension.
That pair is what lets the numeric kernels take `(Val* block, rows, cols, ld)` and never see either
factor class. And `factorization()` has no counterpart at all, because a matrix carries no method
identity while a factor must say whether the solve conjugates.

### Reading the names

An offset array keeps the `Ptr` suffix (`colPtr`, `snodeNodeIdxPtr`, `snodeValPtr`); a lookup does
not (`rowIdx(j)`, `val(j)`, `nodeIdx(jj)`, `val(jj)`). Where a class holds both the whole array and
the per-entity slice, they are overloads of one name: `nodeIdx()` returns the flat buffer and
`nodeIdx(jj)` the address of one supernode's slice within it. See the 2026-07-19 and 2026-07-17
entries in DESIGN_DECISIONS for how those names were settled.
