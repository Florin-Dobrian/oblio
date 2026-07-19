# Ordering Experiment

`OrderEngine` calls a vendored AMD and a vendored MMD, and until now nothing in the port had
looked inside either. Both are hard to read: MMD is translated Fortran with `goto` labels and
1-based arithmetic, AMD is 1800 lines in one function over eleven parallel `int` arrays. This
experiment rebuilds the family one mechanism at a time, so that what each code *does* can be
separated from how it is *written*. Reference / teaching only, **not** part of the main Oblio
build.

Each layer exists twice, as `<name>.py` and `<name>.cpp`, printing the same trace. The Python is
the one to read and the C++ is the one to trust; `make test` checks they have not drifted.

## The layers

Each adds exactly one mechanism to the one before. The right column cites
`archive/sparse_factorization.md`, where the prose lives.

| layer | adds | sections |
|---|---|---|
| `md1` | naive minimum degree, materializing the fill | 5.1, 5.2 |
| `md3` | the quotient graph: cliques instead of fill | 5.3, 5.4 |
| `md4` | supervariables and mass elimination | 5.5, 5.6 |
| `md5` | maintained degrees, refreshed only where they changed | 5.7, 5.8 |
| `md6` | degree buckets, so the minimum is walked to, not scanned | 5.9, 5.10 |
| `mmd` | multiple elimination: a batch of pivots per refresh | 5.11, 5.12 |
| `amd` | approximate degree: a bound instead of a set union | 5.13, 5.14 |

`md2` is an earlier naive variant, kept for reference and not built.

## What the layers show

**`md1` through `md6` return the same ordering.** Every one of them. That is the point of the
first five sections: the heuristic was fixed in `md1` and everything after is implementation, so
the layers can be verified by demanding an identical permutation.

One refinement, because "same ordering" is not quite true across the whole run. Three things can
happen when a layer is added, and all three occur here:

```
                              order        fill        what the change is
md1 -> md3                    same         same        a change of representation
md3 -> md4  (mass elim.)      DIFFERENT    same        a provably free reordering
md4 -> md5 -> md6             same         same        a change of implementation
md6 -> mmd  (multiple elim.)  different     DIFFERENT  a wager
```

So `mmd` is not the first layer to change the permutation. Mass elimination already does, on nine
of twelve test graphs, with identical `nnz(L)` on all twelve. What `mmd` is first to change is the
**fill**.

**`mmd` and `amd` give up different things.** `mmd`'s pivots are always true minimum-degree
vertices; only the tie among equals is broken differently, because a batch evicts what it has
touched. `amd`'s pivots may simply not be minimal, because an overcounted bound can hide the true
minimum. MMD perturbs; AMD can be wrong. Both cost well under a percent of fill, in either
direction, and both are noticeably faster.

## The test graphs

`graph1` is a 4-cycle, `graph2` has six vertices, `graph3` has twelve and is the first whose
ordering is not the identity. `graph4` has eight vertices and fourteen edges and exists for one
reason: **it is the smallest graph we found on which AMD's bound is ever loose.** The bound
overcounts only when a vertex belongs to two elements that overlap outside the new one, which
needs enough eliminations to have made several elements and enough fill for them to intersect.
Checked exhaustively, no connected graph on five or six vertices is ever loose anywhere in its
run, and none in thirty thousand samples on seven. Without `graph4` the `amd` trace would display
the whole algorithm and never once show it approximating.

## Build

```
make        build every C++ prototype
make test   build them, run each, and check the C++ agrees with its Python twin
make clean
```

`make test` compares traces after stripping spaces, brackets, braces, quotes and commas, because
Python renders a dict as `{0: [1,3]}` and the C++ renders `0: {1,3}`. Every number, name and label
still has to match. The alternative was editing fourteen files to unify the spellings, which
seemed the worse trade.

## What is not implemented

The `mmd` and `amd` prototypes are deliberately subsets of the vendored routines. Each file header
carries its own list, and sections 5.11 and 5.13 carry the same lists in prose. In brief:

- **`mmd`** lacks the prepass that numbers degree-0 and degree-1 vertices before the main loop, and
  `mmdupd`'s `q2h` merging of vertices indistinguishable *from each other* rather than from the
  pivot. It also files degrees at a different offset and never uses bucket 0.
- **`amd`** performs its degree update in one pass where `Amd.cpp` uses two, which changes the
  ordering on four of eleven test graphs. It also lacks dense-row handling and the postorder, both
  of which change the output, along with `amd_aat`, `amd_valid`, the `Control`/`Info` interface and
  the workspace compression, none of which do.

## Two bugs this found, both ours

Worth recording, because both were invisible to the checks in place at the time.

**`amd` did not shrink the new element on mass elimination.** When a vertex is mass-eliminated it
joins the pivot's supervariable, so it stops being outside the new element and must stop
contributing to `|L|`; `Amd.cpp` does this at `degme -= nvi`. We computed `|L|` once before the
loop. The effect is nearly invisible: identical results on all four test graphs and on every grid,
surfacing only on a five-vertex bowtie where a bound came out one too large.

**`mmd.cpp` printed display lines its Python twin did not**, left over from the `md6` file it was
derived from. This survived because the verification at the time used a `grep` filter narrow
enough to skip exactly those lines, which is not a test. `make test` exists because of this one:
it compares whole outputs, and it found the drift immediately.

## Related

- `archive/sparse_factorization.md` section 5, the prose, pseudocode and worked examples.
- `src/Mmd.cpp` and `src/Amd.cpp`, the vendored routines these are read against.
- `src/OrderEngine.cpp`, the glue that calls them.
