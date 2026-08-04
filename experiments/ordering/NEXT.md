# Next: why amd2 and Amd2 disagree on a 10x10 grid

A handoff note, written 2026-08-03 at the end of a long session. It says where the ordering
experiment stands, what the next question is, and how far it was chased, so that the next session
can start on the question rather than on rediscovering the reasoning. Delete it once the answer
lands and is recorded in this folder's README.

The same material is item 5 of the ordering questions in `docs/TODO.md`, which is the durable
record; this file is the kickoff and is meant to be thrown away.

## Where things stand

Thirteen layers, thirteen twin agreements, four production agreements. `make test`:

```
md1 md2 md3 md4 md5            the line, one idea per rung
mmd1 mmd2                      the MMD fork, base and extras
amd1 amd2 amd3                 the AMD fork, base, extras, and the vendored remainder
mdm2 mda2 mdam2                the md2 square, three of its four corners
```

**Two things landed this session.** The md2 square, four files isolating exact-against-bounded from
recomputed-against-maintained, described in the README section "Zooming in on md2". And the amd2
split: the old amd2 became amd3, and a new amd2 is amd1 plus the two extras and nothing else, which
made it comparable to production `Amd2` by PERMUTATION rather than by fill. `PORTED_FILL` in the
Makefile is now empty and every ported layer is on the strong oracle.

That split is what exposed the question below. It was invisible while amd2 was checked by fill.

## The question

**amd2 and production's Amd2 produce different permutations on grids of 10x10 and larger.**

```
seven example graphs      permutations agree        make test passes
grid 8x8, grid 9x9        permutations agree
grid 10x10 and larger     permutations DIFFER
amd1 against Amd1         agrees on every grid tried
```

It predates the split. The old amd2 was checked by fill on the same seven graphs, which is weaker in
two ways at once, so this has probably been there since amd2 was written.

## How far it was chased

**It is aggressive absorption, not the hash.** On a 10x10 grid the first difference is at step 4.
Both pick pivot 10, both report supervariable size 1, neither has made a hash merge yet, and the
prototype absorbs one clique where production absorbs none.

**Production's `touchedCliques` is empty at that step and the prototype's is not.** So the two
disagree about which cliques the step touches at all, before any absorption test is applied.

**Which means the graph state has already diverged before any merge of either kind.** The pivots
agree, no hash merge has fired, and every step up to there reports size 1, so no mass elimination
either. Nothing that changes the graph has run except the eliminator itself.

**So start at the eliminator or at absorption's purge**, not at either extra mechanism. The
candidates worth checking first, in order: what `qg.absorb` does to `I[u]` against what the
prototype's list comprehension does; whether the incidence lists differ in ORDER, which the
tie-break makes visible; and whether `eliminate`'s pruning differs in some case the seven small
graphs cannot reach.

**One hypothesis tested and rejected, recorded so it is not repeated.** `cliqueDegree` in
`src/Amd2.cpp` is written once when a clique is formed, and the comment above it lists the
invariants that keep it exact WITHOUT mentioning hash merges. A hash merge does break it: the test
requires `I[u] == I[v]`, so any other clique in that set holds both `u` and `v`, and merging `v`
away leaves its weight counted twice. Maintaining `cliqueDegree` at the merge site changes nothing
on this reproduction, because the divergence appears before any merge. **The staleness is still real
and still unmentioned by that comment**, so it is worth fixing on its own account even though it is
not this bug.

## Reproducing it

Both sides need a per-step line: pivot index, pivot, supervariable size, absorbed count, hash count.
The first line that differs is the one to open up.

**The prototype** filters its own output through a `CounterSink` that keeps only lines matching a
fixed list of prefixes, so a new line has to be added to that list in `main` as well as printed:

```
./amd2_cpp grid 10
```

**Production** needs a `printf` at the end of the elimination loop in a scratch copy of
`src/Amd2.cpp`, driven by a small main that builds the same grid and calls `orderAmd2`. The grid the
prototypes use is `gridGraph(side)` in any of the `.cpp` twins; five-point, row-major, which is what
the scratch main has to match exactly or the comparison is meaningless.

## One bug found and fixed on the way

The new amd2 inherited amd1's `external_degree = len(C[pivot])`, an unweighted count, where hash
detection folds a vertex into a LIVE one so a member of the new clique can stand for several
original vertices. amd3 had the weighted version with a comment saying exactly why. Fixed in both
twins; the prototype's `nnz(L)` now matches production at 8x8 where it did not before. Larger grids
still differ, but that is the permutation divergence above rather than the accounting.

## And the harness is the other half of the finding

**Seven graphs of at most twelve vertices cannot validate a mechanism that only fires on larger
structures.** `make test` comparing prototype against production on a grid or two would have caught
this when amd2 was written, and would cost almost nothing: the prototypes already have a grid mode,
and production needs only the small driver described above. That is the same test-set point at the
top of `docs/TODO.md`, arriving from a new direction.

## After that

- **The fourth corner of the md2 square is written and the ladder above it is not.** mda4 and mda5
  would show the bounded branch under maintained degrees and under buckets, and the README's rule
  says both would use the bound against `C[p]` rather than mda2's pivot-free form.
- **The mark-and-tag counters are `int32` with no guard**, which item 4 of the ordering questions in
  `docs/TODO.md` covers in full. Not urgent, and the fix is decided: keep `int32` and port the
  sweep, which can be an unconditional reset here because our mark array carries nothing but tags.
- **What the hash actually buys is measured and small.** On grids it fires hundreds of times and
  saves one pivot or none, the saving not growing with n. Most of what it finds, mass elimination
  would have found a step or two later. It is not zero: on graph3 it finds a pair mass elimination
  cannot see at all. Recorded in the README's supervariable section, and the caveat that grids are
  the worst case for judging it still stands.
