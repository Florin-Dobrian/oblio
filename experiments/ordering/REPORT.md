# Where Oblio's own orderings stand against the vendored ones

Written 2026-08-03, at the point where MMD1, MMD2, AMD1 and AMD2 first became stable together:
each agrees with its prototype on the seven examples and on grids, the prototypes agree with their
Python twins, and the two defects in the amd2 prototype found that day are fixed. This is a
snapshot for picking the work up later, not a plan. It records what we measured, what we conclude,
and what we do not yet understand.

The question behind it is narrow. Our MMD and AMD were written by reading genmmd and AMD and
porting their ideas, so those two routines are the oracle. Not for the permutation, which is a
heuristic and may legitimately differ, but for the two things a heuristic is judged on: how little
fill it leaves and how fast it runs. Both vendored routines have been in service for decades. Where
we come out worse, the working assumption is not that we made a different trade. It is that we got
something wrong.

## How to read the numbers, before reading any

Two different confidence levels are mixed below and they should not be read the same way.

- **Fill is exact.** Every nnz(L) here comes from Oblio's own symbolic factorization of the emitted
  permutation, including for the vendored routines, so no counter is compared against a different
  counter. It is deterministic and reproduces to the digit.
- **Time is indicative only.** These runs are on the Linux sandbox, which is not a measurement
  platform. Repeating one case three times gave MMD between 12.91 and 16.63 ms and AMD2 between
  23.76 and 28.97 ms, a spread of 20 to 25 percent. Only the large ratios below survive that, and
  every one of them needs re-measuring on alpamayo before it is quoted as a number.

A third caveat has nothing to do with noise. **Every timing and fill figure in
`benchmarks/ordering/README.md` is from a 2D grid**, and the single most important thing in this
report is that 2D and 3D disagree about how good our code is. What follows measures both.

## The measurement

Six orderings, two problem families, three sizes each. Five-point 2D grids at 64, 100 and 140 a
side, and seven-point 3D grids at 20, 26 and 32 a side, the largest of which is n = 32768 with far
more fill and much larger cliques. MMD and AMD are the vendored routines; MMD1, MMD2, AMD1 and AMD2
are ours. Time is ordering only, best of five in 2D and best of three in 3D.

Each of ours against the vendored routine it was ported from, as a time multiple and a fill
difference:

```
2D              MMD1             MMD2             AMD1             AMD2
  64      1.59x +13.4%     1.01x +12.0%     1.41x  +1.1%     1.71x  +2.4%
 100      2.03x +19.8%     1.19x +16.2%     1.19x  -2.2%     1.75x  +3.0%
 140      2.30x +19.4%     1.44x +21.6%     1.83x  -4.1%     2.45x  +2.6%

3D              MMD1             MMD2             AMD1             AMD2
  20      3.79x  +2.4%     1.04x  +2.8%     1.26x  +4.2%     2.52x  +5.3%
  26      5.13x  +7.6%     1.06x  +5.7%     1.58x +13.7%     4.41x  +9.1%
  32      5.69x  +0.1%     0.86x +12.4%     1.45x +11.5%     3.33x  +6.3%
```

The raw figures behind it, time in ms and nnz(L):

```
2D 140x140, n = 19600          3D 32x32x32, n = 32768
  MMD    3.02    412921          MMD   35.25   7898321
  MMD1   6.94    492921          MMD1 200.68   7903225
  MMD2   4.35    501951          MMD2  30.18   8876907
  AMD    3.11    474995          AMD   16.39   7746501
  AMD1   5.69    455472          AMD1  23.81   8639479
  AMD2   7.61    487111          AMD2  54.65   8232585
```

## What is solid, and it is most of it

The foundations hold, and the evidence is that AMD1 is a competitive ordering. In 2D it **beats**
the vendored AMD on fill, by 2.2 percent at 100 a side and 4.1 percent at 140, while being within
a small time multiple. AMD1 is the approximate degree bound and nothing else, so everything
underneath it is carrying its weight:

- **The quotient graph representation.** Adjacency as a flat array, cliques in a bump arena,
  incidence as a vector per vertex. This is a different encoding from the vendored workspace pool
  with its `Iw`, `Pe`, `pfree` and compaction count, and the deliberate decision not to port that
  pool is not costing us. A structure that produced less fill than AMD in 2D is not a structure
  with a representational defect.
- **Element absorption and pruning**, and the conservation lemma `|A[u]| + |I[u]| <= deg(u)` that
  governs them.
- **Mass elimination**, whose partition the rest of the ordering depends on.
- **Maintained degrees in buckets**, rather than recomputed or scanned.
- **Both degree kinds.** The exact degree in the MMD branch and the approximate bound in the AMD
  branch, with the bound's three-way minimum. The bound never falls below the true degree, checked
  now at every step in the prototypes and verified directly in production by comparing against
  `QuotientGraph::reachableWeight`, zero violations on grids of 10, 20, 32 and 50.

So the report is not that the design is wrong. It is that a small number of decisions are wrong,
and they are identifiable.

## Finding 1: 2D flattered us, and all our published numbers are 2D

AMD1 against vendored AMD on fill:

```
2D  64     +1.1%          3D  20      +4.2%
2D 100     -2.2%          3D  26     +13.7%
2D 140     -4.1%          3D  32     +11.5%
```

The same layer that beats the oracle by 4 percent in 2D loses to it by 12 to 14 percent in 3D, and
the 2D advantage improves with size while the 3D deficit appears and then holds. Two conclusions
follow, and the second is the uncomfortable one.

The first is that our test set has been showing us the flattering half of the picture. Everything
in `benchmarks/ordering` is a 2D grid.

The second is that **this gap is in AMD1, which is the bound and nothing else.** It is not in the
late-stage extras, because AMD1 has none. Whatever produces it is in the base: the bound itself,
the tie-break, the filing convention, the degree floor, or the order in which mass elimination and
the refresh interleave. That was not where we expected to be looking.

## Finding 2: our two branches fail in mirror-image places

**In AMD the extras are the problem and the base is fine.** AMD1 is 1.19x to 1.83x of the vendored
routine in 2D and beats it on fill; adding the extras to get AMD2 makes it slower everywhere and
worse on fill in 2D.

**In MMD the base is the problem and the extras are the fix.** MMD1 is 3.79x, 5.13x and 5.69x of
the vendored MMD in 3D, and the gap grows with size. MMD2 adds the six passes genmmd has that MMD1
lacks, and the result is 1.04x, 1.06x and 0.86x. At the largest 3D case MMD2 is **faster than the
vendored routine**. So genmmd's extras are not extras at all: they are load-bearing, exactly as the
experiment README says, and the cost of not having them is a factor of five.

The two branches therefore need opposite work, which is worth knowing before either is opened.

## Finding 3: AMD2's extras are a net loss, and the hash is almost all of it

AMD2 adds aggressive absorption and hash supervariable detection to AMD1. Gating each independently
in the prototype, on 2D grids:

```
grid 64x64          pivots   nnz(L)   degree comps   clique reads   hash merges
  both               3084    68822          22640          29212          1012
  absorption only    3084    67950          25684          32185             0
  neither            3084    67950          25684          32256             0
```

- **Aggressive absorption is worth nothing measurable here.** It fires once per run on a grid. It
  does move the permutation, but the fill is identical to the digit on every 2D and 3D grid tried,
  and it saves 0.4 percent of clique reads.
- **The hash fires constantly and saves no pivots at all.** 1012 merges at 64 a side and the run
  still takes exactly 3084 pivots, and the same holds at 16, 32 and 45 a side. That is the direct
  evidence for the suspicion that mass elimination reaches the same partition a step or two later.
  It does buy 12 percent fewer degree refreshes and 10 percent fewer clique reads.
- **And it costs far more than that.** Gating the hash out of production:

```
                     AMD1     AMD2   AMD2 no hash    the hash's share of AMD2's penalty
2D 140x140           4.51     6.38       5.04              1.34 of 1.87 ms      72%
3D 26x26x26         13.82    34.15      15.44             18.71 of 20.33 ms     92%
```

## Finding 4: but the mechanism is right, so this is ours to fix

Gating the vendored routine's own supervariable detection off, in a scratch copy of `Amd.cpp`:

```
                    with hash        without hash        removing it costs
2D  64x64             1.02 ms          1.05 ms        +2.9% time   +5.0% fill
2D 100x100            2.29 ms          2.30 ms        +0.5%        +1.6%
2D 140x140            4.37 ms          4.39 ms        +0.5%        +0.4%
3D 12x12x12           1.08 ms          1.20 ms       +10.8%        +2.6%
3D 20x20x20           5.11 ms          5.57 ms        +8.9%       +11.0%
3D 26x26x26          11.20 ms         12.67 ms       +13.2%        +8.7%
```

Taking the hash out of the vendored routine makes it slower **and** fills more, in every case, and
most of all in 3D, which is precisely where ours is worst. So hash supervariable detection is a
good mechanism that pays for itself twice over, and the conclusion is not that it is not worth
having. It is that our implementation of it is not yet worth having.

## The one gap we can explain

For the time half of AMD2's hash, the cause is visible in one line of each implementation.

`AMD_2` accumulates its key **inside scan 2**, as `hval += e` and `hval += j`, on lists it is
already walking to compute the degree and compact the pool. The key costs a few additions on data
already in registers.

Ours builds the key in a **separate traversal** of `A[u]` and `I[u]` for every member of `C[p]`,
on top of the bound loop that has just walked both of those same lists for `explicitPart` and the
`outside` sum. We pay a full extra pass over the structure per step for a quantity the vendored
routine gets for free. `C[p]` grows with fill, which is why the penalty is 72 percent of AMD2's
overhead in 2D and 92 percent in 3D.

The indicated repair is to compute the key in the bound loop and store it, leaving the reverse pass
to read a stored key so the tie-break is preserved. Two cautions before that is budgeted. It needs
an array of size n, which is exactly the footprint trade that made AMD1B slower at large n after
being faster at small. And it would not touch the fill half below.

## The gaps we cannot explain yet

These are the open questions, and they are the reason this is a report and not a fix.

1. **Why AMD1 fills 12 to 14 percent more than vendored AMD in 3D while beating it in 2D.** The
   layer is the bound and nothing else, so the cause is in the base rather than the extras. This is
   the largest single unexplained result here.
2. **Why our hash costs fill in 2D when the vendored one saves it.** Removing the vendored hash
   costs 0.4 to 5.0 percent more fill in 2D; adding ours costs 2.4 to 3.0 percent. Same mechanism,
   opposite sign. It may be the interaction with our different tie-breaks. We do not know, and
   nothing in the time analysis accounts for it.
3. **Why MMD2 fills 12 to 22 percent more than vendored MMD.** This is the largest quality gap
   anywhere in the table, it is worst in 2D at 140 a side, and it coexists with MMD2 matching or
   beating the vendored routine on time. A layer that is fast and fills a fifth more than its
   oracle has traded something away, and we did not choose that trade knowingly.

## What to do next, in the order the evidence supports

1. **Re-measure on alpamayo.** Nothing above is a measurement of record. The fill columns will
   reproduce exactly; the time columns will not.
2. **Add 3D grids to `benchmarks/ordering`.** Every conclusion the benchmark currently supports is
   drawn from one problem family, and the other family disagrees with it. This is cheap and it
   changes what everything else is measured against.
3. **Open question 3, the MMD2 fill gap.** Largest quality gap, and MMD2's time is already
   competitive, so there is room to spend.
4. **Open question 1, the AMD1 3D fill gap.** Hardest, because it is in the base, and most valuable
   for the same reason: everything above AMD1 inherits it.
5. **Then the hash key fusion.** It is the best understood of the four and the only one with a
   known repair, which is exactly why it should not go first. It is a constant factor on a layer
   whose fill behavior we do not yet understand.

## What was not measured, and should not be assumed

- **Anything but structured grids.** Two families now instead of one, which is progress, but a
  matrix with genuinely irregular structure remains untested. This is the same test-set item that
  has been open in `docs/TODO.md` throughout.
- **Whether ours and the vendored routines find the same supervariable population.** The experiment
  README records both firing 2488 merges across the small test set. That was not re-verified here,
  and the opposite-signed fill results argue for re-checking it rather than trusting it.
- **The effect of any of this on factorization time.** A previous session established that fill
  differences among orderings that are all roughly good do not propagate into factorization time:
  19 percent of fill spread produced 5 percent of factor spread with no visible correlation. The
  gaps in this report are of that magnitude, so their end-to-end value is unknown and may be small.
  What makes them worth chasing is not the milliseconds. It is that each one is a signal that we
  read something in the vendored code and implemented it differently than we thought we did.
