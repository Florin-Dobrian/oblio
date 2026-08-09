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

**REVERSED ON FILL, 2026-08-08, and it was a defect rather than a mechanism.** The alignment method said to
treat this finding as unresolved until the permutations matched, and they now do. `Amd2` and
`Amd2B` filed a supervariable one bucket too high per vertex a hash merge absorbed, the bound
having subtracted the vertex's own weight before the merge that grows it. Corrected, AMD2's fill
is 11900 at 32 a side, 199386 at 100 and 444191 at 140, so it beats AMD1 at every size by 1 to 3
percent where this section records it as worse. The extras were not costing fill; the filing was,
and the hash was being charged for it, which is exactly the attribution alignment was supposed to
buy.

**The time half stands.** Nothing here makes the hash cheaper, and its 72 to 92 percent share of
AMD2's overhead is unaffected. The gating measurement in "The one gap we can explain" is still
the right one.

The numbers below are kept as the record of the run that produced them.


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

## Where the ordering time gap comes from, measured 2026-08-07

Once MMD3 matched genmmd's permutation exactly, fill stopped being a variable and the time gap
became a clean question with a single answer available. This section records that answer, because
the largest part of it is a DELIBERATE DESIGN CHOICE rather than a defect, and the point of
writing it down is to be able to say where the difference comes from without re-deriving it.

**The counters say there is no algorithmic difference left.** Instrumenting genmmd with the same
counters the prototype prints, on square grids from 32 to 140 a side:

```
                 mmd3      genmmd
degree updates   equal      equal      1575, 5515, 12644, 23861
pair merges      equal      equal        236,  988,  2446,  4826
outmatched       equal      equal        595, 2088,  4635,  8469
ncsub            equal      equal
bucket probes      289        193      the only counter that differs, and it is O(n) either way
```

Same permutation, same fill, same work. So the remaining time is implementation and nothing else.

### The width of our arrays is the largest single contributor

`QuotientGraph` stores positions and measures as `std::size_t`, eight bytes; genmmd stores
everything in `int`, four. **Narrowing one of ours to test that is what we tried first and it
proves nothing**: `mWeight` alone moved the total by nothing at 140x140 and about four percent at
400x400, unreadable against noise, because it is one of five eight-byte arrays read in the same
loop.

So the experiment was inverted: widen the ORACLE. A scratch copy of genmmd with every working
array retyped from `int` to `std::int64_t`, verified to return the identical permutation, doing
byte-for-byte the same work and differing only in width. (`int64_t` and not `size_t`, because
genmmd stores negatives as sentinels, `invp[mn] = -num` and `bwd[nd] = -maxint`.)

```
                          140x140     400x400
genmmd, int                 1.000x      1.000x
genmmd, int64_t             1.170x      1.255x
our MMD3                    1.373x      1.619x
```

**Width alone costs 17 to 26 percent.** Confirmed from our own side afterwards, narrowing the four
measures that are bounded by `n` and leaving the two positions alone:

```
                                    140x140     400x400
MMD3, size_t                          1.373x      1.619x
MMD3, four measures narrowed          1.159x      1.464x
```

Sixteen percent at 140x140 and ten at 400x400, from four member declarations. The two experiments
agree: genmmd widened is 1.17x of itself, ours narrowed is 1.16x of genmmd. That agreement is the
check that the mechanism was identified rather than a number being fitted.

**It is not being taken, and that is the design trade.** `docs/CODING_RULES.md` defines a position
as an offset into a vector which "measures, so it is never negative, never `NIL`, and free to
exceed 2^31", and a measure is the same kind of quantity. Nothing in `QuotientGraph` violates that
rule; the finding CONFLICTS with it. Keeping `std::size_t` for sizes and positions is a decision
about what the code MEANS, not an oversight, and the price of it is now known rather than
suspected. The subsection after next follows the question to where it actually leads, which is not
the ordering at all but whether `A` and `L` can share one encoding.

**And the one operational conclusion: this cannot be tested piecemeal.** One array at a time sits
below the noise floor and will read as "no significant difference" every time it is tried, which
has happened more than once. Only the whole set moves the number.

**The six arrays are not one case but three**, which matters because only one of them is hard:

```
bounded by n         mAdjacencySize  mIncidenceSize  mCliqueSize  mWeight   int32 always safe
bounded by nnz(A)    mSourcePtr                                        int32 if nnz(A) is capped
bounded by nnz(L)    mCliquePtr                                             needs the width today
```

The first four can never exceed 2^31, `n` being an `int32` by construction, so narrowing them
costs only the rule's wording. `mSourcePtr` indexes an arena laid out once at construction and
never grown, so it follows `nnz(A)`. Only `mCliquePtr` follows the factor, and the next subsection
is about why that is a choice rather than a necessity.

### Why AMD does not have this problem, and what it costs to not have it

The width question is really a storage question, and following it one step further explains why
`mCliquePtr` is the one array that genuinely needs eight bytes.

**AMD's arena cannot overflow, and that is a theorem rather than a margin.** Its workspace is
allocated once at `slen = nzaat + nzaat/5 + 7n`, roughly `1.2 * nnz(A + A')`, and never grows.
When `pfree` reaches `iwlen` it COMPACTS, reclaiming the space of absorbed elements and counting
the compactions in `ncmpa`. That is sound because of the conservation bound in section 5.3 of
`archive/sparse_factorization.md`: `|A_i| + |C_i| <= deg(i)`, since every elimination that reaches
`i` replaces at least one source with the new clique and nothing manufactures a source, and every
later mechanism either destroys sources or empties a list. So the live quotient graph never
exceeds what the input occupied. The 20 percent elbow room buys fewer compactions, not
correctness. **AMD refuses only at its entry guard, on `nnz(A)`, never at runtime.**

**Ours grows to `nnz(L)` because we never reclaim, not because the cliques are that large.** The
same theorem applies to us: the LIVE cliques fit in `nnz(A)`. `mCliqueArena` appends every clique
ever formed and lets dead ones lie, so the difference between `nnz(A)` and `nnz(L)` is entirely
garbage. At 400x400 that is 4.86 million entries against roughly 800 thousand, six times the
memory, all of it dead.

**Which means the width and the reclaim are one decision, not two.** `mCliquePtr` needs the width
only because the arena tracks the factor. With compaction it would be bounded by `nnz(A)` like
`mSourcePtr`, and a cap on the matrix would cover it. **You cannot take the narrowing without
taking the reclaim**: the compaction is not a memory optimization there, it is what makes the
narrower type sound.

The experiment README lists AMD's `iwlen`/`pfree`/`ncmpa` machinery under "what is deliberately
excluded", as a consequence of packing state into reusable arrays rather than a feature of the
ordering. That is true of the mechanism and wrong about the consequence: it is exactly what bounds
AMD's storage by the matrix instead of the factor.

**But AMD's position is not available to us, and it is worth being exact about why.** AMD never
factors. It sees only `A`, so `int32` throughout costs it nothing: one guard at the entry on a
quantity the caller handed it, and the conservation bound guarantees that nothing downstream
exceeds it. It never has to represent an OUTPUT. So AMD is not evidence about what a solver should
do; it is evidence about what a routine that stops before the factor exists can do.

**`nnz(A)` is an input and `nnz(L)` is an output, and that is the whole difficulty.** `A` can be
capped once, before anything starts, and the caller learns immediately. `L` cannot. With a static
factorization the symbolic phase at least predicts it, so one test before allocation would serve.
With PIVOTING it is not predictable at all: a delayed column grows its parent's front, so the
bound would have to be re-checked every time a column is delayed. That check would live in the
delay path, which is the most delicate code in the numeric factorization, and it would fire after
real work had been done on an input that was accepted. Aborting with a partial factor or unwinding
are both poor answers to a question that should have been settled at the door.

**So there are two consistent positions, not three, and the trade is consistency against speed.**

```
size_t positions everywhere    A and L alike, no cap on either, no exception anywhere.
                               Costs the width in every hot loop: 17 to 26 percent here.
int32 positions everywhere     Requires capping nnz(L) as well as nnz(A), which means a
                               running check through the pivoting path, and a failure mode
                               that arrives mid-factorization.
```

Encoding `L`'s positions differently from `A`'s is the third thing one could do, and it is not a
boundary between phases but a SINGULARITY: `A` and `L` are the same kind of object, a sparse
matrix, so a different encoding for one of them means every piece of code touching both has to
know which it is holding. The libraries that ship do not do this. SuiteSparse carries ONE source
with a width parameter, the `using Int = int32_t` alias visible in `private/Amd.cpp`, and the
caller selects 32 or 64 at build time; MUMPS offers the same choice. `A` and `L` share a type in
either configuration, so neither has the exception. Note that this does not answer the question,
it relocates it: a 32-bit build still has to cap `nnz(L)` or overflow, and the decision moves to
whoever configures the build.

**Oblio holds the first, deliberately, as of 2026-08-07.** Consistency over the measured 17 to 26
percent. The full reasoning is a statement about Oblio's integer model rather than about the
ordering, so it lives in `docs/DESIGN_DECISIONS.md`, recorded 2026-08-08: the split is by
DIMENSION, one dimensional quantities bounded below 2^31 and two dimensional ones held in 64
bits, which is why the band between a 26 GB and a 52 GB factor that troubles a 32-bit build
never arises here. What belongs in this report is the price, and the price is the measurement
above.

**One narrowing remains available without touching that decision**, and it is where most of the
measured gain was: the four MEASURES, `mAdjacencySize`, `mIncidenceSize`, `mCliqueSize` and
`mWeight`, are bounded by `n` rather than by any `nnz`, so they are not positions at all and could
narrow while every position stays `std::size_t`. Measured on its own that configuration ran 1.159x
of genmmd at 140x140 against 1.373x, so it is most of the width gain for none of the consistency
cost. It is not being taken now, but it is the cheap half if the question returns.

**And it points at a missing category in the type rules, 2026-08-07.** `docs/CODING_RULES.md` has
two kinds of integer, an index that names an entity and may be `NIL`, and a position that offsets
into a vector and is "free to exceed 2^31". Everything above falls between them. The distinction
that actually holds is DIMENSIONAL:

```
kind       what it is                          bounded by      type
index      one dimension, may be NIL           n <= 2^31       std::int32_t
count      one dimension, never NIL            n <= 2^32       the missing one
position   into an n x n object                may exceed      std::size_t
```

**A count is bounded by a SIDE, a position by an AREA.** That is a fact about the object rather
than a convention, so unlike "free to exceed 2^31" it needs no defending. It also makes the two
boundaries land in different places, which is right: an index gives up a bit to `NIL` and so caps
at 2^31, while a count needs none and reaches 2^32, comfortably above any representable `n`.

The taxonomy and the measurement agree about where the line falls, which is the reason to trust
it. The four arrays that carried the 17 to 26 percent are all COUNTS by this reading;
`mSourcePtr` and `mCliquePtr` are positions and stay wide. Front size and update size in
`SymFactor` are counts too, and are `std::size_t` today for the same reason.

**There is no 32-bit `std::size_t`**, by design: it is defined to be able to size any object. So a
count would be `std::uint32_t`, or an alias over it, which has the additional merit of staying a
distinct type from `std::int32_t` so counts and indices cannot silently interchange.

**Both branches measured together, alpamayo, 2026-08-08.** `make width` in this folder builds a
second copy of each vendored routine with its integer type widened to `int64_t` and times the two
against each other in one run, with AMD's fill and MMD's permutation asserted equal on every row so
the work is known to be identical.

**Read this as a RANGE and not as digits.** Two runs of the identical binary disagreed by a factor
of two on MMD at 140 a side, 14.5 percent against 6.6, and both reported the WIDER build as faster
at the smallest sizes. The harness now picks its repeat count from a timed probe so every row runs
for about the same wall time, which was the cause; before that, the small rows were four
milliseconds of measurement and the large ones seven orderings. What survives across runs:

```
                     AMD          MMD
n up to ~10000    within noise, both signs seen
n = 19600         10 to 19%    7 to 15%
n = 40000+        14 to 33%   13 to 37%
```

**Width is a shared effect and it explains none of the branch gap.** The two branches pay within
noise of each other at every size, and MMD3 runs at 1.2x its vendored routine where AMD3 runs at
2.55x. Whatever separates them, it is not the width of the arrays.

**Three claims of mine died here and all were one error.** I measured this on the Linux sandbox
first and read AMD's penalty as larger and steeper than MMD's, projecting that a count sweep would
take AMD3 from 2.55x to about 2.1x. On alpamayo the penalty at 140 is around 10 to 19 percent, not
the 19 to 31 the sandbox showed, so the sweep is worth a few percent there and not half the gap.
The "AMD's is worse than MMD's" comparison was made across two harnesses run at different times,
which measures the harnesses as well; running them together reverses it. And I then quoted single
digits off one run of a harness too coarse to support them.

**Comparing a fresh number against a recorded one is not a comparison**, and a number is not a
measurement until the harness has been shown to reproduce it. Both are worth more than the result.

The count sweep remains worth doing, for the reason it was always worth doing: a count is not a
position and the code should say so. It is not a performance item at the sizes that matter here.

**Two harness notes, because each cost a misreading.** The first routine timed at the smallest size
carried the core's clock ramp and reported 0.19 ms against the second routine's 0.07 on the same
work; taking the minimum does not remove that, since the ramp outlasts the loop at n = 1024, so one
untimed call per routine per row is discarded. And the repeat count is chosen from a timed probe
targeting 400 ms per row, which is 1596 repeats at n = 1024 and 20 at n = 160000, where a fixed
count gave 25 and 7.

**A second instance, 2026-08-08, and it is smaller and sharper than the width one.** The amd hash
pass walks a vertex's incidence list with an `std::size_t` counter and stamps cliques at
`c + static_cast<std::int32_t>(size)`. Neither quantity is two dimensional: the counter is bounded
by `deg(v)` and so by `n`, and `size` is the matrix order. Both are counts. With a count type the
line would read `mark[incidenceV[i] + size]` with no cast at all, and the counter would be the
width of what it holds.

The cast was named as `cliqueStamp`, one crossing per driver rather than four, which is worth doing
whatever type the counts eventually get. **The counters are simply wrong and are owed a sweep**:
`std::size_t` is for a position into an area, these are bounded by a side, and until a count type
exists the honest spelling is `std::int32_t` with the bound cast in the condition.

**One argument against that sweep was made here and withdrawn, and it is recorded because it was
nearly convincing.** `benchmarks/ordering/README.md` measures hoisting these loops' bounds at 2.3
percent worse for AMD2, since they exit on the first mismatch. That is a fact about HOISTING THE
LOAD and says nothing about the counter's type: an `std::int32_t` counter with the call left in the
condition performs the same load and the same comparison at a different width, and the cast folds
into the compare. Reading a measurement of one change as evidence about a different change is the
same error as the entry-4 sentence that ruled the amd branch out by an argument that reached past
where it was checked.

That is the useful shape of the finding: **the missing category shows up first as casts nobody can
place and loop counters nobody can type**, well before it shows up as time. The width measurement
was 17 to 26 percent and is the reason to care; this is the reason it will keep being noticed.

**One caution, and it is already live.** Unsigned counts wrap on subtraction, and the ordering
subtracts constantly. `Mmd3.cpp` has

```cpp
degrees[u] = std::max<std::size_t>(degree - qg.weight(u) + 1, 1);
```

which is safe by invariant, `dg0` including `u`'s own weight, but nothing guards it and the
`std::max` actively hides it: if the subtraction ever went negative it would wrap to an enormous
value and the clamp would select it, filing the vertex past the end of the bucket array. At 64
bits a wrapped value is 1.8e19 and would fault immediately; at 32 bits it is 4e9, which is far
more likely to look like a plausible degree. So narrowing counts is not only a width change. The
subtractions want a signed temporary at the point of use, or the invariant wants asserting.

**For reference, the sizes involved.** With `int32` row indices a matrix could in principle hold
`2^62` entries, `n` being up to `2^31` with every column dense. AMD refuses above
`2^31 / sizeof(Int)`, about 537 million, which is a 6.4 GB real matrix. `nnz = 2^31` is a 25.8 GB
real matrix or a 42.9 GB complex one; as a FACTOR, `2^31` entries is reached by a 2D grid at
`n = 5.6e7` or a 3D grid at 130 a side, both of which fit on a large workstation.

### What the rest of the gap is

At 400x400, `1.619 / 1.255 = 1.29x` remains after width. In descending order of what we know:

1. **Vector against flat arena.** Our clique members live in `mCliqueArena` behind `mCliquePtr`
   and `mCliqueSize`; genmmd's live in `adjncy` in a linked segment structure it walks in place.
   Not measured separately.
2. **Two arrays for one test.** Our membership test is `mMark[v] != mTag && (!live ||
   mEliminated[v] == 0)`, two loads from two arrays. genmmd's is `marker[nb] < tag`, ONE load,
   because its `marker` carries a permanent `maxint` sentinel for dead vertices, so a single
   ordered comparison decides membership and liveness together. AMD does the same fusion by the
   sign of `Nv[i]`. We are the only one of the three that separates them, and it is the price of
   our unconditional tag sweep, which is only possible because our mark array has no sentinel.
   Fusing would trade a load in the hot loop for a selective sweep and a derived ceiling, which is
   genmmd's exact bargain.
3. **Setup.** `QuotientGraph`'s constructor against `mmdint`, 194 ms against 57 in the
   Instruments trace of 2026-08-07, the worst ratio at 3.4x though only four percent of the
   run. This has NO independent fix: it initializes 61 bytes per vertex against genmmd's 28,
   and narrowing the four measures makes it 15 to 20 percent faster on its own, so it is the
   width item again rather than a separate one. Replacing its `reserve` plus `push_back` with
   a sized write and an index was tried and is SLOWER, because `resize` value-initializes the
   whole arena and then it is overwritten; the capacity check being removed is cheaper than the
   zero-fill being added.

### One thing that was fixed rather than documented

`QuotientGraph::orderAscending`, written on 2026-08-07 for MMD3's numbering, cost 244 ms of a
4.94 s profile where genmmd's `mmdint` and `mmdnum` together cost 116 ms, while doing strictly
less work. It allocated four arrays of size `n` and made six passes. Rewritten to one scratch array
and two passes, with the root's cursor and the member's root marker sharing that array by sign, it
is 60 ms, level with `mmdnum`. That was 14 percent of the whole gap over vendored MMD, from one
function, and it is the reason the numbers above are quoted against 1.373x rather than 1.43x.

## The gaps we cannot explain yet

These are the open questions, and they are the reason this is a report and not a fix. Two entries
that stood here on 2026-08-03 have since been answered and are recorded above rather than below:
the MMD fill gap, which was tie-break convention plus one defect and is now zero, and the largest
part of the MMD time gap, which is array width.

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
3. **Open question 3, the MMD2 fill gap: CLOSED 2026-08-07.** It was tie-break convention plus one
   defect, found by building MMD3 and aligning it to genmmd one divergence at a time. MMD3's fill
   is now genmmd's exactly at every size and MMD2's gap fell from about 20 percent to 7. See the
   experiment README's mmd3 section for the ledger and the wrong turns.
4. **Open question 2, why our hash costs fill where the vendored one saves it: CLOSED 2026-08-08.**
   It was not the hash. `Amd2` and `Amd2B` filed a supervariable one bucket too high per vertex a
   hash merge absorbed, the bound having subtracted the vertex's own weight before the merge that
   grows it. Corrected, AMD2 beats AMD1 at every grid size and the vendored AMD at the two larger
   ones. Found by building AMD3 and aligning it to `AMD_2` one divergence at a time; see the
   experiment README's amd3 section. Finding 3 reverses on fill and stands on time.
5. **Open question 1, the AMD1 3D fill gap.** Hardest, because it is in the base, and most valuable
   for the same reason: everything above AMD1 inherits it. Note that AMD2's fill is no longer
   evidence about it in either direction, the defect above having moved every AMD2 figure.
6. **Then the hash key fusion.** It is the best understood of the remaining items and the only one
   with a known repair. Its fill behavior is now understood, which was the reason given for not
   taking it first; what is left is the constant factor, and the time half of finding 3 is
   untouched by anything done since.
7. **And a new one, from the alignment.** AMD3 fills more than the corrected AMD2 on grids, 474995
   against 444191 at 140 a side, so our tie-break now beats AMD's by 6.5 percent where aligning
   MMD improved ours. Whether that survives outside 2D grids is item 2 above, and it is the second
   data point for the LIFO-against-FIFO question below.

## Open question: is LIFO actually better, or is genmmd merely good?

Raised by the mmd alignment, unanswered, and now with a second data point from amd. It was written
down in the handover note between the two sessions and moved here when that note was deleted.

**The mmd evidence.** The tie-break entries alone took mmd2's fill gap from about 20 percent to
about 10, and entry 5's defect alone took it from 20 to 7. So the tie-break is worth roughly 7 to 10
percent of fill on its own, which is a real effect and not noise. But converging on genmmd would
lower fill whether or not LIFO is inherently better, so the measurement does not separate the two.

**The amd evidence, and it points the OTHER WAY.** `amd3` is aligned, so its fill is the vendored
routine's: 474995 at 140 a side, where the corrected `Amd2` reaches 444191. On grids our tie-break
beats AMD's by 6.5 percent. Aligning MMD improved our fill; aligning AMD costs it.

**So the two branches disagree**, which is what makes the question worth answering rather than
assuming. Either LIFO is not inherently better and genmmd is simply a good ordering, or the two
vendored routines differ in ways the tie-break label does not capture.

**The experiment is cheap and is worth running once, from an ALIGNED baseline.** Take mmd3, flip ONE
walk back to FIFO, measure fill, restore it, repeat for each of the four; then the same for amd3's
six entries. From a matched baseline that isolates each walk's contribution and answers whether the
direction itself carries quality. Doing it before alignment would have been meaningless, which is
why it was not done.

**And it wants 3D grids, not 2D.** Both figures above are from square grids, which is one problem
family and the flattering one. The 3D result in this report is the reason to distrust a 2D-only
conclusion about fill.

That minimum degree is tie-break sensitive is well established, and is why the vendored routines are
used as oracles rather than reimplemented. Whether LIFO against FIFO specifically is characterized
in the literature has NOT been checked and should not be asserted either way.

## Leads and observations parked, 2026-08-06

None of these is measured and none may matter. They are here so that they are not lost, since
each came out of reading the vendored codes closely and each could plausibly bear on a gap above.

**1. We touch four arrays where AMD touches one, and it is a candidate for AMD1's time gap.**
Both our AMD layers already compute `outside[c]`, the count of a clique's members lying outside
the new clique, by subtraction rather than by walking: `outside[c] = cliqueDegree[c] - weightU`
on first sighting and `outside[c] -= weightU` after. That is AMD's scan 1, line for line. What
differs is where the state lives. Per touched clique we use

- `outside[c]`, the running count,
- `cliqueDegree[c]`, which AMD does not need because `Degree[e]` already holds it,
- `mark[c]`, to answer "first sighting this iteration",
- `touchedCliques`, built by push_back and walked again at the end of the iteration only to zero
  `outside`.

AMD's `W[e]` does all four jobs at once by storing the count OFFSET ABOVE the current stamp
`wflg`: the value above `wflg` is the count, seeding from `Degree[e]` needs no second array,
"below `wflg`" is the first-sighting test, and advancing `wflg` clears every clique at once with
no list and no second pass. So we pay roughly three extra cache lines per touched clique plus a
full second pass per iteration that AMD never makes.

Worth measuring, not worth believing yet. Reasoning about cache traffic has a poor record on this
project: cachegrind was wrong three times in one afternoon in three directions, and two other
performance leads raised the same day this was written measured as 291-to-1 irrelevant and as
right-mechanism-wrong-cause. What would settle it is counting touched cliques per iteration and
profiling the share taken by `touchedCliques` and the clearing pass.

It cannot explain any fill gap. Encoding does not change which pivot is chosen.

**2. Our own mark and tag have no overflow guard, and both vendored routines do.** `QuotientGraph`
holds `mMark` and `mTag`, and the drivers `Mmd2`, `Amd1`, `Amd1B`, `Amd2` and `Amd2B` hold their
own, all unguarded. Measured on 2D and 3D grids, the stamp advances about 15n over a run, so an
`int32` wrap needs n near 140 million and nothing is at risk at any size we run. The reason it is
still worth recording is that our `mMark` carries no permanent sentinel, so a wrap would give a
stale MATCH rather than a collision: silent, data-dependent, and benign-looking.

**The prototypes now have one, 2026-08-06**, at a `TAG_CEILING` of `2^30 - 1`, in all thirteen
layers, with a `tag sweeps` counter as the witness that it stays inert. Production still does not,
and the ceiling is a pragmatic placeholder rather than a derivation. The experiment README's "The
tag guard" section carries the rule for where a check may land, which is the part that transfers;
`docs/TODO.md` item 4 carries what is left.

**3. The two vendored routines set the ceiling 256 times apart, and only one of them derives it.**
AMD computes `wbig = Int_MAX_VAL - n`, with the header stating the rule outright: `wflg` may not
reach it, and the `- n` is headroom because the largest value ever stored is `wflg + |Le \ Lme|`
and that set size is bounded by n. So AMD's ceiling is exactly what the constraint demands, and it
is a variable so the same source is correct for `int32_t` and `int64_t`.

genmmd takes the ceiling as a caller argument, and OUR wrapper passes the historical `8388607`,
which is `2^23 - 1`, in `private/Mmd.cpp` line 202. It is the only occurrence in the repo. The
constraint there is the same in kind, since `mmdupd` computes `mt = tag + md0` before testing it,
but it would permit anything up to about `INT_MAX - n`. Why 1985 chose a value 256 times lower is
not recoverable from the code: the original Fortran comments did not survive the port. The
best-fitting guess is SPARSPAK's shared workspace with integer and real views of one array, where
a value passing through a `REAL` slot must be exactly representable in a 24-bit mantissa, which
puts the largest safe integer at `2^24`. That is a guess and should be marked as one. Netlib's
SPARSPAK source would settle it in one line.

Consequence worth knowing: with `8388607` and a stamp advancing at about 15n, the vendored MMD's
sweep becomes reachable around n of half a million. Measured, it fires zero times at every size we
benchmark, so it is not a hidden term in the MMD timing comparison. The wrapper is ours, not
vendored, so the value is our choice if that ever changes. It must stay below the maximum by at
least `n`; raising it to `INT_MAX` would defeat the guard rather than disable it, since the
addition would wrap before the test could catch it.

**4. AMD's `W` array carries three meanings at once, and we have not evaluated whether to.** For
an element: zero means absorbed, permanently; at or above `wflg` the excess is the set size; below
`wflg` but nonzero means live and not yet seen this scan. For a variable: equal to `wflg` means in
the pattern. Hence `clear_flag` collapses everything nonzero to 1 rather than to 0, so absorption
survives the sweep, and restarts `wflg` at 2. genmmd does a smaller version of the same thing,
getting two logical marks from one array by offsetting the stamp, `mt = tag + md0`. Ours holds
stamps only, with `mCliqueSize[c] = 0` as a separate absorption sentinel and `outside[c]` as a
separate counter. This is the same observation as item 1 seen from the other side, and the
trade is not obviously in either direction: theirs is denser, ours cannot suffer this class of
aliasing bug and is far easier to reason about.

**5. A stale claim in the experiment README.** It says amd1 obtains `|C[c] - C[p]|` by walking each
touched clique's members while amd2 obtains it by subtraction. Production `Amd1.cpp` does the
subtraction, same as `Amd2.cpp`. Either the prototype and production have diverged here or the
sentence predates a change to one of them. Worth checking before that paragraph is trusted again,
since a 3.7x work ratio is quoted from it.

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
