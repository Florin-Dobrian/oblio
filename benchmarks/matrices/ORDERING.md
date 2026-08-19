# Oblio: Ordering

This report asks what the ORDERING STEP ALONE costs and what space it takes, on 246 matrices from
the SuiteSparse Matrix Collection that Oblio's authors did not choose and did not generate.
`ACCURACY.md` beside it asks how good the answers are and `PERFORMANCE.md` what a whole solve
costs; this one isolates the first phase, where no values are read and no definiteness is needed,
so the input set is everything on disk rather than the subsets those two require.

**BOTH BRANCHES ARE HERE, 2026-08-19.** `MMD3` against `MMD` and `AMD3` against `AMD`, over
the same 246 matrices, in two runs of the same program: `make mmdorder` and `make amdorder`.

The comparison is unusually clean on both. **Each of ours returns the vendored routine's exact
permutation on every matrix in the set**, checked per matrix inside the run rather than assumed, so
there is no fill to trade against time: one `nnz(L)` column serves both routines in each table, and
time is the only thing that differs.

That check is new and it earned itself. It found FOUR divergences on the amd side, one of them a
correctness bug in supervariable detection that had been silently costing fill; see
`docs/NEXT.md` (2026-08-18). The mmd side needed nothing.

**The short answers, one line each.**

| | vendored | ours | ratio | ours slower on |
|---|---:|---:|---:|---:|
| mmd | 75.9 s | 52.9 s | **0.70** | 156 of 237 |
| amd | 3.52 s | 4.12 s | **1.17** | 162 of 242 |

Those two rows are not comparable to each other and the reasons are in the next section. Within
each, "slower on most matrices and faster overall" is the mmd shape and reconciling it is most of
what this report is about; the amd branch is slower on both counts and honestly so.

## The one number that reframes the rest

**MMD costs 21.6 times what AMD costs on this set, for 1.6 per cent more fill.** 75.9
seconds against 3.52, and 18.41 billion factor entries against 18.11 billion.

That is a property of the two ALGORITHMS rather than of either implementation, and it dwarfs every
implementation ratio in this report. Minimum degree recomputes exact degrees; approximate minimum
degree bounds them and never looks at a clique twice. On a set like this the bound costs almost
nothing in quality and saves an order of magnitude in time.

It also means the two halves of this report are asking different questions. The mmd table is about
a 76-second budget where a 30 per cent saving is 23 seconds; the amd table is about a 3.5-second
budget where a 17 per cent loss is 0.60 seconds. Read them accordingly.

### Which wins, counted matrix by matrix

The 20.9x is a total and the totals are the giant matrices. Counting instead how often each branch
wins, one vote per matrix, all 246:

| | mmd wins | amd wins | equal |
|---|---:|---:|---:|
| **time, vendored**, `MMD` against `AMD` | 86 | **158** | 2 |
| **time, ours**, `MMD3` against `AMD3` | 50 | **196** | 0 |
| **fill** | 92 | **136** | 18 |

**Amd wins on every axis and by a clear margin, and the margin is wider for our pair than for the
vendored one**, 196 to 50 against 158 to 86. That is the two implementation results seen from
another angle: `MMD3` beats `MMD` mostly on the giant matrices, which are few, while `AMD3` loses to
`AMD` fairly evenly, so on a one-vote-per-matrix count our amd pulls further ahead.

**The fill row is one column, not two.** Each of ours returns its reference's exact permutation, so
`MMD3` and `MMD` fill identically and so do `AMD3` and `AMD`; the 92 against 136 is a statement
about the two ALGORITHMS. And it is closer than the 1.6 per cent total suggests: the median matrix
has `mmd / amd` fill of 1.003, effectively a tie, with p10 0.966 and p90 1.056 and a long tail
either way, 0.688 at best for mmd and 2.293 at worst. Eighteen matrices tie exactly.

## The environment

| | |
|---|---|
| Machine | Apple M4, 32 GB, macOS 26.6.1 |
| Compiler | Apple Clang (Xcode 26.6), C++17, `-O3 -DNDEBUG` |
| References | `MMD` is genmmd, the Fortran multiple-minimum-degree of George and Liu, translated |
| | `AMD` is `amd_order` from SuiteSparse AMD 3.3.4, Amestoy, Davis and Duff |
| Naming | `MMD` and `AMD` are those two throughout; `MMD3` and `AMD3` are ours. The routine |
| | names appear only in this table, so every comparison below reads the same way |
| Timing | best of N after a warm-up, both routines through the same helper with the same N |

Both routines are timed by the same path on purpose. A column timed one way against a column timed
another differed 2.4 percent on 2026-08-10 and the difference was briefly read as a result.

## What the columns mean

```
n   m   tril(A)   A+I   cC   pC   nnz(L)   cC/tril   cC/nnzL   pC/cC   MMD ms   MMD3 ms   MMD3/MMD
```

`n` and `m` are vertices and edges, `m` counted from the pattern rather than assumed, a Matrix
Market file not necessarily carrying every diagonal entry. **`tril(A) = n + m`** is A in its stored
form. **`A+I = 2m`** is what the quotient graph's source pool holds, each edge appearing in both
endpoints' lists. **`nnz(L)`** includes the diagonal.

**THE TWO CLIQUE FIGURES ANSWER DIFFERENT QUESTIONS, and the words are chosen to keep them apart.**

- **`cC`, CUMULATIVE.** Every entry the arena ever handed out. Nothing is reclaimed in this layout,
  so it holds dead cliques and the members that contractions dropped out of live ones. It is what
  the flat clique store actually costs.
- **`pC`, PEAK LIVE MEMBERS.** The most that were alive at any one instant. It is what a CHUNKED
  store, one allocation per clique, would ask the allocator for at its worst moment, and it is
  PAYLOAD ONLY: no per-clique header, no allocator rounding, no capacity above size.

`pC` is a property of the ALGORITHM rather than of the layout. Two drivers running the same method
form the same cliques and lose the same members at the same moments whatever their storage, so the
figure is identical for `Mmd3`, `Mmd3B` and `Mmd3C`, and likewise for `Amd3` and `Amd3B`. Neither
vendored routine can report one. Only `cC` is ours.

Three ratios. `cC/tril(A)` says whether the arena tracked the INPUT. `cC/nnz(L)` says how much the
compression bought. **`pC/cC` is the one that decides anything**: the fraction a chunked store would
hold at its worst, and so the saving that scheme would buy.

**Read `pC/cC` conditioned on `cC/nnzL`, or it will mislead.** Where the arena is one per cent of
the factor, clique storage is noise beside what the factorization will need and no clique scheme is
worth building for it whatever `pC/cC` says. Where the arena is half the factor, it decides.

### Each table carries its B sibling too, and only its time

The last two columns of each table are the branch's B sibling: `MMD3B ms` and `MMD3B/MMD` on the
mmd side, `AMD3B ms` and `AMD3B/AMD` on the amd side. A sibling is its branch's driver on the
VENDORED ROUTINE'S OWN CLIQUE STORAGE rather than on our arena, so **it is the layout comparison**:
same algorithm, same encodings, same C++, different clique store. The section "The layout question"
below sets the three stores out, adds the chunked one that none of them is, and reads the result.

**Nothing else about a sibling is printed, and that is deliberate.** Its clique storage is a
different quantity from `cC` in both cases, a compacted workspace on one side and chained dead
segments on the other, so there is nothing for that column to compare. Everything else about it
MUST equal its driver's, the two being the same algorithm, so the order, the fill and `pC` are
CHECKED rather than shown and a mismatch is flagged at the end of the row: `AMD3B order differs`,
`MMD3B pC differs`, and so on.

**That last check is the sharp one, and it is why the columns exist.** `pC` is a property of the
algorithm and not of the layout, so two drivers that agree on the permutation can still be caught
doing different work. Order and fill compare the ANSWER, and a B layer exists to differ in
MECHANISM while agreeing on the answer, so until this check there was nothing watching the thing it
varies. It found two defects in `Amd3B`'s mid-walk garbage collector within an hour of existing,
neither of which moved a permutation or a fill figure and neither of which the digest, the vendored
acceptance checks, `test_order` or the sanitizers could see; see `docs/NEXT.md` (2026-08-19). Every
marker was clear across the 246 on both branches afterwards.

On the amd table the same columns carry `AMD ms`, `AMD3 ms` and `AMD3/AMD`, with one caveat that has
no counterpart on the mmd side: **the AMD clock is on the whole vendored call**, which is what a
caller
pays and includes `AMD_aat` and `AMD_postorder`, neither of which we do. A phase split puts those
plus argument validation at about 13 per cent of the vendored column, so against the comparable
region alone the amd ratio is nearer 1.3 than the 1.15 the table shows. Everything except the clock
comes from the hooked pre-postorder copy, since `AMD3` returns that order. MMD does no postorder
at all, so its table needs no such copy and carries no such caveat.

## Where the time actually goes

The set spans five orders of magnitude in ordering time, so a count of wins and losses says almost
nothing. Grouped by how large the arena grew, `cC/tril(A)`. **The AGGREGATED ratio is the band's
total ours divided by its total vendored**, so the big matrices dominate it; the MEDIAN ratio
weights every matrix in the band equally. The two disagree wherever a band mixes sizes, and that
disagreement is the report's main finding on the mmd side:

**mmd, `MMD` against `MMD3`**

| `cC/tril(A)` | matrices | share of total time | `MMD3/MMD` | `MMD3B/MMD` |
|---|---|---|---|---|
| 0, no off-diagonal at all | 14 | 0.0% | 2.24 | 2.20 |
| 0 to 0.5 | 52 | 0.3% | 0.79 | 0.88 |
| 0.5 to 1.0 | 64 | 0.4% | 0.99 | 0.94 |
| 1.0 to 2.0 | 79 | 1.7% | 0.83 | 0.94 |
| 2.0 to 5.0 | 25 | 1.6% | **0.73** | 1.13 |
| 5.0 and above | 12 | **96.0%** | **0.69** | 1.09 |

**The twelve hardest matrices are 96 per cent of all the work, and we are 31 per cent faster on
them.** The 156 matrices where we are slower cost **28.5 milliseconds in total, all of them
together**; the 81 where we are faster save **23.1 seconds**. The largest single loss anywhere is
2.75 ms, on `GHS_indef/bloweybl`.

**AND `MMD3B` LOSES EXACTLY WHERE `MMD3` WINS**, 1.09 against 0.69 in the band holding 96 per cent
of the time. That is chaining's price and it is the mirror of the amd result; the layout section
below reads both.

**amd, `AMD` against `AMD3`**

| `cC/tril(A)` | matrices | share of total time | `AMD3/AMD` | `AMD3B/AMD` |
|---|---|---|---|---|
| 0, no off-diagonal at all | 13 | 0.0% | 0.77 | 1.25 |
| 0 to 0.5 | 71 | 2.2% | 0.87 | **0.78** |
| 0.5 to 1.0 | 69 | 5.4% | 1.09 | **0.99** |
| 1.0 to 2.0 | 72 | 6.2% | 1.14 | **1.05** |
| 2.0 to 5.0 | 19 | **73.9%** | 1.18 | **1.12** |
| 5.0 and above | 2 | 12.3% | 1.20 | 1.20 |

Both ratios are aggregated, and `AMD3B` is there because it answers the layout question; the
section after the top-ten tables reads it.

**The amd branch has no such reversal, and the shape is the opposite one.** We are faster on the
small end, where the dense-row rule and the empty-row prepass now do what `AMD` does and our
container layer has little to be slow at, and slower on the large end where 86 per cent of the time
is. Read against the comparable region rather than the whole vendored call, the aggregated figures
move up by roughly a seventh.

**The ten most expensive matrices, each branch.** They are 95 per cent of the mmd total and 85 per
cent of the amd one.

| mmd | `MMD` | `MMD3/MMD` | `MMD3B/MMD` | `cC/tril` | `pC/cC` |
|---|---:|---:|---:|---:|---:|
| FlowIPM22/uni_chimera_i5 | 15.77 s | **0.40** | 1.12 | 52.71 | 0.01 |
| PARSEC/Ga41As41H72 | 13.72 s | 0.84 | 1.04 | 5.78 | 0.02 |
| PARSEC/Si87H76 | 11.89 s | 0.83 | 1.05 | 7.87 | 0.02 |
| JGD_Trefethen/Trefethen_20000 | 10.37 s | 0.75 | 1.12 | 7.27 | 0.07 |
| FlowIPM22/uni_chimera_i1 | 10.04 s | **0.59** | 1.07 | 48.92 | 0.01 |
| JGD_Trefethen/Trefethen_20000b | 6.50 s | 0.85 | 1.16 | 24.14 | 0.03 |
| Pajek/Reuters911 | 1.33 s | 0.80 | 1.14 | 10.07 | 0.03 |
| DIMACS10/vsp_befref_fxm_2_4_air02 | 0.99 s | 0.67 | 1.10 | 9.58 | 0.08 |
| Schenk/nlpkkt80 | 0.74 s | 0.79 | **0.91** | 1.42 | 0.61 |
| Andrews/Andrews | 0.73 s | 0.66 | 1.06 | 5.65 | 0.06 |

| amd | `AMD` | `AMD3/AMD` | `AMD3B/AMD` | `cC/tril` | `pC/cC` |
|---|---:|---:|---:|---:|---:|
| Schenk/nlpkkt80 | 0.94 s | 1.17 | 1.18 | 2.18 | 0.38 |
| PARSEC/Ga41As41H72 | 0.76 s | 1.19 | **1.09** | 2.05 | 0.05 |
| PARSEC/Si87H76 | 0.59 s | 1.23 | **1.11** | 2.91 | 0.06 |
| FlowIPM22/uni_chimera_i1 | 0.37 s | 1.21 | 1.20 | 15.16 | 0.03 |
| FlowIPM22/uni_chimera_i5 | 0.07 s | 1.11 | 1.11 | 4.37 | 0.09 |
| Andrews/Andrews | 0.07 s | 1.14 | 1.09 | 3.85 | 0.09 |
| JGD_Trefethen/Trefethen_20000b | 0.06 s | 1.18 | 1.19 | 6.08 | 0.10 |
| FlowIPM22/uni_chimera_i2 | 0.06 s | 1.04 | **0.99** | 2.40 | 0.16 |
| JGD_Trefethen/Trefethen_20000 | 0.04 s | 1.10 | **1.02** | 3.27 | 0.15 |
| Rothberg/cfd1 | 0.03 s | 0.95 | 0.93 | 0.81 | 0.24 |

**The two lists are almost the same matrices at wildly different costs.** `uni_chimera_i5` is
15.8 seconds of `MMD` and 0.07 of `AMD`, a factor of 230. That is the 21.6x above, seen one
row at a time.

## The layout question: three clique stores, measured

**WHAT THE TWO ACTUALLY ARE, because the ratio means nothing without it.** Four axes, and they
agree on only one of them:

| | clique space | shape | bounded | reclamation |
|---|---|---|---|---|
| `AMD3`, `MMD3` | separate arena | FLAT | no, it grows | none, ever |
| `AMD3B`, `AMD` | unified with A and I | FLAT | yes, plus a fifth | compaction |
| `MMD3B`, `MMD` | dead segments, chained | FLAT | yes, exactly nnz | reuse by chaining |
| a chunked store | one allocation per clique | CHUNKED | no | free, on death |

So our arena buys simplicity with memory: no collector to write, no bound to respect, and a store
that on this set reaches 98.1 million entries on the amd side and 199.9 on the mmd, 15 and 53 times
`tril(A)` on their worst rows. The two vendored schemes buy a hard bound with machinery, and with
DIFFERENT machinery, which is the whole point of having both siblings. `AMD` compacts a workspace
sized at a fifth over the pattern; `MMD` chains a clique through the dead segments of the cliques
it absorbed and needs no spare space at all.

**FLAT IS WHY COMPACTION IS THE ONLY WAY TO RECLAIM, and that is the axis the fourth row is on.**
In a flat store a clique must be one contiguous block, so the space a dead one leaves is a hole
that nothing can use until everything live is slid down over it; the collector is a consequence of
the shape, not a separate design choice. Give each clique its own allocation and reclamation costs
nothing at all: a dead clique is freed where it lies, no bound is needed, and no collector exists.
What is paid instead is the allocator, a header per clique and rounding, on every clique ever
formed.

**That fourth row is NOT BUILT**, and this report is where the case for it would come from. The
`pC/cC` column says what it would save, and the section on the arena below reads it: over the amd
set a chunked store would hold 0.19 of what the flat arena does, and on the mmd side 0.10, with the
worst rows at 0.01. `docs/DESIGN_DECISIONS.md` (2026-08-18) sets out the two axes all four rows sit
on, memory against machinery, and names the four quotient graph classes;
`experiments/ordering/README.md` has the mechanism.

**The comparison below is therefore not "is compaction expensive".** It is the whole bargain at
once: the collector's cost, minus whatever the unified layout gains by keeping a vertex's cliques
and its lists in one array, against an arena that pays nothing to reclaim because it never does.
Over the 246:

| | total | against its vendored | against our arena |
|---|---:|---:|---:|
| `AMD` | 3.52 s | | |
| `AMD3`, arena | 4.12 s | 1.17 | |
| `AMD3B`, compacted pool | 3.92 s | **1.11** | **0.95** |
| `MMD` | 75.9 s | | |
| `MMD3`, arena | 52.9 s | 0.70 | |
| `MMD3B`, chained | 82.5 s | **1.09** | **1.56** |

**THE TWO SIBLINGS SAY OPPOSITE THINGS, AND THAT IS THE RESULT.** Compaction is nearly free:
`AMD3B` beats the arena by 5 per cent over the set and on 224 of 242 matrices, median ratio 0.90,
winning in every band that has work in it. Chaining is not: `MMD3B` costs 56 per cent over the
arena and wins on only 60 of 237.

**And chaining loses precisely where it matters.** In the band holding 96 per cent of the mmd time,
`MMD3` reads 0.69 against `MMD` and `MMD3B` reads 1.09, a 31 per cent win turned into a 9 per cent
loss. The grid ladders understate this badly, reading 1.1 to 1.3, because grids have short cliques
and few links to follow.

**The reason is in the two mechanisms and it is not a surprise.** A compaction is a rare sweep: at
a fifth of elbow room it runs about once for a whole ordering, so its cost is amortised to nearly
nothing. A chain is a link test on EVERY READ of every clique, forever, and no amount of headroom
reduces it because chaining exists precisely to need none. That asymmetry was written down as
reasoning in `experiments/ordering/README.md` before it was measured; this is the measurement.

**Read both as directions and not as coefficients.** Each sibling differs from its driver in more
than the clique store: `Amd3B` also carries `AMD`'s run order, the three-move rotation and
positions in the pool walks; `Mmd3B` carries genmmd's list conventions. The clique store is the
largest of those differences and the one each file exists to price, but it is not the only one, and
a clean measurement would be the same code with only the placement rule swapped. See
`experiments/ordering/README.md` on the two axes, memory against machinery.

**The bands with no work in them say nothing.** The purely diagonal matrices read 1.25 and 2.20 for
the two siblings; there is no ordering to do there and the figure is the cost of standing up a
store that is then never used.

## Where we lose, and why it does not matter

**The 14 matrices with no off-diagonal entries at all.** `Boeing/bcsstm39`, `HB/bcsstm25`,
`Cunningham/m3plates` and eleven others are purely diagonal: the ordering has nothing to do and the
permutation is the identity. On the mmd side the ratio is 2.15 and the difference is a few hundred
microseconds for the whole group. There is no result here. A routine that does nothing is measuring
its own entry and exit, and ours has a `std::vector` layer around a quotient graph it then never
uses.

**On the amd side that band now reads 0.75, ours faster**, because the empty-row prepass added on
2026-08-18 numbers those vertices where they stand instead of filing and popping them from a degree
bucket. It is the same non-result in the other direction and should be read the same way.

**The nearly-trivial rest.** `Bai/mhd3200b` at 1.78, `HB/bcsstm25` at 2.10, `GHS_indef/spmsrtls` at
1.79: sub-millisecond orderings on matrices with almost no fill. The pattern across the whole set
is that our fixed overhead, a container layer measured at 17 to 20 percent plus construction of a
quotient graph that a trivial matrix barely touches, is a constant, and a constant is most visible
where the work is smallest.

**That is the honest shape of it and it is not being minimised.** If the question were "order ten
thousand tiny matrices", `MMD` would win and we would need to answer for it. It is not the
question this benchmark asks.

## The arena: what it costs, and what a chunked store would cost instead

`cC` is what our storage costs, and the reason it is worth watching is that ours is the one scheme
here NOT bounded from the input: `MMD` chains dead segments and stays inside `O(n + m)`, `AMD`
compacts a fixed workspace, and we allocate cliques as they form and never reclaim the space. See
`docs/DESIGN_DECISIONS.md` (2026-08-18) for the two axes that frames, and why the trade is
deliberate.

**The compression is often enormous.** `PARSEC/Ga41As41H72` factors to 7.08 billion entries and the
mmd arena holds 54.2 million, `cC/nnz(L) = 0.008`. `Schenk/nlpkkt80`, a million rows, factors to
3.68 billion and the arena is 21.2 million, 0.006. A clique holds the update part of a
supervariable's column and nothing else, and mass elimination gives one clique per supervariable
rather than one per column, so what we store tracks the input far more closely than the factor.

**And sometimes it is not.** `Pajek/Reuters911` has `cC/nnz(L) = 0.723` on the mmd side: the arena
is nearly three-quarters of the factor and 10.07 times `tril(A)`. `SNAP/email-Enron` is 0.512.
`FlowIPM22/uni_chimera_i5` reaches `cC/tril(A) = 52.71`, an arena fifty times the input it came
from. These are the cases the guarantee argument is about: social graphs and interior-point systems
where supervariables barely form, the compression fails, and our space is neither small nor
predictable.

**The striking thing is that those are also where the mmd branch is fastest.** `uni_chimera_i5` at
0.43, `uni_chimera_i1` at 0.59, `Reuters911` at 0.81, `vsp_befref` at 0.71. Spending memory to avoid
chaining is exactly the trade being made, and this set says it pays where it costs the most.

**AMD'S ARENA IS HALF THE SIZE FOR THE SAME MATRICES**, 98.1 million entries against mmd's 199.9
across the set, and the worst cases are far milder: `uni_chimera_i5` is 4.37 times `tril(A)` where
mmd is 52.71, `Reuters911` is 1.16 against 10.07. Fewer, larger cliques and aggressive absorption
are why.

### What `pC` says, and it is the new question

`pC/cC` is the fraction a CHUNKED clique store would hold at its worst moment. Over the whole set it
is **0.10 on the mmd side and 0.19 on the amd side**; per matrix the medians are 0.26 and 0.20.

**The extremes are where it matters and they are all on the mmd side.** `uni_chimera_i1` and
`uni_chimera_i5` read 0.01, `Reuters911` and `email-Enron` 0.03, `Trefethen_20000b` 0.03. Those are
exactly the matrices where `cC/tril(A)` is 10 to 50. So the arena's worst behaviour is almost
entirely dead space: a store that freed a clique when it died would hold **one to three per cent**
of what the flat one holds, on the rows that hurt.

**But read it against `cC/nnzL` before concluding anything.** On those same rows the arena is 1 to 9
per cent of the factor, so clique storage was never the binding constraint there; the factorization
was. The rows where clique storage genuinely competes with the factor, `Reuters911` at 0.723 and
`email-Enron` at 0.512, are also rows with `pC/cC` of 0.03, which is the case FOR building it. The
rows where `pC/cC` is high, `nlpkkt80` at 0.61 and `LF10000` at 0.73, are rows where little would be
saved.

That is the whole argument for a chunked clique store, and it is now measured rather than assumed.

## Reading the whole table

Three habits, all learned by getting them wrong.

**Do not average the ratio column.** The mean over 246 rows is dominated by matrices that take
microseconds. Read the aggregated ratio, which is total against total, or a median within a band.

**A ratio hides which side moved.** Where a row looks anomalous, convert to time per vertex or per
edge before believing it. That failure cost a full day on the grid benchmark, where a pattern read
as ours turned out to be an artifact of the vendored routine.

**Do not compare a number in the mmd table against one in the amd table except deliberately.** The
denominators are different routines doing different amounts of work, and the amd one includes two
phases we do not run. Cross-branch statements in this report are confined to the section near the
top that makes them on purpose.

## What this report does not establish

**It is two orderings against two references, and nothing else.** Nested dissection is absent from
this tree, and on the largest matrices here it is what a production solver would actually use.

**It says almost nothing about ordering QUALITY.** Each of ours reproduces its reference's
permutation exactly, so within a branch fill is identical by construction. The one quality
statement this set supports is the cross-branch one near the top: amd fills less on 136 matrices,
mmd on 92, 18 tie, and the median matrix is a tie to within 0.3 per cent. `PERFORMANCE.md` looks at
what that does to a whole solve.

**The AMD clock is not the amd branch's cost.** It is the whole vendored call, about 13 per cent of
which
is work we do not do. Every amd ratio here is therefore optimistic by roughly that much, and the
figure to quote for the C++ and container cost is not in this report; `benchmarks/ordering`'s
`make phases2d` has the split.

**And it says nothing about what a full solve costs.** Ordering is one phase and usually a small
one; see `PERFORMANCE.md` for the split.

## Reproducing this

```
cd benchmarks/matrices
python3 ssget.py            # fetch the collection subset into data/
make mmdorder               # the mmd table, all 246 files
make amdorder               # the amd table, the same 246
```

The private vendored sources are needed for the vendored columns and, on the amd side, for the
hooked pre-postorder copy the fill column comes from. Without them either table still runs and
prints ours alone.

**A row ending in `order differs` or `raw order differs` means our permutation and the vendored one
parted on that matrix**, and the shared `nnz(L)` column has stopped being both routines'. Both
tables were clean on 2026-08-19. That check is what found the amd divergences; before it existed,
the equality was asserted in this file's opening and verified only on generated shapes.
