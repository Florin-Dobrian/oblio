# Oblio: Ordering

This report asks what the ORDERING STEP ALONE costs and what space it takes, on 246 matrices from
the SuiteSparse Matrix Collection that Oblio's authors did not choose and did not generate.
`ACCURACY.md` beside it asks how good the answers are and `PERFORMANCE.md` what a whole solve costs;
this one isolates the first phase, where no values are read and no definiteness is needed, so the
input set is everything on disk rather than the subsets those two require.

**MEASURED 2026-08-21.** Every figure below comes from two runs of one program on that day,
`make mmdorder` and `make amdorder`, over the same 246 matrices. This report is a snapshot: it says
what is true of the code as it stands, not how it got there.

**BOTH BRANCHES ARE HERE.** `MmdFlat` against `MmdVendored` and `AmdFlat` against `AmdVendored`,
with each branch's alternative-store siblings alongside.

**A DRIVER IS ITS BRANCH PLUS THE CLIQUE STORE IT RUNS ON.** `MmdVendored` and `AmdVendored` are
the two reference codes; `MmdFlat` and `AmdFlat` are ours on our own arena; `MmdCompacted` and
`AmdCompacted` are ours on `AMD_2`'s compacted workspace; `MmdChained` is ours on genmmd's chained
segments. Per-matrix tables print three-letter tags, `Vnd Flt Chn Com`, with a legend line above
them; this report uses the full names.

**THREE LAYOUT SIBLINGS RIDE ALONG.** `MmdChained` and `MmdCompacted` on the mmd side,
`AmdCompacted` on the amd side. Each is its branch's driver on a different clique store, so the
sibling columns ARE the layout comparison. Every driver here compiles its quotient graph into its
own translation unit, so a comparison between two of ours is between two things built the same way;
that was not true before 2026-08-19 and figures from before then are not comparable to these.

The comparison is unusually clean on both. **Each of ours returns the vendored routine's exact
permutation on every matrix in the set**, checked per matrix inside the run rather than assumed, so
there is no fill to trade against time: one `nnz(L)` column serves both routines in each table, and
time is the only thing that differs.

That check is new and it earned itself. It found FOUR divergences on the amd side, one of them a
correctness bug in supervariable detection that had been silently costing fill; see `docs/NEXT.md`
(2026-08-18). The mmd side needed nothing. **In these runs not one marker fired**, on either branch
or any of the three siblings.

**The short answers, one line each.**

| | vendored | ours | ratio | ours slower on |
|---|---:|---:|---:|---:|
| mmd | 75.5 s | 49.1 s | **0.65** | 126 of 228 |
| amd | 3.55 s | 4.02 s | **1.13** | 142 of 233 |

The default is neither of those two. `AmdCompacted` orders the set in **3.96 s**, and it is what a
caller who sets nothing now gets; see `include/oblio/OrderEngine.h`.

Those two rows are not comparable to each other and the reasons are in the next section. Within
each, "slower on most matrices and faster overall" is the mmd shape and reconciling it is most of
what this report is about; the amd branch is slower on both counts and honestly so.

## The one number that reframes the rest

**MMD costs 21.3 times what AMD costs on this set, for 1.6 per cent more fill.** 75.5 seconds
against 3.55, and 18.41 billion factor entries against 18.11 billion. Read the direction carefully:
mmd is the one that costs more, on BOTH axes: it is not a trade of time against quality.

That is a property of the two ALGORITHMS rather than of either implementation, and it dwarfs every
implementation ratio in this report. Minimum degree recomputes exact degrees; approximate minimum
degree bounds them and never looks at a clique twice. On a set like this the bound costs almost
nothing in quality and saves an order of magnitude in time.

It also means the two halves of this report are asking different questions. The mmd table is about a
75-second budget where a 35 per cent saving is 26 seconds; the amd table is about a 3.6-second
budget where a 13 per cent loss is 0.47 seconds. Read them accordingly.

### Which wins, counted matrix by matrix

The 21.3x is a total and the totals are the giant matrices. Counting instead how often each branch
wins, one vote per matrix, all 246:

| | mmd wins | amd wins | equal |
|---|---:|---:|---:|
| **time, vendored**, `MmdVendored` against `AmdVendored` | 86 | **159** | 1 |
| **time, ours**, `MmdFlat` against `AmdFlat` | 52 | **190** | 4 |
| **time, the defaults**, `MmdFlat` against `AmdCompacted` | 41 | **201** | 4 |
| **fill** | 92 | **136** | 18 |

**Amd wins on every axis and by a clear margin, and the margin is wider for our pair than for the
vendored one**, 190 to 52 against 159 to 86. That is the two implementation results seen from
another angle: `MmdFlat` beats `MmdVendored` mostly on the giant matrices, which are few, while
`AmdFlat` loses to `AmdVendored` fairly evenly, so on a one-vote-per-matrix count our amd pulls
further ahead. Against `AmdCompacted` it is 201 to 41.

**The fill row is one column, not two.** Each of ours returns its reference's exact permutation, so
every mmd driver fills identically and so does every amd one; the 92 against 136 is a statement
about the two ALGORITHMS. And it is closer than the 1.6 per cent total suggests: the median matrix
has `mmd / amd` fill of 1.003, effectively a tie, with p10 0.967 and p90 1.057 and a long tail
either way, 0.688 at best for mmd and 2.293 at worst. Eighteen matrices tie exactly.

## The environment

| | |
|---|---|
| Machine | MacBook Pro (Mac16,1), Apple M4, 32 GB, macOS 26.6.2 |
| Compiler | Apple Clang (Xcode 26.6), C++17, `-O3 -DNDEBUG` |
| References | `MmdVendored` is genmmd, the multiple-minimum-degree of George and Liu, translated |
| | `AmdVendored` is `amd_order` from SuiteSparse AMD 3.3.4, Amestoy, Davis and Duff |
| Naming | `MmdVendored` and `AmdVendored` are those two throughout; the other five are ours. |
| | The routine names appear only in this table, so every comparison reads the same way |
| Timing | best of N after a warm-up, both routines through the same helper with the same N |

Both routines are timed by the same path on purpose. A column timed one way against a column timed
another differed 2.4 percent on 2026-08-10 and the difference was briefly read as a result.

## The summary tables

Five statistics per ratio, and they answer different questions. **median** and **mean** weight every
matrix equally, so a microsecond matrix counts as much as a fifteen-second one. **p10** and **p90**
are the tenth and ninetieth percentiles and say how wide the spread is. **total** is the sum of one
column's milliseconds over the whole set divided by the sum of the other's, which is what ordering
all 246 would actually cost, so it is dominated by the expensive ones.

Where the median and the total disagree, the winner differs between the cheap matrices and the
costly ones. On the mmd side they disagree sharply and that is this report's main finding.

A row is dropped from a ratio when either side is under 0.05 ms, a ratio of two readings that small
being noise; `n` is how many of the 246 survive that.

**mmd**

| | median | mean | p10 | p90 | total | n |
|---|---:|---:|---:|---:|---:|---:|
| `MmdFlat / MmdVendored` | 1.042 | 1.053 | 0.776 | 1.340 | **0.651** | 228 |
| `MmdChained / MmdVendored` | 1.219 | 1.251 | 1.069 | 1.437 | 1.103 | 228 |
| `MmdCompacted / MmdVendored` | 1.038 | 1.047 | 0.782 | 1.331 | **0.660** | 228 |
| `MmdCompacted / MmdFlat` | **0.998** | 0.997 | 0.946 | 1.045 | **1.014** | 234 |
| `MmdChained / MmdFlat` | **1.154** | 1.221 | 1.042 | 1.398 | 1.695 | 234 |

**amd**

| | median | mean | p10 | p90 | total | n |
|---|---:|---:|---:|---:|---:|---:|
| `AmdFlat / AmdVendored` | 1.037 | 1.051 | 0.766 | 1.329 | 1.132 | 233 |
| `AmdCompacted / AmdVendored` | **0.996** | 1.004 | 0.743 | 1.266 | 1.115 | 233 |
| `AmdCompacted / AmdFlat` | **0.957** | 0.956 | 0.892 | 1.015 | **0.986** | 235 |

**Four things to read off them.**

**The mmd median and total disagree completely**, 1.042 against 0.651. We lose slightly on most
matrices and win by a third on the expensive ones. The amd branch does not do this, 1.037 against
1.132, so it loses a little everywhere and slightly more where it costs.

**`MmdCompacted / MmdFlat` IS A WASH.** 0.998 at the median, 1.014 in total, p10 to p90 spanning
0.946 to 1.045. On the mmd branch the clique store makes no measurable difference to time in either
direction, and the spread is wider than the effect, so a single matrix says nothing about it.

**`AmdCompacted / AmdFlat` IS A REAL WIN AND A SMALL ONE.** 0.957 at the median, winning on 188 of
235 matrices, with p10 to p90 of 0.892 to 1.015. Its total is 0.986, so the aggregate saving over
the whole set is under two per cent while the per-matrix win is nearer four: the count and the
median carry this result, not the total.

**The p10 to p90 spread is wide against the vendored routines and narrow between our own drivers.**
0.776 to 1.340 for `MmdFlat / MmdVendored` against 0.946 to 1.045 for `MmdCompacted / MmdFlat`.
Comparing two of our drivers is a much sharper instrument than comparing ours to a reference, which
mixes in every difference between two independently written codes, most of which nobody has
enumerated. The translation-unit work of 2026-08-19 reached the same conclusion from the other
direction: aligning our drivers with each other made our own ratios comparable and moved the ratios
against the references not at all. The narrow rows are the ones to reason from.

**`AmdCompacted / AmdVendored` reaches 0.996 at the median.** The compacted store takes the amd
branch to parity with `AMD_2` on a typical matrix, from 1.037. It sits four thousandths below 1.0,
which is parity and not a win.

## What the columns mean

```
n   m   tril(A)   A+I   cC   pC   nnz(L)   cC/tril   cC/nnzL   pC/cC   Vnd ms   Flt ms   Flt/Vnd
```

**THE DRIVER COLUMNS ARE TAGGED AND THE PROGRAM PRINTS A LEGEND LINE ABOVE THE TABLE.** `Vnd` the
vendored routine, `Flt` our arena, `Chn` the chained store, `Com` the compacted one, and on the amd
side `Raw` the hooked pre-postorder copy. A run is mmd or amd and never both, so the branch is said
once in that line rather than in every heading. The tags are the 2026-08-21 naming entry's and they
live here, in the same legend as `cC` and `pC`; prose and the summary tables above use the full
names.

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
figure is identical for `MmdFlat`, `MmdChained` and `MmdCompacted`, and likewise for `AmdFlat`
and `AmdCompacted`. Neither
vendored routine can report one. Only `cC` is ours.

Three ratios. `cC/tril(A)` says whether the arena tracked the INPUT. `cC/nnz(L)` says how much the
compression bought. **`pC/cC` is the one that decides anything**: the fraction a chunked store would
hold at its worst, and so the saving that scheme would buy.

**AND TWO COMPACTION COUNTS AT THE END OF EACH ROW.** `Vnd cmp` on the amd table is `AMD_2`'s own
`Info[AMD_NCMPA]`, read from an extra untimed call; `Com cmp` on either table is ours. They are
directly comparable because the three workspaces are the same size: `amd_order` allocates `nzaat +
nzaat/5 + 7n` and carves six per-vertex arrays off the front, leaving `iwlen = nzaat + nzaat/5 + n`
for the pool, and that is what both of our compacted siblings resize to. The `n` is `AMD_2`'s
enforced minimum elbow room and the `nzaat/5` is its recommendation; with neither the algorithm
still runs, and very slowly. `MmdChained` chains and never compacts, so it has no such column, and
neither does genmmd; the mmd table prints a dash where the amd table prints `Vnd cmp`.

**Read `pC/cC` conditioned on `cC/nnzL`, or it will mislead.** Where the arena is one per cent of
the factor, clique storage is noise beside what the factorization will need and no clique scheme is
worth building for it whatever `pC/cC` says. Where the arena is half the factor, it decides.

### Each table carries its siblings too, and only their times

The mmd table ends with `Chn ms`, `Chn/Vnd`, `Com ms`, `Com/Vnd` and `Com cmp`; the amd table with
`Com ms`, `Com/Vnd`, `Vnd cmp` and `Com cmp`. A sibling is its branch's driver on a
DIFFERENT CLIQUE STORE rather than on our arena, so **the sibling columns are the layout
comparison**: same algorithm, same encodings, same C++, different clique store. The section "The
layout question" below sets the stores out and reads the result.

**Nothing else about a sibling is printed except its compaction count, and that is deliberate.** Its
clique storage is a different quantity from `cC`, a compacted workspace in two cases and chained
dead segments in the third, so there is nothing for that column to compare. Everything else about it
MUST equal its driver's, the two being the same algorithm, so the order, the fill and `pC` are
CHECKED rather than shown and a mismatch is flagged at the end of the row: `AmdCompacted order
differs`, `MmdCompacted pC differs`, and so on. Those flags carry the full name rather than the tag,
being prose that fires rarely rather than a column.

**That last check is the sharp one, and it is why the columns exist.** `pC` is a property of the
algorithm and not of the layout, so two drivers that agree on the permutation can still be caught
doing different work. Order and fill compare the ANSWER, and a sibling exists to differ in MECHANISM
while agreeing on the answer, so until this check there was nothing watching the thing it varies. It
found two defects in `AmdCompacted`'s mid-walk collector within an hour of existing, neither of
which moved a permutation or a fill figure and neither of which the digest, the vendored acceptance
checks, `test_order` or the sanitizers could see; see `docs/NEXT.md` (2026-08-19).

**Every marker was clear across the 246 on both branches**, on all three siblings.

On the amd table there is one caveat with no counterpart on the mmd side: **the vendored clock is on
the whole vendored call**, which is what a caller pays and includes `AMD_aat` and `AMD_postorder`,
neither of which we do. A phase split puts those plus argument validation at about 13 per cent of
the vendored column, so against the comparable region alone the amd ratio is nearer 1.28 than the
1.13 the table shows. Everything except the clock comes from the hooked pre-postorder copy, since
`AmdFlat` returns that order. MMD does no postorder at all, so its table needs no such copy and
carries no such caveat.

## Where the time actually goes

The set spans five orders of magnitude in ordering time, so a count of wins and losses says almost
nothing. Grouped by how large the arena grew, `cC/tril(A)`. **The ratios below are AGGREGATED**, the
band's total ours divided by its total vendored, so the big matrices dominate them.

**mmd, against `MmdVendored`**

| `cC/tril(A)` | matrices | share of total time | `MmdFlat` | `MmdChained` | `MmdCompacted` |
|---|---:|---:|---:|---:|---:|
| 0, no off-diagonal at all | 14 | 0.0% | 2.17 | 2.22 | 2.20 |
| 0 to 0.5 | 53 | 0.3% | 0.79 | 1.02 | 0.79 |
| 0.5 to 1.0 | 63 | 0.4% | 0.91 | 1.12 | 0.91 |
| 1.0 to 2.0 | 81 | 2.2% | 0.75 | 1.10 | 0.73 |
| 2.0 to 5.0 | 23 | 1.1% | 0.76 | 1.08 | 0.78 |
| 5.0 and above | 12 | **96.0%** | **0.65** | 1.10 | **0.66** |

**The twelve hardest matrices are 96 per cent of all the work, and we are 35 per cent faster on
them.** The 144 matrices where we are slower cost **20.3 milliseconds in total, all of them
together**; the 101 where we are faster save **26.4 seconds**. The largest single loss anywhere is
2.17 ms, on `Schmid/thermal1`.

**AND `MmdChained` LOSES EXACTLY WHERE `MmdFlat` WINS**, 1.10 against 0.65 in the band holding 96
per cent of the time. That is chaining's price and it is the one store result that is not in doubt.

**`MmdCompacted` DOES NOT SEPARATE FROM `MmdFlat` IN ANY BAND.** It is ahead in the 1.0-to-2.0
band, 0.73 against 0.75, level in two, and behind by a hundredth in the two that carry the work,
0.78 against 0.76 and 0.66 against 0.65. Hundredths either way across six bands is what no effect
looks like.

**amd, against `AmdVendored`**

| `cC/tril(A)` | matrices | share of total time | `AmdFlat` | `AmdCompacted` |
|---|---:|---:|---:|---:|
| 0, no off-diagonal at all | 13 | 0.0% | 0.79 | **0.77** |
| 0 to 0.5 | 72 | 2.2% | 0.86 | **0.82** |
| 0.5 to 1.0 | 68 | 5.4% | 1.07 | **1.03** |
| 1.0 to 2.0 | 73 | 6.2% | 1.13 | **1.10** |
| 2.0 to 5.0 | 18 | **73.9%** | 1.14 | **1.12** |
| 5.0 and above | 2 | 12.3% | 1.19 | 1.20 |

**The amd store leads in five bands of six**, and the exception is the two-matrix top band, which
carries 12 per cent of the time and where the two are within a hundredth. That is the shape the mmd
side does not have.

**The amd branch has no such reversal, and the shape is the opposite one.** We are faster on the
small end, where the dense-row rule and the empty-row prepass now do what `AmdVendored` does and our
container layer has little to be slow at, and slower on the large end where 86 per cent of the time
is. Read against the comparable region rather than the whole vendored call, the aggregated figures
move up by roughly a seventh.

**The ten most expensive matrices, each branch.** They are 95 per cent of the mmd total and 85 per
cent of the amd one.

| mmd | `MmdVendored` | `MmdFlat` | `MmdChained` | `MmdCompacted` | `cC/tril` | `pC/cC` |
|---|---:|---:|---:|---:|---:|---:|
| FlowIPM22/uni_chimera_i5 | 15.91 s | **0.38** | 1.11 | **0.38** | 52.71 | 0.01 |
| PARSEC/Ga41As41H72 | 13.55 s | 0.79 | 1.09 | **0.78** | 5.78 | 0.02 |
| PARSEC/Si87H76 | 11.71 s | **0.77** | 1.07 | 0.78 | 7.87 | 0.02 |
| JGD_Trefethen/Trefethen_20000 | 10.29 s | **0.69** | 1.13 | 0.73 | 7.27 | 0.07 |
| FlowIPM22/uni_chimera_i1 | 10.00 s | **0.55** | 1.07 | 0.56 | 48.92 | 0.01 |
| JGD_Trefethen/Trefethen_20000b | 6.45 s | **0.80** | 1.17 | 0.83 | 24.14 | 0.03 |
| Pajek/Reuters911 | 1.28 s | **0.77** | 1.15 | 0.81 | 10.07 | 0.03 |
| DIMACS10/vsp_befref_fxm_2_4_air02 | 0.99 s | **0.64** | 1.15 | 0.67 | 9.58 | 0.08 |
| Andrews/Andrews | 0.78 s | 0.60 | 1.03 | 0.61 | 5.65 | 0.06 |
| Schenk/nlpkkt80 | 0.71 s | 0.81 | 1.05 | **0.78** | 1.42 | 0.61 |

| amd | `AmdVendored` | `AmdFlat` | `AmdCompacted` | `cC/tril` | `pC/cC` |
|---|---:|---:|---:|---:|---:|
| Schenk/nlpkkt80 | 0.96 s | **1.10** | 1.18 | 2.18 | 0.38 |
| PARSEC/Ga41As41H72 | 0.77 s | 1.15 | **1.06** | 2.05 | 0.05 |
| PARSEC/Si87H76 | 0.58 s | 1.21 | **1.11** | 2.91 | 0.06 |
| FlowIPM22/uni_chimera_i1 | 0.37 s | **1.19** | 1.20 | 15.16 | 0.03 |
| Andrews/Andrews | 0.07 s | 1.06 | **1.05** | 3.85 | 0.09 |
| FlowIPM22/uni_chimera_i5 | 0.07 s | **1.12** | 1.17 | 4.37 | 0.09 |
| JGD_Trefethen/Trefethen_20000b | 0.06 s | **1.18** | 1.19 | 6.08 | 0.10 |
| FlowIPM22/uni_chimera_i2 | 0.05 s | **1.02** | 1.08 | 2.40 | 0.16 |
| JGD_Trefethen/Trefethen_20000 | 0.04 s | **1.01** | 1.02 | 3.27 | 0.15 |
| Rothberg/cfd1 | 0.03 s | **0.96** | 0.98 | 0.81 | 0.24 |

**The two lists are almost the same matrices at wildly different costs.** `uni_chimera_i5` is 15.9
seconds of `MmdVendored` and 0.07 of `AmdVendored`, a factor of 220. That is the 21.3x above, seen
one row at a time.

**CHAINING WINS NOWHERE IN EITHER LIST**, and over the whole set `MmdChained` beats `MmdFlat` on 4
matrices of 234, all four under a millisecond of vendored time.

**AND THE TOP TEN DOES NOT SEPARATE THE TWO MMD STORES EITHER.** `MmdCompacted` is ahead on 3 rows
of the 10 and behind on 6, by hundredths. On the amd side it is ahead on 4 of 10 and the rows it
loses are the cheap ones, which is the aggregate result seen row by row.

## The layout question: four clique stores, three of them measured

**WHAT THEY ACTUALLY ARE, because the ratios mean nothing without it.** Four axes, and they agree on
only one of them:

| | clique space | shape | bounded | reclamation |
|---|---|---|---|---|
| `AmdFlat`, `MmdFlat` | separate arena | FLAT | no, it grows | none, ever |
| the two compacted, `AmdVendored` | unified with A and I | FLAT | yes, plus a fifth | compaction |
| `MmdChained`, `MmdVendored` | dead segments | FLAT | yes, exactly nnz | reuse by chaining |
| a chunked store | one allocation per clique | CHUNKED | no | free, on death |

So our arena buys simplicity with memory: no collector to write, no bound to respect, and a store
that on this set reaches 98.1 million entries on the amd side and 199.9 on the mmd, 4 and 53 times
`tril(A)` on their worst rows. The two vendored schemes buy a hard bound with machinery, and with
DIFFERENT machinery, which is the whole point of having siblings on both branches. `AmdVendored`
compacts a workspace sized at a fifth over the pattern; `MmdVendored` chains a clique through the
dead segments of the
cliques it absorbed and needs no spare space at all.

**FLAT IS WHY COMPACTION IS THE ONLY WAY TO RECLAIM, and that is the axis the fourth row is on.** In
a flat store a clique must be one contiguous block, so the space a dead one leaves is a hole that
nothing can use until everything live is slid down over it; the collector is a consequence of the
shape, not a separate design choice. Give each clique its own allocation and reclamation costs
nothing at all: a dead clique is freed where it lies, no bound is needed, and no collector exists.
What is paid instead is the allocator, a header per clique and rounding, on every clique ever
formed.

**That fourth row is NOT BUILT**, and this report is where the case for it would come from. The
`pC/cC` column says what it would save, and the section on the arena below reads it: over the amd
set a chunked store would hold 0.19 of what the flat arena does, and on the mmd side 0.10, with the
worst rows at 0.01. `docs/DESIGN_DECISIONS.md` (2026-08-18) sets out the two axes all four rows sit
on, memory against machinery; `experiments/ordering/README.md` has the mechanism.

**The comparison below is therefore not "is compaction expensive".** It is the whole bargain at
once: the collector's cost, minus whatever the unified layout gains by keeping a vertex's cliques
and its lists in one array, against an arena that pays nothing to reclaim because it never does.
Over the 246:

| | total | against its vendored | against our arena |
|---|---:|---:|---:|
| `AmdVendored` | 3.55 s | | |
| `AmdFlat`, arena | 4.02 s | 1.13 | |
| `AmdCompacted`, compacted | 3.96 s | **1.12** | **0.99** |
| `MmdVendored` | 75.5 s | | |
| `MmdFlat`, arena | 49.1 s | 0.65 | |
| `MmdChained`, chained | 83.3 s | 1.10 | **1.70** |
| `MmdCompacted`, compacted | 49.9 s | 0.66 | **1.01** |

**CHAINING LOSES AND IT LOSES BADLY. COMPACTION COSTS NOTHING ON EITHER BRANCH AND WINS A LITTLE
ON AMD.**

- `AmdCompacted` beats the arena on **188 of 235** matrices, median **0.957**, and the two are
  within **1.4 per cent** over the whole set.
- `MmdCompacted` beats the arena on **118 of 234**, a coin flip, median **0.998**, and is **1.4
  per cent** behind over the set.
- `MmdChained` costs **70 per cent** over the arena and wins on **4 of 234**.

**ON TIME ALONE THE COMPACTED STORE IS FREE, AND THAT IS ENOUGH, BECAUSE THE SECOND ARGUMENT IS
NOT ABOUT TIME.** A compacted store is bounded from the input; our arena is not, and cannot be,
since it grows with a quantity that depends on the ordering being computed. Given a machine you can
say whether A fits. You cannot say whether the arena will. Two of these matrices grow it past fifty
times `tril(A)`.

**SO THE COMPACTED PAIR ARE OBLIO'S MAIN DRIVERS**, `AmdCompacted` the default and `MmdCompacted`
the mmd one, and this report is what settles it. The bound is the reason and the timings are what
make it affordable: a bound bought at a real cost would be a trade, and this one is bought at none.
`MmdFlat` and `AmdFlat` stay, as the unbounded pair that the bounded ones must equal entry for
entry, which is what makes either of them checkable. `MmdChained` stays for the same reason and
because its result is the sharpest thing here: the other way of bounding the store is measurably
worse and this says by how much.

**Chaining loses precisely where it matters.** In the band holding 96 per cent of the mmd time,
`MmdFlat` reads 0.65 against `MmdVendored` and `MmdChained` reads 1.10. The grid ladders understate
this, reading 1.1 to 1.3, because grids have short cliques and few links to follow.

**The reason is in the two mechanisms and it is not a surprise.** A compaction is a rare sweep: at a
fifth of elbow room it runs about once for a whole ordering, so its cost is amortised to nearly
nothing. A chain is a link test on EVERY READ of every clique, forever, and no amount of headroom
reduces it because chaining exists precisely to need none. That asymmetry was written down as
reasoning in `experiments/ordering/README.md` before it was measured; this is the measurement.

### How often the pool is actually compacted

**`AmdCompacted` MATCHES `AMD_2` ON ALL 246 MATRICES, COMPACTION FOR COMPACTION.** Not a single row
differs, from the 122 that never compact to `uni_chimera_i1` at ten. Same workspace, same collection
schedule, on inputs neither code was tuned for. It is the sharpest faithfulness check in this
report: the permutation says the two agree on the ANSWER, and this says they agree on when the
storage ran out, which nothing else here can see.

**COMPACTION IS RARE, and that is why it is nearly free.** The median matrix compacts ONCE on both
branches, and about half never compact at all: 122 of 246 on the amd side, 100 on the mmd. Against
`MmdChained`, which pays a link test on every read of every clique for the whole run, this is the
asymmetry the timings show, seen directly rather than inferred.

| compactions | amd | mmd |
|---|---:|---:|
| none | 122 | 100 |
| one | 113 | 116 |
| more than one | 11 | 30 |
| total over the set | 150 | 268 |
| worst single matrix | 10 | 31 |

**MMD FILLS THE SAME WORKSPACE FASTER, and this is the figure that was open.** `MmdCompacted`'s
header has carried the question since the file was written: the two branches fill the pool at
different rates, amd absorbing cliques aggressively where mmd does not, so `AMD_2`'s elbow room
might not suit mmd's cliques. It suits them, but less well. Mmd compacts more often on 51 matrices,
fewer on 2, the same on 193, and its total over the set is 268 against 150.

**The gap is concentrated exactly where the arena is worst.** By `cC/tril(A)` band, the mean count
per matrix:

| `cC/tril(A)` | matrices, amd | amd | matrices, mmd | mmd |
|---|---:|---:|---:|---:|
| below 0.5, the diagonal-only rows included | 85 | 0.00 | 67 | 0.00 |
| 0.5 to 1.0 | 68 | 0.50 | 63 | 0.48 |
| 1.0 to 2.0 | 73 | 0.97 | 81 | 1.00 |
| 2.0 to 5.0 | 18 | 1.72 | 23 | 1.78 |
| 5.0 and above | 2 | 7.00 | 12 | 9.67 |

Below `cC/tril(A)` of 2 the two branches are indistinguishable and the count is under one. The whole
difference is in the top band, and it is the same band that holds 96 per cent of the mmd time: the
twelve costliest mmd matrices account for 112 of the 268 compactions.

**`uni_chimera_i5` is the extreme and it is instructive.** Mmd compacts 31 times where amd compacts
4, on a matrix where the mmd arena reaches 52.7 times `tril(A)` and the amd one 4.4. Fewer, larger
cliques and aggressive absorption keep amd inside the same headroom; mmd's many small cliques do
not. Yet `MmdCompacted` ties `MmdFlat` for fastest on that row at 0.38 against `MmdVendored`, so 31
sweeps of a pool that size still cost no more than an arena that never reclaims.

**THE ARGUMENT ABOVE IS ABOUT COST. THERE IS A SECOND ONE, ABOUT THE BOUND, and it is not in any
column here.** Compaction is the only one of the three stores that is both bounded from the input
and FLAT, so that a clique is read as one contiguous block. Chaining is bounded too and pays for it
on every read, the clique being spread across the dead segments it was threaded through. Our arena
reads a block like compaction does and cannot be bounded at all, since it grows with a quantity that
depends on the ordering being computed. That is a property rather than a measurement and it does not
move with the hardware.

**Read the sibling ratios as directions and not as coefficients.** Each sibling differs from its
driver in more than the clique store: `Amd3B` and `Mmd3C` also carry `AMD_2`'s pool walks in
positions and its mid-walk collector, and `MmdChained` carries genmmd's list conventions. The clique
store is the largest of those differences and the one each file exists to price, but it is not the
only one, and a clean measurement would be the same code with only the placement rule swapped.

**The bands with no work in them say nothing.** The purely diagonal matrices read 2.19 to 2.23 for
all three mmd columns; there is no ordering to do there and the figure is the cost of standing up a
store that is then never used.

## Where we lose, and why it does not matter

**The 14 matrices with no off-diagonal entries at all.** `Boeing/bcsstm39`, `HB/bcsstm25`,
`Cunningham/m3plates` and eleven others are purely diagonal: the ordering has nothing to do and the
permutation is the identity. On the mmd side the ratio is 2.20 and the difference is a few hundred
microseconds for the whole group. There is no result here. A routine that does nothing is measuring
its own entry and exit, and ours has a `std::vector` layer around a quotient graph it then never
uses.

**On the amd side that band reads 0.82, ours faster**, because the empty-row prepass added on
2026-08-18 numbers those vertices where they stand instead of filing and popping them from a degree
bucket. It is the same non-result in the other direction and should be read the same way.

**The nearly-trivial rest.** Sub-millisecond orderings on matrices with almost no fill, where the
ratio runs past 1.8. The pattern across the whole set is that our fixed overhead, a container layer
plus construction of a quotient graph that a trivial matrix barely touches, is a constant, and a
constant is most visible where the work is smallest.

**That is the honest shape of it and it is not being minimised.** If the question were "order ten
thousand tiny matrices", `MmdVendored` would win and we would need to answer for it. It is not the
question this benchmark asks.

## The arena: what it costs, and what a chunked store would cost instead

`cC` is what our storage costs, and the reason it is worth watching is that ours is the one scheme
here NOT bounded from the input: `MmdVendored` chains dead segments and stays inside `O(n + m)`,
`AmdVendored`
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
0.38, `uni_chimera_i1` at 0.55, `Reuters911` at 0.77, `vsp_befref` at 0.64. Spending memory to avoid
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

**But read it against `cC/nnzL` before concluding anything.** On the `uni_chimera` and `Trefethen`
rows the arena is 8 to 9 per cent of the factor, so clique storage was never the binding constraint
there; the factorization was. The rows where clique storage genuinely competes with the factor,
`Reuters911` at 0.723 and `email-Enron` at 0.512, are also rows with `pC/cC` of 0.03, which is the
case FOR building it. The rows where `pC/cC` is high, `nlpkkt80` at 0.61 and `LF10000` at 0.73, are
rows where little would be saved.

That is the whole argument for a chunked clique store, and it is now measured rather than assumed.

## Reading the whole table

Four habits, all learned by getting them wrong.

**Do not average the ratio column.** The mean over 246 rows is dominated by matrices that take
microseconds. Read the total, which is total against total, or a median within a band.

**A ratio hides which side moved.** Where a row looks anomalous, convert to time per vertex or per
edge before believing it. That failure cost a full day on the grid benchmark, where a pattern read
as ours turned out to be an artifact of the vendored routine.

**Prefer a ratio between two of ours to a ratio against a reference.** The first varies one thing;
the second varies everything that differs between two independently written codes. The p10 to p90
spreads in the summary tables are three times wider on the second kind.

**Do not compare a number in the mmd table against one in the amd table except deliberately.** The
denominators are different routines doing different amounts of work, and the amd one includes two
phases we do not run. Cross-branch statements in this report are confined to the section near the
top that makes them on purpose.

## What this report does not establish

**It is two orderings against two references, and nothing else.** Nested dissection is absent from
this tree, and on the largest matrices here it is what a production solver would actually use.

**It says almost nothing about ordering QUALITY.** Each of ours reproduces its reference's
permutation exactly, so within a branch fill is identical by construction. The one quality statement
this set supports is the cross-branch one near the top: amd fills less on 136 matrices, mmd on 92,
18 tie, and the median matrix is a tie to within 0.3 per cent. `PERFORMANCE.md` looks at what that
does to a whole solve.

**The `AmdVendored` clock is not the amd branch's cost.** It is the whole vendored call, about 13
per cent of which is work we do not do. Every amd ratio here is therefore optimistic by roughly that
much, and the figure to quote for the C++ and container cost is not in this report;
`benchmarks/ordering`'s `make phases2d` has the split.

**A ratio against a vendored routine is a figure to watch rather than one to reason from.** Our
drivers and the references differ in more than the thing being measured, and those differences have
not been enumerated. The sibling ratios, which vary one thing between two builds of our own code,
are what this report's conclusions rest on.

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

**A row ending in `order differs`, `fill differs` or `pC differs` means a sibling parted from its
driver on that matrix**, and for `order differs` the shared `nnz(L)` column has stopped being both
routines'. Both tables are clean, on all three siblings. That check is what found the
amd divergences; before it existed, the equality was asserted in this file's opening and verified
only on generated shapes.