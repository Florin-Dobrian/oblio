# Oblio: Ordering

This report asks what the ORDERING STEP ALONE costs and what space it takes, on 246 matrices from
the SuiteSparse Matrix Collection that Oblio's authors did not choose and did not generate.
`ACCURACY.md` beside it asks how good the answers are and `PERFORMANCE.md` what a whole solve costs;
this one isolates the first phase, where no values are read and no definiteness is needed, so the
input set is everything on disk rather than the subsets those two require.

**MEASURED 2026-08-23.** Every figure below comes from two runs of one program on that day,
`make mmdorder` and `make amdorder`, over the same 246 matrices. This report is a snapshot: it says
what is true of the code as it stands, not how it got there.

**BOTH BRANCHES ARE HERE.** `MmdFlat` against `MmdCorrected` and `AmdFlat` against `AmdVendored`,
with each branch's alternative-store siblings alongside.

**THE MMD REFERENCE IS `MmdCorrected` AND NOT `MmdVendored`, SINCE 2026-08-23.** genmmd files a
vertex under its degree in `mmdint` and under its degree PLUS ONE in `mmdupd`, so two scales are
live in one bucket array and a vertex the refresh has touched is penalised by one against a vertex
no pivot has reached: the minimum it selects is not always the minimum. `private/MmdCorrected.cpp`
is genmmd with that repaired and is what our three mmd drivers reproduce exactly.
`private/MmdVendored.cpp` is frozen beside it and appears in one column only, its FILL, which the
section below reads. Nothing is timed against it and no driver is compared to it.

**A DRIVER IS ITS BRANCH PLUS THE CLIQUE STORE IT RUNS ON.** `MmdCorrected` and `AmdVendored` are
the two reference codes; `MmdFlat` and `AmdFlat` are ours on our own arena; `MmdCompacted` and
`AmdCompacted` are ours on `AMD_2`'s compacted workspace; `MmdChained` is ours on genmmd's chained
segments. Per-matrix tables print three-letter tags, `Cor` or `Vnd` for the reference and `Flt Chn
Com` for the stores, with a legend line above them; this report uses the full names.

**THREE LAYOUT SIBLINGS RIDE ALONG.** `MmdChained` and `MmdCompacted` on the mmd side,
`AmdCompacted` on the amd side. Each is its branch's driver on a different clique store, so the
sibling columns ARE the layout comparison. Every driver here compiles its quotient graph into its
own translation unit, so a comparison between two of ours is between two things built the same way;
that was not true before 2026-08-19 and figures from before then are not comparable to these.

The comparison is unusually clean on both. **Each of ours returns its reference's exact permutation
on every matrix in the set**, checked per matrix inside the run rather than assumed, so there is no
fill to trade against time: one `nnz(L)` column serves both routines in each table, and time is the
only thing that differs.

That check is new and it earned itself. It found FOUR divergences on the amd side, one of them a
correctness bug in supervariable detection that had been silently costing fill; see `notes/NEXT.md`
(2026-08-18). The mmd side needed nothing. **In these runs not one marker fired**, on either branch
or any of the three siblings.

**The short answers, one line each.**

| | reference | ours | ratio | ours slower on |
|---|---:|---:|---:|---:|
| mmd | 75.2 s | 51.0 s | **0.67** | 131 of 228 |
| amd | 3.59 s | 4.04 s | **1.13** | 135 of 234 |

The default is neither of those two. `AmdCompacted` orders the set in **4.04 s**, and it is what a
caller who sets nothing now gets; see `include/oblio/OrderEngine.h`.

Those two rows are not comparable to each other and the reasons are in the next section. Within
each, "slower on most matrices and faster overall" is the mmd shape and reconciling it is most of
what this report is about; the amd branch is slower on both counts and honestly so.

## The one number that reframes the rest

**MMD costs 21 times what AMD costs on this set, for 3.5 per cent more fill.** 75.2 seconds against
3.59, and 18.74 billion factor entries against 18.11 billion. Read the direction carefully: mmd is
the one that costs more, on BOTH axes: it is not a trade of time against quality.

That is a property of the two ALGORITHMS rather than of either implementation, and it dwarfs every
implementation ratio in this report. Minimum degree recomputes exact degrees; approximate minimum
degree bounds them and never looks at a clique twice. On a set like this the bound costs almost
nothing in quality and saves an order of magnitude in time.

It also means the two halves of this report are asking different questions. The mmd table is about a
75-second budget where a 33 per cent saving is 24 seconds; the amd table is about a 3.6-second
budget where a 13 per cent loss is 0.45 seconds. Read them accordingly.

### Which wins, counted matrix by matrix

The 21x is a total and the totals are the giant matrices. Counting instead how often each branch
wins, one vote per matrix, all 246:

| | mmd wins | amd wins | equal |
|---|---:|---:|---:|
| **time, references**, `MmdCorrected` against `AmdVendored` | 89 | **157** | 0 |
| **time, ours**, `MmdFlat` against `AmdFlat` | 56 | **189** | 1 |
| **time, the defaults**, `MmdCompacted` against `AmdCompacted` | 50 | **193** | 3 |
| **fill** | 93 | **131** | 22 |

**Amd wins on every axis and by a clear margin, and the margin is wider for our pair than for the
reference one**, 189 to 56 against 157 to 89. That is the two implementation results seen from
another angle: `MmdFlat` beats `MmdCorrected` mostly on the giant matrices, which are few, while
`AmdFlat` loses to `AmdVendored` fairly evenly, so on a one-vote-per-matrix count our amd pulls
further ahead.

**The fill row is one column, not two.** Each of ours returns its reference's exact permutation, so
every mmd driver fills identically and so does every amd one; the 93 against 131 is a statement
about the two ALGORITHMS. And it is far closer than the 3.5 per cent total suggests: the median
matrix has `mmd / amd` fill of **1.000**, an exact tie, with p10 0.969 and p90 1.059 and a long tail
either way, 0.688 at best for mmd and 2.293 at worst. Twenty-two matrices tie exactly. The total and
the median disagree because three PARSEC-class matrices carry billions of entries each; see the
section below on what the correction cost.

## What the corrected rule cost, and it is three matrices

The mmd reference changed on 2026-08-23 and this set is the first real measurement of what that is
worth. genmmd files on two scales and the corrected rule files at the true degree everywhere; the
`Vnd nnz(L)` and `Cor/Vnd` columns on the mmd table are the frozen routine's fill and the ratio to
ours. It is the only thing `MmdVendored` is used for here.

**THE TOTAL SAYS THE CORRECTION COSTS 1.8 PER CENT OF FILL.** 18.744 billion entries against
18.408, a net **+336 million**.

**THE COUNT SAYS THE OPPOSITE, AND THE COUNT IS THE HONEST ONE.**

| | |
|---|---|
| corrected fills MORE on | 83 matrices |
| corrected fills LESS on | **92 matrices** |
| identical on | 71 matrices |
| `Cor/Vnd` median | **1.0000** |
| mean, p10, p90 | 0.9992, 0.9703, 1.0319 |
| worst either way | 0.8236 and 1.3053 |

**The two disagree because three matrices carry the total.**

| | corrected | frozen | change |
|---|---:|---:|---:|
| `PARSEC/Si87H76` | 6.059 B | 5.680 B | **+380 M** |
| `Schenk/nlpkkt80` | 3.486 B | 3.679 B | **-193 M** |
| `PARSEC/Ga41As41H72` | 7.224 B | 7.078 B | **+145 M** |
| the other 243 together | | | +4 M |

Those three are **98.9 per cent** of the net change and the top ten are 99.8. `Si87H76` alone moves
more than the whole net.

**SO THE MEASUREMENT DOES NOT SUPPORT "THE MIXED SCALE IS A BETTER HEURISTIC".** On a typical matrix
the two rules are indistinguishable, the median is an exact tie, and slightly more matrices improve
than worsen. What the total reflects is a handful of PARSEC-class matrices where minimum degree is
already producing billions of entries and any perturbation moves hundreds of millions. Those are the
rows this report elsewhere names as the case for nested dissection, which is absent from this tree.

**And the bias it removes was systematic rather than random.** A refreshed vertex always sat one
bucket high against one no pivot had reached, so genmmd consistently preferred the stale vertex
among near-equals. That is a real difference in behaviour; this set says it is worth nothing on a
typical matrix and is a coin flip in sign. Whether the tail behaviour is a heuristic worth having is
not settled here and would need the large matrices studied individually.

## The environment

| | |
|---|---|
| Machine | MacBook Pro (Mac16,1), Apple M4, 32 GB, macOS 26.6.2 |
| Compiler | Apple Clang (Xcode 26.6), C++17, `-O3 -DNDEBUG` |
| References | `MmdCorrected` is genmmd, the multiple-minimum-degree of George and Liu, with its |
| | degree scale repaired; `MmdVendored` is the same code frozen, used for fill only |
| | `AmdVendored` is `amd_order` from SuiteSparse AMD 3.3.4, Amestoy, Davis and Duff |
| Naming | `MmdCorrected` and `AmdVendored` are the references throughout; five are ours. |
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
| `MmdFlat / MmdCorrected` | 1.043 | 1.061 | 0.786 | 1.361 | **0.678** | 228 |
| `MmdChained / MmdCorrected` | 1.237 | 1.271 | 1.067 | 1.501 | 1.113 | 228 |
| `MmdCompacted / MmdCorrected` | 1.039 | 1.058 | 0.784 | 1.332 | **0.676** | 228 |
| `MmdCompacted / MmdFlat` | **0.996** | 0.999 | 0.966 | 1.040 | **0.997** | 234 |
| `MmdChained / MmdFlat` | **1.172** | 1.229 | 1.046 | 1.411 | 1.641 | 234 |

**amd**

| | median | mean | p10 | p90 | total | n |
|---|---:|---:|---:|---:|---:|---:|
| `AmdFlat / AmdVendored` | 1.022 | 1.036 | 0.759 | 1.308 | 1.125 | 234 |
| `AmdCompacted / AmdVendored` | 1.010 | 1.010 | 0.733 | 1.273 | 1.126 | 234 |
| `AmdCompacted / AmdFlat` | **0.977** | 0.975 | 0.909 | 1.034 | **1.001** | 235 |

**Four things to read off them.**

**The mmd median and total disagree completely**, 1.043 against 0.678. We lose slightly on most
matrices and win by a third on the expensive ones. The amd branch does not do this, 1.022 against
1.125, so it loses a little everywhere and slightly more where it costs.

**`MmdCompacted / MmdFlat` IS A WASH.** 0.996 at the median, 0.997 in total, p10 to p90 spanning
0.966 to 1.040. On the mmd branch the clique store makes no measurable difference to time in either
direction, and the spread is wider than the effect, so a single matrix says nothing about it. It
read 0.998 and 1.014 in the 2026-08-21 run, on the other side of 1.0 and by about as much, which is
the size of the run-to-run noise on this ratio.

**`AmdCompacted / AmdFlat` IS A REAL WIN AND A SMALL ONE.** 0.977 at the median, winning on 158 of
235 matrices, with p10 to p90 of 0.909 to 1.034. Its total is 1.001, so over the whole set the two
are level while the per-matrix win is two to three per cent: the count and the median carry this
result, not the total.

**THE AMD ORDERING DID NOT CHANGE BETWEEN THE TWO RUNS AND ITS RATIOS STILL MOVED.**
`AmdCompacted / AmdVendored` read 0.996 at the median on 2026-08-21 and reads 1.010 here;
`AmdFlat / AmdVendored` moved from 1.037 to 1.022. Fill reproduces to the entry, 18.11 billion both
times, which is the control that says the ANSWER did not move.

Be careful with the word unchanged, though: the shared `Buckets` encoding was rewritten on
2026-08-23, unbiasing the predecessor, moving both sentinels and dropping the head array from n + 1
to n, and both amd drivers use that class. So the amd movement is some mixture of that and the
measurement's own spread, and this run cannot separate them. What it does establish is the size:
one to two hundredths, on a branch whose output is bit-identical. Read a hundredth in these tables
as noise either way.

**The p10 to p90 spread is wide against the references and narrow between our own drivers.**
0.786 to 1.361 for `MmdFlat / MmdCorrected` against 0.966 to 1.040 for `MmdCompacted / MmdFlat`.
Comparing two of our drivers is a much sharper instrument than comparing ours to a reference, which
mixes in every difference between two independently written codes, most of which nobody has
enumerated. The translation-unit work of 2026-08-19 reached the same conclusion from the other
direction: aligning our drivers with each other made our own ratios comparable and moved the ratios
against the references not at all. The narrow rows are the ones to reason from.

**`AmdCompacted / AmdVendored` reaches 1.010 at the median.** The compacted store takes the amd
branch to within one per cent of `AMD_2` on a typical matrix, from 1.022. That is parity to the
precision these runs support, and not a win.

## What the columns mean

```
n   m   tril(A)   A+I   cC   pC   nnz(L)   cC/tril   cC/nnzL   pC/cC   Cor ms   Flt ms   Flt/Cor
```

**THE DRIVER COLUMNS ARE TAGGED AND THE PROGRAM PRINTS A LEGEND LINE ABOVE THE TABLE.** `Cor` the
mmd reference and `Vnd` the amd one, `Flt` our arena, `Chn` the chained store, `Com` the compacted
one, and on the amd side `Raw` the hooked pre-postorder copy. A run is mmd or amd and never both, so
the branch is said once in that line rather than in every heading. The tags are the 2026-08-21
naming entry's and they live here, in the same legend as `cC` and `pC`; prose and the summary tables
above use the full names.

**AND THE MMD TABLE ENDS WITH TWO COLUMNS THAT ARE NOT A DRIVER.** `Vnd nnz(L)` is the frozen
`MmdVendored`'s fill and `Cor/Vnd` the ratio to ours, above 1 where the corrected rule fills more.
That routine is called once per matrix outside the clock and nothing else is taken from it.

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
figure is identical for `MmdFlat`, `MmdChained` and `MmdCompacted`, and likewise for `AmdFlat` and
`AmdCompacted`. No reference routine can report one, and `cC` exists only on the flat store: the
compacted pool is sized at construction and a chained clique lives in its own pivot's dead segment,
so on neither is clique storage a quantity that varies.

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
neither does the mmd reference; the mmd table prints a dash where the amd table prints `Vnd cmp`.

**Read `pC/cC` conditioned on `cC/nnzL`, or it will mislead.** Where the arena is one per cent of
the factor, clique storage is noise beside what the factorization will need and no clique scheme is
worth building for it whatever `pC/cC` says. Where the arena is half the factor, it decides.

### Each table carries its siblings too, and only their times

The mmd table ends with `Chn ms`, `Chn/Cor`, `Com ms`, `Com/Cor`, `Com cmp` and the two fill
columns; the amd table with `Com ms`, `Com/Vnd`, `Vnd cmp` and `Com cmp`. A sibling is its branch's
driver on a DIFFERENT CLIQUE STORE rather than on our arena, so **the sibling columns are the layout
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
checks, `test_order` or the sanitizers could see; see `notes/NEXT.md` (2026-08-19).

**Every marker was clear across the 246 on both branches**, on all three siblings.

On the amd table there is one caveat with no counterpart on the mmd side: **the reference clock is
on the whole vendored call**, which is what a caller pays and includes `AMD_aat` and
`AMD_postorder`, neither of which we do. A phase split puts those plus argument validation at about
13 per cent of that column, so against the comparable region alone the amd ratio is nearer 1.28 than
the 1.13 the table shows. Everything except the clock comes from the hooked pre-postorder copy,
since `AmdFlat` returns that order. MMD does no postorder at all, so its table needs no such copy
and carries no such caveat.

## Where the time actually goes

The set spans five orders of magnitude in ordering time, so a count of wins and losses says almost
nothing. Grouped by how large the arena grew, `cC/tril(A)`. **The ratios below are AGGREGATED**, the
band's total ours divided by its total vendored, so the big matrices dominate them.

**mmd, against `MmdCorrected`**

| `cC/tril(A)` | matrices | share of total time | `MmdFlat` | `MmdChained` | `MmdCompacted` |
|---|---:|---:|---:|---:|---:|
| 0, no off-diagonal at all | 14 | 0.0% | 2.14 | 2.24 | 2.25 |
| 0 to 0.5 | 54 | 0.3% | 0.76 | 1.01 | 0.76 |
| 0.5 to 1.0 | 60 | 0.4% | 0.89 | 1.12 | 0.89 |
| 1.0 to 2.0 | 80 | 1.7% | 0.79 | 1.05 | 0.78 |
| 2.0 to 5.0 | 26 | 1.6% | 0.72 | 1.12 | 0.72 |
| 5.0 and above | 12 | **96.0%** | **0.67** | 1.11 | **0.67** |

**The twelve hardest matrices are 96 per cent of all the work, and we are 33 per cent faster on
them.** The 149 matrices where we are slower cost **19.5 milliseconds in total, all of them
together**; the 97 where we are faster save **24.2 seconds**. The largest single loss anywhere is
1.11 ms, on `UTEP/Dubcova2`.

**AND `MmdChained` LOSES EXACTLY WHERE `MmdFlat` WINS**, 1.11 against 0.67 in the band holding 96
per cent of the time. That is chaining's price and it is the one store result that is not in doubt.

**`MmdCompacted` DOES NOT SEPARATE FROM `MmdFlat` IN ANY BAND.** It is level in four bands to the
hundredth, including both that carry the work, ahead by a hundredth in one and behind by eleven
thousandths in the diagonal-only band where there is no work at all. That is what no effect looks
like.

**amd, against `AmdVendored`**

| `cC/tril(A)` | matrices | share of total time | `AmdFlat` | `AmdCompacted` |
|---|---:|---:|---:|---:|
| 0, no off-diagonal at all | 13 | 0.0% | **0.77** | 0.80 |
| 0 to 0.5 | 71 | 2.2% | 0.85 | **0.81** |
| 0.5 to 1.0 | 69 | 5.4% | 1.07 | **1.03** |
| 1.0 to 2.0 | 72 | 6.2% | 1.11 | **1.10** |
| 2.0 to 5.0 | 19 | **73.8%** | 1.13 | 1.13 |
| 5.0 and above | 2 | 12.4% | 1.19 | 1.20 |

**The amd store leads or ties in five bands of six**, and where it trails it does so in bands with
little or no work. That is the shape the mmd side does not have, but read it against the run-to-run
spread above: these are hundredths on an unchanged branch.

**The amd branch has no such reversal, and the shape is the opposite one.** We are faster on the
small end, where the dense-row rule and the empty-row prepass now do what `AmdVendored` does and our
container layer has little to be slow at, and slower on the large end where 86 per cent of the time
is. Read against the comparable region rather than the whole vendored call, the aggregated figures
move up by roughly a seventh.

**The ten most expensive matrices, each branch.** They are 95 per cent of the mmd total and 85 per
cent of the amd one.

| mmd | `MmdCorrected` | `MmdFlat` | `MmdChained` | `MmdCompacted` | `cC/tril` | `pC/cC` |
|---|---:|---:|---:|---:|---:|---:|
| FlowIPM22/uni_chimera_i5 | 15.17 s | **0.41** | 1.14 | **0.41** | 53.79 | 0.01 |
| PARSEC/Ga41As41H72 | 13.60 s | 0.80 | 1.07 | **0.79** | 5.63 | 0.02 |
| PARSEC/Si87H76 | 11.81 s | 0.79 | 1.07 | **0.78** | 8.26 | 0.02 |
| FlowIPM22/uni_chimera_i1 | 10.59 s | **0.57** | 1.09 | 0.58 | 53.18 | 0.01 |
| JGD_Trefethen/Trefethen_20000 | 10.02 s | **0.75** | 1.15 | **0.75** | 7.27 | 0.07 |
| JGD_Trefethen/Trefethen_20000b | 6.33 s | **0.84** | 1.19 | **0.84** | 24.14 | 0.03 |
| Pajek/Reuters911 | 1.35 s | **0.78** | 1.15 | **0.78** | 10.07 | 0.03 |
| DIMACS10/vsp_befref_fxm_2_4_air02 | 0.93 s | **0.68** | 1.16 | **0.68** | 9.35 | 0.08 |
| Andrews/Andrews | 0.77 s | **0.62** | 1.06 | **0.62** | 6.13 | 0.06 |
| Schenk/nlpkkt80 | 0.72 s | 0.77 | 1.04 | **0.75** | 1.42 | 0.61 |

| amd | `AmdVendored` | `AmdFlat` | `AmdCompacted` | `cC/tril` | `pC/cC` |
|---|---:|---:|---:|---:|---:|
| Schenk/nlpkkt80 | 0.96 s | **1.09** | 1.18 | 2.18 | 0.38 |
| PARSEC/Ga41As41H72 | 0.78 s | 1.14 | **1.07** | 2.05 | 0.05 |
| PARSEC/Si87H76 | 0.58 s | 1.21 | **1.16** | 2.91 | 0.06 |
| FlowIPM22/uni_chimera_i1 | 0.38 s | 1.21 | 1.21 | 15.16 | 0.03 |
| FlowIPM22/uni_chimera_i5 | 0.08 s | 1.10 | 1.10 | 4.37 | 0.09 |
| Andrews/Andrews | 0.07 s | **1.11** | 1.12 | 3.85 | 0.09 |
| JGD_Trefethen/Trefethen_20000b | 0.07 s | **1.09** | 1.12 | 6.08 | 0.10 |
| FlowIPM22/uni_chimera_i2 | 0.06 s | **0.97** | 1.01 | 2.40 | 0.16 |
| JGD_Trefethen/Trefethen_20000 | 0.04 s | **1.00** | 1.02 | 3.27 | 0.15 |
| Rothberg/cfd1 | 0.03 s | 0.96 | 0.96 | 0.81 | 0.24 |

**The two lists are almost the same matrices at wildly different costs.** `uni_chimera_i5` is 15.2
seconds of `MmdCorrected` and 0.08 of `AmdVendored`, a factor of 190. That is the 21x above, seen
one row at a time.

**CHAINING WINS NOWHERE IN EITHER LIST**, and over the whole set `MmdChained` beats `MmdFlat` on 8
matrices of 234, all of them cheap rows.

**AND THE TOP TEN DOES NOT SEPARATE THE TWO MMD STORES EITHER.** `MmdCompacted` ties `MmdFlat` on 6
rows of the 10 and is ahead on 3 by a hundredth. On the amd side it is ahead on 2 of 10 and level on
3, the rows it loses being the cheap ones.

## The layout question: four clique stores, three of them measured

**WHAT THEY ACTUALLY ARE, because the ratios mean nothing without it.** Four axes, and they agree on
only one of them:

| | clique space | shape | bounded | reclamation |
|---|---|---|---|---|
| `AmdFlat`, `MmdFlat` | separate arena | FLAT | no, it grows | none, ever |
| the two compacted, `AmdVendored` | unified with A and I | FLAT | yes, plus a fifth | compaction |
| `MmdChained`, `MmdCorrected` | dead segments | FLAT | yes, exactly nnz | reuse by chaining |
| a chunked store | one allocation per clique | CHUNKED | no | free, on death |

So our arena buys simplicity with memory: no collector to write, no bound to respect, and a store
that on this set reaches 98.1 million entries on the amd side and 203.5 on the mmd, 4 and 54 times
`tril(A)` on their worst rows. The two reference schemes buy a hard bound with machinery, and with
DIFFERENT machinery, which is the whole point of having siblings on both branches. `AmdVendored`
compacts a workspace sized at a fifth over the pattern; the mmd reference chains a clique through
the dead segments of the cliques it absorbed and needs no spare space at all.

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
worst rows at 0.01. `notes/DESIGN_DECISIONS.md` (2026-08-18) sets out the two axes all four rows sit
on, memory against machinery; `experiments/ordering/README.md` has the mechanism.

**The comparison below is therefore not "is compaction expensive".** It is the whole bargain at
once: the collector's cost, minus whatever the unified layout gains by keeping a vertex's cliques
and its lists in one array, against an arena that pays nothing to reclaim because it never does.
Over the 246:

| | total | against its reference | against our arena |
|---|---:|---:|---:|
| `AmdVendored` | 3.59 s | | |
| `AmdFlat`, arena | 4.04 s | 1.13 | |
| `AmdCompacted`, compacted | 4.04 s | **1.13** | **1.00** |
| `MmdCorrected` | 75.2 s | | |
| `MmdFlat`, arena | 51.0 s | 0.68 | |
| `MmdChained`, chained | 83.7 s | 1.11 | **1.64** |
| `MmdCompacted`, compacted | 50.8 s | 0.68 | **1.00** |


**CHAINING LOSES AND IT LOSES BADLY. COMPACTION COSTS NOTHING ON EITHER BRANCH AND WINS A LITTLE
ON AMD.**

- `AmdCompacted` beats the arena on **158 of 235** matrices, median **0.977**, and the two are
  level over the whole set, 1.001.
- `MmdCompacted` beats the arena on **135 of 234**, median **0.996**, and is level over the set at
  0.997. The 2026-08-21 run had this at 0.998 and 1.014: a wash in both runs.
- `MmdChained` costs **64 per cent** over the arena and wins on **8 of 234**.

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
`MmdFlat` reads 0.67 against `MmdCorrected` and `MmdChained` reads 1.11. The grid ladders understate
this, reading 1.1 to 1.3, because grids have short cliques and few links to follow.

**The reason is in the two mechanisms and it is not a surprise.** A compaction is a rare sweep: at a
fifth of elbow room it runs about once for a whole ordering, so its cost is amortised to nearly
nothing. A chain is a link test on EVERY READ of every clique, forever, and no amount of headroom
reduces it because chaining exists precisely to need none. That asymmetry was written down as
reasoning in `experiments/ordering/README.md` before it was measured; this is the measurement.

### How often the pool is actually compacted

**`AmdCompacted` MATCHES `AMD_2` ON ALL 246 MATRICES, COMPACTION FOR COMPACTION.** Not a single row
differs, from the 122 that never compact to `uni_chimera_i1` at ten. It held on 2026-08-21 and holds
again here. Same workspace, same collection
schedule, on inputs neither code was tuned for. It is the sharpest faithfulness check in this
report: the permutation says the two agree on the ANSWER, and this says they agree on when the
storage ran out, which nothing else here can see.

**COMPACTION IS RARE, and that is why it is nearly free.** The median matrix compacts ONCE on both
branches, and about half never compact at all: 122 of 246 on the amd side, 98 on the mmd. Against
`MmdChained`, which pays a link test on every read of every clique for the whole run, this is the
asymmetry the timings show, seen directly rather than inferred.

| compactions | amd | mmd |
|---|---:|---:|
| none | 122 | 98 |
| one | 113 | 119 |
| more than one | 11 | 29 |
| total over the set | 150 | 273 |
| worst single matrix | 10 | 33 |

**MMD FILLS THE SAME WORKSPACE FASTER, and this is the figure that was open.** `MmdCompacted`'s
header has carried the question since the file was written: the two branches fill the pool at
different rates, amd absorbing cliques aggressively where mmd does not, so `AMD_2`'s elbow room
might not suit mmd's cliques. It suits them, but less well. Mmd compacts more often on 52 matrices,
fewer on 2, the same on 192, and its total over the set is 273 against 150.

**The gap is concentrated exactly where the arena is worst.** By `cC/tril(A)` band, the mean count
per matrix:

| `cC/tril(A)` | matrices, amd | amd | matrices, mmd | mmd |
|---|---:|---:|---:|---:|
| below 0.5, the diagonal-only rows included | 84 | 0.00 | 68 | 0.00 |
| 0.5 to 1.0 | 69 | 0.49 | 60 | 0.50 |
| 1.0 to 2.0 | 72 | 0.97 | 80 | 1.00 |
| 2.0 to 5.0 | 19 | 1.68 | 26 | 1.65 |
| 5.0 and above | 2 | 7.00 | 12 | 10.00 |

Below `cC/tril(A)` of 2 the two branches are indistinguishable and the count is under one. The whole
difference is in the top band, and it is the same band that holds 96 per cent of the mmd time: the
twelve costliest mmd matrices account for 120 of the 273 compactions.

**`uni_chimera_i5` is the extreme and it is instructive.** Mmd compacts 33 times where amd compacts
4, on a matrix where the mmd arena reaches 53.8 times `tril(A)` and the amd one 4.4. Fewer, larger
cliques and aggressive absorption keep amd inside the same headroom; mmd's many small cliques do
not. Yet `MmdCompacted` ties `MmdFlat` for fastest on that row at 0.41 against `MmdCorrected`, so 33
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

**The bands with no work in them say nothing.** The purely diagonal matrices read 2.14 to 2.25 for
all three mmd columns; there is no ordering to do there and the figure is the cost of standing up a
store that is then never used.

## Where we lose, and why it does not matter

**The 14 matrices with no off-diagonal entries at all.** `Boeing/bcsstm39`, `HB/bcsstm25`,
`Cunningham/m3plates` and eleven others are purely diagonal: the ordering has nothing to do and the
permutation is the identity. On the mmd side the ratio is 2.14 and the difference is a few hundred
microseconds for the whole group. There is no result here. A routine that does nothing is measuring
its own entry and exit, and ours has a `std::vector` layer around a quotient graph it then never
uses.

**On the amd side that band reads 0.77, ours faster**, because the empty-row prepass added on
2026-08-18 numbers those vertices where they stand instead of filing and popping them from a degree
bucket. It is the same non-result in the other direction and should be read the same way.

**The nearly-trivial rest.** Sub-millisecond orderings on matrices with almost no fill, where the
ratio runs past 1.8. The pattern across the whole set is that our fixed overhead, a container layer
plus construction of a quotient graph that a trivial matrix barely touches, is a constant, and a
constant is most visible where the work is smallest.

**That is the honest shape of it and it is not being minimised.** If the question were "order ten
thousand tiny matrices", `MmdCorrected` would win and we would need to answer for it. It is not the
question this benchmark asks.

## The arena: what it costs, and what a chunked store would cost instead

`cC` is what our storage costs, and the reason it is worth watching is that ours is the one scheme
here NOT bounded from the input: the mmd reference chains dead segments and stays inside `O(n + m)`,
`AmdVendored` compacts a fixed workspace, and we allocate cliques as they form and never reclaim the
space. See
`notes/DESIGN_DECISIONS.md` (2026-08-18) for the two axes that frames, and why the trade is
deliberate.

**The compression is often enormous.** `PARSEC/Ga41As41H72` factors to 7.22 billion entries and the
mmd arena holds 52.8 million, `cC/nnz(L) = 0.007`. `Schenk/nlpkkt80`, a million rows, factors to
3.49 billion and the arena is 21.1 million, 0.006. A clique holds the update part of a
supervariable's column and nothing else, and mass elimination gives one clique per supervariable
rather than one per column, so what we store tracks the input far more closely than the factor.

**And sometimes it is not.** `Pajek/Reuters911` has `cC/nnz(L) = 0.726` on the mmd side: the arena
is nearly three-quarters of the factor and 10.07 times `tril(A)`. `SNAP/email-Enron` is 0.501.
`FlowIPM22/uni_chimera_i5` reaches `cC/tril(A) = 53.79`, an arena fifty times the input it came
from. These are the cases the guarantee argument is about: social graphs and interior-point systems
where supervariables barely form, the compression fails, and our space is neither small nor
predictable.

**The striking thing is that those are also where the mmd branch is fastest.** `uni_chimera_i5` at
0.41, `uni_chimera_i1` at 0.57, `Reuters911` at 0.78, `vsp_befref` at 0.68. Spending memory to avoid
chaining is exactly the trade being made, and this set says it pays where it costs the most.

**AMD'S ARENA IS LESS THAN HALF THE SIZE FOR THE SAME MATRICES**, 98.1 million entries against mmd's
203.5 across the set, and the worst cases are far milder: `uni_chimera_i5` is 4.37 times `tril(A)`
where mmd is 53.79, `Reuters911` is 1.16 against 10.07. Fewer, larger cliques and aggressive
absorption are why.

### What `pC` says, and it is the new question

`pC/cC` is the fraction a CHUNKED clique store would hold at its worst moment. Over the whole set it
is **0.10 on the mmd side and 0.19 on the amd side**; per matrix the medians are 0.25 and 0.20.

**The extremes are where it matters and they are all on the mmd side.** `uni_chimera_i1` and
`uni_chimera_i5` read 0.01, `Reuters911`, `email-Enron` and `Trefethen_20000b` 0.03. Those are
exactly the matrices where `cC/tril(A)` is 7 to 54. So the arena's worst behaviour is almost
entirely dead space: a store that freed a clique when it died would hold **one to three per cent**
of what the flat one holds, on the rows that hurt.

**But read it against `cC/nnzL` before concluding anything.** On the `uni_chimera` and `Trefethen`
rows the arena is 3 to 9 per cent of the factor, so clique storage was never the binding constraint
there; the factorization was. The rows where clique storage genuinely competes with the factor,
`Reuters911` at 0.726 and `email-Enron` at 0.501, are also rows with `pC/cC` of 0.03, which is the
case FOR building it. The rows where `pC/cC` is high, `Lourakis/bundle1` at 0.88 and `nlpkkt80` at
0.61, are rows where little would be saved.

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
spreads in the summary tables are seven times wider on the second kind.

**And do not read a hundredth as a result.** The amd branch produced a bit-identical ordering in the
2026-08-21 and 2026-08-23 runs and its ratios still moved by one to two hundredths. Some of that is
the shared bucket encoding, which was rewritten between them; the rest is the measurement.
Differences of that size in these tables are not results.

**Do not compare a number in the mmd table against one in the amd table except deliberately.** The
denominators are different routines doing different amounts of work, and the amd one includes two
phases we do not run. Cross-branch statements in this report are confined to the section near the
top that makes them on purpose.

## What this report does not establish

**It is two orderings against two references, and nothing else.** Nested dissection is absent from
this tree, and on the largest matrices here it is what a production solver would actually use.

**It says almost nothing about ordering QUALITY.** Each of ours reproduces its reference's
permutation exactly, so within a branch fill is identical by construction. Two quality statements
this set does support. The cross-branch one near the top: amd fills less on 131 matrices, mmd on 93,
22 tie, and the median matrix is an exact tie. And the within-branch one, the frozen genmmd against
the corrected rule: 92 to 83 in the corrected rule's favour by count, a median of exactly 1.0000,
and a total that goes the other way on the strength of three matrices. `PERFORMANCE.md` looks at
what fill does to a whole solve.

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

The private sources are needed for the reference columns, for the mmd table's frozen fill column,
and, on the amd side, for the hooked pre-postorder copy the fill column comes from. Without them
either table still runs and prints ours alone.

**A row ending in `order differs`, `fill differs` or `pC differs` means a sibling parted from its
driver on that matrix**, and for `order differs` the shared `nnz(L)` column has stopped being both
routines'. Both tables are clean, on all three siblings. That check is what found the
amd divergences; before it existed, the equality was asserted in this file's opening and verified
only on generated shapes.