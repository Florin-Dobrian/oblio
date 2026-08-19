# Oblio: Ordering

This report asks what the ORDERING STEP ALONE costs and what space it takes, on 246 matrices from
the SuiteSparse Matrix Collection that Oblio's authors did not choose and did not generate.
`ACCURACY.md` beside it asks how good the answers are and `PERFORMANCE.md` what a whole solve
costs; this one isolates the first phase, where no values are read and no definiteness is needed,
so the input set is everything on disk rather than the subsets those two require.

The comparison is unusually clean. Oblio's `MMD3` returns **the same permutation as genmmd** on
every matrix here, which `make mmdmatrices` in `experiments/ordering` checks independently. So
there is no fill to trade against time: one `nnz(L)` column serves both routines, and time is the
only thing that differs.

**The short answer, in one line: 74.7 seconds of genmmd against 52.7 of ours, a ratio of 0.71, and
Oblio is slower on 163 of the 246 matrices.** Those two facts are not in tension, and reconciling
them is most of what this report is about.

## The environment

| | |
|---|---|
| Machine | Apple M4, 32 GB, macOS 26.6.1 |
| Compiler | Apple Clang (Xcode 26.6), C++17, `-O3 -DNDEBUG` |
| Reference | genmmd, the Fortran multiple-minimum-degree of George and Liu, translated |
| Timing | best of N after a warm-up, both routines through the same helper with the same N |

Both routines are timed by the same path on purpose. A column timed one way against a column timed
another differed 2.4 percent on 2026-08-10 and the difference was briefly read as a result.

## What the columns mean

```
n   m   tril(A)   A+I   C   nnz(L)   C/tril   C/nnzL   MMD ms   MMD3 ms   MMD3/MMD
```

`n` and `m` are vertices and edges, `m` counted from the pattern rather than assumed, a Matrix
Market file not necessarily carrying every diagonal entry. **`tril(A) = n + m`** is A in its stored
form. **`A+I = 2m`** is what the quotient graph's source pool holds, each edge appearing in both
endpoints' lists. **`C`** is the clique arena, reported by the ordering itself: a size rather than a
capacity, and, this arena never shrinking, also its peak. **`nnz(L)`** includes the diagonal.

The two ratios are the point of the space half. `C/tril(A)` says whether the arena tracked the
INPUT; `C/nnz(L)` says how much the compression bought.

## Where the time actually goes

The set spans five orders of magnitude in ordering time, from 2 microseconds to 15 seconds, so a
count of wins and losses says almost nothing. Grouped by how large the arena grew:

| `C/tril(A)` | matrices | median ratio | share of total time | pooled ratio |
|---|---|---|---|---|
| 0, no off-diagonal at all | 14 | 2.13 | 0.0% | 2.14 |
| 0 to 0.5 | 52 | 1.19 | 0.3% | 0.82 |
| 0.5 to 1.0 | 62 | 1.14 | 0.4% | 0.97 |
| 1.0 to 2.0 | 81 | 1.05 | 1.8% | 0.78 |
| 2.0 to 5.0 | 25 | 0.87 | 1.6% | 0.75 |
| 5.0 and above | 12 | 0.76 | **95.9%** | **0.70** |

**The twelve hardest matrices are 96 percent of all the work, and we are 30 percent faster on
them.** The 163 matrices where we are slower cost **30.9 milliseconds in total, all of them
together**; the 83 where we are faster save **22.0 seconds**. The largest single loss anywhere in
the set is 1.96 ms, on `Schmid/thermal1`.

The ten most expensive matrices are 95 percent of the total:

| matrix | genmmd | ratio | C/tril | C/nnz(L) |
|---|---:|---:|---:|---:|
| FlowIPM22/uni_chimera_i5 | 15.4 s | **0.41** | 52.71 | 0.090 |
| PARSEC/Ga41As41H72 | 13.7 s | 0.85 | 5.78 | 0.008 |
| PARSEC/Si87H76 | 11.9 s | 0.83 | 7.87 | 0.008 |
| JGD_Trefethen/Trefethen_20000 | 10.1 s | 0.77 | 7.27 | 0.024 |
| FlowIPM22/uni_chimera_i1 | 9.7 s | **0.60** | 48.92 | 0.025 |
| JGD_Trefethen/Trefethen_20000b | 6.3 s | 0.87 | 24.14 | 0.084 |
| Pajek/Reuters911 | 1.3 s | 0.81 | 10.07 | 0.723 |
| DIMACS10/vsp_befref_fxm_2_4_air02 | 0.93 s | 0.69 | 9.58 | 0.417 |
| Schenk/nlpkkt80 | 0.80 s | 0.73 | 1.42 | 0.006 |
| Andrews/Andrews | 0.70 s | 0.70 | 5.65 | 0.019 |

Pooled ratio over those ten: **0.70**.

## Where we lose, and why it does not matter

**The 14 matrices with no off-diagonal entries at all.** `Boeing/bcsstm39`, `HB/bcsstm25`,
`Cunningham/m3plates` and eleven others are purely diagonal: the ordering has nothing to do and the
permutation is the identity. genmmd takes 0.259 ms across all fourteen and we take 0.553, so the
ratio is 2.14 and the difference is **294 microseconds for the whole group**. There is no result
here. A routine that does nothing is measuring its own entry and exit, and ours has a
`std::vector` layer around a quotient graph it then never uses.

**The nearly-trivial rest.** `Bai/mhd3200b` at 1.78, `HB/bcsstm25` at 2.10, `GHS_indef/spmsrtls` at
1.79: sub-millisecond orderings on matrices with almost no fill. The pattern across the whole set
is that our fixed overhead, a container layer measured at 17 to 20 percent plus construction of a
quotient graph that a trivial matrix barely touches, is a constant, and a constant is most visible
where the work is smallest.

**That is the honest shape of it and it is not being minimised.** If the question were "order ten
thousand tiny matrices", genmmd would win and we would need to answer for it. It is not the
question this benchmark asks.

## The arena, and two matrices that show its extremes

`C` is what our storage costs, and the reason it is worth watching is that ours is the one scheme
here NOT bounded from the input: genmmd chains dead segments and stays inside `O(n + m)`; we
allocate cliques as they form and never reclaim the space. See `docs/DESIGN_DECISIONS.md`
(2026-08-16) for why that trade is deliberate and what it buys.

**The compression is often enormous.** `PARSEC/Ga41As41H72` factors to 7.08 billion entries and our
arena holds 54.2 million, `C/nnz(L) = 0.008`. `Schenk/nlpkkt80`, a million rows, factors to 3.68
billion and the arena is 21.2 million, 0.006. A clique holds the update part of a supervariable's
column and nothing else, and mass elimination gives one clique per supervariable rather than one
per column, so what we store tracks the input far more closely than the factor.

**And sometimes it is not.** `Pajek/Reuters911` has `C/nnz(L) = 0.723`: the arena is nearly
three-quarters of the factor, and 10.07 times `tril(A)`. `SNAP/email-Enron` is 0.512.
`FlowIPM22/uni_chimera_i5` reaches `C/tril(A) = 52.71`, an arena fifty times the size of the input
it came from. These are the cases the guarantee argument is about: social graphs and interior-point
systems where supervariables barely form, the compression fails, and our space is neither small nor
predictable.

**The striking thing is that those are also where we are fastest.** `uni_chimera_i5` at 0.41,
`uni_chimera_i1` at 0.60, `Reuters911` at 0.81, `vsp_befref` at 0.69. Spending memory to avoid
chaining is exactly the trade being made, and this set says the trade pays where it costs the most.

## Reading the whole table

Two habits, both learned by getting them wrong.

**Do not average the ratio column.** The mean over 246 rows is dominated by matrices that take
microseconds. Weight by time, or read the pooled ratio within a band.

**A ratio hides which side moved.** Where a row looks anomalous, convert to time per vertex or per
edge before believing it. That failure cost a full day on the grid benchmark, where a pattern read
as ours turned out to be an artifact of the vendored routine.

## What this report does not establish

**It is one ordering against one reference, on one branch.** `AMD3` against the vendored AMD
belongs in this table and is not in it yet; the amd branch is currently under active change and
adding it now would date the report within a week.

**It says nothing about ordering QUALITY.** `MMD3` reproduces genmmd's permutation exactly, so fill
is identical by construction, and the interesting question of whether a different ordering would
fill less is not asked here. `PERFORMANCE.md` compares MMD against AMD on that axis.

**And it says nothing about what a full solve costs.** Ordering is one phase and usually a small
one; see `PERFORMANCE.md` for the split.

## Reproducing this

```
cd benchmarks/matrices
python3 ssget.py            # fetch the collection subset into data/
make mmdorder               # this table, all 246 files
```

The private vendored sources are needed for the genmmd columns. Without them the table still runs
and prints ours alone.
