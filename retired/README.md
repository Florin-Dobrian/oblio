# Retired: the earlier ladder layers

`Mmd1`, `Mmd2`, `Amd1` and `Amd2`, with their headers. Retired 2026-08-21, still here, out of the
build and out of the `Ordering` enum.

**WHAT THEY WERE.** Ours, over the shared `QuotientGraph`, each carrying a base algorithm without
its reference's later refinements, so each ordered differently and none was a drop-in replacement
for a vendored routine:

| | what it had |
|---|---|
| `Mmd1` | the batch alone: multiple elimination, mass elimination, natural absorption |
| `Mmd2` | plus the prepass, the q2h refresh, pairwise merging and outmatching |
| `Amd1` | the approximate bound alone |
| `Amd2` | plus aggressive absorption and hash detection |

`Mmd3` and `Amd3` are the tops of those two ladders and are what remains in production. Each
reproduces its reference's permutation exactly on every matrix the benchmarks cover, which is what
makes them worth keeping and what these four were never trying to be.

**WHY THEY WERE RETIRED, and it was not that they were wrong.** The flat `QuotientGraph` was
serving six drivers across THREE list-order conventions: these four on the defaults, `Mmd3` on
genmmd's reversed incidence walk, `Amd3` on `AMD_2`'s cliques-first order. That third convention is
ours and predates both vendored ones, and its existence blocked aligning the two that matter: the
class had to select a convention with a flag where the compacted class uses a suffixed pair per
branch. With these four gone the flat class serves two drivers and two conventions, the same as the
compacted one, and the two can be brought into line. See `docs/QUOTIENT_GRAPH_USAGE.md`.

**THEY WILL NOT COMPILE FOR LONG.** `97f4bc6` is the last commit where they build. They depend on
`QuotientGraph`'s default convention and on `mVendoredListOrder`, both of which are on their way
out, so treat these files as a record of what the layers were rather than as code to resurrect
unchanged.

**AND THE SAME ALGORITHMS ARE ALIVE ELSEWHERE.** `experiments/ordering/` carries the whole ladder as
working C++ and Python twins, `md1` through `md5`, `mmd1` to `mmd3`, `amd1` to `amd4`, each checked
against the other and against production by `make test` in that directory. That is where the
pedagogy lives and where it is maintained. These four were a second copy of it in production form.

**WHAT WENT WITH THEM.** Four enumerators, `MMD1`, `MMD2`, `AMD1` and `AMD2`, so `Ordering` is now
five values. Sixteen validity assertions in `test_order` and twelve more across the other suites,
which is why the suite reports 251 where it reported 279. The digest went from nine drivers to
five, 365 digests where it recorded 657, and `.digest-baseline` was re-recorded. Tests and examples
that named `Ordering::MMD2` because they wanted *an* ordering now name `Ordering::MMD3`.

**IF THEY COME BACK.** The question they would answer is what each mechanism costs in isolation,
which is a real question and is currently answered on grids only, in `experiments/ordering`. Doing
it on real matrices would need these four building again, against whatever the quotient graph looks
like then.
