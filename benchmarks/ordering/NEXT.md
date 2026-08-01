# Next: the incidence lists off std::vector

A handoff note, written 2026-08-01 at the end of a long session. It says where the ordering work
stands, what the next change is, and why, so that the next session can start on the change rather
than on rediscovering the reasoning. Delete it once the change lands and its result is recorded in
this folder's README.

## Where things stand

Seven ordering methods. `Natural`, the vendored `MMD` and `AMD`, and ours: `MMD1`, `MMD2`, `AMD1`,
`AMD2`, all over the shared `QuotientGraph`. 213 assertions pass. `experiments/ordering`'s
`make test` checks each of ours against the prototype it came from, by permutation for `mmd1`,
`mmd2` and `amd1`, and by nnz(L) for `amd2`, which skips the postorder its prototype ends with.

alpamayo, 140x140 grid, milliseconds:

```
MMD   1.16    MMD1  3.46 (2.98x)   MMD2  2.28 (1.97x)
AMD   1.27    AMD1  2.42 (1.90x)   AMD2  3.81 (3.00x)
```

Three of ours in the 1.9 to 2.0 band. AMD2 is the outlier: on a grid its mechanisms fire exactly as
often as the vendored routine's (1 clique absorbed, 2488 hash merges, measured both sides) while
costing us the hash pass, so it pays for something it does not get back here.

## The next change

**Get `I[u]` off `std::vector`.** A Time Profiler trace of `amd1` at 140x140, 7.08 s over 3000
orderings:

```
eliminate                      3.35 s   47.3%
  push_back                    1.12 s   15.9%   <- the largest single line in the program
  reachableSet                 666 ms    9.4%
  eliminate self               637 ms    9.0%
  iterator increment           320 ms
  insert / resize              335 ms
orderAmd1 self                 1.19 s   16.8%
constructor                    644 ms    9.1%
Buckets file/unfile/refile     686 ms    9.7%
destructor                     259 ms    3.7%
```

There is exactly one `push_back` in that hot path: rebuilding `I[u]` in `eliminate`, which compacts
in place and then appends the pivot, once per reached vertex per elimination, on a vector that may
reallocate. `amd_2` never appends, because a vertex's adjacency and element list share one list and
the new element is written during the same pass that rewrites it.

So the incidence lists want to be written rather than appended to. `A[u]` is already one flat array,
since it only shrinks; `C[c]` is a bump arena, since its size is known when it is formed. `I[u]` is
the one that grows, which is why it was left alone, and is now the measured target rather than a
plausible one. Expect something on the order of 15 percent, with `insert` and `resize` likely to
follow.

Design notes that are already settled and should not be relitigated: the arena's offsets must be
indices rather than pointers, so growth invalidates nothing; nothing needs reclaiming, since the
totals here are bounded by the factor's pattern, which has to be resident anyway; and AMD's
copy-to-the-end with compaction exists to honor a caller-supplied fixed workspace, which we do not
promise.

## What NOT to try, all measured this session

Each of these was plausible, was tried, and returned nothing. `../README.md` has the detail.

- **Fusing the driver's passes.** AMD1 walks `A[u]` and `I[u]` three times per elimination against
  the vendored scan's one, which looks like the whole gap and is not: fusing them measured 1.44x
  against 1.45x, nothing. Fusion removes loop setup and a re-fetch of warm data, not per-element
  work. **Do not restructure the driver-eliminator seam on the strength of the visit count.**
- **A flat weight array** for locality: no difference (it is in the tree now for a different
  reason, since the supervariable members became a chain).
- **Reserving the incidence buffers** up front: moved allocations rather than removing them.
- **Narrowing the per-vertex arrays** to `int32_t`: about 2 percent, inside the noise.
- **Writing the adjacency through a pointer** instead of `push_back` at construction: slower,
  `resize` zero-fills.

## How to measure it

`../README.md` is the method. In short: `make run` for the timing table, and Instruments for
anything that claims to be faster.

```
make profile
xcrun xctrace record --template 'Time Profiler' --launch -- ./order_profile_cpp amd1 140 3000
xcrun xctrace record --template 'CPU Counters' --launch -- ./order_profile_cpp amd1 140 3000
```

Time Profiler for where, CPU Counters for why: its useful-cycle percentage separates doing more
work from doing it less efficiently, which is the distinction that settled the MMD question. Both
branches' remaining gap is work rather than stalling, so a change that does not remove work should
not be expected to pay.

Correctness is the harness: `experiments/ordering`'s `make test` must still report all four of ours
agreeing with their prototypes, and the suite must still be at 213. A change to `QuotientGraph`
touches every driver, so both checks matter more than usual here.

## After that

- **The remaining 9 percent in the constructor** is mostly first touch of freshly allocated memory,
  which no loop rewriting reaches. It shrinks only by needing less memory, and moving `I[u]` off
  per-vertex vectors removes 19600 vector headers from it at n = 19600.
- **`Buckets` at 10 percent** has not been looked at at all.
- **Widening the test set.** Every number in this folder comes from grid Laplacians. The
  experiment's own open items say the test set is thin, and the AMD2 result in particular, where
  mechanisms cost more than they return, is exactly the kind of finding that might reverse on a
  problem with real supernodal structure.
- **Dense row handling** is still absent from all of ours, which the vendored AMD has. Left out for
  want of evidence rather than on principle: no matrix measured so far has a dense row.
