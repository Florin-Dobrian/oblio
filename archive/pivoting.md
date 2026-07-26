# Pivoting

Transcriptions of two algorithms from Ashcraft, Grimes and Lewis, *Accurate Symmetric Indefinite
Linear Equation Solvers*, SIAM J. Matrix Anal. Appl. 20 (1998), 513-561.

The two are independent and answer different questions. Figure 3.3 is an *acceptance test*: given a
candidate pivot, may we use it? Figure 3.4 is a *search*: which candidates are offered, in what
order, and what happens when one is refused? Either test can be dropped into either search, and the
paper's own experiments run this search with three different tests. Keeping them apart is the point
of setting them down separately here.

Neither is Oblio's code. These are the sources, on their own, so that the port can be read against
them.

**Status, 2026-07-26.** The two proposals below have since shipped, so the sections headed "the
nonroot loop" and "the root loop" transcribe code that no longer exists: the kernel was split into
`factorDynamicRootSupernode` and `factorDynamicNonRootSupernode`, the non-root one was reshaped to
Figure 3.4 and lost its symmetric-maximum clause, and the root one was replaced by Figure 2.4. Those
transcriptions are kept because the comparison sections argue against them, and because the
before-and-after is the whole record of why the change was made.

## Notation

The frontal matrix is partitioned

```
A = [ A11  A21^H ]
    [ A21  A22   ]
```

where `A11` holds the *fully assembled* columns, the only ones this front may eliminate, and `A22`
is not available at all. Matrices are indexed from `1`, and `A(i,j)` is an entry.

`gamma(j)` is the largest magnitude among the off-diagonal entries of column `j`, taken over the
**whole** column, `A11` and `A21` alike. A pivot *partner* `q`, by contrast, is always drawn from
`A11`, because an entry of `A22` is not assembled and cannot be pivoted on. That asymmetry is
load-bearing and is the one thing easiest to lose in transcription: the entries of `L` include the
rows lying in `A21`, so a maximum taken over `A11` alone would not bound them.

`alpha` is the pivot threshold, written `alpha-hat` in the paper to distinguish it from the dense
Bunch-Kaufman `alpha = (1 + sqrt(17)) / 8`. It is a free parameter in `(0, 1]`, and the paper
recommends `0.01`. Accepted pivots satisfy `|l_ij| <= 1 / alpha`.

The paper is real symmetric throughout, so a transposed entry is the same entry and the determinant
below is written with `A(q,1)^2` rather than with a conjugate.

## Figure 3.3, an explicit bounding sparse pivot strategy

Stated for column `1` of the reduced matrix, with `q` a partner column already chosen from `A11`.

```
explicit_bounding_test(A, q) -> pivot decision for column 1:

    if gamma(1) == 0:                                   # nothing below the diagonal to eliminate
        return nothing necessary for column 1

    if |A(1,1)| >= alpha * gamma(1):                    # 1x1 on this column's own diagonal
        return 1x1 at 1

    if |A(q,q)| >= alpha * gamma(q):                    # 1x1 on the partner's diagonal
        return 1x1 at q

    det = A(1,1) * A(q,q) - A(q,1)^2

    if max( |A(q,q)| * gamma(1) + |A(q,1)| * gamma(q),  # numerators of the bounds on the two
            |A(1,1)| * gamma(q) + |A(q,1)| * gamma(1)   #   columns of L this block would produce
          ) <= |det| / alpha:
        return 2x2 at 1, q

    return no pivot found                               # repeat the search from the next column
```

The `2 x 2` test is the whole of the paper's contribution on this axis. It is a direct bound on the
entries of `L`, computable from five scalars, and it depends on no `1 x 1` case having failed first.

## Figure 3.4, standard pivot ordering

A refinement of the search in Duff and Reid's MA27. The acceptance test is left open, which is why
the line below names no test.

```
standard_pivot_ordering(A11) -> pivots taken, and columns postponed:

    queue = i1 < i2 < ... < iu,                         # eligible by structure alone
            then iu+1 < iu+2 < ... < iv                 #   postponed at earlier fronts, at the rear

    repeat:
        pivot found = false

        repeat:
            j = the index at the front of the queue
            q = the index of an off-diagonal entry of largest
                  magnitude in column j of A11
            remove j from the queue

            if gamma(j) == 0:                           # j is consumed either way
                nothing necessary for column j

            else if |A(j,j)| >= alpha * gamma(j):
                pivot found = true
                use A(j,j) as a 1x1 pivot

            else if D(j,q) is acceptable as a 2x2 pivot:
                pivot found = true
                remove q from the middle of the queue
                use D(j,q) as a 2x2 pivot

            else:
                add j to the rear of the queue

        until pivot found
           or every column in the queue has been tested since the last successful pivot

    until the queue is empty
       or no pivot was found

where   D(j,q) = [ A(j,j)  A(q,j) ]
                 [ A(q,j)  A(q,q) ]
```

Two columns are examined per candidate and only one `2 x 2` block per column is ever tried, so of
the `m (m - 1) / 2` blocks available in an `A11` of `m` columns, this search tests `m - 1`. All `m`
columns are postponed if none of those `2 m - 1` candidates is accepted.

The rear of the initial queue is where the columns delayed from earlier fronts sit, so the fresh
ones are offered first. The paper calls this the initial skip over previously postponed columns and
reports it effective in both of its orderings.

## Two notes on reading Figure 3.4

The figure's own wording, "an off-diagonal entry of largest magnitude `gamma_j` in column `j` of
`A11`", reads as though `gamma(j)` were the maximum over `A11` only. Section 3.1, which defines the
quantities, and section 3.5, which speaks of the largest entry in a column of `A`, both say
otherwise: the maximum is over the whole column and only the *partner index* is restricted to
`A11`. The transcription above follows the definitions, not the figure's shorthand, because the
alternative reading would not bound the entries of `L` in the update rows.

The `gamma(j) == 0` branch does not set `pivot found`, so an isolated column consumes a queue slot
without counting as progress. Figure 3.6, the exhaustive search of section 3.5, sets `pivot found =
true` in the same branch. The difference is visible in the two figures and the paper does not remark
on it.

## The three scanned quantities

A candidate is judged from three magnitudes, and they differ on two independent axes: which column
is scanned, and how far down. Confusing them is easy, because two of the three are maxima of the
same column.

```
                paper              column scanned   how far down        argmax kept    used by
--------------  -----------------  ---------------  ------------------  -------------  ---------
gammaJ          gamma_1, gamma_j   the candidate j  full block height   no             1x1 and 2x2
frontGammaJ     |a_q1|, |a_qj,j|   the candidate j  front columns only  yes, it is q   2x2 only
gammaQ          gamma_q            the partner q    full block height   no             2x2 only
```

The full height in the first and third rows is what bounds `L`, since a pivot writes factor entries
into the update rows as well as the front, so a maximum stopping at the front would not cover them.
The restriction in the second row is what keeps a partner pivotable: an update row has no column and
no diagonal entry, so it can be measured but not eliminated. That the paper writes the second as an
entry rather than as a maximum, `a_q1` rather than a restricted gamma, is a presentational
difference and not a different quantity; `q` is the argmax of that restricted scan, so the entry at
`(q, j)` has that maximum as its magnitude.

Two asymmetries follow. There is no `frontGammaQ` and no use for one, because the partner is chosen
by scanning the candidate's column, so column `q` is only ever measured and never searched. And
`gammaQ` is confined to the `2 x 2` test in Figure 3.4 and in Oblio, but not in Figure 3.3, whose
second `1 x 1` test reads `gamma_q` outside any `2 x 2` branch.

## Figures 2.4 and 2.5, the dense algorithms

Two figures from the paper's dense half, here because a root front is a dense problem. A root has no
`A22`, so `A11` is the whole matrix, nothing is held back, and section 3's central difficulty, that
an acceptable pivot may not exist in the assembled part, cannot arise.

The property that matters at a root is that neither algorithm can refuse. Both search by chasing the
column maximum, `i <- r` and repeat, and the chase terminates because the off-diagonal magnitude
never decreases along it and strictly increases except in the case that stops it. So every path out
of these figures accepts a pivot. A sparse front may find nothing acceptable and postpone every
column; a root has no parent to postpone to, and needs exactly an algorithm that always succeeds.

`alpha` here is the dense threshold, `(1 + sqrt(17)) / 8`, about `0.6404`, not the sparse `alpha`
of Figures 3.3 and 3.4. Both algorithms bound the factor by `max{1/alpha, 1/(1-alpha)}` and the
condition number of a `2 x 2` block by `(1+alpha)/(1-alpha)`.

```
bounded_bunch_kaufman(A) -> pivot decision:

    if gamma(1) == 0:                                  # nothing below the diagonal to eliminate
        return nothing necessary for column 1

    if |A(1,1)| >= alpha * gamma(1):                   # 1x1 on the first column's own diagonal
        return 1x1 at 1

    i = 1
    repeat:
        r = the index of the first entry of maximum magnitude in column i

        if |A(r,r)| >= alpha * gamma(r):               # 1x1 on the column the chase reached
            return 1x1 at r

        if gamma(i) == gamma(r):                       # the maximum is mutual: a local maximum
            return 2x2 at i, r

        i = r                                          # gamma(r) > gamma(i), so this terminates
```

Figure 2.5 differs in two places. It starts from the largest diagonal in the whole reduced matrix
rather than from the first column, and having done so it drops the `1 x 1` test from the chase,
because the best candidate for a `1 x 1` has already been tried.

```
fast_bunch_parlett(A) -> pivot decision:

    s = the index of the diagonal entry of largest magnitude

    if gamma(s) == 0:
        return nothing necessary for column s

    if |A(s,s)| >= alpha * gamma(s):                   # 1x1 on the largest diagonal
        return 1x1 at s

    i = s
    repeat:
        r = the index of the first entry of maximum magnitude in column i

        if gamma(i) == gamma(r):
            return 2x2 at i, r

        i = r
```

## Choosing between them for a root

Accuracy does not separate them. Section 2.9 reports no appreciable difference and says the choice
should be made on speed in a particular setting, both algorithms carrying the same bounds on the
factor and on the block condition number.

Speed separates them by blocking, and the split reverses. Unblocked at order 400, Table 2.7 gives
17.4 Mflops for Bunch-Kaufman, 16.8 for bounded Bunch-Kaufman and 11.8 for fast Bunch-Parlett.
Blocked, fast Bunch-Parlett wins, because its search usually terminates in a single step and so
discards almost none of the precomputed columns, where the bounded algorithm's extra searches
discard many.

Oblio's front is unblocked. `factorDynamicSupernode` takes one pivot at a time and applies the
rank-1 or rank-2 update immediately, so no column is ever computed and thrown away, which is the
whole basis of the fast Bunch-Parlett advantage. That places us on the side of the table where the
bounded algorithm wins. Two smaller points agree. Fast Bunch-Parlett is `1 x 1`-heavy, 97 percent of
columns against 44 in Table 2.2, and we have a rank-2 kernel that would go unused. And its opening
cost, a full scan for the largest diagonal, buys nothing when nothing is discarded.

One objection to fast Bunch-Parlett does not apply here, worth recording so it is not raised later:
the paper's footnote 7 notes it must maintain the reduced matrix's diagonal separately at `O(n^2)`
extra work, because LAPACK holds that matrix implicitly. Our front is an explicit dense block, so
the diagonal is simply there.

So bounded Bunch-Kaufman for roots, with fast Bunch-Parlett reconsidered if a root front ever proves
hard enough to warrant measuring.

## Oblio notation

Both loops below are transcribed from `src/NumFactorEngine.cpp` in the code's own names, with no
attempt yet to bring them into line with the notation above. That alignment is the next step, and
holding the two vocabularies side by side is what makes it possible to see what has to move.

The names that matter. `jj` is the supernode being factored. `j` is the position of the next pivot
column in its front, so it is the elimination frontier and it advances by one or two as pivots are
accepted. `pivotList` is the candidate queue, holding *global* column indices; `lk1` and `lk2` are
such indices and `k1`, `k2` are their *positions* in `jj`'s front. `jjVal` is the front's value
block, column-major with leading dimension `jjNumNodeIdx`, and `at(r, c)` is the position of row `r`
of column `c` within it. `jjPreFactorFrontSize` is the number of fully assembled columns and
`jjNumNodeIdx` is the full block height, front columns plus update rows. `threshold` is
`mPivotThreshold`.

The two passes are selected by `mUpdateSize[jj]`, which is to say by whether the supernode is a
root, and a supernode runs one or the other and never both.

## Oblio, the nonroot loop

`mUpdateSize[jj] > 0`, so a parent exists and a refused column can be delayed to it.

```
factorDynamicSupernode(nf, jj), pass 2:

    pivotList = the front's global column indices, in front order
    j         = 0

    while pivotList is not empty:
        pivotFound = false
        trials     = pivotList.size()

        while trials > 0:
            lk1 = pivotList.pop_front()
            k1  = gblToLcl[lk1]

            # lk1's largest off-diagonal magnitude: its row to the left, its column below
            max1 = max of |jjVal[at(k1, i)]|  over i = j .. k1-1
                      and |jjVal[at(i, k1)]|  over i = k1+1 .. jjNumNodeIdx-1
            diagonal1 = jjVal[at(k1, k1)]

            if max1 == 0:                                   # isolated column, nothing to eliminate
                pivotFound = true
                swap(jj, j, k1) if j != k1
                rank-- if |diagonal1| == 0
                pivotType[lk1] = 1
                j += 1
                break

            if |diagonal1| > 0 and |diagonal1| >= threshold * max1:             # accept 1x1
                pivotFound = true
                applyPivot1x1(nf, jj, j, k1, lk1, ...)
                j += 1
                break

            else if pivotList is not empty:                 # try a 2x2 with lk1 and a front partner
                # the partner scan: max1's scan again, stopped at the end of the front
                max, k2 = max of |jjVal[at(k1, i)]|  over i = j .. k1-1
                             and |jjVal[at(i, k1)]|  over i = k1+1 .. jjPreFactorFrontSize-1,
                          with k2 the position attaining it
                lk2  = jjNodeIdx[k2]
                max2 = max of |jjVal[at(k2, i)]|  over i = j .. k2-1
                          and |jjVal[at(i, k2)]|  over i = k2+1 .. jjNumNodeIdx-1

                d      = readPivotBlock2x2(jjVal, ld, k1, k2)
                maxmax = max( |d.d22| * max1 + max * max2,
                              |d.d11| * max2 + max * max1 )

                if (max == max1 and max == max2 and max != 0)                   # accept 2x2
                   or (|d.det| > 0 and |d.det| >= threshold * maxmax):
                    pivotFound = true
                    pivotList.remove(lk2)
                    applyPivot2x2(nf, jj, j, k1, k2, lk1, lk2, ...)
                    j += 2
                    break
                else:
                    pivotList.push_back(lk1)                # delay lk1

            else:
                pivotList.push_back(lk1)                    # no partner available, delay lk1

            trials -= 1

        if not pivotFound:
            break

    # whatever is left in pivotList is delayed to the parent, and frontSize contracts by that many
```

The first clause of that `2 x 2` condition, `max == max1 and max == max2 and max != 0`, is not in
Figure 3.3. It was carried from 0.9 and removed on 2026-07-26, so the transcription above is of the
code as it stood before that. It reached for the Bunch-Parlett condition of AGL section 2.3, which
requires the off-diagonal to be a local maximum in both columns *and* both diagonals to be at most
`alpha` times it. Reaching the clause establishes only the first of those for the candidate, from
the `1 x 1` test that failed; the partner's diagonal is never tested, and without it the bound on
`L` does not hold. A six by six banded matrix with a zero diagonal beside one of 5.5e5 drives it:
the clause accepted a block the growth bound rejects, `L` took an entry of 1.98e5 and the residual
came out at 1.17e-11 against 4.8e-16 once the clause was gone. The root loop below keeps the same
predicate, where the chase supplies the missing condition through control flow rather than as a
conjunct.

## Oblio, the root loop

`mUpdateSize[jj] == 0`, so there is no parent and a refused column has nowhere to go. Note that
`jjNumNodeIdx == jjPreFactorFrontSize` here, so a scan stopped at the end of the front and one run
to the full block height are the same scan, and the loop keeps only one of them.

```
factorDynamicSupernode(nf, jj), pass 1:

    pivotList = the front's global column indices, in front order
    j         = 0

    while pivotList is not empty:
        pivotFound = false
        trials     = pivotList.size()

        while trials > 0:
            lk1 = pivotList.pop_front()
            k1  = gblToLcl[lk1]

            if pivotList is empty:                          # last candidate standing, forced 1x1
                pivotFound = true
                rank-- if |jjVal[at(k1, k1)]| == 0
                pivotType[lk1] = 1
                j += 1
                break

            # lk1's largest off-diagonal magnitude: its row to the left, its column below
            max1, k2 = max of |jjVal[at(k1, i)]|  over i = j .. k1-1
                          and |jjVal[at(i, k1)]|  over i = k1+1 .. jjNumNodeIdx-1,
                       with k2 the position attaining it
            diagonal1 = jjVal[at(k1, k1)]

            if max1 == 0:                                   # isolated column, nothing to eliminate
                pivotFound = true
                swap(jj, j, k1) if j != k1
                rank-- if |diagonal1| == 0
                pivotType[lk1] = 1
                j += 1
                break

            if |diagonal1| > 0 and |diagonal1| >= threshold * max1:             # accept 1x1
                pivotFound = true
                applyPivot1x1(nf, jj, j, k1, lk1, ...)
                j += 1
                break

            else:                                           # try a 2x2 with lk1 and its max partner
                lk2  = jjNodeIdx[k2]
                max2 = max of |jjVal[at(k2, i)]|  over i = j .. k2-1
                          and |jjVal[at(i, k2)]|  over i = k2+1 .. jjNumNodeIdx-1

                if max1 == max2:                            # accept 2x2, on the magnitudes alone
                    pivotFound = true
                    pivotList.remove(lk2)
                    applyPivot2x2(nf, jj, j, k1, k2, lk1, lk2, ...)
                    j += 2
                    break
                else:
                    pivotList.push_back(lk1)

            trials -= 1

        if not pivotFound:
            break
```

Two differences in shape are worth having in view before the notation moves. The root loop has a
forced `1 x 1` at the top, taken when `lk1` is the only candidate left, and it performs no
elimination and no swap, because that column is already at position `j` with nothing below it. And
where the nonroot loop runs three scans, `max1` over the full height, `max` over the front, and
`max2`, the root loop runs two: `k2` is the position attaining `max1` itself, so the partner is
whatever the full scan found rather than a separately chosen front column.

## Figure 3.4 against the nonroot loop

The two agree on everything that decides an outcome: the queue, the order of the tests, `1 x 1`
before `2 x 2`, a refused column sent to the rear, and an inner loop that gives up once every
candidate has been tried since the last acceptance. What differs is shape, and the differences are
of a kind that could be removed without changing a single decision.

**The partner is chosen late.** Figure 3.4 computes `qj` in the same breath as `gamma_j`, right
after the candidate is popped, and from then on the candidate and its partner are both in hand. The
nonroot loop computes `max1` at the top but leaves the partner scan inside the `2 x 2` branch, so
`k2` does not exist until two tests have already been decided.

**One column is walked twice.** Because the partner scan is a separate scan, it re-reads entries
that `max1`'s scan has already touched: the row segment is walked twice in full, and the column
segment twice down to the end of the front. A single pass produces `max1`, `max` and `k2` together,
because the front-restricted maximum is a prefix of the full-height one. That is not an
approximation, it is the same arithmetic in one traversal instead of two.

**The `2 x 2` branch is guarded by the queue rather than by the scan.** The nonroot loop asks
whether `pivotList` is still nonempty, where Figure 3.4 simply has a partner in hand. The two
conditions coincide: `k2` is set by the first entry either scan range examines, so `k2 < 0` exactly
when both ranges are empty, which happens only when `k1 == j == jjPreFactorFrontSize-1`, which is
exactly when every other front column has already been eliminated and the queue is empty after the
pop. Asking the scan is the more local of the two questions and needs nothing outside the candidate.

**The delay appears twice.** Figure 3.4 has one chain of four outcomes and one place that pushes to
the rear. The nonroot loop has an `if` for the isolated column, then a chain whose `2 x 2` case is
nested one level deeper, so the same `push_back(lk1)` is written at two sites and `trials -= 1`
serves both.

**The acceptance test is inline.** Figure 3.4 writes `if D(j,q) is acceptable as a 2 x 2 pivot`,
a hook, which is why the same search accommodates three different tests in the paper's experiments.
The nonroot loop spells `max2`, `readPivotBlock2x2`, `maxmax` and a two-part condition into the loop
body, and that is most of what makes the loop long enough to obscure its shape.

## Oblio, a proposed nonroot loop

Structure only, measured against Figure 3.4. **The restructure was behavior-preserving, exactly**:
the same pivots, in the same order, from the same fronts, which the pinned pivot and delay counts
confirmed when it landed. What the acceptance test contains is a separate question, settled
separately afterwards.

The names follow the figure. `lj` and `lq` are the candidate and its partner as global column
indices, which is what Figure 3.4 calls `j` and `q_j`; the subscript is unnecessary here because
there is one candidate in hand at a time and `lq` is by construction the partner of that one. `j`
and `q` are their positions in the front, which the figure has no need of and an implementation
cannot do without. `nextPivot` is the frontier, the count of front columns already factored and
equally the position the next accepted pivot is swapped down to. The figure has no counterpart for
it either, because it works on a matrix already permuted and leaves `P A P^T` to record silently
what `nextPivot` and `swap` record explicitly.

```
factorDynamicSupernode(nf, jj), pass 2, proposed:

    pivotList = the front's global column indices, in front order
    nextPivot = 0

    while pivotList is not empty:
        pivotFound = false
        trials     = pivotList.size()

        while trials > 0:
            lj = pivotList.pop_front()
            j  = gblToLcl[lj]

            gammaJ, frontGammaJ, q = scanPivotColumn(jj, nextPivot, j)
            diagonalJ              = jjVal[at(j, j)]

            if gammaJ == 0:                            # isolated column, nothing to eliminate
                pivotFound = true
                swap(jj, nextPivot, j) if nextPivot != j
                rank-- if |diagonalJ| == 0
                pivotType[lj] = 1
                nextPivot += 1
                break

            else if |diagonalJ| > 0 and |diagonalJ| >= threshold * gammaJ:
                pivotFound = true
                applyPivot1x1(nf, jj, nextPivot, j, lj, ...)
                nextPivot += 1
                break

            else if q >= 0 and acceptPivot2x2(nf, jj, nextPivot, j, q, gammaJ, frontGammaJ):
                pivotFound = true
                lq = jjNodeIdx[q]
                pivotList.remove(lq)
                applyPivot2x2(nf, jj, nextPivot, j, q, lj, lq, ...)
                nextPivot += 2
                break

            else:
                pivotList.push_back(lj)                    # delay lj
                trials -= 1

        if not pivotFound:
            break

    # whatever is left in pivotList is delayed to the parent, and frontSize contracts by that many
```

with two helpers lifted out of the body. The scan, which walks each entry of `lj`'s line once:

```
scanPivotColumn(jj, nextPivot, j) -> gammaJ, frontGammaJ, q:

    frontGammaJ = -1;  q = -1

    for i = nextPivot .. j-1:                              # lj's row, all of it front columns
        if frontGammaJ < |jjVal[at(j, i)]|:  q = i;  frontGammaJ = |jjVal[at(j, i)]|

    for i = j+1 .. jjPreFactorFrontSize-1:                 # its column, front part
        if frontGammaJ < |jjVal[at(i, j)]|:  q = i;  frontGammaJ = |jjVal[at(i, j)]|

    gammaJ = frontGammaJ                                   # the front maximum is a prefix of the
                                                           #   full-height one, so continue from it
    for i = jjPreFactorFrontSize .. jjNumNodeIdx-1:        # its column, update rows
        if gammaJ < |jjVal[at(i, j)]|:  gammaJ = |jjVal[at(i, j)]|

    return gammaJ, frontGammaJ, q
```

and the acceptance test, which is the hook Figure 3.4 leaves open:

```
acceptPivot2x2(nf, jj, nextPivot, j, q, gammaJ, frontGammaJ) -> bool:

    gammaQ = max of |jjVal[at(q, i)]|  over i = nextPivot .. q-1
                and |jjVal[at(i, q)]|  over i = q+1 .. jjNumNodeIdx-1

    d      = readPivotBlock2x2(jjVal, ld, j, q)
    maxmax = max( |d.d22| * gammaJ + frontGammaJ * gammaQ,
                  |d.d11| * gammaQ + frontGammaJ * gammaJ )

    return |d.det| > 0 and |d.det| >= threshold * maxmax
```

The test above is Figure 3.3's fourth branch and nothing else. The transcription of the current loop
higher up still shows a second clause, `max == max1 and max == max2 and max != 0`, carried from 0.9
and accepted alongside the growth bound. It was removed on 2026-07-26, and the two-line note under
that transcription says why: it reached for the Bunch-Parlett condition of section 2.3, which needs
both diagonals bounded as well as the mutual maximum, and only the candidate's diagonal is known
small where the clause sits.

On efficiency, the proposal is ahead where it matters and level elsewhere. The fused scan replaces
two traversals of `lj`'s line with one whenever a `2 x 2` is considered, which is the expensive
path, and the three loops stay branch-free because the front part and the update part are separate
loops rather than one loop with a test in it. Where the `1 x 1` succeeds the fused scan does track
`q` and a second running maximum that the current code would not have computed, but over entries it
is loading anyway, so the cost is a comparison and a store and no extra memory traffic. Nothing is
computed earlier than it is needed except that: `gammaQ` and the `2 x 2` block are still read only
inside the test, exactly as now, and `lq` is still looked up only on acceptance.

The scan also makes the relation between the two structural. `frontGammaJ` is computed first and
`gammaJ` continues from it rather than starting over, so `frontGammaJ <= gammaJ` follows from the
loop bounds rather than standing as an invariant to be argued across two separate scans.

Three things in the loop have no counterpart in Figure 3.4, and all three are implementation detail
rather than shape. The `q >= 0` conjunct covers the case the figure does not show, where the front
holds no partner to pair with. The isolated-column arm sets `pivotFound`, where the figure leaves it
alone. And the figure's nested `else { if acceptable ... else rear }` is flattened into a fourth arm
of the one chain.

`frontGammaJ` follows `nodeIdx` and `frontNodeIdx`, and `frontSize`: the same quantity as `gammaJ`,
narrowed to the front. It is the paper's `|a_qj,j|`, and the two descriptions agree because `q` is
the argmax of the front-restricted scan, so the entry at `(q, j)` has that maximum as its magnitude
by definition.

One name is still open. `maxmax` is the quantity the determinant is tested against, and its name was
formed from `max1` and `max2`, which no longer exist. The code's comment calls it the growth bound,
which is the Duff-Reid motivation; AGL's is that the two terms are the numerators of the bounds on
the entries of `L`, so following the paper the name would say `L` rather than growth.

The `2 x 2` block needs nothing, since `readPivotBlock2x2` already returns `d11`, `d22`, `d21`,
`d12` and `det`. That is `D`, the object the paper also writes as `D_{k,k+1}`, and it is the right
letter because it names what those four entries become.

## Oblio, a proposed root loop

**This one is not behavior-preserving**, unlike the nonroot proposal. It changes which pivots are
chosen, and the pinned tier-1 pivot counts will move. That is the point of it: the current root
acceptance is the local-maximum condition with the diagonal conditions dropped, and dropping them is
what leaves the entries of the factor unbounded.

The shape is Figure 2.4 rather than Figure 3.4, and the queue goes with it. A root cannot delay a
column, so `pivotList`, `trials`, the push to the rear and the outer sweep are all vestigial there:
nothing is ever postponed, and the chase reaches a pivot from any starting column. What is left is
take the next unfactored column, chase, accept, advance.

`scanPivotColumn` serves unchanged. At a root `jjNumNodeIdx == jjPreFactorFrontSize`, so its third
loop is empty and `frontGammaJ == gammaJ` identically; the root ignores the second output and uses
the third, `q`, which is the chase's `r`. Figure 2.4's `i` and `r` are our `j` and `q`.

```
factorRootSupernode(nf, jj), proposed, bounded Bunch-Kaufman:

    nextPivot = 0

    while nextPivot < jjPreFactorFrontSize:
        j = nextPivot
        lj = jjNodeIdx[j]
        gammaJ, _, q = scanPivotColumn(jj, nextPivot, j)

        if gammaJ == 0:                                # isolated column, nothing to eliminate
            rank-- if |jjVal[at(j, j)]| == 0
            pivotType[lj] = 1
            nextPivot += 1
            continue

        if |jjVal[at(j, j)]| >= alphaRoot * gammaJ:     # 1x1 on this column's own diagonal
            applyPivot1x1(nf, jj, nextPivot, j, lj, ...)
            nextPivot += 1
            continue

        repeat:                                        # the chase
            gammaQ, _, qNext = scanPivotColumn(jj, nextPivot, q)
            lq = jjNodeIdx[q]

            if |jjVal[at(q, q)]| >= alphaRoot * gammaQ: # 1x1 on the column the chase reached
                applyPivot1x1(nf, jj, nextPivot, q, lq, ...)
                nextPivot += 1
                break

            if gammaJ == gammaQ:                       # mutual maximum: accept the 2x2
                applyPivot2x2(nf, jj, nextPivot, j, q, lj, lq, ...)
                nextPivot += 2
                break

            j = q;  lj = lq;  gammaJ = gammaQ;  q = qNext
```

What the change buys, in the order the guarantees fall out.

No delays, structurally. There is no branch that fails to produce a pivot, so a root front is always
factored completely, and it is not a probabilistic claim: `gamma` never decreases along the chase
and strictly increases except in the case that stops it, so no column is visited twice.

The `2 x 2` block is nonsingular by construction. The chase reaches its `2 x 2` only after both
diagonal tests have failed, `|A(j,j)| < alphaRoot * gammaJ` from how it arrived and
`|A(q,q)| < alphaRoot * gammaQ` from the branch above, while `|A(q,j)| == gammaJ == gammaQ`. So
`|det| > gammaJ^2 (1 - alphaRoot^2)`, bounded away from zero. Those two failed tests are exactly the
guard the current `gammaJ == gammaQ` condition lacks, and they are recovered here not as a patch but
because the chase cannot reach the acceptance without them.

The forced `1 x 1` disappears. The current loop accepts its last remaining candidate whatever it
looks like, which is where the open question about dividing by zero comes from. Here there is no
last-candidate case, because acceptance does not depend on candidates remaining.

And the block is well conditioned, `kappa(D) < (1 + alphaRoot) / (1 - alphaRoot)`, which is the
condition under which the paper allows explicit inversion. Cramer's rule in `applyPivot2x2` is
therefore sound at roots, though not at nonroots.

Two things this needs that the nonroot path does not. `alphaRoot` is the dense threshold,
`(1 + sqrt(17)) / 8`, not `mPivotThreshold`. The strong value is free at a root, since the front is
dense, a symmetric permutation of it stays dense, and no column can be delayed, so there is no fill
to trade against and the only cost is comparisons. But it does mean two thresholds in one
factorization, and the name above is a placeholder. And the scan is now called once per chase step
rather than once per candidate, so the worst case is `O(m)` traversals per pivot; the paper's
Appendix C bounds the expected number of column searches at about `e`, under 2.72, and Table 2.6
measures 2.44 on random matrices of order 400.
