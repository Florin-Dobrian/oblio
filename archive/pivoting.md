# Pivoting

A record of the pivoting strategies in Ashcraft, Grimes and Lewis, *Accurate Symmetric Indefinite
Linear Equation Solvers*, SIAM J. Matrix Anal. Appl. 20 (1998), 513-561, and of what Oblio takes
from them.

Four algorithms matter here, and they answer two different questions. An *acceptance test* asks,
given a candidate pivot, whether it may be used: Figure 3.3 is one. A *search* asks which candidates
are offered, in what order, and what happens when one is refused: Figures 2.4, 2.5 and 3.4 are
searches. The two axes are independent, and the paper's own experiments run one search with three
different tests, which is the reason for keeping them apart throughout.

The dense algorithms come first here, against the paper's own order, because they are the simpler
case and because Oblio's root kernel is one of them. The sparse ones follow, then the material that
belongs to neither, then Oblio's two loops.

**Oblio uses both halves.** Root supernodes run bounded Bunch-Kaufman, Figure 2.4, because a root
cannot delay a column and needs a search that cannot refuse. Non-root supernodes run Figure 3.4's
search with Figure 3.3's test, and delay what they cannot pivot.

## The paper's notation

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


## Dense algorithms

### Figures 2.4 and 2.5

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

### Choosing between them for a root

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


## Sparse algorithms

### Figure 3.3, an explicit bounding sparse pivot strategy

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

### Figure 3.4, standard pivot ordering

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

### Two notes on reading Figure 3.4

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

### The three scanned quantities

A candidate is judged from three magnitudes, and they differ on two independent axes: which column
is scanned, and how far down. Confusing them is easy, because two of the three are maxima of the
same column.

```
                paper              column scanned   how far down        argmax kept    used by
--------------  -----------------  ---------------  ------------------  -------------  ---------
jGamma          gamma_1, gamma_j   the candidate j  full block height   no             1x1 and 2x2
jFrontGamma     |a_q1|, |a_qj,j|   the candidate j  front columns only  yes, it is q   2x2 only
qGamma          gamma_q            the partner q    full block height   no             2x2 only
```

The full height in the first and third rows is what bounds `L`, since a pivot writes factor entries
into the update rows as well as the front, so a maximum stopping at the front would not cover them.
The restriction in the second row is what keeps a partner pivotable: an update row has no column and
no diagonal entry, so it can be measured but not eliminated. That the paper writes the second as an
entry rather than as a maximum, `a_q1` rather than a restricted gamma, is a presentational
difference and not a different quantity; `q` is the argmax of that restricted scan, so the entry at
`(q, j)` has that maximum as its magnitude.

Two asymmetries follow. There is no `qFrontGamma` and no use for one, because the partner is chosen
by scanning the candidate's column, so column `q` is only ever measured and never searched. And
`qGamma` is confined to the `2 x 2` test in Figure 3.4 and in Oblio, but not in Figure 3.3, whose
second `1 x 1` test reads `gamma_q` outside any `2 x 2` branch.


## Extras

### Solving the 2 x 2 systems, and why complete pivoting is fine here

A separate axis from selection, and easy to conflate with it. Once a `2 x 2` block is accepted, two
places solve a system against it: the factorization forms the pair of `L` columns as `[l1 l2] D =
[t1 t2]` for every row below the block, and the solve later applies `D` to the right-hand side. The
paper is emphatic about how, and its statements look contradictory until the regimes are separated.

**Unbounded `L`.** Bunch-Kaufman does not bound the entries of `L`, section 2.2 being the
demonstration. Stability there rests on large entries of `L` always being paired with small entries
of `D`, so that their products stay bounded, and Higham's analysis needs the `2 x 2` solve to be
**componentwise** backward stable to preserve that pairing. The scaled Cramer's rule of LINPACK and
LAPACK qualifies, and so does Gaussian elimination with partial pivoting. Complete pivoting and the
SVD do not, and this is where the paper's alarming sentence comes from: "Implementations using
complete pivoting are unstable, while partial pivoting is provably stable." Appendix A derives the
SVD case, exhibiting a reduced matrix computed as `O(1) - O(10^t) - O(10^t)` with the cancellation
eating the result. The complete-pivoting case is asserted alongside it rather than derived.

**Bounded `L`.** Every other algorithm in the paper bounds it, and the requirement then drops to
**normwise** backward stability, which is the hypothesis Appendix B's proof is built on. Partial
pivoting, complete pivoting and QR all qualify. This is why AGL can write, with no contradiction,
that their own sparse codes use Gaussian elimination with complete pivoting: their scheme bounds
`L`, so the stronger condition is not required.

So the regimes, and which algorithm sits in which:

```
                             bounds L        needs                which solves are admissible
---------------------------  --------------  -------------------  -------------------------------
Bunch-Kaufman                no              componentwise        scaled Cramer, partial pivoting
bounded BK, fast BP (2.4-5)  yes             normwise             any of them
Duff-Reid / AGL 3.3          yes             normwise             any of them
```

**Both of Oblio's paths are in the bounded row.** Non-roots bound `L` by `1 / alpha` through the
Figure 3.3 test; roots bound it by `max{1/alpha, 1/(1-alpha)}` through bounded Bunch-Kaufman, that
bound being the whole content of the word "bounded" in its name. Appendix A's warning therefore
does not reach either path, and partial and complete pivoting are both valid choices for us.
Nothing in the stability argument would decide between them.

One asymmetry survives the regimes, and it concerns Cramer's rule rather than the pivoting. Cramer
is admissible under a *bounded condition number* of `D`, which is a stronger property than a
bounded `L` and which the two paths do not share. The chase supplies it at a root, where reaching a
`2 x 2` means both diagonals have already failed their own tests, giving `kappa(D) < (1 + alpha)/(1
- alpha)` by the Gerschgorin argument of section 2.4. Figure 3.3 does not supply it at a non-root,
and page 29 says so outright: the test bounds `L` but not the conditioning, which is exactly why
that page rules Cramer out for this family.

Oblio used Cramer in `factor2x2` and partial-pivoted LU in `diagonalDynamic` until 2026-07-26,
which is the arrangement the paper warns against, unstably formed factor and stably applied one.
Both now use partial-pivoted LU. Partial is not required over complete; it was chosen so the two
halves of the same system match, and because it remains valid should some later path ever admit an
unbounded `L`.

### Why the sparse search must keep refused candidates, and the dense one need not

Figures 2.4 and 2.5 hold no candidate list. They take a column, chase, and accept. Figures 3.4 and
3.6 both maintain a queue, rotate refused columns to the rear, and carry a termination condition
counting how many have been tried since the last acceptance. That difference is not a matter of
taste, and it is worth following precisely, because it is the deepest structural consequence of
sparsity in this whole area. What the difference does *not* settle is the data structure, which the
section after this one takes up separately.

**The dense chase terminates because one scan produces both of its outputs.** In Figure 2.4, `r` is
where column `i`'s largest off-diagonal sits and `gamma_i` is that entry's magnitude. Same scan,
same entry. The monotonicity is then one line: the entry `a_ri` has magnitude `gamma_i` and it lies
in column `r`, so it is a lower bound on column `r`'s own maximum, giving `gamma_r >= gamma_i`. The
chase never decreases, so it never revisits a column, so it terminates, and it terminates on
`gamma_i == gamma_r`, which is a local maximum and acceptable by construction. **No path out of the
figure refuses**, and an algorithm that cannot refuse has nothing to remember: the frontier alone
suffices.

**Sparsity splits the argmax from the bound, and the proof goes with it.** The partner must be a
column of `A11`, since an update row has no column and no diagonal, but the bound has to run the
full column height, since the entries of `L` this pivot writes reach into the update rows. So the
two come from different ranges: `q` is the argmax of `jFrontGamma` while the quantity bounding `L`
is `jGamma`. What survives of the chain is

```
qGamma >= |a_qj| = jFrontGamma        and        jFrontGamma <= jGamma
```

with nothing relating `qGamma` to `jGamma`. The maximum can fall as the chase steps, so the chase
can cycle.

Neither repair works, which is what makes this a barrier rather than an inconvenience. Chase on
`jFrontGamma` and it terminates, since that quantity is a genuine maximum over a fixed set, but it
terminates at a maximum local to the front only, and that bounds nothing in the update rows, which
is the case AGL make in section 3.3 against carrying the dense solution over. Chase on `jGamma` and
the column holding it may lie in `A21`, where there is no column to jump to.

**So refusal becomes possible, and refusal is temporary.** If a refused column were refused for
good, the loop would delay it and move on: one pass, a list of the delayed, nothing to remember.
But every acceptance applies a trailing update, so the values in every remaining column change, and
a column refused a moment ago may pass afterwards. Something must therefore keep the candidates
still in play and offer them again, and must know when to stop offering: when every one has been
tried since the last acceptance, no acceptance is possible without new values and none are coming.

### What the queue is actually for

Keeping candidates in play does not require a queue, and it is worth separating what the structure
is needed for from what it merely happens to provide.

**Not termination, and not correctness.** The unfactored columns are always the contiguous range
`[nextPivot, frontSize)`, because every acceptance swaps its columns down to the frontier. So the
same loop can be written on indices alone:

```
while nextPivot < frontSize:
    accepted = false
    for j = nextPivot .. frontSize-1:
        scan j, test 1x1, test 2x2 with its partner q
        if accepted: swap down, nextPivot += 1 or 2, break
    if not accepted: break

delay [nextPivot, frontSize)
```

A full pass with no acceptance is the same stopping condition `trials` expresses, and the bound on
each pass is `frontSize` rather than a counter. Any order is correct, since the acceptance test
bounds `L` whichever candidate passes it.

**What the queue provides is a scheduling policy.** A refused column goes to the rear, so it is not
offered again until every other candidate has had a turn. The index form restarts at `nextPivot`
after each acceptance, so a column refused a moment ago is among the first tried again. Neither is
obviously better, and the arguments run in opposite directions: a column that has just failed will
probably fail again after one rank-1 or rank-2 update, so retrying it at once tends to waste a scan;
but an acceptance changes the trailing block, and the columns nearest the frontier are the ones it
changes most, so the recently refused may be exactly the ones worth retrying. AGL do not settle it.
Both their figures rotate, and section 3.5 argues about which *pairs* to test rather than about
revisit order.

**A note on the initial order.** Figure 3.4 initializes its queue to `i1 < ... < iu` then
`iu+1 < ... < iv`, structurally eligible columns first and previously postponed ones at the rear,
and AGL report that skip effective in both their orderings. Oblio does not have it: delayed columns
are placed at the left of the parent's front by `assembleDelay`, with the parent's own columns
shifted right, and `pivotList` is built by walking the index set in front order. So the postponed
columns are offered first. The index form would inherit exactly the same order, so this is a
property of the storage layout rather than of the container.

**Cost.** `push_back` and `pop_front` are `O(1)`, but `remove(lq)` for the partner of an accepted
`2 x 2` is a linear scan, `O(m)`. Against an elimination that is `O(m h)` for the same pivot that is
a `1/h` share, and `h` is the large dimension, so the list's real cost is more likely its allocation
and pointer chasing than that scan. If it ever mattered, the cheap fix is not the textbook one: a
vector of list iterators would make removal `O(1)` at the price of another array to maintain, where
lazy deletion needs no new structure at all. `gblToLcl` already tracks each column's position and
`swap` keeps it current, so `gblToLcl[lj] < nextPivot` is exactly "already eliminated as someone's
partner", and a popped entry failing that test can simply be skipped, provided it is not counted
against `trials`.

So the honest summary is that the root kernel needs no queue because bounded Bunch-Kaufman cannot
refuse, while the non-root kernel needs to keep refused candidates in play because its refusals are
provisional, and the queue is one of at least two ways to do that. Which way is better is an
experiment, wanting a counter for scans or sweeps that the suite does not currently have, and the
comparison would only be fair with the cheap removal in place.


## Oblio

### Oblio notation

Two kernels, `factorDynamicRootSupernode` and `factorDynamicNonRootSupernode` in
`src/NumFactorEngine.cpp`, chosen by the caller on `parent[jj] == NIL`. A supernode runs one or the
other and never both.

The names, which follow the paper's letters wherever the paper has one:

```
jj                      the supernode being factored; doubled, as all supernode names are
nextPivot               the elimination frontier: how many front columns are done, and equally the
                        position the next accepted pivot is swapped down to. The figures have no
                        counterpart, working on an already permuted matrix
lj, lq, lr              global column indices, the paper's j, q_j and r
j, q, r                 their positions in jj's front, which the figures never need
jGamma, qGamma, rGamma  the corresponding column maxima, over the full block height
jFrontGamma             the same maximum for column j, stopped at the end of the front; q is where
                        it sits
jjVal                   the front's value block, column-major with leading dimension jjNumNodeIdx
at(r, c)                the position of row r of column c within it
jjPreFactorFrontSize    the number of fully assembled columns, non-root; the qualifier is there
                        because the epilogue contracts the front by whatever it delays
jjFrontSize             the same count at a root, where nothing is delayed and so nothing contracts
jjNumNodeIdx            the full block height, front columns plus update rows
threshold               mPivotThreshold, tunable, non-root only
alphaRoot               (1 + sqrt(17)) / 8, a compile-time constant, root only
pivotList, trials       the candidate queue and its sweep counter, non-root only
```

The root partner is `r` and the non-root partner is `q`, which is the paper's own distinction:
section 3.1 introduces `q` for the partner *restricted to* `A11`, a restriction that exists only
because the unrestricted one may lie in `A21`. The dense figures have no such restriction and write
`r` from Figure 2.1 onward, and a root has no `A21`. The candidate stays `j` where Figures 2.4 and
2.5 write `i`, since `i` is a row index by house convention and a pivot eliminates a column.

### The root loop

Bounded Bunch-Kaufman, Figure 2.4. Landed 2026-07-26, replacing a weakened copy of the non-root
pass: that one forced its last remaining candidate as a `1 x 1` whatever it looked like, and
accepted a `2 x 2` on the local-maximum condition with the diagonal conditions dropped, which is
what leaves the entries of the factor unbounded. The change moved the pinned pivot counts, as
intended, and brought a regression case with it.

The shape is Figure 2.4 rather than Figure 3.4, and the queue goes with it. A root cannot delay a
column, so `pivotList`, `trials`, the push to the rear and the outer sweep are all vestigial there:
nothing is ever postponed, and the chase reaches a pivot from any starting column. What is left is
take the next unfactored column, chase, accept, advance.

The partner is `r`, not `q`. `q` is the sparse letter: AGL introduce it in section 3.1 for the
partner *restricted to* `A11`, precisely because the unrestricted one may lie in `A21` and be
unavailable. The dense figures have no such restriction and call it `r` throughout, from Figure 2.1
onward. A root has no `A21`, so the restriction is vacuous and the letter should be the dense one;
writing `q` here would import a distinction the algorithm does not make.

`scanPivotColumn` serves unchanged. At a root `jjNumNodeIdx == jjPreFactorFrontSize`, so its third
loop is empty and `jFrontGamma == jGamma` identically. Its third output is the position attaining
the front-restricted maximum, which at a root is the maximum outright, so the root binds it to `r`.
Figure 2.4's `r` is ours; its candidate `i` is our `j`, since `i` is a row index by house
convention.

```
factorDynamicRootSupernode(nf, jj), bounded Bunch-Kaufman:

    nextPivot = 0

    while nextPivot < jjPreFactorFrontSize:
        j = nextPivot
        lj = jjNodeIdx[j]
        jGamma, _, r = scanPivotColumn(jj, nextPivot, j)

        if jGamma == 0:                                # isolated column, nothing to eliminate
            rank-- if |jjVal[at(j, j)]| == 0
            pivotType[lj] = 1
            nextPivot += 1
            continue

        if |jjVal[at(j, j)]| >= alphaRoot * jGamma:     # 1x1 on this column's own diagonal
            factor1x1(nf, jj, nextPivot, j, lj, ...)
            nextPivot += 1
            continue

        repeat:                                        # the chase
            rGamma, _, rNext = scanPivotColumn(jj, nextPivot, r)
            lr = jjNodeIdx[r]

            if |jjVal[at(r, r)]| >= alphaRoot * rGamma: # 1x1 on the column the chase reached
                factor1x1(nf, jj, nextPivot, r, lr, ...)
                nextPivot += 1
                break

            if jGamma == rGamma:                       # mutual maximum: accept the 2x2
                factor2x2(nf, jj, nextPivot, j, r, lj, lr, ...)
                nextPivot += 2
                break

            j = r;  lj = lr;  jGamma = rGamma;  r = rNext
```

What the change buys, in the order the guarantees fall out.

No delays, structurally. There is no branch that fails to produce a pivot, so a root front is always
factored completely, and it is not a probabilistic claim: `gamma` never decreases along the chase
and strictly increases except in the case that stops it, so no column is visited twice.

The `2 x 2` block is nonsingular by construction. The chase reaches its `2 x 2` only after both
diagonal tests have failed, `|A(j,j)| < alphaRoot * jGamma` from how it arrived and
`|A(r,r)| < alphaRoot * rGamma` from the branch above, while `|A(r,j)| == jGamma == rGamma`. So
`|det| > jGamma^2 (1 - alphaRoot^2)`, bounded away from zero. Those two failed tests are exactly the
guard the condition this replaced lacked, that one accepting a `2 x 2` on the equality of the two
maxima alone, and they are recovered here not as a patch but because the chase cannot reach the
acceptance without them.

The forced `1 x 1` disappears. The loop this replaced accepted its last remaining candidate
whatever it looks like, which is where the open question about dividing by zero comes from. Here
there is no last-candidate case, because acceptance does not depend on candidates remaining.

And the block is well conditioned, `kappa(D) < (1 + alphaRoot) / (1 - alphaRoot)`, which is the
condition under which the paper allows explicit inversion. Cramer's rule in `factor2x2` is
therefore sound at roots, though not at nonroots.

Two things this needs that the nonroot path does not. `alphaRoot` is the dense threshold,
`(1 + sqrt(17)) / 8`, not `mPivotThreshold`. The strong value is free at a root, since the front is
dense, a symmetric permutation of it stays dense, and no column can be delayed, so there is no fill
to trade against and the only cost is comparisons. But it does mean two thresholds in one
factorization, and the name above is a placeholder. And the scan is now called once per chase step
rather than once per candidate, so the worst case is `O(m)` traversals per pivot; the paper's
Appendix C bounds the expected number of column searches at about `e`, under 2.72, and Table 2.6
measures 2.44 on random matrices of order 400.

### The nonroot loop

Figure 3.4's search with Figure 3.3's test. The reshaping to this form landed 2026-07-26 and was
**behavior-preserving, exactly**: the same pivots, in the same order, from the same fronts, which
the pinned pivot and delay counts confirmed. What the acceptance test contains was a separate
question, settled separately afterwards, when the symmetric-maximum clause came out.

The names follow the figure. `lj` and `lq` are the candidate and its partner as global column
indices, which is what Figure 3.4 calls `j` and `q_j`; the subscript is unnecessary here because
there is one candidate in hand at a time and `lq` is by construction the partner of that one. `j`
and `q` are their positions in the front, which the figure has no need of and an implementation
cannot do without. `nextPivot` is the frontier, the count of front columns already factored and
equally the position the next accepted pivot is swapped down to. The figure has no counterpart for
it either, because it works on a matrix already permuted and leaves `P A P^T` to record silently
what `nextPivot` and `swap` record explicitly.

```
factorDynamicNonRootSupernode(nf, jj):

    pivotList = the front's global column indices, in front order
    nextPivot = 0

    while pivotList is not empty:
        pivotFound = false
        trials     = pivotList.size()

        while trials > 0:
            lj = pivotList.pop_front()
            j  = gblToLcl[lj]

            jGamma, jFrontGamma, q = scanPivotColumn(jj, nextPivot, j)
            jDiagonal              = jjVal[at(j, j)]

            if jGamma == 0:                            # isolated column, nothing to eliminate
                pivotFound = true
                swap(jj, nextPivot, j) if nextPivot != j
                rank-- if |jDiagonal| == 0
                pivotType[lj] = 1
                nextPivot += 1
                break

            else if |jDiagonal| > 0 and |jDiagonal| >= threshold * jGamma:
                pivotFound = true
                factor1x1(nf, jj, nextPivot, j, lj, ...)
                nextPivot += 1
                break

            else if q >= 0 and acceptPivot2x2(nf, jj, nextPivot, j, q, jGamma, jFrontGamma):
                pivotFound = true
                lq = jjNodeIdx[q]
                pivotList.remove(lq)
                factor2x2(nf, jj, nextPivot, j, q, lj, lq, ...)
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
scanPivotColumn(jj, nextPivot, j) -> jGamma, jFrontGamma, q:

    jFrontGamma = -1;  q = -1

    for i = nextPivot .. j-1:                              # lj's row, all of it front columns
        if jFrontGamma < |jjVal[at(j, i)]|:  q = i;  jFrontGamma = |jjVal[at(j, i)]|

    for i = j+1 .. jjPreFactorFrontSize-1:                 # its column, front part
        if jFrontGamma < |jjVal[at(i, j)]|:  q = i;  jFrontGamma = |jjVal[at(i, j)]|

    jGamma = jFrontGamma                                   # the front maximum is a prefix of the
                                                           #   full-height one, so continue from it
    for i = jjPreFactorFrontSize .. jjNumNodeIdx-1:        # its column, update rows
        if jGamma < |jjVal[at(i, j)]|:  jGamma = |jjVal[at(i, j)]|

    return jGamma, jFrontGamma, q
```

and the acceptance test, which is the hook Figure 3.4 leaves open:

```
acceptPivot2x2(nf, jj, nextPivot, j, q, jGamma, jFrontGamma) -> bool:

    qGamma = max of |jjVal[at(q, i)]|  over i = nextPivot .. q-1
                and |jjVal[at(i, q)]|  over i = q+1 .. jjNumNodeIdx-1

    d      = readPivotBlock2x2(jjVal, ld, j, q)
    maxmax = max( |d.d22| * jGamma + jFrontGamma * qGamma,
                  |d.d11| * qGamma + jFrontGamma * jGamma )

    return |d.det| > 0 and |d.det| >= threshold * maxmax
```

The test above is Figure 3.3's fourth branch and nothing else. Until 2026-07-26 it carried a second
clause alongside the growth bound, the equality of the candidate's two maxima with the partner's,
inherited from 0.9. That clause reached for the Bunch-Parlett condition of section 2.3, which needs
both diagonals bounded as well as the mutual maximum, and only the candidate's diagonal is known
small where the clause sat, so it admitted blocks whose factor entries are unbounded. Removing it
rather than adding the missing conjunct was the right repair, because with the conjunct restored the
clause accepts nothing the growth bound does not already accept, for any threshold at or below one
half.

On efficiency, this shape is ahead where it matters and level elsewhere. The fused scan replaces
two traversals of `lj`'s line with one whenever a `2 x 2` is considered, which is the expensive
path, and the three loops stay branch-free because the front part and the update part are separate
loops rather than one loop with a test in it. Where the `1 x 1` succeeds the fused scan does track
`q` and a second running maximum that the current code would not have computed, but over entries it
is loading anyway, so the cost is a comparison and a store and no extra memory traffic. Nothing is
computed earlier than it is needed except that: `qGamma` and the `2 x 2` block are still read only
inside the test, exactly as now, and `lq` is still looked up only on acceptance.

The scan also makes the relation between the two structural. `jFrontGamma` is computed first and
`jGamma` continues from it rather than starting over, so `jFrontGamma <= jGamma` follows from the
loop bounds rather than standing as an invariant to be argued across two separate scans.

Three things in the loop have no counterpart in Figure 3.4, and all three are implementation detail
rather than shape. The `q >= 0` conjunct covers the case the figure does not show, where the front
holds no partner to pair with. The isolated-column arm sets `pivotFound`, where the figure leaves it
alone. And the figure's nested `else { if acceptable ... else rear }` is flattened into a fourth arm
of the one chain.

`jFrontGamma` follows `nodeIdx` and `frontNodeIdx`, and `frontSize`: the same quantity as `jGamma`,
narrowed to the front. It is the paper's `|a_qj,j|`, and the two descriptions agree because `q` is
the argmax of the front-restricted scan, so the entry at `(q, j)` has that maximum as its magnitude
by definition.

One name is still open. `maxmax` is the quantity the determinant is tested against, and its name is
inherited from a pair of variables that no longer exist. The code's comment calls it the growth
bound, which is the Duff-Reid motivation; AGL's is that the two terms are the numerators of the
bounds on the entries of `L`, so following the paper the name would say `L` rather than growth.

The `2 x 2` block needs nothing, since `readPivotBlock2x2` already returns `d11`, `d22`, `d21`,
`d12` and `det`. That is `D`, the object the paper also writes as `D_{k,k+1}`, and it is the right
letter because it names what those four entries become.

### What the reshaping to Figure 3.4 changed

Kept as the record of why the non-root loop has the shape it has. The loop described below is the
one that stood before 2026-07-26, and the five differences named are what the reshaping removed;
none of them changed a single pivot decision, which is why the pinned counts were the check.

The two agreed on everything that decides an outcome: the queue, the order of the tests, `1 x 1`
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
