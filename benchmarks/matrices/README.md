# Matrices Benchmark

**`ACCURACY.md` and `PERFORMANCE.md` beside this file are the reports.** They are written for
someone evaluating Oblio: how good the answers are, and what they cost. This README is the working
notes behind both, including what went wrong on the way, the two measurement mistakes that had to
be corrected before the accuracy results could be trusted, and the two ordering investigations the
performance run set off.

The pipeline on matrices nobody generated. Every figure in this tree until now was measured on a
grid Laplacian or on one of seven hand-built examples, so every claim drawn from them is a claim
about those. This folder is where real ones get in, from the SuiteSparse Matrix Collection.

**A benchmark, not an experiment**, on the same terms as `../ordering` and `../pipeline`: it links
`../../src` directly and is expected to keep compiling as the tree moves, which is why `make`
builds it and only the run targets measure.

**And it is named for its input rather than for what it measures**, unlike its two siblings, which
are named for what they compare: one phase against itself, and the phases against each other. What
distinguishes this folder is the matrices.

```
make             build the three programs
make accuracy    the accuracy pass over accuracy_candidates.txt
make performance the timing pass over performance_candidates.txt
make read        the reader alone over data/, an inventory of what is usable
make clean
make help        print this list, including the fetch commands
```

**One driver per report, both reading the same files.** `matrix_accuracy_cpp` answers how good the
answers are; `matrix_performance_cpp` answers what they cost, across four orderings and three
traversals with Cholesky held still. `mmread_cpp` inventories the pool and factors nothing.

## Getting the matrices

They are not in the repository. They are large, they are somebody else's, and they are
reproducible from one script, which is what is committed instead.

```
./ssget.py list --per-kind 8 --values --max-nnz 500000 > accuracy_candidates.txt
./ssget.py fetch accuracy_candidates.txt
make accuracy
```

`list` downloads only the collection's index, one line per matrix, and prints what matches the
filter. It is cheap to re-run until the candidate set looks right, and nothing moves until
`fetch`. The output of `list` is valid input to `fetch` unedited, so the loop is print, delete the
rows we do not want, fetch. Matrices land in `../../data/<Group>/<name>.mtx`, the collection's own
layout, so a path says where a file came from and two groups may share a name without colliding.
`data/` is gitignored in full.

### The candidates files

**A candidates file is the INPUT to a run, not its output.** It records which matrices a report was
built from, which is what makes that report reproducible on another machine, and it is regenerated
rather than edited. Two are committed, one per report, and the commands that produced them are:

```
./ssget.py list --per-kind 8 --values --max-nnz 500000 > accuracy_candidates.txt
./ssget.py list --posdef --max-nnz 2000000            > performance_candidates.txt
```

A difference in one of these files is meaningful only if the command that produced it changed, so
the commands belong here beside them.

**They select differently because the reports ask different questions.** The accuracy set samples
by problem kind, since structural variety is what generated grids lack. The performance set takes
**every** positive definite matrix in range without sampling, `--posdef` narrowing the pool to 89
by itself, and sampling a pool that small would be worse than taking all of it.

**And `data/` holds the union**, plus whatever else has been fetched by hand, which is the point of
keeping the pool separate from the run: neither list has to be curated when the other changes, and
a matrix fetched for one report costs the other nothing.

### A list can name anything the collection has

`fetch` looks names up in the **full index**, not in the filtered one, so a list may name a matrix
no filter here would select. `extras.txt` is committed as the third such list, beside the two
candidates files, and holds three of them:

```
./ssget.py fetch extras.txt
```

```
PARSEC/Si87H76
PARSEC/Ga41As41H72
Schenk/nlpkkt80
```

**It is a hand-written list rather than generated**, which is the difference from the two
candidates files: no command produced it, so nothing regenerates it and editing it is the way it
changes. Add to it, or write another beside it, whenever a set is wanted that no filter describes.

That `fetch` accepts such a list is deliberate and it was not always so. **The filter's job is to
PRODUCE candidates; once a list exists, the list IS the selection.** An earlier version looked
names up in the filtered set and answered `not in the filtered index`, which made a recommended
set from somebody else unusable without first widening a filter to admit matrices nobody wanted to
select on. Only a name the collection does not have is refused now.

**Nothing checks whether the drivers can use what arrives, and nothing needs to.** Each driver
screens what it reads: `matrix_accuracy_cpp` by field type, by structural rank and by predicted
fill, with a fork guard behind that; `matrix_performance_cpp` by definiteness and by predicted
fill. And neither reads `data/` directly, both taking their own candidates file, so a matrix
fetched by hand appears in a run only if a list names it.

The three in `extras.txt` are the case that prompted this. They are standard hard matrices from
the sparse direct literature and all three are far beyond what either report covers: `Si87H76` is
n = 240369 with 10.7 million nonzeros, and Oblio's analysis predicts **5.68 billion entries of fill
under MMD3**, which is 42.3 GB of values. The accuracy driver reported that and declined, which is
the cap working: the analysis completed, since the symbolic factor stores index sets rather than
entries, and only the numeric phase would have allocated. `PARSEC` matrices fill this way under any
minimum degree ordering; what they want is nested dissection, which Oblio does not have.

**The filter is narrow and matches what the reader accepts**: real, square, and numerically
symmetric. A matrix stored `general` whose values happen to be symmetric is excluded, and so is
anything complex. Widening any of that is a change in two places at once.

**Sampling by kind is the point of `--per-kind`.** The `kind` column is the only account of
structure the index gives us, and structural variety is the whole reason for coming here: a
hundred matrices that are all `structural problem` would be one family again in a new costume. It
takes the cheapest N of each kind, which biases small within every kind. That is what a
correctness pass wants and is not what a timing ladder will want later.

**`--max-nnz` is what bounds the download**, where `n` bounds nothing useful. The unfiltered
candidate pool at 1000 <= n <= 100000 was 719 matrices and 718 million nonzeros, which as Matrix
Market text is on the order of ten gigabytes; `--per-kind 2 --max-nnz 500000` took that to 60
matrices and 4.7 million nonzeros. It does NOT bound the run: see the fill cap section above.

**`--values` excludes binary matrices, and it works because binary means pattern here.** A binary
matrix has values that are all 1, so the collection stores its pattern and nothing else, and the
file comes back with Matrix Market field type `pattern`, two indices per line and no third column.
The index's `isReal` flag does not predict this, meaning only "not complex". The index's
`isBinary` flag predicts it exactly, on two independent samples: **19 of 19 pattern files matched
the binary rows and no others, then 41 of 41.** So `--values` makes `--per-kind` count matrices the
accuracy run can use, and `--per-kind 6 --values` produced 99 of them where reaching the same
number without it would have meant downloading around 160 and discarding 60. Leave it off when a
pattern is all a study needs, which is true of any ordering question.

**It narrows the structural variety, though, and that is the cost.** The kind count fell from 35
at `--per-kind 2` to 28 at `--per-kind 8 --values`, because the graph kinds vanish entirely: a
graph is stored as a pattern. So the accuracy set is structurally narrower than the pool it came
from, and the pattern files hold variety the valued ones cannot.

**`fetch` reports `have` against `got`.** Widening a filter re-runs the whole selection, so the
summary says how much was already on disk: `99 of 99 in .../data: 66 already there, 33
downloaded`. Growing the set never looks like starting over, and the fetch is idempotent.

### What the collection's index does and does not mark

Two properties were checked against the index directly, and one of the answers is a trap worth
recording.

**`posdef` is marked and usable on the real side.** With our filter plus `posdef == 1`, the pool in
the 1000 to 100000 range under 500000 nonzeros is 89 matrices, which is the whole performance set.

**Numerical symmetry is computed WITHOUT the conjugate**, so it measures `A == A.'` rather than
`A == A'`. For a complex Hermitian matrix with nonzero imaginary parts every off-diagonal
disagrees, and the score collapses: `Bai/qc324` and `Bai/qc2534` are Hermitian and the index scores
them at **exactly 0** while their pattern symmetry is 1. So a filter on `nsym == 1` selects complex
**symmetric** matrices and silently excludes every complex Hermitian one, which would look like the
collection having none rather than the filter hiding them.

**Nothing in the index distinguishes Hermitian from merely pattern-symmetric.** Both give
`psym == 1` with `nsym` below 1. Only the Matrix Market banner says `hermitian`, which is why
`fetch` printing every file's banner is worth its line of code: the index narrows and the file
decides.

**And complex positive definiteness is effectively unmarked.** Of the 23 complex square matrices
with `psym == 1`, one is flagged `posdef`, and that one, `Bai/mhd1280b`, has `nsym == 1` and so is
complex symmetric rather than Hermitian. **There is no complex Cholesky set to be had from real
matrices**, which is why complex Cholesky belongs in a scaling report on synthetic Hermitian grids,
where definiteness comes free from the construction.

### The pool and the run are separate selections

Worth stating because it is easy to collapse the two and then fight the filter.

**`ssget.py` builds a pool.** It selects on what the collection's index knows: size, nnz, symmetry,
binary, kind. A pool can also come from a curated list somebody else recommends, which has no
filter behind it at all.

**The driver selects over that pool**, on criteria the index cannot supply, and today that is the
shell glob plus `--max-fill`. The two need not agree and should not be made to: `data/` holding
the 41 pattern files from an earlier pass costs the accuracy run nothing, they are skipped by
field type, and they are exactly what an ordering study will want later.

## The size columns, and the fill cap that has to exist

```
                                      nnz(L) ------------------
matrix                 n    nnz(A)  predicted     actual   ratio
Marini/eurqsa.mtx   7245     48263     156243    3101560   19.9x
Cote/vibrobox.mtx  12328    342828    2332570    7113191    3.0x
Test/spd_grid.mtx    400      1920       3648       3648    1.0x
```

**`predicted` is what the analysis said and `actual` is what the factorization held, and the ratio
between them is the price of the pivoting.** They differ ONLY under dynamic pivoting: a delayed
column widens its parent's front, so the gap is exactly what the delays cost, and the `actual`
column is read from the dynamic run for that reason.

A statically pivoted factor reads 1.0x, which is every row with no `Nd` in its note, so the column
is self-checking. `eurqsa` is the case worth seeing: 40281 delayed columns take the factor from
156243 entries to 3101560, **twenty times the prediction**. That is the price of the pivoting that
gave it a backward error of 6.6e-19 where static LDL, refusing to move, got 3.8e-02.

### The cap, and why nothing before the analysis could have set it

`GHS_indef/bloweya`, n = 30004, killed a run outright: SIGKILL from the operating system's
out-of-memory killer, not a fault in the code.

**`ssget.py --max-nnz` cannot prevent that, because it bounds `nnz(A)` and the cost is `nnz(L)`.**
The two differ by orders of magnitude and nothing in the collection's index predicts the second:
`PARSEC/benzene` reaches 13.9 million entries of fill from 8219 rows. There is exactly one moment
when the number is known and nothing has been allocated on it yet, which is after `analyze` and
before any factorization, and that is where the check goes:

```
./matrix_accuracy_cpp --max-fill=2e8 ../../data/*/*.mtx
```

The default is 50 million, generous enough for everything fetched so far. A row over it prints its
predicted fill and an estimate of what the values alone would take, then moves on like any other
skip, and the accounting at the end counts it separately from the refusals.

**Eight bytes per entry is a floor rather than an estimate.** It counts the value array and nothing
else: no row indices, no update matrices, and no allowance for a delayed column widening a front.

**And the cap is checked against the PREDICTION, which is the one thing it cannot help with.**
`Oberwolfach/LFAT5000` predicts 67463 entries and factors into 15678721, **232 times** what the cap
saw. `eurqsa` is 19.9x and `bloweya` ran off the scale entirely. So the cap bounds what the
analysis says rather than what the factorization will hold, and on a row like that it cannot see
the cost coming. Nothing before the factorization can, which is why the fork exists as well.

## What delays cost, and the two matrices that show it

`GHS_indef/bloweya`, a Cahn-Hilliard problem at n = 30004, killed two runs outright with SIGKILL
from the operating system's out-of-memory killer. It is worth writing down in full because the
cause is not a defect and is not the benchmark's business either.

**The analysis is innocent.** 30002 supernodes, largest front **3**, predicted fill 110007
entries. Under a megabyte of values. So the fill cap could not have fired and would not have
helped.

**Static LDL factors it exactly as predicted**, 110007 entries. **Dynamic LDL exceeds 3.6 GB** and
no limit was found that it fits in.

**The shape is a chain, and 20003 of its 30004 diagonal entries are structurally zero:**

```
                n      supernodes   trees   height
bloweya      30004         30002        1     5004
bloweybq     10001          9997        1     4999
eurqsa        7245          6425        3       37
sherman1      1000           854      318       22
```

So almost every pivot candidate is one that dynamic pivoting must delay or pair, and **delays
cascade up a chain**: a delayed column widens its parent, the wider parent has more columns it
cannot use, and over a depth of 5004 the fronts grow without bound. The same effect is visible
and bounded elsewhere in the table, `eurqsa` going from 156243 predicted entries to 3101560 actual
and `vibrobox` from 2332570 to 7113191. This is that effect with nothing to stop it.

**`bloweybq` is the control that makes this an argument rather than a guess.** Same author, same
problem family, the same chain shape at height 4999, but positive definite: zero delays, and it
factors in 39997 entries. The chain alone is harmless. The chain plus a diagonal that forces
delays is not.

**Nothing here is wrong.** The factorization is doing what it is defined to do, and the question
it raises is whether delay growth should be bounded on a deep tree, which is a question for the
strategy rather than for this folder. Ashcraft, Grimes and Lewis is where Oblio's pivoting comes
from and bounding this is part of what that paper is about, so the finding is recorded here and
belongs upstream of it.

### The same mechanism, survived: LFAT5000

`Oberwolfach/LFAT5000`, n = 19994, is `bloweya`'s cascade caught in the act rather than off the
end of the scale:

```
matrix                     n   nnz(A)  predicted     actual   ratio   delayed
Oberwolfach/LFAT5000   19994    79966      67463   15678721  232.4x   3128750
```

**232 times the predicted fill**, from 3.1 million delayed columns. It finishes because the
cascade lands at 15.7 million entries rather than running out of memory, and every factorization
reaches a backward error near 1e-19, so the answer is right. The cost is the whole of the finding.

**And it is a sharper case than `bloweya`, because the matrix is positive definite.** 19994
positive eigenvalues, zero negative, Cholesky succeeds. So **every pivot was acceptable and 3.1
million columns were delayed anyway.** On `bloweya` the delays are forced, 20003 of 30004 diagonal
entries being structurally zero. Here nothing forces them. Something in the pivot search is
delaying when it does not have to.

**The two controls now sit on both sides of it:**

```
                     shape          definite   delays        ratio
bloweybq             chain, h 4999  yes        0             1.0x
LFAT5000             chain-like     yes        3128750     232.4x
bloweya              chain, h 5004  no         off the end   -
```

`bloweybq` says a deep chain alone is harmless. `LFAT5000` says definiteness alone is not enough
to prevent the cascade. So neither the tree shape nor the pivot quality explains it on its own,
and the question is what the pivot search does on a chain that it does not do on a bushy tree.
That is a question for the strategy, upstream of this folder, and Ashcraft, Grimes and Lewis is
where Oblio's pivoting comes from.

**One consequence for the fill cap.** It is checked against the PREDICTION, and this factor was
232 times it. So `--max-fill` bounds what the analysis says, not what the factorization will hold,
and on a row like this one it cannot see the cost coming. Nothing before the factorization can.

### Why a try/catch does not save the run, and what does

**SIGKILL cannot be caught.** The kernel destroys the process without unwinding, so no handler
runs and no destructor fires. `std::bad_alloc` is catchable and the driver already catches it, but
it is thrown only when an allocation FAILS, and under overcommit the allocation succeeds and the
kill arrives later when the pages are touched.

Under a `ulimit -v` on Linux the failure does arrive as `std::bad_alloc`, and the driver's existing
catch turns it into a clean skip. macOS generally does not enforce `RLIMIT_AS`, which is why
alpamayo got SIGKILL where the sandbox got an exception. So the portable answer is **one child
process per matrix**: the parent forks, the child does the row, and a death by signal becomes a
labeled row rather than the end of the run. `setrlimit` inside the child is worth having beside
it, since where the kernel honors it the message is much better than "killed".

**A skip list would also work and is worse.** It records a decision about one matrix on one day,
goes stale silently, and the next such matrix, arriving with somebody's recommended set, takes the
run down again. `fork` and `waitpid` never need maintaining and catch the next case automatically.

**So the driver runs one child process per matrix.** The parent forks, waits, and reads the
child's findings off a pipe as a single line, since a child cannot append to the parent's vectors.
The human-readable row goes to stdout, which the child inherits, so it lands exactly where it
would have without the fork. A child that exits normally reported everything; a child killed by a
signal reported nothing, and the parent finishes the half-written row:

```
GHS_indef/bloweya.mtx     30004    170012     110007  KILLED by signal 9, out of memory
Test/spd_grid.mtx           400      1920       3648       3648    1.0x  ...
```

The run continues, and the accounting counts the killed rows separately. The size columns flush
before any factorization starts, which is what makes the half-written row worth finishing: it
already carries `n`, `nnz(A)` and the predicted fill.

**`--max-memory-gb=G` is the better path where the kernel honors it, and macOS does not.** Measured
rather than assumed: `--max-memory-gb=8` on alpamayo still produced `KILLED by signal 9`, so
`RLIMIT_AS` is not enforced there and the fork is what holds. The option stays for Linux, where it
sets the limit in the child, the allocation fails, `std::bad_alloc` is thrown, and the existing
catch turns it into an ordinary cell note. That is a far better message than a signal, and it also
isolates the cause:

```
GHS_indef/bloweya.mtx  30004  170012  110007  -  -  ... 1.0e-18 1.0e-06 10001p  ... - std::bad_alloc
```

Static LDL completes with 10001 perturbed pivots and only the dynamic factorization runs out,
which is the finding above stated by the driver rather than by hand. It is **off by default**, a
guard that is silently not applied being worse than no guard at all.

**The cost of the fork is a process per matrix and nothing else measurable.** The local set runs in
the same time it did before.

## The two measures, and why the residual alone is not one

Each cell of the table holds two numbers off the same residual vector, in infinity norms
throughout:

```
bwd = ||Ax - b|| / (||A|| ||x|| + ||b||)        the backward error
res = ||Ax - b|| / ||b||                        the relative residual
```

**The backward error is the verdict.** It says we solved a nearby problem, `(A + E)x = b` with
`||E||` small against `||A||`, and for a stable factorization it sits at machine precision
whatever the conditioning. **The relative residual is not a verdict**, because it differs from
the backward error by a factor of `(||A|| ||x|| + ||b||) / ||b||`, which reaches the condition
number when `x` is large relative to `b`. Real matrices span a range of conditioning that no grid
in this tree ever did, so on this table a large `res` beside a tiny `bwd` is a statement about the
matrix and not about us.

**The note beside each pair carries what a number cannot:**

```
Np    pivots the static LDL had to perturb, having no way to pivot
Nd    columns the dynamic LDL delayed to a parent
```

### How to read the pair

**Read `bwd` first, always.** It is the verdict, and it is the only column that means the same
thing on every row whatever the matrix. The residual is not a second opinion on it: it answers a
different question, and the pair is read as two axes rather than as a check and a cross-check.

- **`bwd` small, `res` small.** Nothing to see: a well conditioned matrix, solved.
- **`bwd` small, `res` large.** We solved it correctly and either the matrix is ill conditioned or
  the system has no solution. `eurqsa` and `bloweybq` below are the two shapes. The ratio between
  the columns is NOT a clean measure of conditioning, since it also grows with the scaling of `A`;
  the summary section says why that matters.
- **`bwd` large.** The only case that is about us rather than about the matrix, and it needs
  three exclusions checked before it is one. `Np` in the note says static LDL perturbed, which is
  the method doing what it is documented to do, so a large `bwd` there is expected. A nonzero
  `zero` in the inertia columns says the matrix is numerically singular, so the system is probably
  not solvable at all. `Nd` says columns were delayed. **A large `bwd` with none of the three
  would be the genuine alarm**, and there is no such row today.
- **`bwd` large, `res` small** is close to impossible, `res` being the larger of the two whenever
  `||A|| ||x||` exceeds `||b||`.

**And `bwd` is not comparable BETWEEN methods when one of them perturbs**, which is the least
obvious thing on this page. On `saylr3` static LDL reads 1.9e-15 against dynamic's 8.9e-05, and
the static figure is the less useful of the two: it got there by perturbing the two zero pivots
and returning an `x` of 2e14, and a huge `x` with a bounded residual satisfies the backward-error
definition trivially. Dynamic returned a bounded `x` and an honest failure to solve a system that
has no solution. Compare a column against machine precision, not against the column beside it.

So the residual never rescues a bad `bwd` and never condemns a good one. What it adds is the
second axis: `bwd` says whether we did our job, and the gap between the two says how hard the
matrix was.

**And `bwd` gets no threshold, deliberately.** A static LDL that perturbs is EXPECTED to land far
from machine precision, so one constant cannot be checked against three factorizations that carry
different contracts. The note column is what says which contract a row is under.

### The second case, and the free control beside it

`HB/saylr3` and `HB/sherman1` are the same 10x10x10 grid problem from two different authors, and
they are identical in n, in nnz, in nnz(L) at 9872, in diagonal range and in `||A||`. One is clean
and one is not:

```
                     bwd        res     note      pos   neg  zero
saylr3    static   1.9e-15    2.0e+00    2p
          dynamic  8.9e-05    2.0e+00               0   998     2
sherman1  both     8.1e-17    1.8e-12               0  1000     0
```

**`saylr3` is numerically singular, rank 998 of 1000**, and three independent witnesses say so:
the inertia reports two zero eigenvalues, static LDL perturbs exactly two pivots, and a consistent
right-hand side removes the problem completely. Handed `b = A x`, dynamic LDL on `saylr3` returns
**9.31e-17, identical to `sherman1`**. So `b` all ones is simply not in the range, and the 2.0
residual in both columns is that and nothing else.

**This is the case the structural label cannot catch by construction.** Structural rank is a full
1000, the matching passes it, and dynamic LDL then neither perturbs nor delays because accepting a
zero pivot and losing rank is exactly what it should do. Any alarm criterion built on "no
perturbations and no delays" fires here, which is why the inertia columns are printed on every
row.

**And `sherman1` is negative definite**, 0 positive against 1000 negative, which is why Cholesky
refuses it. `-A` would factor by Cholesky and `A` will not.

**And the pair is a free control worth keeping.** Same pattern, same symbolic factor, different
values: whatever separates them is in the numerics rather than in the ordering, the forest or the
symbolic factorization, and no work was needed to isolate that.

### The case that forced the change

`Marini/eurqsa` is the case that forced the change and it is worth keeping as the illustration.
It is a KKT system from economic time series reconciliation, structurally nonsingular, with 2121
zero diagonal entries, a diagonal of 1 to 2 where present, and off-diagonals reaching 1.19e6:

```
                       bwd        res      max|x|   perturbed
static LDL          3.8e-02    6.9e+78     6.7e+94      2079
dynamic LDL         6.6e-19    2.3e+15     1.7e+17         0
```

**Dynamic LDL is backward stable here at 1e-18**, and the enormous residual is the amplification
factor and nothing else. On alpamayo the static column read `nan` rather than a number, which is
the same event seen through an overflow: `x` reaches where squaring it gives infinity, and
`inf - inf` is a nan. Static LDL's 3.8e-02 is meanwhile the honest figure and it is genuinely
poor, but it is documented behavior rather than a defect: 2079 perturbed pivots means it factored
a matrix differing from `A` by a few percent, which is what a factorization that cannot pivot has
to do on this input.

**And the matrix is numerically near-singular while structurally full**, which is exactly the case
the structural label cannot catch. Handed a consistent right-hand side with `||xExact|| = 229`,
the computed solution differs from it by 5.6e17 while the residual stays at 7.7e6, so the error is
almost all a near-null-space direction and the smallest singular value is around 1.4e-11 against a
norm of 6.4e6. Symmetric equilibration takes `||A||` from 6.4e6 to 11 and improves matters by
orders of magnitude without fixing them.

## The inertia, which is of A and not of a factorization

The last three columns are the eigenvalues of `A` counted by sign, and they close the table
rather than sitting beside `n` for a mechanical reason: they are read from a factorization, so
they are not known until one has run.

**Read from `D`, not from `L`.** `L` is unit lower triangular and its eigenvalues are all 1, so it
carries nothing. `A = L D L^H` is a congruence, a congruence preserves the signs of the
eigenvalues, which is Sylvester's law of inertia, and so counting them in `D` counts them in `A`
without forming an eigenvalue.

**One set of columns and not three, because it is a property of A**, so all three factorizations
should agree. **They do not always, and the exception is the reason the printed one is taken from
the DYNAMIC factorization**: a static LDL that perturbed reports the inertia of the matrix it
actually factored, which is a nearby one. On `eurqsa` static reads 83 zero eigenvalues where
dynamic reads 8, and on `saylr3` static reads none where dynamic reads 2. Dynamic does not
perturb, so its answer is about `A`.

**A nonzero `zero` gets the word SINGULAR at the end of the row**, the same word the
structurally singular rows carry, so the two read alike at a glance:

```
HB/saylr3.mtx      1000    9872  ...  8.9e-05  2.0e+00        0    998      2  SINGULAR
Newman/netscience.mtx    1589    7073  SINGULAR, structural rank 1424 of 1589
```

**They are not established alike, and the difference is the whole reason both exist.** The
structural one is proved from the pattern before anything runs, and that row is declined. This one
is found by the factorization, so the row is solved and reported: it is a property of the outcome
where structural rank is a property of the input, and only the input can be gated on.

It is also the one class the structural label cannot catch: full structural rank, no
perfect-matching failure, and singular anyway. And it is a marker rather than a proof, for a
reason stated at the declaration in `DirectSolver.h`: for a singular `A` a zero eigenvalue lands
on whichever side rounding puts it, so the split is least reliable exactly where it is most
interesting.

**And the columns cross-check the Cholesky column for free.** Cholesky succeeds exactly when `A`
is positive definite, so it should refuse on every row whose `pos` is below `n` and on no other.
`sherman1` is the case worth seeing: 0 positive and 1000 negative, so it is negative definite,
`-A` would factor by Cholesky and `A` will not, and the refusal is correct rather than a
limitation.

## The summary at the end

Three blocks close the run, and the split between what they can prove and what they can only
point at is the whole design.

**The first block accounts for every file**, so nothing disappears between the glob and the
classification and the counts can be checked against each other:

```
of the 160 files read:
  pattern, no values      41   nothing to solve, and nothing lost: an
                               ordering study would want exactly these
  with values            119
    structurally singular  12   declined before factoring, listed above
    killed                  1   the child did not survive the matrix, so
                               the row above says how far it got
    solved                106
```

**The classification is exact and has no threshold in it.** The inertia gives it outright:

```
of the 106 solved, by the inertia of A:
  positive definite       30   every eigenvalue strictly positive
  negative definite       18   every eigenvalue strictly negative; Cholesky
                               must refuse these, and -A would factor
  indefinite              56   both signs present, none zero
  numerically singular     2   at least one exactly zero, whatever the other
                               signs. Counted ONLY here, not in the three
                               rows above, the zero test coming first
```

**The four classes are disjoint and the zero test comes first**, since a matrix carrying a zero
eigenvalue is singular whatever the other two counts say. So the definiteness rows see only the
matrices with no zero at all, and a singular matrix that is otherwise indefinite is not counted
among the indefinite ones. The four therefore sum to the solved count.

**The Cholesky check is free and it is a real assertion.** Cholesky succeeds exactly on a positive
definite matrix, so its column and the inertia columns are two accounts of the same fact and must
agree on every row. The run says so or names the rows where they do not. Nothing else in this
folder cross-checks one part of the library against another.

**The third block ranks rather than selects.** It lists the largest residuals over the whole set,
taking the best of the three factorizations, since a row where one refuses and another solves
cleanly is a solved row. A ranking needs no cutoff invented for it, and the rows at the top are
the rows where we could not produce a usable answer whatever the rest of the table looks like.

```
the 8 largest residuals, taking the best of the three factorizations on
EACH MEASURE SEPARATELY, which may be two different factorizations:
matrix                               best res   best bwd  perturbed    zero
Marini/eurqsa.mtx                     2.0e+16    5.8e-18       2084       8
Guettel/TEM27623.mtx                  1.1e+04    3.5e-18         79       0
GHS_indef/sit100.mtx                  2.6e+03    9.6e-15       1063       0
Oberwolfach/t2dal_bci.mtx             1.4e+03    1.1e-16          0       0
HB/plat1919.mtx                       1.5e+02    7.3e-17          1       0
```

`best res` says whether ANY factorization got the residual down and `best bwd` whether ANY was
backward stable. A large residual beside a `bwd` at machine precision is the matrix and not us: no
`x` reproduces `b`, either because none exists or because `||x||` is so large that forming `Ax`
loses `b` entirely, and `zero` says which of those is proved. **A large residual beside a large
`bwd` and no perturbations would be ours, and there is no such row today.**

**The two columns may come from different factorizations, and that is the point.** Taking them as
a pair went wrong twice, in opposite directions, which is what settled it:

- **Choosing by residual**, on `GHS_indef/sit100`, static LDL's 2.6e+03 edged dynamic's 3.2e+03 and
  dragged its backward error of 1.3e-02 along, against dynamic's 9.6e-15.
- **Choosing by backward error**, on `HB/plat1919`, Cholesky's 7.3e-17 edged dynamic's 9.9e-17 and
  dragged its residual of 1.0e+03 along, against dynamic's 1.5e+02.

Both times a rounding-level difference in one column decided the other one. The two measures
answer different questions and neither answer is improved by insisting one factorization supply
both, so decoupling them removes the tie-break entirely.

### The third case, which needs no marker at all

`GHS_indef/bloweybq` is the row that ranks third and carries nothing: `zero` is 0, no
perturbations, full structural rank, and `bwd` at 2.1e-19.

```
b = ones    dynamic   bwd 2.1e-19   res 8.8e-01   max|x| 8.3e+14
b = A x     dynamic   bwd 1.8e-16   res 9.1e-16   max|x| 4.0e+00
```

**It is not singular, it is conditioned at about 1e15.** Given a consistent right-hand side all
three factorizations solve it and recover `x_exact` to 1.6e-3 relative. Given `b` all ones the
solution is legitimately of size 1e15, and at that magnitude `Ax` cannot be formed accurately
enough to reproduce `b`, so the residual near 1 is a cancellation artifact rather than a failure to
solve. The estimate off the consistent run puts the smallest singular value near 1e-11 against a
norm of 5e3.

Two things about it are worth not misreading. Its group name is `GHS_indef` and the matrix is
positive definite, all 10001 eigenvalues, which is why Cholesky succeeds on it: the group is a
collection for indefinite solvers rather than a claim about each member. And its own title in the
file is "matrix for which early version of MA57 fails", so it is a known-hard case that we handle
correctly rather than a case we handle badly.

**So the two invisible cases are different and only one of them wants a marker.** `saylr3` is
invisible singularity, which the `zero` column catches. This is invisible conditioning, which
nothing marks and nothing should: the answer is correct and the question was ill-posed.

## The four programs

**`mmread_cpp`** reads and reports, and does not factor anything. It is the inventory: what each
file is, what it converted to, and what was refused. Run it first on a new set.

**`matrix_ordering_cpp`** answers what the ORDERING STEP ALONE costs, genmmd against `Mmd3`, and
is the newest of the four, 2026-08-15. `make ordering`, over all of `data/` rather than a
candidates file: an ordering needs neither values nor definiteness, so what it wants is everything
on disk, which is 246 files against the performance report's 107.

It is a different question from the ordering column in `matrix_performance_cpp`'s table, which
prices ordering AGAINST the rest of a solve. This one prices two implementations of the same
ordering against each other, and the comparison is unusually clean: `make mmdmatrices` in
`experiments/ordering` shows the two return the SAME permutation on every one of these matrices, so
there is no fill to trade against time. One `nnz(L)` column serves both, and time is the only thing
that differs. Both routines go through the same timing helper with the same warm-up and the same
repeat count, which is not fussiness: a column timed by one path against a column timed by another
differed 2.4 percent on 2026-08-10 and the difference was read as a result.

**`matrix_accuracy_cpp`** answers how good the answers are. One ordering, MMD3, and the
left-looking traversal, both held fixed on purpose: the question is whether we compute correctly on
matrices nobody generated, not which ordering is best. It drops `nnz(A)` from its table, which
`mmread_cpp` already reports, to make room for the two measures. `ACCURACY.md` is its report.

**`matrix_performance_cpp`** answers what they cost. Four orderings, the vendored MMD and AMD
beside our MMD3 and AMD3, and three traversals, with **Cholesky held still**. That choice buys
three things at once: the factor's structure is exactly what the analysis predicted, so nnz(L) is a
property of the ordering alone; no column is delayed, so the fill cap is exact rather than a lower
bound and no fork guard is needed; and the numeric phase does the same arithmetic under every
traversal, so a difference between them is a difference in scheduling rather than in work. Its
phase split and timing protocol are `../pipeline`'s, so rows can be read across the two folders,
and the vendored pair is a live oracle in every row: MMD and MMD3 must agree on nnz(L), as must AMD
and AMD3.

All three take file names, so a subset needs no option:

```
./matrix_accuracy_cpp ../../data/HB/*.mtx
./matrix_accuracy_cpp --max-fill=2e8 ../../data/*/*.mtx
./matrix_performance_cpp --repeats=1 ../../data/HB/*.mtx
```

`matrix_performance_cpp` needs `../../private` for its MMD and AMD rows; without it those two
refuse and the other two still run. It is the only thing in this folder that does.

## What the reader accepts, and what that turned out to matter for

```
object      matrix          the format also covers dense arrays under `array`
format      coordinate      the sparse spelling
field       real, integer   both become double
            pattern         ON REQUEST ONLY, and only `matrix_ordering_cpp` asks. A pattern file
                            carries structure and no values, two indices per line, which is useless
                            for a residual and exactly an ordering's input. The default refuses
                            them so a caller who needs numbers cannot get a matrix of ones by
                            accident. Added 2026-08-15; it is most of why the ordering set is twice
                            the size of the other two.
symmetry    symmetric       and the file then stores ONE TRIANGLE, which we mirror
```

**`integer` is a spelling, not a meaning.** It says the file writes `3` rather than `3.0`, and
several matrices here are legitimately in it: Gset's weights are +1 and -1, and Trefethen's
entries are built from primes. Refusing it would have cost five matrices for no reason.

**`pattern` files carry no values at all**, two indices per line, so they are useless for a
residual. They are exactly what a future ordering-only check would want, which is why they are
fetched and kept rather than filtered out.

**The reader mirrors, and it checks the triangle.** The format says a `symmetric` file holds the
lower triangle and the collection's files do, but nothing in a file enforces it, and an upper
entry in a file we then mirror would produce a duplicate that the conversion silently SUMS. So a
malformed file would become a wrong number rather than an error. The check is one comparison per
entry and it is the only validation anywhere in this tree today.

**The conversion is `examples/example_matrix.cpp`'s `fromTriplets`**, whose own comment calls it
the conversion worth stealing and states the three preconditions a `SparseMatrix` has: full
storage, a structurally present diagonal, and sorted rows with no duplicates. It inserts a
structural zero on any diagonal nothing landed on, which is what the symbolic factorization needs
and is ordinary input for LDL. So the requirement that a matrix arrive with a diagonal is not one
we have to impose. It is the readable form rather than the fast one, one map node per nonzero: a
30000 by 500000 file reads in 1.4 seconds, which is fine here and is a counting sort if the
matrices grow.

### The field types, 2026-08-10

The first fetch of 60 matrices, filtered on the index alone, came back:

```
36 real       usable
19 pattern    no values, refused
 5 integer    usable
```

**So the index's `isReal` flag means "not complex" rather than "the file has values"**, and only
the file can tell the two apart. This surfaced at fetch time rather than in the driver because
`fetch` prints each file's Matrix Market banner as it lands, which is worth keeping.

## Singular matrices are labeled, not hidden and not solved

**The test is the structural rank**, the size of a maximum matching between columns and rows. A
permutation contributing to the determinant needs a nonzero in every position, which is exactly a
perfect matching in the bipartite graph of the pattern; where none exists every term of the
expansion vanishes and the matrix is singular for almost any values. So a rank below `n` is a
proof of singularity taken from the pattern, with nothing factored. Those rows print it in place
of the residual columns, and again as a group at the end:

```
Newman/netscience.mtx    1589    7073    SINGULAR, structural rank 1424 of 1589

matrix                        n   s-rank  deficit    empty
Newman/netscience.mtx      1589     1424      165      128
```

`empty` counts columns holding no nonzero at all. It is kept beside the rank because it names the
common case here, the isolated vertex, and because it is a subset of the deficit: an empty column
can never be matched.

**The matching is Hopcroft-Karp**, imported from the `combinatorial-suite` package and carried in
`Matching.h` with three mechanical changes and no change to any algorithm. That header also
records why general (non-bipartite) matching is the wrong tool, and what the rest of
Dulmage-Mendelsohn is for, which is reducibility rather than singularity.

### Why they are not solved, worked through on netscience

`Newman/netscience` is a collaboration network, n = 1589, and 128 of its columns are empty:
authors who published alone, stored with no diagonal, so after conversion their column carries a
structural zero and nothing else. Its first pass looked like a defect in dynamic LDL, reading
1.25e+33 against static LDL's 4.02e-01.

It is not. **With `b` all ones the system provably has no solution.** For each of those 128 rows
`(Ax)_i` is zero whatever `x` is, while `b_i` is one, so no vector in existence gets the residual
below

```
sqrt(128 / 1589) = 0.2838
```

Static LDL returned 0.3104, which is essentially that floor. The dynamic figure varied by
thirty-three orders of magnitude between two machines running the same code, and that is the tell
rather than a symptom: when no solution exists, `x` is decided by arbitrary tiny pivots and its
magnitude is noise.

**The matching puts the deficiency at 165**, against the 128 the empty columns account for and the
171 zero eigenvalues the inertia reported after factoring. Three numbers arrived at three ways,
and they bracket as they should: 128 provable from single columns, 165 provable from the pattern,
and 171 found numerically, the last being at or below the structural bound plus whatever the
values add.

**Handed a consistent right-hand side it comes out right.** Taking a non-constant `x_exact` and
forming `b = A x_exact`, which lies in the range of A by construction and so has a solution:

```
                b = ones                        b = A x
static     3.10e-01   max|x| 1.0e+14      3.94e-03   max|x| 8.5e+00
dynamic    2.19e+00   max|x| 2.7e+18      2.17e-16   max|x| 2.3e+02
```

Dynamic LDL solves this rank-deficient matrix to 2.2e-16. Static's 1e+14 is its perturbation of
1e-14 showing through, exactly as documented. (Measured in the Linux sandbox, which is not a
timing platform; residuals are not timings and the two machines agree on everything here except
the meaningless column.)

The degenerate shapes underneath were checked separately on matrices small enough to verify by
hand: a single isolated vertex, a 2x2 with no diagonal, and mixtures of the two. The 2x2 takes a
2x2 pivot and gives an exact answer; the singular ones report the correct inertia and a residual
of 1. Nothing misbehaves.

### What the label does not catch, and why that is safe

**The structural half is now complete**, and only the numerical half is left: a matrix that has a
perfect matching, so a determinant that does not vanish identically, but whose values make it
singular anyway. Nothing short of factoring will find that, and a matrix singular by a hair is
indistinguishable from an ill-conditioned one in any case, which is a property of the problem
rather than a limitation here.

**That incompleteness would be fatal under one policy and is tolerable under this one.** If we
excluded singular matrices in order to trust the table, an undetected one would sit there looking
like a defect. What actually protects the table is the pair: label what we can prove, and be ready
to switch to a consistent right-hand side, which gives every matrix a solvable system whatever its
rank. The label is then an explanation rather than a gate.

**The empty-column test that came first was a real gap, not a smaller version of this one.** It
finds only the isolated vertex. A star graph, every vertex joined to one center and to nothing
else, has no empty column at all and a structural rank of 2 whatever its size; it went through
unlabeled until the matching landed. On `netscience` the two differ by 37, the matching finding a
deficit of 165 where empty columns account for 128.

**One wrong proposal is recorded in `Matching.h`** rather than here, since that is where someone
would reach for it: an odd component with no diagonal is NOT singular, the permutation not having
to be an involution, and a triangle refutes it.

**And the collection does not mark singularity.** `ssstats.csv` carries size, nnz, the two
symmetry measures, `posdef` and `kind`, and nothing about rank. One class is marked exactly:
`posdef == 1` implies nonsingular, and it is also the narrowest possible set, excluding every
indefinite matrix, which is where the interesting work is.

### What to change when singular matrices come back, and it is the right-hand side

Written down now because the reasoning is fresh and the change is small, and because a reader
arriving later will otherwise re-derive it from a table that no longer has the rows in it.

**The change is `b`.** Today it is all ones, which is safe only because the rows that need
better are skipped. Restoring them means forming the right-hand side from a chosen solution:

```
for i:  xExact[i] = <something not constant>
        b = A * xExact          MultiplyEngine::compute
```

Four things to get right, and each of them is a mistake this folder has either made or nearly
made:

- **`xExact` must not be constant.** With `xExact` all ones, `b` is the vector of row sums, which
  is exactly zero for any matrix whose rows sum to zero. Laplacians are the obvious family and
  they are common in this collection. A relative residual then divides by zero.
- **Guard `||b||` anyway.** Non-constant makes the degenerate case unlikely rather than
  impossible, since `xExact` can still fall near the null space. A matrix whose `||b||` is
  negligible against `||A|| ||xExact||` should be labeled, exactly as a rank deficit is, rather
  than divided by.
- **Do not check the forward error.** On a singular matrix the solution is not unique, so a
  correct answer may differ from `xExact` by any null-space vector while having a residual at
  machine precision. That is the method working, not failing. The residual is the thing to
  report and `||x - xExact||` is not a check.
- **Know that it is a slightly easier test.** `b` formed this way inherits A's scaling, where an
  arbitrary `b` does not. That is a fair criticism of it as a stress test and not a problem for a
  correctness pass, which is what this folder is.

**What it buys is more than the labeled rows.** Every matrix gets a solvable system whatever its
rank, so the singular ones the empty-column test cannot see stop producing figures that read as
defects. That is why this change and the label are a pair rather than alternatives: the label
explains a row, and the consistent right-hand side is what makes the rest of the table trustworthy
without a complete detector behind it.

The line to change is in `matrix_accuracy.cpp`, where `b` is filled, and its header comment points
back here.

## Why some matrices order slowly, 2026-08-11

The first performance run turned up ordering times that looked wrong: `GHS_indef/bloweybq` at 148
ms under AMD3 where the vendored AMD took 0.38, `Lourakis/bundle1` at 121 against 1.93,
`Mulvey/finan512` at 347 against 6.42. **Two different causes, and neither is a defect.**

### One is a missing feature, and its remedy is documented in the vendored source

**A single dense row or column makes minimum degree quadratic.** `private/Amd.cpp` says so in as
many words at the top of the file: "the presence of a dense row/column can increase the ordering
time by up to O(n^2), unless they are removed prior to ordering". Its remedy is a threshold,
`max(16, 10 * sqrt(n))` by default, above which a row is called dense, pulled out before ordering
and placed last in the output.

**Oblio implements no such rule**, in `Amd3.cpp` or in `QuotientGraph`, and neither does MMD3 or
the vendored genmmd. So all three carry a vertex adjacent to everything through every degree
update, and the vendored AMD does not.

Removing every column above AMD's own threshold, alpamayo, milliseconds:

```
                     n  threshold  dense    MMD             MMD3            AMD           AMD3
bloweybq         10001       1000      1    70.70 -> 0.83  145.37 -> 1.31  1.36 -> 1.35  470.32 ->
1.54
bundle1          10581       1028    252   259.34 -> 10.98 322.25 -> 15.93 7.47 -> 2.46  107.30 ->
11.12
```

The vendored AMD barely moves, having already done this internally; everything else falls by one
to two orders of magnitude. **`bloweybq` has exactly one column of degree 10000 and 9992 columns of
degree 5**, so a single vertex accounts for a factor of 300.

**Fill is unaffected**, 39996 against 39997 across all four on `bloweybq`, so this is time and not
quality. A dense-row threshold in `QuotientGraph` would fix all six ordering drivers at once and is
its own piece of work.

### The other is minimum degree's own weakness, and it is not ours

`Mulvey/finan512` has **no dense column at all**: maximum degree 54 against a threshold of 2734,
and removing nothing changes nothing, 815 ms against AMD's 28. It also fills 6.5 million against
AMD's 2.8 million, **losing on both axes at once**.

Its degree histogram says why: 512 columns at 54, 512 at 51, 512 at 22, 1024 at 20, 24064 at 6.
That is a nested block structure with massive degree ties, which is exactly where exact minimum
degree spends its time and where AMD's approximate degree does not. **Our MMD3 is nearly twice as
fast as the vendored genmmd here**, 452 against 815, so this one is not ours to answer for.

**Two lessons for reading any ordering timing.** A minimum degree ordering's cost depends on the
degree distribution far more than on n or nnz, so a mean over a set says little about any row in
it. And the two families fail differently: AMD is protected against the dense case by a feature and
exposed to nothing else here, where MMD is exposed to both.

## What the runs showed, 2026-08-10

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate.** MMD3, left-looking throughout. Four
runs, at 60, 107, 140 and 160 files. The 160-file run is the current record:

```
of the 160 files read:
  pattern, no values      41
  with values            119
    structurally singular  12
    killed                  1
    solved                106

  positive definite       30
  negative definite       18
  indefinite              56
  numerically singular     2
```

**Widening the set widens the mix rather than deepening it.** Indefinite went 14, 28, 46, 56 across
the four runs and positive definite 17, 17, 25, 30, which is what `--per-kind` was meant to do.

**But the accuracy set is now structurally NARROWER than the pool.** 35 kinds at `--per-kind 2`
against 28 at `--per-kind 8`, because `--values` removes the graph kinds entirely, a graph being
stored as a pattern. The 41 pattern files carry variety the valued ones cannot, which is an
argument for any ordering study using the whole pool rather than this subset.

**No alarm row at any size.** Every large residual sits beside a `bwd` at machine precision. At 106
solved the largest are `eurqsa` at 2.0e+16 with `bwd` 5.8e-18, `TEM27623` at 1.1e+04 with 3.5e-18,
`sit100` at 2.6e+03 with 9.6e-15 and `t2dal_bci` at 1.4e+03 with 1.1e-16.

**Cholesky and the inertia agree on every row of every run**, 35, 58, 88, then 106. That check has
held across four independent samples and it is the only place one part of the library is tested
against another.

**And the largest-residual block has no cliff**, which corrects what the first run suggested. At 35
solved it read as three rows and then four orders of magnitude. It is one outlier and a continuum
of ill-conditioned matrices, which is why the block is a ranking: a cutoff chosen against the first
run would have been wrong by the second.

**The price of dynamic pivoting is visible per row.** The ratio column is 1.0x across most of the
table and reaches 232.4x on `LFAT5000`, 19.9x on `eurqsa`, 3.5x on `c-66b`, 3.2x on
`spaceStation_10`, 3.1x on `c-23`. Every 1.0x row has no `Nd` in its note, which is the column
checking itself.

**Three things worth keeping from the earlier runs**, none of them a defect:

- **Cholesky refuses on most of the set**, which is correct. Roughly a quarter of these matrices
  have a structurally absent diagonal, and a graph adjacency matrix is not positive definite.
- **Dynamic LDL is doing exactly what it exists for.** `Gset/G32` goes from 1.21e+00 under static
  LDL to 1.34e-13 under dynamic; `Schenk_IBMNA/c-18` from 5.69e-09 to 3.57e-15; `ML_Graph/mice_10NN`
  from 1.65e-02 to 1.84e-15. All three have no diagonal at all.
- **What looked like conditioning was singularity.** The first run recorded `HB/saylr3` at 8.9e-02
  beside `HB/sherman1`, identical in pattern and nnz, at 9e-14, and read the gap as conditioning.
  It is not: `saylr3` is numerically singular by 2, and the section above works it through.

**Free oracles keep appearing without being arranged.** `Nemeth/nemeth02` through `nemeth09` are
consecutive members of one sequence, and the first seven come out at nnz(L) 227391 exactly while
`nemeth09` differs at 227568, its nnz(A) differing too. `Schenk_IBMNA/c-66` and `c-66b` share a
pattern, `Bindel/ted_B` against `ted_B_unscaled` is the same matrix scaled, and `QY/case9` and
`TSOPF/TSOPF_FS_b9_c6` are the same matrix under two names, agreeing to every digit in every
column. `Oberwolfach/t2dal`, `t2dal_a`, `t2dal_bci` and `t2dal_e` are one descriptor system in four
parts, three sharing a pattern at 117365 fill and the mass matrix diagonal at 4257. A structural
disagreement inside any of these would be a defect by construction.

**And structural mass matrices keep turning up rank deficient.** `m3plates`, `bcsstm13` with 762
empty columns, then `bcsstm38` with 2833. That is a property of how mass matrices are built rather
than a curiosity, and this folder has now found several without looking for them.

**The `GHS_indef/aug*` family is structurally singular as a family**, `aug2d`, `aug2dc` and `aug3d`
all declined, with deficits of 9800, 10200 and 12636 and **no empty columns at all**. That is the
matching earning its place: the empty-column test would have passed every one of them and their
residuals would have been noise in the table.

## What the performance run showed, 2026-08-11

**alpamayo (Apple Silicon), macOS, Apple Clang, Accelerate.** Cholesky, real, best of three after a
warm-up, four orderings and three traversals over the positive definite set.

```
of the 114 files read: 107 measured, 1 not positive definite, 5 over the fill cap, 1 skipped

geometric mean relative to the best ordering on each matrix:
  order      nnz(L)      order    analyze     factLL
  MMD         1.032      1.566      1.335      1.081
  MMD3        1.032      1.821      1.416      1.074
  AMD         1.020      1.126      1.045      1.073
  AMD3        1.020      1.529      1.282      1.040

traversals at MMD3, relative to the best traversal on each matrix:
  left-looking   1.364
  right-looking  1.388
  multifrontal   1.028
```

**Multifrontal wins, and by a wide margin**: 1.028 against 1.364 and 1.388. That is the clearest
result in the table and it held at both set sizes, 87 matrices and then 107.

**The fill difference between orderings is one to three percent.** On grids `../ordering` measures
up to 13 percent. That gap is a property of grids, where nearly every live vertex has the same
degree and the tie-break decides almost every pick; on real structure the ordering choice barely
moves the fill.

**AMD fills slightly less than MMD here, 1.020 against 1.032**, which is the opposite of what square
grids say. Both differences are small enough that the honest statement is that neither wins.

**Analysis is where the spread lives**, 1.045 to 1.416, and the section above explains most of it:
the dense-row rule the vendored AMD has and nothing else does.

**The fill cap fired five times**, including `FlowIPM22/uni_chimera_i1` at a predicted 1.2
**billion**
entries, about 10 GB of values. Without it that would have been another `bloweya`.

**And the vendored pair is an oracle in every row.** MMD and MMD3 agreed on nnz(L) on all 107
matrices. **AMD and AMD3 did not**, differing on a minority of rows, once substantially:
`HB/bcsstk08` at 31153 against 29922, where ours fills 4 percent less. The acceptance tests in
`experiments/ordering` have AMD3 reproducing the vendored raw order exactly on all 38 cases, but
those are grids and random patterns. This is the first evidence from real structure and it is
unexplained; it is recorded here rather than chased.

## What the ordering run showed, 2026-08-15

The first `make ordering`, all 246 matrices, `Mmd3` against genmmd on the ordering step alone.

**Where there is work to do, we win, and the margin grows with the work.** 0.83 to 0.85x on the two
PARSEC giants, 0.81x on `Schenk/nlpkkt80`, 0.60x and 0.40x on the `FlowIPM22/uni_chimera` pair,
0.45x on `AG-Monien/se`, 0.53x on `Mulvey/finan512`, 0.57x on `Lourakis/bundle1`. The largest of
these is `uni_chimera_i5`, 6.3 seconds against 15.6.

**Where there is nothing to do, we lose, and the pattern is sharp.** The worst five rows are all
matrices with NO OFF-DIAGONAL ENTRIES, `nnz(A) = n` and `nnz(L) = n`: `Boeing/bcsstm39`,
`Cunningham/m3plates`, `HB/bcsstm25`, `Oberwolfach/t3dl_e` and `Oberwolfach/t2dal_e`, reading 2.03
to 2.28x. The tier below, 1.5 to 1.9x, is the same thing weaker: `Bai/mhd3200b` at 3.4 entries of L
per column, `Oberwolfach/LFAT5000` at 3.4, `GHS_indef/linverse` at 4.5.

**So the shape is a higher per-vertex constant and a lower per-unit-of-work cost**, and on a pure
diagonal the constant is the entire run. One term of it was found and removed the same day: the
prepass collected the degree-1 bucket into a vector and then walked it, where genmmd reads each
successor before unfiling and needs no list. That was worth 8.6 percent of a pure-diagonal ordering
and 0.2 percent of a grid, and it moved those five rows from 2.5 to 2.8x down to 2.0 to 2.3x.

**What remains of that constant is CONSTRUCTION**, and it is measured: with no elimination work at
all, an `Mmd3` ordering is roughly a third `QuotientGraph` construction and a sixth
`orderAscending`. Construction allocates and initializes about ten size-n arrays where genmmd
allocates five plus its 1-based copies, which is the array-count finding of that morning moved into
the constructor, invisible on a grid because real work amortizes it. `docs/NEXT.md` item 8 carries
it as an experiment rather than as an optimization, and the reason is in these numbers: the
matrices in question order in tenths of a millisecond and need no ordering at all.

**Two figures worth keeping from the fill column.** `PARSEC/Si87H76` produces 5,679,875,732 entries
of fill under MMD3, confirming to three digits a prediction this tree had only extrapolated from
grids. And there is a worse case nothing had recorded: `FlowIPM22/uni_chimera_i1`, n = 100000 and
nnz(A) = 1100592, produces 1,179,373,506, which is 1072x nnz(A) from a matrix an eighth of
`Si87H76`'s size. Both are arguments for nested dissection, and `uni_chimera_i1` is the cheaper one
to experiment on.

## What this folder still needs

- **A consistent right-hand side**, `b = A x_exact` with a guard on `||b||`, which is what lets
  the labeled rows back into the table and protects it against the singular matrices the label
  cannot see. What to change and what to watch for is written out above, under "What to change
  when singular matrices come back".
- **`DirectSolver` does not forward `rank()`, and the `Nz` note is standing in for it.** It
  forwards `numPerturbations`, `numDelayedColumns`, the two pivot counts and `inertia`, but not
  the numerical rank, which `NumFactorDynamic` maintains and which names this situation directly
  rather than through a count whose own header warns it is least reliable on singular input. That
  is a library change rather than a benchmark one, and it is the first thing that would improve
  this table.
- **A dense-row threshold**, which is the single largest ordering-time item and the one with a
  documented remedy. It belongs in `QuotientGraph` so that all six ordering drivers gain it at
  once. See "Why some matrices order slowly" above.
- **Why AMD3 and the vendored AMD disagree on fill** on a minority of real matrices, where the 38
  acceptance cases in `experiments/ordering` show exact agreement. Ours fills less where they
  differ, so it is not urgent, but it is a divergence the acceptance tests cannot see.
- **The scaling report**, the third of the three: synthetic 2D and 3D grids, Cholesky throughout,
  covering real and complex Hermitian, which real matrices cannot supply, since the collection
  marks one complex matrix in 23 as positive definite and that one is complex symmetric.
- **Complex matrices.** None are fetched. They are worth having, complex Hermitian dynamic LDL
  being the one part of the numeric code that is an extension rather than a port, so its only
  oracles today are the residual and reconstruction on generated input. Fetching them needs a
  separate branch in the filter rather than a relaxed threshold, for the reason in "What the
  collection's index does and does not mark" above. The reader refuses `complex` today, so any
  that land in `data/` skip harmlessly.
- **The other orderings and traversals.** One ordering and one traversal here, deliberately: they
  are the subject of the performance report rather than a gap in this one.
- **The 41 pattern files.** They have no values and are exactly what widening
  `experiments/ordering`'s `make amdorder` to real patterns would want. They also carry the graph
  kinds the accuracy set no longer has.
- **`fetch` should look names up in the FULL index, not the filtered one.** It currently reports
  `not in the filtered index` for a matrix outside the size range or the symmetry test, which is
  wrong once a list exists: the filter's job is to produce candidates, and a list handed to
  `fetch` IS the selection. This blocks using a set somebody else recommends.
- **Driver-side selection beyond the glob and the fill cap.** `nnz(A)` is enough as an indicator
  for now. The point is that the run chooses over the pool rather than the pool being curated per
  run.
- **More matrices.** 160 files is one sample of one filter, and the kinds are saturating: 66
  valued at `--per-kind 4`, 99 at 6, 119 at 8. Growing further means raising `--max-nnz` rather
  than `--per-kind`, since several kinds are thin only because their cheap end sits above the cap.
  The larger end of the collection is untouched.
- **Why dynamic pivoting cascades on a chain**, which is the largest open question this folder has
  raised and the only one that is about the library rather than about the benchmark. `LFAT5000`
  delays 3.1 million columns on a POSITIVE DEFINITE matrix where no pivot needed delaying, and
  `bloweya` does the same past the end of memory. The section above has the controls and the
  argument; the question belongs with the pivoting strategy.

**What is deliberately NOT owed.** The block triangular form. It is the other half of
Dulmage-Mendelsohn and it answers reducibility rather than singularity, and on symmetric input it
degenerates: a symmetric permutation cannot produce a nontrivial triangular form, so the blocks
are the connected components and a depth-first pass gives them. `Matching.h` carries the argument.

**What the matching costs**, since it now runs on every matrix before anything else: 85 ms on a
30000 by 970000 pattern, 1 ms at n = 10000, and it took the whole pass over a 100-grid from 68 ms
to 96. Cheap against a factorization, and worth watching if the set grows toward the large end.

**One thing deliberately NOT here: a skip list.** It would unblock a run in five minutes and it
records a decision about one matrix on one day, going stale silently while the next such matrix
takes the run down again. If one is ever added as an interim, every entry carries a reason and the
run prints them as a block, so nothing passes over quietly.

**Record the machine, the date and the compiler beside any numbers kept here.** A timing without
them cannot be compared against anything later, and on this folder's evidence a residual on a
singular system cannot be compared against anything at all.
