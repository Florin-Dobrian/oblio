# TODO

Work we intend to do that is not a porting question. The porting ledger tracks whether a 0.9 unit
has been carried over and checked, which is a question with an end; this file holds everything that
outlives it, including things 0.9 never did and that we want anyway.

Two rules keep the split clean. If the item is "0.9 has this and we do not yet", it belongs in
PORTING_LEDGER under Owed. If the item is "neither reference has this, or we want it done
differently", it belongs here. Where an item straddles the two, it lives here and the ledger points
at it, because this is the file meant for browsing.

Each entry says what, why it matters, and what makes it urgent. Order within a section is rough
priority, but nothing here is scheduled.

## Robustness

### Validate the input matrix

Nothing checks that `A` is a well-formed input, and three distinct properties are assumed silently.
Today we generate our own test matrices, so the assumption holds by construction; it stops holding
the moment anyone else calls the solver.

**The diagonal must be structurally present**, even where it is numerically zero. Symbolic
factorization builds a column's index set from A's column structure, so a column with no stored
diagonal never enters its own index set, and the resulting supernode index sets are wrong. This is
not a corner case: a zero diagonal is completely ordinary in indefinite problems, which are exactly
what dynamic LDL exists for. A KKT matrix whose zero block is stored sparsely, dropping the zeros,
hits it immediately. Found on 2026-07-19 while validating dynamic LDL, where it produced both wrong
answers and spurious refusals before the cause was clear.

**No duplicate entries.** `assembleFromA` and `assembleDelay` assign rather than accumulate,
which is correct when each position has one writer and silently keeps only the last value if a
column stores the same row twice. 0.9 accumulates in the matrix case and so tolerates duplicates.
The header already documents CSC with sorted, and by implication unique, row indices, so this is an
unchecked precondition rather than a disagreement.

**Symmetry, and which kind.** Complex Cholesky requires a Hermitian matrix and complex LDL^T
requires a complex-symmetric one, and neither is verified. The failure is silent: `zpotrf` reads
only the lower triangle and assumes the upper is its conjugate, so a complex symmetric matrix is
factored as though it were the Hermitian matrix agreeing with its lower triangle, succeeds, and
returns a plausible wrong answer. The agreed fix is two flags computed on `SparseMatrix` in one
construction pass, `isHermitian()` and `isSymmetric()`, with the engine requiring the first for
Cholesky. For `double` the two coincide; for complex a matrix may be one, the other, or neither.
This entry supersedes the Hermitian item under Owed in PORTING_LEDGER, which is the same fix seen
from the porting side.

The three want one validation pass and one place to report from, so they should be done together
rather than one at a time. Trigger: before anyone outside this effort runs the solver.

## Capability

### Complex Hermitian dynamic LDL: DONE, 2026-07-19

Recorded here as done rather than deleted, because it is the one factorization in the library with
no reference behind it and that is worth being able to find later. 0.9's complex LDL is symmetric
only, so nothing was transcribed; the oracles are the residual in `test_pipeline` and reconstruction
of `L D L^H` on dense fronts.

Complex `DynamicLDLT` turned out to need no kernel change at all, 0.9's complex `factorDynamicLDL_`
differs from its real one in six lines, all declaring the pivot magnitudes real rather than scalar,
which this port had done from the start, so only the dispatch guard widened.

`DynamicLDLH` needed four changes, one place each thanks to the pivot-body merge landing first:
the conjugate in `readPivotBlock2x2`, conjugated `L` where the `D L^H` rows are formed in both
eliminations, and `forceReal` on the diagonal.

**And one that was not in the plan, in `swap`.** Two of its loops move values across the diagonal
and needed conjugating, which was anticipated; the one that was not is the entry *between* the two
swapped positions, `block[at(k_, j_)]`. It is its own reflection under the permutation, so no loop
touches it, and for a symmetric factor it is invariant, which is why 0.9 leaves it alone and why
the omission is invisible until the factorization is Hermitian. Left out, the factor reconstructs
the *conjugate* of the matrix in the affected rows, with no crash, no NaN, and no symptom at all
when the swaps happen to be the identity. Roughly a third of random Hermitian matrices failed.

The lesson worth keeping: **a symmetric-to-Hermitian port is not a search for `conj` opportunities
in the arithmetic.** The arithmetic was nearly all fine. The bug was in a data-movement verb, in an
entry whose symmetry made it a no-op for every case the reference ever ran.

### Dynamic multifrontal: DONE, 2026-07-21

Recorded as done rather than deleted, because it was the last and hardest traversal, the one where
delayed columns meet the update stack, and 0.9 needed three attempts at its factor kernel before the
live one settled. The port needed none of that drama: the multifrontal factor kernel turned out to be
the same computation as our `factorDynamicSupernode`, reused unchanged. The driver is the left-looking
dynamic skeleton with the `std::vector<UpdateMatrix>` stack in place of the pull queue; assembling a
child does both halves of the assembly (`assembleDelay` for its delayed columns, `assembleUpdateMatrix`
for its contribution block), and the block is formed by the new `updateDynamicUpdateMatrix`
(`formDynamicUpper` + `gemmLower`). Verified by residual at all three tiers, real and complex,
symmetric and Hermitian; tier 1 matches left-looking's counts exactly. See PORTING_LEDGER for the
mechanism. With this the whole order-by-factorization-by-traversal matrix resolves to a residual.

### Static factorization into dynamic storage through multifrontal

The one cell left returning false without an assertion behind it. `NumFactorEngine::compute` into
`NumFactorDynamic` dispatches Cholesky and static LDL to left- and right-looking but not to
multifrontal, where it still returns false. This is a convenience path, not a missing capability: the
same static factor is reached through flat `NumFactorStatic` storage, which multifrontal already
serves. `factorStaticMultifrontal` is templated on the factor type and the other two static traversals
already instantiate for dynamic storage, so wiring this is likely a one-line dispatch change plus a
test; the reason to leave it for now is that nothing needs a statically pivoted factor in
delayed-column storage. Trigger: if a caller ever wants one factor object that can hold either.

### Multiple right-hand sides### Multiple right-hand sides

`Vector<Val>` carries one right-hand side and the solve is scalar, which is the right call for one
column: there is no level-3 BLAS to be had. Many right-hand sides make a supernode a matrix
operation, and then gathering the rows into a dense block, calling TRSM and GEMM, and scattering
back pays for itself. 0.9 has this as its `MultipleVector` path, so the algorithm is available; what
is missing on our side is `DenseMatrix`, which the ledger lists as a unit and which nothing has
needed yet.

### The update stack is not a true LIFO, and may never need to be

Our multifrontal update stack is a `std::vector<UpdateMatrix>` indexed by supernode, with
arbitrary-slot allocate and discard, not the classical push/pop arena moved by a single stack pointer.
DESIGN_DECISIONS (2026-07-22) records the full reasoning. This entry began as a note tracking that
departure as a debt to be repaid; the posture below is weaker than that, and deliberately so.

**The arrangement we have is already the pure-parallelism design.** One slot per supernode, filled
when the supernode is reached and freed when its parent assembles it, is disjoint by construction: no
shared stack pointer, no arena to size, no per-thread bookkeeping, and concurrent branches touching
memory that cannot overlap. Nothing about it needs changing to run branches in parallel. So if the
storage the update matrices occupy is not a binding constraint, this item is closed rather than
deferred, and what is recorded elsewhere as a departure from 0.9 is simply the right answer for a
concurrent in-core solver.

**Whether it binds is a measurement, and it is the only thing gating this entry.** The number wanted is
the update stack's share of peak memory, against the size of L and against the total. A few percent and
the item is dead and the current arrangement final; a third of the peak on a 3D problem and it is not.
That is cheap instrumentation, and it is the same discipline the OpenMP experiment applied to the
threading question, where measuring turned a documented instruction into a corrected no-op.

**The reason to keep the question open at all is that memory, not speed, usually decides how large a
problem a sparse direct solver can take.** That is what out-of-core multifrontal was invented for, and
0.9's abstract `UpdateStack` with its out-of-core concrete `UpdateStackDynamic` is the reference shape
for it. The factor L normally dominates the total, but the update stack is the part that is pure
scratch: it holds no results and is gone when the factorization ends. So if memory binds, the stack is
the first place to look, precisely because shrinking it costs nothing in output.

**Out of core is the one regime where the LIFO discipline is required rather than merely tidy, and
the reason is the shape of a file.** A file is a linear sequence of bytes that can be appended to and
truncated at its end cheaply and at neither of those cheaply anywhere else. A stack maps onto that
exactly: pushing is an append, popping is a truncate, and the file's one end is the stack's one end.
Our per-supernode slots do not map onto it at all, since freeing a block in the middle of the run
would mean freeing bytes in the middle of the file, and the ways out of that are all bad. Compacting
rewrites the tail on every free. Leaving holes and tracking them is a free-space allocator on disk,
written by hand. One file per update matrix avoids both and is overkill, thousands of files with the
open, close and metadata cost of each, for blocks that live a few steps.

So out of core does not merely prefer an arena; it is the case where nothing else works, and the
requirement comes from the medium rather than from any property of the factorization. The postorder
is what makes the mapping available, since nesting live intervals is exactly the condition under
which a one-ended structure suffices, and that in turn is why the discipline collapses under
concurrency: interleaved branches destroy the nesting, so per-thread stacks each spilling their own
file is the shape an out-of-core parallel solver would need.

One consequence is worth keeping: the resident part of a file-backed stack is always a suffix, the
most recently pushed blocks, so deciding how much to hold in memory is sizing a single buffer rather
than choosing among blocks. Sequential transfers come along as well, but that is the append pattern
restated rather than a second argument for it, and contiguity in core is a different question, the
fragmentation one below, on which our scattered per-supernode allocations do not deliver anything.

**The shortest argument against a stack in core needs no measurement at all: the hole.**
`updateMatrix[K]` is allocated while its children are still live, so it necessarily sits above them,
and freeing them opens a gap beneath it. In an arena there are only two ways out, moving K's block
down into the gap, which is a copy of `U(K)` bytes per supernode, or leaving the gap and abandoning
the bump pointer, at which point it is not a stack. With separate allocations there is no gap to
close, because there was no adjacency to preserve in the first place.

The cost of that compaction is what settles the question, because it is bought with different things
in the two regimes. In core it buys contiguity, which we have argued is worth little, while the peak
does not move at all: paying a copy per supernode for nothing. Out of core the same copy is what
keeps the file a file, so it buys the ability to run at all. The compaction is therefore the visible
symptom of a discipline applied outside the regime that pays for it.

**The three regimes, then, and they do not agree.** In core and serial, random-access slots are fine:
the peak is identical to a stack's and the child ordering is what pays. In core and parallel, the
slots are actively better, disjoint by construction with no shared stack pointer, and a single global
arena is incoherent once branches interleave. Out of core, the stack is not an optimization but a
requirement, a file having only one cheap end. So the arena is not a general good being deferred;
it is specifically the out-of-core enabler, plus a fragmentation hedge in core, which is narrower
than a plain reading of this entry suggests.

**If it does bind, the target differs between the serial and the parallel case.** Serially the answer
is unambiguous: one arena, full LIFO discipline, the smallest peak the traversal admits. Under
concurrency the peak rises for a reason that has nothing to do with storage management, since N
branches in flight hold N live sets where serial execution holds one, and per-thread arenas add a
second and smaller cost because each must be sized for its own worst case and cannot lend slack to a
neighbor. That increase is intrinsic to parallel multifrontal rather than a flaw in the arrangement,
and it is the same effect that makes MUMPS's memory grow with rank count. A single global arena stops
being an option once branches interleave, since the push and pop sequence is then not LIFO at all and
one stack pointer would have to serialize the branches to remain a stack. So the shape is full LIFO
serially, per-thread LIFO in parallel, and a bit more total storage as the price of the concurrency.

**Two things keep that increase modest.** What gets multiplied is the small blocks and not the large
ones: tree parallelism is available near the leaves, where the forest is wide and contribution blocks
are small, while near the root, where blocks are enormous, the tree has narrowed to a few fronts and
there is little concurrency left to exploit. The multiplier therefore applies where blocks are cheap,
and the expensive blocks near the root stay singular. And the multiplier is a knob rather than a
constant, since the number of concurrently active branches bounds the increase directly: capping
active tasks caps the peak, which lets the tradeoff be tuned at run time rather than fixed at design
time. A shared pool handing fixed chunks to per-thread arenas would recover some of the remaining
slack without reintroducing a shared stack pointer.

**The lever that actually lowers the peak is child ordering, and it is orthogonal to all of the
above.** Processing a parent's children in decreasing order of contribution-block size lowers the
stack peak, serially and per thread alike, whatever the storage arrangement, because it is the order
of the pushes and not their management that determines how high the stack climbs. That is 0.9's
`sortForOptimalMultifrontal`, recorded as unported in PORTING_LEDGER, where the stated reason for
deferring it (not needed until the numeric factorization exists) no longer holds now that multifrontal
is complete for every factorization. If the goal is to keep the extra storage small, it is the
cheapest real win available and it does not touch the LIFO question at all.

**The trap, for anyone who does reintroduce stack discipline: reverse the assembly order but not the delay
prepend.** The assembly is commutative and free to reorder; the delay prepend is structural, since it
writes kk's index set and so fixes the front's column order that the pivot sequence reads and that the
pinned tier-1 counts depend on. A LIFO refactor is right for the assembly and wrong for the delays, and
the two cannot be swept together.

**One reason to pool survives even if this entry stays shut, and it is not about memory.** Per-call
allocate and free at leaf granularity, from many threads at once, is allocator contention, which is a
throughput problem rather than a footprint one. The forest parallelism entry below already flags it for
`updateStaticUpdateMatrix`'s per-call `upper` buffer, and update matrices would face the same pressure.
Per-thread pools could therefore earn their place on speed grounds alone, the same mechanism as an
arena with the opposite justification. The two should not be conflated: a decision to pool the update
matrices may be right while the LIFO question stays permanently closed.

Trigger: a measurement showing the update stack to be a large share of peak memory. Absent that,
closed rather than deferred.

### Experiment: update-matrix storage against factor storage, and how much the traversal changes it

The measurement the entry above turns on, written down so it can be picked up cold. It belongs in
`experiments/`, alongside `ordering/` and `openmp/`.

**It is entirely symbolic, which is what makes it cheap.** Both quantities come from the symbolic
factor: a factor block is `(frontSize + updateSize) * frontSize` and an update matrix is
`updateSize * updateSize`, full square rather than triangular, and the lifetimes come from the tree.
So the driver needs `SymFactorEngine` and nothing below it, no numeric factorization, no BLAS, no
right-hand side, and it runs on problems far larger than we would care to factor.

**The live-set rule is fixed by the code and should be simulated, not assumed.**
`factorStaticMultifrontal` allocates `stack[K]` when K is reached, then assembles each child and
discards it one at a time, so every child block is live when the parent allocates its own. Writing
`U(X)` for the size of X's update matrix, `updateSize(X)` squared, and letting I and J range over K's
children with `I < J` meaning I is assembled before J,

```
peak(K) = max( max_J [ sum_{I<J} U(I) + peak(J) ],     term 1, order-dependent
               sum_J U(J) + U(K) )                     term 2, order-invariant
```

**Term 1 is the running term: the worst moment while a child's subtree is still working.** When J is
running, two things are resident. Its earlier siblings have finished and parked their update
matrices, which is `sum_{I<J} U(I)`, and J's own subtree is at its own worst instant, which is
`peak(J)`, recursively. The maximum is over which child is running. This is the only place the child
order enters, through `sum_{I<J}`: change the order and each child faces a different pile of parked
matrices, so a child that needs much should run while little is parked.

**Term 2 is the assembly term: the single instant at K's own step.** Every child has finished, so all
of their update matrices are parked and none is freed yet, and K allocates its own on top of them.
Not a maximum over anything, and no order can touch it, being a total over all the children with
addition commutative. It is also the only place `U(K)` appears, the join where every child is
simultaneously resident, and the non-LIFO moment, since `U(K)` is allocated above children that then
die beneath it. The recursion bottoms out at a leaf, which has no children, so `peak(leaf) =
U(leaf)`.

**Why a sort is enough.** Minimizing term 1 ranges over all orderings of the children a priori, but an
exchange argument collapses it to a single key. Take two adjacent children X and Y, and let S be the
sum of U over everything assembled before them. Their two positions contribute

```
P1  (X first) = max( peak(X),  U(X) + peak(Y) )
P2  (Y first) = max( peak(Y),  U(Y) + peak(X) )
```

with S dropped, appearing in both, and everything after the pair untouched, the running sum growing
by `U(X) + U(Y)` either way. The claim is that `peak(X) - U(X) >= peak(Y) - U(Y)` implies
`P1 <= P2`, so X first is never worse.

One term of P2 dominates both terms of P1. Write `T = U(Y) + peak(X)`, the second argument of P2, so
that `P2 >= T` by definition of the maximum. Then

```
peak(X)        <= T     since T = U(Y) + peak(X) and U(Y) >= 0
U(X) + peak(Y) <= T     this is the hypothesis rearranged, U(Y) + peak(X) >= U(X) + peak(Y)
```

Both arguments of P1 are at most T, so `P1 <= T <= P2`. Any out-of-order adjacent pair can therefore
be swapped without loss, the key `peak - U` induces a total order, and sorting by it decreasing is
optimal with no search. That is Liu's rule, and 0.9's key exactly.

Against the worked example below, where A should go first: `T = U(B) + peak(A) = 310`, `P1 = max(210,
160) = 210`, `P2 = max(150, 310) = 310`, and the chain `210 <= 310 <= 310` holds.

**What it costs.** `O(m log m)` per node of m children, and since the child counts sum to at most the
supernode count, `O(snodeSize log m)` overall with m tiny in practice. It runs once at forest
construction, is purely structural and so is reused by every numeric factorization sharing the
pattern, and the `peak` values it needs come free in the same upward pass, children being visited
before parents in label order.

**The greedy is exact because the substructure is.** `peak(J)` is a single number summarizing
everything below J, and a parent needs nothing else from inside J: the internal ordering of J's own
subtree does not change how J and its siblings should be sequenced. So minimizing each child
independently and then ordering them at the parent yields the global optimum, and `peak(root)` is the
whole factorization's peak rather than a local quantity.

**A worked example, since the recursion is easier to trust once it has been run by hand.** Take the
forest

```
        K   U(K)=0        (root)
       / \
      A    B              U(A)=10    U(B)=100
      |    |
     A1    B1             U(A1)=200  U(B1)=50
```

with subtree A costing 10 units of time, subtree B costing 6, and K itself 1.

Serially, bottom up. A leaf has no children, so `peak = U`, giving `peak(A1) = 200` and
`peak(B1) = 50`. A and B have one child each, so there is no ordering decision at either:

```
peak(A) = max( peak(A1), U(A1) + U(A) ) = max( 200, 210 ) = 210
peak(B) = max( peak(B1), U(B1) + U(B) ) = max(  50, 150 ) = 150
```

At K there are two children, so the sort applies. The keys are `peak - U`, which is 200 for A and 50
for B, so decreasing key puts A first:

```
term 1 = max( peak(A), U(A) + peak(B) ) = max( 210, 10 + 150 ) = 210
term 2 = U(A) + U(B) + U(K)             = 10 + 100 + 0         = 110
peak(K) = 210
```

Taking B first instead gives `term 1 = max( 150, 100 + 210 ) = 310`, so the ordering is worth 100
here, and the mechanism is visible: A needs a great deal and leaves almost nothing, so it should run
while nothing is parked, which is exactly what a large `peak - U` says.

Now choose a schedule: A and B run concurrently at K, everything else sequential. K evaluates the sum
form:

```
term 1 = peak(A) + peak(B)  = 210 + 150 = 360
term 2 = U(A) + U(B) + U(K) = 110
peak(K) = 360
```

```
                     peak     time
serial, sorted        210       17        (10 + 6 + 1)
serial, unsorted      310       17
parallel, A || B      360       11        (max(10, 6) + 1)
```

1.71 times the memory for 1.55 times the speed, and under this schedule **the sort buys nothing at
K**: with the children concurrent there is no `sum_{I<J}` for it to act on, and both orders give 360.
The ordering still does its work at every node that remains sequenced, which here is inside A and
inside B.

**One structural fact the example makes visible: a parent whose children are all leaves has no
ordering decision at all.** A leaf has `peak = U`, so its key `peak - U` is zero, every key ties, and
the arithmetic agrees, term 1 collapsing to `sum_J U(J)`, which is order-invariant, and term 2 being
that plus `U(K)`. The ordering only does work where the children are themselves internal nodes with
depth beneath them, which is why the banded matrices in the table below show a saving of zero while
the grids show 16 to 38 percent.


An update matrix is `updateSize(K)` squared, and `updateSize` is fixed by the symbolic factorization
and never written afterwards: `mUpdateSize` is read-only throughout `NumFactorEngine`, `expandNodeIdx`
only lengthens the index array while the new entries become *front* columns, and `resetVal` sizes a
block as `frontSize * (frontSize + updateSize)` with `updateSize` read. The invariant recorded
elsewhere, that block height `frontSize + delaySize + updateSize` is preserved because
`factorDynamicSupernode` decrements `frontSize` by exactly what it sets `delaySize` to, is the same
fact seen from the other side. **Delayed columns move between `frontSize` and `delaySize` and never
touch `updateSize`.**

Three consequences follow, and they are all useful.

The stack peak is *identical* for dynamic and static LDL on the same forest, not approximately but to
the element, so every measurement in this entry applies unchanged to dynamic multifrontal.

The child ordering is exactly as valid for dynamic LDL as for static, and for the same reason: Liu's
rule runs on symbolic `updateSize` at forest-construction time, and nothing pivoting does can
invalidate its inputs. No re-derivation is needed for the dynamic case.

And the update stack is predictable before the values are known, even under dynamic pivoting, which
is the surprising one, pivoting being what normally destroys predictability. What does grow with
delays is the factor, since `NumFactorDynamic`'s blocks widen as fronts absorb delayed columns. So
`|L|` grows while the stack does not, and the stack's *share* of peak memory falls as pivoting gets
harder: the 105 percent of `|L|` measured on 3D is the static figure, and heavy delaying would put it
lower.

**Four peaks from the same simulation.** Serial in natural sibling order, which is what runs today;
serial in the order that minimizes the peak, which is Liu's rule, children sorted by decreasing
`peak(c) - B(c)`, and therefore the payoff `sortForOptimalMultifrontal` would buy; parallel with N
branches live for N of 2, 4, 8, 16; and parallel unbounded as the ceiling. Report each against
`numVal()` as a percentage, and the spread from the optimal serial peak to the parallel one, which is
the number the entry above needs.

**The test set decides the answer more than the code does**, since the ratio is known to be
problem-dependent. A 2D grid (thin fronts, tall tree, ratio expected small), a 3D grid (fat fronts
near the root, expected largest), a band, and something tree-like and very sparse. Sweep n on each,
because whether the ratio grows with n or stays flat is the part that settles the question.

**One modeling choice to make honestly.** The N-branch peak depends on which branches happen to be
live together, so there is no single true number. Compute a deterministic greedy schedule and the
unbounded ceiling, and report both rather than presenting one as the answer.

**A first look, 2026-07-23.** Run from a scratch probe outside the repo rather than the experiment
above, so treat the numbers as indicative and not as the measurement. It builds the symbolic factor
through `SymFactorEngine` and simulates the live-set rule transcribed from `factorStaticMultifrontal`,
comparing only the two serial child orders: forward, which is first-to-last as the code assembles today,
and reverse, which is last-to-first. Sizes in elements.

```
AMD ordering            n       |L|        forward            reverse           rev/fwd
grid2D 20x20           400      4106       1845   ( 44.9%)    1218   ( 29.7%)    0.66
grid2D 40x40          1600     23558       8196   ( 34.8%)    6900   ( 29.3%)    0.84
grid2D 80x80          6400    133261      42735   ( 32.1%)   26425   ( 19.8%)    0.62
grid3D 8x8x8           512     13229      16500   (124.7%)   13328   (100.7%)    0.81
grid3D 12x12x12       1728     91199      97339   (106.7%)   68102   ( 74.7%)    0.70
grid3D 16x16x16       4096    340513     357398   (105.0%)  255433   ( 75.0%)    0.72
band n=2000 bw=5      2000     11995         75   (  0.6%)      75   (  0.6%)    1.00
band n=2000 bw=20     2000     41980       1200   (  2.9%)    1200   (  2.9%)    1.00
```

**The stack is not a rounding error, and on 3D it exceeds the factor.** At 105 to 125 percent of `|L|`
under AMD, peak working memory for a 3D problem is roughly double what the factor alone suggests, and
the stack is the half that is pure scratch. On 2D it is 20 to 45 percent and drifting down as n grows;
on banded problems it is under 3 percent and can be ignored. So the question this entry gates is
answered differently per problem class, which is the reason the test set was specified before the code.

**Reversing the child order is worth 16 to 38 percent on the grids, and costs nothing.** Every grid
case improves, `rev/fwd` between 0.62 and 0.84, and last-to-first beats the first-to-last the code
assembles today. That is not even Liu's order, merely the reverse, so it is a lower bound on what child
ordering could buy and it corroborates `sortForOptimalMultifrontal` being worth porting.

**Under Natural ordering the child order changes nothing at all**, `rev/fwd` exactly 1.000 in every
case. Natural on these graphs produces essentially chain-like forests, so there are no sibling choices
to reorder. The ordering lever exists only once AMD builds genuinely branching trees, which is worth
knowing before anyone measures the lever on the wrong ordering and concludes it is absent.

**A separate lever falls out of the formula, and it is larger than either of the above.** The update
matrix is stored as a full `u * u` square while being symmetric, so keeping only the lower triangle
would halve the stack outright, independent of ordering and of LIFO. The full square is presumably
deliberate, the same rectangle-for-BLAS tradeoff `NumFactorStatic` documents for its own zeros, but at
105 percent of `|L|` on 3D the factor of two is the biggest single number in this entry and deserves
weighing rather than assuming.

**What the probe does not cover**, and what the experiment proper still owes: the parallel dimension
is entirely unmeasured, only the two serial orders were compared; Liu's optimal order was not
computed, so the true serial minimum is lower than the reverse column; and the live-set rule was
transcribed from the driver by reading it rather than instrumented, so the simulation should be
checked against a real run before any of these numbers are relied on.

**Update, 2026-07-23: both halves are ported, and the saving is realized.**
`sortForOptimalMultifrontal` and `labelDepthFirst` are now `ElmForestEngine` methods, called as a
pair behind `setOptimizeMultifrontal` and off by default, exactly the sequence 0.9 runs. The sort
chooses each supernode's child order by Liu's key; the relabeling makes the labels follow that order
as a postorder, which is what the drivers need, since they loop in increasing label order and so a
contribution block is live from its own supernode to its parent in the numbering rather than in the
links. Either alone is inert; the pair is not.

Measured on the same probe, the peak the multifrontal driver actually reaches:

```
                       |L|         option off          option on          saving
grid2D 20x20           4106       1845   ( 44.9%)     1218   ( 29.7%)     34.0%
grid2D 40x40          23558       8196   ( 34.8%)     6900   ( 29.3%)     15.8%
grid2D 80x80         133261      42735   ( 32.1%)    26425   ( 19.8%)     38.2%
grid3D 8x8x8          13229      16500   (124.7%)    13328   (100.7%)     19.2%
grid3D 12x12x12       91199      97339   (106.7%)    68102   ( 74.7%)     30.0%
grid3D 16x16x16      340513     357398   (105.0%)   255433   ( 75.0%)     28.5%
band n=2000 bw=5      11995         75   (  0.6%)       75   (  0.6%)      0.0%
band n=2000 bw=20     41980       1200   (  2.9%)     1200   (  2.9%)      0.0%
```

So the 3D case falls from above `|L|` to three quarters of it, and the bands do not move, having no
sibling choice to make. The multifrontal drivers were not touched.

**What this buys is allocator-agnostic, which is why it was worth doing first.** The saving is a
smaller live set, and a smaller live set is smaller under every arrangement that might follow: an
arena, per-thread arenas, or the per-supernode slots we have. Nothing decided later about parallelism
or about storage can un-earn it, which is what makes it safe to bank before those questions are
settled.

**But Liu's rule minimizes the serial peak, and that should not be quietly assumed to carry over.**
The recursion it optimizes assumes a parent's children run one after another, so that a finished
child's block sits waiting while its later siblings run. Under concurrency children may run at the
same time instead, and the quantity being minimized is no longer the same quantity. The ordering
installed here is a sound starting point and certainly better than an arbitrary one, but it is not
proven optimal for the parallel case, and the child ordering is a free parameter there too. Whether
it wants revisiting is part of the unmeasured parallel dimension, not something this measurement
settles.

**One consequence of the option to keep in view, unrelated to the peak.** The forest is shared by all
three traversals, so reordering and relabeling change the delayed-column prepend order for dynamic
LDL, which fixes a front's column order and therefore the pivot sequence. Turning the option on will
give a different and equally valid dynamic factorization with different delay and pivot counts, which
is why the pinned tier-1 counts hold only with it off. That is the same structural-not-commutative
property the stack entry above records for the delay prepend, and it is the reason the option is off
by default rather than simply always on for multifrontal.

**What remains for the experiment proper**, now that the serial ordering question is settled: the
parallel dimension is still entirely unmeasured, and the triangular-storage question below is
untouched and larger than everything measured here.

**Reversing the assembly order changes nothing, measured 2026-07-23.** Assembling a parent's children last
to first rather than first to last is the order a true LIFO stack pops in, so it is the obvious cheap
alternative to any arena work. It was tried in an instrumented scratch copy of
`factorStaticMultifrontal`, on the six grids above, with the child-ordering option both off and on.
The peak is identical in all twelve configurations, to the element.

It cannot be otherwise, and the reason is the shape of the driver rather than anything numerical. The
peak is reached at `stack[kk].allocate(...)`, the instant kk's own block joins its children's, all of
which are still live; reversing the assembly merely permutes frees that happen after that instant. The
assembly loop discards the same set either way.

**The same shape is why the pattern is not stack-realizable, which is the more useful half of the
result.** Last to first is a genuine LIFO order for the children, but kk's block was allocated
*before* the assembly and therefore sits above them, so assembling in reverse pops the children in LIFO
order among themselves while a stack pointer still cannot free beneath its own top. A true stack
would need `stack[kk]` allocated after all children are assembled, which this design forbids:
`assembleUpdateMatrix` needs kk's block already present as the destination for the child entries
falling in kk's update rows. So the LIFO *order* is available and buys nothing, and the LIFO
*discipline* is not reachable without changing what the assembly writes into.

**The instrumentation also validates the simulation.** These numbers came from inside the real
driver, not from the model, and they match the simulated table above to the element, with the option
off and on. That closes the caveat recorded above, that the live-set rule had been transcribed by
reading the driver rather than instrumented.

**Implications for parallelism, which are the reason this negative result is worth keeping.**

The moment that sets the peak is a join: a parent allocates while every child's block is still
resident, so all of a supernode's children must be complete and simultaneously in memory before it
can begin. That is inherent to the current shape, and under branch concurrency it happens in N places
at once, which is the multiplier the stack entry above describes.

That assembly order does not affect the peak is a freedom for a scheduler rather than a curiosity. A
parent may assemble its children in whatever order they finish, which is what branch parallelism will
hand it, without paying anything in memory and without needing a canonical order. Note the contrast
with the delay prepend for dynamic LDL, which is structural and not free to reorder; the static assembly
genuinely is.

The assemblies of one parent are done by whichever task owns that parent, so branch parallelism races on
nothing: children write their own blocks, the parent reads them and writes its own. Parallelizing the
assemblies *within* one parent is the case that would race, several tasks accumulating into one block, and
would need per-thread partial blocks or atomics. That is node parallelism, a different axis, and not
what the forest work is about.

Finally, the shape suggests an option that only becomes visible once the allocate and free pattern is
written down. A parent could allocate its block early and assemble each child as it completes, rather
than waiting for all of them, freeing each child's block at once. Its index set is known from the
symbolic factorization, so nothing prevents it. That trades kk's own block being live for longer
against its children never accumulating, giving `B(kk) + max over children of their subtree peak`
where the present shape gives `sum of all children's blocks + B(kk)`. Which is smaller depends on the
node: a wide parent with many small children favors assembly on arrival, a narrow one with a single
deep child favors the present shape. It is a genuine third option alongside the arena and the
ordering, unmeasured, and it interacts with parallelism directly, since assembly on arrival is exactly
what a task that wakes on each child's completion would naturally do.

### Forest (tree) parallelism on Apple Silicon

ARCHITECTURE's parallelism note ends on the observation that Oblio gets node parallelism for free
through Accelerate, but that tree parallelism, factoring independent forest branches on separate
cores, is Oblio's to build. This is that item. The node side needs nothing; this is the tree side.

The how is worked out in DESIGN_DECISIONS (2026-07-22). The split that matters is portable versus
architecture-targeted, and portability wins by default, since Oblio runs anywhere a BLAS does. Portable
and recommended first is OpenMP: a compiler feature rather than a dependency, its `task depend` pragmas
map onto the tree DAG, and they degrade to valid serial code when OpenMP is off, so the any-CPU
property is untouched (the one snag is Apple Clang shipping no runtime, so macOS needs libomp). The C++
task-graph runtimes (Taskflow, oneTBB) are the portable, nicer-ergonomics option at the cost of a
dependency. Architecture-targeted stacks (GCD on Apple, a vendor equivalent elsewhere) come only when
platform tuning like P-core placement is worth giving up portability. Not MPI, which is portable but
process-based and memory-heavy on a single shared-memory node, where threads are the right model; not
pthreads, which means hand-rolling the work-stealing DAG; not the GPU, which is the wrong shape for
many small irregular tasks.

**Measured 2026-07-23 on alpamayo (M4); see `experiments/openmp/` for the tables and method.** Two
paragraphs of the above were written from documentation and are now corrected by measurement.

The cap is a no-op here, and worse, it may be unimplementable. `VECLIB_MAXIMUM_THREADS=1` changes
nothing at any order, and neither does `OMP_NUM_THREADS=1`. The second is the informative one: if
Accelerate threaded *through OpenMP*, capping OpenMP would have capped it, so either the library does
not thread at these sizes or it threads through something honoring neither variable, with Grand
Central Dispatch the obvious candidate. If the latter, then when fronts eventually run as concurrent
tasks each front's BLAS call could spawn its own work with no lever to stop it, and the instruction
this entry gives, cap Accelerate to one thread per front, could not be carried out at all. One lever
remains untried, Accelerate's own `BLASSetThreading` called from code rather than set in the
environment; checking Apple's documentation for it is the cheapest next step and would settle the
threading question at the same time, since a GF/s column that moves at all proves the threading was
there.

**The split between library threading and own threading is a machine-dependent decision, not a
design to commit to.** A solver that assumes a sequential BLAS and threads the dense work itself is
tuned for x86, where each core has its own vector units, and wasted here, where the threads would
funnel into one shared matrix unit. A solver that assumes a threaded BLAS is tuned for here and
leaves x86 cores idle. This is why the MUMPS work reaches for a performance model rather than a fixed
rule, and it argues against baking either assumption into Oblio.

The window is far narrower than "the leaves" suggests. Running two whole `dgemm`s side by side rather
than back to back gains about 1.5x at order 32, 1.14x at 64, and nothing from 128 up. A hand-written
triple loop under identical scheduling gains 1.9x across that whole range, which is the control
proving the cores are real and the shortfall belongs to the kernel. So the dense side has a payoff
window of roughly order 32 to 64 and no more. The mechanism is not fully resolved: either Accelerate
already spreads a call across cores, so the serial baseline was never serial, or the cores are free
and all funnel into a matrix unit shared per cluster that one call already saturates. Both may hold.
For this decision the disjunction suffices, since either way a second core adds nothing above order 64.

**Where the opportunity actually is.** Dense kernels run at about 448 GF/s against the hand loop's 21,
so they are already at the hardware ceiling and cannot get faster, while everything else in a
multifrontal pass is our own code: the assembly and `gblToLcl` scatter, the update stack, and for
dynamic LDL the threshold pivot search, which is a scan with no BLAS in it. That work uses per-core
resources, which is precisely what the hand loop's steady 1.9x says two cores deliver. As the kernels
stay pinned and problems grow, the structural share only rises. **The number that decides this item is
therefore the structural-to-dense time ratio inside `factorMultifrontal`, and nothing else needs
measuring first.** Timers around the BLAS calls against everything else would answer it; bucketing by
supernode size afterwards sharpens it into the front-size histogram the same instrumentation gives.

**Why multifrontal and not the other two.** MF's dependency set is exactly the parent-child edge: a
front is ready when its children are done, which is a dependency the tree already encodes and which
`task depend` expresses natively. LL gathers from a scattered set of descendants and RL pushes to a
scattered set of ancestors, so both have a dependency graph much denser than the tree even though
their traversal order is tree-consistent. The memory argument is stronger still: MF accumulates into
a private contribution block per supernode and the branches meet only at the assembly into the
shared parent, whereas concurrent branches under RL would accumulate into the same ancestor columns
and need locks or per-thread privatization. MF is the traversal whose parallel form needs no
synchronization its serial form did not already have.

**What `factorStaticMultifrontal` would need.** Four changes, none large, listed so the work can be
picked up cold.

The driver loop has the wrong shape. It is a flat loop over supernodes in increasing order, relying
on the postorder numbering to have children done before parents. A `#pragma omp for` over it is
simply incorrect, since it expresses that children come earlier and not that `kk` waits for them.
Either a recursive walk that spawns a task per child and waits before assembling, or `task depend`
keyed on the stack slots. The recursive form is natural and needs no `depend` clauses.

`gblToLcl` is a genuine race. One `std::vector<std::int32_t>` sized to the matrix order, filled by
`setGlobalToLocal` and emptied by `clearGlobalToLocal` around every supernode. It must become
per-thread, at `O(n)` per thread, which is cheap at these thread counts.

The early exits cannot stay. Two `return false` paths sit inside the loop, for `assembleFromA`
failure and for non-positive-definite Cholesky, and a return cannot leave an OpenMP structured block.
They become an atomic flag, read by sibling tasks if they should stop rather than finish.

Task granularity needs a cutoff, since spawning a task per order-32 front spends more on scheduling
than on arithmetic. Relatedly, `updateStaticUpdateMatrix` allocates a local `upper` buffer per call,
which at leaf granularity becomes allocator traffic from every thread at once and may want a
per-thread buffer.

**One structure is already right, and it is the one recorded above as a wart.** `stack` is a
`std::vector<UpdateMatrix<Val>>` with one slot per supernode, written by `kk` and read by its parent,
so concurrent branches touch disjoint slots. The preceding entry records this non-LIFO shape as a
departure from 0.9, conditional on whether the scratch ever proves to be a binding constraint; for
parallelism it is exactly right, since a *single global* LIFO
would be shared mutable state that every concurrent branch would contend for, and would stop being a
LIFO at all once branches interleaved. Per-branch arenas are the form that keeps both properties.
Anyone reintroducing stack discipline should read that entry and this one together.

All of the above is the static path. Dynamic multifrontal is harder, because a delayed column
migrates from a child into its parent's index set, which is real cross-front coupling rather than the
clean parent-child handoff the static case has.

**Order of attack, by effort against certainty.** Take the free parallelism first, which means
leaving dense kernels to whatever the BLAS does, since the measurement says there is nothing left to
take there on this hardware. Then do the work no library can do, the structural loops, which pays on
any machine because it uses per-core resources rather than a shared unit and does not depend on which
BLAS is linked. Branch parallelism comes with that rather than before it: it is worth about 1.5x on
dense fronts below order 64 and nothing above, but it needs the same restructuring the structural
work needs, a driver expressed as a tree walk with per-thread scratch, so the two arrive together.
Getting smarter than that, a performance model choosing per region of the tree, is a later and much
larger project.

#### Memory under concurrency, and what the child ordering does and does not do

**The peak formula has one form, and concurrency changes one term of it.** With the notation of the
experiment entry above, and I and J ranging over K's children,

```
peak(K) = max( TERM1, sum_J U(J) + U(K) )

TERM1 = max_J [ sum_{I<J} U(I) + peak(J) ]      children sequenced
      = sum_J peak(J)                           children concurrent
```

Term 2 is untouched: it was already order-invariant, and it was already the join where every child is
simultaneously resident. Term 1 goes from a maximum to a sum, which is the entire memory cost of
concurrency, and `sum_{I<J}` disappears with it. **So at a node whose children run concurrently, the
child ordering buys exactly nothing**, there being no "before" for it to act on. The two levers are
disjoint per node: ordering matters where children are sequenced, and only there.

Note this is per node rather than per region. A node whose children are sequenced evaluates the
maximum form and a node whose children run together evaluates the sum form, and mixtures sit between,
three children with the first two concurrent giving `max( peak(J1) + peak(J2), U(J1) + U(J2) +
peak(J3) )`. `peak(K)` still summarizes everything below K in one number either way, so the recursion
still composes and a parent still needs nothing else from inside a child.

**The trade is the ordinary one, space for time, and it should be stated plainly.** Serial execution
gives the minimum aggregate peak, and Liu's rule over the whole forest attains it exactly. Running
branches concurrently raises the aggregate, because `sum_J peak(J)` exceeds `max_J peak(J)`. The extra
memory is the concurrency; there is no arrangement in which parallel execution matches the serial
minimum. What is bought with it is wall-clock time.

**A correction worth recording, because it was made in discussion and is easy to make again.** The
serial minimum is not a ceiling in any enforcing sense. Liu's rule is an ordering and enforces
nothing: nothing prevents two subtrees from running concurrently and exceeding the serial peak, and
if they do, the program simply uses more memory. Treating the computed minimum as a constraint that
forbids concurrency is a category error, and it leads to elaborate reasoning about which parallelism
is "allowed" when in fact none of it is gated. Any real bound on memory would have to come from
admission control at run time, refusing to start a branch when the budget is short, and that is a
mechanism nobody has built here.

**What the per-subtree optimization actually gives, then, is honesty and information.** It minimizes
each subtree's own peak, which is free and never harmful, and along the way it computes `peak(K)` for
every supernode: what the run would cost sequentially, and the per-subtree figures a memory-aware
scheduler would need. Without it one can still parallelize freely and simply never know what is being
spent. With it, the numbers exist. It blocks nothing and informs everything, which is the right way
to hold it.

**The numbers are currently computed and thrown away.** `sortForOptimalMultifrontal` fills a local
`maximumStorage` array over all supernodes and lets it die at the end of the call. Keeping it on the
forest would give a scheduler its inputs directly. The live total at any instant is not only the
active subtrees but `sum of peak(J) over active subtrees`, plus `sum of U(I) over siblings finished
and not yet assembled`, plus whatever ancestors hold, all of which is arithmetic on that same array.

**Ordering and schedule are one joint decision, not two independent ones, and this is the practical
bind.** Liu's key is `peak(J) - U(J)`, and `peak(J)` depends on whether J's own children were run
concurrently, so the sort cannot be computed before the schedule is known, and a schedule cannot be
costed before the sort. The schedule must be fixed first, and then the recursion is well defined and
the sort correct for it. Our current code does exactly this silently: it assumes everything is
sequenced, evaluates the maximum form throughout, and sorts accordingly, which is right for the serial
driver and is what was measured. Introduce concurrency without revisiting it and the ordering is not
wrong, it is simply optimal for an execution that is not the one happening.

**Which leaves a fork that any future work here has to pick.** Either a static plan, deciding in
advance which subtrees run in parallel, evaluating the sum form at those nodes and the maximum form
elsewhere, and sorting under those numbers, which is coherent and exact but needs a scheduler that
honors the plan; or dynamic scheduling, keeping the serial ordering as a heuristic and treating
memory as a run-time constraint, which predicts nothing and needs nothing predicted. The MUMPS
performance model is the first, work-stealing runtimes are the second, and the choice determines
whether the ordering numbers above are a plan or a heuristic.

Requires multifrontal, whose subtrees are already independent tasks with a per-subtree stack. Trigger:
when the structural share of `factorMultifrontal` is measured and found large, which supersedes the
earlier trigger of serial leaf-front time, since the dense leaf window is now known to be narrow.

Reference: J.-Y. L'Excellent and W. M. Sid-Lakhdar, "A study of shared-memory parallelism in a
multifrontal solver", *Parallel Computing*, 2014, already cited in ARCHITECTURE. It measures threaded
BLAS and solver OpenMP as separate contributors rather than lumping them, proposes the performance
model mentioned above, and finds the payoff turns on the ratio of large fronts to small ones, which
is the same front-size question named here.

## Structure

Two extractions are available in the dynamic factorization code, recorded here in enough detail to
be picked up cold. Both are refactors of working, tested code: 143 assertions cover this path,
including tier 1 and tier 2 matrices that exercise both pivot selections and both traversals, so
either can be done with a real safety net rather than by eye.

Measured 2026-07-19, code lines with comments stripped: `factorDynamicSupernode` 230 against 0.9's
313, `updateDynamicUpdateBlock` 49 against 61. The port is already about a quarter shorter than the
reference, almost all of it from dropping the error-macro scaffolding and the raw-pointer preamble
rather than from restructuring. What follows is the restructuring that has not been done.

### Merge the duplicated pivot bodies: DONE, 2026-07-19

The eliminations are now `applyPivot1x1` and `applyPivot2x2`, shared by both selection loops, and
the 2x2 block is read once in `readPivotBlock2x2`. 230 code lines became 204 across three functions,
under 147 passing assertions throughout.

The non-goal held: **the selection loops stayed separate.** They are two algorithms rather than one
with flags, differing in three ways, no forced 1x1 in pass 2, a partner scan bounded by the front,
and a Bunch-Kaufman determinant test in place of `max1 == max2`. Merging those behind parameters
would save lines and cost the reader the ability to see what each pass does. That remains the
recommendation if anyone revisits it.

Two things surfaced that the duplication had hidden, which is the usual argument for removing
duplication and is worth recording as evidence for it. **Pass 1 never reads the 2x2 block**: it
decides on magnitudes alone, and the compiler pointed this out with an unused-variable warning the
moment one shared body served both. And `readPivotBlock2x2` is the single place the symmetry of D is
decided, `d12 = d21` being the symmetric statement, which makes complex `LDL^H` a change in one
function rather than four.

### Extract the front expansion shared by all three dynamic drivers

Not previously recorded. `factorDynamicLeftLooking`, `factorDynamicRightLooking`, and now
`factorDynamicMultifrontal` expand a parent's index set with about 16 identical lines: extend the
index array by the children's total, shift the existing indices right by that amount, then prepend
each child's delayed globals in sibling order. Verified identical apart from local variable names.
The third copy landed with dynamic multifrontal, which raises the payoff and the case for doing this.

They differ only in the verb that closes the block, `resetVal` for left-looking and multifrontal
(both discard, nothing being in the front yet) and `expandVal` for right-looking (which preserves the
values already accumulated). That is the one place the traversals genuinely part company on delayed
columns.

So the extraction is clean and the signature is narrow, something like

```
std::int32_t expandFrontIndexSet(NumFactorDynamic<Val>& nf, std::int32_t parent,
                                 const std::vector<std::int32_t>& firstChild,
                                 const std::vector<std::int32_t>& nextSibling) const;
```

returning the number of delayed columns assembled, with each driver calling its own val verb
afterwards (`resetVal` or `expandVal`). That also has a documentation benefit beyond the line count:
it puts the drivers' real difference on one visible line each, instead of at the bottom of a block
that otherwise reads identically in all three. (Name uses the current `expand` verb, matching
`expandNodeIdx` / `expandVal`; an earlier draft of this entry said `grow`, before that rename.)

Lower risk and smaller than the pivot merge, and independent of it. Either can be done first.

### The max1 == 0 swap branch: leave duplicated, do not extract

`factorDynamicSupernode` has the same five-line `max1 == 0` branch in both pivot passes (mark pivot
found, `swap` if `j_ != k1_`, set `mPivotType` to 1, advance `j_`, break). It reads as a candidate
for extraction alongside the front expansion, but it is not one, and this note records the decision
so it is not re-proposed. The five lines are welded to two different scan structures: pass 1 does a
two-loop partner scan, pass 2 a one-loop scan, and the branch's `break` and `++j_` act on *that
pass's* loop and cursor. Lifting it would mean either threading "what to do next" signals back out
to the caller, or extracting the whole pivot-accept region, which would bury the
pass-1-versus-pass-2 distinction that is the actual content there. Five visible lines repeated is
the cheaper reading than a helper that has to route control flow. Unlike the front expansion, the
duplication here is the lesser evil.

### Left-looking A-assembly ordering: whether to align its loop skeleton with right-looking

Open design question, deferred for thought, not a bug. The two dynamic drivers currently order the
per-supernode work differently, and the difference is the storage-level shadow of pull-versus-push.
Three orderings are on the table.

Right-looking, as it stands. A separate prepass assembles A into every front first; then the main
loop, per supernode jj, expands jj from its children (`expandVal`, which *preserves* what is already
there), contracts children and assembles their delays, factors jj, and pushes jj's updates up into
ancestors. A is in place before the loop because a push lands on an ancestor not yet reached, and
`assembleFromA` assigns, so it cannot run after a push without overwriting it.

Left-looking, as it stands. No prepass. Per supernode kk, the loop expands kk from its children
(`resetVal`, which *discards*: nothing has been written into kk yet, so the widened block is simply
re-zeroed), assembles A into the now-expanded front at the delayed-column offset, contracts children
and assembles their delays, pulls updates from descendants, and factors kk last. Assembling A after
the expand is what lets the expand discard, and discarding is cheaper than the shift `expandVal`
does: A is placed directly at its final offset rather than written low and moved.

Left-looking, the proposed aligned form. Move `assembleFromA` to the first line of the `kk` loop, so
every driver reads `assemble A, expand, delays, updates-or-factor` in one skeleton. This is not a
pure reorder: assembling A before the expand means the expand can no longer discard, so `resetVal`
must become `expandVal` (preserve and shift), exactly right-looking's verb. The cost is real if
small: left-looking would assemble A into the unexpanded front and then shift those values
during the expand, rather than placing them once at the final offset. It also erases a distinction
the code and `docs/ARCHITECTURE.md` currently draw deliberately: left-looking fills a front from
scratch each turn (`resetVal`), right-looking carries a front's accumulated state through expansions
(`expandVal`). Note this does not remove the prepass asymmetry: left-looking would still have no
prepass, since a front's shape still is not known until its descendants are factored. It removes
only the `resetVal`-versus-`expandVal` asymmetry.

The trade is uniform loop skeleton against a genuine (if small) efficiency and a meaningful
distinction. If taken, the change must switch `resetVal` to `expandVal` in left-looking, re-verify
that A lands identically under assemble-then-shift versus expand-then-assemble-at-offset, and
rewrite the `resetVal`/`expandVal` header comments in `NumFactorDynamic.h` (currently at the two
verbs) and the lazy-assembly paragraph in `ARCHITECTURE.md`, which would no longer be accurate.

### Align, and perhaps unify, the two update kernels

`updateStaticUpdateBlock` and `updateDynamicUpdateBlock` compute the same thing, `block -= L21 D L21^H`,
over the same geometry. Both walk jj's front (`f` columns) against the update rows the ancestor
receives (`L21`, starting at `offset`), form an upper factor `U := D L21^H` into an `f` by `width`
scratch, then close with `gemmLower` for the symmetric square and a plain `gemm` for the rectangle
below. The closing two multiplies are already identical. The geometry matches once the names are
read correctly: static's `f` is the whole front and its stride is `frontSize + updateSize`;
dynamic's `f` is the post-factor front (reduced by the delayed columns, which no longer drive
updates) and its stride is `frontSize + delaySize + updateSize`, the one extra term being the
delayed columns' rows, which a delayed column keeps even though it left as a column.

Near term, align the coding style so the two read as twins. The prologue is where they diverge
cosmetically: static pulls `frontSize`/`numNodeIdx` as named `size_t` locals through accessors and
casts once to `int`; dynamic goes straight to `int` off raw members (`nf.mFrontSize[jj]`,
`nf.mDelaySize[jj]`, `nf.mUpdateSize[jj]`) and is inconsistent with itself, using accessors for
`val`/`nodeIdx` but raw members for the sizes. Bring dynamic to static's form: named `size_t` locals
via accessors, then the cast. Keep the differences that are essential, not cosmetic: dynamic's
`f == 0` early return (every column delayed, no pivot to update with), and its `withHermitian`/`idx`
locals, which the pivot walk needs. Leave the kernel free of the pre/post/expand vocabulary; inside
it jj is simply whatever the caller passed, so bare `frontSize` matches static and does not leak the
expand-block naming.

The one substantive difference is forming `U`. Static calls `formUpper`, which walks D as a pure
diagonal, one `d[dp]` per column (`BlasLapack.cpp`), correct because static's D is diagonal
(identity for Cholesky, a plain diagonal for static LDL). Dynamic inlines its own form-upper loop
because its D is block diagonal: `mPivotType` marks each column 1x1 or 2x2, the 1x1 branch mirrors
`formUpper` line for line, and the 2x2 branch couples two columns through four D entries, which
`formUpper` cannot express. So dynamic's loop is `formUpper` generalized to block-diagonal D.

Longer term, decide whether to truly unify. A pivot-type-aware `formUpper` that handled 2x2 blocks
would let dynamic call out instead of inlining, after which the two kernels would differ only in
where `f` and the stride come from, and could plausibly become one templated kernel. The cost is
coupling the BLAS-adjacent layer to pivot types, and losing the inline pivot walk, which is readable
where it sits. Align first, then take this decision with the symmetric versions side by side.

### Pass the factor to the assemble functions: done for assembleUpdateBlock, open for assembleFromA

**`assembleUpdateBlock` now takes `(jjKkUpdateBlock, nf, kk, gblToLcl)`**, done 2026-07-23, so two of
the three assemble functions resolve their own destination and only `assembleFromA` still takes a raw
`(std::size_t numNodeIdx, Val* block)` pair. The change was mechanical and the reasoning below held
up, so it is kept here as the record for finishing the job on `assembleFromA`.

It works because at every call site the destination `kk` is still unfactored when the assemble runs
(the factor call comes later, or `kk` is a not-yet-reached ancestor in the push case), so
`delaySize(kk) == 0` and the destination derives cleanly: `kkVal = nf.val(kk)` and
`kkNumNodeIdx = nf.frontSize(kk) + nf.updateSize(kk)`, which is exactly what the callers passed. All
four call sites were checked to be deriving those two values that way already.

The appeal is readability at the call site and consistency with `assembleDelay`, `assembleUpdateMatrix`
and the factor and update functions, which all take `nf`. There is **no execution-time cost**: at
`-O3` the accessors inline to the same loads the raw-pointer form used, and re-deriving
`kkNumNodeIdx` per call is an add the optimizer absorbs.

The costs, confirmed rather than predicted. Passing `nf` turned `template<class Val>` into
`template<class Val, class Factor>`, doubling the instantiations, since the static drivers run these
into both `NumFactorStatic` and `NumFactorDynamic`. The call-site cleanup was smaller than it looks,
as expected: the drivers bracket the scatter with `setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, ...)`
and `clearGlobalToLocal(...)`, so `kkNumNodeIdx` and `kkNodeIdx` stay; only the `kkVal` locals went,
and only in the two right-looking drivers, where the compiler flagged them unused. And it gives up
the storage-blind raw-block shape, which is the same principle that keeps the numeric kernels blind
to static-versus-dynamic storage (the storage-options experiment, cited in `ARCHITECTURE.md`).

One incidental gain worth noting: deriving `nf.val(kk)` inside the function re-reads the pointer per
call rather than relying on a hoisted one, which retires the comment in the dynamic left-looking
driver explaining why the hoisted `kkVal` stayed valid across its loop.

### Measure left-looking against right-looking

The two traversals compute the same factor by the same kernels and differ only in bookkeeping, so
the comparison is clean, and neither side is obviously cheaper. **The reasoning lives in "The life of
an update" in docs/ARCHITECTURE.md**, what the per-pair costs are, why the unit is pairs rather than
supernodes, and what could be done about left-looking's queues if it loses. That belongs there
rather than here: it describes how the thing works, and stays true whether or not this item is ever
picked up.

What is owed here is only the measurement. It wants matrices large enough to separate an allocation
from a linear sweep, and it should report per-pair costs rather than totals, so the two sides are
comparable. If left-looking measures badly, try the flat-array queues described in ARCHITECTURE
before concluding anything about the traversal itself.

Worth doing alongside any other timing work rather than on its own, and worth doing before either
traversal is chosen as a default for anything.

### Amalgamation tie-break: front size versus list position

`compressThreshold` in ElmForestEngine breaks a fill tie between two children by taking the one with
the larger front, and only then by list position. The larger-front rule is a greedy heuristic,
"absorb the wide children before the parent's front grows and prices them out", and section 4.5 of
docs/sparse_factorization.md explains it. It is locally tempting but has no global guarantee: each
absorption widens the shared front, so a locally good pick changes the state every later pick is
priced against and can foreclose a better sequence. The objective it optimizes, block quality, is
itself machine-dependent and fuzzy, so there is no clean optimum being approximated.

That raises a simpler alternative: drop the front-size level and break fill ties on list position
alone. It loses an unproven heuristic and gains a shorter, more honest rule, one that does not imply
a rigor the greedy algorithm lacks. It does not remove arbitrariness, only moves it from front size
to list order, but both are conventions, not theorems.

Which is better is empirical, and the two variants differ by one line, so the experiment is clean:
run both over a matrix suite, holding everything else fixed, and measure. Worth reporting how often
the tie-break even fires (it needs an exact fill tie *and* differing front sizes, likely rare), how
often the partition then differs, and, the thing that actually matters, the factor flops or timing
of the resulting blocks on the target machine. If the deltas are in the noise, positional-only wins
on simplicity; if front-size-first measurably helps, it earns its complexity. Fits the experiments/
prototype-first pattern: an isolated harness that runs both and reports the comparison, kept out of
the main code until the conclusion is drawn.

## Testing

### Examples are built but never run

`make` compiles everything in `examples/`, and nothing ever executes it or checks what it prints. So
an example can go stale silently, and one did: `examples/pipeline.cpp` fixed its storage at
`NumFactorStatic` and therefore reported every dynamic LDL cell as "not implemented" for as long as
it took someone to run it by hand and notice. The code was right; the example was lying about it.

An example is documentation that compiles, which is most of its value, but compiling is not the
claim it makes. The cheap fix is to run each example under `make test` and require exit status zero,
which catches a crash but not a wrong number. The honest fix is for examples with deterministic
output to have that output checked. Worth deciding which, rather than drifting into the first
because it is easier.

### Assert the LDL perturbation branch

The branch runs and nothing checks it, which PORTING_LEDGER now records with measurements. The
recipe is known: `StaticLDLT` on `bandIndefinite(40, 3, 0.3, 7)` under Natural ordering perturbs one
pivot and lands at 1.4e-03.

Two assertions, and the second is the one with content: that `numPerturbations()` is nonzero, and
that the factor reconstructs a matrix differing from `A` by about the perturbation. A residual
assertion alone would just record that the answer is worse, which says nothing about whether the
right thing was done.

Worth pairing with a note in the specification that a poor residual from a static factorization on
an indefinite matrix is expected behavior rather than a defect. That confusion has now cost time
twice in one session.

### A specification the suite can be regenerated from

The suite is currently defined by its own source: what is tested is whatever the test files happen
to assert. That is enough to catch regressions and not enough to notice a gap, and it makes the
combinations we deliberately do not support indistinguishable from the ones we forgot.

What we want is a document that states, per combination of ordering, factorization, traversal and
scalar type, whether it is supported and what is asserted about it, plus the catalog of test
matrices and what each is for. Tests are then written from the spec rather than the spec being
recovered from the tests, and a combination that is not yet implemented is a cell that asserts a
refusal rather than a hole.

Being worked on now; this entry stands until the spec exists and the suite matches it.
