# Decisions

Durable record of structural choices, newest first. Each entry: date, decision,
why. This is the file to open after a gap to reconstruct the project's shape.

---

## Why this file exists

**In a solver, the bottleneck was never writing the code.** A Cholesky is a day's work. A
supernodal Cholesky is a week's. What takes years is deciding *how the pieces should fit*: where
the factor's values live, whether the engines hold state or the data does, how symmetry is
expressed, what the API promises, which phase owns which fact. Those are the decisions Oblio was
made of, and they were the hard part in the late 1990s.

**They are expensive because they propagate.** Choose the storage and you have chosen the API;
choose the API and you have chosen what every engine can see; and by the time a mistake shows up
it is behind ten thousand lines that assume it. A wrong turn is not a wrong function, it is a
redesign. So each decision carries risk, so you deliberate, so they queue, so the project slows.
That is the real reason a library like this takes a decade rather than a year.

**And the tooling of the time made it worse.** C++98 with shaky compiler support: templates you
could not trust, no move semantics, containers that might or might not be there, and no cheap way
to *try* an abstraction and measure what it cost. The design space was constrained and the
instruments for exploring it were bad. So you chose carefully, once, and lived with it.

**What has changed is not that decisions became easier. It is that exploring one became
disposable.**

A design question can now be turned into a measurement in an hour (`experiments/storage-options`
answered flat-versus-VV with a number and a symbol table). A ported kernel can be tested standalone
before it touches anything (0.9's LDL was reconstructing `A` before the engine had ever seen it). A
clever solution can be built, looked at, called a hack, and rewritten before lunch (the
identity-offsets array in `SymFactorEngine`, which existed for about twenty minutes). A claim can
be made, and superseded the same day when it turns out to have been true only of the case that
prompted it ("symmetry is determined, not chosen", which LDL falsified within hours).

**So the thing to optimize for is not choosing correctly the first time. It is making a wrong
choice cheap to discover and cheap to reverse.**

That reframes what infrastructure is *for*:

- **Per-phase oracles.** Each phase is checkable against something that shares no code with it, so
  a mistake is localized rather than inferred.
- **An end-to-end residual.** `||Ax - b|| / ||b||` in one number, across six phases, tells you the
  pipeline is *consistent*, which no per-phase test can.
- **`experiments/`.** A design argument that can be settled by measurement should be, and the
  measurement should be kept, so nobody re-litigates it from memory.
- **This file.** A decision that was reversed is more instructive than one that was right, and a
  reversal that is not written down will simply be made again. Entries here supersede rather than
  overwrite, and say what was wrong and why it looked right.
- **A willingness to throw work away.** The cost of a wrong turn is now measured in minutes, which
  means the correct response to "that is a hack" is to do it again, not to defend it.

**And the human's role changes rather than shrinks.** Code is fast; *judgment* is not. Which claim
is load-bearing and which is decoration, whether an argument is sound or merely plausible, whether
an abstraction is elegant or merely clever, these do not get cheaper, and they are what decides
whether the fast part was worth doing. Nearly every real correction in this project came from
someone looking at a plausible answer and saying "no, that is not quite it".

The best example is the smallest. Cholesky spent years as an open question here, "in real it is
`CC^T`; in complex, can I have both `CC^T` and `CC^H`?", and it would not close, because it was
the wrong question. `CC^T` in complex does not exist: positive definiteness requires `x* A x` to
be real, which requires Hermitian. Once that is on the table, the design collapses to a sentence:
**Cholesky is `CC^H`, always, and in real that *is* `CC^T`.** No option, no flag, no forbidden
combination to reject. The answer was not hard. Asking the right question was.

---

## 2026-07-15, An in-header throwing constructor slowed a hot loop in the same translation unit

A concrete, measured finding from the storage-options experiment, worth recording because it is
counterintuitive and it informs where numfact puts its code.

**Symptom.** After the matrix constructors gained the nnz/dimension guard (a `throw`), the static
direct `multiply()` slowed from ~1.03x to ~1.20x of the hand-written baseline on the M4 (Accelerate,
`g++ -O3`). Consistent across many runs, so not noise. The dynamic and baseline rows were unaffected.

**Diagnosis.** The multiply loop's own source had not changed, and the constructor is never called
anywhere near it. What changed is that the constructor, defined *inline in the header*, became a
potentially-throwing body visible in the translation unit that compiles
`multiply<SparseMatrixStatic>` (`MultiplyEngine.cpp`, which includes the matrix header for the
explicit instantiations). The exception path it introduced perturbed the optimizer's treatment of
the hot loop in that unit, even though the loop neither throws nor constructs anything.

**Confirmation, both directions.** Commenting out the `throw` restored 1.03x. Moving the whole
constructor into a `.cpp` (so the throwing body is no longer in the header, hence no longer in
`MultiplyEngine.cpp`'s unit) also restored 1.03x, with the guard intact. Object-file check: with the
constructor in its own `.cpp`, `MultiplyEngine.o` contains no `length_error` or `__cxa_throw`
machinery at all; the exception path is gone from that unit.

**Fix.** Both matrix constructors moved from their headers into `SparseMatrixStatic.cpp` /
`SparseMatrixDynamic.cpp`; the headers declare only. Accessors and the non-throwing mutators
(`setValues`, `setColumn`, which return false rather than throw) stay inline, since only the
constructor carries an exception path. This mirrors the main-code `SparseMatrix`, whose constructor
is already in `SparseMatrix.cpp`, so the experiment now matches the real tree here.

**The lesson, and why it is a third reason.** The usual "small bodies in headers, large bodies in the
`.cpp`" rule (CLAUDE.md, and the explicit-instantiation entry) is argued from inlining and build
time. This adds a distinct, optimization-quality reason: an exception path anywhere in a translation
unit can degrade the codegen of unrelated hot code in that same unit, at a cost that is not small
(17% here). So a throwing body is "heavy" for header purposes even when it is textually short. The
rule for numfact follows directly: keep throwing and heavy bodies out of the translation units that
compile the numeric kernels, since the factorization and solve loops are exactly the hot code this
would silently tax.

---

## 2026-07-15, The nnz cap is an ordering constraint, not a representability one; A is capped, L is not

A follow-on to the index-types entry (2026-07-09), sharpening one thing it left blurred: it speaks
of "the ~2.1-billion index cap" as if a single ceiling covered both indices and nnz. There are two
ceilings, with different origins, equal today only by coincidence.

**Two separate ceilings.** The *dimension/index* ceiling holds every id to `2^31 - 1`, because ids
are `std::int32_t`. But note the dimension `n` is capped one below what index representability alone
would allow, and for a distinct reason (see the next paragraph). The *nnz(A)*
ceiling is not that. nnz is a count, stored as `std::size_t` everywhere (`colPtr`, block offsets,
every position), so nothing internal caps it at `2^31`. What caps it is the ordering handoff: the
vendored AMD/MMD are the `int`-based build, and A's pattern reaches them as `int` arrays, `Ap` (with
`Ap[n] = nnz`) and `Ai`. So `nnz(A)` must fit `int`, `<= 2^31 - 1`, only to be handed to the orderer.
The two ceilings coincide at `INT32_MAX` solely because both are "largest value that fits a signed
32-bit int," one as an id, one as a count in `Ap[n]`. That coincidence is why one constant,
`MAX_IDX = INT32_MAX = 2^31 - 1`, does both jobs. It is a coincidence, not a shared fact.

**Why the dimension cap is `2^31 - 1`, not `2^31`.** Pure index representability would allow
`n = 2^31`: an `n x n` matrix then has largest index `n - 1 = 2^31 - 1`, which fits `int32`. What
removes that last value is the *loop counter*, not the index. Entity loops run
`for (std::int32_t j = 0; j < static_cast<std::int32_t>(size); ++j)`, and at `size = 2^31` the cast
is `INT32_MIN` (negative), so the loop runs zero times, silently; an `int32` counter also cannot reach
`2^31` without overflowing on the final increment. So `n <= 2^31 - 1` is the largest dimension an
`int32` counter can walk, and the dimension cap binds through the counter, not the index. The cost is
one unusable value: the index `INT32_MAX` is representable but never reached (the largest valid index
is `n - 1 = 2^31 - 2`), sacrificed to keep the `int32` entity-loop convention sound. This is why the
guarded constructors reject `size > MAX_IDX` (not `>=`), and it is the same ceiling nnz(A) meets from
the other direction.

**int32 indices and size_t offsets are a self-consistent CSC pairing.** A fully dense matrix at int32
dimension has `nnz = n^2 <= (2^31)^2 = 2^62`, and `size_t` holds `2^64`, so `colPtr[n] = nnz` is
always representable. The condition is `offset_bits >= 2 * index_bits`, and `2 * 32 = 64` fits
exactly (int64 indices would need 128-bit offsets and would *not* fit, which is a further reason int32
is the comfortable resting point). Consequence: no storage or representability guardrail is needed on
A at all; a dense A is representable, and memory is the only real limit. The dense int32 ceiling is
~48 EiB (real), a formality no machine reaches this century. The only guardrail that exists is the
ordering one.

**Where the guard lives, and why only on the matrix.** The `SparseMatrix` constructor throws
`std::length_error` when `mSize > MAX_IDX || nnz > MAX_IDX`. It belongs to the matrix because A is
what gets ordered. It is a limit borrowed from a third-party interface, held at the one seam between
the `int`-based orderer and our `size_t` world, not a limit of our own storage.

**L is not capped, and the factor must not import the cap.** L, and the dynamic factor during
elimination, is never handed to an orderer, so the nnz ceiling does not apply. Its offsets are
`size_t`, its row ids are int32 (same id space, so L inherits the *dimension* bound but not the *nnz*
bound), and it grows to memory. The factor's constructor carries no nnz cap and throws nothing. The
one discipline: nothing on the numeric or solve path may narrow an L count or offset to `int` out of
symmetry with A. The correct design for the factor is the *absence* of the guard. L cannot approach
the dense ceiling regardless: it is triangular (`nnz(L) <= n^2/2`), and a dense block within it is
bounded by `dimension^2` with per-side dimensions bounded by `n <= 2^31 - 1`, exactly what BLAS's
`int` m/n/k accept, so no block dimension can overflow the `int` BLAS wants. That is foreclosed by the
id type, not by any check we write.

**Liftability.** The nnz(A) cap is cleanly removable: SuiteSparse ships a 64-bit
(`SuiteSparse_long`) AMD build, so lifting it is a link-and-widen of the ordering interface, not an
algorithm change, and the constructor's nnz guard then becomes a clean deletion (indices stay int32;
only the orderer's integer width changes). The reason AMD/MMD store nnz in an `int` at all is not
carelessness but a uniform-width choice: the orderer builds and mutates a quotient graph whose ids and
offsets share one integer type, which is cache-friendly for a pointer-chasing kernel, and it ships two
instantiations (`int` and `_long`) rather than mixing widths. So the "better design" is the `_long`
build, a link choice, not a rewrite. Keeping our cap localized and honestly labeled (A-side,
ordering-driven) is what makes its eventual removal a deletion rather than an archaeology problem.

Recorded now because the storage-options experiment's matrix constructors just gained this guard
(mirroring `SparseMatrix`), and it is the model for numfact: the matrix keeps the cap, the factor
analog drops it.

---

## 2026-07-14, Bulk versus direct access: direct wins, and bulk was never worth its one advantage

A consumer that reads a CSC-style object can do it two ways. *Direct*: hold the object and ask it
for one column's (or one supernode's) pointer and length at the moment of use, through the storage's
own lookup (`rowIdx` / `val` / `colSize` on the matrix, `blockPtr` on the factor). *Bulk*:
extract every pointer up front into three plain arrays, then run a kernel that takes nothing but
those arrays, so the object is out of the picture and its dimensions must be extracted alongside.

Two questions decide between them, and they have different answers.

**Is bulk safe? Only for a consumer that does not mutate the object during its own sweep.** A
read-only sweep (matvec, the solve, both read a frozen structure) may be bulk: nothing grows, so no
extracted pointer goes stale. A mutating sweep must be direct: the numeric factorization grows a
front under itself (a delayed pivot reallocates an ancestor's buffer), so a pointer extracted up
front dangles, silently. This is the hard constraint, and it puts `NumFactorEngine` on direct by
necessity rather than preference.

**Is bulk worth it, where it is safe? No.** Bulk's single advantage is one storage-blind compiled
kernel, verifiable in the symbol table, against direct's one instantiation per storage. For us that
trades a real cost for nothing we need:

- The instantiation it saves is cheap. We have a small fixed set of storages, and the arithmetic
  kernels stay shared regardless (they take raw pointers); only the traversal, which is bookkeeping,
  monomorphizes. One extra copy of the bookkeeping per storage is noise.
- It is slower. Bulk streams three extra arrays that direct never touches, plus an `O(n)` extraction
  pass, and on a memory-bound sweep that traffic shows at full price. Measured, the matvec: bulk runs
  1.05x to 1.07x of hand-written flat, while direct matches hand-written flat.
- It is more API and more hazard. The extractor is method surface a reader must understand, and it
  carries an invalidation warning (a stale extracted pointer) that exists only because something was
  extracted. Direct deletes the method, the warning, and the paragraph explaining the warning, all
  at once.

So the rule is direct, and the storage's lookup is the whole interface: a fact about the layout,
answered by the class that owns it, called by a consumer templated on the storage type. No consumer
carries an extractor of its own, which was the repetition the lookup-versus-view corollary set out to
kill and then, by leaving `columnPointers` on the engine, half-kept. That is corrected in the same
entry: `columnPointers` was a bulk copy of the lookups, not a genuine view, and once the lookups
moved onto the matrix it was removed and the matvec reshaped to direct access.

| engine | access | why |
|---|---|---|
| NumFactorEngine | direct | must be: the factor grows under it, so a bulk pointer would dangle |
| SolveEngine | direct | frozen factor, so bulk would be safe, but direct is faster and smaller |
| MultiplyEngine | direct | frozen matrix, same reason |

One thing stays genuinely open, and only one: the deferred multi-RHS solve. There the kernel is
BLAS-3, `O(n^3)` work on `O(n^2)` data, so an extraction's cost amortizes over heavy block work
rather than showing at full price as it does on the memory-bound single-RHS sweeps. If bulk ever
earns its keep it is there, and even then the arrays would be built once from the same storage
lookups, in one shared helper, never a method per engine. The no-per-consumer-extractor rule holds
regardless of how that measures.

This aligns the matrix and the factor on one pattern, which is the point: whether an object is static
or dynamic because the user changes it (the matrix) or because the algorithm changes it (the factor),
the consumer sees one lookup interface and templates over the storage. See "an object offers what its
storage makes cheap" for the lookup-versus-view rule this completes, and the flat-versus-VV entry for
the storage taxonomy.

## 2026-07-14, An object offers what its storage makes cheap, and nothing else

A third rule in the same family as the two below, and it came out of naming the storage-options
experiment's classes for their *purpose* (`SparseMatrixStatic`, `SparseMatrixDynamic`) rather than
their layout (`Csc`, `Vv`).

**The two classes do not have the same API, and that is the design.**

| | static (flat) | dynamic (vector of vectors) |
|---|---|---|
| `setValues`, same structure, new numbers | **yes**, cheap: nothing moves | **yes**, cheap |
| `setColumn`, one column's structure | **absent by design** | **yes**, cheap: the column owns its buffer |
| restructure | build a new one | `setColumn` |

**`setValues` is the mutation a solver actually performs most often**, a Newton iteration, a time
step: same pattern, new numbers, refactorize. The flat layout is perfectly happy with it, so both
classes have it. This is also why "we need a mutable matrix" is a weaker argument for VV than it
first appears: the common case needs no VV at all.

**Amendment (2026-07-14): `setValues` landed at column granularity.** The experiment settled on
`setValues(std::int32_t j, const std::vector<double>& val)`, one column, with an *identical
signature on both classes*: static overwrites the contiguous run `mVal[colPtr[j]..]`, dynamic
overwrites `mVal[j]` in place, and neither invalidates a pointer. Setting every value is a loop over
columns, the write-twin of reading every column through the accessors. This strengthens the point
rather than changing it: the cheap operation is now not merely present on both classes but the *same
call* on both, and the asymmetry stays exactly where it belongs, on `setColumn` (structural), which
remains dynamic-only.

**`setColumn` is absent from the flat one, and its absence is the point.** Changing a column's
*structure* there means shifting every later column: `O(nnz)`, not `O(column)`. An API that looks
cheap and is secretly linear in the whole matrix is a trap, and the caller who writes it in a loop
will not find out until their program crawls. **Refusing to offer it is not a limitation; it is
telling the truth about the storage.**

And it puts the decision where it belongs. The caller knows whether they are changing one column or
rebuilding, and can pick the object that suits: *want a column replaced for you, use the dynamic
one; want to shift data around a flat buffer, do it yourself, this class will not pretend it is
cheap.*

**This is the API-side argument for having no common base class**, and it is stronger than the
performance one already recorded. A shared interface would force one of two lies:

- `setColumn` on the flat matrix, pretending an `O(nnz)` shift is a column operation, or
- `setColumn` on **neither**, crippling the dynamic one to match its sibling's weakness.

So the asymmetry between `NumFactorStatic` and `NumFactorDynamic` is not a wart to be tidied away
when dynamic LDL lands. It is what the two storages **are**, and the interfaces should say so.

**The rule: an object offers what its storage makes cheap, and nothing else.** Which is the same
instinct as the two rules below, seen from the API rather than from the conventions: *do not name a
thing you cannot honour.*

**A second corollary: a lookup belongs to the storage; a view belongs to the consumer.**

This one took a wrong turn first, which is worth recording because the wrong version *looked*
right.

The storage-options experiment puts `columnPointers` on the **engine**, not on the matrix, and that
is correct. The natural generalization is "access functions live in the algorithm class", and on
that reasoning `blockOf` was first written as a private helper on `NumFactorEngine` and again on
`SolveEngine`. Two copies of one expression, which is exactly the duplication the rule below
forbids.

The generalization was too broad. The two things are not the same kind of thing:

| | `columnPointers` | `blockPtr` |
|---|---|---|
| what it does | **materializes** three new arrays | **computes an address** in the existing storage |
| cost | `O(n)`, called once | `O(1)`, called in the hot loop |
| ownership | somebody owns the arrays | nothing is owned |
| shaped by | the traversal that wants it | the storage that holds it |

`columnPointers` is a **view**: a new structure, built in the shape one algorithm wants. A different
algorithm would want a different one, so putting it on the matrix would grow the matrix a method
per consumer, each committing it to a format chosen for somebody else. **A view belongs to the
consumer.**

`blockPtr` is a **lookup**: it answers "where does supernode kk's block live", and the answer is a
fact about the layout, not a shape chosen by anyone. The layout (one flat buffer with offsets, or
one vector per supernode) is the factor's own business, and no consumer should have to restate it.
**A lookup belongs to the storage.**

So `blockPtr` now lives on `NumFactorStatic` and `NumFactorDynamic`, private, with the engines as
friends. One definition per storage rather than one per engine, and the engines cannot tell the two
factors apart. It still inlines to nothing (no `blockPtr` symbol survives in either object file), so
the move is free.

**Amendment (2026-07-14, later the same day): the view was not a view, and the right answer was no
extractor at all.** The rule above is right and `blockPtr` is right, but it mis-cast
`columnPointers`. Calling it a view and leaving it on the engine treated it as a genuine
algorithm-shaped structure, which it is not: the three arrays it builds are a bulk copy of the
per-column lookups the matrix already answers (`rowIdx` / `val` / `colSize`, the matrix-side
twin of `blockPtr`, added to the storage after this entry was written). Once those lookups exist a
consumer reads a column directly, at the point of use, and the extractor is pure redundancy: an
`O(n)` up-front copy of what the storage already holds, carrying an invalidation hazard that exists
only because something was extracted. Its one property, a single storage-blind compiled kernel, has
no value to us: a small fixed set of storages, so direct costs one instantiation per storage, and
bulk measures a few percent slower on the memory-bound matvec besides. So `columnPointers` was
removed and the matvec reshaped to direct access, three lookups per matrix and one templated
`multiply` that calls them, no extractor. What stands: a lookup belongs to the storage; a genuine
view (a transpose, an algorithm-specific reordering) still belongs to the consumer, but
`columnPointers` was never one. See the bulk-versus-direct entry above.

**And a corollary that dynamic LDL will live or die by: structural mutation invalidates every
pointer previously extracted; value mutation does not.**

```
setValues   does NOT invalidate.  The buffer stays put; only its contents change.
setColumn   DOES invalidate.      The column's buffer is replaced; anything into it dangles.
```

The experiment demonstrates it rather than asserting it (`testInvalidation` extracts the pointers,
mutates, and observes which still point where they did). The rule holds in both storages, and it is
exactly the rule the dynamic factor needs: a delayed pivot grows an ancestor's front, which
reallocates its buffer, which dangles every pointer previously taken into it. **Silently.**

```cpp
eng.blockPointers(f, block);        // extracted once
for (kk) {
    ... factor kk ...
    f.mVal[pp].resize(bigger);      // a delayed column grows ancestor pp
    ...                             // block[pp] is now DANGLING
}
```

Nothing in C++ enforces this. The remedy is to **fetch a supernode's block pointer at the moment of
use rather than up front**, which is one indirection and which `storage-options` measured at
essentially nothing. The alternative, re-extracting after every growth, is more code and more to
forget.

This is the one thing the experiment does *not* rehearse, since its structures change only between
runs of the algorithm, never during one. Worth knowing before writing dynamic LDL rather than
after.

## 2026-07-14, Two rules for conventions: one predicate, and a sparse choice matrix

Two ideas that keep recurring, worth stating once rather than rediscovering.

### One predicate, in one place, applied everywhere

**When several places must agree about a convention, they must not each decide it.** Give them one
named thing to ask, and let it be the only place the decision exists.

**Example 1, and it is the reason for the rule: `hermitian()`.** The factorization and the solve
must both know whether the factor conjugates. 10.12 lets each site decide for itself:

```cpp
// in the factorization
SYRK('L','N', ...)              // decides: no conjugate
GEMM('N','T', ...)              // decides: no conjugate, again

// in the solve, a different file
y[lc] -= y[lr] * val[lij];      // decides: no conjugate, a third time
```

Nowhere is it written down that a Cholesky factor is Hermitian. Each site *inferred* the
convention from context, and each inherited LDL's (symmetric) habit. So Cholesky is wrong in
**all three**, independently, and fixing one would not fix the others. Nothing reveals the
disagreement, because there is nothing to disagree *with*.

Ours has one predicate:

```cpp
static bool hermitian(Factorization f) {
    return f == Factorization::Cholesky
        || f == Factorization::StaticLDLH
        || f == Factorization::DynamicLDLH;
}
```

The factorization asks it; the solve asks it. Add a factorization type and both follow. Get the
predicate wrong and *everything* breaks loudly, which is far better than one site being quietly
wrong.

**Example 2, and it shows the other half of the rule: `exactPatterns()`.**

Symbolic factorization may read one front column per supernode when the supernodes' columns share
a pattern exactly, and must read them all when they do not. Who decides?

Not the caller. A `bool useExactPatterns` parameter on `SymFactorEngine` would put the decision in
the hands of whoever wires the phases together, and that person **cannot know the answer**:
amalgamation at threshold zero merges only free merges, so it *may or may not* store a zero
depending on the matrix. Threshold alone does not determine it. Only the forest engine, which did
the merging, knows what it actually did.

So the forest **records** it (`ElmForest::exactPatterns()`, set false by `compressThreshold` only
when a merge actually pays fill), and symbolic factorization **reads** it. One place decides, at
the moment it can, and everyone downstream asks rather than infers.

**The two examples differ in an instructive way.**

| | `hermitian()` | `exactPatterns()` |
|---|---|---|
| what it is | a pure function of a setting | a fact about what happened at runtime |
| where it lives | computed on demand, anywhere | **recorded** by the code that knows |
| why | the answer is determined by the input | the answer is not determined by any input |

So the discipline is one rule with two shapes: **if the convention is a function of the inputs,
make it a predicate; if it depends on what happened, record it at the moment it happens. Either
way, decide once and let everyone else ask.** The failure mode both avoid is the same: a
convention re-derived at each use, drifting apart, with nothing to compare against.

`Blas<Val>::conjTrans` is the same rule at the smallest scale: one trait, not `'T'` and `'C'` typed
out at each call site.

### The choice matrix is sparse. Enumerate its nonzeros, not the product.

The obvious modelling of "what do we factor" is a product:

```
(Cholesky, LDL)  x  (transpose, conjugate-transpose)   ->  four combinations
```

and then one of the four (**complex symmetric Cholesky**) has to be *forbidden*, because positive
definiteness is meaningless for it. A runtime rejection, guarding a combination the API invited
the caller to ask for.

**We did not do that.** The enum lists only the combinations that exist:

```cpp
enum class Factorization {
    Cholesky,      // always A = LL^H. There is no LL^T variant to name.
    StaticLDLT, StaticLDLH,
    DynamicLDLT, DynamicLDLH
};
```

Cholesky carries no transpose suffix because **it has no choice to offer**: it is always the
conjugate transpose, and over the reals that *is* the plain transpose. LDL carries one because it
genuinely has two forms. The forbidden cell is not rejected; **it does not exist**.

The cost is one letter of naming convention. The gain is that a whole class of runtime error is
unrepresentable, which is the same principle as `BlasLapack`'s operation-named wrappers: **make the
wrong thing unwriteable, not merely refused.**

Worth naming the general shape, since it will recur: a design space is rarely a full product. Model
the combinations that mean something and leave the rest unnameable. In a sparse-matrix library the
metaphor is right there: **the choice matrix is sparse, so store its nonzeros.**

## 2026-07-13, The solve, and the first test that checks the pipeline rather than a phase

**Every test before this one checks a phase against an oracle.** The forest against a
recomputation, the symbolic factor against a dense pattern, the numeric factor against a dense
Cholesky or by reconstruction. Each says *this phase computed what it should*. **None says the
phases compose.**

```
|| A x - b ||  /  || b ||
```

says that, in one number, through ordering, elimination forest, symbolic factorization, numeric
factorization, triangular solve and sparse matvec. Six phases. It is where two phases disagreeing
about a convention would show, an ordering, a conjugate, an index base, and nothing else we have
would catch that. It comes out at 3e-16 for every factorization, both traversals, both scalar
types.

The right-hand side is manufactured from a known solution (`b := A x`), so the test needs no
reference solver. And it checks the **residual**, not the distance to the manufactured `x`: those
differ by the conditioning of `A`, and the residual is the honest thing to require of a direct
solver.

**`MultiplyEngine` exists for this**, and only for this. It is fifteen lines and it is what makes
the check possible.

**One right-hand side, and the scalar solve that follows from it.** 0.9 has two vector classes and
the split is principled: with a *single* right-hand side there is no level-3 BLAS to be had, so its
`SingleVector` solve is scalar and works directly on the vector through indirect indexing. With
*many*, a supernode's rows become a dense block and the solve becomes TRSM and GEMM, which is worth
the gather and scatter that packing demands. We take the first. The multi-column path is real and
worth adding; it is a performance path, not a correctness one.

**And the conjugate in the backward pass is where 10.12 repeats its bug.** Its backward solve is

```cpp
y[col] -= y[row] * val[...];        // no conjugate: this applies L^-T, not L^-H
```

which is right for its complex-symmetric LDL and **wrong for its Cholesky**, exactly as its `SYRK`
is. Ours conjugates when the factorization does, using the *same* `hermitian()` predicate the
factorization uses, so it is one rule stated once rather than two rules that must be kept in step.
That is the pattern to hold on to: **when two places must agree about a convention, give them one
predicate, not two copies of a decision.**

## 2026-07-13, Static LDL: three kernels BLAS does not have, and why the traversals did not change

**LAPACK has no unpivoted LDL.** `?sytrf` is Bunch-Kaufman, which pivots, and pivoting is exactly
what a *static* factorization refuses to do. So the kernel is ours. 0.9 wrote it; we port it.

**Three kernels, all from 0.9, all absent from BLAS:**

- **`ldl`** (0.9's `OBLIO_POTRF2`). Recursive: split the block in half, factor the leading part,
  solve for the off-diagonal, form its upper counterpart, take the Schur complement, recurse.
- **`formUpper`** (`OBLIO_COMPUTE_U`). `U := D L^T` (or `D L^H`), into the upper triangle.
- **`gemmLower`** (`OBLIO_GEMM`). `A -= L U`, filling only the lower triangle, because the product
  is symmetric. BLAS has nothing for this: `syrk` computes `A A^T`, and there is no "`A B` where
  the product is known symmetric".

**A nice factoring in 0.9 worth naming.** `OBLIO_POTRF1` and `OBLIO_POTRF2` are the *same
algorithm*, differing only at `n == 1`: the first fails on a non-positive pivot, the second
replaces a tiny one and counts it. So 0.9's `NO_LAPACK` "Cholesky" is really an LDL that refuses
to proceed unless `D > 0`, which is legitimate (`LDL^T` with positive `D` *is* Cholesky, with
`L_chol = L_ldl sqrt(D)`; same factorization, different storage). We need only the second, since
LAPACK gives us a real Cholesky.

**The storage, which the statistics predicted.** In an LDL block:

```
the diagonal        holds D          where L's implicit 1s would be
the lower triangle  holds L          unit lower triangular
the upper triangle  holds U = D L^T  which Cholesky leaves as explicit zeros
```

So an LDL block uses the whole rectangle. `U` is not redundant: the recursion needs `D L^T` in two
places (to solve for the next `L`, and to form the Schur complement), so it is computed once and
kept. We reasoned to exactly this layout from 0.9's `numberOfEntries` versus
`numberOfAllocatedEntries` counts, before seeing the code.

**Perturbation is not a refinement; it is the only recourse.** A static factorization cannot pivot,
so a pivot too small to divide by has no remedy but replacement. We then factor a matrix slightly
different from the one we were given, which is a real thing to have done, so the count is reported
on the factor (`NumFactorStatic::numPerturbations`) rather than hidden. It belongs on the factor,
not the engine: it is a property of what was computed, not of the thing that computed it.

**And the update is not a rank-k operation, which is the deep difference from Cholesky.**

```
Cholesky:  T -= L21 L21^H          one HERK. Exactly what a rank-k routine computes.
LDL:       T -= L21 D L21^H        no BLAS routine at all: the D in the middle rules it out.
```

So LDL forms `U := D L21^H` into a **scratch block**, then multiplies. That scratch is the whole
price of the `D`, and it is why Cholesky's update is one call and LDL's is three.

**The traversals did not change by a line, and that was the design working.** `factorSupernode` and
`updateSupernode` dispatch on the factorization type internally, so left-looking and right-looking
decide *when* to factor and *when* to update, never *how*. Adding LDL touched those two functions
and one `switch`. That separation is worth protecting when multifrontal arrives: a traversal is a
schedule, not an arithmetic.

**Complex Hermitian LDL (`StaticLDLH`) is an extension, not a port.** 0.9's complex LDL is
symmetric only, which is why its complex LDL correctly uses `SYRK` and `'T'`. We support both, and
the `T`/`H` distinction runs through all three kernels as a single `bool hermitian`, which for
`double` is a no-op in every branch. That is the honest expression of "real symmetric and real
Hermitian are the same case".

**Checked by reconstruction, which is a better oracle than a second implementation.** Cholesky is
compared against an independently written dense Cholesky. LDL is checked by multiplying the factor
back out: `L D L^H == P A P^T`. That needs no second implementation, validates `D`, `L`, the
storage layout and the supernodal assembly in one statement, and works unchanged across all three
symmetries, which a dense oracle would not.

## 2026-07-13, Symmetry is part of the factorization, not a setting beside it

*(Supersedes the "symmetry is determined, not chosen" half of the entry below, written earlier
today. That claim was true for Cholesky and false the moment LDL arrived. Recording the
correction rather than editing the original, because the reason it was wrong is worth keeping.)*

**What changed.** The earlier entry argued that `Val` and the factorization type fix the symmetry
between them, so symmetry never needs to be named:

```
real    + anything   ->  symmetric
complex + Cholesky   ->  Hermitian
complex + LDL        ->  symmetric
```

The third line is where it breaks. It is what **0.9** does, not what LDL *is*. A complex matrix
may be symmetric (`A = A^T`, so `D` comes out complex) or Hermitian (`A = A^H`, so `D` comes out
real), and both are legitimate, useful factorizations. `LDL^T` and `LDL^H` are different
computations, and nothing in `Val` or in "static LDL" says which one is wanted. **Symmetry becomes
a genuine choice.**

**The choice goes into the factorization type, not beside it.**

```cpp
enum class Factorization {
    Cholesky,      // A = LL^H  (LL^T for real). Positive definite, hence Hermitian.
    StaticLDLT,    // A = LDL^T          complex: D complex
    StaticLDLH,    // A = LDL^H          complex: D real
    DynamicLDLT,
    DynamicLDLH
};
```

**Why here and not as a `Symmetry` flag beside it**, which was the obvious alternative: the same
principle that shaped `BlasLapack`, **make the wrong thing unwriteable.** A separate flag would let
a caller ask for `Cholesky` + `Symmetric` + complex, which we would then reject at runtime. In the
enum above there is no such value. The combination is not *forbidden*; it does not *exist*. That is
strictly stronger, and it costs one letter in a name.

A flag would also be **inert for real**, where symmetric and Hermitian coincide. An API parameter
that is sometimes meaningless is a small lie, and worth avoiding for the price of an enumerator.

The cost, stated honestly: the enum grows, and `T`/`H` is a suffix convention a reader must learn.
Both seem cheap against a runtime rejection that exists only to guard a combination we could
simply not have named.

**And complex Hermitian LDL is an extension, not a port.** 0.9 does only the symmetric one, which
is why its complex LDL correctly uses `SYRK` and `'T'`. That its *Cholesky* also uses them is the
bug we have already recorded; the two are related, since Cholesky is the one place 0.9 needed the
Hermitian convention and did not reach for it.

## 2026-07-13, The factorization space, and a BLAS layer that names operations rather than routines

Entering the numeric phase, three questions had to be settled together: which combinations of
scalar type, symmetry and factorization we support, how the code selects among them, and what
shape the BLAS wrapper takes. They turn out to be one question.

**The space, and it collapses.** Three axes, `Val` in {real, complex}, symmetry in
{symmetric, Hermitian}, factorization in {`CC^T`, `LDL^T`}:

| `Val` | symmetry | `CC^T` | `LDL^T` |
|---|---|---|---|
| real | symmetric **is** Hermitian | yes, SPD | yes, indefinite |
| complex | Hermitian | yes, HPD (`zpotrf`) | yes; 0.9 does not do it, we will |
| complex | symmetric | **forbidden** | yes, the standard case |

**For real the symmetry axis does not exist.** Conjugation is the identity, so `A^H = A^T` and
the two conditions are the same condition. That row is one row wearing two hats.

**Complex symmetric Cholesky is forbidden, and not merely unimplemented.** Positive definiteness
means `x* A x > 0`, which requires that quantity to be *real*, which happens for all `x` exactly
when `A` is Hermitian. For a complex *symmetric* `A`, `x^T A x` is a complex number and the
inequality does not even typecheck. Concretely: `A = [[0, 1], [1, 0]]` is symmetric and
nonsingular, and `a11 = 0` kills the very first square root. Cholesky has no pivoting to recover
with, because not needing pivoting is the entire point of Cholesky, and that guarantee comes from
positive definiteness. LAPACK has no complex-symmetric Cholesky for precisely this reason. So the
API rejects it, as a hard error rather than a to-do.

The stable factorization for complex symmetric matrices is `LDL^T` with 2x2 pivots, which is
already on the plan. Complex symmetric `CC^T` would be an unstable duplicate of a thing we are
building anyway.

**Complex Hermitian `LDL^T` (`A = LDL^H`, `D` real) is a gap in 0.9, not a gap in the
mathematics.** It is a perfectly good factorization, indefinite Hermitian, and we intend to
support it. 0.9 simply never wrote it. Worth distinguishing sharply from the cell above: one is a
mathematical impossibility and must be rejected forever, the other is work not yet done. They
must not report the same error.

**So symmetry is determined, not chosen, and there is no third parameter.** Given `Val` and the
factorization type, everything else follows:

```
real    + anything   ->  symmetric  ->  'T',  syrk
complex + CC^T       ->  Hermitian  ->  'C',  herk
complex + LDL^T      ->  symmetric  ->  'T',  syrk
```

`Val` is the template parameter and `FactorType` is an engine setting (10.12 already has the
enum: `eCC`, `eStaticLDL`, `eDynamicLDL`). Between them the transpose character and the rank-k
routine are fixed. No symmetry flag, no extra argument.

**The BLAS layer names operations, not routines. This is the part 0.9 gets wrong.** A wrapper is
needed regardless: BLAS is a Fortran interface, everything by pointer, and without a layer every
call site carries `&uplo, &trans, &n, &k, &alpha, ...` plus a branch to choose `d` or `z`. 0.9
has such a layer, overloaded inline wrappers on `Real*` and `Complex*`, and that much is right.

What is wrong is *what it wraps*. 0.9 wraps BLAS **routine by routine**: `SYRK`, `GEMM`, `TRSM`,
`POTRF`, leaving the caller to choose `'T'` versus `'C'` and `SYRK` versus `HERK`. The convention
therefore leaks into the engine, and the engine gets it wrong: **0.9's complex Cholesky calls
`SYRK` and `TRSM('T')` and `GEMM('N','T')`, all of which are the complex *symmetric* pattern,
while `POTRF` maps to `zpotrf_`, which is Hermitian.** There is no `HERK` anywhere in 0.9's
`BlasLapack`. For Hermitian `A = LL^H` the update must be `L21 L21^H` and the off-diagonal solve
must be against `L11^H`; using `'T'` is correct only when `L11` is real. Almost certainly never
exercised on a genuinely complex Hermitian matrix.

Our layer therefore exposes an operation whose *meaning* is fixed and lets the type pick the
routine:

```cpp
// "A times A-conjugate-transpose", whatever that means for this Val.
void herk(char uplo, char trans, ...);   // double -> dsyrk_ ;  complex -> zherk_

template<class Val> struct Blas;
template<> struct Blas<double>               { static constexpr char conjTrans = 'T'; };
template<> struct Blas<std::complex<double>> { static constexpr char conjTrans = 'C'; };
```

so the Cholesky kernel is one piece of code, correct for both:

```cpp
potrf('L', f, block, ld);
trsm ('R', 'L', Blas<Val>::conjTrans, 'N', u, f, ...);
herk ('L', 'N', ...);
gemm ('N', Blas<Val>::conjTrans, ...);
```

No branch, no `if constexpr`, and **0.9's bug becomes unwriteable**: the engine never names
`SYRK` or `HERK`, so it cannot pick the wrong one. `syrk` and a literal `'T'` remain available in
the header for `LDL^T`, where plain transpose is what is wanted, and there the *algorithm* asks
for them explicitly rather than inheriting them by accident.

**Storage: this is what the storage-options experiment was for.** Yesterday's study
(`experiments/storage-options/`) established that one compiled algorithm serves both a flat CSC
storage and a vector-of-vectors, through nothing but a pointer array, and that the abstraction
costs about nothing (1.07x flat, 1.10x packed VV, one `multiply` symbol in the binary). It was
run against exactly this moment. `NumFactor` uses both:

- **static** (`CC^T`, static `LDL^T`): **flat**. Symbolic has already sized every block, nothing
  grows, so one buffer with per-supernode offsets. 0.9 does the same (`FactorsStatic` allocates
  one array and points into it).
- **dynamic** (dynamic `LDL^T`): **VV**. Delayed pivoting grows a front at runtime by an amount
  symbolic never predicted, and the growth is local. 0.9 does the same (`FactorsDynamic`
  allocates one array per supernode).

**Static and dynamic, not flat and VV**, in the naming. The layout is a *consequence* of
mutability, not the thing being chosen, and 0.9 names its classes the same way. Flat-versus-VV
describes the bytes; static-versus-dynamic describes why.

**Status.** Settled: the table above, symmetry determined rather than chosen, operation-named
BLAS, static-flat and dynamic-VV. The objects are `NumFactor` and `NumFactorEngine`, the engine
taking a factorization type (`CC^T` first, then `LDL^T` static and dynamic) and a traversal
(left-looking, right-looking, multifrontal). Cholesky first, left- and right-looking first.

## 2026-07-12, Amalgamation: a second compression, and the third bug in 10.12

**Two compressions, and they are not the same algorithm with a knob.** Fundamental compression
contracts *paths*: it merges a supernode with its child when that child is an **only** child
and shares its pattern. Amalgamation contracts *stars*: it merges a supernode with any number
of its children, and will pay for the privilege in explicitly stored zeros, up to a budget.
They are orthogonal settings in the engine (`Supernodes::Fundamental` and an optional
threshold), as they are in both references, and both may run, fundamental first.

**Absent is not zero.** The threshold is a `std::optional<std::size_t>`, because "do not
amalgamate" and "amalgamate but pay nothing" are different instructions. At threshold zero it
still merges, and it merges strictly more than fundamental compression does.

**Which is the interesting fact: fundamental supernodes are not maximal.** The zero-fill
condition and the fundamental pattern condition are the *same* test. The only difference is
the only-child requirement, which exists to make supernodes paths, and hence unique, not
because merging would cost anything. Drop it and free merges appear. The smallest case is a
three-column star: fundamental gives 3 supernodes, amalgamation at threshold zero gives 2, at
no cost. On the grid of the notes it takes 7 supernodes to 6, column 3 joining the separator
for nothing.

**Uniqueness is what is being traded away, and it is worth stating plainly.** Where two
children could each merge for free, only one can (absorbing the first widens the front, which
prices out the second), so the algorithm must **break ties**, and a tie-break is a convention
rather than a theorem. 0.9's rule is: least fill, then largest front, then first in the child
list. That last clause is arbitrariness made deterministic, which is exactly what a canonical
algorithm never needs. Fundamental supernodes are unique because they refuse to make the
choice; amalgamation is not, because it makes it.

**The third 10.12 bug, and the most instructive.** 0.9 updates a parent's front size after it
absorbs children:

```
frontSizeArray[kk] += frontIncrement;
```

10.12 has that line, **commented out**, with the explanatory comment above it left intact:

```
// update the front size of supernode s2.
//numFrntIdxsArr[s2] += frntInc;
```

It matters. Parents are processed in increasing label order, and a child's label is below its
parent's, so by the time a parent is reached its children **have already been parents
themselves** and may have absorbed children of their own. Pricing a merge needs the child's
*current* front size. 10.12 reads the stale one, understating both the fill and the resulting
block, silently. And the cause is visible in its own source: the array is declared `const`, so
the update would not compile, and someone commented out the statement rather than fix the
constness. We restore it.

That is now three bugs found in 10.12 (the sibling-link copy in `SymbolicEngine`, the wrong
front-size operand in the fundamental merge test, and this), and two of them left a corpse in
place: a commented-out line whose comment still promises the behavior. The pattern is worth
naming, since it will recur: **10.12 transcribes 0.9's prose faithfully and its code
approximately.** Read the comments for intent; verify the code against 0.9.

**The unsigned subtraction is safe by theorem, and says so.** The merge cost contains
`|indexSet(K)| - |update(J)|`, which on unsigned types would wrap if it went negative. It
cannot: by the containment theorem `update(J)` is a **subset** of `K`'s index set, so this is
a set-difference size. We name it (`zerosPerCol`) and cite the theorem, rather than leaving a
subtraction that happens not to wrap. Unlike the fundamental merge test, it cannot be
rearranged into an addition, because the running fill total genuinely needs the difference.

**Testing what survives tie-breaking.** Since the partition is not canonical, the tests assert
only what is invariant under the tie-break: at threshold zero the stored-zero count is exactly
zero; the factor's true nonzeros are all still present; the supernode count never rises, and a
larger budget never raises it; the links, height and topological labelling stay valid. The
specific partition is deliberately not asserted, except on the star and the grid, where it is
forced.

## 2026-07-12, Friendship is a write grant, and reading needs no friend

**Every argument an engine takes falls into one of three cases, and the third is the only one
that is subtle.**

**Written: friend, and pass the object.** An engine is declared `friend` by exactly the object
it fills, and by no other. `ElmForestEngine` writes `ElmForest`, `SymFactorEngine` writes
`SymFactor`, `OrderEngine` writes `Permutation`. That is the complete list, and friendship exists
for no other purpose. Having granted it, pass the object rather than its fields: the engine can
reach them regardless, so enumerating them in the signature restricts nothing and only makes it
long. This is the subject of the entry above.

**Read: pass the object, through its public API, with no friendship.** We considered granting
the engines read friendship on `SparseMatrix` and `Permutation`, on the theory that it would
compact the signatures further. It compacts nothing. Every read already goes through an
accessor returning a `const&`, so there is no copy to avoid and no access to gain: friendship
would only let us write `A.mColPtr` where we write `A.colPtr()`, at the price of new friends on
objects whose headers say plainly that they have no writer. `SymFactorEngine` reads a dozen
fields of `ElmForest` and is not its friend, which is exactly right.

**Read, and only part of the object is needed: take that part.** `SparseMatrix` is the one
object an engine does not simply take whole, and the reason is not access but need. Ordering,
the elimination forest and the symbolic factorization are graph algorithms: they read a
sparsity pattern and never touch a value. So they take one. `SparseMatrix` offers two overloads,
one taking the matrix and one taking `colPtr` and `rowIdx`; the second is the implementation,
the first a one-line adapter over it. See the entry below.

**Friendship was never the constraint, and it is worth being precise about that.** The
field-taking signatures are not there because access is blocked, it is not. They are there
because a structural algorithm has no business asking for values. No amount of friendship would
change that, and granting it would not have shortened a single signature.

## 2026-07-12, Engine helpers take whole objects, or exactly the part they need

**Passing an object's pieces to a function that is already its friend restricts nothing.**
The engine can reach every field regardless; unpacking the arrays into the signature does not
narrow that access, it only spreads it across nine parameters. We had been writing helpers
that took the forest apart, on the theory that a signature listing `const&` inputs and `&`
outputs documents the data flow. It does not, and the cost of pretending otherwise was real.

**The signature could not carry the contract that actually matters.** `compressFundamental`
leaves `mFirstChild`, `mLastChild`, the sibling links and the roots **stale**: they still
describe the nodal forest and are wrong for the compressed one until rebuilt. A `const&`
parameter says "not written", which a reader takes as "still valid", which is exactly
backwards. The read set, the write set and the stale set always lived in a comment. We were
paying nine parameters for documentation we never received.

**Meanwhile the parameters were manufacturing a bug class.** `finalizeLinks` took five
`std::vector<std::int32_t>&` arguments in a row. Transposing `firstChild` and `lastChild`
compiles silently and builds a subtly wrong forest. Same type, same reference-ness, no help
from the compiler. `finalizeLinks(f)` cannot express that error at all.

**So: helpers take the object.** Across `ElmForestEngine` this took the private helpers from
38 parameters to 11 (`finalizeLinks` 9 to 1, `compressFundamental` 8 to 1, `computeHeight` 5
to 1). The gain is not brevity. It is that `compute()` now shows the shape of the algorithm,
parent links, child links, sizes, compress, child links again, height, instead of the
plumbing, and a structurally interesting fact like `finalizeLinks` running twice is legible
at a glance rather than buried in argument lists.

**The read/write/stale sets move to the comment, where they can be said properly.** That is
not a loss. A signature can express "I do not write this"; only prose can express "I leave
this stale, and you must rebuild it before anything else reads it".

**Three cases, and the third dissolves rather than trading off.**

- **Objects we write through friendship** (`ElmForest`, `SymFactor`): pass whole. This is what
  the entry is about, and it is also 0.9's shape, whose `compress_()` is a member taking no
  arguments at all. Our array-passing came from following 10.12.
- **Objects we only read, and take whole** (`Permutation`): pass whole. No friendship needed,
  the public accessors suffice, and it is free. This collapses `oldToNew, newToOld` into one
  parameter with no downside whatever.
- **Objects we only read, and of which we need only a part** (`SparseMatrix`): take the part,
  and adapt at the boundary. See below.

**Pass only the structural part of a matrix when only structural work is done.** Ordering, the
elimination forest and the symbolic factorization are graph algorithms. They read a sparsity
pattern; they never touch a value. So the implementation takes a pattern, and the overload
taking a matrix is a one-line adapter over it:

```
template<class Val>
bool compute(const SparseMatrix<Val>& A, const Permutation& p, ElmForest& f) const
{ return compute(A.colPtr(), A.rowIdx(), p, f); }         // adapter, one line

bool compute(const std::vector<std::size_t>&  colPtr,
             const std::vector<std::int32_t>& rowIdx,
             const Permutation& p, ElmForest& f) const;   // the engine
```

**This is about honest dependencies, not about C++.** A function's parameters should say what
it actually consumes. A structural algorithm that demands a matrix is lying about what it needs,
and it forces a caller holding only a graph to fabricate numbers to satisfy a signature that
will ignore them. Both overloads are public for that reason: the lower one is not an internal
shortcut, it is the **graph interface**. The whole structural pipeline,
`OrderEngine -> ElmForestEngine -> SymFactorEngine`, runs on a bare graph today: no `SparseMatrix`,
no scalar type, no numbers. The rule would be right in any language, and we should not restate
it as a fact about templates.

**The template mechanics are a consequence, not the motive.** In C++ the structural overload
happens to be non-templated, so it is compiled once rather than once per scalar type, and the
adapter is a forwarding line that inlines away. Pleasant, and worth nothing on its own. We first
wrote this entry as though `Val` were the reason, which inverted cause and effect and made an
honest interface look like a workaround for a language wart. It is not: the pattern-taking
overload would earn its place even if every type in the codebase were concrete.

**Adapt once, at the public boundary.** The adapter belongs on the entry point and nowhere
else. We first tried it on the private helpers too, and those adapters immediately became dead
code: once `compute` has unpacked the matrix, every helper below it already holds the pattern,
and a second layer has no callers. All three engines now have exactly one adapter each, on
`compute`.

This is not a new idea in the code, only a newly named one: `compute` already pulled `colPtr`
and `rowIdx` out of `A` and handed them to the helpers. We had simply not seen that the same
move, applied at the entry point, makes the pattern a public capability rather than an internal
convenience.

**`SparsityPattern` is a packaging question, not a design one.** `SparseMatrix` is two things
under one name: a pattern, and values indexed by it. Since the interface already passes the
pattern, introducing a type to name it would only replace two array parameters with one. The
interface is already right; whether the pattern travels as a named type or as two arrays is a
matter of convenience. Recorded as an improvement available, not a debt owed.

**Entry points are named `compute`.** One verb across every engine, since an engine's job is
to derive a fact from its inputs. `OrderEngine::order` was renamed to match. 0.9 calls them
all `run`, which is uniform but says nothing; 10.12 uses `ComputeElmForest` and
`ComputeSymbolic`, right about the verb but appending a noun the class already carries
(`ElmForestEngine::ComputeElmForest` stutters). `ElmForestEngine::compute(A, p, f)` says what
is computed three times over: in the class, in the arguments, and in the output.

## 2026-07-11, Choice objects are constructible; derived-fact objects are engine-filled

**The pipeline holds two kinds of object, and they deserve opposite rules.** A permutation
is a choice. Given a matrix there are n! valid permutations, and which one is wanted is
policy: AMD, MMD, nested dissection, the problem's own numbering, an ordering read back from
a file, an ordering composed from two others. An elimination forest is not a choice. Given a
matrix and a permutation there is exactly one correct forest, one correct symbolic
factorization, one correct numeric factor. They are derived facts, and any value other than
the derived one is simply wrong.

**So: choice objects are freely constructible, derived-fact objects are engine-filled only.**
`SparseMatrix` and `Permutation` are inputs in the sense that matters, they encode decisions
the caller is entitled to make, so they take a full public interface for building and setting
them. `ElmForest`, `SymFactor` and the numeric factor to come are written only by their engine,
through friend access, and expose read-only accessors. A caller-supplied elimination forest
is meaningless at best and silently wrong at worst, so nothing is gained by permitting one,
while a whole class of bugs is prevented by forbidding it.

**Input versus output is the wrong test, and would have misled us here.** A permutation is
`OrderEngine`'s output and `ElmForestEngine`'s input, so that framing gives no answer. The
question that does give an answer is whether the object is determined by what precedes it.
Determinacy, not position in the pipeline.

**A consequence worth stating: an ordering engine is optional.** `OrderEngine` is a
convenience for computing a good permutation; it is not the sole authority on what a
permutation may be, which is why `setOldToNew`/`setNewToOld` are public API rather than a
testing backdoor. This is not an innovation, 0.9 had `set`, `get`, `read`, `write`, `compose`
and the grid orderings, all public. We had merely under-ported the class, and the missing
setter first showed up as an inability to test the permutation maps against a known answer.
The friend grant to `OrderEngine` survives, but only as an optimization (it can skip
re-checking a bijection it just built), not as the mechanism by which a permutation gets
filled.

## 2026-07-11, Symbolic factorization: 10.12's design, 0.9's behavior

**The source-of-truth rule got its first real test here, and both halves of it fired.** The
union recurrence at the heart of symbolic factorization is identical in 0.9 and 10.12: the
same two-part union (the sparsity patterns of a supernode's own front columns, then the
update indices of its children), the same two skip tests, the same marker array. There was
no behavioral choice to make. Every difference between the two references was a matter of
shape, and in one case a matter of one of them being wrong.

**We took 10.12's naming, and it is a real improvement.** `s1`/`s2` for a child and parent
supernode, `lc`/`lr` for a local (factor-order) column and row, `ac`/`ar` for their
counterparts in the original matrix. 0.9 uses `jj`/`kk`/`i`/`di_`, which is hard to read and
harder to review. The 10.12 vocabulary also lines up with the `lc1`/`lc2` names already in
our `computeParent`, so the port reads consistently across engines.

**Front size is computed by counting the map, not by filling ones.** 10.12's
`rComputeNumIdxs` fills the front sizes with 1 unconditionally. We count how many columns
map to each supernode. Both run only on the nodal forest, where the map is the identity, so
both give all ones and the choice is stylistic: counting derives the value from the map
rather than asserting it. We originally justified this as generality, on the grounds that it
would stay correct after compression. That was wrong, and worth recording as such:
compression derives the merged front sizes itself, and the rest of the function is
nodal-only anyway (its update-size walk indexes by column, which coincides with supernode
only while supernodes are trivial), so it is never re-run afterwards. The map-count buys no
capability. This is also not a departure from 0.9, whose symbolic factorization recomputes
front sizes from the map with exactly this loop; we merely compute them one stage earlier,
on the forest, where they belong as an attribute.

**SymFactor stores its index sets flat, as 0.9 does, not as a vector of vectors.** This
follows directly from the flat-vs-VV decision below: the symbolic factor is written once
into a structure whose size the forest already knows (`frontSize + updateSize` per
supernode), so a flat buffer with per-supernode offsets is the right shape and stays
comparable to 0.9 buffer-for-buffer. 10.12 uses a vector of vectors here and we do not
follow it. One modernization within the flat layout: 0.9's `pointerToIndex` holds absolute
pointers into the index array, while we hold `std::size_t` offsets. Same layout, but the
offsets bracket each supernode's block as `[s]`/`[s + 1]` and match the convention already
used by `SparseMatrix`.

**SymFactor copies three links, not five.** The forest is doubly linked (parent, first and
last child, next and previous sibling), because `computeHeight` walks it backward. Symbolic
factorization only ever walks a child list forward, so it needs parent, first child and next
sibling, and that is exactly what both 0.9 and 10.12 store in their symbolic object. We
follow. The other two are cheap to add later if numeric factorization wants a backward walk.

**10.12's symbolic engine has a bug, and finding it is the point.** It builds its sibling
vector with `std::copy(elmForest.mLstChldVec.begin(), ..., nxtSblgVec.begin())`, copying the
last-child links into the next-sibling links. The factorization itself still comes out right,
because the recurrence reads the sibling links straight off the forest rather than off the
copy, but the sibling links stored in the resulting symbolic object are corrupt, and those
are what numeric factorization would later traverse. This is the concrete reason the rule is
"favor 10.12's design, verify against 0.9's behavior" and not "follow 10.12". We took 0.9's
correct three-link copy.

**The generality for compression is written but currently dead.** The recurrence unions the
patterns of every front column of a supernode, not just the first, because threshold-based
compression groups columns whose patterns are merely similar rather than identical. While
supernodes are trivial every supernode has exactly one front column, so that loop runs once
and the generality is untested. It is deliberately kept (both references keep it) and will
get its first real exercise when fundamental-supernode compression lands.

**Verification is against an independent oracle, not against ourselves.** The tests compare
every supernode's index set to a dense simulation of Cholesky fill (eliminate a column, make
its subdiagonal rows pairwise adjacent), which shares no code path with the forest or the
symbolic factorization. Agreement is therefore evidence rather than tautology. It runs under
both the natural ordering and AMD, so the permutation maps are exercised too, and it
incidentally reconfirms the forest's update sizes on the same matrices.

## 2026-07-11, Flat vs vector-of-vectors storage, and the descriptor view that spans both

*Partly measured (2026-07-12, `experiments/storage-options/`): the abstraction cost and the
layout cost are now numbers, not guesses, and one of this entry's original claims is superseded
by them. What remains provisional is the dynamic-numfactor storage, still recorded ahead of having
0.9's numfactor code, and the flatten-or-not question the measurement reopened. See Status.*

**Not an API question, a data-structure question.** Flat storage (one contiguous buffer
plus a per-column/per-supernode offset array, CSC style) versus vector-of-vectors (VV, one
inner `std::vector` per column/supernode) is not about the public interface. A clean API
could wrap either; that part is cosmetic. What matters is the friend algorithms: they reach
into the internals directly for speed (the friend-access decision), so the storage layout
is part of each algorithm's contract. Change the layout and the friends adapt; they cannot
go through a high-level per-element accessor in a hot loop.

**The dividing line is mutability.** VV earns its keep only when the structure is mutable
and edits must stay local: inserting into one column without reflowing every downstream
offset in a flat buffer. Flat wins everywhere the structure is write-once with sizes known
up front: one allocation, contiguous streaming, cheap offsets to hand to BLAS. The cost of
flat is that each column's size must be known before filling; the cost of VV is scattered
heap allocations (a cache miss at every column boundary) plus per-vector overhead.

**Symfactor and static numfactor are flat.** Symbolic never grows: the elimination forest
supplies `frontSize + updateSize` per supernode, so the flat buffer is sized exactly and
each supernode's block is written once with a local cursor. This is what 0.9 does:
`pointerToIndex` from the known sizes, then the union loop fills. A mutable symfactor would be
overkill; nobody patches a symbolic factor in place, one recomputes it after the matrix
structure changes. Static numeric factorization is likewise write-once into a known
structure. Both stay directly comparable to 0.9 buffer-for-buffer.

**Dynamic numfactor is why VV matters.** Dynamic LDL pivoting delays an unstable pivot and
passes it up to an ancestor, so that ancestor's front acquires columns symbolic never
predicted; its index set and value block grow at runtime by an amount unknowable until the
numerics run. That is genuine mutability, arriving through the numerics rather than through
matrix edits, and it is local, independent growth (one ancestor grows, its siblings do
not), which is precisely what VV does cheaply. Flat for dynamic needs contortions:
preallocate every front to its worst case (mostly wasted memory), or reflow the buffer on
each delay (O(nnz) per delay). Flat for dynamic is overkill; VV is the natural structure.

**A hybrid can keep the stored object flat.** The dynamic structure is unpredictable only
*during* elimination; once numfactor finishes, every delay has resolved and the factor is a
fixed structure with known counts, so it is flattenable, we just cannot size the flat
buffer until the numerics have run. The option on the table: dynamic (growable) working
storage inside the numfactor engine, then one flatten pass into the persistent flat factor.
The persistent factor is then uniformly flat whether it came from static or dynamic
numfactor, so the solve phase and the verification against 0.9 are identical either way, and
dynamism becomes an implementation detail of the engine's scratch. The copy is O(nnz(L)),
one pass, noise against numfactor's O(flops), and we factor once and solve many. Two honest
caveats: the flat factor is sized from the numfactor *result*, not from symbolic (delayed
pivoting adds fill, so symbolic stays the static lower bound and the dynamic flat factor is
a separate, post-hoc-sized object); and the flatten canonicalizes, it is copy-and-order
into the fixed front-then-update sorted layout the solve expects, not a pure memcpy. An
incremental flatten, supernode by supernode, releasing each dynamic front as it is copied,
avoids a peak that holds both representations at once.

**The descriptor view unifies the code across both layouts.** A VV is, in essence, an array
of (base pointer, length) pairs, one per column/supernode; flat presents the identical
shape by handing out pointers into its one buffer plus the counts. Write each friend against
that view (`{std::int32_t* ptr; std::size_t len}` for an index run, `{std::int32_t* rowIdx;
Val* val; std::size_t len}` for a matrix column) and the algorithm does not know which
layout produced the base pointer. 0.9 is already in this shape: `pointerToIndex[kk]` is an
absolute pointer into the flat index buffer and the consumer walks `(base, len)`. This is
the friend-access decision generalized from "member access" to "a layout-agnostic view over
members," and it applies to any CSC-style object, the matrix included.

**The view must not be a per-element virtual call. It need not be a template either.**
*(Written before the measurement; the second half of this paragraph is what changed.)* Our
worry was that unifying the layouts would reintroduce the non-inlinable per-element cost the
friend-access experiment clocked at ~6x, so we planned a compile-time abstraction: one algorithm
source, `template<class Storage>`, instantiated per layout. Two monomorphizations, not one binary
branching on layout inside the loop.

The measurement says even that is more machinery than the problem needs. Hoist the descriptors
out of the loop entirely, into three plain arrays (`rowIdxPtr[j]`, `valPtr[j]`, `len[j]`), and
the algorithm becomes an ordinary non-template function taking nothing but pointers. **One**
compiled kernel, no instantiation per layout, no virtual call, and it matches hand-written CSC.
The per-element indirection we feared never arises, because the indirection happens once per
column, not once per entry. So the rule survives in its negative half (never a virtual accessor
in the hot loop) and is superseded in its positive half (a template is sufficient but not
necessary).

The write side unifies too: because sizes are known up front for the flat cases, both flat and a
pre-sized VV inner accept the same `base[cursor++] = i` fill.

**One VV-only hazard: growth invalidates descriptors.** When a VV inner vector reallocates
on growth its base pointer moves, so any cached `{ptr, len}` goes stale. For flat and
symfactor this never happens (read-only after fill), so the view is stable. For dynamic the
rule is grow, then fetch the descriptor, then stream; never hold a descriptor across a
growth event. A `reserve()` per front to a symbolic-size-plus-margin reduces reallocations
but cannot guarantee zero, so the fetch-after-growth discipline still governs. Delays land
at assembly and the BLAS streaming comes after, so the ordering is natural.

**Measured, in `experiments/storage-options/`, and the result is stronger than predicted.**
Two sparse matrix classes (CSC and VV) holding identical content, one `MultiplyEngine`, and a
sparse matvec. Each class fills the same three arrays, `rowIdxPtr[j]`, `valPtr[j]`, `len[j]`.
CSC points into its single buffer; VV reads each inner vector's `data()`. The multiply takes
nothing but those arrays, and on an M4:

```
hand-written CSC (baseline)     1.362 ms
multiply(), CSC pointers        1.454 ms   1.07x
multiply(), VV pointers         1.499 ms   1.10x
multiply(), VV scattered        8.723 ms   6.41x
```

Bit-identical results. Two findings, and the first supersedes what this entry originally
claimed.

**The abstraction is free, and it needs no template.** We had assumed the unifier would have to
be a compile-time template, monomorphized per layout, to avoid a virtual call's cost. It does
not. There is exactly **one compiled multiply**, verified in the symbol table, and its signature
names neither matrix class:

```
T MultiplyEngine::multiply(unsigned long, int const* const*, double const* const*,
                           unsigned long const*, double const*, double*) const
```

It cannot tell CSC from VV because by the time it runs there is nothing left to tell apart. The
storage question and the algorithm question are **separable**: the layout decides where the
pointers come from, and nothing else. Two layouts therefore cost us zero duplicated kernels, not
"two monomorphizations" as this entry first supposed.

**The interface costs the same on both layouts**, 1.07x and 1.10x, a three percent spread. What
costs is *locality*: the two VV rows are the same class with the same content, differing only in
the order their inner vectors were allocated. **A flat buffer guarantees consecutive columns are
adjacent; a vector of vectors only ever borrows that from the allocator.**

**Two caveats on the 6.41x, both important.** It is *constructed*, shuffled allocation plus
interleaved spacers, engineered to remove the allocator's help entirely, so read 1.10x as VV's
structural cost and 6.41x as an upper bound. And it is hardware-dependent: the same code measures
8.87x on a machine with smaller caches. More importantly, **this kernel is the harshest possible
setting for a cache miss**: a sparse matvec does about two flops per element loaded, so a miss
shows at full price with nothing to hide behind.

**Which layout for which object, with reasons rather than symmetry.** Having both layouts is now
free, so the temptation is to offer both everywhere. Most of that would be storage nobody uses.

- **`A`: CSC.** It is built once and read forever, so VV buys it nothing but scatter. A VV `A`
  is defensible only for *structural* mutation, and that case is weaker than it looks. Value-only
  mutation (same pattern, new numbers, refactorize) simply overwrites `mVal` in place and needs
  no VV at all. Structural mutation invalidates the ordering, the forest, the symbolic
  factorization and the factor's size, so it forces the whole analysis phase again, against which
  rebuilding `A` in CSC is one `O(nnz)` pass and therefore noise. Which is also why one *batches*
  structural changes: each forces re-analysis anyway, so nobody applies them one at a time. And
  incremental assembly, the one place VV genuinely helps, has a better answer already: triplets,
  then convert once, with no per-column allocation at all.
- **`SymFactor`: CSC.** Write-once into a size the forest already knows
  (`frontSize + updateSize` per supernode). The textbook case for flat. VV buys it nothing.
- **Static numfactor: CSC**, for the same reason: write-once into a structure symbolic has sized.
- **Dynamic numfactor: VV.** Delayed pivoting grows a front at runtime by an amount symbolic never
  predicted, and the growth is local. This is the one place the algorithm genuinely mutates.

**Flatten-or-not is now a real fork, not a foregone conclusion.** This entry originally treated
the flatten-to-flat hybrid as settled. At 1.10x for a packed VV, "do not flatten; let the solve
read the VV factor" is defensible, and it saves an `O(nnz(L))` copy plus the peak memory of
holding both representations at once. Against that: the solve *streams* the factor, which is
exactly the cache-hostile shape the experiment measures, so it is the phase most likely to want
flat. Both options are now open, and the experiment is what opened them.

**The measurement we owe, and it is load-bearing.** The matvec is a toy for this purpose: about
two flops per element loaded, so a cache miss shows at full price. Numeric factorization is the
opposite, a dense front handed to BLAS level 3, `O(n^3)` arithmetic on `O(n^2)` data, where the
same miss amortizes over far more work. **We have predicted that VV-during-elimination is
therefore affordable, and we have not measured it.** That prediction is what the entire dynamic
design rests on, and it deserves its own study against a realistic front, not a sparse column
scatter. Until then, treat "VV is affordable in numfactor" as a hypothesis, not a result.

(Note also that the ~6x virtual-call figure this entry cites was carried over from the
friend-access experiment and is *not* tested by storage-options, which has no virtual path. The
question turned out not to arise, since the unifier needed no polymorphism at all.)

**Status.** Settled: CSC for `A`, for symfactor, and for static numfactor. Settled: VV for dynamic
numfactor's working storage. Open: whether the persistent factor is flattened, and whether a
VV-all-the-way numfactor is viable, both of which want the BLAS-3 measurement above.

**Confirmed against 0.9 (2026-07-13).** This entry predicted the split and flagged it as
unverified, pending sight of 0.9's numeric code. That code has now arrived, and 0.9 does exactly
what we reasoned it must:

- **`FactorsStatic`** sums `numberOfAllocatedEntries` over every supernode, allocates **one
  buffer**, and points `pointerToEntry[jj]` into it. **Flat.**
- **`FactorsDynamic`** runs `for (jj) pointerToEntryArray[jj] = new Real[jjEntrySize];`. **One
  allocation per supernode. A vector of vectors in all but name.**

So the split is a *port*, not a modernization: 0.9 already separates the two storages by
mutability, for the same reason we would. The prediction and the oracle agree, and the entry's
provisional status on this point is closed. What remains open is the *hybrid* question, whether
the dynamic factor is flattened afterwards, which 0.9 answers only implicitly (it keeps a
`FactorsDynamic` object, and `SolveEngine` reads it through the same abstract base), and the
BLAS-3 measurement that would tell us whether flattening is worth it.

**Amendment (2026-07-14): the taxonomy stands; the descriptor view in the title did not.** The
storage split above (CSC for the matrix, the symbolic factor, and the static numeric factor; a
vector of vectors for the dynamic factor) is the settled, 0.9-confirmed core and is unchanged. What
did not survive is the mechanism this entry pairs with it, the "descriptor view that spans both": a
bulk extractor that materializes `rowIdxPtr[j]` / `valPtr[j]` / `len[j]` arrays and feeds one
storage-blind compiled multiply. That was superseded by **direct access**. Each storage now exposes
per-column accessors (`rowIdx(j)` / `val(j)` / `colSize(j)`, one pointer or size per call), a
consumer templates over the storage and calls them at the moment of use, and there is no extractor
at all, so nothing is materialized, owned, or left to go stale. Direct is also faster (it streams no
extra arrays) and it is what the numeric factorization must use regardless, since a growing dynamic
factor would dangle any pointer extracted up front. So the array-valued names in this entry describe
the retired bulk form; the current interface is the three per-column accessors. See the
bulk-versus-direct entry above, and the storage-options README.

## 2026-07-09, Index types: `std::int32_t` IDs (NIL = -1), `std::size_t` offsets

Two kinds of integer, two types:

- **IDs**, a value that *names* a vertex/row/column/supernode, and may need a "none"
  marker, are **`std::int32_t`**, with sentinel **`NIL = -1`** (`constexpr std::int32_t`).
  E.g. `SparseMatrix::rowIdx`, the permutation maps, and the forest's parent / child /
  sibling / supernode-map arrays.
- **Offsets / counts / sizes**, row-pointers, `nnz`, dimensions, anything that indexes
  into or measures, are **`std::size_t`**. E.g. `SparseMatrix::colPtr`, `size()`.

**Why signed int32 for IDs.** The forcing function was the sentinel. The forest needs a
"no parent"/"no child" marker; on an unsigned `std::size_t` that can only be the max
value, spelled `static_cast<std::size_t>(-1)`, defined behaviour but an ugly wraparound
smell, and easy to misuse in arithmetic. A signed `std::int32_t` gives a clean, obvious
`-1`. This also matches (a) the graph code, which already uses `int32_t` vertex IDs with
`static const int NIL = -1;` so companion arrays like `mate` hold `-1` naturally, and (b)
the vendored AMD/MMD, which are `int`-based, so `rowIdx` as `int32_t` largely removes the
`size_t→int` conversion at that boundary. The cost is a ~2.1-billion index cap (int32 vs
size_t's range). Accepted deliberately: cleaner and more agile, and Oblio isn't targeting
matrices past 2^31 structural indices for now. If that changes, widen the ID type to
`std::int64_t` in one place.

**Why `size_t` stays for positions.** A position is never negative and can legitimately exceed
2^31 even when the index *count* does not: a row-pointer is a position into an `nnz`-length
array, and `nnz` outgrows `n`. So positions keep the full unsigned range and never carry a
sentinel. This mirrors the graph's `GeneralGraph` exactly: `idx` (the row-pointer) is `size_t`,
`adj` (the neighbour indices) is `int32_t`.

(CODING_RULES now calls these **indices** and **positions** rather than IDs and offsets, and
bans the phrase "index into a vector": `index` names a matrix or graph entity, and reusing it
for array access collides with that meaning exactly where the distinction matters. Same
distinction, sharper words.)

**No signed/unsigned friction, because loop counters stay `std::size_t`.** g2_csr
demonstrates the discipline: counters that enumerate positions are `size_t` (so they
compare against `.size()` cleanly), and `int32_t` appears only as *stored values* that
may be `NIL`, never as a loop variable. `size_t` is a safe superset of the non-negative
`int32_t` range, so viewing an index as `size_t` in a bounded, non-negative loop loses
nothing. The two types reconcile with explicit casts at exactly two crossings:
`static_cast<std::int32_t>(counter)` when storing a counter into an ID array (the
narrowing where the 2^31 cap lives), and `static_cast<std::size_t>(id)` when an ID
subscripts an array (widening; guard against `NIL` first if it could be a sentinel). Casts
are few and mark the ID↔offset boundary.

Spelling: `std::int32_t` from `<cstdint>`, `std::size_t` from `<cstddef>`, std-qualified,
C++ headers, same rule as every other stdlib type. Applied first to `SparseMatrix`
(`rowIdx`→`int32_t`, `colPtr`→`size_t`); `Permutation` (maps→`int32_t`) and `ElmForest`
(parent/etc.→`int32_t` with a shared `NIL`) follow. The graph code uses bare `int32_t`;
it's the code to bring in line later, like the matching codebase was for `size_t`.

**Later note (2026-07-15):** the "~2.1-billion index cap" above is really two separate ceilings that
happen to meet at `INT32_MAX`, a dimension/index cap (representability, intrinsic to `int32`) and an
nnz(A) cap (the `int`-based ordering interface, external and A-only). nnz(L) is not capped. See the
ordering-constraint entry for the distinction.

## 2026-07-09, Friend grants write access; reads are public

The engine↔data access rule:

- **`friend` = write access.** An engine befriends exactly the data class(es) it
  *produces/mutates*. There is no public mutation API for those internals, only the
  producing engine reaches the private members. E.g. `OrderEngine` writes
  `Permutation`; `ElmForestEngine` writes `ElmForest`. Friendship is declared by the
  data class (`friend class FooEngine;`), not by the engine.
- **Reads are public, everywhere, including hot paths.** Engines, tests, and users
  all read via the public const API. `friend` is *not* needed for reading.
- **This corrects an earlier overstatement** (the friend/BLAS entry below implied
  friend was needed for hot-path reads). The `experiments/friend-access/` study
  measured a *per-element, cross-translation-unit accessor call*, which can't inline,
  so it blocks vectorization, against direct friend member access; friend won ~6×.
  But that gap comes from calling the accessor *per element*, not from using a public
  accessor at all. Two public-accessor patterns match friend's performance exactly:
  (a) hand BLAS the block via a single `.data()` call, BLAS then owns the O(n³) loop,
  so the one call is free; (b) for a hand-written element loop, bind the returned
  container once (`const auto& v = A.val();`) and loop over *that*, one non-inlined
  call total, then a vectorizable loop over contiguous memory. So the hot-path
  discipline is **"bind-once / pass the pointer," not "be a friend."**
- **Consequence for `SparseMatrix`:** A is input and nothing writes it (its
  construction path is TBD), so it needs *no* friends. Its earlier
  `friend class OrderEngine` is removed; `OrderEngine` reads A via
  `colPtr()`/`rowIdx()`/`size()` and remains a friend only of `Permutation` (which
  it writes). `ElmForestEngine` already followed this (friend of `ElmForest` only).
- **Numeric data classes** (`Factors`, …) still befriend their producing engines for
  *writes*; their hot-path reads also go through public accessors with the
  bind-once/pointer discipline.
- **Exposure stance (pragmatic, not purist):** exposing internals read-only is fine;
  we don't design curated "won't-break" read APIs up front. A representation change
  already forces editing the friends (writers); public read exposure adds only the
  *tests* to that blast radius, cheap, and tests *should* feel a rep change. Curate a
  narrower public API later only if a structure's representation proves unstable. For
  the canonical structures we have (CSC, etree `parent[]`), the representation is the
  settled standard form, so exposing it by reference is low-risk.

## 2026-07-09, Sparse matrix storage: flat CSC, stored FULLY (both triangles)

`Matrix` (input A) stores its structure and values as flat **compressed sparse column
(CSC)**: three contiguous `std::vector`s, `colPtr` (size+1), `rowIdx` (nnz), `val`
(nnz), row indices sorted ascending per column. A symmetric matrix is stored **fully
(both triangles)**, each column holding its complete neighbour list plus the diagonal.

Storage layout (flat CSC vs vector-of-vectors) and triangle (full vs lower) are two
separate decisions:

**Layout, flat CSC**, across the three generations:
- **0.9** stored CSC via four manually `new`/`delete`d `Array*` pointers, the manual
  memory management the port removes.
- **10.12** modernized to `std::vector` but chose **vector-of-vectors** (one inner
  vector per column), RAII but the wrong layout: columns scatter across the heap.
- **PoC / port** use flat CSC (`mColPtr`/`mRowIdx`/`mVal`), vectors *and* contiguous,
  satisfying the "contiguous storage to BLAS via `.data()`" invariant.

**Triangle, full, not lower.** Both 0.9 and 10.12 store A **fully**: 0.9's
`getNumberOfNonzeroEntries() = size + 2*numOffDiagonals` (each off-diagonal in both
triangles) and its "storing A within the structure of A+Aᵀ"; 10.12 has `SymmetrizeStrc`
and its etree reads full per-column neighbour lists. The **PoC diverged**, storing the
lower triangle only ("stored as lower triangle in CSC"). That divergence is what forced
every structural consumer to expand lower→full first (the MMD path in `OrderEngine`, the
etree in `ElmForestEngine`), and the etree bug where lower-triangle input silently
produced an empty tree until expansion was added. **The port matches the oracle: A is
stored fully.** Consequences: structural phases (ordering, elimination forest, symbolic)
read each column's neighbours directly with no expansion (the etree's diagonal
self-skips via `lc1 < lc2`; MMD just strips the diagonal; AMD ignores it); it's the
faithful port (lower-triangle was a rewrite of the data structure); and it's the natural
substrate for a future **unsymmetric extension**, factor the symmetrized structure
A+Aᵀ while carrying asymmetric values. Cost: ~2× off-diagonal storage for A, and the
numeric phase carries a redundant triangle. Accepted, matching 0.9/10.12; A is the input
and is far smaller than the factors, where the real memory lives.

Open for the port: the PoC exposed this as a public `struct` with a `fromCOO` builder
and a weak `isValid`. The modern `SparseMatrix` keeps the flat-CSC layout but is a
`class` with `friend` engines and a structural interface; 0.9 is the oracle for the
COO→CSC assembly details (zero-diagonal insertion, duplicate merging, symmetrization),
10.12 shows which operations the solver actually calls.

## 2026-07-09, Two layers of modernization: rules prevent, clang-tidy catches

The coding rules and `.clang-tidy` are complementary layers catching different
failures at different times, not redundant work:
- **Coding rules** (CODING_RULES + CLAUDE invariants), preventive and broad. They
  shape code as it's written and cover what no tool can judge: port-verbatim
  discipline, friend→BLAS, when to split a header, `.cpp`, `mFoo`, `std::size_t`.
- **`.clang-tidy` `modernize-*`**, a mechanical safety net, narrow but certain. It
  can't reason about intent, but within scope it catches the idiom slips that get
  through (a stray `NULL`, a `typedef`) and can auto-fix them.

So a `NULL` written by mistake is exactly what the rules *say* and the tool
*guarantees*, belt and braces, each doing what the other can't.

`modernize-*` is aligned with this project's purpose at the concept level: porting 0.9
(late-90s C++) forward *is* turning `typedef`→`using`, `NULL`→`nullptr`, raw loops→
range-for, and so on across old code, exactly what that checkset does. So it can do
part of the mechanical modernization *for* you on each ported file (`--fix`), leaving
your attention on the algorithmic faithfulness no tool can verify. Hence the per-unit
workflow (in CLAUDE.md Process): port faithfully → `clang-tidy --fix` → verify vs 0.9.

## 2026-07-08, `experiments/` convention (runnable design studies)

`experiments/<name>/` holds self-contained, runnable studies that establish or
validate a coding standard before it is applied in the main tree. Each is its own
folder with its own sources, `Makefile`, and `README.md`; builds standalone
(`make test`); and is reference/teaching material, **not** part of the main Oblio
build. Executables carry the `_cpp` suffix and are gitignored.

Distinct from its two neighbors:
- `archive/`, frozen history (superseded PoC devlog, 0.9-analysis notes, old
  harnesses). Not maintained, not built.
- `examples/`, usage samples showing how to *call* the library (`examples/basic.cpp`).

An experiment answers a design question with code you can run and measure, then feeds
a decision here. Current studies: `template-instantiation/` (how to instantiate the
`Val` template, implicit vs plain/guarded explicit) and `friend-access/` (public API
vs `friend`-direct access, with timing). Experiments use the already-settled standards
(guarded explicit, `.cpp`, `mFoo`, `Oblio` namespace), so they double as worked
references for those standards.

## 2026-07-08, Numeric hot path: `friend` access, then BLAS (carried from 0.9)

The 0.9 design for numeric work, which the port preserves: **an engine reaches the
data's contiguous block via `friend`, then hands that raw block to BLAS** wherever
BLAS applies (gemv/gemm/syrk/trsm/potrf, via Accelerate on macOS). Not a new choice,
this two-step (`friend` → BLAS) is how 0.9 does dense numerics.

This *is* supernodal numeric factorization: supernodes are dense blocks embedded in
the sparse structure, and the numeric phase is a long sequence of dense BLAS calls on
them, `syrk`/`gemm` for Schur-complement updates, `potrf`/`getrf` to factor the pivot
block, `trsm` for the off-diagonal solve, repeated per supernode across the whole
elimination. `FactorEngine` reaches each supernode's contiguous storage via `friend`
and passes the pointer straight to BLAS, no copy, thousands of times per factorization.
So `friend` isn't an optimization detail; it's the access mechanism the entire numeric
phase is built on. The `experiments/friend-access/` mat-vec is the single-block toy of
this pattern.

Important distinction, **the supernode blocks live in the factors, not in A.** `A`
(the input `Matrix`) is never handed to BLAS block-by-block; it's *read*, its structure
by the ordering and symbolic phases, its values once when they're scattered into the
factor. The dense blocks that BLAS operates on are created during factorization and
stored in `Factors`. So the `friend`→BLAS hot path is specifically a `Factors` /
`FactorEngine` story. `A`'s storage (CSC) is chosen for a different reason: cheap,
cache-friendly *sequential column traversal* by the structural phases, plus being the
standard interchange format (what AMD/MMD expect). Both `A` and the factors favor flat,
contiguous storage over 10.12's vector-of-vectors, but for `A` the reason is streaming
structural reads, and for the factors it's the contiguous block BLAS needs.

Data classes (`Matrix`, `Vector`, `Factors`, `Symbolic`) expose a public,
bounds-checked API (`operator()`, `operator[]`) for reads, by all callers, on hot
paths too (see the "friend = write access" entry above: reads are public; the hot-path
discipline is bind-once / pass `.data()` to BLAS, not friendship). Engines befriend a
data class only to *write* it. 0.9 grants `FactorEngine` friend access into `Matrix`,
`Factors`, `Symbolic`; in the port that friendship is retained only where the engine
writes (the factor storage), and reads go through the public API.

Why (performance): the public-operator path is one non-inlined, cross-translation-unit
call per element (data-class body in its `.cpp`, loop in the engine's), which blocks
vectorization. Direct `friend` access fetches the raw block pointer once and walks
contiguous memory, which vectorizes, measured ~6× on Apple Silicon (M4/AppleClang),
~3× on x86/g++, over the API path. But the raw block should then go to **BLAS**, not a
hand loop: on the M4, Accelerate's `dgemv` ran the 2000×2000 mat-vec ~6× faster than
the vectorized hand loop (and ~36× over the API path), because it breaks past a single
core's bandwidth ceiling (multithreading + prefetch/blocking) that a hand loop can't.
The advantage only grows for the compute-bound O(n³) kernels (`dgemm`/`syrk`/`trsm`/
`dpotrf`) the factorization leans on. So the hand loop is a fallback/baseline; where
BLAS applies, use it.

Why `friend`, not public getters: `friend` is *tighter* encapsulation, not looser, it
grants access to exactly the named engine classes, where a public `data()`/getter
exposes internals to the whole program. It honestly encodes that a data class and its
engines are one subsystem split for organization, not two modules talking through a
narrow API. Deliberate pragmatic choice over OO-purist accessors, and it's faster.

Consequence for porting: `friend` couples a data class to its engines, so the natural
port/verify unit is the *cluster* (e.g. `Factors` + `FactorEngine`), not the data class
in isolation, the friend boundary sets porting granularity.

Measured note: the gap is structural (non-inlined calls / vectorization), not assertion
overhead, toggling bounds-check asserts barely moves it. See `experiments/friend-access/`.

## 2026-07-08, Source extension: `.cpp` (headers stay `.h`)

Switch Oblio source files from `.cc` to `.cpp`; headers remain `.h`. All
extensions (`.cc`, `.cpp`, `.cxx`, `.C`) are identical to the compiler, so this is
convention, not correctness.

Why `.cpp`: cross-project consistency. The matching codebase uses `.cpp` (alongside
`.rs`, `.py`), one extension per language across the ecosystem, so Oblio being
`.cc` made it the odd one out. `.cpp` is also the more common choice in the wider
world (`.cc` is mainly the Google-style corner). Done now because the tree is still
scaffold + PoC with no ported units yet, so the rename is at its cheapest.

Blast radius (all mechanical, filename-level, no semantics): `git mv *.cc → *.cpp`
across `src/`, `tests/`, `examples/`, `archive/`, `experiments/`; the CMake source and
executable lists; the manual build glob in CLAUDE.md / README (`src/*.cpp`); `.clang-format`;
the example `Makefile`s and the build-command comments in the example files. Headers and
object files are untouched, `#include`s point at `.h`, and `Foo.o` derives from the
source basename regardless of extension, so no `.o` reference changes. This makes
the rename strictly safer than the `exp`→`ext` one, which touched `#include`s.

## 2026-07-08, Explicit instantiation over header-only templates (rationale)

Decision (already active in CLAUDE.md; this entry records *why*, which otherwise
lives only in the `experiments/template-instantiation/` example comments): Val-dependent classes keep a single
`Val` template, but their definitions live in `.cpp` files with explicit
instantiation for the supported scalar types, and headers carry declarations plus
`extern template`. Fuller treatments: `archive/oblio_modernization_notes.md` §"Why
explicit instantiation still works" (the Val-surface table) and
`archive/oblio-new-devlog.md` Session 3 (adoption + the link-failure proof). The
`_tpl` / `_ext` example files are the compact head-to-head (see naming below).

Mental model: a template is a recipe, not code, it generates code once a type is
plugged in. Header-only templates plug in *late and everywhere*: every translation
unit that includes the header re-runs the recipe for each type it uses, and the
linker discards the duplicates (N files × 2 types → the same bodies compiled 2N
times). That was the real cost of 0.9's header-heavy templating, not the templates,
but the repeated late instantiation. Explicit instantiation is "one template,
applied early, once per type": the recipe runs exactly twice, in one `.cpp`, at
library-build time, and every other file links the existing result instead of
re-running the recipe. (Two mechanisms achieve that "instead of re-running",
declaration-only headers and/or `extern template`, see History below; they are
not the same feature and have different dates.) Generality is preserved (adding `float` is one
early application), but instantiation collapses from scattered-2N to centralized-2.

Key framing: build cost that scales with the number of scalar types is *incidental*
to the C++ compilation model, not *inherent* to supporting real and complex.
Nothing about "a matrix can hold real or complex" requires recompiling matrix code
in every includer, that only happens because header-only defers instantiation to
include time. Explicit instantiation removes the accident and keeps the capability.

The tradeoff, and why it's nearly free here: explicit instantiation gives up
instantiating *arbitrary* types at use sites, a consumer can't spin up
`Matrix<long double>` unless that line exists in the `.cpp`. For a maximally-generic
header library (Eigen) that's a real loss. For Oblio it isn't, because the scalar
world is closed and tiny: a type only makes sense if a dense BLAS/LAPACK kernel
exists for it, which bounds the space to BLAS's four, `float`, `double`,
`complex<float>`, `complex<double>` (s/d/c/z). We use two, might add the others.
Enumerating even all four is a handful of lines, far below the build cost of keeping
them implicit. So the one thing explicit instantiation costs is something a
closed-world numerical solver doesn't want anyway.

Bonus, and it matters for a verification-focused port: the `template class Foo<...>;`
lines *are* the list of supported types, in one place. Support becomes a declared,
reviewable fact rather than an emergent property of whatever anyone happened to
instantiate, and adding a type is a deliberate act that forces confronting whether
kernels and tests exist for it, exactly the gap the appendix flagged when complex
was "a new code path with zero test coverage."

Three mechanisms (not two, this frame makes the history click). Three distinct
tools, three natures, three dates:
- **Forcing**, `template class Foo<double>;`, emit all of `Foo<double>` in this
  TU. C++98.
- **Suppressing**, `extern template class Foo<double>;`, do *not* implicitly
  instantiate here; link it from elsewhere. C++11 (GCC extension earlier).
- **Definition-hiding**, not a keyword but a code-organization move: put member
  bodies in a `.cpp`, leave declarations in the header. A TU that can't *see* a body
  can't implicitly instantiate it. Works in every era.

Definition-hiding is the hinge, and it pairs with forcing. The three configurations:

| Case | Available | Approach | Build cost |
|---|---|---|---|
| 1 | neither forcing nor suppressing | Inclusion model: all definitions in headers, every TU re-instantiates what it uses, linker merges duplicates. No way to move bodies out and still get symbols. | High (~2N), unavoidable |
| 2 | forcing only (C++98) | Definition-hiding + forcing: bodies to `.cpp`, declaration-only headers, `template class Foo<double>;` in the `.cpp`. Other TUs see declarations only → can't implicitly instantiate → link the forced symbols. | Low |
| 3 | forcing + suppressing (C++11) | Case 2 still works and stays the choice; *additionally* you may keep bodies in headers and use `extern template` to suppress re-instantiation, needed only when definitions must stay header-visible. | Low |

Key insight: the big jump is **1 → 2, not 2 → 3**. Forcing is what unlocks the whole
technique (hide a definition, still guarantee the symbol). Suppressing is the
incremental step that only adds a second route for the case where you insist on
header-visible definitions. If you move bodies to the `.cpp` (Oblio does), you never
need it.

Precondition for all of it: an **enumerable** type set. For genuinely arbitrary
types you can't force or suppress anything (you don't know the list), so you're in
Case 1 regardless of language version. Forcing/suppressing are tools for closed type
sets, which Oblio's is (BLAS s/d/c/z), the same fact that makes the tradeoff above
nearly free.

Dates and 0.9: forcing is C++98, suppressing is C++11, so when 0.9 was written (late
90s) only Case 1 was portably reachable, header-only, inclusion model. That was the
*correct* choice for the era, and not because suppressing was missing (Case 2 doesn't
need it) but because template separate-compilation was a portability minefield then
(the `export` saga, inconsistent two-phase lookup, compilers disagreeing on
inclusion-vs-separation). Header-only-everything was the safe default. The modern
refactor applies matured portability; it does not correct a 0.9 error.

Where Oblio sits: current `ext` code is **Case 3**, bodies in `.cpp` (declaration-only
headers) *plus* `extern template`. But because the headers are already
declaration-only, the build win is really Case 2's (definition-hiding + forcing); the
`extern template` lines suppress nothing here (no visible header body to instantiate),
so they are documentation, not mechanism, a header annotation of intent, latent
unless a body is later (wrongly) added to a header. See the naming note below for the
full implicit / plain explicit / guarded explicit framing. So Oblio's pattern was achievable in C++98; C++11 was
not strictly required.

Template-instantiation example, naming (one algorithm, dense mat-vec, built
three ways). File names are `<Class><Variant>` with no separator, e.g.
`MatrixImplicit.h`, `MatrixPlainExplicit.{h,cpp}`, `MatrixGuardedExplicit.{h,cpp}`
(same for `Vector`, `MultiplyEngine`). The three variants:
- `Implicit`, body in the header; instantiated implicitly per TU. Stands in for
  what 0.9 effectively was.
- `PlainExplicit`, bodies in the `.cpp`, header signatures only (explicit
  instantiation, forcing only).
- `GuardedExplicit`, plain explicit + `extern template` in the header. The
  pattern used in the real tree.

Conceptual framing (two axes):

- **Axis 1, where the body lives (implicit vs explicit).** *implicit* = body in
  the header (`.h`), instantiated implicitly per translation unit → `Implicit`.
  *explicit* = body outside the header, forced in the `.cpp`; header carries
  signatures only → `PlainExplicit` and `GuardedExplicit`.
- **Axis 2, applies only within explicit; guarded vs plain.** *plain explicit* =
  `.cpp` bodies, declaration-only header, nothing more → `PlainExplicit`.
  *guarded explicit* = same, plus `extern template` in the header → `GuardedExplicit`.

  The three named layers:
  - **implicit**, body in `.h` (`Implicit`)
  - **plain explicit**, bodies in `.cpp`, header signatures only (`PlainExplicit`)
  - **guarded explicit**, same as plain explicit, plus the guard (`GuardedExplicit`)

  Important: the guard is NOT a sub-kind of "more correct" explicit, plain explicit
  is a complete, valid design. `extern template` only ever acts on *visible* header
  bodies (the implicit failure mode); in a declaration-only design there is nothing
  for it to suppress, so it is pure documentation, a header annotation reading
  "instantiated elsewhere," aimed back at the implicit branch it guards against. It
  gains mechanical effect only if someone reintroduces a header body (which the
  "definitions live in `.cpp`" invariant in CLAUDE.md forbids). So "guarded" (a
  reminder/guard), not "suppressed" or "enforced", in this design it suppresses
  nothing and enforces nothing; the invariant does the enforcing.

All three are built and tested together via the example's `Makefile` (`make test`)
against one shared source, `test_multiply.cpp`, they must produce identical
results, and the plain explicit and guarded explicit variants share the same link-failure behaviour when
their `.cpp` files are omitted (empirical confirmation that with declaration-only
headers `extern template` suppresses nothing, it is documentation, not mechanism).
Selector macros: `OBLIO_TI_IMPLICIT` / `OBLIO_TI_PLAIN_EXPLICIT` /
`OBLIO_TI_GUARDED_EXPLICIT`.

Naming history: suffixes were once `_tpl`/`_exp`/`_ext`, where `_exp` inaccurately
labeled the extern-template variant. Renamed in two steps, first `_exp`→`_ext` to
free `_exp` for the genuine forcing-only variant, then all three to the conceptual
`Implicit`/`PlainExplicit`/`GuardedExplicit` once the two-axis framing settled.

## 2026-07-08, Matrix naming: explicit `SparseMatrix` / `DenseMatrix`

Rename the sparse matrix type from `Matrix` to **`SparseMatrix`**; keep
**`DenseMatrix`**. Both are plain concrete types with **no shared base class**.

Why: the old `Matrix` (implicitly sparse) + `DenseMatrix` (explicitly marked) is a
half-committed convention, it privileges sparse as the unmarked default but marks
dense as the exception, so a reader must already know "unmarked = sparse" to parse
it. Going fully explicit removes that. It also makes matrix naming consistent with
the graph naming used elsewhere (`GeneralGraph` / `BipartiteGraph`, no bare
`Graph`). Note the typed-library mainstream (Eigen, scipy) actually defaults the
*other* way (unmarked = dense, `SparseMatrix` marked), so `Matrix` = sparse
actively misleads anyone arriving from there. The sparse-first precedent
(CSparse/SuiteSparse `cs`) is the one respectable counter-argument, but explicit
wins on cross-codebase consistency and reader-independence.

No speculative base `Matrix` interface. If something ever needs to consume sparse
and dense polymorphically, add the interface then, on top of the two concrete
types, not before a caller forces it.

Scope: this is a rename (wrapping/API), not an algorithm change, it's on the
port-and-modernize track, not the rewrite track. Do it as one deliberate
mechanical pass **before** porting proper, since it's cross-cutting (every `friend`
decl, FactorEngine, SolveEngine, tests name `Matrix`) and only gets more expensive
as code solidifies around the name. **Oracle mapping: 0.9 `Matrix` ↔ modern
`SparseMatrix`**, record this so output comparisons against 0.9 stay unambiguous.

## 2026-07-08, Minimal abstraction; containers are the structure that matters

Design stance for the port: concrete types, minimal OO ceremony, `std::vector` as
the spine. Don't build base classes, single-implementation interfaces, or
inheritance ahead of an actual polymorphic caller (YAGNI). A direct solver is a
pipeline of concrete transforms, not a class hierarchy.

One guard against a tempting over-correction: "AI makes code cheap to generate, so
structure matters less" is only half true. Generation got cheap; **verification did
not**, and this project is the proof (every PoC bug was cheap to write, expensive
to trust, and fixed only by slow comparison against the 0.9 oracle). Structure's
real job was never to save typing; it's to keep code readable and checkable and to
localize where a bug can hide. That matters *more* when code is machine-generated,
because the bottleneck shifts onto review. So: drop OO flavor that only served
human writing time; keep the structure that serves verification (clear module
boundaries, one concern per unit, the friend seams that let a supernode block be
diffed against 0.9).

This is why "proper containers everywhere" (`std::vector` over `Array`) is not a
style preference but the load-bearing invariant: a vector carries its own size (no
drifting length variables), self-frees (no leak/double-free surface), bounds-checks
under sanitizers, and hands clean pointers to BLAS via `.data()`. Each removes a
bug class you'd otherwise verify the absence of by hand. For this codebase, the
container discipline *is* the architecture that matters.

## 2026-07-07, Align with standard project files; adopt clang tooling

The doc set maps to established conventions rather than being bespoke:
CODING_RULES.md ≈ a style/conventions guide, DESIGN_DECISIONS.md ≈ a lightweight ADR log
(single-file variant of the Nygard/ADR pattern), CLAUDE.md ≈ the agent-instructions
file (AGENTS.md is the emerging cross-tool equivalent, if we ever run more than
Claude Code, keep content in AGENTS.md and make CLAUDE.md a one-line `@AGENTS.md`).
PORTING_LEDGER.md stays bespoke, it's specific to a port, no standard analog.

Renamed CPP_RULES.md → **CODING_RULES.md** and made it language-general (Rust is a
likely future scope), with per-language sections.

Adopted **`.clang-format`** and **`.clang-tidy`** so mechanical rules (`nullptr`,
`using`, `enum class`, `= delete`, narrowing) are tool-enforced instead of written
as prose an agent may skip. `.clang-tidy`'s `modernize-*` family does part of the
port mechanically. The prose rules doc now holds only judgment calls tools can't
check.

Stubbed ahead of need (deliberate, not yet load-bearing): CONTRIBUTING.md and
CHANGELOG.md ("Keep a Changelog" format). CONTRIBUTING fills out on going public
(and becomes the canonical build/test source then); CHANGELOG gets its first real
entry at the first tagged release, which requires settling the tree's version
identity (Oblio 11 vs a fresh 1.0), tied to the working-tree question below. Still
outstanding: LICENSE (bundled AMD is BSD-3-clause, so a project license decision is
eventually needed).

Tension worth remembering: canonical ADR is many small numbered files, which suits
humans browsing history but fights the "always in context" goal, a growing pile
can't all stay loaded. The single-file log + distilled always-on layer is a
deliberate adaptation of the standard, not an oversight.

## 2026-07-07, What stays always in context vs on demand

Goal: coding rules and active design constraints present in every session, without
letting context grow unbounded.

Mechanism: only CLAUDE.md and its `@`-imports load every session; everything else
is read on demand. `@import` costs the same context as inlining, it just keeps
the file separate and editable.

Split:
- **Always-on** (in CLAUDE.md or imported): invariants (inline), conventions
  (`@docs/CODING_RULES.md`), and a distilled *Active design constraints* summary.
- **On demand**: the full DECISIONS log (this file, grows over time, so importing
  it would erode context), PORTING_LEDGER (porting-specific; read after a gap),
  README, archive history.

Rule of thumb: the always-on set is the distilled essence; narrative, dates,
history, and open questions stay here. Keep CLAUDE.md well under ~200 lines, past
that, adherence drops. If the always-on set grows, distill harder rather than
importing the growing logs.

## 2026-07-07, Invariants live in CLAUDE.md, not CODING_RULES.md

The C++ **invariants** (port-verbatim, `std::vector` default, no signed→unsigned
index slips, `.data()` to BLAS) are written directly in CLAUDE.md. Only the
**conventions** (style preferences) stay in CODING_RULES.md.

Why: Claude Code auto-loads CLAUDE.md every session (directory walk up to repo
root, concatenated) but does **not** auto-load other files, a pointer from
CLAUDE.md to CODING_RULES.md does not pull its contents in. So a rule meant to be
active every session must physically live in CLAUDE.md, or it isn't loaded until
Claude happens to read the file it's in. Invariants are always-on; conventions
are fine to read on demand.

Do not "consolidate" the invariants back into CODING_RULES.md, that silently
disables them.

## 2026-07-07, Documentation structure

Established this doc set: CLAUDE.md (operating contract + index),
PORTING_LEDGER.md (per-unit status), CODING_RULES.md (conventions a linter can't
enforce), DESIGN_DECISIONS.md (this log), plus `.clang-format` / `.clang-tidy` for
tool-enforced mechanical style. The existing md files, devlog, modernization
notes, appendix, README, are kept as history/rationale for now; decide later
whether to fence or retire them.

Why: after a multi-month gap the project needs to be reconstructable from three
things (history, layout, decisions) without re-deriving them. CLAUDE.md holds
the operational index (it's what Claude Code loads each session); the reasoning
lives here so CLAUDE.md stays lean. The two live in different files on purpose and
should not duplicate each other.

(This entry documents a list that includes DESIGN_DECISIONS.md itself. The self-reference
is intentional, not an oversight.)

## 2026-07-07, OPEN: which tree is the working copy

Unresolved; resolve before the first port step and record the outcome here.

Candidates:
- **(a) Fresh port from 0.9** into a new tree, treating the PoC the way 10.12 is
  treated, learn from, don't build on. "Replacing Array gradually" only applies
  to a tree that still has `Array` (0.9 or 10.12), which points here.
- **(b) Continue the PoC tree**, where Array→vector is already done, in which
  case the remaining work is different (coverage, finishing, cleanup), and this
  is not really an Array migration.

The stated plan ("port carefully from 0.9, one function at a time, replace Array
gradually") reads as (a). Confirm.

## 2026-03-07 (PoC), Choices carried from the proof of concept

Recorded for continuity. The PoC was exploratory, so revisit each on its merits
rather than inheriting it unquestioned.

- **One `Val` template** instead of 0.9's separate `*Real.h` / `*Complex.h` file
  pairs. Cholesky treated as Hermitian input; LDL^T as complex-symmetric.
- **`std::vector` storage** instead of the hand-rolled `Array`.
- **Explicit template instantiation** for `double` and `std::complex<double>`,
  headers declare, `.cpp` files define and instantiate. Faster builds; `float` /
  `long double` remain one line away.
- **Namespaced `include/oblio/` headers**, declarations only.
- **Flat `src/`**, all sources, including `Mmd.cpp` / `Amd.cpp`, directly in `src/`.
