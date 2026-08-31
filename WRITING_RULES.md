# Writing Rules

Conventions for prose and documentation, not enforced by tooling. For code, the
mechanical style layer is `.clang-format` and `.clang-tidy`; this file is the same
softer layer for narrative, holding conventions for consistency across the prose
deliverables: README, CHANGELOG, CONTRIBUTING, DESIGN_DECISIONS, PORTING_LEDGER,
the devlog, and educational docs such as `sparse_factorization.md`.

**Invariants** (must not break) live in CLAUDE.md, always loaded. This file is the
softer layer: conventions for consistency, not correctness. It is the prose
counterpart to CODING_RULES.md.

## Prose conventions (all documents)

- **No em-dashes in prose. Use commas.** Never the em-dash character, in any
  document, at any time. A comma does the job; where a stronger break is wanted, a
  colon, parentheses, or a new sentence. This is a hard rule, not a preference.
- **ASCII only.** Every character in every document is plain ASCII: prose,
  diagrams, matrices, pseudocode, and math alike. No smart quotes, en-dashes or
  em-dashes, arrows, Greek letters, or other Unicode. Write `->` for arrows,
  `alpha` for the symbol, `<=` for inequalities.
- **American spelling.** `neighbor`, not `neighbour`; likewise `behavior`,
  `analyze`, `center`, and so on.
- **Itemize what is enumerable.** Parallel, countable things get a list: options
  being weighed, alternatives, ordered steps, distinct cases, the fields of a
  struct. Prose buries exactly what a reader is meant to compare, so a paragraph
  listing three tradeoffs in sequence is worse than three bullets. Argument still
  reads as argument: why one option wins, or how a conclusion follows, is
  reasoning rather than an inventory and stays in prose. The test is whether the
  items are parallel and countable, not whether the passage is long.
- **Bold and headers, sparingly.** Bold for genuine emphasis or defined terms, not
  decoration. Headers only where they aid navigation, not on every paragraph.
- **Package name is `Oblio`** (capital O) in prose and documentation, matching the
  C++ namespace.
- **First person plural: always `we`.** Address the reader as `we`, not `you` and
  not `I`. "We compute the forest first" rather than "you compute" or "I compute".
  Use `our` for the possessive. This keeps a single, collaborative voice across all
  documents.

## Terminology

**Scope note: this section is the one part of this file that also governs code.** Names and
comments must use the same words as the prose, or the two drift apart and a reader has to hold
a translation table. So the vocabulary below applies to prose, to code comments, and to
identifiers alike. Everything else in this file is prose-only.

**Numeric factorization, not numerical.** The phases are **symbolic** and **numeric**. Both
adjectives are legitimate English, and the distinction is real: *numeric* means "consisting of
numbers", *numerical* means "relating to the study or manipulation of numbers". The
factorization consists of numbers, as against the symbolic one, which consists of structure. It
is not a study of numbers. The field is near-unanimous here (Davis, SuiteSparse, MUMPS,
UMFPACK, PARDISO all say "numeric factorization"), so this is a term of art and not a matter of
taste.

Keep *numerical* for the mathematics, where it is correct: **numerical stability**, **numerical
pivoting**, **numerical rank**, **numerical analysis**. Those really are about the manipulation
of numbers.

The pair is not parallel, which is why only one half of it ever feels uncertain: *symbolic* has
no competing form (there is no "symbolical"), so the ambiguity is entirely on the numeric side.
Knowing that is half of not second-guessing it.

Short forms: `SymFactor` / `symfactor` and `NumFactor` / `numfactor`.

**Front is the noun, frontal the adjective.** A supernode has a *front* (its columns), so the
field is `frontSize` and the prose says "front indices". A function may use the adjective, hence
`gatherFrontalIndices`. The split is the ordinary one: the front is a thing, a frontal index is
a kind of index. See CODING_RULES for the naming half of this.

**Case carries the level: lowercase for a column, uppercase for a supernode.** In prose and in
pseudocode, `j` and `k` are single columns and `J` and `K` are supernodes. This is the ordinary
convention of numerical linear algebra, where a block of a partitioned matrix is `A_{IJ}` and the
uppercase subscript is an *index set*. A supernode is exactly that: a set of columns. So the
notation is doing its usual job, not being borrowed for a new one.

No collision with `A` and `L`, which are matrices: nobody will read `Idx(K)` as indexing a
matrix, and the context is never in doubt.

Code cannot use case this way (`J` and `K` are poor variable names), so **the doc's uppercase is
the code's doubling**: `K` in prose is `kk` in code, `J` is `jj`. Both say "the supernode", one
by case and one by repetition, and each is the natural device for its medium. See CODING_RULES
for the doubling rule.

**Node and supernode, column and row, index.** These name the same object seen three ways, and
all three are legitimate; the structure is symmetric, so an index is an index, and row `k`,
column `k` and vertex `k` are the same `k`. Prefer the vocabulary of the context: matrix words
(row, column, `colPtr`, `rowIdx`) when discussing `A` or `L` as matrices, graph words (node,
supernode, parent, child) when discussing the elimination forest, and the neutral **index** and
**snode** where the two meet, which is most of the symbolic factorization. `snode` is deliberately
neutral: a node is neither row nor column (the structure is symmetric), so a supernode of nodes
commits to neither.

**Clique, not element.** An eliminated pivot's fill-in structure is a **clique**. AMD and genmmd
both call it an *element*, and every ordering paper in that lineage follows them, so the word is
in every source we read and it leaks into our own writing constantly. It is still not our word.
Clique says what the object IS, a set of vertices now mutually adjacent, and it is the word the
rest of this tree already uses: `QuotientGraphFlat::clique`, `cliqueSize`, `cliqueDegree`, the clique
arena, and the counter "clique-member visits".

The one place *element* is correct is inside a quotation of code that uses it. When the prose is
walking through `mmdupd` or `AMD_2`, naming their `ehead` chain or their `Elen` array, their word
has to stand or the explanation stops matching the source. Make it visibly theirs, and say once
per document that element is their name for our clique. Everywhere else, including comments and
identifiers in our own files, it is clique.

*Element* in the ordinary English sense is untouched: "a d-element list", "an element of `A[u]`",
`std::min_element`. That is a different word and it collides with nothing.

**The Cholesky factor is `C`, the LDL factor is `L`.** Cholesky is `A = CC^H` (`CC^T` in real);
LDL is `A = LDL^H` (or `LDL^T`). The letter carries the distinction the two turn on: `C` holds its
own diagonal, whereas `L` is *unit* lower triangular with the diagonal pulled out into `D`. Writing
both as `L` erases that contrast; "for Cholesky the diagonal holds the factor's own, for LDL it
holds `D`" says itself once it reads as `C` against `L`. This is already the notation of the design
notes (`CC^H` throughout); the rule only makes it bind on comments too. The one exception is the
shared numeric-kernel blocks `L11`, `L21`: a single kernel computes them for both factorizations,
so they are generic and stay `L`.

**Default to `^H` for the conjugate transpose; use `^T` only where the factorization does not
conjugate.** On a real matrix `^H` is `^T`, since conjugation does nothing there, so `^H` is the
general form that specializes correctly and needs no "(`^T` in real)" parenthetical. Write `A = CC^H`,
`C21 = A21 C11^-H`, and so on. Cholesky is always Hermitian, so it is `^H` throughout with no
exception. LDL has two variants and the notation must track which: the Hermitian variant is
`A = LDL^H` and the complex-symmetric variant is `A = LDL^T`, which deliberately does not conjugate,
so `^H` there would be wrong. This mirrors the code, where `hermitian(factorization())` and
`maybeConjugate(...)` select between the two. So `^H` is the default and the real case is silent
under it; `^T` is reserved for the complex-symmetric path, where it is load-bearing.

**Static and dynamic mean pivoting on an algorithm, storage on a container.** The two words name
two different axes, and the noun they qualify says which one is meant, so nothing has to be spelled
out. On an algorithm object, static and dynamic are about *pivoting*:
`NumFactorEngine::factorStaticLeftLooking` runs the factorizations whose pivots are fixed by the
symbolic structure (Cholesky, static LDL), and `factorDynamicNonRootSupernode` runs the one that
chooses pivots while the arithmetic runs. On a data object they are about *storage*: `NumFactorStatic` holds
flat buffers that cannot grow, `NumFactorDynamic` holds one vector per supernode so that a front
can. A field describing a storage choice says so in its own name (`mUsesDynamicStorage`).

The axes are aligned but not identical, and the asymmetry is worth knowing: static pivoting runs in
either storage, while dynamic pivoting requires the dynamic one, because delaying a column grows a
front. `dynamicPivoting()` in `Types.h` is where that rule is stated once.

**Assembly. Not folding, and not extend-add.** The operation that adds a child's contribution into
its parent is **assembly**. Write *assembled into*, *assembles*, *the assembly step*, *assembly
order*. Every function that performs it is already named for it, `assembleFromA`, `assembleUpdateBlock`,
`assembleUpdateMatrix`, `assembleDelay`, so the prose and the identifiers say one word between them.

Two alternatives are excluded. *Folded in*, *folds*, *the fold* is an informal coinage with no
standing anywhere. *Extend-add* does have standing, being the usual term in the multifrontal
literature, but it names a mechanism (extend the index set, add the values) where we want the
operation, and it reads as jargon beside four functions that already say `assemble`. Where the
sparse-to-dense routing itself is the point, name that directly: the *scatter*, or `gblToLcl`.

*Fold* stays correct in its ordinary English sense of absorbing one thing into another where no
assembly is meant, as in Cholesky folding its diagonal into `C`, or a prepass being folded into a
driver to keep one shape. The rule is about the operation, not about the word.

**Spell the term out, and give the acronym once where the concept is defined.** The prose says
*default member initializer*, not NSDMI, every time; the acronym appears exactly once, in
parentheses at the rule in `CODING_RULES.md`. Two reasons, and the second is the one that gets
forgotten. An acronym is opaque to a reader who has not met it, and expanding it every time is
noise, so the full term as the standing form is simply easier to read. But a term that appears
nowhere in its short form is also **unfindable**: someone who knows the concept by its acronym
searches for it, gets nothing, and concludes the tree does not cover it. The single parenthetical
is what makes a search from either direction land on the rule.

This applies to jargon acronyms, the ones naming a language or library mechanism. Acronyms that
are the field's ordinary vocabulary, BLAS, CSC, AMD, MMD, are used bare and need no expansion:
the audience for a sparse direct solver has them already. The test is whether a competent reader
of this codebase would have to look it up.

Two existing gaps, left as they are for now: `OBLIO_NOTES_FROM_POLYGLOT.md` uses NRVO three times
and `DESIGN_DECISIONS.md` uses RAII once, neither expanded anywhere.

## A code comment says what the code does. Everything else goes in a document

**A comment states what the code does and what a reader must not break. Why it is that way, what
it was before, what alternatives were tried, and what any of it measured all belong in a document.**
The documents carry a date, a machine, a build and a method, and a later dated entry supersedes
them where a reader will see it. A comment carries none of that: it is a claim about one afternoon,
sitting in a file nobody re-reads to check whether it still holds, and the next reader cannot
falsify it because they have no idea what produced it. It is not a weak version of a document
entry; it is noise wearing a document's clothes.

**FOUR KINDS THAT DO NOT BELONG, and the fourth is the one that keeps coming back:**

- **Measurements.** Times, ratios, speedups, instruction counts, cache figures, compaction counts,
  allocation sizes, percentages of anything including how often a branch is taken, the size of the
  input a bug was found on, how long a run took. "Measured at n = ...", "reads 1 from 8 to 401 a
  side", "worth about 5 per cent".
- **History.** What the code was before, what commit changed it, what the previous arrangement
  cost. A version-control system already holds this and holds it accurately.
- **Alternatives tried and rejected.** The version that was built and reverted, the approach that
  turned out slower, the design that was considered and dropped.
- **Justification for the arrangement.** Why this file is laid out this way, why the class sits
  where it does, what the choice buys. This is the one that looks harmless because it contains no
  numbers, and it is the same defect: it is an argument about a decision, and arguments about
  decisions live in `notes/DESIGN_DECISIONS.md`.

**WHAT IS OFTEN HIDING BEHIND ALL FOUR IS AN INVARIANT, AND THAT DOES BELONG IN THE CODE.** The
useful content of "the collector fires once per ordering at the shipped headroom" is that the
mid-walk path is not exercised by an ordinary run, and a reader of that code needs to know it. The
useful content of "we made these header-only so the comparison is fair" is that every out-of-class
definition must be `inline` or the next one added will not link. Write the property, drop the
argument, and let the dated entry hold the reasoning:

```
// Every out-of-class definition below is `inline`, and so is any variable defined here: several
// translation units include this. See CODING_RULES.md.
```

not

```
// A driver has to be in the same translation unit as its quotient graph or a comparison against
// the vendored routine is not apples to apples, both vendored files being single units...
```

**TWO TESTS, and a comment has to pass both:**

- **Would it have to change if someone reran the benchmark on different hardware, or if the
  decision were revisited?** Then it is a measurement or an argument, and it goes in a document.
- **Does it tell a reader something they need in order not to break the code?** If not, it does
  not earn its place regardless of whether it is true.

**And the reference goes the other way.** A comment may point at `notes/DESIGN_DECISIONS.md` for
the reasoning, which is a free reference by the rule below, and must still read completely without
it.

## References, and why the rule differs by artifact

Only one case is restricted: **a reference from a coding artifact to another coding artifact.**
Everything else is free.

**In any `.md` file, anywhere in the tree, link freely.** No restriction at all. Cross-references
are how the record works and there is no reason to ration them. Markdown is a special category:
the set is small, it is under constant scrutiny, it is reorganized as a deliberate act in which
every reference can be brought along at once, and a stale link is caught the next time anyone
reads the paragraph. **Pointing from code INTO a `.md` file is equally free**, for the same
reason: the target is a document whose job is to be pointed at.

**Between coding artifacts, keep it tight.** Sources, headers and Makefiles get a stricter rule
for one reason: nobody re-reads a comment to check whether its pointer still resolves. A reference
in code is never audited, so it has to be one that cannot go quietly wrong.

- **Reference along tight connections, not loose ones.** A reference is safe when something already
  binds the two files, so it cannot rot in silence: a Makefile naming a source it compiles, a rule
  naming its own prerequisite, a header naming the unit that implements it. If that file moves, the
  build or the compile breaks immediately and loudly, and the breakage is what keeps the comment
  honest. A reference is unsafe when nothing binds them, and one Makefile in an experiment pointing
  at another Makefile in a benchmark is the clearest case: they share no dependency, either may be
  renamed, moved or deleted with everything still building, and the pointer dies without a sound.
- **A comment is complete where it stands.** Even along a tight connection, the reference adds
  DEPTH and never carries the explanation. Someone who never opens the other file must be able to
  tell what the line does and why it is there. `See X: without this the binary goes stale` is the
  shape to avoid, since it leaves a sentence with its subject in another file. The test is to
  delete the reference and read what remains: if it no longer explains, the explanation was in the
  wrong place. Never cite a line number.

## Dated entries (DESIGN_DECISIONS, CHANGELOG, the devlog)

**Read the date before writing it.** An entry is stamped with the date it is *written*, and that
sounds too obvious to need a rule until you notice how it goes wrong: not by getting the date
wrong, but by **never looking it up**. The date carries forward by habit, from the last entry, from
the last session, from whatever was true when the work began.

Two failures, both of which happened here:

- **A session that runs past midnight.** Work started on the 13th and continued into the 14th; five
  entries were stamped 07-13 and one of them was written the next morning. Nothing felt wrong at
  the time, because nothing had changed.
- **A gap between writing and committing.** An entry written on the 14th and dated from memory as
  the 13th, because the previous entries said 13.

**So: look the date up, per entry, at the moment of writing.** Not once at the start of a session,
not inherited from the entry above.

**And the record can be checked, which is better than trusting it.** The git history knows when an
entry actually landed:

```
git show <commit> -- notes/DESIGN_DECISIONS.md | grep "^+## "
```

lists the entry headings a commit introduced. If a heading's date and its commit's date disagree,
the heading is wrong. Worth running when a decision matters enough to be cited later, which most of
them are.

## Document-specific style

Conventions internal to a single document stay with that document, not here. For
example the 1-indexed math, the `j < k` index roles, the lower-triangle `A`
accesses, and the pseudocode comment style of `sparse_factorization.md` are that
document's own style. This file holds only what applies across all prose.

### README build targets

The README's Build section lists the whole-project targets a user reaches for (`make`, `make all`,
`make test`, `make tests`, `make examples`, `make clean`). Single-unit and inner-loop targets,
`make objs` and the per-file `%_cpp` / `example_%_cpp` rules, stay in the Makefile's own header
comment, which is the exhaustive list; the README is the curated subset. The test: a target that
acts on a whole category (all tests, all examples, the whole build) is user-facing and belongs in
both; one that compiles a single unit or just checks that the core builds is a contributor
convenience and stays in the Makefile. So updating the Makefile obliges a README edit only when a
whole-category verb changes, which is why adding `objs` did not.

`make help` is named on the `make` line rather than given a row of its own. It is not a
whole-category verb and does nothing the default goal has not already done, so a row would be the
only entry in that block adding a name rather than a capability. It is still worth naming, because
the README is read on GitHub by people who never run anything, and `help` is the name they would
try first.
