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
- **Minimal formatting.** Prefer flowing prose to lists. Use bold sparingly, for
  genuine emphasis or defined terms, not decoration. Add headers only where they
  aid navigation, not on every paragraph.
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

**The Cholesky factor is `C`, the LDL factor is `L`.** Cholesky is `A = CC^H` (`CC^T` in real);
LDL is `A = LDL^H` (or `LDL^T`). The letter carries the distinction the two turn on: `C` holds its
own diagonal, whereas `L` is *unit* lower triangular with the diagonal pulled out into `D`. Writing
both as `L` erases that contrast; "for Cholesky the diagonal holds the factor's own, for LDL it
holds `D`" says itself once it reads as `C` against `L`. This is already the notation of the design
notes (`CC^H` throughout); the rule only makes it bind on comments too. The one exception is the
shared numeric-kernel blocks `L11`, `L21`: a single kernel computes them for both factorizations,
so they are generic and stay `L`.

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
git show <commit> -- docs/DESIGN_DECISIONS.md | grep "^+## "
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

The README's Build section lists the whole-project targets a user reaches for (`make`,
`make test`, `make tests`, `make examples`, `make clean`). Single-unit and inner-loop targets,
`make objs` and the per-file `%_cpp` / `example_%_cpp` rules, stay in the Makefile's own header
comment, which is the exhaustive list; the README is the curated subset. The test: a target that
acts on a whole category (all tests, all examples, the whole build) is user-facing and belongs in
both; one that compiles a single unit or just checks that the core builds is a contributor
convenience and stays in the Makefile. So updating the Makefile obliges a README edit only when a
whole-category verb changes, which is why adding `objs` did not.
