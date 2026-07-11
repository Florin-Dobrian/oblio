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

## Document-specific style

Conventions internal to a single document stay with that document, not here. For
example the 1-indexed math, the `j < k` index roles, the lower-triangle `A`
accesses, and the pseudocode comment style of `sparse_factorization.md` are that
document's own style. This file holds only what applies across all prose.
