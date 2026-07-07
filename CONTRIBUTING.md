# Contributing

> **Stub.** Stubbed ahead of need — this is a solo, pre-release refactor with no
> outside contributors yet. Flesh this out when the repo goes public or gains a
> collaborator. At that point, make this file the **single source of truth** for
> build/test/style so it doesn't fork from the docs below.

For now, the working docs are canonical; this file only points at them:

- **Build & test** — see CLAUDE.md (`## Build`, `## Tooling`).
- **Coding rules** — invariants in CLAUDE.md; conventions in docs/CODING_RULES.md;
  mechanical style enforced by `.clang-format` / `.clang-tidy`.
- **Design rationale & history** — docs/DESIGN_DECISIONS.md, and archive/ for the
  development history.
- **Porting workflow** — one unit at a time, verified against oblio 0.9; status
  tracked in docs/PORTING_LEDGER.md.

## When this goes public (checklist for later)

- [ ] Move the actual build/test commands here as the canonical copy.
- [ ] Add PR / review expectations.
- [ ] Add a `LICENSE` and state it here (note: bundled AMD is BSD-3-clause).
- [ ] Point contributors at the ledger + 0.9-oracle discipline.
