# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning
intent is [Semantic Versioning](https://semver.org/).

> **Stub.** No release of the modern tree has been cut yet, so there are no
> versioned entries. Start real entries at the first tagged version.
>
> **Open — version identity of the modern tree.** Undecided whether this continues
> the old line (e.g. oblio 11) or restarts (e.g. 1.0). The first release entry
> can't be written until this is settled; it ties to the working-tree question in
> DESIGN_DECISIONS.md. Record the choice there and seed the first entry here.

## [Unreleased]

Work in progress on the modern refactor. Session-by-session detail lives in
archive/oblio-new-devlog.md; per-unit port status in PORTING_LEDGER.md. Summarize
here only when cutting a release.

### Lineage (historical, for context)

- **oblio 0.9** — original, complete, correct sparse direct solver (circa
  2003–2005). Reference and oracle.
- **oblio 10.12** — partial, unfinished refactor of 0.9.
- **modern refactor** (this tree) — port of 0.9 into modern C++17: `Array` →
  `std::vector`, modernized idioms, one function at a time.
