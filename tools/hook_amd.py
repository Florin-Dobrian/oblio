#!/usr/bin/env python3
"""Generate a hooked copy of the vendored AMD, in one of two modes.

    python3 tools/hook_amd.py private/AmdVendored.cpp <out>.cpp                 the RAW ORDER
    python3 tools/hook_amd.py private/AmdVendored.cpp <out>.cpp --mode=timed    the PHASE TIMES

Both are additive: neither changes what the routine computes, and each is checked against the
unhooked routine by its consumer before any number is read off it.

THE TWO MODES ANSWER DIFFERENT QUESTIONS, and the second exists because the first cannot be timed.

  raw     What order does AMD_2 eliminate in? The copy accumulates supervariable membership and
          emits the pivot sequence with member order, which is the acceptance test for our amd
          alignment. It carries a vector per vertex and a push per member, so it is an ORACLE and
          not a stopwatch: timing it would measure the bookkeeping.

  timed   How much of `amd_order` is work we also do? The copy takes six timestamps at phase
          boundaries and nothing else, so it is a stopwatch and not an oracle: it reports no order
          and reproduces the shipped routine's output exactly.

WHY THE SECOND MODE IS NEEDED AT ALL. Every ratio this tree quotes for an ordering is against
`amd_order` whole, and two of its phases have no counterpart on our side. `AMD_aat` forms the
pattern of A+A' because AMD takes a one-sided matrix; ours arrives full-symmetric with the diagonal
present, which is what a `SparseMatrix` holds. `AMD_postorder` relabels the assembly tree, which
`ElmForestEngine` redoes later with real front and update sizes on the exact supernodal tree. So a
figure like `AmdFlat 1.82x` is measured against a denominator carrying work we deliberately do
not do, and `benchmarks/ordering/README.md` has recorded one estimate of that, 8 percent at 140 a side,
without ever measuring it per size. This is the instrument that measures it per size.

WHY THE PHASES ARE TIMED RATHER THAN THE POSTORDER REMOVED, which is the first thing to reach for
and does not work. `AMD_postorder` runs INSIDE `AMD_2`, and the block after it builds the output
permutation out of the ranks it writes into `W`. Delete the call and `W` is unwritten, the returned
permutation is meaningless, and both of the checks that make a generated oracle trustworthy go with
it: the permutation compared entry for entry against the unhooked routine, and `Info[AMD_LNZ]`
compared against the fill this benchmark already records. A timestamp costs one call at a loop
boundary and keeps both.

THE PHASES, and what each is:

  0  valid    AMD_valid, amd_preprocess where the input is jumbled, and the Len and Pinv vectors.
  1  aat      AMD_aat, forming the pattern of A+A'.        WE DO NOT DO THIS.
  2  build    the S workspace, and AMD_1's construction of Iw and Pe from that pattern.
  3  core     AMD_2 from entry to the end of its main loop.
  4  post     the assembly-tree path compression, AMD_postorder, and the output permutation.
                                                            WE DO NOT DO THIS.

`build + core` is the comparable region: it is the vendored routine turning a caller's pattern into
its working structure and then ordering it, which is what `orderAmdFlat` does from `QuotientGraph`'s
constructor to its last pivot.

WHERE THIS LIVES, AND WHY NOT IN THE STUDY THAT FIRST NEEDED IT. It was written for
`experiments/ordering` and moved here on 2026-08-09 when `benchmarks/ordering` wanted the same
oracle. It could not simply be reached across: `benchmarks/README.md` states that the studies under
`experiments/` are frozen and that nothing in benchmarks may depend on their contents staying
current, so a benchmark rule invoking a script inside one would have broken the rule that keeps the
two kinds of folder apart. Copying it into both was the other option and is worse: two copies of a
generator whose whole job is to fail loudly when the vendored source moves, with nothing comparing
them. Promoting it makes the dependency a TOOL that both consumers use, which is a direction that
does not cross the boundary at all.

WHY A GENERATOR AND NOT A CHECKED-IN COPY. `private/` holds the vendored routine as SuiteSparse
ships it, and its whole value is being exactly that. A copy carrying our edits would be a third
thing, drifting from the original with nothing to notice, and a test comparing against a stale
oracle is worse than no test. So the copy is generated from whatever `private/AmdVendored.cpp`
currently says, gitignored, and removed by `make clean` -- the same arrangement as the int64
copies the width
study uses.

WHAT THE HOOK IS FOR. `AMD_2` does not emit the elimination order. `amd_order` returns `Perm`,
which `AMD_1` has already relabeled by `AMD_postorder`, and that postorder is a heuristic tidy of
an assembly tree that Amd.cpp's own header says is "not guaranteed to be the precise supernodal
elimination tree". Oblio replaces it with Liu's rule on the exact supernodal tree, so the two
output vectors differ by construction. What we need is the order AMD_2 would emit if it stopped at
the end of its main loop: the pivot sequence together with the member order inside each
supervariable, which is the whole algorithm.

AND IT CANNOT BE RECOVERED WITHOUT THIS, which is worth recording because the cheaper route is the
first thing anyone looks for. `AMD_1` builds its permutation from `Pe`, `Nv` and the postorder
ranks alone, and the elimination order is not among what `AMD_2` leaves behind: the assembly tree
says which element absorbed which, and the order pivots were taken in is a linear extension of that
tree rather than a property of it. Nor is there a knob: `Control` carries exactly two, `AMD_DENSE`
and `AMD_AGGRESSIVE`, and the postorder is unconditional.

WHAT IT CHANGES. Nothing about behaviour. Two containers and three statements:

  PB_members[i]   what supervariable i currently stands for, seeded as {i} and grown at the hash
                  merge, beside `Nv [i] += Nv [j]`
  PB_raw          the raw order, taking a pivot's members once per elimination at FINALIZE THE NEW
                  ELEMENT, which is before AMD_postorder ever runs, AND taking empty variables
                  where the initialization numbers them: a vertex with no off-diagonal entry never
                  forms an element and so never reaches the finalize marker

and the entry point renamed so both the original and the hooked copy can be linked at once.

EVERY ANCHOR IS ASSERTED. If the vendored source moves under us the generation FAILS LOUDLY rather
than producing a copy with the hook in the wrong place, which would compare against a plausible
wrong answer. That is the failure mode this whole arrangement exists to avoid.
"""

import sys

DECLS_RAW = """#include <vector>

/* ---- added by tools/hook_amd.py --mode=raw; not part of the vendored source ---- */
static std::vector<std::vector<int> > PB_members;   /* members[i], what supervariable i stands for */
static std::vector<int>               PB_raw;       /* the raw elimination order */
static std::vector<int>               PB_dense;     /* the rows the dense rule set aside */
static void PB_reset (int n)
{
    PB_members.assign (n, std::vector<int> ()) ;
    for (int i = 0 ; i < n ; i++) PB_members [i].push_back (i) ;
    PB_raw.clear () ;
    PB_dense.clear () ;
}
/* ----------------------------------------------------------------------------------- */

"""

# The timed mode's declarations. Six timestamps and nothing else: no container, no per-vertex
# state, no allocation, so the copy computes exactly what the shipped routine computes and its
# consumer can check that it does.
#
# A CLOCK READ IS A COMPILER BARRIER, which is the honest caveat and the reason the consumer prints
# the hooked total beside the unhooked one. Each of the five sits at a phase boundary, between two
# loops rather than inside one, so there is nothing across it for the optimizer to have been moving;
# but that is an argument, and the two totals side by side are the measurement. Five reads against
# an ordering of tens of microseconds is on the order of a hundred nanoseconds either way.
#
# `PB_mark(k)` closes phase k and opens the next, so the phases partition the run and the caller
# needs one stamp per boundary rather than a pair per phase. `PB_reset` zeroes them, so a call that
# returns early, on n = 0 or on invalid input, reports zeros rather than a stale reading from the
# previous ordering.
DECLS_TIMED = """#include <chrono>

/* ---- added by tools/hook_amd.py --mode=timed; not part of the vendored source ---- */
static double PB_ms [5] ;      /* 0 valid, 1 aat, 2 build, 3 core, 4 post; milliseconds */
static std::chrono::steady_clock::time_point PB_at ;

static void PB_mark (int phase)
{
    const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now () ;
    PB_ms [phase] = std::chrono::duration<double, std::milli> (now - PB_at).count () ;
    PB_at = now ;
}

static void PB_reset (void)
{
    for (int k = 0 ; k < 5 ; k++) PB_ms [k] = 0.0 ;
    PB_at = std::chrono::steady_clock::now () ;
}

extern "C" void amd_phase_ms (double out [5])
{
    for (int k = 0 ; k < 5 ; k++) out [k] = PB_ms [k] ;
}
/* ----------------------------------------------------------------------------------- */

"""


def need(text, anchor, what):
    if text.count(anchor) != 1:
        sys.exit(f"hook_amd.py: {what} -- found {text.count(anchor)} times, expected 1.\n"
                 f"  The vendored source has moved. Fix the anchor here rather than\n"
                 f"  loosening it: a hook in the wrong place compares against a plausible\n"
                 f"  wrong answer, which is worse than a build failure.\n"
                 f"  anchor: {anchor!r}")


def transform_raw(s):
    # 1. the hash merge, where a supervariable absorbs another
    #    j is folded into i, so i's list grows and j's is EMPTIED. Clearing matters: without it a
    #    member reached through a chain of merges is copied more than once.
    merge = "Nv [i] += Nv [j] ;"
    need(s, merge, "the hash merge moved")
    s = s.replace(merge, merge + "\n\t\t\t    PB_members [i].insert (PB_members [i].end (),\n"
                                 "\t\t\t\tPB_members [j].begin (), PB_members [j].end ()) ;\n"
                                 "\t\t\t    PB_members [j].clear () ;")

    # 2. MASS ELIMINATION, the second merge path and the easy one to miss. Scan 2 folds a variable
    #    whose degree has fallen to the pivot's straight into `me`, with `Nv [i] = 0`, and it never
    #    goes through the hash merge above. Missing it left the raw order short by three entries on
    #    most grids and six on one, which is exactly the irregular deficit that says a whole path
    #    is unhandled rather than an off-by-one.
    mass = "\t\tNv [i] = 0 ;\n\t\tElen [i] = EMPTY ;"
    need(s, mass, "the mass-elimination branch moved")
    s = s.replace(mass, "\t\tPB_members [me].insert (PB_members [me].end (),\n"
                        "\t\t    PB_members [i].begin (), PB_members [i].end ()) ;\n"
                        "\t\tPB_members [i].clear () ;\n" + mass)

    # 3. EMPTY VARIABLES, numbered in the initialization before the main loop ever runs. A vertex
    #    with no off-diagonal entry never forms an element and so never reaches the finalize
    #    marker below; it is simply numbered where it stands. Missing these is what made a first
    #    version of this hook come up short by two to six entries on every grid.
    empty = "\t    Elen [i] = FLIP (1) ;\n\t    nel++ ;"
    need(s, empty, "the empty-variable branch moved")
    s = s.replace(empty, empty + "\n\t    PB_raw.push_back (i) ;")

    # 3b. DENSE VARIABLES, which AMD_2 does not eliminate at all: `Nv [i] = 0`, `Pe [i] = EMPTY`,
    #    and the comment there says they take no part in the postorder. So they never reach the
    #    finalize marker and the raw order came up short by `ndense` on any matrix with a hub,
    #    which is most social and power-law graphs in the collection.
    #
    #    THEY ARE COLLECTED HERE AND APPENDED AT THE END, not pushed where they are found, because
    #    that is where AMD_2's own output assembly puts them: `Next [i] = nel++` under "This is a
    #    dense unordered variable, with no parent. Place it last in the output order", over i
    #    ascending. This branch also runs over i ascending, so collecting here gives that order.
    dense = "\t    ndense++ ;\n\t    Nv [i] = 0 ;\t\t/* do not postorder this node */"
    need(s, dense, "the dense-variable branch moved")
    s = s.replace(dense, dense + "\n\t    PB_dense.push_back (i) ;")

    # 4. once per pivot, at the finalize marker, before AMD_postorder runs
    fin = "/* FINALIZE THE NEW ELEMENT */"
    need(s, fin, "the finalize marker moved")
    cut = s.index("\n", s.index("*/", s.index(fin) + len(fin))) + 1
    s = s[:cut] + "\n\tfor (size_t q = 0 ; q < PB_members [me].size () ; q++)\n" \
                  "\t    PB_raw.push_back (PB_members [me][q]) ;\n" + s[cut:]

    # 5. the entry point, renamed so both copies can be linked, and wrapped to hand the order back
    entry = 'extern "C" int amd_order('
    need(s, entry, "the entry point moved")
    s = s.replace(entry, 'extern "C" int amd_order_raw(')

    call = "    return AMD_order(n, Ap, Ai, P, Control, Info);"
    need(s, call, "the AMD_order call inside the entry point moved")
    s = s.replace(call,
                  "    PB_reset (n) ;\n"
                  "    const int pb_status = AMD_order(n, Ap, Ai, P, Control, Info);\n"
                  "    PB_raw.insert (PB_raw.end (), PB_dense.begin (), PB_dense.end ()) ;\n"
                  "    { extern void pbRawOrder (const int*, int) ;\n"
                  "      pbRawOrder (PB_raw.empty () ? 0 : &PB_raw[0], (int) PB_raw.size ()) ; }\n"
                  "    return pb_status;")

    return DECLS_RAW + s


def transform_timed(s):
    # 1. AMD_aat, bracketed. The mark before it closes `valid`, which is everything AMD_order does
    #    up to here: the validity check, amd_preprocess where the input is jumbled, and the Len and
    #    Pinv vectors. The mark after it closes `aat` itself.
    aat = "    size_t nzaat = AMD_aat(n, Cp, Ci, Len.data(), P, Info);"
    need(s, aat, "the AMD_aat call moved")
    s = s.replace(aat, "    PB_mark (0) ;\n" + aat + "\n    PB_mark (1) ;")

    # 2. the AMD_2 call, whose mark closes `build`: the S workspace allocated in AMD_order, and
    #    AMD_1's construction of Iw and Pe out of the A+A' pattern. That is the phase our own
    #    QuotientGraph constructor corresponds to.
    two = "    AMD_2 (n, Pe, Iw, Len, iwlen, pfree,"
    need(s, two, "the AMD_2 call moved")
    s = s.replace(two, "    PB_mark (2) ;\n" + two)

    # 3. the end of AMD_2's main loop, which is where the path compression for the assembly tree
    #    begins. Everything from here is for AMD_postorder and the output permutation, and Oblio
    #    does none of it: ElmForestEngine postorders the exact supernodal tree later, with real
    #    front and update sizes.
    #
    #    The mark goes AFTER the banner rather than before it, so that the phase boundary sits with
    #    the code it opens rather than between two comment lines. Same idiom as the finalize marker
    #    in the raw transform above: find the closing `*/` of the banner line below the anchor, then
    #    the newline after it.
    compress = "/* compress the paths of the variables */"
    need(s, compress, "the path-compression marker moved")
    cut = s.index("\n", s.index("*/", s.index(compress) + len(compress))) + 1
    s = s[:cut] + "\n    PB_mark (3) ;\n" + s[cut:]

    # 4. the entry point, renamed so the shipped routine, the raw copy and this one can all be
    #    linked at once, and wrapped so the phases are zeroed before the call and `post` is closed
    #    after it. Nothing between the last mark and here but AMD_order's own two-line tail.
    entry = 'extern "C" int amd_order('
    need(s, entry, "the entry point moved")
    s = s.replace(entry, 'extern "C" int amd_order_timed(')

    call = "    return AMD_order(n, Ap, Ai, P, Control, Info);"
    need(s, call, "the AMD_order call inside the entry point moved")
    s = s.replace(call,
                  "    PB_reset () ;\n"
                  "    const int pb_status = AMD_order(n, Ap, Ai, P, Control, Info);\n"
                  "    PB_mark (4) ;\n"
                  "    return pb_status;")

    return DECLS_TIMED + s


MODES = {"raw": transform_raw, "timed": transform_timed}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    mode = "raw"
    for f in flags:
        if f.startswith("--mode="):
            mode = f[len("--mode="):]
        else:
            sys.exit(f"hook_amd.py: unknown option {f!r}")
    if len(args) != 2 or mode not in MODES:
        sys.exit("usage: hook_amd.py <vendored Amd.cpp> <output> [--mode=raw|timed]")

    src, dst = args
    s = MODES[mode](open(src).read())
    open(dst, "w").write(s)
    print(f"hook_amd.py: wrote {dst} from {src}, mode {mode}")


if __name__ == "__main__":
    main()
