#!/usr/bin/env python3
"""Generate a hooked copy of the vendored AMD that emits its RAW elimination order.

    python3 hook_amd.py ../../private/Amd.cpp amd_raw.cpp

WHY A GENERATOR AND NOT A CHECKED-IN COPY. `private/` holds the vendored routine as SuiteSparse
ships it, and its whole value is being exactly that. A copy carrying our edits would be a third
thing, drifting from the original with nothing to notice, and a test comparing against a stale
oracle is worse than no test. So the copy is generated from whatever `private/Amd.cpp` currently
says, gitignored, and removed by `make clean` -- the same arrangement as the int64 copies the width
study uses.

WHAT THE HOOK IS FOR. `AMD_2` does not emit the elimination order. `amd_order` returns `Perm`,
which `AMD_1` has already relabeled by `AMD_postorder`, and that postorder is a heuristic tidy of
an assembly tree that Amd.cpp's own header says is "not guaranteed to be the precise supernodal
elimination tree". Oblio replaces it with Liu's rule on the exact supernodal tree, so the two
output vectors differ by construction. What we need is the order AMD_2 would emit if it stopped at
the end of its main loop: the pivot sequence together with the member order inside each
supervariable, which is the whole algorithm.

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

DECLS = """#include <vector>

/* ---- added by experiments/ordering/hook_amd.py; not part of the vendored source ---- */
static std::vector<std::vector<int> > PB_members;   /* members[i], what supervariable i stands for */
static std::vector<int>               PB_raw;       /* the raw elimination order */
static void PB_reset (int n)
{
    PB_members.assign (n, std::vector<int> ()) ;
    for (int i = 0 ; i < n ; i++) PB_members [i].push_back (i) ;
    PB_raw.clear () ;
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


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: hook_amd.py <vendored Amd.cpp> <output>")
    src, dst = sys.argv[1], sys.argv[2]
    s = open(src).read()

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
                  "    { extern void pbRawOrder (const int*, int) ;\n"
                  "      pbRawOrder (PB_raw.empty () ? 0 : &PB_raw[0], (int) PB_raw.size ()) ; }\n"
                  "    return pb_status;")

    open(dst, "w").write(DECLS + s)
    print(f"hook_amd.py: wrote {dst} from {src}")


if __name__ == "__main__":
    main()
