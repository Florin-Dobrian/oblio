# ordering experiment - the minimum degree family, built up one mechanism at a time.
# Reference/teaching only; not part of the main Oblio build.
#
#:   make          build every C++ prototype
#:   make test     build them, run each, check the C++ agrees with its Python twin, and check
#:                 each ported layer agrees with the production driver extracted from it
#:   make production   the production drivers on these graphs, linked from ../../src
#:   make width        what the integer width costs, measured on the vendored AMD by retyping it
#:                     to int64_t and timing the two against each other. Needs private/.
#:   make aligned      both acceptance tests below, which is the one word for "is either
#:                     ordering still what the vendored routine computes". Needs private/.
#:   make amdorder     production Amd3 against the vendored AMD's RAW elimination order, which is
#:                     the acceptance test for the amd alignment. Four shapes: the seven
#:                     examples, 2D grids, 3D grids, random patterns. Needs private/.
#:   make mmdorder     production Mmd3 against genmmd's elimination order, the same four shapes.
#:                     The mmd counterpart, and it needs no hook: genmmd emits that order and
#:                     does no postorder. Needs private/.
#:   make vendored     the vendored MMD and AMD on the same graphs, for reference. This is the
#:                 only target here that needs private/, and it reports rather than fails without
#:                 it; OBLIO_PUBLIC=1 makes it behave as though private/ were absent
#:   make clean
#:   make help     print this list
#:
#: Prefix any of these with OBLIO_PUBLIC=1 to build as everyone else does, without the vendored
#: orderings; see docs/DESIGN_DECISIONS.md. The same word works in the repo root and in both
#: benchmarks. Here it reaches only `make vendored`, that being the only target that needs them.
#
# Each layer exists twice, as <name>.py and <name>.cpp, and the two are written to
# print the same trace. That is the point of the pairing: the Python is the one to
# read and the C++ is the one to trust, and `make test` checks they have not drifted.
# Executables carry the _cpp language suffix (coding convention) and are gitignored.
#
# The layers, each adding exactly one mechanism to the one before:
#
#   md1   naive minimum degree, materializing the fill              5.1, 5.2
#   md2   the quotient graph: cliques instead of fill               5.3, 5.4
#   md3   supervariables and mass elimination                       5.5, 5.6
#   md4   maintained degrees, refreshed only where they changed     5.7, 5.8
#   md5   degree buckets, so the minimum is walked to not scanned   5.9, 5.10
#   mmd1  multiple elimination: a batch of pivots per refresh       5.11, 5.12
#   mmd2  the rest of genmmd, one pass at a time                    5.11, 5.12
#   amd   approximate degree: a bound instead of a set union        5.13
#
# md1 through md5 all return the SAME ordering; mmd1 and amd do not, deliberately.
# See archive/sparse_factorization.md section 5 for the sections cited above.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra

# Only the production check links Oblio, whose sources call BLAS. Accelerate on macOS, the system
# BLAS elsewhere.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    BLAS_LIBS ?= -framework Accelerate
else
    BLAS_LIBS ?= -lblas -llapack -lm
endif

# The interpreter that runs the Python twins. Bare python3 by default, for two
# reasons. It is a PATH name rather than a particular tool, so the prerequisite
# is only that some Python 3 is present, where `uv run` needs uv installed
# first and then resolves a Python of its own. And it is the same binary PyCharm
# is pointed at, so the editor and `make test` cannot come apart: one file,
# moved by one upgrade.
#
# `uv run` is a fine override, and it resolves elsewhere rather than badly. uv
# prefers its own managed interpreters under ~/.local/share/uv/python, which are
# on a separate upgrade cycle from Homebrew's. Measured on alpamayo and
# quitaraju, 2026-07-31: uv's 3.14.3 against Homebrew's 3.14.6, and after
# `uv python upgrade` both reading 3.14.6 while still being different builds,
# Jul 23 / Clang 22.1.3 against Jun 10 / Clang 21.0.0. Nothing in these layers
# is sensitive to that, being stdlib only, with no `hash()` reaching any trace
# and amd2's key arithmetic modulo n.
#
#   make test PYTHON="uv run"
PYTHON   ?= python3

# The main ladder, then the three md2 alternatives. Those four together are the square described
# in the README section "Zooming in on md2": exact against bounded on one axis, recomputed against
# maintained on the other. They are not rungs and nothing is extracted from them, so they never
# join PORTED below; their oracle is each other and their own twins.
LAYERS = md1 md2 md3 md4 md5 mmd1 mmd2 mmd3 amd1 amd2 amd3 amd4 mdm2 mda2 mdam2
BINS   = $(addsuffix _cpp,$(LAYERS))

# The layers that have been pulled into the main tree, and so have a production driver to be
# checked against. Grows as the extraction proceeds; empty would simply skip that check.
# Layers checked by permutation, and layers checked by fill. PORTED_FILL is empty and the weaker
# check exists only in case a future layer needs it: amd2 was in it until the postorder moved to
# what is now amd4, and every ported layer now matches production entry for entry.
#
# amd4 is not here at all and is not expected to be. Its four items are input conditioning, the
# postorder and a control interface, none of which Oblio has or wants, so there is no production
# counterpart to compare it against. Its oracle is its own C++ twin. That is also why the
# prototype's matrix1 has nothing to compare against on this side: it is the one example given as
# a matrix rather than a graph, unsymmetric with duplicates and an unsorted column, and it exists
# to exercise amd4Aat and amd4Preprocess.
#
# The digit is not the chain here: amd4 forks sideways from amd2 and is temporary, renamed out of
# the way so that amd3 can be the layer aligned to the vendored routine, which is what mmd3
# already is on the other branch. Its own header says so and says when it goes.
PORTED      = mmd1 mmd2 mmd3 amd1 amd2 amd3
PORTED_FILL =


# The grid sides each ported layer is also compared on, beyond the seven examples. The examples
# are at most twelve vertices, which is too small to fire a mechanism that needs real structure:
# two defects in amd2 left every one of them byte for byte identical while the ordering was wrong
# on any grid of 10 a side or more. So the examples check that the mechanisms are right in
# miniature and these check that they are still right at a size where they all fire.
#
# Two sides rather than one, and small ones. 10 is the smallest that caught the amd2 pair, so it
# is the regression case; 20 is four times the vertices, enough that mass elimination, aggressive
# absorption and the hash all fire repeatedly. Anything larger belongs in benchmarks/ordering,
# which measures rather than compares.
#
# What it costs, on the sandbox with everything already built: 0.46 s for the example checks
# alone, 0.90 s with the prototype-against-production grids added, 7.8 s with the twin grids as
# well. Nearly all of the last step is Python at side 20, and nearly all of that is the amd
# layers' exact-degree instrumentation, which computes per refreshed vertex the very union the
# bound exists to avoid. At side 10 alone the whole run is 1.2 s. So if the wait ever grates, the
# edit is to give the twin loop below its own sides rather than to drop the check: the two loops
# share this variable today only because they happened to want the same numbers.
GRID_SIDES  = 10 20

# The layers that have a grid mode, in both twins, and so can be compared on one. Every layer has
# it as of 2026-08-07, so this is LAYERS and is kept as its own variable only so the twin grid loop
# can be narrowed without touching the example loop. It used to be mmd1 upward, on the reasoning
# that a grid is unreadable as a trace and pointless as an example so the md ladder had no use for
# one. The first half is still true and no longer costs anything, since a run above SHOW_THRESHOLD
# prints no trace at all; the second half was answering the wrong question, because the grid is not
# there to be read, it is there to be diffed, and the md layers had no check at any size the seven
# examples cannot reach. amd4 is here although it is not in PORTED: it has no production
# counterpart to be checked against, which makes its twin the only oracle it has, and that is a
# reason to exercise it harder rather than less.
GRID_LAYERS = $(LAYERS)

# The production sources those drivers need, linked directly rather than copied: the point of this
# check is to notice when production and prototype come apart, and a copy would hide exactly that.
# The fill check runs Oblio's symbolic factorization, so the whole library is linked.
#
# No vendored code here. production.cpp calls Mmd1, Mmd2, Amd1 and Amd2 directly and never
# OrderEngine, so nothing in this link reaches Ordering::MMD or Ordering::AMD, and the two vendored
# files under ../../private are neither needed nor compiled. They used to be, but only because
# OrderEngine.cpp once declared their entry points unconditionally and the linker therefore
# demanded them.
PROD_SRCS   = $(wildcard ../../src/*.cpp)

# The headers too. Without them a change to include/oblio/QuotientGraph.h alone leaves
# production_cpp stale, and `make test` then compares this session's prototype against LAST
# session's production and reports agreement. See benchmarks/ordering/Makefile.
PROD_HDRS   = $(wildcard ../../include/oblio/*.h)
#
# VENDOR_SRCS is spelled the same here, in ../../Makefile and in both benchmark Makefiles, so
# `make <target> VENDOR_SRCS=` builds the published configuration from any of them.
# OBLIO_PUBLIC=1 builds as everyone else does, ignoring private/ even when it is there. Same word
# in every Makefile here that links the library. In this directory it reaches only the `vendored`
# target, that being the only thing here that needs those two files: `production` links the library
# but calls Mmd1, Mmd2, Amd1 and Amd2 directly and never OrderEngine, so it compiles and links
# identically either way.
ifdef OBLIO_PUBLIC
  VENDOR_SRCS =
else
  VENDOR_SRCS = $(wildcard ../../private/Amd.cpp ../../private/Mmd.cpp)
endif

.PHONY: all test production vendored width amdorder mmdorder aligned clean help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: $(BINS)

%_cpp: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Run each C++ prototype, then diff its trace against the Python twin's.
#
# The two print the same CONTENT but not the same punctuation: Python renders a
# dict as {0: [1,3], 1: [0,2]} and the C++ renders 0: {1,3}, 1: {0,2}, so the
# braces land in different places. Rather than churn fourteen files to make the
# spellings identical, the comparison strips spaces, brackets, braces, quotes and commas
# and requires the remainder to match. Every number, name and label still has to
# agree, which is what the check is for.
#
# Temp files rather than process substitution, so this runs under /bin/sh.
test: all production_cpp
	@rc=0; \
	for L in $(LAYERS); do \
	  printf '%-6s ' "$$L"; \
	  $(PYTHON) $$L.py | tr -d " []{},'" > .$$L.py.out; \
	  ./$${L}_cpp   | tr -d " []{},'" > .$$L.cpp.out; \
	  if diff .$$L.py.out .$$L.cpp.out > /dev/null; then \
	    echo "py and cpp agree"; \
	  else \
	    echo "DIFFER"; rc=1; \
	  fi; \
	  rm -f .$$L.py.out .$$L.cpp.out; \
	done; \
	for L in $(GRID_LAYERS); do \
	  for S in $(GRID_SIDES); do \
	    printf '%-6s ' "$$L"; \
	    $(PYTHON) $$L.py grid $$S | tr -d " []{},'" > .$$L.py.out; \
	    ./$${L}_cpp   grid $$S | tr -d " []{},'" > .$$L.cpp.out; \
	    if diff .$$L.py.out .$$L.cpp.out > /dev/null; then \
	      echo "py and cpp agree on grid $$S"; \
	    else \
	      echo "DIFFER on grid $$S"; rc=1; \
	    fi; \
	    rm -f .$$L.py.out .$$L.cpp.out; \
	  done; \
	done; \
	for L in $(PORTED); do \
	  printf '%-6s ' "$$L"; \
	  ./$${L}_cpp     | grep -E '^===|^order:' > .$$L.proto.out; \
	  ./production_cpp $$L | grep -E '^===|^order:' > .$$L.prod.out; \
	  if diff .$$L.proto.out .$$L.prod.out > /dev/null; then \
	    echo "prototype and production agree"; \
	  else \
	    echo "DIFFER"; rc=1; \
	  fi; \
	  rm -f .$$L.proto.out .$$L.prod.out; \
	done; \
	for L in $(PORTED); do \
	  for S in $(GRID_SIDES); do \
	    printf '%-6s ' "$$L"; \
	    ./$${L}_cpp     grid $$S | grep -E '^order:' > .$$L.proto.out; \
	    ./production_cpp $$L grid $$S | grep -E '^order:' > .$$L.prod.out; \
	    if diff .$$L.proto.out .$$L.prod.out > /dev/null; then \
	      echo "prototype and production agree on grid $$S"; \
	    else \
	      echo "DIFFER on grid $$S"; rc=1; \
	    fi; \
	    rm -f .$$L.proto.out .$$L.prod.out; \
	  done; \
	done; \
	for L in $(PORTED_FILL); do \
	  printf '%-6s ' "$$L"; \
	  ./$${L}_cpp | sed '/=== matrix1 ===/,$$d' | grep -E '^===|^nnz\(L\)' > .$$L.proto.out; \
	  ./production_cpp $$L | grep -E '^===|^nnz\(L\)' > .$$L.prod.out; \
	  if diff .$$L.proto.out .$$L.prod.out > /dev/null; then \
	    echo "prototype and production agree on fill"; \
	  else \
	    echo "DIFFER"; rc=1; \
	  fi; \
	  rm -f .$$L.proto.out .$$L.prod.out; \
	done; \
	exit $$rc

# Oblio's own ordering drivers, run on the same graphs as the prototypes. The mirror image of
# vendored_cpp: that one runs what we are reading against, this one runs what we extracted.
production_cpp: production.cpp graphs.h $(PROD_SRCS) $(PROD_HDRS)
	$(CXX) $(CXXFLAGS) -I../../include -DOBLIO_BLAS_UNDERSCORE production.cpp \
	  $(PROD_SRCS) $(BLAS_LIBS) -o $@

production: production_cpp

# The vendored routines, built separately: they are not layers and have no Python twin. They are
# compiled straight from private/ at the repo root rather than from a copy here, which is both one
# file less to keep in step and the only place they live now: private/ is not tracked and not
# published. Warnings are suppressed for them: they are not ours to clean up.
#
# The target reports rather than fails when the directory is absent, so a clone of the public tree
# can still build every layer.
ifeq ($(VENDOR_SRCS),)
vendored:
	@echo "vendored: skipped, private/ is not present"
width:
	@echo "width: skipped, private/ is not present"
# The two alignment checks say so too. They used to fall through with no rule at all, so make
# answered "Nothing to be done for 'amdorder'" and exited zero, which is a check silently
# declining to run: the same failure mode as grid mode's dead whitelist and order_profile's
# missing else, both recorded in the README. A target that cannot run must say why.
amdorder:
	@echo "amdorder: skipped, private/ is not present"
mmdorder:
	@echo "mmdorder: skipped, private/ is not present"
aligned:
	@echo "aligned: skipped, private/ is not present"
else
vendored_cpp: vendored.cpp graphs.h $(VENDOR_SRCS)
	$(CXX) -std=c++17 -O3 -w vendored.cpp $(VENDOR_SRCS) -o $@

vendored: vendored_cpp

# The width study. A SECOND copy of the vendored AMD is generated here with its one `using Int`
# line retyped to int64_t and its entry point renamed, so both widths can be linked into one
# program and timed against each other. Generated rather than committed: it is a mechanical
# transform of a file that is not ours and not tracked, and a stale copy of private/ living in the
# experiment would be worse than no copy at all.
amd_int64.cpp: ../../private/Amd.cpp
	@sed -e 's/using Int = int32_t;/#include <vector>\nusing Int = int64_t;/' \
	     -e 's/^extern "C" int amd_order(/extern "C" int amd_order64(/' $< > $@.tmp
	@python3 -c "import sys;\
p='$@.tmp';s=open(p).read();\
old='    return AMD_order(n, Ap, Ai, P, Control, Info);';\
new='    std::vector<Int> ap(n + 1), ai(Ap[n]), p_(n);\n\
    for (int32_t k = 0; k <= n; ++k) ap[k] = Ap[k];\n\
    for (int32_t k = 0; k < Ap[n]; ++k) ai[k] = Ai[k];\n\
    const int status = AMD_order(static_cast<Int>(n), ap.data(), ai.data(), p_.data(), Control, Info);\n\
    for (int32_t k = 0; k < n; ++k) P[k] = static_cast<int32_t>(p_[k]);\n\
    return status;';\
assert s.count(old)==1, 'the entry point moved; width.cpp needs updating';\
open('$@','w').write(s.replace(old,new))"
	@rm -f $@.tmp

# genmmd has no `using Int` line: it is plain `int` throughout, all of them the arithmetic type,
# so its transform is a blanket rename to a typedef rather than a one-line change. Mechanical for
# that file, and it would not be for one mixing widths.
mmd_int64.cpp: ../../private/Mmd.cpp
	@python3 -c "import re;\
s=open('$<').read();\
s=re.sub(r'\\bint\\b','MmdInt',s);\
assert s.count('void mmd_order64_inner(')==0 and s.count('void mmd_order(')==1, 'the entry point moved';\
s=s.replace('void mmd_order(','void mmd_order64_inner(',1);\
open('$@','w').write('#include <cstdint>\n#include <vector>\nusing MmdInt = int64_t;\n'+s+'''\n\
extern \"C\" void mmd_order64(std::int32_t n, const std::int32_t colPtr[], const std::int32_t rowIdx[],\n\
                            std::int32_t perm[], std::int32_t invp[])\n\
{\n\
    if (n <= 0) return;\n\
    const std::int32_t nnz = colPtr[n];\n\
    std::vector<MmdInt> cp(n + 1), ri(nnz), p(n), ip(n);\n\
    for (std::int32_t k = 0; k <= n; ++k) cp[k] = colPtr[k];\n\
    for (std::int32_t k = 0; k < nnz; ++k) ri[k] = rowIdx[k];\n\
    mmd_order64_inner(n, cp.data(), ri.data(), p.data(), ip.data());\n\
    for (std::int32_t k = 0; k < n; ++k) { perm[k] = static_cast<std::int32_t>(p[k]);\n\
                                           invp[k] = static_cast<std::int32_t>(ip[k]); }\n\
}\n''')"

# ---- the two alignment checks -----------------------------------------------------------------
#
# `make amdorder` checks production Amd3 against the vendored AMD's RAW elimination order: the
# pivot sequence together with the member order inside each supervariable, which is the whole
# algorithm. It is the acceptance test for the amd alignment and the README's method section
# explains why it, and not amd_order's output vector, is the object to compare.
#
# AMD_2 does not emit that order. amd_order returns Perm, which AMD_1 has already relabeled by
# AMD_postorder, and that postorder is a heuristic tidy of an approximate assembly tree which
# ElmForest replaces with Liu's rule on the exact supernodal tree. So the vendored OUTPUT differs
# from ours by construction and always should.
#
# A hooked copy is GENERATED here rather than kept, exactly as the int64 copies above are, and for
# the same reason: a checked-in copy of vendored code with our edits in it would drift from
# private/ silently, and a test comparing against a stale oracle is worse than no test. Generated,
# gitignored, removed by `make clean`, and rebuilt from whatever private/Amd.cpp currently says.
#
# The hook is additive: two containers and three lines, no behaviour change. PB_members[i]
# accumulates what supervariable i stands for, growing at the hash merge; PB_raw takes a pivot's
# members once per elimination, at FINALIZE THE NEW ELEMENT, before AMD_postorder ever runs.
#
# It runs FOUR SHAPES, not one shape at many sizes: the seven examples, 2D grids, 3D grids and
# random patterns. The 2D-only version of this list was green while Amd3 carried a stale clique
# degree that a 3D grid at 16 a side finds, which is what the widening is for. It takes under a
# second, so there is no reason to run less of it.
#
# `graphs.h` is a prerequisite so that editing it rebuilds, and `$(filter-out %.h,$^)` keeps it
# out of the compiler's argument list, which is not optional and not tidiness. A header named as
# a source is a header to PRECOMPILE, so the command then generates two outputs and `-o` is no
# longer meaningful: Apple Clang says "cannot specify -o when generating multiple output files"
# and refuses, while g++ accepts it silently and builds. So this breaks on macOS and passes on
# Linux, which is the compiler asymmetry CLAUDE.md's Build section names, in the direction that
# catches you out less often.
amd_raw.cpp: ../../private/Amd.cpp hook_amd.py
	@python3 hook_amd.py $< $@

amdorder_cpp: amdorder.cpp graphs.h $(VENDOR_SRCS) amd_raw.cpp ../../src/Amd3.cpp ../../src/QuotientGraph.cpp ../../src/Types.cpp
	$(CXX) -std=c++17 -O3 -DNDEBUG -w -I../../include $(filter-out %.h,$^) -o $@

amdorder: amdorder_cpp
	./amdorder_cpp

# The mmd half, and it needs none of the machinery above. genmmd emits the order it eliminates in
# and does no postorder, so `mmd_order`'s own output vector is the object to compare: no hook, no
# generated copy, no anchors, no Control array. MMD GETS THE ORDER DIRECTLY AND AMD GETS IT THROUGH
# A HOOK, which is genmmd's asymmetry rather than ours. The two targets are named for the BRANCH
# they check and not for that difference, which is written down in both drivers instead: naming one
# for how it works and the other for what it checks made a matched pair read as two different kinds
# of thing. This target was `raworder` until 2026-08-09 for exactly that reason.
mmdorder_cpp: mmdorder.cpp graphs.h $(VENDOR_SRCS) ../../src/Mmd3.cpp ../../src/QuotientGraph.cpp ../../src/Types.cpp
	$(CXX) -std=c++17 -O3 -DNDEBUG -w -I../../include $(filter-out %.h,$^) -o $@

mmdorder: mmdorder_cpp
	./mmdorder_cpp

# One word for "is either ordering still aligned", which is the question a reader has. The same
# shape as `make scale` over `scale-mmd` and `scale-amd`.
aligned: amdorder mmdorder

width_cpp: width.cpp $(VENDOR_SRCS) amd_int64.cpp mmd_int64.cpp
	$(CXX) -std=c++17 -O3 -DNDEBUG -w width.cpp $(VENDOR_SRCS) amd_int64.cpp mmd_int64.cpp -o $@

width: width_cpp
	./width_cpp
endif

clean:
	rm -f $(BINS) vendored_cpp production_cpp width_cpp amdorder_cpp mmdorder_cpp \
	      amd_int64.cpp mmd_int64.cpp amd_raw.cpp Amd.o Mmd.o
