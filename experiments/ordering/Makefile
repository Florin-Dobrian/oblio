# ordering experiment - the minimum degree family, built up one mechanism at a time.
# Reference/teaching only; not part of the main Oblio build.
#
#:   make          build every C++ prototype
#:   make test     build them, run each, check the C++ agrees with its Python twin, and check
#:                 each ported layer agrees with the production driver extracted from it
#:   make production   the production drivers on these graphs, linked from ../../src
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
LAYERS = md1 md2 md3 md4 md5 mmd1 mmd2 amd1 amd2 amd3 mdm2 mda2 mdam2
BINS   = $(addsuffix _cpp,$(LAYERS))

# The layers that have been pulled into the main tree, and so have a production driver to be
# checked against. Grows as the extraction proceeds; empty would simply skip that check.
# Layers checked by permutation, and layers checked by fill. PORTED_FILL is empty and the weaker
# check exists only in case a future layer needs it: amd2 was in it until the postorder moved to
# amd3, and every ported layer now matches production entry for entry.
#
# amd3 is not here at all and is not expected to be. Its four items are input conditioning, the
# postorder and a control interface, none of which Oblio has or wants, so there is no production
# counterpart to compare it against. Its oracle is its own C++ twin. That is also why the
# prototype's matrix1 has nothing to compare against on this side: it is the one example given as
# a matrix rather than a graph, unsymmetric with duplicates and an unsorted column, and it exists
# to exercise amd3Aat and amd3Preprocess.
PORTED      = mmd1 mmd2 amd1 amd2
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
# examples cannot reach. amd3 is here although it is not in PORTED: it has no production
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

.PHONY: all test production vendored clean help

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
production_cpp: production.cpp $(PROD_SRCS)
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
else
vendored_cpp: vendored.cpp $(VENDOR_SRCS)
	$(CXX) -std=c++17 -O3 -w vendored.cpp $(VENDOR_SRCS) -o $@

vendored: vendored_cpp
endif

clean:
	rm -f $(BINS) vendored_cpp production_cpp Amd.o Mmd.o
