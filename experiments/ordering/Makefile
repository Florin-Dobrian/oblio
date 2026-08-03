# ordering experiment - the minimum degree family, built up one mechanism at a time.
# Reference/teaching only; not part of the main Oblio build.
#
#   make          build every C++ prototype
#   make test     build them, run each, check the C++ agrees with its Python twin, and check
#                 each ported layer agrees with the production driver extracted from it
#   make production   the production drivers on these graphs, linked from ../../src
#   make clean
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
LAYERS = md1 md2 md3 md4 md5 mmd1 mmd2 amd1 amd2 mdm2 mda2 mdam2
BINS   = $(addsuffix _cpp,$(LAYERS))

# The layers that have been pulled into the main tree, and so have a production driver to be
# checked against. Grows as the extraction proceeds; empty would simply skip that check.
# Layers checked by permutation, and layers checked by fill. AMD2 is in the second group
# because production skips the postorder the prototype ends with, so the permutations
# legitimately differ while the fill does not.
#
# The prototype's matrix1 is dropped from the comparison. It is the one example given as a matrix
# rather than a graph, unsymmetric with duplicates and an unsorted column, and it exists to
# exercise amd2Aat and amd2Preprocess. Production has no input conditioning on purpose, Oblio's
# matrices being clean by construction, so there is nothing on this side to compare it against.
PORTED      = mmd1 mmd2 amd1
PORTED_FILL = amd2

# The production sources those drivers need, linked directly rather than copied. A copy is right
# for the vendored routines, which are not ours to edit, and wrong here: the point of this check
# is to notice when production and prototype come apart.
# The fill check runs Oblio's symbolic factorization, so the whole library is linked. The
# two vendored orderers compile with warnings off, as everywhere else.
PROD_VENDOR = ../../src/Amd.cpp ../../src/Mmd.cpp
PROD_SRCS   = $(filter-out $(PROD_VENDOR),$(wildcard ../../src/*.cpp))

.PHONY: all test production vendored clean

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
production_cpp: production.cpp $(PROD_SRCS) $(PROD_VENDOR)
	$(CXX) $(CXXFLAGS) -I../../include -DOBLIO_BLAS_UNDERSCORE -w -c $(PROD_VENDOR)
	$(CXX) $(CXXFLAGS) -I../../include -DOBLIO_BLAS_UNDERSCORE production.cpp $(PROD_SRCS) \
	  Amd.o Mmd.o $(BLAS_LIBS) -o $@

production: production_cpp

# The vendored routines, built separately: they are not layers and have no Python
# twin. vendored/vendored_mmd.cpp and vendored/vendored_amd.cpp are copies of
# src/Mmd.cpp and src/Amd.cpp, never edited. The lowercase prefixed names matter:
# macOS is case-insensitive by default, so Amd.cpp and amd.cpp would be one path.
# Warnings are suppressed for these two only: they are not ours to clean up.
vendored_cpp: vendored.cpp vendored/vendored_mmd.cpp vendored/vendored_amd.cpp
	$(CXX) -std=c++17 -O3 -w $^ -o $@

vendored: vendored_cpp

clean:
	rm -f $(BINS) vendored_cpp production_cpp Amd.o Mmd.o
