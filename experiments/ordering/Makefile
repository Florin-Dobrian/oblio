# ordering experiment - the minimum degree family, built up one mechanism at a time.
# Reference/teaching only; not part of the main Oblio build.
#
#   make        build every C++ prototype
#   make test   build them, run each, and check the C++ agrees with its Python twin
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

LAYERS = md1 md2 md3 md4 md5 mmd1 mmd2 amd1 amd2
BINS   = $(addsuffix _cpp,$(LAYERS))

.PHONY: all test vendored clean

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
test: all
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
	exit $$rc

# The vendored routines, built separately: they are not layers and have no Python
# twin. vendored/vendored_mmd.cpp and vendored/vendored_amd.cpp are copies of
# src/Mmd.cpp and src/Amd.cpp, never edited. The lowercase prefixed names matter:
# macOS is case-insensitive by default, so Amd.cpp and amd.cpp would be one path.
# Warnings are suppressed for these two only: they are not ours to clean up.
vendored_cpp: vendored.cpp vendored/vendored_mmd.cpp vendored/vendored_amd.cpp
	$(CXX) -std=c++17 -O3 -w $^ -o $@

vendored: vendored_cpp

clean:
	rm -f $(BINS) vendored_cpp
