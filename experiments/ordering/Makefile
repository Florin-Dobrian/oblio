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
#   md3   the quotient graph: cliques instead of fill               5.3, 5.4
#   md4   supervariables and mass elimination                       5.5, 5.6
#   md5   maintained degrees, refreshed only where they changed     5.7, 5.8
#   md6   degree buckets, so the minimum is walked to not scanned   5.9, 5.10
#   mmd   multiple elimination: a batch of pivots per refresh       5.11, 5.12
#   amd   approximate degree: a bound instead of a set union        5.13
#
# md1 through md6 all return the SAME ordering; mmd and amd do not, deliberately.
# See archive/sparse_factorization.md section 5 for the sections cited above.
#
# md2 is an earlier naive variant, kept for reference and not built.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra

LAYERS = md1 md3 md4 md5 md6 mmd amd
BINS   = $(addsuffix _cpp,$(LAYERS))

.PHONY: all test clean

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
	  python3 $$L.py | tr -d " []{},'" > .$$L.py.out; \
	  ./$${L}_cpp   | tr -d " []{},'" > .$$L.cpp.out; \
	  if diff .$$L.py.out .$$L.cpp.out > /dev/null; then \
	    echo "py and cpp agree"; \
	  else \
	    echo "DIFFER"; rc=1; \
	  fi; \
	  rm -f .$$L.py.out .$$L.cpp.out; \
	done; \
	exit $$rc

clean:
	rm -f $(BINS)
