# matrices benchmark - accuracy of the pipeline on real matrices, from the SuiteSparse
# Matrix Collection. Nothing here is timed.
#
#:   make        build both programs
#:   make read   build and run the reader alone over data/, an inventory of what is usable
#:   make run    build and run the residual pass over everything in data/
#:   make clean
#:   make help   print this list
#:
#: The matrices are NOT in the repository. Get them first, from this directory. This is the
#: command that produced the committed candidates.txt:
#:   ./ssget.py list --per-kind 8 --values --max-nnz 500000 > candidates.txt
#:   ./ssget.py fetch candidates.txt
#: They land in ../../data/<Group>/<name>.mtx, which is gitignored. `list` downloads nothing but
#: the collection's index, so it is cheap to re-run until the candidate set looks right, and
#: `fetch` skips what is already on disk.
#:
#: ./mmread_cpp ../../data/*/*.mtx                     what each file is, and what is refused
#: ./matrix_residual_cpp ../../data/HB/*.mtx           any subset, the shell choosing it
#: ./matrix_residual_cpp --max-fill=2e8 ../../data/*/*.mtx    raise the predicted-fill cap
#:
#: ACCURACY.md is the report; README.md is the working notes behind it.
#
# Not part of the main Oblio build, and unlike experiments/ this is not a frozen study: a
# benchmark runs against the current tree forever, so it links ../../src directly and must keep
# compiling as the code moves. Building it is cheap and running it is not, which is why `make`
# builds and only the run targets measure.
#
# NDEBUG is on for the same reason it is in the sibling folders, so that what runs here is what a
# caller would run.
#
# THIS FOLDER HAS NO VENDORED-ORDERING BRANCH, where ../ordering and ../pipeline both do. Nothing
# here selects Ordering::MMD or Ordering::AMD, so private/ changes nothing about what it builds or
# prints, and the OBLIO_PUBLIC switch would be a mechanism with no effect. Add both together if a
# vendored ordering is ever wanted as a column.
#
# BLAS: Accelerate on macOS, the system BLAS on Linux. mmread does not factor anything and could
# link without it, but it shares the library sources with the driver that does, so one rule covers
# both rather than two rules differing in a way nobody would remember.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra -I../../include -DOBLIO_BLAS_UNDERSCORE

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    BLAS_LIBS ?= -framework Accelerate
else
    BLAS_LIBS ?= -lblas -llapack -lm
endif

# The whole library, globbed rather than listed: naming the few units a solve needs would be a
# second inventory to keep in step by hand.
OBLIO_SRCS = $(wildcard ../../src/*.cpp)

# The headers are prerequisites too. Leaving them out was a real defect in ../ordering: a change
# to a header alone left the run targets reporting the OLD binary, silently, every source being
# older than the target.
OBLIO_HDRS = $(wildcard ../../include/oblio/*.h)

MATRICES = ../../data/*/*.mtx

# The run targets pass no options, so `run` uses the driver's own defaults: a predicted-fill
# cap of 50 million entries and no address-space limit on the child processes. Both are
# command-line flags; call the binary directly to change them.

# `make` builds, as the comment above says it does. Without this, the default goal is the first
# rule in the file, which is `help`, and a bare `make` would print the usage block and build
# nothing. ../ordering and ../pipeline have the same shape and therefore the same behavior, which
# is worth knowing when reading CLAUDE.md's clone checklist.
.DEFAULT_GOAL := all

.PHONY: all read run clean help

help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: mmread_cpp matrix_residual_cpp

mmread_cpp: mmread.cpp MatrixMarket.h $(OBLIO_SRCS) $(OBLIO_HDRS)
	$(CXX) $(CXXFLAGS) $< $(OBLIO_SRCS) $(BLAS_LIBS) -o $@

matrix_residual_cpp: matrix_residual.cpp MatrixMarket.h $(OBLIO_SRCS) $(OBLIO_HDRS)
	$(CXX) $(CXXFLAGS) $< $(OBLIO_SRCS) $(BLAS_LIBS) -o $@

# Both run targets take the whole of data/ through the shell rather than through the program, so
# any subset is a matter of typing a different glob and needs no option here.
read: mmread_cpp
	@./mmread_cpp $(MATRICES)

run: matrix_residual_cpp
	@./matrix_residual_cpp $(MATRICES)

clean:
	rm -f mmread_cpp matrix_residual_cpp
	rm -rf *.dSYM
