# storage-options experiment - one multiply function over two sparse storage formats.
# Reference/teaching only; not part of the main Oblio build.
#
#:   make        print this list (the default goal; `make help` is the same thing)
#:   make all    build the test executable
#:   make test   build and run it
#:   make clean
#:   make help   print this list
#
# The executable carries the _cpp language suffix (coding convention) and is gitignored.
#
# NDEBUG is on deliberately: this experiment measures performance, so asserts are compiled
# out (release build). -O3 matters: the claim is that reaching a column through a pointer
# array costs nothing against hand-written CSC, and at -O0 that would measure nothing.
#
# No BLAS: the kernel is a sparse column scatter, not a dense block operation. (The
# friend-access experiment is the one that hands a raw block to BLAS.)

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra

SRCS = MultiplyEngine.cpp SparseMatrixStatic.cpp SparseMatrixDynamic.cpp
HDRS = SparseMatrixStatic.h SparseMatrixDynamic.h MultiplyEngine.h

.PHONY: all test clean help

# Bare `make` prints the target list rather than building. Stated explicitly because the default
# goal is otherwise the FIRST target in the file, so moving a rule would silently change it and the
# header above would stop being true with nothing to catch it. Every Makefile in this tree answers
# `make` the same way; see docs/CODING_RULES.md.
.DEFAULT_GOAL := help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: test_multiply_cpp

test_multiply_cpp: test_multiply.cpp $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) test_multiply.cpp $(SRCS) -o $@

test: all
	./test_multiply_cpp

# The *.dSYM line removes macOS debug-symbol bundles, which are directories, hence -rf rather
# than -f. Nothing here emits one under the committed flags, and that is not the test: CXXFLAGS
# is overridable, `make CXXFLAGS="... -g"` produces one on macOS, and clean removes what a build
# in this directory CAN produce.
clean:
	rm -f test_multiply_cpp
	rm -rf *.dSYM
