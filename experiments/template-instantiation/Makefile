# Template-instantiation example - Matrix/Vector/MultiplyEngine across the three variants.
# Reference/teaching only; not part of the main Oblio build.
#
#:   make         print this list (the default goal; `make help` is the same thing)
#:   make all     build all three variants
#:   make test    build and run all three (they must print identical results)
#:   make clean
#:   make help    print this list
#
# Executables carry the _cpp language suffix (coding convention), so *_cpp is
# gitignored and never committed.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra

# .cpp sources to link per explicit variant. (Implicit is header-only, it links
# nothing but the test itself, so it has no source list.)
PLAIN_EXPLICIT_SRCS   = MatrixPlainExplicit.cpp VectorPlainExplicit.cpp MultiplyEnginePlainExplicit.cpp
GUARDED_EXPLICIT_SRCS = MatrixGuardedExplicit.cpp VectorGuardedExplicit.cpp MultiplyEngineGuardedExplicit.cpp

.PHONY: all test clean help

# Bare `make` prints the target list rather than building. Stated explicitly because the default
# goal is otherwise the FIRST target in the file, so moving a rule would silently change it and the
# header above would stop being true with nothing to catch it. Every Makefile in this tree answers
# `make` the same way; see CODING_RULES.md.
.DEFAULT_GOAL := help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: test_multiply_implicit_cpp test_multiply_plain_explicit_cpp test_multiply_guarded_explicit_cpp

# Implicit - body in header; nothing to link but the test
test_multiply_implicit_cpp: test_multiply.cpp \
                            MatrixImplicit.h VectorImplicit.h MultiplyEngineImplicit.h
	$(CXX) $(CXXFLAGS) -DOBLIO_TI_IMPLICIT test_multiply.cpp -o $@

# Plain explicit - bodies in .cpp, header signatures only
test_multiply_plain_explicit_cpp: test_multiply.cpp $(PLAIN_EXPLICIT_SRCS) \
                                  MatrixPlainExplicit.h VectorPlainExplicit.h MultiplyEnginePlainExplicit.h
	$(CXX) $(CXXFLAGS) -DOBLIO_TI_PLAIN_EXPLICIT test_multiply.cpp $(PLAIN_EXPLICIT_SRCS) -o $@

# Guarded explicit - plain explicit + extern template guard
test_multiply_guarded_explicit_cpp: test_multiply.cpp $(GUARDED_EXPLICIT_SRCS) \
                                    MatrixGuardedExplicit.h VectorGuardedExplicit.h MultiplyEngineGuardedExplicit.h
	$(CXX) $(CXXFLAGS) -DOBLIO_TI_GUARDED_EXPLICIT test_multiply.cpp $(GUARDED_EXPLICIT_SRCS) -o $@

test: all
	./test_multiply_implicit_cpp
	./test_multiply_plain_explicit_cpp
	./test_multiply_guarded_explicit_cpp

# The *.dSYM line removes macOS debug-symbol bundles, which are directories, hence -rf rather
# than -f. Nothing here emits one under the committed flags, and that is not the test: CXXFLAGS
# is overridable, `make CXXFLAGS="... -g"` produces one on macOS, and clean removes what a build
# in this directory CAN produce.
clean:
	rm -f test_multiply_implicit_cpp test_multiply_plain_explicit_cpp test_multiply_guarded_explicit_cpp
	rm -rf *.dSYM
