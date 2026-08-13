# Two independent dense products, run back to back and then side by side.
#
#:   make            print this list (the default goal; `make help` is the same thing)
#:   make all        build all four binaries
#:   make test       build them and run the comparison
#:   make example    build and run the example pair alone
#:   make clean
#:   make help       print this list
#
# Two binaries from one source, which is the point of the second one:
#
#   test_parallel_cpp          built with OpenMP
#   test_parallel_serial_cpp   built without it, and still correct
#
# The serial build is the degradation claim in docs/DESIGN_DECISIONS.md made
# checkable. It also needs -Wno-unknown-pragmas: the pragmas are correctly ignored,
# but -Wall reports every one of them, and the tree is otherwise -Wall -Wextra clean.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -DNDEBUG

# Same BLAS selection as the main tree's Makefile.
CXXFLAGS += -DOBLIO_BLAS_UNDERSCORE

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  BLAS_LIBS = -framework Accelerate
else
  BLAS_LIBS = -llapack -lblas
endif

# OpenMP is a compiler feature everywhere except on Apple's Clang, which ships no
# runtime, so macOS needs libomp from Homebrew and the include and library paths that
# come with it. Set OMP_PREFIX if Homebrew lives somewhere other than /opt/homebrew.
ifeq ($(UNAME),Darwin)
  OMP_PREFIX ?= /opt/homebrew/opt/libomp
  OMP_FLAGS   = -Xpreprocessor -fopenmp -I$(OMP_PREFIX)/include
  OMP_LIBS    = -L$(OMP_PREFIX)/lib -lomp
else
  OMP_FLAGS   = -fopenmp
  OMP_LIBS    = -fopenmp
endif

# The cap that keeps a self-threading BLAS from oversubscribing the cores when two
# fronts run at once. Accelerate reads VECLIB_MAXIMUM_THREADS; OpenBLAS reads its own.
ifeq ($(UNAME),Darwin)
  BLAS_CAP = VECLIB_MAXIMUM_THREADS=1
else
  BLAS_CAP = OPENBLAS_NUM_THREADS=1
endif

SRCS = Gemm.cpp test_parallel.cpp

all: test_parallel_cpp test_parallel_serial_cpp example_parallel_cpp example_parallel_serial_cpp

test_parallel_cpp: $(SRCS) Gemm.h
	$(CXX) $(CXXFLAGS) $(OMP_FLAGS) $(SRCS) $(OMP_LIBS) $(BLAS_LIBS) -o $@

test_parallel_serial_cpp: $(SRCS) Gemm.h
	$(CXX) $(CXXFLAGS) -Wno-unknown-pragmas $(SRCS) $(BLAS_LIBS) -o $@

test: all
	@echo "=== serial build, no OpenMP: correct, and a speedup of one ==="
	@./test_parallel_serial_cpp
	@echo
	@echo "=== OpenMP, BLAS uncapped: two fronts against a self-threading BLAS ==="
	@./test_parallel_cpp
	@echo
	@echo "=== OpenMP, BLAS capped to one thread: the honest one-against-two ==="
	@$(BLAS_CAP) ./test_parallel_cpp

# The teaching example, built and run separately from the measurement above. No BLAS: it
# needs only the compiler, which is the point, since it is about the pragmas rather than
# about this machine. Built both ways for the same reason the test is.
example_parallel_cpp: example_parallel.cpp
	$(CXX) $(CXXFLAGS) $(OMP_FLAGS) example_parallel.cpp $(OMP_LIBS) -o $@

example_parallel_serial_cpp: example_parallel.cpp
	$(CXX) $(CXXFLAGS) -Wno-unknown-pragmas example_parallel.cpp -o $@

example: example_parallel_cpp example_parallel_serial_cpp
	@echo "=== without OpenMP: pragmas ignored, answers still correct ==="
	@./example_parallel_serial_cpp
	@echo
	@echo "=== with OpenMP ==="
	@./example_parallel_cpp

# The *.dSYM line removes macOS debug-symbol bundles, which are directories, hence -rf rather
# than -f. Nothing here emits one under the committed flags, and that is not the test: CXXFLAGS
# is overridable, `make CXXFLAGS="... -g"` produces one on macOS, and clean removes what a build
# in this directory CAN produce.
clean:
	rm -f test_parallel_cpp test_parallel_serial_cpp
	rm -f example_parallel_cpp example_parallel_serial_cpp
	rm -rf *.dSYM

.PHONY: all test example clean help

# Bare `make` prints the target list rather than building. Stated explicitly because the default
# goal is otherwise the FIRST target in the file, so moving a rule would silently change it and the
# header above would stop being true with nothing to catch it. Every Makefile in this tree answers
# `make` the same way; see docs/CODING_RULES.md.
.DEFAULT_GOAL := help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-
