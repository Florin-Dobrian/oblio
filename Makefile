# Builds the Oblio units and every test in tests/. A simple alternative to the
# CMake build (CMakeLists.txt); both compile the same src/ and tests/.
#
#:   make            print this list (the default goal; `make help` is the same thing)
#:   make all        build everything (tests and examples)
#:   make tests      build the test executables
#:   make test       build and run the tests, then run the examples for exit status
#:   make examples   build the example programs (examples/*.cpp)
#:   make objs       compile the library sources to .o only (a fast core compile check)
#:   make <name>_cpp build one test or example (e.g. make test_order_cpp,
#:                   make example_basic_cpp)
#:   make clean
#:   make help       print this list
#
#: Prefix any of these with OBLIO_PUBLIC=1 to build as everyone else does, without the
#: vendored orderings; see docs/DESIGN_DECISIONS.md. The same works in benchmarks/*/ and
#: experiments/ordering/.
#:
#
# -Iinclude points at include/, so #include "oblio/X.h" resolves to the project
# headers. Executables carry the _cpp suffix (coding convention) and are gitignored.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -Iinclude

# BLAS and LAPACK, for the numeric factorization.
#
# OBLIO_BLAS_UNDERSCORE selects the trailing-underscore Fortran symbol convention
# (dpotrf_ rather than dpotrf), which is what Accelerate and the reference BLAS both
# use. See include/oblio/BlasLapack.h.
CXXFLAGS += -DOBLIO_BLAS_UNDERSCORE

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  BLAS_LIBS = -framework Accelerate
else
  BLAS_LIBS = -llapack -lblas
endif

# Oblio units, full warnings apply.
#
# Mmd3B IS IN THIS LIST AND IS TEMPORARY, 2026-08-16. It is the one remaining B layer, Mmd3 on the
# vendored clique storage scheme, and its whole obligation is to reproduce Mmd3's permutation entry
# for entry. tests/test_order.cpp asserts exactly that, so it has to link into the test binaries,
# and a check that never links is not a check. It leaves this list with the file.
#
# The two amd B layers retired on the same day were never here: they were reached only through the
# benchmark's own glob of src/*.cpp, which is why their pair check, though written, never ran under
# `make test` at all. That is worth knowing before another B layer is added.
#
# NO COMMENTS INSIDE THE LIST BELOW. It is one assignment continued with backslashes, so a comment
# line between entries is swallowed by the continuation and make reports a missing separator.
OBLIO_SRCS = \
  src/Types.cpp \
  src/SparseMatrix.cpp \
  src/Permutation.cpp \
  src/OrderEngine.cpp \
  src/QuotientGraph.cpp \
  src/Mmd1.cpp \
  src/Mmd2.cpp \
  src/Mmd3.cpp \
  src/Mmd3B.cpp \
  src/Amd1.cpp \
  src/Amd2.cpp \
  src/Amd3.cpp \
  src/ElmForestEngine.cpp \
  src/SymFactorEngine.cpp \
  src/BlasLapack.cpp \
  src/UpdateBlock.cpp \
  src/UpdateMatrix.cpp \
  src/NumFactorStatic.cpp \
  src/NumFactorDynamic.cpp \
  src/NumFactorEngine.cpp \
  src/Vector.cpp \
  src/MultiplyEngine.cpp \
  src/SolveEngine.cpp \
  src/DirectSolver.cpp

# Vendored ordering codes (SuiteSparse AMD, Sparspak MMD), copied verbatim, not maintained here,
# so warnings are silenced (-w) for these two files only.
#
# They live in private/, which is not tracked, and are optional: the build detects them rather than
# requiring them. Present, Ordering::MMD and Ordering::AMD work as they always have. Absent, the
# library still builds and those two enumerators refuse, everything else being unaffected. Nothing
# is switched by hand and there is no flag to remember; a tree that has the directory behaves one
# way and a clone behaves the other. See docs/DESIGN_DECISIONS.md.
# OBLIO_PUBLIC=1 builds as everyone else does, ignoring private/ even when it is there. It is the
# same word in every Makefile in this repo that links the library, so one habit covers all of them.
ifdef OBLIO_PUBLIC
  VENDOR_SRCS =
else
  VENDOR_SRCS = $(wildcard private/Amd.cpp private/Mmd.cpp)
endif

ifneq ($(VENDOR_SRCS),)
  CXXFLAGS += -DOBLIO_VENDORED_ORDERINGS
endif

LIB_SRCS = $(OBLIO_SRCS) $(VENDOR_SRCS)
LIB_HDRS = $(wildcard include/oblio/*.h)

# Which configuration the build is in. The two differ by a macro, so objects and binaries from one
# are wrong for the other, and no timestamp would say so: the sources are identical either way.
# The $(shell) runs while this file is read and deletes the stamp when the recorded mode is not the
# current one, which makes it out of date, which rebuilds what depends on it. Switching therefore
# rebuilds and repeating does not, with no clean in between.
BUILD_MODE := $(if $(VENDOR_SRCS),vendored,public)
$(shell [ "`cat .build-mode 2>/dev/null`" = "$(BUILD_MODE)" ] || rm -f .build-mode)

.build-mode:
	@echo "$(BUILD_MODE)" > $@

# Compile each source to an object so per-file flags can differ. Oblio objects get
# full warnings; vendored objects get -w.
OBLIO_OBJS  = $(OBLIO_SRCS:.cpp=.o)
VENDOR_OBJS = $(VENDOR_SRCS:.cpp=.o)
LIB_OBJS    = $(OBLIO_OBJS) $(VENDOR_OBJS)

# One executable per tests/*.cpp file, named <stem>_cpp.
TEST_SRCS = $(wildcard tests/*.cpp)
TEST_BINS = $(patsubst tests/%.cpp,%_cpp,$(TEST_SRCS))

# One executable per examples/*.cpp file, named <stem>_cpp, exactly as for the tests. The example_
# prefix is carried by the source file itself (examples/example_basic.cpp), the way tests/ carries
# test_, so this rule adds only the _cpp suffix and both directories follow one convention.
EXAMPLE_SRCS = $(wildcard examples/*.cpp)
EXAMPLE_BINS = $(patsubst examples/%.cpp,%_cpp,$(EXAMPLE_SRCS))

.PHONY: all tests test examples objs clean help

# Bare `make` prints the target list rather than building, which is a decision and not an
# accident of ordering: a newcomer's first command in a repository root should say what is on
# offer, and the build everyone actually wants here is `make test` rather than `make all`.
# Stated explicitly because the default goal is otherwise the FIRST target in the file, so
# moving `help` below `all` would silently change it, and the header above would stop being
# true with nothing to catch it.
.DEFAULT_GOAL := help

# Print the target list from this file's header. The lines there are marked with #: so there is one
# source of truth: the header comment is the help text, and neither can drift from the other.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: tests examples

tests: $(TEST_BINS)

$(OBLIO_OBJS): %.o: %.cpp .build-mode $(LIB_HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(VENDOR_OBJS): %.o: %.cpp .build-mode
	$(CXX) $(CXXFLAGS) -w -c $< -o $@

%_cpp: tests/%.cpp $(LIB_OBJS) $(LIB_HDRS) $(wildcard tests/*.h)
	$(CXX) $(CXXFLAGS) tests/$*.cpp $(LIB_OBJS) $(BLAS_LIBS) -o $@

# Two pattern rules share the %_cpp target, one per source directory. Make picks the one whose
# prerequisite exists, and the stems cannot collide, tests/ being all test_* and examples/ all
# example_*.
%_cpp: examples/%.cpp $(LIB_OBJS) $(LIB_HDRS)
	$(CXX) $(CXXFLAGS) examples/$*.cpp $(LIB_OBJS) $(BLAS_LIBS) -o $@

# The suites first, each printing its own PASS/FAIL lines, then the examples. An example is
# documentation that compiles, and compiling is not the claim it makes: one of them once reported
# every dynamic cell as unimplemented long after they worked, and nothing caught it because nothing
# ran it. Running them here catches a crash, a hard refusal and a stale build, and nothing about
# whether the numbers are right; their output is discarded so it cannot drown the suites. Checking
# the numbers is the open half, in docs/TODO.md.
test: tests examples
	@pass=0; total=0; suites=0; \
	for t in $(TEST_BINS); do \
	  echo "== $$t =="; \
	  out=`./$$t` || { echo "$$out"; exit 1; }; \
	  echo "$$out"; echo; \
	  count=`echo "$$out" | grep -oE '[0-9]+/[0-9]+ passed' | tail -1`; \
	  p=`echo "$$count" | cut -d/ -f1`; \
	  n=`echo "$$count" | cut -d/ -f2 | cut -d' ' -f1`; \
	  pass=`expr $$pass + $$p`; total=`expr $$total + $$n`; \
	  suites=`expr $$suites + 1`; \
	done; \
	echo "== examples =="; \
	runs=0; \
	for e in $(EXAMPLE_BINS); do \
	  if ./$$e > /dev/null 2>&1; then echo "  PASS  $$e"; runs=`expr $$runs + 1`; \
	  else echo "  FAIL  $$e (exit status)"; exit 1; fi; \
	done; \
	echo; \
	echo "== total =="; \
	echo "  $$pass/$$total assertions across $$suites suites, $$runs examples run"; \
	if [ "$$pass" != "$$total" ]; then exit 1; fi
	@echo

examples: $(EXAMPLE_BINS)

# The library sources compiled to objects, nothing linked: a fast check that the core still builds.
objs: $(LIB_OBJS)

# private/*.o is named directly rather than through $(LIB_OBJS), which is empty of it under
# OBLIO_PUBLIC: clean should remove everything a build can leave behind regardless of which mode
# it is invoked in, so that the switch is never something to remember here.
# The *.dSYM line removes macOS debug-symbol bundles, which are directories, hence -rf rather
# than -f. Nothing here emits one under the committed flags, and that is not the test: CXXFLAGS
# is overridable, `make CXXFLAGS="... -g"` produces one on macOS, and clean removes what a build
# in this directory CAN produce.
clean:
	rm -f $(TEST_BINS) $(EXAMPLE_BINS) $(OBLIO_OBJS) private/Amd.o private/Mmd.o .build-mode
	rm -rf *.dSYM
