# pipeline benchmark - where a solve's time goes, across all its phases.
#
#:   make          build both benchmarks
#:   make run2d    where a solve's time goes, the default square-grid sizes, every ordering
#:   make run3d    the same on cubic grids, which are a DIFFERENT PROBLEM FAMILY
#:   make scale2d  how the cost GROWS, a longer square ladder, every factorization
#:   make scale3d  the same, cubic
#:   make clean
#:   make help     print this list
#:
#: ./pipeline_timing_cpp 200            any square sides
#: ./pipeline_timing_cpp 3d 26          any cubic sides
#: ./pipeline_scaling_cpp 2d 100 200    any square sides
#: ./pipeline_scaling_cpp 3d 16 20      any cubic sides
#:
#: run2d and run3d use ordering/'s ladders exactly, 32, 64, 100, 140 on squares and 6, 12, 16,
#: 20, 26 on cubes, so a row here can be read beside a row there without converting sizes. The
#: SCALE ladders are longer and are their own, 100 to 600 and 16 to 48, and their matrices carry a
#: dominance margin the short ladders' do not, so rows are NOT comparable across the two. There is
#: no bare `make run` in either folder: the family is always named, because a claim measured on
#: one of them is a claim about that one.
#:
#: Prefix any of these with OBLIO_PUBLIC=1 to build as everyone else does, without the vendored
#: orderings; see docs/DESIGN_DECISIONS.md. The same word works in the repo root, both benchmarks
#: and experiments/ordering.
#
# The sibling folder, ordering/, measures one phase against itself. This one measures the phases
# against each other, which is what says whether optimizing the first is worth an afternoon.
#
# Not part of the main Oblio build, and unlike experiments/ this is not a frozen study: a
# benchmark runs against the current tree forever, so it links ../../src directly and must keep
# compiling as the code moves. Building it is cheap and running it is not, which is why `make`
# builds and only `make run2d` and `make run3d` measure.
#
# NDEBUG is on deliberately, as in the experiments that measure: the claim is about what a solve
# costs, and at -O0 with asserts in it would measure something else.
#
# BLAS: Accelerate on macOS, the system BLAS on Linux. Unlike the ordering benchmark, this one
# factors, so the BLAS is in the measurement rather than only on the link line.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra -I../../include -DOBLIO_BLAS_UNDERSCORE

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    BLAS_LIBS ?= -framework Accelerate
else
    BLAS_LIBS ?= -lblas -llapack -lm
endif

# The whole library, which this benchmark genuinely needs end to end. Globbed rather than listed,
# so a new unit never has to be added here by hand.
# The two vendored ordering codes are separate for two reasons: they are not ours to clean up, so
# they compile with warnings off exactly as the root Makefile does it, and they are not published,
# so they live in ../../private and are optional. Absent, this benchmark still builds and any row
# that asks for Ordering::MMD or Ordering::AMD reports a refusal.
OBLIO_SRCS  = $(wildcard ../../src/*.cpp)

# See benchmarks/ordering/Makefile: without this a header-only change is invisible to make and the
# benchmark reports a stale binary as if it were current.
OBLIO_HDRS  = $(wildcard ../../include/oblio/*.h)
# OBLIO_PUBLIC=1 builds as everyone else does, ignoring private/ even when it is there. Same word
# in every Makefile here that links the library, so one habit covers all of them.
ifdef OBLIO_PUBLIC
  VENDOR_SRCS =
else
  VENDOR_SRCS = $(wildcard ../../private/Amd.cpp ../../private/Mmd.cpp)
endif
# Which configuration this directory is built in. The two differ by a macro, so a binary from one
# is wrong for the other and no timestamp would say so. The $(shell) deletes the stamp when the
# recorded mode is not the current one, which rebuilds what depends on it; switching therefore
# rebuilds and repeating does not.
BUILD_MODE := $(if $(VENDOR_SRCS),vendored,public)
$(shell [ "`cat .build-mode 2>/dev/null`" = "$(BUILD_MODE)" ] || rm -f .build-mode)

.build-mode:
	@echo "$(BUILD_MODE)" > $@

ifneq ($(VENDOR_SRCS),)
  VENDOR_DEFS = -DOBLIO_VENDORED_ORDERINGS
  VENDOR_OBJS = Amd.o Mmd.o
endif

# `make` builds, as the header says it does. Without this the default goal is the first rule in
# the file, which is `help`, and a bare `make` would print the usage block and build nothing.
# CLAUDE.md's pre-release clone checklist runs `make` here and expects a build.
.DEFAULT_GOAL := all

.PHONY: all run2d run3d scale2d scale3d clean help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: pipeline_timing_cpp pipeline_scaling_cpp

pipeline_timing_cpp: pipeline_timing.cpp $(OBLIO_SRCS) $(OBLIO_HDRS) $(VENDOR_SRCS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) $(VENDOR_DEFS) $< $(OBLIO_SRCS) $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

# The scaling driver, which measures growth rather than a share: every factorization Oblio has,
# real and complex Hermitian and complex symmetric, over a longer ladder. See its header.
pipeline_scaling_cpp: pipeline_scaling.cpp $(OBLIO_SRCS) $(OBLIO_HDRS) $(VENDOR_SRCS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) $(VENDOR_DEFS) $< $(OBLIO_SRCS) $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

run2d: pipeline_timing_cpp
	./pipeline_timing_cpp 2d

run3d: pipeline_timing_cpp
	./pipeline_timing_cpp 3d

scale2d: pipeline_scaling_cpp
	./pipeline_scaling_cpp 2d

scale3d: pipeline_scaling_cpp
	./pipeline_scaling_cpp 3d

clean:
	rm -f .build-mode pipeline_timing_cpp pipeline_scaling_cpp Amd.o Mmd.o
