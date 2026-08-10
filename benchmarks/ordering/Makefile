# ordering benchmark - what each ordering method costs, in time and in fill.
#
#:   make        build the benchmark
#:   make run2d    build and run it, the default square-grid sizes, every method
#:   make run3d    the same on cubic grids, which are a DIFFERENT PROBLEM FAMILY
#:   make scale2d  a wider ladder, one branch at a time, with the gap columns
#:   make scale3d  the same, cubic
#:                 Also scale-mmd-2d, scale-amd-2d, scale-mmd-3d, scale-amd-3d for one branch
#:   make clean
#:   make help   print this list
#:
#: ./order_profile_cpp amd3 140 200      one method in a loop, for a sampling profiler
#: ./order_profile_cpp amd3 3d 26 200    the same on a CUBIC grid: profile the family the
#:                                       question is about, since every profile taken here
#:                                       before 2026-08-09 was square and the amd hash pass
#:                                       looked diffuse there while doing eight times the work
#:                                       per pivot on a cube
#:
#: Prefix any of these with OBLIO_PUBLIC=1 to build as everyone else does, without the vendored
#: orderings; see docs/DESIGN_DECISIONS.md. The same word works in the repo root, both benchmarks
#: and experiments/ordering.
#
# Not part of the main Oblio build, and unlike experiments/ this is not a frozen study: a
# benchmark runs against the current tree forever, so it links ../../src directly and must keep
# compiling as the code moves. Building it is cheap and running it is not, which is why `make`
# builds and only the run targets measure.
#
# NDEBUG is on deliberately, as in the experiments that measure: the claim is about what an
# ordering costs, and at -O0 with asserts in it would measure something else.
#
# BLAS: Accelerate on macOS, the system BLAS on Linux. Nothing here calls BLAS, but the Oblio
# sources it links do, so the link line needs it.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra -I../../include -DOBLIO_BLAS_UNDERSCORE

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    BLAS_LIBS ?= -framework Accelerate
else
    BLAS_LIBS ?= -lblas -llapack -lm
endif

# The whole library. An ordering needs only a few units, but naming them here would be a second
# inventory to keep in step by hand, which is the failure the root build already fixed by globbing.
# The two vendored ordering codes are separate for two reasons: they are not ours to clean up, so
# they compile with warnings off exactly as the root Makefile does it, and they are not published,
# so they live in ../../private and are optional. Absent, this benchmark still builds and any row
# that asks for Ordering::MMD or Ordering::AMD reports a refusal.
OBLIO_SRCS  = $(wildcard ../../src/*.cpp)

# The headers are prerequisites too, and leaving them out was a real defect: a change to
# include/oblio/QuotientGraph.h alone left the run targets reporting the OLD binary, silently,
# because
# every source was older than the target. The library's own Makefile has always carried LIB_HDRS
# for this reason; the benchmarks link ../../src directly and had never been given the same.
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
  # The hooked copy, for the AMDraw column. Present only where private/ is, exactly as the vendored
  # objects are, and the column disappears from the table with them.
  RAW_SRCS    = amd_raw.cpp
  RAW_DEFS    = -DOBLIO_AMD_RAW
endif

.PHONY: all run2d run3d scale2d scale3d scale-mmd-2d scale-amd-2d scale-mmd-3d \
        scale-amd-3d profile clean help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: order_timing_cpp

# The hooked copy of the vendored AMD, for the AMDraw column. GENERATED from whatever
# private/Amd.cpp currently says, gitignored, removed by `make clean`, so no edited copy of vendored
# code can drift from the original unnoticed. The generator lives at the repo root rather than in
# experiments/ordering, which also uses it: ../README.md says nothing here may depend on a frozen
# study's contents, and copying it would leave two versions of a tool whose whole job is to fail
# loudly when the vendored source moves.
amd_raw.cpp: ../../private/Amd.cpp ../../tools/hook_amd.py
	@python3 ../../tools/hook_amd.py $< $@

order_timing_cpp: order_timing.cpp $(OBLIO_SRCS) $(OBLIO_HDRS) $(VENDOR_SRCS) $(RAW_SRCS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) $(VENDOR_DEFS) $(RAW_DEFS) $< $(OBLIO_SRCS) $(RAW_SRCS) $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

# TWO AXES, AND EVERY TARGET NAMES BOTH. `run` is what each method costs at four sizes; `scale` is
# how far ours sit from the vendored routine of the same lineage and whether that distance moves
# with n, over a ladder wide enough to tell a constant factor from a trend, one branch at a time
# because ten methods at 400 a side is a long wait for columns nobody is reading. `2d` and `3d` are
# the problem family.
#
# THE CUBIC HALF IS NOT A WIDER LADDER. Fill grows far faster, cliques are larger, mass elimination
# fires far more often, and the two families do not agree about how good our orderings are: our
# tie-break fills 6.5 percent below the vendored AMD on square grids and alternates sign on cubes.
# Separate runs and separate tables, never extra columns, because a row is a matrix and a column
# read down across two families is the error the second family exists to prevent.
#
# The targets were `run` and `scale` until 2026-08-09, when the cubic half arrived and left one
# shape named and the other implied. A default that is invisible in the name is how every figure
# recorded here came to be 2D without anyone deciding that it should be.
run2d: all
	./order_timing_cpp 2d

scale2d: scale-mmd-2d scale-amd-2d

scale-mmd-2d: all
	./order_timing_cpp mmd 2d

scale-amd-2d: all
	./order_timing_cpp amd 2d

run3d: all
	./order_timing_cpp 3d

scale3d: scale-mmd-3d scale-amd-3d

scale-mmd-3d: all
	./order_timing_cpp mmd 3d

scale-amd-3d: all
	./order_timing_cpp amd 3d

# Built separately, with symbols and without NDEBUG's effect on inlining left to chance: a profile
# has to describe the code that actually runs. See the header of order_profile.cpp for how to
# attach a profiler.
profile: order_profile_cpp

order_profile_cpp: order_profile.cpp $(OBLIO_SRCS) $(OBLIO_HDRS) $(VENDOR_SRCS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -g -fno-omit-frame-pointer -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) -g -fno-omit-frame-pointer $(VENDOR_DEFS) order_profile.cpp $(OBLIO_SRCS) \
	  $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

clean:
	rm -f .build-mode order_timing_cpp order_profile_cpp Amd.o Mmd.o amd_raw.cpp
