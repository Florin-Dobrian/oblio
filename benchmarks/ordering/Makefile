# ordering benchmark - what each ordering method costs, in time and in fill.
#
#:   make        build the benchmark
#:   make run    build and run it, default grid sizes
#:   make clean
#:   make help   print this list
#:
#: Prefix any of these with OBLIO_PUBLIC=1 to build as everyone else does, without the vendored
#: orderings; see docs/DESIGN_DECISIONS.md. The same word works in the repo root, both benchmarks
#: and experiments/ordering.
#
# Not part of the main Oblio build, and unlike experiments/ this is not a frozen study: a
# benchmark runs against the current tree forever, so it links ../../src directly and must keep
# compiling as the code moves. Building it is cheap and running it is not, which is why `make`
# builds and only `make run` measures.
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

.PHONY: all run profile clean help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: order_timing_cpp

order_timing_cpp: order_timing.cpp $(OBLIO_SRCS) $(VENDOR_SRCS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) $(VENDOR_DEFS) $< $(OBLIO_SRCS) $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

run: all
	./order_timing_cpp

# Built separately, with symbols and without NDEBUG's effect on inlining left to chance: a profile
# has to describe the code that actually runs. See the header of order_profile.cpp for how to
# attach a profiler.
profile: order_profile_cpp

order_profile_cpp: order_profile.cpp $(OBLIO_SRCS) $(VENDOR_SRCS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -g -fno-omit-frame-pointer -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) -g -fno-omit-frame-pointer $(VENDOR_DEFS) order_profile.cpp $(OBLIO_SRCS) \
	  $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

clean:
	rm -f .build-mode order_timing_cpp order_profile_cpp Amd.o Mmd.o
