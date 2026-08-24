# ordering benchmark - what each ordering method costs, in time and in fill.
#
#:   make          print this list (the default goal; `make help` is the same thing)
#:   make all      build the benchmark
#:   make run2d    build and run it, the default square-grid sizes, every method
#:   make run3d    the same on cubic grids, which are a DIFFERENT PROBLEM FAMILY
#:   make scale2d  a wider ladder, one branch at a time, with the gap columns
#:   make scale3d  the same, cubic
#:                 Also scale-mmd-2d, scale-amd-2d, scale-mmd-3d, scale-amd-3d for one branch
#:   make digest-record   record a permutation digest for every driver over 73 small grids
#:   make digest          compare against it. "Did anything move?", asked of all nine drivers
#:                 at once, in half a second. For re-encodings, which must leave every
#:                 permutation untouched. NOT a substitute for `make amdorder` and
#:                 `make mmdorder`, which compare against the vendored routines themselves;
#:                 this only says "same as last time". The baseline is gitignored.
#:   make phases2d what the vendored AMD's own time is MADE OF, per size: AMD_aat and
#:                 AMD_postorder are work Oblio does not do, so every ratio in the tables
#:                 above is against a denominator larger than the comparable region. This
#:                 measures it rather than estimating it, and says whether the share is flat
#:                 in n. Needs private/.
#:   make phases3d the same, cubic
#:   make profile  build order_profile_cpp, the binary the invocations below run: one method in
#:                 a loop, with symbols, for attaching a sampling profiler
#:   make clean
#:   make help     print this list
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
# that asks for Ordering::MmdVendored or Ordering::AmdVendored reports a refusal.
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
  VENDOR_SRCS = $(wildcard ../../private/AmdVendored.cpp ../../private/MmdVendored.cpp \
                          ../../private/MmdCorrected.cpp)
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
  # DERIVED FROM VENDOR_SRCS RATHER THAN SPELLED OUT. `-c` on several sources writes each
  # object to the basename in the CURRENT directory, so these must track the source names
  # exactly; they were spelled out until 2026-08-21 and a rename of the sources broke the
  # link, the compile step having quietly produced objects under the new names.
  VENDOR_OBJS = $(notdir $(VENDOR_SRCS:.cpp=.o))
  # The hooked copy, for the AMDraw column. Present only where private/ is, exactly as the vendored
  # objects are, and the column disappears from the table with them.
  RAW_SRCS    = amd_raw.cpp
  RAW_DEFS    = -DOBLIO_AMD_RAW
  # The second hooked copy, for the phases table. A DIFFERENT INSTRUMENT rather than a second use
  # of the first: the raw copy carries a vector per vertex, so it is an oracle and cannot be timed,
  # and this one carries five timestamps, so it can be timed and reports no order. Both are
  # generated from the same source by the same tool, and both are gitignored and removed by clean.
  TIMED_SRCS  = amd_timed.cpp
  TIMED_DEFS  = -DOBLIO_AMD_TIMED
endif

.PHONY: all run2d run3d scale2d scale3d scale-mmd-2d scale-amd-2d scale-mmd-3d \
        scale-amd-3d phases2d phases3d digest digest-record profile clean help

# Bare `make` prints the target list rather than building. Stated explicitly because the default
# goal is otherwise the FIRST target in the file, so moving a rule would silently change it and the
# header above would stop being true with nothing to catch it. Every Makefile in this tree answers
# `make` the same way; see docs/CODING_RULES.md.
.DEFAULT_GOAL := help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: order_timing_cpp

# The hooked copy of the vendored AMD, for the AMDraw column. GENERATED from whatever
# private/AmdVendored.cpp currently says, gitignored, removed by `make clean`, so no edited copy
# of vendored code can drift from the original unnoticed. The generator lives at the repo root
# rather than in
# experiments/ordering, which also uses it: ../README.md says nothing here may depend on a frozen
# study's contents, and copying it would leave two versions of a tool whose whole job is to fail
# loudly when the vendored source moves.
amd_raw.cpp: ../../private/AmdVendored.cpp ../../tools/hook_amd.py
	@python3 ../../tools/hook_amd.py $< $@

# The timed copy, for the phases table. Same generator and the same discipline; what differs is
# what it is for. Its own header sets out the five phases and why the postorder is TIMED rather
# than removed: AMD_postorder runs inside AMD_2 and the output permutation is built from the ranks
# it writes, so deleting the call would leave the returned permutation meaningless and take both of
# the checks with it.
amd_timed.cpp: ../../private/AmdVendored.cpp ../../tools/hook_amd.py
	@python3 ../../tools/hook_amd.py $< $@ --mode=timed

order_timing_cpp: order_timing.cpp $(OBLIO_SRCS) $(OBLIO_HDRS) $(VENDOR_SRCS) $(RAW_SRCS) \
                  $(TIMED_SRCS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) $(VENDOR_DEFS) $(RAW_DEFS) $(TIMED_DEFS) $< $(OBLIO_SRCS) $(RAW_SRCS) \
	  $(TIMED_SRCS) $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

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

# A THIRD QUESTION AND SO A THIRD TABLE. `run` and `scale` both ask what one of our orderings costs
# against the vendored routine of its lineage; this asks what that BASELINE is made of, which is a
# question about the denominator every one of those ratios divides by. Two of the vendored routine's
# phases have no counterpart here, `AMD_aat` because our pattern arrives full-symmetric and
# `AMD_postorder` because `ElmForestEngine` postorders the exact supernodal tree later, so a figure
# like `AmdFlat 1.82x` is measured against more work than the comparison is about. This folder has
# carried one estimate of the correction since 2026-08-01, 8 percent from a single profile at 140 a
# side; this measures it per size, which is the part that matters, since a correction that FALLS
# with n would mean part of the amd branch's 2D growth is the baseline shrinking rather than us.
#
# Amd only, and the mode word says so on its own: `phases` implies the amd branch, there being no
# mmd counterpart. genmmd has neither an aat nor a postorder, so `mmd_order`'s time is already the
# comparable region and MmdFlat's ratios need no correction at all.
phases2d: all
	./order_timing_cpp phases 2d

phases3d: all
	./order_timing_cpp phases 3d

# Built separately, with symbols and without NDEBUG's effect on inlining left to chance: a profile
# has to describe the code that actually runs. See the header of order_profile.cpp for how to
# attach a profiler.
profile: order_profile_cpp

# THE DIGEST HARNESS. It links only our own drivers, so it needs neither the vendored sources nor
# the hooked copies and builds in a published tree exactly as it does here. See its header for what
# it is for and, more importantly, for what it is not.
order_digest_cpp: order_digest.cpp $(OBLIO_SRCS) $(OBLIO_HDRS)
	$(CXX) $(CXXFLAGS) $< $(OBLIO_SRCS) $(BLAS_LIBS) -o $@

# RECORD THEN COMPARE. Record at the start of a session, compare after each change. A baseline that
# is regenerated to make a failure go away is worse than none, so this one is gitignored and
# short-lived by design: the durable external check is `make amdorder` and `make mmdorder`, which
# cannot be quietly re-blessed.
digest-record: order_digest_cpp
	./order_digest_cpp record

digest: order_digest_cpp
	./order_digest_cpp

order_profile_cpp: order_profile.cpp $(OBLIO_SRCS) $(OBLIO_HDRS) $(VENDOR_SRCS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -g -fno-omit-frame-pointer -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) -g -fno-omit-frame-pointer $(VENDOR_DEFS) order_profile.cpp $(OBLIO_SRCS) \
	  $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

# `make profile` links with -g, and on macOS Apple clang then runs dsymutil, which drops a .dSYM
# debug bundle beside the executable. It is a DIRECTORY, so `rm -f` cannot remove it and
# .gitignore's *_cpp does not match it either (the name does not end there); both are handled, here
# with -rf and there with an explicit *.dSYM/ entry. Every clean rule in the tree carries the -rf
# line, because CXXFLAGS is overridable everywhere and any directory produces a bundle under
# `make CXXFLAGS="... -g"`, which is what debugging looks like.
#
# The comment lives ABOVE the rule and not inside it: make echoes every recipe line before running
# it, and a # there is a shell comment rather than a make one, so an explanation written in the
# recipe gets printed in full every time someone types `make clean`.
clean:
	rm -f .build-mode order_timing_cpp order_profile_cpp order_digest_cpp .digest-baseline
	rm -f AmdVendored.o MmdVendored.o MmdCorrected.o amd_raw.cpp amd_timed.cpp
	rm -rf *.dSYM
