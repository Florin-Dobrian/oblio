# matrices benchmark - accuracy and performance of the pipeline on real matrices, from the
# SuiteSparse Matrix Collection. One driver per report, sharing the reader.
#
#:   make             print this list (the default goal; `make help` is the same thing)
#:   make all         build the four programs
#:   make accuracy    the accuracy pass over accuracy_candidates.txt
#:   make performance the timing pass over performance_candidates.txt
#:   make mmdorder    the ORDERING STEP ALONE, genmmd against MmdFlat, over all of data/. A
#:                    different question from the performance table's ordering column and so a
#:                    different set: every matrix whose PATTERN can be read, values or not,
#:                    positive definite or not, which is roughly twice as many. The two return
#:                    the same permutation on every matrix experiments/ordering's `make
#:                    mmdmatrices` checks, so the fill column is both of theirs and time is the
#:                    only difference between them. NAMED FOR ITS BRANCH because the amd half is
#:                    coming and will be `make amdorder` beside it. Not to be confused with
#:                    experiments/ordering's targets of the same names, which check that our
#:                    orderings match the vendored ones rather than timing them.
#:   make amdorder    the same for the amd branch, amd_order against AmdFlat, over the same set.
#:                    The CLOCK is on amd_order whole, which is what a caller pays and includes
#:                    AMD_aat and AMD_postorder, neither of which we do; every other column comes
#:                    from the hooked PRE-POSTORDER copy, since AmdFlat returns that order.
#:   make read        the reader alone over data/, an inventory of what is usable
#:   make clean
#:   make help        print this list
#:
#: Prefix any of these with OBLIO_PUBLIC=1 to build without the vendored orderings, as everyone
#: else does; the performance table's MMD and AMD rows then report a refusal.
#:
#: The matrices are NOT in the repository. Get them first, from this directory. These are the
#: commands that produced the two committed candidates files, one per report:
#:   ./ssget.py list --per-kind 8 --values --max-nnz 500000 > accuracy_candidates.txt
#:   ./ssget.py list --posdef --max-nnz 500000              > performance_candidates.txt
#:   ./ssget.py fetch accuracy_candidates.txt
#: They land in ../../data/<Group>/<name>.mtx, which is gitignored. `list` downloads nothing but
#: the collection's index, so it is cheap to re-run until the candidate set looks right, and
#: `fetch` skips what is already on disk.
#:
#: ./mmread_cpp ../../data/*/*.mtx                       what each file is, and what is refused
#: ./matrix_accuracy_cpp ../../data/HB/*.mtx             any subset, the shell choosing it
#: ./matrix_accuracy_cpp --max-fill=2e8 ../../data/*/*.mtx      raise the predicted-fill cap
#: ./matrix_performance_cpp --repeats=1 ../../data/HB/*.mtx     a quick timing pass
#:
#: ACCURACY.md and PERFORMANCE.md are the reports; README.md is the working notes behind them.
#
# Not part of the main Oblio build, and unlike experiments/ this is not a frozen study: a
# benchmark runs against the current tree forever, so it links ../../src directly and must keep
# compiling as the code moves. Building it is cheap and running it is not, which is why `make`
# builds and only the run targets measure.
#
# NDEBUG is on for the same reason it is in the sibling folders, so that what runs here is what a
# caller would run.
#
# THE VENDORED ORDERINGS, which matrix_performance names and matrix_accuracy does not. They are
# not ours to clean up, so they compile with warnings off exactly as the root Makefile does, and
# they are not published, so they live in ../../private and are optional. Absent, this folder still
# builds and the MMD and AMD rows of the performance table report a refusal.
#
# OBLIO_PUBLIC=1 builds as everyone else does, ignoring private/ even when it is there. Same word
# in every Makefile here that links the library, so one habit covers all of them.
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

ifdef OBLIO_PUBLIC
  VENDOR_SRCS =
else
  VENDOR_SRCS = $(wildcard ../../private/AmdVendored.cpp ../../private/MmdVendored.cpp \
                          ../../private/MmdCorrected.cpp)
endif

ifneq ($(VENDOR_SRCS),)
  VENDOR_DEFS = -DOBLIO_VENDORED_ORDERINGS
  # DERIVED FROM VENDOR_SRCS RATHER THAN SPELLED OUT. `-c` on several sources writes each
  # object to the basename in the CURRENT directory, so these must track the source names
  # exactly; they were spelled out until 2026-08-21 and a rename of the sources broke the
  # link, the compile step having quietly produced objects under the new names.
  VENDOR_OBJS = $(notdir $(VENDOR_SRCS:.cpp=.o))
  # The hooked copy, for `make amdorder`. It reports the elimination order BEFORE AMD_postorder,
  # which is the order AmdFlat returns, so the fill column can be both routines' as it is on the mmd
  # side. Present only where private/ is, exactly as the vendored objects are. GENERATED from
  # whatever private/AmdVendored.cpp currently says, gitignored, removed by clean, so no edited
  # copy of vendored code can drift from the original unnoticed. Same generator and same
  # discipline as ../ordering; see its Makefile for why the tool lives at the repo root.
  RAW_SRCS    = amd_raw.cpp
endif

# Which configuration this directory is built in. The two differ by a macro, so a binary from one
# is wrong for the other and no timestamp would say so. The $(shell) deletes the stamp when the
# recorded mode is not the current one, which rebuilds what depends on it; switching therefore
# rebuilds and repeating does not.
BUILD_MODE := $(if $(VENDOR_SRCS),vendored,public)
$(shell [ "`cat .build-mode 2>/dev/null`" = "$(BUILD_MODE)" ] || rm -f .build-mode)

.build-mode:
	@echo "$(BUILD_MODE)" > $@

# Everything in the pool, which is what `read` wants: the reader's job is to inventory what is on
# disk, not to serve a report. The two report targets take their candidates file instead, since
# data/ holds the union of both sets.
MATRICES = ../../data/*/*.mtx

# The run targets pass no options, so `run` uses the driver's own defaults: a predicted-fill
# cap of 50 million entries and no address-space limit on the child processes. Both are
# command-line flags; call the binary directly to change them.

# Bare `make` prints the target list rather than building. Stated explicitly because the default
# goal is otherwise the FIRST target in the file, so moving a rule would silently change it and the
# header above would stop being true with nothing to catch it. Every Makefile in this tree answers
# `make` the same way; see CODING_RULES.md.
.DEFAULT_GOAL := help

.PHONY: all read accuracy performance mmdorder amdorder clean help

# Print the target list from this file's header. Those lines are marked with #: so there is one
# source of truth: the header comment is the help text.
help:
	@grep '^#:' $(firstword $(MAKEFILE_LIST)) | cut -c3-

all: mmread_cpp matrix_accuracy_cpp matrix_performance_cpp matrix_ordering_cpp

mmread_cpp: mmread.cpp MatrixMarket.h $(OBLIO_SRCS) $(OBLIO_HDRS)
	$(CXX) $(CXXFLAGS) $< $(OBLIO_SRCS) $(BLAS_LIBS) -o $@

matrix_accuracy_cpp: matrix_accuracy.cpp MatrixMarket.h Matching.h $(OBLIO_SRCS) $(OBLIO_HDRS)
	$(CXX) $(CXXFLAGS) $< $(OBLIO_SRCS) $(BLAS_LIBS) -o $@

# The vendored sources compile separately with -w, exactly as ../ordering does it: they are not
# ours to clean up, and their warnings would bury ours.
matrix_performance_cpp: matrix_performance.cpp MatrixMarket.h $(OBLIO_SRCS) $(VENDOR_SRCS) \
                        $(OBLIO_HDRS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) $(VENDOR_DEFS) $< $(OBLIO_SRCS) $(VENDOR_OBJS) $(BLAS_LIBS) -o $@

amd_raw.cpp: ../../private/AmdVendored.cpp ../../tools/hook_amd.py
	@python3 ../../tools/hook_amd.py $< $@

matrix_ordering_cpp: matrix_ordering.cpp MatrixMarket.h $(OBLIO_SRCS) $(VENDOR_SRCS) \
                     $(RAW_SRCS) $(OBLIO_HDRS) .build-mode
	$(if $(VENDOR_SRCS),$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS))
	$(CXX) $(CXXFLAGS) $(VENDOR_DEFS) $< $(OBLIO_SRCS) $(RAW_SRCS) $(VENDOR_OBJS) \
	  $(BLAS_LIBS) -o $@

# Both run targets take the whole of data/ through the shell rather than through the program, so
# any subset is a matter of typing a different glob and needs no option here.
read: mmread_cpp
	@./mmread_cpp $(MATRICES)

# Over ALL of data/ rather than a candidates file, and that is the difference from the two report
# targets below. They each read their own list because the pool holds the union of two sets and
# neither report wants the other's matrices; an ordering has no such constraint, needing neither
# values nor definiteness, so what it wants is everything on disk.
mmdorder: matrix_ordering_cpp
	@./matrix_ordering_cpp mmd $(MATRICES)

# The amd half, over the same set. The clock is on `amd_order` whole and everything else comes
# from the hooked pre-postorder copy; see the driver's own note for why.
amdorder: matrix_ordering_cpp
	@./matrix_ordering_cpp amd $(MATRICES)

# The two report targets read their own candidates file rather than globbing data/, because the
# pool holds the union of both sets and neither report wants the other's matrices. A candidates
# line is `  Group/Name  n  nnz ...`, so the first field is the path minus the directory and the
# extension.
#
# THE FILTER MATCHES THE SLASH rather than excluding the comment character, and that is not a
# style choice: `#` starts a comment in a Makefile even inside $(shell ...), so an awk program
# containing /^#/ is truncated there and the call is left unterminated. Every data line has a
# Group/Name with a slash and every header and total line begins with `#`, so the slash separates
# them exactly.
ACCURACY_SET    = $(shell awk '$$1 ~ /\// {print "../../data/" $$1 ".mtx"}' accuracy_candidates.txt)
PERFORMANCE_SET = $(shell awk '$$1 ~ /\// {print "../../data/" $$1 ".mtx"}' performance_candidates.txt)

accuracy: matrix_accuracy_cpp
	@./matrix_accuracy_cpp $(ACCURACY_SET)

performance: matrix_performance_cpp
	@./matrix_performance_cpp $(PERFORMANCE_SET)

# The *.dSYM line removes macOS debug-symbol bundles, which are directories, hence -rf rather
# than -f. Nothing here emits one under the committed flags, and that is not the test: CXXFLAGS
# is overridable, `make CXXFLAGS="... -g"` produces one on macOS, and clean removes what a build
# in this directory CAN produce.
clean:
	rm -f amd_raw.cpp \
	      mmread_cpp matrix_accuracy_cpp matrix_performance_cpp matrix_ordering_cpp \
	      AmdVendored.o MmdVendored.o MmdCorrected.o .build-mode
	rm -rf *.dSYM
