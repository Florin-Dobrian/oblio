# pipeline benchmark - where a solve's time goes, across all its phases.
#
#   make        build the benchmark
#   make run    build and run it, default grid sizes
#   make clean
#
# The sibling folder, ordering/, measures one phase against itself. This one measures the phases
# against each other, which is what says whether optimizing the first is worth an afternoon.
#
# Not part of the main Oblio build, and unlike experiments/ this is not a frozen study: a
# benchmark runs against the current tree forever, so it links ../../src directly and must keep
# compiling as the code moves. Building it is cheap and running it is not, which is why `make`
# builds and only `make run` measures.
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
# The two vendored ordering codes are split back out for one reason only: they are not ours to
# clean up, so they compile with warnings off, exactly as the root Makefile does it.
VENDOR_SRCS = ../../src/Amd.cpp ../../src/Mmd.cpp
OBLIO_SRCS  = $(filter-out $(VENDOR_SRCS),$(wildcard ../../src/*.cpp))

.PHONY: all run clean

all: pipeline_timing_cpp

pipeline_timing_cpp: pipeline_timing.cpp $(OBLIO_SRCS) $(VENDOR_SRCS)
	$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS)
	$(CXX) $(CXXFLAGS) $< $(OBLIO_SRCS) Amd.o Mmd.o $(BLAS_LIBS) -o $@

run: all
	./pipeline_timing_cpp

clean:
	rm -f pipeline_timing_cpp Amd.o Mmd.o
