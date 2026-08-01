# ordering benchmark - what each ordering method costs, in time and in fill.
#
#   make        build the benchmark
#   make run    build and run it, default grid sizes
#   make clean
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
# The two vendored ordering codes are split back out for one reason only: they are not ours to
# clean up, so they compile with warnings off, exactly as the root Makefile does it.
VENDOR_SRCS = ../../src/Amd.cpp ../../src/Mmd.cpp
OBLIO_SRCS  = $(filter-out $(VENDOR_SRCS),$(wildcard ../../src/*.cpp))

.PHONY: all run profile clean

all: order_timing_cpp

order_timing_cpp: order_timing.cpp $(OBLIO_SRCS) $(VENDOR_SRCS)
	$(CXX) $(CXXFLAGS) -w -c $(VENDOR_SRCS)
	$(CXX) $(CXXFLAGS) $< $(OBLIO_SRCS) Amd.o Mmd.o $(BLAS_LIBS) -o $@

run: all
	./order_timing_cpp

# Built separately, with symbols and without NDEBUG's effect on inlining left to chance: a profile
# has to describe the code that actually runs. See the header of order_profile.cpp for how to
# attach a profiler.
profile: order_profile_cpp

order_profile_cpp: order_profile.cpp $(OBLIO_SRCS) $(VENDOR_SRCS)
	$(CXX) $(CXXFLAGS) -g -fno-omit-frame-pointer -w -c $(VENDOR_SRCS)
	$(CXX) $(CXXFLAGS) -g -fno-omit-frame-pointer order_profile.cpp $(OBLIO_SRCS) \
	  Amd.o Mmd.o $(BLAS_LIBS) -o $@

clean:
	rm -f order_timing_cpp order_profile_cpp Amd.o Mmd.o
