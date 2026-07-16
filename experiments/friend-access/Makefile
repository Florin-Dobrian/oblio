# friend-access experiment - matvec via public API vs friend direct access vs BLAS.
# Reference/teaching only; not part of the main Oblio build.
#
#   make        build the test executable
#   make test   build and run it
#   make clean
#
# The executable carries the _cpp language suffix (coding convention) and is gitignored.
#
# NDEBUG is on deliberately: this experiment measures performance, so asserts are
# compiled out (release build). The gap it shows is then the structural cost of the
# API path (non-inlined cross-TU calls, no vectorization), not assertion overhead.
# (The template-instantiation experiment omits NDEBUG, its goal is different.)
#
# BLAS: Accelerate on macOS, the system BLAS on Linux. multiplyWithBlas calls the
# Fortran gemv (dgemv_/zgemv_), which both provide.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    BLAS_LIBS ?= -framework Accelerate
else
    BLAS_LIBS ?= -lblas
endif

SRCS = Matrix.cpp Vector.cpp MultiplyEngine.cpp

.PHONY: all test clean

all: test_multiply_cpp

test_multiply_cpp: test_multiply.cpp $(SRCS) Matrix.h Vector.h MultiplyEngine.h
	$(CXX) $(CXXFLAGS) test_multiply.cpp $(SRCS) $(BLAS_LIBS) -o $@

test: all
	./test_multiply_cpp

clean:
	rm -f test_multiply_cpp
