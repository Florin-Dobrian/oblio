# storage-options experiment - one multiply function over two sparse storage formats.
# Reference/teaching only; not part of the main Oblio build.
#
#   make        build the test executable
#   make test   build and run it
#   make clean
#
# The executable carries the _cpp language suffix (coding convention) and is gitignored.
#
# NDEBUG is on deliberately: this experiment measures performance, so asserts are compiled
# out (release build). -O3 matters: the claim is that reaching a column through a pointer
# array costs nothing against hand-written CSC, and at -O0 that would measure nothing.
#
# No BLAS: the kernel is a sparse column scatter, not a dense block operation. (The
# friend-access experiment is the one that hands a raw block to BLAS.)

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra

SRCS = MultiplyEngine.cpp SparseMatrixStatic.cpp SparseMatrixDynamic.cpp
HDRS = SparseMatrixStatic.h SparseMatrixDynamic.h MultiplyEngine.h

.PHONY: all test clean

all: test_multiply_cpp

test_multiply_cpp: test_multiply.cpp $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) test_multiply.cpp $(SRCS) -o $@

test: all
	./test_multiply_cpp

clean:
	rm -f test_multiply_cpp
