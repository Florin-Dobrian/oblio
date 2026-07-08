# friend-access experiment — matvec via public API vs friend direct access.
# Reference/teaching only; not part of the main Oblio build.
#
#   make        build the test executable
#   make test   build and run it
#   make clean
#
# The executable carries the _cpp language suffix (coding convention) and is gitignored.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra

SRCS = Matrix.cpp Vector.cpp MultiplyEngine.cpp

.PHONY: all test clean

all: test_multiply_cpp

test_multiply_cpp: test_multiply.cpp $(SRCS) Matrix.h Vector.h MultiplyEngine.h
	$(CXX) $(CXXFLAGS) test_multiply.cpp $(SRCS) -o $@

test: all
	./test_multiply_cpp

clean:
	rm -f test_multiply_cpp
