# Template-instantiation example — Matrix/Vector/MultiplyEngine across the three variants.
# Reference/teaching only; not part of the main Oblio build.
#
#   make         build all three variants
#   make test    build and run all three (they must print identical results)
#   make clean
#
# Executables carry the _cpp language suffix (coding convention), so *_cpp is
# gitignored and never committed.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra

# .cpp sources to link per explicit variant. (Implicit is header-only — it links
# nothing but the test itself, so it has no source list.)
PLAIN_EXPLICIT_SRCS   = MatrixPlainExplicit.cpp VectorPlainExplicit.cpp MultiplyEnginePlainExplicit.cpp
GUARDED_EXPLICIT_SRCS = MatrixGuardedExplicit.cpp VectorGuardedExplicit.cpp MultiplyEngineGuardedExplicit.cpp

.PHONY: all test clean

all: test_multiply_implicit_cpp test_multiply_plain_explicit_cpp test_multiply_guarded_explicit_cpp

# Implicit — body in header; nothing to link but the test
test_multiply_implicit_cpp: test_multiply.cpp \
                            MatrixImplicit.h VectorImplicit.h MultiplyEngineImplicit.h
	$(CXX) $(CXXFLAGS) -DOBLIO_TI_IMPLICIT test_multiply.cpp -o $@

# Plain explicit — bodies in .cpp, header signatures only
test_multiply_plain_explicit_cpp: test_multiply.cpp $(PLAIN_EXPLICIT_SRCS) \
                                  MatrixPlainExplicit.h VectorPlainExplicit.h MultiplyEnginePlainExplicit.h
	$(CXX) $(CXXFLAGS) -DOBLIO_TI_PLAIN_EXPLICIT test_multiply.cpp $(PLAIN_EXPLICIT_SRCS) -o $@

# Guarded explicit — plain explicit + extern template guard
test_multiply_guarded_explicit_cpp: test_multiply.cpp $(GUARDED_EXPLICIT_SRCS) \
                                    MatrixGuardedExplicit.h VectorGuardedExplicit.h MultiplyEngineGuardedExplicit.h
	$(CXX) $(CXXFLAGS) -DOBLIO_TI_GUARDED_EXPLICIT test_multiply.cpp $(GUARDED_EXPLICIT_SRCS) -o $@

test: all
	./test_multiply_implicit_cpp
	./test_multiply_plain_explicit_cpp
	./test_multiply_guarded_explicit_cpp

clean:
	rm -f test_multiply_implicit_cpp test_multiply_plain_explicit_cpp test_multiply_guarded_explicit_cpp
