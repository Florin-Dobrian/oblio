# Spike build — Matrix/Vector/MultiplyEngine across the three compilation cases.
# Reference/teaching only; not part of the main oblio build.
#
#   make         build all three variants
#   make test    build and run all three (they must print identical results)
#   make clean

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

EXP_SRCS = Matrix_exp.cpp Vector_exp.cpp MultiplyEngine_exp.cpp
EXT_SRCS = Matrix_ext.cpp Vector_ext.cpp MultiplyEngine_ext.cpp

.PHONY: all test clean

all: test_multiply_tpl test_multiply_exp test_multiply_ext

# Case 1 — header-only, nothing to link
test_multiply_tpl: test_multiply.cpp Matrix_tpl.h Vector_tpl.h MultiplyEngine_tpl.h
	$(CXX) $(CXXFLAGS) -DOBLIO_SPIKE_TPL test_multiply.cpp -o $@

# Case 2 — explicit instantiation, forcing only; must link the three _exp.cpp
test_multiply_exp: test_multiply.cpp $(EXP_SRCS) \
                   Matrix_exp.h Vector_exp.h MultiplyEngine_exp.h
	$(CXX) $(CXXFLAGS) -DOBLIO_SPIKE_EXP test_multiply.cpp $(EXP_SRCS) -o $@

# Case 3 — explicit instantiation + extern template; must link the three _ext.cpp
test_multiply_ext: test_multiply.cpp $(EXT_SRCS) \
                   Matrix_ext.h Vector_ext.h MultiplyEngine_ext.h
	$(CXX) $(CXXFLAGS) -DOBLIO_SPIKE_EXT test_multiply.cpp $(EXT_SRCS) -o $@

test: all
	./test_multiply_tpl
	./test_multiply_exp
	./test_multiply_ext

clean:
	rm -f test_multiply_tpl test_multiply_exp test_multiply_ext
