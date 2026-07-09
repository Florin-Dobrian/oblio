# Parallel port tree — builds the new units and every test in tests/.
# The PoC build (top-level CMake) is untouched; this compiles port/ in isolation.
#
#   make            build all test executables
#   make test       build and run them all
#   make <name>_cpp build one test (e.g. make test_order_cpp)
#   make clean
#
# -Iinclude points at port/include, so #include "oblio/X.h" resolves to the NEW
# headers, never the PoC ones one level up. Executables carry the _cpp suffix
# (coding convention) and are gitignored.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -Iinclude

# Oblio units — full warnings apply.
OBLIO_SRCS = \
  src/SparseMatrix.cpp \
  src/Permutation.cpp \
  src/OrderEngine.cpp

# Vendored ordering codes (SuiteSparse AMD, Sparspak MMD) — copied verbatim, not
# maintained here, so warnings are silenced (-w) for these two files only.
VENDOR_SRCS = \
  src/Amd.cpp \
  src/Mmd.cpp

LIB_SRCS = $(OBLIO_SRCS) $(VENDOR_SRCS)
LIB_HDRS = $(wildcard include/oblio/*.h)

# Compile each source to an object so per-file flags can differ. Oblio objects get
# full warnings; vendored objects get -w.
OBLIO_OBJS  = $(OBLIO_SRCS:.cpp=.o)
VENDOR_OBJS = $(VENDOR_SRCS:.cpp=.o)
LIB_OBJS    = $(OBLIO_OBJS) $(VENDOR_OBJS)

# One executable per tests/*.cpp file, named <stem>_cpp.
TEST_SRCS = $(wildcard tests/*.cpp)
TEST_BINS = $(patsubst tests/%.cpp,%_cpp,$(TEST_SRCS))

.PHONY: all test clean

all: $(TEST_BINS)

$(OBLIO_OBJS): %.o: %.cpp $(LIB_HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(VENDOR_OBJS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -w -c $< -o $@

%_cpp: tests/%.cpp $(LIB_OBJS) $(LIB_HDRS)
	$(CXX) $(CXXFLAGS) tests/$*.cpp $(LIB_OBJS) -o $@

test: all
	@for t in $(TEST_BINS); do echo "== $$t =="; ./$$t || exit 1; echo; done

clean:
	rm -f $(TEST_BINS) $(LIB_OBJS)
