# Builds the Oblio units and every test in tests/. A simple alternative to the
# CMake build (CMakeLists.txt); both compile the same src/ and tests/.
#
#   make            build all test executables
#   make test       build and run them all
#   make <name>_cpp build one test (e.g. make test_order_cpp)
#   make clean
#
# -Iinclude points at include/, so #include "oblio/X.h" resolves to the project
# headers. Executables carry the _cpp suffix (coding convention) and are gitignored.

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -Iinclude

# BLAS and LAPACK, for the numeric factorization.
#
# OBLIO_BLAS_UNDERSCORE selects the trailing-underscore Fortran symbol convention
# (dpotrf_ rather than dpotrf), which is what Accelerate and the reference BLAS both
# use. See include/oblio/BlasLapack.h.
CXXFLAGS += -DOBLIO_BLAS_UNDERSCORE

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  BLAS_LIBS = -framework Accelerate
else
  BLAS_LIBS = -llapack -lblas
endif

# Oblio units, full warnings apply.
OBLIO_SRCS = \
  src/Types.cpp \
  src/SparseMatrix.cpp \
  src/Permutation.cpp \
  src/OrderEngine.cpp \
  src/ElmForestEngine.cpp \
  src/SymFactorEngine.cpp \
  src/BlasLapack.cpp \
  src/UpdateBlock.cpp \
  src/NumFactorStatic.cpp \
  src/NumFactorDynamic.cpp \
  src/NumFactorEngine.cpp \
  src/Vector.cpp \
  src/MultiplyEngine.cpp \
  src/SolveEngine.cpp

# Vendored ordering codes (SuiteSparse AMD, Sparspak MMD), copied verbatim, not
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

%_cpp: tests/%.cpp $(LIB_OBJS) $(LIB_HDRS) $(wildcard tests/*.h)
	$(CXX) $(CXXFLAGS) tests/$*.cpp $(LIB_OBJS) $(BLAS_LIBS) -o $@

test: all
	@for t in $(TEST_BINS); do echo "== $$t =="; ./$$t || exit 1; echo; done

clean:
	rm -f $(TEST_BINS) $(LIB_OBJS)
