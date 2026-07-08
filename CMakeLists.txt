cmake_minimum_required(VERSION 3.14)
project(oblio VERSION 0.1.0 LANGUAGES C CXX)
# C is enabled (alongside CXX) because the BLAS underscore probe below uses
# check_c_source_compiles — a bare extern avoids C++ name-mangling of dpotrf_.
# NOTE: VERSION is a placeholder — reconcile with the version-identity decision
# (oblio 11 vs a fresh 1.0) recorded in docs/DESIGN_DECISIONS.md before first release.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----- BLAS / LAPACK ---------------------------------------------------------
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Detect whether BLAS symbols have trailing underscore (dpotrf_ vs dpotrf).
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
check_c_source_compiles("
    extern void dpotrf_(char*,int*,double*,int*,int*);
    int main(){ char u='L'; int n=1,lda=1,info=0; double a=1.0;
                dpotrf_(&u,&n,&a,&lda,&info); return 0; }
" OBLIO_BLAS_UNDERSCORE)

# ----- Library ---------------------------------------------------------------
# Source list matches REPO_MAP src/*.cc. Amd.cc and Mmd.cc hold the ordering
# implementations that OrderEngine.cc dispatches to — omitting them link-fails
# on mmd_order / amd_order.
add_library(oblio STATIC
    src/OrderEngine.cc
    src/Mmd.cc
    src/Amd.cc
    src/Matrix.cc
    src/Vector.cc
    src/DenseMatrix.cc
    src/Factors.cc
    src/FactorEngine.cc
    src/SolveEngine.cc
    src/OblioEngine.cc
)

# Users #include "oblio/OblioEngine.h" (namespaced). PUBLIC so tests/examples
# that link oblio inherit the include path — no separate global include needed.
target_include_directories(oblio PUBLIC include)
target_link_libraries(oblio PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(OBLIO_BLAS_UNDERSCORE)
    target_compile_definitions(oblio PUBLIC OBLIO_BLAS_UNDERSCORE)
endif()

# ----- Tests -----------------------------------------------------------------
enable_testing()

# The four suites in tests/ (real/complex × smoke/extended). Each is a standalone
# executable linking oblio.
foreach(t test_smoke_real test_smoke_complex test_extended_real test_extended_complex)
    add_executable(${t} tests/${t}.cc)
    target_link_libraries(${t} PRIVATE oblio)
    add_test(NAME ${t} COMMAND ${t})
endforeach()

# ----- Examples --------------------------------------------------------------
add_executable(example_basic examples/basic.cc)
target_link_libraries(example_basic PRIVATE oblio)
