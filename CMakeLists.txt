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

# Detect whether the Fortran symbols have a trailing underscore (dpotrf_ vs
# dpotrf). The probe calls dpotrf_, which is a LAPACK symbol, so it must link
# LAPACK, not just BLAS: with a combined library (Accelerate) either works, but
# with split reference libs dpotrf_ lives in liblapack, and linking BLAS alone
# fails the probe and wrongly reports no underscore. BLAS and LAPACK share the
# convention, so this one result drives the macro for both. Underscore is the
# common case (Accelerate and the reference libs); the else branch is the rare
# no-underscore build.
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
check_c_source_compiles("
    extern void dpotrf_(char*,int*,double*,int*,int*);
    int main(){ char u='L'; int n=1,lda=1,info=0; double a=1.0;
                dpotrf_(&u,&n,&a,&lda,&info); return 0; }
" OBLIO_BLAS_UNDERSCORE)
unset(CMAKE_REQUIRED_LIBRARIES)

# ----- Library ---------------------------------------------------------------
# Source list matches REPO_MAP src/*.cpp. Amd.cpp and Mmd.cpp hold the ordering
# implementations that OrderEngine.cpp dispatches to — omitting them link-fails
# on mmd_order / amd_order.
add_library(oblio STATIC
    src/Amd.cpp
    src/BlasLapack.cpp
    src/ElmForestEngine.cpp
    src/Mmd.cpp
    src/MultiplyEngine.cpp
    src/NumFactorDynamic.cpp
    src/NumFactorEngine.cpp
    src/NumFactorStatic.cpp
    src/OrderEngine.cpp
    src/Permutation.cpp
    src/SolveEngine.cpp
    src/SparseMatrix.cpp
    src/SymFactorEngine.cpp
    src/Types.cpp
    src/UpdateBlock.cpp
    src/Vector.cpp
)

# Users #include the "oblio/..." headers (namespaced). PUBLIC so tests/examples
# that link oblio inherit the include path — no separate global include needed.
target_include_directories(oblio PUBLIC include)
target_link_libraries(oblio PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(OBLIO_BLAS_UNDERSCORE)
    target_compile_definitions(oblio PUBLIC OBLIO_BLAS_UNDERSCORE)
endif()

# ----- Tests -----------------------------------------------------------------
enable_testing()

# The suites in tests/. Each is a standalone executable linking oblio.
foreach(t smoke test_forest test_numfactor test_order test_permutation test_solve test_symfactor)
    add_executable(${t} tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE oblio)
    add_test(NAME ${t} COMMAND ${t})
endforeach()

# ----- Examples --------------------------------------------------------------
# No example targets yet. examples/basic.cc sketches the intended OblioEngine driver
# API, which is not ported, so it does not compile against the current library. Wire
# example executables up here once the driver and real examples exist.
