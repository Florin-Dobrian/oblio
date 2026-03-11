cmake_minimum_required(VERSION 3.14)
project(oblio VERSION 0.1.0 LANGUAGES CXX)

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

# ----- Include path ----------------------------------------------------------
# Users #include "oblio/OblioEngine.h" (namespaced)
include_directories(include)

# ----- Library ---------------------------------------------------------------
add_library(oblio STATIC
    src/OrderEngine.cc
    src/Matrix.cc
    src/Vector.cc
    src/DenseMatrix.cc
    src/Factors.cc
    src/FactorEngine.cc
    src/SolveEngine.cc
    src/OblioEngine.cc
)

target_include_directories(oblio PUBLIC include)
target_link_libraries(oblio PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(OBLIO_BLAS_UNDERSCORE)
    target_compile_definitions(oblio PUBLIC OBLIO_BLAS_UNDERSCORE)
endif()

# ----- Tests -----------------------------------------------------------------
enable_testing()

add_executable(test_smoke tests/test_smoke.cc)
target_link_libraries(test_smoke PRIVATE oblio)
add_test(NAME smoke COMMAND test_smoke)

# ----- Examples --------------------------------------------------------------
add_executable(example_basic examples/basic.cc)
target_link_libraries(example_basic PRIVATE oblio)
