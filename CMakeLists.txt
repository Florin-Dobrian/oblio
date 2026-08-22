cmake_minimum_required(VERSION 3.14)
project(oblio VERSION 0.1.0 LANGUAGES C CXX)
# C is enabled (alongside CXX) because the BLAS underscore probe below uses
# check_c_source_compiles, a bare extern avoids C++ name-mangling of dpotrf_.
# NOTE: VERSION is a placeholder, reconcile with the version-identity decision
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
# Source list matches src/*.cpp. The two vendored ordering codes are not
# here: they live in private/, which is not part of the published tree, and are
# added below only when that directory exists.
#
# THIS LIST HAD DRIFTED BADLY AND WAS REPAIRED 2026-08-17. It named src/Amd1B.cpp and
# src/Amd2B.cpp, which were deleted on 2026-08-16, and omitted src/MmdFlat.cpp, which is the DEFAULT
# ordering, along with MmdChained, AmdCompacted and MmdCompacted. So this build had been failing to
# configure for a day and nobody noticed, the Makefile being what the development loop runs.
#
# NOTHING DETECTS THE DRIFT, because the two builds share no source list and neither reads the
# other. Adding or removing a file in src/ means editing BOTH this list and OBLIO_SRCS in the
# Makefile. The Makefile's own comment beside that list says the same thing from the other side.
# CMake is the release path, so running it occasionally is the only thing that catches this; see
# README.md, "Which to use".
add_library(oblio STATIC
    src/AmdFlat.cpp
    src/AmdCompacted.cpp
    src/BlasLapack.cpp
    src/DirectSolver.cpp
    src/ElmForestEngine.cpp
    src/MmdFlat.cpp
    src/MmdChained.cpp
    src/MmdCompacted.cpp
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
    src/UpdateMatrix.cpp
    src/Vector.cpp
)

# Users #include the "oblio/..." headers (namespaced). PUBLIC so tests/examples
# that link oblio inherit the include path, so no separate global include is needed.
# The vendored orderings, optional. Present, Ordering::MmdVendored and Ordering::AmdVendored
# work; absent, the library still builds and those two enumerators refuse. Detected rather than
# configured, so a
# tree that has the directory needs no flag and a clone needs no edit. Warnings are off for these
# two files only, since they are copied verbatim and not maintained here.
file(GLOB OBLIO_VENDOR_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/private/AmdVendored.cpp
                            ${CMAKE_CURRENT_SOURCE_DIR}/private/MmdVendored.cpp)
if(OBLIO_VENDOR_SRCS)
    target_sources(oblio PRIVATE ${OBLIO_VENDOR_SRCS})
    set_source_files_properties(${OBLIO_VENDOR_SRCS} PROPERTIES COMPILE_OPTIONS "-w")
    target_compile_definitions(oblio PUBLIC OBLIO_VENDORED_ORDERINGS)
    message(STATUS "Vendored orderings: found in private/, MMD and AMD enabled")
else()
    message(STATUS "Vendored orderings: private/ absent, MMD and AMD will refuse")
endif()

target_include_directories(oblio PUBLIC include)
target_link_libraries(oblio PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(OBLIO_BLAS_UNDERSCORE)
    target_compile_definitions(oblio PUBLIC OBLIO_BLAS_UNDERSCORE)
endif()

# ----- Tests -----------------------------------------------------------------
enable_testing()

# The suites in tests/. Each is a standalone executable linking oblio.
#
# Globbed rather than listed, deliberately, and against the usual CMake advice.
# The Makefile discovers them with $(wildcard tests/*.cpp), and an enumerated
# list here has to be kept in step with it by hand, which is a step nothing
# detects: both builds succeed, both report passing, and the shorter list is
# simply a smaller claim. That is not hypothetical, test_pipeline was added and
# not listed, so ctest ran seven suites of eight and omitted the largest.
# CONFIGURE_DEPENDS makes the build re-glob when the directory changes, so a new
# tests/*.cpp file needs no edit here, exactly as it needs none in the Makefile.
file(GLOB TEST_SRCS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
foreach(src ${TEST_SRCS})
    get_filename_component(t ${src} NAME_WE)
    add_executable(${t} ${src})
    target_link_libraries(${t} PRIVATE oblio)
    add_test(NAME ${t} COMMAND ${t})
endforeach()

# ----- Examples --------------------------------------------------------------
# Each builds a standalone driver linking oblio. Globbed for the same reason as
# the tests. All three compile against the current library: example_basic.cpp
# drives DirectSolver, and the two example_pipeline_* files wire the engines by
# hand, one per scalar type. The sources carry the example_ prefix, as tests/
# carry test_, so the target names come out prefixed here with no edit. Note that
# none of them is run by anything, here or in the Makefile, which is its own item
# in docs/TODO.md.
file(GLOB EXAMPLE_SRCS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.cpp)
foreach(src ${EXAMPLE_SRCS})
    get_filename_component(e ${src} NAME_WE)
    add_executable(${e} ${src})
    target_link_libraries(${e} PRIVATE oblio)
endforeach()
