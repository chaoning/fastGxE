cmake_minimum_required(VERSION 3.16)
project(fastgxeProject)
enable_testing()

option(FASTGXE_FULL_STATIC_LINK "Link executables as fully static binaries when static runtimes are available" ON)

set(CMAKE_CXX_COMPILER "/home/chaoning/intel/oneapi/compiler/2025.2/bin/icpx")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qopenmp")

# Set external library paths
set(PROJECTROOT ${CMAKE_SOURCE_DIR})
set(COMPILERROOT /home/chaoning/intel/oneapi/compiler/2025.2/)
set(MKLROOT /home/chaoning/intel/oneapi/mkl/2025.2/)
set(GSLROOT ${PROJECTROOT}/external/gsl-2.7/)
set(EIGENROOT ${PROJECTROOT}/external/eigen-5.0.0/)
set(LBFGSppROOT ${PROJECTROOT}/external/LBFGSpp/include/)
set(BOOSTROOT ${PROJECTROOT}/external/boost_1_90_0/)
set(UTILROOT ${PROJECTROOT}/utils/)
set(GMATRIXROOT ${PROJECTROOT}/gmatrix/)
set(SPDROOT ${PROJECTROOT}/external/)
set(CLIROOT ${PROJECTROOT}/external/)

if(FASTGXE_FULL_STATIC_LINK)
    set(CMAKE_SKIP_RPATH TRUE)
    add_link_options(-static -static-libstdc++ -static-libgcc)
    set(FASTGXE_OPENMP_LIB ${COMPILERROOT}/lib/libiomp5.a)
else()
    set(CMAKE_BUILD_RPATH "${COMPILERROOT}/lib")
    set(CMAKE_INSTALL_RPATH "${COMPILERROOT}/lib")
    set(FASTGXE_OPENMP_LIB -liomp5)
endif()

# Include directories
include_directories(
    ${MKLROOT}/include
    ${GSLROOT}/include
    ${EIGENROOT}
    ${EIGENROOT}/unsupported
    ${BOOSTROOT}
    ${UTILROOT}
    ${GMATRIXROOT}
    ${LBFGSppROOT}
    ${SPDROOT}
    ${CLIROOT}
)

# Source files
set(SOURCE_FILES
    utils/chi2score.cpp
    utils/acat_utils.cpp
    utils/EigenMatrix_utils.cpp
    utils/geno.cpp
    utils/iterator_utils.cpp
    utils/partitionLowerTriangle.cpp
    utils/phen.cpp
    utils/random_sampling.cpp
    utils/someCommonFun.cpp
    utils/string_utils.cpp
    gmatrix/gmatrix.cpp
    gmatrix/processGRM.cpp
    fastgxe/fastgxe.cpp
    fastgxe/mom.cpp
    fastgxe/mmsusie.cpp
    fastgxe/main.cpp
)

set(FASTGXE_MKL_LIBS
    -Wl,--start-group
    ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
    ${MKLROOT}/lib/intel64/libmkl_intel_thread.a
    ${MKLROOT}/lib/intel64/libmkl_core.a
    -Wl,--end-group
)

set(FASTGXE_GSL_LIBS
    ${GSLROOT}/lib/libgsl.a
    ${GSLROOT}/lib/libgslcblas.a
)

set(FASTGXE_RUNTIME_LIBS
    ${FASTGXE_OPENMP_LIB}
    -lpthread
    -lm
    -ldl
    -lrt
)

# Define the executable
add_executable(fastgxe ${SOURCE_FILES})

# Link Libraries
target_link_libraries(fastgxe
    ${FASTGXE_MKL_LIBS}
    ${FASTGXE_GSL_LIBS}
    ${FASTGXE_RUNTIME_LIBS}
)

add_executable(chi2score_test
    test/chi2score_test.cpp
    utils/chi2score.cpp
)

target_link_libraries(chi2score_test
    ${FASTGXE_GSL_LIBS}
    ${FASTGXE_RUNTIME_LIBS}
)

add_test(NAME chi2score_test COMMAND chi2score_test)

add_executable(eigenmatrix_utils_example
    test/eigenmatrix_utils_example.cpp
    utils/EigenMatrix_utils.cpp
)

target_link_libraries(eigenmatrix_utils_example
    ${FASTGXE_MKL_LIBS}
    ${FASTGXE_RUNTIME_LIBS}
)

add_test(NAME eigenmatrix_utils_example COMMAND eigenmatrix_utils_example)

add_executable(iterator_utils_test
    test/iterator_utils_test.cpp
    utils/iterator_utils.cpp
)

add_test(NAME iterator_utils_test COMMAND iterator_utils_test)

add_executable(string_utils_test
    test/string_utils_test.cpp
    utils/string_utils.cpp
)

add_test(NAME string_utils_test COMMAND string_utils_test)

add_executable(acat_utils_test
    test/acat_utils_test.cpp
    utils/acat_utils.cpp
)

add_test(NAME acat_utils_test COMMAND acat_utils_test)

add_executable(random_sampling_example
    test/random_sampling_example.cpp
    utils/random_sampling.cpp
)

add_test(NAME random_sampling_example COMMAND random_sampling_example)

add_executable(partition_lower_triangle_example
    test/partition_lower_triangle_example.cpp
    utils/partitionLowerTriangle.cpp
)

add_executable(geno_test
    test/geno_test.cpp
    utils/geno.cpp
    utils/string_utils.cpp
)

target_link_libraries(geno_test
    ${FASTGXE_MKL_LIBS}
    ${FASTGXE_RUNTIME_LIBS}
)

add_test(
    NAME geno_test
    COMMAND ${CMAKE_COMMAND} -E env
        OMP_NUM_THREADS=1
        KMP_AFFINITY=disabled
        $<TARGET_FILE:geno_test>
)

add_executable(phen_workflow_example
    test/phen_workflow_example.cpp
    utils/phen.cpp
    utils/iterator_utils.cpp
    utils/string_utils.cpp
    utils/EigenMatrix_utils.cpp
)

target_link_libraries(phen_workflow_example
    ${FASTGXE_MKL_LIBS}
    ${FASTGXE_RUNTIME_LIBS}
)

add_test(
    NAME phen_workflow_example
    COMMAND ${CMAKE_COMMAND} -E env
        OMP_NUM_THREADS=1
        KMP_AFFINITY=disabled
        $<TARGET_FILE:phen_workflow_example>
)

add_executable(grm_binary_roundtrip_test
    test/grm_binary_roundtrip_test.cpp
    gmatrix/gmatrix.cpp
    gmatrix/processGRM.cpp
    utils/geno.cpp
    utils/string_utils.cpp
    utils/EigenMatrix_utils.cpp
    utils/partitionLowerTriangle.cpp
    utils/iterator_utils.cpp
)

target_link_libraries(grm_binary_roundtrip_test
    ${FASTGXE_MKL_LIBS}
    ${FASTGXE_RUNTIME_LIBS}
)

add_test(
    NAME grm_binary_roundtrip_test
    COMMAND ${CMAKE_COMMAND} -E env
        OMP_NUM_THREADS=1
        KMP_AFFINITY=disabled
        $<TARGET_FILE:grm_binary_roundtrip_test>
)

add_executable(gmatrix_run_smoke_test
    test/gmatrix_run_smoke_test.cpp
    gmatrix/gmatrix.cpp
    gmatrix/processGRM.cpp
    utils/geno.cpp
    utils/string_utils.cpp
    utils/EigenMatrix_utils.cpp
    utils/partitionLowerTriangle.cpp
    utils/iterator_utils.cpp
)

target_link_libraries(gmatrix_run_smoke_test
    ${FASTGXE_MKL_LIBS}
    ${FASTGXE_RUNTIME_LIBS}
)

add_test(
    NAME gmatrix_run_smoke_test
    COMMAND ${CMAKE_COMMAND} -E env
        OMP_NUM_THREADS=1
        KMP_AFFINITY=disabled
        $<TARGET_FILE:gmatrix_run_smoke_test>
)
