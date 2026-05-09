cmake_minimum_required(VERSION 3.16)
project(fastgxeProject)

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

