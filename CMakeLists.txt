cmake_minimum_required(VERSION 3.16)
project(fastgxeProject)

set(CMAKE_CXX_COMPILER "/home/chaoning/intel/oneapi/compiler/2025.2/bin/icpx")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qopenmp")

# Set external library paths (modify accordingly)
set(MKLROOT /home/chaoning/intel/oneapi/mkl/2025.2/)
set(GSLROOT /home/chaoning/software/fastgxe/external/gsl-2.7/)
set(EIGENROOT /home/chaoning/software/fastgxe/external/eigen-3.4.0/)
set(LBFGSppROOT /home/chaoning/software/fastgxe/external/LBFGSpp/include/)
set(BOOSTROOT /home/chaoning/software/fastgxe/external/boost_1_89_0/)
set(UTILROOT /home/chaoning/software/fastgxe/utils/)
set(GMATRIXROOT /home/chaoning/software/fastgxe/gmatrix/)
set(SPDROOT /home/chaoning/software/fastgxe/external/)
set(CLIROOT /home/chaoning/software/fastgxe/external/)

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
    utils/EigenMatrix_utils.cpp
    utils/geno.cpp
    utils/iterator_utils.cpp
    utils/partitionLowerTriangle.cpp
    utils/phen.cpp
    utils/random_num.cpp
    utils/someCommonFun.cpp
    utils/string_utils.cpp
    gmatrix/gmatrix.cpp
    gmatrix/processGRM.cpp
    fastgxe/fastgxe.cpp
    fastgxe/mom.cpp
    fastgxe/mmsusie.cpp
    fastgxe/main.cpp
)

# Define the executable
add_executable(fastgxe ${SOURCE_FILES})

# Link Libraries
target_link_libraries(fastgxe
    -Wl,--start-group
    ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
    ${MKLROOT}/lib/intel64/libmkl_intel_thread.a
    ${MKLROOT}/lib/intel64/libmkl_core.a
    -Wl,--end-group
    ${GSLROOT}/lib/libgsl.a
    ${GSLROOT}/lib/libgslcblas.a
    -liomp5 -lpthread -lm -ldl
)
