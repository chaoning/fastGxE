cmake_minimum_required(VERSION 3.16)
project(GmatrixProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

# Set external library paths (modify accordingly)
set(MKLROOT /net/mulan/home/chaon/intel/oneapi/mkl/2024.1/)
set(GSLROOT /net/mulan/home/chaon/software/gsl-2.7.1/)
set(EIGENROOT /net/mulan/home/chaon/software/eigen-3.4.0/)
set(BOOSTROOT /net/mulan/home/chaon/software/boost_1_81_0/)
set(LBFGSppROOT /net/mulan/home/chaon/WORK/GMAT/external/LBFGSpp/include/)
set(UTILROOT /net/zootopia/disk1/chaon/WORK/GMAT/utils/)
set(SPDROOT /net/zootopia/disk1/chaon/WORK/GMAT/external/)
set(CLIROOT /net/zootopia/disk1/chaon/WORK/GMAT/external/)

# Include directories
include_directories(
    ${MKLROOT}/include
    ${GSLROOT}/include
    ${EIGENROOT}
    ${EIGENROOT}/unsupported
    ${BOOSTROOT}
    ${UTILROOT}
    ${LBFGSppROOT}
    ${SPDROOT}
    ${CLIROOT}
)

# Source files
set(SOURCE_FILES
    ../utils/string_utils.cpp
    ../utils/partitionLowerTriangle.cpp
    ../utils/iterator_utils.cpp
    ../utils/EigenMatrix_utils.cpp
    ../geno/geno.cpp
    ../gmatrix/gmatrix.cpp
    ../gmatrix/processGRM.cpp
    main.cpp
)

# Define the executable
add_executable(gmatrix ${SOURCE_FILES})

# Link Libraries
target_link_libraries(gmatrix
    -Wl,--start-group
    ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
    ${MKLROOT}/lib/intel64/libmkl_intel_thread.a
    ${MKLROOT}/lib/intel64/libmkl_core.a
    -Wl,--end-group
    ${GSLROOT}/lib/libgsl.a
    ${GSLROOT}/lib/libgslcblas.a
    -liomp5 -lpthread -lm -ldl
)
