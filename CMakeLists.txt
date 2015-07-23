cmake_minimum_required(VERSION 2.8)
project(PadeOps)

enable_language(Fortran)

set( CMAKE_MODULE_PATH ${PadeOps_SOURCE_DIR} )

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(FFTW_PATH $ENV{FFTW_PATH})

# Where to look for FFTW
set(FFTW_LIBRARY_PATH "${FFTW_PATH}/lib")
set(FFTW_INCLUDE_PATH "${FFTW_PATH}/include")

find_package(MPI REQUIRED)

message(STATUS "MPI Fortran compiler is: ${MPI_Fortran_COMPILER}" )
message(STATUS "Fortran compiler ID is: ${CMAKE_Fortran_COMPILER_ID}" )

if ( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -ipo -xhost -mcmodel=\"medium\" -i8 -qopt-report=2 -qopt-report-phase=vec")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -pg -check all -fpe0 -warn -traceback -debug extended")
    endif()
endif()

# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lfftw3")

# Add source file directories
add_subdirectory(src)
add_subdirectory(tests)
