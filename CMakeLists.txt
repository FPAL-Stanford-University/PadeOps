cmake_minimum_required(VERSION 2.8)
project(PadeOps)

enable_language(Fortran)

set( CMAKE_MODULE_PATH ${PadeOps_SOURCE_DIR} )

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Where to look for FFTW
set(FFTW_PATH $ENV{FFTW_PATH})
set(FFTW_LIBRARY_PATH "${FFTW_PATH}/lib")
set(FFTW_INCLUDE_PATH "${FFTW_PATH}/include")

# Where to look for DECOMP&FFT
set(DECOMP_PATH $ENV{DECOMP_PATH})
set(DECOMP_LIBRARY_PATH "${DECOMP_PATH}/lib")
set(DECOMP_INCLUDE_PATH "${DECOMP_PATH}/include")

set(CMAKE_Fortran_COMPILER_ID $ENV{COMPILER_ID})
if ( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -xhost -dynamic -qopt-report=2 -qopt-report-phase=vec")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -pg -check all -fpe0 -warn -traceback -debug extended")
    endif()
elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wconversion -Wextra -Waliasing -ffree-form -ffree-line-length-none -ffast-math -march=native -funroll-loops -fno-protect-parens")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "-Og -g -pg -ffree-form -ffree-line-length-none -fbounds-check -ffpe-trap=zero,overflow,underflow -Wall -Wconversion -Wextra -Waliasing -Wsurprising -fbacktrace")
    endif()
endif()

# Add source file directories
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(problems)
