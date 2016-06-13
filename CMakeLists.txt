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

# Where to look for Lib_VTK_IO
set(VTK_IO_PATH $ENV{VTK_IO_PATH})
set(VTK_IO_LIBRARY_PATH "${VTK_IO_PATH}/lib")
set(VTK_IO_INCLUDE_PATH "${VTK_IO_PATH}/modules")

set(CMAKE_Fortran_COMPILER_ID $ENV{COMPILER_ID})

# Intel compiler suite
if ( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O5 -mkl -warn all -xhost -dynamic -qopt-report=2 -qopt-report-phase=vec -qoverride-limits")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -pg -mkl -check all -check noarg-temp-created -fpe0 -warn -traceback -debug extended -assume realloc_lhs -fstack-protector -assume protect_parens -implicitnone")
    endif()

# GNU compiler on OSX
elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "GNU_OSX")
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack -Wall -Wconversion -Wextra -Waliasing -ffree-form -ffree-line-length-none -ffast-math -funroll-loops -fno-protect-parens -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "-Og -g -pg -llapack -ffree-form -ffree-line-length-none -fbounds-check -ffpe-trap=zero,overflow,underflow -Wall -Wconversion -Wextra -Waliasing -Wsurprising -fbacktrace -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
    endif()

# Standard GNU compilers
elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack -Wall -Wconversion -Wextra -Waliasing -ffree-form -ffree-line-length-none -ffast-math -march=native -funroll-loops -fno-protect-parens")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "-Og -g -pg -llapack -ffree-form -ffree-line-length-none -fcheck=all -fbounds-check -ffpe-trap=zero,overflow,underflow -Wall -Wconversion -Wextra -Waliasing -Wsurprising -fbacktrace")
    endif()
endif()

# Add source file directories
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(problems)
