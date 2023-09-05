cmake_minimum_required(VERSION 3.10)
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

# Where to look for HDF5
set(HDF5_PATH $ENV{HDF5_PATH})
set(HDF5_LIBRARY_PATH "${HDF5_PATH}/lib")
set(HDF5_INCLUDE_PATH "${HDF5_PATH}/include")

# Where to look for fftpack
set(FFTPACK_PATH $ENV{FFTPACK_PATH})
set(FFTPACK_LIBRARY_PATH "${FFTPACK_PATH}/lib")
set(FFTPACK_INCLUDE_PATH "${FFTPACK_PATH}/include")


set(CMAKE_Fortran_COMPILER_ID $ENV{COMPILER_ID})



# Intel compiler suite
if ( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )	
	if ($ENV{ARCH_OPT_FLAG})
	    set(OPTFLAG "-xhost")
	else()
	    set(OPTFLAG $ENV{ARCH_OPT_FLAG})
	endif()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -traceback -mkl -heap-arrays 1024 -warn all ${OPTFLAG} -dynamic -qopt-report=2 -qopt-report-phase=vec -qopenmp")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -traceback -g -mkl -heap-arrays 1024 -check all -check noarg-temp-created -fpe0 -warn -traceback -debug extended -assume realloc_lhs -fstack-protector -assume protect_parens -implicitnone")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Profile" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -traceback -mkl -heap-arrays 1024 -warn all ${OPTFLAG} -dynamic -qopt-report=2 -qopt-report-phase=vec -qopenmp -pg -g")
    elseif ( CMAKE_BUILD_TYPE MATCHES "O2" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -traceback -mkl -heap-arrays 1024 -warn all ${OPTFLAG} -dynamic -qopt-report=2 -qopt-report-phase=vec -qopenmp")
    endif()

# GNU compiler on OSX
elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "GNU_OSX")
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wconversion -Wextra -Waliasing -fopenmp -ffree-form -ffree-line-length-none -ffast-math -funroll-loops -fno-protect-parens -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "-Og -g -fbacktrace -pg -llapack -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -ffpe-trap=zero,overflow -Wall -Wconversion -Wextra -Waliasing -Wsurprising -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
    endif()

# Standard GNU compilers
elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
	if ($ENV{ARCH_OPT_FLAG})
	    set(OPTFLAG "-march=native")
	else()
	    set(OPTFLAG $ENV{ARCH_OPT_FLAG})
	endif()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wconversion -Wextra -Waliasing -ffree-form -ffree-line-length-none -ffast-math ${OPTFLAG} -funroll-loops -fno-protect-parens -fopenmp -fallow-argument-mismatch")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "-Og -g -fbacktrace -pg -llapack -ffree-form -ffree-line-length-none -fcheck=all -fbounds-check -ffpe-trap=zero,overflow -Wall -Wconversion -Wextra -Waliasing -Wsurprising")
    endif()

 elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wconversion -Wextra -Waliasing -ffree-form -ffree-line-length-none -ffast-math -march=native -funroll-loops -fno-protect-parens -fopenmp")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "-Og -g -fbacktrace -pg -llapack -ffree-form -ffree-line-length-none -fcheck=all -fbounds-check -ffpe-trap=zero,overflow -Wall -Wconversion -Wextra -Waliasing -Wsurprising")
    endif()
 
 elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "IBM")
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O5 -qsuffix=cpp=f90 -qxlf2003=polymorphic")
    elseif ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set(CMAKE_Fortran_FLAGS "-O0 -qsuffix=cpp=f90 -qxlf2003=polymorphic")
    endif()
 
endif()

# Add source file directories
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(problems)
