file(GLOB utilities_source_files utilities/*.F90)
file(GLOB derivatives_source_files derivatives/*.F90)
file(GLOB filters_source_files filters/*.F90)
file(GLOB io_source_files io/*.F90)
file(GLOB grid_source_files grid/*.F90)
file(GLOB compressible_source_files compressible/*.F90)
file(GLOB solid_source_files solid/*.F90)
file(GLOB poisson_source_files incompressible/poisson/*.F90)
file(GLOB incompressible_source_files incompressible/*.F90)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{FORTRUAN_STDLIB_PATH}")

find_package(fortran_stdlib REQUIRED)
if ( fortran_stdlib_FOUND )
  message(STATUS "fortran_stdlib FOUND ${fortran_stdlib}")
else()
  message(FATAL "fortran_stdlib NOT found")
endif()

# Where to look for FFTW
set(FFTW_LIBRARY_PATH "${FFTW_PATH}/lib")
set(FFTW_INCLUDE_PATH "${FFTW_PATH}/include")

# Include directories
include_directories(${MPI_INCLUDE_PATH} ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${VTK_IO_INCLUDE_PATH} ${HDF5_INCLUDE_PATH})

# Link directories
link_directories(${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${VTK_IO_LIBRARY_PATH} ${HDF5_LIBRARY_PATH})

# Create a library called PadeOps
add_library(PadeOps STATIC ${utilities_source_files} ${derivatives_source_files} ${filters_source_files} ${io_source_files} )
if ( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
   target_link_libraries(PadeOps fftw3 2decomp_fft ${VTK_IO_LIBRARY_PATH}/libVTK_IO.a ${HDF5_LIBRARY_PATH}/libhdf5hl_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5_hl.a ${HDF5_LIBRARY_PATH}/libhdf5_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5.a ${MPI_LIBRARIES} -lz -ldl -lm -Wl,-rpath -Wl,${HDF5_LIBRARY_PATH})
else()
   target_link_libraries(PadeOps fftw3 2decomp_fft blas lapack ${VTK_IO_LIBRARY_PATH}/libVTK_IO.a ${HDF5_LIBRARY_PATH}/libhdf5hl_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5_hl.a ${HDF5_LIBRARY_PATH}/libhdf5_fortran.a ${HDF5_LIBRARY_PATH}/libhdf5.a ${MPI_LIBRARIES} -lz -ldl -lm -Wl,-rpath -Wl,${HDF5_LIBRARY_PATH})
endif()

if ( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
   target_link_libraries(PadeOps fftw3 2decomp_fft ${VTK_IO_LIBRARY_PATH}/libVTK_IO.a ${MPI_LIBRARIES})
elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   target_link_libraries(PadeOps fftw3 2decomp_fft blas lapack ${VTK_IO_LIBRARY_PATH}/libVTK_IO.a ${MPI_LIBRARIES})
elseif ( CMAKE_Fortran_COMPILER_ID MATCHES "GNU_OSX")
   target_link_libraries(PadeOps fftw3 2decomp_fft blas lapack ${VTK_IO_LIBRARY_PATH}/libVTK_IO.a ${MPI_LIBRARIES})
endif()
   
if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(PadeOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

# Create a library called CompressibleOps
add_library(CompressibleOps STATIC ${grid_source_files} ${compressible_source_files})
target_link_libraries(CompressibleOps PadeOps ${MPI_LIBRARIES})

if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(CompressibleOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

# Create a library called SolidOps
add_library(SolidOps STATIC ${grid_source_files} ${solid_source_files})
target_link_libraries(SolidOps CompressibleOps ${MPI_LIBRARIES})

if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(SolidOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

# Create a library called IncompressibleOps
add_library(IncompressibleOps STATIC ${grid_source_files} ${poisson_source_files} ${incompressible_source_files})
target_link_libraries(IncompressibleOps PUBLIC 
  PadeOps 
  ${MPI_LIBRARIES}
  fortran_stdlib::fortran_stdlib
)

if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(IncompressibleOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()
