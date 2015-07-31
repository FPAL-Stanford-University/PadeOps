file(GLOB utilities_source_files utilities/*.F90)
file(GLOB derivatives_source_files derivatives/*.F90)
file(GLOB filters_source_files filters/*.F90)
file(GLOB poisson_source_files poisson/*.F90)

# Where to look for FFTW
set(FFTW_LIBRARY_PATH "${FFTW_PATH}/lib")
set(FFTW_INCLUDE_PATH "${FFTW_PATH}/include")

# Include directories
include_directories(${MPI_INCLUDE_PATH} ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH})

# Link directories
link_directories(${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH})

# Create a library called PadeOps
add_library(PadeOps SHARED ${utilities_source_files} ${derivatives_source_files} ${filters_source_files} ${poisson_source_files})
target_link_libraries(PadeOps ${MPI_LIBRARIES})

if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(PadeOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()
