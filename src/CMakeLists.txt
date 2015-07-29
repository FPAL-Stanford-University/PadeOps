file(GLOB utilities_source_files utilities/*.F90)
file(GLOB poisson_source_files poisson/*.F90)
# file(GLOB derivatives_source_files derivatives/*.F90)

# Where to look for FFTW
set(FFTW_LIBRARY_PATH "${FFTW_PATH}/lib")
set(FFTW_INCLUDE_PATH "${FFTW_PATH}/include")

# Include directories
<<<<<<< HEAD
include_directories(${FFTW_INCLUDE_PATH})
=======
include_directories(${MPI_INCLUDE_PATH} ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH})
>>>>>>> 44f523fa8667c0fb162c2d56941e559e105a90f8

# Link directories
link_directories(${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH})

# Create a library called PadeOps
<<<<<<< HEAD
add_library(PadeOps SHARED ${utilities_source_files} ${derivatives_source_files})
# target_link_libraries(PadeOps ${MPI_LIBRARIES})
=======
add_library(PadeOps SHARED ${utilities_source_files} ${poisson_source_files} ${derivatives_source_files})
target_link_libraries(PadeOps ${MPI_LIBRARIES})
>>>>>>> 44f523fa8667c0fb162c2d56941e559e105a90f8

if (MPI_Fortran_COMPILER_FLAGS)
    set_target_properties(PadeOps PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
endif()

# if (MPI_LINK_FLAGS)
#     set_target_properties(PadeOps PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
# endif()
