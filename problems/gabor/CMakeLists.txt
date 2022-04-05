#file(GLOB GaborProblems_source_files *.F90)
file(GLOB GaborProblems_source_files *.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )
# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the incompressible problem executables
foreach ( GaborProblemFile ${GaborProblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/gabor/" "" GaborProblemFile ${GaborProblemFile} )
    string( REPLACE ".F90" "" GaborProblemexec ${GaborProblemFile} )
  
    add_executable( ${GaborProblemexec} ${GaborProblemFile} )
    
    target_link_libraries( ${GaborProblemexec} IncompressibleOps 2decomp_fft fftw3)

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${GaborProblemexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
    
    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/incompressible/${GaborProblemexec}_files)
    #    file(COPY ${GaborProblemexec}_files/*.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/incompressible/${GaborProblemexec}_files) 
    file(GLOB GaborProblems_input_files ${GaborProblemexec}_files/*.dat)
    file(COPY ${GaborProblems_input_files} DESTINATION ${PadeOps_BINARY_DIR}/problems/incompressible/${GaborProblemexec}_files) 
endforeach()
