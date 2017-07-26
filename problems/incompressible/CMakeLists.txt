file(GLOB iproblems_source_files *.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${FFTPACK_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )
# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${FFTPACK_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the incompressible problem executables
foreach ( iproblemfile ${iproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/incompressible/" "" iproblemfile ${iproblemfile} )
    string( REPLACE ".F90" "" iproblemexec ${iproblemfile} )
  
    add_executable( ${iproblemexec} ${iproblemfile} )
    
    target_link_libraries( ${iproblemexec} IncompressibleOps 2decomp_fft fftpack fftw3)

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${iproblemexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
    
    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files)
    #    file(COPY ${iproblemexec}_files/*.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files) 
    file(GLOB iproblems_input_files ${iproblemexec}_files/*.dat)
    file(COPY ${iproblems_input_files} DESTINATION ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files) 
endforeach()
