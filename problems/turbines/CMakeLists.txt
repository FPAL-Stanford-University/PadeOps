file(GLOB turbproblems_source_files *.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )

# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the turbines problem executables
foreach ( turbproblemfile ${turbproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/turbines/" "" turbproblemfile ${turbproblemfile} )
    string( REPLACE ".F90" "" turbproblemexec ${turbproblemfile} )
  
    add_executable( ${turbproblemexec} ${turbproblemfile} )
    
    target_link_libraries( ${turbproblemexec} IncompressibleOps 2decomp_fft fftw3)

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${turbproblemexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
    
    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/turbines/${turbproblemexec}_files)
    #    file(COPY ${turbproblemexec}_files/*.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/turbines/${turbproblemexec}_files) 
    file(GLOB turbproblems_input_files ${turbproblemexec}_files/*.dat)
    file(COPY ${turbproblems_input_files} DESTINATION ${PadeOps_BINARY_DIR}/problems/turbines/${turbproblemexec}_files) 
endforeach()
