file(GLOB sproblems_source_files *.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )

# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the solids problem executables
foreach ( sproblemfile ${sproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/solids/" "" sproblemfile ${sproblemfile} )
    string( REPLACE ".F90" "" sproblemexec ${sproblemfile} )
  
    add_executable( ${sproblemexec} ${sproblemfile} )
    
    target_link_libraries( ${sproblemexec} PadeOps )

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${sproblemexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()

    # file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/solids/${sproblemexec}_files)
    # file(COPY ${sproblemexec}/input.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/solids/${sproblemexec}_files/input.dat)
    
endforeach()
