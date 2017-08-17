file(GLOB problems_source_files *.F90)
file(GLOB cproblems_source_files compressible/*.F90)
file(GLOB iproblems_source_files incompressible/*.F90)
file(GLOB turbproblems_source_files turbines/*.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )

# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the problem executables
foreach ( problemfile ${problems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/" "" problemfile ${problemfile} )
    string( REPLACE ".F90" "" problemexec ${problemfile} )
  
    add_executable( ${problemexec} ${problemfile} )
    
    target_link_libraries( ${problemexec} PadeOps )

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${problemexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
    
endforeach()

add_subdirectory(compressible)
add_subdirectory(incompressible)
add_subdirectory(turbines)
add_subdirectory(postprocessing)
add_subdirectory(solids)
