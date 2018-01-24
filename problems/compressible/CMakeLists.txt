file(GLOB cproblems_source_files *.F90)

# Include directories
include_directories(${MPI_INCLUDE_PATH} ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${VTK_IO_INCLUDE_PATH} ${HDF5_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src )

# Link directories
link_directories(${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${VTK_IO_LIBRARY_PATH} ${HDF5_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the compressible problem executables
foreach ( cproblemfile ${cproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/compressible/" "" cproblemfile ${cproblemfile} )
    string( REPLACE ".F90" "" cproblemexec ${cproblemfile} )
  
    add_executable( ${cproblemexec} ${cproblemfile} )
    
    target_link_libraries( ${cproblemexec} CompressibleOps )

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${cproblemexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()

    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/compressible/${cproblemexec}_files)
    file(COPY ${cproblemexec}_files/input.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/compressible/${cproblemexec}_files) 
endforeach()
