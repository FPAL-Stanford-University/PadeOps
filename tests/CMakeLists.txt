file(GLOB tests_source_files *.F90)

# set(tests_source_files "test_cd10.F90;test_cd06.F90")

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${VTK_IO_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src )

# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${VTK_IO_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src "/soft/libraries/alcf/current/xl/LAPACK64/lib" )

# Create the test executables
foreach ( testfile ${tests_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/tests/" "" testfile ${testfile} )
    string( REPLACE ".F90" "" testexec ${testfile} )
  
    add_executable( ${testexec} ${testfile} )
    
    target_link_libraries( ${testexec} IncompressibleOps CompressibleOps 2decomp_fft fftw3 lapack ${VTK_IO_LIBRARY_PATH}/libVTK_IO.a )

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${testexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
    
    # if (MPI_LINK_FLAGS)
    #     set_target_properties(${testexec} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    # endif()

endforeach()
