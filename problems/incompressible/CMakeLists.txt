file(GLOB iproblems_source_files *.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )
# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

#set(GLIBC_DIR "/ocean/projects/phy230038p/awu2/glibc_install/glibc-2.29")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--dynamic-linker=${GLIBC_DIR}/lib/ld-linux-x86-64.so.2")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${GLIBC_DIR}/lib -L${FFTW_LIBRARY_PATH} -L${DECOMP_LIBRARY_PATH} -L${PadeOps_BINARY_DIR}/src -L/lib64")

# new include and link dirs
#include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src ${GLIBC_DIR}/include)
#link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src ${GLIBC_DIR}/lib)

#set(CMAKE_INSTALL_RPATH "${GLIBC_DIR}/lib:${FFTW_LIBRARY_PATH}:${DECOMP_INCLUDE_PATH}:${PadeOps_BINARY_DIR}/src:/lib64")
#set(CMAKE_BUILD_RPATH "${GLIBC_DIR}/lib:${FFTW_LIBRARY_PATH}:${DECOMP_INCLUDE_PATH}:${PadeOps_BINARY_DIR}/src:/lib64")

# Create the incompressible problem executables
foreach ( iproblemfile ${iproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/incompressible/" "" iproblemfile ${iproblemfile} )
    string( REPLACE ".F90" "" iproblemexec ${iproblemfile} )
  
    add_executable( ${iproblemexec} ${iproblemfile} )
    
    target_link_libraries( ${iproblemexec} PUBLIC IncompressibleOps 2decomp_fft fftw3 tf_load pytorch_load ${MPI_LIBRARIES})
    #target_link_libraries( ${iproblemexec} IncompressibleOps 2decomp_fft fftw3 tf_load pytorch_load ${MPI_LIBRARIES})
    #target_link_libraries( ${iproblemexec} PRIVATE pytorch_load)
    #target_link_libraries( ${iproblemexec} IncompressibleOps 2decomp_fft fftw3 tf_load ${MPI_LIBRARIES})
    #target_link_libraries( ${iproblemexec} IncompressibleOps 2decomp_fft fftw3 ${MPI_LIBRARIES} tf_load)
    #target_link_libraries( ${iproblemexec} IncompressibleOps 2decomp_fft fftw3 ${MPI_LIBRARIES})


    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${iproblemexec} PROPERTIES LINKER_LANGUAGE Fortran COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    else()
        set_target_properties(${iproblemexec} PROPERTIES LINKER_LANGUAGE Fortran)
    endif()
    
    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files)
    #    file(COPY ${iproblemexec}_files/*.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files) 
    file(GLOB iproblems_input_files ${iproblemexec}_files/*.dat)
    file(COPY ${iproblems_input_files} DESTINATION ${PadeOps_BINARY_DIR}/problems/incompressible/${iproblemexec}_files) 
endforeach()
