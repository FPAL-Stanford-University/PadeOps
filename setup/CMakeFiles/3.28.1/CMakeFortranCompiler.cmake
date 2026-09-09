set(CMAKE_Fortran_COMPILER "/opt/packages/oneapi/v2023.2.0/mpi/2021.10.0/bin/mpiifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "19.1.3.20200925")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_TAPI "CMAKE_TAPI-NOTFOUND")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95;f03;F03;f08;F08)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
set(CMAKE_Fortran_LINKER_DEPFILE_SUPPORTED )
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/packages/oneapi/v2023.2.0/mpi/2021.10.0/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/ipp/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/pstl/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/tbb/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/daal/include;/opt/packages/oneapi/v2023.2.0/mkl/2023.2.0/include;/opt/packages/oneapi/v2023.2.0/ipp/2021.9.0/include;/opt/packages/oneapi/v2023.2.0/ippcp/2021.8.0/include;/opt/packages/oneapi/v2023.2.0/dpl/2022.2.0/linux/include;/opt/packages/oneapi/v2023.2.0/dpcpp-ct/2023.2.0/include;/opt/packages/oneapi/v2023.2.0/dnnl/2023.2.0/cpu_dpcpp_gpu_dpcpp/include;/opt/packages/oneapi/v2023.2.0/dev-utilities/2021.10.0/include;/opt/packages/oneapi/v2023.2.0/dal/2023.2.0/include;/opt/packages/oneapi/v2023.2.0/ccl/2021.10.0/include/cpu_gpu_dpcpp;/opt/packages/oneapi/v2023.2.0/tbb/2021.10.0/include;/jet/packages/intel/compilers_and_libraries_2020.4.304/linux/compiler/include/intel64;/jet/packages/intel/compilers_and_libraries_2020.4.304/linux/compiler/include/icc;/jet/packages/intel/compilers_and_libraries_2020.4.304/linux/compiler/include;/usr/local/include;/usr/lib/gcc/x86_64-redhat-linux/8/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/intel/clck/2019.5/lib/intel64;/opt/intel/compilers_and_libraries_2020.4.304/linux/ipp/lib/intel64;/opt/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin;/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin;/opt/intel/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64/gcc4.8;/opt/intel/compilers_and_libraries_2020.4.304/linux/daal/lib/intel64_lin;/opt/packages/oneapi/v2023.2.0/mpi/2021.10.0/libfabric/lib;/opt/packages/oneapi/v2023.2.0/mkl/2023.2.0/lib/intel64;/opt/packages/oneapi/v2023.2.0/ipp/2021.9.0/lib/intel64;/opt/packages/oneapi/v2023.2.0/ippcp/2021.8.0/lib/intel64;/opt/packages/oneapi/v2023.2.0/dnnl/2023.2.0/cpu_dpcpp_gpu_dpcpp/lib;/opt/packages/oneapi/v2023.2.0/dal/2023.2.0/lib/intel64;/opt/packages/oneapi/v2023.2.0/ccl/2021.10.0/lib/cpu_gpu_dpcpp;/opt/packages/oneapi/v2023.2.0/compiler/2023.2.1/linux/lib;/opt/packages/oneapi/v2023.2.0/tbb/2021.10.0/lib/intel64/gcc4.8;/jet/packages/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin;/usr/lib/gcc/x86_64-redhat-linux/8;/usr/lib64;/lib64;/usr/lib;/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
