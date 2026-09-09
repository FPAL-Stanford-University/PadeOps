set(CMAKE_C_COMPILER "/opt/packages/oneapi/v2023.2.0/mpi/2021.10.0/bin/mpiicc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "Intel")
set(CMAKE_C_COMPILER_VERSION "19.1.3.20200925")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_EXTENSIONS_COMPUTED_DEFAULT "ON")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")
set(CMAKE_C17_COMPILE_FEATURES "")
set(CMAKE_C23_COMPILE_FEATURES "")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "GNU")
set(CMAKE_C_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_C_SIMULATE_VERSION "8.5.0")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_MT "")
set(CMAKE_TAPI "CMAKE_TAPI-NOTFOUND")
set(CMAKE_COMPILER_IS_GNUCC )
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)

set(CMAKE_C_COMPILER_ENV_VAR "CC")

set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)
set(CMAKE_C_LINKER_DEPFILE_SUPPORTED )

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/opt/packages/oneapi/v2023.2.0/mpi/2021.10.0/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/pstl/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/ipp/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/tbb/include;/opt/intel/compilers_and_libraries_2020.4.304/linux/daal/include;/opt/packages/oneapi/v2023.2.0/mkl/2023.2.0/include;/opt/packages/oneapi/v2023.2.0/ipp/2021.9.0/include;/opt/packages/oneapi/v2023.2.0/ippcp/2021.8.0/include;/opt/packages/oneapi/v2023.2.0/dpl/2022.2.0/linux/include;/opt/packages/oneapi/v2023.2.0/dpcpp-ct/2023.2.0/include;/opt/packages/oneapi/v2023.2.0/dnnl/2023.2.0/cpu_dpcpp_gpu_dpcpp/include;/opt/packages/oneapi/v2023.2.0/dev-utilities/2021.10.0/include;/opt/packages/oneapi/v2023.2.0/dal/2023.2.0/include;/opt/packages/oneapi/v2023.2.0/ccl/2021.10.0/include/cpu_gpu_dpcpp;/opt/packages/oneapi/v2023.2.0/tbb/2021.10.0/include;/jet/packages/intel/compilers_and_libraries_2020.4.304/linux/compiler/include/intel64;/jet/packages/intel/compilers_and_libraries_2020.4.304/linux/compiler/include/icc;/jet/packages/intel/compilers_and_libraries_2020.4.304/linux/compiler/include;/usr/local/include;/usr/lib/gcc/x86_64-redhat-linux/8/include;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "mpifort;mpi;dl;rt;pthread;imf;svml;irng;m;ipgo;decimal;cilkrts;stdc++;gcc;gcc_s;irc;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/opt/packages/oneapi/v2023.2.0/mpi/2021.10.0/lib/release;/opt/packages/oneapi/v2023.2.0/mpi/2021.10.0/lib;/opt/intel/clck/2019.5/lib/intel64;/opt/intel/compilers_and_libraries_2020.4.304/linux/ipp/lib/intel64;/opt/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin;/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin;/opt/intel/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64/gcc4.8;/opt/intel/compilers_and_libraries_2020.4.304/linux/daal/lib/intel64_lin;/opt/packages/oneapi/v2023.2.0/mpi/2021.10.0/libfabric/lib;/opt/packages/oneapi/v2023.2.0/mkl/2023.2.0/lib/intel64;/opt/packages/oneapi/v2023.2.0/ipp/2021.9.0/lib/intel64;/opt/packages/oneapi/v2023.2.0/ippcp/2021.8.0/lib/intel64;/opt/packages/oneapi/v2023.2.0/dnnl/2023.2.0/cpu_dpcpp_gpu_dpcpp/lib;/opt/packages/oneapi/v2023.2.0/dal/2023.2.0/lib/intel64;/opt/packages/oneapi/v2023.2.0/ccl/2021.10.0/lib/cpu_gpu_dpcpp;/opt/packages/oneapi/v2023.2.0/compiler/2023.2.1/linux/lib;/opt/packages/oneapi/v2023.2.0/tbb/2021.10.0/lib/intel64/gcc4.8;/jet/packages/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin;/usr/lib/gcc/x86_64-redhat-linux/8;/usr/lib64;/lib64;/usr/lib;/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
