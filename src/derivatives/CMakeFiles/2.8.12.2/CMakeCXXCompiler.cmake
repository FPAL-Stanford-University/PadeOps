set(CMAKE_CXX_COMPILER "/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin_ohpc/mpicxx")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "Intel")
set(CMAKE_CXX_COMPILER_VERSION "18.0.0.20180210")
set(CMAKE_CXX_PLATFORM_ID "Linux")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_COMPILER_IS_GNUCXX )
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;CPP)
set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()




set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "mpicxx;mpifort;mpi;mpigi;dl;rt;pthread;imf;svml;irng;stdc++;m;ipgo;decimal;cilkrts;stdc++;irc;svml;c;irc_s;dl;c")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/mpi/intel64/lib/release_mt;/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/mpi/intel64/lib;/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/ipp/lib/intel64;/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/compiler/lib/intel64_lin;/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/mkl/lib/intel64_lin;/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/tbb/lib/intel64/gcc4.1;/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/daal/lib/intel64_lin;/opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/tbb/lib/intel64_lin/gcc4.4;/usr/lib/gcc/x86_64-redhat-linux/4.8.5;/usr/lib64;/lib64;/usr/lib;/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")



