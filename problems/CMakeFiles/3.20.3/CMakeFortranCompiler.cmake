set(CMAKE_Fortran_COMPILER "/share/software/user/restricted/impi/2018/impi/2018.0.128/bin64/mpiifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "18.0.0.20170811")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
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





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/share/software/user/restricted/impi/2018/impi/2018.0.128/intel64/include;/share/software/user/open/python/2.7.13/include/python2.7;/share/software/user/open/python/2.7.13/include;/share/software/user/open/sqlite/3.18.0/include;/share/software/user/open/tcltk/8.6.6/include;/share/software/user/open/xz/5.2.3/include;/share/software/user/open/zlib/1.2.11/include;/share/software/user/open/curl/7.54.0/include;/share/software/user/restricted/imkl/2018/mkl/include;/share/software/user/restricted/impi/2018/impi/2018.0.128/include64;/share/software/user/restricted/icc/2018/include;/share/software/user/restricted/ifort/2018/include;/share/software/user/open/libressl/2.5.3/include;/share/software/user/restricted/ifort/2018/compilers_and_libraries_2018.0.128/linux/compiler/include/intel64;/share/software/user/restricted/ifort/2018/compilers_and_libraries_2018.0.128/linux/compiler/include;/usr/local/include;/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "mpifort;mpi;mpigi;dl;rt;pthread;ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/share/software/user/restricted/impi/2018/impi/2018.0.128/intel64/lib/release_mt;/share/software/user/restricted/impi/2018/impi/2018.0.128/intel64/lib;/share/software/user/open/python/2.7.13/lib;/share/software/user/open/sqlite/3.18.0/lib;/share/software/user/open/tcltk/8.6.6/lib;/share/software/user/open/xz/5.2.3/lib;/share/software/user/open/zlib/1.2.11/lib;/share/software/user/open/curl/7.54.0/lib;/share/software/user/restricted/imkl/2018/mkl/lib/intel64;/share/software/user/restricted/impi/2018/impi/2018.0.128/lib64;/share/software/user/restricted/icc/2018/lib/intel64;/share/software/user/restricted/ifort/2018/lib/intel64;/share/software/user/open/libressl/2.5.3/lib;/share/software/user/restricted/ifort/2018/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64_lin;/usr/lib/gcc/x86_64-redhat-linux/4.8.5;/usr/lib64;/lib64;/usr/lib;/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
