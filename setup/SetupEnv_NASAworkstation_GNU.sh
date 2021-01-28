#!/bin/bash

CWD=`pwd`
export COMPILER_ID=GNU
export FC="mpif90 -L/u/wk/aghate/OpenBLAS/build/lib"
export CC=mpicc
export CXX=mpic++
export FFTW_PATH=/u/wk/aghate/Codes/PadeOps/dependencies/fftw-3.3.5
export DECOMP_PATH=/u/wk/aghate/Codes/PadeOps/dependencies/2decomp_fft
export VTK_IO_PATH=/u/wk/aghate/Codes/PadeOps/dependencies/Lib_VTK_IO/build
export HDF5_PATH=/u/wk/aghate/Codes/PadeOps/dependencies/hdf5-1.8.18
export ARCH_OPT_FLAG="-march=skylake-avx512"
