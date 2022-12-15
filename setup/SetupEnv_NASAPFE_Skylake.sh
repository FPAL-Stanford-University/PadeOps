#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpif90 
export CC="icc -lmpi"
export CXX="icc -lmpi"
export FFTW_PATH=/nobackup/aghate/PadeOps/dependencies/fftw-3.3.10
export DECOMP_PATH=/nobackup/aghate/PadeOps/dependencies/2decomp_fft
export VTK_IO_PATH=/nobackup/aghate/PadeOps/dependencies/Lib_VTK_IO/build
export HDF5_PATH=/nobackup/aghate/PadeOps/dependencies/hdf5-1.8.18
export FFTPACK_PATH=/nobackup/aghate/PadeOps/dependencies/fftpack
export ARCH_OPT_FLAG="-xCORE-AVX512 -qopt-zmm-usage=high -qopt-report=5"
