#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpif90 
export CC="icc -lmpi"
export CXX="icc -lmpi"
export FFTW_PATH=$CWD/dependencies/fftw-3.3.10
export DECOMP_PATH=$CWD/dependencies/2decomp_fft
export VTK_IO_PATH=$CWD/dependencies/Lib_VTK_IO/build
export HDF5_PATH=$CWD/dependencies/hdf5-1.8.18
export FFTPACK_PATH=$CWD/dependencies/fftpack
export ARCH_OPT_FLAG="-axCORE-AVX2 -qopt-zmm-usage=high"
