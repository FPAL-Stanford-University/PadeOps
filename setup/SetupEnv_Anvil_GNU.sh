#!/bin/bash

CWD=`pwd`
export FC="mpif90 -L/home/x-ryanhass/OpenBLAS/build/lib"
export COMPILER_ID=GNU
export FFTW_PATH=/home/x-ryanhass/codes/PadeOps/dependencies/GNUcompiled/fftw-3.3.10
export DECOMP_PATH=/home/x-ryanhass/codes/PadeOps/dependencies/GNUcompiled/2decomp_fft
export VTK_IO_PATH=/home/x-ryanhass/codes/PadeOps/dependencies/GNUcompiled/Lib_VTK_IO/build
export HDF5_PATH=/home/x-ryanhass/codes/PadeOps/dependencies/GNUcompiled/hdf5-1.8.18
export FFTPACK_PATH=/home/x-ryanhass/codes/PadeOps/dependencies/GNUcompiled/fftpack
