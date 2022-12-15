#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpif90 
export CC="icc -lmpi"
export CXX=mpicxx
export FFTW_PATH=/u/wk/aghate/fftw-3.3.8
export DECOMP_PATH=/u/wk/aghate/PadeOps/dependencies/2decomp_fft
export VTK_IO_PATH=/u/wk/aghate/PadeOps/dependencies/Lib_VTK_IO/build
export HDF5_PATH=/u/wk/aghate/hdf5-1.8.18
