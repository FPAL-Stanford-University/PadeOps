#!/bin/bash

CWD=`pwd`
export COMPILER_ID=GNU
export FC=mpif90
export CC=mpicc
export CXX=mpic++
export FFTW_PATH=/u/wk/aghate/fftw-3.3.8
export DECOMP_PATH=/u/wk/aghate/PadeOps/dependencies/2decomp_fft
export VTK_IO_PATH=/u/wk/aghate/PadeOps/dependencies/Lib_VTK_IO/build
export HDF5_PATH=/u/wk/aghate/hdf5-1.8.18
