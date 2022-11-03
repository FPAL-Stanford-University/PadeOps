#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpiifort 
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=/home/x-aditya90/PadeOps/dependencies/fftw-3.3.10
export DECOMP_PATH=/home/x-aditya90/PadeOps/dependencies/2decomp_fft
export VTK_IO_PATH=/home/x-aditya90/PadeOps/dependencies/Lib_VTK_IO/build
export HDF5_PATH=/home/x-aditya90/PadeOps/dependencies/hdf5-1.8.18
export FFTPACK_PATH=/home/x-aditya90/PadeOps/dependencies/fftpack
export ARCH_OPT_FLAG="-axCORE-AVX2 -qopt-zmm-usage=high"
