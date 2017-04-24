#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpif90
export CC=mpicc
export CXX=mpic++
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.5-bigmem
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft-bigmem
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO-bigmem/build
export HDF5_PATH=${CWD}/dependencies/hdf5-1.8.18-bigmem
