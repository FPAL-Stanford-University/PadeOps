#!/bin/bash

CWD=`pwd`
export COMPILER_ID=GNU
export FC=mpif90
export CC=mpicc
export CXX=mpicxx
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.4_gnu
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft_gnu
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO_gnu/build
export HDF5_PATH=${CWD}/dependencies/hdf5-1.8.18_gnu
