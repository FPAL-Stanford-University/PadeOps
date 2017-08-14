#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=/share/software/user/open/fftw/3.3.6/
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
export HDF5_PATH=/share/software/user/open/hdf5/1.10.0p1/
export FFTPACK_PATH=${CWD}/dependencies/fftpack
