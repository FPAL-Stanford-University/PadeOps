#!/bin/bash

CWD=`pwd`

export COMPILER_ID=GNU
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.4-GNU
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft-GNU
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO-GNU/build
