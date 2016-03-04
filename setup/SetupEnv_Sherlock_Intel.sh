#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.4
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
