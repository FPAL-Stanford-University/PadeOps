#!/bin/bash

CWD=`pwd`
export COMPILER_ID=GNU
export FC=mpif90
export CC=mpicc
export CXX=mpicxx
export FFTW_PATH=/soft/libraries/alcf/current/gcc/FFTW3
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft-GNU
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
