#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=${CWD}/dependencies/sherlock2/fftw-3.3.5
export DECOMP_PATH=${CWD}/dependencies/sherlock2/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/sherlock2/Lib_VTK_IO/build
export HDF5_PATH=${CWD}/dependencies/sherlock2/hdf5-1.8.18
