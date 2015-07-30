#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpif90
export CC=mpicc
export CXX=mpicxx
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.4
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
