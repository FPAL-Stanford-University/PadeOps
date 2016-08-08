#!/bin/bash

CWD=`pwd`
export COMPILER_ID=IBM
export FC=mpixlf2008_r
export CC=mpixlc_r
export CXX=mpixlcxx_r
export FFTW_PATH=/soft/libraries/alcf/current/xl/FFTW3
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
