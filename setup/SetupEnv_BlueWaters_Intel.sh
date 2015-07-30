#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=ftn
export CC=cc
export CXX=CC
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.4
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
