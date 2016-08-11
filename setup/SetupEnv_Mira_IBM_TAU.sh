#!/bin/bash

CWD=`pwd`
export COMPILER_ID=IBM
export FC=tau_f90.sh
export CC=tau_cc.sh
export CXX=tau_cxx.sh
export FFTW_PATH=/soft/libraries/alcf/current/xl/FFTW3
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
