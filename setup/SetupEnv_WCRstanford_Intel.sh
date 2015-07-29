#!/bin/bash

CWD=`pwd`

export FC=ifort
export FFTW_PATH=`dirname ${CWD}`/fftw-3.3.4
export DECOMP_PATH=`dirname ${CWD}`/2decomp_fft
