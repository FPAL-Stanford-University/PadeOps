#!/bin/bash

vpkg_devrequire intel
vpkg_devrequire intel-oneapi
vpkg_devrequire git
vpkg_devrequire cmake
vpkg_devrequire openmpi
CWD=/lustre/xg-phy230025/users/3108/PadeOps/PadeOps
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=${CWD}/dependencies/fftw-3.3.5
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
export HDF5_PATH=${CWD}/dependencies/hdf5-1.8.18
export FFTPACK_PATH=${CWD}/dependencies/fftpack
#export ARCH_OPT_FLAG="-xCORE-AVX2 -axMIC-AVX512"
export ARCH_OPT_FLAG="-xHost"
