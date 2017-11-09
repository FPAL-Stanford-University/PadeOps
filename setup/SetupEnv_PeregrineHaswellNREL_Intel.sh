#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F90=ifort
export I_MPI_F77=ifort
export FFTW_PATH=${CWD}/dependencies/haswell/fftw-3.3.5
export DECOMP_PATH=${CWD}/dependencies/haswell/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/haswell/Lib_VTK_IO/build
export HDF5_PATH=${CWD}/dependencies/hdf5-1.8.18
export ARCH_OPT_FLAG="-xCORE-AVX2"
