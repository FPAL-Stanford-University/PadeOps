#!/bin/bash

export COMPILER_ID=PGI
export FC=/opt/pgi/linux86-64/2016/mpi/openmpi/bin/mpif90
export CC=/opt/pgi/linux86-64/2016/mpi/openmpi/bin/mpicc
export CXX=/opt/pgi/linux86-64/2016/mpi/openmpi/bin/mpic++
export FFTW_PATH=/opt/fftw-3.3.4-PGI
export DECOMP_PATH=/opt/2decomp_fft-PGI
export VTK_IO_PATH=/home/akshays/Codes/PadeOps/dependencies/Lib_VTK_IO_PGI/build
