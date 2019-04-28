#!/bin/bash

export COMPILER_ID=Intel
export FC=mpiifort
export CC=mpiicc
export CXX=mpiicpc
export FFTW_PATH=/work/04076/tg833754/stampede2/PadeOps_Shared/padeops_dependencies/fftw-3.3.5
export DECOMP_PATH=/work/04076/tg833754/stampede2/PadeOps_Shared/padeops_dependencies/2decomp_fft
export VTK_IO_PATH=/work/04076/tg833754/stampede2/PadeOps_Shared/padeops_dependencies/Lib_VTK_IO/build
export HDF5_PATH=/work/04076/tg833754/stampede2/PadeOps_Shared/padeops_dependencies/hdf5-1.8.18
export FFTPACK_PATH=/work/04076/tg833754/stampede2/PadeOps_Shared/fftpack
export ARCH_OPT_FLAG="-xCOMMON-AVX512 -axCORE-AVX512,MIC-AVX512 -qopt-zmm-usage=high"
