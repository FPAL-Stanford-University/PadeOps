#!/bin/bash

#module load intel-compiler
#module load intel-icc
#module load intel-oneapi
#module load intel/20.4
#module load gcc
module load lazygit
module load openmpi/4.0.2-intel20.4
CWD=/ocean/projects/phy230038p/awu2/PadeOps/PadeOps
source /ocean/projects/phy230038p/awu2/2022.3.0.8767/setvars.sh --force
export COMPILER_ID=Intel
export FC=/ocean/projects/phy230038p/awu2/2022.3.0.8767/mpi/2021.7.0/bin/mpiifort
export CC=/ocean/projects/phy230038p/awu2/2022.3.0.8767/mpi/2021.7.0/bin/mpiicc
export CXX=/ocean/projects/phy230038p/awu2/2022.3.0.8767/mpi/2021.7.0/bin/mpiicpc

#export COMPILER_ID=GNU
#export FC=/ocean/projects/phy230038p/awu2/gcc_install/GCC-12.2.0/bin/gfortran
#export CC=/ocean/projects/phy230038p/awu2/gcc_install/GCC-12.2.0/bin/gcc
#export CXX=/ocean/projects/phy230038p/awu2/gcc_install/GCC-12.2.0/bin/g++

export FFTW_PATH=${CWD}/dependencies/fftw-3.3.5
export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
export HDF5_PATH=${CWD}/dependencies/hdf5-1.8.18
export FFTPACK_PATH=${CWD}/dependencies/fftpack
#export ARCH_OPT_FLAG="-xCORE-AVX2 -axMIC-AVX512"
export ARCH_OPT_FLAG="-xHost"
