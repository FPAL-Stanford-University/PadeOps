module purge 
module load cmake/3.20.2 gcc/9.1.0 mvapich2/2.3.7
export COMPILER_ID=GNU
export FC="mpif90 -L/work2/06632/ryanhass/stampede2/software/PadeOps_dependencies/GNUcompiled/OpenBLAS-0.3.21"
export CC=mpicc
export CXX=mpicxx
export FFTW_PATH=/work2/06632/ryanhass/stampede2/software/PadeOps_dependencies/GNUcompiled/fftw-3.3.10
export HDF5_PATH=/work2/06632/ryanhass/stampede2/software/PadeOps_dependencies/GNUcompiled/hdf5-1.8.18
export FFTPACK_PATH=/work2/06632/ryanhass/stampede2/software/PadeOps_dependencies/GNUcompiled/fftpack
export VTK_IO_PATH=/work2/06632/ryanhass/stampede2/software/PadeOps_dependencies/GNUcompiled/Lib_VTK_IO/build
export DECOMP_PATH=/work2/06632/ryanhass/stampede2/software/PadeOps_dependencies/GNUcompiled/2decomp_fft
