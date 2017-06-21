#!/bin/sh
# This configuration file was taken originally from the mpi4py project
# <http://mpi4py.scipy.org/>, and then modified for Julia

set -e
set -x

os=`uname`

case "$os" in
    Linux)
        wget --no-check-certificate https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.0-patch1/src/hdf5-1.10.0-patch1.tar
        if [ ! -d "hdf5-1.10.0-patch1" ]; then
            tar -xf hdf5-1.10.0-patch1.tar
        fi
        cd hdf5-1.10.0-patch1
        FC=mpif90 CC=mpicc CXX=mpic++ RUNPARALLEL='mpiexec -n $${NPROCS:=4}' sh ./configure --enable-fortran --enable-parallel --prefix=$HOME/HDF5 > /dev/null
        make -j > /dev/null
        make install > /dev/null
        ;;

    *)
        echo "Unknown operating system: $os"
        exit 1
        ;;
esac
