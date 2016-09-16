#!/bin/sh
# This configuration file was taken originally from the mpi4py project
# <http://mpi4py.scipy.org/>, and then modified for Julia

set -e
set -x

os=`uname`

case "$os" in
    Linux)
        cd dependencies
        if [[ ! -d "Lib_VTK_IO" ]]; then
            tar -zxf Lib_VTK_IO.tar.gz
        fi
        cd Lib_VTK_IO
        if [[ ! -d "build" ]]; then
            mkdir build
        fi
        cd build
        export CC=mpicc
        export CXX=mpicxx
        export FC=mpif90
        cmake ..
        make -j > /dev/null
        ;;

    *)
        echo "Unknown operating system: $os"
        exit 1
        ;;
esac
