#!/bin/sh
# This configuration file was taken originally from the mpi4py project
# <http://mpi4py.scipy.org/>, and then modified for Julia

set -e
set -x

os=`uname`

case "$os" in
    Linux)
        cd dependencies
        if [ ! -d "fftw-3.3.4" ]; then
            tar -zxf fftw-3.3.4.tar.gz
        fi
        cd fftw-3.3.4
        sh ./configure --prefix=$HOME/FFTW F77=gfortran-6 MPICC=mpicc > /dev/null
        make -j > /dev/null
        sudo make install > /dev/null
        ;;

    *)
        echo "Unknown operating system: $os"
        exit 1
        ;;
esac
