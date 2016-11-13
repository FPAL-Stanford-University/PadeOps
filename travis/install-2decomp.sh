#!/bin/sh
# This configuration file was taken originally from the mpi4py project
# <http://mpi4py.scipy.org/>, and then modified for Julia

set -e
set -x

os=`uname`

case "$os" in
    Linux)
        cd dependencies
        if [ ! -d "2decomp_fft" ]; then
            tar -zxf 2decomp_fft-1.5.847.tar.gz
        fi
        cd 2decomp_fft
        cp ../../travis/2decomp_fft_Makefile.inc src/Makefile.inc
        make > /dev/null
        ;;

    *)
        echo "Unknown operating system: $os"
        exit 1
        ;;
esac
