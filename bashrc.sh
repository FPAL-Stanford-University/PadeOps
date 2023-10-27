# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User Specific Aliases and Functions
alias cdp="cd ~/Codes/PadeOps"

###Needed for PadeOps
module load ifort/2018 icc/2018 impi/2018 imkl/2018 cmake

#visit
ml apps/visit/3.1.4
# Modules to load 
if [[ "$SHERLOCK" == "1" ]]; then
    echo "In Sherlock 1"
elif [[ "$SHERLOCK" == "2" ]]; then
    echo "In Sherlock 2"
else
    echo "Uh-oh, not sure where we are..." 
fi

#Load modules for postprocessing
module load viz
module load py-numpy
module load py-matplotlib

function setup_padeops {
	CWD=$HOME/Codes/PadeOps

	export COMPILER_ID=Intel
	export FC=mpiifort
	export CC=mpiicc
	export CXX=mpiicpc
	export I_MPI_CC=icc
	export I_MPI_CXX=icpc
	export I_MPI_F90=ifort
	export I_MPI_F77=ifort
	export FFTW_PATH=${CWD}/dependencies/fftw-3.3.5
	export DECOMP_PATH=${CWD}/dependencies/2decomp_fft
	export VTK_IO_PATH=${CWD}/dependencies/Lib_VTK_IO/build
	export HDF5_PATH=${CWD}/dependencies/hdf5-1.8.18
	export OMP_PROC_BIND=spread
	export OMP_PLACES=threads

}
