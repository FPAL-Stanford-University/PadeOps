# README #

This   README would normally document whatever steps are necessary to get your application up and running.

### PadeOps ###

* Quick summary : Hybrid OpenMP/MPI derivative operators using Compact Difference (6th and 10th Order) and Spectral (Fourier and Chebyshev) Methods to solve PDEs.
* Version : 0.1
* Tutorials : TBD/Incomplete

### How do I get set up? ###

* Summary of set up
* Configuration :
    * The Fortran compiler (FC) (preferably even the CC and CXX variables) need to be set to the desired MPI fortran compiler. Also the FFTW library path (FFTW_PATH) and 2DECOMP&FFT library path (DECOMP_PATH) need to be set.
    * Examples of this are in the setup folder. For your system, either use one of the SetupEnv_<Machine>\_<CompilerID>.sh files or copy the closest one to your own SetupEnv_<Machine>\_<CompilerID>.sh. Then, from the main directory (PadeOps), run `source setup/SetupEnv_<Machine>_<CompilerID>.sh` to set the correct environment variables.

    * Now, build the required dependencies:
         * FFTW
             * Extract the fftw-<version>.tar.gz file in the dependencies folder. Change directory using `cd fftw-<version>`. Set the envireonment variables `F77` and `MPICC` to the correct fortran and MPI C compilers. Configure the build using `./configure --prefix=<current directory> --enable-avx`. Then build the library using `make; make install`. At the end of this, there should be a folder in this directory called `lib` and a folder called `include` with the static library file and the include files respectively.
         * 2DECOMP&FFT
             * Extract the 2decomp_fft.tar.gz file in the dependencies folder. Change directory using `2decomp_fft/src`. Now, in the file Makefile.inc.x86, change the FFT=<MKL> (line 25) to FFT=fftw_f03. In line 32, change the FFTW_PATH variable to where you installed FFTW in the previous step. In line 57, set the F90 variable to your MPI Fortran compiler (like mpif90). In line 71, change the CC variable to your MPI C compiler. This next step is optional and might be required on some systems. Remove the -lfftw3f flag in line 77 if required for compilation. Make a symbolic link to `Makefile.inc` using `ln -s Makefile.inc.x86 Makefile.inc`. Move to the directory above using `cd ..`. Then build the library using `make`. At the end of this, there should be a folder in this directory called `lib` and a folder called `include` with the static library file and the include files respectively. Now set the DECOMP_PATH variable to the current directory in your SetupEnv.sh script.
         * Lib_VTK_IO
              * Extract the Lib_VTK_IO.tar.gz file and `cd` into the created directory. Make a build directory and move to it using `mkdir build; cd build`. Then build the library using `cmake ..; make`. Now, set the VTK_IO_PATH to the current directory in your SetupEnv.sh script. 
          * HDF5
              * Extract the hdf5-1.8.18.tar.gz file and `cd` into the directory created. Configure the build using `CC=</path/to/mpicc> FC=</path/to/mpif90> CXX=</path/to/mpic++> ./configure --enable-parallel --enable-fortran --enable-production --prefix=<current directory>`. Build HDF5 using `make; make install`. Set the `HDF5_PATH` variable in the `SetupEnv_<MACHINE>_<COMPILER>.sh` script to the directory that you built HDF5 in.

    * To build the code, run the following commands:
~~~
           mkdir build
           cd build
           cmake ..
           make
~~~
* Dependencies :
    * CMake 2.8 or above
    * MPI Library (https://www.mpich.org/)
    * FFTW (www.fftw.org/)
    * 2DECOMP&FFT (http://www.2decomp.org/)
    * Parallel HDF5 (https://www.hdfgroup.org/)
* Database configuration : TBD
* How to run tests : TBD
* Deployment instructions : TBD

### Contribution guidelines ###

* To merge a branch "dev" to master, use the following commands:
```bash
git merge --no-commit dev
git checkout .travis.yml
git commit -m "merge dev into master"
```

### Who do I talk to? ###

* Akshay Subramaniam or Aditya Ghate
