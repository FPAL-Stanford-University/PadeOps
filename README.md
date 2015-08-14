# README #

This README would normally document whatever steps are necessary to get your application up and running.

### PadeOps ###

* Quick summary : Hybrid OpenMP/MPI derivative operators using Compact Difference (6th and 10th Order) and Spectral (Fourier and Chebyshev) Methods to solve PDEs.
* Version : 0.1
* Tutorials : TBD/Incomplete

### How do I get set up? ###

* Summary of set up
* Configuration :
    * The Fortran compiler (FC) (preferably even the CC and CXX variables) need to be set to the desired MPI fortran compiler. Also the FFTW library path (FFTW_PATH) and 2DECOMP&FFT library path (DECOMP_PATH) need to be set.
    * Examples of this are in the setup folder. For your system, either use one of the SetupEnv_*.sh files or copy the closest one to SetupEnv_<Machine>_<CompilerID>.sh. Then, from the main directory (PadeOps), run
~~~
           source setup/SetupEnv_<Machine>_<CompilerID>.sh
~~~
      to set the correct environment variables.
    * To build the code, run the following commands:
~~~
           mkdir build
           cd build
           cmake ..
           make
~~~
* Dependencies :
    * MPI Library (https://www.mpich.org/)
    * FFTW (www.fftw.org/)
    * 2DECOMP&FFT (http://www.2decomp.org/)
* Database configuration : TBD
* How to run tests : TBD
* Deployment instructions : TBD

### Contribution guidelines ###

* TBD

### Who do I talk to? ###

* Akshay Subramaniam or Aditya Ghate