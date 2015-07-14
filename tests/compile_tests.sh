# Compile test_cd10.F90
ifort -warn all -O3 -qopenmp -ipo -mcmodel="medium" -i8 -qopt-report=2 -qopt-report-phase=vec ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/timer.F90 ../src/utilities/cd10.F90 test_cd10.F90 -o test_cd10

# Compile test_cd06.F90
# ifort -warn all ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/cd06.F90 test_cd06.F90 -o test_cd06

# Compile test_ffts.F90
# ifort -warn all -g -pg -O3 -mkl -mcmodel="medium" -i8 -ipo -qopt-report=2 -qopt-report-phase=vec ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/timer.F90 ../src/utilities/ffts.F90 test_ffts.F90 -o test_ffts

# Compile test_dcts.F90
# ifort -warn all -mkl ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/dcts.F90 test_dcts.F90 -o test_dcts

# Compile test_derivatives.F90
#ifort -warn all -g -O3 -mkl -ipo -mcmodel="medium" -i8 -qopt-report=2 -qopt-report-phase=vec ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/timer.F90 ../src/utilities/exits.F90 ../src/utilities/cd10.F90 ../src/utilities/cd06.F90 ../src/utilities/ffts.F90 ../src/utilities/dcts.F90 ../src/derivatives/Derivatives.F90 test_derivatives.F90 -o test_derivatives
#gfortran -g -O0 -pg -fdefault-integer-8 -fbounds-check -fno-automatic -Wall -Wextra -Wconversion -fbacktrace -I/afs/ir/users/a/k/akshays/FFTW_INSTALL/include -L/afs/ir/users/a/k/akshays/FFTW_INSTALL/lib ../src/utilities/kind_parameters.F90 ../src/utilities/constants.F90 ../src/utilities/timer.F90 ../src/utilities/exits.F90 ../src/utilities/cd10.F90 ../src/utilities/cd06.F90 ../src/utilities/ffts.F90 ../src/utilities/dcts.F90 ../src/derivatives/Derivatives.F90 test_derivatives.F90 -o test_derivatives -lfftw3
