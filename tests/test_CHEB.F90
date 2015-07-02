#include "../src/utilities/kind_parameters.F90"
#include "../src/utilities/constants.F90"
#include "../src/utilities/dcts.F90"
#include "../src/utilities/timer.F90"

module test_module
contains
#include "../src/derivatives/CD10.F90"
end module

program test_FOUR

    use kind_parameters, only: rkind
    use dctstuff, only: dcts
    use constants, only: one,two,pi,imi,half
    use test_module
    use timer, only: tic, toc
    implicit none

    type(dcts) :: my1dDCT
    integer, parameter :: nx=17,ny=1,nz=1
    real(rkind), dimension(nx,9) :: LU
    real(rkind), dimension(nx) :: tmp
    real(rkind) :: d=1.0_rkind, a=0.5_rkind, c=0.5_rkind, e=0.05_rkind, ff=0.05_rkind
    real(rkind), dimension(nx) :: x,f,kx,fp,df, df_CD10
    real(rkind) :: onebydx, sigma

    integer :: i,j,ierr

    
   
    ! Initialize f as a periodic gaussian
    do i=1,nx
        x(i) = real((i-1),rkind)
        f(i) = real(i,rkind) 
    end do

    ! Exact derivative
    ! fp = 

    ! Initialize fft
    ierr = my1dDCT % init(nx)


    fp = my1dDCT % dct(f)
    print*, fp
    print*, my1dDCT % idct( fp )
end program
