#include "../src/utilities/kind_parameters.F90"
#include "../src/utilities/constants.F90"
#include "../src/utilities/ffts.F90"

module test_module
contains
#include "../src/derivatives/CD10.F90"
end module

program test_FOUR

    use kind_parameters, only: rkind
    use ffts, only: FFTW
    use constants, only: one,two,pi,imi,half
    use test_module
    implicit none

    type(FFTW) :: my1dFFT
    integer, parameter :: nx=1024
    real(rkind), dimension(nx,9) :: LU
    real(rkind), dimension(nx) :: tmp
    real(rkind) :: d=1.0_rkind, a=0.5_rkind, c=0.5_rkind, e=0.05_rkind, ff=0.05_rkind
    real(rkind), dimension(nx) :: x,f,kx,fp,df, df_CD10
    real(rkind) :: dx, sigma

    integer :: i,ierr

    dx = two*pi/real(nx,rkind)
    
    call ComputeLU10(LU,nx,e,a,d,c,ff)
    sigma = half
    do i=1,nx
        x(i) = real((i-1),rkind)*dx
        f(i) = exp( -(x(i) - pi)**2/(2*sigma**2) )
    end do

    fp = -(one/(2*sigma**2))*two*(x - pi)*f
    print*, f(nx),fp(nx)

    ierr = my1dFFT % init(nx,dx)

    df = my1dFFT % ifft( imi * my1dFFT % k() *  my1dFFT % fft(f) )
    tmp = CD10D1RHS(f,nx,.TRUE.)
    call SolveLU10(LU,tmp,nx)
    df_CD10 = tmp / dx
    print*, "Maximum error (FOUR) = ", MAXVAL(ABS(fp - df))
    print*, "Maximum error (CD10) = ", MAXVAL(ABS(fp - df_CD10))
end program
