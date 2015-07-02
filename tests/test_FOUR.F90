#include "../src/utilities/kind_parameters.F90"
#include "../src/utilities/constants.F90"
#include "../src/utilities/ffts.F90"
#include "../src/utilities/timer.F90"

module test_module
contains
#include "../src/derivatives/CD10.F90"
end module

program test_FOUR

    use kind_parameters, only: rkind
    use fftstuff, only: ffts
    use constants, only: one,two,pi,imi,half
    use test_module
    use timer, only: tic, toc
    implicit none

    type(ffts) :: my1dFFT
    integer, parameter :: nx=2048,ny=150,nz=150
    real(rkind), dimension(nx,9) :: LU
    real(rkind), dimension(nx) :: tmp
    real(rkind) :: d=1.0_rkind, a=0.5_rkind, c=0.5_rkind, e=0.05_rkind, ff=0.05_rkind
    real(rkind), dimension(nx,ny,nz) :: x,f,kx,fp,df, df_CD10
    real(rkind) :: dx, onebydx, sigma

    integer :: i,j,ierr

    dx = two*pi/real(nx,rkind)
    onebydx = one / dx
    
    ! Compute LU for CD10
    call ComputeLU10(LU,nx,e,a,d,c,ff)
   
    ! Initialize f as a periodic gaussian
    sigma = half
    do i=1,nx
        x(i,:,:) = real((i-1),rkind)*dx
        f(i,:,:) = exp( -(x(i,:,:) - pi)**2/(2*sigma**2) )
    end do

    ! Exact derivative
    fp = -(one/(2*sigma**2))*two*(x - pi)*f
    print*, f(nx,1,1),fp(nx,1,1)

    ! Initialize fft
    ierr = my1dFFT % init(nx,dx)

    ! Get fft derivative
    call tic()
    do i=1,ny
      do j=1,nz
        df(:,i,j) = my1dFFT % ifft( imi * my1dFFT % k() *  my1dFFT % fft(f(:,i,j)) )
      end do
    end do
    call toc()
    
    ! Get CD10 derivative
    call tic()
    do i=1,ny
      do j=1,nz
        !tmp = CD10D1RHS( f(:,i,j),nx,.TRUE. )
        df_CD10(:,i,j) = CD10D1RHS( f(:,i,j),nx,.TRUE. ) * onebydx
        call SolveLU10(LU,df_CD10(:,i,j),nx)
        df_CD10(:,i,j) = df_CD10(:,i,j) ! * onebydx
      end do
    end do
    call toc()

    print*, "Maximum error (FOUR) = ", MAXVAL(ABS(fp - df))
    print*, "Maximum error (CD10) = ", MAXVAL(ABS(fp - df_CD10))
end program
