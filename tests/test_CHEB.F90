#include "../src/utilities/kind_parameters.F90"
#include "../src/utilities/constants.F90"
#include "../src/utilities/dcts.F90"
#include "../src/utilities/ffts.F90"
#include "../src/utilities/timer.F90"

module test_module
contains
#include "../src/derivatives/CD10.F90"
end module

program test_FOUR

    use kind_parameters, only: rkind
    use dctstuff, only: dcts
    use fftstuff, only: ffts
    use constants, only: one,two,pi,imi,half
    use test_module
    use timer, only: tic, toc
    implicit none

    type(dcts) :: my1dDCT
    type(ffts) :: my1dFFT
    integer, parameter :: nx=2049,ny=50,nz=50
    real(rkind), dimension(nx,9) :: LU
    real(rkind), dimension(nx) :: tmp
    real(rkind) :: d=1.0_rkind, a=0.5_rkind, c=0.5_rkind, e=0.05_rkind, ff=0.05_rkind
    real(rkind), dimension(nx,ny,nz) :: x,f,kx,fp,df, df_CD10
    complex(rkind), dimension(nx,ny,nz) :: fpp
    real(rkind) :: onebydx, sigma

    integer :: i,j,ierr
   
    ! Initialize f as a periodic gaussian
    do i=1,nx
        x(i,:,:) = - cos( real((i-1),rkind) * pi / real((nx-1),rkind) )
        f(i,:,:) = x(i,:,:) 
    end do

    ! Initialize dct
    ierr = my1dDCT % init(nx)
    ierr = my1dFFT % init(nx-1,one)

    call tic()
    do i=1,ny
      do j=1,nz
        fp(:,i,j) = my1dDCT % dct(f(:,i,j))
        df(:,i,j) = my1dDCT % idct(fp(:,i,j))
      end do
    end do
    call toc()
    print*, "Maximum error = ", MAXVAL(ABS(df - f))

  !  print*, fp(:,1,1)
    call tic() 
    do i=1,ny
      do j=1,nz
        fpp(1:(nx-1)/2 + 1,i,j) = my1dFFT % fft(f(:,i,j), .TRUE.)
        df(:,i,j) = my1dFFT % ifft(fpp(:,i,j))
      end do
    end do
    call toc()

    print*, "Maximum error = ", MAXVAL(ABS(df(1:nx-1,:,:) - f(1:nx-1,:,:)))

end program
