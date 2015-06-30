#include "../src/utilities/kind_parameters.F90"
#include "../src/derivatives/CD10.F90"
#include "../src/utilities/timer.F90"

program test_LU10

    use kind_parameters, only: rkind
    use timer, only: tic, toc
    implicit none

    integer, parameter :: ny=8192,nz=100,nx=100
    real(rkind), dimension(ny,9) :: LU
    real(rkind), dimension(nx,ny,nz) :: y
    real(rkind), dimension(ny) :: tmp
    real(rkind) :: d=1.0_rkind, a=0.5_rkind, c=0.5_rkind, e=0.05_rkind, f=0.05_rkind
    integer :: i,j
    double precision :: t0,t1

    call tic()
    call ComputeLU10(LU,ny,e,a,d,c,f)
    call toc()

    call tic()
    do i = 1,nz
        y(:,i,:) = real(i,rkind)
    end do
    call toc()

    !allocate( tmp(nx) ) 

    call tic()
    do i=1,nx
      do j=1,nz
        tmp = y(i,:,j)
        call SolveLU10(LU,tmp,ny)
      end do
    end do
    call toc()

    ! do i = 1,nx
    !     print ('(13F8.4)'), LU(i,:), y(i,:,:)
    ! end do
    !deallocate( tmp ) 

end program

