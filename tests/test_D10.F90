#include "../src/utilities/timer.F90"
#include "../src/utilities/kind_parameters.F90"

module test_module
contains
#include "../src/derivatives/CD10.F90"
end module

program test_LU10

    use kind_parameters, only: rkind
    use timer, only: tic, toc
    use test_module
    implicit none

    integer, parameter :: nx=256,ny=100,nz=100
    real(rkind), dimension(nx,9) :: LU
    real(rkind), dimension(nx,ny,nz) :: y,yp,dy
    real(rkind), dimension(nx) :: tmp
    real(rkind) :: dx
    real(rkind) :: d=1.0_rkind, a=0.5_rkind, c=0.5_rkind, e=0.05_rkind, f=0.05_rkind
    integer :: i,j
    double precision :: t0,t1

    dx = 8._rkind*atan(1._rkind)/real(nx,rkind)

    call tic()
    call ComputeLU10(LU,nx,e,a,d,c,f)
    call toc()

    call tic()
    do i = 1,nx
        y(i,:,:) = sin( 10._rkind*real((i-1),rkind)*dx)
        yp(i,:,:) = 10._rkind*cos( 10._rkind*real((i-1),rkind)*dx)
    end do
    call toc()

    !print*, "nx = ", nx

    call tic()
    do i=1,ny
      do j=1,nz
        tmp = CD10D1RHS(y(:,i,j),nx,.TRUE.)
        call SolveLU10(LU,tmp,nx)
        dy(:,i,j) = tmp / dx
      end do
    end do
    call toc()

    yp = dy - yp

    ! do i = 1,nx
    !     print ('(1F8.4)'), y(i,1,1)
    ! end do
    
    !print*, "Maximum error = ", MAXVAL(ABS(yp))

end program

