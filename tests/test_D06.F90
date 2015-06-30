#include "../src/utilities/timer.F90"
#include "../src/utilities/kind_parameters.F90"

module test_module
contains
#include "../src/derivatives/CD06.F90"
end module

program test_LU06

    use kind_parameters, only: rkind
    use timer, only: tic, toc
    use test_module
    implicit none

    integer, parameter :: nx=2048,ny=100,nz=100
    real(rkind), dimension(nx,5) :: LU
    real(rkind), dimension(nx,ny,nz) :: y,yp,dy
    real(rkind), dimension(nx) :: tmp
    real(rkind) :: dx
    real(rkind) :: d=1.0_rkind, a=1.0_rkind/3.0_rkind, b=1.0_rkind/3.0_rkind
    integer :: i,j
    double precision :: t0,t1

    dx = 8._rkind*atan(1._rkind)/real(nx,rkind)

    call tic()
    call ComputeLU06(LU,nx,b,d,a)
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
        tmp = CD06D1RHS(y(:,i,j),nx,.TRUE.)
        call SolveLU06(LU,tmp,nx)
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

