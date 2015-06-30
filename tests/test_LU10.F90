#include '../src/utilities/kind_parameters.F90'
#include '../src/derivatives/CD10.F90'

program test_LU10

    use kind_parameters, only: rkind
    implicit none

    integer, parameter :: N=1000000
    integer, parameter :: nvec=4
    real(rkind), dimension(N,9) :: LU
    real(rkind), dimension(N,nvec)   :: y
    real(rkind) :: d=1.0_rkind, a=0.5_rkind, c=0.5_rkind, e=0.05_rkind, f=0.05_rkind
    integer :: i,j

    call ComputeLU10(LU,N,e,a,d,c,f)

    do i = 1,N
        y(i,:) = real(i,rkind)
    end do

    do i=1,10
      !do j=1,nvec
        call SolveLU10(LU,y,N,nvec)
      !end do
    end do

    ! do i = 1,N
    !     print ('(13F8.4)'), LU(i,:), y(i,:)
    ! end do

end program

