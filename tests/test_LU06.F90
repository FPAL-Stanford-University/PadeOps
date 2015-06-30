#include '../src/utilities/kind_parameters.F90'
#include '../src/derivatives/CD06.F90'

program test_LU06

    use kind_parameters, only: rkind
    implicit none

    integer, parameter :: N=1000000
    integer, parameter :: nvec=4
    real(rkind), dimension(N,5) :: LU
    real(rkind), dimension(N,nvec) :: k
    real(rkind) :: d=1.0_rkind, a=1.0_rkind/3.0_rkind, b=1.0_rkind/3.0_rkind
    integer :: i,j

    call ComputeLU06(LU,N,b,d,a)

    do i=1,N
        k(i,:) = real(i,rkind)
    end do

    do i=1,10
      do j=1,nvec
        call SolveLU06(LU,k(:,i),N,1)
      end do
    end do

    ! do i = 1,N
    !     print ('(9F10.4)'), LU(i,:), k(i,:)
    ! end do
    
end program

