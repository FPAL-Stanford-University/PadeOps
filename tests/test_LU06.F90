include '../src/utilities/kind_parameters.F90'
include '../src/derivatives/CD06.F90'

program test_LU06

    use kind_parameters, only: rkind
    implicit none

    integer, parameter :: N=1000000
    real(rkind), dimension(N,5) :: LU
    real(rkind) :: d=1.0_rkind, a=1.0_rkind/3.0_rkind, b=1.0_rkind/3.0_rkind
    integer :: i

    call ComputeLU06(LU,N,b,d,a)

    ! do i = 1,N
    !     print ('(5F8.4)'), LU(i,:)
    ! end do
    

end program

