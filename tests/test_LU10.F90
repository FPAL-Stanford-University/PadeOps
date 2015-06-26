include '../src/utilities/kind_parameters.F90'
include '../src/derivatives/CD10.F90'

program test_LU10

    use kind_parameters, only: rkind
    implicit none

    integer, parameter :: N=100000
    real(rkind), dimension(N,9) :: LU
    real(rkind) :: d=1.0_rkind, a=0.5_rkind, c=0.5_rkind, e=0.05_rkind, f=0.05_rkind
    integer :: i

    call ComputeLU10(LU,N,e,a,d,c,f)

    ! do i = 1,N
    !     print ('(9F8.4)'), LU(i,:)
    ! end do
    

end program

