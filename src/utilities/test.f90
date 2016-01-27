program test
    implicit none
    real(kind=8) :: A(20,3)

    A = 2.d0
    call dynamic_inout(A)



end program test



subroutine dynamic_inout(A)
    implicit none
    real(kind=8), intent(in) :: A(:,:)
    integer :: n

    n = size(A,1)
    print*, n

    !A(1,1) = 1.d0



end subroutine dynamic_inout
