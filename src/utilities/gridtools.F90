module gridtools
    use kind_parameters, only: rkind
    use constants, only: one, ten 
    implicit none

contains

    pure subroutine logspace(x,x_min,x_max,n)
        integer, intent(in) :: n                    ! Desired size
        real(rkind), dimension(n), intent(out) :: x ! Output array
        real(rkind), intent(in) :: x_min, x_max     ! Left and Right bounds in POWERS OF 10

        real(rkind) :: step
        integer :: i
        real(rkind), dimension(n) :: xtemp

        if (n .le. 1) then
            x = ten**(x_min)
            return 
        end if 
        step = (x_max - x_min)/(real(n,rkind)- one)
        xtemp(1) = x_min

        do i = 2,n
            xtemp(i) = xtemp(i-1) + step
        end do

         x = (ten)**(xtemp)

    end subroutine 



end module 
