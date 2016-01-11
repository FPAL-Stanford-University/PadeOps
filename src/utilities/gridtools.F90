module gridtools
    use kind_parameters, only: rkind
    use constants, only: one, ten
    use decomp_2d, only: decomp_info
    use exits,     only: GracefulExit
    implicit none

    interface destroy_buffs
        module procedure destroy_buffs_real, destroy_buffs_complex 
    end interface
    
    interface alloc_buffs
        module procedure alloc_buffs_real, alloc_buffs_complex 
    end interface

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

    subroutine destroy_buffs_real(buff)
        real(rkind), dimension(:,:,:,:),allocatable, intent(inout) :: buff

        if (allocated(buff)) deallocate(buff)
    end subroutine

    subroutine destroy_buffs_complex(buff)
        complex(rkind), dimension(:,:,:,:),allocatable, intent(inout) :: buff

        if (allocated(buff)) deallocate(buff)
    end subroutine

    subroutine alloc_buffs_real(buff,vars,dir,decomp)
        character(len=1), intent(in) :: dir
        class(decomp_info), intent(in) :: decomp
        integer, intent(in) :: vars
        real(rkind), dimension(:,:,:,:), allocatable, intent(out) :: buff

        if (allocated(buff)) deallocate(buff)

        select case (dir)
        case ("x")
            allocate(buff(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),vars))
        case("y")
            allocate(buff(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),vars))
        case("z")
            allocate(buff(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3),vars))
        case default
            call GracefulExit("Incorrect direction selected in ALLOC_BUFFS_REAL subroutine", 13123)
        end select

    end subroutine

    subroutine alloc_buffs_complex(buff,vars,dir,decomp)
        character(len=1), intent(in) :: dir
        class(decomp_info), intent(in) :: decomp
        integer, intent(in) :: vars
        complex(rkind), dimension(:,:,:,:), allocatable, intent(out) :: buff

        if (allocated(buff)) deallocate(buff)

        select case (dir)
        case ("x")
            allocate(buff(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),vars))
        case("y")
            allocate(buff(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),vars))
        case("z")
            allocate(buff(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3),vars))
        case default
            call GracefulExit("Incorrect direction selected in ALLOC_BUFFS_REAL subroutine", 13123)
        end select

    end subroutine

end module 
