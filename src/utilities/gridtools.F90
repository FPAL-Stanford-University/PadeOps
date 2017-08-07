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
    
   pure function mytrapz(x,f) result(intf)
        real(rkind), dimension(:), intent(in) :: x, f
        real(rkind) :: intf
        real(rkind) :: dx

        dx = x(2) - x(1)
        intf = sum(f)*dx
    end function 

    pure function logspace(x_min,x_max,n) result(x)
        integer, intent(in) :: n                    ! Desired size
        real(rkind), dimension(n) :: x ! Output array
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

    end function

    pure function linspace(x_min,x_max,n) result(x)
        integer, intent(in) :: n
        real(rkind), dimension(n) :: x
        real(rkind), intent(in) :: x_min, x_max
        
        real(rkind) :: step
        integer :: i

        if (n .le. 1) then
            x = x_min
            return
        end if 
        step = (x_max - x_min)/(real(n,rkind) - one)
        
        x(1) = x_min
        do i = 2,n
            x(i) = x(i-1) + step
        end do 

    end function
    
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
        type(decomp_info), intent(in) :: decomp
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
        type(decomp_info), intent(in) :: decomp
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

    pure subroutine upsample_Periodic_1d(vecin, vecout, nx)
        integer, intent(in) :: nx
        real(rkind), intent(in) , dimension(nx) :: vecin
        real(rkind), intent(out), dimension(2*nx) :: vecout
        integer :: i

        do i = 1,nx
            vecout(2*i - 1) = vecin(i)
        end do 

        do i = 2,2*nx-2,2
            vecout(i) = 0.5d0*(vecout(i-1) + vecout(i+1))
        end do 
        vecout(2*nx) = 0.5d0*(vecout(2*nx-1) + vecout(1))
    end subroutine
end module 
