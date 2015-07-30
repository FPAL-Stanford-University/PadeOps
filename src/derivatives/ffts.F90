module fftstuff
    use kind_parameters, only: rkind
    use constants, only: zero,one,half,two,pi,imi

    implicit none
    private
    public :: ffts

    include "fftw3.f"


    type ffts
        private 

        integer :: n
        real(rkind) :: onebyn

        real(rkind) :: dx
        integer(kind=8) :: plan_fwd
        integer(kind=8) :: plan_bwd

        real(rkind), dimension(:,:,:), allocatable :: k1d

        integer :: split

        character(len=1) :: dir = "x"

        integer :: n2, n3

        logical :: initialized = .false.
        
        integer :: fftw_plan = FFTW_MEASURE

        contains

            procedure :: init
            procedure :: destroy
            procedure :: dd1
            procedure :: dd2
            procedure :: dd3

            procedure  :: fftx
            procedure  :: ifftx
            procedure  :: ffty
            procedure  :: iffty
            procedure  :: fftz
            procedure  :: ifftz


    end type


contains

    function init(this, n_, dir_, n2_, n3_, dx_, exhaustive) result(ierr)
        class(ffts), intent(inout) :: this
        integer, intent(in) :: n_, n2_, n3_
        character(len=1), intent(in) :: dir_
        integer :: ierr
        logical, optional :: exhaustive
        real(rkind), intent(in) :: dx_
        real(rkind), dimension(:,:,:), allocatable :: in_arr
        complex(rkind), dimension(:,:,:), allocatable :: out_arr
        real(rkind), dimension(:,:), allocatable :: in_arr_2d
        complex(rkind), dimension(:,:), allocatable :: out_arr_2d

        real(rkind), dimension(:), allocatable :: k_tmp
        integer :: i, j, k
        

        this%dir = dir_
        this%n = n_
        this%n2 = n2_
        this%n3 = n3_
        this%dx = dx_
        this%onebyn = one/real(this%n,rkind)

        this%split = n_/2 + 1
        
        if (present(exhaustive)) then
            if (exhaustive) this%fftw_plan = FFTW_EXHAUSTIVE
        end if 

        if (allocated(this%k1d)) deallocate(this%k1d)
        allocate(k_tmp(this%split))
        k_tmp = GetWaveNums(this%n,this%dx)
        
        
        select case (this%dir)
        case ("x")
            allocate (this%k1d(this%split,this%n2,this%n3))
            do k = 1,this%n3
                do j = 1,this%n2
                    this%k1d(:,j,k) = k_tmp
                end do 
            end do
            allocate (in_arr (this%n    , this%n2, this%n3))
            allocate (out_arr(this%split, this%n2, this%n3))

            call dfftw_plan_many_dft_r2c(this%plan_fwd, 1, this%n, &
                this%n2*this%n3, in_arr, this%n, 1, &
                this%n, out_arr, this%split, 1, this%split, &
                this%fftw_plan)

            call dfftw_plan_many_dft_c2r(this%plan_bwd, 1, this%n, &
                this%n2*this%n3, out_arr, this%split, 1, &
                this%split, in_arr, this%n, 1, this%n, &
                this%fftw_plan)

             deallocate (in_arr, out_arr)

        case ("y")
            allocate (this%k1d(this%n2,this%split,this%n3))
            do k = 1,this%n3
                do i = 1,this%n2
                    this%k1d(i,:,k) = k_tmp
                end do 
            end do
    
            allocate (in_arr_2d(this%n2, this%n))
            allocate (out_arr_2d(this%n2, this%split))

            call dfftw_plan_many_dft_r2c(this%plan_fwd, 1, this%n, &
                this%n2, in_arr_2d, this%n, this%n2, 1, out_arr_2d, this%split, &
                this%n2, 1, this%fftw_plan)

            call dfftw_plan_many_dft_c2r(this%plan_bwd, 1, this%n, &
                this%n2, out_arr_2d, this%split, this%n2, 1, in_arr_2d, this%n, &
                this%n2, 1, this%fftw_plan)

            deallocate (in_arr_2d, out_arr_2d) 
        
        case ("z")
            allocate (this%k1d(this%n2,this%n3,this%n))
            do i = 1,this%n2
                do j = 1,this%n3
                    this%k1d(i,j,:) = k_tmp
                end do 
            end do
            allocate (in_arr (this%n2, this%n3, this%n ))
            allocate (out_arr(this%n2, this%n3, this%split))

            call dfftw_plan_many_dft_r2c(this%plan_fwd, 1, this%n, &
                this%n2*this%n3, in_arr, this%n, this%n2*this%n3, &
                1, out_arr, this%split, this%n2*this%n3, 1, &
                this%fftw_plan)

            call dfftw_plan_many_dft_c2r(this%plan_bwd, 1, this%n, &
                this%n2*this%n3, out_arr, this%split, this%n2*this%n3, &
                1, out_arr, this%n, this%n2*this%n3, 1, &
                this%fftw_plan)
             
            deallocate (in_arr, out_arr)
        case default 
           ierr = 20
            return  
        end select 

        this%initialized = .true. 
        ierr = 0
    end function 

    subroutine dd1(this,f, df)
        class(ffts), intent(in) :: this
        real(rkind), dimension(this%n,this%n2, this%n3), intent(in) :: f 
        real(rkind), dimension(this%n,this%n2, this%n3), intent(out) :: df 
        complex(rkind), dimension(this%split,this%n2,this%n3) :: f_hat

        if ((this%dir .eq. "x") .and. (this%n > 1)) then
           call this%fftx(f,f_hat)

           f_hat = imi*this%k1d*f_hat
           f_hat(this%split,:,:) = zero

           call this%ifftx(f_hat,df)
        else 
            df = zero
        end if 
    end subroutine 


    subroutine dd2(this,f, df)
        class(ffts), intent(in) :: this
        real(rkind), dimension(this%n2,this%n, this%n3), intent(in) :: f 
        real(rkind), dimension(this%n2,this%n, this%n3), intent(out) :: df 
        complex(rkind), dimension(this%n2,this%split) :: f_hat
        integer :: k

        if (( this%dir .eq. "y") .and. (this%n > 1)) then
            do k = 1,this%n3
                call this%ffty(f(:,:,k),f_hat)
                f_hat = imi*this%k1d(:,:,k)*f_hat
                f_hat(:,this%split) = zero
                call this%iffty(f_hat, df(:,:,k))
            end do 
        else
            df = zero
        end if 
    
    end subroutine 

    subroutine dd3(this,f, df)
        class(ffts), intent(in) :: this
        real(rkind), dimension(this%n2,this%n3, this%n), intent(in) :: f 
        real(rkind), dimension(this%n2,this%n3, this%n), intent(out) :: df 
        complex(rkind), dimension(this%n2,this%n3,this%split) :: f_hat

        if ((this%dir .eq. "z") .and. (this%n > 1)) then
           call this%fftz(f,f_hat)

           f_hat = imi*this%k1d*f_hat
           f_hat(:,:,this%split) = zero

           call this%ifftz(f_hat,df)
        else 
            df = zero
        end if 
    end subroutine 

    subroutine fftx(this,in_arr, out_arr)
        class(ffts), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in) :: in_arr
        complex(rkind), dimension(:,:,:), intent(out) :: out_arr

        call dfftw_execute_dft_r2c(this%plan_fwd,in_arr,out_arr)

    end subroutine
   
    subroutine ffty(this,in_arr, out_arr)
        class(ffts), intent(in) :: this
        real(rkind), dimension(:,:), intent(in) :: in_arr
        complex(rkind), dimension(:,:), intent(out) :: out_arr

        call dfftw_execute_dft_r2c(this%plan_fwd,in_arr,out_arr)

    end subroutine
    
    subroutine fftz(this,in_arr, out_arr)
        class(ffts), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in) :: in_arr
        complex(rkind), dimension(:,:,:), intent(out) :: out_arr

        call dfftw_execute_dft_r2c(this%plan_fwd,in_arr,out_arr)

    end subroutine


    subroutine ifftx(this,in_arr, out_arr)
        class(ffts), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(out) :: out_arr
        complex(rkind), dimension(:,:,:), intent(in) :: in_arr

        call dfftw_execute_dft_c2r(this%plan_bwd,in_arr,out_arr)

        out_arr = out_arr*this%onebyn
    end subroutine

    subroutine iffty(this,in_arr, out_arr)
        class(ffts), intent(in) :: this
        real(rkind), dimension(:,:), intent(out) :: out_arr
        complex(rkind), dimension(:,:), intent(in) :: in_arr

        call dfftw_execute_dft_c2r(this%plan_bwd,in_arr,out_arr)

        out_arr = out_arr*this%onebyn
    
    end subroutine
    
    subroutine ifftz(this,in_arr, out_arr)
        class(ffts), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(out) :: out_arr
        complex(rkind), dimension(:,:,:), intent(in) :: in_arr

        call dfftw_execute_dft_c2r(this%plan_bwd,in_arr,out_arr)

        out_arr = out_arr*this%onebyn
    
    end subroutine

    pure function GetWaveNums(nx,dx) result(k_small)

        integer, intent(in) :: nx
        real(rkind), intent(in) :: dx
        real(rkind), dimension(nx) :: k
        real(rkind), dimension(nx) :: k_small

        integer :: i,dummy

        dummy = nx - MOD(nx,2)

        do i = 1,nx
            k(i) = ( -pi + (i-1)*two*pi/real(dummy,rkind) ) / dx
        end do

        k = ifftshift(k)

        k_small = k(1:(nx/2 + 1))
    end function

    pure function ifftshift(k) result(kshift)

        real(rkind), dimension(:), intent(in) :: k
        real(rkind), dimension(SIZE(k)) :: kshift
        integer :: n

        n = SIZE(k)

        select case ( MOD(n,2) )
        case (0)
            kshift(1:n/2) = k(n/2+1:n)
            kshift(n/2+1:n) = k(1:n/2)
        case (1)
            kshift(1:(n+1)/2) = k((n+1)/2:n)
            kshift((n+1)/2+1:n) = k(1:(n-1)/2)
        end select

    end function

    subroutine destroy(this) 
        class(ffts), intent(inout) :: this

        if (allocated(this%k1d)) deallocate (this%k1d)
        call dfftw_destroy_plan ( this%plan_fwd )
        call dfftw_destroy_plan ( this%plan_bwd )

        this%initialized = .false.
    end subroutine 

end module 
