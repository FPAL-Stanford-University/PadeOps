module dct_1d_stuff

    use kind_parameters, only : rkind
    use constants, only : one, two 
    implicit none 
    private
    public :: dct1d

    include "fftw3.f"

    type dct1d
        private
        integer :: n, n1, n2
        integer(kind=8) :: plan_fwdDCT
        integer(kind=8) :: plan_bwdDCT
        integer(kind=8) :: plan_fwdDCT_ip 
        integer(kind=8) :: plan_bwdDCT_ip
        logical :: ortho = .false. 
        real(rkind) :: normfact, orthofactAll, orthofactFirst, orthofactAll_inv, orthofactFirst_inv

        contains 
            procedure :: init
            procedure :: destroy
            procedure, private :: dctx_oop
            procedure, private :: idctx_oop
            procedure, private :: dctx_ip
            procedure, private :: idctx_ip
            generic :: dctx => dctx_oop, dctx_ip
            generic :: idctx => idctx_oop, idctx_ip
    end type

contains
    subroutine init(this, n, n1, n2, ortho)
        class(dct1d), intent(inout) :: this
        integer, intent(in) :: n, n1, n2
        real(rkind), dimension(:), allocatable :: arr_in, arr_out
        logical :: ortho

        this%n = n; this%n2 = n2; this%n1 = n1
        this%normfact = one/(two*real(n,rkind))
        this%ortho = ortho

        allocate(arr_in(n), arr_out(n))
        call dfftw_plan_r2r_1d( this%plan_fwdDCT, n, arr_in, arr_out, FFTW_REDFT10, FFTW_EXHAUSTIVE )   
        call dfftw_plan_r2r_1d( this%plan_bwdDCT, n, arr_in, arr_out, FFTW_REDFT01, FFTW_EXHAUSTIVE )   

        call dfftw_plan_r2r_1d( this%plan_fwdDCT_ip, n, arr_in, arr_in, FFTW_REDFT10, FFTW_EXHAUSTIVE )   
        call dfftw_plan_r2r_1d( this%plan_bwdDCT_ip, n, arr_in, arr_in, FFTW_REDFT01, FFTW_EXHAUSTIVE )   
        
        this%orthofactAll = sqrt(one/(two*real(n,rkind)))
        this%orthofactFirst = sqrt(one/two)
        this%orthofactAll_inv = this%normfact*one/sqrt(one/(two*real(n,rkind)))
        this%orthofactFirst_inv = one/sqrt(one/two)
    end subroutine


    subroutine destroy(this)

        class( dct1d ), intent(inout) :: this

        call dfftw_destroy_plan ( this%plan_fwdDCT )
        call dfftw_destroy_plan ( this%plan_bwdDCT )


    end subroutine

    subroutine dctx_oop(this, input, output)
        class( dct1d ), intent(in) :: this
        real(rkind), dimension(this%n,this%n1,this%n2), intent(in) :: input
        real(rkind), dimension(this%n,this%n1,this%n2), intent(out) :: output
        integer :: j, k

        select case(this%ortho)
        case(.true.)
            do k = 1,this%n2
                do j = 1,this%n1
                    call dfftw_execute_r2r( this%plan_fwdDCT, input(:,j,k), output(:,j,k) )
                    output(:,j,k) = this%orthoFactAll*output(:,j,k)
                    output(1,j,k) = this%orthoFactFirst*output(1,j,k)
                end do 
            end do 
        case(.false.)
            do k = 1,this%n2
                do j = 1,this%n1
                    call dfftw_execute_r2r( this%plan_fwdDCT, input(:,j,k), output(:,j,k) )
                end do 
            end do 
        end select

    end subroutine 

    subroutine idctx_oop(this, input, output)
        class( dct1d ), intent(in) :: this
        real(rkind), dimension(this%n,this%n1,this%n2), intent(in) :: input
        real(rkind), dimension(this%n,this%n1,this%n2), intent(out) :: output
        real(rkind), dimension(this%n) :: tmp
        integer :: j, k


        select case(this%ortho)
        case(.true.)
        do k = 1,this%n2
            do j = 1,this%n1
                tmp = input(:,j,k)*this%orthofactAll_inv
                tmp(1) = tmp(1)*this%orthoFactFirst_inv
                call dfftw_execute_r2r( this%plan_bwdDCT, tmp, output(:,j,k) )
            end do 
        end do 
        case(.false.)
            do k = 1,this%n2
                do j = 1,this%n1
                    tmp = input(:,j,k)*this%normfact
                    call dfftw_execute_r2r( this%plan_bwdDCT, tmp, output(:,j,k) )
                end do 
            end do 
            !output = this%normfact*output
        end select

    end subroutine

    subroutine dctx_ip(this, input)
        class( dct1d ), intent(in) :: this
        real(rkind), dimension(this%n,this%n1,this%n2), intent(inout) :: input
        integer :: j, k

        select case(this%ortho)
        case(.true.)
            do k = 1,this%n2
                do j = 1,this%n1
                    call dfftw_execute_r2r( this%plan_fwdDCT, input(:,j,k), input(:,j,k) )
                    input(1,j,k) = this%orthoFactFirst*input(1,j,k)
                    input(:,j,k) = this%orthoFactAll*input(:,j,k)
                end do 
            end do 
        case(.false.)
            do k = 1,this%n2
                do j = 1,this%n1
                    call dfftw_execute_r2r( this%plan_fwdDCT, input(:,j,k), input(:,j,k) )
                end do 
            end do 
        end select

    end subroutine 

    subroutine idctx_ip(this, input)
        class( dct1d ), intent(in) :: this
        real(rkind), dimension(this%n,this%n1,this%n2), intent(inout) :: input
        real(rkind), dimension(this%n) :: tmp
        integer :: j, k


        select case(this%ortho)
        case(.true.)
        do k = 1,this%n2
            do j = 1,this%n1
                input(:,j,k) = input(:,j,k)*this%orthofactAll_inv
                input(1,j,k) = input(1,j,k)*this%orthoFactFirst_inv
                call dfftw_execute_r2r( this%plan_bwdDCT, input(:,j,k), input(:,j,k) )
            end do 
        end do 
        case(.false.)
            do k = 1,this%n2
                do j = 1,this%n1
                    input(:,j,k) = input(:,j,k)*this%normfact
                    call dfftw_execute_r2r( this%plan_bwdDCT, input(:,j,k), input(:,j,k) )
                end do 
            end do 
            !input = this%normfact*input
        end select

    end subroutine

end module 
