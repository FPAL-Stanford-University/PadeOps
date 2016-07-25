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

        real(rkind) :: normfact, orthofactAll, orthofactFirst

        contains 
            procedure :: init
            procedure :: destroy
            procedure :: dctx
            procedure :: idctx
    end type

contains
    subroutine init(this, n, n1, n2)
        class(dct1d), intent(inout) :: this
        integer, intent(in) :: n, n1, n2
        real(rkind), dimension(:), allocatable :: arr_in, arr_out
         
        this%n = n; this%n2 = n2; this%n1 = n1
        this%normfact = one/(two*real(n,rkind))
       
        allocate(arr_in(n), arr_out(n))
        call dfftw_plan_r2r_1d( this%plan_fwdDCT, n, arr_in, arr_out, FFTW_REDFT10, FFTW_EXHAUSTIVE )   
        call dfftw_plan_r2r_1d( this%plan_bwdDCT, n, arr_in, arr_out, FFTW_REDFT01, FFTW_EXHAUSTIVE )   

        this%orthofactAll = sqrt(one/(two*real(n,rkind)))
        this%orthofactFirst = sqrt(one/two)
    end subroutine


    subroutine destroy(this)

        class( dct1d ), intent(inout) :: this

        call dfftw_destroy_plan ( this%plan_fwdDCT )
        call dfftw_destroy_plan ( this%plan_bwdDCT )


    end subroutine


    subroutine dctx(this, input, output)
        class( dct1d ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: input
        real(rkind), dimension(this%n), intent(out) :: output


        call dfftw_execute_r2r( this%plan_fwdDCT, input, output )
        output = this%orthoFactAll*output
        output(1) = this%orthoFactFirst*output(1)

    end subroutine 

    subroutine idctx(this, input, output)
        class( dct1d ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: input
        real(rkind), dimension(this%n), intent(out) :: output
        real(rkind), dimension(this%n) :: tmp

        tmp = input/this%orthofactAll
        tmp(1) = tmp(1)/this%orthoFactFirst
        
        call dfftw_execute_r2r( this%plan_bwdDCT, tmp, output )
        
        output = output*this%normfact

    end subroutine 
end module 
