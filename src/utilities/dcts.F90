module dctstuff

    use kind_parameters, only: rkind
    use constants, only: zero,one,half,two,pi,imi
    implicit none

    private
    public :: dcts

    include "fftw3.F90"

    type dcts
   
        private

        integer             :: n
        integer(kind=8)     :: plan_fwd, plan_bwd
        
        real(rkind)         :: onebynxm1
        
        logical             :: initialized = .false.

        contains

        procedure           :: init                                     ! Initialize type
       
        procedure, private :: dct_flipped
        procedure, private :: dct_unflipped
        generic            :: dct => dct_flipped, dct_unflipped         ! Compute Forward Transform 
        
        procedure, private :: idct_flipped
        procedure, private :: idct_unflipped
        generic            :: idct => idct_flipped, idct_unflipped      ! Compute Backward Transform 
    
    end type
        
contains
    
    function init(this, n_) result(ierr)
   
        class( dcts ), intent(inout)    :: this
        integer, intent(in)             :: n_
        integer                         :: ierr

        real(rkind), dimension(n_) :: arr_in, arr_out

        this%n = n_

        this%onebynxm1 = one / ( n_ - one )

        call dfftw_plan_r2r_1d( this%plan_fwd, n_, arr_in, arr_out, FFTW_REDFT00, FFTW_EXHAUSTIVE )   
        call dfftw_plan_r2r_1d( this%plan_bwd, n_, arr_in, arr_out, FFTW_REDFT00, FFTW_EXHAUSTIVE )   

        ierr = 0

    end function     
   
    function dct_unflipped(this, f, flipflag) result(fhat)
    
        class( dcts ), intent(in)                   :: this
        real(rkind), dimension(this%n), intent(in)  :: f
        logical, intent(in)                         :: flipflag
        real(rkind), dimension(this%n)              :: fhat

        call dfftw_execute_r2r( this%plan_fwd, f, fhat )
        fhat(1) = half*fhat(1)
        fhat = fhat*this%onebynxm1
        fhat(this%n) = half*fhat(this%n)

    end function 

    function dct_flipped(this, f) result(fhat)
   
        class( dcts ), intent(in)                   :: this
        real(rkind), dimension(this%n), intent(in)  :: f
        real(rkind), dimension(this%n)              :: fhat

        integer                                     :: i
        real(rkind)                                 :: dummy

        call dfftw_execute_r2r( this%plan_fwd, f, fhat )
        fhat(1) = half*fhat(1)
        fhat = fhat*this%onebynxm1
        fhat(this%n) = half*fhat(this%n)

        do i = 1,this%n/2
            dummy = fhat(i)
            fhat(i) = fhat(this%n + 1 - i)
            fhat(this%n + 1 - i) = dummy
        end do

    end function 

    function idct_unflipped(this, fhat, flipflag) result(f)

        class( dcts ), intent(in)                     :: this
        real(rkind), dimension(this%n), intent(inout) :: fhat
        real(rkind), dimension(this%n)                :: f
        logical, intent(in)                           :: flipflag

        fhat(2:this%n - 1) = half*fhat(2:this%n - 1)

        call dfftw_execute_r2r( this%plan_bwd, fhat, f )
    
    end function 

    function idct_flipped(this, fhat) result(f)
    
        class( dcts ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: fhat
        real(rkind), dimension(this%n) :: f

        integer :: i
        real(rkind) :: dummy

        fhat(2:this%n - 1) = half*fhat(2:this%n - 1)

        do i = 1,this%n/2
            dummy = fhat(i)
            fhat(i) = fhat(this%n + 1 - i)
            fhat(this%n + 1 - i) = dummy
        end do

        call dfftw_execute_r2r( this%plan_bwd, fhat, f )
    
    end function 

end module 
