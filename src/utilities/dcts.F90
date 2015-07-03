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
       
        procedure            :: dct                                     ! Compute Forward Transform 
        
        procedure            :: idct                                      ! Compute Backward Transform 
    
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
   
    function dct(this, f) result(fhat)
   
        class( dcts ), intent(in)                       :: this
        real(rkind), dimension(this%n), intent(in)   :: f
        real(rkind), dimension(this%n)                  :: fhat
        real(rkind), parameter :: minusOne = -one

        integer                                     :: i
        real(rkind)                                 :: dummy


        call dfftw_execute_r2r( this%plan_fwd, f, fhat )
        
        fhat(this%n) = half*fhat(this%n)
        
        do i = this%n,2,-2
            fhat(i) = fhat(i)*this%onebynxm1
            fhat(i-1) = minusOne*fhat(i-1)*this%onebynxm1
        end do 
        
        if (i == 1) then
            fhat(1) = fhat(1)*this%onebynxm1
        end if

        fhat(1) = half*fhat(1)
  
    end function 

    function idct(this, fhat) result(f)
    
        class( dcts ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: fhat
        real(rkind), dimension(this%n) :: f
        real(rkind), parameter :: minusOne = -one
        integer :: i
        real(rkind) :: dummy

        fhat(2:this%n - 1) = half*fhat(2:this%n - 1)

        do i = this%n,2,-2
            fhat(i-1) = minusOne*fhat(i-1)
        end do 

        call dfftw_execute_r2r( this%plan_bwd, fhat, f )
    
    end function 

end module 
