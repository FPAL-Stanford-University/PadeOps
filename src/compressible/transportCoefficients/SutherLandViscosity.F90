module SutherLandViscosityMod

    use kind_parameters,   only: rkind
    use constants,         only: one
    use ShearViscosityMod, only: shearViscosity

    implicit none

    type, extends(shearViscosity) :: sutherlandViscosity

        real(rkind) :: mu_ref = one
        real(rkind) :: T_ref  = one
        real(rkind) :: n = 1.5_rkind
        real(rkind) :: S = 110.4_rkind 

    contains

        procedure :: get_mu

    end type

    interface sutherlandViscosity
        module procedure init
    end interface

contains

    function init(mu_ref, T_ref, n, S) result(this)
        type(sutherlandViscosity) :: this
        real(rkind),           intent(in) :: mu_ref, T_ref, n, S

        this%mu_ref = mu_ref
        this%T_ref  = T_ref
        this%n      = n
        this%S      = S
    end function

    pure subroutine get_mu(this, T, mu)
        class(sutherlandViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T
        real(rkind), dimension(:,:,:), intent(out) :: mu

        mu = this%mu_ref * ( T/this%T_ref )**this%n * ((this%T_ref + this%S) / (T + this%S))
    
    end subroutine

end module
