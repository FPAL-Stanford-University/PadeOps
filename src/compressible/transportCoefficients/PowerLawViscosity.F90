module PowerLawViscosityMod

    use kind_parameters,   only: rkind
    use constants,         only: one
    use ShearViscosityMod, only: shearViscosity

    implicit none

    type, extends(shearViscosity) :: powerLawViscosity

        real(rkind) :: mu_ref = one
        real(rkind) :: T_ref  = one
        real(rkind) :: n = 0.75_rkind

    contains

        procedure :: get_mu

    end type

    interface powerLawViscosity
        module procedure init
    end interface

contains

    function init(mu_ref, T_ref, n) result(this)
        type(powerLawViscosity) :: this
        real(rkind),           intent(in) :: mu_ref, T_ref, n

        this%mu_ref = mu_ref
        this%T_ref  = T_ref
        this%n      = n
    end function

    pure subroutine get_mu(this, T, simtime, mu)
        class(powerLawViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T
        real(rkind),                   intent(in)  :: simtime
        real(rkind), dimension(:,:,:), intent(out) :: mu

        mu = this%mu_ref * ( T/this%T_ref )**this%n

    end subroutine

end module
