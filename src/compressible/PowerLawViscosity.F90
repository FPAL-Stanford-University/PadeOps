module PowerLawViscosityMod

    use kind_parameters, only: rkind
    use constants,       only: zero,one

    implicit none

    type :: powerLawViscosity

        real(rkind) :: mu_ref = one
        real(rkind) :: T_ref  = one
        real(rkind) :: n = 0.75_rkind
        real(rkind) :: Pr = 0.75_rkind
        real(rkind) :: Sc = one

    contains

        procedure :: get_mu
        procedure :: get_beta
        procedure :: get_kappa
        procedure :: get_diff

    end type

    interface powerLawViscosity
        module procedure init
    end interface

contains

    function init(mu_ref, T_ref, n, Pr, Sc) result(this)
        type(powerLawViscosity) :: this
        real(rkind),           intent(in) :: mu_ref, T_ref, n, Pr
        real(rkind), optional, intent(in) :: Sc

        this%mu_ref = mu_ref
        this%T_ref  = T_ref
        this%n      = n
        this%Pr     = Pr

        if (present(Sc)) this%Sc = Sc

    end function

    pure subroutine get_mu(this, T, mu)
        class(powerLawViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T
        real(rkind), dimension(:,:,:), intent(out) :: mu

        mu = this%mu_ref * ( T/this%T_ref )**this%n

    end subroutine

    pure subroutine get_beta(this, T, beta)
        class(powerLawViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T
        real(rkind), dimension(:,:,:), intent(out) :: beta

        beta = zero

    end subroutine

    pure subroutine get_kappa(this, Cp, mu, kappa)
        class(powerLawViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:),           intent(in)  :: Cp, mu
        real(rkind), dimension(:,:,:),           intent(out) :: kappa

        kappa = Cp * mu / this%Pr

    end subroutine

    pure subroutine get_diff(this, rho, mu, diff)
        class(powerLawViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:),           intent(in)  :: rho, mu
        real(rkind), dimension(:,:,:),           intent(out) :: diff

        diff = mu / (rho * this%Sc)

    end subroutine

end module
