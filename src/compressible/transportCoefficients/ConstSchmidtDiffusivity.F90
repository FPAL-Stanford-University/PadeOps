module ConstSchmidtDiffusivityMod

    use kind_parameters,    only: rkind
    use constants,          only: zero,one
    use MassDiffusivityMod, only: massDiffusivity

    implicit none

    type, extends(massDiffusivity) :: constSchmidtDiffusivity

        ! Compute constant diffusivity for all species
        ! based on a fixed Schmidt number.
        ! Intended mainly for use with 2 species

        real(rkind) :: mu_ref  = one
        real(rkind) :: rho_ref = one
        real(rkind) :: Sc      = one

    contains

        procedure :: get_diff

    end type

    interface constSchmidtDiffusivity
        module procedure init
    end interface

contains

    function init(mu_ref, rho_ref, Sc) result(this)
        type(constSchmidtDiffusivity) :: this
        real(rkind),       intent(in) :: mu_ref, rho_ref, Sc

        this%mu_ref  = mu_ref
        this%rho_ref = rho_ref
        this%Sc      = Sc

    end function

    pure subroutine get_diff(this, p, T, Xs, diff)
        class(constSchmidtDiffusivity),  intent(in)  :: this
        real(rkind), dimension(:,:,:),   intent(in)  :: p, T
        real(rkind), dimension(:,:,:,:), intent(in)  :: Xs
        real(rkind), dimension(:,:,:,:), intent(out) :: diff

        diff = this%mu_ref / (this%rho_ref * this%Sc)

    end subroutine

end module
