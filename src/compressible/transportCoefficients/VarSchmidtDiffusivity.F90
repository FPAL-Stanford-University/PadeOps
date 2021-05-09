module VarSchmidtDiffusivityMod

    use kind_parameters,    only: rkind
    use constants,          only: zero,one
    use MassDiffusivityMod, only: massDiffusivity

    implicit none

    type, extends(massDiffusivity) :: varSchmidtDiffusivity

        ! Compute constant diffusivity for all species
        ! based on a fixed Schmidt number.
        ! Intended mainly for use with 2 species

        real(rkind) :: mu_ref  = one
        real(rkind) :: rho_ref = one
        real(rkind) :: Sc      = one
        real(rkind) :: T_ref = 1.0_rkind

    contains

        procedure :: get_diff

    end type

    interface varSchmidtDiffusivity
        module procedure init
    end interface

contains

    function init(mu_ref, rho_ref, Sc, T_ref) result(this)
        type(varSchmidtDiffusivity) :: this
        real(rkind),       intent(in) :: mu_ref, rho_ref, Sc, T_ref

        this%mu_ref  = mu_ref
        this%rho_ref = rho_ref
        this%Sc      = Sc
        this%T_ref   = T_ref

    end function

    pure subroutine get_diff(this, p, T, Xs, simtime, diff)
        class(varSchmidtDiffusivity),    intent(in)  :: this
        real(rkind), dimension(:,:,:),   intent(in)  :: p, T
        real(rkind), dimension(:,:,:,:), intent(in)  :: Xs
        real(rkind),                     intent(in)  :: simtime
        real(rkind), dimension(:,:,:,:), intent(out) :: diff

        ! assuming 2 species for now
        diff(:,:,:,1) = (T/this%T_ref)**1.50_rkind * this%mu_ref / (this%rho_ref * this%Sc)
        diff(:,:,:,2) = (T/this%T_ref)**1.50_rkind * this%mu_ref / (this%rho_ref * this%Sc)

    end subroutine

end module
