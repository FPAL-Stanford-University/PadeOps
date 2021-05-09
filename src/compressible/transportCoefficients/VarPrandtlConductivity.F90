module VarPrandtlConductivityMod

    use kind_parameters,        only: rkind
    use ThermalConductivityMod, only: thermalConductivity

    implicit none

    type, extends(thermalConductivity) :: varPrandtlConductivity

        real(rkind) :: RePr = 0.75_rkind
        real(rkind) :: kappa_ref = 1.0_rkind
        real(rkind) :: T_ref = 1.0_rkind

    contains

        procedure :: get_kappa

    end type

    interface varPrandtlConductivity
        module procedure init
    end interface

contains

    function init(RePr, kappa_ref, T_ref) result(this)
        type(varPrandtlConductivity) :: this
        real(rkind),           intent(in) :: RePr, kappa_ref, T_ref

        this%RePr       = RePr
        this%kappa_ref  = kappa_ref
        this%T_ref      = T_ref

    end function

    pure subroutine get_kappa(this, T, Cp, mu, simtime, kappa)
        class(varPrandtlConductivity),           intent(in) :: this
        real(rkind), dimension(:,:,:),           intent(in)  :: T, Cp, mu
        real(rkind),                             intent(in)  :: simtime
        real(rkind), dimension(:,:,:),           intent(out) :: kappa

        kappa = this%kappa_ref*sqrt(T/this%T_ref) / this%RePr

    end subroutine

end module
