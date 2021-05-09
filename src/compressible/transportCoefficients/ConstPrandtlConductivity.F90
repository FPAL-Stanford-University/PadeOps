module ConstPrandtlConductivityMod

    use kind_parameters,        only: rkind
    use ThermalConductivityMod, only: thermalConductivity

    implicit none

    type, extends(thermalConductivity) :: constPrandtlConductivity

        real(rkind) :: Pr = 0.75_rkind

    contains

        procedure :: get_kappa

    end type

    interface constPrandtlConductivity
        module procedure init
    end interface

contains

    function init(Pr) result(this)
        type(constPrandtlConductivity) :: this
        real(rkind),           intent(in) :: Pr

        this%Pr     = Pr

    end function

    pure subroutine get_kappa(this, T, Cp, mu, simtime, kappa)
        class(constPrandtlConductivity), intent(in) :: this
        real(rkind), dimension(:,:,:),           intent(in)  :: T, Cp, mu
        real(rkind),                             intent(in)  :: simtime
        real(rkind), dimension(:,:,:),           intent(out) :: kappa

        kappa = Cp * mu / this%Pr

    end subroutine

end module
