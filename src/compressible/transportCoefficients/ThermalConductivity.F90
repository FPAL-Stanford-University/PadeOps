module ThermalConductivityMod

    use kind_parameters, only: rkind

    implicit none

    type, abstract :: thermalConductivity

    contains

        procedure(get_kappa_interface), deferred :: get_kappa

    end type

    abstract interface

        pure subroutine get_kappa_interface(this, T, Cp, mu, simtime, kappa)
            import :: thermalConductivity
            import :: rkind
            class(thermalConductivity),    intent(in)  :: this
            real(rkind), dimension(:,:,:), intent(in)  :: T, Cp, mu
            real(rkind),                   intent(in)  :: simtime
            real(rkind), dimension(:,:,:), intent(out) :: kappa
        end subroutine

    end interface

contains

end module
