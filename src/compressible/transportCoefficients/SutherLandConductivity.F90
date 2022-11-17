module SutherLandConductivityMod

    use kind_parameters,        only: rkind
    use constants,              only: one
    use ThermalConductivityMod, only: thermalConductivity

    implicit none

    type, extends(thermalConductivity) :: sutherlandConductivity

        real(rkind) :: T_ref = one
        real(rkind) :: Pr    = 0.7_rkind
        real(rkind) :: n     = 1.5_rkind
        real(rkind) :: Sk    = 194.0_rkind

    contains

        procedure :: get_kappa

    end type

    interface sutherlandConductivity
        module procedure init
    end interface

contains

    function init(Pr,T_ref,n,Sk) result(this)
        type(sutherlandConductivity)      :: this
        real(rkind),           intent(in) :: Pr, T_ref, n, Sk

        this%Pr     = Pr
        this%T_ref  = T_ref
        this%n      = n
        this%Sk     = Sk

    end function

    pure subroutine get_kappa(this, T, Cp, mu, kappa)
        class(sutherlandConductivity), intent(in) :: this
        real(rkind), dimension(:,:,:),           intent(in)  :: T, Cp, mu
        real(rkind), dimension(:,:,:),           intent(out) :: kappa

        kappa = (Cp * mu / this%Pr) * ( T/this%T_ref )**this%n * ((this%T_ref +this%Sk) / (T + this%Sk)) 

    end subroutine

end module
