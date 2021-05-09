module VariableConductivityMod

    use kind_parameters,        only: rkind
    use constants,              only: zero, one
    use ThermalConductivityMod, only: thermalConductivity

    implicit none

    type, extends(thermalConductivity) :: variableConductivity

        real(rkind) :: kap_ref = one
        real(rkind) :: T_ref  = one
        real(rkind) :: n = 0.75_rkind
        real(rkind) :: t_start = 100_rkind
        real(rkind) :: t_width = one
        integer :: cond_type = 0

    contains

        procedure :: get_kappa

    end type

    interface variableConductivity
        module procedure init
    end interface

contains

    function init(kap_ref, T_ref, n, t_start, t_width, cond_type) result(this)
        type(variableConductivity) :: this
        real(rkind),           intent(in) :: kap_ref, T_ref, n, t_start, t_width
        integer,               intent(in) :: cond_type

        this%kap_ref = kap_ref
        this%T_ref  = T_ref
        this%n      = n
        this%t_start = t_start
        this%t_width = t_width
        this%cond_type = cond_type
    end function

    pure subroutine get_kappa(this, T, Cp, mu, simtime, kappa)
        class(variableConductivity), intent(in) :: this
        real(rkind), dimension(:,:,:),           intent(in)  :: T, Cp, mu
        real(rkind),                             intent(in)  :: simtime
        real(rkind), dimension(:,:,:),           intent(out) :: kappa

        real(rkind) :: tmp

        select case(this%cond_type)
        case(1)
            kappa = this%kap_ref * ( T/this%T_ref )**this%n
        case(2)
            tmp = 0.5_rkind*(1.0_rkind + erf((simtime-this%t_start)/this%t_width))
            kappa = this%kap_ref*(1.0_rkind-tmp) + tmp * this%kap_ref * ( T/this%T_ref )**this%n
        case(3)
            if (simtime < this%t_start-this%t_width/2.0_rkind) then 
                tmp = zero
            else if (simtime > this%t_start+this%t_width/2.0_rkind) then
                tmp = one
            else
                tmp = one/this%t_width*(simtime-this%t_start+this%t_width/2.0_rkind)
            end if
            kappa = this%kap_ref*(1.0_rkind-tmp) + tmp * this%kap_ref * ( T/this%T_ref )**this%n
        case default
            ! or case 0
            kappa = this%kap_ref
        end select

    end subroutine

end module
