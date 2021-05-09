module VariableViscosityMod

    use kind_parameters,   only: rkind
    use constants,         only: zero, one
    use ShearViscosityMod, only: shearViscosity

    implicit none

    type, extends(shearViscosity) :: variableViscosity

        real(rkind) :: mu_ref = one
        real(rkind) :: T_ref  = one
        real(rkind) :: n = 0.75_rkind
        real(rkind) :: t_start = 100_rkind
        real(rkind) :: t_width = one
        integer :: visc_type = 0

    contains

        procedure :: get_mu

    end type

    interface variableViscosity
        module procedure init
    end interface

contains

    function init(mu_ref, T_ref, n, t_start, t_width, visc_type) result(this)
        type(variableViscosity) :: this
        real(rkind),           intent(in) :: mu_ref, T_ref, n, t_start, t_width
        integer,               intent(in) :: visc_type

        this%mu_ref = mu_ref
        this%T_ref  = T_ref
        this%n      = n
        this%t_start = t_start
        this%t_width = t_width
        this%visc_type = visc_type
    end function

    pure subroutine get_mu(this, T, simtime, mu)
        class(variableViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T
        real(rkind),                   intent(in)  :: simtime
        real(rkind), dimension(:,:,:), intent(out) :: mu

        real(rkind) :: tmp

        select case(this%visc_type)
        case(1)
            mu = this%mu_ref * ( T/this%T_ref )**this%n
        case(2)
            tmp = 0.5_rkind*(1.0_rkind + erf((simtime-this%t_start)/this%t_width))
            mu = this%mu_ref*(1.0_rkind-tmp) + tmp * this%mu_ref * ( T/this%T_ref )**this%n
        case(3)
            if (simtime < this%t_start-this%t_width/2.0_rkind) then
                tmp = zero
            else if (simtime > this%t_start+this%t_width/2.0_rkind) then
                tmp = one
            else
                tmp = one/this%t_width*(simtime-this%t_start+this%t_width/2.0_rkind)
            end if
            mu = this%mu_ref*(1.0_rkind-tmp) + tmp * this%mu_ref * ( T/this%T_ref )**this%n
        case default
            ! or case 0
            mu = this%mu_ref
        end select
    end subroutine

end module
