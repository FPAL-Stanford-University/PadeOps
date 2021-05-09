module VariableDiffusivityMod

    use kind_parameters,    only: rkind
    use constants,          only: zero,one
    use MassDiffusivityMod, only: massDiffusivity

    implicit none

    type, extends(massDiffusivity) :: variableDiffusivity

        ! Compute constant diffusivity for all species
        ! based on a fixed Schmidt number.
        ! Intended mainly for use with 2 species

        real(rkind) :: diff_ref = one
        real(rkind) :: T_ref  = one
        real(rkind) :: n = 0.75_rkind
        real(rkind) :: t_start = 100_rkind
        real(rkind) :: t_width = one
        integer :: diff_type = 0

    contains

        procedure :: get_diff

    end type

    interface variableDiffusivity
        module procedure init
    end interface

contains

    function init(diff_ref, T_ref, n, t_start, t_width, diff_type) result(this)
        type(variableDiffusivity) :: this
        real(rkind),       intent(in) :: diff_ref, T_ref, n, t_start, t_width
        integer,           intent(in) :: diff_type

        this%diff_ref = diff_ref
        this%T_ref  = T_ref
        this%n      = n
        this%t_start = t_start
        this%t_width = t_width
        this%diff_type = diff_type
    end function

    pure subroutine get_diff(this, p, T, Xs, simtime, diff)
        class(variableDiffusivity),  intent(in)  :: this
        real(rkind), dimension(:,:,:),   intent(in)  :: p, T
        real(rkind), dimension(:,:,:,:), intent(in)  :: Xs
        real(rkind),                     intent(in)  :: simtime
        real(rkind), dimension(:,:,:,:), intent(out) :: diff

        real(rkind) :: tmp

        select case(this%diff_type)
        case(1)
            diff(:,:,:,1) = this%diff_ref * ( T/this%T_ref )**this%n
            diff(:,:,:,2) = this%diff_ref * ( T/this%T_ref )**this%n
        case(2)
            tmp = 0.5_rkind*(1.0_rkind + erf((simtime-this%t_start)/this%t_width))
            diff(:,:,:,1) = this%diff_ref*(1.0_rkind-tmp) + tmp * this%diff_ref * ( T/this%T_ref )**this%n
            diff(:,:,:,2) = this%diff_ref*(1.0_rkind-tmp) + tmp * this%diff_ref * ( T/this%T_ref )**this%n
        case(3)
            if (simtime < this%t_start-this%t_width/2.0_rkind) then
                tmp = zero
            else if (simtime > this%t_start+this%t_width/2.0_rkind) then
                tmp = one
            else
                tmp = one/this%t_width*(simtime-this%t_start+this%t_width/2.0_rkind)
            end if
            diff(:,:,:,1) = this%diff_ref*(1.0_rkind-tmp) + tmp * this%diff_ref * ( T/this%T_ref )**this%n
            diff(:,:,:,2) = this%diff_ref*(1.0_rkind-tmp) + tmp * this%diff_ref * ( T/this%T_ref )**this%n
        case default
            ! or case 0
            diff = this%diff_ref
        end select
    end subroutine

end module
