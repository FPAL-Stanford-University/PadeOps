module SolidMixtureMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,one
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use EOSMod,          only: eos
    use StiffGasEOS,     only: stiffgas
    use Sep1SolidEOS,    only: sep1solid
    use SolidMod,        only: solid

    implicit none

    type :: solid_mixture

        integer :: ns
        integer :: nxp, nyp, nzp
        type(solid), dimension(:), allocatable :: material

    contains

        procedure :: set_material
        procedure :: relaxPressure
        final     :: destroy

    end type

    interface solid_mixture
        module procedure init
    end interface

contains

    function init(decomp,der,ns) result(this)
        type(solid_mixture)               :: this
        type(decomp_info),  intent(in)    :: decomp
        type(derivatives),  intent(in)    :: der
        integer,            intent(in)    :: ns

        type(solid), allocatable :: dummy
        integer :: i

        this%ns = ns
        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        ! Allocate array of solid objects (Use a dummy to avoid memory leaks)
        allocate(dummy, source=solid(decomp,der))
        if (allocated(this%material)) deallocate(this%material)
        allocate(this%material(this%ns), source=dummy)
        deallocate(dummy)

    end function

    pure elemental subroutine destroy(this)
        type(solid_mixture), intent(inout)  :: this
        integer :: i

        ! Deallocate array of solids (Destructor of solid should take care of everything else)
        if (allocated(this%material)) deallocate(this%material)
    end subroutine

    subroutine set_material(this, imat, hydro, elastic)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: imat
        class(stiffgas ),     intent(in)    :: hydro
        class(sep1solid),     intent(in)    :: elastic

        if ((imat .GT. this%ns) .OR. (imat .LE. 0)) call GracefulExit("Cannot set material with index greater than the number of species.",4534)

        if (allocated(this%material(imat)%hydro)) deallocate(this%material(imat)%hydro)
        allocate( this%material(imat)%hydro, source=hydro )
        
        if (allocated(this%material(imat)%elastic)) deallocate(this%material(imat)%elastic)
        allocate( this%material(imat)%elastic, source=elastic )
    end subroutine

end module
