module MixtureEOSMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,one
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use EOSMod,          only: eos
    use IdealGasEOS,     only: idealgas

    type :: material_eos
        class(idealgas), allocatable :: mat
    end type

    type :: mixture

        integer :: ns
        integer :: nxp, nyp, nzp
        type(material_eos), dimension(:), allocatable :: material
        real(rkind), dimension(:,:,:), allocatable :: gam, Rgas, Cp, Cv

    contains

        procedure :: init
        procedure :: set_material
        procedure :: update
        procedure :: get_p
        procedure :: get_T
        procedure :: get_e_from_p
        procedure :: get_sos
        final     :: destroy

    end type

    ! interface mixture
    !     module procedure init
    ! end interface

contains

    ! function init(decomp,ns) result(this)
    subroutine init(this,decomp,ns)
        class(mixture),                                intent(inout) :: this
        type(decomp_info),                             intent(in)    :: decomp
        integer,                                       intent(in)    :: ns

        this%ns  = ns
        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        if (allocated(this%material)) deallocate(this%material); allocate(this%material(this%ns))

        if (this%ns .GT. 1) then
            if (allocated(this%gam )) deallocate(this%gam ); allocate(this%gam (this%nxp,this%nyp,this%nzp))
            if (allocated(this%Rgas)) deallocate(this%Rgas); allocate(this%Rgas(this%nxp,this%nyp,this%nzp))
            if (allocated(this%Cv  )) deallocate(this%Cv  ); allocate(this%Cv  (this%nxp,this%nyp,this%nzp))
            if (allocated(this%Cp  )) deallocate(this%Cp  ); allocate(this%Cp  (this%nxp,this%nyp,this%nzp))
        end if
    end subroutine
    ! end function

    subroutine set_material(this, imat, mat)
        class(mixture),  intent(inout) :: this
        integer,         intent(in)    :: imat
        class(idealgas), intent(in)    :: mat

        if ((imat .GT. this%ns) .OR. (imat .LE. 0)) call GracefulExit("Cannot set material with index greater than the number of species.",4534)

        if (allocated(this%material(imat)%mat)) deallocate(this%material(imat)%mat)
        allocate( this%material(imat)%mat, source=mat )
    end subroutine

    pure subroutine update(this,Ys)
        class(mixture), target,                                     intent(inout) :: this 
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(in)    :: Ys
        
        class(idealgas), pointer :: mat
        integer :: i
        
        if (this%ns .GT. 1) then
            this%Rgas = zero
            this%Cp = zero
            do i = 1,this%ns
                mat => this%material(i)%mat
                this%Rgas = this%Rgas + (mat%Rgas) * Ys(:,:,:,i)
                this%Cp = this%Cp + (mat%gam * mat%Rgas * mat%onebygam_m1 ) * Ys(:,:,:,i) ! Cp_i = gam/(gam-1) * Rgas
            end do

            this%gam = this%Cp / (this%Cp - this%Rgas)
            this%Cv = this%Cp / this%gam
        end if

    end subroutine

    pure subroutine get_p(this,rho,e,p)
        class(mixture),                                             intent(in)  :: this 
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: p

        select case(this%ns)
        case(1)
            call this%material(1)%mat%get_p(rho,e,p)
        case default
            p = (this%gam-one)*rho*e
        end select

    end subroutine

    pure subroutine get_T(this,e,T)
        class(mixture),                                             intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: T

        select case(this%ns)
        case(1)
            call this%material(1)%mat%get_T(e,T)
        case default
            T = e / this%Cv
        end select

    end subroutine

    pure subroutine get_e_from_p(this,rho,p,e)
        class(mixture),                                             intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,p
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: e

        select case(this%ns)
        case(1)
            call this%material(1)%mat%get_e_from_p(rho,p,e)
        case default
            e = p / ((this%gam - one)*rho)
        end select

    end subroutine

    pure subroutine get_sos(this,rho,p,sos)
        class(mixture),                                             intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,p
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: sos

        select case(this%ns)
        case(1)
            call this%material(1)%mat%get_sos(rho,p,sos)
        case default
            sos = sqrt(this%gam*p/rho)
        end select

    end subroutine

    subroutine destroy(this)
        type(mixture), intent(inout)  :: this
        integer :: i

        do i = 1,this%ns
            if (allocated(this%material(i)%mat)) deallocate(this%material(i)%mat)
        end do

        if (allocated(this%gam )) deallocate(this%gam )
        if (allocated(this%Rgas)) deallocate(this%Rgas)
        if (allocated(this%Cv  )) deallocate(this%Cv  )
        if (allocated(this%Cp  )) deallocate(this%Cp  )

        deallocate(this%material)
    end subroutine


end module
