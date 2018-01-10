module MixtureEOSMod

    use kind_parameters,        only: rkind,clen
    use constants,              only: zero,one
    use decomp_2d,              only: decomp_info
    use DerivativesMod,         only: derivatives
    use FiltersMod,             only: filters
    use exits,                  only: GracefulExit
    use EOSMod,                 only: eos
    use IdealGasEOS,            only: idealgas
    use ShearViscosityMod,      only: shearViscosity
    use BulkViscosityMod,       only: bulkViscosity
    use ThermalConductivityMod, only: thermalConductivity
    use MassDiffusivityMod,     only: massDiffusivity

    implicit none

    type :: material_eos
        class(idealgas), allocatable :: mat
        class(shearViscosity),      allocatable :: shearvisc
        class(bulkViscosity),       allocatable :: bulkvisc
        class(thermalConductivity), allocatable :: thermcond
    end type

    type :: mixture

        logical :: inviscid = .true.

        integer :: ns
        integer :: nxp, nyp, nzp
        type(material_eos), dimension(:), allocatable :: material
        class(massDiffusivity),           allocatable :: massdiff
        real(rkind), dimension(:,:,:), allocatable :: gam, Rgas, Cp, Cv

    contains

        procedure :: init
        procedure :: set_material
        procedure :: set_massdiffusivity
        procedure :: check_initialization
        procedure :: update
        procedure :: get_p
        procedure :: get_T
        procedure :: get_e_from_p
        procedure :: get_sos
        procedure :: get_transport_properties
        final     :: destroy

    end type

    ! interface mixture
    !     module procedure init
    ! end interface

contains

    ! function init(decomp,ns) result(this)
    subroutine init(this,decomp,ns,inviscid)
        class(mixture),                                intent(inout) :: this
        type(decomp_info),                             intent(in)    :: decomp
        integer,                                       intent(in)    :: ns
        logical, optional,                             intent(in)    :: inviscid

        this%ns  = ns
        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        if (present(inviscid)) then
            this%inviscid = inviscid
        else
            this%inviscid = .true.
        end if

        if (allocated(this%material)) deallocate(this%material); allocate(this%material(this%ns))

        if (this%ns .GT. 1) then
            if (allocated(this%gam )) deallocate(this%gam ); allocate(this%gam (this%nxp,this%nyp,this%nzp))
            if (allocated(this%Rgas)) deallocate(this%Rgas); allocate(this%Rgas(this%nxp,this%nyp,this%nzp))
            if (allocated(this%Cv  )) deallocate(this%Cv  ); allocate(this%Cv  (this%nxp,this%nyp,this%nzp))
            if (allocated(this%Cp  )) deallocate(this%Cp  ); allocate(this%Cp  (this%nxp,this%nyp,this%nzp))
        end if
    end subroutine
    ! end function

    subroutine set_material(this, imat, mat, shearvisc, bulkvisc, thermcond)
        class(mixture),                       intent(inout) :: this
        integer,                              intent(in)    :: imat
        class(idealgas),                      intent(in)    :: mat
        class(shearViscosity),      optional, intent(in)    :: shearvisc
        class(bulkViscosity),       optional, intent(in)    :: bulkvisc
        class(thermalConductivity), optional, intent(in)    :: thermcond
        ! class(shearViscosity),                intent(in)    :: shearvisc
        ! class(bulkViscosity),                 intent(in)    :: bulkvisc
        ! class(thermalConductivity),           intent(in)    :: thermcond

        if ((imat .GT. this%ns) .OR. (imat .LE. 0)) call GracefulExit("Cannot set material with index greater than the number of species.",4534)

        ! Allocate and set the material eos object
        if (allocated(this%material(imat)%mat)) deallocate(this%material(imat)%mat)
        allocate( this%material(imat)%mat, source=mat )

        if (.not. this%inviscid) then
            if (.not. present(shearvisc)) call GracefulExit("Cannot run viscous simulation without &
                                          &setting shear viscosity transport property object", 4534)

            ! Allocate and set the material transport property object
            if (allocated(this%material(imat)%shearvisc)) deallocate(this%material(imat)%shearvisc)
            allocate( this%material(imat)%shearvisc, source=shearvisc )

            if (.not. present(bulkvisc)) call GracefulExit("Cannot run viscous simulation without &
                                         &setting bulk viscosity transport property object", 4534)

            ! Allocate and set the material transport property object
            if (allocated(this%material(imat)%bulkvisc)) deallocate(this%material(imat)%bulkvisc)
            allocate( this%material(imat)%bulkvisc, source=bulkvisc )

            if (.not. present(thermcond)) call GracefulExit("Cannot run viscous simulation without &
                                          &setting thermal conductivity transport property object", 4534)

            ! Allocate and set the material transport property object
            if (allocated(this%material(imat)%thermcond)) deallocate(this%material(imat)%thermcond)
            allocate( this%material(imat)%thermcond, source=thermcond )

        end if
    end subroutine

    subroutine set_massdiffusivity(this, massdiff)
        class(mixture),         intent(inout) :: this
        class(massDiffusivity), intent(in)    :: massdiff

        if (.not. this%inviscid) then
            if (this%ns > 1) then
                ! Allocate and set the material transport property object
                if (allocated(this%massdiff)) deallocate(this%massdiff)
                allocate( this%massdiff, source=massdiff )
            end if
        end if
    end subroutine

    subroutine check_initialization(this)
        class(mixture), intent(in) :: this

        integer :: imat

        if (.not. this%inviscid) then
            do imat = 1,this%ns
                if (.not. allocated(this%material(imat)%shearvisc)) call GracefulExit("Cannot run viscous simulation without &
                                                                    &setting shear viscosity transport property object", 4534)

                if (.not. allocated(this%material(imat)%bulkvisc)) call GracefulExit("Cannot run viscous simulation without &
                                                                   &setting bulk viscosity transport property object", 4534)

                if (.not. allocated(this%material(imat)%thermcond)) call GracefulExit("Cannot run viscous simulation without &
                                                                    &setting thermal conductivity transport property object", 4534)
            end do

            if (this%ns > 1) then
                if (.not. allocated(this%massdiff)) call GracefulExit("Cannot run viscous multispecies simulation without &
                                                    &setting mass diffusivity transport property object", 4534)
            end if
        end if
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

    pure subroutine get_transport_properties(this, p, T, Ys, mu, bulk, kappa, diff)
        ! class(mixture), target,                                     intent(in)  :: this
        class(mixture),                                             intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: p, T
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(in)  :: Ys
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: mu, bulk, kappa
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(out) :: diff

        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: den, tmp, mu_i, Cp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns) :: Xs

        select case(this%inviscid)
        case(.true.)
            mu = zero
            bulk = zero
            kappa = zero
            diff = zero
        case default
            select case(this%ns)
            case(1)
                ! mat  => this%material(1)%mat
                ! visc => this%material(1)%visc
                Cp = this%material(1)%mat%gam * this%material(1)%mat%Rgas * this%material(1)%mat%onebygam_m1 ! Cp
                call this%material(1)%shearvisc%get_mu(T, mu)
                call this%material(1)%bulkvisc%get_beta(T, mu, bulk)
                call this%material(1)%thermcond%get_kappa(T, Cp, mu, kappa)
                diff = zero
            case default
                den = zero
                mu = zero
                bulk = zero
                kappa = zero
                do i=1,this%ns
                    ! mat  => this%material(i)%mat
                    ! visc => this%material(i)%visc

                    Cp = this%material(i)%mat%gam * this%material(i)%mat%Rgas * this%material(i)%mat%onebygam_m1 ! Cp
                    den = den + Ys(:,:,:,i)*sqrt(this%material(i)%mat%Rgas)

                    call this%material(i)%shearvisc%get_mu(T, mu_i)
                    mu = mu + mu_i*Ys(:,:,:,i)*sqrt(this%material(i)%mat%Rgas)

                    call this%material(i)%thermcond%get_kappa(T, Cp, mu_i, tmp)
                    kappa = kappa + tmp*Ys(:,:,:,i)*sqrt(this%material(i)%mat%Rgas)

                    call this%material(i)%bulkvisc%get_beta(T, mu_i, tmp)
                    bulk = bulk + tmp*Ys(:,:,:,i)*sqrt(this%material(i)%mat%Rgas)

                end do
                mu = mu / den
                bulk = bulk / den
                kappa = kappa / den

                do i = 1,this%ns
                    Xs(:,:,:,i) = this%material(i)%mat%Rgas * Ys(:,:,:,i) / this%Rgas
                end do
                call this%massdiff%get_diff(p, T, Xs, diff)
            end select
        end select
    end subroutine

    subroutine destroy(this)
        type(mixture), intent(inout)  :: this
        integer :: i

        do i = 1,this%ns
            if (allocated(this%material(i)%mat))     deallocate(this%material(i)%mat)
            if (allocated(this%material(i)%shearvisc)) deallocate(this%material(i)%shearvisc)
            if (allocated(this%material(i)%bulkvisc))  deallocate(this%material(i)%bulkvisc)
            if (allocated(this%material(i)%thermcond)) deallocate(this%material(i)%thermcond)
        end do
        if (allocated(this%massdiff)) deallocate(this%massdiff)

        if (allocated(this%gam )) deallocate(this%gam )
        if (allocated(this%Rgas)) deallocate(this%Rgas)
        if (allocated(this%Cv  )) deallocate(this%Cv  )
        if (allocated(this%Cp  )) deallocate(this%Cp  )

        deallocate(this%material)
    end subroutine


end module
