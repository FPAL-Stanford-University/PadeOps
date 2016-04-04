module MixtureEOSMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,one
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use EOSMod,          only: eos
    use IdealGasEOS,     only: idealgas

    integer, parameter :: max_materials = 99

    type :: material_eos
        class(idealgas), allocatable :: mat
    end type

    type, abstract :: mixture

        integer :: ns
        integer :: nxp, nyp, nzp
        type(material_eos), dimension(:), allocatable :: material
        real(rkind), dimension(:,:,:), allocatable :: gam, Rgas, Cp, Cv

    contains

        procedure :: init
        procedure :: update
        procedure :: get_p
        procedure :: get_T
        procedure :: get_e_from_p
        procedure :: get_sos
        procedure :: destroy

    end type

contains

    subroutine init(this,decomp,ns,gam,Rgas)
        class(mixture),                                intent(inout) :: this
        type(decomp_info),                             intent(in)    :: decomp
        integer,                                       intent(in)    :: ns
        real(rkind), dimension(max_materials),         intent(in)    :: gam, Rgas

        integer :: i

        if (ns .GT. max_materials) call GracefulExit('ns exceeds max_materials.', 3457)

        this%ns = ns
        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        allocate(this%material(this%ns))

        do i = 1,this%ns
            !if (allocated(this%material(i)%mat)) deallocate(this%material(i)%mat); 
            allocate(idealgas :: this%material(i)%mat)
            call this%material(i)%mat%init(gam(i), Rgas(i))
        end do

        if (allocated(this%gam )) deallocate(this%gam ); allocate(this%gam (this%nxp,this%nyp,this%nzp))
        if (allocated(this%Rgas)) deallocate(this%Rgas); allocate(this%Rgas(this%nxp,this%nyp,this%nzp))
        if (allocated(this%Cv  )) deallocate(this%Cv  ); allocate(this%Cv  (this%nxp,this%nyp,this%nzp))
        if (allocated(this%Cp  )) deallocate(this%Cp  ); allocate(this%Cp  (this%nxp,this%nyp,this%nzp))

    end subroutine

    pure subroutine update(this,Ys)
        class(mixture), target,                                     intent(inout) :: this 
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(in)    :: Ys
        
        class(idealgas), pointer :: mat
        integer :: i
        
        this%Rgas = zero
        this%Cp = zero
        do i = 1,this%ns
            mat => this%material(i)%mat
            this%Rgas = this%Rgas + (mat%Rgas) * Ys(:,:,:,i)
            this%Cp = this%Cp + (mat%gam * mat%Rgas * mat%onebygam_m1 ) * Ys(:,:,:,i) ! Cp_i = gam/(gam-1) * Rgas
        end do

        this%gam = this%Cp / (this%Cp - this%Rgas)
        this%Cv = this%Cp / this%gam

    end subroutine

    pure subroutine get_p(this,rho,e,p)
        class(mixture),                                             intent(in)  :: this 
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: p

        p = (this%gam-one)*rho*e

    end subroutine

    pure subroutine get_T(this,e,T)
        class(mixture),                                             intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: T

        T = e / this%Cv
    end subroutine

    pure subroutine get_e_from_p(this,rho,p,e)
        class(mixture),                                             intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,p
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: e

        e = p / ((this%gam - one)*rho)
    end subroutine

    pure subroutine get_sos(this,rho,p,sos)
        class(mixture),                                             intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,p
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: sos

        sos = sqrt(this%gam*p/rho)
    end subroutine

    subroutine destroy(this)
        class(mixture), intent(inout)  :: this
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
