module Sep1SolidEOSMod

    use kind_parameters, only: rkind,clen
    use constants,       only: one
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use AbstractEOSMod,  only: abstracteos
    use StiffGasEOS,     only: stiffgas
    use Sep1SolidEOS,    only: sep1solid

    implicit none

    type, extends(abstracteos) :: sep1solideos_abstract

        type(stiffgas) , allocatable :: hydro
        type(sep1solid), allocatable :: elastic

    contains

        procedure :: get_e_from_rho_g_T
        procedure :: get_p_devstress_T_sos2
        final     :: destroy

    end type

    interface sep1solideos_abstract
        module procedure init
    end interface

contains

    function init(gam,Rgas,PInf,rho0,mu,yield,tau0,usegTg) result(this)
        type(sep1solideos_abstract)             :: this
        real(rkind),                 intent(in) :: gam, Rgas, PInf, rho0, mu, yield, tau0
        logical,                     intent(in) :: usegTg

        this%rho0 = rho0
        this%usegTg = usegTg

        allocate(this%hydro,source=stiffgas(gam,Rgas,PInf))
        allocate(this%elastic,source=sep1solid(this%rho0,mu,yield,tau0))
    end function

    pure elemental subroutine destroy(this)
        type(sep1solideos_abstract), intent(inout) :: this

        if (allocated(this%hydro)  ) deallocate(this%hydro)
        if (allocated(this%elastic)) deallocate(this%elastic)
    end subroutine

    ! Subroutine to get the internal energy from the material density, inverse deformation gradient tensor (of Finger tensor)
    ! and temperature
    subroutine get_e_from_rho_g_T(this,rho,g,T,e)
        class(sep1solideos_abstract),                                   intent(in)  :: this
        real(rkind), dimension(:,:,:,:),                       intent(in)  :: g
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)), intent(in)  :: rho,T
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)), intent(out) :: e
        
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)  ) :: trG, trG2, detG, eelastic
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6) :: finger, fingersq
    
        call this%elastic%get_finger(g,finger,fingersq,trG,trG2,detG,this%usegTg)
        call this%elastic%get_eelastic(trG,trG2,detG,eelastic)

        call this%hydro%get_e_from_T(rho,T,e)

        e = e + eelastic
    end subroutine

    subroutine get_p_devstress_T_sos2(this,g,rho,e,p,T,devstress,sos_sq)
        class(sep1solideos_abstract),                                     intent(in)  :: this
        real(rkind), dimension(:,:,:,:),                         intent(in)  :: g
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)),   intent(in)  :: rho, e
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)),   intent(out) :: p, T, sos_sq
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6), intent(out) :: devstress
    
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)  ) :: trG, trG2, detG, eelastic
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6) :: finger, fingersq
    
        call this%elastic%get_finger(g,finger,fingersq,trG,trG2,detG,this%usegTg)
        call this%elastic%get_eelastic(trG,trG2,detG,eelastic)

        eelastic = e - eelastic ! Contains ehydro now

        call this%hydro%get_p(rho,eelastic,p)
        call this%hydro%get_T(eelastic,T,rho)

        call this%elastic%get_devstress(finger,fingersq,trG,trG2,detG,devstress)

        call this%hydro%get_sos2(rho,p,sos_sq)
        call this%elastic%get_sos2(rho,sos_sq) ! Add elastic component to sos^2
    end subroutine

end module
