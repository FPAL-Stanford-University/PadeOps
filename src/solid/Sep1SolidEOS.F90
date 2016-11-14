module Sep1SolidEOSMod

    use kind_parameters,      only: rkind,clen
    use constants,            only: one
    use decomp_2d,            only: decomp_info
    use DerivativesMod,       only: derivatives
    use FiltersMod,           only: filters
    use exits,                only: GracefulExit
    use AbstractEOSMod,       only: abstracteos
    use StiffGasEOS,          only: stiffgas
    use Sep1Solid_elasticMod, only: sep1solid_elastic

    implicit none

    type, extends(abstracteos) :: sep1solideos

        type(stiffgas)          :: hydro
        type(sep1solid_elastic) :: elastic

    contains

        procedure :: get_e_from_rho_g_T
        procedure :: get_p_devstress_T_sos2
        procedure :: get_pT_derivatives_wrt_energyVF
        procedure :: get_pT_from_energyVF
        final     :: destroy

    end type

    interface sep1solideos
        module procedure init
    end interface

contains

    function init(gam,Rgas,PInf,rho0,mu,yield,tau0,usegTg) result(this)
        type(sep1solideos)             :: this
        real(rkind),                 intent(in) :: gam, Rgas, PInf, rho0, mu, yield, tau0
        logical,                     intent(in) :: usegTg

        this%rho0 = rho0
        this%usegTg = usegTg

        call this%hydro%init(gam,Rgas,PInf)
        call this%elastic%init(this%rho0,mu,yield,tau0)
    end function

    pure elemental subroutine destroy(this)
        type(sep1solideos), intent(inout) :: this

        ! if (allocated(this%hydro)  ) deallocate(this%hydro)
        ! if (allocated(this%elastic)) deallocate(this%elastic)
    end subroutine

    ! Subroutine to get the internal energy from the material density, inverse deformation gradient tensor (of Finger tensor)
    ! and temperature
    subroutine get_e_from_rho_g_T(this,rho,g,T,e)
        class(sep1solideos),                                   intent(in)  :: this
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
        class(sep1solideos),                                     intent(in)  :: this
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

    subroutine get_pT_from_energyVF(this, VF0, g0, energy, VF, p, T)
        class(sep1solideos),       intent(in)  :: this
        real(rkind), dimension(9), intent(in)  :: g0
        real(rkind),               intent(in)  :: VF0, VF, energy
        real(rkind),               intent(out) :: p, T
        
        real(rkind) :: trG, trG2, detG, eelastic
        real(rkind) :: GG11, GG12, GG13, GG22, GG23, GG33

        if (this%usegTg) then
            GG11 = g0(1); GG12 = g0(2); GG13 = g0(3);
                          GG22 = g0(4); GG23 = g0(5);
                                        GG33 = g0(6);
        else
            GG11 = g0(1)*g0(1) + g0(4)*g0(4) + g0(7)*g0(7)
            GG12 = g0(1)*g0(2) + g0(4)*g0(5) + g0(7)*g0(8)
            GG13 = g0(1)*g0(3) + g0(4)*g0(6) + g0(7)*g0(9)
            GG22 = g0(2)*g0(2) + g0(5)*g0(5) + g0(8)*g0(8)
            GG23 = g0(2)*g0(3) + g0(5)*g0(6) + g0(8)*g0(9)
            GG33 = g0(3)*g0(3) + g0(6)*g0(6) + g0(9)*g0(9)
        end if
        
        call this%elastic%get_invariants_elemental(GG11,GG12,GG13,GG22,GG23,GG33,trG,trG2,detG)
        call this%elastic%get_eelastic(trG,trG2,detG,eelastic)

        detG = sqrt(detG) ! This is really det(g0) now

        p = (this%hydro%gam - one)*this%rho0*VF0*detG/VF*(energy-eelastic) &
          - this%hydro%gam * this%hydro%PInf
        
        T = (energy - eelastic - this%hydro%PInf*VF/(this%rho0*VF0*detG))/this%hydro%Cv

        !write(*,*) 'VF0 = ', VF0
        !write(*,*) 'VF = ', VF
        !write(*,'(a,9(e19.12,1x))') 'g0 = ', g0
        !write(*,*) 'energy = ', energy
        !write(*,*) 'p = ', p
        !write(*,*) 'T = ', T

    end subroutine

    subroutine get_pT_derivatives_wrt_energyVF(this, VF0, g0, energy, VF, dpde, dpdVF, dTde, dTdVF)
        class(sep1solideos),       intent(in)  :: this
        real(rkind), dimension(9), intent(in)  :: g0
        real(rkind),               intent(in)  :: VF0, VF, energy
        real(rkind),               intent(out) :: dpde, dpdVF, dTde, dTdVF
        
        real(rkind) :: trG, trG2, detG, eelastic
        real(rkind) :: GG11, GG12, GG13, GG22, GG23, GG33
        
        if (this%usegTg) then
            GG11 = g0(1); GG12 = g0(2); GG13 = g0(3);
                          GG22 = g0(4); GG23 = g0(5);
                                        GG33 = g0(6);
        else
            GG11 = g0(1)*g0(1) + g0(4)*g0(4) + g0(7)*g0(7)
            GG12 = g0(1)*g0(2) + g0(4)*g0(5) + g0(7)*g0(8)
            GG13 = g0(1)*g0(3) + g0(4)*g0(6) + g0(7)*g0(9)
            GG22 = g0(2)*g0(2) + g0(5)*g0(5) + g0(8)*g0(8)
            GG23 = g0(2)*g0(3) + g0(5)*g0(6) + g0(8)*g0(9)
            GG33 = g0(3)*g0(3) + g0(6)*g0(6) + g0(9)*g0(9)
        end if

        call this%elastic%get_invariants_elemental(GG11,GG12,GG13,GG22,GG23,GG33,trG,trG2,detG)
        call this%elastic%get_eelastic(trG,trG2,detG,eelastic)

        detG = sqrt(detG) ! This is really det(g0) now

        dpde = (this%hydro%gam - one)*this%rho0*VF0*detG / VF

        dpdVF = - (this%hydro%gam - one)*this%rho0*VF0*detG*(energy - eelastic) / VF**2

        dTde = one / this%hydro%Cv

        dTdVF = -this%hydro%PInf / (this%hydro%Cv*this%rho0*VF0*detG)
    end subroutine

end module
