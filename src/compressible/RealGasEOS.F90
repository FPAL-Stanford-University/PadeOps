module RealGasEOS

    use kind_parameters, only: rkind
    use constants,       only: half,one
    use EOSMod_Real,     only: eos_real
    use decomp_2d,       only: decomp_info
    
    implicit none

    real(rkind) :: Ru = 8.3145 ! J/mol/K    Universal gas constant 
    
    ! Fitting for internal energy and temps
    real(rkind), dimension(7) :: an = (/2.356773520E+00,\
       8.984596770E-03,-7.123562690E-06,2.459190220E-09,\
       -1.436995480E-13,-4.837196970E+04,9.901052220E+00/)
    
    type, extends(eos_real) :: realgas

        real(rkind) :: Tc ! Critical temp [K]
        real(rkind) :: Pc ! Critical pressure [Pa]
        real(rkind) :: omega ! acentricity
        real(rkind) :: Rgas ! gas constant for this species
        real(rkind) :: MW ! molecular weight 
        real(rkind) :: PR_A, PR_B, PR_kappa ! peng-robinson constants
        real(rkind) :: b0,b1,b2 ! used to calculate K1 
        real(rkind), dimension(7) :: an ! fitting for energy 

        ! Some other buffers and     
        real(rkind), dimension(:,:,:), allocatable :: vm ! 1/rho
        real(rkind), dimension(:,:,:), allocatable :: Tr ! Tr = T/Tc
        real(rkind), dimension(:,:,:), allocatable :: Ti ! for temp iteration
        real(rkind), dimension(:,:,:), allocatable :: ei ! for temp iteration
        real(rkind), dimension(:,:,:), allocatable :: rhoi!for temp iteration
        real(rkind), dimension(:,:,:), allocatable :: alphaT ! thermal expansivity
        real(rkind), dimension(:,:,:), allocatable :: betaT  ! isothermal compressivity
        real(rkind), dimension(:,:,:), allocatable :: cp,cv,gam
        real(rkind), dimension(:,:,:), allocatable :: tmp1,tmp2,tmp3 ! buffers
    
    contains

        ! procedure :: init
        procedure :: get_p
        !procedure :: get_e_from_p
        procedure :: get_e_from_T
        procedure :: get_T
        procedure :: get_sos
        procedure :: get_enthalpy
        procedure :: update_props

    end type

    interface realgas
        module procedure init
    end interface

contains

    function init(decomp,MW,Tc,Pc,omega,T0,rho0) result(this)
        type(decomp_info), intent(in)   :: decomp
        real(rkind), intent(in) :: MW,Tc,Pc,omega
        real(rkind), intent(in) :: T0,rho0 ! initial T,rho for iters
        type(realgas) :: this
        integer :: nxp,nyp,nzp

        ! Properties
        this%Tc = Tc ! Critical temp [K]
        this%Pc = Pc ! Critical pressure [Pa]
        this%omega = omega ! acentricity
        this%Rgas = Ru/MW ! gas constant for this species
        this%MW = MW

        ! PR coefficients
        this%PR_A = 0.457236*(Ru*this%Tc)**2/this%Pc 
        this%PR_B = 0.0778*Ru*this%Tc/this%Pc
        this%PR_kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
        
        ! Used to calculate K1
        this%b0 = 8**(-0.5)/this%PR_B 
        this%b1 = (1.d0 - 2**0.5)*this%PR_B
        this%b2 = (1.d0 + 2**0.5)*this%PR_B
       
        ! Fitting for internal energy and temps
        this%an = (/2.356773520E+00,\
           8.984596770E-03,-7.123562690E-06,2.459190220E-09,\
           -1.436995480E-13,-4.837196970E+04,9.901052220E+00/)
        
        ! Allocate necessary fields
        nxp = decomp%ysz(1)
        nyp = decomp%ysz(2)
        nzp = decomp%ysz(3)
        if (allocated(this%vm))   deallocate(this%vm);   allocate(this%vm  (nxp,nyp,nzp))
        if (allocated(this%Tr))   deallocate(this%Tr);   allocate(this%Tr  (nxp,nyp,nzp))
        if (allocated(this%Ti))   deallocate(this%Ti);   allocate(this%Ti  (nxp,nyp,nzp))
        if (allocated(this%ei))   deallocate(this%ei);   allocate(this%ei  (nxp,nyp,nzp))
        if (allocated(this%rhoi)) deallocate(this%rhoi); allocate(this%rhoi(nxp,nyp,nzp))
        if (allocated(this%cp))   deallocate(this%cp);   allocate(this%cp  (nxp,nyp,nzp))
        if (allocated(this%cv))   deallocate(this%cv);   allocate(this%cv  (nxp,nyp,nzp))
        if (allocated(this%gam))  deallocate(this%gam);  allocate(this%gam (nxp,nyp,nzp))
        if (allocated(this%tmp1)) deallocate(this%tmp1); allocate(this%tmp1(nxp,nyp,nzp))
        if (allocated(this%tmp2)) deallocate(this%tmp2); allocate(this%tmp2(nxp,nyp,nzp))
        if (allocated(this%tmp3)) deallocate(this%tmp3); allocate(this%tmp3(nxp,nyp,nzp))

        ! Initialize density, temp and energy iteration values
        ! eventually overwritted with T,rho from previous timestep
        this%Ti = T0
        this%rhoi = rho0
    end function

    pure subroutine get_p(this,rho,T,p)
        class(realgas), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,T
        real(rkind), dimension(:,:,:), intent(out) :: p

        this%vm = this%MW/rho
        this%Tr = T/this%Tc
        this%tmp1 = this%PR_A * (1 + this%PR_kappa*(1-this%Tr**0.5))**2 ! alpha*a
        p = Ru*T/(this%vm-this%PR_B) - this%tmp1 / (this%vm**2 + 2*this%vm*this%PR_B - this%PR_B**2)
    end subroutine

    !pure subroutine get_e_from_p(this,rho,p,e)
    pure subroutine get_e_from_T(this,rho,T,e)
        class(realgas), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,T
        real(rkind), dimension(:,:,:), intent(out) :: e
        real(rkind), dimension(7) :: an
        real(rkind) :: b0,b1,b2

        ! ideal enthalpy h_ideal (store in e for now)
        an = this%an
        e = this%Rgas*T* (an(1) + an(2)*T/2d0 + an(3)*T**2d0/3d0 + \
            an(4)*T**3d0/4d0 + an(5)*T**4d0/5d0 + an(6)/T )
        
        ! Calculate K1
        this%vm = 1.d0/rho ! Use 1/rho instead of MW/rho since R=Ru/MW 
        this%tmp1 = b0*log( (this%vm+b1)/(this%vm+b2) ) !K1

        ! alpha*a and d(a*alpha)/dT
        this%Tr = T/this%Tc
        this%tmp2 = this%PR_A * (1 + this%PR_kappa*(1-this%Tr**0.5))**2 
        this%tmp3 = -1d0/T*this%tmp2 * \
            ( this%PR_kappa*this%Tr**0.5/(1+this%PR_kappa*(1-this%Tr**0.5)) )
        
        ! e = h_ideal - R*T + K1*[(a*alpha)-T*d(a*alpha)/dT] 
        e = e - this%Rgas*T + (this%tmp2-T*this%tmp3)*this%tmp1
    end subroutine

    !pure subroutine get_T(this,e,T)
    subroutine get_T(this,rho,e,T)
        class(realgas), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,e
        real(rkind), dimension(:,:,:), intent(out) :: T
        integer :: iter, max_iters=100, max_err=1d-6

        ! iteration guesses rhoi and Ti from previous calc
        T = this%Ti

        ! Try to match internal energy
        do iter = 1,max_iters 
            call get_e_from_T(this,rho,T,this%ei)
            if (maxval(abs((this%ei-e)/e)).lt.max_err) then
                this%rhoi = rho
                this%Ti = T
                exit 
            else 
                ! Update cp,cp,alphaT,betaT and p (store in tmp3)
                call update_props(this,rho,T)
                call get_p(this,rho,T,this%tmp1)

                ! Update temperature guess
                ! Ti = Ti + dv*dTdv + de*dTde
                this%tmp2 = rho/this%alphaT*(1d0-this%cp/this%cv) + this%tmp1/this%cv !dTdv
                this%tmp3 = this%MW/this%Cv !dTde
                T = T + ((1d0/rho-1d0/this%rhoi)*this%tmp2 + (e-this%ei)*this%tmp3)
            endif
            ! TODO: check if tmp3 should be MW/Cv or 1/Cv
        enddo

        if (maxval(abs((this%ei-e)/e)) .ge. max_err) then
            print *, 'Could not converge in get_T from e in 500 steps'
            pause
        endif
    end subroutine

    pure subroutine update_props(this,rho,T)
        class(realgas), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,T
        real(rkind), dimension(7) :: an
        real(rkind) :: c

        ! Updates cp,cv,alphaT,betaT

        ! Cp = dh/dT at constant p
        !    = dh_ideal/dT - R - T*d^2(a*alpha)/dT^2

        ! Ideal enthalphy derivative dh_ideal/dT (store in cp for now)
        this%cp = this%Rgas*(this%an(1) + this%an(2)*T + this%an(3)*T**2.d0 + this%an(4)*T**3.d0 + this%an(5)*T**4.d0 )
       
        ! Second derivative of d2(a*alpha)/dT^2 store in tmp2
        c = this%PR_kappa
        this%tmp2 = 0.457236*this%Rgas**2/T/2*c*(1+c)*this%Tc/this%Pc * sqrt(this%Tc/T)
        this%cp = this%cp - this%Rgas -T*this%tmp2
        
        ! Store tmp1 = d(alpha)/dT, alpha as in wikipedia
        !       tmp2 = dp/dT, constant vm
        ! Store tmp1 = a*alpha
        !       tmp3 = dp/dvm, constant T
        this%vm = 1d0/rho
        this%tmp1 = 2d0*( 1+c*(1-(T/this%Tc)**0.5)) * (-c/2/this%Tc**0.5/T**0.5)
        this%tmp2 = this%Rgas/(this%vm-this%PR_B) - \
            this%PR_A/(this%vm**2+2*this%vm*this%PR_B-this%PR_B**2) * this%tmp1
        this%tmp1 = this%PR_A * (1 + c*(1-(T/this%Tc)**0.5))**2 ! alpha*a
        this%tmp3 = -this%Rgas*T/(this%vm-this%PR_B)**2 + (2*this%vm + 2*this%PR_B)* \
            this%tmp1/(this%vm**2+2*this%vm*this%PR_B-this%PR_B**2)**2
             
        ! Thermal expansion coefficient
        ! alphaT = - dpdT / ( vm*dpdvm ) 
        this%alphaT = - this%tmp2 * rho / this%tmp3
        
        ! Isothermal compressibility 
        ! betaT = - 1 / ( vm*dpdT ) 
        this%betaT  = -rho / this%tmp3

        ! Specific heat constant volume
        ! Cv = Cp - alphaT**2 * vm * T / betaT
        this%cv = this%cp - this%alphaT**2 * T / rho / this%betaT
        this%gam = this%cp/this%cv
    end subroutine

    pure subroutine get_sos(this,rho,T,sos)
        class(realgas), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,T
        real(rkind), dimension(:,:,:), intent(out) :: sos
        
        call update_props(this,rho,T)
        sos = (rho*(this%betaT-T*this%alphaT**2/rho/this%Cp))**-0.5
    end subroutine

    !pure subroutine get_enthalpy(this,T,h)
    pure subroutine get_enthalpy(this,rho,p,T,h)
        class(realgas), intent(inout) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,p,T
        real(rkind), dimension(:,:,:), intent(out) :: h
        
        call get_e_from_T(this,rho,T,this%tmp1)
        h = this%tmp1 + p/rho
    end subroutine

end module
