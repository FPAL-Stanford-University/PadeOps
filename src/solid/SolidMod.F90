module SolidMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,one,two,epssmall,three
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use EOSMod,          only: eos
    use StiffGasEOS,     only: stiffgas
    use Sep1SolidEOS,    only: sep1solid

    implicit none

    type :: solid
        integer :: nxp, nyp, nzp

        logical :: plast = .FALSE.
        logical :: explPlast = .FALSE.
        logical :: PTeqb = .TRUE.
        logical :: use_gTg = .FALSE.

        class(stiffgas ), allocatable :: hydro
        class(sep1solid), allocatable :: elastic

        type(decomp_info), pointer :: decomp
        type(derivatives), pointer :: der
        type(filters),     pointer :: fil

        real(rkind), dimension(:,:,:), allocatable :: Ys
        real(rkind), dimension(:,:,:), allocatable :: VF
        real(rkind), dimension(:,:,:), allocatable :: eh
        real(rkind), dimension(:,:,:), allocatable :: eel

        real(rkind), dimension(:,:,:,:), allocatable :: g
        real(rkind), dimension(:,:,:),   pointer     :: g11
        real(rkind), dimension(:,:,:),   pointer     :: g12
        real(rkind), dimension(:,:,:),   pointer     :: g13
        real(rkind), dimension(:,:,:),   pointer     :: g21
        real(rkind), dimension(:,:,:),   pointer     :: g22
        real(rkind), dimension(:,:,:),   pointer     :: g23
        real(rkind), dimension(:,:,:),   pointer     :: g31
        real(rkind), dimension(:,:,:),   pointer     :: g32
        real(rkind), dimension(:,:,:),   pointer     :: g33
        
        real(rkind), dimension(:,:,:,:), allocatable :: devstress
        real(rkind), dimension(:,:,:),   pointer     :: sxx
        real(rkind), dimension(:,:,:),   pointer     :: sxy
        real(rkind), dimension(:,:,:),   pointer     :: sxz
        real(rkind), dimension(:,:,:),   pointer     :: syy
        real(rkind), dimension(:,:,:),   pointer     :: syz
        real(rkind), dimension(:,:,:),   pointer     :: szz
        
        real(rkind), dimension(:,:,:),   allocatable :: p
        real(rkind), dimension(:,:,:),   allocatable :: T

        ! species-specific artificial properties
        real(rkind), dimension(:,:,:),   allocatable :: kap
        real(rkind), dimension(:,:,:,:), allocatable :: qi
        real(rkind), dimension(:,:,:),   allocatable :: diff
        real(rkind), dimension(:,:,:,:), allocatable :: Ji

        ! species-specific conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: consrv

        ! work arrays
        real(rkind), dimension(:,:,:,:), allocatable :: Qtmpg
        real(rkind), dimension(:,:,:),   allocatable :: QtmpYs
        real(rkind), dimension(:,:,:),   allocatable :: Qtmpeh
        real(rkind), dimension(:,:,:),   allocatable :: QtmpVF

        real(rkind), dimension(:,:,:), allocatable :: modDevSigma
    contains

        procedure :: init
        procedure :: getRHS_g
        procedure :: getRHS_gTg
        procedure :: getRHS_Ys
        procedure :: getRHS_eh
        procedure :: getRHS_VF
        procedure :: update_g
        procedure :: update_gTg
        procedure :: update_Ys
        procedure :: update_eh
        procedure :: update_VF
        procedure :: getPhysicalProperties
        procedure :: getPlasticSources
        procedure :: get_p_from_ehydro
        procedure :: get_ehydroT_from_p
        procedure :: get_eelastic_devstress
        procedure :: get_conserved
        procedure :: get_primitive
        procedure :: getSpeciesDensity
        procedure :: getSpeciesDensity_from_g
        procedure :: get_enthalpy
        procedure :: checkNaN
        procedure :: filter
        final     :: destroy

    end type

    ! interface solid
    !     module procedure init
    ! end interface

    ! hooks for sources

    interface hook_material_g_source
        subroutine hook_material_g_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
            import :: rkind
            import :: decomp_info
            import :: stiffgas
            import :: sep1solid
            type(decomp_info),               intent(in)    :: decomp
            type(stiffgas),                  intent(in)    :: hydro
            type(sep1solid),                 intent(in)    :: elastic
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
            real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
            real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
        end subroutine
    end interface

    interface hook_material_mass_source
        subroutine hook_material_mass_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
            import :: rkind
            import :: decomp_info
            import :: stiffgas
            import :: sep1solid
            type(decomp_info),               intent(in)    :: decomp
            type(stiffgas),                  intent(in)    :: hydro
            type(sep1solid),                 intent(in)    :: elastic
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
            real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
            real(rkind), dimension(:,:,:),   intent(inout) :: rhs
        end subroutine
    end interface

    interface hook_material_VF_source
        subroutine hook_material_VF_source(decomp,hydro,elastic,x,y,z,tsim,u,v,w,Ys,VF,p,rhs)
            import :: rkind
            import :: decomp_info
            import :: stiffgas
            import :: sep1solid
            type(decomp_info),               intent(in)    :: decomp
            type(stiffgas),                  intent(in)    :: hydro
            type(sep1solid),                 intent(in)    :: elastic
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
            real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w,Ys,VF,p
            real(rkind), dimension(:,:,:),   intent(inout) :: rhs
        end subroutine
    end interface

    interface hook_material_energy_source
        subroutine hook_material_energy_source(decomp,hydro,elastic,x,y,z,tsim,rho,u,v,w,Ys,VF,p,rhs)
            import :: rkind
            import :: decomp_info
            import :: stiffgas
            import :: sep1solid
            type(decomp_info),               intent(in)    :: decomp
            type(stiffgas),                  intent(in)    :: hydro
            type(sep1solid),                 intent(in)    :: elastic
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:),   intent(in)    :: x,y,z
            real(rkind), dimension(:,:,:),   intent(in)    :: rho,u,v,w,Ys,VF,p
            real(rkind), dimension(:,:,:),   intent(inout) :: rhs
        end subroutine
    end interface

contains

    !function init(decomp,der,fil,hydro,elastic) result(this)
    subroutine init(this,decomp,der,fil,PTeqb,use_gTg)
        class(solid), target, intent(inout) :: this
        type(decomp_info), target, intent(in) :: decomp
        type(derivatives), target, intent(in) :: der
        type(filters),     target, intent(in) :: fil
        logical, intent(in) :: PTeqb
        logical, intent(in) :: use_gTg

        this%decomp => decomp
        this%der => der
        this%fil => fil
       
        this%PTeqb = PTeqb

        this%use_gTg = use_gTg

        ! Assume everything is in Y decomposition
        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        if (allocated(this%hydro)) deallocate(this%hydro)
        allocate( this%hydro )
        
        if (allocated(this%elastic)) deallocate(this%elastic)
        allocate( this%elastic )
        
        ! Allocate material massfraction
        if( allocated( this%Ys ) ) deallocate( this%Ys )
        allocate( this%Ys(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material volume fraction
        if( allocated( this%VF ) ) deallocate( this%VF )
        allocate( this%VF(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material hydrodynamic energy
        if( allocated( this%eh ) ) deallocate( this%eh )
        allocate( this%eh(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material elastic energy
        if( allocated( this%eel ) ) deallocate( this%eel )
        allocate( this%eel(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material inverse elastic deformation gradients and associate pointers
        if( allocated( this%g ) ) deallocate( this%g )
        allocate( this%g(this%nxp,this%nyp,this%nzp,9) )
        this%g11 => this%g(:,:,:,1)   
        this%g12 => this%g(:,:,:,2)   
        this%g13 => this%g(:,:,:,3)   
        this%g21 => this%g(:,:,:,4)   
        this%g22 => this%g(:,:,:,5)   
        this%g23 => this%g(:,:,:,6)   
        this%g31 => this%g(:,:,:,7)   
        this%g32 => this%g(:,:,:,8)   
        this%g33 => this%g(:,:,:,9)   

        ! Allocate material deviatoric stress array
        if( allocated( this%devstress ) ) deallocate( this%devstress )
        allocate( this%devstress(this%nxp,this%nyp,this%nzp,6) )
        this%sxx  => this%devstress(:,:,:,1)   
        this%sxy  => this%devstress(:,:,:,2)   
        this%sxz  => this%devstress(:,:,:,3)   
        this%syy  => this%devstress(:,:,:,4)   
        this%syz  => this%devstress(:,:,:,5)   
        this%szz  => this%devstress(:,:,:,6)   
        
        ! Allocate mrray to store contraction of deviatoric stress
        if( allocated( this%modDevSigma ) ) deallocate( this%modDevSigma )
        allocate( this%modDevSigma(this%nxp,this%nyp,this%nzp) )
        
        ! Allocate material pressure array
        if( allocated( this%p ) ) deallocate( this%p )
        allocate( this%p(this%nxp,this%nyp,this%nzp) )

        ! Allocate material temperature array
        if( allocated( this%T ) ) deallocate( this%T )
        allocate( this%T(this%nxp,this%nyp,this%nzp) )

        ! Allocate material kappa array
        if( allocated( this%kap ) ) deallocate( this%kap )
        allocate( this%kap(this%nxp,this%nyp,this%nzp) )

        ! Allocate material diffusive flux
        if( allocated( this%qi ) ) deallocate( this%qi )
        allocate( this%qi(this%nxp,this%nyp,this%nzp,3) )
        
        ! Allocate material diffusivity array
        if( allocated( this%diff ) ) deallocate( this%diff )
        allocate( this%diff(this%nxp,this%nyp,this%nzp) )

        ! Allocate material diffusive flux
        if( allocated( this%Ji ) ) deallocate( this%Ji )
        allocate( this%Ji(this%nxp,this%nyp,this%nzp,3) )
        
        ! Allocate material conserved variables
        if( allocated( this%consrv ) ) deallocate( this%consrv )
        if(this%PTeqb) then
            allocate( this%consrv(this%nxp,this%nyp,this%nzp,1) )
        else
            allocate( this%consrv(this%nxp,this%nyp,this%nzp,2) )
        endif

        ! Allocate work arrays
        ! g tensor equation
        if( allocated( this%Qtmpg ) ) deallocate( this%Qtmpg )
        allocate( this%Qtmpg(this%nxp,this%nyp,this%nzp,9) )

        ! Ys equation
        if( allocated( this%QtmpYs ) ) deallocate( this%QtmpYs )
        allocate( this%QtmpYs(this%nxp,this%nyp,this%nzp) )

        if(.NOT. this%PTeqb) then
            ! eh equation
            if( allocated( this%Qtmpeh ) ) deallocate( this%Qtmpeh )
            allocate( this%Qtmpeh(this%nxp,this%nyp,this%nzp) )

            ! VF equation
            if( allocated( this%QtmpVF ) ) deallocate( this%QtmpVF )
            allocate( this%QtmpVF(this%nxp,this%nyp,this%nzp) )
        endif
    !end function
    end subroutine

    pure elemental subroutine destroy(this)
        type(solid), intent(inout) :: this

        ! First deallocate all the arrays
        if( allocated( this%QtmpVF ) ) deallocate( this%QtmpVF )
        if( allocated( this%Qtmpeh ) ) deallocate( this%Qtmpeh )
        if( allocated( this%QtmpYs ) ) deallocate( this%QtmpYs )
        if( allocated( this%Qtmpg  ) ) deallocate( this%Qtmpg )
        if( allocated( this%consrv ) ) deallocate( this%consrv )
        
        if( allocated( this%Ji )   ) deallocate( this%Ji )
        if( allocated( this%diff ) ) deallocate( this%diff )
        if( allocated( this%kap )  ) deallocate( this%kap )
        if( allocated( this%T )    ) deallocate( this%T )
        if( allocated( this%p )    ) deallocate( this%p )
        if( allocated( this%modDevSigma ) ) deallocate( this%modDevSigma )

        nullify( this%sxx ); nullify( this%sxy ); nullify( this%sxz )
                             nullify( this%syy ); nullify( this%syz )
                                                  nullify( this%szz )
        if( allocated( this%devstress ) ) deallocate( this%devstress )

        nullify( this%g11 ); nullify( this%g12 ); nullify( this%g13 )
        nullify( this%g21 ); nullify( this%g22 ); nullify( this%g23 )
        nullify( this%g31 ); nullify( this%g32 ); nullify( this%g33 )
        if( allocated( this%g )         ) deallocate( this%g )

        if( allocated( this%eel ) ) deallocate( this%eel )
        if( allocated( this%eh )  ) deallocate( this%eh )
        if( allocated( this%VF )  ) deallocate( this%VF )
        if( allocated( this%Ys )  ) deallocate( this%Ys )

        ! Now deallocate the EOS objects
        if ( allocated(this%hydro)   ) deallocate(this%hydro)
        if ( allocated(this%elastic) ) deallocate(this%elastic)

        nullify( this%fil    )
        nullify( this%der    )
        nullify( this%decomp )
    end subroutine

    subroutine getPhysicalProperties(this)
        class(solid), intent(inout) :: this

        this%diff = zero
        this%kap = zero

    end subroutine

    subroutine getSpeciesDensity_from_g(this,rho)
        class(solid), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(out)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: detg

        detg = this%g11*(this%g22*this%g33-this%g23*this%g32) &
             - this%g12*(this%g21*this%g33-this%g31*this%g23) &
             + this%g13*(this%g21*this%g32-this%g31*this%g22)

        if (this%use_gTg) then
            detg = sqrt(detg)
        end if

        rho = detg*this%elastic%rho0

    end subroutine

    subroutine update_g(this,isub,dt,rho,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg  ! RHS for g tensor equation
        real(rkind) :: max_modDevSigma

        call this%getRHS_g(rho,u,v,w,dt,rhsg,x_bc,y_bc,z_bc)
        call hook_material_g_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsg)

        ! advance sub-step
        if(isub==1) this%Qtmpg = zero                   ! not really needed, since RK45_A(1) = 0
        this%Qtmpg  = dt*rhsg + RK45_A(isub)*this%Qtmpg
        this%g = this%g  + RK45_B(isub)*this%Qtmpg

        ! Now project g tensor to SPD space
        call this%elastic%make_tensor_SPD(this%g)

        if(this%plast) then
            if (.NOT. this%explPlast) then
                ! Effect plastic deformations
                call this%get_eelastic_devstress()
                this%modDevSigma = this%sxx*this%sxx + this%syy*this%syy + this%szz*this%szz + &
                        two*( this%sxy*this%sxy + this%sxz*this%sxz + this%syz*this%syz )
                max_modDevSigma = sqrt(3.0d0/2.0d0*maxval(this%modDevSigma))
                !if(max_modDevSigma > this%elastic%yield) then
                !  write(*,'(a,2(e19.12,1x))') 'Entering plasticity. Stresses = ', max_modDevSigma, this%elastic%yield
                !endif
                if(this%PTeqb) this%kap = sqrt(two/three*this%modDevSigma)

                call this%elastic%plastic_deformation(this%g, this%use_gTg)

                call this%get_eelastic_devstress()
                this%modDevSigma = this%sxx*this%sxx + this%syy*this%syy + this%szz*this%szz + &
                        two*( this%sxy*this%sxy + this%sxz*this%sxz + this%syz*this%syz )
                max_modDevSigma = sqrt(3.0d0/2.0d0*maxval(this%modDevSigma))
                !if(max_modDevSigma > this%elastic%yield) then
                !  write(*,'(a,2(e19.12,1x))') 'Exiting plasticity. Stresses = ', max_modDevSigma, this%elastic%yield
                !endif
                
            end if
        end if

    end subroutine

    subroutine getRHS_g(this,rho,u,v,w,dt,rhsg,x_bc,y_bc,z_bc)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl
        use reductions, only: P_MAXLOC
        class(solid),                                         intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg
        real(rkind), parameter :: etafac = one/6._rkind

        real(rkind) :: pmax
        integer :: imax, jmax, kmax, rmax

        ! Symmetry and anti-symmetry properties of g are assumed as below
        ! In x g_{ij}: [S A A; A S S; A S S]
        ! In y g_{ij}: [S A S; A S A; S A S]
        ! In z g_{ij}: [S S A; S S A; A A S]

        rhsg = zero

        detg = this%g11*(this%g22*this%g33-this%g23*this%g32) &
             - this%g12*(this%g21*this%g33-this%g31*this%g23) &
             + this%g13*(this%g21*this%g32-this%g31*this%g22)

        ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
        tmp = (rho*this%Ys + this%elastic%rho0*detg*epssmall)/(this%VF + epssmall)   
        ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF
        penalty = etafac*( tmp/detg/this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
        if (this%elastic%mu < eps) penalty = zero

        tmp = -u*this%g11-v*this%g12-w*this%g13
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,1),rhsg(:,:,:,2),rhsg(:,:,:,3),-x_bc, y_bc, z_bc)
        
        call curl(this%decomp, this%der, this%g11, this%g12, this%g13, curlg, -x_bc, y_bc, z_bc)
        ! call P_MAXLOC( abs(penalty*this%g11), pmax, imax, jmax, kmax, rmax)
        ! if (nrank == rmax) then
        !     print*, "Maximum penalty term = ", pmax
        !     print*, "        curl    term = ", v(imax,jmax,kmax)*curlg(imax,jmax,kmax,3) - w(imax,jmax,kmax)*curlg(imax,jmax,kmax,2)
        !     print*, "        g11          = ", this%g11(imax,jmax,kmax)
        !     print*, "        VF           = ", this%VF(imax,jmax,kmax)
        ! end if
        ! call P_MAXLOC( abs(v*curlg(:,:,:,3) - w*curlg(:,:,:,2)), pmax, imax, jmax, kmax, rmax)
        ! if (nrank == rmax) then
        !     print*, "Maximum curl    term = ", pmax
        !     print*, "        penalty term = ", penalty(imax,jmax,kmax)*this%g11(imax,jmax,kmax)
        !     print*, "        g11          = ", this%g11(imax,jmax,kmax)
        !     print*, "        VF           = ", this%VF(imax,jmax,kmax)
        !     print*, ""
        ! end if
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g11
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g12
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g13
 
        tmp = -u*this%g21-v*this%g22-w*this%g23
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,4),rhsg(:,:,:,5),rhsg(:,:,:,6), x_bc,-y_bc, z_bc)   
        
        call curl(this%decomp, this%der, this%g21, this%g22, this%g23, curlg, x_bc, -y_bc, z_bc)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g21
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g22
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g23
 
        tmp = -u*this%g31-v*this%g32-w*this%g33
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,7),rhsg(:,:,:,8),rhsg(:,:,:,9), x_bc, y_bc,-z_bc)

        call curl(this%decomp, this%der, this%g31, this%g32, this%g33, curlg, x_bc, y_bc, -z_bc)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g31
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g32
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g33

        if (this%plast) then
            if(this%explPlast) then
                call this%getPlasticSources(detg,rhsg)
            end if
        end if

    end subroutine

    subroutine update_gTg(this,isub,dt,rho,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg  ! RHS for gTg tensor equation

        call this%getRHS_gTg(rho,u,v,w,dt,rhsg,x_bc,y_bc,z_bc)
        call hook_material_g_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsg)

        ! advance sub-step
        if(isub==1) this%Qtmpg = zero                   ! not really needed, since RK45_A(1) = 0
        this%Qtmpg  = dt*rhsg + RK45_A(isub)*this%Qtmpg
        this%g = this%g  + RK45_B(isub)*this%Qtmpg

        if(this%plast) then
            if (.NOT. this%explPlast) then
                call this%elastic%plastic_deformation(this%g, this%use_gTg)
            end if
        end if

    end subroutine

    subroutine getRHS_gTg(this,rho,u,v,w,dt,rhsg,x_bc,y_bc,z_bc)
        use operators, only: gradient, curl
        class(solid),                                         intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradG
        real(rkind), parameter :: etafac = one/6._rkind

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        call gradient(this%decomp, this%der, u, dudx, dudy, dudz, -x_bc,  y_bc,  z_bc)
        call gradient(this%decomp, this%der, v, dvdx, dvdy, dvdz,  x_bc, -y_bc,  z_bc)
        call gradient(this%decomp, this%der, w, dwdx, dwdy, dwdz,  x_bc,  y_bc, -z_bc)

        ! Symmetry and anti-symmetry properties of gTg are assumed as below (same as g)
        ! In x gTg_{ij}: [S A A; A S S; A S S]
        ! In y gTg_{ij}: [S A S; A S A; S A S]
        ! In z gTg_{ij}: [S S A; S S A; A A S]

        rhsg = zero

        detG = this%G11*(this%G22*this%G33-this%G23*this%G32) &
             - this%G12*(this%G21*this%G33-this%G31*this%G23) &
             + this%G13*(this%G21*this%G32-this%G31*this%G22)

        ! penalty = etafac*(this%rho/sqrt(detG)/this%rho0-one)
        
        detg = sqrt(detg)
        tmp = (rho*this%Ys + this%elastic%rho0*detg*epssmall)/(this%VF + epssmall) ! species density 
        penalty = etafac*( tmp/detg/this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density

        call gradient(this%decomp, this%der, this%G11, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,1) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
                        - (dudx*this%G11 + dvdx*this%G21 + dwdx*this%G31) &
                        - (this%G11*dudx + this%G12*dvdx + this%G13*dwdx) + penalty*this%G11

        call gradient(this%decomp, this%der, this%G12, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3),-x_bc,-y_bc, z_bc)
        rhsg(:,:,:,2) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
                        - (dudx*this%G12 + dvdx*this%G22 + dwdx*this%G32) &
                        - (this%G11*dudy + this%G12*dvdy + this%G13*dwdy) + penalty*this%G12

        call gradient(this%decomp, this%der, this%G13, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3),-x_bc, y_bc,-z_bc)
        rhsg(:,:,:,3) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
                        - (dudx*this%G13 + dvdx*this%G23 + dwdx*this%G33) &
                        - (this%G11*dudz + this%G12*dvdz + this%G13*dwdz) + penalty*this%G13

        rhsg(:,:,:,4) = rhsg(:,:,:,2)  ! Since symmetric

        call gradient(this%decomp, this%der, this%G22, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,5) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
                        - (dudy*this%G12 + dvdy*this%G22 + dwdy*this%G32) &
                        - (this%G21*dudy + this%G22*dvdy + this%G23*dwdy) + penalty*this%G22

        call gradient(this%decomp, this%der, this%G23, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3), x_bc,-y_bc,-z_bc)
        rhsg(:,:,:,6) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
                        - (dudy*this%G13 + dvdy*this%G23 + dwdy*this%G33) &
                        - (this%G21*dudz + this%G22*dvdz + this%G23*dwdz) + penalty*this%G23

        rhsg(:,:,:,7) = rhsg(:,:,:,3)  ! Since symmetric
        rhsg(:,:,:,8) = rhsg(:,:,:,6)  ! Since symmetric

        call gradient(this%decomp, this%der, this%G33, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,9) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
                        - (dudz*this%G13 + dvdz*this%G23 + dwdz*this%G33) &
                        - (this%G31*dudz + this%G32*dvdz + this%G33*dwdz) + penalty*this%G33

        if (this%plast) then
            if (this%explPlast) then
                call this%getPlasticSources(detg,rhsg)
            end if
        end if

    end subroutine

    subroutine getPlasticSources(this,detg,rhsg)
        use constants, only: twothird
        class(solid),                                         intent(in)    :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)    :: detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(inout) :: rhsg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: invtaurel

        ! Get S'S'
        invtaurel = this%sxx*this%sxx + two*this%sxy*this%sxy + two*this%sxz*this%sxz &
                                      +     this%syy*this%syy + two*this%syz*this%syz &
                                                              +     this%szz*this%szz

        ! 1/tau_rel
        invtaurel = (one / this%elastic%tau0) * ( invtaurel - (twothird)*this%elastic%yield**2 ) / this%elastic%mu**2
        where (invtaurel .LE. zero)
            invtaurel = zero
        end where
        invtaurel = invtaurel / (two * this%elastic%mu * detg)

        if (this%use_gTg) then
            invtaurel = two*invtaurel  ! Factor of 2 for gTg implementation
        end if

        ! Add (1/tau_rel)*g*S to the rhsg (explicit plastic source terms)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + invtaurel * ( this%g11*this%sxx + this%g12*this%sxy + this%g13*this%sxz ) ! g11 
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + invtaurel * ( this%g11*this%sxy + this%g12*this%syy + this%g13*this%syz ) ! g12 
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + invtaurel * ( this%g11*this%sxz + this%g12*this%syz + this%g13*this%szz ) ! g13 
 
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + invtaurel * ( this%g21*this%sxx + this%g22*this%sxy + this%g23*this%sxz ) ! g21 
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + invtaurel * ( this%g21*this%sxy + this%g22*this%syy + this%g23*this%syz ) ! g22 
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + invtaurel * ( this%g21*this%sxz + this%g22*this%syz + this%g23*this%szz ) ! g23 

        rhsg(:,:,:,7) = rhsg(:,:,:,7) + invtaurel * ( this%g31*this%sxx + this%g32*this%sxy + this%g33*this%sxz ) ! g31 
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + invtaurel * ( this%g31*this%sxy + this%g32*this%syy + this%g33*this%syz ) ! g32 
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + invtaurel * ( this%g31*this%sxz + this%g32*this%syz + this%g33*this%szz ) ! g33 

    end subroutine

    subroutine update_Ys(this,isub,dt,rho,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhsYs  ! RHS for mass fraction equation

        call this%getRHS_Ys(rho,u,v,w,rhsYs,x_bc,y_bc,z_bc)
        call hook_material_mass_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsYs)

        ! advance sub-step
        if(isub==1) this%QtmpYs = zero                   ! not really needed, since RK45_A(1) = 0
        this%QtmpYs  = dt*rhsYs + RK45_A(isub)*this%QtmpYs
        this%consrv(:,:,:,1) = this%consrv(:,:,:,1)  + RK45_B(isub)*this%QtmpYs

    end subroutine

    subroutine getRHS_Ys(this,rho,u,v,w,rhsYs,x_bc,y_bc,z_bc)
        use operators, only: divergence
        class(solid),                                         intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(out) :: rhsYs
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp1, tmp2, tmp3

        rhsYs = -rho*this%Ys
        tmp1 = rhsYs*u - this%Ji(:,:,:,1)   
        tmp2 = rhsYs*v - this%Ji(:,:,:,2)
        tmp3 = rhsYs*w - this%Ji(:,:,:,3)

        call divergence(this%decomp,this%der,tmp1,tmp2,tmp3,rhsYs,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric

    end subroutine

    subroutine update_eh(this,isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,x_bc,y_bc,z_bc)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,divu,viscwork
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhseh  ! RHS for eh equation

        if(this%PTeqb) then
            call GracefulExit("update_eh shouldn't be called with PTeqb. Exiting.",4809)
        endif

        call this%getRHS_eh(rho,u,v,w,divu,viscwork,rhseh,x_bc,y_bc,z_bc)
        call hook_material_energy_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhseh)

        ! advance sub-step
        if(isub==1) this%Qtmpeh = zero                   ! not really needed, since RK45_A(1) = 0
        this%Qtmpeh  = dt*rhseh + RK45_A(isub)*this%Qtmpeh
        this%consrv(:,:,:,2) = this%consrv(:,:,:,2) + RK45_B(isub)*this%Qtmpeh

    end subroutine

    subroutine getRHS_eh(this,rho,u,v,w,divu,viscwork,rhseh,x_bc,y_bc,z_bc)
        use operators, only: gradient, divergence
        class(solid),                                       intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: divu,viscwork
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhseh
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp1, tmp2, tmp3

        ! artificial conductivity and inter-species enthalpy diffusion fluxes
        tmp1 = -this%qi(:,:,:,1) ! * this%VF 
        tmp2 = -this%qi(:,:,:,2) ! * this%VF 
        tmp3 = -this%qi(:,:,:,3) ! * this%VF 

        ! add convective fluxes
        rhseh = -rho*this%Ys*this%eh
        tmp1 = tmp1 + rhseh*u
        tmp2 = tmp2 + rhseh*v
        tmp3 = tmp3 + rhseh*w

        ! Take divergence of fluxes
        call divergence(this%decomp,this%der,tmp1,tmp2,tmp3,rhseh,-x_bc,-y_bc,-z_bc)     ! energy has to be anti-symmetric

        ! Add pressure and viscous work terms
        rhseh = rhseh - this%VF * (this%p*divu + viscwork)  ! full viscous stress tensor here so equation is exact in the stiffened gas limit

    end subroutine

    subroutine update_VF(this,other,isub,dt,rho,u,v,w,x,y,z,tsim,divu,x_bc,y_bc,z_bc)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        class(solid), intent(in)    :: other
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,divu
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhsVF  ! RHS for mass fraction equation

        if(this%PTeqb) then
            call GracefulExit("update_VF shouldn't be called with PTeqb. Exiting.",4809)
        endif

        call this%getRHS_VF(other,rho,u,v,w,divu,rhsVF,x_bc,y_bc,z_bc)
        call hook_material_VF_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,u,v,w,this%Ys,this%VF,this%p,rhsVF)

        ! advance sub-step
        if(isub==1) this%QtmpVF = zero                   ! not really needed, since RK45_A(1) = 0
        this%QtmpVF  = dt*rhsVF + RK45_A(isub)*this%QtmpVF
        this%VF = this%VF  + RK45_B(isub)*this%QtmpVF

    end subroutine

    subroutine getRHS_VF(this,other,rho,u,v,w,divu,rhsVF,x_bc,y_bc,z_bc)
        use operators, only: gradient, divergence
        use constants, only: one
        class(solid),                                       intent(in)  :: this
        class(solid),                                       intent(in)  :: other
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w,divu
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhsVF
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp1, tmp2, tmp3, rhocsq1, rhocsq2

        ! any penalty term to keep VF between 0 and 1???
        call this%getSpeciesDensity(rho,tmp1)
        call this%hydro%get_sos2(tmp1,this%p,tmp2)
        call this%elastic%get_sos2(tmp1,tmp2)
        rhocsq1 = tmp1*tmp2

        call other%getSpeciesDensity(rho,tmp1)
        call other%hydro%get_sos2(tmp1,other%p,tmp3)
        call other%elastic%get_sos2(tmp1,tmp3)
        rhocsq2 = tmp1*tmp3

        rhsVF = other%VF*tmp2 + this%VF*tmp3

        call this%get_enthalpy(tmp2)
        call other%get_enthalpy(tmp3)

        call divergence(this%decomp,this%der,-this%Ji(:,:,:,1),-this%Ji(:,:,:,2),-this%Ji(:,:,:,3),tmp1,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric

        rhsVF = -(rhocsq1 - rhocsq2)*divu*this%VF*other%VF + (rhsVF + this%VF*(other%hydro%gam-one)*(tmp2-tmp3))*tmp1

        rhsVF = rhsVF / (other%VF*rhocsq1 + this%VF*rhocsq2)

        call gradient(this%decomp,this%der,-this%VF,tmp1,tmp2,tmp3,x_bc,y_bc,z_bc)
        rhsVF = rhsVF + u*tmp1 + v*tmp2 + w*tmp3

    end subroutine

    subroutine filter(this, iflag, x_bc, y_bc, z_bc)
        use operators, only: filter3D
        class(solid),  intent(inout) :: this
        integer,       intent(in)    :: iflag
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        ! filter g
        call filter3D(this%decomp, this%fil, this%g11, iflag, x_bc, y_bc, z_bc)
        call filter3D(this%decomp, this%fil, this%g12, iflag,-x_bc,-y_bc, z_bc)
        call filter3D(this%decomp, this%fil, this%g13, iflag,-x_bc, y_bc,-z_bc)
        call filter3D(this%decomp, this%fil, this%g21, iflag,-x_bc,-y_bc, z_bc)
        call filter3D(this%decomp, this%fil, this%g22, iflag, x_bc, y_bc, z_bc)
        call filter3D(this%decomp, this%fil, this%g23, iflag, x_bc,-y_bc,-z_bc)
        call filter3D(this%decomp, this%fil, this%g31, iflag,-x_bc, y_bc,-z_bc)
        call filter3D(this%decomp, this%fil, this%g32, iflag, x_bc,-y_bc,-z_bc)
        call filter3D(this%decomp, this%fil, this%g33, iflag, x_bc, y_bc, z_bc)

        ! filter Ys
        call filter3D(this%decomp, this%fil, this%consrv(:,:,:,1), iflag, x_bc, y_bc, z_bc)

        if(.NOT. this%PTeqb) then
            ! filter eh
            call filter3D(this%decomp, this%fil, this%consrv(:,:,:,2), iflag,x_bc, y_bc, z_bc)

            ! filter VF
            call filter3D(this%decomp, this%fil, this%VF, iflag, x_bc, y_bc,z_bc)
        endif

    end subroutine

    subroutine getSpeciesDensity(this,rho,rhom)
        class(solid), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhom

        ! Get detg in rhom
        rhom = this%g11*(this%g22*this%g33-this%g23*this%g32) &
             - this%g12*(this%g21*this%g33-this%g31*this%g23) &
             + this%g13*(this%g21*this%g32-this%g31*this%g22)

        if (this%use_gTg) then
            rhom = sqrt(rhom)
        end if

        ! Get rhom = rho*Ys/VF (Additional terms to give correct limiting behaviour when Ys and VF tend to 0)
        rhom = (rho*this%Ys + this%elastic%rho0*rhom*epssmall)/(this%VF + epssmall)   

    end subroutine

    ! computes p from ehydro
    subroutine get_p_from_ehydro(this, rho)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)              :: rhom

        call this%getSpeciesDensity(rho,rhom)
        call this%hydro%get_p( rhom, this%eh, this%p )
        ! call this%hydro%get_p( this%Ys*rho/(this%VF+epssmall), this%eh, this%p )

    end subroutine

    ! computes ehydro from p; and T from ehydro
    subroutine get_ehydroT_from_p(this, rho)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)              :: rhom

        call this%getSpeciesDensity(rho,rhom)
        call this%hydro%get_e_from_p( rhom, this%p, this%eh )
        call this%hydro%get_T(this%eh, this%T, rhom)
        ! call this%hydro%get_e_from_p( this%Ys*rho/(this%VF+epssmall), this%p, this%eh )
        ! call this%hydro%get_T(this%eh, this%T, this%Ys*rho/(this%VF+epssmall))

    end subroutine

    subroutine get_enthalpy(this,enthalpy)
        class(solid), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: enthalpy

        call this%hydro%get_enthalpy(this%T,enthalpy)
    end subroutine

    subroutine get_eelastic_devstress(this)
        class(solid), intent(inout) :: this

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6)  :: finger,fingersq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)    :: trG, trG2, detG

        call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG,this%use_gTg)
        call this%elastic%get_eelastic(trG,trG2,detG,this%eel)
        call this%elastic%get_devstress(finger, fingersq, trG, trG2, detG, this%devstress)

    end subroutine

    pure subroutine get_conserved(this,rho)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho

        this%consrv(:,:,:,1) = rho * this%Ys
        if(.NOT. this%PTeqb) this%consrv(:,:,:,2) = this%consrv(:,:,:,1) * this%eh

    end subroutine

    subroutine get_primitive(this,rho)
        use operators, only : gradient
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)                :: rhom

        this%Ys = this%consrv(:,:,:,1) / rho
        if(.NOT. this%PTeqb) then
            this%eh = this%consrv(:,:,:,2) / this%consrv(:,:,:,1)

            call this%getSpeciesDensity(rho,rhom)
            call this%hydro%get_T(this%eh, this%T, rhom)
            ! call this%hydro%get_T(this%eh, this%T, this%consrv(:,:,:,1)/(this%VF+epssmall))
        endif

        ! Get gradients of Ys and put in Ji for subsequent use
        call gradient(this%decomp,this%der,this%Ys,this%Ji(:,:,:,1),this%Ji(:,:,:,2),this%Ji(:,:,:,3))

    end subroutine

    subroutine checkNaN(this,imat)
        use exits, only : nancheck
        class(solid), intent(in) :: this
        integer, intent(in) :: imat
        character(len=clen) :: charout

        integer :: i,j,k,l

        if ( nancheck(this%g,i,j,k,l) ) then
            write(charout,'(A,5(I0,A))') "NaN encountered in material ", imat,&
                         &" at (",i,", ",j,", ",k,", ",l,") of g"
            call GracefulExit(charout,4809)
        end if
        if ( nancheck(this%VF,i,j,k) ) then
            write(charout,'(A,4(I0,A))') "NaN encountered in material ", imat,&
                         &" at (",i,", ",j,", ",k,") of VF"
            call GracefulExit(charout,4809)
        end if
        if ( nancheck(this%consrv,i,j,k,l) ) then
            write(charout,'(A,5(I0,A))') "NaN encountered in material ", imat,&
                         &" at (",i,", ",j,", ",k,", ",l,") of consrv"
            call GracefulExit(charout,4809)
        end if

    end subroutine

end module
