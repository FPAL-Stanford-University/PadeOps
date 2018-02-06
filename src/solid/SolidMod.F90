module SolidMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,one,two,epssmall,three,half
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

        logical :: sliding = .FALSE.
        logical :: plast = .FALSE.
        logical :: explPlast = .FALSE.
        logical :: PTeqb = .TRUE., pEqb = .FALSE., pRelax = .FALSE., updateEtot = .TRUE., includeSources = .FALSE.
        logical :: use_gTg = .FALSE.

        class(stiffgas ), allocatable :: hydro
        class(sep1solid), allocatable :: elastic

        type(decomp_info), pointer :: decomp
        type(derivatives), pointer :: der
        type(filters),     pointer :: fil
        type(filters),     pointer :: gfil

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
        
        real(rkind), dimension(:,:,:),   allocatable :: rhom
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
        !procedure :: get_eelastic_devstress_mixture
        procedure :: get_conserved
        procedure :: get_primitive
        procedure :: getSpeciesDensity
        procedure :: getSpeciesDensity_from_g
        procedure :: get_enthalpy
        procedure :: checkNaN
        procedure :: filter
        procedure :: sliding_deformation
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
    subroutine init(this,decomp,der,fil,gfil,PTeqb,pEqb,pRelax,use_gTg,updateEtot)
        class(solid), target, intent(inout) :: this
        type(decomp_info), target, intent(in) :: decomp
        type(derivatives), target, intent(in) :: der
        type(filters),     target, intent(in) :: fil, gfil
        logical, intent(in) :: PTeqb,pEqb,pRelax,updateEtot
        logical, intent(in) :: use_gTg

        this%decomp => decomp
        this%der  => der
        this%fil  => fil
        this%gfil => gfil
       
        this%PTeqb  = PTeqb
        this%pEqb   = pEqb
        this%pRelax = pRelax

        this%use_gTg = use_gTg
        this%updateEtot  = updateEtot

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
        
        ! Allocate material density
        if( allocated( this%rhom ) ) deallocate( this%rhom )
        allocate( this%rhom(this%nxp,this%nyp,this%nzp) )
        
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
        if(this%PTeqb .or. this%pEqb) then
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

        if(this%pEqb) then
            ! VF equation
            if( allocated( this%QtmpVF ) ) deallocate( this%QtmpVF )
            allocate( this%QtmpVF(this%nxp,this%nyp,this%nzp) )
        endif

        if(this%pRelax) then
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

        if( allocated( this%rhom) ) deallocate( this%rhom )
        if( allocated( this%eel ) ) deallocate( this%eel )
        if( allocated( this%eh )  ) deallocate( this%eh )
        if( allocated( this%VF )  ) deallocate( this%VF )
        if( allocated( this%Ys )  ) deallocate( this%Ys )

        ! Now deallocate the EOS objects
        if ( allocated(this%hydro)   ) deallocate(this%hydro)
        if ( allocated(this%elastic) ) deallocate(this%elastic)

        nullify( this%fil    )
        nullify( this%gfil   )
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

    subroutine update_g(this,isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,rho0mix,mumix,yieldmix,solidVF)
        use constants,  only: eps
        use RKCoeffs,   only: RK45_A,RK45_B
        use reductions, only: P_MAXVAL
        use operators,  only: gradient, filter3D
        use exits,      only : nancheck
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in), optional  :: rho0mix, mumix, yieldmix, solidVF

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg   ! RHS for g tensor equation
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: normal ! Interface normal for sliding treatment
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: mask, auxVF
        real(rkind), parameter :: auxVF_exponent = 0.5_rkind
        real(rkind) :: max_modDevSigma
        integer :: i,j,k,l
        character(len=clen) :: charout
        integer, parameter :: mask_exponent = 40

        ! ! Set normal to be radially outwards for a cylindrical interface
        ! theta = atan2(y,x)
        ! normal(:,:,:,1) = cos(theta); normal(:,:,:,2) = sin(theta); normal(:,:,:,3) = zero

        ! Set normal to be the maximal gradient direction of the volume fraction
        ! Use auxiliary function as in Shukla, Pantano and Freund (JCP 2010)
        auxVF = this%VF
        where (auxVF < zero)
            auxVF = zero
        end where
        where (auxVF > one)
            auxVF = one
        end where
        call filter3D(this%decomp, this%gfil, auxVF, 1)
        auxVF = (auxVF**auxVF_exponent) / ( auxVF**auxVF_exponent + (one - auxVF)**auxVF_exponent )
        this%kap = auxVF
        call gradient(this%decomp, this%der, auxVF, normal(:,:,:,1), normal(:,:,:,2), normal(:,:,:,3), x_bc, y_bc, z_bc)
        ! Normalize to magnitude 1
        mask = sqrt( normal(:,:,:,1)*normal(:,:,:,1) + normal(:,:,:,2)*normal(:,:,:,2) + normal(:,:,:,3)*normal(:,:,:,3) )
        do l = 1,3
            do k = 1,this%nzp
                do j = 1,this%nyp
                    do i = 1,this%nxp
                        if (mask(i,j,k) > eps) then
                            normal(i,j,k,l) = normal(i,j,k,l) / mask(i,j,k)
                        else
                            normal(i,j,k,:) = [one, zero, zero] ! This is arbitrary but should not affect anything
                        end if
                    end do
                end do
            end do
        end do
        if ( nancheck(normal,i,j,k,l) ) then
            write(charout,'(A,4(I0,A))') "NaN encountered in interface normal at (",i,", ",j,", ",k,", ",l,")"
            call GracefulExit(charout,4809)
        end if

        ! Get the mask function to apply sliding treatment
        if (this%sliding) then
            mask = this%VF*(one - this%VF)
            mask = mask/P_MAXVAL(mask)
            mask = one - (one - mask)**mask_exponent
        else
            mask = zero
        end if

        if(present(rho0mix)) then
            call this%getRHS_g(rho,u,v,w,dt,src,normal,mask,rhsg,x_bc,y_bc,z_bc,rho0mix,solidVF)
        else
            call this%getRHS_g(rho,u,v,w,dt,src,normal,mask,rhsg,x_bc,y_bc,z_bc)
        endif
        call hook_material_g_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsg)

        ! advance sub-step
        if(isub==1) this%Qtmpg = zero                   ! not really needed, since RK45_A(1) = 0
        this%Qtmpg  = dt*rhsg + RK45_A(isub)*this%Qtmpg
        !print *, 'Before : ', this%g11(89,1,1), rhsg(89,1,1,1)
        this%g = this%g  + RK45_B(isub)*this%Qtmpg
        !print *, 'After 1: ', this%g11(89,1,1)

        ! Now project g tensor to SPD space
        call this%elastic%make_tensor_SPD(this%g)
        !print *, 'After 2: ', this%g11(89,1,1)

        ! Sliding treatment using plasticity (Using zero yield everywhere for now)
        if (this%sliding) call this%sliding_deformation(normal, mask)

        if(this%plast) then
            if (.NOT. this%explPlast) then
                ! Effect plastic deformations
                if(present(rho0mix) .and. present(mumix)) then
                    call this%get_eelastic_devstress(rho0mix,mumix)
                else
                    call this%get_eelastic_devstress()
                endif
                this%modDevSigma = this%sxx*this%sxx + this%syy*this%syy + this%szz*this%szz + &
                        two*( this%sxy*this%sxy + this%sxz*this%sxz + this%syz*this%syz )
                max_modDevSigma = sqrt(3.0d0/2.0d0*maxval(this%modDevSigma))
                !if(max_modDevSigma > this%elastic%yield) then
                !  write(*,'(a,2(e19.12,1x))') 'Entering plasticity. Stresses = ', max_modDevSigma, this%elastic%yield
                !endif
                if(this%PTeqb) this%kap = sqrt(two/three*this%modDevSigma) !-- is this temporary storage ?? -- NSG

                if(present(yieldmix) .and. present(mumix)) then
                    call this%elastic%plastic_deformation(this%g, this%use_gTg, mumix, yieldmix)
                else
                    call this%elastic%plastic_deformation(this%g, this%use_gTg)
                endif

                if(present(rho0mix) .and. present(mumix)) then
                    call this%get_eelastic_devstress(rho0mix,mumix)
                else
                    call this%get_eelastic_devstress()
                endif
                this%modDevSigma = this%sxx*this%sxx + this%syy*this%syy + this%szz*this%szz + &
                        two*( this%sxy*this%sxy + this%sxz*this%sxz + this%syz*this%syz )
                max_modDevSigma = sqrt(3.0d0/2.0d0*maxval(this%modDevSigma))
                !if(max_modDevSigma > this%elastic%yield) then
                !  write(*,'(a,2(e19.12,1x))') 'Exiting plasticity. Stresses = ', max_modDevSigma, this%elastic%yield
                !endif
                !write(*,*) 'modDevSig: ', max_modDevSigma, this%elastic%yield
                
            end if
        end if
        !print *, 'After 3: ', this%g11(89,1,1)

    end subroutine

    subroutine getRHS_g(this,rho,u,v,w,dt,src,normal,mask,rhsg,x_bc,y_bc,z_bc,rho0mix,solidVF)
        use decomp_2d, only: nrank
        use constants, only: eps, pi
        use operators, only: gradient, curl
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src,mask
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3), intent(in)  :: normal
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix, solidVF

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dutdx,dutdy,dutdz,dvtdx,dvtdy,dvtdz,dwtdx,dwtdy,dwtdz

        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg, ut, vt, wt, rad, theta
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg
        real(rkind), parameter :: etafac = one/6._rkind

        real(rkind) :: dx, dy, x, y
        ! real(rkind), parameter :: theta = 45._rkind * pi / 180._rkind

        real(rkind) :: pmax
        integer :: imax, jmax, kmax, rmax
        integer :: i,j,k

        ! Symmetry and anti-symmetry properties of g are assumed as below
        ! In x g_{ij}: [S A A; A S S; A S S]
        ! In y g_{ij}: [S A S; A S A; S A S]
        ! In z g_{ij}: [S S A; S S A; A A S]

        rhsg = zero

        ! if (this%sliding) then
        !     mask = this%VF*(one - this%VF)
        !     mask = mask/P_MAXVAL(mask)
        !     mask = one - (one - mask)**mask_exponent
        ! else
        !     mask = zero
        ! end if
        ! mask = zero

        ! Set normal to be aligned with the x direction
        ! normal(:,:,:,1) = one; normal(:,:,:,2) = zero; normal(:,:,:,3) = zero

        ! Set normal to be normal to the specified oblique interface
        ! normal(:,:,:,1) = cos(theta); normal(:,:,:,2) =-sin(theta); normal(:,:,:,3) = zero

        ! Set normal to be radially outwards for a cylindrical interface
        ! dx = one/real(this%decomp%xsz(1)-1,rkind)
        ! dy = dx
        ! do k=1,this%decomp%ysz(3)
        !     do j=1,this%decomp%ysz(2)
        !         do i=1,this%decomp%ysz(1)
        !             x = - half + real( this%decomp%yst(1) - 1 + i - 1, rkind ) * dx
        !             y = - half + real( this%decomp%yst(2) - 1 + j - 1, rkind ) * dy
        !             rad(i,j,k) = sqrt(x**2 + y**2)
        !             theta(i,j,k) = atan2(y,x)
        !         end do
        !     end do
        ! end do
        ! normal(:,:,:,1) = cos(theta); normal(:,:,:,2) = sin(theta); normal(:,:,:,3) = zero

        ! ! Set normal to be the maximal volume fraction gradient direction
        ! call gradient(this%decomp, this%der, this%VF, normal(:,:,:,1), normal(:,:,:,2), normal(:,:,:,3), x_bc, y_bc, z_bc)
        ! ! Normalize to magnitude 1
        ! tmp = sqrt( normal(:,:,:,1)*normal(:,:,:,1) + normal(:,:,:,2)*normal(:,:,:,2) + normal(:,:,:,3)*normal(:,:,:,3) )
        ! do i = 1,3
        !     normal(:,:,:,i) = normal(:,:,:,i) / tmp
        ! end do

        ! print *, "mask = ", mask(39,83,1)
        ! print *, "normal = ", normal(39,83,1,:)

        ! Magnitude of normal velocity
        wt = u*normal(:,:,:,1) + v*normal(:,:,:,2) + w*normal(:,:,:,3)

        ! print *, "normal velocity = ", wt(39,83,1)

        ! Get tangential velocity
        ut = u - normal(:,:,:,1)*wt
        vt = v - normal(:,:,:,2)*wt
        wt = w - normal(:,:,:,3)*wt

        ! print *, "tangential velocity = ", sqrt(ut(39,83,1)**2 + vt(39,83,1)**2 + wt(39,83,1)**2)

        dutdx => duidxj(:,:,:,1); dutdy => duidxj(:,:,:,2); dutdz => duidxj(:,:,:,3);
        dvtdx => duidxj(:,:,:,4); dvtdy => duidxj(:,:,:,5); dvtdz => duidxj(:,:,:,6);
        dwtdx => duidxj(:,:,:,7); dwtdy => duidxj(:,:,:,8); dwtdz => duidxj(:,:,:,9);
        
        call gradient(this%decomp, this%der, ut, dutdx, dutdy, dutdz, -x_bc,  y_bc,  z_bc)
        call gradient(this%decomp, this%der, vt, dvtdx, dvtdy, dvtdz,  x_bc, -y_bc,  z_bc)
        call gradient(this%decomp, this%der, wt, dwtdx, dwtdy, dwtdz,  x_bc,  y_bc, -z_bc)

        ! Get the normal derivative of the tangential velocity and recast them in 3D coordinates as a tensor
        tmp = dutdx*normal(:,:,:,1) + dutdy*normal(:,:,:,2) + dutdz*normal(:,:,:,3)
        dutdx = tmp*normal(:,:,:,1); dutdy = tmp*normal(:,:,:,2); dutdz = tmp*normal(:,:,:,3); 

        tmp = dvtdx*normal(:,:,:,1) + dvtdy*normal(:,:,:,2) + dvtdz*normal(:,:,:,3)
        dvtdx = tmp*normal(:,:,:,1); dvtdy = tmp*normal(:,:,:,2); dvtdz = tmp*normal(:,:,:,3); 

        tmp = dwtdx*normal(:,:,:,1) + dwtdy*normal(:,:,:,2) + dwtdz*normal(:,:,:,3)
        dwtdx = tmp*normal(:,:,:,1); dwtdy = tmp*normal(:,:,:,2); dwtdz = tmp*normal(:,:,:,3); 

        detg = this%g11*(this%g22*this%g33-this%g23*this%g32) &
             - this%g12*(this%g21*this%g33-this%g31*this%g23) &
             + this%g13*(this%g21*this%g32-this%g31*this%g22)

        ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
        tmp = (rho*this%Ys + this%elastic%rho0*detg*epssmall)/(this%VF + epssmall)   
        ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF

        if(present(rho0mix)) then
            penalty = etafac*(rho/rho0mix/detg - one)/dt
        else
            penalty = etafac*( tmp/detg/this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
        endif
        if(this%pRelax) penalty = this%VF*etafac*( tmp/detg/this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density -- change2

        if (this%elastic%mu < eps) penalty = zero

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty = penalty - src/this%VF
        endif

        tmp = -u*this%g11-v*this%g12-w*this%g13
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,1),rhsg(:,:,:,2),rhsg(:,:,:,3),-x_bc, y_bc, z_bc)
        !do i = 1, size(u,1)
        !  print '(4(e19.12,1x))', u(i,1,1), this%g11(i,1,1), tmp(i,1,1), rhsg(i,1,1,1)
        !enddo
        !print *, 'rhsg 1 : ', rhsg(89,1,1,1)
        
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
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g11 + mask*(this%g11*dutdx + this%g12*dvtdx + this%g13*dwtdx)
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g12 + mask*(this%g11*dutdy + this%g12*dvtdy + this%g13*dwtdy)
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g13 + mask*(this%g11*dutdz + this%g12*dvtdz + this%g13*dwtdz)
        !print *, 'rhsg 2 : ', rhsg(89,1,1,1)
        !print *, '------ : ', v(89,1,1), curlg(89,1,1,3)
        !print *, '------ : ', w(89,1,1), curlg(89,1,1,2)
        !print *, '------ : ', penalty(89,1,1), this%g11(89,1,1)
        !print *, '------ : ', mask(89,1,1)
 
        tmp = -u*this%g21-v*this%g22-w*this%g23
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,4),rhsg(:,:,:,5),rhsg(:,:,:,6), x_bc,-y_bc, z_bc)   
        
        call curl(this%decomp, this%der, this%g21, this%g22, this%g23, curlg, x_bc, -y_bc, z_bc)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g21 + mask*(this%g21*dutdx + this%g22*dvtdx + this%g23*dwtdx)
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g22 + mask*(this%g21*dutdy + this%g22*dvtdy + this%g23*dwtdy)
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g23 + mask*(this%g21*dutdz + this%g22*dvtdz + this%g23*dwtdz)
 
        tmp = -u*this%g31-v*this%g32-w*this%g33
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,7),rhsg(:,:,:,8),rhsg(:,:,:,9), x_bc, y_bc,-z_bc)

        call curl(this%decomp, this%der, this%g31, this%g32, this%g33, curlg, x_bc, y_bc, -z_bc)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) + v*curlg(:,:,:,3) - w*curlg(:,:,:,2) + penalty*this%g31 + mask*(this%g31*dutdx + this%g32*dvtdx + this%g33*dwtdx)
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + w*curlg(:,:,:,1) - u*curlg(:,:,:,3) + penalty*this%g32 + mask*(this%g31*dutdy + this%g32*dvtdy + this%g33*dwtdy)
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + u*curlg(:,:,:,2) - v*curlg(:,:,:,1) + penalty*this%g33 + mask*(this%g31*dutdz + this%g32*dvtdz + this%g33*dwtdz)

        if (this%plast) then
            if(this%explPlast) then
                call this%getPlasticSources(detg,rhsg)
            end if
        end if

        if(present(solidVF)) then
           do i = 1, 9
               rhsg(:,:,:,i) = rhsg(:,:,:,i) * solidVF
           enddo
        endif
        !print *, 'rhsg 3 : ', rhsg(89,1,1,1)

    end subroutine

    subroutine update_gTg(this,isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,rho0mix,mumix,yieldmix,solidVF)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in), optional  :: rho0mix, mumix, yieldmix, solidVF

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg  ! RHS for gTg tensor equation

        if(present(rho0mix) .and. present(solidVF)) then
            call this%getRHS_gTg(rho,u,v,w,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix,solidVF)
        else
            call this%getRHS_gTg(rho,u,v,w,dt,src,rhsg,x_bc,y_bc,z_bc)
        endif
        call hook_material_g_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsg)

        ! advance sub-step
        if(isub==1) this%Qtmpg = zero                   ! not really needed, since RK45_A(1) = 0
        this%Qtmpg  = dt*rhsg + RK45_A(isub)*this%Qtmpg
        this%g = this%g  + RK45_B(isub)*this%Qtmpg

        if(this%plast) then
            if (.NOT. this%explPlast) then
                if(present(yieldmix) .and. present(mumix)) then
                    call this%elastic%plastic_deformation(this%g, this%use_gTg, mumix, yieldmix)
                else
                    call this%elastic%plastic_deformation(this%g, this%use_gTg)
                endif
            end if
        end if

    end subroutine

    subroutine getRHS_gTg(this,rho,u,v,w,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix,solidVF)
        use operators, only: gradient, curl
        class(solid),                                         intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix, solidVF

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradG
        real(rkind), parameter :: etafac = one/6._rkind
        integer :: i

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
        if(present(rho0mix)) then
            penalty = etafac*( rho/rho0mix/detg-one)/dt ! Penalty term to keep g consistent with species density
        else
            penalty = etafac*( tmp/detg/this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
        endif

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

        if(present(solidVF)) then
           do i = 1, 9
               rhsg(:,:,:,i) = rhsg(:,:,:,i) * solidVF
           enddo
        endif

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
!print *, 'rhs Ys:', rhsYs(89,1,1), rho(89,1,1), this%Ys(89,1,1)
!print *, 'Ji:    ', this%Ji(89,1,1,1), this%Ji(89,1,1,2), this%Ji(89,1,1,3)
!print *, 'cns Ys:', this%consrv(89,1,1,1)
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

    subroutine update_eh(this,isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,src,taustar,x_bc,y_bc,z_bc)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,divu,viscwork,src
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(in) :: taustar
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhseh  ! RHS for eh equation

        if(.not. this%pRelax) then
            call GracefulExit("update_eh shouldn't be called without pRelax. Exiting.",4809)
        endif

        call this%getRHS_eh(rho,u,v,w,divu,viscwork,src,taustar,rhseh,x_bc,y_bc,z_bc)
        call hook_material_energy_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhseh)

        ! advance sub-step
        if(isub==1) this%Qtmpeh = zero                   ! not really needed, since RK45_A(1) = 0
        this%Qtmpeh  = dt*rhseh + RK45_A(isub)*this%Qtmpeh
        this%consrv(:,:,:,2) = this%consrv(:,:,:,2) + RK45_B(isub)*this%Qtmpeh

    end subroutine

    subroutine getRHS_eh(this,rho,u,v,w,divu,viscwork,src,taustar,rhseh,x_bc,y_bc,z_bc)
        use operators, only: gradient, divergence
        class(solid),                                       intent(in)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: divu,viscwork,src
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(in) :: taustar
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhseh
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp1, tmp2, tmp3

        ! artificial conductivity and inter-species enthalpy diffusion fluxes
        tmp1 = -this%qi(:,:,:,1) ! * this%VF 
        tmp2 = -this%qi(:,:,:,2) ! * this%VF 
        tmp3 = -this%qi(:,:,:,3) ! * this%VF 

        ! add convective fluxes
        if(this%updateEtot) then
            rhseh = -rho*this%Ys*(this%eh + this%eel + half*(u*u + v*v + w*w))
        else
            rhseh = -rho*this%Ys*this%eh
        endif
        tmp1 = tmp1 + rhseh*u
        tmp2 = tmp2 + rhseh*v
        tmp3 = tmp3 + rhseh*w

        if(this%updateEtot) then
            tmp1 = tmp1 + this%VF*((taustar(:,:,:,1) + this%sxx - this%p)*u + (taustar(:,:,:,2) + this%sxy         )*v + (taustar(:,:,:,3) + this%sxz         )*w)
            tmp2 = tmp2 + this%VF*((taustar(:,:,:,2) + this%sxy         )*u + (taustar(:,:,:,4) + this%syy - this%p)*v + (taustar(:,:,:,5) + this%syz         )*w)
            tmp3 = tmp3 + this%VF*((taustar(:,:,:,3) + this%sxz         )*u + (taustar(:,:,:,5) + this%syz         )*v + (taustar(:,:,:,6) + this%szz - this%p)*w)
        endif

        ! Take divergence of fluxes
        call divergence(this%decomp,this%der,tmp1,tmp2,tmp3,rhseh,-x_bc,-y_bc,-z_bc)     ! energy has to be anti-symmetric

        if(this%updateEtot) then
            ! Add source
            if(this%includeSources) rhseh = rhseh + src
        else
            ! Add pressure and viscous work terms
            rhseh = rhseh - this%VF * (this%p*divu + viscwork)  ! full viscous stress tensor here so equation is exact in the stiffened gas limit
        endif

    end subroutine

    subroutine update_VF(this,other,isub,dt,rho,u,v,w,x,y,z,tsim,divu,src,x_bc,y_bc,z_bc)
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        class(solid), intent(in)    :: other
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,divu,src
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhsVF  ! RHS for mass fraction equation

        if(this%PTeqb) then
            call GracefulExit("update_VF shouldn't be called with PTeqb. Exiting.",4809)
        endif

        call this%getRHS_VF(other,rho,u,v,w,divu,src,rhsVF,x_bc,y_bc,z_bc)
        call hook_material_VF_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,u,v,w,this%Ys,this%VF,this%p,rhsVF)

        ! advance sub-step
        if(isub==1) this%QtmpVF = zero                   ! not really needed, since RK45_A(1) = 0
        this%QtmpVF  = dt*rhsVF + RK45_A(isub)*this%QtmpVF
        this%VF = this%VF  + RK45_B(isub)*this%QtmpVF

    end subroutine

    subroutine getRHS_VF(this,other,rho,u,v,w,divu,src,rhsVF,x_bc,y_bc,z_bc)
        use operators, only: gradient, divergence
        use constants, only: one
        class(solid),                                       intent(in)  :: this
        class(solid),                                       intent(in)  :: other
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w,divu,src
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhsVF
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp1, tmp2, tmp3, rhocsq1, rhocsq2

        ! Add C/rhom to Fsource
        !--call this%getSpeciesDensity(rho,tmp1)  !-- use this%rhom
        call divergence(this%decomp,this%der,this%Ji(:,:,:,1),this%Ji(:,:,:,2),this%Ji(:,:,:,3),rhsVF,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric
        if(this%pEqb) then
            rhsVF = src - rhsVF/this%rhom
        else
            rhsVF = -rhsVF/this%rhom
        endif
        if(.not. this%includeSources) rhsVF = zero

        call gradient(this%decomp,this%der,-this%VF,tmp1,tmp2,tmp3,x_bc,y_bc,z_bc)

        rhsVF = rhsVF + u*tmp1 + v*tmp2 + w*tmp3

    end subroutine

    subroutine filter(this, iflag, x_bc, y_bc, z_bc, fil)
        use operators, only: filter3D
        class(solid),  intent(inout) :: this
        integer,       intent(in)    :: iflag
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        type(filters), target, optional, intent(in) :: fil

        type(filters), pointer :: fil_

        if (present(fil)) then
            fil_ => fil
        else
            fil_ => this%fil
        end if

        ! filter g
        call filter3D(this%decomp, fil_, this%g11, iflag, x_bc, y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%g12, iflag,-x_bc,-y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%g13, iflag,-x_bc, y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%g21, iflag,-x_bc,-y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%g22, iflag, x_bc, y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%g23, iflag, x_bc,-y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%g31, iflag,-x_bc, y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%g32, iflag, x_bc,-y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%g33, iflag, x_bc, y_bc, z_bc)

        ! filter Ys
        call filter3D(this%decomp, fil_, this%consrv(:,:,:,1), iflag, x_bc, y_bc, z_bc)

        if(this%pEqb) then
            ! filter VF
            call filter3D(this%decomp, fil_, this%VF, iflag, x_bc, y_bc,z_bc)
        endif

        if(this%pRelax) then
            ! filter eh
            call filter3D(this%decomp, fil_, this%consrv(:,:,:,2), iflag,x_bc, y_bc, z_bc)

            ! filter VF
            call filter3D(this%decomp, fil_, this%VF, iflag, x_bc, y_bc,z_bc)
        endif

    end subroutine

    subroutine getSpeciesDensity(this,rho,rhom)
        class(solid), intent(inout) :: this
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

        this%rhom = rhom

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

    !subroutine get_eelastic_devstress_mixture(this,rho0mix, mumix)
    !    class(solid), intent(inout) :: this
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho0mix, mumix

    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp,6)  :: finger,fingersq
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp)    :: trG, trG2, detG

    !    call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG,this%use_gTg)
    !    call this%elastic%get_eelastic(trG,trG2,detG,this%eel,rho0mix,mumix)
    !    call this%elastic%get_devstress(finger, fingersq, trG, trG2, detG, this%devstress,rho0mix,mumix)

    !end subroutine

    subroutine get_eelastic_devstress(this,rho0mix,mumix)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional  :: rho0mix, mumix

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6)  :: finger,fingersq
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)    :: trG, trG2, detG

        call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG,this%use_gTg)
        if(present(rho0mix) .and. present(mumix)) then
          call this%elastic%get_eelastic(trG,trG2,detG,this%eel,rho0mix,mumix)
          call this%elastic%get_devstress(finger, fingersq, trG, trG2, detG, this%devstress, rho0mix, mumix)
        else
          call this%elastic%get_eelastic(trG,trG2,detG,this%eel)
          call this%elastic%get_devstress(finger, fingersq, trG, trG2, detG, this%devstress)
        endif

    end subroutine

    pure subroutine get_conserved(this,rho,u,v,w)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w

        this%consrv(:,:,:,1) = rho * this%Ys
        if(this%pRelax) then
            if(this%updateEtot) then
                this%consrv(:,:,:,2) = this%consrv(:,:,:,1) * (this%eh + this%eel + half*(u*u+v*v+w*w))
            else
                this%consrv(:,:,:,2) = this%consrv(:,:,:,1) * this%eh
            endif
        endif

    end subroutine

    subroutine get_primitive(this,rho,u,v,w)
        use operators, only : gradient
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)                :: rhom

        this%Ys = this%consrv(:,:,:,1) / rho
        if(this%pRelax) then
            if(this%updateEtot) then
                call this%get_eelastic_devstress()
                this%eh = this%consrv(:,:,:,2) / this%consrv(:,:,:,1) - this%eel - half*(u*u+v*v+w*w)
            else
                this%eh = this%consrv(:,:,:,2) / this%consrv(:,:,:,1)
            endif

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

    subroutine sliding_deformation(this, normal, mask)
        use kind_parameters, only: clen
        use constants,       only: zero, eps, third, twothird, pi
        use decomp_2d,       only: nrank
        use operators,       only: filter3D
        use exits,           only: GracefulExit, nancheck, message
        class(solid),                                         intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3), intent(in)    :: normal ! Interface normal for sliding treatment
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)    :: mask   ! Mask to avoid sliding away from interface
        real(rkind) :: yield = zero

        real(rkind), dimension(3,3) :: g, u, vt, gradf, gradf_new, sigma, sigma_tilde, v_tilde
        real(rkind), dimension(3)   :: sval, beta, Sa, f, f1, f2, dbeta, beta_new, dbeta_new
        real(rkind) :: sqrt_om, betasum, Sabymu_sq, ycrit, C0, t
        real(rkind) :: tol = real(1.D-12,rkind), residual, residual_new
        integer :: i,j,k
        integer, dimension(1) :: min_norm
        integer :: iters
        integer, parameter :: niters = 500
        integer :: lwork, info
        integer, dimension(3) :: ipiv
        integer :: nxp, nyp, nzp
        character(len=clen) :: charout

        real(rkind), dimension(:), allocatable :: svdwork     ! Work array for SVD stuff

        real(rkind) :: dx, dy, x, y, rad, theta = 45._rkind * pi / 180._rkind

        ! print *, "In sliding_deformation"

        ! where ((this%VF > 0.01) .and. (this%VF < 0.99))
        !     mask = one
        ! elsewhere
        !     mask = zero
        ! end where
        ! call filter3D(this%decomp, this%gfil, mask, 1)
        ! mask = this%VF*(one - this%VF)
        ! mask = mask/maxval(mask)
        ! mask = one - (one - mask)**mask_exponent

        allocate(svdwork(1))
        if (this%use_gTg) then
            ! Get optimal lwork
            lwork = -1
            call dsyev('V', 'U', 3, G, 3, sval, svdwork, lwork, info)
            lwork = svdwork(1)
        else
            ! Get optimal lwork
            lwork = -1
            call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, svdwork, lwork, info)
            lwork = svdwork(1)
        end if

        if (lwork .GT. size(svdwork)) then
            deallocate(svdwork); allocate(svdwork(lwork))
        end if

        if ( nancheck(this%g) ) then
            call message("NaN found in g during plastic relaxation.")
        end if

        dx = one/real(this%decomp%xsz(1)-1,rkind)
        dy = dx

        do k = 1,this%nzp
            do j = 1,this%nyp
                do i = 1,this%nxp

                    if (abs(mask(i,j,k)) > real(1.0D-2,rkind)) then
                        ! yield = zero
                        yield = (one-mask(i,j,k))*this%elastic%yield

                        g(1,1) = this%g(i,j,k,1); g(1,2) = this%g(i,j,k,2); g(1,3) = this%g(i,j,k,3)
                        g(2,1) = this%g(i,j,k,4); g(2,2) = this%g(i,j,k,5); g(2,3) = this%g(i,j,k,6)
                        g(3,1) = this%g(i,j,k,7); g(3,2) = this%g(i,j,k,8); g(3,3) = this%g(i,j,k,9)

                        !!!! HACK !!!
                        ! Hard code v_tilde for now as the cartesian basis
                        ! Should eventually set v_tilde(:,1) to the interface normal and the rest as orthogonal directions

                        ! 1D
                        ! v_tilde(1,1) = one;  v_tilde(1,2) = zero; v_tilde(1,3) = zero
                        ! v_tilde(2,1) = zero; v_tilde(2,2) = one;  v_tilde(2,3) = zero
                        ! v_tilde(3,1) = zero; v_tilde(3,2) = zero; v_tilde(3,3) = one

                        ! 2D oblique
                        ! v_tilde(1,1) = cos(theta); v_tilde(1,2) = sin(theta); v_tilde(1,3) = zero
                        ! v_tilde(2,1) =-sin(theta); v_tilde(2,2) = cos(theta); v_tilde(2,3) = zero
                        ! v_tilde(3,1) = zero;       v_tilde(3,2) = zero;       v_tilde(3,3) = one

                        ! 2D circular
                        ! x = - half + real( this%decomp%yst(1) - 1 + i - 1, rkind ) * dx
                        ! y = - half + real( this%decomp%yst(2) - 1 + j - 1, rkind ) * dy
                        ! rad = sqrt(x**2 + y**2)
                        ! theta = atan2(y,x)
                        ! this%kap(i,j,k) = theta
                        ! v_tilde(1,1) = cos(theta); v_tilde(1,2) =-sin(theta); v_tilde(1,3) = zero
                        ! v_tilde(2,1) = sin(theta); v_tilde(2,2) = cos(theta); v_tilde(2,3) = zero
                        ! v_tilde(3,1) = zero;       v_tilde(3,2) = zero;       v_tilde(3,3) = one

                        ! General interface normal
                        v_tilde(:,1) = normal(i,j,k,:) ! Assume normalized

                        ! Get a vector not parallel to v_tilde(:,1)
                        min_norm = minloc(abs(normal(i,j,k,:)))
                        v_tilde(:,3) = zero; v_tilde(min_norm(1),3) = one 
                        ! if ((i == 44) .and. (j == 2)) then
                        !     print *, "normal = ", normal(i,j,k,:)
                        !     print *, "2nd vector = ", v_tilde(:,3)
                        ! end if

                        ! Get an othogonal vector to v_tilde(:,1) through a crossproduct with a non-parallel vector
                        v_tilde(1,2) = v_tilde(2,1)*v_tilde(3,3) - v_tilde(3,1)*v_tilde(2,3)
                        v_tilde(2,2) = v_tilde(3,1)*v_tilde(1,3) - v_tilde(1,1)*v_tilde(3,3)
                        v_tilde(3,2) = v_tilde(1,1)*v_tilde(2,3) - v_tilde(2,1)*v_tilde(1,3)
                        v_tilde(:,2) = v_tilde(:,2)/sqrt(v_tilde(1,2)**2 + v_tilde(2,2)**2 + v_tilde(3,2)**2) ! Normalize

                        ! Get third orthogonal vector through a crossproduct of the first two
                        v_tilde(1,3) = v_tilde(2,1)*v_tilde(3,2) - v_tilde(3,1)*v_tilde(2,2)
                        v_tilde(2,3) = v_tilde(3,1)*v_tilde(1,2) - v_tilde(1,1)*v_tilde(3,2)
                        v_tilde(3,3) = v_tilde(1,1)*v_tilde(2,2) - v_tilde(2,1)*v_tilde(1,2)
                        v_tilde(:,3) = v_tilde(:,3)/sqrt(v_tilde(1,3)**2 + v_tilde(2,3)**2 + v_tilde(3,3)**2) ! Normalize

                        ! if ((i == 44) .and. (j == 2)) then
                        !     print *, "v_tilde:"
                        !     print *, "   ",  v_tilde(1,1), v_tilde(1,2), v_tilde(1,3)
                        !     print *, "   ",  v_tilde(2,1), v_tilde(2,2), v_tilde(2,3)
                        !     print *, "   ",  v_tilde(3,1), v_tilde(3,2), v_tilde(3,3)
                        ! end if

                        if (this%use_gTg) then
                            ! Get eigenvalues and eigenvectors of G
                            call dsyev('V', 'U', 3, G, 3, sval, svdwork, lwork, info)
                            if(info .ne. 0) print '(A,I6,A)', 'proc ', nrank, ': Problem with DSYEV. Please check.'
                            sval = sqrt(sval)  ! Get singular values of g
                        else
                            ! Get SVD of g
                            call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, svdwork, lwork, info)
                            if(info .ne. 0) then
                                write(charout, '(A,I0,A)') 'proc ', nrank, ': Problem with SVD. Please check.'
                                call GracefulExit(charout,3475)
                            end if
                        end if

                        sqrt_om = sval(1)*sval(2)*sval(3)
                        beta = sval**two / sqrt_om**(two/three)

                        betasum = sum( beta*(beta-one) ) / three
                        Sa = -this%elastic%mu*sqrt_om * ( beta*(beta-one) - betasum )

                        ! Get the actual stress
                        sigma = zero
                        sigma(1,1) = Sa(1); sigma(2,2) = Sa(2); sigma(3,3) = Sa(3)
                        ! sigma = matmul(sigma, transpose(vt))
                        ! sigma = matmul(vt, sigma)
                        sigma = matmul(sigma, vt)
                        sigma = matmul(transpose(vt), sigma)

                        ! if ((i == 44) .and. (j == 2)) then
                        !     print *, "V:"
                        !     print *, "   ",  vt(1,1), vt(2,1), vt(3,1)
                        !     print *, "   ",  vt(1,2), vt(2,2), vt(3,2)
                        !     print *, "   ",  vt(1,3), vt(2,3), vt(3,3)
                        ! end if

                        ! if ((i == 44) .and. (j == 2)) then
                        !     print *, "sigma:"
                        !     print *, "   ",  sigma(1,1), sigma(1,2), sigma(1,3)
                        !     print *, "   ",  sigma(2,1), sigma(2,2), sigma(2,3)
                        !     print *, "   ",  sigma(3,1), sigma(3,2), sigma(3,3)
                        ! end if

                        ! Get new stress
                        sigma_tilde = zero
                        sigma_tilde(1,1) = sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,1) ) )
                        sigma_tilde(2,2) = sum( v_tilde(:,2) * matmul( sigma, v_tilde(:,2) ) )
                        sigma_tilde(2,3) = sum( v_tilde(:,3) * matmul( sigma, v_tilde(:,2) ) ) ! Need this off-diagonal term
                                                                                                 ! since v_tilde(:,2:3) are not
                                                                                                 ! the singular vectors but just
                                                                                                 ! the basis for the interface
                                                                                                 ! tangential singular vectors
                        sigma_tilde(3,2) = sigma_tilde(2,3) ! Symmetric
                        sigma_tilde(3,3) = sum( v_tilde(:,3) * matmul( sigma, v_tilde(:,3) ) )

                        ! sigma_tilde = matmul(sigma_tilde, v_tilde)
                        ! sigma_tilde = matmul(transpose(v_tilde), sigma_tilde)
                        sigma_tilde = matmul(sigma_tilde, transpose(v_tilde))
                        sigma_tilde = matmul(v_tilde, sigma_tilde)

                        ! if ((i == 21) .and. (j == 64)) then
                        !     print *, "theta = ", theta
                        !     print *, "sigma_tilde:"
                        !     print *, "   ",  sigma_tilde(1,1), sigma_tilde(1,2), sigma_tilde(1,3)
                        !     print *, "   ",  sigma_tilde(2,1), sigma_tilde(2,2), sigma_tilde(2,3)
                        !     print *, "   ",  sigma_tilde(3,1), sigma_tilde(3,2), sigma_tilde(3,3)
                        !     print *, "tangential: ", sum( v_tilde(:,1) * matmul( sigma_tilde, v_tilde(:,2) ) ), &
                        !                              sum( v_tilde(:,1) * matmul( sigma_tilde, v_tilde(:,3) ) )
                        ! end if

                        ! Try to make it a smoother transition to sliding
                        ! sigma_tilde = mask(i,j,k)*sigma_tilde + (one - mask(i,j,k))*sigma

                        ! Get eigenvalues and eigenvectors of sigma
                        call dsyev('V', 'U', 3, sigma_tilde, 3, Sa, svdwork, lwork, info)
                        if(info .ne. 0) then
                            print '(A,I6,A)', 'proc ', nrank, ': Problem with DSYEV. Please check.'
                            print *, "sigma:"
                            print *, "   ",  sigma(1,1), sigma(1,2), sigma(1,3)
                            print *, "   ",  sigma(2,1), sigma(2,2), sigma(2,3)
                            print *, "   ",  sigma(3,1), sigma(3,2), sigma(3,3)
                            print *, "tangential: ", sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,2) ) ), &
                                                     sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,3) ) )
                        end if

                        ! if ((i == 108) .and. (j == 2)) then
                        !     print *, "v_tilde_new:"
                        !     print *, "   ",  sigma_tilde(1,1), sigma_tilde(1,2), sigma_tilde(1,3)
                        !     print *, "   ",  sigma_tilde(2,1), sigma_tilde(2,2), sigma_tilde(2,3)
                        !     print *, "   ",  sigma_tilde(3,1), sigma_tilde(3,2), sigma_tilde(3,3)
                        ! end if

                        ! Now get new beta
                        f = Sa / (this%elastic%mu*sqrt_om); f(3) = beta(1)*beta(2)*beta(3)     ! New function value (target to attain)
                        
                        betasum = sum( beta*(beta-one) ) / three
                        f1 = -( beta*(beta-one) - betasum ); f1(3) = beta(1)*beta(2)*beta(3)   ! Original function value

                        ! Get newton step
                        gradf(1,1) = -twothird*(two*beta(1)-one); gradf(1,2) =     third*(two*beta(2)-one); gradf(1,3) = third*(two*beta(3)-one)
                        gradf(2,1) =     third*(two*beta(1)-one); gradf(2,2) = -twothird*(two*beta(2)-one); gradf(2,3) = third*(two*beta(3)-one)
                        gradf(3,1) = beta(2)*beta(3);             gradf(3,2) = beta(3)*beta(1);             gradf(3,3) = beta(1)*beta(2)

                        dbeta = (f-f1)
                        call dgesv(3, 1, gradf, 3, ipiv, dbeta, 3, info)
                   
                        ! Compute residual
                        residual = -sum( (f1-f)*dbeta )                                    ! lambda**2
                        iters = 0
                        t = 1._rkind
                        do while ( (iters < niters) .AND. (abs(residual) .GT. tol) )
                            ! Backtracking line search
                            t = 1._rkind
                            beta_new = beta + t * dbeta

                            ! Get new residual
                            gradf_new(1,1) = -twothird*(two*beta_new(1)-one);
                            gradf_new(1,2) =     third*(two*beta_new(2)-one);
                            gradf_new(1,3) = third*(two*beta_new(3)-one)

                            gradf_new(2,1) =     third*(two*beta_new(1)-one);
                            gradf_new(2,2) = -twothird*(two*beta_new(2)-one);
                            gradf_new(2,3) = third*(two*beta_new(3)-one)

                            gradf_new(3,1) = beta_new(2)*beta_new(3);
                            gradf_new(3,2) = beta_new(3)*beta_new(1);
                            gradf_new(3,3) = beta_new(1)*beta_new(2)

                            betasum = sum( beta_new*(beta_new-one) ) / three
                            f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)
                            dbeta_new = (f-f2)
                            call dgesv(3, 1, gradf_new, 3, ipiv, dbeta_new, 3, info)
                            residual_new = -sum( (f2-f)*dbeta_new )                                    ! lambda**2

                            do while ( (abs(residual_new) .GE. abs(residual)) .AND. (t > eps) )
                                if (iters .GT. (niters - 10)) then
                                    print '(A,I0,3(A,ES15.5))', 'iters = ', iters, ', t = ', t, ', residual_new = ', residual_new, ', residual = ', residual
                                end if

                                t = half*t
                                beta_new = beta + t * dbeta

                                gradf_new(1,1) = -twothird*(two*beta_new(1)-one);
                                gradf_new(1,2) =     third*(two*beta_new(2)-one);
                                gradf_new(1,3) = third*(two*beta_new(3)-one)

                                gradf_new(2,1) =     third*(two*beta_new(1)-one);
                                gradf_new(2,2) = -twothird*(two*beta_new(2)-one);
                                gradf_new(2,3) = third*(two*beta_new(3)-one)

                                gradf_new(3,1) = beta_new(2)*beta_new(3);
                                gradf_new(3,2) = beta_new(3)*beta_new(1);
                                gradf_new(3,3) = beta_new(1)*beta_new(2)

                                betasum = sum( beta_new*(beta_new-one) ) / three
                                f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)

                                dbeta_new = (f-f2)
                                call dgesv(3, 1, gradf_new, 3, ipiv, dbeta_new, 3, info)
                                residual_new = -sum( (f2-f)*dbeta_new )                                    ! lambda**2
                            end do
                            beta = beta_new
                            f1 = f2
                            dbeta = dbeta_new
                            residual = residual_new

                            iters = iters + 1
                            if (t <= eps) then
                                print '(A)', 'Newton solve in sliding_deformation did not converge'
                                exit
                            end if
                        end do
                        if ((iters >= niters) .OR. (t <= eps)) then
                            write(charout,'(4(A,I0))') 'Newton solve in sliding_deformation did not converge at index ',i,',',j,',',k,' of process ',nrank
                            print '(A)', charout
                            print '(A,3(X,I0))', 'sty = ', this%decomp%yst
                            print '(A)', 'g = '
                            print '(4X,3(ES15.5))', this%g(i,j,k,1), this%g(i,j,k,2), this%g(i,j,k,3)
                            print '(4X,3(ES15.5))', this%g(i,j,k,4), this%g(i,j,k,5), this%g(i,j,k,6)
                            print '(4X,3(ES15.5))', this%g(i,j,k,7), this%g(i,j,k,8), this%g(i,j,k,9)
                            print '(A,ES15.5)', '( ||S||^2 - (2/3) sigma_Y^2 )/mu^2 = ', ycrit

                            print '(A,ES15.5)', 'Relaxation, t = ', t
                            print '(A,ES15.5)', 'Residual = ', residual
                            call GracefulExit(charout,6382)
                        end if

                        ! Then get new svals
                        sval = sqrt(beta) * sqrt_om**(one/three)

                        if (this%use_gTg) then
                            sval = sval*sval ! New eigenvalues of G
                            
                            ! Get g = v*sval*vt
                            u = sigma_tilde; vt = transpose(u)
                            vt(1,:) = vt(1,:)*sval(1); vt(2,:) = vt(2,:)*sval(2); vt(3,:) = vt(3,:)*sval(3)  ! eigval*vt
                            G = MATMUL(u,vt) ! v*eigval*vt
                        else
                            ! Get g = u*sval*vt
                            u = sigma_tilde; vt = transpose(u)
                            vt(1,:) = vt(1,:)*sval(1); vt(2,:) = vt(2,:)*sval(2); vt(3,:) = vt(3,:)*sval(3)  ! sval*vt
                            g = MATMUL(u,vt) ! u*sval*vt
                        end if


                        ! Try to make it a smoother transition to sliding
                        this%g(i,j,k,1) = mask(i,j,k)*g(1,1) + (one - mask(i,j,k))*this%g(i,j,k,1)
                        this%g(i,j,k,2) = mask(i,j,k)*g(1,2) + (one - mask(i,j,k))*this%g(i,j,k,2)
                        this%g(i,j,k,3) = mask(i,j,k)*g(1,3) + (one - mask(i,j,k))*this%g(i,j,k,3)
                        this%g(i,j,k,4) = mask(i,j,k)*g(2,1) + (one - mask(i,j,k))*this%g(i,j,k,4)
                        this%g(i,j,k,5) = mask(i,j,k)*g(2,2) + (one - mask(i,j,k))*this%g(i,j,k,5)
                        this%g(i,j,k,6) = mask(i,j,k)*g(2,3) + (one - mask(i,j,k))*this%g(i,j,k,6)
                        this%g(i,j,k,7) = mask(i,j,k)*g(3,1) + (one - mask(i,j,k))*this%g(i,j,k,7)
                        this%g(i,j,k,8) = mask(i,j,k)*g(3,2) + (one - mask(i,j,k))*this%g(i,j,k,8)
                        this%g(i,j,k,9) = mask(i,j,k)*g(3,3) + (one - mask(i,j,k))*this%g(i,j,k,9)

                        ! this%g(i,j,k,1) = g(1,1); this%g(i,j,k,2) = g(1,2); this%g(i,j,k,3) = g(1,3)
                        ! this%g(i,j,k,4) = g(2,1); this%g(i,j,k,5) = g(2,2); this%g(i,j,k,6) = g(2,3)
                        ! this%g(i,j,k,7) = g(3,1); this%g(i,j,k,8) = g(3,2); this%g(i,j,k,9) = g(3,3)

                    end if
                end do
            end do
        end do

        ! ! Get devstress for debugging
        ! call this%get_eelastic_devstress()
        ! sigma(1,1) = this%sxx(21,64,1); sigma(1,2) = this%sxy(21,64,1); sigma(1,3) = this%sxz(21,64,1);
        ! sigma(2,1) = this%sxy(21,64,1); sigma(2,2) = this%syy(21,64,1); sigma(2,3) = this%syz(21,64,1);
        ! sigma(3,1) = this%sxz(21,64,1); sigma(3,2) = this%syz(21,64,1); sigma(3,3) = this%szz(21,64,1);
        ! print *, "New sigma:"
        ! print *, "   ",  sigma(1,1), sigma(1,2), sigma(1,3)
        ! print *, "   ",  sigma(2,1), sigma(2,2), sigma(2,3)
        ! print *, "   ",  sigma(3,1), sigma(3,2), sigma(3,3)

        ! ! 2D circular
        ! x = - half + real( this%decomp%yst(1) - 1 + 21 - 1, rkind ) * dx
        ! y = - half + real( this%decomp%yst(2) - 1 + 64 - 1, rkind ) * dy
        ! rad = sqrt(x**2 + y**2)
        ! theta = atan2(y,x)

        ! v_tilde(1,1) = cos(theta); v_tilde(1,2) =-sin(theta); v_tilde(1,3) = zero
        ! v_tilde(2,1) = sin(theta); v_tilde(2,2) = cos(theta); v_tilde(2,3) = zero
        ! v_tilde(3,1) = zero;       v_tilde(3,2) = zero;       v_tilde(3,3) = one

        ! print *, "tangential: ", sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,2) ) ), &
        !                          sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,3) ) )


        deallocate(svdwork)
    end subroutine

end module
