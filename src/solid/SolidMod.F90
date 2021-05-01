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
        integer :: nxp, nyp, nzp, ns

        logical :: plast = .FALSE.
        logical :: explPlast = .FALSE.
        logical :: PTeqb = .TRUE., pEqb = .FALSE., pRelax = .FALSE., updateEtot = .TRUE., includeSources = .FALSE.
        logical :: use_gTg = .FALSE.,useOneG = .FALSE.,intSharp = .FALSE.,intSharp_spf = .TRUE.,intSharp_ufv = .TRUE.,intSharp_d02 = .TRUE.,strainHard = .TRUE.,cnsrv_g = .FALSE.,cnsrv_gt = .FALSE.,cnsrv_gp = .FALSE.,cnsrv_pe = .FALSE., intSharp_cpg_west = .FALSE.

        class(stiffgas ), allocatable :: hydro
        class(sep1solid), allocatable :: elastic

        type(decomp_info), pointer :: decomp
        type(derivatives), pointer :: der,derD02
        type(filters),     pointer :: fil
        type(filters),     pointer :: gfil

        real(rkind), dimension(:,:,:), allocatable :: Ys
        real(rkind), dimension(:,:,:), allocatable :: VF
        real(rkind), dimension(:,:,:), allocatable :: eh
        real(rkind), dimension(:,:,:), allocatable :: eel

        real(rkind), dimension(:,:,:,:), allocatable :: g,g_t,g_p,rg,rg_t,rg_p
        real(rkind), dimension(:,:,:),   allocatable :: e_p,e_pp,pe,rpe
        real(rkind), dimension(:,:,:),   allocatable :: curl_e,curl_t,curl_p,det_e,det_t,det_p
        real(rkind), dimension(:,:,:),   pointer     :: g11,gt11,gp11
        real(rkind), dimension(:,:,:),   pointer     :: g12,gt12,gp12
        real(rkind), dimension(:,:,:),   pointer     :: g13,gt13,gp13
        real(rkind), dimension(:,:,:),   pointer     :: g21,gt21,gp21
        real(rkind), dimension(:,:,:),   pointer     :: g22,gt22,gp22
        real(rkind), dimension(:,:,:),   pointer     :: g23,gt23,gp23
        real(rkind), dimension(:,:,:),   pointer     :: g31,gt31,gp31
        real(rkind), dimension(:,:,:),   pointer     :: g32,gt32,gp32
        real(rkind), dimension(:,:,:),   pointer     :: g33,gt33,gp33

        real(rkind), dimension(:,:,:),   pointer     :: rg11,rgt11,rgp11
        real(rkind), dimension(:,:,:),   pointer     :: rg12,rgt12,rgp12
        real(rkind), dimension(:,:,:),   pointer     :: rg13,rgt13,rgp13
        real(rkind), dimension(:,:,:),   pointer     :: rg21,rgt21,rgp21
        real(rkind), dimension(:,:,:),   pointer     :: rg22,rgt22,rgp22
        real(rkind), dimension(:,:,:),   pointer     :: rg23,rgt23,rgp23
        real(rkind), dimension(:,:,:),   pointer     :: rg31,rgt31,rgp31
        real(rkind), dimension(:,:,:),   pointer     :: rg32,rgt32,rgp32
        real(rkind), dimension(:,:,:),   pointer     :: rg33,rgt33,rgp33
        
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
        real(rkind), dimension(:,:,:),   allocatable :: diff_g
        real(rkind), dimension(:,:,:),   allocatable :: diff_gt
        real(rkind), dimension(:,:,:),   allocatable :: diff_gp
        real(rkind), dimension(:,:,:),   allocatable :: diff_pe
        real(rkind), dimension(:,:,:,:), allocatable :: Ji

        ! species-specific variables for interface sharpening
        real(rkind), dimension(:,:,:,:), allocatable :: intSharp_a,intSharp_R,intSharp_aDiff,intSharp_RDiff
        real(rkind), dimension(:,:,:),   allocatable :: intSharp_aFV,intSharp_RFV
        real(rkind), dimension(:),       allocatable :: intSharp_ysc
        real(rkind) :: intSharp_cut
        real(rkind), dimension(:,:,:,:,:),allocatable :: intSharp_rg,intSharp_rgDiff,intSharp_rgt,intSharp_rgtDiff,intSharp_rgp,intSharp_rgpDiff
        real(rkind), dimension(:,:,:,:), allocatable ::  intSharp_rgFV,intSharp_rgtFV,intSharp_rgpFV,intSharp_gFV,intSharp_gtFV,intSharp_gpFV

        ! species-specific conserved variables
        real(rkind), dimension(:,:,:,:), allocatable :: consrv

        ! work arrays
        real(rkind), dimension(:,:,:,:), allocatable :: Qtmpg,Qtmpg_t,Qtmpg_p
        real(rkind), dimension(:,:,:),   allocatable :: QtmpYs
        real(rkind), dimension(:,:,:),   allocatable :: Qtmpeh
        real(rkind), dimension(:,:,:),   allocatable :: QtmpVF
        real(rkind), dimension(:,:,:),   allocatable :: Qtmppe
        real(rkind), dimension(:,:,:), allocatable :: modDevSigma
    contains

        procedure :: init
        procedure :: getRHS_g
        procedure :: getRHS_rg
        procedure :: getRHS_gt
        procedure :: getRHS_rgt
        procedure :: getRHS_gp
        procedure :: getRHS_rgp

        procedure :: getRHS_gTg
        procedure :: getRHS_rgTg
        procedure :: getRHS_gtTgt
        procedure :: getRHS_rgtTgt
        procedure :: getRHS_gpTgp
        procedure :: getRHS_rgpTgp

        procedure :: getRHS_pe
        procedure :: getRHS_rpe
        !procedure :: getRHS_gTg
        procedure :: getRHS_Ys
        procedure :: getRHS_Ys_intSharp
        procedure :: getRHS_eh
        procedure :: getRHS_VF
        procedure :: update_g
        !procedure :: update_gTg
        procedure :: update_Ys
        procedure :: update_eh
        procedure :: update_VF
        procedure :: getPhysicalProperties
        procedure :: getPlasticSources
        procedure :: implicit_plastic
        procedure :: get_p_from_ehydro
        procedure :: get_ehydroT_from_p
        procedure :: get_eelastic_devstress
        !procedure :: get_eelastic_devstress_mixture
        procedure :: getModifiedModulii
        procedure :: get_conserved
        procedure :: get_conserved_g
        procedure :: get_primitive
        procedure :: get_primitive_g
        procedure :: getSpeciesDensity
        procedure :: getSpeciesDensity_from_g
        procedure :: get_enthalpy
        procedure :: checkNaN
        procedure :: filter
        procedure :: filter_g
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
    subroutine init(this,decomp,der,derD02,fil,gfil,PTeqb,pEqb,pRelax,use_gTg,useOneG,intSharp,intSharp_spf,intSharp_ufv,intSharp_d02,intSharp_cut,intSharp_cpg_west,updateEtot,strainHard,cnsrv_g,cnsrv_gt,cnsrv_gp,cnsrv_pe,ns)
        class(solid), target, intent(inout) :: this
        type(decomp_info), target, intent(in) :: decomp
        type(derivatives), target, intent(in) :: der,derD02
        type(filters),     target, intent(in) :: fil, gfil
        logical, intent(in) :: PTeqb,pEqb,pRelax,updateEtot
        logical, intent(in) :: use_gTg,useOneG,intSharp,intSharp_spf,intSharp_ufv,intSharp_d02,intSharp_cpg_west,strainHard,cnsrv_g,cnsrv_gt,cnsrv_gp,cnsrv_pe
        integer, intent(in) :: ns
        real(rkind), intent(in) :: intSharp_cut

        this%decomp => decomp
        this%der  => der
        this%derD02  => derD02
        this%fil  => fil
        this%gfil => gfil
       
        this%PTeqb  = PTeqb
        this%pEqb   = pEqb
        this%pRelax = pRelax

        this%use_gTg = use_gTg
        this%useOneG = useOneG
        this%intSharp = intSharp
        this%intSharp_spf = intSharp_spf
        this%intSharp_ufv = intSharp_ufv
        this%intSharp_d02 = intSharp_d02
        this%intSharp_cut = intSharp_cut
        this%intSharp_cpg_west = intSharp_cpg_west
        this%updateEtot  = updateEtot

        this%strainHard = strainHard
        this%cnsrv_g  = cnsrv_g
        this%cnsrv_gt = cnsrv_gt
        this%cnsrv_gp = cnsrv_gp
        this%cnsrv_pe = cnsrv_pe

        this%ns = ns

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
        
        ! Allocate material elastic inverse deformation gradients and associate pointers
        if( allocated( this%g ) ) deallocate( this%g)
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

        ! Allocate material plastic inverse deformation gradients and associate pointers
        if( allocated( this%g_t ) ) deallocate( this%g_t)
        allocate( this%g_t(this%nxp,this%nyp,this%nzp,9) )
        this%gt11 => this%g_t(:,:,:,1)   
        this%gt12 => this%g_t(:,:,:,2)   
        this%gt13 => this%g_t(:,:,:,3)   
        this%gt21 => this%g_t(:,:,:,4)   
        this%gt22 => this%g_t(:,:,:,5)   
        this%gt23 => this%g_t(:,:,:,6)   
        this%gt31 => this%g_t(:,:,:,7)   
        this%gt32 => this%g_t(:,:,:,8)   
        this%gt33 => this%g_t(:,:,:,9)   

        if( allocated( this%g_p ) ) deallocate( this%g_p)
        allocate( this%g_p(this%nxp,this%nyp,this%nzp,9) )
        this%gp11 => this%g_p(:,:,:,1)   
        this%gp12 => this%g_p(:,:,:,2)   
        this%gp13 => this%g_p(:,:,:,3)   
        this%gp21 => this%g_p(:,:,:,4)   
        this%gp22 => this%g_p(:,:,:,5)   
        this%gp23 => this%g_p(:,:,:,6)   
        this%gp31 => this%g_p(:,:,:,7)   
        this%gp32 => this%g_p(:,:,:,8)   
        this%gp33 => this%g_p(:,:,:,9)   

        if( allocated( this%e_p ) ) deallocate( this%e_p) !integrated plastic Eulerian-Almansi strain
        allocate( this%e_p(this%nxp,this%nyp,this%nzp) )
        this%e_p=0.0d0

        if( allocated( this%e_pp ) ) deallocate( this%e_pp) !integrated plastic Eulerian-Almansi strain
        allocate( this%e_pp(this%nxp,this%nyp,this%nzp) )
        this%e_pp=0.0d0

        if( allocated( this%pe ) ) deallocate( this%pe) !integrated plastic Eulerian-Almansi strain
        allocate( this%pe(this%nxp,this%nyp,this%nzp) )
        this%pe=0.0d0

        if( allocated( this%curl_e ) ) deallocate( this%curl_e) !curl compatibility rhs
        allocate( this%curl_e(this%nxp,this%nyp,this%nzp) )
        this%curl_e=0.0d0

        if( allocated( this%curl_t ) ) deallocate( this%curl_t) !curl compatibility rhs
        allocate( this%curl_t(this%nxp,this%nyp,this%nzp) )
        this%curl_t=0.0d0

        if( allocated( this%curl_p ) ) deallocate( this%curl_p) !curl compatibility rhs
        allocate( this%curl_p(this%nxp,this%nyp,this%nzp) )
        this%curl_p=0.0d0

        if( allocated( this%det_e ) ) deallocate( this%det_e) !det compatibility rhs
        allocate( this%det_e(this%nxp,this%nyp,this%nzp) )
        this%det_e=0.0d0

        if( allocated( this%det_t ) ) deallocate( this%det_t) !det compatibility rhs
        allocate( this%det_t(this%nxp,this%nyp,this%nzp) )
        this%det_t=0.0d0

        if( allocated( this%det_p ) ) deallocate( this%det_p) !det compatibility rhs
        allocate( this%det_p(this%nxp,this%nyp,this%nzp) )
        this%det_p=0.0d0


        !conservative form quantities
        ! Allocate material elastic inverse deformation gradients and associate pointers
        if( allocated( this%rg ) ) deallocate( this%rg)
        allocate( this%rg(this%nxp,this%nyp,this%nzp,9) )
        this%rg11 => this%rg(:,:,:,1)   
        this%rg12 => this%rg(:,:,:,2)   
        this%rg13 => this%rg(:,:,:,3)   
        this%rg21 => this%rg(:,:,:,4)   
        this%rg22 => this%rg(:,:,:,5)   
        this%rg23 => this%rg(:,:,:,6)   
        this%rg31 => this%rg(:,:,:,7)   
        this%rg32 => this%rg(:,:,:,8)   
        this%rg33 => this%rg(:,:,:,9)   

        ! Allocate material plastic inverse deformation gradients and associate pointers
        if( allocated( this%rg_t ) ) deallocate( this%rg_t)
        allocate( this%rg_t(this%nxp,this%nyp,this%nzp,9) )
        this%rgt11 => this%rg_t(:,:,:,1)   
        this%rgt12 => this%rg_t(:,:,:,2)   
        this%rgt13 => this%rg_t(:,:,:,3)   
        this%rgt21 => this%rg_t(:,:,:,4)   
        this%rgt22 => this%rg_t(:,:,:,5)   
        this%rgt23 => this%rg_t(:,:,:,6)   
        this%rgt31 => this%rg_t(:,:,:,7)   
        this%rgt32 => this%rg_t(:,:,:,8)   
        this%rgt33 => this%rg_t(:,:,:,9)   

        if( allocated( this%rg_p ) ) deallocate( this%rg_p)
        allocate( this%rg_p(this%nxp,this%nyp,this%nzp,9) )
        this%rgp11 => this%rg_p(:,:,:,1)   
        this%rgp12 => this%rg_p(:,:,:,2)   
        this%rgp13 => this%rg_p(:,:,:,3)   
        this%rgp21 => this%rg_p(:,:,:,4)   
        this%rgp22 => this%rg_p(:,:,:,5)   
        this%rgp23 => this%rg_p(:,:,:,6)   
        this%rgp31 => this%rg_p(:,:,:,7)   
        this%rgp32 => this%rg_p(:,:,:,8)   
        this%rgp33 => this%rg_p(:,:,:,9)   

        if( allocated( this%rpe ) ) deallocate( this%rpe)
        allocate( this%rpe(this%nxp,this%nyp,this%nzp) )
        this%rpe=0.0d0


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

        ! Allocate g_e diffusivity array
        if( allocated( this%diff_g ) ) deallocate( this%diff_g )
        allocate( this%diff_g(this%nxp,this%nyp,this%nzp) )

        ! Allocate g_t diffusivity array
        if( allocated( this%diff_gt ) ) deallocate( this%diff_gt )
        allocate( this%diff_gt(this%nxp,this%nyp,this%nzp) )

        ! Allocate g_p diffusivity array
        if( allocated( this%diff_gp ) ) deallocate( this%diff_gp )
        allocate( this%diff_gp(this%nxp,this%nyp,this%nzp) )

        ! Allocate pe diffusivity array
        if( allocated( this%diff_pe ) ) deallocate( this%diff_pe )
        allocate( this%diff_pe(this%nxp,this%nyp,this%nzp) )

        ! Allocate material diffusive flux
        if( allocated( this%Ji ) ) deallocate( this%Ji )
        allocate( this%Ji(this%nxp,this%nyp,this%nzp,3) )

        ! Allocate interface sharpening volume fraction flux
        if( allocated( this%intSharp_a ) ) deallocate( this%intSharp_a )
        allocate( this%intSharp_a(this%nxp,this%nyp,this%nzp,3) )

        ! Allocate interface sharpening mass fraction flux
        if( allocated( this%intSharp_R ) ) deallocate( this%intSharp_R )
        allocate( this%intSharp_R(this%nxp,this%nyp,this%nzp,3) )

        ! Allocate interface sharpening volume fraction diffusion flux
        if( allocated( this%intSharp_aDiff ) ) deallocate( this%intSharp_aDiff )
        allocate( this%intSharp_aDiff(this%nxp,this%nyp,this%nzp,3) )

        ! Allocate interface sharpening mass fraction diffusion flux
        if( allocated( this%intSharp_RDiff ) ) deallocate( this%intSharp_RDiff )
        allocate( this%intSharp_RDiff(this%nxp,this%nyp,this%nzp,3) )

        ! Allocate interface sharpening volume fraction FV flux
        if( allocated( this%intSharp_aFV ) ) deallocate( this%intSharp_aFV )
        allocate( this%intSharp_aFV(this%nxp,this%nyp,this%nzp) )

        ! Allocate interface sharpening mass fraction FV flux
        if( allocated( this%intSharp_RFV ) ) deallocate( this%intSharp_RFV )
        allocate( this%intSharp_RFV(this%nxp,this%nyp,this%nzp) )

        ! Allocate interface sharpening rho*g flux
        if( allocated( this%intSharp_rg ) ) deallocate( this%intSharp_rg )
        allocate( this%intSharp_rg(this%nxp,this%nyp,this%nzp,9,3) )
        ! Allocate interface sharpening rho*g diffusion flux
        if( allocated( this%intSharp_rgDiff ) ) deallocate( this%intSharp_rgDiff )
        allocate( this%intSharp_rgDiff(this%nxp,this%nyp,this%nzp,9,3) )
        ! Allocate interface sharpening rho*g FV flux
        if( allocated( this%intSharp_rgFV ) ) deallocate( this%intSharp_rgFV )
        allocate( this%intSharp_rgFV(this%nxp,this%nyp,this%nzp,9) )
        ! Allocate interface sharpening g FV flux
        if( allocated( this%intSharp_gFV ) ) deallocate( this%intSharp_gFV )
        allocate( this%intSharp_gFV(this%nxp,this%nyp,this%nzp,9) )


        ! Allocate interface sharpening rho*gt flux
        if( allocated( this%intSharp_rgt ) ) deallocate( this%intSharp_rgt )
        allocate( this%intSharp_rgt(this%nxp,this%nyp,this%nzp,9,3) )
        ! Allocate interface sharpening rho*gt diffusion flux
        if( allocated( this%intSharp_rgtDiff ) ) deallocate( this%intSharp_rgtDiff )
        allocate( this%intSharp_rgtDiff(this%nxp,this%nyp,this%nzp,9,3) )
        ! Allocate interface sharpening rho*gt FV flux
        if( allocated( this%intSharp_rgtFV ) ) deallocate( this%intSharp_rgtFV )
        allocate( this%intSharp_rgtFV(this%nxp,this%nyp,this%nzp,9) )
        ! Allocate interface sharpening gt FV flux
        if( allocated( this%intSharp_gtFV ) ) deallocate( this%intSharp_gtFV )
        allocate( this%intSharp_gtFV(this%nxp,this%nyp,this%nzp,9) )

        ! Allocate interface sharpening rho*gp flux
        if( allocated( this%intSharp_rgp ) ) deallocate( this%intSharp_rgp )
        allocate( this%intSharp_rgp(this%nxp,this%nyp,this%nzp,9,3) )
        ! Allocate interface sharpening rho*gp diffusion flux
        if( allocated( this%intSharp_rgpDiff ) ) deallocate( this%intSharp_rgpDiff )
        allocate( this%intSharp_rgpDiff(this%nxp,this%nyp,this%nzp,9,3) )
        ! Allocate interface sharpening rho*gp FV flux
        if( allocated( this%intSharp_rgpFV ) ) deallocate( this%intSharp_rgpFV )
        allocate( this%intSharp_rgpFV(this%nxp,this%nyp,this%nzp,9) )
        ! Allocate interface sharpening gp FV flux
        if( allocated( this%intSharp_gpFV ) ) deallocate( this%intSharp_gpFV )
        allocate( this%intSharp_gpFV(this%nxp,this%nyp,this%nzp,9) )


        ! Allocate interface sharpening Ys bound
        if( allocated( this%intSharp_ysc ) ) deallocate( this%intSharp_ysc )
        allocate( this%intSharp_ysc(2) )
        this%intSharp_ysc(1) = this%intSharp_cut*1.D-3 !default values -- ok up to density ratio 1000
        this%intSharp_ysc(2) = one-this%intSharp_cut*1.D-3

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

        ! g_t tensor equation
        if( allocated( this%Qtmpg_t ) ) deallocate( this%Qtmpg_t )
        allocate( this%Qtmpg_t(this%nxp,this%nyp,this%nzp,9) )
        
        ! g_p tensor equation
        if( allocated( this%Qtmpg_p ) ) deallocate( this%Qtmpg_p )
        allocate( this%Qtmpg_p(this%nxp,this%nyp,this%nzp,9) )
        
        ! Ys equation
        if( allocated( this%QtmpYs ) ) deallocate( this%QtmpYs )
        allocate( this%QtmpYs(this%nxp,this%nyp,this%nzp) )

        ! pe equation
        if( allocated( this%Qtmppe ) ) deallocate( this%Qtmppe )
        allocate( this%Qtmppe(this%nxp,this%nyp,this%nzp) )

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
        if( allocated( this%Qtmpg_t  ) ) deallocate( this%Qtmpg_t )
        if( allocated( this%Qtmpg_p  ) ) deallocate( this%Qtmpg_p )
        if( allocated( this%Qtmppe  ) ) deallocate( this%Qtmppe )
        if( allocated( this%consrv ) ) deallocate( this%consrv )
        
        if( allocated( this%intSharp_a )   ) deallocate( this%intSharp_a )
        if( allocated( this%intSharp_R )   ) deallocate( this%intSharp_R )
        if( allocated( this%intSharp_aDiff )   ) deallocate( this%intSharp_aDiff )
        if( allocated( this%intSharp_RDiff )   ) deallocate( this%intSharp_RDiff )
        if( allocated( this%intSharp_aFV )   ) deallocate( this%intSharp_aFV )
        if( allocated( this%intSharp_RFV )   ) deallocate( this%intSharp_RFV )
        if( allocated( this%intSharp_rg )   ) deallocate( this%intSharp_rg )
        if( allocated( this%intSharp_rgDiff )   ) deallocate( this%intSharp_rgDiff )
        if( allocated( this%intSharp_rgFV )   ) deallocate( this%intSharp_rgFV )
        if( allocated( this%intSharp_gFV )   ) deallocate( this%intSharp_gFV )
        if( allocated( this%intSharp_rgt )   ) deallocate( this%intSharp_rgt )
        if( allocated( this%intSharp_rgtDiff )   ) deallocate( this%intSharp_rgtDiff )
        if( allocated( this%intSharp_rgtFV )   ) deallocate( this%intSharp_rgtFV )
        if( allocated( this%intSharp_gtFV )   ) deallocate( this%intSharp_gtFV )
        if( allocated( this%intSharp_rgp )   ) deallocate( this%intSharp_rgp )
        if( allocated( this%intSharp_rgpDiff )   ) deallocate( this%intSharp_rgpDiff )
        if( allocated( this%intSharp_rgpFV )   ) deallocate( this%intSharp_rgpFV )
        if( allocated( this%intSharp_gpFV )   ) deallocate( this%intSharp_gpFV )
        if( allocated( this%intSharp_ysc )   ) deallocate( this%intSharp_ysc )

        if( allocated( this%Ji )   ) deallocate( this%Ji )
        if( allocated( this%diff ) ) deallocate( this%diff )
        if( allocated( this%diff_g ) ) deallocate( this%diff_g )
        if( allocated( this%diff_gt ) ) deallocate( this%diff_gt )
        if( allocated( this%diff_gp ) ) deallocate( this%diff_gp )
        if( allocated( this%diff_pe ) ) deallocate( this%diff_pe )
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

        nullify( this%gt11 ); nullify( this%gt12 ); nullify( this%gt13 )
        nullify( this%gt21 ); nullify( this%gt22 ); nullify( this%gt23 )
        nullify( this%gt31 ); nullify( this%gt32 ); nullify( this%gt33 )
        nullify( this%gp11 ); nullify( this%gp12 ); nullify( this%gp13 )
        nullify( this%gp21 ); nullify( this%gp22 ); nullify( this%gp23 )
        nullify( this%gp31 ); nullify( this%gp32 ); nullify( this%gp33 )
        if( allocated( this%g_t )         ) deallocate( this%g_t )
        if( allocated( this%g_p )         ) deallocate( this%g_p )
        if( allocated( this%e_p )         ) deallocate( this%e_p )
        if( allocated( this%e_pp )         ) deallocate( this%e_pp )
        if( allocated( this%pe )         ) deallocate( this%pe )
        if( allocated( this%curl_e )         ) deallocate( this%curl_e )
        if( allocated( this%curl_t )         ) deallocate( this%curl_t )
        if( allocated( this%curl_p )         ) deallocate( this%curl_p )
        if( allocated( this%det_e )         ) deallocate( this%det_e )
        if( allocated( this%det_t )         ) deallocate( this%det_t )
        if( allocated( this%det_p )         ) deallocate( this%det_p )

        if( allocated( this%rg )         ) deallocate( this%rg )
        if( allocated( this%rg_t )         ) deallocate( this%rg_t )
        if( allocated( this%rg_p )         ) deallocate( this%rg_p )
        if( allocated( this%rpe )         ) deallocate( this%rpe )
        nullify( this%rg11 ); nullify( this%rg12 ); nullify( this%rg13 )
        nullify( this%rg21 ); nullify( this%rg22 ); nullify( this%rg23 )
        nullify( this%rg31 ); nullify( this%rg32 ); nullify( this%rg33 )
        nullify( this%rgt11 ); nullify( this%rgt12 ); nullify( this%rgt13 )
        nullify( this%rgt21 ); nullify( this%rgt22 ); nullify( this%rgt23 )
        nullify( this%rgt31 ); nullify( this%rgt32 ); nullify( this%rgt33 )
        nullify( this%rgp11 ); nullify( this%rgp12 ); nullify( this%rgp13 )
        nullify( this%rgp21 ); nullify( this%rgp22 ); nullify( this%rgp23 )
        nullify( this%rgp31 ); nullify( this%rgp32 ); nullify( this%rgp33 )



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
        nullify( this%derD02    )
        nullify( this%decomp )
    end subroutine

    subroutine getPhysicalProperties(this)
        class(solid), intent(inout) :: this

        this%diff = zero
        this%diff_g = zero
        this%diff_gt = zero
        this%diff_gp = zero
        this%diff_pe = zero
        this%kap = zero

    end subroutine

    subroutine getSpeciesDensity_from_g(this,rho)
        class(solid), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(out)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: detg

        detg = this%g11*(this%g22*this%g33-this%g23*this%g32) &
             - this%g12*(this%g21*this%g33-this%g31*this%g23) &
             + this%g13*(this%g21*this%g32-this%g31*this%g22)

        if (this%use_gTg.and.(.not.this%strainHard)) then
            detg = sqrt(detg)
        end if

        rho = detg*this%elastic%rho0

    end subroutine

    subroutine update_g(this,isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,rho0mix,mumix,yieldmix,solidVF)
        use constants,  only: eps
        use RKCoeffs,   only: RK45_A,RK45_B
        use reductions, only: P_MAXVAL,P_MINVAL
        use operators,  only: gradient, filter3D
        use exits,      only : nancheck
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9)  :: duidxj
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in), optional  :: rho0mix, mumix, yieldmix, solidVF

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg,rhsgt,rhsgp   ! RHS for g tensor equation
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: rhspe, tmpt, detgt
        real(rkind) :: max_modDevSigma
        integer :: i,j,k,l
        character(len=clen) :: charout

        rhsg = zero
        rhsgt = zero
        rhsgp = zero
        rhspe = zero

        call gradient(this%decomp,this%der,u, duidxj(:,:,:,1), duidxj(:,:,:,2), duidxj(:,:,:,3), -x_bc,  y_bc,  z_bc)
        call gradient(this%decomp,this%der,v, duidxj(:,:,:,4), duidxj(:,:,:,5), duidxj(:,:,:,6),  x_bc, -y_bc,  z_bc)
        call gradient(this%decomp,this%der,w, duidxj(:,:,:,7), duidxj(:,:,:,8), duidxj(:,:,:,9),  x_bc,  y_bc, -z_bc)


        if(this%strainHard) then

           if(this%cnsrv_g) then
              if(this%useOneG) then
                 call this%getRHS_rg(rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
              else
                 call this%getRHS_rg(rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc)
              endif
           else
              if(this%useOneG) then
                 call this%getRHS_g(rho,u,v,w,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
              else
                 call this%getRHS_g(rho,u,v,w,dt,src,rhsg,x_bc,y_bc,z_bc)
              endif
           endif

           if(this%use_gTg) then !gtTgt and gpTgp

              if(this%cnsrv_gt) then
                 if(this%useOneG) then
                    call this%getRHS_rgtTgt(rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_rgtTgt(rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc)
                 endif
              else
                 if(this%useOneG) then
                    call this%getRHS_gtTgt(rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_gtTgt(rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc)
                 endif
              endif

              if(this%cnsrv_gp) then
                 if(this%useOneG) then
                    call this%getRHS_rgpTgp(rho,u,v,w,duidxj,dt,src,rhsgp,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_rgpTgp(rho,u,v,w,duidxj,dt,src,rhsgp,x_bc,y_bc,z_bc)
                 endif
              else
                 if(this%useOneG) then
                    call this%getRHS_gpTgp(rho,u,v,w,duidxj,dt,src,rhsgp,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_gpTgp(rho,u,v,w,duidxj,dt,src,rhsgp,x_bc,y_bc,z_bc)
                 endif
              endif

           else

              if(this%cnsrv_gt) then
                 if(this%useOneG) then
                    call this%getRHS_rgt(rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_rgt(rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc)
                 endif
              else
                 if(this%useOneG) then
                    call this%getRHS_gt(rho,u,v,w,dt,src,rhsgt,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_gt(rho,u,v,w,dt,src,rhsgt,x_bc,y_bc,z_bc)
                 endif
              endif

              if(this%cnsrv_gp) then
                 if(this%useOneG) then
                    call this%getRHS_rgp(rho,u,v,w,dt,src,rhsgp,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_rgp(rho,u,v,w,dt,src,rhsgp,x_bc,y_bc,z_bc)
                 endif
              else
                 if(this%useOneG) then
                    call this%getRHS_gp(rho,u,v,w,dt,src,rhsgp,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_gp(rho,u,v,w,dt,src,rhsgp,x_bc,y_bc,z_bc)
                 endif
              endif

           endif

           if(this%cnsrv_pe) then
              if(this%useOneG) then
                 call this%getRHS_rpe(rho,u,v,w,dt,src,rhspe,x_bc,y_bc,z_bc,rho0mix)
              else
                 call this%getRHS_rpe(rho,u,v,w,dt,src,rhspe,x_bc,y_bc,z_bc)
              endif
           else
              if(this%useOneG) then
                 call this%getRHS_pe(rho,u,v,w,dt,src,rhspe,x_bc,y_bc,z_bc,rho0mix)
              else
                 call this%getRHS_pe(rho,u,v,w,dt,src,rhspe,x_bc,y_bc,z_bc)
              endif
           endif

        else

           !if no strain hardening only update g^e

           if(this%use_gTg) then

              if(this%cnsrv_g) then
                 if(this%useOneG) then
                    call this%getRHS_rgTg(rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_rgTg(rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc)
                 endif
              else
                 if(this%useOneG) then
                    call this%getRHS_gTg(rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_gTg(rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc)
                 endif
              endif

           else

              if(this%cnsrv_g) then
                 if(this%useOneG) then
                    call this%getRHS_rg(rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_rg(rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc)
                 endif
              else
                 if(this%useOneG) then
                    call this%getRHS_g(rho,u,v,w,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
                 else
                    call this%getRHS_g(rho,u,v,w,dt,src,rhsg,x_bc,y_bc,z_bc)
                 endif
              endif

           endif

        endif




        if(present(solidVF)) then
           do i = 1, 9
               rhsg(:,:,:,i) = rhsg(:,:,:,i) * solidVF
               rhsgt(:,:,:,i) = rhsgt(:,:,:,i) * solidVF
               rhsgp(:,:,:,i) = rhsgp(:,:,:,i) * solidVF
           enddo
           rhspe = rhspe * solidVF
        endif


        call hook_material_g_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsg)


        ! advance sub-step: g
        if(isub==1) this%Qtmpg = zero                   
        this%Qtmpg  = dt*rhsg + RK45_A(isub)*this%Qtmpg
        this%rg = this%rg  + RK45_B(isub)*this%Qtmpg

        ! advance sub-step: g_t
        if(isub==1) this%Qtmpg_t = zero  
        this%Qtmpg_t  = dt*rhsgt + RK45_A(isub)*this%Qtmpg_t
        this%rg_t = this%rg_t  + RK45_B(isub)*this%Qtmpg_t

        ! advance sub-step: g_p
        if(isub==1) this%Qtmpg_p = zero  
        this%Qtmpg_p  = dt*rhsgp + RK45_A(isub)*this%Qtmpg_p
        this%rg_p = this%rg_p  + RK45_B(isub)*this%Qtmpg_p

        ! advance sub-step: pe
        if(isub==1) this%Qtmppe = zero 
        this%Qtmppe  = dt*rhspe + RK45_A(isub)*this%Qtmppe
        this%rpe = this%rpe  + RK45_B(isub)*this%Qtmppe


        ! ! Now project tensors to SPD space -- mca: do not do this -- g not symmetric for > 1-D
        ! call this%elastic%make_tensor_SPD(this%rg)
        ! call this%elastic%make_tensor_SPD(this%rg_t)
        ! call this%elastic%make_tensor_SPD(this%rg_p)


        ! instead you can use this -- but it is redundant b/c RHS forced symmetric
        ! !enforce symmetry for gTg
        ! if(this%use_gTg) then
        !    if(this%strainHard) then
        !       this%rg_t(:,:,:,4) = this%rg_t(:,:,:,2)
        !       this%rg_t(:,:,:,7) = this%rg_t(:,:,:,3)
        !       this%rg_t(:,:,:,8) = this%rg_t(:,:,:,6)

        !       this%rg_p(:,:,:,4) = this%rg_p(:,:,:,2)
        !       this%rg_p(:,:,:,7) = this%rg_p(:,:,:,3)
        !       this%rg_p(:,:,:,8) = this%rg_p(:,:,:,6)

        !       call this%elastic%make_tensor_SPD(this%rg_t)
        !       call this%elastic%make_tensor_SPD(this%rg_p)

        !    else
        !       this%rg(:,:,:,4) = this%rg(:,:,:,2)
        !       this%rg(:,:,:,7) = this%rg(:,:,:,3)
        !       this%rg(:,:,:,8) = this%rg(:,:,:,6)

        !       call this%elastic%make_tensor_SPD(this%rg)
        !    endif
        ! endif

    end subroutine


    subroutine getRHS_g(this,rho,u,v,w,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        ! Symmetry and anti-symmetry properties of g are assumed as below
        ! In x g_{ij}: [S A A; A S S; A S S]
        ! In y g_{ij}: [S A S; A S A; S A S]
        ! In z g_{ij}: [S S A; S S A; A A S]


        rhsg = zero

        detg  = this%g11*(this%g22*this%g33-this%g23*this%g32) &
              - this%g12*(this%g21*this%g33-this%g31*this%g23) &
              + this%g13*(this%g21*this%g32-this%g31*this%g22)

        ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
        tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
        ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF

        if(present(rho0mix)) then
            penalty  = this%elastic%eta_det_ge*(rho/rho0mix/detg  - one)/dt
            this%det_e = (rho/rho0mix/detg  - one) !diagnostic
        else
            penalty  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
            this%det_e = (tmp /detg /this%elastic%rho0-one) !diagnostic
        endif
        if(this%pRelax) then
           penalty  = this%VF*this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density -- change2
        endif

        where (this%elastic%mu .LT. eps)
            penalty = zero
            this%det_e = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty  = penalty  - src/this%VF
        endif

        !Transport terms
        tmp = -u*this%g11-v*this%g12-w*this%g13
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,1),rhsg(:,:,:,2),rhsg(:,:,:,3),-x_bc, y_bc, z_bc)
        call curl(this%decomp, this%der, this%g11, this%g12, this%g13, curlg, -x_bc, y_bc, z_bc)
        this%curl_e = curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

        !LAD terms
        call gradient(this%decomp,this%der,this%g11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%g11

        !LAD terms
        call gradient(this%decomp,this%der,this%g12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%g12

        !LAD terms
        call gradient(this%decomp,this%der,this%g13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%g13

        !curl diffusion: Miller and Colella -- mca: leads to instability, check implementation
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,2) = rhsg(:,:,:,2) - this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,3) = rhsg(:,:,:,3) - this%elastic%diff_c_ge/dt * diffg(:,:,:,1)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) - this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + this%elastic%diff_c_ge/dt * diffg(:,:,:,1)


        !Transport terms
        tmp = -u*this%g21-v*this%g22-w*this%g23
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,4),rhsg(:,:,:,5),rhsg(:,:,:,6), x_bc,-y_bc, z_bc)
        call curl(this%decomp, this%der, this%g21, this%g22, this%g23, curlg, x_bc, -y_bc, z_bc)
        this%curl_e = this%curl_e + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

        !LAD terms
        call gradient(this%decomp,this%der,this%g21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc,-y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%g21

        !LAD terms
        call gradient(this%decomp,this%der,this%g22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%g22

        !LAD terms
        call gradient(this%decomp,this%der,this%g23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%g23


        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,5) = rhsg(:,:,:,5) - this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,6) = rhsg(:,:,:,6) - this%elastic%diff_c_ge/dt * diffg(:,:,:,1)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) - this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + this%elastic%diff_c_ge/dt * diffg(:,:,:,1)


        !Transport terms
        tmp = -u*this%g31-v*this%g32-w*this%g33
        call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,7),rhsg(:,:,:,8),rhsg(:,:,:,9), x_bc, y_bc,-z_bc)
        call curl(this%decomp, this%der, this%g31, this%g32, this%g33, curlg, x_bc, y_bc, -z_bc)
        this%curl_e = this%curl_e + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2
        this%curl_e = sqrt(this%curl_e)

        !LAD terms
        call gradient(this%decomp,this%der,this%g31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%g31

        !LAD terms
        call gradient(this%decomp,this%der,this%g32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%g32

        !LAD terms
        call gradient(this%decomp,this%der,this%g33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%g33

        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,8) = rhsg(:,:,:,8) - this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,9) = rhsg(:,:,:,9) - this%elastic%diff_c_ge/dt * diffg(:,:,:,1)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) + this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) - this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + this%elastic%diff_c_ge/dt * diffg(:,:,:,1)

        !RHS update for explicit plastic terms
        if (this%plast) then
           if(this%explPlast) then
               call this%getPlasticSources(detg,rhsg)
           end if
        end if

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_cpg_west) then
               !TODO: finish implementation of new sharpening term
               CONTINUE     
           else
           endif
           if(this%intSharp_spf) then
              do i=1,9
                 rhsg(:,:,:,i) = rhsg(:,:,:,i) + this%intSharp_rg(:,:,:,i,1)/rho !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,1,1)/rho,this%intSharp_rgDiff(:,:,:,1,2)/rho,this%intSharp_rgDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,2,1)/rho,this%intSharp_rgDiff(:,:,:,2,2)/rho,this%intSharp_rgDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,3,1)/rho,this%intSharp_rgDiff(:,:,:,3,2)/rho,this%intSharp_rgDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,4,1)/rho,this%intSharp_rgDiff(:,:,:,4,2)/rho,this%intSharp_rgDiff(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,5,1)/rho,this%intSharp_rgDiff(:,:,:,5,2)/rho,this%intSharp_rgDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,6,1)/rho,this%intSharp_rgDiff(:,:,:,6,2)/rho,this%intSharp_rgDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,7,1)/rho,this%intSharp_rgDiff(:,:,:,7,2)/rho,this%intSharp_rgDiff(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,8,1)/rho,this%intSharp_rgDiff(:,:,:,8,2)/rho,this%intSharp_rgDiff(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,9,1)/rho,this%intSharp_rgDiff(:,:,:,9,2)/rho,this%intSharp_rgDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,1,1)/rho,this%intSharp_rg(:,:,:,1,2)/rho,this%intSharp_rg(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,2,1)/rho,this%intSharp_rg(:,:,:,2,2)/rho,this%intSharp_rg(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,3,1)/rho,this%intSharp_rg(:,:,:,3,2)/rho,this%intSharp_rg(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,4,1)/rho,this%intSharp_rg(:,:,:,4,2)/rho,this%intSharp_rg(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,5,1)/rho,this%intSharp_rg(:,:,:,5,2)/rho,this%intSharp_rg(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,6,1)/rho,this%intSharp_rg(:,:,:,6,2)/rho,this%intSharp_rg(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,7,1)/rho,this%intSharp_rg(:,:,:,7,2)/rho,this%intSharp_rg(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,8,1)/rho,this%intSharp_rg(:,:,:,8,2)/rho,this%intSharp_rg(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,9,1)/rho,this%intSharp_rg(:,:,:,9,2)/rho,this%intSharp_rg(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,1,1)/rho,this%intSharp_rgDiff(:,:,:,1,2)/rho,this%intSharp_rgDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,2,1)/rho,this%intSharp_rgDiff(:,:,:,2,2)/rho,this%intSharp_rgDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,3,1)/rho,this%intSharp_rgDiff(:,:,:,3,2)/rho,this%intSharp_rgDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,4,1)/rho,this%intSharp_rgDiff(:,:,:,4,2)/rho,this%intSharp_rgDiff(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,5,1)/rho,this%intSharp_rgDiff(:,:,:,5,2)/rho,this%intSharp_rgDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,6,1)/rho,this%intSharp_rgDiff(:,:,:,6,2)/rho,this%intSharp_rgDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,7,1)/rho,this%intSharp_rgDiff(:,:,:,7,2)/rho,this%intSharp_rgDiff(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,8,1)/rho,this%intSharp_rgDiff(:,:,:,8,2)/rho,this%intSharp_rgDiff(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,9,1)/rho,this%intSharp_rgDiff(:,:,:,9,2)/rho,this%intSharp_rgDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
              !FV terms
              rhsg = rhsg + this%intSharp_gFV
           endif
        endif

    end subroutine


    subroutine getRHS_rg(this,rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg, diff_rg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        ! Symmetry and anti-symmetry properties of g are assumed as below
        ! In x g_{ij}: [S A A; A S S; A S S]
        ! In y g_{ij}: [S A S; A S A; S A S]
        ! In z g_{ij}: [S S A; S S A; A A S]

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        ! call gradient(this%decomp, this%der, u, dudx, dudy, dudz, -x_bc,  y_bc,  z_bc)
        ! call gradient(this%decomp, this%der, v, dvdx, dvdy, dvdz,  x_bc, -y_bc,  z_bc)
        ! call gradient(this%decomp, this%der, w, dwdx, dwdy, dwdz,  x_bc,  y_bc, -z_bc)



        rhsg = zero
        diff_rg = rho*this%diff_g

        detg  = this%g11*(this%g22*this%g33-this%g23*this%g32) &
              - this%g12*(this%g21*this%g33-this%g31*this%g23) &
              + this%g13*(this%g21*this%g32-this%g31*this%g22)

        ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
        tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
        ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF

        if(present(rho0mix)) then
            penalty  = this%elastic%eta_det_ge*(rho/rho0mix/detg  - one)/dt
            this%det_e = (rho/rho0mix/detg  - one) !diagnostic
        else
            penalty  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
            this%det_e = (tmp /detg /this%elastic%rho0-one) !diagnostic
        endif
        if(this%pRelax) then
           penalty  = this%VF*this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density -- change2
        endif

        where (this%elastic%mu .LT. eps)
            penalty = zero
            this%det_e = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty  = penalty  - src/this%VF
        endif


        !Curl for diagnostic
        call curl(this%decomp, this%der, this%g11, this%g12, this%g13, curlg, -x_bc, y_bc, z_bc)
        this%curl_e = curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg11,-v*this%rg11,-w*this%rg11,rhsg(:,:,:,1),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) - this%rg11*dudx - this%rg12*dvdx - this%rg13*dwdx
        !LAD terms
        call gradient(this%decomp,this%der,this%g11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp + penalty*this%rg11

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg12,-v*this%rg12,-w*this%rg12,rhsg(:,:,:,2),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,2) = rhsg(:,:,:,2) - this%rg11*dudy - this%rg12*dvdy - this%rg13*dwdy
        !LAD terms
        call gradient(this%decomp,this%der,this%g12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp + penalty*this%rg12

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg13,-v*this%rg13,-w*this%rg13,rhsg(:,:,:,3),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,3) = rhsg(:,:,:,3) - this%rg11*dudz - this%rg12*dvdz - this%rg13*dwdz
        !LAD terms
        call gradient(this%decomp,this%der,this%g13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp + penalty*this%rg13

        !curl diffusion: Miller and Colella -- mca: leads to instability, check implementation
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,2) = rhsg(:,:,:,2) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,3)
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,3) = rhsg(:,:,:,3) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,1)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,2)
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,1)


        !Curl for diagnostic
        call curl(this%decomp, this%der, this%g21, this%g22, this%g23, curlg, x_bc, -y_bc, z_bc)
        this%curl_e = this%curl_e + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg21,-v*this%rg21,-w*this%rg21,rhsg(:,:,:,4),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) - this%rg21*dudx - this%rg22*dvdx - this%rg23*dwdx
        !LAD terms
        call gradient(this%decomp,this%der,this%g21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp + penalty*this%rg21

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg22,-v*this%rg22,-w*this%rg22,rhsg(:,:,:,5),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,5) = rhsg(:,:,:,5) - this%rg21*dudy - this%rg22*dvdy - this%rg23*dwdy
        !LAD terms
        call gradient(this%decomp,this%der,this%g22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp + penalty*this%rg22

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg23,-v*this%rg23,-w*this%rg23,rhsg(:,:,:,6),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,6) = rhsg(:,:,:,6) - this%rg21*dudz - this%rg22*dvdz - this%rg23*dwdz
        !LAD terms
        call gradient(this%decomp,this%der,this%g23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp + penalty*this%rg23


        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,5) = rhsg(:,:,:,5) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,3)
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,6) = rhsg(:,:,:,6) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,1)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,4) = rhsg(:,:,:,4) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,2)
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,1)


        !Curl for diagnostic
        call curl(this%decomp, this%der, this%g31, this%g32, this%g33, curlg, x_bc, y_bc, -z_bc)
        this%curl_e = this%curl_e + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2
        this%curl_e = sqrt(this%curl_e)

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg31,-v*this%rg31,-w*this%rg31,rhsg(:,:,:,7),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) - this%rg31*dudx - this%rg32*dvdx - this%rg33*dwdx
        !LAD terms
        call gradient(this%decomp,this%der,this%g31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp + penalty*this%rg31

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg32,-v*this%rg32,-w*this%rg32,rhsg(:,:,:,8),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,8) = rhsg(:,:,:,8) - this%rg31*dudy - this%rg32*dvdy - this%rg33*dwdy
        !LAD terms
        call gradient(this%decomp,this%der,this%g32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp + penalty*this%rg32

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg33,-v*this%rg33,-w*this%rg33,rhsg(:,:,:,9),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,9) = rhsg(:,:,:,9) - this%rg31*dudz - this%rg32*dvdz - this%rg33*dwdz
        !LAD terms
        call gradient(this%decomp,this%der,this%g33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp + penalty*this%rg33

        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,8) = rhsg(:,:,:,8) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,3)
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,9) = rhsg(:,:,:,9) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,1)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsg(:,:,:,7) = rhsg(:,:,:,7) - this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,2)
        rhsg(:,:,:,8) = rhsg(:,:,:,8) + this%elastic%diff_c_ge/dt * rho*diffg(:,:,:,1)

        !RHS update for explicit plastic terms
        if (this%plast) then
           if(this%explPlast) then
               call this%getPlasticSources(detg,rhsg)
           end if
        end if

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_spf) then
              rhsg = rhsg + this%intSharp_rg(:,:,:,:,1) !ignore components 2 and 3 when not in divergence form
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,1,1),this%intSharp_rgDiff(:,:,:,1,2),this%intSharp_rgDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,2,1),this%intSharp_rgDiff(:,:,:,2,2),this%intSharp_rgDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,3,1),this%intSharp_rgDiff(:,:,:,3,2),this%intSharp_rgDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,4,1),this%intSharp_rgDiff(:,:,:,4,2),this%intSharp_rgDiff(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,5,1),this%intSharp_rgDiff(:,:,:,5,2),this%intSharp_rgDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,6,1),this%intSharp_rgDiff(:,:,:,6,2),this%intSharp_rgDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,7,1),this%intSharp_rgDiff(:,:,:,7,2),this%intSharp_rgDiff(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,8,1),this%intSharp_rgDiff(:,:,:,8,2),this%intSharp_rgDiff(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,9,1),this%intSharp_rgDiff(:,:,:,9,2),this%intSharp_rgDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,1,1),this%intSharp_rg(:,:,:,1,2),this%intSharp_rg(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,2,1),this%intSharp_rg(:,:,:,2,2),this%intSharp_rg(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,3,1),this%intSharp_rg(:,:,:,3,2),this%intSharp_rg(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,4,1),this%intSharp_rg(:,:,:,4,2),this%intSharp_rg(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,5,1),this%intSharp_rg(:,:,:,5,2),this%intSharp_rg(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,6,1),this%intSharp_rg(:,:,:,6,2),this%intSharp_rg(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,7,1),this%intSharp_rg(:,:,:,7,2),this%intSharp_rg(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,8,1),this%intSharp_rg(:,:,:,8,2),this%intSharp_rg(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,9,1),this%intSharp_rg(:,:,:,9,2),this%intSharp_rg(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,1,1),this%intSharp_rgDiff(:,:,:,1,2),this%intSharp_rgDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,2,1),this%intSharp_rgDiff(:,:,:,2,2),this%intSharp_rgDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,3,1),this%intSharp_rgDiff(:,:,:,3,2),this%intSharp_rgDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,4,1),this%intSharp_rgDiff(:,:,:,4,2),this%intSharp_rgDiff(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,5,1),this%intSharp_rgDiff(:,:,:,5,2),this%intSharp_rgDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,6,1),this%intSharp_rgDiff(:,:,:,6,2),this%intSharp_rgDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,7,1),this%intSharp_rgDiff(:,:,:,7,2),this%intSharp_rgDiff(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,8,1),this%intSharp_rgDiff(:,:,:,8,2),this%intSharp_rgDiff(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,9,1),this%intSharp_rgDiff(:,:,:,9,2),this%intSharp_rgDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
              !FV terms
              rhsg = rhsg + this%intSharp_rgFV
           endif
        endif

    end subroutine


    subroutine getRHS_gt(this,rho,u,v,w,dt,src,rhsgt,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsgt
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detgt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        ! Symmetry and anti-symmetry properties of g are assumed as below
        ! In x g_{ij}: [S A A; A S S; A S S]
        ! In y g_{ij}: [S A S; A S A; S A S]
        ! In z g_{ij}: [S S A; S S A; A A S]

        rhsgt = zero

        !g_t rhs
        detgt = this%gt11*(this%gt22*this%gt33-this%gt23*this%gt32) &
             - this%gt12*(this%gt21*this%gt33-this%gt31*this%gt23) &
             + this%gt13*(this%gt21*this%gt32-this%gt31*this%gt22)

        tmp = (rho*this%Ys + this%elastic%rho0*detgt*epssmall)/(this%VF + epssmall)   

        if(present(rho0mix)) then
           penalty = this%elastic%eta_det_gt*(rho/rho0mix/detgt - one)/dt
           this%det_t = (rho/rho0mix/detgt - one) !diagnostic
        else
           penalty = this%elastic%eta_det_gt*( tmp/detgt/this%elastic%rho0-one)/dt
           this%det_t = (tmp/detgt/this%elastic%rho0 - one)
        endif
        if(this%pRelax) then
           penalty = this%VF*this%elastic%eta_det_gt*( tmp/detgt/this%elastic%rho0-one)/dt
        endif

        where (this%elastic%mu .LT. eps)
           penalty = zero
           this%det_t = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
           ! add Fsource term to penalty 
           penalty = penalty - src/this%VF
        endif


        !Transport terms
        tmp = -u*this%gt11-v*this%gt12-w*this%gt13
        call gradient(this%decomp,this%der,tmp,rhsgt(:,:,:,1),rhsgt(:,:,:,2),rhsgt(:,:,:,3),-x_bc, y_bc, z_bc)
        call curl(this%decomp, this%der, this%gt11, this%gt12, this%gt13, curlg, -x_bc, y_bc, z_bc)
        this%curl_t = curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

        !LAD terms
        call gradient(this%decomp,this%der,this%gt11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%gt11

        !LAD terms
        call gradient(this%decomp,this%der,this%gt12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%gt12 

        !LAD terms
        call gradient(this%decomp,this%der,this%gt13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%gt13 

        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) - this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) - this%elastic%diff_c_gt/dt * diffg(:,:,:,1)
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) - this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + this%elastic%diff_c_gt/dt * diffg(:,:,:,1)


        !Transport terms
        tmp = -u*this%gt21-v*this%gt22-w*this%gt23
        call gradient(this%decomp,this%der,tmp,rhsgt(:,:,:,4),rhsgt(:,:,:,5),rhsgt(:,:,:,6), x_bc,-y_bc, z_bc)   
        call curl(this%decomp, this%der, this%gt21, this%gt22, this%gt23, curlg, x_bc, -y_bc, z_bc)
        this%curl_t = this%curl_t + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

        !LAD terms
        call gradient(this%decomp,this%der,this%gt21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), -x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, x_bc, y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%gt21 

        !LAD terms
        call gradient(this%decomp,this%der,this%gt22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, -x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%gt22 

        !LAD terms
        call gradient(this%decomp,this%der,this%gt23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, -x_bc, y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%gt23 

        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) - this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) - this%elastic%diff_c_gt/dt * diffg(:,:,:,1)
        rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,4) = rhsgt(:,:,:,4) - this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + this%elastic%diff_c_gt/dt * diffg(:,:,:,1)


        !Transport terms
        tmp = -u*this%gt31-v*this%gt32-w*this%gt33
        call gradient(this%decomp,this%der,tmp,rhsgt(:,:,:,7),rhsgt(:,:,:,8),rhsgt(:,:,:,9), x_bc, y_bc,-z_bc)
        call curl(this%decomp, this%der, this%gt31, this%gt32, this%gt33, curlg, x_bc, y_bc, -z_bc)
        this%curl_t = this%curl_t + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2
        this%curl_t = sqrt(this%curl_t)

        !LAD terms
        call gradient(this%decomp,this%der,this%gt31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), -x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, x_bc, -y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%gt31 

        !LAD terms
        call gradient(this%decomp,this%der,this%gt32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, -x_bc, y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%gt32 

        !LAD terms
        call gradient(this%decomp,this%der,this%gt33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, -x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%gt33 

        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,8) = rhsgt(:,:,:,8) - this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) - this%elastic%diff_c_gt/dt * diffg(:,:,:,1)
        rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,7) = rhsgt(:,:,:,7) - this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
        rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + this%elastic%diff_c_gt/dt * diffg(:,:,:,1)


        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_cpg_west) then
               !TODO: finish implementation of new sharpening term
               CONTINUE     
           else
           endif
           if(this%intSharp_spf) then
              do i=1,9
                 rhsgt(:,:,:,i) = rhsgt(:,:,:,i) + this%intSharp_rgt(:,:,:,i,1)/rho !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,1,1)/rho,this%intSharp_rgtDiff(:,:,:,1,2)/rho,this%intSharp_rgtDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,2,1)/rho,this%intSharp_rgtDiff(:,:,:,2,2)/rho,this%intSharp_rgtDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,3,1)/rho,this%intSharp_rgtDiff(:,:,:,3,2)/rho,this%intSharp_rgtDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,4,1)/rho,this%intSharp_rgtDiff(:,:,:,4,2)/rho,this%intSharp_rgtDiff(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,5,1)/rho,this%intSharp_rgtDiff(:,:,:,5,2)/rho,this%intSharp_rgtDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,6,1)/rho,this%intSharp_rgtDiff(:,:,:,6,2)/rho,this%intSharp_rgtDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,7,1)/rho,this%intSharp_rgtDiff(:,:,:,7,2)/rho,this%intSharp_rgtDiff(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,8,1)/rho,this%intSharp_rgtDiff(:,:,:,8,2)/rho,this%intSharp_rgtDiff(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,9,1)/rho,this%intSharp_rgtDiff(:,:,:,9,2)/rho,this%intSharp_rgtDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,1,1)/rho,this%intSharp_rgt(:,:,:,1,2)/rho,this%intSharp_rgt(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,2,1)/rho,this%intSharp_rgt(:,:,:,2,2)/rho,this%intSharp_rgt(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,3,1)/rho,this%intSharp_rgt(:,:,:,3,2)/rho,this%intSharp_rgt(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,4,1)/rho,this%intSharp_rgt(:,:,:,4,2)/rho,this%intSharp_rgt(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,5,1)/rho,this%intSharp_rgt(:,:,:,5,2)/rho,this%intSharp_rgt(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,6,1)/rho,this%intSharp_rgt(:,:,:,6,2)/rho,this%intSharp_rgt(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,7,1)/rho,this%intSharp_rgt(:,:,:,7,2)/rho,this%intSharp_rgt(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,8,1)/rho,this%intSharp_rgt(:,:,:,8,2)/rho,this%intSharp_rgt(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,9,1)/rho,this%intSharp_rgt(:,:,:,9,2)/rho,this%intSharp_rgt(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,1,1)/rho,this%intSharp_rgtDiff(:,:,:,1,2)/rho,this%intSharp_rgtDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,2,1)/rho,this%intSharp_rgtDiff(:,:,:,2,2)/rho,this%intSharp_rgtDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,3,1)/rho,this%intSharp_rgtDiff(:,:,:,3,2)/rho,this%intSharp_rgtDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,4,1)/rho,this%intSharp_rgtDiff(:,:,:,4,2)/rho,this%intSharp_rgtDiff(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,5,1)/rho,this%intSharp_rgtDiff(:,:,:,5,2)/rho,this%intSharp_rgtDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,6,1)/rho,this%intSharp_rgtDiff(:,:,:,6,2)/rho,this%intSharp_rgtDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,7,1)/rho,this%intSharp_rgtDiff(:,:,:,7,2)/rho,this%intSharp_rgtDiff(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,8,1)/rho,this%intSharp_rgtDiff(:,:,:,8,2)/rho,this%intSharp_rgtDiff(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,9,1)/rho,this%intSharp_rgtDiff(:,:,:,9,2)/rho,this%intSharp_rgtDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
              !FV terms
              rhsgt = rhsgt + this%intSharp_gtFV
           endif
        endif
           
    end subroutine


    subroutine getRHS_rgt(this,rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsgt
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detgt, diff_rgt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        ! Symmetry and anti-symmetry properties of g are assumed as below
        ! In x g_{ij}: [S A A; A S S; A S S]
        ! In y g_{ij}: [S A S; A S A; S A S]
        ! In z g_{ij}: [S S A; S S A; A A S]

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        ! call gradient(this%decomp, this%der, u, dudx, dudy, dudz, -x_bc,  y_bc,  z_bc)
        ! call gradient(this%decomp, this%der, v, dvdx, dvdy, dvdz,  x_bc, -y_bc,  z_bc)
        ! call gradient(this%decomp, this%der, w, dwdx, dwdy, dwdz,  x_bc,  y_bc, -z_bc)



        rhsgt = zero
        diff_rgt = rho*this%diff_gt

        !g_t rhs
        detgt = this%gt11*(this%gt22*this%gt33-this%gt23*this%gt32) &
             - this%gt12*(this%gt21*this%gt33-this%gt31*this%gt23) &
             + this%gt13*(this%gt21*this%gt32-this%gt31*this%gt22)

        tmp = (rho*this%Ys + this%elastic%rho0*detgt*epssmall)/(this%VF + epssmall)   

        if(present(rho0mix)) then
           penalty = this%elastic%eta_det_gt*(rho/rho0mix/detgt - one)/dt
           this%det_t = (rho/rho0mix/detgt - one) !diagnostic
        else
           penalty = this%elastic%eta_det_gt*( tmp/detgt/this%elastic%rho0-one)/dt
           this%det_t = (tmp/detgt/this%elastic%rho0 - one)
        endif
        if(this%pRelax) then
           penalty = this%VF*this%elastic%eta_det_gt*( tmp/detgt/this%elastic%rho0-one)/dt
        endif

        where (this%elastic%mu .LT. eps)
           penalty = zero
           this%det_t = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
           ! add Fsource term to penalty 
           penalty = penalty - src/this%VF
        endif



        !Curl for diagnostic
        call curl(this%decomp, this%der, this%gt11, this%gt12, this%gt13, curlg, -x_bc, y_bc, z_bc)
        this%curl_t = curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt11,-v*this%rgt11,-w*this%rgt11,rhsgt(:,:,:,1),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) - this%rgt11*dudx - this%rgt12*dvdx - this%rgt13*dwdx
        !LAD terms
        call gradient(this%decomp,this%der,this%gt11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp + penalty*this%rgt11

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt12,-v*this%rgt12,-w*this%rgt12,rhsgt(:,:,:,2),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) - this%rgt11*dudy - this%rgt12*dvdy - this%rgt13*dwdy
        !LAD terms
        call gradient(this%decomp,this%der,this%gt12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp + penalty*this%rgt12

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt13,-v*this%rgt13,-w*this%rgt13,rhsgt(:,:,:,3),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) - this%rgt11*dudz - this%rgt12*dvdz - this%rgt13*dwdz
        !LAD terms
        call gradient(this%decomp,this%der,this%gt13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp + penalty*this%rgt13

        !curl diffusion: Miller and Colella -- mca: leads to instability, check implementation
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,3)
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,1)
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,2)
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,1)


        !Curl for diagnostic
        call curl(this%decomp, this%der, this%gt21, this%gt22, this%gt23, curlg, x_bc, -y_bc, z_bc)
        this%curl_t = this%curl_t + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt21,-v*this%rgt21,-w*this%rgt21,rhsgt(:,:,:,4),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,4) = rhsgt(:,:,:,4) - this%rgt21*dudx - this%rgt22*dvdx - this%rgt23*dwdx
        !LAD terms
        call gradient(this%decomp,this%der,this%gt21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp + penalty*this%rgt21

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt22,-v*this%rgt22,-w*this%rgt22,rhsgt(:,:,:,5),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) - this%rgt21*dudy - this%rgt22*dvdy - this%rgt23*dwdy
        !LAD terms
        call gradient(this%decomp,this%der,this%gt22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp + penalty*this%rgt22

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt23,-v*this%rgt23,-w*this%rgt23,rhsgt(:,:,:,6),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) - this%rgt21*dudz - this%rgt22*dvdz - this%rgt23*dwdz
        !LAD terms
        call gradient(this%decomp,this%der,this%gt23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp + penalty*this%rgt23


        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,3)
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,1)
        rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,4) = rhsgt(:,:,:,4) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,2)
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,1)


        !Curl for diagnostic
        call curl(this%decomp, this%der, this%gt31, this%gt32, this%gt33, curlg, x_bc, y_bc, -z_bc)
        this%curl_t = this%curl_t + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2
        this%curl_t = sqrt(this%curl_t)

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt31,-v*this%rgt31,-w*this%rgt31,rhsgt(:,:,:,7),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,7) = rhsgt(:,:,:,7) - this%rgt31*dudx - this%rgt32*dvdx - this%rgt33*dwdx
        !LAD terms
        call gradient(this%decomp,this%der,this%gt31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp + penalty*this%rgt31

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt32,-v*this%rgt32,-w*this%rgt32,rhsgt(:,:,:,8),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,8) = rhsgt(:,:,:,8) - this%rgt31*dudy - this%rgt32*dvdy - this%rgt33*dwdy
        !LAD terms
        call gradient(this%decomp,this%der,this%gt32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp + penalty*this%rgt32

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt33,-v*this%rgt33,-w*this%rgt33,rhsgt(:,:,:,9),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) - this%rgt31*dudz - this%rgt32*dvdz - this%rgt33*dwdz
        !LAD terms
        call gradient(this%decomp,this%der,this%gt33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp + penalty*this%rgt33

        !curl diffusion: Miller and Colella
        call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,8) = rhsgt(:,:,:,8) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,3)
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,2)
        call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,1)
        rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,3)
        call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgt(:,:,:,7) = rhsgt(:,:,:,7) - this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,2)
        rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + this%elastic%diff_c_gt/dt * rho*diffg(:,:,:,1)

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_spf) then
              rhsgt = rhsgt + this%intSharp_rgt(:,:,:,:,1) !ignore components 2 and 3 when not in divergence form
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,1,1),this%intSharp_rgtDiff(:,:,:,1,2),this%intSharp_rgtDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,2,1),this%intSharp_rgtDiff(:,:,:,2,2),this%intSharp_rgtDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,3,1),this%intSharp_rgtDiff(:,:,:,3,2),this%intSharp_rgtDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,4,1),this%intSharp_rgtDiff(:,:,:,4,2),this%intSharp_rgtDiff(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,5,1),this%intSharp_rgtDiff(:,:,:,5,2),this%intSharp_rgtDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,6,1),this%intSharp_rgtDiff(:,:,:,6,2),this%intSharp_rgtDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,7,1),this%intSharp_rgtDiff(:,:,:,7,2),this%intSharp_rgtDiff(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,8,1),this%intSharp_rgtDiff(:,:,:,8,2),this%intSharp_rgtDiff(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,9,1),this%intSharp_rgtDiff(:,:,:,9,2),this%intSharp_rgtDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,1,1),this%intSharp_rgt(:,:,:,1,2),this%intSharp_rgt(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,2,1),this%intSharp_rgt(:,:,:,2,2),this%intSharp_rgt(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,3,1),this%intSharp_rgt(:,:,:,3,2),this%intSharp_rgt(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,4,1),this%intSharp_rgt(:,:,:,4,2),this%intSharp_rgt(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,5,1),this%intSharp_rgt(:,:,:,5,2),this%intSharp_rgt(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,6,1),this%intSharp_rgt(:,:,:,6,2),this%intSharp_rgt(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,7,1),this%intSharp_rgt(:,:,:,7,2),this%intSharp_rgt(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,8,1),this%intSharp_rgt(:,:,:,8,2),this%intSharp_rgt(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,9,1),this%intSharp_rgt(:,:,:,9,2),this%intSharp_rgt(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,1,1),this%intSharp_rgtDiff(:,:,:,1,2),this%intSharp_rgtDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,2,1),this%intSharp_rgtDiff(:,:,:,2,2),this%intSharp_rgtDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,3,1),this%intSharp_rgtDiff(:,:,:,3,2),this%intSharp_rgtDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,4,1),this%intSharp_rgtDiff(:,:,:,4,2),this%intSharp_rgtDiff(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,5,1),this%intSharp_rgtDiff(:,:,:,5,2),this%intSharp_rgtDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,6,1),this%intSharp_rgtDiff(:,:,:,6,2),this%intSharp_rgtDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,7,1),this%intSharp_rgtDiff(:,:,:,7,2),this%intSharp_rgtDiff(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,8,1),this%intSharp_rgtDiff(:,:,:,8,2),this%intSharp_rgtDiff(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,9,1),this%intSharp_rgtDiff(:,:,:,9,2),this%intSharp_rgtDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
              !FV terms
              rhsgt = rhsgt + this%intSharp_rgtFV
           endif
        endif
           
    end subroutine


    subroutine getRHS_gp(this,rho,u,v,w,dt,src,rhsgp,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsgp
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detgp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        ! Symmetry and anti-symmetry properties of g are assumed as below
        ! In x g_{ij}: [S A A; A S S; A S S]
        ! In y g_{ij}: [S A S; A S A; S A S]
        ! In z g_{ij}: [S S A; S S A; A A S]

        rhsgp = zero

        !g_p rhs
        detgp = this%gp11*(this%gp22*this%gp33-this%gp23*this%gp32) &
             - this%gp12*(this%gp21*this%gp33-this%gp31*this%gp23) &
             + this%gp13*(this%gp21*this%gp32-this%gp31*this%gp22)

        ! !use this to set g^p correction in terms of g^e and g^t corrections
        ! tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
        ! if(present(rho0mix)) then
        !    penalty_e  = this%elastic%eta_det_ge*(rho/rho0mix/detg)
        ! else
        !    penalty_e  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0)
        ! endif
        ! penalty = (penalty_e*(this%elastic%eta_det_gp/detgp - one) + this%elastic%eta_det_ge*(1.0-this%elastic%eta_det_gp))/dt

        !independent g^p correction
        penalty = this%elastic%eta_det_gp*(one/detgp - one)/dt

        this%det_p = one/detgp - one

        where (this%elastic%mu .LT. eps)
           penalty = zero
           this%det_p = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
           ! add Fsource term to penalty 
           penalty = penalty - src/this%VF !mca check???
        endif

        this%curl_p = 0 !free variable to track compatability condition -- not set up
             

        !Transport terms
        call gradient(this%decomp,this%der,this%gp11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsgp(:,:,:,1) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp + penalty*this%gp11 

        !Transport terms
        call gradient(this%decomp,this%der,this%gp12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        rhsgp(:,:,:,2) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp + penalty*this%gp12

        !Transport terms
        call gradient(this%decomp,this%der,this%gp13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        rhsgp(:,:,:,3) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp + penalty*this%gp13 

        !Transport terms
        call gradient(this%decomp,this%der,this%gp21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), -x_bc, -y_bc, z_bc)
        rhsgp(:,:,:,4) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, x_bc, y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp + penalty*this%gp21

        !Transport terms
        call gradient(this%decomp,this%der,this%gp22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
        rhsgp(:,:,:,5) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, -x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp + penalty*this%gp22 

        !Transport terms
        call gradient(this%decomp,this%der,this%gp23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, -y_bc, -z_bc)
        rhsgp(:,:,:,6) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, -x_bc, y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp + penalty*this%gp23 

        !Transport terms
        call gradient(this%decomp,this%der,this%gp31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), -x_bc, y_bc, -z_bc)
        rhsgp(:,:,:,7) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, x_bc, -y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp + penalty*this%gp31 

        !Transport terms
        call gradient(this%decomp,this%der,this%gp32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, -y_bc, -z_bc)
        rhsgp(:,:,:,8) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, -x_bc, y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp + penalty*this%gp32 

        !Transport terms
        call gradient(this%decomp,this%der,this%gp33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
        !LAD terms
        rhsgp(:,:,:,9) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, -x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp + penalty*this%gp33



        !RHS update for explicit plastic terms
        if (this%plast) then
           if(this%explPlast) then
               print*,"explict plastic terms not yet implemented for g^p"
               stop
           end if
        end if

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_cpg_west) then
               !TODO: finish implementation of new sharpening term
               CONTINUE     
           else
           endif
           if(this%intSharp_spf) then
              do i=1,9
                 rhsgp(:,:,:,i) = rhsgp(:,:,:,i) + this%intSharp_rgp(:,:,:,i,1)/rho !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,1,1)/rho,this%intSharp_rgpDiff(:,:,:,1,2)/rho,this%intSharp_rgpDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,2,1)/rho,this%intSharp_rgpDiff(:,:,:,2,2)/rho,this%intSharp_rgpDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,3,1)/rho,this%intSharp_rgpDiff(:,:,:,3,2)/rho,this%intSharp_rgpDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,4,1)/rho,this%intSharp_rgpDiff(:,:,:,4,2)/rho,this%intSharp_rgpDiff(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,5,1)/rho,this%intSharp_rgpDiff(:,:,:,5,2)/rho,this%intSharp_rgpDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,6,1)/rho,this%intSharp_rgpDiff(:,:,:,6,2)/rho,this%intSharp_rgpDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,7,1)/rho,this%intSharp_rgpDiff(:,:,:,7,2)/rho,this%intSharp_rgpDiff(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,8,1)/rho,this%intSharp_rgpDiff(:,:,:,8,2)/rho,this%intSharp_rgpDiff(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,9,1)/rho,this%intSharp_rgpDiff(:,:,:,9,2)/rho,this%intSharp_rgpDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,1,1)/rho,this%intSharp_rgp(:,:,:,1,2)/rho,this%intSharp_rgp(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,2,1)/rho,this%intSharp_rgp(:,:,:,2,2)/rho,this%intSharp_rgp(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,3,1)/rho,this%intSharp_rgp(:,:,:,3,2)/rho,this%intSharp_rgp(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,4,1)/rho,this%intSharp_rgp(:,:,:,4,2)/rho,this%intSharp_rgp(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,5,1)/rho,this%intSharp_rgp(:,:,:,5,2)/rho,this%intSharp_rgp(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,6,1)/rho,this%intSharp_rgp(:,:,:,6,2)/rho,this%intSharp_rgp(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,7,1)/rho,this%intSharp_rgp(:,:,:,7,2)/rho,this%intSharp_rgp(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,8,1)/rho,this%intSharp_rgp(:,:,:,8,2)/rho,this%intSharp_rgp(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,9,1)/rho,this%intSharp_rgp(:,:,:,9,2)/rho,this%intSharp_rgp(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,1,1)/rho,this%intSharp_rgpDiff(:,:,:,1,2)/rho,this%intSharp_rgpDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,2,1)/rho,this%intSharp_rgpDiff(:,:,:,2,2)/rho,this%intSharp_rgpDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,3,1)/rho,this%intSharp_rgpDiff(:,:,:,3,2)/rho,this%intSharp_rgpDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,4,1)/rho,this%intSharp_rgpDiff(:,:,:,4,2)/rho,this%intSharp_rgpDiff(:,:,:,4,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,5,1)/rho,this%intSharp_rgpDiff(:,:,:,5,2)/rho,this%intSharp_rgpDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,6,1)/rho,this%intSharp_rgpDiff(:,:,:,6,2)/rho,this%intSharp_rgpDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,7,1)/rho,this%intSharp_rgpDiff(:,:,:,7,2)/rho,this%intSharp_rgpDiff(:,:,:,7,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,8,1)/rho,this%intSharp_rgpDiff(:,:,:,8,2)/rho,this%intSharp_rgpDiff(:,:,:,8,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,9,1)/rho,this%intSharp_rgpDiff(:,:,:,9,2)/rho,this%intSharp_rgpDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
              !FV terms
              rhsgp = rhsgp + this%intSharp_gpFV
           endif
        endif

    end subroutine


    subroutine getRHS_rgp(this,rho,u,v,w,dt,src,rhsgp,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsgp
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detgp, diff_rgp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        ! Symmetry and anti-symmetry properties of g are assumed as below
        ! In x g_{ij}: [S A A; A S S; A S S]
        ! In y g_{ij}: [S A S; A S A; S A S]
        ! In z g_{ij}: [S S A; S S A; A A S]

        rhsgp = zero
        diff_rgp = rho*this%diff_gp

        !g_p rhs
        detgp = this%gp11*(this%gp22*this%gp33-this%gp23*this%gp32) &
             - this%gp12*(this%gp21*this%gp33-this%gp31*this%gp23) &
             + this%gp13*(this%gp21*this%gp32-this%gp31*this%gp22)

        ! !use this to set g^p correction in terms of g^e and g^t corrections
        ! tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
        ! if(present(rho0mix)) then
        !    penalty_e  = this%elastic%eta_det_ge*(rho/rho0mix/detg)
        ! else
        !    penalty_e  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0)
        ! endif
        ! penalty = (penalty_e*(this%elastic%eta_det_gp/detgp - one) + this%elastic%eta_det_ge*(1.0-this%elastic%eta_det_gp))/dt

        !independent g^p correction
        penalty = this%elastic%eta_det_gp*(one/detgp - one)/dt

        this%det_p = one/detgp - one

        where (this%elastic%mu .LT. eps)
           penalty = zero
           this%det_p = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
           ! add Fsource term to penalty 
           penalty = penalty - src/this%VF !mca check???
        endif

        this%curl_p = 0 !free variable to track compatability condition -- not set up
             

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp11,-v*this%rgp11,-w*this%rgp11,rhsgp(:,:,:,1),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp + penalty*this%rgp11 

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp12,-v*this%rgp12,-w*this%rgp12,rhsgp(:,:,:,2),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp + penalty*this%rgp12

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp13,-v*this%rgp13,-w*this%rgp13,rhsgp(:,:,:,3),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp + penalty*this%rgp13 

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp21,-v*this%rgp21,-w*this%rgp21,rhsgp(:,:,:,4),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), -x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp + penalty*this%rgp21

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp22,-v*this%rgp22,-w*this%rgp22,rhsgp(:,:,:,5),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp + penalty*this%rgp22 

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp23,-v*this%rgp23,-w*this%rgp23,rhsgp(:,:,:,6),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp + penalty*this%rgp23 

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp31,-v*this%rgp31,-w*this%rgp31,rhsgp(:,:,:,7),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), -x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp + penalty*this%rgp31 

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp32,-v*this%rgp32,-w*this%rgp32,rhsgp(:,:,:,8),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp + penalty*this%rgp32 

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp33,-v*this%rgp33,-w*this%rgp33,rhsgp(:,:,:,9),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp + penalty*this%rgp33



        !RHS update for explicit plastic terms
        if (this%plast) then
           if(this%explPlast) then
               print*,"explict plastic terms not yet implemented for g^p"
               stop
           end if
        end if

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_spf) then
              rhsgp = rhsgp + this%intSharp_rgp(:,:,:,:,1) !ignore components 2 and 3 when not in divergence form
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,1,1),this%intSharp_rgpDiff(:,:,:,1,2),this%intSharp_rgpDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,2,1),this%intSharp_rgpDiff(:,:,:,2,2),this%intSharp_rgpDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,3,1),this%intSharp_rgpDiff(:,:,:,3,2),this%intSharp_rgpDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,4,1),this%intSharp_rgpDiff(:,:,:,4,2),this%intSharp_rgpDiff(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,5,1),this%intSharp_rgpDiff(:,:,:,5,2),this%intSharp_rgpDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,6,1),this%intSharp_rgpDiff(:,:,:,6,2),this%intSharp_rgpDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,7,1),this%intSharp_rgpDiff(:,:,:,7,2),this%intSharp_rgpDiff(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,8,1),this%intSharp_rgpDiff(:,:,:,8,2),this%intSharp_rgpDiff(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,9,1),this%intSharp_rgpDiff(:,:,:,9,2),this%intSharp_rgpDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,1,1),this%intSharp_rgp(:,:,:,1,2),this%intSharp_rgp(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,2,1),this%intSharp_rgp(:,:,:,2,2),this%intSharp_rgp(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,3,1),this%intSharp_rgp(:,:,:,3,2),this%intSharp_rgp(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,4,1),this%intSharp_rgp(:,:,:,4,2),this%intSharp_rgp(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,5,1),this%intSharp_rgp(:,:,:,5,2),this%intSharp_rgp(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,6,1),this%intSharp_rgp(:,:,:,6,2),this%intSharp_rgp(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,7,1),this%intSharp_rgp(:,:,:,7,2),this%intSharp_rgp(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,8,1),this%intSharp_rgp(:,:,:,8,2),this%intSharp_rgp(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,9,1),this%intSharp_rgp(:,:,:,9,2),this%intSharp_rgp(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,1,1),this%intSharp_rgpDiff(:,:,:,1,2),this%intSharp_rgpDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,2,1),this%intSharp_rgpDiff(:,:,:,2,2),this%intSharp_rgpDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,3,1),this%intSharp_rgpDiff(:,:,:,3,2),this%intSharp_rgpDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,4,1),this%intSharp_rgpDiff(:,:,:,4,2),this%intSharp_rgpDiff(:,:,:,4,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,5,1),this%intSharp_rgpDiff(:,:,:,5,2),this%intSharp_rgpDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,6,1),this%intSharp_rgpDiff(:,:,:,6,2),this%intSharp_rgpDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,7,1),this%intSharp_rgpDiff(:,:,:,7,2),this%intSharp_rgpDiff(:,:,:,7,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,8,1),this%intSharp_rgpDiff(:,:,:,8,2),this%intSharp_rgpDiff(:,:,:,8,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,9,1),this%intSharp_rgpDiff(:,:,:,9,2),this%intSharp_rgpDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
              !FV terms
              rhsgp = rhsgp + this%intSharp_rgpFV
           endif
        endif

    end subroutine


    subroutine getRHS_pe(this,rho,u,v,w,dt,src,rhspe,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhspe
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: tmp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: LADg
        integer :: i,j,k,l

        rhspe = zero

        !pe rhs
        call gradient(this%decomp,this%der,this%pe,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
        rhspe = - u*LADg(:,:,:,1) - v*LADg(:,:,:,2) - w*LADg(:,:,:,3)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_pe*LADg(:,:,:,1),this%diff_pe*LADg(:,:,:,2),this%diff_pe*LADg(:,:,:,3),tmp, -x_bc, -y_bc, -z_bc)
        rhspe = rhspe + tmp
        
        !rhspe = rhspe / this%T

        if (this%plast) then
           if(this%explPlast) then
               print*,"explict plastic terms not yet implemented for pe"
               stop
           end if
        end if

    end subroutine


    subroutine getRHS_rpe(this,rho,u,v,w,dt,src,rhspe,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhspe
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: tmp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: LADg
        integer :: i,j,k,l

        rhspe = zero

        !pe rhs
        call divergence(this%decomp,this%der,-rho*u*this%pe,-rho*v*this%pe,-rho*w*this%pe,rhspe,-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%pe,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,rho*this%diff_pe*LADg(:,:,:,1),rho*this%diff_pe*LADg(:,:,:,2),rho*this%diff_pe*LADg(:,:,:,3),tmp, -x_bc, -y_bc, -z_bc)
        rhspe = rhspe + tmp

        rhspe = rhspe / this%T

        if (this%plast) then
           if(this%explPlast) then
               print*,"explict plastic terms not yet implemented for pe"
               stop
           end if
        end if

    end subroutine

    subroutine getRHS_gTg(this,rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        ! Symmetry and anti-symmetry properties of gTg are assumed as below (same as g)
        ! In x gTg_{ij}: [S A A; A S S; A S S]
        ! In y gTg_{ij}: [S A S; A S A; S A S]
        ! In z gTg_{ij}: [S S A; S S A; A A S]
    
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        

        rhsg = zero


        detg  = this%g11*(this%g22*this%g33-this%g23*this%g32) &
              - this%g12*(this%g21*this%g33-this%g31*this%g23) &
              + this%g13*(this%g21*this%g32-this%g31*this%g22)
        detg = sqrt(detg)

        ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
        tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
        ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF

        if(present(rho0mix)) then
            penalty  = this%elastic%eta_det_ge*(rho/rho0mix/detg  - one)/dt
            this%det_e = (rho/rho0mix/detg  - one) !diagnostic
        else
            penalty  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
            this%det_e = (tmp /detg /this%elastic%rho0-one) !diagnostic
        endif
        if(this%pRelax) then
           penalty  = this%VF*this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density -- change2
        endif

        where (this%elastic%mu .LT. eps)
            penalty = zero
            this%det_e = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty  = penalty  - src/this%VF
        endif


        !Curl can be used for diagnostic variable
        this%curl_e = zero

        !Transport terms
        call gradient(this%decomp,this%der,this%g11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsg(:,:,:,1) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - 2.0*(this%g11*dudx + this%g12*dvdx + this%g13*dwdx)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp + penalty*this%g11

        !Transport terms
        call gradient(this%decomp,this%der,this%g12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        rhsg(:,:,:,2) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - (this%g11*dudy + this%g12*dvdy + this%g13*dwdy) - (this%g21*dudx + this%g22*dvdx + this%g23*dwdx)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp + penalty*this%g12

        !Transport terms
        call gradient(this%decomp,this%der,this%g13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        rhsg(:,:,:,3) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - (this%g11*dudz + this%g12*dvdz + this%g13*dwdz) - (this%g31*dudx + this%g32*dvdx + this%g33*dwdx)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp + penalty*this%g13


        !Symmetric update
        rhsg(:,:,:,4) = rhsg(:,:,:,2)

        !Transport terms
        call gradient(this%decomp,this%der,this%g22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsg(:,:,:,5) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - 2.0*(this%g21*dudy + this%g22*dvdy + this%g23*dwdy)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp + penalty*this%g22

        !Transport terms
        call gradient(this%decomp,this%der,this%g23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,6) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - (this%g21*dudz + this%g22*dvdz + this%g23*dwdz) - (this%g31*dudy + this%g32*dvdy + this%g33*dwdy)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp + penalty*this%g23

        !Symmetric update
        rhsg(:,:,:,7) = rhsg(:,:,:,3)

        !Symmetric update
        rhsg(:,:,:,8) = rhsg(:,:,:,6)

        !Transport terms
        call gradient(this%decomp,this%der,this%g33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsg(:,:,:,9) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - 2.0*(this%g31*dudz + this%g32*dvdz + this%g33*dwdz)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp + penalty*this%g33


        !RHS update for explicit plastic terms
        if (this%plast) then
           if(this%explPlast) then
              print*,"check this"
              stop
               !call this%getPlasticSources(detg,rhsg)
           end if
        end if

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_cpg_west) then
               !TODO: finish implementation of new sharpening term
               CONTINUE     
           else
           endif
           if(this%intSharp_spf) then
              do i=1,9
                 rhsg(:,:,:,i) = rhsg(:,:,:,i) + this%intSharp_rg(:,:,:,i,1)/rho !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,1,1)/rho,this%intSharp_rgDiff(:,:,:,1,2)/rho,this%intSharp_rgDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,2,1)/rho,this%intSharp_rgDiff(:,:,:,2,2)/rho,this%intSharp_rgDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,3,1)/rho,this%intSharp_rgDiff(:,:,:,3,2)/rho,this%intSharp_rgDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,5,1)/rho,this%intSharp_rgDiff(:,:,:,5,2)/rho,this%intSharp_rgDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,6,1)/rho,this%intSharp_rgDiff(:,:,:,6,2)/rho,this%intSharp_rgDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,9,1)/rho,this%intSharp_rgDiff(:,:,:,9,2)/rho,this%intSharp_rgDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,1,1)/rho,this%intSharp_rg(:,:,:,1,2)/rho,this%intSharp_rg(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,2,1)/rho,this%intSharp_rg(:,:,:,2,2)/rho,this%intSharp_rg(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,3,1)/rho,this%intSharp_rg(:,:,:,3,2)/rho,this%intSharp_rg(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,5,1)/rho,this%intSharp_rg(:,:,:,5,2)/rho,this%intSharp_rg(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,6,1)/rho,this%intSharp_rg(:,:,:,6,2)/rho,this%intSharp_rg(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,9,1)/rho,this%intSharp_rg(:,:,:,9,2)/rho,this%intSharp_rg(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,1,1)/rho,this%intSharp_rgDiff(:,:,:,1,2)/rho,this%intSharp_rgDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,2,1)/rho,this%intSharp_rgDiff(:,:,:,2,2)/rho,this%intSharp_rgDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,3,1)/rho,this%intSharp_rgDiff(:,:,:,3,2)/rho,this%intSharp_rgDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,5,1)/rho,this%intSharp_rgDiff(:,:,:,5,2)/rho,this%intSharp_rgDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,6,1)/rho,this%intSharp_rgDiff(:,:,:,6,2)/rho,this%intSharp_rgDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,9,1)/rho,this%intSharp_rgDiff(:,:,:,9,2)/rho,this%intSharp_rgDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
              !FV terms
              rhsg = rhsg + this%intSharp_gFV

              rhsg(:,:,:,4) = rhsg(:,:,:,2)
              rhsg(:,:,:,7) = rhsg(:,:,:,3)
              rhsg(:,:,:,8) = rhsg(:,:,:,6)
           endif
        endif

    end subroutine


    subroutine getRHS_rgTg(this,rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg, diff_rg
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        ! Symmetry and anti-symmetry properties of gTg are assumed as below (same as g)
        ! In x gTg_{ij}: [S A A; A S S; A S S]
        ! In y gTg_{ij}: [S A S; A S A; S A S]
        ! In z gTg_{ij}: [S S A; S S A; A A S]
    
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        

        rhsg = zero
        diff_rg = rho*this%diff_g

        detg  = this%g11*(this%g22*this%g33-this%g23*this%g32) &
              - this%g12*(this%g21*this%g33-this%g31*this%g23) &
              + this%g13*(this%g21*this%g32-this%g31*this%g22)
        detg = sqrt(detg)

        ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
        tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
        ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF

        if(present(rho0mix)) then
            penalty  = this%elastic%eta_det_ge*(rho/rho0mix/detg  - one)/dt
            this%det_e = (rho/rho0mix/detg  - one) !diagnostic
        else
            penalty  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
            this%det_e = (tmp /detg /this%elastic%rho0-one) !diagnostic
        endif
        if(this%pRelax) then
           penalty  = this%VF*this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density -- change2
        endif

        where (this%elastic%mu .LT. eps)
            penalty = zero
            this%det_e = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty  = penalty  - src/this%VF
        endif


        !Curl can be used for diagnostic variable
        this%curl_e = zero

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg11,-v*this%rg11,-w*this%rg11,rhsg(:,:,:,1),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,1) = rhsg(:,:,:,1) - 2.0*(this%rg11*dudx + this%rg12*dvdx + this%rg13*dwdx)
        !LAD terms
        call gradient(this%decomp,this%der,this%g11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp + penalty*this%rg11

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg12,-v*this%rg12,-w*this%rg12,rhsg(:,:,:,2),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,2) = rhsg(:,:,:,2) - (this%rg11*dudy + this%rg12*dvdy + this%rg13*dwdy) - (this%rg21*dudx + this%rg22*dvdx + this%rg23*dwdx)
        !LAD terms
        call gradient(this%decomp,this%der,this%g12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp + penalty*this%rg12

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg13,-v*this%rg13,-w*this%rg13,rhsg(:,:,:,3),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,3) = rhsg(:,:,:,3) - (this%rg11*dudz + this%rg12*dvdz + this%rg13*dwdz) - (this%rg31*dudx + this%rg32*dvdx + this%rg33*dwdx)
        !LAD terms
        call gradient(this%decomp,this%der,this%g13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp + penalty*this%rg13


        !Symmetric update
        rhsg(:,:,:,4) = rhsg(:,:,:,2)

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg22,-v*this%rg22,-w*this%rg22,rhsg(:,:,:,5),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,5) = rhsg(:,:,:,5) - 2.0*(this%rg21*dudy + this%rg22*dvdy + this%rg23*dwdy)
        !LAD terms
        call gradient(this%decomp,this%der,this%g22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp + penalty*this%rg22

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg23,-v*this%rg23,-w*this%rg23,rhsg(:,:,:,6),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,6) = rhsg(:,:,:,6) - (this%rg21*dudz + this%rg22*dvdz + this%rg23*dwdz) - (this%rg31*dudy + this%rg32*dvdy + this%rg33*dwdy)
        !LAD terms
        call gradient(this%decomp,this%der,this%g23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp + penalty*this%rg23

        !Symmetric update
        rhsg(:,:,:,7) = rhsg(:,:,:,3)

        !Symmetric update
        rhsg(:,:,:,8) = rhsg(:,:,:,6)

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rg33,-v*this%rg33,-w*this%rg33,rhsg(:,:,:,9),-x_bc, -y_bc, -z_bc)
        rhsg(:,:,:,9) = rhsg(:,:,:,9) - 2.0*(this%rg31*dudz + this%rg32*dvdz + this%rg33*dwdz)
        !LAD terms
        call gradient(this%decomp,this%der,this%g33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rg*LADg(:,:,:,1),diff_rg*LADg(:,:,:,2),diff_rg*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp + penalty*this%rg33


        !RHS update for explicit plastic terms
        if (this%plast) then
           if(this%explPlast) then
              print*,"check this"
              stop
               !call this%getPlasticSources(detg,rhsg)
           end if
        end if

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_spf) then
              do i=1,9
                 rhsg(:,:,:,i) = rhsg(:,:,:,i) + this%intSharp_rg(:,:,:,i,1) !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,1,1),this%intSharp_rgDiff(:,:,:,1,2),this%intSharp_rgDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,2,1),this%intSharp_rgDiff(:,:,:,2,2),this%intSharp_rgDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,3,1),this%intSharp_rgDiff(:,:,:,3,2),this%intSharp_rgDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,5,1),this%intSharp_rgDiff(:,:,:,5,2),this%intSharp_rgDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,6,1),this%intSharp_rgDiff(:,:,:,6,2),this%intSharp_rgDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,9,1),this%intSharp_rgDiff(:,:,:,9,2),this%intSharp_rgDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,1,1),this%intSharp_rg(:,:,:,1,2),this%intSharp_rg(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,2,1),this%intSharp_rg(:,:,:,2,2),this%intSharp_rg(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,3,1),this%intSharp_rg(:,:,:,3,2),this%intSharp_rg(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,5,1),this%intSharp_rg(:,:,:,5,2),this%intSharp_rg(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,6,1),this%intSharp_rg(:,:,:,6,2),this%intSharp_rg(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rg(:,:,:,9,1),this%intSharp_rg(:,:,:,9,2),this%intSharp_rg(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,1,1),this%intSharp_rgDiff(:,:,:,1,2),this%intSharp_rgDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,2,1),this%intSharp_rgDiff(:,:,:,2,2),this%intSharp_rgDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,3,1),this%intSharp_rgDiff(:,:,:,3,2),this%intSharp_rgDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,5,1),this%intSharp_rgDiff(:,:,:,5,2),this%intSharp_rgDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,6,1),this%intSharp_rgDiff(:,:,:,6,2),this%intSharp_rgDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgDiff(:,:,:,9,1),this%intSharp_rgDiff(:,:,:,9,2),this%intSharp_rgDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp
              
              !FV terms
              rhsg = rhsg + this%intSharp_rgFV

              rhsg(:,:,:,4) = rhsg(:,:,:,2)
              rhsg(:,:,:,7) = rhsg(:,:,:,3)
              rhsg(:,:,:,8) = rhsg(:,:,:,6)
           endif
        endif

    end subroutine


    subroutine getRHS_gtTgt(this,rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL,P_MINVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsgt
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detgt, diff_g
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz


        ! Symmetry and anti-symmetry properties of gTg are assumed as below (same as g)
        ! In x gTg_{ij}: [S A A; A S S; A S S]
        ! In y gTg_{ij}: [S A S; A S A; S A S]
        ! In z gTg_{ij}: [S S A; S S A; A A S]


        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        

        rhsgt = zero


        detgt  = this%gt11*(this%gt22*this%gt33-this%gt23*this%gt32) &
               - this%gt12*(this%gt21*this%gt33-this%gt31*this%gt23) &
               + this%gt13*(this%gt21*this%gt32-this%gt31*this%gt22)
        detgt = sqrt(detgt)

        ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
        tmp  = (rho*this%Ys + this%elastic%rho0*detgt *epssmall)/(this%VF + epssmall)   
        ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF

        if(present(rho0mix)) then
            penalty  = this%elastic%eta_det_gt*(rho/rho0mix/detgt  - one)/dt
            this%det_t = (rho/rho0mix/detgt  - one) !diagnostic
        else
            penalty  = this%elastic%eta_det_gt*( tmp /detgt /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
            this%det_t = (tmp /detgt /this%elastic%rho0-one) !diagnostic
        endif
        if(this%pRelax) then
           penalty  = this%VF*this%elastic%eta_det_gt*( tmp /detgt /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density -- change2
        endif

        where (this%elastic%mu .LT. eps)
            penalty = zero
            this%det_t = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty  = penalty  - src/this%VF
        endif


        !Curl can be used for diagnostic variable
        this%curl_t = zero

        !Transport terms
        call gradient(this%decomp,this%der,this%gt11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsgt(:,:,:,1) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - 2.0*(this%gt11*dudx + this%gt12*dvdx + this%gt13*dwdx)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp + penalty*this%gt11

        !Transport terms
        call gradient(this%decomp,this%der,this%gt12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        rhsgt(:,:,:,2) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - (this%gt11*dudy + this%gt12*dvdy + this%gt13*dwdy) - (this%gt21*dudx + this%gt22*dvdx + this%gt23*dwdx)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp + penalty*this%gt12

        !Transport terms
        call gradient(this%decomp,this%der,this%gt13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        rhsgt(:,:,:,3) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - (this%gt11*dudz + this%gt12*dvdz + this%gt13*dwdz) - (this%gt31*dudx + this%gt32*dvdx + this%gt33*dwdx)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp + penalty*this%gt13


        !Symmetric update
        rhsgt(:,:,:,4) = rhsgt(:,:,:,2)

        !Transport terms
        call gradient(this%decomp,this%der,this%gt22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsgt(:,:,:,5) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - 2.0*(this%gt21*dudy + this%gt22*dvdy + this%gt23*dwdy)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp + penalty*this%gt22

        !Transport terms
        call gradient(this%decomp,this%der,this%gt23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,6) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - (this%gt21*dudz + this%gt22*dvdz + this%gt23*dwdz) - (this%gt31*dudy + this%gt32*dvdy + this%gt33*dwdy)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp + penalty*this%gt23

        !Symmetric update
        rhsgt(:,:,:,7) = rhsgt(:,:,:,3)

        !Symmetric update
        rhsgt(:,:,:,8) = rhsgt(:,:,:,6)

        !Transport terms
        call gradient(this%decomp,this%der,this%gt33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsgt(:,:,:,9) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3)) - 2.0*(this%gt31*dudz + this%gt32*dvdz + this%gt33*dwdz)
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp + penalty*this%gt33

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_cpg_west) then
               !TODO: finish implementation of new sharpening term
               CONTINUE     
           else
           endif
           if(this%intSharp_spf) then
              do i=1,9
                 rhsgt(:,:,:,i) = rhsgt(:,:,:,i) + this%intSharp_rgt(:,:,:,i,1)/rho !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,1,1)/rho,this%intSharp_rgtDiff(:,:,:,1,2)/rho,this%intSharp_rgtDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,2,1)/rho,this%intSharp_rgtDiff(:,:,:,2,2)/rho,this%intSharp_rgtDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,3,1)/rho,this%intSharp_rgtDiff(:,:,:,3,2)/rho,this%intSharp_rgtDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,5,1)/rho,this%intSharp_rgtDiff(:,:,:,5,2)/rho,this%intSharp_rgtDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,6,1)/rho,this%intSharp_rgtDiff(:,:,:,6,2)/rho,this%intSharp_rgtDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,9,1)/rho,this%intSharp_rgtDiff(:,:,:,9,2)/rho,this%intSharp_rgtDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,1,1)/rho,this%intSharp_rgt(:,:,:,1,2)/rho,this%intSharp_rgt(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,2,1)/rho,this%intSharp_rgt(:,:,:,2,2)/rho,this%intSharp_rgt(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,3,1)/rho,this%intSharp_rgt(:,:,:,3,2)/rho,this%intSharp_rgt(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,5,1)/rho,this%intSharp_rgt(:,:,:,5,2)/rho,this%intSharp_rgt(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,6,1)/rho,this%intSharp_rgt(:,:,:,6,2)/rho,this%intSharp_rgt(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,9,1)/rho,this%intSharp_rgt(:,:,:,9,2)/rho,this%intSharp_rgt(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,1,1)/rho,this%intSharp_rgtDiff(:,:,:,1,2)/rho,this%intSharp_rgtDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,2,1)/rho,this%intSharp_rgtDiff(:,:,:,2,2)/rho,this%intSharp_rgtDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,3,1)/rho,this%intSharp_rgtDiff(:,:,:,3,2)/rho,this%intSharp_rgtDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,5,1)/rho,this%intSharp_rgtDiff(:,:,:,5,2)/rho,this%intSharp_rgtDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,6,1)/rho,this%intSharp_rgtDiff(:,:,:,6,2)/rho,this%intSharp_rgtDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,9,1)/rho,this%intSharp_rgtDiff(:,:,:,9,2)/rho,this%intSharp_rgtDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
              !FV terms
              rhsgt = rhsgt + this%intSharp_gtFV

              rhsgt(:,:,:,4) = rhsgt(:,:,:,2)
              rhsgt(:,:,:,7) = rhsgt(:,:,:,3)
              rhsgt(:,:,:,8) = rhsgt(:,:,:,6)
           endif
        endif

    end subroutine


    subroutine getRHS_rgtTgt(this,rho,u,v,w,duidxj,dt,src,rhsgt,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsgt
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detgt, diff_rgt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        ! Symmetry and anti-symmetry properties of gTg are assumed as below (same as g)
        ! In x gTg_{ij}: [S A A; A S S; A S S]
        ! In y gTg_{ij}: [S A S; A S A; S A S]
        ! In z gTg_{ij}: [S S A; S S A; A A S]
    
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        

        rhsgt = zero
        diff_rgt = rho*this%diff_gt

        detgt  = this%gt11*(this%gt22*this%gt33-this%gt23*this%gt32) &
              - this%gt12*(this%gt21*this%gt33-this%gt31*this%gt23) &
              + this%gt13*(this%gt21*this%gt32-this%gt31*this%gt22)
        detgt = sqrt(detgt)

        ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
        tmp  = (rho*this%Ys + this%elastic%rho0*detgt *epssmall)/(this%VF + epssmall)   
        ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF

        if(present(rho0mix)) then
            penalty  = this%elastic%eta_det_gt*(rho/rho0mix/detgt  - one)/dt
            this%det_t = (rho/rho0mix/detgt  - one) !diagnostic
        else
            penalty  = this%elastic%eta_det_gt*( tmp /detgt /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
            this%det_t = (tmp /detgt /this%elastic%rho0-one) !diagnostic
        endif
        if(this%pRelax) then
           penalty  = this%VF*this%elastic%eta_det_gt*( tmp /detgt /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density -- change2
        endif

        where (this%elastic%mu .LT. eps)
            penalty = zero
            this%det_t = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty  = penalty  - src/this%VF
        endif


        !Curl can be used for diagnostic variable
        this%curl_t = zero

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt11,-v*this%rgt11,-w*this%rgt11,rhsgt(:,:,:,1),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) - 2.0*(this%rgt11*dudx + this%rgt12*dvdx + this%rgt13*dwdx)
        !LAD terms
        call gradient(this%decomp,this%der,this%gt11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp + penalty*this%rgt11

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt12,-v*this%rgt12,-w*this%rgt12,rhsgt(:,:,:,2),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) - (this%rgt11*dudy + this%rgt12*dvdy + this%rgt13*dwdy) - (this%rgt21*dudx + this%rgt22*dvdx + this%rgt23*dwdx)
        !LAD terms
        call gradient(this%decomp,this%der,this%gt12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp + penalty*this%rgt12

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt13,-v*this%rgt13,-w*this%rgt13,rhsgt(:,:,:,3),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) - (this%rgt11*dudz + this%rgt12*dvdz + this%rgt13*dwdz) - (this%rgt31*dudx + this%rgt32*dvdx + this%rgt33*dwdx)
        !LAD terms
        call gradient(this%decomp,this%der,this%gt13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp + penalty*this%rgt13


        !Symmetric update
        rhsgt(:,:,:,4) = rhsgt(:,:,:,2)

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt22,-v*this%rgt22,-w*this%rgt22,rhsgt(:,:,:,5),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) - 2.0*(this%rgt21*dudy + this%rgt22*dvdy + this%rgt23*dwdy)
        !LAD terms
        call gradient(this%decomp,this%der,this%gt22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp + penalty*this%rgt22

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt23,-v*this%rgt23,-w*this%rgt23,rhsgt(:,:,:,6),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) - (this%rgt21*dudz + this%rgt22*dvdz + this%rgt23*dwdz) - (this%rgt31*dudy + this%rgt32*dvdy + this%rgt33*dwdy)
        !LAD terms
        call gradient(this%decomp,this%der,this%gt23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp + penalty*this%rgt23

        !Symmetric update
        rhsgt(:,:,:,7) = rhsgt(:,:,:,3)

        !Symmetric update
        rhsgt(:,:,:,8) = rhsgt(:,:,:,6)

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgt33,-v*this%rgt33,-w*this%rgt33,rhsgt(:,:,:,9),-x_bc, -y_bc, -z_bc)
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) - 2.0*(this%rgt31*dudz + this%rgt32*dvdz + this%rgt33*dwdz)
        !LAD terms
        call gradient(this%decomp,this%der,this%gt33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgt*LADg(:,:,:,1),diff_rgt*LADg(:,:,:,2),diff_rgt*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp + penalty*this%rgt33

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_spf) then
              do i=1,9
                 rhsgt(:,:,:,i) = rhsgt(:,:,:,i) + this%intSharp_rgt(:,:,:,i,1) !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,1,1),this%intSharp_rgtDiff(:,:,:,1,2),this%intSharp_rgtDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,2,1),this%intSharp_rgtDiff(:,:,:,2,2),this%intSharp_rgtDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,3,1),this%intSharp_rgtDiff(:,:,:,3,2),this%intSharp_rgtDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,5,1),this%intSharp_rgtDiff(:,:,:,5,2),this%intSharp_rgtDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,6,1),this%intSharp_rgtDiff(:,:,:,6,2),this%intSharp_rgtDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,9,1),this%intSharp_rgtDiff(:,:,:,9,2),this%intSharp_rgtDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,1,1),this%intSharp_rgt(:,:,:,1,2),this%intSharp_rgt(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,2,1),this%intSharp_rgt(:,:,:,2,2),this%intSharp_rgt(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,3,1),this%intSharp_rgt(:,:,:,3,2),this%intSharp_rgt(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,5,1),this%intSharp_rgt(:,:,:,5,2),this%intSharp_rgt(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,6,1),this%intSharp_rgt(:,:,:,6,2),this%intSharp_rgt(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgt(:,:,:,9,1),this%intSharp_rgt(:,:,:,9,2),this%intSharp_rgt(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,1,1),this%intSharp_rgtDiff(:,:,:,1,2),this%intSharp_rgtDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,2,1),this%intSharp_rgtDiff(:,:,:,2,2),this%intSharp_rgtDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,3,1),this%intSharp_rgtDiff(:,:,:,3,2),this%intSharp_rgtDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,5,1),this%intSharp_rgtDiff(:,:,:,5,2),this%intSharp_rgtDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,6,1),this%intSharp_rgtDiff(:,:,:,6,2),this%intSharp_rgtDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgtDiff(:,:,:,9,1),this%intSharp_rgtDiff(:,:,:,9,2),this%intSharp_rgtDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp
              
              !FV terms
              rhsgt = rhsgt + this%intSharp_rgtFV

              rhsgt(:,:,:,4) = rhsgt(:,:,:,2)
              rhsgt(:,:,:,7) = rhsgt(:,:,:,3)
              rhsgt(:,:,:,8) = rhsgt(:,:,:,6)
           endif
        endif

    end subroutine


    subroutine getRHS_gpTgp(this,rho,u,v,w,duidxj,dt,src,rhsgp,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsgp
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detgp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        ! Symmetry and anti-symmetry properties of gTg are assumed as below (same as g)
        ! In x gTg_{ij}: [S A A; A S S; A S S]
        ! In y gTg_{ij}: [S A S; A S A; S A S]
        ! In z gTg_{ij}: [S S A; S S A; A A S]
    
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        

        rhsgp = zero


        detgp  = this%gp11*(this%gp22*this%gp33-this%gp23*this%gp32) &
              - this%gp12*(this%gp21*this%gp33-this%gp31*this%gp23) &
              + this%gp13*(this%gp21*this%gp32-this%gp31*this%gp22)
        ! !debug
        ! do i=1,this%nxp
        !    do j=1,this%nyp
        !       do k=1,this%nzp
        !          if (detgp(i,j,k) .le. 10e-16) then
        !             print*,i,j,k,detgp(i,j,k),this%gp11(i,j,k),this%gp12(i,j,k),this%gp13(i,j,k),this%gp21(i,j,k),this%gp22(i,j,k),this%gp23(i,j,k),this%gp31(i,j,k),this%gp32(i,j,k),this%gp33(i,j,k)
        !          endif
        !       enddo
        !    enddo
        ! enddo
        ! !end

        detgp = sqrt(detgp)


        penalty = this%elastic%eta_det_gp*(one/detgp - one)/dt
        this%det_p = one/detgp - one

        where (this%elastic%mu .LT. eps)
            penalty = zero
            this%det_p = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty  = penalty  - src/this%VF
        endif


        !Curl can be used for diagnostic variable
        this%curl_p = zero

        !Transport terms
        call gradient(this%decomp,this%der,this%gp11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsgp(:,:,:,1) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3))
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp + penalty*this%gp11

        !Transport terms
        call gradient(this%decomp,this%der,this%gp12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        rhsgp(:,:,:,2) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3))
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp + penalty*this%gp12

        !Transport terms
        call gradient(this%decomp,this%der,this%gp13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        rhsgp(:,:,:,3) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3))
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp + penalty*this%gp13


        !Symmetric update
        rhsgp(:,:,:,4) = rhsgp(:,:,:,2)

        !Transport terms
        call gradient(this%decomp,this%der,this%gp22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsgp(:,:,:,5) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3))
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp + penalty*this%gp22

        !Transport terms
        call gradient(this%decomp,this%der,this%gp23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        rhsgp(:,:,:,6) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3))
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp + penalty*this%gp23

        !Symmetric update
        rhsgp(:,:,:,7) = rhsgp(:,:,:,3)

        !Symmetric update
        rhsgp(:,:,:,8) = rhsgp(:,:,:,6)

        !Transport terms
        call gradient(this%decomp,this%der,this%gp33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        rhsgp(:,:,:,9) = - (u*LADg(:,:,:,1) + v*LADg(:,:,:,2) + w*LADg(:,:,:,3))
        !LAD terms
        call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp + penalty*this%gp33


        !RHS update for explicit plastic terms
        if (this%plast) then
           if(this%explPlast) then
              print*,"check this"
              stop
               !call this%getPlasticSources(detgp,rhsgp)
           end if
        end if

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_cpg_west) then
               !TODO: finish implementation of new sharpening term
               CONTINUE     
           else
           endif
           if(this%intSharp_spf) then
              do i=1,9
                 rhsgp(:,:,:,i) = rhsgp(:,:,:,i) + this%intSharp_rgp(:,:,:,i,1)/rho !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,1,1)/rho,this%intSharp_rgpDiff(:,:,:,1,2)/rho,this%intSharp_rgpDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,2,1)/rho,this%intSharp_rgpDiff(:,:,:,2,2)/rho,this%intSharp_rgpDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,3,1)/rho,this%intSharp_rgpDiff(:,:,:,3,2)/rho,this%intSharp_rgpDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,5,1)/rho,this%intSharp_rgpDiff(:,:,:,5,2)/rho,this%intSharp_rgpDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,6,1)/rho,this%intSharp_rgpDiff(:,:,:,6,2)/rho,this%intSharp_rgpDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,9,1)/rho,this%intSharp_rgpDiff(:,:,:,9,2)/rho,this%intSharp_rgpDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,1,1)/rho,this%intSharp_rgp(:,:,:,1,2)/rho,this%intSharp_rgp(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,2,1)/rho,this%intSharp_rgp(:,:,:,2,2)/rho,this%intSharp_rgp(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,3,1)/rho,this%intSharp_rgp(:,:,:,3,2)/rho,this%intSharp_rgp(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,5,1)/rho,this%intSharp_rgp(:,:,:,5,2)/rho,this%intSharp_rgp(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,6,1)/rho,this%intSharp_rgp(:,:,:,6,2)/rho,this%intSharp_rgp(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,9,1)/rho,this%intSharp_rgp(:,:,:,9,2)/rho,this%intSharp_rgp(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,1,1)/rho,this%intSharp_rgpDiff(:,:,:,1,2)/rho,this%intSharp_rgpDiff(:,:,:,1,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,2,1)/rho,this%intSharp_rgpDiff(:,:,:,2,2)/rho,this%intSharp_rgpDiff(:,:,:,2,3)/rho,tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,3,1)/rho,this%intSharp_rgpDiff(:,:,:,3,2)/rho,this%intSharp_rgpDiff(:,:,:,3,3)/rho,tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,5,1)/rho,this%intSharp_rgpDiff(:,:,:,5,2)/rho,this%intSharp_rgpDiff(:,:,:,5,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,6,1)/rho,this%intSharp_rgpDiff(:,:,:,6,2)/rho,this%intSharp_rgpDiff(:,:,:,6,3)/rho,tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,9,1)/rho,this%intSharp_rgpDiff(:,:,:,9,2)/rho,this%intSharp_rgpDiff(:,:,:,9,3)/rho,tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
              !FV terms
              rhsgp = rhsgp + this%intSharp_gpFV

              rhsgp(:,:,:,4) = rhsgp(:,:,:,2)
              rhsgp(:,:,:,7) = rhsgp(:,:,:,3)
              rhsgp(:,:,:,8) = rhsgp(:,:,:,6)
           endif
        endif

    end subroutine


    subroutine getRHS_rgpTgp(this,rho,u,v,w,duidxj,dt,src,rhsgp,x_bc,y_bc,z_bc,rho0mix)
        use decomp_2d, only: nrank
        use constants, only: eps
        use operators, only: gradient, curl, divergence
        use reductions, only: P_MAXVAL
        class(solid),                                         intent(inout)  :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
        real(rkind),                                          intent(in)  :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsgp
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detgp, diff_rgp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
        integer :: i,j,k,l

        real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        ! Symmetry and anti-symmetry properties of gTg are assumed as below (same as g)
        ! In x gTg_{ij}: [S A A; A S S; A S S]
        ! In y gTg_{ij}: [S A S; A S A; S A S]
        ! In z gTg_{ij}: [S S A; S S A; A A S]
    
        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        

        rhsgp = zero
        diff_rgp = rho*this%diff_gp

        detgp  = this%gp11*(this%gp22*this%gp33-this%gp23*this%gp32) &
              - this%gp12*(this%gp21*this%gp33-this%gp31*this%gp23) &
              + this%gp13*(this%gp21*this%gp32-this%gp31*this%gp22)
        detgp = sqrt(detgp)

        penalty = this%elastic%eta_det_gp*(one/detgp - one)/dt
        this%det_p = one/detgp - one

        where (this%elastic%mu .LT. eps)
            penalty = zero
            this%det_p = zero
        end where

        if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
            ! add Fsource term to penalty 
            penalty  = penalty  - src/this%VF
        endif


        !Curl can be used for diagnostic variable
        this%curl_p = zero

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp11,-v*this%rgp11,-w*this%rgp11,rhsgp(:,:,:,1),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp + penalty*this%rgp11

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp12,-v*this%rgp12,-w*this%rgp12,rhsgp(:,:,:,2),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, -y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,x_bc, y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp + penalty*this%rgp12

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp13,-v*this%rgp13,-w*this%rgp13,rhsgp(:,:,:,3),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),-x_bc, y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,x_bc, -y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp + penalty*this%rgp13


        !Symmetric update
        rhsgp(:,:,:,4) = rhsgp(:,:,:,2)

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp22,-v*this%rgp22,-w*this%rgp22,rhsgp(:,:,:,5),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp + penalty*this%rgp22

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp23,-v*this%rgp23,-w*this%rgp23,rhsgp(:,:,:,6),-x_bc, -y_bc, -z_bc)
        !LAD terms
        call gradient(this%decomp,this%der,this%gp23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, -y_bc, -z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, y_bc, z_bc)
        !RHS update
        rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp + penalty*this%rgp23

        !Symmetric update
        rhsgp(:,:,:,7) = rhsgp(:,:,:,3)

        !Symmetric update
        rhsgp(:,:,:,8) = rhsgp(:,:,:,6)

        !Transport terms
        call divergence(this%decomp,this%der,-u*this%rgp33,-v*this%rgp33,-w*this%rgp33,rhsgp(:,:,:,9),-x_bc, -y_bc, -z_bc)
         !LAD terms
        call gradient(this%decomp,this%der,this%gp33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
        call divergence(this%decomp,this%der,diff_rgp*LADg(:,:,:,1),diff_rgp*LADg(:,:,:,2),diff_rgp*LADg(:,:,:,3),tmp,-x_bc, -y_bc, -z_bc)
        !RHS update
        rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp + penalty*this%rgp33


        !RHS update for explicit plastic terms
        if (this%plast) then
           if(this%explPlast) then
              print*,"check this"
              stop
               !call this%getPlasticSources(detgp,rhsgp)
           end if
        end if

        !RHS update for interface sharpening terms -- once we settle on a version -- don't repeat the divergence calulations -- work it into the fluxes above
        if (this%intSharp) then
           if(this%intSharp_spf) then
              do i=1,9
                 rhsgp(:,:,:,i) = rhsgp(:,:,:,i) + this%intSharp_rgp(:,:,:,i,1) !ignore components 2 and 3 when not in divergence form
              enddo
              
              !high order VF bounds diffusion terms
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,1,1),this%intSharp_rgpDiff(:,:,:,1,2),this%intSharp_rgpDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,2,1),this%intSharp_rgpDiff(:,:,:,2,2),this%intSharp_rgpDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,3,1),this%intSharp_rgpDiff(:,:,:,3,2),this%intSharp_rgpDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,5,1),this%intSharp_rgpDiff(:,:,:,5,2),this%intSharp_rgpDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,6,1),this%intSharp_rgpDiff(:,:,:,6,2),this%intSharp_rgpDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,9,1),this%intSharp_rgpDiff(:,:,:,9,2),this%intSharp_rgpDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
           else
              !low order terms
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,1,1),this%intSharp_rgp(:,:,:,1,2),this%intSharp_rgp(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,2,1),this%intSharp_rgp(:,:,:,2,2),this%intSharp_rgp(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,3,1),this%intSharp_rgp(:,:,:,3,2),this%intSharp_rgp(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,5,1),this%intSharp_rgp(:,:,:,5,2),this%intSharp_rgp(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,6,1),this%intSharp_rgp(:,:,:,6,2),this%intSharp_rgp(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%derD02,this%intSharp_rgp(:,:,:,9,1),this%intSharp_rgp(:,:,:,9,2),this%intSharp_rgp(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
              !high order terms
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,1,1),this%intSharp_rgpDiff(:,:,:,1,2),this%intSharp_rgpDiff(:,:,:,1,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,2,1),this%intSharp_rgpDiff(:,:,:,2,2),this%intSharp_rgpDiff(:,:,:,2,3),tmp, x_bc, y_bc,-z_bc)
              rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,3,1),this%intSharp_rgpDiff(:,:,:,3,2),this%intSharp_rgpDiff(:,:,:,3,3),tmp, x_bc,-y_bc, z_bc)
              rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,5,1),this%intSharp_rgpDiff(:,:,:,5,2),this%intSharp_rgpDiff(:,:,:,5,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,6,1),this%intSharp_rgpDiff(:,:,:,6,2),this%intSharp_rgpDiff(:,:,:,6,3),tmp,-x_bc, y_bc, z_bc)
              rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp
              call divergence(this%decomp,this%der,this%intSharp_rgpDiff(:,:,:,9,1),this%intSharp_rgpDiff(:,:,:,9,2),this%intSharp_rgpDiff(:,:,:,9,3),tmp,-x_bc,-y_bc,-z_bc)
              rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp
              
              !FV terms
              rhsgp = rhsgp + this%intSharp_rgpFV

              rhsgp(:,:,:,4) = rhsgp(:,:,:,2)
              rhsgp(:,:,:,7) = rhsgp(:,:,:,3)
              rhsgp(:,:,:,8) = rhsgp(:,:,:,6)
           endif
        endif

    end subroutine


    ! subroutine getRHS_original(this,rho,u,v,w,dt,src,rhsg,rhsgt,rhsgp,rhspe,x_bc,y_bc,z_bc,rho0mix)
    !     use decomp_2d, only: nrank
    !     use constants, only: eps
    !     use operators, only: gradient, curl, divergence
    !     use reductions, only: P_MAXVAL
    !     class(solid),                                         intent(inout)  :: this
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
    !     real(rkind),                                          intent(in)  :: dt
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg,rhsgt,rhsgp
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhspe
    !     integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, penalty_e, tmp, detg, detgt, detgp
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: curlg,LADg,diffg
    !     integer :: i,j,k,l

    !     ! Symmetry and anti-symmetry properties of g are assumed as below
    !     ! In x g_{ij}: [S A A; A S S; A S S]
    !     ! In y g_{ij}: [S A S; A S A; S A S]
    !     ! In z g_{ij}: [S S A; S S A; A A S]


    !     rhsg = zero
    !     rhsgt = zero
    !     rhsgp = zero
    !     rhspe = zero

    !     detg  = this%g11*(this%g22*this%g33-this%g23*this%g32) &
    !           - this%g12*(this%g21*this%g33-this%g31*this%g23) &
    !           + this%g13*(this%g21*this%g32-this%g31*this%g22)

    !     ! Get the species density = rho*Y/VF (additional terms to give correct limiting behaviour as Ys and VF tend to 0)
    !     tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
    !     ! tmp = rho*this%Ys/(this%VF + epssmall)   ! Get the species density = rho*Y/VF

    !     if(present(rho0mix)) then
    !         penalty  = this%elastic%eta_det_ge*(rho/rho0mix/detg  - one)/dt
    !         this%det_e = (rho/rho0mix/detg  - one) !diagnostic
    !     else
    !         penalty  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
    !         this%det_e = (tmp /detg /this%elastic%rho0-one) !diagnostic
    !     endif
    !     if(this%pRelax) then
    !        penalty  = this%VF*this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density -- change2
    !     endif

    !     where (this%elastic%mu .LT. eps)
    !         penalty = zero
    !         this%det_e = zero
    !     end where

    !     if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
    !         ! add Fsource term to penalty 
    !         penalty  = penalty  - src/this%VF
    !     endif

    !     !Transport terms
    !     tmp = -u*this%g11-v*this%g12-w*this%g13
    !     call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,1),rhsg(:,:,:,2),rhsg(:,:,:,3),-x_bc, y_bc, z_bc)
    !     call curl(this%decomp, this%der, this%g11, this%g12, this%g13, curlg, -x_bc, y_bc, z_bc)
    !     this%curl_e = curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,1) = rhsg(:,:,:,1) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%g11

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,2) = rhsg(:,:,:,2) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%g12

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,3) = rhsg(:,:,:,3) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%g13

    !     !curl diffusion: Miller and Colella
    !     call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,2) = rhsg(:,:,:,2) - this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
    !     rhsg(:,:,:,3) = rhsg(:,:,:,3) + this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,3) = rhsg(:,:,:,3) - this%elastic%diff_c_ge/dt * diffg(:,:,:,1)
    !     rhsg(:,:,:,1) = rhsg(:,:,:,1) + this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,1) = rhsg(:,:,:,1) - this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
    !     rhsg(:,:,:,2) = rhsg(:,:,:,2) + this%elastic%diff_c_ge/dt * diffg(:,:,:,1)


    !     !Transport terms
    !     tmp = -u*this%g21-v*this%g22-w*this%g23
    !     call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,4),rhsg(:,:,:,5),rhsg(:,:,:,6), x_bc,-y_bc, z_bc)
    !     call curl(this%decomp, this%der, this%g21, this%g22, this%g23, curlg, x_bc, -y_bc, z_bc)
    !     this%curl_e = this%curl_e + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,4) = rhsg(:,:,:,4) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%g21

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,5) = rhsg(:,:,:,5) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%g22

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,6) = rhsg(:,:,:,6) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%g23


    !     !curl diffusion: Miller and Colella
    !     call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,5) = rhsg(:,:,:,5) - this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
    !     rhsg(:,:,:,6) = rhsg(:,:,:,6) + this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,6) = rhsg(:,:,:,6) - this%elastic%diff_c_ge/dt * diffg(:,:,:,1)
    !     rhsg(:,:,:,4) = rhsg(:,:,:,4) + this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,4) = rhsg(:,:,:,4) - this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
    !     rhsg(:,:,:,5) = rhsg(:,:,:,5) + this%elastic%diff_c_ge/dt * diffg(:,:,:,1)


    !     !Transport terms
    !     tmp = -u*this%g31-v*this%g32-w*this%g33
    !     call gradient(this%decomp,this%der,tmp,rhsg(:,:,:,7),rhsg(:,:,:,8),rhsg(:,:,:,9), x_bc, y_bc,-z_bc)
    !     call curl(this%decomp, this%der, this%g31, this%g32, this%g33, curlg, x_bc, y_bc, -z_bc)
    !     this%curl_e = this%curl_e + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2
    !     this%curl_e = sqrt(this%curl_e)

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,7) = rhsg(:,:,:,7) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%g31

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,8) = rhsg(:,:,:,8) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%g32

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%g33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_g*LADg(:,:,:,1),this%diff_g*LADg(:,:,:,2),this%diff_g*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsg(:,:,:,9) = rhsg(:,:,:,9) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%g33

    !     !curl diffusion: Miller and Colella
    !     call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,8) = rhsg(:,:,:,8) - this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
    !     rhsg(:,:,:,9) = rhsg(:,:,:,9) + this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,9) = rhsg(:,:,:,9) - this%elastic%diff_c_ge/dt * diffg(:,:,:,1)
    !     rhsg(:,:,:,7) = rhsg(:,:,:,7) + this%elastic%diff_c_ge/dt * diffg(:,:,:,3)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,7) = rhsg(:,:,:,7) - this%elastic%diff_c_ge/dt * diffg(:,:,:,2)
    !     rhsg(:,:,:,8) = rhsg(:,:,:,8) + this%elastic%diff_c_ge/dt * diffg(:,:,:,1)




    !     !pe rhs
    !     call divergence(this%decomp,this%der,-u*this%pe,-v*this%pe,-w*this%pe,rhspe,x_bc, y_bc, z_bc)
    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%pe,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_pe*LADg(:,:,:,1),this%diff_pe*LADg(:,:,:,2),this%diff_pe*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     rhspe = rhspe + tmp


    !     !g_t rhs
    !     detgt = this%gt11*(this%gt22*this%gt33-this%gt23*this%gt32) &
    !          - this%gt12*(this%gt21*this%gt33-this%gt31*this%gt23) &
    !          + this%gt13*(this%gt21*this%gt32-this%gt31*this%gt22)

    !     tmp = (rho*this%Ys + this%elastic%rho0*detgt*epssmall)/(this%VF + epssmall)   

    !     if(present(rho0mix)) then
    !        penalty = this%elastic%eta_det_gt*(rho/rho0mix/detgt - one)/dt
    !        this%det_t = (rho/rho0mix/detgt - one) !diagnostic
    !     else
    !        penalty = this%elastic%eta_det_gt*( tmp/detgt/this%elastic%rho0-one)/dt
    !        this%det_t = (tmp/detgt/this%elastic%rho0 - one)
    !     endif
    !     if(this%pRelax) then
    !        penalty = this%VF*this%elastic%eta_det_gt*( tmp/detgt/this%elastic%rho0-one)/dt
    !     endif

    !     where (this%elastic%mu .LT. eps)
    !        penalty = zero
    !        this%det_t = zero
    !     end where

    !     if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
    !        ! add Fsource term to penalty 
    !        penalty = penalty - src/this%VF
    !     endif


    !     !Transport terms
    !     tmp = -u*this%gt11-v*this%gt12-w*this%gt13
    !     call gradient(this%decomp,this%der,tmp,rhsgt(:,:,:,1),rhsgt(:,:,:,2),rhsgt(:,:,:,3),-x_bc, y_bc, z_bc)
    !     call curl(this%decomp, this%der, this%gt11, this%gt12, this%gt13, curlg, -x_bc, y_bc, z_bc)
    !     this%curl_t = curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%gt11

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%gt12 

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%gt13 

    !     !curl diffusion: Miller and Colella
    !     call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,2) = rhsgt(:,:,:,2) - this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
    !     rhsgt(:,:,:,3) = rhsgt(:,:,:,3) + this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,3) = rhsgt(:,:,:,3) - this%elastic%diff_c_gt/dt * diffg(:,:,:,1)
    !     rhsgt(:,:,:,1) = rhsgt(:,:,:,1) + this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,1) = rhsgt(:,:,:,1) - this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
    !     rhsgt(:,:,:,2) = rhsgt(:,:,:,2) + this%elastic%diff_c_gt/dt * diffg(:,:,:,1)


    !     !Transport terms
    !     tmp = -u*this%gt21-v*this%gt22-w*this%gt23
    !     call gradient(this%decomp,this%der,tmp,rhsgt(:,:,:,4),rhsgt(:,:,:,5),rhsgt(:,:,:,6), x_bc,-y_bc, z_bc)   
    !     call curl(this%decomp, this%der, this%gt21, this%gt22, this%gt23, curlg, x_bc, -y_bc, z_bc)
    !     this%curl_t = this%curl_t + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%gt21 

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%gt22 

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%gt23 

    !     !curl diffusion: Miller and Colella
    !     call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,5) = rhsgt(:,:,:,5) - this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
    !     rhsgt(:,:,:,6) = rhsgt(:,:,:,6) + this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,6) = rhsgt(:,:,:,6) - this%elastic%diff_c_gt/dt * diffg(:,:,:,1)
    !     rhsgt(:,:,:,4) = rhsgt(:,:,:,4) + this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,4) = rhsgt(:,:,:,4) - this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
    !     rhsgt(:,:,:,5) = rhsgt(:,:,:,5) + this%elastic%diff_c_gt/dt * diffg(:,:,:,1)


    !     !Transport terms
    !     tmp = -u*this%gt31-v*this%gt32-w*this%gt33
    !     call gradient(this%decomp,this%der,tmp,rhsgt(:,:,:,7),rhsgt(:,:,:,8),rhsgt(:,:,:,9), x_bc, y_bc,-z_bc)
    !     call curl(this%decomp, this%der, this%gt31, this%gt32, this%gt33, curlg, x_bc, y_bc, -z_bc)
    !     this%curl_t = this%curl_t + curlg(:,:,:,1)**2 + curlg(:,:,:,2)**2 + curlg(:,:,:,3)**2
    !     this%curl_t = sqrt(this%curl_t)

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + tmp + (v*curlg(:,:,:,3) - w*curlg(:,:,:,2)) + penalty*this%gt31 

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + tmp + (w*curlg(:,:,:,1) - u*curlg(:,:,:,3)) + penalty*this%gt32 

    !     !LAD terms
    !     call gradient(this%decomp,this%der,this%gt33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     call divergence(this%decomp,this%der,this%diff_gt*LADg(:,:,:,1),this%diff_gt*LADg(:,:,:,2),this%diff_gt*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + tmp + (u*curlg(:,:,:,2) - v*curlg(:,:,:,1)) + penalty*this%gt33 

    !     !curl diffusion: Miller and Colella
    !     call gradient(this%decomp,this%der,curlg(:,:,:,1),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,8) = rhsgt(:,:,:,8) - this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
    !     rhsgt(:,:,:,9) = rhsgt(:,:,:,9) + this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,2),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,9) = rhsgt(:,:,:,9) - this%elastic%diff_c_gt/dt * diffg(:,:,:,1)
    !     rhsgt(:,:,:,7) = rhsgt(:,:,:,7) + this%elastic%diff_c_gt/dt * diffg(:,:,:,3)
    !     call gradient(this%decomp,this%der,curlg(:,:,:,3),diffg(:,:,:,1),diffg(:,:,:,2),diffg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgt(:,:,:,7) = rhsgt(:,:,:,7) - this%elastic%diff_c_gt/dt * diffg(:,:,:,2)
    !     rhsgt(:,:,:,8) = rhsgt(:,:,:,8) + this%elastic%diff_c_gt/dt * diffg(:,:,:,1)


    !     !g_p rhs
    !     detgp = this%gp11*(this%gp22*this%gp33-this%gp23*this%gp32) &
    !          - this%gp12*(this%gp21*this%gp33-this%gp31*this%gp23) &
    !          + this%gp13*(this%gp21*this%gp32-this%gp31*this%gp22)

    !     ! !use this to set g^p correction in terms of g^e and g^t corrections
    !     ! tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall)   
    !     ! if(present(rho0mix)) then
    !     !    penalty_e  = this%elastic%eta_det_ge*(rho/rho0mix/detg)
    !     ! else
    !     !    penalty_e  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0)
    !     ! endif
    !     ! penalty = (penalty_e*(this%elastic%eta_det_gp/detgp - one) + this%elastic%eta_det_ge*(1.0-this%elastic%eta_det_gp))/dt

    !     !independent g^p correction
    !     penalty = this%elastic%eta_det_gp*(one/detgp - one)/dt

    !     this%det_p = one/detgp - one

    !     where (this%elastic%mu .LT. eps)
    !        penalty = zero
    !        this%det_p = zero
    !     end where

    !     if(this%pEqb) then  !--actually, these source terms should be included for PTeqb as well -- NSG
    !        ! add Fsource term to penalty 
    !        penalty = penalty - src/this%VF !mca check???
    !     endif

    !     this%curl_p = 0 !free variable to track compatability condition -- not set up
             

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp11,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     rhsgp(:,:,:,1) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     !LAD terms
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,1) = rhsgp(:,:,:,1) + tmp + penalty*this%gp11 

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp12,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     rhsgp(:,:,:,2) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     !LAD terms
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,2) = rhsgp(:,:,:,2) + tmp + penalty*this%gp12

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp13,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3),x_bc, y_bc, z_bc)
    !     rhsgp(:,:,:,3) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     !LAD terms
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp,x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,3) = rhsgp(:,:,:,3) + tmp + penalty*this%gp13 

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp21,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgp(:,:,:,4) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     !LAD terms
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,4) = rhsgp(:,:,:,4) + tmp + penalty*this%gp21

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp22,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgp(:,:,:,5) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     !LAD terms
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,5) = rhsgp(:,:,:,5) + tmp + penalty*this%gp22 

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp23,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgp(:,:,:,6) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     !LAD terms
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,6) = rhsgp(:,:,:,6) + tmp + penalty*this%gp23 

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp31,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgp(:,:,:,7) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     !LAD terms
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,7) = rhsgp(:,:,:,7) + tmp + penalty*this%gp31 

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp32,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsgp(:,:,:,8) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     !LAD terms
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,8) = rhsgp(:,:,:,8) + tmp + penalty*this%gp32 

    !     !Transport terms
    !     call gradient(this%decomp,this%der,this%gp33,LADg(:,:,:,1),LADg(:,:,:,2),LADg(:,:,:,3), x_bc, y_bc, z_bc)
    !     !LAD terms
    !     rhsgp(:,:,:,9) = -u*LADg(:,:,:,1)-v*LADg(:,:,:,2)-w*LADg(:,:,:,3)
    !     call divergence(this%decomp,this%der,this%diff_gp*LADg(:,:,:,1),this%diff_gp*LADg(:,:,:,2),this%diff_gp*LADg(:,:,:,3),tmp, x_bc, y_bc, z_bc)
    !     !RHS update
    !     rhsgp(:,:,:,9) = rhsgp(:,:,:,9) + tmp + penalty*this%gp33

            
    !     if (this%plast) then
    !        tmp = (rho*this%Ys + this%elastic%rho0*detg*epssmall)/(this%VF + epssmall)
    !        if(present(rho0mix)) then
    !           tmp = rho/rho0mix
    !        else
    !           tmp = tmp/this%elastic%rho0
    !        endif
    !        call getModifiedModulii(this,tmp)

    !        if(this%explPlast) then
    !            call this%getPlasticSources(detg,rhsg)
    !        end if
    !     end if

    ! end subroutine



    subroutine implicit_plastic(this,rho,rho0mix,mumix,yieldmix)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in), optional  :: rho0mix, mumix, yieldmix
        real(rkind) :: max_modDevSigma


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

                if(present(yieldmix) .and. present(mumix) .and. present(rho0mix)) then
                    call this%elastic%plastic_deformation(this%g, this%g_p, this%pe, rho, this%T, this%sxx, this%sxy, this%sxz, this%syy, this%syz, this%szz,this%use_gTg,this%useOneG,this%strainHard,this%cnsrv_g,this%cnsrv_gt,this%cnsrv_gp,this%cnsrv_pe, mumix, yieldmix, rho0mix)
                else
                    call this%elastic%plastic_deformation(this%g, this%g_p, this%pe, rho, this%T, this%sxx, this%sxy, this%sxz, this%syy, this%syz, this%szz,this%use_gTg,this%useOneG,this%strainHard,this%cnsrv_g,this%cnsrv_gt,this%cnsrv_gp,this%cnsrv_pe)
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
                
            endif
        end if
    end subroutine



    ! subroutine update_gTg(this,isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,rho0mix,mumix,yieldmix,solidVF)
    !     use RKCoeffs,   only: RK45_A,RK45_B
    !     class(solid), intent(inout) :: this
    !     integer, intent(in) :: isub
    !     real(rkind), intent(in) :: dt,tsim
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
    !     integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in), optional  :: rho0mix, mumix, yieldmix, solidVF

    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp,9) :: rhsg,rhsgt,rhsgp  ! RHS for gTg tensor equation

    !     if(present(rho0mix) .and. present(solidVF)) then
    !         call this%getRHS_gTg(rho,u,v,w,dt,src,rhsg,rhsgt,rhsgp,x_bc,y_bc,z_bc,rho0mix,solidVF)
    !     else
    !         call this%getRHS_gTg(rho,u,v,w,dt,src,rhsg,rhsgt,rhsgp,x_bc,y_bc,z_bc)
    !     endif
    !     call hook_material_g_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsg)
    !     call hook_material_g_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsgt)

    !     ! advance sub-step g
    !     if(isub==1) this%Qtmpg = zero                   ! not really needed, since RK45_A(1) = 0
    !     this%Qtmpg  = dt*rhsg + RK45_A(isub)*this%Qtmpg
    !     this%g = this%g  + RK45_B(isub)*this%Qtmpg

    !     ! advance sub-step g_t
    !     if(isub==1) this%Qtmpg_t = zero                   ! not really needed, since RK45_A(1) = 0
    !     this%Qtmpg_t  = dt*rhsgt + RK45_A(isub)*this%Qtmpg_t
    !     this%g_t = this%g_t  + RK45_B(isub)*this%Qtmpg_t
        
    !     ! advance sub-step g_p
    !     if(isub==1) this%Qtmpg_p = zero                   ! not really needed, since RK45_A(1) = 0
    !     this%Qtmpg_p  = dt*rhsgp + RK45_A(isub)*this%Qtmpg_p
    !     this%g_p = this%g_p  + RK45_B(isub)*this%Qtmpg_p

    !     if(this%plast) then
    !         if (.NOT. this%explPlast) then
    !             if(present(yieldmix) .and. present(mumix)) then
    !                !call this%elastic%plastic_deformation(this%g, this%use_gTg, mumix, yieldmix)
    !                !call this%elastic%plastic_deformation(this%g, this%pe, this%use_gTg, rho, this%T, tdel, this%sxx, this%sxy, this%sxz, this%syy, this%syz, this%szz, mumix, yieldmix, rho0mix)
    !                print*,'need to fix gTg for plastic deformation'
    !                stop
    !             else
    !                !call this%elastic%plastic_deformation(this%g, this%use_gTg)
    !                !call this%elastic%plastic_deformation(this%g, this%pe, this%use_gTg, rho, this%T, tdel, this%sxx, this%sxy, this%sxz, this%syy, this%syz, this%szz)
    !                print*,'need to fix gTg for plastic deformation'
    !                stop
    !             endif
    !         end if
    !     end if

    ! end subroutine

    ! !original gTg
    ! !subroutine getRHS_gTg(this,rho,u,v,w,dt,src,rhsg,rhsgt,rhsgp,x_bc,y_bc,z_bc,rho0mix,solidVF)
    ! subroutine getRHS_gTg(this,rho,u,v,w,duidxj,dt,src,rhsg,x_bc,y_bc,z_bc,rho0mix,solidVF)
    !     use operators, only: gradient, curl, divergence
    !     !class(solid),                                         intent(in)  :: this
    !     class(solid),                                         intent(inout)  :: this !mca
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w,src
    !     real(rkind),                                          intent(in)  :: dt
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp,9), intent(out) :: rhsg!,rhsgt,rhsgp
    !     integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in), optional :: rho0mix, solidVF

    !     real(rkind), dimension(this%nxp, this%nyp, this%nzp,9), target, intent(in) :: duidxj
    !     real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: penalty, tmp, detg, penaltyt, tmpt, detgt, detgp
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradG,gradGt,gradGp
    !     !real(rkind), parameter :: etafac = one/6._rkind
    !     integer :: i

    !     dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
    !     dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
    !     dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
    !     ! call gradient(this%decomp, this%der, u, dudx, dudy, dudz, -x_bc,  y_bc,  z_bc)
    !     ! call gradient(this%decomp, this%der, v, dvdx, dvdy, dvdz,  x_bc, -y_bc,  z_bc)
    !     ! call gradient(this%decomp, this%der, w, dwdx, dwdy, dwdz,  x_bc,  y_bc, -z_bc)

    !     ! Symmetry and anti-symmetry properties of gTg are assumed as below (same as g)
    !     ! In x gTg_{ij}: [S A A; A S S; A S S]
    !     ! In y gTg_{ij}: [S A S; A S A; S A S]
    !     ! In z gTg_{ij}: [S S A; S S A; A A S]

    !     rhsg = zero
    !     !rhsgt = zero
    !     !rhsgp = zero

    !     detG = this%G11*(this%G22*this%G33-this%G23*this%G32) &
    !          - this%G12*(this%G21*this%G33-this%G31*this%G23) &
    !          + this%G13*(this%G21*this%G32-this%G31*this%G22)

    !     detg = sqrt(detg)
    !     tmp  = (rho*this%Ys + this%elastic%rho0*detg *epssmall)/(this%VF + epssmall) ! species density 

    !     if(present(rho0mix)) then
    !         penalty  = this%elastic%eta_det_ge*( rho/rho0mix/detg -one)/dt ! Penalty term to keep g consistent with species density
    !     else
    !         penalty  = this%elastic%eta_det_ge*( tmp /detg /this%elastic%rho0-one)/dt ! Penalty term to keep g consistent with species density
    !     endif

    !     call gradient(this%decomp, this%der, this%G11, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,1) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
    !                     - (dudx*this%G11 + dvdx*this%G21 + dwdx*this%G31) &
    !                     - (this%G11*dudx + this%G12*dvdx + this%G13*dwdx) + penalty*this%G11

    !     call gradient(this%decomp, this%der, this%G12, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3),-x_bc,-y_bc, z_bc)
    !     rhsg(:,:,:,2) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
    !                     - (dudx*this%G12 + dvdx*this%G22 + dwdx*this%G32) &
    !                     - (this%G11*dudy + this%G12*dvdy + this%G13*dwdy) + penalty*this%G12

    !     call gradient(this%decomp, this%der, this%G13, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3),-x_bc, y_bc,-z_bc)
    !     rhsg(:,:,:,3) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
    !                     - (dudx*this%G13 + dvdx*this%G23 + dwdx*this%G33) &
    !                     - (this%G11*dudz + this%G12*dvdz + this%G13*dwdz) + penalty*this%G13

    !     rhsg(:,:,:,4) = rhsg(:,:,:,2)  ! Since symmetric

    !     call gradient(this%decomp, this%der, this%G22, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,5) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
    !                     - (dudy*this%G12 + dvdy*this%G22 + dwdy*this%G32) &
    !                     - (this%G21*dudy + this%G22*dvdy + this%G23*dwdy) + penalty*this%G22

    !     call gradient(this%decomp, this%der, this%G23, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3), x_bc,-y_bc,-z_bc)
    !     rhsg(:,:,:,6) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
    !                     - (dudy*this%G13 + dvdy*this%G23 + dwdy*this%G33) &
    !                     - (this%G21*dudz + this%G22*dvdz + this%G23*dwdz) + penalty*this%G23

    !     rhsg(:,:,:,7) = rhsg(:,:,:,3)  ! Since symmetric
    !     rhsg(:,:,:,8) = rhsg(:,:,:,6)  ! Since symmetric

    !     call gradient(this%decomp, this%der, this%G33, gradG(:,:,:,1), gradG(:,:,:,2), gradG(:,:,:,3), x_bc, y_bc, z_bc)
    !     rhsg(:,:,:,9) = - (u*gradG(:,:,:,1) + v*gradG(:,:,:,2) + w*gradG(:,:,:,3)) &
    !                     - (dudz*this%G13 + dvdz*this%G23 + dwdz*this%G33) &
    !                     - (this%G31*dudz + this%G32*dvdz + this%G33*dwdz) + penalty*this%G33




    !     ! !g_t
    !     ! detGt = this%Gt11*(this%Gt22*this%Gt33-this%Gt23*this%Gt32) &
    !     !      - this%Gt12*(this%Gt21*this%Gt33-this%Gt31*this%Gt23) &
    !     !      + this%Gt13*(this%Gt21*this%Gt32-this%Gt31*this%Gt22)
        
    !     ! detgt = sqrt(detgt)
    !     ! tmpt = (rho*this%Ys + this%elastic%rho0*detgt*epssmall)/(this%VF + epssmall)
        
    !     ! if(present(rho0mix)) then
    !     !    penaltyt = this%elastic%eta_det_gt*( rho/rho0mix/detgt-one)/dt
    !     ! else
    !     !    penaltyt = this%elastic%eta_det_gt*( tmpt/detgt/this%elastic%rho0-one)/dt
    !     ! endif

    !     ! !plastic g terms
    !     ! call gradient(this%decomp, this%der, this%Gt11, gradGt(:,:,:,1), gradGt(:,:,:,2), gradGt(:,:,:,3), x_bc, y_bc, z_bc)
    !     ! rhsgt(:,:,:,1) = - (u*gradGt(:,:,:,1) + v*gradGt(:,:,:,2) + w*gradGt(:,:,:,3)) &
    !     !      - (dudx*this%Gt11 + dvdx*this%Gt21 + dwdx*this%Gt31) &
    !     !      - (this%Gt11*dudx + this%Gt12*dvdx + this%Gt13*dwdx) + penaltyt*this%Gt11

    !     ! call gradient(this%decomp, this%der, this%Gt12, gradGt(:,:,:,1), gradGt(:,:,:,2), gradGt(:,:,:,3),-x_bc,-y_bc, z_bc)
    !     ! rhsgt(:,:,:,2) = - (u*gradGt(:,:,:,1) + v*gradGt(:,:,:,2) + w*gradGt(:,:,:,3)) &
    !     !      - (dudx*this%Gt12 + dvdx*this%Gt22 + dwdx*this%Gt32) &
    !     !      - (this%Gt11*dudy + this%Gt12*dvdy + this%Gt13*dwdy) + penaltyt*this%Gt12

    !     ! call gradient(this%decomp, this%der, this%Gt13, gradGt(:,:,:,1), gradGt(:,:,:,2), gradGt(:,:,:,3),-x_bc, y_bc,-z_bc)
    !     ! rhsgt(:,:,:,3) = - (u*gradGt(:,:,:,1) + v*gradGt(:,:,:,2) + w*gradGt(:,:,:,3)) &
    !     !      - (dudx*this%Gt13 + dvdx*this%Gt23 + dwdx*this%Gt33) &
    !     !      - (this%Gt11*dudz + this%Gt12*dvdz + this%Gt13*dwdz) + penaltyt*this%Gt13

    !     ! rhsgt(:,:,:,4) = rhsgt(:,:,:,2)  ! Since symmetric

    !     ! call gradient(this%decomp, this%der, this%Gt22, gradGt(:,:,:,1), gradGt(:,:,:,2), gradGt(:,:,:,3), x_bc, y_bc, z_bc)
    !     ! rhsgt(:,:,:,5) = - (u*gradGt(:,:,:,1) + v*gradGt(:,:,:,2) + w*gradGt(:,:,:,3)) &
    !     !      - (dudy*this%Gt12 + dvdy*this%Gt22 + dwdy*this%Gt32) &
    !     !      - (this%Gt21*dudy + this%Gt22*dvdy + this%Gt23*dwdy) + penaltyt*this%Gt22

    !     ! call gradient(this%decomp, this%der, this%Gt23, gradGt(:,:,:,1), gradGt(:,:,:,2), gradGt(:,:,:,3), x_bc,-y_bc,-z_bc)
    !     ! rhsgt(:,:,:,6) = - (u*gradGt(:,:,:,1) + v*gradGt(:,:,:,2) + w*gradGt(:,:,:,3)) &
    !     !      - (dudy*this%Gt13 + dvdy*this%Gt23 + dwdy*this%Gt33) &
    !     !      - (this%Gt21*dudz + this%Gt22*dvdz + this%Gt23*dwdz) + penaltyt*this%Gt23

    !     ! rhsgt(:,:,:,7) = rhsgt(:,:,:,3)  ! Since symmetric
    !     ! rhsgt(:,:,:,8) = rhsgt(:,:,:,6)  ! Since symmetric

    !     ! call gradient(this%decomp, this%der, this%Gt33, gradGt(:,:,:,1), gradGt(:,:,:,2), gradGt(:,:,:,3), x_bc, y_bc, z_bc)
    !     ! rhsgt(:,:,:,9) = - (u*gradGt(:,:,:,1) + v*gradGt(:,:,:,2) + w*gradGt(:,:,:,3)) &
    !     !      - (dudz*this%Gt13 + dvdz*this%Gt23 + dwdz*this%Gt33) &
    !     !      - (this%Gt31*dudz + this%Gt32*dvdz + this%Gt33*dwdz) + penaltyt*this%Gt33


    !     ! !g_p


    !     ! detGp = this%Gp11*(this%Gp22*this%Gp33-this%Gp23*this%Gp32) &
    !     !      - this%Gp12*(this%Gp21*this%Gp33-this%Gp31*this%Gp23) &
    !     !      + this%Gp13*(this%Gp21*this%Gp32-this%Gp31*this%Gp22)

    !     ! detgp = sqrt(detgp)
    !     ! tmpt = (rho*this%Ys + this%elastic%rho0*detgp*epssmall)/(this%VF + epssmall)

    !     ! ! if(present(rho0mix)) then
    !     ! !    penaltyt = this%elastic%eta_det*( rho/rho0mix/detgp-one)/dt
    !     ! ! else
    !     ! !    penaltyt = this%elastic%eta_det*( tmpt/detgp/this%elastic%rho0-one)/dt
    !     ! ! endif

    !     ! !detgp = 1 condition
    !     ! penaltyt = this%elastic%eta_det_gp*(one/detgp-one)/dt

    !     ! !plastic g terms
    !     ! call gradient(this%decomp, this%der, this%Gp11, gradGp(:,:,:,1), gradGp(:,:,:,2), gradGp(:,:,:,3), x_bc, y_bc, z_bc)
    !     ! rhsgp(:,:,:,1) = - (u*gradGp(:,:,:,1) + v*gradGp(:,:,:,2) + w*gradGp(:,:,:,3)) &
    !     !      - (dudx*this%Gp11 + dvdx*this%Gp21 + dwdx*this%Gp31) &
    !     !      - (this%Gp11*dudx + this%Gp12*dvdx + this%Gp13*dwdx) + penaltyt*this%Gp11

    !     ! call gradient(this%decomp, this%der, this%Gp12, gradGp(:,:,:,1), gradGp(:,:,:,2), gradGp(:,:,:,3),-x_bc,-y_bc, z_bc)
    !     ! rhsgp(:,:,:,2) = - (u*gradGp(:,:,:,1) + v*gradGp(:,:,:,2) + w*gradGp(:,:,:,3)) &
    !     !      - (dudx*this%Gp12 + dvdx*this%Gp22 + dwdx*this%Gp32) &
    !     !      - (this%Gp11*dudy + this%Gp12*dvdy + this%Gp13*dwdy) + penaltyt*this%Gp12

    !     ! call gradient(this%decomp, this%der, this%Gp13, gradGp(:,:,:,1), gradGp(:,:,:,2), gradGp(:,:,:,3),-x_bc, y_bc,-z_bc)
    !     ! rhsgp(:,:,:,3) = - (u*gradGp(:,:,:,1) + v*gradGp(:,:,:,2) + w*gradGp(:,:,:,3)) &
    !     !      - (dudx*this%Gp13 + dvdx*this%Gp23 + dwdx*this%Gp33) &
    !     !      - (this%Gp11*dudz + this%Gp12*dvdz + this%Gp13*dwdz) + penaltyt*this%Gp13

    !     ! rhsgp(:,:,:,4) = rhsgp(:,:,:,2)  ! Since symmetric

    !     ! call gradient(this%decomp, this%der, this%Gp22, gradGp(:,:,:,1), gradGp(:,:,:,2), gradGp(:,:,:,3), x_bc, y_bc, z_bc)
    !     ! rhsgp(:,:,:,5) = - (u*gradGp(:,:,:,1) + v*gradGp(:,:,:,2) + w*gradGp(:,:,:,3)) &
    !     !      - (dudy*this%Gp12 + dvdy*this%Gp22 + dwdy*this%Gp32) &
    !     !      - (this%Gp21*dudy + this%Gp22*dvdy + this%Gp23*dwdy) + penaltyt*this%Gp22

    !     ! call gradient(this%decomp, this%der, this%Gp23, gradGp(:,:,:,1), gradGp(:,:,:,2), gradGp(:,:,:,3), x_bc,-y_bc,-z_bc)
    !     ! rhsgp(:,:,:,6) = - (u*gradGp(:,:,:,1) + v*gradGp(:,:,:,2) + w*gradGp(:,:,:,3)) &
    !     !      - (dudy*this%Gp13 + dvdy*this%Gp23 + dwdy*this%Gp33) &
    !     !      - (this%Gp21*dudz + this%Gp22*dvdz + this%Gp23*dwdz) + penaltyt*this%Gp23

    !     ! rhsgp(:,:,:,7) = rhsgp(:,:,:,3)  ! Since symmetric
    !     ! rhsgp(:,:,:,8) = rhsgp(:,:,:,6)  ! Since symmetric

    !     ! call gradient(this%decomp, this%der, this%Gp33, gradGp(:,:,:,1), gradGp(:,:,:,2), gradGp(:,:,:,3), x_bc, y_bc, z_bc)
    !     ! rhsgp(:,:,:,9) = - (u*gradGp(:,:,:,1) + v*gradGp(:,:,:,2) + w*gradGp(:,:,:,3)) &
    !     !      - (dudz*this%Gp13 + dvdz*this%Gp23 + dwdz*this%Gp33) &
    !     !      - (this%Gp31*dudz + this%Gp32*dvdz + this%Gp33*dwdz) + penaltyt*this%Gp33


    !     ! if (this%plast) then
    !     !    !mca: check if this is correct for gTg
    !     !    tmp = (rho*this%Ys + this%elastic%rho0*detg*epssmall)/(this%VF + epssmall)
    !     !    if(present(rho0mix)) then
    !     !       tmp = rho/rho0mix
    !     !    else
    !     !       tmp = tmp/this%elastic%rho0
    !     !    endif

    !     !     call getModifiedModulii(this) !mca
    !     !     if (this%explPlast) then
    !     !         call this%getPlasticSources(detg,rhsg)
    !     !     end if
    !     ! end if

    !     if(present(solidVF)) then
    !        do i = 1, 9
    !            rhsg(:,:,:,i) = rhsg(:,:,:,i) * solidVF
    !            !rhsgt(:,:,:,i) = rhsgt(:,:,:,i) * solidVF
    !            !rhsgp(:,:,:,i) = rhsgp(:,:,:,i) * solidVF
    !        enddo
    !     endif

    ! end subroutine


    subroutine getModifiedModulii(this)
        class(solid),                                         intent(inout)    :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)    :: fm,wrk11,wrk12,wrk13,wrk21,wrk22,wrk23,wrk31,wrk32,wrk33,detg,gp11,gp12,gp13,gp21,gp22,gp23,gp31,gp32,gp33,pmu,fmu,hmu,pyield,fyield,hyield,et
        character(len=clen) :: charout

        wrk11 = this%g22*this%g33 - this%g23*this%g32
        wrk21 = this%g23*this%g31 - this%g21*this%g33
        wrk31 = this%g21*this%g32 - this%g22*this%g31
        wrk12 = this%g13*this%g32 - this%g12*this%g33
        wrk22 = this%g11*this%g33 - this%g13*this%g31 
        wrk32 = this%g12*this%g31 - this%g11*this%g32
        wrk13 = this%g12*this%g23 - this%g13*this%g22
        wrk23 = this%g13*this%g21 - this%g11*this%g23
        wrk33 = this%g11*this%g22 - this%g12*this%g21
        !det(g_e)
        detg = wrk13*this%g31 + wrk23*this%g32 + wrk33*this%g33

        if(this%strainHard) then        
           if(this%use_gTg) then
              !G_p_ij = G_kl *(g_e)^-1_lj *(g_e)^-1_ki  
              gp11 =  (  (this%gt11*wrk11 + this%gt12*wrk21  + this%gt13*wrk31) * wrk11 &
                   + (this%gt21*wrk11 + this%gt22*wrk21  + this%gt23*wrk31) * wrk21 &
                   + (this%gt31*wrk11 + this%gt32*wrk21  + this%gt33*wrk31) * wrk31 ) / detg**2
              gp12 =  (  (this%gt11*wrk12 + this%gt12*wrk22  + this%gt13*wrk32) * wrk11 &
                   + (this%gt21*wrk12 + this%gt22*wrk22  + this%gt23*wrk32) * wrk21 &
                   + (this%gt31*wrk12 + this%gt32*wrk22  + this%gt33*wrk32) * wrk31 ) / detg**2
              gp13 =  (  (this%gt11*wrk13 + this%gt12*wrk23  + this%gt13*wrk33) * wrk11 &
                   + (this%gt21*wrk13 + this%gt22*wrk23  + this%gt23*wrk33) * wrk21 &
                   + (this%gt31*wrk13 + this%gt32*wrk23  + this%gt33*wrk33) * wrk31 ) / detg**2

              gp21 =  (  (this%gt11*wrk11 + this%gt12*wrk21  + this%gt13*wrk31) * wrk12 &
                   + (this%gt21*wrk11 + this%gt22*wrk21  + this%gt23*wrk31) * wrk22 &
                   + (this%gt31*wrk11 + this%gt32*wrk21  + this%gt33*wrk31) * wrk32 ) / detg**2
              gp22 =  (  (this%gt11*wrk12 + this%gt12*wrk22  + this%gt13*wrk32) * wrk12 &
                   + (this%gt21*wrk12 + this%gt22*wrk22  + this%gt23*wrk32) * wrk22 &
                   + (this%gt31*wrk12 + this%gt32*wrk22  + this%gt33*wrk32) * wrk32 ) / detg**2
              gp23 =  (  (this%gt11*wrk13 + this%gt12*wrk23  + this%gt13*wrk33) * wrk12 &
                   + (this%gt21*wrk13 + this%gt22*wrk23  + this%gt23*wrk33) * wrk22 &
                   + (this%gt31*wrk13 + this%gt32*wrk23  + this%gt33*wrk33) * wrk32 ) / detg**2

              gp31 =  (  (this%gt11*wrk11 + this%gt12*wrk21  + this%gt13*wrk31) * wrk13 &
                   + (this%gt21*wrk11 + this%gt22*wrk21  + this%gt23*wrk31) * wrk23 &
                   + (this%gt31*wrk11 + this%gt32*wrk21  + this%gt33*wrk31) * wrk33 ) / detg**2
              gp32 =  (  (this%gt11*wrk12 + this%gt12*wrk22  + this%gt13*wrk32) * wrk13 &
                   + (this%gt21*wrk12 + this%gt22*wrk22  + this%gt23*wrk32) * wrk23 &
                   + (this%gt31*wrk12 + this%gt32*wrk22  + this%gt33*wrk32) * wrk33 ) / detg**2
              gp33 =  (  (this%gt11*wrk13 + this%gt12*wrk23  + this%gt13*wrk33) * wrk13 &
                   + (this%gt21*wrk13 + this%gt22*wrk23  + this%gt23*wrk33) * wrk23 &
                   + (this%gt31*wrk13 + this%gt32*wrk23  + this%gt33*wrk33) * wrk33 ) / detg**2

           else
              !g_p_ij = g_ik *(g_e)^-1_kj
              gp11 = (this%gt11*wrk11 + this%gt12*wrk21 + this%gt13*wrk31) /detg
              gp12 = (this%gt11*wrk12 + this%gt12*wrk22 + this%gt13*wrk32) /detg
              gp13 = (this%gt11*wrk13 + this%gt12*wrk23 + this%gt13*wrk33) /detg
              gp21 = (this%gt21*wrk11 + this%gt22*wrk21 + this%gt23*wrk31) /detg
              gp22 = (this%gt21*wrk12 + this%gt22*wrk22 + this%gt23*wrk32) /detg
              gp23 = (this%gt21*wrk13 + this%gt22*wrk23 + this%gt23*wrk33) /detg
              gp31 = (this%gt31*wrk11 + this%gt32*wrk21 + this%gt33*wrk31) /detg
              gp32 = (this%gt31*wrk12 + this%gt32*wrk22 + this%gt33*wrk32) /detg
              gp33 = (this%gt31*wrk13 + this%gt32*wrk23 + this%gt33*wrk33) /detg
           end if

           !Eulerian-Almansi strain tensor: ea = (I-(g_p)^T.g_p)/2 --- not explicitly calculated
           !strain norm: e_p = sqrt(2/3*ea_ij*ea_ij)
           if (this%use_gTg) then
              this%e_p = sqrt(( (1.0d0 - gp11)**2 + gp12**2 + gp13**2 + gp21**2 + (1.0d0 - gp22)**2 + gp23**2 + gp31**2 + gp32**2 + (1.0d0 - gp33)**2 )/6.0d0)
           else
              this%e_p = sqrt( ( (1.0d0 - gp11**2 - gp21**2 - gp31**2)**2 + (1.0d0 - gp12**2 - gp22**2 - gp32**2)**2 + (1.0d0 - gp13**2 - gp23**2 - gp33**2)**2 )/6.0d0 + ( (-gp11*gp12 - gp21*gp22 - gp31*gp32)**2 + (-gp11*gp13 - gp21*gp23 - gp31*gp33)**2 + (-gp12*gp13 - gp22*gp23 - gp32*gp33)**2 )/3.0d0 )
           endif
           ! print*, maxval(this%e_p)

           if (this%use_gTg) then
              this%e_pp = sqrt(( (1.0d0 - this%gp11)**2 + this%gp12**2 + this%gp13**2 + this%gp21**2 + (1.0d0 - this%gp22)**2 + this%gp23**2 + this%gp31**2 + this%gp32**2 + (1.0d0 - this%gp33)**2 )/6.0d0)
           else
              this%e_pp = sqrt( ( (1.0d0 - this%gp11**2 - this%gp21**2 - this%gp31**2)**2 + (1.0d0 - this%gp12**2 - this%gp22**2 - this%gp32**2)**2 + (1.0d0 - this%gp13**2 - this%gp23**2 - this%gp33**2)**2 )/6.0d0 + ( (-this%gp11*this%gp12 - this%gp21*this%gp22 - this%gp31*this%gp32)**2 + (-this%gp11*this%gp13 - this%gp21*this%gp23 - this%gp31*this%gp33)**2 + (-this%gp12*this%gp13 - this%gp22*this%gp23 - this%gp32*this%gp33)**2 )/3.0d0 )
           endif

        else
           this%e_p  = 0.0
           this%e_pp = 0.0

           if(this%use_gTg) detg = sqrt(detg)
        endif

        ! !difference between g_p and g_t predictions for e_p
        ! if (this%use_gTg) then
        !    this%e_pdiff = sqrt(( (1.0d0 - this%gp11)**2 + this%gp12**2 + this%gp13**2 + this%gp21**2 + (1.0d0 - this%gp22)**2 + this%gp23**2 + this%gp31**2 + this%gp32**2 + (1.0d0 - this%gp33)**2 )/6.0d0) - this%e_p
        ! else
        !    this%e_pdiff = sqrt( ( (1.0d0 - this%gp11**2 - this%gp21**2 - this%gp31**2)**2 + (1.0d0 - this%gp12**2 - this%gp22**2 - this%gp32**2)**2 + (1.0d0 - this%gp13**2 - this%gp23**2 - this%gp33**2)**2 )/6.0d0 + ( (-this%gp11*this%gp12 - this%gp21*this%gp22 - this%gp31*this%gp32)**2 + (-this%gp11*this%gp13 - this%gp21*this%gp23 - this%gp31*this%gp33)**2 + (-this%gp12*this%gp13 - this%gp22*this%gp23 - this%gp32*this%gp33)**2 )/3.0d0 ) - this%e_p
        ! endif
        ! ! print*, maxval(this%e_pdiff)


        ! !Kospall model from LANL PAGOSA manual
        ! this%elastic%kos_b     = 0.0d0
        ! this%elastic%kos_t     = 3.0d0 !zeroed later
        ! this%elastic%kos_h     = 0.0d0
        ! this%elastic%kos_g     = 0.0d0   
        ! this%elastic%kos_m     = 0.0d0
        ! this%elastic%kos_q     = 1.0d0
        ! this%elastic%kos_f     = this%elastic%kos_g
        ! this%elastic%kos_alpha = 1.0d0
        ! this%elastic%kos_beta  = 0.20d0 !0.20 aluminum, 0.54 copper
        ! this%elastic%kos_e     = 0.0d0

        et = this%eh + this%eel

        if(minval(detg).le.1.0e-10) then
           print*,'detg = ',minval(detg)
        endif

        !modified shear modulus
        pmu = 1.0d0 + this%elastic%kos_b*this%p*detg**(-1.0/3.0)
        fmu = (this%T-this%elastic%kos_t)*this%elastic%kos_h
        hmu = exp(-this%elastic%kos_g*et/max(this%elastic%kos_m-et,1.0d-8))
        !comment for uncomment later for full model
        this%elastic%mu = this%elastic%mu0*max(pmu-fmu,0.0d0)*hmu

        !modified yield modulus
        pyield = 1.0d0 + this%elastic%kos_b*this%elastic%kos_q*this%p*detg**(-1.0/3.0)
        fyield = (this%T-this%elastic%kos_t)*this%elastic%kos_h
        hyield = exp(-this%elastic%kos_f*et/max(this%elastic%kos_m-et,1.0d-8))
        !strain hardening choice
        if (this%elastic%kos_sh.eq.1) then
           this%elastic%yield = this%elastic%yield0*(1.0d0+this%elastic%kos_alpha*(this%elastic%kos_e+this%e_p))**this%elastic%kos_beta*max(pyield-fyield,0.0d0)*hyield
           !print*,this%elastic%kos_sh,this%elastic%kos_alpha,this%elastic%kos_beta,maxval(hyield),minval(hyield),this%elastic%yield0
           !stop
           !print*,maxval(this%elastic%yield),maxval(this%e_p)
        elseif (this%elastic%kos_sh.eq.2) then
           this%elastic%yield = this%elastic%yield0*(1.0d0+this%elastic%kos_alpha*(this%elastic%kos_e+this%e_pp))**this%elastic%kos_beta*max(pyield-fyield,0.0d0)*hyield
        elseif (this%elastic%kos_sh.eq.2) then
           this%elastic%yield = this%elastic%yield0*(1.0d0+this%elastic%kos_alpha*(this%elastic%kos_e+this%pe))**this%elastic%kos_beta*max(pyield-fyield,0.0d0)*hyield
        endif
        !multiply Johnson-cook factor for dependence on strain rate (1+c*log(e_p'))




        !simple thermal softening -- set melt_c to zero in input to turn off -- otherwise it modifes above predicted yield and mu to account for temperature based thermal softening
        fm = exp(-this%elastic%melt_c*this%T/max(this%elastic%melt_t-this%T,1.0d-8))
        !fm = 1.0
        !this%elastic%yield = this%elastic%yield0*fm
        !this%elastic%mu = this%elastic%mu0*fm
        this%elastic%yield = this%elastic%yield*fm
        this%elastic%mu = this%elastic%mu*fm

    end subroutine getModifiedModulii


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
           print*,"explicit plastic needs fix for gTg"
           stop
            invtaurel = two*invtaurel  ! Factor of 2 for gTg implementation
        end if

        ! Add (1/tau_rel)*g_e*S to the rhsg (explicit plastic source terms)
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
        use decomp_2d,  only: nrank
        use RKCoeffs,   only: RK45_A,RK45_B
        class(solid), intent(inout) :: this
        integer, intent(in) :: isub
        real(rkind), intent(in) :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho,u,v,w
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhsYs  ! RHS for mass fraction equation

        if(this%intSharp) then        
           call this%getRHS_Ys_intSharp(rho,u,v,w,rhsYs,x_bc,y_bc,z_bc)
        else
           call this%getRHS_Ys(rho,u,v,w,rhsYs,x_bc,y_bc,z_bc)
        endif
        call hook_material_mass_source(this%decomp,this%hydro,this%elastic,x,y,z,tsim,rho,u,v,w,this%Ys,this%VF,this%p,rhsYs)

        ! advance sub-step
        if(isub==1) this%QtmpYs = zero                   ! not really needed, since RK45_A(1) = 0
        this%QtmpYs  = dt*rhsYs + RK45_A(isub)*this%QtmpYs
        this%consrv(:,:,:,1) = this%consrv(:,:,:,1)  + RK45_B(isub)*this%QtmpYs
!print *, 'rhs Ys:', rhsYs(89,1,1), rho(89,1,1), this%Ys(89,1,1)
!print *, 'Ji:    ', this%Ji(89,1,1,1), this%Ji(89,1,1,2), this%Ji(89,1,1,3)
!print *, 'cns Ys:', this%consrv(89,1,1,1)
    end subroutine

    subroutine getRHS_Ys_intSharp(this,rho,u,v,w,rhsYs,x_bc,y_bc,z_bc)
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
        
        !original terms
        call divergence(this%decomp,this%der,tmp1,tmp2,tmp3,rhsYs,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric
        
        if(this%intSharp_spf) then
           rhsYs = rhsYs + this%intSharp_R(:,:,:,1) !ignore components 2 and 3 when not in divergence form

           !high order VF bounds diffusion terms
           call divergence(this%decomp,this%der,this%intSharp_RDiff(:,:,:,1),this%intSharp_RDiff(:,:,:,2),this%intSharp_RDiff(:,:,:,3),tmp1,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric
           rhsYs = rhsYs + tmp1

        else
           !low order terms
           call divergence(this%decomp,this%derD02,this%intSharp_R(:,:,:,1),this%intSharp_R(:,:,:,2),this%intSharp_R(:,:,:,3),tmp1,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric
           rhsYs = rhsYs + tmp1
           
           !high order terms
           call divergence(this%decomp,this%der,this%intSharp_RDiff(:,:,:,1),this%intSharp_RDiff(:,:,:,2),this%intSharp_RDiff(:,:,:,3),tmp1,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric
           rhsYs = rhsYs + tmp1
           
           !FV terms
           rhsYs = rhsYs + this%intSharp_RFV
           
        endif


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

        ! ! filter rg
        ! call filter3D(this%decomp, fil_, this%rg11, iflag, x_bc, y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rg12, iflag,-x_bc,-y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rg13, iflag,-x_bc, y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rg21, iflag,-x_bc,-y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rg22, iflag, x_bc, y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rg23, iflag, x_bc,-y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rg31, iflag,-x_bc, y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rg32, iflag, x_bc,-y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rg33, iflag, x_bc, y_bc, z_bc)

        ! ! filter rgt
        ! call filter3D(this%decomp, fil_, this%rgt11, iflag, x_bc, y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rgt12, iflag,-x_bc,-y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rgt13, iflag,-x_bc, y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rgt21, iflag,-x_bc,-y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rgt22, iflag, x_bc, y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rgt23, iflag, x_bc,-y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rgt31, iflag,-x_bc, y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rgt32, iflag, x_bc,-y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rgt33, iflag, x_bc, y_bc, z_bc)

        ! ! filter rgp
        ! call filter3D(this%decomp, fil_, this%rgp11, iflag, x_bc, y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rgp12, iflag,-x_bc,-y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rgp13, iflag,-x_bc, y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rgp21, iflag,-x_bc,-y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rgp22, iflag, x_bc, y_bc, z_bc)
        ! call filter3D(this%decomp, fil_, this%rgp23, iflag, x_bc,-y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rgp31, iflag,-x_bc, y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rgp32, iflag, x_bc,-y_bc,-z_bc)
        ! call filter3D(this%decomp, fil_, this%rgp33, iflag, x_bc, y_bc, z_bc)

        ! !filter rpe
        ! call filter3D(this%decomp, fil_, this%rpe, iflag, x_bc, y_bc, z_bc)

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


    subroutine filter_g(this, iflag, x_bc, y_bc, z_bc, fil)
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

        ! filter gt
        call filter3D(this%decomp, fil_, this%gt11, iflag, x_bc, y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%gt12, iflag,-x_bc,-y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%gt13, iflag,-x_bc, y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%gt21, iflag,-x_bc,-y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%gt22, iflag, x_bc, y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%gt23, iflag, x_bc,-y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%gt31, iflag,-x_bc, y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%gt32, iflag, x_bc,-y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%gt33, iflag, x_bc, y_bc, z_bc)

        ! filter gp
        call filter3D(this%decomp, fil_, this%gp11, iflag, x_bc, y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%gp12, iflag,-x_bc,-y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%gp13, iflag,-x_bc, y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%gp21, iflag,-x_bc,-y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%gp22, iflag, x_bc, y_bc, z_bc)
        call filter3D(this%decomp, fil_, this%gp23, iflag, x_bc,-y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%gp31, iflag,-x_bc, y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%gp32, iflag, x_bc,-y_bc,-z_bc)
        call filter3D(this%decomp, fil_, this%gp33, iflag, x_bc, y_bc, z_bc)

        !filter pe
        call filter3D(this%decomp, fil_, this%pe, iflag, x_bc, y_bc, z_bc)

    end subroutine



    subroutine getSpeciesDensity(this,rho,rhom)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rhom

        ! Get detg in rhom
        rhom = this%g11*(this%g22*this%g33-this%g23*this%g32) &
             - this%g12*(this%g21*this%g33-this%g31*this%g23) &
             + this%g13*(this%g21*this%g32-this%g31*this%g22)

        if (this%use_gTg.and.(.not.this%strainHard)) then
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

        !call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG,this%use_gTg)
        call this%elastic%get_finger(this%g,finger,fingersq,trG,trG2,detG,((this%use_gTg).and.(.not.this%strainHard)))

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


    subroutine get_primitive_g(this,rho)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho
        integer:: i

        if (this%cnsrv_g) then
           do i=1,9
              this%g(:,:,:,i) = this%rg(:,:,:,i) / rho
           enddo
        else
           this%g = this%rg
        endif
        if (this%cnsrv_gt) then
           do i=1,9
              this%g_t(:,:,:,i) = this%rg_t(:,:,:,i) / rho
           enddo
        else
           this%g_t = this%rg_t
        endif
        if (this%cnsrv_gp) then
           do i=1,9
              this%g_p(:,:,:,i) = this%rg_p(:,:,:,i) / rho
           enddo
        else
           this%g_p = this%rg_p
        endif
        if (this%cnsrv_pe) then
           this%pe = this%rpe / rho
        else
           this%pe = this%rpe
        endif

    end subroutine

    subroutine get_conserved_g(this,rho)
        class(solid), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho
        integer:: i

        if (this%cnsrv_g) then
           do i=1,9
              this%rg(:,:,:,i) = this%g(:,:,:,i) * rho
           enddo
        else
           this%rg = this%g
        endif
        if (this%cnsrv_gt) then
           do i=1,9
              this%rg_t(:,:,:,i) = this%g_t(:,:,:,i) * rho
           enddo
        else
           this%rg_t = this%g_t
        endif
        if (this%cnsrv_gp) then
           do i=1,9
              this%rg_p(:,:,:,i) = this%g_p(:,:,:,i) * rho
           enddo
        else
           this%rg_p = this%g_p
        endif
        if (this%cnsrv_pe) then
           this%rpe = this%pe * rho
        else
           this%rpe = this%pe
        endif

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
            print*,'g NaN',i+this%decomp%yst(1)-1,j+this%decomp%yst(2)-1,k+this%decomp%yst(3)-1,l
            call GracefulExit(charout,4809)
        end if
        if ( nancheck(this%g_t,i,j,k,l) ) then
           write(charout,'(A,5(I0,A))') "NaN encountered in material ", imat,&
                &" at (",i,", ",j,", ",k,", ",l,") of g_t"
           print*,'g_t NaN',i+this%decomp%yst(1)-1,j+this%decomp%yst(2)-1,k+this%decomp%yst(3)-1,l
           call GracefulExit(charout,4809)
        end if
        if ( nancheck(this%g_p,i,j,k,l) ) then
           write(charout,'(A,5(I0,A))') "NaN encountered in material ", imat,&
                &" at (",i,", ",j,", ",k,", ",l,") of g_p"
           print*,'g_p NaN',i+this%decomp%yst(1)-1,j+this%decomp%yst(2)-1,k+this%decomp%yst(3)-1,l
           call GracefulExit(charout,4809)
        end if
        if ( nancheck(this%pe,i,j,k) ) then
           write(charout,'(A,5(I0,A))') "NaN encountered in material ", imat,&
                &" at (",i,", ",j,", ",k,") of pe"
           print*,'pe NaN',i,j,k
           call GracefulExit(charout,4809)
        end if
        if ( nancheck(this%VF,i,j,k) ) then
            write(charout,'(A,4(I0,A))') "NaN encountered in material ", imat,&
                         &" at (",i,", ",j,", ",k,") of VF"
            print*,'VF NaN',i+this%decomp%yst(1)-1,j+this%decomp%yst(2)-1,k+this%decomp%yst(3)-1
            call GracefulExit(charout,4809)
        end if
        if ( nancheck(this%consrv,i,j,k,l) ) then
            write(charout,'(A,5(I0,A))') "NaN encountered in material ", imat,&
                         &" at (",i,", ",j,", ",k,", ",l,") of consrv"
            print*,'consrv NaN',i+this%decomp%yst(1)-1,j+this%decomp%yst(2)-1,k+this%decomp%yst(3)-1,l
            call GracefulExit(charout,4809)
        end if

    end subroutine

    ! subroutine sliding_deformation(this, normal, mask)
    !     use kind_parameters, only: clen
    !     use constants,       only: zero, eps, third, twothird, pi
    !     use decomp_2d,       only: nrank
    !     use operators,       only: filter3D
    !     use exits,           only: GracefulExit, nancheck, message
    !     class(solid),                                         intent(inout) :: this
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp,3), intent(in)    :: normal ! Interface normal for sliding treatment
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)    :: mask   ! Mask to avoid sliding away from interface
    !     real(rkind) :: yield = zero

    !     real(rkind), dimension(3,3) :: g, u, vt, gradf, gradf_new, sigma, sigma_tilde, v_tilde
    !     real(rkind), dimension(3)   :: sval, beta, Sa, f, f1, f2, dbeta, beta_new, dbeta_new
    !     real(rkind) :: sqrt_om, betasum, Sabymu_sq, ycrit, C0, t
    !     real(rkind) :: tol = real(1.D-12,rkind), residual, residual_new
    !     integer :: i,j,k
    !     integer, dimension(1) :: min_norm
    !     integer :: iters
    !     integer, parameter :: niters = 500
    !     integer :: lwork, info
    !     integer, dimension(3) :: ipiv
    !     integer :: nxp, nyp, nzp
    !     character(len=clen) :: charout

    !     real(rkind), dimension(:), allocatable :: svdwork     ! Work array for SVD stuff

    !     real(rkind) :: dx, dy, x, y, rad, theta = 45._rkind * pi / 180._rkind

    !     print*,"g_t not included for sliding deformation: need to add"
    !     stop

    !     ! print *, "In sliding_deformation"

    !     ! where ((this%VF > 0.01) .and. (this%VF < 0.99))
    !     !     mask = one
    !     ! elsewhere
    !     !     mask = zero
    !     ! end where
    !     ! call filter3D(this%decomp, this%gfil, mask, 1)
    !     ! mask = this%VF*(one - this%VF)
    !     ! mask = mask/maxval(mask)
    !     ! mask = one - (one - mask)**mask_exponent

    !     allocate(svdwork(1))
    !     if (this%use_gTg) then
    !         ! Get optimal lwork
    !         lwork = -1
    !         call dsyev('V', 'U', 3, G, 3, sval, svdwork, lwork, info)
    !         lwork = svdwork(1)
    !     else
    !         ! Get optimal lwork
    !         lwork = -1
    !         call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, svdwork, lwork, info)
    !         lwork = svdwork(1)
    !     end if

    !     if (lwork .GT. size(svdwork)) then
    !         deallocate(svdwork); allocate(svdwork(lwork))
    !     end if

    !     if ( nancheck(this%g) ) then
    !         call message("NaN found in g during plastic relaxation.")
    !     end if

    !     dx = one/real(this%decomp%xsz(1)-1,rkind)
    !     dy = dx

    !     do k = 1,this%nzp
    !         do j = 1,this%nyp
    !             do i = 1,this%nxp

    !                 if (abs(mask(i,j,k)) > real(1.0D-2,rkind)) then
    !                     ! yield = zero
    !                     !yield = (one-mask(i,j,k))*this%elastic%yield
    !                     yield = (one-mask(i,j,k))*this%elastic%yield(i,j,k) !mca

    !                     g(1,1) = this%g(i,j,k,1); g(1,2) = this%g(i,j,k,2); g(1,3) = this%g(i,j,k,3)
    !                     g(2,1) = this%g(i,j,k,4); g(2,2) = this%g(i,j,k,5); g(2,3) = this%g(i,j,k,6)
    !                     g(3,1) = this%g(i,j,k,7); g(3,2) = this%g(i,j,k,8); g(3,3) = this%g(i,j,k,9)

    !                     !!!! HACK !!!
    !                     ! Hard code v_tilde for now as the cartesian basis
    !                     ! Should eventually set v_tilde(:,1) to the interface normal and the rest as orthogonal directions

    !                     ! 1D
    !                     ! v_tilde(1,1) = one;  v_tilde(1,2) = zero; v_tilde(1,3) = zero
    !                     ! v_tilde(2,1) = zero; v_tilde(2,2) = one;  v_tilde(2,3) = zero
    !                     ! v_tilde(3,1) = zero; v_tilde(3,2) = zero; v_tilde(3,3) = one

    !                     ! 2D oblique
    !                     ! v_tilde(1,1) = cos(theta); v_tilde(1,2) = sin(theta); v_tilde(1,3) = zero
    !                     ! v_tilde(2,1) =-sin(theta); v_tilde(2,2) = cos(theta); v_tilde(2,3) = zero
    !                     ! v_tilde(3,1) = zero;       v_tilde(3,2) = zero;       v_tilde(3,3) = one

    !                     ! 2D circular
    !                     ! x = - half + real( this%decomp%yst(1) - 1 + i - 1, rkind ) * dx
    !                     ! y = - half + real( this%decomp%yst(2) - 1 + j - 1, rkind ) * dy
    !                     ! rad = sqrt(x**2 + y**2)
    !                     ! theta = atan2(y,x)
    !                     ! this%kap(i,j,k) = theta
    !                     ! v_tilde(1,1) = cos(theta); v_tilde(1,2) =-sin(theta); v_tilde(1,3) = zero
    !                     ! v_tilde(2,1) = sin(theta); v_tilde(2,2) = cos(theta); v_tilde(2,3) = zero
    !                     ! v_tilde(3,1) = zero;       v_tilde(3,2) = zero;       v_tilde(3,3) = one

    !                     ! General interface normal
    !                     v_tilde(:,1) = normal(i,j,k,:) ! Assume normalized

    !                     ! Get a vector not parallel to v_tilde(:,1)
    !                     min_norm = minloc(abs(normal(i,j,k,:)))
    !                     v_tilde(:,3) = zero; v_tilde(min_norm(1),3) = one 
    !                     ! if ((i == 44) .and. (j == 2)) then
    !                     !     print *, "normal = ", normal(i,j,k,:)
    !                     !     print *, "2nd vector = ", v_tilde(:,3)
    !                     ! end if

    !                     ! Get an othogonal vector to v_tilde(:,1) through a crossproduct with a non-parallel vector
    !                     v_tilde(1,2) = v_tilde(2,1)*v_tilde(3,3) - v_tilde(3,1)*v_tilde(2,3)
    !                     v_tilde(2,2) = v_tilde(3,1)*v_tilde(1,3) - v_tilde(1,1)*v_tilde(3,3)
    !                     v_tilde(3,2) = v_tilde(1,1)*v_tilde(2,3) - v_tilde(2,1)*v_tilde(1,3)
    !                     v_tilde(:,2) = v_tilde(:,2)/sqrt(v_tilde(1,2)**2 + v_tilde(2,2)**2 + v_tilde(3,2)**2) ! Normalize

    !                     ! Get third orthogonal vector through a crossproduct of the first two
    !                     v_tilde(1,3) = v_tilde(2,1)*v_tilde(3,2) - v_tilde(3,1)*v_tilde(2,2)
    !                     v_tilde(2,3) = v_tilde(3,1)*v_tilde(1,2) - v_tilde(1,1)*v_tilde(3,2)
    !                     v_tilde(3,3) = v_tilde(1,1)*v_tilde(2,2) - v_tilde(2,1)*v_tilde(1,2)
    !                     v_tilde(:,3) = v_tilde(:,3)/sqrt(v_tilde(1,3)**2 + v_tilde(2,3)**2 + v_tilde(3,3)**2) ! Normalize

    !                     ! if ((i == 44) .and. (j == 2)) then
    !                     !     print *, "v_tilde:"
    !                     !     print *, "   ",  v_tilde(1,1), v_tilde(1,2), v_tilde(1,3)
    !                     !     print *, "   ",  v_tilde(2,1), v_tilde(2,2), v_tilde(2,3)
    !                     !     print *, "   ",  v_tilde(3,1), v_tilde(3,2), v_tilde(3,3)
    !                     ! end if

    !                     if (this%use_gTg) then
    !                         ! Get eigenvalues and eigenvectors of G
    !                         call dsyev('V', 'U', 3, G, 3, sval, svdwork, lwork, info)
    !                         if(info .ne. 0) print '(A,I6,A)', 'proc ', nrank, ': Problem with DSYEV. Please check.'
    !                         sval = sqrt(sval)  ! Get singular values of g
    !                     else
    !                         ! Get SVD of g
    !                         call dgesvd('A', 'A', 3, 3, g, 3, sval, u, 3, vt, 3, svdwork, lwork, info)
    !                         if(info .ne. 0) then
    !                             write(charout, '(A,I0,A)') 'proc ', nrank, ': Problem with SVD. Please check.'
    !                             call GracefulExit(charout,3475)
    !                         end if
    !                     end if

    !                     sqrt_om = sval(1)*sval(2)*sval(3)
    !                     beta = sval**two / sqrt_om**(two/three)

    !                     betasum = sum( beta*(beta-one) ) / three
    !                     !Sa = -this%elastic%mu*sqrt_om * ( beta*(beta-one) - betasum )
    !                     Sa = -this%elastic%mu(i,j,k)*sqrt_om * ( beta*(beta-one) - betasum ) !mca

    !                     ! Get the actual stress
    !                     sigma = zero
    !                     sigma(1,1) = Sa(1); sigma(2,2) = Sa(2); sigma(3,3) = Sa(3)
    !                     ! sigma = matmul(sigma, transpose(vt))
    !                     ! sigma = matmul(vt, sigma)
    !                     sigma = matmul(sigma, vt)
    !                     sigma = matmul(transpose(vt), sigma)

    !                     ! if ((i == 44) .and. (j == 2)) then
    !                     !     print *, "V:"
    !                     !     print *, "   ",  vt(1,1), vt(2,1), vt(3,1)
    !                     !     print *, "   ",  vt(1,2), vt(2,2), vt(3,2)
    !                     !     print *, "   ",  vt(1,3), vt(2,3), vt(3,3)
    !                     ! end if

    !                     ! if ((i == 44) .and. (j == 2)) then
    !                     !     print *, "sigma:"
    !                     !     print *, "   ",  sigma(1,1), sigma(1,2), sigma(1,3)
    !                     !     print *, "   ",  sigma(2,1), sigma(2,2), sigma(2,3)
    !                     !     print *, "   ",  sigma(3,1), sigma(3,2), sigma(3,3)
    !                     ! end if

    !                     ! Get new stress
    !                     sigma_tilde = zero
    !                     sigma_tilde(1,1) = sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,1) ) )
    !                     sigma_tilde(2,2) = sum( v_tilde(:,2) * matmul( sigma, v_tilde(:,2) ) )
    !                     sigma_tilde(2,3) = sum( v_tilde(:,3) * matmul( sigma, v_tilde(:,2) ) ) ! Need this off-diagonal term
    !                                                                                              ! since v_tilde(:,2:3) are not
    !                                                                                              ! the singular vectors but just
    !                                                                                              ! the basis for the interface
    !                                                                                              ! tangential singular vectors
    !                     sigma_tilde(3,2) = sigma_tilde(2,3) ! Symmetric
    !                     sigma_tilde(3,3) = sum( v_tilde(:,3) * matmul( sigma, v_tilde(:,3) ) )

    !                     ! sigma_tilde = matmul(sigma_tilde, v_tilde)
    !                     ! sigma_tilde = matmul(transpose(v_tilde), sigma_tilde)
    !                     sigma_tilde = matmul(sigma_tilde, transpose(v_tilde))
    !                     sigma_tilde = matmul(v_tilde, sigma_tilde)

    !                     ! if ((i == 21) .and. (j == 64)) then
    !                     !     print *, "theta = ", theta
    !                     !     print *, "sigma_tilde:"
    !                     !     print *, "   ",  sigma_tilde(1,1), sigma_tilde(1,2), sigma_tilde(1,3)
    !                     !     print *, "   ",  sigma_tilde(2,1), sigma_tilde(2,2), sigma_tilde(2,3)
    !                     !     print *, "   ",  sigma_tilde(3,1), sigma_tilde(3,2), sigma_tilde(3,3)
    !                     !     print *, "tangential: ", sum( v_tilde(:,1) * matmul( sigma_tilde, v_tilde(:,2) ) ), &
    !                     !                              sum( v_tilde(:,1) * matmul( sigma_tilde, v_tilde(:,3) ) )
    !                     ! end if

    !                     ! Try to make it a smoother transition to sliding
    !                     ! sigma_tilde = mask(i,j,k)*sigma_tilde + (one - mask(i,j,k))*sigma

    !                     ! Get eigenvalues and eigenvectors of sigma
    !                     call dsyev('V', 'U', 3, sigma_tilde, 3, Sa, svdwork, lwork, info)
    !                     if(info .ne. 0) then
    !                         print '(A,I6,A)', 'proc ', nrank, ': Problem with DSYEV. Please check.'
    !                         print *, "sigma:"
    !                         print *, "   ",  sigma(1,1), sigma(1,2), sigma(1,3)
    !                         print *, "   ",  sigma(2,1), sigma(2,2), sigma(2,3)
    !                         print *, "   ",  sigma(3,1), sigma(3,2), sigma(3,3)
    !                         print *, "tangential: ", sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,2) ) ), &
    !                                                  sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,3) ) )
    !                     end if

    !                     ! if ((i == 108) .and. (j == 2)) then
    !                     !     print *, "v_tilde_new:"
    !                     !     print *, "   ",  sigma_tilde(1,1), sigma_tilde(1,2), sigma_tilde(1,3)
    !                     !     print *, "   ",  sigma_tilde(2,1), sigma_tilde(2,2), sigma_tilde(2,3)
    !                     !     print *, "   ",  sigma_tilde(3,1), sigma_tilde(3,2), sigma_tilde(3,3)
    !                     ! end if

    !                     ! Now get new beta
    !                     !f = Sa / (this%elastic%mu*sqrt_om); f(3) = beta(1)*beta(2)*beta(3)     ! New function value (target to attain)
    !                     f = Sa / (this%elastic%mu(i,j,k)*sqrt_om); f(3) = beta(1)*beta(2)*beta(3)     ! New function value (target to attain) !mca
                        
    !                     betasum = sum( beta*(beta-one) ) / three
    !                     f1 = -( beta*(beta-one) - betasum ); f1(3) = beta(1)*beta(2)*beta(3)   ! Original function value

    !                     ! Get newton step
    !                     gradf(1,1) = -twothird*(two*beta(1)-one); gradf(1,2) =     third*(two*beta(2)-one); gradf(1,3) = third*(two*beta(3)-one)
    !                     gradf(2,1) =     third*(two*beta(1)-one); gradf(2,2) = -twothird*(two*beta(2)-one); gradf(2,3) = third*(two*beta(3)-one)
    !                     gradf(3,1) = beta(2)*beta(3);             gradf(3,2) = beta(3)*beta(1);             gradf(3,3) = beta(1)*beta(2)

    !                     dbeta = (f-f1)
    !                     call dgesv(3, 1, gradf, 3, ipiv, dbeta, 3, info)
                   
    !                     ! Compute residual
    !                     residual = -sum( (f1-f)*dbeta )                                    ! lambda**2
    !                     iters = 0
    !                     t = 1._rkind
    !                     do while ( (iters < niters) .AND. (abs(residual) .GT. tol) )
    !                         ! Backtracking line search
    !                         t = 1._rkind

    !                         !original
    !                         beta_new = beta + t * dbeta
    !                         !modified
    !                         !beta_new = beta + t * dbeta * exp( (1.0-real(iters,rkind))/real(niters,rkind)*30.0 )

    !                         ! Get new residual
    !                         gradf_new(1,1) = -twothird*(two*beta_new(1)-one);
    !                         gradf_new(1,2) =     third*(two*beta_new(2)-one);
    !                         gradf_new(1,3) = third*(two*beta_new(3)-one)

    !                         gradf_new(2,1) =     third*(two*beta_new(1)-one);
    !                         gradf_new(2,2) = -twothird*(two*beta_new(2)-one);
    !                         gradf_new(2,3) = third*(two*beta_new(3)-one)

    !                         gradf_new(3,1) = beta_new(2)*beta_new(3);
    !                         gradf_new(3,2) = beta_new(3)*beta_new(1);
    !                         gradf_new(3,3) = beta_new(1)*beta_new(2)

    !                         betasum = sum( beta_new*(beta_new-one) ) / three
    !                         f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)
    !                         dbeta_new = (f-f2)
    !                         call dgesv(3, 1, gradf_new, 3, ipiv, dbeta_new, 3, info)
    !                         residual_new = -sum( (f2-f)*dbeta_new )                                    ! lambda**2

    !                         do while ( (abs(residual_new) .GE. abs(residual)) .AND. (t > eps) )
    !                             if (iters .GT. (niters - 10)) then
    !                                 print '(A,I0,3(A,ES15.5))', 'iters = ', iters, ', t = ', t, ', residual_new = ', residual_new, ', residual = ', residual
    !                             end if

    !                             t = half*t
    !                             beta_new = beta + t * dbeta

    !                             gradf_new(1,1) = -twothird*(two*beta_new(1)-one);
    !                             gradf_new(1,2) =     third*(two*beta_new(2)-one);
    !                             gradf_new(1,3) = third*(two*beta_new(3)-one)

    !                             gradf_new(2,1) =     third*(two*beta_new(1)-one);
    !                             gradf_new(2,2) = -twothird*(two*beta_new(2)-one);
    !                             gradf_new(2,3) = third*(two*beta_new(3)-one)

    !                             gradf_new(3,1) = beta_new(2)*beta_new(3);
    !                             gradf_new(3,2) = beta_new(3)*beta_new(1);
    !                             gradf_new(3,3) = beta_new(1)*beta_new(2)

    !                             betasum = sum( beta_new*(beta_new-one) ) / three
    !                             f2 = -( beta_new*(beta_new-one) - betasum ); f2(3) = beta_new(1)*beta_new(2)*beta_new(3)

    !                             dbeta_new = (f-f2)
    !                             call dgesv(3, 1, gradf_new, 3, ipiv, dbeta_new, 3, info)
    !                             residual_new = -sum( (f2-f)*dbeta_new )                                    ! lambda**2
    !                         end do!
    !                         beta = beta_new
    !                         f1 = f2
    !                         dbeta = dbeta_new
    !                         residual = residual_new

    !                         iters = iters + 1
    !                         if (t <= eps) then
    !                             print '(A)', 'Newton solve in sliding_deformation did not converge'
    !                             exit
    !                         end if
    !                     end do
    !                     if ((iters >= niters) .OR. (t <= eps)) then
    !                         write(charout,'(4(A,I0))') 'Newton solve in sliding_deformation did not converge at index ',i,',',j,',',k,' of process ',nrank
    !                         print '(A)', charout
    !                         print '(A,3(X,I0))', 'sty = ', this%decomp%yst
    !                         print '(A)', 'g = '
    !                         print '(4X,3(ES15.5))', this%g(i,j,k,1), this%g(i,j,k,2), this%g(i,j,k,3)
    !                         print '(4X,3(ES15.5))', this%g(i,j,k,4), this%g(i,j,k,5), this%g(i,j,k,6)
    !                         print '(4X,3(ES15.5))', this%g(i,j,k,7), this%g(i,j,k,8), this%g(i,j,k,9)
    !                         print '(A,ES15.5)', '( ||S||^2 - (2/3) sigma_Y^2 )/mu^2 = ', ycrit

    !                         print '(A,ES15.5)', 'Relaxation, t = ', t
    !                         print '(A,ES15.5)', 'Residual = ', residual
    !                         call GracefulExit(charout,6382)
    !                     end if

    !                     ! Then get new svals
    !                     sval = sqrt(beta) * sqrt_om**(one/three)

    !                     if (this%use_gTg) then
    !                         sval = sval*sval ! New eigenvalues of G
                            
    !                         ! Get g = v*sval*vt
    !                         u = sigma_tilde; vt = transpose(u)
    !                         vt(1,:) = vt(1,:)*sval(1); vt(2,:) = vt(2,:)*sval(2); vt(3,:) = vt(3,:)*sval(3)  ! eigval*vt
    !                         G = MATMUL(u,vt) ! v*eigval*vt
    !                     else
    !                         ! Get g = u*sval*vt
    !                         u = sigma_tilde; vt = transpose(u)
    !                         vt(1,:) = vt(1,:)*sval(1); vt(2,:) = vt(2,:)*sval(2); vt(3,:) = vt(3,:)*sval(3)  ! sval*vt
    !                         g = MATMUL(u,vt) ! u*sval*vt
    !                     end if


    !                     ! Try to make it a smoother transition to sliding
    !                     this%g(i,j,k,1) = mask(i,j,k)*g(1,1) + (one - mask(i,j,k))*this%g(i,j,k,1)
    !                     this%g(i,j,k,2) = mask(i,j,k)*g(1,2) + (one - mask(i,j,k))*this%g(i,j,k,2)
    !                     this%g(i,j,k,3) = mask(i,j,k)*g(1,3) + (one - mask(i,j,k))*this%g(i,j,k,3)
    !                     this%g(i,j,k,4) = mask(i,j,k)*g(2,1) + (one - mask(i,j,k))*this%g(i,j,k,4)
    !                     this%g(i,j,k,5) = mask(i,j,k)*g(2,2) + (one - mask(i,j,k))*this%g(i,j,k,5)
    !                     this%g(i,j,k,6) = mask(i,j,k)*g(2,3) + (one - mask(i,j,k))*this%g(i,j,k,6)
    !                     this%g(i,j,k,7) = mask(i,j,k)*g(3,1) + (one - mask(i,j,k))*this%g(i,j,k,7)
    !                     this%g(i,j,k,8) = mask(i,j,k)*g(3,2) + (one - mask(i,j,k))*this%g(i,j,k,8)
    !                     this%g(i,j,k,9) = mask(i,j,k)*g(3,3) + (one - mask(i,j,k))*this%g(i,j,k,9)

    !                     ! this%g(i,j,k,1) = g(1,1); this%g(i,j,k,2) = g(1,2); this%g(i,j,k,3) = g(1,3)
    !                     ! this%g(i,j,k,4) = g(2,1); this%g(i,j,k,5) = g(2,2); this%g(i,j,k,6) = g(2,3)
    !                     ! this%g(i,j,k,7) = g(3,1); this%g(i,j,k,8) = g(3,2); this%g(i,j,k,9) = g(3,3)

    !                 end if
    !             end do
    !         end do
    !     end do

    !     ! ! Get devstress for debugging
    !     ! call this%get_eelastic_devstress()
    !     ! sigma(1,1) = this%sxx(21,64,1); sigma(1,2) = this%sxy(21,64,1); sigma(1,3) = this%sxz(21,64,1);
    !     ! sigma(2,1) = this%sxy(21,64,1); sigma(2,2) = this%syy(21,64,1); sigma(2,3) = this%syz(21,64,1);
    !     ! sigma(3,1) = this%sxz(21,64,1); sigma(3,2) = this%syz(21,64,1); sigma(3,3) = this%szz(21,64,1);
    !     ! print *, "New sigma:"
    !     ! print *, "   ",  sigma(1,1), sigma(1,2), sigma(1,3)
    !     ! print *, "   ",  sigma(2,1), sigma(2,2), sigma(2,3)
    !     ! print *, "   ",  sigma(3,1), sigma(3,2), sigma(3,3)

    !     ! ! 2D circular
    !     ! x = - half + real( this%decomp%yst(1) - 1 + 21 - 1, rkind ) * dx
    !     ! y = - half + real( this%decomp%yst(2) - 1 + 64 - 1, rkind ) * dy
    !     ! rad = sqrt(x**2 + y**2)
    !     ! theta = atan2(y,x)

    !     ! v_tilde(1,1) = cos(theta); v_tilde(1,2) =-sin(theta); v_tilde(1,3) = zero
    !     ! v_tilde(2,1) = sin(theta); v_tilde(2,2) = cos(theta); v_tilde(2,3) = zero
    !     ! v_tilde(3,1) = zero;       v_tilde(3,2) = zero;       v_tilde(3,3) = one

    !     ! print *, "tangential: ", sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,2) ) ), &
    !     !                          sum( v_tilde(:,1) * matmul( sigma, v_tilde(:,3) ) )


    !     deallocate(svdwork)
    ! end subroutine

end module
