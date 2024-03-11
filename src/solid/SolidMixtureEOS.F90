module SolidMixtureMod

    use kind_parameters,         only: rkind,clen
    use constants,               only: zero,epssmall,eps,one,two,third,half,four
    use decomp_2d,               only: decomp_info
    use DerivativesMod,          only: derivatives
    use DerivativesStaggeredMod, only: derivativesStagg
    use InterpolatorsMod,        only: interpolators
    use FiltersMod,              only: filters
    use LADMod,                  only: ladobject
    use exits,                   only: GracefulExit
    use EOSMod,                  only: eos
    use StiffGasEOS,             only: stiffgas
    use Sep1SolidEOS,            only: sep1solid
    use SolidMod,                only: solid
    use mpi

    implicit none

    type :: solid_mixture
        
        integer :: mpi_rank, mpi_size, ierror, mpi_rank_prev, mpi_rank_next
        integer :: ns
        integer :: nxp, nyp, nzp
        type(solid), dimension(:), allocatable :: material

        type(decomp_info), pointer      :: decomp
        type(derivatives), pointer      :: der,derD04,derD02, derD06,derCD06
        type(interpolators), pointer    :: interpMid
        type(derivativesStagg), pointer :: derStagg
        type(interpolators), pointer    :: interpMid02
        type(derivativesStagg), pointer :: derStaggd02
        type(filters),     pointer      :: fil
        type(filters),     pointer      :: gfil
        type(ladobject),   pointer      :: LAD

        real(rkind), allocatable, dimension(:,:,:) :: rho0mix, mumix, yieldmix, solidVF

        logical :: SOSmodel = .FALSE., usePhiForm = .TRUE.       ! is sound speed given by `equilibrium' model? Alternative is `frozen' model. Check Saurel et al., JCP 2009.
        logical :: PTeqb = .TRUE., pEqb = .FALSE., pRelax = .FALSE., updateEtot = .FALSE., useAkshayForm = .FALSE., twoPhaseLAD = .FALSE., LAD5eqn = .FALSE.
        logical :: use_gTg = .FALSE., useOneG = .FALSE., intSharp = .TRUE., intSharp_cpl = .TRUE., intSharp_cpg = .TRUE., intSharp_cpg_west = .FALSE., intSharp_spf = .FALSE., intSharp_ufv = .TRUE., intSharp_utw = .FALSE., intSharp_d02 = .TRUE., intSharp_msk = .TRUE., intSharp_flt = .FALSE., strainHard = .TRUE., cnsrv_g = .FALSE., cnsrv_gt = .FALSE., cnsrv_gp = .FALSE., cnsrv_pe = .FALSE.

        logical :: LADInt, LADN2F 
        logical :: use_gradXi
        logical :: use_normFV, use_normInt
	logical :: surface_mask !flag to use masking on surface tension calculation
	logical :: use_FV, use_D04, use_Stagg, use_XiLS      !flag to use FV scheme when solving for kappa in surface tension	
	logical :: use_gradphi !flag to use phi formulation of VF in calculating surface tension
	logical :: use_gradVF, energy_surfTen  !flag to use VF in calculating surface tension
        logical :: weightedcurvature !flag for curvature correction
        logical     :: use_surfaceTension   !flag to turn on/off surface tension (in momentum and energy equations)
        real(rkind) :: surfaceTension_coeff !constant coefficient for surface tension
        real(rkind) :: intSharp_gam, intSharp_eps, intSharp_cut, intSharp_dif, intSharp_tnh, intSharp_pfloor, XiLS_eps
        real(rkind), allocatable, dimension(:,:,:,:) :: intSharp_f,intSharp_h,VFboundDiff,intSharp_fDiff,intSharp_hDiff,intSharp_fFV, antiDiffTerm, DiffTerm, LADVF,intSharp_fDiffFV
        real(rkind), allocatable, dimension(:,:,:,:) :: surfaceTension_f,xi
        real(rkind), allocatable, dimension(:,:,:)   :: surfaceTension_fxx, surfaceTension_fyy, surfaceTension_fzz, surfaceTension_fxy,surfaceTension_fxz, surfaceTension_fyz, surfaceTension_pe
        real(rkind), allocatable, dimension(:,:,:)   :: surfaceTension_pe_x,surfaceTension_pe_y, surfaceTension_pe_z
        real(rkind), allocatable, dimension(:,:,:) :: surfaceTension_e,surfaceTension_fxy_x,surfaceTension_fxy_y, surfaceTension_fxz_x,surfaceTension_fxz_z, surfaceTension_fyz_y, surfaceTension_fyz_z
        real(rkind), allocatable, dimension(:,:,:) :: intSharp_hFV,intSharp_kFV,intSharp_hDiffFV,intSharp_kDiffFV, intSharp_pFV
        
        integer, dimension(2) :: x_bc, y_bc, z_bc
        real(rkind), allocatable, dimension(:,:,:)   :: kappa, maskKappa,VF_intx, VF_inty,VF_intz, DerX, DerY,DerZ, ddx_exact, ddy_exact, Pmix 
        real(rkind), allocatable, dimension(:,:,:)   ::  intX_error, intY_error,derX_error, derY_error, intx_exact, inty_exact, lapTest, DivTest, lap_error, div_error
        real(rkind), allocatable, dimension(:,:,:,:) :: norm, normFV,gradp,gradVF,gradxi
        real(rkind), allocatable, dimension(:,:,:,:,:) :: gradVF_FV
	real(rkind), allocatable, dimension(:,:,:)   :: phi
        real(rkind), allocatable, dimension(:,:,:)   :: fmask, antidiff, intDiff
        real(rkind), allocatable, dimension(:,:,:)   :: buffer_send_1, buffer_send_k_1
        real(rkind), allocatable, dimension(:,:,:)   :: buffer_send_2, buffer_send_k_2
        real(rkind), allocatable, dimension(:,:,:)   :: buffer_recieve_1, buffer_recieve_k_1
        real(rkind), allocatable, dimension(:,:,:)   :: buffer_recieve_2, buffer_recieve_k_2
!        real(rkind), dimension(:,:,:), pointer       :: buffs1, buffs2, buffr1, buffr2
!        real(rkind), dimension(:,:,:), pointer       :: buffsk1, buffsk2, buffrk1, buffrk2

        integer, dimension(8) :: MPI_req, MPI_Stats
      


    contains

        procedure :: init
        procedure :: set_material
        procedure :: pressureLiquidGas
        procedure :: relaxPressure
        procedure :: relaxPressure_os
        procedure :: equilibratePressure
        procedure :: updateP_VF
        procedure :: equilibratePressureTemperature
        procedure :: equilibratePressureTemperature_new
        procedure :: equilibrateTemperature
        procedure :: getLAD
        procedure :: calculate_source
        procedure :: update_g
        procedure :: implicit_plastic
        procedure :: Test_Der_NP
        procedure :: Test_Der_N2F
        procedure :: Test_Der
        procedure :: Test_Der_Periodic
        procedure :: update_Ys
        procedure :: update_eh
        procedure :: update_VF
        procedure :: filter
        procedure :: filter_g
        procedure :: get_rho
        procedure :: get_primitive
        procedure :: get_primitive_g
        procedure :: get_conserved
        procedure :: get_conserved_g
        procedure :: get_ehydro_from_p
        procedure :: get_p_from_ehydro
        procedure :: get_rhoYs_from_gVF
        procedure :: get_emix
        procedure :: get_pmix
        procedure :: get_Tmix
        procedure :: getSOS
        procedure :: get_J
        procedure :: get_q
        procedure :: get_intSharp_clean
        procedure :: get_intSharp_clean2
        procedure :: get_intSharp
        procedure :: get_gradp
        procedure :: get_surfaceTensionPE
        procedure :: get_surfaceTensionCnsrv
        procedure :: getYs
        procedure :: CheckP
        procedure :: getLAD_5eqn
        procedure :: get_surfaceTension
        procedure :: interpolateFV !TODO: check BCs on all FV terms    
        procedure :: divergenceFV
        procedure :: interpolateFV_6 !TODO: check BCs on all FV terms    
        procedure :: divergenceFV_6
        procedure :: gradientFV
        procedure :: get_qmix
        procedure :: get_dt
        procedure :: get_eelastic_devstress
        procedure :: get_mixture_properties
        procedure :: checkNaN
        procedure :: fnumden
        procedure :: fnumden_new
        procedure :: rootfind_nr_1d
        procedure :: rootfind_nr_1d_new
        procedure :: thick_calculations
        final     :: destroy

    end type

    !interface solid_mixture
    !    module procedure init
    !end interface

contains

    !function init(decomp,der,fil,LAD,ns) result(this)
    subroutine init(this,decomp,der,derD02,derStagg,derStaggd02,derD06,derCD06,derD04,interpMid,interpMid02,use_Stagg,LADN2F,LADInt,fil,gfil,LAD,ns,PTeqb,pEqb,pRelax,SOSmodel,use_gTg,updateEtot,useAkshayForm,twoPhaseLAD,LAD5eqn,useOneG,intSharp,usePhiForm,intSharp_cpl,intSharp_cpg,intSharp_cpg_west,intSharp_spf,intSharp_ufv,intSharp_utw,intSharp_d02,intSharp_msk,intSharp_flt,intSharp_gam,intSharp_eps,intSharp_cut,intSharp_dif,intSharp_tnh,intSharp_pfloor,use_surfaceTension,use_normFV,use_normInt,use_gradXi, energy_surfTen,use_gradphi, use_gradVF, surfaceTension_coeff, use_FV,use_XiLS, XiLS_eps,use_D04,surface_mask, weightedcurvature, strainHard,cnsrv_g,cnsrv_gt,cnsrv_gp,cnsrv_pe,x_bc,y_bc,z_bc)

        class(solid_mixture), target,    intent(inout) :: this
        type(decomp_info), target,       intent(in)    :: decomp
        type(filters),     target,       intent(in)    :: fil, gfil
        type(derivatives), target,       intent(in)    :: der,derD02,derD04, derD06, derCD06
        type(derivativesStagg), target,  intent(in)    :: derStagg, derStaggd02
        type(interpolators), target,     intent(in)    :: interpMid, interpMid02
        type(ladobject),   target,       intent(in)    :: LAD
        integer,                         intent(in)    :: ns
        logical,                         intent(in)    :: PTeqb,pEqb,pRelax,updateEtot,useAkshayForm, twoPhaseLAD
        logical,                         intent(in)    :: SOSmodel, LAD5eqn
        logical,                         intent(in)    :: use_gTg,useOneG,intSharp,intSharp_cpl,intSharp_cpg,intSharp_cpg_west,intSharp_spf,intSharp_ufv,usePhiForm, intSharp_utw,intSharp_d02,intSharp_msk,intSharp_flt,strainHard,cnsrv_g,cnsrv_gt,cnsrv_gp,cnsrv_pe
        logical,                         intent(in)    :: use_surfaceTension, use_gradphi,use_Stagg, use_gradVF, use_FV, use_D04,surface_mask, weightedcurvature, use_gradXi, energy_surfTen, use_XiLS,LADN2F,LADInt,use_normFV, use_normInt
        real(rkind),                     intent(in)    :: surfaceTension_coeff,XiLS_eps
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        real(rkind) :: intSharp_gam, intSharp_eps, intSharp_cut, intSharp_dif, intSharp_tnh, intSharp_pfloor


        type(solid), allocatable :: dummy
        integer :: i

        if (ns < 1) call GracefulExit("Must have at least 1 species in the problem. Check input file for errors",3457)

        this%PTeqb      = PTeqb
        this%pEqb       = pEqb
        this%pRelax     = pRelax
        this%SOSmodel   = SOSmodel
        this%use_gTg    = use_gTg
        this%updateEtot     = updateEtot
        this%useOneG        = useOneG
        this%intSharp       = intSharp
        this%intSharp_cpl   = intSharp_cpl
        this%usePhiForm     = usePhiForm
        this%useAkshayForm  = useAkshayForm
        this%twoPhaseLAD    = twoPhaseLAD
        this%LAD5eqn        = LAD5eqn
        this%intSharp_cpg      = intSharp_cpg
        this%intSharp_cpg_west = intSharp_cpg_west
        this%intSharp_spf   = intSharp_spf
        this%intSharp_ufv   = intSharp_ufv
        this%intSharp_utw   = intSharp_utw
        this%intSharp_d02   = intSharp_d02
        this%intSharp_msk   = intSharp_msk
        this%intSharp_flt   = intSharp_flt
        this%intSharp_gam   = intSharp_gam
        this%intSharp_eps   = intSharp_eps
        this%intSharp_cut   = intSharp_cut
        this%intSharp_dif   = intSharp_dif
        this%intSharp_tnh   = intSharp_tnh
        this%intSharp_pfloor   = intSharp_pfloor
        this%XiLS_eps          = XiLS_eps
        this%strainHard = strainHard
        this%cnsrv_g  = cnsrv_g
        this%cnsrv_gt = cnsrv_gt
        this%cnsrv_gp = cnsrv_gp
        this%cnsrv_pe = cnsrv_pe

        this%weightedcurvature    = weightedcurvature
	this%surface_mask 	  = surface_mask
	this%use_FV 		  = use_FV
        this%use_XiLS             = use_XiLS
        this%use_D04              = use_D04
        this%use_surfaceTension   = use_surfaceTension
	this%use_gradphi 	  = use_gradphi
        this%energy_surfTen       = energy_surfTen
        this%use_gradXi           = use_gradXi
        this%use_normFV           = use_normFV
        this%use_normInt          = use_normInt
	this%use_gradVF 	  = use_gradVF  
        this%surfaceTension_coeff = surfaceTension_coeff
        this%use_Stagg            = use_Stagg
        this%LADInt               = LADInt
        this%LADN2F               = LADN2F

        this%x_bc = x_bc
        this%y_bc = y_bc
        this%z_bc = z_bc


        this%ns = ns

        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        this%decomp => decomp
        this%der  => der
        this%derD06 => derD06
        this%derD04  => derD04
        this%derD02  => derD02
        this%derCD06 => derCD06
        this%fil  => fil
        this%gfil => gfil
        this%LAD  => LAD
        this%derStagg  => derStagg
        this%interpMid  => interpMid
        this%derStaggd02  => derStaggd02
        this%interpMid02  => interpMid02

        ! Allocate array of solid objects (Use a dummy to avoid memory leaks)
        allocate(dummy)
        call dummy%init(decomp,der,derD02,derD04,derD06,derStagg,derStaggd02,interpMid,interpMid02,this%use_Stagg,this%LADN2F,this%LADInt,fil,gfil,this%PTeqb,this%pEqb,this%pRelax,this%use_gTg,this%useOneG,this%intSharp,this%intSharp_spf,intSharp_ufv,this%intSharp_d02,intSharp_cut,this%intSharp_cpg_west,this%useAkshayForm,this%twoPhaseLAD,this%LAD5eqn,this%updateEtot,this%strainHard,this%cnsrv_g,this%cnsrv_gt,this%cnsrv_gp,this%cnsrv_pe,this%ns, this%x_bc, this%y_bc, this%z_bc)

        if (allocated(this%material)) deallocate(this%material)
        allocate(this%material(this%ns))!, source=dummy)
   

       print *, "called init mod"



     do i=1,this%ns
            call this%material(i)%init(decomp,der,derD02,derD04,derD06,derStagg,derStaggd02,interpMid,interpMid02,this%use_Stagg,this%LADN2F,this%LADInt,fil,gfil,this%PTeqb,this%pEqb,this%pRelax,this%use_gTg,this%useOneG,this%intSharp,this%intSharp_spf,intSharp_ufv,this%intSharp_d02,intSharp_cut,this%intSharp_cpg_west,this%useAkshayForm,this%twoPhaseLAD,this%LAD5eqn,this%updateEtot,this%strainHard,this%cnsrv_g,this%cnsrv_gt,this%cnsrv_gp,this%cnsrv_pe,this%ns, this%x_bc, this%y_bc, this%z_bc)
        end do
        deallocate(dummy)


         print *, "finished call to init mod"

        this%material(1)%Ys = one
        this%material(1)%VF = one
        do i=2,this%ns
            this%material(i)%Ys = zero
            this%material(i)%VF = zero
        end do

        if(allocated(this%rho0mix)) deallocate(this%rho0mix)
        allocate(this%rho0mix(this%nxp, this%nyp, this%nzp))

        if(allocated(this%mumix)) deallocate(this%mumix)
        allocate(this%mumix(this%nxp, this%nyp, this%nzp))

        if(allocated(this%yieldmix)) deallocate(this%yieldmix)
        allocate(this%yieldmix(this%nxp, this%nyp, this%nzp))

        if(allocated(this%solidVF)) deallocate(this%solidVF)
        allocate(this%solidVF(this%nxp, this%nyp, this%nzp))

        if(allocated(this%intSharp_f)) deallocate(this%intSharp_f)
        allocate(this%intSharp_f(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%intSharp_h)) deallocate(this%intSharp_h)
        allocate(this%intSharp_h(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%intSharp_fDiff)) deallocate(this%intSharp_fDiff)
        allocate(this%intSharp_fDiff(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%intSharp_hDiff)) deallocate(this%intSharp_hDiff)
        allocate(this%intSharp_hDiff(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%intSharp_fFV)) deallocate(this%intSharp_fFV)
        allocate(this%intSharp_fFV(this%nxp, this%nyp, this%nzp,3))

        if(allocated(this%intSharp_fDiffFV)) deallocate(this%intSharp_fDiffFV)
        allocate(this%intSharp_fDiffFV(this%nxp, this%nyp, this%nzp,3))

        if(allocated(this%intSharp_hFV)) deallocate(this%intSharp_hFV)
        allocate(this%intSharp_hFV(this%nxp, this%nyp, this%nzp))

        if(allocated(this%intSharp_pFV)) deallocate(this%intSharp_pFV)
        allocate(this%intSharp_pFV(this%nxp, this%nyp, this%nzp))

        if(allocated(this%intSharp_kFV)) deallocate(this%intSharp_kFV)
        allocate(this%intSharp_kFV(this%nxp, this%nyp, this%nzp))
 
        if(allocated(this%intSharp_hDiffFV)) deallocate(this%intSharp_hDiffFV)
        allocate(this%intSharp_hDiffFV(this%nxp, this%nyp, this%nzp))

        if(allocated(this%intSharp_kDiffFV)) deallocate(this%intSharp_kDiffFV)
        allocate(this%intSharp_kDiffFV(this%nxp, this%nyp, this%nzp))

        if(allocated(this%VFboundDiff)) deallocate(this%VFboundDiff)
        allocate(this%VFboundDiff(this%nxp, this%nyp, this%nzp, this%ns))

        if(allocated(this%antiDiffTerm)) deallocate(this%antiDiffTerm)
        allocate(this%antiDiffTerm(this%nxp, this%nyp, this%nzp, this%ns))

        if(allocated(this%DiffTerm)) deallocate(this%DiffTerm)
        allocate(this%DiffTerm(this%nxp, this%nyp, this%nzp, this%ns))

        if(allocated(this%LADVF)) deallocate(this%LADVF)
        allocate(this%LADVF(this%nxp, this%nyp, this%nzp, this%ns))

        if(allocated(this%surfaceTension_f)) deallocate(this%surfaceTension_f)
        allocate(this%surfaceTension_f(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%surfaceTension_e)) deallocate(this%surfaceTension_e)
        allocate(this%surfaceTension_e(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_pe)) deallocate(this%surfaceTension_pe)
        allocate(this%surfaceTension_pe(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_pe_x)) deallocate(this%surfaceTension_pe_x)
        allocate(this%surfaceTension_pe_x(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_pe_y)) deallocate(this%surfaceTension_pe_y)
        allocate(this%surfaceTension_pe_y(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_pe_z)) deallocate(this%surfaceTension_pe_z)
        allocate(this%surfaceTension_pe_z(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fxx)) deallocate(this%surfaceTension_fxx)
        allocate(this%surfaceTension_fxx(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fyy)) deallocate(this%surfaceTension_fyy)
        allocate(this%surfaceTension_fyy(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fzz)) deallocate(this%surfaceTension_fzz)
        allocate(this%surfaceTension_fzz(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fxy)) deallocate(this%surfaceTension_fxy)
        allocate(this%surfaceTension_fxy(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fxz)) deallocate(this%surfaceTension_fxz)
        allocate(this%surfaceTension_fxz(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fyz)) deallocate(this%surfaceTension_fyz)
        allocate(this%surfaceTension_fyz(this%nxp, this%nyp, this%nzp))
         
        if(allocated(this%surfaceTension_fxy_x)) deallocate(this%surfaceTension_fxy_x)
        allocate(this%surfaceTension_fxy_x(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fxz_x)) deallocate(this%surfaceTension_fxz_x)
        allocate(this%surfaceTension_fxz_x(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fyz_y)) deallocate(this%surfaceTension_fyz_y)
        allocate(this%surfaceTension_fyz_y(this%nxp, this%nyp, this%nzp))
 
       if(allocated(this%surfaceTension_fxy_y)) deallocate(this%surfaceTension_fxy_y)
        allocate(this%surfaceTension_fxy_y(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fxz_z)) deallocate(this%surfaceTension_fxz_z)
        allocate(this%surfaceTension_fxz_z(this%nxp, this%nyp, this%nzp))

        if(allocated(this%surfaceTension_fyz_z)) deallocate(this%surfaceTension_fyz_z)
        allocate(this%surfaceTension_fyz_z(this%nxp, this%nyp, this%nzp))

        if(allocated(this%norm)) deallocate(this%norm)
        allocate(this%norm(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%normFV)) deallocate(this%normFV)
        allocate(this%normFV(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%Pmix)) deallocate(this%Pmix)
        allocate(this%Pmix(this%nxp, this%nyp, this%nzp))

        if(allocated(this%gradp)) deallocate(this%gradp)
        allocate(this%gradp(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%gradVF)) deallocate(this%gradVF)
        allocate(this%gradVF(this%nxp, this%nyp, this%nzp, 3))
        
        if(allocated(this%gradVF_FV)) deallocate(this%gradVF_FV)
        allocate(this%gradVF_FV(this%nxp, this%nyp, this%nzp, 3,3))
        
        if(allocated(this%gradxi)) deallocate(this%gradxi)
        allocate(this%gradxi(this%nxp, this%nyp, this%nzp, 3))

        if(allocated(this%kappa)) deallocate(this%kappa)
        allocate(this%kappa(this%nxp, this%nyp, this%nzp))

        if(allocated(this%intX_error)) deallocate(this%intX_error)
        allocate(this%intX_error(this%nxp, this%nyp, this%nzp))

        if(allocated(this%intY_error)) deallocate(this%intY_error)
        allocate(this%intY_error(this%nxp, this%nyp, this%nzp))

        if(allocated(this%intY_error)) deallocate(this%intY_error)
        allocate(this%intY_error(this%nxp, this%nyp, this%nzp))

        if(allocated(this%derX_error)) deallocate(this%derX_error)
        allocate(this%derX_error(this%nxp, this%nyp, this%nzp))

        if(allocated(this%derY_error)) deallocate(this%derY_error)
        allocate(this%derY_error(this%nxp, this%nyp, this%nzp))

        if(allocated(this%VF_intx)) deallocate(this%VF_intx)
        allocate(this%VF_intx(this%nxp, this%nyp, this%nzp))

        if(allocated(this%VF_inty)) deallocate(this%VF_inty)
        allocate(this%VF_inty(this%nxp, this%nyp, this%nzp))

        if(allocated(this%VF_intz)) deallocate(this%VF_intz)
        allocate(this%VF_intz(this%nxp, this%nyp, this%nzp))

        if(allocated(this%DerZ)) deallocate(this%DerZ)
        allocate(this%DerZ(this%nxp, this%nyp, this%nzp))

        if(allocated(this%DivTest)) deallocate(this%DivTest)
        allocate(this%DivTest(this%nxp, this%nyp, this%nzp))

        if(allocated(this%lapTest)) deallocate(this%lapTest)
        allocate(this%lapTest(this%nxp, this%nyp, this%nzp))

        if(allocated(this%div_error)) deallocate(this%div_error)
        allocate(this%div_error(this%nxp, this%nyp, this%nzp))

        if(allocated(this%lap_error)) deallocate(this%lap_error)
        allocate(this%lap_error(this%nxp, this%nyp, this%nzp))

        if(allocated(this%ddx_exact)) deallocate(this%ddx_exact)
        allocate(this%ddx_exact(this%nxp, this%nyp, this%nzp))

        if(allocated(this%ddy_exact)) deallocate(this%ddy_exact)
        allocate(this%ddy_exact(this%nxp, this%nyp, this%nzp))

         if(allocated(this%intx_exact)) deallocate(this%intx_exact)
        allocate(this%intx_exact(this%nxp, this%nyp, this%nzp))

        if(allocated(this%inty_exact)) deallocate(this%inty_exact)
        allocate(this%inty_exact(this%nxp, this%nyp, this%nzp))

        if(allocated(this%DerY)) deallocate(this%DerY)
        allocate(this%DerY(this%nxp, this%nyp, this%nzp))

        if(allocated(this%DerX)) deallocate(this%DerX)
        allocate(this%DerX(this%nxp, this%nyp, this%nzp))

        if(allocated(this%buffer_send_1)) deallocate(this%buffer_send_1)
        allocate(this%buffer_send_1(1,this%nyp, this%nzp))

        
        if(allocated(this%buffer_send_2)) deallocate(this%buffer_send_2)
        allocate(this%buffer_send_2(1,this%nyp,this%nzp))

        
        if(allocated(this%buffer_recieve_1)) deallocate(this%buffer_recieve_1)
        allocate(this%buffer_recieve_1(1,this%nyp,this%nzp))

        
        if(allocated(this%buffer_recieve_2)) deallocate(this%buffer_recieve_2)
        allocate(this%buffer_recieve_2(1,this%nyp,this%nzp))


        if(allocated(this%buffer_send_k_1)) deallocate(this%buffer_send_k_1)
        allocate(this%buffer_send_k_1(1,this%nyp, this%nzp))


        if(allocated(this%buffer_send_k_2)) deallocate(this%buffer_send_k_2)
        allocate(this%buffer_send_k_2(1,this%nyp,this%nzp))


        if(allocated(this%buffer_recieve_k_1)) deallocate(this%buffer_recieve_k_1)
        allocate(this%buffer_recieve_k_1(1,this%nyp,this%nzp))


        if(allocated(this%buffer_recieve_k_2)) deallocate(this%buffer_recieve_k_2)
        allocate(this%buffer_recieve_k_2(1,this%nyp,this%nzp))

        if(allocated(this%maskKappa)) deallocate(this%maskKappa)
        allocate(this%maskKappa(this%nxp, this%nyp, this%nzp))

	if(allocated(this%phi)) deallocate(this%phi)
        allocate(this%phi(this%nxp, this%nyp, this%nzp))


        if(allocated(this%xi)) deallocate(this%xi)
        allocate(this%xi(this%nxp, this%nyp, this%nzp, 3))

	if(allocated(this%fmask)) deallocate(this%fmask)
        allocate(this%fmask(this%nxp, this%nyp, this%nzp))
       
        if(allocated(this%antidiff)) deallocate(this%antidiff)
        allocate(this%antidiff(this%nxp, this%nyp, this%nzp))

        if(allocated(this%intdiff)) deallocate(this%intdiff)
        allocate(this%intdiff(this%nxp, this%nyp, this%nzp))
 
        call MPI_COMM_SIZE(MPI_COMM_WORLD, this%mpi_size, this%ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, this%mpi_rank, this%ierror)

        this%mpi_rank_prev = mod((this%mpi_rank - 1 + this%mpi_size), this%mpi_size)
        this%mpi_rank_next = mod((this%mpi_rank + 1), this%mpi_size)



    end subroutine
    !end function

    pure elemental subroutine destroy(this)
        type(solid_mixture), intent(inout)  :: this

        if(allocated(this%solidVF)) deallocate(this%solidVF)
        if(allocated(this%yieldmix)) deallocate(this%yieldmix)
        if(allocated(this%mumix)) deallocate(this%mumix)
        if(allocated(this%rho0mix)) deallocate(this%rho0mix)
        if(allocated(this%intSharp_f)) deallocate(this%intSharp_f)
        if(allocated(this%intSharp_h)) deallocate(this%intSharp_h)
        if(allocated(this%intSharp_fDiff)) deallocate(this%intSharp_fDiff)
        if(allocated(this%intSharp_hDiff)) deallocate(this%intSharp_hDiff)
        if(allocated(this%intSharp_fFV)) deallocate(this%intSharp_fFV)
        if(allocated(this%intSharp_hFV)) deallocate(this%intSharp_hFV)
        if(allocated(this%intSharp_pFV)) deallocate(this%intSharp_pFV)
        if(allocated(this%intSharp_kFV)) deallocate(this%intSharp_kFV)
        if(allocated(this%intSharp_fDiffFV)) deallocate(this%intSharp_fDiffFV)
        if(allocated(this%intSharp_hDiffFV)) deallocate(this%intSharp_hDiffFV)
        if(allocated(this%intSharp_kDiffFV)) deallocate(this%intSharp_kDiffFV)
        if(allocated(this%VFboundDiff)) deallocate(this%VFboundDiff)
        if(allocated(this%antiDiffTerm)) deallocate(this%antiDiffTerm)
        if(allocated(this%DiffTerm)) deallocate(this%DiffTerm)
        if(allocated(this%LADVF)) deallocate(this%LADVF)

        if(allocated(this%surfaceTension_f)) deallocate(this%surfaceTension_f)
        if(allocated(this%surfaceTension_e)) deallocate(this%surfaceTension_e)
        if(allocated(this%surfaceTension_pe)) deallocate(this%surfaceTension_pe)
        if(allocated(this%surfaceTension_pe_x)) deallocate(this%surfaceTension_pe_x)
        if(allocated(this%surfaceTension_pe_y)) deallocate(this%surfaceTension_pe_y)
        if(allocated(this%surfaceTension_pe_z)) deallocate(this%surfaceTension_pe_z)
        if(allocated(this%surfaceTension_fxx)) deallocate(this%surfaceTension_fxx)
        if(allocated(this%surfaceTension_fyy)) deallocate(this%surfaceTension_fyy)
        if(allocated(this%surfaceTension_fzz)) deallocate(this%surfaceTension_fzz)
        if(allocated(this%surfaceTension_fxy_x)) deallocate(this%surfaceTension_fxy_x)
        if(allocated(this%surfaceTension_fxz_x)) deallocate(this%surfaceTension_fxz_x)
        if(allocated(this%surfaceTension_fyz_y)) deallocate(this%surfaceTension_fyz_y)
        if(allocated(this%surfaceTension_fxy_y)) deallocate(this%surfaceTension_fxy_y)
        if(allocated(this%surfaceTension_fxz_z)) deallocate(this%surfaceTension_fxz_z)
        if(allocated(this%surfaceTension_fyz_z)) deallocate(this%surfaceTension_fyz_z)
        if(allocated(this%surfaceTension_fxy)) deallocate(this%surfaceTension_fxy)
        if(allocated(this%surfaceTension_fxz)) deallocate(this%surfaceTension_fxz) 
        if(allocated(this%surfaceTension_fyz)) deallocate(this%surfaceTension_fyz)

        if(allocated(this%norm)) deallocate(this%norm)
        if(allocated(this%normFV)) deallocate(this%normFV)
        if(allocated(this%gradp)) deallocate(this%gradp)
        if(allocated(this%gradVF)) deallocate(this%gradVF)
        if(allocated(this%gradVF_FV)) deallocate(this%gradVF_FV)
        if(allocated(this%gradxi)) deallocate(this%gradxi)
        if(allocated(this%kappa)) deallocate(this%kappa)
        if(allocated(this%Pmix)) deallocate(this%Pmix)
        if(allocated(this%VF_intx)) deallocate(this%VF_intx)
        if(allocated(this%VF_inty)) deallocate(this%VF_inty)
        if(allocated(this%VF_intz)) deallocate(this%VF_intz)
        if(allocated(this%DerZ)) deallocate(this%DerZ)
        if(allocated(this%lapTest)) deallocate(this%lapTest)
        if(allocated(this%DivTest)) deallocate(this%DivTest)
        if(allocated(this%div_error)) deallocate(this%div_error)
        if(allocated(this%lap_error)) deallocate(this%lap_error)
        if(allocated(this%intx_exact)) deallocate(this%intx_exact)
        if(allocated(this%inty_exact)) deallocate(this%inty_exact)
        if(allocated(this%ddx_exact)) deallocate(this%ddx_exact)
        if(allocated(this%ddy_exact)) deallocate(this%ddy_exact)
        if(allocated(this%DerY)) deallocate(this%DerY)
        if(allocated(this%DerX)) deallocate(this%DerX)
        if(allocated(this%derY_error)) deallocate(this%derY_error)
        if(allocated(this%derX_error)) deallocate(this%derX_error)
        if(allocated(this%intX_error)) deallocate(this%intY_error)
	if(allocated(this%buffer_send_1)) deallocate(this%buffer_send_1)
        if(allocated(this%buffer_send_2)) deallocate(this%buffer_send_2)
        if(allocated(this%buffer_recieve_1)) deallocate(this%buffer_recieve_1)
        if(allocated(this%buffer_recieve_2)) deallocate(this%buffer_recieve_2)
        if(allocated(this%buffer_send_k_1)) deallocate(this%buffer_send_k_1)
        if(allocated(this%buffer_send_k_2)) deallocate(this%buffer_send_k_2)
        if(allocated(this%buffer_recieve_k_1)) deallocate(this%buffer_recieve_k_1)
        if(allocated(this%buffer_recieve_k_2)) deallocate(this%buffer_recieve_k_2)


        if(allocated(this%maskKappa)) deallocate(this%maskKappa)
	if(allocated(this%fmask)) deallocate(this%fmask)
        if(allocated(this%antidiff)) deallocate(this%antidiff)
        if(allocated(this%intdiff)) deallocate(this%intdiff)

        ! Deallocate array of solids (Destructor of solid should take care of everything else)
        if (allocated(this%material)) deallocate(this%material)

        nullify(this%LAD)
        nullify(this%fil)
        nullify(this%gfil)
        nullify(this%der)
        nullify(this%derD02)
        nullify(this%derD04)
        nullify(this%derD06)
        nullify(this%derCD06)
        nullify(this%decomp)
        nullify(this%derStagg)
        nullify(this%interpMid)

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


     subroutine pressureLiquidGas(this,mixRho,mixE,mixP)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix,rhom,e_species, delta_max, delta_min, del_e, delta_e1, delta_e2, rhom2, rhom1,denom
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rho1gambyone,rho2gambyone, diffPInf,  tmp, q, Cv
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: VF1, VF2, A1,A2,mixT
        integer :: imat, i,j,k,a(3), s = 2
        real(rkind) :: gamfac,eps2, eps1 = 1D-8, minVF = 1D-8, thresh, cutoff

      !  print *, "liquidgas"
        Cv = this%material(1)%Ys*this%material(1)%hydro%Cv+this%material(2)%Ys*this%material(2)%hydro%Cv
        
      !  print *, "Cv"
         
        A1 = this%material(1)%Ys*(this%material(1)%hydro%gam-1)*this%material(1)%hydro%Cv/Cv*(mixRho*mixE- this%material(1)%hydro%PInf)
        A2 = this%material(2)%Ys*(this%material(2)%hydro%gam-1)*this%material(2)%hydro%Cv/Cv*(mixRho*mixE-this%material(2)%hydro%PInf)

        mixP = 0.5*( A1 + A2 -(this%material(1)%hydro%PInf+this%material(2)%hydro%PInf)) + sqrt( 0.25*(A2 - A1- (this%material(2)%hydro%PInf - this%material(1)%hydro%PInf))**2 + A1*A2)        

      !  print *, "mixP"
        do k=1,this%nzp
           do j=1,this%nyp
              do i=1,this%nxp
                 

                 if(mixP(i,j,k) < 0) then
       !          print *, "mixP = ", mixP(i,j,k)
                 endif
                 mixT(i,j,k) = zero
                 do imat = 1, this%ns
                    mixT(i,j,k) = mixT(i,j,k) +this%material(imat)%Ys(i,j,k)*this%material(imat)%hydro%Cv * (mixP(i,j,k) +this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf)/(mixP(i,j,k) +this%material(imat)%hydro%PInf) !original
                 enddo
                 mixT(i,j,k) = ehmix(i,j,k)/mixT(i,j,k)

                 do imat = 1, this%ns
                    this%material(imat)%T(i,j,k) = mixT(i,j,k)
                    this%material(imat)%p(i,j,k) = mixP(i,j,k)

                    this%material(imat)%VF(i,j,k)= mixRho(i,j,k)*this%material(imat)%Ys(i,j,k)*(this%material(imat)%hydro%gam-one)* this%material(imat)%hydro%Cv*mixT(i,j,k)/(mixP(i,j,k) +this%material(imat)%hydro%PInf) !original

                    this%material(imat)%eh(i,j,k) = this%material(imat)%hydro%Cv*mixT(i,j,k)*(mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf) / (mixP(i,j,k) + this%material(imat)%hydro%PInf)
                 end do

              enddo
           enddo
        enddo

        !this may need to be modified for ns>2
        tmp = zero
        do imat = 1, this%ns-one
           tmp = tmp + this%material(imat)%VF
        enddo
        this%material(this%ns)%VF = one - tmp

        tmp = zero
        do imat = 1, this%ns-one
           tmp = tmp + this%material(imat)%eh
        enddo
        this%material(this%ns)%eh = ehmix - tmp

    end subroutine

    subroutine relaxPressure_os(this,rho,u,v,w,mixE,dtsim,mixP)
        class(solid_mixture), intent(inout), target :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w,mixE
        real(rkind),                                        intent(in)  :: dtsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: Fsrc,dt
        type(solid), pointer :: mat1, mat2
        integer :: i, tt 
        real(rkind) :: muc = one

        mat1 => this%material(1); mat2 => this%material(2)


        do tt = 1, 1000
            !print *, tt, maxval(abs(mat1%p-mat2%p))
            Fsrc = abs(mat1%VF*mat2%VF)*(mat1%p - mat2%p)/muc
            !if(maxval(abs(Fsrc))*dtsim < 1.0d-10) exit
            if(maxval(abs((mat1%p-mat2%p)/(mat1%p+mat2%p))) < 1.0d-10) exit
            dt = 0.00000002_rkind*(mat1%VF/Fsrc)

            !do i = 1, 9
            !    mat1%g(:,:,:,i) = mat1%g(:,:,:,i) * (one - dt*Fsrc/mat1%VF)
            !    mat2%g(:,:,:,i) = mat2%g(:,:,:,i) * (one - dt*Fsrc/mat2%VF)
            !enddo

            mat1%VF = mat1%VF + Fsrc*dt
            mat2%VF = mat2%VF - Fsrc*dt
            !print *, tt, maxval((mat1%p*Fsrc*dt)), minval(mat1%p*Fsrc*dt)
            !print *, "p1 before = ", mat1%p(100,1,1)

            mat1%consrv(:,:,:,2) = mat1%consrv(:,:,:,2) - dt*mat1%p*Fsrc
            mat2%consrv(:,:,:,2) = mat2%consrv(:,:,:,2) + dt*mat1%p*Fsrc

            call mat1%get_primitive(rho,u,v,w)
            call mat2%get_primitive(rho,u,v,w)
            call this%get_p_from_ehydro(rho)   ! Get species pressures from species hydrodynamic energy 
            !print *, "p1 after = ", mat1%p(100,1,1)

            if (tt==1000) then
                print *, "Pressure relaxation OS did not converge", maxval(abs(Fsrc)), dtsim
            endif
        enddo

        call this%get_pmix(mixP)              ! Get mixture pressure from species pressures
stop
    end subroutine

    subroutine relaxPressure(this,rho,mixE,mixP)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(1)             :: iparams
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf

        integer :: i,j,k,imat,minlc(1),thisind,othrind
        real(rkind) :: refp, peqb, peqb2, minvf

        real(rkind), dimension(1:this%ns) :: fac

        ehmix = mixE
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
        enddo

        ! equilibrate and reset species pressures, reset volume fractions

        vf   => fparams(  1:this%ns)
        gam  => fparams(  this%ns+1:2*this%ns)
        psph => fparams(2*this%ns+1:3*this%ns)
        pinf => fparams(3*this%ns+1:4*this%ns)
        
        do imat=1,this%ns
          fparams(  this%ns+imat) = this%material(imat)%hydro%gam    ! gamma
          fparams(3*this%ns+imat) = this%material(imat)%hydro%PInf   ! PInf
        enddo

        ! Set reference pressure
        refp = maxval(fparams(3*this%ns+1:4*this%ns)) ! max over all PInfs
        if(refp < 1.0D-5) refp = 1.0D0

        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp

            do imat=1,this%ns
                ! set fparams
                fparams(          imat) = this%material(imat)%VF(i,j,k)    ! volume fractions
                fparams(2*this%ns+imat) = this%material(imat)%p(i,j,k)     ! pressure before eqb
                !print '(1(i4,1x),3(e19.12,1x))', i, fparams(2*this%ns+imat), fparams(imat), this%material(1)%VF(i,j,k) + this%material(2)%VF(i,j,k)
            end do

            !minvf = minval(vf(1:this%ns));     minlc = minloc(vf(1:this%ns))
            !if(minvf < zero) then
            !    ! do not solve non-linear problem. assume relaxed pressure is
            !    ! equal to the pressure of the dominant species
            !    thisind = minlc(1);  othrind = mod(thisind, 2) + 1
            !    psph = this%material(othrind)%p(i,j,k)
            !    do imat=1,this%ns
            !        vf(imat) = this%material(imat)%Ys(i,j,k)* this%material(imat)%eh(i,j,k)*rho(i,j,k)*(gam(imat)-one)/ &
            !                   (psph(imat)+gam(imat)*pinf(imat))
            !    end do
            !    peqb = psph(1)
            !    peqb2 = peqb
            !    
            !else
                ! solve non-linear problem for relaxed pressure

                ! set iparams
                iparams(1) = 1     !   used in fnumden; 2 for PTeqb 

                ! scale all pressures by pref
                fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)/refp

                ! set initial guess
                peqb = sum(fparams(1:this%ns)*fparams(2*this%ns+1:3*this%ns))
                ! solve non-linear equation
                call this%rootfind_nr_1d(peqb,fparams,iparams)

                ! rescale all pressures by refp
                fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*refp
                peqb = peqb*refp

                ! update species VF, eh, ...
                fac = (psph + gam*pinf + (gam-one)*peqb)/(gam*(pinf+peqb))
                vf = vf*fac
                !this%material(1:this%ns)%g = this%material(1:this%ns)%g*fac**third        !! --- not clear if this is needed or if it works
                peqb2 = (rho(i,j,k)*ehmix(i,j,k) - sum(vf*gam*pinf/(gam-one))) / sum(vf/(gam-one))
                psph = peqb2
            !endif
                !print *, 'peqb = ', peqb
                !if(peqb2 < zero) then
                !    write(*,*) i, peqb2, peqb
                !    write(*,*) 'rho: ', rho(i,j,k)
                !    write(*,*) 'eh: ', ehmix(i,j,k)
                !    write(*,*) 'vf : ', vf
                !    write(*,*) 'gam: ', gam
                !    write(*,*) 'fac: ', fac
                !    write(*,*) 'Ys : ', this%material(1)%Ys(i,j,k), this%material(2)%Ys(i,j,k)
                !    write(*,*) 'eh : ', this%material(1)%eh(i,j,k), this%material(2)%eh(i,j,k)
                !    write(*,*) '-----------------------------------'
                !    write(*,*) '-----------------------------------'
                !endif

            mixP(i,j,k) = peqb2
            do imat=1,this%ns
              this%material(imat)%VF(i,j,k) = vf(imat)
              this%material(imat)%p(i,j,k) = psph(imat)
            enddo

            

          enddo
         enddo
        enddo
        ! get species energy from species pressure
        do imat = 1, this%ns
            call this%material(imat)%get_ehydroT_from_p(rho)
        end do

    end subroutine


    subroutine updateP_VF(this,mixRho,mixE,mixP)
        use operators, only: divergence,gradient,filter3D
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: mixP

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix,rhom,e_species, delta_max, delta_min, del_e, delta_e1, delta_e2, rhom2, rhom1,denom
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rho1gambyone,rho2gambyone, diffPInf,  tmp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: VF1,VF2,mixrhoE,fac, num,den,pold,pnew, pvalue
        integer :: imat, i,j,k,a(3), s = 2, iter, iter2 = 0
        real(rkind) :: gamfac,eps2=1d-8,  eps1 = 1D-8, minVF = 1D-8, e =1d-100,thresh, cutoff, p_res

        ! print *, '---'
        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.
        ehmix = mixE*mixRho
        mixrhoE = ehmix
        pold = mixP 
        p_res = 1
        iter = 0
        !DO WHILE((p_res .GT. 1d-14) .AND. (iter .LT. 100))
          fac = (pold + this%material(1)%hydro%gam*this%material(1)%hydro%PInf)*this%material(1)%hydro%onebygam_m1 - &
             (pold + this%material(2)%hydro%gam*this%material(2)%hydro%PInf)*this%material(2)%hydro%onebygam_m1
          this%material(1)%VF = (mixrhoE -  (pold + this%material(2)%hydro%gam*this%material(2)%hydro%PInf)*this%material(2)%hydro%onebygam_m1)/fac
        !this%material(2)%VF = (-mixrhoE +  (mixP + this%material(1)%hydro%gam*this%material(1)%hydro%PInf)*this%material(1)%hydro%onebygam_m1)/fac
        !this%material(2)%VF = 1 - this%material(1)%VF
        !call this%filter(1, this%x_bc, this%y_bc, this%z_bc)
          this%material(2)%VF = 1 - this%material(1)%VF
        do imat = 1, this%ns

            this%material(imat)%p = mixP
            call this%material(imat)%get_ehydroT_from_p(mixRho)

        enddo
    end subroutine

    subroutine equilibratePressure(this,mixRho,mixE,mixP)
        use operators, only: divergence,gradient,filter3D
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix,rhom, e_species, delta_max, delta_min, del_e, delta_e1, delta_e2, rhom2, rhom1,denom
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rho1gambyone, rho2gambyone, diffPInf,  tmp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: VF1, VF2,mixrhoE,fac
        integer :: imat, i,j,k,a(3), s = 2
        real(rkind) :: gamfac,eps2=1d-8,  eps1 = 1D-8, minVF = 1D-8, e = 1d-100,thresh, cutoff

        ! print *, '---'
        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.
        ehmix = mixE !*mixRho
        mixrhoE = ehmix
     !   ehmix = ehmix*mixRho
        tmp = zero

        thresh = 1d-16

!        where((this%material(1)%VF .LE. thresh) .OR. (this%material(1)%Ys .LE. thresh*this%material(1)%elastic%rho0) )
!             this%material(1)%VF = thresh
!             this%material(1)%Ys = thresh*this%material(1)%elastic%rho0
        !     this%material(2)%VF = 1-this%material(1)%VF 
        !     this%material(2)%Ys = 1 - this%material(1)%Ys
              
!        endwhere


!        where((this%material(2)%VF .LE. thresh) .OR. (this%material(2)%Ys .LE. thresh/this%material(1)%elastic%rho0) )
!             this%material(2)%VF = thresh
!             this%material(2)%Ys = thresh/this%material(1)%elastic%rho0
         !    this%material(1)%VF = 1-this%material(2)%VF
         !    this%material(1)%Ys = 1 - this%material(2)%Ys

!        endwhere


!        where((this%material(1)%VF .GE. 1-thresh) .OR. (this%material(1)%Ys .GE. 1-thresh*this%material(1)%elastic%rho0) )
!             this%material(1)%VF = 1-thresh
!             this%material(1)%Ys = 1-thresh*this%material(1)%elastic%rho0
!        endwhere


!        where((this%material(2)%VF .GE. 1-thresh) .OR. (this%material(2)%Ys .GE. 1-thresh/this%material(1)%elastic%rho0) )
!             this%material(2)%VF = 1-thresh
!             this%material(2)%Ys = 1-thresh/this%material(1)%elastic%rho0
!        endwhere

!        do imat = 1, this%ns
!         do k=1,this%nzp
!                do j = 1,this%nyp
!                    do i = 1,this%nxp
!                        if ( this%material(imat)%Ys(i,j,k) *this%material(imat)%VF(i,j,k) < zero ) then
!                            print *, "Found negative Ys*VF in material after thresh ", imat," at ", i, j, k
!                            print *, "  Ys = ", this%material(imat)%Ys(i,j,k)
!                            print *, "  VF = ", this%material(imat)%Vf(i,j,k)
!                            ! exit
!                        end if
!                    end do
!                end do
!            end do
!        enddo

       ! cutoff = 1d-5
       ! where( this%material(1)%Ys .LE. cutoff*this%material(1)%elastic%rho0)
       !    this%material(1)%rhom = cutoff
       !    this%material(1)%VF = this%material(1)%Ys*mixRho/this%material(1)%rhom 
       ! endwhere

       ! where( this%material(2)%Ys .LE. cutoff/this%material(1)%elastic%rho0)
       !    this%material(2)%rhom = cutoff
       !    this%material(2)%VF = this%material(2)%Ys*mixRho/this%material(2)%rhom  
       ! endwhere

       ! this%material(2)%VF = 1 - this%material(1)%VF
       ! this%material(2)%Ys = 1 - this%material(1)%Ys

!       do imat = 1, this%ns

!         call this%material(imat)%getSpeciesDensity(mixRho,rhom)

!       enddo     

        do imat = 1, this%ns
        
      
!       ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
        call this%material(imat)%getSpeciesDensity(mixRho,rhom)
        gamfac =  this%material(imat)%hydro%gam * this%material(imat)%hydro%onebygam_m1*this%material(imat)%hydro%PInf

    !    where( this%material(imat)%VF .LE. thresh)

    !        this%material(imat)%VF = thresh

    !    endwhere
      ! where( this%material(imat)%VF .lt. this%intSharp_cut)      
           !  ehmix = ehmix - gamfac*this%material(imat)%VF
           !ehmix = ehmix - gamfac*this%material(imat)%VF/mixRho
           !  tmp = tmp + this%material(imat)%hydro%onebygam_m1 * this%material(imat)%Ys/rhom
      !     ehmix = ehmix - gamfac*(this%intSharp_cut)
      !     tmp = tmp + this%material(imat)%hydro%onebygam_m1*(this%intSharp_cut)

      !  elsewhere( this%material(imat)%VF .gt. 1-this%intSharp_cut)

      !    ehmix = ehmix - gamfac*(1-this%intSharp_cut)
      !    tmp = tmp + this%material(imat)%hydro%onebygam_m1*(1-this%intSharp_cut)

      !  elsewhere

         ehmix = ehmix - gamfac*this%material(imat)%Ys/rhom
        ! tmp = tmp + this%material(imat)%hydro%onebygam_m1 * this%material(imat)%VF
         tmp = tmp + this%material(imat)%hydro%onebygam_m1 * this%material(imat)%Ys/rhom
      !  endwhere
            

            
 
!             do k=1,this%nzp
!                do j = 1,this%nyp
!                    do i = 1,this%nxp
!                        if ( this%material(imat)%Ys(i,j,k) * this%material(imat)%VF(i,j,k) < zero ) then
!                            print *, "Found negative Ys*VF in material in pressure loop ", imat, " at ", i, j, k
!                            print *, "  Ys = ", this%material(imat)%Ys(i,j,k)
!                            print *, "  VF = ", this%material(imat)%Vf(i,j,k)
                            ! exit
!                        end if
!                    end do
!                end do
!            end do

       enddo


       mixP = ehmix/tmp
       
       !call filter3D(this%decomp, this%fil, mixP, 1, this%x_bc, this%y_bc,this%z_bc)
       !fac = (mixP + this%material(1)%hydro%gam*this%material(1)%hydro%PInf)*this%material(1)%hydro%onebygam_m1 - &
       !      (mixP + this%material(2)%hydro%gam*this%material(2)%hydro%PInf)*this%material(2)%hydro%onebygam_m1
       !this%material(1)%VF = (mixrhoE -  (mixP + this%material(2)%hydro%gam*this%material(2)%hydro%PInf)*this%material(2)%hydro%onebygam_m1)/fac 
       !this%material(2)%VF = (-mixrhoE +  (mixP + this%material(1)%hydro%gam*this%material(1)%hydro%PInf)*this%material(1)%hydro%onebygam_m1)/fac       
         ! a = minloc(ehmix)
         ! !a(1) = 122; a(2:3) = 1
         ! print *, "Mat1 Ys, VF = ", this%material(1)%Ys(a(1),a(2),a(3)), this%material(1)%VF(a(1),a(2),a(3))
         ! print *, "Mat2 Ys, VF = ", this%material(2)%Ys(a(1),a(2),a(3)), this%material(2)%VF(a(1),a(2),a(3))
         ! print *, "Mat1,2 T1   = ", this%material(1)%hydro%gam * this%material(1)%hydro%onebygam_m1 * this%material(1)%hydro%PInf   * this%material(1)%VF(a(1),a(2),a(3))/mixRho(a(1),a(2),a(3)),  this%material(2)%hydro%gam * this%material(2)%hydro%onebygam_m1 * this%material(2)%hydro%PInf   * this%material(2)%VF(a(1),a(2),a(3))/mixRho(a(1),a(2),a(3))
         ! print *, "Mat1,2 eel = ",  this%material(1)%eel(a(1),a(2),a(3)),  this%material(2)%eel(a(1),a(2),a(3))
         ! print *, "Mat1,2 g11 = ",  this%material(1)%g11(a(1),a(2),a(3)),  this%material(2)%g11(a(1),a(2),a(3))
         ! print *, "Mat1,2 g22 = ",  this%material(1)%g22(a(1),a(2),a(3)),  this%material(2)%g22(a(1),a(2),a(3))
         ! print *, "Mix eh      = ",  mixE(a(1),a(2),a(3)), mixE(a(1),a(2),a(3)) -   (this%material(1)%Ys(a(1),a(2),a(3)) * this%material(1)%eel(a(1),a(2),a(3)) + this%material(2)%Ys(a(1),a(2),a(3)) * this%material(2)%eel(a(1),a(2),a(3)) )
         ! print *, "ehmix      = ",  ehmix(a(1),a(2),a(3))
         ! print *, "tmp        = ",  tmp(a(1),a(2),a(3))
         ! print *, "Mat2 Ys, VF = ", this%material(2)%Ys(a(1),a(2),a(3)), this%material(2)%VF(a(1),a(2),a(3))
         ! print *, "Mat1,2 T1   = ", this%material(1)%hydro%gam * this%material(1)%hydro%onebygam_m1 * this%material(1)%hydro%PInf   * this%material(1)%VF(a(1),a(2),a(3))/mixRho(a(1),a(2),a(3)),  this%material(2)%hydro%gam * this%material(2)%hydro%onebygam_m1 * this%material(2)%hydro%PInf   * this%material(2)%VF(a(1),a(2),a(3))/mixRho(a(1),a(2),a(3))
         ! print *, "Mat1,2 eel = ",  this%material(1)%eel(a(1),a(2),a(3)),  this%material(2)%eel(a(1),a(2),a(3))
         ! print *, "Mat1,2 g11 = ",  this%material(1)%g11(a(1),a(2),a(3)),  this%material(2)%g11(a(1),a(2),a(3))
         ! print *, "Mat1,2 g22 = ",  this%material(1)%g22(a(1),a(2),a(3)),  this%material(2)%g22(a(1),a(2),a(3))
         ! print *, "Mix eh      = ",  mixE(a(1),a(2),a(3)), mixE(a(1),a(2),a(3)) -   (this%material(1)%Ys(a(1),a(2),a(3)) * this%material(1)%eel(a(1),a(2),a(3)) + this%material(2)%Ys(a(1),a(2),a(3)) * this%material(2)%eel(a(1),a(2),a(3)) )
         ! print *, "ehmix      = ",  ehmix(a(1),a(2),a(3))
         ! print *, "tmp        = ",  tmp(a(1),a(2),a(3))
         ! print *, "mixP       = ",  mixP(a(1),a(2),a(3))
         
         !print*, "ehmix: ", minval(ehmix), maxval(ehmix)
         !print*, "tmp  : ", minval(tmp), maxval(tmp)
         !print*, "mixP : ", minval(mixP), maxval(mixP)
      
         !if (minval(mixP) < zero) then
         !    print *, "  loc: ", minloc(mixP)
         !end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !! Zoe's scheme for 2 species!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!         !! Retrieve material densities
!         call this%material(1)%getSpeciesDensity(mixRho,rhom1)
!         call this%material(2)%getSpeciesDensity(mixRho,rhom2)
    
!         rho2gambyone = rhom2*(this%material(2)%hydro%gam -1)
!         rho1gambyone = rhom1*(this%material(1)%hydro%gam -1)

!         denom = 0
         
        ! eps1 = this%intSharp_cut
          
!        where(this%material(1)%Ys .lt. eps)
!           VF1 = abs(this%material(1)%Ys)
!        elsewhere(this%material(1)%Ys .gt. 1+eps)
!           VF1 = 1 - abs(this%material(1)%Ys -1)
!        elsewhere
!           VF1 = this%material(1)%Ys
!        endwhere

        ! eps2 = this%intSharp_cut
!         where(this%material(2)%Ys .lt. eps)
!           VF2 = abs(this%material(2)%Ys)
!         elsewhere(this%material(2)%Ys .gt. 1+eps)
!           VF2 = 1 - abs(this%material(2)%Ys - 1)
!         elsewhere
!           VF2 = this%material(2)%Ys
!         endwhere
        
!         denom = VF1*rho2gambyone + VF2*rho1gambyone

!         diffPInf = this%material(1)%hydro%gam*this%material(1)%hydro%PInf - this%material(2)%hydro%gam*this%material(2)%hydro%PInf

         !! Solve for species energys (these formulas are from Conservation of
         !Energy where emix = summation Yi*eh_i & from isobaric assumption p1 = p2 = ei*rho_i*(gam_i-1)
        
      !  e_species = this%material(1)%eh

       !  where(this%material(2)%Ys .lt. eps)
       !     this%material(1)%eh =(rho2gambyone*ehmix + diffPInf*eps)/denom
       !  elsewhere(this%material(2)%Ys .gt. (1))
       !     this%material(1)%eh =(rho2gambyone*ehmix + diffPInf*(1))/denom
       !  elsewhere
       !     this%material(1)%eh = (rho2gambyone*ehmix +  diffPInf*VF2)/denom
       ! endwhere
     
         !del_e = this%material(1)%eh-e_species
         !delta_max = (s - 1)*e_species
         !delta_min = (1/s - 1)*e_species

         ! Use limitors on species energy
         !delta_e1 =  min( del_e, delta_max )
         !delta_e2 = max( delta_e1, delta_min )

         !this%material(1)%eh = e_species + delta_e2

         !e_species = this%material(2)%eh

         !where(this%material(1)%Ys .lt. eps)

        !   this%material(2)%eh =(rho1gambyone*ehmix - diffPInf*eps)/denom
  
        ! elsewhere(this%material(1)%Ys .gt. (1))

        !   this%material(2)%eh =(rho1gambyone*ehmix - diffPInf*(1))/denom

        ! elsewhere

        ! this%material(2)%eh = (rho1gambyone*ehmix -diffPInf*VF1)/denom

       !  endwhere
        ! del_e = this%material(2)%eh-e_species
        ! delta_max = (s - 1)*e_species
        ! delta_min = (1/s - 1)*e_species

         ! Use limitors on species energy
        ! delta_e1 =  min( del_e, delta_max )
        ! delta_e2 = max( delta_e1, delta_min )

        ! this%material(2)%eh = e_species + delta_e2

              
        ! mixP = this%material(1)%eh*rhom1*(this%material(1)%hydro%gam -1)-this%material(1)%hydro%gam*this%material(1)%hydro%PInf

         do imat = 1, this%ns

             this%material(imat)%p = mixP
             call this%material(imat)%get_ehydroT_from_p(mixRho)

         enddo

         !call this%material(1)%hydro%get_T(this%material(1)%eh,this%material(1)%T, rhom1)
         !call this%material(2)%hydro%get_T(this%material(2)%eh,this%material(2)%T, rhom2)
!rint *, " call T "
    end subroutine

    subroutine equilibratePressureTemperature(this,mixRho,mixE,mixP,mixT,isub)
        use reductions, only : P_MINVAL,P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP, mixT

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(4)             :: iparams
        integer, intent(in)               :: isub
        ! real(rkind), dimension(:), pointer :: vf, psph, pinf

        integer :: i,j,k,imat,icount
        real(rkind) :: maxp, peqb,rhomin,pmin,pmax, minVF = 1d-6 !, pest, pdiffmax
        
        ! rhomin = 1.0/eps
        ! pmin = 1.0/eps
        ! pmax = zero
        ! do imat = 1, this%ns
        !    rhomin = min(P_MINVAL(this%material(imat)%consrv(:,:,:,1)),rhomin)
        !    pmin = min(P_MINVAL(this%material(imat)%p),pmin)
        !    pmax = max(P_MAXVAL(this%material(imat)%p),pmax)
        ! enddo
        

        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.
        ehmix = mixE
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
        enddo

        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp
            ! set fparams
            fparams(1) = mixRho(i,j,k)*ehmix(i,j,k)

            ! set iparams
            iparams(1) = 2
            iparams(2) = i; iparams(3) = j; iparams(4) = k;

            maxp = zero; peqb = zero
            do imat=1,this%ns
              !! determine max over all PInfs
              !maxp = maxval(maxp, this%material(imat)%hydro%PInf)

              ! set initial guess
              peqb = peqb + this%material(imat)%VF(i,j,k)*this%material(imat)%p(i,j,k)
            end do
            !pest = peqb

            ! solve non-linear equation
            call this%rootfind_nr_1d(peqb,fparams,iparams) !old
            !call this%rootfind_nr_1d_new(peqb,fparams,iparams,pmin,pmax,icount,isub) !now coupled with bisection method for improved stability --- fixes problem in fnumden -> fnumden_new when negative mass fraction !new

            !pdiffmax = max(dabs(pest-peqb),pdiffmax)

            !! rescale all pressures by maxp
            !fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*maxp
            mixP(i,j,k) = peqb !*maxp


            !add here
            !if (peqb.le.pfloor) then !pressure iteration did not converge
            ! call setVFpressureEq
            ! endif
          enddo
         enddo
        enddo
!print *, 'pafter: ', this%material(1)%p(89,1,1), this%material(2)%p(89,1,1)

        mixT = zero
        do i = 1, this%ns
          mixT = mixT + this%material(i)%Ys*this%material(i)%hydro%Cv * &
                (mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf)/(mixP + this%material(i)%hydro%PInf)
        enddo
        mixT = ehmix/mixT

        do i = 1, this%ns
            this%material(i)%T = mixT
            this%material(i)%p = mixP
            this%material(i)%VF = mixRho*this%material(i)%Ys*(this%material(i)%hydro%gam-one)* &
                                         this%material(i)%hydro%Cv*mixT/(mixP + this%material(i)%hydro%PInf)

            this%material(i)%eh = this%material(i)%hydro%Cv*mixT*(mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf) / &
                                                                 (mixP + this%material(i)%hydro%PInf)
        end do

        !add here
        !do ijk
        !do ns
        !if (VF .lt. cut) then
        ! call setVFpressureEq
        !endif
        
    end subroutine

    subroutine equilibratePressureTemperature_new(this,mixRho,mixE,mixP,mixT,isub,nsubs)
        use reductions, only : P_MINVAL,P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP, mixT

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix,tmp
        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(4)             :: iparams
        integer, intent(in)               :: isub,nsubs
        ! real(rkind), dimension(:), pointer :: vf, psph, pinf

        integer :: i,j,k,imat,icount = zero,icount2 = zero
        real(rkind) :: maxp, peqb,rhomin,pmin,pmax !, pest, pdiffmax
        
        rhomin = 1.0/eps
        pmin = 1.0/eps
        pmax = zero
        do imat = 1, this%ns
           rhomin = min(P_MINVAL(this%material(imat)%consrv(:,:,:,1)),rhomin)
           pmin = min(P_MINVAL(this%material(imat)%p),pmin)
           pmax = max(P_MAXVAL(this%material(imat)%p),pmax)
        enddo
        

        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.
        ehmix = mixE
        do imat = 1, this%ns
           ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel !original
           !ehmix = ehmix - min(max(this%material(imat)%Ys,zero),one) * this%material(imat)%eel  !new
        enddo
!print *, 'pinpiut: ', mixE(89,1,1), ehmix(89,1,1), mixRho(89,1,1)
!print *, 'VF     : ', this%material(1)%VF(89,1,1), this%material(2)%VF(89,1,1)
!print *, 'pstart : ', this%material(1)%p(89,1,1), this%material(2)%p(89,1,1)

        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp
            ! set fparams
            fparams(1) = mixRho(i,j,k)*ehmix(i,j,k)

            ! set iparams
            iparams(1) = 2
            iparams(2) = i; iparams(3) = j; iparams(4) = k;

            maxp = zero; peqb = zero
            do imat=1,this%ns
              !! determine max over all PInfs
              !maxp = maxval(maxp, this%material(imat)%hydro%PInf)

              ! set initial guess
               peqb = peqb + this%material(imat)%VF(i,j,k)*this%material(imat)%p(i,j,k) !original
               !peqb = peqb + min(max(this%material(imat)%VF(i,j,k),zero),one)*this%material(imat)%p(i,j,k)  !new

               ! !debug
               ! if(this%material(imat)%p(i,j,k).le.eps) then
               !    print*,'negative p 1 ',this%material(imat)%p(i,j,k),imat,i,j,k
               ! endif
               ! !end 
            end do
            !pest = peqb

            ! !debug
            ! if(peqb.le.eps) then
            !    print*,'negative p 2 ',peqb,imat,i,j,k
            ! endif
            ! !end

            ! solve non-linear equation
            !call this%rootfind_nr_1d(peqb,fparams,iparams) !original
            call this%rootfind_nr_1d_new(peqb,fparams,iparams,pmin,pmax,icount,icount2,isub,nsubs) !now coupled with bisection method for improved stability --- fixes problem in fnumden -> fnumden_new when negative mass fraction
            
            !debug
            if(peqb.le.eps) then
               print*,'negative p 3 ',peqb,imat,i,j,k
            endif
            !end

            !pdiffmax = max(dabs(pest-peqb),pdiffmax)

            !! rescale all pressures by maxp
            !fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*maxp
            mixP(i,j,k) = peqb !*maxp

          enddo
         enddo
        enddo
!print *, 'pafter: ', this%material(1)%p(89,1,1), this%material(2)%p(89,1,1)
        
        
        do k=1,this%nzp
           do j=1,this%nyp
              do i=1,this%nxp

                 mixT(i,j,k) = zero
                 do imat = 1, this%ns
                    mixT(i,j,k) = mixT(i,j,k) + this%material(imat)%Ys(i,j,k)*this%material(imat)%hydro%Cv * (mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !original
                    !mixT(i,j,k) = mixT(i,j,k) + min(max(this%material(imat)%Ys(i,j,k),zero),one)*this%material(imat)%hydro%Cv * (mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !new
                    !mixT(i,j,k) = mixT(i,j,k) + min(max(this%material(imat)%Ys(i,j,k),zero+1.D-6),one-1.D-6)*this%material(imat)%hydro%Cv * (mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !new
                 enddo
                 mixT(i,j,k) = ehmix(i,j,k)/mixT(i,j,k)

                 do imat = 1, this%ns
                    this%material(imat)%T(i,j,k) = mixT(i,j,k)
                    this%material(imat)%p(i,j,k) = mixP(i,j,k)

                    this%material(imat)%VF(i,j,k) = mixRho(i,j,k)*this%material(imat)%Ys(i,j,k)*(this%material(imat)%hydro%gam-one) * this%material(imat)%hydro%Cv*mixT(i,j,k)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !original
                    !this%material(imat)%VF(i,j,k) = mixRho(i,j,k)*min(max(this%material(imat)%Ys(i,j,k),zero),one)*(this%material(imat)%hydro%gam-one) * this%material(imat)%hydro%Cv*mixT(i,j,k)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !new
                    !this%material(imat)%VF(i,j,k) = mixRho(i,j,k)*min(max(this%material(imat)%Ys(i,j,k),zero+1.D-6),one-1.D-6)*(this%material(imat)%hydro%gam-one) * this%material(imat)%hydro%Cv*mixT(i,j,k)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !new

                    this%material(imat)%eh(i,j,k) = this%material(imat)%hydro%Cv*mixT(i,j,k)*(mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf) / (mixP(i,j,k) + this%material(imat)%hydro%PInf)
                 end do

              enddo
           enddo
        enddo
        
        !this may need to be modified for ns>2
        tmp = zero
        do imat = 1, this%ns-one
           tmp = tmp + this%material(imat)%VF
        enddo
        this%material(this%ns)%VF = one - tmp

        tmp = zero
        do imat = 1, this%ns-one
           tmp = tmp + this%material(imat)%eh
        enddo
        this%material(this%ns)%eh = ehmix - tmp

        ! ehmixNew = zero
        ! do imat = 1, this%ns
        !    ! !call this%material(imat)%get_ehydroT_from_p(rho)
        !    ! e = (p + this%gam*this%PInf) * this%onebygam_m1 / rho
        !    ! T =  (e - this%PInf/rho)/this%Cv
        !    ehmixNew = ehmixNew + (mixP + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf) / ( (this%material(imat)%hydro%gam-one) * mixRho)

        !             this%material(imat)%eh(i,j,k) = this%material(imat)%hydro%Cv*mixT(i,j,k)*(mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf) / (mixP(i,j,k) + this%material(imat)%hydro%PInf)

        !             this%material(imat)%hydro%Cv*mixT(i,j,k) * (this%material(imat)%hydro%gam-one) * mixRho  =  (mixP(i,j,k) + this%material(imat)%hydro%PInf)



        ! end do

        

    end subroutine


subroutine equilibrateTemperature(this,mixRho,mixE,mixP,mixT,isub, nsubs)
        use reductions, only : P_MINVAL,P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP, mixT

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix, tmp
        real(rkind), dimension(4*this%ns), target :: fparams
        real(rkind), dimension(this%ns) :: e_species, T_species, p_species 
        real(rkind) :: Teqb, Peqb2, T, numer, denom, eps_1 = 1D0, eps_2 = 1D0, vf_cutoff= 0
        integer, dimension(4)             :: iparams
        integer, intent(in)               :: isub, nsubs

        integer :: i,j,k,imat,icount = zero, icount2 = zero
        real(rkind) :: maxp, peqb,rhomin,pmin,pmax, delta_max, delta_min, delta_e1, delta_e2, s=2, VF, e_species_n1

        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.

        
        rhomin = 1.0/eps
        pmin = 1.0/eps
        pmax = zero
        do imat = 1, this%ns
           rhomin = min(P_MINVAL(this%material(imat)%consrv(:,:,:,1)),rhomin)
           pmin = min(P_MINVAL(this%material(imat)%p),pmin)
           pmax = max(P_MAXVAL(this%material(imat)%p),pmax)
        enddo



        ehmix = mixE
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
        enddo


        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp
           this%material(2)%rho(i,j,k) = this%material(2)%elastic%rho0

           VF = mixRho(i,j,k)*this%material(2)%Ys(i,j,k)/this%material(2)%rho(i,j,k)

           if ( ( VF .LE. vf_cutoff) .OR. (VF .GE. (1-vf_cutoff)) ) then
           fparams(1) = mixRho(i,j,k)*ehmix(i,j,k)

            ! set iparams
            iparams(1) = 2
            iparams(2) = i; iparams(3) = j; iparams(4) = k;

            maxp = zero; peqb = zero
            do imat=1,this%ns
              !! determine max over all PInfs
              !maxp = maxval(maxp, this%material(imat)%hydro%PInf)

              ! set initial guess
               peqb = peqb + this%material(imat)%VF(i,j,k)*this%material(imat)%p(i,j,k) !original
               !peqb = peqb + min(max(this%material(imat)%VF(i,j,k),zero),one)*this%material(imat)%p(i,j,k)  !new

               ! !debug
               ! if(this%material(imat)%p(i,j,k).le.eps) then
               !    print*,'negative p 1 ',this%material(imat)%p(i,j,k),imat,i,j,k
               ! endif
               ! !end 
            end do
            !pest = peqb

            ! !debug
            ! if(peqb.le.eps) then
            !    print*,'negative p 2 ',peqb,imat,i,j,k
            ! endif
            ! !end

            ! solve non-linear equation
            !call this%rootfind_nr_1d(peqb,fparams,iparams) !original
            call this%rootfind_nr_1d_new(peqb,fparams,iparams,pmin,pmax,icount,icount2,isub,nsubs) !now coupled with bisection method for improved stability --- fixes problem in fnumden -> fnumden_new when negative mass fraction
            
            !debug
            if(peqb.le.eps) then
               print*,'negative p 3 ',peqb,imat,i,j,k

            endif
            !end

           if (ehmix(i,j,k) .LE. 0) then
           print *, ehmix(i,j,k)
           print *, 'ehmix in PTEQB'
           end if



            !! rescale all pressures by maxp
            !fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*maxp
            mixP(i,j,k) = peqb !*maxp
                 mixT(i,j,k) = zero
                 do imat = 1, this%ns
                    mixT(i,j,k) = mixT(i,j,k) + this%material(imat)%Ys(i,j,k)*this%material(imat)%hydro%Cv * (mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !original
                    !mixT(i,j,k) = mixT(i,j,k) + min(max(this%material(imat)%Ys(i,j,k),zero),one)*this%material(imat)%hydro%Cv * (mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !new
                    !mixT(i,j,k) = mixT(i,j,k) + min(max(this%material(imat)%Ys(i,j,k),zero+1.D-6),one-1.D-6)*this%material(imat)%hydro%Cv * (mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !new
                 enddo
                 mixT(i,j,k) = ehmix(i,j,k)/mixT(i,j,k)

                 do imat = 1, this%ns
                    this%material(imat)%T(i,j,k) = mixT(i,j,k)
                    this%material(imat)%p(i,j,k) = mixP(i,j,k)

                    this%material(imat)%VF(i,j,k) = mixRho(i,j,k)*this%material(imat)%Ys(i,j,k)*(this%material(imat)%hydro%gam-one) * this%material(imat)%hydro%Cv*mixT(i,j,k)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !original
                    !this%material(imat)%VF(i,j,k) = mixRho(i,j,k)*min(max(this%material(imat)%Ys(i,j,k),zero),one)*(this%material(imat)%hydro%gam-one) * this%material(imat)%hydro%Cv*mixT(i,j,k)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !new
                    !this%material(imat)%VF(i,j,k) = mixRho(i,j,k)*min(max(this%material(imat)%Ys(i,j,k),zero+1.D-6),one-1.D-6)*(this%material(imat)%hydro%gam-one) * this%material(imat)%hydro%Cv*mixT(i,j,k)/(mixP(i,j,k) + this%material(imat)%hydro%PInf) !new

                    this%material(imat)%eh(i,j,k) = this%material(imat)%hydro%Cv*mixT(i,j,k)*(mixP(i,j,k) + this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf) / (mixP(i,j,k) + this%material(imat)%hydro%PInf)
                 end do

                tmp = zero
                do imat = 1, this%ns-one
                  tmp = tmp + this%material(imat)%VF
                enddo
                this%material(this%ns)%VF = one - tmp

                tmp = zero
                do imat = 1, this%ns-one
                   tmp = tmp + this%material(imat)%eh
                enddo
                this%material(this%ns)%eh = ehmix - tmp

           else
           !! NEW incompressible phase EOS code
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
           !! set reference density
           this%material(2)%rho(i,j,k) = this%material(2)%elastic%rho0

!           print *,'Teqb routine'

           !! Find VF from mass conservation
           this%material(2)%VF(i,j,k) = mixRho(i,j,k)*this%material(2)%Ys(i,j,k)/this%material(2)%rho(i,j,k)
           this%material(1)%VF(i,j,k) = 1 - this%material(2)%VF(i,j,k)
            
           !! Find density of material 1
           this%material(1)%rho(i,j,k) = mixRho(i,j,k)*this%material(1)%Ys(i,j,k)/this%material(1)%VF(i,j,k)

                      do imat = 1, this%ns


              e_species(imat) = this%material(imat)%eh(i,j,k)



              T_species(imat) = this%material(imat)%T(i,j,k)

               if (((e_species(1) .LT. 0) .OR. (e_species(2) .LT. 0))) then

                print *, e_species(imat)
                print *, 'especie sinitial'

                print *, T_species(imat)
                print *, 'T_species initial'
                end if

           enddo




          if ((e_species(1) .LT. 0) .OR. (e_species(2) .LT. 0)) then
           print *, mixRho(i,j,k)
           print *, 'mixRho'


           print *, this%material(1)%p(i,j,k)
           print *, 'pressure1'

           print *, this%material(2)%p(i,j,k)
           print *, 'pressure2'



          print *, this%material(1)%rho(i,j,k)
           print *, 'rho gas'

           print *, ehmix(i,j,k)
           print *, 'ehmix'

           print *, this%material(1)%Ys(i,j,k)
           print *, 'Ys gas'


           print *, this%material(2)%Ys(i,j,k)
           print *, 'Ys liquid'

           print *, this%material(2)%VF(i,j,k)
           print *, 'VF liquid'

          end if
           !!! Iteratively solve for temperature equilibrium

           !  get initial guess for species energy

            eps_1 = 1
            eps_2 = 1

          !!!!! Iterate until convergence is reached !!!!!!
          do while ( (eps_1 .gt. 10d-5) .OR. (eps_2 .gt. 10d-5) )


               
               !find summations
               numer = 0
               denom = 0

               do imat=1,this%ns

                    !!Temperature floor!!
                    if (T_species(imat) .LE. 0) then
                        T_species(imat) = 10D-12
                    endif

                    if (e_species(imat) .LE. 0) then
                        e_species(imat) = 10D-12
                    end if

                    numer = numer + this%material(imat)%Ys(i,j,k)*(e_species(imat) -  T_species(imat)*(this%material(imat)%hydro%Cv) )

                    denom = denom + this%material(imat)%Ys(i,j,k)*(this%material(imat)%hydro%Cv)

               end do
 

               !Calculate Tn+1
               Teqb  = (ehmix(i,j,k) - numer)/denom
                !print *, Teqb
                !print *, 'teqb'


               ! Update convergence criterion
               eps_1 = abs((this%material(1)%Ys(i,j,k)*(Teqb - T_species(1) ) ) / Teqb)
               eps_2 = abs((this%material(2)%Ys(i,j,k)*(Teqb - T_species(2) ) ) / Teqb)

 !               print *, eps_1
 !               print *, 'eps1'

 !               print *, eps_2
 !               print *, 'eps2'


               !print *, '---Material 1---'
               !write(*,'(3(a,e12.5))') 'rho_1 = ', this%material(1)%rho(i,j,k), ', Tspecies  = ', T_species(1), ', e_species = ', e_species(1)
               !write(*,'(3(a,e12.5))') 'VF  = ', this%material(1)%VF(i,j,k)
               !print *, '---Material 2---'
               !write(*,'(3(a,e12.5))') 'rho_2 = ', this%material(2)%rho(i,j,k), ', Tspecies  = ', T_species(2), ', e_species = ', e_species(2)
               !write(*,'(3(a,e12.5))') 'VF  = ', this%material(1)%VF(i,j,k)
               !write(*,'(3(a,e12.5))') 'Teqb  = ', Teqb,'ehmix  = ', ehmix(i,j,k)
               !write(*,'(3(a,e12.5))') 'numerator  = ', numer, ', denominator  = ', denom
         
!               if ( (eps_1 .GT. 10D-5) .OR. (eps_2 .GT. 10D-5) ) then
               ! Update species energy
               do imat = 1, this%ns

                     e_species_n1 = ( Teqb - T_species(imat) )*(this%material(imat)%hydro%Cv)
     
                     !!! delta max and min needed for limitors
                     delta_max = (s - 1)*e_species(imat)
                     delta_min = (1/s - 1)*e_species(imat)

                     !!! Use limitors on species energy
                     delta_e1 =  min( e_species_n1, delta_max )
                     delta_e2 = max( delta_e1, delta_min )

                     !! species energy update

                     e_species(imat) = e_species(imat) + delta_e2

                     T_species(imat) = (e_species(imat) - (this%material(imat)%hydro%PInf/this%material(imat)%rho(i,j,k))) /this%material(imat)%hydro%Cv

                    ! print *, T_species(imat)
                    ! print *, 'tspecies'
                    ! print *, e_species(imat)
                    ! print *, 'especies'

                             

                 enddo
 !               endif

              end do


               !! Set temperature, energy, and pressure
               Peqb2 = 0
               do imat = 1, this%ns

                    this%material(imat)%T(i,j,k) = Teqb
                    this%material(imat)%eh(i,j,k) = e_species(imat) !this%material(imat)%hydro%Cv * Teqb + this%material(imat)%hydro%PInf/this%material(imat)%rho(i,j,k)

                    p_species(imat) =  (this%material(imat)%hydro%gam-one)*this%material(imat)%rho(i,j,k)*this%material(imat)%eh(i,j,k) &
                                       - this%material(imat)%hydro%gam*this%material(imat)%hydro%PInf

                    Peqb2 = Peqb2 + this%material(imat)%VF(i,j,k)*p_species(imat)
               
               if ((e_species(1) .LT. 0) .OR. (e_species(2) .LT. 0)) then        
               print *, e_species(imat)
               print *, 'final especies'

               print *, p_species(imat)
               print *, 'pspecies'

               print *, T_species(imat)
               print *, 'Tspecies'
               endif

              enddo
              
              mixP(i,j,k) = Peqb2
              mixT(i,j,k) = Teqb
                    
 

              if ((e_species(1) .LT. 0) .OR. (e_species(2) .LT. 0)) then
              print *, Teqb
              print *, 'teqb'

              print *, Peqb2

              print *, 'peqb'
              endif

             do imat = 1, this%ns
                 this%material(imat)%p(i,j,k) = p_species(imat)
             enddo
             endif


          enddo
         enddo
        enddo
 

           ! 1) calculate the species densities and volume fractions
           ! (conservation of mass)
           ! reference density is stored in this%material(imat)%elastic%rho0
           ! this%material(2)%rho = this%material(2)%elastic%rho0

           ! 2) iteratively solve for equilibrium temperature
           ! 2a) set initial guess for temperaures and energies based on
           ! previous timestep
           ! 2b) calculate partial derivatives (analytical based Stiff EOS)
           ! 2c) calculate update of T_eqb until it converges
                 ! (also update species energies)
                 ! (may be necessary to use a limiter)

           ! 3) calculate species pressures analytical based on
           ! densities and temperature
           !
           ! 4) calculate mixture pressure based on the volume fractions
           ! and species pressures

!           print *, 'end routine'
!           stop
    end subroutine


    subroutine get_dt(this, rho, delta, dtkap, dtDiff, dtDiff_g, dtDiff_gt, dtDiff_gp, dtplast)
        use reductions, only: P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho
        real(rkind), intent(in)  :: delta
        real(rkind), intent(out) :: dtkap, dtDiff, dtDiff_g, dtDiff_gt, dtDiff_gp, dtplast

        integer :: imat

        dtkap = one/epssmall; dtDiff = dtkap; dtplast = dtDiff; dtDiff_g = dtDiff; dtDiff_gt = dtDiff_g; dtDiff_gp = dtDiff_gt

        do imat = 1, this%ns
          if (.NOT. this%PTeqb) then
              dtkap =   min(dtkap,   one / ( (P_MAXVAL( this%material(imat)%kap*this%material(imat)%T/(rho* delta**4)))**(third) + eps))
          end if
          dtDiff =  min(dtDiff,  one / ( (P_MAXVAL( this%material(imat)%diff/delta**2) + eps)) )
          dtDiff_g =  min(dtDiff_g,  one / ( (P_MAXVAL( this%material(imat)%diff_g/delta**2) + eps)) )
          dtDiff_gt =  min(dtDiff_gt,  one / ( (P_MAXVAL( this%material(imat)%diff_gt/delta**2) + eps)) )
          dtDiff_gp =  min(dtDiff_gp,  one / ( (P_MAXVAL( this%material(imat)%diff_gp/delta**2) + eps)) )
          dtplast = min(dtplast, this%material(imat)%elastic%tau0)
        enddo

        ! For now disable plastic time step limit by setting a large value
        dtplast = real(1.0D32,rkind)

    end subroutine

    ! Subroutine to get species art. conductivities and diffusivities
    subroutine getLAD(this,rho,p,e,u,v,w,duidxj,sos,use_gTg,strainHard,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc,tfloor)
        use operators, only: divergence,gradient,filter3D
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,e,sos,p,u,v,w  ! Mixture density and speed of sound
        real(rkind),dimension(this%nxp, this%nyp,this%nzp,9), target, intent(in)  :: duidxj
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical, intent(in) :: periodicx, periodicy, periodicz
        real(rkind), intent(in) :: tfloor
        real(rkind), dimension(2) :: minYs
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradrYs
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: umag
        integer :: i

        logical, intent(in) :: use_gTg,strainHard
        minYs(1) = this%material(1)%elastic%rho0/this%material(2)%elastic%rho0*this%intSharp_cut
        minYs(2) = this%material(2)%elastic%rho0/this%material(1)%elastic%rho0*this%intSharp_cut
        umag = u*u + v*v + w*w

        do i = 1,this%ns
            call this%material(i)%getPhysicalProperties()
            if (.NOT. this%PTeqb) then
                ! Artificial conductivity
                !call this%LAD%get_conductivity(rho, this%material(i)%eh, this%material(i)%T, sos, &
                !                                    this%material(i)%kap, x_bc, y_bc, z_bc)
                !call this%LAD%get_conductivity(rho,p, this%material(i)%Ys*this%material(i)%eh, e, sos, this%material(i)%kap, x_bc, y_bc, z_bc,tfloor)
            end if


            if( .NOT. this%twoPhaseLAD) then
                ! Artificial diffusivity (grad(Ys) is stored in Ji at this stage) 
                !print *, 'bef LAD:', this%material(1)%Ji(89,1,1,1), sos(89,1,1), this%material(1)%diff(89,1,1)
                call this%LAD%get_diffusivity(this%material(i)%Ys, this%material(i)%Ji(:,:,:,1), &
                                          this%material(i)%Ji(:,:,:,2), this%material(i)%Ji(:,:,:,3), &
                                          sos, this%material(i)%diff, x_bc, y_bc, z_bc)
            else
                call gradient(this%decomp,this%der,rho*this%material(i)%Ys,gradrYs(:,:,:,1),gradrYs(:,:,:,2),gradrYs(:,:,:,3))
                call this%LAD%get_diffusivity_5eqn(rho,this%material(i)%VF,rho*this%material(i)%Ys,gradrYs(:,:,:,1), gradrYs(:,:,:,2),gradrYs(:,:,:,3),umag,duidxj,minYs(i),this%intSharp_cut,sos,this%material(i)%adiff,this%material(i)%rhodiff,this%material(i)%fd,x_bc,y_bc, z_bc)
               !call this%LAD%get_diffusivity_Steve(rho,this%material(i)%VF,rho*this%material(i)%Ys,minYs(i),this%intSharp_cut,sos,this%material(i)%adiff,this%material(i)%rhodiff,x_bc,y_bc,z_bc)
                !call this%LAD%get_diffusivity_Ys(rho,this%material(i)%VF,rho*this%material(i)%Ys,minYs(i),this%intSharp_cut,sos,this%material(i)%Ysdiff,x_bc,y_bc,z_bc)
                this%material(i)%Ysdiff = 0.0
                !call this%LAD%get_diffusivity_N2F(rho,this%material(i)%VF,rho*this%material(i)%Ys,minYs(i),this%intSharp_cut,sos,this%material(i)%adiff_stagg,this%material(i)%rhodiff_stagg, periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
                !call this%LAD%get_diffusivity_Aslani(rho,this%material(i)%VF,rho*this%material(i)%Ys,minYs(i),this%intSharp_cut,sos,this%material(i)%adiff,this%material(i)%rhodiff,x_bc,y_bc,z_bc)
                ! call this%LAD%get_diffusivity_Bounds(rho,this%material(i)%VF,rho*this%material(i)%Ys,minYs(i),this%intSharp_cut,sos,this%material(i)%outdiff,x_bc,y_bc,z_bc)


            endif


!print *, 'aft LAD:', this%material(1)%Ji(89,1,1,1), sos(89,1,1), this%material(1)%diff(89,1,1)!artificial diffusivity for g and g_t equations

            call this%LAD%get_diff_pe(this%material(i)%pe, sos, this%material(i)%diff_pe, x_bc, y_bc, z_bc)
            call this%LAD%get_diff_g(this%material(i)%g, this%material(i)%g_t, this%material(i)%g_p, sos, this%material(i)%diff_g, this%material(i)%diff_gt, this%material(i)%diff_gp, use_gTg, strainHard,x_bc, y_bc, z_bc)

        end do

    end subroutine

    subroutine get_p_from_ehydro(this, rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho

        integer :: imat

        ! get species pressure from species energy
        do imat = 1, this%ns
            call this%material(imat)%get_p_from_ehydro(rho)
        enddo

    end subroutine

    subroutine get_ehydro_from_p(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho

        integer :: imat

        ! get species energy from species pressure
        do imat = 1, this%ns
            call this%material(imat)%get_ehydroT_from_p(rho)             ! computes species ehydro and T from species p
        enddo

    end subroutine

    subroutine get_Tmix(this,T)
        class(solid_mixture), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: T  ! Mixture temperature

        integer :: i

        T = zero
        do i = 1,this%ns
            T = T + this%material(i)%VF * this%material(i)%T  ! Volume fraction weighted sum
        end do
    end subroutine

    subroutine get_mixture_properties(this)
        class(solid_mixture), intent(inout) :: this

        integer :: imat

        ! compute rho0mix, mumix, yieldmix as VF weighted sums
        this%rho0mix  = this%material(1)%elastic%rho0  * this%material(1)%VF
        !this%mumix    = this%material(1)%elastic%mu/this%material(1)%elastic%rho0 * this%material(1)%VF
        this%mumix    = this%material(1)%elastic%mu * this%material(1)%VF
        this%yieldmix = this%material(1)%elastic%yield * this%material(1)%VF
        do imat = 2, this%ns
            this%rho0mix  = this%rho0mix  + this%material(imat)%elastic%rho0 * this%material(imat)%VF
            !this%mumix    = this%mumix    + this%material(imat)%elastic%mu/this%material(imat)%elastic%rho0 * this%material(imat)%VF
            this%mumix    = this%mumix    + this%material(imat)%elastic%mu * this%material(imat)%VF
            this%yieldmix = this%yieldmix + this%material(imat)%elastic%yield * this%material(imat)%VF
        enddo
        !this%mumix = this%mumix * this%rho0mix

        this%solidVF = zero
        do imat = 1, this%ns
           !if(this%material(imat)%elastic%mu > eps) this%solidVF = this%solidVF + this%material(imat)%VF
           if(this%material(imat)%elastic%mu0 > eps) this%solidVF = this%solidVF + this%material(imat)%VF !mca
        enddo

    end subroutine

    subroutine get_eelastic_devstress(this,devstress)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress

        integer :: imat


        !calculate softening and hardening changes to to yield strength and shear modulus
        do imat = 1, this%ns
           call this%material(imat)%getModifiedModulii()
        enddo
        if(this%useOneG) then
           call get_mixture_properties(this)
        endif
           

        devstress = zero

        if(this%useOneG) then
            call this%material(1)%get_eelastic_devstress(this%rho0mix,this%mumix)
            !print '(a,4(e21.14,1x))', 'mu/rho0_m', this%mumix(89,1,1), this%rho0mix(89,1,1), this%mumix(89,1,1)/this%rho0mix(89,1,1)
            !print '(a,4(e21.14,1x))', 'mu/rho0_1', this%material(1)%elastic%mu, this%material(1)%elastic%rho0, this%material(1)%elastic%mu/this%material(1)%elastic%rho0, this%material(1)%VF(89,1,1)
            !print '(a,4(e21.14,1x))', 'mu/rho0_2', this%material(2)%elastic%mu, this%material(2)%elastic%rho0, this%material(2)%elastic%mu/this%material(2)%elastic%rho0, this%material(2)%VF(89,1,1)
            do imat = 2, this%ns
              this%material(imat)%eel = this%material(1)%eel
              !this%material(imat)%devstress = this%material(1)%devstress
            enddo
            devstress = this%material(1)%devstress
            !print '(a,4(e21.14,1x))', 'mixt sxx', this%material(1)%devstress(89,1,1,1)
            !print '(a,4(e21.14,1x))', 'mixt eel', this%material(1)%eel(89,1,1), this%material(2)%eel(89,1,1)
            !print '(a,4(e21.14,1x))', 'mixt Ys ', this%material(1)%Ys(89,1,1), this%material(2)%Ys(89,1,1)
            !print '(a,4(e21.14,1x))', 'tot  mu ', this%mumix(89,1,1)/this%rho0mix(89,1,1), this%mumix(89,1,1), this%rho0mix(89,1,1)
            !print '(a,4(e21.14,1x))', 'eel 1   ', this%material(1)%Ys(89,1,1)*this%material(1)%eel(89,1,1)
            !print '(a,4(e21.14,1x))', 'eel 2   ', this%material(2)%Ys(89,1,1)*this%material(2)%eel(89,1,1)
            !print '(a,4(e21.14,1x))', 'tot  eel', this%material(1)%Ys(89,1,1)*this%material(1)%eel(89,1,1) + this%material(2)%Ys(89,1,1)*this%material(2)%eel(89,1,1)
            !print '(a,9(e21.14,1x))', 'mix1 g  ', this%material(1)%g(89,1,1,:)
            !print '(a,9(e21.14,1x))', 'mix2 g  ', this%material(2)%g(89,1,1,:)
        else
            do imat = 1, this%ns
              call this%material(imat)%get_eelastic_devstress()
              ! print *, "Material ", imat, " sigma:"
              ! print *, "   ",  this%material(imat)%devstress(200,1,1,1), this%material(imat)%devstress(200,1,1,2), this%material(imat)%devstress(200,1,1,3)
              ! print *, "   ",  this%material(imat)%devstress(200,1,1,2), this%material(imat)%devstress(200,1,1,4), this%material(imat)%devstress(200,1,1,5)
              ! print *, "   ",  this%material(imat)%devstress(200,1,1,3), this%material(imat)%devstress(200,1,1,5), this%material(imat)%devstress(200,1,1,6)
              devstress(:,:,:,1) = devstress(:,:,:,1) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,1)
              devstress(:,:,:,2) = devstress(:,:,:,2) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,2)
              devstress(:,:,:,3) = devstress(:,:,:,3) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,3)
              devstress(:,:,:,4) = devstress(:,:,:,4) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,4)
              devstress(:,:,:,5) = devstress(:,:,:,5) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,5)
              devstress(:,:,:,6) = devstress(:,:,:,6) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,6)
            end do
            !print '(a,4(e21.14,1x))', 'mult sxx', devstress(89,1,1,1)
            !print '(a,4(e21.14,1x))', 'mult eel', this%material(1)%eel(89,1,1), this%material(2)%eel(89,1,1)
            !print '(a,4(e21.14,1x))', 'mult Ys ', this%material(1)%Ys(89,1,1), this%material(2)%Ys(89,1,1)
            !print '(a,4(e21.14,1x))', 'tot  mu1', (this%material(1)%VF(89,1,1)*this%material(1)%elastic%mu + this%material(2)%VF(89,1,1)*this%material(2)%elastic%mu)/(this%material(1)%VF(89,1,1)*this%material(1)%elastic%rho0 + this%material(2)%VF(89,1,1)*this%material(2)%elastic%rho0), (this%material(1)%VF(89,1,1)*this%material(1)%elastic%mu + this%material(2)%VF(89,1,1)*this%material(2)%elastic%mu), (this%material(1)%VF(89,1,1)*this%material(1)%elastic%rho0 + this%material(2)%VF(89,1,1)*this%material(2)%elastic%rho0)
            !print '(a,4(e21.14,1x))', 'tot  mu2', (this%material(1)%Ys(89,1,1)*this%material(1)%elastic%mu/this%material(1)%elastic%rho0 + this%material(2)%Ys(89,1,1)*this%material(2)%elastic%mu/this%material(2)%elastic%rho0)
            !print '(a,4(e21.14,1x))', 'tot  eel', this%material(1)%Ys(89,1,1)*this%material(1)%eel(89,1,1) + this%material(2)%Ys(89,1,1)*this%material(2)%eel(89,1,1)
            !print '(a,9(e21.14,1x))', 'mat1 g  ', this%material(1)%g(89,1,1,:)
            !print '(a,9(e21.14,1x))', 'mat2 g  ', this%material(2)%g(89,1,1,:)
            !print '(a,4(e21.14,1x))', 'mat1 sxx', this%material(1)%devstress(89,1,1,1), this%material(1)%VF(89,1,1)
            !print '(a,4(e21.14,1x))', 'mat2 sxx', this%material(2)%devstress(89,1,1,1), this%material(2)%VF(89,1,1)
        endif

    end subroutine

    subroutine get_J(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho  ! Mixture density

        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: sumJx,sumJy,sumJz
!print *, 'diff Ji:', this%material(1)%Ji(89,1,1,1), this%material(1)%diff(89,1,1)
        sumJx = zero; sumJy = zero; sumJz = zero;
        ! Get diff*gradYs (gradYs are in Ji's)
        do i=1,this%ns
            this%material(i)%Ji(:,:,:,1) = this%material(i)%Ji(:,:,:,1)*this%material(i)%diff
            sumJx = sumJx + this%material(i)%Ji(:,:,:,1)

            this%material(i)%Ji(:,:,:,2) = this%material(i)%Ji(:,:,:,2)*this%material(i)%diff
            sumJy = sumJy + this%material(i)%Ji(:,:,:,2)

            this%material(i)%Ji(:,:,:,3) = this%material(i)%Ji(:,:,:,3)*this%material(i)%diff
            sumJz = sumJz + this%material(i)%Ji(:,:,:,3)
        end do

        ! Correct Ji's so this sum becomes zero (No net diffusive flux)
        do i=1,this%ns
            this%material(i)%Ji(:,:,:,1) = -rho*( this%material(i)%Ji(:,:,:,1) - this%material(i)%Ys*sumJx )
            this%material(i)%Ji(:,:,:,2) = -rho*( this%material(i)%Ji(:,:,:,2) - this%material(i)%Ys*sumJy )
            this%material(i)%Ji(:,:,:,3) = -rho*( this%material(i)%Ji(:,:,:,3) - this%material(i)%Ys*sumJz )
        end do

    end subroutine

    pure subroutine get_conserved(this,rho,u,v,w)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%get_conserved(rho,u,v,w)
        end do

    end subroutine

    subroutine get_primitive(this,rho,u,v,w,e,devstress,p,sos)
        use operators, only: divergence,gradient,filter3D
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: u,v,w,e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: p, sos

        real(rkind) :: minYs1,minVF = 1d-6, minYs2
        integer :: imat
       
        minVF = 1d-6

        if(this%PTeqb) then
        do imat = 1, this%ns
          call this%material(imat)%get_primitive(rho,u,v,w)
        end do
        endif

        if(this%pEqb) then


         do imat = 1, this%ns
            call this%material(imat)%get_primitive(rho,u,v,w)
         enddo

        ! minYs2 = minVF/this%material(1)%elastic%rho0

        ! where( (this%material(2)%VF .LE. minVF) .OR. (this%material(2)%Ys .LE. minYs2))

        !   this%material(2)%rhom = this%material(2)%elastic%rho0
          ! this%material(2)%VF = this%material(2)%consrv(:,:,:,1) / this%material(2)%elastic%rho0

        ! endwhere

        ! this%material(1)%VF = 1 - this%material(2)%VF
      
        ! minYs1 = minVF*this%material(1)%elastic%rho0

        ! where( (this%material(1)%VF .LE. minVF) .OR. (this%material(1)%Ys .LE. minYs1))

        !   this%material(1)%rhom = this%material(1)%elastic%rho0
          ! this%material(1)%VF = this%material(1)%consrv(:,:,:,1) /this%material(1)%elastic%rho0
           

        ! endwhere

         

       !  this%material(2)%Ys = this%material(2)%consrv(:,:,:,1) / rho
       !  minYs2 = minVF/this%material(1)%elastic%rho0

       !  where((this%material(2)%Ys .LE. minYs2) .OR. (this%material(2)%VF .LE. minVF))

       !    this%material(2)%Ys = minYs2
       !    this%material(2)%VF = minVF
        
       !  endwhere

       !  this%material(2)%rhom  = (this%material(2)%consrv(:,:,:,1) + this%material(2)%elastic%rho0*epssmall)/(this%material(2)%VF + epssmall)

         !this%material(1)%VF = 1 - this%material(2)%VF
       !  this%material(1)%Ys = this%material(1)%consrv(:,:,:,1) / rho
       !  minYs1 = minVF*this%material(1)%elastic%rho0

       !  where((this%material(1)%Ys .LE. minYs1) .OR. (this%material(1)%VF .LE. minVF))

       !    this%material(1)%Ys = minYs1
       !    this%material(1)%VF = minVF
  
       !  endwhere

       !  this%material(1)%rhom  = (this%material(1)%consrv(:,:,:,1) +this%material(1)%elastic%rho0*epssmall)/(this%material(1)%VF + epssmall)
       !  rho = this%material(1)%rhom*this%material(1)%VF + this%material(2)%rhom*this%material(2)%VF
    
       !  do imat = 1, this%ns
       !    call gradient(this%decomp,this%der,this%material(imat)%Ys,this%material(imat)%Ji(:,:,:,1),this%material(imat)%Ji(:,:,:,2),this%material(imat)%Ji(:,:,:,3))
       !  enddo

        endif

  
        call this%get_eelastic_devstress(devstress)   ! Get species elastic energies, and mixture and species devstress
        if(this%ns == 1) then
          this%material(1)%eh = e - this%material(1)%eel ! Since eh equation is not exact and this is a better alternative for 1 species
        endif

        if(this%PTeqb) then    ! --- why is this condition needed here? this provides initial guess for pressure even in pRelax case -- NSG
            call this%get_p_from_ehydro(rho)   ! Get species pressures from species hydrodynamic energy 
            call this%get_pmix(p)              ! Get mixture pressure from species pressures
        endif
        
        call this%getSOS(rho,p,sos)

    end subroutine

    subroutine get_primitive_g(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho
        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%get_primitive_g(rho)
        end do
      end subroutine get_primitive_g

    subroutine get_conserved_g(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho
        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%get_conserved_g(rho)
        end do
      end subroutine get_conserved_g

    subroutine get_rho(this,rho)
        use reductions, only: P_MAXVAL, P_MINVAL
        use decomp_2d,  only: nrank
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rho
        real(rkind) :: rhomin

        integer :: imat

        rho = zero
        do imat = 1, this%ns

           !rhomin = P_MINVAL(this%material(imat)%consrv(:,:,:,1))
           !if (nrank.eq.0) print*,imat,rhomin
        

          rho = rho + this%material(imat)%consrv(:,:,:,1)
        end do
!print *, '--rho--', rho(89,1,1), this%material(1)%consrv(89,1,1,1), this%material(2)%consrv(89,1,1,1)
    end subroutine

    subroutine get_rhoYs_from_gVF(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhom

        integer :: imat

        !f(this%pEqb) then

         ! rho = zero
          !do imat = 1, this%ns
          !  call this%material(imat)%getSpeciesDensity(rhom)
          !  rho = rho + this%material(imat)%VF * rhom
         ! end do

         ! do imat = 1,this%ns
         !   call this%material(imat)%getSpeciesDensity(rhom)
         !   this%material(imat)%Ys = this%material(imat)%VF * rhom / rho
         ! end do

       ! else

          rho = zero
          do imat = 1, this%ns
            call this%material(imat)%getSpeciesDensity_from_g(rhom)
            rho = rho + this%material(imat)%VF * rhom
          end do

          do imat = 1,this%ns
            call this%material(imat)%getSpeciesDensity_from_g(rhom)
            this%material(imat)%Ys = this%material(imat)%VF * rhom / rho
          end do

      ! endif
        
    end subroutine

    subroutine get_q(this,x_bc,y_bc,z_bc)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        class(solid_mixture), intent(inout) :: this
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: i
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: tmp1_in_x, tmp2_in_x
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: tmp1_in_y
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: tmp1_in_z, tmp2_in_z

        do i = 1,this%ns
            if(.NOT. this%PTeqb) then
                ! Step 1: Get qy
                !call this%der%ddy(this%material(i)%T,tmp1_in_y,y_bc(1),y_bc(2))
                !this%material(i)%qi(:,:,:,2) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y
                call this%der%ddy(this%material(i)%Ys*this%material(i)%eh,tmp1_in_y,y_bc(1),y_bc(2)) ! --change2
                this%material(i)%qi(:,:,:,2) = -this%material(i)%kap*tmp1_in_y

                ! Step 2: Get qx
                !call transpose_y_to_x(this%material(i)%T,tmp1_in_x,this%decomp)
                !call this%der%ddx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
                !call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
                !this%material(i)%qi(:,:,:,1) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y
                call transpose_y_to_x(this%material(i)%Ys*this%material(i)%eh,tmp1_in_x,this%decomp) ! --change2
                call this%der%ddx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
                call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
                this%material(i)%qi(:,:,:,1) = -this%material(i)%kap*tmp1_in_y

                ! Step 3: Get qz
                !call transpose_y_to_z(this%material(i)%T,tmp1_in_z,this%decomp)
                !call this%der%ddz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
                !call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
                !this%material(i)%qi(:,:,:,3) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y
                call transpose_y_to_z(this%material(i)%Ys*this%material(i)%eh,tmp1_in_z,this%decomp) ! --change2
                call this%der%ddz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
                call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
                this%material(i)%qi(:,:,:,3) = -this%material(i)%kap*tmp1_in_y
            else
                this%material(i)%qi(:,:,:,:) = zero     ! calculated in sgrid
            endif

            if( .NOT. this%twoPhaseLAD) then
            ! If multispecies, add the inter-species enthalpy flux
               if (this%ns .GT. 1) then
                   call this%material(i)%get_enthalpy(tmp1_in_y)
                   this%material(i)%qi(:,:,:,1) = this%material(i)%qi(:,:,:,1) + ( tmp1_in_y * this%material(i)%Ji(:,:,:,1) )
                   this%material(i)%qi(:,:,:,2) = this%material(i)%qi(:,:,:,2) + ( tmp1_in_y * this%material(i)%Ji(:,:,:,2) )
                   this%material(i)%qi(:,:,:,3) = this%material(i)%qi(:,:,:,3) + ( tmp1_in_y * this%material(i)%Ji(:,:,:,3) )
               end if
            end if
        end do

        ! Done
    end subroutine


    subroutine getLAD_5eqn(this,rho,p,e,Frho,Fenergy,Fp,x_bc,y_bc,z_bc,dx,dy,dz,periodicx,periodicy,periodicz)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        use operators, only: divergence,gradient,filter3D, interpolateFV,interpolateFV_x, interpolateFV_y, interpolateFV_z,gradFV_N2Fx,gradFV_N2Fy,gradFV_N2Fz
        use constants,       only: zero,epssmall,eps,one,two,third,half
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3), intent(out) :: Frho, Fenergy, Fp
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: hi
        real(rkind), intent(in) :: dx,dy,dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,p,e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradRYs, gradVF,gradRYs_int, gradVF_int, rhodiff_int, p_int,adiff_int, hi_int, rhom_int,rho_int,gradYs_int,outdiff_int, Ysdiff_int, VF_int, esum, gamsum, e_int, Ys_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: dRYdx_x,dRYdy_y,dRdz_z,dVFdx_x,dVFdy_y,dVFdz_z, rhom
         
        integer :: i, imat
        logical :: periodicx,periodicy,periodicz

        Frho = 0.0
        Fenergy = 0.0
        Fp      = 0.0
        if( this%LADInt ) then
          do i = 1,this%ns

            call this%material(i)%get_enthalpy(hi)
            call this%material(i)%getSpeciesDensity(rho, rhom)
            call gradient(this%decomp,this%der,this%material(i)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))
            call gradient(this%decomp,this%der,rho*this%material(i)%Ys,gradRYs(:,:,:,1),gradRYs(:,:,:,2),gradRYs(:,:,:,3))

            call interpolateFV_x(this%decomp,this%interpMid,gradVF(:,:,:,1),gradVF_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV_y(this%decomp,this%interpMid,gradVF(:,:,:,2),gradVF_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV_z(this%decomp,this%interpMid,gradVF(:,:,:,3),gradVF_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

            call interpolateFV_x(this%decomp,this%interpMid,gradRYs(:,:,:,1),gradRYs_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV_y(this%decomp,this%interpMid,gradRYs(:,:,:,2),gradRYs_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV_z(this%decomp,this%interpMid,gradRYs(:,:,:,3),gradRYs_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

          
            call interpolateFV(this%decomp,this%interpMid,this%material(i)%rhodiff,rhodiff_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,p,p_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,rhom*hi,hi_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
 
            do imat = 1,3
             Frho(:,:,:,imat)    = Frho(:,:,:,imat) + rhodiff_int(:,:,:,imat)*gradRYs_int(:,:,:,imat) !   + this%material(i)%rhodiff*gradRYs(:,:,:,imat)
             !Fenergy(:,:,:,imat)    = Fenergy(:,:,:,imat) + hi*this%material(i)%rhodiff*gradRYs(:,:,:,imat)
             !Fenergy(:,:,:,imat) = Fenergy(:,:,:,imat) + rhodiff_int(:,:,:,imat)*gradVF_int(:,:,:,imat)*( (p_int(:,:,:,imat) + this%material(i)%hydro%Pinf*this%material(i)%hydro%gam)*this%material(i)%hydro%onebygam_m1)
             Fenergy(:,:,:,imat) = Fenergy(:,:,:,imat) + rhodiff_int(:,:,:,imat)*gradVF_int(:,:,:,imat)*hi_int(:,:,:,imat)
            enddo
          !this%material(i)%rhodiff*gradVF(:,:,:,imat)*( (p + this%material(i)%hydro%Pinf*this%material(i)%hydro%gam)*this%material(i)%hydro%onebygam_m1
          enddo 

        elseif(this%LADN2F) then

          do i = 1,this%ns

            call this%material(i)%get_enthalpy(hi)
            call this%material(i)%getSpeciesDensity(rho, rhom)
            call gradient(this%decomp,this%der,this%material(i)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))
            call gradient(this%decomp,this%der,rho*this%material(i)%Ys,gradRYs(:,:,:,1),gradRYs(:,:,:,2),gradRYs(:,:,:,3))

            call gradFV_N2Fx(this%decomp,this%derStagg,this%material(i)%VF,gradVF_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call gradFV_N2Fy(this%decomp,this%derStagg,this%material(i)%VF,gradVF_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call gradFV_N2Fz(this%decomp,this%derStagg,this%material(i)%VF,gradVF_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

            call gradFV_N2Fx(this%decomp,this%derStagg,rho*this%material(i)%Ys,gradRYs_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call gradFV_N2Fy(this%decomp,this%derStagg,rho*this%material(i)%Ys,gradRYs_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call gradFV_N2Fz(this%decomp,this%derStagg,rho*this%material(i)%Ys,gradRYs_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

            call gradFV_N2Fx(this%decomp,this%derStagg,this%material(i)%Ys,gradYs_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call gradFV_N2Fy(this%decomp,this%derStagg,this%material(i)%Ys,gradYs_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call gradFV_N2Fz(this%decomp,this%derStagg,this%material(i)%Ys,gradYs_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

            call interpolateFV(this%decomp,this%interpMid,this%material(i)%Ysdiff,Ysdiff_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,this%material(i)%adiff,adiff_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,this%material(i)%rhodiff,rhodiff_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,this%material(i)%outdiff,outdiff_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,p,p_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,rho,rho_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,e,e_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,rhom,rhom_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,hi,hi_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,this%material(i)%VF,VF_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
            call interpolateFV(this%decomp,this%interpMid,this%material(i)%Ys,Ys_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc) 
            do imat = 1,3
             Frho(:,:,:,imat)    = Frho(:,:,:,imat) + ( rhodiff_int(:,:,:,imat)*gradRYs_int(:,:,:,imat)) !rho_int(:,:,:,imat)*Ysdiff_int(:,:,:,imat)*gradYs_int(:,:,:,imat)
 
             Fenergy(:,:,:,imat) = Fenergy(:,:,:,imat) + (adiff_int(:,:,:,imat)*gradVF_int(:,:,:,imat))*((this%material(i)%hydro%gam*p_int(:,:,:, imat) + this%material(i)%hydro%gam*this%material(i)%hydro%Pinf)*this%material(i)%hydro%onebygam_m1 )

            enddo
          !this%material(i)%rhodiff*gradVF(:,:,:,imat)*( (p +
          !this%material(i)%hydro%Pinf*this%material(i)%hydro%gam)*this%material(i)%hydro%onebygam_m1
          enddo


        else

           do i = 1,this%ns

            call this%material(i)%get_enthalpy(hi)
            call gradient(this%decomp,this%der,this%material(i)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))
            call gradient(this%decomp,this%der,rho*this%material(i)%Ys,gradRYs(:,:,:,1),gradRYs(:,:,:,2),gradRYs(:,:,:,3))
        


            do imat = 1,3
             Frho(:,:,:,imat)    = Frho(:,:,:,imat) + this%material(i)%rhodiff*gradRYs(:,:,:,imat)
             Fenergy(:,:,:,imat) = Fenergy(:,:,:,imat) + this%material(i)%rhodiff*gradVF(:,:,:,imat)*( (p +this%material(i)%hydro%Pinf*this%material(i)%hydro%gam)*this%material(i)%hydro%onebygam_m1 )
            enddo
          enddo



        endif 

    end subroutine

    subroutine get_intSharp_clean(this,rho,x_bc,y_bc,z_bc,dx,dy,dz,periodicx,periodicy,periodicz,u,v,w)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        use operators, only: divergence,gradient,filter3D
        use constants,       only: zero,epssmall,eps,one,two,third,half
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), intent(in) :: dx,dy,dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: norm,gradRhoYs,gradVF,gradVFdiff,fv_f,fv_h,tmp4, gradphi
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns) :: rhoi,VFbound,hi,spf_a,spf_r
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: phiint, gradxi,VFint,uFVint,vFVint,wFVint,gFVint,gtFVint,gpFVint,rhoFVint, kernel
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,this%ns) :: antiDiffFVint,rhoiFVint,hiFVint,rhoiVFint_local
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,3) :: NMint,gradVF_FV,gradVFint
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp,VF_fil,antiDiff,mask,RhoYsbound,filt,antiDiffFV,fmask,tanhmask,mask2,maskDiff,spf_f,spf_h,GVFmag,GVFmagT,antiDiffT

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns) :: J_i,VF_RHS_i, Kij_coeff_i
        real(rkind) :: intSharp_alp = 0.1, r= 0.5, nmask = 40, intSharp_adm =1.0D-1,e = 1d-100, intSharp_exp = -1.0D0,gradDiff,md1,md2 !, intSharp_tnh = 0.1!1.0D-2

        integer :: i,j,ii,jj,kk,iflag = one,im,jm,km,k
        logical :: useTiwari = .FALSE., useRhoYsbound = .FALSE., useTotalRho =.FALSE., useRhoLocal = .FALSE.
        logical :: periodicx,periodicy,periodicz

        do i = 1, this%ns

           !! gradient for out of bounds diffusion
           call gradient(this%decomp,this%der,this%material(i)%VF,gradVFdiff(:,:,:,1),gradVFdiff(:,:,:,2),gradVFdiff(:,:,:,3))

           if(this%usePhiForm) then

               !! Construct xi and account for outbounds behaviour
               where( this%material(i)%VF .GE. 1-this%intSharp_cut)
                    this%xi(:,:,:,i) = this%intSharp_eps*log( (1-2*this%intSharp_cut + e )/ (e))*(1/(1-2*this%intSharp_cut))
               elsewhere( this%material(i)%VF .LE. this%intSharp_cut )                
                    this%xi(:,:,:,i) = this%intSharp_eps*log( ( e )/ (1-2*this%intSharp_cut + e))*(1/(1-2*this%intSharp_cut))
               elsewhere
                    this%xi(:,:,:,i) = this%intSharp_eps*(1/(1-2*this%intSharp_cut))*log( ( this%material(i)%VF - this%intSharp_cut + e )/ (1 - this%intSharp_cut - this%material(i)%VF + e) )
               endwhere



               if(this%intSharp_d02) then    
                     call gradient(this%decomp,this%derD02,this%xi(:,:,:,i),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3)) !low order
               else        
                     call gradient(this%decomp,this%derD06,this%xi(:,:,:,i),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3))             
               endif


           else

              if(this%intSharp_d02)then
                    call gradient(this%decomp,this%derD02,this%material(i)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))!low order
              else
                    call gradient(this%decomp,this%derD06,this%material(i)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))
              endif

           endif


           if (this%usePhiForm) then
           !magnitude of surface vector

              call filter3D(this%decomp, this%gfil, gradxi(:,:,:,1), iflag, x_bc,y_bc,z_bc)
              call filter3D(this%decomp, this%gfil, gradxi(:,:,:,2), iflag,x_bc,y_bc,z_bc)
              call filter3D(this%decomp, this%gfil, gradxi(:,:,:,3), iflag,x_bc,y_bc,z_bc)

              GVFmag = sqrt( gradxi(:,:,:,1)**two + gradxi(:,:,:,2)**two +  gradxi(:,:,:,3)**two )

              !surface normal
              where (GVFmag < eps)
                 norm(:,:,:,1) = zero
                 norm(:,:,:,2) = zero
                 norm(:,:,:,3) = zero
             elsewhere
                 norm(:,:,:,1) = gradxi(:,:,:,1) / GVFmag
                 norm(:,:,:,2) = gradxi(:,:,:,2) / GVFmag
                 norm(:,:,:,3) = gradxi(:,:,:,3) / GVFmag
              endwhere

          else
             
              GVFmag = sqrt( gradVF(:,:,:,1)**two + gradVF(:,:,:,2)**two +  gradVF(:,:,:,3)**two )

              !surface normal
              where (GVFmag < eps)
                 norm(:,:,:,1) = zero
                 norm(:,:,:,2) = zero
                 norm(:,:,:,3) = zero
             elsewhere
                 norm(:,:,:,1) = gradVF(:,:,:,1) / GVFmag
                 norm(:,:,:,2) = gradVF(:,:,:,2) / GVFmag
                 norm(:,:,:,3) = gradVF(:,:,:,3) / GVFmag
              endwhere

 
          endif

          call interpolateFV(this,this%material(i)%VF,VFint,periodicx,periodicy,periodicz,this%x_bc, this%y_bc, this%z_bc)
         !TODO: make sure these BCS for surface normal are correct
          call interpolateFV(this,norm(:,:,:,1),NMint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc, this%y_bc, this%z_bc)
          call interpolateFV(this,norm(:,:,:,2),NMint(:,:,:,:,2),periodicx,periodicy,periodicz, this%x_bc,-this%y_bc, this%z_bc)
          call interpolateFV(this,norm(:,:,:,3),NMint(:,:,:,:,3),periodicx,periodicy,periodicz, this%x_bc, this%y_bc,-this%z_bc)

          call interpolateFV(this,u,uFVint,periodicx,periodicy,periodicz,-this%x_bc, this%y_bc,this%z_bc)
          call interpolateFV(this,v,vFVint,periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
          call interpolateFV(this,w,wFVint,periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)
          call interpolateFV(this,rho,rhoFVint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
      
          !! TO DO: COMPARE RHOLOCAL AND RHO CONST 
          if(useRhoLocal) then !use local component density --- recommended
              call interpolateFV(this,rhoi(:,:,:,i),rhoiFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
          else
              call interpolateFV(this,rhoi(:,:,:,i)*zero+this%material(i)%elastic%rho0,rhoiFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
          endif
          

          call interpolateFV(this,gradVF(:,:,:,1),gradVFint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc,this%y_bc, this%z_bc)
          call interpolateFV(this,gradVF(:,:,:,2),gradVFint(:,:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
          call interpolateFV(this,gradVF(:,:,:,3),gradVFint(:,:,:,:,3),periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)

          !calculate antiDiffFVint term
          antiDiffFVint(:,:,:,1,i) = -this%intSharp_gam *(VFint(:,:,:,1)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,1))*NMint(:,:,:,1,1)
          antiDiffFVint(:,:,:,2,i) = -this%intSharp_gam *(VFint(:,:,:,2)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,2))*NMint(:,:,:,2,2)
          antiDiffFVint(:,:,:,3,i) = -this%intSharp_gam *(VFint(:,:,:,3)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,3))*NMint(:,:,:,3,3)



          if (this%usePhiForm) then
              call interpolateFV(this,this%xi,phiint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              antiDiffFVint(:,:,:,1,i) = -this%intSharp_gam*(0.25*(1 -(tanh((1-2*this%intSharp_cut)* phiint(:,:,:,1)/(2 *this%intSharp_eps)))**2) & 
                                             -0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,1)/(2 *this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,1,1)
              antiDiffFVint(:,:,:,2,i) = -this%intSharp_gam*(0.25*(1 -(tanh((1-2*this%intSharp_cut)* phiint(:,:,:,2)/(2 *this%intSharp_eps)))**2) &
                                             -0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,2)/(2 *this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,2,2)
              antiDiffFVint(:,:,:,3,i) = -this%intSharp_gam*(0.25*(1 -(tanh((1-2*this%intSharp_cut)* phiint(:,:,:,3)/(2 *this%intSharp_eps)))**2) &
                                             -0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,3)/(2 *this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,3,3)
              call divergenceFV(this,antiDiffFVint(:,:,:,:,i),this%antidiff,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc) 
          endif

          call gradientFV(this,this%material(i)%VF,gradVF_FV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
          antiDiffFVint(:,:,:,1,i) = antiDiffFVint(:,:,:,1,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,1,1))
          antiDiffFVint(:,:,:,2,i) = antiDiffFVint(:,:,:,2,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,2,2))
          antiDiffFVint(:,:,:,3,i) = antiDiffFVint(:,:,:,3,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,3,3))

          if( this%intSharp_msk ) then

            
              !FV mask term for antidiffusion
              where(VFint(:,:,:,1).lt.this%intSharp_cut)
                 antiDiffFVint(:,:,:,1,i) = zero
              elsewhere(VFint(:,:,:,1).gt.one-this%intSharp_cut)
                 antiDiffFVint(:,:,:,1,i) = zero
             endwhere
              where(VFint(:,:,:,2).lt.this%intSharp_cut)
                 antiDiffFVint(:,:,:,2,i) = zero
              elsewhere(VFint(:,:,:,2).gt.one-this%intSharp_cut)
                 antiDiffFVint(:,:,:,2,i) = zero
              endwhere
              where(VFint(:,:,:,3).lt.this%intSharp_cut)
                 antiDiffFVint(:,:,:,3,i) = zero
              elsewhere(VFint(:,:,:,3).gt.one-this%intSharp_cut)
                 antiDiffFVint(:,:,:,3,i) = zero
              endwhere



          endif

          if(i.eq.this%ns) then
             antiDiffFVint(:,:,:,:,this%ns) = zero
             do j=1,this%ns-1
                antiDiffFVint(:,:,:,:,this%ns) = antiDiffFVint(:,:,:,:,this%ns) - antiDiffFVint(:,:,:,:,j)
             enddo
          endif

          !compute divergence and calculate RHS terms
          call divergenceFV(this,antiDiffFVint(:,:,:,:,i),this%material(i)%intSharp_aFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
          call divergenceFV(this,rhoiFVint(:,:,:,:,i)*antiDiffFVint(:,:,:,:,i),this%material(i)%intSharp_RFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

           this%VFboundDiff(:,:,:,i) = zero
           maskDiff = zero
           md1 = this%intSharp_cut**half !need to set a larger cutoff toeffectively blend FV and FD
           do ii=1,this%nxp
              do jj=1,this%nyp
                 do kk=1,this%nzp
                    this%VFboundDiff(ii,jj,kk,i) = max(this%VFboundDiff(ii,jj,kk,i), (md1-this%material(i)%VF(ii,jj,kk))/(md1)) !VF < sqrt(cut)
                    this%VFboundDiff(ii,jj,kk,i) = max(this%VFboundDiff(ii,jj,kk,i),(this%material(i)%VF(ii,jj,kk)-(one-md1))/(md1)) !VF > 1-sqrt(cut)
                 enddo
              enddo
           enddo

          this%material(i)%intSharp_aDiff = zero

          if(this%intSharp_msk) then
                 this%VFboundDiff(:,:,:,i) = this%intSharp_dif*this%VFboundDiff(:,:,:,i)
                 mask2 = mask
              else
                 this%VFboundDiff(:,:,:,i) = zero
                 mask2 = one !too diffusive for FD !original

             do j=1,3

                this%material(i)%intSharp_aDiff(:,:,:,j) = this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)
                !Gaussian filter to smooth out-of-bounds diffusion --
                !similar to LAD
                call filter3D(this%decomp,this%gfil,this%material(i)%intSharp_aDiff(:,:,:,j), iflag, x_bc, y_bc,z_bc)

             enddo
           endif
       enddo

       

       tmp4 = zero
       do j = 1,this%ns-1
            tmp4(:,:,:,1) = tmp4(:,:,:,1) + this%material(j)%intSharp_aDiff(:,:,:,1)
            tmp4(:,:,:,2) = tmp4(:,:,:,2) + this%material(j)%intSharp_aDiff(:,:,:,2) 
            tmp4(:,:,:,3) = tmp4(:,:,:,3) + this%material(j)%intSharp_aDiff(:,:,:,3)
       enddo
       this%material(this%ns)%intSharp_aDiff(:,:,:,1) = -tmp4(:,:,:,1)
       this%material(this%ns)%intSharp_aDiff(:,:,:,2) = -tmp4(:,:,:,2)
       this%material(this%ns)%intSharp_aDiff(:,:,:,3) = -tmp4(:,:,:,3)

       do i = 1, this%ns
           this%material(i)%intSharp_RDiff(:,:,:,1) = rhoi(:,:,:,i)*this%material(i)%intSharp_aDiff(:,:,:,1)
           this%material(i)%intSharp_RDiff(:,:,:,2) = rhoi(:,:,:,i)*this%material(i)%intSharp_aDiff(:,:,:,2)
           this%material(i)%intSharp_RDiff(:,:,:,3) = rhoi(:,:,:,i)*this%material(i)%intSharp_aDiff(:,:,:,3)
       enddo
 

       this%intSharp_kFV   = zero
       this%intSharp_hFV   = zero
       this%intSharp_pFV   = zero
       this%intSharp_hDiff = zero
       !enthalpy term is nonzero if coupling is turned on
       if(this%intSharp_cpl) then
              
              do i = 1,this%ns
                  
                  call this%material(i)%get_enthalpy(hi(:,:,:,i))
                  !low order FD terms
                  this%intSharp_h(:,:,:,1) = this%intSharp_h(:,:,:,1) + hi(:,:,:,i) * this%material(i)%intSharp_R(:,:,:,1)
                  this%intSharp_h(:,:,:,2) = this%intSharp_h(:,:,:,2) + hi(:,:,:,i) * this%material(i)%intSharp_R(:,:,:,2)
                  this%intSharp_h(:,:,:,3) = this%intSharp_h(:,:,:,3) + hi(:,:,:,i) * this%material(i)%intSharp_R(:,:,:,3)

                  !high order FD terms
                  this%intSharp_hDiff(:,:,:,1) = this%intSharp_hDiff(:,:,:,1) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,1)
                  this%intSharp_hDiff(:,:,:,2) = this%intSharp_hDiff(:,:,:,2) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,2)
                  this%intSharp_hDiff(:,:,:,3) = this%intSharp_hDiff(:,:,:,3) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,3) 

                  call interpolateFV(this,hi(:,:,:,i),hiFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              enddo

              !finish FV enthalpy term
              fv_h = zero
              do i=1,this%ns
                  fv_h = fv_h + rhoiFVint(:,:,:,:,i)*antiDiffFVint(:,:,:,:,i)*(hiFVint(:,:,:,:,i)+half*sqrt(uFVint**two+vFVint**two+wFVint**two))
              enddo

              call divergenceFV(this,fv_h,this%intSharp_hFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

         endif

         this%intSharp_fFV = zero
         this%intSharp_fDiff = zero
         !momentum term is nonzero if coupling is turned on
         if(this%intSharp_cpl) then
                 fv_f = zero
                 !FV term
                 do i = 1,this%ns
                       fv_f = fv_f + rhoiFVint(:,:,:,:,i)*antiDiffFVint(:,:,:,:,i)
                 enddo

                 call divergenceFV(this,fv_f*uFVint,this%intSharp_fFV(:,:,:,1),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                 call divergenceFV(this,fv_f*vFVint,this%intSharp_fFV(:,:,:,2),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                 call divergenceFV(this,fv_f*wFVint,this%intSharp_fFV(:,:,:,3),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

                do i = 1, this%ns

                  this%intSharp_fDiff(:,:,:,1) = this%intSharp_fDiff(:,:,:,1) + this%material(i)%intSharp_RDiff(:,:,:,1)
                  this%intSharp_fDiff(:,:,:,2) = this%intSharp_fDiff(:,:,:,2) + this%material(i)%intSharp_RDiff(:,:,:,2)
                  this%intSharp_fDiff(:,:,:,3) = this%intSharp_fDiff(:,:,:,3) + this%material(i)%intSharp_RDiff(:,:,:,3)

               enddo

        endif

                   !kinematic fluxes
        do i=1,this%ns
              this%material(i)%intSharp_rg = zero
              this%material(i)%intSharp_rgt = zero
              this%material(i)%intSharp_rgp = zero
              this%material(i)%intSharp_rgDiff = zero
              this%material(i)%intSharp_rgtDiff = zero
              this%material(i)%intSharp_rgpDiff = zero
              this%material(i)%intSharp_rgFV  = zero
              this%material(i)%intSharp_rgtFV = zero
              this%material(i)%intSharp_rgpFV = zero
              this%material(i)%intSharp_gFV  = zero
              this%material(i)%intSharp_gtFV = zero
              this%material(i)%intSharp_gpFV = zero
              this%material(i)%intSharp_a    = zero
              this%material(i)%intSharp_R    = zero
              this%intSharp_h    = zero
              this%intSharp_f    = zero
         enddo



 


    end subroutine

    subroutine get_intSharp_clean2(this,rho,x_bc,y_bc,z_bc,dx,dy,dz,periodicx,periodicy,periodicz,u,v,w)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        use operators, only: divergence,gradient,filter3D,interpolateFV_x,interpolateFV_y, interpolateFV_z, gradFV_N2Fx, gradFV_N2Fy,gradFV_N2Fz
        use constants,       only: zero,epssmall,eps,one,two,third,half
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), intent(in) :: dx,dy,dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,u,v,w
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: norm,gradRhoYs,gradVF,gradVFdiff,fv_f,fv_h,tmp4, gradphi, fv_k, Db_int,hDiff,kDiff,uDiff
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns) :: rhoi,VFbound,hi,spf_a,spf_r
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: phiint,gradxi,VFint,uFVint,vFVint,wFVint,gFVint,gtFVint,gpFVint,rhoFVint, kerne,gradFVFV, gradFV_N2F
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: uFVint_6,vFVint_6,wFVint_6,rhoFVint_6
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,this%ns) :: antiDiffFVint,rhoiFVint,hiFVint, rhoiFVint_local,pFVint, intDiff,hiFVint_6, rhoiFVint_6, pFVint_6
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,3) :: NMint,gradVF_FV,gradVFint, gradXi_FV
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp,VF_fil,antiDiff,mask,RhoYsbound,filt,antiDiffFV,fmask,tanhmask,mask2,maskDiff,spf_f,spf_h,GVFmag,GVFmagT,antiDiffT,rhom, Db
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: gradVF_x,gradVF_y, gradVF_z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns) :: J_i,VF_RHS_i, Kij_coeff_i
        real(rkind) :: intSharp_alp = 0.1, r= 0.5, nmask = 40, intSharp_adm =1.0D-1,e = 1d-100, intSharp_exp = -1.0D0,gradDiff,md1,md2 !, intSharp_tnh = 0.1
!1.0D-2
        integer :: i,j,ii,jj,kk,iflag = one,im,jm,km,k
        logical :: useTiwari = .FALSE., useRhoYsbound = .FALSE., useTotalRho = .FALSE.
        logical :: periodicx,periodicy,periodicz, useGradPsi = .FALSE., useRhoLocal = .FALSE., useHighOrder = .TRUE.,useYSbound = .TRUE., useNewSPF = .TRUE., useNewSPFfull = .FALSE.

            do i=1,this%ns
              this%material(i)%intSharp_a = zero
              this%material(i)%intSharp_aDiff = zero
              this%material(i)%intSharp_R = zero
              this%material(i)%intSharp_RDiff = zero
              this%material(i)%intSharp_RFV = zero
              this%material(i)%intSharp_rg = zero
              this%material(i)%intSharp_rgt = zero
              this%material(i)%intSharp_rgp = zero
              this%material(i)%intSharp_rgDiff = zero
              this%material(i)%intSharp_rgtDiff = zero
              this%material(i)%intSharp_rgpDiff = zero
              this%material(i)%intSharp_rgFV = zero
              this%material(i)%intSharp_rgtFV = zero
              this%material(i)%intSharp_rgpFV = zero
              this%material(i)%intSharp_gFV = zero
              this%material(i)%intSharp_gtFV = zero
              this%material(i)%intSharp_gpFV = zero
            enddo
          
            do i = 1,this%ns
               !component density
               rhoi(:,:,:,i) = rho*this%material(i)%Ys/this%material(i)%VF 
            enddo
            
            do i = 1,this%ns
             !! gradVF derivatives
             call gradient(this%decomp,this%derD02,this%material(i)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))
             call gradient(this%decomp,this%derD06,this%material(i)%VF,gradVFdiff(:,:,:,1),gradVFdiff(:,:,:,2),gradVFdiff(:,:,:,3))

             !! Construct Psi
             where( this%material(i)%VF .GE. 1-this%intSharp_cut)
                 this%xi(:,:,:,i) = this%intSharp_eps*log( (1-2*this%intSharp_cut + e )/ (e))*(1/(1-2*this%intSharp_cut))

             elsewhere( this%material(i)%VF .LE. this%intSharp_cut )
                 this%xi(:,:,:,i) = this%intSharp_eps*log( ( e )/ (1-2*this%intSharp_cut + e))*(1/(1-2*this%intSharp_cut))
             elsewhere
                 this%xi(:,:,:,i) = this%intSharp_eps*(1/(1-2*this%intSharp_cut))*log( ( this%material(i)%VF - this%intSharp_cut + e )/ (1 - this%intSharp_cut - this%material(i)%VF + e) )

             endwhere

             call gradientFV(this,this%xi,gradXi_FV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

             !! Psi Grad
             if(this%intSharp_d02) then
               call gradient(this%decomp,this%derD02,this%xi(:,:,:,i),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3))
!low order
             else
               call gradient(this%decomp,this%derCD06,this%xi(:,:,:,i),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3))
             endif

             !! Psi Norm
             
             call filter3D(this%decomp, this%gfil, gradxi(:,:,:,1), iflag, x_bc,y_bc,z_bc)
             call filter3D(this%decomp, this%gfil, gradxi(:,:,:,2), iflag,x_bc,y_bc,z_bc)
             call filter3D(this%decomp, this%gfil, gradxi(:,:,:,3), iflag,x_bc,y_bc,z_bc)

             GVFmag = sqrt( gradxi(:,:,:,1)**two + gradxi(:,:,:,2)**two + gradxi(:,:,:,3)**two )

             fmask = 1 - (1 - 4*this%material(i)%VF*(1-this%material(i)%VF))**nmask

             !surface normal
             where (GVFmag < eps)
                norm(:,:,:,1) = zero
                norm(:,:,:,2) = zero
                norm(:,:,:,3) = zero
             elsewhere
                norm(:,:,:,1) = gradxi(:,:,:,1) / GVFmag
                norm(:,:,:,2) = gradxi(:,:,:,2) / GVFmag
                norm(:,:,:,3) = gradxi(:,:,:,3) / GVFmag
             endwhere

             !interpolate nodes to faces: ( i, j, k ) -> ( i+1/2, j+1/2, k+1/2 )
              call interpolateFV(this,this%material(i)%VF,VFint,periodicx,periodicy,periodicz,this%x_bc, this%y_bc, this%z_bc)
              !TODO: make sure these BCS for surface normal are correct

              call interpolateFV(this,norm(:,:,:,1),NMint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc,this%y_bc, this%z_bc)
              call interpolateFV(this,norm(:,:,:,2),NMint(:,:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
              call interpolateFV(this,norm(:,:,:,3),NMint(:,:,:,:,3),periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)


              call interpolateFV(this,u,uFVint,periodicx,periodicy,periodicz,-this%x_bc, this%y_bc,this%z_bc)
              call interpolateFV(this,v,vFVint,periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
              call interpolateFV(this,w,wFVint,periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)
              call interpolateFV(this,rho,rhoFVint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

              call interpolateFV_6(this,u,uFVint_6,periodicx,periodicy,periodicz,-this%x_bc,this%y_bc,this%z_bc)
              call interpolateFV_6(this,v,vFVint_6,periodicx,periodicy,periodicz,this%x_bc,-this%y_bc,this%z_bc)
              call interpolateFV_6(this,w,wFVint_6,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,-this%z_bc)
              call interpolateFV_6(this,rho,rhoFVint_6,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

              if(useRhoLocal) then !use local component density --- recommended
                 call interpolateFV(this,rhoi(:,:,:,i),rhoiFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              else
                 call interpolateFV(this,rhoi(:,:,:,i)*zero+this%material(i)%elastic%rho0,rhoiFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              endif

              call interpolateFV(this,rhoi(:,:,:,i),rhoiFVint_local(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              if(this%intSharp_d02) then !which FD VF gradient to interpolate
                 !not used
                 !TODO: make sure these BCS for gradVF are correct
                 call interpolateFV(this,gradVF(:,:,:,1),gradVFint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc,this%y_bc, this%z_bc)
                 call interpolateFV(this,gradVF(:,:,:,2),gradVFint(:,:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
                 call interpolateFV(this,gradVF(:,:,:,3),gradVFint(:,:,:,:,3),periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)

              else
                 !used in high order version
                 !TODO: make sure these BCS for gradVF are correct


                 call interpolateFV_6(this,gradVFdiff(:,:,:,1),gradVFint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc,this%y_bc, this%z_bc)
                 call interpolateFV_6(this,gradVFdiff(:,:,:,2),gradVFint(:,:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
                 call interpolateFV_6(this,gradVFdiff(:,:,:,3),gradVFint(:,:,:,:,3),periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)
              endif


              this%normFV(:,:,:,1) = NMint(:,:,:,1,1)
              this%normFV(:,:,:,2) = NMint(:,:,:,2,2)
              this%normFV(:,:,:,3) = NMint(:,:,:,3,3)

              !calculate antiDiffFVint term
              antiDiffFVint(:,:,:,1,i) = -this%intSharp_gam * (VFint(:,:,:,1)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,1))*NMint(:,:,:,1,1)
              antiDiffFVint(:,:,:,2,i) = -this%intSharp_gam * (VFint(:,:,:,2)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,2))*NMint(:,:,:,2,2)
              antiDiffFVint(:,:,:,3,i) = -this%intSharp_gam * (VFint(:,:,:,3)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,3))*NMint(:,:,:,3,3)
              this%antidiff = antiDiffFVint(:,:,:,1,i)


!              call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,1, i),
!              iflag, x_bc,y_bc,z_bc)
!              call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,2, i),
!              iflag, x_bc,y_bc,z_bc)
!              call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,3, i),
!              iflag, x_bc,y_bc,z_bc)


              if (this%usePhiForm) then
                 call interpolateFV(this,this%xi,phiint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                 antiDiffFVint(:,:,:,1,i) = -this%intSharp_gam*(0.25*(1 - (tanh((1-2*this%intSharp_cut)* phiint(:,:,:,1)/(2 *this%intSharp_eps)))**2)-0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,1)/(2 *this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,1,1)
                 antiDiffFVint(:,:,:,2,i) = -this%intSharp_gam*(0.25*(1 - (tanh((1-2*this%intSharp_cut)* phiint(:,:,:,2)/(2 *this%intSharp_eps)))**2)-0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,2)/(2 *this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,2,2)
                 antiDiffFVint(:,:,:,3,i) = -this%intSharp_gam*(0.25*(1 -(tanh((1-2*this%intSharp_cut)* phiint(:,:,:,3)/(2 *this%intSharp_eps)))**2)-0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,3)/(2 *this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,3,3)


                 call divergenceFV(this,antiDiffFVint(:,:,:,:,i),this%antidiff,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

              endif

               !! This is HIGH ORDER
              call gradFV_N2Fx(this%decomp,this%derStagg,this%material(i)%VF,gradFV_N2F(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
              call gradFV_N2Fy(this%decomp,this%derStagg,this%material(i)%VF,gradFV_N2F(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
              call gradFV_N2Fz(this%decomp,this%derStagg,this%material(i)%VF,gradFV_N2F(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

                 !call
                 !divergenceFV(this,gradFVFV,this%intdiff,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

              antiDiffFVint(:,:,:,1,i) = antiDiffFVint(:,:,:,1,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,1))
              antiDiffFVint(:,:,:,2,i) = antiDiffFVint(:,:,:,2,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,2))
              antiDiffFVint(:,:,:,3,i) = antiDiffFVint(:,:,:,3,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,3))

             if(this%intSharp_msk) then
                where(VFint(:,:,:,1).lt.this%intSharp_cut)
                   antiDiffFVint(:,:,:,1,i) = zero
                elsewhere(VFint(:,:,:,1).gt.one-this%intSharp_cut)
                   antiDiffFVint(:,:,:,1,i) = zero
                endwhere
                where(VFint(:,:,:,2).lt.this%intSharp_cut)
                   antiDiffFVint(:,:,:,2,i) = zero
                elsewhere(VFint(:,:,:,2).gt.one-this%intSharp_cut)
                   antiDiffFVint(:,:,:,2,i) = zero
                endwhere
                where(VFint(:,:,:,3).lt.this%intSharp_cut)
                   antiDiffFVint(:,:,:,3,i) = zero
                elsewhere(VFint(:,:,:,3).gt.one-this%intSharp_cut)
                   antiDiffFVint(:,:,:,3,i) = zero
                endwhere
             endif 


             if(i.eq.this%ns) then
                 antiDiffFVint(:,:,:,:,this%ns) = zero
                 do j=1,this%ns-1
                    antiDiffFVint(:,:,:,:,this%ns) = antiDiffFVint(:,:,:,:,this%ns) - antiDiffFVint(:,:,:,:,j)
                 enddo
              endif

              Db = maxval(1 - this%material(i)%VF / (this%intSharp_cut)**(0.5) )
              call filter3D(this%decomp,this%gfil,Db,iflag,x_bc, y_bc,z_bc)
              call interpolateFV_6(this,Db,Db_int,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,-this%z_bc)             
              intDiff(:,:,:,:,i) = 0
 
              if(this%intSharp_msk) then
                where(VFint(:,:,:,1).gt.this%intSharp_cut)
                   Db_int(:,:,:,1) = zero
                elsewhere(VFint(:,:,:,1).lt.one-this%intSharp_cut)
                   Db_int(:,:,:,1) = zero
                endwhere
                where(VFint(:,:,:,2).gt.this%intSharp_cut)
                   Db_int(:,:,:,2)  = zero
                elsewhere(VFint(:,:,:,2).lt.one-this%intSharp_cut)
                   Db_int(:,:,:,2) = zero
                endwhere
                where(VFint(:,:,:,3).gt.this%intSharp_cut)
                   Db_int(:,:,:,3) = zero
                elsewhere(VFint(:,:,:,3).lt.one-this%intSharp_cut)
                   Db_int(:,:,:,3) = zero
                endwhere

                intDiff(:,:,:,:,i) = this%intSharp_eps*this%intSharp_dif*Db_int*gradFV_N2F

             endif
              
              !compute divergence and calculate RHS terms
              call divergenceFV(this,antiDiffFVint(:,:,:,:,i),this%material(i)%intSharp_aFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV(this,rhoiFVint(:,:,:,:,i)*antiDiffFVint(:,:,:,:,i),this%material(i)%intSharp_RFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV_6(this,intDiff(:,:,:,:,i),this%material(i)%intSharp_aDiffFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV_6(this,rhoiFVint_6(:,:,:,:,i)*intDiff(:,:,:,:,i),this%material(i)%intSharp_RDiffFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

            enddo
          

              fv_f = zero
              uDiff = zero
              do i = 1,this%ns
                 !FV term
                 fv_f = fv_f +rhoiFVint(:,:,:,:,i)*antiDiffFVint(:,:,:,:,i)
                 uDiff = uDiff + rhoiFVint_6(:,:,:,:,i)*intDiff(:,:,:,:,i)
              enddo

              call divergenceFV(this,fv_f*uFVint,this%intSharp_fFV(:,:,:,1),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV(this,fv_f*vFVint,this%intSharp_fFV(:,:,:,2),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV(this,fv_f*wFVint,this%intSharp_fFV(:,:,:,3),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
             
              call divergenceFV_6(this,uDiff*uFVint_6,this%intSharp_fDiffFV(:,:,:,1),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV_6(this,uDiff*vFVint_6,this%intSharp_fDiffFV(:,:,:,2),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV_6(this,uDiff*wFVint_6,this%intSharp_fDiffFV(:,:,:,3),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
 
               do i = 1,this%ns
                    if(this%PTeqb) then
                         call this%material(i)%get_enthalpy(hi(:,:,:,i))
                    else
                         call this%material(i)%getSpeciesDensity(rho,rhom)

                         hi(:,:,:,i) = this%material(i)%hydro%gam*(this%material(i)%p + this%material(i)%hydro%PInf) * this%material(i)%hydro%onebygam_m1

                    endif

                       call interpolateFV(this,this%material(i)%p,pFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                       hiFVint(:,:,:,:,i) = (this%material(i)%hydro%gam*pFVint(:,:,:,:,i) + this%material(i)%hydro%gam*this%material(i)%hydro%PInf)*this%material(i)%hydro%onebygam_m1
                       call interpolateFV_6(this,this%material(i)%p,pFVint_6(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                       hiFVint_6(:,:,:,:,i) = (this%material(i)%hydro%gam*pFVint_6(:,:,:,:,i) + this%material(i)%hydro%gam*this%material(i)%hydro%PInf)*this%material(i)%hydro%onebygam_m1

              enddo
             
              fv_h = zero
              fv_k = zero
              kDiff = zero
              hDiff = zero
              do i=1,this%ns

                 if(this%PTeqb) then
                   fv_h = fv_h + antiDiffFVint(:,:,:,:,i)*(rhoiFVint_local(:,:,:,:,i)*hiFVint(:,:,:,:,i))
                   fv_k = fv_k + antiDiffFVint(:,:,:,:,i)*rhoiFVint(:,:,:,:,i)*half*(uFVint**two+vFVint**two+wFVint**two)
                 else
                   fv_h = fv_h + antiDiffFVint(:,:,:,:,i)*(hiFVint(:,:,:,:,i))
                   fv_k = fv_k + antiDiffFVint(:,:,:,:,i)*rhoiFVint(:,:,:,:,i)*half*(uFVint**two+vFVint**two+wFVint**two)
                   hDiff = hDiff + intDiff(:,:,:,:,i)*hiFVint_6(:,:,:,:,i)
                   kDiff = kDiff + intDiff(:,:,:,:,i)*half*(uFVint_6**two+vFVint_6**two+wFVint_6**two)
                 endif
              enddo

              call divergenceFV(this,fv_h,this%intSharp_hFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV(this,fv_k,this%intSharp_kFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              
              call divergenceFV_6(this,hDiff,this%intSharp_hDiffFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV_6(this,kDiff,this%intSharp_kDiffFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

    end subroutine
    subroutine get_intSharp(this,rho,x_bc,y_bc,z_bc,dx,dy,dz,periodicx,periodicy,periodicz,u,v,w)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        use operators, only: divergence,gradient,filter3D, interpolateFV_x,interpolateFV_y, interpolateFV_z, gradFV_N2Fx, gradFV_N2Fy, gradFV_N2Fz
        use constants,       only: zero,epssmall,eps,one,two,third,half
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        real(rkind), intent(in) :: dx,dy,dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,u,v,w

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: norm,gradRhoYs,gradVF,gradVFdiff,fv_f,fv_h,tmp4, gradphi, fv_k, fv_p
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns) :: rhoi,VFbound,hi,spf_a,spf_r
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: phiint, gradxi,VFint,uFVint,vFVint,wFVint,gFVint,gtFVint,gpFVint,rhoFVint, kerne, gradFVFV, gradFV_N2F
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,this%ns) :: antiDiffFVint,rhoiFVint,hiFVint, rhoiFVint_local,pFVint
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,3) :: NMint,gradVF_FV,gradVFint, gradXi_FV
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: tmp,VF_fil,antiDiff,mask,RhoYsbound,filt,antiDiffFV,fmask,tanhmask,mask2,maskDiff,spf_f,spf_h,GVFmag,GVFmagT,antiDiffT, rhom
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: gradVF_x,gradVF_y, gradVF_z

        real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns) :: J_i, VF_RHS_i, Kij_coeff_i

        real(rkind) :: intSharp_alp = 0.1, r= 0.5, nmask = 40, intSharp_adm = 1.0D-1,e = 1d-100, intSharp_exp = -1.0D0,gradDiff,md1,md2 !, intSharp_tnh = 0.1 !1.0D-2

        integer :: i,j,ii,jj,kk,iflag = one,im,jm,km,k
        logical :: useTiwari = .FALSE., useRhoYsbound = .FALSE., useTotalRho = .FALSE.
        logical :: periodicx,periodicy,periodicz, useGradPsi = .FALSE., useRhoLocal = .FALSE., useHighOrder = .TRUE.,useYSbound = .TRUE., useNewSPF = .TRUE., useNewSPFfull = .FALSE.!, useTiwari = .TRUE.

        this%xi = 0.0d0
        gradxi  = 0.0d0

        !correct intSharp_tnh if negative
        if(this%intSharp_tnh.le.zero) then
           this%intSharp_tnh = one/eps
        endif



        !print method selection warnings
        if(this%intSharp_ufv.AND.this%intSharp_spf) then
           call GracefulExit("Cannot yet use finite volume and Shukla formulations together...",3457)
        endif



        !loop through species
        do i = 1,this%ns

           !component density
           rhoi(:,:,:,i) = rho*this%material(i)%Ys/this%material(i)%VF


           !calculate surface vector
           if(useGradPsi) then !use modified function -- Shukla, Pantano, Freund JCP 2010 version -- this give problems for out of bounds volume fraction -- not recommended to use
              !bounded volume fraction field
              do ii=1,this%nxp
                 do jj=1,this%nyp
                    do kk=1,this%nzp
                       VFbound(ii,jj,kk,i) = min(max(this%material(i)%VF(ii,jj,kk),zero),one)
                       !VFbound(ii,jj,kk,i) = min(max(this%material(i)%VF(ii,jj,kk),this%intSharp_cut),one-this%intSharp_cut)
                    enddo
                 enddo
              enddo

              !gradient of modified function
              call gradient(this%decomp,this%der,VFbound(:,:,:,i)**intSharp_alp/(VFbound(:,:,:,i)**intSharp_alp+(one-VFbound(:,:,:,i))**intSharp_alp),gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3)) !grad Psi -- only high order implemented
              !factor to transform modifed function to surface vector
              tmp = one/intSharp_alp * (VFbound(:,:,:,i)*(one-VFbound(:,:,:,i)))**(one-intSharp_alp) * (VFbound(:,:,:,i)**intSharp_alp+(one-VFbound(:,:,:,i))**intSharp_alp)**two !grad Phi/grad Psi

              !surface vector -- VF gradient field -- grad Phi
              gradVF(:,:,:,1) = gradVF(:,:,:,1) * tmp
              gradVF(:,:,:,2) = gradVF(:,:,:,2) * tmp
              gradVF(:,:,:,3) = gradVF(:,:,:,3) * tmp

           else if(this%usePhiForm) then

             VF_fil = this%material(i)%VF

             call filter3D(this%decomp, this%fil, VF_fil,iflag, x_bc,y_bc,z_bc)

             !! Construct xi and account for outbounds behaviour
             where( this%material(i)%VF .GE. 1-this%intSharp_cut)
                 this%xi(:,:,:,i) = this%intSharp_eps*log( ( 1-2*this%intSharp_cut + e )/ (e))*(1/(1-2*this%intSharp_cut))

             elsewhere( this%material(i)%VF .LE. this%intSharp_cut )
                 this%xi(:,:,:,i) = this%intSharp_eps*log( ( e )/ (1 -2*this%intSharp_cut + e))*(1/(1-2*this%intSharp_cut))
             elsewhere
                 this%xi(:,:,:,i) = this%intSharp_eps*(1/(1-2*this%intSharp_cut))*log( ( this%material(i)%VF - this%intSharp_cut + e )/ (1 - this%intSharp_cut - this%material(i)%VF + e) )

             endwhere
             call gradientFV(this,this%xi,gradXi_FV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

             !!!!!! UNCOMMENT FOR 2ND ORDER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if(this%intSharp_d02) then       
               call gradient(this%decomp,this%derD02,this%xi(:,:,:,i),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3)) !low order
             else
               call gradient(this%decomp,this%derCD06,this%xi(:,:,:,i),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3))
             endif
                           

    
           else !use volume fraction

              where( this%material(i)%VF .LE. this%intSharp_cut)
         
                 this%phi = this%intSharp_cut**r / ( this%intSharp_cut**r +(1-this%intSharp_cut)**r)

              elsewhere( this%material(i)%VF .GE. 1)

                 this%phi = 1

              elsewhere
                 this%phi = (this%material(i)%VF**r)/(this%material(i)%VF**r +(1-this%material(i)%VF)**r)


              endwhere


              call gradient(this%decomp,this%derD02,this%phi,gradphi(:,:,:,1),gradphi(:,:,:,2),gradphi(:,:,:,3)) !low order
              call gradient(this%decomp,this%derD02,this%material(i)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3)) !low order
              call gradient(this%decomp,this%derD06,this%material(i)%VF,gradVFdiff(:,:,:,1),gradVFdiff(:,:,:,2),gradVFdiff(:,:,:,3)) 


             ! call gradientFV(this,this%material(i)%VF,gradVFdiff,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)


           endif

           !calculate surface normal
           if (.NOT.this%intSharp_d02) then !low order FD
              !magnitude of surface vector
              GVFmag = sqrt( gradVFdiff(:,:,:,1)**two + gradVFdiff(:,:,:,2)**two + gradVFdiff(:,:,:,3)**two )

              call filter3D(this%decomp, this%gfil, gradVFdiff(:,:,:,1),iflag,x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%gfil, gradVFdiff(:,:,:,2),iflag,x_bc,y_bc,z_bc)
              call filter3D(this%decomp, this%gfil, gradVFdiff(:,:,:,3),iflag,x_bc,y_bc,z_bc)

              !surface normal
              where (GVFmag < eps)
                 norm(:,:,:,1) = zero
                 norm(:,:,:,2) = zero
                 norm(:,:,:,3) = zero
              elsewhere
                 norm(:,:,:,1) = gradVFdiff(:,:,:,1) / GVFmag
                 norm(:,:,:,2) = gradVFdiff(:,:,:,2) / GVFmag
                 norm(:,:,:,3) = gradVFdiff(:,:,:,3) / GVFmag
              endwhere

          


          else !high order FD -- high order norm with FV is stable
              !magnitude of surface vector
            


              GVFmag = sqrt( gradVF(:,:,:,1)**two + gradVF(:,:,:,2)**two + gradVF(:,:,:,3)**two )

              !surface normal
              where (GVFmag < eps)
                 norm(:,:,:,1) = zero
                 norm(:,:,:,2) = zero
                 norm(:,:,:,3) = zero
              elsewhere
                 norm(:,:,:,1) = gradVF(:,:,:,1) / GVFmag
                 norm(:,:,:,2) = gradVF(:,:,:,2) / GVFmag
                 norm(:,:,:,3) = gradVF(:,:,:,3) / GVFmag
              endwhere


           endif
           
         if (this%usePhiForm) then
           !magnitude of surface vector

          
  
         call filter3D(this%decomp, this%gfil, gradxi(:,:,:,1), iflag, x_bc, y_bc,z_bc)
         call filter3D(this%decomp, this%gfil, gradxi(:,:,:,2), iflag, x_bc,y_bc,z_bc)
         call filter3D(this%decomp, this%gfil, gradxi(:,:,:,3), iflag, x_bc,y_bc,z_bc)

          GVFmag = sqrt( gradxi(:,:,:,1)**two + gradxi(:,:,:,2)**two + gradxi(:,:,:,3)**two )

          fmask = 1 - (1 - 4*this%material(i)%VF*(1-this%material(i)%VF))**nmask

           !surface normal
           where (GVFmag < eps)
              norm(:,:,:,1) = zero
              norm(:,:,:,2) = zero
              norm(:,:,:,3) = zero
           elsewhere
              norm(:,:,:,1) = gradxi(:,:,:,1) / GVFmag
              norm(:,:,:,2) = gradxi(:,:,:,2) / GVFmag
              norm(:,:,:,3) = gradxi(:,:,:,3) / GVFmag
           endwhere


           
        endif

           ! !framework to filter surface normal if needed --- normal is not conservative so best not to?
           ! if(this%intSharp_flt) then
           !    call filter3D(this%decomp, this%gfil, norm(:,:,:,1), iflag, x_bc, y_bc,z_bc)
           !    call filter3D(this%decomp, this%gfil, norm(:,:,:,2), iflag, x_bc, y_bc,z_bc)
           !    call filter3D(this%decomp, this%gfil, norm(:,:,:,3), iflag, x_bc, y_bc,z_bc)
           ! endif

            this%norm = norm


           !2nd order finite volume discretization
           if(this%intSharp_ufv) then

              !interpolate nodes to faces: ( i, j, k ) -> ( i+1/2, j+1/2, k+1/2 )
              call interpolateFV(this,this%material(i)%VF,VFint,periodicx,periodicy,periodicz,      this%x_bc, this%y_bc, this%z_bc)
              !TODO: make sure these BCS for surface normal are correct
             
              call interpolateFV(this,norm(:,:,:,1),NMint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc, this%y_bc, this%z_bc)
              call interpolateFV(this,norm(:,:,:,2),NMint(:,:,:,:,2),periodicx,periodicy,periodicz, this%x_bc,-this%y_bc, this%z_bc)
              call interpolateFV(this,norm(:,:,:,3),NMint(:,:,:,:,3),periodicx,periodicy,periodicz, this%x_bc, this%y_bc,-this%z_bc)
             

              call interpolateFV(this,u,uFVint,periodicx,periodicy,periodicz,-this%x_bc, this%y_bc, this%z_bc)
              call interpolateFV(this,v,vFVint,periodicx,periodicy,periodicz, this%x_bc,-this%y_bc, this%z_bc)
              call interpolateFV(this,w,wFVint,periodicx,periodicy,periodicz, this%x_bc, this%y_bc,-this%z_bc)
              call interpolateFV(this,rho,rhoFVint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              if(useRhoLocal) then !use local component density --- recommended
                 call interpolateFV(this,rhoi(:,:,:,i),rhoiFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              else
                 call interpolateFV(this,rhoi(:,:,:,i)*zero+this%material(i)%elastic%rho0,rhoiFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              endif

              call interpolateFV(this,rhoi(:,:,:,i),rhoiFVint_local(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              if(this%intSharp_d02) then !which FD VF gradient to interpolate
                 !not used
                 !TODO: make sure these BCS for gradVF are correct
                 call interpolateFV(this,gradVF(:,:,:,1),gradVFint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc, this%y_bc, this%z_bc)
                 call interpolateFV(this,gradVF(:,:,:,2),gradVFint(:,:,:,:,2),periodicx,periodicy,periodicz, this%x_bc,-this%y_bc, this%z_bc) 
                 call interpolateFV(this,gradVF(:,:,:,3),gradVFint(:,:,:,:,3),periodicx,periodicy,periodicz, this%x_bc, this%y_bc,-this%z_bc) 
    
	      else
                 !used in high order version
                 !TODO: make sure these BCS for gradVF are correct
                 

                 call interpolateFV(this,gradVFdiff(:,:,:,1),gradVFint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc, this%y_bc, this%z_bc)
                 call interpolateFV(this,gradVFdiff(:,:,:,2),gradVFint(:,:,:,:,2),periodicx,periodicy,periodicz, this%x_bc,-this%y_bc, this%z_bc) 
                 call interpolateFV(this,gradVFdiff(:,:,:,3),gradVFint(:,:,:,:,3),periodicx,periodicy,periodicz, this%x_bc, this%y_bc,-this%z_bc) 
              endif

              
              this%normFV(:,:,:,1) = NMint(:,:,:,1,1)
              this%normFV(:,:,:,2) = NMint(:,:,:,2,2)
              this%normFV(:,:,:,3) = NMint(:,:,:,3,3)

            

              !calculate antiDiffFVint term
              antiDiffFVint(:,:,:,1,i) = -this%intSharp_gam * (VFint(:,:,:,1)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,1))*NMint(:,:,:,1,1)
              antiDiffFVint(:,:,:,2,i) = -this%intSharp_gam * (VFint(:,:,:,2)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,2))*NMint(:,:,:,2,2)
              antiDiffFVint(:,:,:,3,i) = -this%intSharp_gam * (VFint(:,:,:,3)-this%intSharp_cut)*(one-this%intSharp_cut-VFint(:,:,:,3))*NMint(:,:,:,3,3)

              this%antidiff = antiDiffFVint(:,:,:,1,i)
 

!              call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,1, i), iflag, x_bc,y_bc,z_bc)
!              call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,2, i), iflag, x_bc,y_bc,z_bc)
!              call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,3, i), iflag, x_bc,y_bc,z_bc)


              if (this%usePhiForm) then 
                 call interpolateFV(this,this%xi,phiint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc, this%z_bc)
                 antiDiffFVint(:,:,:,1,i) = -this%intSharp_gam*(0.25*(1 - (tanh((1-2*this%intSharp_cut)* phiint(:,:,:,1)/(2 * this%intSharp_eps)))**2)-0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,1)/(2 * this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,1,1)
                 antiDiffFVint(:,:,:,2,i) = -this%intSharp_gam*(0.25*(1 - (tanh((1-2*this%intSharp_cut)* phiint(:,:,:,2)/(2 * this%intSharp_eps)))**2)-0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,2)/(2 * this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,2,2)
                 antiDiffFVint(:,:,:,3,i) = -this%intSharp_gam*(0.25*(1 - (tanh((1-2*this%intSharp_cut)* phiint(:,:,:,3)/(2 * this%intSharp_eps)))**2)-0.5*(1+tanh((1-2*this%intSharp)*phiint(:,:,:,3)/(2 * this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut)*NMint(:,:,:,3,3)


                 call divergenceFV(this,antiDiffFVint(:,:,:,:,i),this%antidiff,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

              endif

             !call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,1, i),iflag, x_bc,y_bc,z_bc)
             !call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,2, i),iflag, x_bc,y_bc,z_bc)
             !call filter3D(this%decomp, this%fil, antiDiffFVint(:,:,:,3, i),iflag, x_bc,y_bc,z_bc)



              if (.NOT.this%intSharp_msk) then
                 ! ! Do not do this for FV
                 ! !Mask blending function matches FD but now on faces
                 ! antiDiffFVint(:,:,:,:,i) = antiDiffFVint(:,:,:,:,i) * tanh( ( ((VFint-this%intSharp_cut)*(one-VFint-this%intSharp_cut)) / (this%intSharp_cut/this%intSharp_tnh**two) )**two )
              else

                 !FV diffusion terms
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !call gradientFV(this,this%material(i)%VF,gradVF_FV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                 !gradFVFV(:,:,:,1) = gradVF_FV(:,:,:,1,1)
                 !gradFVFV(:,:,:,2) = gradVF_FV(:,:,:,2,2)
                 !gradFVFV(:,:,:,3) = gradVF_FV(:,:,:,3,3)
 
                 call gradFV_N2Fx(this%decomp,this%derStaggd02,this%material(i)%VF,gradFV_N2F(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
                 call gradFV_N2Fy(this%decomp,this%derStaggd02,this%material(i)%VF,gradFV_N2F(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
                 call gradFV_N2Fz(this%decomp,this%derStaggd02,this%material(i)%VF,gradFV_N2F(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

                 !call divergenceFV(this,gradFVFV,this%intdiff,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

                ! antiDiffFVint(:,:,:,1,i) = antiDiffFVint(:,:,:,1,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,1))
                ! antiDiffFVint(:,:,:,2,i) = antiDiffFVint(:,:,:,2,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,2))
                ! antiDiffFVint(:,:,:,3,i) = antiDiffFVint(:,:,:,3,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,3))

                 !antiDiffFVint(:,:,:,1,i) = antiDiffFVint(:,:,:,1,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,1,1))
                 !antiDiffFVint(:,:,:,2,i) = antiDiffFVint(:,:,:,2,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,2,2))
                 !antiDiffFVint(:,:,:,3,i) = antiDiffFVint(:,:,:,3,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,3,3))
                 !this%gradVF_FV =  gradVF_FV



              endif

              !FV mask term for antidiffusion
              where(VFint(:,:,:,1).lt.this%intSharp_cut)
                 antiDiffFVint(:,:,:,1,i) = zero
              elsewhere(VFint(:,:,:,1).gt.one-this%intSharp_cut)
                 antiDiffFVint(:,:,:,1,i) = zero
             endwhere
              where(VFint(:,:,:,2).lt.this%intSharp_cut)
                 antiDiffFVint(:,:,:,2,i) = zero
              elsewhere(VFint(:,:,:,2).gt.one-this%intSharp_cut)
                 antiDiffFVint(:,:,:,2,i) = zero
              endwhere
              where(VFint(:,:,:,3).lt.this%intSharp_cut)
                 antiDiffFVint(:,:,:,3,i) = zero
              elsewhere(VFint(:,:,:,3).gt.one-this%intSharp_cut)
                 antiDiffFVint(:,:,:,3,i) = zero
              endwhere
            


              if (.NOT.this%intSharp_msk) then
                 !FV diffusion terms
                 !!!!!!!!!!!!!!!!!!!!!!!!!!! UNCOMMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !call gradientFV(this,this%material(i)%VF,gradVF_FV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                 !antiDiffFVint(:,:,:,1,i) = antiDiffFVint(:,:,:,1,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,1,1))
                 !antiDiffFVint(:,:,:,2,i) = antiDiffFVint(:,:,:,2,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,2,2))
                 !antiDiffFVint(:,:,:,3,i) = antiDiffFVint(:,:,:,3,i) + this%intSharp_gam * (this%intSharp_eps * gradVF_FV(:,:,:,3,3))

                 call gradFV_N2Fx(this%decomp,this%derStagg,this%material(i)%VF,gradFV_N2F(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
                 call gradFV_N2Fy(this%decomp,this%derStagg,this%material(i)%VF,gradFV_N2F(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
                 call gradFV_N2Fz(this%decomp,this%derStagg,this%material(i)%VF,gradFV_N2F(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

                 !call
                 !divergenceFV(this,gradFVFV,this%intdiff,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

                 antiDiffFVint(:,:,:,1,i) = antiDiffFVint(:,:,:,1,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,1))
                 antiDiffFVint(:,:,:,2,i) = antiDiffFVint(:,:,:,2,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,2))
                 antiDiffFVint(:,:,:,3,i) = antiDiffFVint(:,:,:,3,i) + this%intSharp_gam * (this%intSharp_eps * gradFV_N2F(:,:,:,3))

              endif

              !ENFORCE antiDiffFVint sum to zero over i species
              if(i.eq.this%ns) then
                 antiDiffFVint(:,:,:,:,this%ns) = zero
                 do j=1,this%ns-1
                    antiDiffFVint(:,:,:,:,this%ns) = antiDiffFVint(:,:,:,:,this%ns) - antiDiffFVint(:,:,:,:,j)
                 enddo
              endif

              !compute divergence and calculate RHS terms
              call divergenceFV(this,antiDiffFVint(:,:,:,:,i),this%material(i)%intSharp_aFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
              call divergenceFV(this,rhoiFVint(:,:,:,:,i)*antiDiffFVint(:,:,:,:,i),this%material(i)%intSharp_RFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

           else
              !set FV terms to zero if not using
              this%material(i)%intSharp_aFV = zero
              this%material(i)%intSharp_RFV = zero

           endif

           !antidiffusion term
           antiDiff = (this%material(i)%VF-this%intSharp_cut)*(one-this%material(i)%VF-this%intSharp_cut)

           !mask for interface region
           mask = one
           where(this%material(i)%VF.lt.this%intSharp_cut)
              mask = zero
           elsewhere(this%material(i)%VF.gt.one-this%intSharp_cut)
              mask = zero
           endwhere
           

        !   if(this%usePhiForm) then

        !      mask = one
        !      where(this%xi .le.  this%intSharp_eps*log( ( e )/ (1-2*this%intSharp_cut+e))*(1/(1-2*this%intSharp_cut)) )

        !        mask = zero

        !      elsewhere(this%xi .ge. this%intSharp_eps*log( (1-2*this%intSharp_cut + e )/ (e))*(1/(1-2*this%intSharp_cut)) )

        !        mask = zero

        !      endwhere
                    
        !    endif



           this%VFboundDiff(:,:,:,i) = zero
           maskDiff = zero
           ! if(this%intSharp_ufv) then
           md1 = this%intSharp_cut**half !need to set a larger cutoff toeffectively blend FV and FD
           ! else
           !    md1 = this%intSharp_cut
           ! endif

           do ii=1,this%nxp
              do jj=1,this%nyp
                 do kk=1,this%nzp
                    ! this%VFboundDiff(ii,jj,kk,i) =
                    ! max(this%VFboundDiff(ii,jj,kk,i),
                    ! (this%intSharp_cut-this%material(i)%VF(ii,jj,kk))/this%intSharp_cut)
                    ! ! VF < cut
                    ! this%VFboundDiff(ii,jj,kk,i) =
                    ! max(this%VFboundDiff(ii,jj,kk,i),
                    ! (this%material(i)%VF(ii,jj,kk)-(one-this%intSharp_cut))/this%intSharp_cut)
                    ! !VF > 1-cut
                    this%VFboundDiff(ii,jj,kk,i) =max(this%VFboundDiff(ii,jj,kk,i), (md1-this%material(i)%VF(ii,jj,kk))/(md1)) !VF < sqrt(cut)
                    this%VFboundDiff(ii,jj,kk,i) =max(this%VFboundDiff(ii,jj,kk,i),(this%material(i)%VF(ii,jj,kk)-(one-md1))/(md1))  !VF > 1-sqrt(cut)
                    useYSbound = .FALSE.



         !          if(this%usePhiForm) then

          !            this%VFboundDiff(ii,jj,kk,i) = max(this%VFboundDiff(ii,jj,kk,i), (md1-0.5*(1-this%intSharp_cut)*(1+tanh( (1-2*this%intSharp_cut)*this%xi(ii,jj,kk) /this%intSharp_eps) ))/(md1))
!VF < sqrt(cut)
         !             this%VFboundDiff(ii,jj,kk,i) = max(this%VFboundDiff(ii,jj,kk,i),(0.5*(1-this%intSharp_cut)*(1+tanh( (1-2*this%intSharp_cut)*this%xi(ii,jj,kk) /this%intSharp_eps) )-(one-md1))/(md1))
!VF > 1-sqrt(cut)

          !         endif


                    if(useYSbound) then !possible useful for high density ratioow to get bounds?
                       ! this%VFboundDiff(ii,jj,kk,i) =
                       ! max(this%VFboundDiff(ii,jj,kk,i),
                       ! (this%material(i)%intSharp_ysc(1)-this%material(i)%Ys(ii,jj,kk))/this%material(i)%intSharp_ysc(1))
                       ! ! Ys < YScut(i,1)
                       ! this%VFboundDiff(ii,jj,kk,i) =
                       ! max(this%VFboundDiff(ii,jj,kk,i),
                       ! ((one-this%material(i)%intSharp_ysc(2))-(one-this%material(i)%Ys(ii,jj,kk)))/(one-this%material(i)%intSharp_ysc(2)))
                       ! !Ys > 1-Yscut(i,2)
                       ! this%VFboundDiff(ii,jj,kk,i) =
                       ! max(this%VFboundDiff(ii,jj,kk,i),
                       ! (zero-this%material(i)%Ys(ii,jj,kk))/this%material(i)%intSharp_ysc(1))
                       ! ! Ys < YScut(i,1)
                       ! this%VFboundDiff(ii,jj,kk,i) =
                       ! max(this%VFboundDiff(ii,jj,kk,i),
                       ! ((one-one)-(one-this%material(i)%Ys(ii,jj,kk)))/(one-this%material(i)%intSharp_ysc(2)))
                       ! !Ys > 1-Yscut(i,2)
                       this%VFboundDiff(ii,jj,kk,i) =max(this%VFboundDiff(ii,jj,kk,i),(zero-this%material(i)%Ys(ii,jj,kk))/this%intSharp_cut) ! Ys < YScut(i,1)
                       this%VFboundDiff(ii,jj,kk,i) =max(this%VFboundDiff(ii,jj,kk,i),((one-one)-(one-this%material(i)%Ys(ii,jj,kk)))/this%intSharp_cut)  !Ys >1-Yscut(i,2)

                    endif
                 enddo
              enddo
           enddo


           !blending function between masked and unmasked -- matches version
           !used on faces for FV -- tnh controls transition width
           !tanhmask = tanh( (
           !((this%material(i)%VF-this%intSharp_cut)*(one-this%material(i)%VF-this%intSharp_cut))
           !/ (this%intSharp_cut/this%intSharp_tnh**two) )**two ) !old
           tanhmask = tanh( (((this%material(i)%VF-this%intSharp_cut)*(one-this%material(i)%VF-this%intSharp_cut))/ this%intSharp_tnh )**two ) !new

           if(this%intSharp_spf) then 
              antiDiff = antiDiff * mask * tanhmask
           else
              antiDiff = antiDiff * mask
           endif


           !Shukla, Pantano, Freund JCP 2010 version -- not in divergence form
           if(this%intSharp_spf) then 
              this%material(i)%intSharp_a = zero !components 2 and 3 are zerohen not in divergence form -- they are not used in logic for RHS update
              this%material(i)%intSharp_R = zero !components 2 and 3 are zerowhen not in divergence form -- they are not used in logic for RHS update
              
              !if(.not.useTiwari) then
              if(.not.this%intSharp_utw) then
                 !RHS gradient
                 if(this%intSharp_msk) then
                    !original
                    spf_a(:,:,:,i) = this%intSharp_eps*GVFmag*mask*tanhmask-antiDiff !GVFmag is gradVF magnitude
                    
                    !spf_a(:,:,:,i) = this%intSharp_eps*GVFmag*mask-antiDiff
                    !!GVFmag is gradVF magnitude
                    
                    
                    ! !spf_a(:,:,:,i) = this%intSharp_eps*GVFmag*mask-antiDiff
                    ! !GVFmag is gradVF magnitude
                    
                    ! !new -- for VFboundDiff
                    ! this%VFboundDiff(:,:,:,i) =
                    ! this%intSharp_dif*this%VFboundDiff(:,:,:,i)
                    ! call filter3D(this%decomp, this%gfil,
                    ! this%VFboundDiff(:,:,:,i), iflag, x_bc, y_bc,z_bc) 
                    ! !call
                    ! gradient(this%decomp,this%der,this%material(i)%VF*this%VFboundDiff(:,:,:,i),gradVFdiff(:,:,:,1),gradVFdiff(:,:,:,2),gradVFdiff(:,:,:,3))
                    ! !high order
                    ! !GVFmag = sqrt( gradVFdiff(:,:,:,1)**two +
                    ! gradVFdiff(:,:,:,2)**two + gradVFdiff(:,:,:,3)**two )
                    
                    ! spf_a(:,:,:,i) = spf_a(:,:,:,i) +
                    ! this%intSharp_eps*GVFmag*this%VFboundDiff(:,:,:,i)
                 else
                    spf_a(:,:,:,i) = this%intSharp_eps*GVFmag-antiDiff
                 endif
                 
                 !filter argument of RHS gradient
                 if(this%intSharp_flt) then
                    call filter3D(this%decomp, this%fil, spf_a(:,:,:,i), iflag,x_bc, y_bc,z_bc)
                 endif
                 
                 
                 !update terms for FV --- SPF version
                 if (this%intSharp_d02)  then
                    call gradient(this%decomp,this%derD02,spf_a(:,:,:,i),tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
                 else
                    call gradient(this%decomp,this%der,spf_a(:,:,:,i),tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
                 endif
                 !update VF term -- only a used for this update -- aDiff used
                 !for divergence form diffusion
                 this%material(i)%intSharp_a(:,:,:,1) = this%intSharp_gam *(norm(:,:,:,1)*tmp4(:,:,:,1) + norm(:,:,:,2)*tmp4(:,:,:,2) +norm(:,:,:,3)*tmp4(:,:,:,3))
                 spf_r(:,:,:,i) = spf_a(:,:,:,i)*rhoi(:,:,:,i)
                 
                 !put rhoi in gradient instead of assuming constant
                 if(useNewSPF) then
                    
                    if(this%intSharp_d02)  then
                       call gradient(this%decomp,this%derD02,spf_r(:,:,:,i),tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
                    else
                       call gradient(this%decomp,this%der,spf_r(:,:,:,i),tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
                    endif
                    this%material(i)%intSharp_R(:,:,:,1) = this%intSharp_gam *(norm(:,:,:,1)*tmp4(:,:,:,1) + norm(:,:,:,2)*tmp4(:,:,:,2) +norm(:,:,:,3)*tmp4(:,:,:,3))

                 endif
              endif


              !Use Tiwari Freund Pantano, JCP 2013 version
              !if(useTiwari.and.useNewSPF) then
              if(this%intSharp_utw.and.useNewSPF) then
                 spf_r(:,:,:,i) = zero ! can't do full couple with Tiwari -- notactually in gradient form
                 
                 if (this%intSharp_d02) then
                    call gradient(this%decomp,this%derD02,rho*this%material(i)%Ys,gradRhoYs(:,:,:,1),gradRhoYs(:,:,:,2),gradRhoYs(:,:,:,3))
                 else
                    call gradient(this%decomp,this%der,rho*this%material(i)%Ys,gradRhoYs(:,:,:,1),gradRhoYs(:,:,:,2),gradRhoYs(:,:,:,3))
                 endif

                 tmp = this%intSharp_eps*(norm(:,:,:,1)*gradRhoYs(:,:,:,1) + norm(:,:,:,2)*gradRhoYs(:,:,:,2) + norm(:,:,:,3)*gradRhoYs(:,:,:,3)) 

                 if (this%intSharp_d02) then
                    call gradient(this%decomp,this%derD02,tmp,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
                 else
                    call gradient(this%decomp,this%der,tmp,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
                 endif

                 if(this%intSharp_msk) then
                    do j=1,3
                       tmp4(:,:,:,j) = ( tmp4(:,:,:,j) -(one-two*this%material(i)%VF)*gradRhoYs(:,:,:,j) ) * mask*tanhmask
                    enddo
                 else
                    do j=1,3
                       tmp4(:,:,:,j) =   tmp4(:,:,:,j) -(one-two*this%material(i)%VF)*gradRhoYs(:,:,:,j)   * mask*tanhmask
                    enddo
                 endif

                 this%material(i)%intSharp_R(:,:,:,1) = this%intSharp_gam *(norm(:,:,:,1)*tmp4(:,:,:,1) + norm(:,:,:,2)*tmp4(:,:,:,2) +norm(:,:,:,3)*tmp4(:,:,:,3)) 
                 
                 if(this%intSharp_flt) then
                    call filter3D(this%decomp, this%fil,this%material(i)%intSharp_R(:,:,:,1), iflag, x_bc, y_bc,z_bc)
                 endif
                 
              endif


              !out of bounds diffusion is divergence form always
              if(this%intSharp_msk) then
                 !Divergence form -- high order --- high order --- FD diffusion
                 !away from interface
                 do j=1,3
                    !this%material(i)%intSharp_aDiff(:,:,:,j) =
                    !this%intSharp_gam*this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)


                    !this%material(i)%intSharp_aDiff(:,:,:,j) =
                    !this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)
                    !!tested
                    this%material(i)%intSharp_aDiff(:,:,:,j) =this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)*this%intSharp_dif
!new -- same as divergence version
                    !Gaussian filter to smooth out-of-bounds diffusion --
                    !similar to LAD
                    call filter3D(this%decomp, this%gfil,this%material(i)%intSharp_aDiff(:,:,:,j), iflag, x_bc, y_bc,z_bc)
                 enddo
              else
                 this%material(i)%intSharp_aDiff = zero
              endif
              
              
              

              ! old Tiwari
              ! ! The basic framework is here, but it is now deprecated -- may
              ! have some advantages in not requiering VF field -- but if the
              ! pressure-temperature equilibration stability problems are fixed,
              ! the other version would be equivalent if not better since
              ! compoent rho is not assumed constant across interfaces in the
              ! other version

              ! if(useTiwari) then

              !    if(useRhoYsbound) then                  
              !       !anti diffusion term and mask --- for intSharp_cut > 0
              !       RhoYsbound = rho*this%material(i)%Ys ! =
              !       rhoi(:,:,:,i)*VFbound(:,:,:,i)
              !       where (this%material(i)%VF .lt. this%intSharp_cut)
              !          RhoYsbound = rhoi(:,:,:,i) * this%intSharp_cut
              !       elsewhere (this%material(i)%VF .gt. one-this%intSharp_cut)
              !          RhoYsbound = rhoi(:,:,:,i) * (one-this%intSharp_cut)
              !       endwhere

              !       filt = RhoYsbound
              !    else
              !       filt = rho*this%material(i)%Ys
              !    endif

              !    if(this%intSharp_flt) then
              !       call filter3D(this%decomp, this%fil, filt, iflag, x_bc,
              !       y_bc,z_bc)
              !    endif

              !    call
              !    gradient(this%decomp,this%der,filt,gradRhoYs(:,:,:,1),gradRhoYs(:,:,:,2),gradRhoYs(:,:,:,3))

              !    filt = this%intSharp_eps * ( norm(:,:,:,1)*gradRhoYs(:,:,:,1)
              !    + norm(:,:,:,2)*gradRhoYs(:,:,:,2) +
              !    norm(:,:,:,3)*gradRhoYs(:,:,:,3) )

              !    if(this%intSharp_flt) then
              !       call filter3D(this%decomp, this%fil, filt, iflag, x_bc,
              !       y_bc,z_bc)
              !    endif

              !    call
              !    gradient(this%decomp,this%der,filt,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))

              !    gradVF(:,:,:,1) = gradVF(:,:,:,1) -
              !    (one-two*VFbound(:,:,:,i))/(one-this%intSharp_eps) *
              !    gradRhoYs(:,:,:,1)
              !    gradVF(:,:,:,2) = gradVF(:,:,:,2) -
              !    (one-two*VFbound(:,:,:,i))/(one-this%intSharp_eps) *
              !    gradRhoYs(:,:,:,2)
              !    gradVF(:,:,:,3) = gradVF(:,:,:,3) -
              !    (one-two*VFbound(:,:,:,i))/(one-this%intSharp_eps) *
              !    gradRhoYs(:,:,:,3)

              !    filt = this%intSharp_gam * (norm(:,:,:,1)*gradVF(:,:,:,1) +
              !    norm(:,:,:,2)*gradVF(:,:,:,2) +
              !    norm(:,:,:,3)*gradVF(:,:,:,3)) * mask *
              !    tanh((antiDiff/this%intSharp_tnh)**two)

              !    if(this%intSharp_flt) then
              !       call filter3D(this%decomp, this%fil, filt, iflag, x_bc,
              !       y_bc,z_bc)
              !    endif

              !    this%material(i)%intSharp_R(:,:,:,1) = filt

              !    !set components 2 and 3 to zero when not in divergence form
              !    -- not used
              !    this%material(i)%intSharp_R(:,:,:,2) = zero
              !    this%material(i)%intSharp_R(:,:,:,3) = zero


              ! endif

           else !Divergence form of Jain-Mani-Moin

              this%material(i)%intSharp_a = zero
              this%material(i)%intSharp_aDiff = zero

              if(this%intSharp_msk) then
                 this%VFboundDiff(:,:,:,i) =this%intSharp_dif*this%VFboundDiff(:,:,:,i)  
                 mask2 = mask
               !  mask2 = mask*tanhmask      !original        
              else
                 this%VFboundDiff(:,:,:,i) = zero
                 mask2 = one !too diffusive for FD !original
                 !mask2 = mask
              endif

          !form low and high order FD diffusion terms for VF
              if(this%intSharp_ufv) then

                 !high order --- high order --- FD diffusion away from interface
                 this%material(i)%intSharp_a = zero
                 do j=1,3
                    !this%material(i)%intSharp_aDiff(:,:,:,j) =
                    !this%intSharp_gam*this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)

                    this%material(i)%intSharp_aDiff(:,:,:,j) =this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)

                    !Gaussian filter to smooth out-of-bounds diffusion --
                    !similar to LAD
                    call filter3D(this%decomp, this%gfil,this%material(i)%intSharp_aDiff(:,:,:,j), iflag, x_bc, y_bc,z_bc)
                 enddo


              else

                 if(.NOT.this%intSharp_d02) then !high order FD
                    !high order --- high order --- FD diffusion away from
                    !interface
                    this%material(i)%intSharp_a = zero
                    do j=1,3
                       !this%material(i)%intSharp_aDiff(:,:,:,j) =
                       !this%intSharp_gam*this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)
                       this%material(i)%intSharp_aDiff(:,:,:,j) =this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)
                       !Gaussian filter to smooth out-of-bounds diffusion --
                       !similar to LAD
                       call filter3D(this%decomp, this%gfil,this%material(i)%intSharp_aDiff(:,:,:,j), iflag, x_bc, y_bc,z_bc)
                    enddo


                    if(.NOT.useHighOrder) then ! low order for interface terms
                       !low order -- interface sharpening and coupled diffusion
                       !-- this is unstable -- high order gradVF and norm are
                       !the problem
                       do j=1,3
                          this%material(i)%intSharp_a(:,:,:,j) =this%material(i)%intSharp_a(:,:,:,j) + this%intSharp_gam*(this%intSharp_eps*gradVFdiff(:,:,:,j)*mask2 - antiDiff *norm(:,:,:,j))
                       enddo
                    else !high order for interface terms --- this is unstable
                       !high order -- interface sharpening and coupled diffusion
                       do j=1,3
                          this%material(i)%intSharp_aDiff(:,:,:,j) =this%material(i)%intSharp_aDiff(:,:,:,j) +this%intSharp_gam*(this%intSharp_eps*gradVFdiff(:,:,:,j)*mask2 - antiDiff *norm(:,:,:,j))
                       enddo
                    endif

                 else !low order FD
                    !low order -- diffusion away from interface
                    do j=1,3
                       !this%material(i)%intSharp_a(:,:,:,j) =
                       !this%intSharp_gam*this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)
                       this%material(i)%intSharp_a(:,:,:,j) =this%intSharp_eps*gradVFdiff(:,:,:,j)*this%VFboundDiff(:,:,:,i)
                       !Gaussian filter to smooth out-of-bounds diffusion --
                       !similar to LAD
                       call filter3D(this%decomp, this%gfil,this%material(i)%intSharp_a(:,:,:,j), iflag, x_bc, y_bc,z_bc)
                    enddo

                    !low order -- interface sharpening and coupled diffusion
                    do j=1,3
                       this%material(i)%intSharp_a(:,:,:,j) =this%material(i)%intSharp_a(:,:,:,j) +this%intSharp_gam*(this%intSharp_eps*gradVF(:,:,:,j)*mask2 - antiDiff *norm(:,:,:,j))
                    enddo
                 endif

              endif
           endif


           !filter conservative fluxes and updates
           if(this%intSharp_flt) then
              call filter3D(this%decomp, this%fil,this%material(i)%intSharp_a(:,:,:,1), iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil,this%material(i)%intSharp_a(:,:,:,2), iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil,this%material(i)%intSharp_a(:,:,:,3), iflag, x_bc, y_bc,z_bc)

              call filter3D(this%decomp, this%fil,this%material(i)%intSharp_aDiff(:,:,:,1), iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil,this%material(i)%intSharp_aDiff(:,:,:,2), iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil,this%material(i)%intSharp_aDiff(:,:,:,3), iflag, x_bc, y_bc,z_bc)

              call filter3D(this%decomp, this%fil,this%material(i)%intSharp_aFV, iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil,this%material(i)%intSharp_RFV, iflag, x_bc, y_bc,z_bc)
           endif

        enddo


        if(this%intSharp_spf.and.useNewSPF.and.useNewSPFfull) then
!useNewSPFfull leads to interface instability in coupling terms
           !if(useTiwari) then
           if(this%intSharp_utw) then
              print*, "can't use Tiwari here"
              stop
           endif
        !   !enfoce surface normals sum to zero --- correct intSharp_a for
        !   material ns --- this should be modifed for ns > 2
        !   !i.e. conserve intSharp_a intSharp__r

        !   ! do i = 1,this%ns
        !   !    this%material(i)%intSharp_RDiff(:,:,:,1) = rhoi(:,:,:,i) *
        !   this%material(i)%intSharp_aDiff(:,:,:,1) !this is still used in
        !   divergence form for OOB diffusion
        !   !    this%material(i)%intSharp_RDiff(:,:,:,2) = rhoi(:,:,:,i) *
        !   this%material(i)%intSharp_aDiff(:,:,:,2)
        !   !    this%material(i)%intSharp_RDiff(:,:,:,3) = rhoi(:,:,:,i) *
        !   this%material(i)%intSharp_aDiff(:,:,:,3)
        !   ! enddo



        !   !enfoce surface normals sum to zero --- correct intSharp_a for
        !   material ns --- this should be modifed for ns > 2

        !   !out-of-bounds diffusion terms
        !   tmp4 = zero
        !   do j = 1,this%ns-1
        !      tmp4(:,:,:,1) = tmp4(:,:,:,1) +
        !      this%material(j)%intSharp_aDiff(:,:,:,1)
        !      tmp4(:,:,:,2) = tmp4(:,:,:,2) +
        !      this%material(j)%intSharp_aDiff(:,:,:,2)
        !      tmp4(:,:,:,3) = tmp4(:,:,:,3) +
        !      this%material(j)%intSharp_aDiff(:,:,:,3)
        !   enddo
        !   this%material(this%ns)%intSharp_aDiff(:,:,:,1) = -tmp4(:,:,:,1)
        !   this%material(this%ns)%intSharp_aDiff(:,:,:,2) = -tmp4(:,:,:,2)
        !   this%material(this%ns)%intSharp_aDiff(:,:,:,3) = -tmp4(:,:,:,3)

        !   !mass flux --  rho*Ys flux --  R_i
        !   do i = 1,this%ns
        !      if(useRhoLocal) then !use local component density --- recommended
        !         this%material(i)%intSharp_RDiff(:,:,:,1) = rhoi(:,:,:,i) *
        !         this%material(i)%intSharp_aDiff(:,:,:,1)
        !         this%material(i)%intSharp_RDiff(:,:,:,2) = rhoi(:,:,:,i) *
        !         this%material(i)%intSharp_aDiff(:,:,:,2)
        !         this%material(i)%intSharp_RDiff(:,:,:,3) = rhoi(:,:,:,i) *
        !         this%material(i)%intSharp_aDiff(:,:,:,3)
        !      else
        !         this%material(i)%intSharp_RDiff(:,:,:,1) =
        !         this%material(i)%elastic%rho0 *
        !         this%material(i)%intSharp_aDiff(:,:,:,1)
        !         this%material(i)%intSharp_RDiff(:,:,:,2) =
        !         this%material(i)%elastic%rho0 *
        !         this%material(i)%intSharp_aDiff(:,:,:,2)
        !         this%material(i)%intSharp_RDiff(:,:,:,3) =
        !         this%material(i)%elastic%rho0 *
        !         this%material(i)%intSharp_aDiff(:,:,:,3)
        !      endif
        !   enddo


        !   ! momentum flux -- f_i
        !   this%intSharp_f = zero
        !   this%intSharp_fDiff = zero

        !   !momentum term is nonzero if coupling is turned on
        !   if(this%intSharp_cpl) then
        !      !spf update term

        !      !net mass flux
        !      spf_f = zero
        !      do i = 1,this%ns
        !         spf_f = spf_f + spf_r(:,:,:,i)
        !      enddo

        !      ! u
        !      if (this%intSharp_d02)  then
        !         call
        !         gradient(this%decomp,this%derD02,spf_f*u,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
        !      else
        !         call
        !         gradient(this%decomp,this%der,spf_f*u,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
        !      endif
        !      this%intSharp_f(:,:,:,1) = this%intSharp_gam *
        !      (norm(:,:,:,1)*tmp4(:,:,:,1) + norm(:,:,:,2)*tmp4(:,:,:,2) +
        !      norm(:,:,:,3)*tmp4(:,:,:,3))

        !      ! v
        !      if (this%intSharp_d02)  then
        !         call
        !         gradient(this%decomp,this%derD02,spf_f*v,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
        !      else
        !         call
        !         gradient(this%decomp,this%der,spf_f*v,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
        !      endif
        !      this%intSharp_f(:,:,:,2) = this%intSharp_gam *
        !      (norm(:,:,:,1)*tmp4(:,:,:,1) + norm(:,:,:,2)*tmp4(:,:,:,2) +
        !      norm(:,:,:,3)*tmp4(:,:,:,3))

        !      ! w
        !      if (this%intSharp_d02)  then
        !         call
        !         gradient(this%decomp,this%derD02,spf_f*w,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
        !      else
        !         call
        !         gradient(this%decomp,this%der,spf_f*w,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
        !      endif
        !      this%intSharp_f(:,:,:,3) = this%intSharp_gam *
        !      (norm(:,:,:,1)*tmp4(:,:,:,1) + norm(:,:,:,2)*tmp4(:,:,:,2) +
        !      norm(:,:,:,3)*tmp4(:,:,:,3))


        !      do i = 1,this%ns
        !         !high order FD terms
        !         this%intSharp_fDiff(:,:,:,1) = this%intSharp_fDiff(:,:,:,1) +
        !         this%material(i)%intSharp_RDiff(:,:,:,1)
        !         this%intSharp_fDiff(:,:,:,2) = this%intSharp_fDiff(:,:,:,2) +
        !         this%material(i)%intSharp_RDiff(:,:,:,2)  
        !         this%intSharp_fDiff(:,:,:,3) = this%intSharp_fDiff(:,:,:,3) +
        !         this%material(i)%intSharp_RDiff(:,:,:,3)
        !      enddo
        !   endif


        !   !enthalpy flux -- (intSharp_h)_j = sum_i ( rho_i h_i (a_i)_j)
        !   this%intSharp_h = zero
        !   this%intSharp_hDiff = zero

        !   !enthalpy term is nonzero if coupling is turned on
        !   if(this%intSharp_cpl) then
        !      !spf update term 
        !      spf_h = zero
        !      do i = 1,this%ns
        !         call this%material(i)%get_enthalpy(hi(:,:,:,i))
        !         spf_h = spf_h + spf_r(:,:,:,i)*(hi(:,:,:,i) + half*(u**two +
        !         v**two + w**two) )
        !      enddo
        !      if (this%intSharp_d02)  then
        !         call
        !         gradient(this%decomp,this%derD02,spf_h,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
        !      else
        !         call
        !         gradient(this%decomp,this%der,spf_h,tmp4(:,:,:,1),tmp4(:,:,:,2),tmp4(:,:,:,3))
        !      endif
        !      this%intSharp_h(:,:,:,1) = this%intSharp_gam *
        !      (norm(:,:,:,1)*tmp4(:,:,:,1) + norm(:,:,:,2)*tmp4(:,:,:,2) +
        !      norm(:,:,:,3)*tmp4(:,:,:,3))
        !      this%intSharp_h(:,:,:,2) = zero
        !      this%intSharp_h(:,:,:,3) = zero

        !      do i = 1,this%ns
        !         !high order FD terms
        !         this%intSharp_hDiff(:,:,:,1) = this%intSharp_hDiff(:,:,:,1) +
        !         hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,1)
        !         this%intSharp_hDiff(:,:,:,2) = this%intSharp_hDiff(:,:,:,2) +
        !         hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,2)
        !         this%intSharp_hDiff(:,:,:,3) = this%intSharp_hDiff(:,:,:,3) +
        !         hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,3)
        !      enddo

        !   endif

        !   !fix momentum and energy above
        !   !add kinematic below

        else ! i.e. NOT:  if(this%intSharp_spf.and.useNewSPF.and.useNewSPFfull)  


           !enfoce surface normals sum to zero --- correct intSharp_a for
           !material ns --- this should be modifed for ns > 2

           !sharpening terms
           tmp4 = zero
           do j = 1,this%ns-1
              tmp4(:,:,:,1) = tmp4(:,:,:,1) +this%material(j)%intSharp_a(:,:,:,1)
              tmp4(:,:,:,2) = tmp4(:,:,:,2) +this%material(j)%intSharp_a(:,:,:,2)
              tmp4(:,:,:,3) = tmp4(:,:,:,3) +this%material(j)%intSharp_a(:,:,:,3)
           enddo
           this%material(this%ns)%intSharp_a(:,:,:,1) = -tmp4(:,:,:,1)
           this%material(this%ns)%intSharp_a(:,:,:,2) = -tmp4(:,:,:,2)
           this%material(this%ns)%intSharp_a(:,:,:,3) = -tmp4(:,:,:,3)

           !out-of-bounds diffusion terms
           tmp4 = zero
           do j = 1,this%ns-1
              tmp4(:,:,:,1) = tmp4(:,:,:,1) +this%material(j)%intSharp_aDiff(:,:,:,1)
              tmp4(:,:,:,2) = tmp4(:,:,:,2) +this%material(j)%intSharp_aDiff(:,:,:,2)
              tmp4(:,:,:,3) = tmp4(:,:,:,3) +this%material(j)%intSharp_aDiff(:,:,:,3)
           enddo
           this%material(this%ns)%intSharp_aDiff(:,:,:,1) = -tmp4(:,:,:,1)
           this%material(this%ns)%intSharp_aDiff(:,:,:,2) = -tmp4(:,:,:,2)
           this%material(this%ns)%intSharp_aDiff(:,:,:,3) = -tmp4(:,:,:,3)

           !FV not needed since antiDiffFVint is consistent



           !mass flux --  rho*Ys flux --  R_i
           do i = 1,this%ns
              if(useRhoLocal) then !use local component density --- recommended
                 if(this%intSharp_spf.and.useNewSPF) then
                    !already calculated
                 else
                    this%material(i)%intSharp_R(:,:,:,1) = rhoi(:,:,:,i) *this%material(i)%intSharp_a(:,:,:,1)
                    this%material(i)%intSharp_R(:,:,:,2) = rhoi(:,:,:,i) *this%material(i)%intSharp_a(:,:,:,2)
                    this%material(i)%intSharp_R(:,:,:,3) = rhoi(:,:,:,i) *this%material(i)%intSharp_a(:,:,:,3)
                 endif

                 this%material(i)%intSharp_RDiff(:,:,:,1) = rhoi(:,:,:,i) *this%material(i)%intSharp_aDiff(:,:,:,1)
                 this%material(i)%intSharp_RDiff(:,:,:,2) = rhoi(:,:,:,i) *this%material(i)%intSharp_aDiff(:,:,:,2)
                 this%material(i)%intSharp_RDiff(:,:,:,3) = rhoi(:,:,:,i) *this%material(i)%intSharp_aDiff(:,:,:,3)

              else
                 if(this%intSharp_spf.and.useNewSPF) then
                    !already calculated
                 else
                    this%material(i)%intSharp_R(:,:,:,1) =this%material(i)%elastic%rho0 * this%material(i)%intSharp_a(:,:,:,1)
                    this%material(i)%intSharp_R(:,:,:,2) =this%material(i)%elastic%rho0 * this%material(i)%intSharp_a(:,:,:,2)
                    this%material(i)%intSharp_R(:,:,:,3) =this%material(i)%elastic%rho0 * this%material(i)%intSharp_a(:,:,:,3)
                 endif

                 this%material(i)%intSharp_RDiff(:,:,:,1) =this%material(i)%elastic%rho0 * this%material(i)%intSharp_aDiff(:,:,:,1)
                 this%material(i)%intSharp_RDiff(:,:,:,2) =this%material(i)%elastic%rho0 * this%material(i)%intSharp_aDiff(:,:,:,2)
                 this%material(i)%intSharp_RDiff(:,:,:,3) =this%material(i)%elastic%rho0 * this%material(i)%intSharp_aDiff(:,:,:,3)

              endif

           enddo


           ! momentum flux -- f_i
           this%intSharp_f = zero
           this%intSharp_fDiff = zero
           this%intSharp_fFV = zero

           !momentum term is nonzero if coupling is turned on
           if(this%intSharp_cpl) then
              if(this%intSharp_spf) then 
                 spf_f = zero
                 !spf update term
                 do i = 1,this%ns
                    !this%intSharp_f(:,:,:,1) = this%intSharp_f(:,:,:,1) +
                    !this%material(i)%intSharp_R(:,:,:,1)
                    spf_f = spf_f + this%material(i)%intSharp_R(:,:,:,1)
                 enddo
                 this%intSharp_f(:,:,:,1) = spf_f * u 
                 this%intSharp_f(:,:,:,2) = spf_f * v
                 this%intSharp_f(:,:,:,3) = spf_f * w
                 !this%intSharp_f(:,:,:,2) = zero
                 !this%intSharp_f(:,:,:,3) = zero


                 do i = 1,this%ns
                    !high order FD terms
                    this%intSharp_fDiff(:,:,:,1) = this%intSharp_fDiff(:,:,:,1)+ this%material(i)%intSharp_RDiff(:,:,:,1)
                    this%intSharp_fDiff(:,:,:,2) = this%intSharp_fDiff(:,:,:,2)+ this%material(i)%intSharp_RDiff(:,:,:,2)  
                    this%intSharp_fDiff(:,:,:,3) = this%intSharp_fDiff(:,:,:,3)+ this%material(i)%intSharp_RDiff(:,:,:,3)
                 enddo

              else !i.e. if divergence method 
                 fv_f = zero
                 do i = 1,this%ns
                    !low order FD terms
                    this%intSharp_f(:,:,:,1) = this%intSharp_f(:,:,:,1) +this%material(i)%intSharp_R(:,:,:,1)
                    this%intSharp_f(:,:,:,2) = this%intSharp_f(:,:,:,2) +this%material(i)%intSharp_R(:,:,:,2)  
                    this%intSharp_f(:,:,:,3) = this%intSharp_f(:,:,:,3) +this%material(i)%intSharp_R(:,:,:,3)

                    !high order FD terms
                    this%intSharp_fDiff(:,:,:,1) = this%intSharp_fDiff(:,:,:,1)+ this%material(i)%intSharp_RDiff(:,:,:,1)
                    this%intSharp_fDiff(:,:,:,2) = this%intSharp_fDiff(:,:,:,2)+ this%material(i)%intSharp_RDiff(:,:,:,2)  
                    this%intSharp_fDiff(:,:,:,3) = this%intSharp_fDiff(:,:,:,3)+ this%material(i)%intSharp_RDiff(:,:,:,3)


                    !FV term
                    if(this%intSharp_ufv) then
                       fv_f = fv_f +rhoiFVint(:,:,:,:,i)*antiDiffFVint(:,:,:,:,i) 
                    endif
                 enddo

                 !FV terms
                 if(this%intSharp_ufv) then

                    call divergenceFV(this,fv_f*uFVint,this%intSharp_fFV(:,:,:,1),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                    call divergenceFV(this,fv_f*vFVint,this%intSharp_fFV(:,:,:,2),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                    call divergenceFV(this,fv_f*wFVint,this%intSharp_fFV(:,:,:,3),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                 endif

              endif
           endif


           !enthalpy flux -- (intSharp_h)_j = sum_i ( rho_i h_i (a_i)_j)
           this%intSharp_h = zero
           this%intSharp_hDiff = zero
           this%intSharp_hFV = zero
           this%intSharp_pFV = zero
           !enthalpy term is nonzero if coupling is turned on
           if(this%intSharp_cpl) then
              if(this%intSharp_spf) then
                 !spf update term 
                 do i = 1,this%ns
                    call this%material(i)%get_enthalpy(hi(:,:,:,i))
                    this%intSharp_h(:,:,:,1) = this%intSharp_h(:,:,:,1) + hi(:,:,:,i) * this%material(i)%intSharp_R(:,:,:,1)
                 enddo
                 this%intSharp_h(:,:,:,1) = this%intSharp_h(:,:,:,1) + spf_f * half*(u**two + v**two + w**two) !new form
                 this%intSharp_h(:,:,:,2) = zero
                 this%intSharp_h(:,:,:,3) = zero

                 do i = 1,this%ns
                    !high order FD terms
                    this%intSharp_hDiff(:,:,:,1) = this%intSharp_hDiff(:,:,:,1) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,1)
                    this%intSharp_hDiff(:,:,:,2) = this%intSharp_hDiff(:,:,:,2) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,2)
                    this%intSharp_hDiff(:,:,:,3) = this%intSharp_hDiff(:,:,:,3) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,3)
                 enddo

              else

                 do i = 1,this%ns
                    if(this%PTeqb) then
                         call this%material(i)%get_enthalpy(hi(:,:,:,i))
                    else
                          call this%material(i)%getSpeciesDensity(rho,rhom)

                          hi(:,:,:,i) = this%material(i)%hydro%gam*(this%material(i)%p + this%material(i)%hydro%PInf) * this%material(i)%hydro%onebygam_m1 

                    endif
                    !low order FD terms
                    this%intSharp_h(:,:,:,1) = this%intSharp_h(:,:,:,1) + hi(:,:,:,i) * this%material(i)%intSharp_R(:,:,:,1)
                    this%intSharp_h(:,:,:,2) = this%intSharp_h(:,:,:,2) + hi(:,:,:,i) * this%material(i)%intSharp_R(:,:,:,2)
                    this%intSharp_h(:,:,:,3) = this%intSharp_h(:,:,:,3) + hi(:,:,:,i) * this%material(i)%intSharp_R(:,:,:,3)

                    !high order FD terms
                    this%intSharp_hDiff(:,:,:,1) = this%intSharp_hDiff(:,:,:,1) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,1)
                    this%intSharp_hDiff(:,:,:,2) = this%intSharp_hDiff(:,:,:,2) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,2)
                    this%intSharp_hDiff(:,:,:,3) = this%intSharp_hDiff(:,:,:,3) + hi(:,:,:,i) * this%material(i)%intSharp_RDiff(:,:,:,3)

                    !FV term
                    if(this%intSharp_ufv) then                                
                     ! call interpolateFV(this,hi(:,:,:,i),hiFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                       call interpolateFV(this,this%material(i)%p,pFVint(:,:,:,:,i),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                       hiFVint(:,:,:,:,i) = (this%material(i)%hydro%gam*pFVint(:,:,:,:,i) + this%material(i)%hydro%gam*this%material(i)%hydro%PInf) *this%material(i)%hydro%onebygam_m1 
                    endif
                 enddo

                 !finish FV enthalpy term
                 if(this%intSharp_ufv) then
                    fv_h = zero
                    fv_k = zero
                    fv_p = zero

     
                    do i=1,this%ns
      
                       if(this%PTeqb) then
                           fv_h = fv_h + antiDiffFVint(:,:,:,:,i)*(rhoiFVint_local(:,:,:,:,i)*hiFVint(:,:,:,:,i))
                           fv_k = fv_k + antiDiffFVint(:,:,:,:,i)*rhoiFVint(:,:,:,:,i)*half*(uFVint**two+vFVint**two+wFVint**two)
                       else
                           fv_h = fv_h + antiDiffFVint(:,:,:,:,i)*(hiFVint(:,:,:,:,i)) 
                           fv_k = fv_k + antiDiffFVint(:,:,:,:,i)*rhoiFVint(:,:,:,:,i)*half*(uFVint**two+vFVint**two+wFVint**two)
                           fv_p = fv_p + antiDiffFVint(:,:,:,:,i)*(pFVint(:,:,:,:,i)) 
                       endif
                    enddo

                    call divergenceFV(this,fv_p,this%intSharp_pFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                    call divergenceFV(this,fv_h,this%intSharp_hFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                    call divergenceFV(this,fv_k,this%intSharp_kFV,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                 endif

              endif
           endif


           !kinematic fluxes
           do i=1,this%ns
              this%material(i)%intSharp_rg = zero
              this%material(i)%intSharp_rgt = zero
              this%material(i)%intSharp_rgp = zero
              this%material(i)%intSharp_rgDiff = zero
              this%material(i)%intSharp_rgtDiff = zero
              this%material(i)%intSharp_rgpDiff = zero
              this%material(i)%intSharp_rgFV = zero
              this%material(i)%intSharp_rgtFV = zero
              this%material(i)%intSharp_rgpFV = zero
              this%material(i)%intSharp_gFV = zero
              this%material(i)%intSharp_gtFV = zero
              this%material(i)%intSharp_gpFV = zero
           enddo


           if(this%intSharp_cpl.AND.this%intSharp_cpg) then
               if(this%intSharp_cpg_west) then !new implementation based on Jacob's derivation

                   ! Compute divergence of sharpening mass flux (low order)
                   call divergence(this%decomp,this%derD02,this%intSharp_f(:,:,:,1),this%intSharp_f(:,:,:,2),this%intSharp_f(:,:,:,3),tmp4(:,:,:,1),-x_bc,-y_bc,-z_bc) ! mass fraction equation is anti-symmetric (because a symmetry boundary implies no mass flux across it)
                   ! Compute divergence of sharpening mass flux (high order)
                   call divergence(this%decomp,this%der,this%intSharp_fDiff(:,:,:,1),this%intSharp_fDiff(:,:,:,2),this%intSharp_fDiff(:,:,:,3),tmp4(:,:,:,2),-x_bc,-y_bc,-z_bc) ! mass fraction equation is anti-symmetric (because a symmetry boundary implies no mass flux across it)
                   do i=1,this%ns
                      do j=1,9
                         if(this%intSharp_spf) then
                            !low order FD terms
                            this%material(i)%intSharp_rg (:,:,:,j,1) = this%material(i)%intSharp_rg (:,:,:,j,1) + third * spf_f * this%material(i)%g  (:,:,:,j)
                            this%material(i)%intSharp_rgt(:,:,:,j,1) = this%material(i)%intSharp_rgt(:,:,:,j,1) + third * spf_f * this%material(i)%g_t(:,:,:,j)
                         else
                            !low order FD terms
                            this%material(i)%intSharp_rg (:,:,:,j,1) = this%material(i)%intSharp_rg (:,:,:,j,1) + third * tmp4(:,:,:,1)   * this%material(i)%g  (:,:,:,j)
                            this%material(i)%intSharp_rgt(:,:,:,j,1) = this%material(i)%intSharp_rgt(:,:,:,j,1) + third * tmp4(:,:,:,1)   * this%material(i)%g_t(:,:,:,j)
                         endif

                         !high order FD terms
                         this%material(i)%intSharp_rgDiff (:,:,:,j,1) = this%material(i)%intSharp_rgDiff (:,:,:,j,1) + third * tmp4(:,:,:,2)   * this%material(i)%g  (:,:,:,j)
                         this%material(i)%intSharp_rgtDiff(:,:,:,j,1) = this%material(i)%intSharp_rgtDiff(:,:,:,j,1) + third * tmp4(:,:,:,2)   * this%material(i)%g_t(:,:,:,j)
                      enddo
                   enddo

                   !FV terms
                   if(this%intSharp_ufv) then
                      call divergenceFV(this,fv_f,tmp,dx,dy,dz,periodicx,periodicy,periodicz,-this%x_bc,-this%y_bc,-this%z_bc)
                      do i=1,this%ns
                         do j=1,9
                            this%material(i)%intSharp_rgFV (:,:,:,j) = this%material(i)%intSharp_rgFV (:,:,:,j) + third * tmp * this%material(i)%g  (:,:,:,j)
                            this%material(i)%intSharp_rgtFV(:,:,:,j) = this%material(i)%intSharp_rgtFV(:,:,:,j) + third * tmp * this%material(i)%g_t(:,:,:,j)
                         enddo
                      enddo
                   endif

               else !implementation used in CTR brief

                   do i=1,this%ns
                      do j=1,9
                         do k=1,3
                            if(this%intSharp_spf) then
                               !low order FD terms
                               this%material(i)%intSharp_rg (:,:,:,j,k) = this%material(i)%intSharp_rg (:,:,:,j,k) + spf_f*this%material(i)%g  (:,:,:,j)
                               this%material(i)%intSharp_rgt(:,:,:,j,k) = this%material(i)%intSharp_rgt(:,:,:,j,k) + spf_f*this%material(i)%g_t(:,:,:,j)
                               this%material(i)%intSharp_rgp(:,:,:,j,k) = this%material(i)%intSharp_rgp(:,:,:,j,k) + spf_f*this%material(i)%g_p(:,:,:,j)
                            else
                               !low order FD terms
                               this%material(i)%intSharp_rg (:,:,:,j,k) = this%material(i)%intSharp_rg (:,:,:,j,k) + this%intSharp_f(:,:,:,k)*this%material(i)%g  (:,:,:,j)
                               this%material(i)%intSharp_rgt(:,:,:,j,k) = this%material(i)%intSharp_rgt(:,:,:,j,k) + this%intSharp_f(:,:,:,k)*this%material(i)%g_t(:,:,:,j)
                               this%material(i)%intSharp_rgp(:,:,:,j,k) = this%material(i)%intSharp_rgp(:,:,:,j,k) + this%intSharp_f(:,:,:,k)*this%material(i)%g_p(:,:,:,j)
                            endif

                            !high order FD terms
                            this%material(i)%intSharp_rgDiff (:,:,:,j,k) = this%material(i)%intSharp_rgDiff (:,:,:,j,k) + this%intSharp_fDiff(:,:,:,k)*this%material(i)%g  (:,:,:,j)
                            this%material(i)%intSharp_rgtDiff(:,:,:,j,k) = this%material(i)%intSharp_rgtDiff(:,:,:,j,k) + this%intSharp_fDiff(:,:,:,k)*this%material(i)%g_t(:,:,:,j)
                            this%material(i)%intSharp_rgpDiff(:,:,:,j,k) = this%material(i)%intSharp_rgpDiff(:,:,:,j,k) + this%intSharp_fDiff(:,:,:,k)*this%material(i)%g_p(:,:,:,j)

                         enddo
                      enddo
                   enddo


                   if(this%intSharp_spf) then 
                      do i=1,this%ns
                         do j=1,9
                            do k=2,3 !zero components 2 and 3 for gradient form
                               this%material(i)%intSharp_rg (:,:,:,j,k) = zero
                               this%material(i)%intSharp_rgt(:,:,:,j,k) = zero
                               this%material(i)%intSharp_rgp(:,:,:,j,k) = zero
                            enddo
                         enddo
                      enddo
                   endif


                   !FV terms
                   if(this%intSharp_ufv) then
                      do i=1,this%ns
                         do j=1,9
                            !interpolate in loop to reduce storage
                            call interpolateFV(this,this%material(i)%g  (:,:,:,j), gFVint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                            call interpolateFV(this,this%material(i)%g_t(:,:,:,j),gtFVint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                            call interpolateFV(this,this%material(i)%g_p(:,:,:,j),gpFVint,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

                            !conservative form
                            call divergenceFV(this,fv_f* gFVint,this%material(i)%intSharp_rgFV (:,:,:,j),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                            call divergenceFV(this,fv_f*gtFVint,this%material(i)%intSharp_rgtFV(:,:,:,j),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                            call divergenceFV(this,fv_f*gpFVint,this%material(i)%intSharp_rgpFV(:,:,:,j),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

                            !non-conservative form
                            call divergenceFV(this,fv_f/rhoFVint* gFVint,this%material(i)%intSharp_gFV (:,:,:,j),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                            call divergenceFV(this,fv_f/rhoFVint*gtFVint,this%material(i)%intSharp_gtFV(:,:,:,j),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                            call divergenceFV(this,fv_f/rhoFVint*gpFVint,this%material(i)%intSharp_gpFV(:,:,:,j),dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
                         enddo
                      enddo
                   endif
               endif


           endif

        endif

        !filter remaining conservative fluxes and updates -- need to check bc symmetry flags
        if(this%intSharp_flt) then 
           do i=1,this%ns
              call filter3D(this%decomp, this%fil, this%material(i)%intSharp_R(:,:,:,1), iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil, this%material(i)%intSharp_R(:,:,:,2), iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil, this%material(i)%intSharp_R(:,:,:,3), iflag, x_bc, y_bc,z_bc)

              call filter3D(this%decomp, this%fil, this%material(i)%intSharp_RDiff(:,:,:,1), iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil, this%material(i)%intSharp_RDiff(:,:,:,2), iflag, x_bc, y_bc,z_bc)
              call filter3D(this%decomp, this%fil, this%material(i)%intSharp_RDiff(:,:,:,3), iflag, x_bc, y_bc,z_bc)
           enddo

           call filter3D(this%decomp, this%fil, this%intSharp_f(:,:,:,1), iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_f(:,:,:,2), iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_f(:,:,:,3), iflag, x_bc, y_bc,z_bc)

           call filter3D(this%decomp, this%fil, this%intSharp_fDiff(:,:,:,1), iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_fDiff(:,:,:,2), iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_fDiff(:,:,:,3), iflag, x_bc, y_bc,z_bc)

           call filter3D(this%decomp, this%fil, this%intSharp_h(:,:,:,1), iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_h(:,:,:,2), iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_h(:,:,:,3), iflag, x_bc, y_bc,z_bc)

           call filter3D(this%decomp, this%fil, this%intSharp_hDiff(:,:,:,1), iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_hDiff(:,:,:,2), iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_hDiff(:,:,:,3), iflag, x_bc, y_bc,z_bc)

           call filter3D(this%decomp, this%fil, this%intSharp_fFV, iflag, x_bc, y_bc,z_bc)
           call filter3D(this%decomp, this%fil, this%intSharp_hFV, iflag, x_bc, y_bc,z_bc)

           do i=1,this%ns
              do j=1,9
                 do k=1,3
                    call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rg(:,:,:,j,k), iflag, x_bc, y_bc,z_bc)
                    call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rgt(:,:,:,j,k), iflag, x_bc, y_bc,z_bc)
                    call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rgp(:,:,:,j,k), iflag, x_bc, y_bc,z_bc)

                    call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rgDiff(:,:,:,j,k), iflag, x_bc, y_bc,z_bc)
                    call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rgtDiff(:,:,:,j,k), iflag, x_bc, y_bc,z_bc)
                    call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rgpDiff(:,:,:,j,k), iflag, x_bc, y_bc,z_bc)
                 enddo

                 call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rgFV(:,:,:,j), iflag, x_bc, y_bc,z_bc)
                 call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rgtFV(:,:,:,j), iflag, x_bc, y_bc,z_bc)
                 call filter3D(this%decomp, this%fil, this%material(i)%intSharp_rgpFV(:,:,:,j), iflag, x_bc, y_bc,z_bc)

              enddo
           enddo

        endif

        

 end subroutine get_intSharp

 subroutine get_gradp(this,rho,x_bc,y_bc,z_bc,dx,dy,dz,periodicx,periodicy,periodicz,u,v,w)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z,transpose_z_to_y
        use operators, only: divergence,gradient,filter3D, laplacian,interpolateFV, divergenceFV
        use constants,       only: zero,epssmall,eps,one,two,third,half
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture),                               intent(inout) :: this
        integer, dimension(2),                              intent(in) :: x_bc,y_bc, z_bc
        real(rkind),                                        intent(in) :: dx,dy,dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w
        logical,                                            intent(in) :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)  :: lapVF,VFmag,tanhmask, GVFmag, GPHImag, mask2, updatedKappa, weight, kappaSum,phi, xi, mu,d2vfdx2,d2vfdy2,d2vfdz2
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradVF, gradphi,gradxi, gradVFk, p_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,3) :: NMint,gradVF_FV,gradVFint
        integer :: iflag = one
        real(rkind) :: r = 0.5D0, nmask = 40, minVF = 1D-6, tmask = 0.2d0, e =1D-100
        !TODO: add additional arrays to be used locally in calculation of
        !surface tension force
        integer :: i,j,k,n
        integer, dimension(2) :: coords

        call interpolateFV_6(this,this%material(1)%p,p_int,periodicx,periodicy,periodicz,-this%x_bc,this%y_bc,this%z_bc)
        call divergenceFV_6(this,p_int,this%gradp,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
        !call gradient(this%decomp,this%der,this%material(1)%p,this%gradp(:,:,:,1),this%gradp(:,:,:,2),this%gradp(:,:,:,3))

       
    end subroutine

  subroutine get_surfaceTensionPE(this,rho,x_bc,y_bc,z_bc,dx,dy,dz,periodicx,periodicy,periodicz,u,v,w)
        use decomp_2d, only: transpose_y_to_x,transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        use operators, only: divergence,gradient,filter3D,laplacian,interpolateFV_x, interpolateFV_y, interpolateFV_z
        use constants,       only: zero,epssmall,eps,one,two,third,half
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture),                               intent(inout) :: this
        integer, dimension(2),                              intent(in) :: x_bc,y_bc, z_bc
        real(rkind),                                        intent(in) :: dx,dy,dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w
        logical,                                            intent(in) :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)  :: lapVF,VFmag,tanhmask,GVFmag, GYSmag, GPHImag, mask2, updatedKappa, weight,kappaSum,phi, xi, mu,d2vfdx2,d2vfdy2,d2vfdz2
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradYs,gradphi,gradxi, gradVFk, gradYs_x, gradYs_y, gradYs_z,GVFmag_int,gradVF
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,3) :: NMint,gradVF_FV,gradVF_int
        integer :: iflag = one
        real(rkind) :: r = 0.4D0, nmask = 40, minVF = 1D-6, tmask = 0.2d0, e = 1D-100
        !TODO: add additional arrays to be used locally in calculation of
        !surface tension force
        integer :: i,j,k,n
        integer, dimension(2) :: coords


        if  (this%ns.ne.2) then
            call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
        endif

        if(.NOT. this%use_Stagg) then
            call gradient(this%decomp,this%der,this%material(1)%VF,gradYs(:,:,:,1),gradYs(:,:,:,2),gradYs(:,:,:,3))
            GYSmag = sqrt( gradYs(:,:,:,1)**two +gradYs(:,:,:,2)**two+gradYs(:,:,:,3)**two ) 
            this%surfaceTension_pe = this%surfaceTension_coeff*GYSmag
        else

           call gradient(this%decomp,this%derD06,this%material(1)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))            
           call interpolateFV_x(this%decomp,this%interpMid,gradVF(:,:,:,1),gradVF_int(:,:,:,1,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,gradVF(:,:,:,2),gradVF_int(:,:,:,2,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,gradVF(:,:,:,3),gradVF_int(:,:,:,3,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,gradVF(:,:,:,1),gradVF_int(:,:,:,1,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,gradVF(:,:,:,2),gradVF_int(:,:,:,2,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,gradVF(:,:,:,3),gradVF_int(:,:,:,3,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,gradVF(:,:,:,1),gradVF_int(:,:,:,1,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,gradVF(:,:,:,2),gradVF_int(:,:,:,2,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,gradVF(:,:,:,3),gradVF_int(:,:,:,3,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

           this%surfaceTension_pe_x = this%surfaceTension_coeff*sqrt( gradVF_int(:,:,:,1,1)**2 + gradVF_int(:,:,:,2,1)**2 + gradVF_int(:,:,:,3,1)**2 )
           this%surfaceTension_pe_y = this%surfaceTension_coeff*sqrt( gradVF_int(:,:,:,1,2)**2 + gradVF_int(:,:,:,2,2)**2 + gradVF_int(:,:,:,3,2)**2 )
           this%surfaceTension_pe_z = this%surfaceTension_coeff*sqrt( gradVF_int(:,:,:,1,3)**2 + gradVF_int(:,:,:,2,3)**2 + gradVF_int(:,:,:,3,3)**2 )

        end if 
        !high order derivative

    end subroutine

   subroutine get_surfaceTensionCnsrv(this,rho,x_bc,y_bc,z_bc,dx,dy,dz,periodicx,periodicy,periodicz,u,v,w)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        use operators, only: divergence,gradient,filter3D,laplacian,interpolateFV_x, interpolateFV_y, interpolateFV_z
        use constants,       only: zero,epssmall,eps,one,two,third,half
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture),                               intent(inout) :: this
        integer, dimension(2),                              intent(in) :: x_bc,y_bc, z_bc
        real(rkind),                                        intent(in) :: dx,dy,dz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w
        logical,                                            intent(in) :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)  :: lapVF,VFmag,tanhmask, GYSmag,GVFmag, GPHImag,Gximag, mask2, updatedKappa, weight, kappaSum,phi, xi, mu,d2vfdx2,d2vfdy2,d2vfdz2,s1,s2,s3
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)  :: GVFmag_x,GVFmag_y, GVFmag_z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradYs, gradphi,gradxi, gradVFk, gradYs_x, gradYs_y, gradYs_z, gradVF, GVFmag_int
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,3) :: NMint,gradVF_FV,gradVF_int
        integer :: iflag = one
        real(rkind) :: r = 0.5D0, nmask = 40, minVF = 1D-6, tmask = 0.2d0, e = 1D-100
        !TODO: add additional arrays to be used locally in calculation of
        !surface tension force
        integer :: i,j,k,n
        integer, dimension(2) :: coords


        if  (this%ns.ne.2) then
            call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
        endif
        

        if( .NOT. this%use_Stagg) then

            call gradient(this%decomp,this%der,this%material(1)%VF,gradYs(:,:,:,1),gradYs(:,:,:,2),gradYs(:,:,:,3))

            !high order derivative
            this%gradVF = gradYs

        
            call gradient(this%decomp,this%der,this%material(1)%p,this%gradp(:,:,:,1),this%gradp(:,:,:,2),this%gradp(:,:,:,3))
            GYSmag = sqrt( gradYs(:,:,:,1)**two + gradYs(:,:,:,2)**two+gradYs(:,:,:,3)**two )

            where( GYSmag .LE. eps)
               this%surfaceTension_fxx = this%surfaceTension_coeff*( GYSmag  )
               this%surfaceTension_fyy = this%surfaceTension_coeff*( GYSmag  )
               this%surfaceTension_fzz = this%surfaceTension_coeff*( GYSmag  )
               this%surfaceTension_fxy = this%surfaceTension_coeff*( zero )
               this%surfaceTension_fxz = this%surfaceTension_coeff*( zero )
               this%surfaceTension_fyz = this%surfaceTension_coeff*( zero )
           elsewhere
               this%surfaceTension_fxx = this%surfaceTension_coeff*( GYSmag - gradYs(:,:,:,1)*gradYs(:,:,:,1)/GYSmag )
               this%surfaceTension_fyy = this%surfaceTension_coeff*( GYSmag - gradYs(:,:,:,2)*gradYs(:,:,:,2)/GYSmag )
               this%surfaceTension_fzz = this%surfaceTension_coeff*( GYSmag - gradYs(:,:,:,3)*gradYs(:,:,:,3)/GYSmag )
               this%surfaceTension_fxy = this%surfaceTension_coeff*(  -gradYs(:,:,:,2)*gradYs(:,:,:,1)/GYSmag )
               this%surfaceTension_fxz = this%surfaceTension_coeff*(  -gradYs(:,:,:,1)*gradYs(:,:,:,3)/GYSmag )
               this%surfaceTension_fyz = this%surfaceTension_coeff*(  -gradYs(:,:,:,2)*gradYs(:,:,:,3)/GYSmag )
           endwhere

           call divergence(this%decomp,this%der,this%surfaceTension_fxx,this%surfaceTension_fxy,this%surfaceTension_fxz,this%surfaceTension_f(:,:,:,1),x_bc,y_bc,z_bc)
           call divergence(this%decomp,this%der,this%surfaceTension_fxy,this%surfaceTension_fyy,this%surfaceTension_fyz,this%surfaceTension_f(:,:,:,2),x_bc,y_bc,z_bc)
           call divergence(this%decomp,this%der,this%surfaceTension_fxz,this%surfaceTension_fyz,this%surfaceTension_fzz,this%surfaceTension_f(:,:,:,3),x_bc,y_bc,z_bc)

           s1 = u*this%surfaceTension_fxx + v*this%surfaceTension_fxy + w*this%surfaceTension_fxz
           s2 = u*this%surfaceTension_fxy + v*this%surfaceTension_fyy +w*this%surfaceTension_fyz
           s3 = u*this%surfaceTension_fxz + v*this%surfaceTension_fyz +w*this%surfaceTension_fzz
      
           call divergence(this%decomp,this%der,s1, s2, s3,this%surfaceTension_e,x_bc,y_bc,z_bc)

        else

           call gradient(this%decomp,this%derD06,this%material(1)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))
           call interpolateFV_x(this%decomp,this%interpMid,gradVF(:,:,:,1),gradVF_int(:,:,:,1,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,gradVF(:,:,:,2),gradVF_int(:,:,:,2,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,gradVF(:,:,:,3),gradVF_int(:,:,:,3,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,gradVF(:,:,:,1),gradVF_int(:,:,:,1,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,gradVF(:,:,:,2),gradVF_int(:,:,:,2,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,gradVF(:,:,:,3),gradVF_int(:,:,:,3,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,gradVF(:,:,:,1),gradVF_int(:,:,:,1,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,gradVF(:,:,:,2),gradVF_int(:,:,:,2,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,gradVF(:,:,:,3),gradVF_int(:,:,:,3,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

           GVFmag_x = sqrt( gradVF_int(:,:,:,1,1)**2 + gradVF_int(:,:,:,2,1)**2 + gradVF_int(:,:,:,3,1)**2 )
           GVFmag_y = sqrt( gradVF_int(:,:,:,1,2)**2 + gradVF_int(:,:,:,2,2)**2 + gradVF_int(:,:,:,3,2)**2 )
           GVFmag_z = sqrt( gradVF_int(:,:,:,1,3)**2 + gradVF_int(:,:,:,2,3)**2 + gradVF_int(:,:,:,3,3)**2 )

           where( GVFmag_x .LE. eps)
               this%surfaceTension_fxx   = this%surfaceTension_coeff*(GVFmag_x  )
               this%surfaceTension_fxy_x = this%surfaceTension_coeff*( zero )
               this%surfaceTension_fxz_x = this%surfaceTension_coeff*( zero )
           elsewhere
               this%surfaceTension_fxx   = this%surfaceTension_coeff*(GVFmag_x - gradVF_int(:,:,:,1,1)*gradVF_int(:,:,:,1,1))
               this%surfaceTension_fxy_x = this%surfaceTension_coeff*( - gradVF_int(:,:,:,1,1)*gradVF_int(:,:,:,2,1))
               this%surfaceTension_fxz_x = this%surfaceTension_coeff*( - gradVF_int(:,:,:,1,1)*gradVF_int(:,:,:,3,1))
           endwhere

           where( GVFmag_y .LE. eps)
               this%surfaceTension_fyy   = this%surfaceTension_coeff*(GVFmag_y  )
               this%surfaceTension_fxy_y = this%surfaceTension_coeff*( zero )
               this%surfaceTension_fyz_y = this%surfaceTension_coeff*( zero )
           elsewhere
               this%surfaceTension_fyy   = this%surfaceTension_coeff*(GVFmag_y - gradVF_int(:,:,:,2,2)*gradVF_int(:,:,:,2,2))
               this%surfaceTension_fxy_y = this%surfaceTension_coeff*( - gradVF_int(:,:,:,1,2)*gradVF_int(:,:,:,2,2))
               this%surfaceTension_fyz_y = this%surfaceTension_coeff*( - gradVF_int(:,:,:,2,2)*gradVF_int(:,:,:,3,2))
           
           endwhere

           where( GVFmag_z .LE. eps)
               this%surfaceTension_fzz   = this%surfaceTension_coeff*(GVFmag_z  )
               this%surfaceTension_fxz_z = this%surfaceTension_coeff*( zero )
               this%surfaceTension_fyz_z = this%surfaceTension_coeff*( zero )
           elsewhere
               this%surfaceTension_fzz   = this%surfaceTension_coeff*(GVFmag_z - gradVF_int(:,:,:,3,3)*gradVF_int(:,:,:,3,3))
               this%surfaceTension_fxz_z = this%surfaceTension_coeff*( - gradVF_int(:,:,:,1,3)*gradVF_int(:,:,:,3,3))
               this%surfaceTension_fyz_z = this%surfaceTension_coeff*( - gradVF_int(:,:,:,2,3)*gradVF_int(:,:,:,3,3))
           endwhere




!           call
!           interpolateFV_y(this%decomp,this%interpMid,GVFmag,GVFmag_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
!           call
!           interpolateFV_z(this%decomp,this%interpMid,GVFmag,GVFmag_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

            

!           where( this%material(1)%VF .GE. 1-this%intSharp_cut)
!                 this%xi(:,:,:,1) = this%intSharp_eps*log( (1-2*this%intSharp_cut+ e )/ (e))*(1/(1-2*this%intSharp_cut))
!           elsewhere( this%material(1)%VF .LE. this%intSharp_cut )
!                 this%xi(:,:,:,1)  = this%intSharp_eps*log( ( e) / (1-2*this%intSharp_cut + e))*(1/(1-2*this%intSharp_cut))
!           elsewhere
!                 this%xi(:,:,:,1)  = this%intSharp_eps*(1/(1-2*this%intSharp_cut))*log( (this%material(1)%VF-this%intSharp_cut + e )/ (1 - this%intSharp_cut -this%material(1)%VF + e) )
!           endwhere

!           call gradient(this%decomp,this%derD06,this%xi(:,:,:,1),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3))

!           call filter3D(this%decomp, this%gfil,gradxi(:,:,:,1),iflag,x_bc,y_bc,z_bc)
!           call filter3D(this%decomp, this%gfil,gradxi(:,:,:,2),iflag,x_bc,y_bc,z_bc)
!           call filter3D(this%decomp, this%gfil,gradxi(:,:,:,3),iflag,x_bc,y_bc,z_bc)

!            this%gradxi = gradxi

           !magnitude of surface vector
!           Gximag = sqrt( gradxi(:,:,:,1)**two + gradxi(:,:,:,2)**two+gradxi(:,:,:,3)**two )
!           GVFmag = sqrt( gradVF(:,:,:,1)**two +gradVF(:,:,:,2)**two+gradVF(:,:,:,3)**two )

           !surface normal
!           where (Gximag < eps)
!              this%norm(:,:,:,1) = zero
!              this%norm(:,:,:,2) = zero
!              this%norm(:,:,:,3) = zero
!           elsewhere
!              this%norm(:,:,:,1) = gradxi(:,:,:,1) / Gximag
!              this%norm(:,:,:,2) = gradxi(:,:,:,2) / Gximag
!              this%norm(:,:,:,3) = gradxi(:,:,:,3) / Gximag
!           endwhere

!           call interpolateFV_6(this,this%norm(:,:,:,1),NMint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc,this%y_bc, this%z_bc)
!           call interpolateFV_6(this,this%norm(:,:,:,2),NMint(:,:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
!           call interpolateFV_6(this,this%norm(:,:,:,3),NMint(:,:,:,:,3),periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)
!           this%normFV(:,:,:,1) = NMint(:,:,:,1,1)        
!           this%normFV(:,:,:,2) = NMint(:,:,:,2,2)
!           this%normFV(:,:,:,3) = NMint(:,:,:,3,3)

!           call interpolateFV_x(this%decomp,this%interpMid,GVFmag,GVFmag_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
!           call interpolateFV_y(this%decomp,this%interpMid,GVFmag,GVFmag_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
!           call interpolateFV_z(this%decomp,this%interpMid,GVFmag,GVFmag_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
   
!           where( GVFmag_int(:,:,:,1) .LE. eps)
!               this%surfaceTension_fxx   = this%surfaceTension_coeff*( GVFmag_int(:,:,:,1)  )
!               this%surfaceTension_fxy_x = this%surfaceTension_coeff*( zero )
!               this%surfaceTension_fxz_x = this%surfaceTension_coeff*( zero )
!           elsewhere( GVFmag_int(:,:,:,2) .LE. eps)
!               this%surfaceTension_fyy   = this%surfaceTension_coeff*( GVFmag_int(:,:,:,2)  )
!               this%surfaceTension_fxy_y = this%surfaceTension_coeff*( zero )
!               this%surfaceTension_fyz_y = this%surfaceTension_coeff*( zero )
!           elsewhere( GVFmag_int(:,:,:,3) .LE. eps)
!               this%surfaceTension_fzz   = this%surfaceTension_coeff*( GVFmag_int(:,:,:,3)  )
!               this%surfaceTension_fxz_z = this%surfaceTension_coeff*( zero )
!               this%surfaceTension_fyz_z = this%surfaceTension_coeff*( zero )

!           elsewhere
!               this%surfaceTension_fxx = this%surfaceTension_coeff*( GVFmag_int(:,:,:,1)*( 1 - NMint(:,:,:,1,1)*NMint(:,:,:,1,1)) )
!               this%surfaceTension_fyy = this%surfaceTension_coeff*( GVFmag_int(:,:,:,2)*( 1 - NMint(:,:,:,2,2)*NMint(:,:,:,2,2)) )
!               this%surfaceTension_fzz = this%surfaceTension_coeff*( GVFmag_int(:,:,:,3)*( 1 - NMint(:,:,:,3,3)*NMint(:,:,:,3,3)) )
!               this%surfaceTension_fxy_x = this%surfaceTension_coeff*( -NMint(:,:,:,2,1)*Nhis%normFV(:,:,:,1)*GVFmag_int(:,:,:,1) )
!               this%surfaceTension_fxz_x = this%surfaceTension_coeff*( -this%normFV(:,:,:,1)*this%normFV(:,:,:,3)*GVFmag_int(:,:,:,1) )
!               this%surfaceTension_fyz_y = this%surfaceTension_coeff*( -this%normFV(:,:,:,2)*this%normFV(:,:,:,3)*GVFmag_int(:,:,:,2) )
!               this%surfaceTension_fxy_y = this%surfaceTension_coeff*( -this%normFV(:,:,:,2)*this%normFV(:,:,:,1)*GVFmag_int(:,:,:,2) )
!               this%surfaceTension_fxz_z = this%surfaceTension_coeff*( -this%normFV(:,:,:,1)*this%normFV(:,:,:,3)*GVFmag_int(:,:,:,3) )
!               this%surfaceTension_fyz_z = this%surfaceTension_coeff*( -this%normFV(:,:,:,2)*this%normFV(:,:,:,3)*GVFmag_int(:,:,:,3) )

!           endwhere
 


        endif

    end subroutine 

    subroutine CheckP(this,rho,p,e)
       use operators, only: divergence,gradient,filter3D
       class(solid_mixture), intent(inout) :: this
       real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,p,e
       real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: num, dum, gamfac
       integer :: i

       dum = zero
       num = rho*e
       do i = 1,this%ns
           
          gamfac = this%material(i)%VF*this%material(i)%hydro%gam*this%material(i)%hydro%Pinf*this%material(i)%hydro%onebygam_m1
          num    = num - gamfac
          dum = dum + this%material(i)%VF*this%material(i)%hydro%onebygam_m1
          
       enddo 

       this%Pmix = num/dum

    end subroutine 

    subroutine getYs(this,rho)
        use operators, only: divergence,gradient,filter3D
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(inout) :: rho
        real(rkind) :: minYs1,minVF = 1d-6, minYs2
        integer :: imat

        minVF = 1d-6

        if(this%PTeqb) then
        do imat = 1, this%ns
          call this%material(imat)%getYs(rho)
        end do
        endif

   end subroutine 

    subroutine get_surfaceTension(this,rho,x_bc,y_bc,z_bc,dx,dy,dz,periodicx,periodicy,periodicz,u,v,w)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        use operators, only:divergence,gradient,filter3D,laplacian,interpolateFV_x,interpolateFV_y,interpolateFV_z,gradFV_x, gradFV_y, gradFV_z,gradFV_N2Fx, gradFV_N2Fy, gradFV_N2Fz
        use constants,       only: zero,epssmall,eps,one,two,third,half, pi
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture),                               intent(inout) :: this
        integer, dimension(2),                              intent(in) :: x_bc, y_bc, z_bc
        real(rkind),                                        intent(in) :: dx,dy,dz
	real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w
        logical,                                            intent(in) :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)  :: lapVF,udiv, divuphi,VFmag,tanhmask, GVFmag, GPHImag, mask2, updatedKappa, weight, kappaSum, phi, xi, mu,d2vfdx2,d2vfdy2,d2vfdz2,divu,divphiu,dirac,H
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradVF, gradphi, gradxi, gradVFk, p_int, VF_int, gradH,u_int, uphi_int, gradFV
	real(rkind), dimension(this%nxp,this%nyp,this%nzp,3,3) :: NMint,gradVF_FV,gradVFint
	integer :: iflag = one
	real(rkind) :: r = 0.4D0, nmask = 40, minVF = 1D-6, tmask = 0.2d0, e = 1D-100 
	!TODO: add additional arrays to be used locally in calculation of surface tension force
        integer :: i,j,k,n
        integer, dimension(2) :: coords

        if (this%ns.ne.2) then
            call GracefulExit("Surface tension is not defined for single-species, and not implemented for more than 2 species",4634)
        endif

        !initialize surface tension force and energy source to zero everywhere
        this%surfaceTension_f = 0.0d0
        this%surfaceTension_e = 0.0d0
        !this%kappa            = 0.0d0
        updatedKappa           = 0.0d0
        weight                 = 0.0d0
        kappaSum               = 0.0d0
        gradVFk                = 0.0d0
        gradVF                 = 0.0d0
        this%norm              = 0.0d0
        this%normFV            = 0.0d0
        this%xi(:,:,:,1)       = 0.0d0
        this%fmask             = 0.0d0
        GVFmag                 = 0.0d0
        gradxi                 = 0.0d0
        NMint                  = 0.0d0
        this%gradp             = 0.0d0
        !TODO: Compute curvature (use gradient and divergence operators)
        !for example, to take the gradient of the volume fraction of species # 1

       if( this%use_Stagg) then
           call interpolateFV_x(this%decomp,this%interpMid,this%material(1)%VF,VF_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,this%material(1)%VF,VF_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,this%material(1)%VF,VF_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call gradFV_x(this%decomp,this%derStagg,VF_int(:,:,:,1),gradVF(:,:,:,1),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_y(this%decomp,this%derStagg,VF_int(:,:,:,2),gradVF(:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
           call interpolateFV_x(this%decomp,this%interpMid,this%material(1)%p,p_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_y(this%decomp,this%interpMid,this%material(1)%p,p_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call interpolateFV_z(this%decomp,this%interpMid,this%material(1)%p,p_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
           call gradFV_x(this%decomp,this%derStagg,p_int(:,:,:,1),this%gradp(:,:,:,1),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
           call gradFV_y(this%decomp,this%derStagg,p_int(:,:,:,2),this%gradp(:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

!           call
!           gradFV_z(this%decomp,this%derStagg,VF_int(:,:,:,3),gradVF(:,:,:,3),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%
           call gradFV_z(this%decomp,this%derStagg,VF_int(:,:,:,3),gradVF(:,:,:,3),periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

       else
           call gradient(this%decomp,this%der,this%material(1)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3)) !high order derivative
       endif
       this%gradVF = gradVF
   !   call gradient(this%decomp,this%der,this%material(1)%p,this%gradp(:,:,:,1),this%gradp(:,:,:,2),this%gradp(:,:,:,3))
 
	if (this%use_gradVF) then

                if(this%use_FV) then

                   call gradient(this%decomp,this%derD02,this%material(1)%VF,gradVFk(:,:,:,1),gradVFk(:,:,:,2),gradVFk(:,:,:,3)) !high order derivative
                  ! call gradientFV(this,this%material(1)%VF,gradVF_FV,dx,dy,dz,periodicx,periodicy,periodicz, this%x_bc, this%y_bc, this%z_bc)

                   this%gradVF_FV = gradVF_FV
                else if (this%use_D04) then
                   call gradient(this%decomp,this%derD06,this%material(1)%VF,gradVFk(:,:,:,1),gradVFk(:,:,:,2),gradVFk(:,:,:,3)) !high order derivative
                else
                   call gradient(this%decomp,this%der,this%material(1)%VF,gradVFk(:,:,:,1),gradVFk(:,:,:,2),gradVFk(:,:,:,3)) !high order derivative
                endif



               ! if(this%use_FV) then
                  !magnitude of surface vector
               !   GVFmag = sqrt( gradVF_FV(:,:,:,1,1)**two + gradVF_FV(:,:,:,2,2)**two +gradVF_FV(:,:,:,3,3)**two )
 
               !    where (GVFmag < eps)
               !     this%norm(:,:,:,1) = zero
               !     this%norm(:,:,:,2) = zero
               !     this%norm(:,:,:,3) = zero
               !   elsewhere
               !     this%norm(:,:,:,1) = gradVF_FV(:,:,:,1,1) / GVFmag
               !     this%norm(:,:,:,2) = gradVF_FV(:,:,:,2,2) / GVFmag
               !     this%norm(:,:,:,3) = gradVF_FV(:,:,:,3,3) / GVFmag
               !   endwhere

               ! else 
                 !magnitude of surface vector
                 GVFmag = sqrt( gradVFk(:,:,:,1)**two + gradVFk(:,:,:,2)**two +gradVFk(:,:,:,3)**two )

            	!surface normal
           	  where (GVFmag < eps)
                    this%norm(:,:,:,1) = zero
                    this%norm(:,:,:,2) = zero
                    this%norm(:,:,:,3) = zero
              	  elsewhere
                    this%norm(:,:,:,1) = gradVFk(:,:,:,1) / GVFmag
                    this%norm(:,:,:,2) = gradVFk(:,:,:,2) / GVFmag
                    this%norm(:,:,:,3) = gradVFk(:,:,:,3) / GVFmag
                  endwhere
              !  endif


                this%fmask =(this%material(1)%VF-this%intSharp_cut)*(one-this%intSharp_cut-this%material(1)%VF)
	
		if (this%use_FV) then
             
	         call interpolateFV(this,this%norm(:,:,:,1),NMint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc, this%y_bc, this%z_bc)
		 call interpolateFV(this,this%norm(:,:,:,2),NMint(:,:,:,:,2),periodicx,periodicy,periodicz, this%x_bc,-this%y_bc, this%z_bc)
              	 call interpolateFV(this,this%norm(:,:,:,3),NMint(:,:,:,:,3),periodicx,periodicy,periodicz, this%x_bc, this%y_bc,-this%z_bc)
                 this%normFV(:,:,:,1) = NMint(:,:,:,1,1)	
                 this%normFV(:,:,:,2) = NMint(:,:,:,2,2)
                 this%normFV(:,:,:,3) = NMint(:,:,:,3,3)

                 call divergenceFV(this,this%normFV, this%kappa,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

                !  call divergence(this%decomp,this%derD02,this%norm(:,:,:,1),this%norm(:,:,:,2),this%norm(:,:,:,3),this%kappa,x_bc,y_bc,z_bc)
                else if (this%use_D04) then
                  call divergence(this%decomp,this%derD04,this%norm(:,:,:,1),this%norm(:,:,:,2),this%norm(:,:,:,3),this%kappa,x_bc,y_bc,z_bc)
                else		!kappa, divergence of surface normal	
		  call divergence(this%decomp,this%der,this%norm(:,:,:,1),this%norm(:,:,:,2),this%norm(:,:,:,3),this%kappa,x_bc,y_bc,z_bc)	
		endif
	endif

	if (this%use_gradXi) then

               
          where( this%material(1)%VF .GE. 1-this%intSharp_cut)
                 this%xi(:,:,:,1) = this%intSharp_eps*log( ( 1-2*this%intSharp_cut+ e )/ (e))*(1/(1-2*this%intSharp_cut))
          elsewhere( this%material(1)%VF .LE. this%intSharp_cut )
                 this%xi(:,:,:,1)  = this%intSharp_eps*log( ( e) / (1 -2*this%intSharp_cut + e))*(1/(1-2*this%intSharp_cut))
          elsewhere
                 this%xi(:,:,:,1)  = this%intSharp_eps*(1/(1-2*this%intSharp_cut))*log( (this%material(1)%VF -this%intSharp_cut + e )/ (1 - this%intSharp_cut -this%material(1)%VF + e) )
          endwhere

          if(this%use_FV) then
           ! call gradientFV(this,this%xi(:,:,:,1) ,gradVF_FV,dx,dy,dz,periodicx,periodicy, periodicz, this%x_bc, this%y_bc, this%z_bc)
            call gradient(this%decomp,this%derCD06,this%xi(:,:,:,1),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3))
            this%gradVF_FV = gradVF_FV
          else if(this%use_D04) then
            call gradient(this%decomp,this%derD06,this%xi(:,:,:,1),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3))
          else
            call gradient(this%decomp,this%der,this%xi(:,:,:,1),gradxi(:,:,:,1),gradxi(:,:,:,2),gradxi(:,:,:,3))
          endif
          
            !call filter3D(this%decomp, this%gfil, gradxi(:,:,:,1),iflag,x_bc,y_bc,z_bc)
            !call filter3D(this%decomp, this%gfil, gradxi(:,:,:,2),iflag,x_bc,y_bc,z_bc)
            !call filter3D(this%decomp, this%gfil, gradxi(:,:,:,3),iflag,x_bc,y_bc,z_bc)

            this%gradxi = gradxi
 
            !magnitude of surface vector
            GVFmag = sqrt( gradxi(:,:,:,1)**two + gradxi(:,:,:,2)**two +gradxi(:,:,:,3)**two )

            !surface normal
            where (GVFmag < eps)
              this%norm(:,:,:,1) = zero
              this%norm(:,:,:,2) = zero
              this%norm(:,:,:,3) = zero
            elsewhere
              this%norm(:,:,:,1) = gradxi(:,:,:,1) / GVFmag
              this%norm(:,:,:,2) = gradxi(:,:,:,2) / GVFmag
              this%norm(:,:,:,3) = gradxi(:,:,:,3) / GVFmag
            endwhere
        !  endif
        
         ! this%fmask = 0.25*(1 - (tanh((1-2*this%intSharp_cut)*this%xi(:,:,:,1)/(2 *this%intSharp_eps)))**2)-0.5*(1+tanh((1-2*this%intSharp)*this%xi(:,:,:,1)/(2 *this%intSharp_eps)))*this%intSharp_cut+this%intSharp_cut

          this%fmask =(this%material(1)%VF-this%intSharp_cut)*(one-this%intSharp_cut-this%material(1)%VF)

         if(this%use_FV) then

            if(this%use_normInt) then
               call interpolateFV_6(this,this%norm(:,:,:,1),NMint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc,this%y_bc, this%z_bc)
               call interpolateFV_6(this,this%norm(:,:,:,2),NMint(:,:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
               call interpolateFV_6(this,this%norm(:,:,:,3),NMint(:,:,:,:,3),periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)
               this%normFV(:,:,:,1) = NMint(:,:,:,1,1)        
               this%normFV(:,:,:,2) = NMint(:,:,:,2,2)
               this%normFV(:,:,:,3) = NMint(:,:,:,3,3)
               call divergenceFV_6(this,this%normFV,this%kappa,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
            
            elseif(this%use_normFV) then

              call gradFV_N2Fx(this%decomp,this%derStagg,this%xi,gradxi(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
              call gradFV_N2Fy(this%decomp,this%derStagg,this%xi,gradxi(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
              call gradFV_N2Fz(this%decomp,this%derStagg,this%xi,gradxi(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

              !call filter3D(this%decomp, this%gfil,gradxi(:,:,:,1),iflag,x_bc,y_bc,z_bc)
              !call filter3D(this%decomp, this%gfil,gradxi(:,:,:,2),iflag,x_bc,y_bc,z_bc)
              !call filter3D(this%decomp, this%gfil,gradxi(:,:,:,3),iflag,x_bc,y_bc,z_bc)

              this%gradxi = gradxi

              !magnitude of surface vector
              GVFmag = sqrt( gradxi(:,:,:,1)**two + gradxi(:,:,:,2)**two+gradxi(:,:,:,3)**two )

              !surface normal
              where (GVFmag < eps)
                this%norm(:,:,:,1) = zero
                this%norm(:,:,:,2) = zero
                this%norm(:,:,:,3) = zero
              elsewhere
                this%norm(:,:,:,1) = gradxi(:,:,:,1) / GVFmag
                this%norm(:,:,:,2) = gradxi(:,:,:,2) / GVFmag
                this%norm(:,:,:,3) = gradxi(:,:,:,3) / GVFmag
              endwhere

              call divergenceFV_6(this,this%norm,this%kappa,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)
            endif
         else
            call divergence(this%decomp,this%der,this%norm(:,:,:,1),this%norm(:,:,:,2),this%norm(:,:,:,3),this%kappa,x_bc,y_bc,z_bc)
         endif

        endif

        if (this%use_gradphi) then


                where( (this%material(1)%VF .GT. minVF) .AND. (this%material(1)%VF .LT. 1) )
                         this%phi = (this%material(1)%VF**r)/(this%material(1)%VF**r + (1-this%material(1)%VF)**r)
                elsewhere( (this%material(1)%VF .GE. 1) )
                        this%phi = 1
                elsewhere( this%material(1)%VF .LE. minVF)

                        this%phi = minVF**r / ( minVF**r + (1-minVF)**r)
                endwhere

        	call gradient(this%decomp,this%derD02,this%phi,gradphi(:,:,:,1),gradphi(:,:,:,2),gradphi(:,:,:,3)) !high order derivative


         	!calculate surface normal

              !magnitude of surface vector
              GPHImag = sqrt( gradphi(:,:,:,1)**two + gradphi(:,:,:,2)**two + gradphi(:,:,:,3)**two )


              !surface normal
              where (GPHImag < eps)
                 this%norm(:,:,:,1) = eps
                 this%norm(:,:,:,2) = eps
                 this%norm(:,:,:,3) = eps
              elsewhere
                 this%norm(:,:,:,1) = gradphi(:,:,:,1) / GPHImag
                 this%norm(:,:,:,2) = gradphi(:,:,:,2) / GPHImag
                 this%norm(:,:,:,3) = gradphi(:,:,:,3) / GPHImag
              endwhere

                if (this%use_FV) then

                  call interpolateFV(this,this%norm(:,:,:,1),NMint(:,:,:,:,1),periodicx,periodicy,periodicz,-x_bc, y_bc, z_bc)
                  call interpolateFV(this,this%norm(:,:,:,2),NMint(:,:,:,:,2),periodicx,periodicy,periodicz, x_bc,-y_bc, z_bc)
                  call interpolateFV(this,this%norm(:,:,:,3),NMint(:,:,:,:,3),periodicx,periodicy,periodicz, x_bc, y_bc,-z_bc)
                  call interpolateFV(this,this%norm(:,:,:,1),NMint(:,:,:,:,1),periodicx,periodicy,periodicz,-this%x_bc,this%y_bc, this%z_bc)
                  call interpolateFV(this,this%norm(:,:,:,2),NMint(:,:,:,:,2),periodicx,periodicy,periodicz,this%x_bc,-this%y_bc, this%z_bc)
                  call interpolateFV(this,this%norm(:,:,:,3),NMint(:,:,:,:,3),periodicx,periodicy,periodicz,this%x_bc, this%y_bc,-this%z_bc)
                  this%normFV(:,:,:,1) = NMint(:,:,:,1,1)        
                  this%normFV(:,:,:,2) = NMint(:,:,:,2,2)
                  this%normFV(:,:,:,3) = NMint(:,:,:,3,3)

                  call divergenceFV(this,this%normFV,this%kappa,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

                else

                !kappa, divergence of surface normal
                !TODO: double check BCs for div(norm)... do they follow symmetry, or are they opposite?
                call divergence(this%decomp,this%derD02,this%norm(:,:,:,1),this%norm(:,:,:,2),this%norm(:,:,:,3),this%kappa,x_bc,y_bc,z_bc)
                !call filter3D(this%decomp, this%fil, this%kappa, iflag, x_bc, y_bc, z_bc)

                endif

		
	endif


      if (this%weightedcurvature) then

        do n = 1,5
!   print*, this%mpi_rank
!           this%MPI_req = MPI_REQUEST_NULL

            call MPI_Irecv(this%buffer_recieve_1, this%nyp, MPI_real, this%mpi_rank_prev, 1, MPI_Comm_world, this%MPI_req(1), this%ierror)


            call MPI_Irecv(this%buffer_recieve_2, this%nyp, MPI_real, this%mpi_rank_next, 2, MPI_Comm_world, this%MPI_req(2), this%ierror)


            call MPI_Irecv(this%buffer_recieve_k_1, this%nyp, MPI_real, this%mpi_rank_prev, 3, MPI_Comm_world, this%MPI_req(3), this%ierror)


            call MPI_Irecv(this%buffer_recieve_k_2, this%nyp, MPI_real, this%mpi_rank_next, 4, MPI_Comm_world, this%MPI_req(4), this%ierror)
        


            this%buffer_send_1 = this%material(1)%VF(1:1,:,:);
            call  MPI_Issend(this%buffer_send_1,this%nyp, MPI_real, this%mpi_rank_prev, 2, MPI_Comm_world, this%MPI_req(5), this%ierror)

            this%buffer_send_2 = this%material(1)%VF(this%nxp:this%nxp,:, :);
            call  MPI_Issend(this%buffer_send_2,this%nyp, MPI_real, this%mpi_rank_next, 1, MPI_Comm_world, this%MPI_req(6), this%ierror)

            this%buffer_send_k_1 = this%kappa(1:1,:,:);
            call  MPI_Issend(this%buffer_send_k_1,this%nyp, MPI_real, this%mpi_rank_prev, 4, MPI_Comm_world, this%MPI_req(7), this%ierror)

            this%buffer_send_k_2 = this%kappa(this%nxp:this%nxp,:, :);
            call  MPI_Issend(this%buffer_send_k_2,this%nyp, MPI_real, this%mpi_rank_next, 3, MPI_Comm_world, this%MPI_req(8), this%ierror)

            call MPI_Waitall(8, this%MPI_req, this%MPI_Stats, this%ierror)


                do k = 1, this%nzp
                  do j = 2, (this%nyp-1)
                     i = 1


                      weight(i,j,k) = (this%material(1)%VF(i,j,k)*(1-this%material(1)%VF(i,j,k)))**2 + &
                        (this%material(1)%VF(i+1,j,k)*(1-this%material(1)%VF(i+1,j,k)))**2 + &
                        (this%material(1)%VF(i,j+1,k)*(1-this%material(1)%VF(i,j+1,k)))**2 + &
                        (this%buffer_recieve_1(1,j,k)*(1-this%buffer_recieve_1(1,j,k)))**2 + &
                        (this%buffer_recieve_1(1,j+1,k)*(1-this%buffer_recieve_1(1,j+1,k)))**2 +  &
                        (this%buffer_recieve_1(1,j-1,k)*(1-this%buffer_recieve_1(1,j-1,k)))**2 + &
                        (this%material(1)%VF(i+1,j+1,k)*(1-this%material(1)%VF(i+1,j+1,k)))**2 + &
                        (this%material(1)%VF(i+1,j-1,k)*(1-this%material(1)%VF(i+1,j-1,k)))**2 + &
                        (this%material(1)%VF(i,j-1,k)*(1-this%material(1)%VF(i,j-1,k)))**2


                       kappaSum(i,j,k) = this%kappa(i,j,k)*(this%material(1)%VF(i,j,k)*(1-this%material(1)%VF(i,j,k)))**2 + &
                          this%kappa(i+1,j,k)*(this%material(1)%VF(i+1,j,k)*(1-this%material(1)%VF(i+1,j,k)))**2 + &
                          this%kappa(i+1,j-1,k)*(this%material(1)%VF(i+1,j-1,k)*(1-this%material(1)%VF(i+1,j-1,k)))**2 + &
                          this%buffer_recieve_k_1(1,j,k)*(this%buffer_recieve_1(1,j,k)*(1-this%buffer_recieve_1(1,j,k)))**2 + &
                          this%buffer_recieve_k_1(1,j+1,k)*(this%buffer_recieve_1(1,j+1,k)*(1-this%buffer_recieve_1(1,j+1,k)))**2 +  &
                          this%buffer_recieve_k_1(1,j-1,k)*(this%buffer_recieve_1(1,j-1,k)*(1-this%buffer_recieve_1(1,j-1,k)))**2 + &
                          this%kappa(i,j+1,k)*(this%material(1)%VF(i,j+1,k)*(1-this%material(1)%VF(i,j+1,k)))**2 + &
                          this%kappa(i,j-1,k)*(this%material(1)%VF(i,j-1,k)*(1-this%material(1)%VF(i,j-1,k)))**2 + &
                          this%kappa(i+1,j+1,k)*(this%material(1)%VF(i+1,j+1,k)*(1-this%material(1)%VF(i+1,j+1,k)))**2

                 enddo
              enddo

               do k = 1, this%nzp
                  do j = 2, (this%nyp-1)
                     i = this%nxp


                      weight(i,j,k) = (this%material(1)%VF(i,j,k)*(1-this%material(1)%VF(i,j,k)))**2 + &
                        (this%material(1)%VF(i,j-1,k)*(1-this%material(1)%VF(i,j-1,k)))**2 + &
                        (this%material(1)%VF(i,j+1,k)*(1-this%material(1)%VF(i,j+1,k)))**2 + &
                        (this%buffer_recieve_2(1,j,k)*(1-this%buffer_recieve_2(1,j,k)))**2 + &
                        (this%buffer_recieve_2(1,j+1,k)*(1-this%buffer_recieve_2(1,j+1,k)))**2 +  &
                        (this%buffer_recieve_2(1,j-1,k)*(1-this%buffer_recieve_2(1,j-1,k)))**2 + &
                        (this%material(1)%VF(i-1,j-1,k)*(1-this%material(1)%VF(i-1,j-1,k)))**2 + &
                        (this%material(1)%VF(i-1,j+1,k)*(1-this%material(1)%VF(i-1,j+1,k)))**2 + &
                        (this%material(1)%VF(i-1,j,k)*(1-this%material(1)%VF(i-1,j,k)))**2


                       kappaSum(i,j,k) = this%kappa(i,j,k)*(this%material(1)%VF(i,j,k)*(1-this%material(1)%VF(i,j,k)))**2 + &
                          this%kappa(i,j+1,k)*(this%material(1)%VF(i,j+1,k)*(1-this%material(1)%VF(i,j+1,k)))**2 + &
                          this%kappa(i,j-1,k)*(this%material(1)%VF(i,j-1,k)*(1-this%material(1)%VF(i,j-1,k)))**2 + &
                          this%buffer_recieve_k_2(1,j,k)*(this%buffer_recieve_2(1,j,k)*(1-this%buffer_recieve_2(1,j,k)))**2 + &
                          this%buffer_recieve_k_2(1,j+1,k)*(this%buffer_recieve_2(1,j+1,k)*(1-this%buffer_recieve_2(1,j+1,k)))**2 +  &
                          this%buffer_recieve_k_2(1,j-1,k)*(this%buffer_recieve_2(1,j-1,k)*(1-this%buffer_recieve_2(1,j-1,k)))**2 + &
                          this%kappa(i-1,j,k)*(this%material(1)%VF(i-1,j,k)*(1-this%material(1)%VF(i-1,j,k)))**2 + &
                          this%kappa(i-1,j-1,k)*(this%material(1)%VF(i-1,j-1,k)*(1-this%material(1)%VF(i-1,j-1,k)))**2 + &
                          this%kappa(i-1,j+1,k)*(this%material(1)%VF(i-1,j+1,k)*(1-this%material(1)%VF(i-1,j+1,k)))**2

                 enddo
              enddo


        
        do k = 1, this%nzp
          do j = 2, (this%nyp-1)
            do i = 2, (this%nxp-1)
                            
                              
               weight(i,j,k) = (this%material(1)%VF(i,j,k)*(1-this%material(1)%VF(i,j,k)))**2 + &
                        (this%material(1)%VF(i+1,j,k)*(1-this%material(1)%VF(i+1,j,k)))**2 + &
                        (this%material(1)%VF(i-1,j,k)*(1-this%material(1)%VF(i-1,j,k)))**2 + &
                        (this%material(1)%VF(i,j-1,k)*(1-this%material(1)%VF(i,j-1,k)))**2 + &
                        (this%material(1)%VF(i+1,j-1,k)*(1-this%material(1)%VF(i+1,j-1,k)))**2 +  &
                        (this%material(1)%VF(i-1,j-1,k)*(1-this%material(1)%VF(i-1,j-1,k)))**2 + &
                        (this%material(1)%VF(i,j+1,k)*(1-this%material(1)%VF(i,j+1,k)))**2 + &
                        (this%material(1)%VF(i-1,j+1,k)*(1-this%material(1)%VF(i-1,j+1,k)))**2 + &
                        (this%material(1)%VF(i+1,j+1,k)*(1-this%material(1)%VF(i+1,j+1,k)))**2 
               
                
               kappaSum(i,j,k) = this%kappa(i,j,k)*(this%material(1)%VF(i,j,k)*(1-this%material(1)%VF(i,j,k)))**2 + &
                          this%kappa(i+1,j,k)*(this%material(1)%VF(i+1,j,k)*(1-this%material(1)%VF(i+1,j,k)))**2 + &
                          this%kappa(i-1,j,k)*(this%material(1)%VF(i-1,j,k)*(1-this%material(1)%VF(i-1,j,k)))**2 + &
                          this%kappa(i,j-1,k)*(this%material(1)%VF(i,j-1,k)*(1-this%material(1)%VF(i,j-1,k)))**2 + &
                          this%kappa(i+1,j-1,k)*(this%material(1)%VF(i+1,j-1,k)*(1-this%material(1)%VF(i+1,j-1,k)))**2 +  &
                          this%kappa(i-1,j-1,k)*(this%material(1)%VF(i-1,j-1,k)*(1-this%material(1)%VF(i-1,j-1,k)))**2 + &
                          this%kappa(i,j+1,k)*(this%material(1)%VF(i,j+1,k)*(1-this%material(1)%VF(i,j+1,k)))**2 + &
                          this%kappa(i-1,j+1,k)*(this%material(1)%VF(i-1,j+1,k)*(1-this%material(1)%VF(i-1,j+1,k)))**2 + &
                          this%kappa(i+1,j+1,k)*(this%material(1)%VF(i+1,j+1,k)*(1-this%material(1)%VF(i+1,j+1,k)))**2


             enddo
           enddo
         enddo
         
         
               where (weight .LE. eps) 

                   updatedKappa = eps

                elsewhere

                    updatedKappa = kappaSum/weight

               endwhere

          this%kappa = updatedKappa         
        enddo
        

       endif

       
     !  where( (this%material(1)%VF .LE. 1d-6) .OR. (this%material(1)%VF .GE. (1-1d-6)))

     !     this%kappa = 0

     !  endwhere

   !  this%kappa = this%kappa*tanh((this%material(1)%VF*(1-this%material(1)%VF)/0.01)**2)

    if (this%surface_mask) then
	
	 !his%fmask = 1 - (1 - 4*this%material(1)%VF*(1-this%material(1)%VF))**nmask 
	 !this%fmask = tanh( ( ((this%material(1)%VF-minVF)*(one-this%material(1)%VF-minVF)) / tmask)**two )

       ! where (this%material(1)%VF .GT. (one-1d-5))
       !         mask2 = 0
       ! elsewhere (this%material(1)%VF .LE. 1d-5)
       !         mask2 = 0
       ! elsewhere
       !         mask2 = one
       ! endwhere

        !his%maskKappa = this%fmask*this%kappa*mask2
	!TODO: Compute surface tension force and store in this%surfaceTension_f
        this%surfaceTension_f(:,:,:,1) = -6*this%surfaceTension_coeff*this%fmask*this%kappa*gradVF(:,:,:,1)
        this%surfaceTension_f(:,:,:,2) = -6*this%surfaceTension_coeff*this%fmask*this%kappa*gradVF(:,:,:,2)
        this%surfaceTension_f(:,:,:,3) = -6*this%surfaceTension_coeff*this%fmask*this%kappa*gradVF(:,:,:,3)
          
        this%maskKappa = this%fmask*this%kappa 

	else

       this%fmask = 1 - (1 -4*this%material(1)%VF*(1-this%material(1)%VF))**nmask
       where( this%fmask .LE. eps)
          this%fmask = 0
       endwhere
       ! tanhmask =tanh((((this%material(1)%VF-this%intSharp_cut)*(one-this%material(1)%VF-this%intSharp_cut))/this%intSharp_tnh)**two )

     !  this%kappa = this%fmask*this%kappa
        !TODO: Compute surface tension force and store in this%surfaceTension_f
        !this%kappa = one/0.4
	this%surfaceTension_f(:,:,:,1) = -this%surfaceTension_coeff*this%kappa*gradVF(:,:,:,1)
	this%surfaceTension_f(:,:,:,2) = -this%surfaceTension_coeff*this%kappa*gradVF(:,:,:,2)
        this%surfaceTension_f(:,:,:,3) = -this%surfaceTension_coeff*this%kappa*gradVF(:,:,:,3)

	endif

        if(this%energy_surfTen) then

     !      if(this%use_FV) then
     !         call laplacian(this%decomp,this%derD02,this%material(1)%VF, lapVF,[0,0],[0,0],[0,0])
     !         call gradient(this%decomp,this%der,this%material(1)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))        
          ! else if(this%use_D04) then
          !    call laplacian(this%decomp,this%derD04,this%material(1)%VF, lapVF,[0,0],[0,0],[0,0])
     !      else
              call laplacian(this%decomp,this%derCD06,this%material(1)%VF, lapVF, [0,0],[0,0],[0,0]) 
     !      end if
           mu =6*this%surfaceTension_coeff/this%intSharp_eps*this%material(1)%VF*(1-this%material(1)%VF)*(1-2*this%material(1)%VF) &
            - 6*this%surfaceTension_coeff*this%intSharp_eps*lapVF

           this%surfaceTension_f(:,:,:,1) = mu*gradVF(:,:,:,1)
           this%surfaceTension_f(:,:,:,2) = mu*gradVF(:,:,:,2)
           this%surfaceTension_f(:,:,:,3) = mu*gradVF(:,:,:,3)
           this%fmask = lapVF
           this%kappa = mu
        endif

     !  do i = 1, 3

     !    this%surfaceTension_f(:,:,:,i) =this%fmask*this%surfaceTension_f(:,:,:,i)

     !  enddo

        !TODO: Use this%surfaceTension_f to compute this%surfaceTension_e
     this%surfaceTension_e = u*this%surfaceTension_f(:,:,:,1) +v*this%surfaceTension_f(:,:,:,2) +w*this%surfaceTension_f(:,:,:,3) 

     ! call divergence(this%decomp,this%der,u,v,w,divu,x_bc,y_bc,z_bc)
     ! call divergence(this%decomp,this%der,u*this%material(1)%VF,v*this%material(1)%VF,w*this%material(1)%VF,divphiu,x_bc,y_bc,z_bc)
    
     !if( this%use_Stagg) then
     !   call interpolateFV_x(this%decomp,this%interpMid,u,u_int(:,:,:,1),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
     !   call interpolateFV_y(this%decomp,this%interpMid,v,u_int(:,:,:,2),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
     !   call interpolateFV_z(this%decomp,this%interpMid,w,u_int(:,:,:,3),periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         
     !   call divergenceFV_6(this,u_int,udiv,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)      
     !   call divergenceFV_6(this,u_int*VF_int,divuphi,dx,dy,dz,periodicx,periodicy,periodicz,this%x_bc,this%y_bc,this%z_bc)

     !    if( .NOT. this%energy_surfTen) then
     !       if(.NOT. this%surface_mask) then
     !          this%surfaceTension_e = -this%surfaceTension_coeff*this%kappa*(divuphi - this%material(1)%VF*divu)
     !       else

     !           this%surfaceTension_e = -6*this%surfaceTension_coeff*this%fmask*this%kappa*(divuphi - this%material(1)%VF*udiv)
       
     !       endif
     !   else

           ! this%surfaceTension_e = mu*(divuphi - this%material(1)%VF*divu)
     !      this%surfaceTension_e = u*this%surfaceTension_f(:,:,:,1) +v*this%surfaceTension_f(:,:,:,2) +w*this%surfaceTension_f(:,:,:,3)

     !   endif

     !endif
!     call filter3D(this%decomp,this%fil,this%surfaceTension_e,iflag,x_bc,y_bc,z_bc)
!     call filter3D(this%decomp, this%fil,this%surfaceTension_f(:,:,:,1),iflag,x_bc,y_bc,z_bc)
!     call filter3D(this%decomp, this%fil,this%surfaceTension_f(:,:,:,2),iflag,x_bc,y_bc,z_bc)
!     call filter3D(this%decomp, this%fil,this%surfaceTension_f(:,:,:,3),iflag,x_bc,y_bc,z_bc)

!     call filter3D(this%decomp,this%fil,this%gradp(:,:,:,1),iflag,x_bc,y_bc,z_bc)
!     call filter3D(this%decomp,this%fil,this%gradp(:,:,:,2),iflag,x_bc,y_bc,z_bc)
!     call filter3D(this%decomp,this%fil,this%gradp(:,:,:,3),iflag,x_bc,y_bc,z_bc)

        
    end subroutine get_surfaceTension

    subroutine interpolateFV(this,nodes,faces,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)  
        !interpolates from Nodes to faces for finite volume treatment of terms
        !in interface advection
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: nodes
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 3), intent(out) :: faces
        logical, intent(in) :: periodicx,periodicy,periodicz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xbuf,xint
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: zbuf,zint
        integer :: i,j,k
        
        faces = 0.0
        ! i+1/2 faces
        if(this%decomp%xsz(1).gt.one) then
        call transpose_y_to_x(nodes,xbuf,this%decomp)
        call this%interpMid02 % iN2Fx(xbuf,xint,x_bc(1),x_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
        !call this%interpMid % iN2Fx(xbuf,xint,x_bc(1),x_bc(2))
        call transpose_x_to_y(xint,faces(:,:,:,1),this%decomp)
        endif

        ! j+1/2 faces
        if(this%decomp%ysz(2).gt.one) then
        call this%interpMid02 % iN2Fy(nodes,faces(:,:,:,2),y_bc(1),y_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
        !call this%interpMid % iN2Fy(nodes,faces(:,:,:,2),y_bc(1),y_bc(2)) 
        endif

        ! k+1/2 faces
        if(this%decomp%zsz(3).gt.one) then
        call transpose_y_to_z(nodes,zbuf,this%decomp)
        call this%interpMid02 % iN2Fz(zbuf,zint,z_bc(1),z_bc(2)) !TODO: add BCs (only correct if interface is away from boundary
        !call this%interpMid % iN2Fz(zbuf,zint,z_bc(1),z_bc(2))
        call transpose_z_to_y(zint,faces(:,:,:,3),this%decomp)
        endif

    end subroutine interpolateFV


    subroutine interpolateFV_6(this,nodes,faces,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        !interpolates from Nodes to faces for finite volume treatment of terms
        !in interface advection
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: nodes
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 3), intent(out) :: faces
        logical, intent(in) :: periodicx,periodicy,periodicz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        real(rkind),dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xbuf,xint
        real(rkind),dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: zbuf,zint
        integer :: i,j,k
        faces  = 0.0
        ! i+1/2 faces
        if(this%decomp%xsz(1).gt.one) then
        call transpose_y_to_x(nodes,xbuf,this%decomp)
        call this%interpMid % iN2Fx(xbuf,xint,x_bc(1),x_bc(2)) !TODO: add BCs(only correct if interface is away from boundary)
        call transpose_x_to_y(xint,faces(:,:,:,1),this%decomp)
        endif

        ! j+1/2 faces
        if(this%decomp%ysz(2).gt.one) then
        call this%interpMid % iN2Fy(nodes,faces(:,:,:,2),y_bc(1),y_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
        endif

        ! k+1/2 faces
        if(this%decomp%zsz(3).gt.one) then
        call transpose_y_to_z(nodes,zbuf,this%decomp)
        call this%interpMid % iN2Fz(zbuf,zint,z_bc(1),z_bc(2)) !TODO: add BCs(only correct if interface is away from boundary)
        call transpose_z_to_y(zint,faces(:,:,:,3),this%decomp)
        endif
    end subroutine interpolateFV_6


    subroutine divergenceFV(this,faces,nodes,dx,dy,dz,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)  
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 3), intent(in) :: faces
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out) :: nodes
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: tmp
        logical, intent(in) :: periodicx,periodicy,periodicz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        real(rkind), intent(in) :: dx,dy,dz
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xbuf
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xdiv
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ybuf
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ydiv
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: zbuf
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: zdiv
        integer :: i,j,k

        nodes = zero


        ! i nodes
        if(this%decomp%xsz(1).gt.one) then
           call transpose_y_to_x(faces(:,:,:,1),xbuf,this%decomp)
           call this%derStaggd02 % ddxF2N(xbuf,xdiv,x_bc(1),x_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
           !call this%derStagg % ddxF2N(xbuf,xdiv,x_bc(1),x_bc(2))
           call transpose_x_to_y(xdiv,tmp,this%decomp)
           nodes = nodes + tmp
        endif


        ! j nodes
        if(this%decomp%ysz(2).gt.one) then
           call this%derStaggd02 % ddyF2N(faces(:,:,:,2),ydiv,y_bc(1),y_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
           !call this%derStagg % ddyF2N(faces(:,:,:,2),ydiv,y_bc(1),y_bc(2))
           nodes = nodes + ydiv
        endif

        ! k nodes
        if(this%decomp%zsz(3).gt.one) then
           call transpose_y_to_z(faces(:,:,:,3),zbuf,this%decomp)
           call this%derStaggd02 % ddzF2N(zbuf,zdiv,z_bc(1),z_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
           !call this%derStagg % ddzF2N(zbuf,zdiv,z_bc(1),z_bc(2)) 
           call transpose_z_to_y(zdiv,tmp,this%decomp)
           nodes = nodes + tmp
        endif

    end subroutine divergenceFV

  subroutine divergenceFV_6(this,faces,nodes,dx,dy,dz,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y,transpose_y_to_z, transpose_z_to_y
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp, 3), intent(in) :: faces
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out) :: nodes
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: tmp
        logical, intent(in) :: periodicx,periodicy,periodicz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        real(rkind), intent(in) :: dx,dy,dz
        real(rkind),dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xbuf
        real(rkind),dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xdiv
        real(rkind),dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ybuf
        real(rkind),dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ydiv
        real(rkind),dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: zbuf
        real(rkind),dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: zdiv
        integer :: i,j,k

        nodes = zero


        ! i nodes
        if(this%decomp%xsz(1).gt.one) then
           call transpose_y_to_x(faces(:,:,:,1),xbuf,this%decomp)
           call this%derStagg % ddxF2N(xbuf,xdiv,x_bc(1),x_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
           call transpose_x_to_y(xdiv,tmp,this%decomp)
           nodes = nodes + tmp
        endif


        ! j nodes
        if(this%decomp%ysz(2).gt.one) then
           call this%derStagg % ddyF2N(faces(:,:,:,2),ydiv,y_bc(1),y_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
           nodes = nodes + ydiv
        endif

        ! k nodes
        if(this%decomp%zsz(3).gt.one) then
           call transpose_y_to_z(faces(:,:,:,3),zbuf,this%decomp)
           call this%derStagg % ddzF2N(zbuf,zdiv,z_bc(1),z_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
           call transpose_z_to_y(zdiv,tmp,this%decomp)
           nodes = nodes + tmp
        endif

    end subroutine divergenceFV_6

    subroutine gradientFV(this,nodes,faces,dx,dy,dz,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)  
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        use exits,      only: nancheck
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in) :: nodes
        real(rkind), dimension(this%nxp, this%nyp, this%nzp,3,3), intent(out) :: faces
        real(rkind), dimension(this%nxp, this%nyp, this%nzp) :: tmp

        real(rkind), dimension(this%decomp%xsz(1), this%decomp%ysz(2), this%decomp%zsz(3),3,3) :: facesG !GLOBAL -- take the memory hit for simplicity?
        logical, intent(in) :: periodicx,periodicy,periodicz
        integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        real(rkind), intent(in) :: dx,dy,dz
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xbuf
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: xgrad
        !real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),3) :: ybuf
        !real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: ygrad
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: zbuf
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: zgrad
        integer :: i,j,k,is,js,ks,il,jl,kl


        !reduced version -- do not need cross derivatives
        faces = zero

        ! i+1/2 faces
        if(this%decomp%xsz(1).gt.one) then
        call transpose_y_to_x(nodes,xbuf,this%decomp)
        call this%derStaggd02 % ddxN2F(xbuf,xgrad,x_bc(1),x_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
        call transpose_x_to_y(xgrad,faces(:,:,:,1,1),this%decomp)
        endif

        ! j+1/2 faces
        if(this%decomp%ysz(2).gt.one) then
        call this%derStaggd02 % ddyN2F(nodes,faces(:,:,:,2,2),y_bc(1),y_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
        endif

        ! k+1/2 faces
        if(this%decomp%zsz(3).gt.one) then
        call transpose_y_to_z(nodes,zbuf,this%decomp)
        call this%derStaggd02 % ddzN2F(zbuf,zgrad,z_bc(1),z_bc(2)) !TODO: add BCs (only correct if interface is away from boundary)
        call transpose_z_to_y(zgrad,faces(:,:,:,3,3),this%decomp)
        endif

    end subroutine gradientFV


    subroutine get_qmix(this,qx,qy,qz)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: qx,qy,qz

        integer :: i

        qx = zero; qy = zero; qz = zero

        do i = 1,this%ns
            qx = qx + this%material(i)%qi(:,:,:,1)
            qy = qy + this%material(i)%qi(:,:,:,2)
            qz = qz + this%material(i)%qi(:,:,:,3)
        end do

    end subroutine

    subroutine filter(this, iflag,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: iflag
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%filter(iflag,x_bc,y_bc,z_bc)
        end do

       if(this%ns.eq.2) then
           this%material(2)%VF = one - this%material(1)%VF
        elseif(this%ns.gt.2) then
           print*,"fix needed?, SolidMixtureEOS VF filter"
           stop
        endif

    end subroutine

    subroutine filter_g(this, iflag,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: iflag
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%filter_g(iflag,x_bc,y_bc,z_bc)
        end do
        
    end subroutine

    subroutine calculate_source(this,rho,divu,u,v,w,p,Fsource,x_bc,y_bc,z_bc)
        use operators, only: divergence,gradient
        class(solid_mixture), intent(in), target :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,divu,u,v,w,p
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: Fsource
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
      
        type(solid), pointer :: mat1, mat2
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)  :: tmp1,tmp2,tmp3,rhocsq1,rhocsq2

            mat1 => this%material(1); mat2 => this%material(2)

        if(this%pEqb) then

            !call mat1%getSpeciesDensity(rho,tmp1)
            call mat1%hydro%get_sos2(mat1%rhom,mat1%p,tmp2)
            call mat1%elastic%get_sos2(mat1%rhom,tmp2)
            rhocsq1 = mat1%rhom*tmp2

            !call mat2%getSpeciesDensity(rho,tmp1)
            call mat2%hydro%get_sos2(mat2%rhom,mat2%p,tmp3)
            call mat2%elastic%get_sos2(mat2%rhom,tmp3)
            rhocsq2 = mat2%rhom*tmp3

            Fsource = mat1%VF*tmp3*(one - mat2%rhom/mat1%rhom)

            call mat1%get_enthalpy(tmp2)
            call mat2%get_enthalpy(tmp3)

            call divergence(this%decomp,this%der,-mat1%Ji(:,:,:,1),-mat1%Ji(:,:,:,2),-mat1%Ji(:,:,:,3),tmp1,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric

            Fsource = -(rhocsq1 - rhocsq2)*divu*mat1%VF*mat2%VF + (Fsource + mat1%VF*(mat2%hydro%gam-one)*(tmp2-tmp3))*tmp1

            Fsource = Fsource / (mat2%VF*rhocsq1 + mat1%VF*rhocsq2)

        elseif(this%pRelax) then

            call divergence(this%decomp,this%der,mat1%Ji(:,:,:,1),mat1%Ji(:,:,:,2),mat1%Ji(:,:,:,3),tmp1,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric
            Fsource = -tmp1*(mat1%eh + mat1%eel + half*(u*u+v*v+w*w))

            call gradient(this%decomp,this%der,mat2%VF,tmp1,tmp2,tmp3)
            Fsource = Fsource - p*(u*tmp1 + v*tmp2 + w*tmp3)

            call gradient(this%decomp,this%der,p,tmp1,tmp2,tmp3)
            Fsource = Fsource - mat1%VF*mat2%VF/rho*(mat1%rhom-mat2%rhom)*(u*tmp1 + v*tmp2 + w*tmp3)

        endif

    end subroutine 

    subroutine update_g(this,isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc)
        use operators, only: filter3D
        use reductions, only: P_SUM
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w,src
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        if(this%useOneG) then
            !if (this%use_gTg) then
            !    call this%material(1)%update_gTg(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,this%rho0mix,this%mumix,this%yieldmix,this%solidVF)
            !else
                !call this%material(1)%update_g(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,this%rho0mix,this%mumix,this%yieldmix,this%solidVF)
                 call this%material(1)%update_g(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,this%rho0mix,this%mumix,this%yieldmix,this%solidVF)
            !end if
            do imat = 2, this%ns
                this%material(imat)%g = this%material(1)%g
            enddo
        else
            do imat = 1, this%ns
                !if (this%use_gTg) then
                !    call this%material(imat)%update_gTg(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc)
                !else
                   call this%material(imat)%update_g(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc)
                !end if
            end do
        endif

    end subroutine

    subroutine implicit_plastic(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho
        integer :: imat

        if(this%useOneG) then
           call this%material(1)%implicit_plastic(rho,this%rho0mix,this%mumix,this%yieldmix)
           
           do imat = 2, this%ns
              call this%material(1)%implicit_plastic(rho,this%rho0mix,this%mumix,this%yieldmix)
           enddo
        else
           do imat = 1, this%ns
              call this%material(imat)%implicit_plastic(rho)
           end do
        endif
        
    end subroutine implicit_plastic


    subroutine Test_Der(this,x,y,z,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        use decomp_2d, only: transpose_y_to_x,transpose_x_to_y,transpose_y_to_z,transpose_z_to_y
        use operators, only: gradFV_x, gradFV_y,gradFV_z, interpolateFV_x,interpolateFV_y, interpolateFV_z, divergenceFV, laplacian, interpolateFV
        use constants,       only: zero,epssmall,eps,one,two,third,half, pi
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) ::x_half,y_half,ddx,ddy
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: VF_int
        integer :: imat
        real(rkind) :: dx, dy

        dx = x(2,1,1)-x(1,1,1)
        dy = y(1,2,1) -y(1,1,1)
        x_half = x+0.5*dx;
        y_half = y+0.5*dy;


       !call interpolateFV_x(this%decomp,this%interpMid,this%material(1)%VF,this%VF_intx,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
       !call interpolateFV_y(this%decomp,this%interpMid,this%material(1)%VF,this%VF_inty,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
       !call interpolateFV_z(this%decomp,this%interpMid,this%material(1)%VF,this%VF_intz,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        
        call interpolateFV(this%decomp,this%interpMid,this%material(1)%VF,VF_int,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

        this%VF_intx = VF_int(:,:,:,1)
        this%VF_inty = VF_int(:,:,:,2)
        this%VF_intz = VF_int(:,:,:,3)
 
        this%ddx_exact  =   -2*sin(2*x)              !-2*sin(2.0*y)*cos(2.0*x) ! ( 2.0) * cos(2.0*x) *cos(4.0*y) - 3*sin(x)
        this%ddy_exact  =   -2*sin(2*y)              !-2*sin(2.0*x)*cos(2.0*y) !(-4.0) * sin(2.0*x) *sin(4.0*y) + 5*cos(y)
        this%intx_exact =   cos(2*x_half) + cos(2*y) !cos(2.0*x_half)*cos(2.0*y) !sin(2.0*x_half)*cos(4.0*y) + 5*sin(y) + 3*cos(x_half)
        this%inty_exact =   cos(2*x) + cos(2*y_half) !cos(2.0*x)*cos(2.0*y_half) !sin(2.0*x)*cos(4.0*y_half) +5*sin(y_half) + 3*cos(x)
        this%intX_error =   abs(this%intx_exact - this%VF_intx) !abs( sin(2.0*x_half)*cos(4.0*y) +  5*sin(y) + 3*cos(x_half) - this%VF_intx)
        this%intY_error =   abs(this%inty_exact - this%VF_inty) !abs( sin(2.0*x)*cos(4.0*y_half) +5*sin(y_half) + 3*cos(x) - this%VF_intY) 

        call gradFV_x(this%decomp,this%derStagg,this%VF_intx,this%DerX,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call gradFV_y(this%decomp,this%derStagg,this%VF_inty,this%DerY,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call gradFV_z(this%decomp,this%derStagg,this%VF_intz,this%DerZ,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        print *, this%DerZ(1,1,1)

        call divergenceFV(this%decomp,this%derStagg,this%VF_intx, this%VF_inty, this%VF_intz,this%DivTest,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        print *, this%DivTest(1,1,1)
        call laplacian(this%decomp,this%derCD06, this%material(1)%VF, this%lapTest, x_bc,y_bc, z_bc)
        
        this%lap_error = abs( -20*sin(2*x)*cos(4*y)-3*cos(x) - 5*sin(y) - this%lapTest)
        this%div_error  = abs((this%ddx_exact + this%ddy_exact) - this%DivTest)
        this%derX_error = abs(this%DerX - this%ddx_exact)
        this%derY_error = abs(this%DerY - this%ddy_exact)
    end subroutine

    subroutine Test_Der_N2F(this,x,y,z,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        use decomp_2d, only: transpose_y_to_x,transpose_x_to_y,transpose_y_to_z,transpose_z_to_y
        use operators, only: gradFV_N2Fx, gradFV_N2Fy,gradFV_N2Fz
        use constants,       only: zero,epssmall,eps,one,two,third,half, pi
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: x_half,y_half,ddx,ddy, testfunc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: VF_int
        integer :: imat
        real(rkind) :: dx, dy

        dx = x(2,1,1)-x(1,1,1)
        dy = y(1,2,1) -y(1,1,1)
        x_half = x+0.5*dx;
        y_half = y+0.5*dy;

        this%ddx_exact = ( 2.0) * cos(2.0*x_half) *cos(4.0*y) -3*sin(x_half)
        this%ddy_exact =  (-4.0) * sin(2.0*x) *sin(4.0*y_half) + 5*cos(y_half)

        testfunc = sin(2.0*x)*cos(4.0*y) + 5*sin(y) + 3*cos(x)
        this%lap_error = testfunc

        call gradFV_N2Fx(this%decomp,this%derStagg,testfunc,this%DerX,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call gradFV_N2Fy(this%decomp,this%derStagg,testfunc,this%DerY,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call gradFV_N2Fz(this%decomp,this%derStagg,testfunc,this%DerZ,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
         
        this%derX_error = abs(this%DerX - this%ddx_exact)
        this%derY_error = abs(this%DerY - this%ddy_exact)


    end subroutine
  
    subroutine Test_Der_NP(this,x,y,z,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        use decomp_2d, only: transpose_y_to_x,transpose_x_to_y,transpose_y_to_z,transpose_z_to_y
        use operators, only: gradFV_N2Fx, gradFV_N2Fy,gradFV_N2Fz
        use constants,       only: zero,epssmall,eps,one,two,third,half, pi
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: x_half,y_half,ddx,ddy, testfunc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: VF_int
        integer :: imat
        real(rkind) :: dx, dy

        dx = x(2,1,1)-x(1,1,1)
        dy = y(1,2,1) -y(1,1,1)
        x_half = x+0.5*dx;
        y_half = y+0.5*dy;

        this%ddx_exact = ( 3.0) * cos(3.0*x_half) *cos(5.0*y) -3*sin(x_half)
        this%ddy_exact =  (-5.0) * sin(3.0*x) *sin(5.0*y_half) + 5*cos(y_half)
        this%ddx_exact(this%nxp,:,:) = 0
        this%ddy_exact(:,this%nyp,:) = 0       
        testfunc = sin(3.0*x)*cos(5.0*y) + 5*sin(y) + 3*cos(x)
        this%lap_error = testfunc
 
        call gradFV_N2Fx(this%decomp,this%derStagg,testfunc,this%DerX,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call gradFV_N2Fy(this%decomp,this%derStagg,testfunc,this%DerY,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call gradFV_N2Fz(this%decomp,this%derStagg,testfunc,this%DerZ,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

        this%derX_error = abs(this%DerX - this%ddx_exact)
        this%derY_error = abs(this%DerY - this%ddy_exact)


    end subroutine
 
    subroutine Test_Der_Periodic(this,x,y,z,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        use decomp_2d, only: transpose_y_to_x,transpose_x_to_y,transpose_y_to_z,transpose_z_to_y
        use operators, only: gradFV_x, gradFV_y,gradFV_z,interpolateFV_x,interpolateFV_y, interpolateFV_z, divergenceFV, laplacian, gradient
        use constants,       only: zero,epssmall,eps,one,two,third,half, pi
        use exits,           only: GracefulExit
        use reductions, only : P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical :: periodicx,periodicy,periodicz
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: x_half,y_half,ddx,ddy, testfuncY, testfuncX,intFunc
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: TestF,DerF
        integer :: imat
        real(rkind) :: dx, dy

        dx = x(2,1,1)-x(1,1,1)
        dy = y(1,2,1) -y(1,1,1)
        x_half = x+0.5*dx;
        y_half = y+0.5*dy;
        print *, "in Test"
        !testfuncX = sin(5*y) + sin(7*x_half)
        !testfuncY = sin(5*y_half) + sin(7*x)
        intFunc = y*exp(-x) + x*exp(-y)                                  !sin(5*y) + sin(7*x)
        testfuncX = y*exp(-x_half) + x_half*exp(-y)                      !sin(5*y) + sin(7*x_half)
        testfuncY = y_half*exp(-x) + x*exp(-y_half)
 
        TestF(:,:,:,1) = testfuncX
        TestF(:,:,:,2) = testfuncY
        TestF(:,:,:,3) = 0.0
        call interpolateFV_x(this%decomp,this%interpMid,intFunc,this%VF_intx,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call interpolateFV_y(this%decomp,this%interpMid,intFunc,this%VF_inty,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call interpolateFV_z(this%decomp,this%interpMid,intFunc,this%VF_intz,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)

        this%ddx_exact  = -y*exp(-x) + exp(-y)  !7.0*cos(7.0*x) ! ( 2.0) * cos(2.0*x) *cos(4.0*y) - 3*sin(x)
        this%ddy_exact  = -x*exp(-y) + exp(-x)  !5.0*cos(5.0*y)  !(-4.0) * sin(2.0*x) *sin(4.0*y) + 5*cos(y)
        this%intx_exact =  testfuncX !sin(2.0*x_half)*cos(4.0*y) + 5*sin(y) + 3*cos(x_half)
        this%inty_exact =  testfuncY !sin(2.0*x)*cos(4.0*y_half) +5*sin(y_half) + 3*cos(x)
        this%intX_error =   abs(this%intx_exact - this%VF_intx) !abs( sin(2.0*x_half)*cos(4.0*y) +  5*sin(y) + 3*cos(x_half) - this%VF_intx)
        this%intY_error =   abs(this%inty_exact - this%VF_inty) !abs( sin(2.0*x)*cos(4.0*y_half) +5*sin(y_half) + 3*cos(x) - this%VF_intY)
        call gradFV_x(this%decomp,this%derStagg,testfuncX,this%DerX,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        print *, "gradVF x"
        call gradFV_y(this%decomp,this%derStagg,testfuncY,this%DerY,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        print *, "gradVF_y"
        call gradFV_z(this%decomp,this%derStagg,testfuncX,this%DerZ,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        !call gradient(this%decomp, this%der, TestF, this%DerX,this%DerY,this%DerZ, x_bc, y_bc, z_bc)

        call divergenceFV(this%decomp,this%derStagg,this%VF_intx, this%VF_inty,this%VF_intz,this%DivTest,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        call laplacian(this%decomp,this%derCD06, this%material(1)%VF,this%lapTest, x_bc,y_bc, z_bc)
        
        this%div_error  = abs((this%ddx_exact + this%ddy_exact) - this%DivTest)
        this%derX_error = abs(this%DerX - this%ddx_exact)
        this%derY_error = abs(this%DerY - this%ddy_exact)


    end subroutine


    subroutine update_Ys(this,isub,dt,rho,u,v,w,x,y,z,tsim,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical :: periodicx,periodicy,periodicz
        integer :: imat
        do imat = 1, this%ns
          call this%material(imat)%update_Ys(isub,dt,rho,u,v,w,x,y,z,tsim,periodicx, periodicy, periodicz,x_bc,y_bc,z_bc)
        end do
    end subroutine

    subroutine get_pmix(this,p)
        class(solid_mixture), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: p  ! Mixture pressure

        integer :: i

        p = zero
        do i = 1,this%ns
            p = p + this%material(i)%VF * this%material(i)%p  ! Volume fraction weighted sum
        end do

    end subroutine

    subroutine get_emix(this,e)
        class(solid_mixture), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: e  ! Mixture internal energy

        integer :: i

        e = zero
        do i = 1,this%ns
            e = e + this%material(i)%Ys * ( this%material(i)%eh + this%material(i)%eel )  ! Mass fraction weighted sum
        end do

    end subroutine

    subroutine update_eh(this,isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,src,taustar,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,u,v,w,divu,viscwork,src
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(in) :: taustar
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        !do imat = 1, this%ns
        !  call this%material(imat)%update_eh(isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,src,taustar,x_bc,y_bc,z_bc)
        !end do
        call this%material(1)%update_eh(isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork, src,taustar,x_bc,y_bc,z_bc)
        call this%material(2)%update_eh(isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,-src,taustar,x_bc,y_bc,z_bc)

    end subroutine

    subroutine getSOS(this,rho,p,sos)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho, p  ! Mixture density and pressure
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: sos     ! Mixture speed of sound

        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhom, sosm,gam,pgam,gam_m,pinf_m
!print *, '----In mix%getSOS----'
        sos = zero
        do i = 1,this%ns
            call this%material(i)%getSpeciesDensity(rho,rhom)
            call this%material(i)%hydro%get_sos2(rhom,p,sosm)
            !i:qf(.not. this%pEqb) then
                call this%material(i)%elastic%get_sos2(rhom,sosm)
                !print *, 'mat  dens', i, rhom(89,1,1), p(89,1,1), this%material(i)%Ys(89,1,1), sosm(89,1,1)
            !endif
            if(this%SOSmodel) then
                ! equilibrium model
                sos = sos + this%material(i)%VF/(rhom*sosm)
            else
                ! frozen model (details in Saurel et al. 2009)
                sos = sos + this%material(i)%Ys*sosm
            endif
        end do

        if(this%SOSmodel) then
            sos = one / (sqrt(rho*sos) + epssmall)
        else
            sos = sqrt(abs(sos))
        endif

        if(this%pEqb) then

           sos = zero
           gam = zero
           pgam = zero
           do i = 1, this%ns
              gam = gam + this%material(i)%VF*this%material(i)%hydro%onebygam_m1
              pgam = pgam +this%material(i)%VF*this%material(i)%hydro%onebygam_m1*this%material(i)%hydro%gam*this%material(i)%hydro%PInf
           !   sos = sos + this%material(i)%Ys*sosm
           enddo
           !sos = sqrt(abs(sos)) 
           !gam_m = 1/gam - 1
           !pinf_m = pgam /gam
           sos = sqrt(((gam+1)*p + pgam)/(gam*rho)) 

        end if
        !if(this%useOneG) then
        !    sosm = sos*sos
        !    call this%material(1)%elastic%get_sos2_mixture(rho,this%mumix,sosm)
        !    sos = sqrt(sosm)
        !    print *, 'mixt dens', rho(89,1,1), this%mumix(89,1,1)
        !endif
!print *, '----Exiting mix%getSOS----'
    end subroutine

    subroutine update_VF(this,isub,dt,rho,u,v,w,x,y,z,tsim,divu,src,periodicx, periodicy, periodicz,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w,divu,src
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
        logical :: periodicx,periodicy,periodicz
        integer :: imat

        if (this%ns > 2) call GracefulExit("Figure out 2 materials first!",4356)

        !do imat = 1, this%ns
        !  call this%material(imat)%update_VF(this%material( 2-mod(imat+1,2) ),isub,dt,rho,u,v,w,x,y,z,tsim,divu,src,x_bc,y_bc,z_bc)
        !end do
        call this%material(1)%update_VF(this%material(1),isub,dt,rho,u,v,w,x,y,z,tsim,divu,src,periodicx,periodicy,periodicz,x_bc,y_bc,z_bc)
        !call this%material(2)%update_VF(this%material(1),isub,dt,rho,u,v,w,x,y,z,tsim,divu,-src,x_bc,y_bc,z_bc)

        this%material(2)%VF = one - this%material(1)%VF

    end subroutine

    subroutine checkNaN(this)
        class(solid_mixture), intent(in)   :: this
        integer :: imat

        do imat=1,this%ns
            call this%material(imat)%checkNaN(imat)
        end do

    end subroutine

    subroutine fnumden(this,pf,fparams,iparams,num,den)
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(in)               :: pf
        real(rkind), intent(in), dimension(:), target :: fparams
        integer, intent(in), dimension(:)     :: iparams
        real(rkind), intent(out)              :: num, den

        integer :: im, i, j, k
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf
        real(rkind) :: fac, gm1, YCv, pinfloc, rhoE

        if(iparams(1) == 1) then
            ! relaxPressure
            vf   => fparams(  1:this%ns)
            gam  => fparams(  this%ns+1:2*this%ns)
            psph => fparams(2*this%ns+1:3*this%ns)
            pinf => fparams(3*this%ns+1:4*this%ns)
            
            num = zero; den = zero;
            do im = 1, this%ns
              fac = vf(im)/gam(im)/(pinf(im)+pf)
              num = num + fac*(psph(im)-pf)
              !den = den + fac*(psph(im)-pinf(im)-two*pf)/(pinf(im)+pf)
              den = den - fac*(psph(im)+pinf(im))/(pinf(im)+pf)
              !write(*,*) im, psph(im), pf!-pinf(im)-two*pf, -(psph(im)+pinf(im))
            enddo
            nullify(vf,gam,psph,pinf)
        elseif(iparams(1)==2) then
            ! equilibratePressureTemperature

            i = iparams(2); j = iparams(3); k = iparams(4)
            num = zero; den = zero
            do im = 1, this%ns
              gm1 = this%material(im)%hydro%gam-one
              YCv = this%material(im)%Ys(i,j,k) * this%material(im)%hydro%Cv
              pinfloc = this%material(im)%hydro%PInf
              rhoE = fparams(1)

              num = num + YCv*(gm1*rhoE-(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
              den = den + YCv*gm1*(pinfloc-rhoE) / (pf+pinfloc)**2
            enddo

        endif

    end subroutine


    subroutine fnumden_new(this,pf,fparams,iparams,num,den)
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(in)               :: pf
        real(rkind), intent(in), dimension(:), target :: fparams
        integer, intent(in), dimension(:)     :: iparams
        real(rkind), intent(out)              :: num, den

        integer :: im, i, j, k
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf
        real(rkind) :: fac, gm1, YCv, pinfloc, rhoE

        if(iparams(1) == 1) then
            ! relaxPressure
            vf   => fparams(  1:this%ns)
            gam  => fparams(  this%ns+1:2*this%ns)
            psph => fparams(2*this%ns+1:3*this%ns)
            pinf => fparams(3*this%ns+1:4*this%ns)
            
            num = zero; den = zero;
            do im = 1, this%ns
              fac = vf(im)/gam(im)/(pinf(im)+pf)
              num = num + fac*(psph(im)-pf)
              !den = den + fac*(psph(im)-pinf(im)-two*pf)/(pinf(im)+pf)
              den = den - fac*(psph(im)+pinf(im))/(pinf(im)+pf)
              !write(*,*) im, psph(im), pf!-pinf(im)-two*pf, -(psph(im)+pinf(im))
            enddo
            nullify(vf,gam,psph,pinf)
        elseif(iparams(1)==2) then
            ! equilibratePressureTemperature

            i = iparams(2); j = iparams(3); k = iparams(4)
            num = zero; den = zero
            do im = 1, this%ns
              gm1 = this%material(im)%hydro%gam-one

              YCv = this%material(im)%Ys(i,j,k) * this%material(im)%hydro%Cv !original
              !YCv = min(max(this%material(im)%Ys(i,j,k),zero),one) * this%material(im)%hydro%Cv !new
              !YCv = min(max(this%material(im)%Ys(i,j,k),zero+1.D-6),one-1.D-6) * this%material(im)%hydro%Cv !new
              !YCv = abs(this%material(im)%Ys(i,j,k)) * this%material(im)%hydro%Cv !new



              pinfloc = this%material(im)%hydro%PInf
              rhoE = fparams(1)

              num = num + YCv*(gm1*rhoE-(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
              den = den + YCv*gm1*(pinfloc-rhoE) / (pf+pinfloc)**2

              !new -- den not needed for bisection
              !num = num + this%material(im)%hydro%Cv*( abs(this%material(im)%Ys(i,j,k)) *gm1*rhoE - this%material(im)%Ys(i,j,k) *(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
              !num = num + this%material(im)%hydro%Cv*( this%material(im)%Ys(i,j,k) *gm1*rhoE - abs(this%material(im)%Ys(i,j,k)) *(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
              !num = num + this%material(im)%hydro%Cv*( this%material(im)%Ys(i,j,k) *gm1*rhoE - min(max(this%material(im)%Ys(i,j,k),zero),one) *(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
              !num = num + this%material(im)%hydro%Cv*( min(max(this%material(im)%Ys(i,j,k),zero),one) *gm1*rhoE - this%material(im)%Ys(i,j,k) *(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
              !num = num + this%material(im)%hydro%Cv*( abs(this%material(im)%Ys(i,j,k)) *gm1*rhoE - min(max(this%material(im)%Ys(i,j,k),zero),one) *(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
              !num = num + this%material(im)%hydro%Cv*( min(max(this%material(im)%Ys(i,j,k),zero),one) *gm1*rhoE - abs(this%material(im)%Ys(i,j,k)) *(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
     
              !num = num + this%material(im)%hydro%Cv*( this%material(im)%Ys(i,j,k) *gm1*rhoE - min(max(this%material(im)%Ys(i,j,k),zero),one) *pf - this%material(im)%Ys(i,j,k) *this%material(im)%hydro%gam*pinfloc)/(pf+pinfloc)  --ok
              !num = num + this%material(im)%hydro%Cv*( min(max(this%material(im)%Ys(i,j,k),zero),one) *gm1*rhoE - this%material(im)%Ys(i,j,k) *pf - min(max(this%material(im)%Ys(i,j,k),zero),one) *this%material(im)%hydro%gam*pinfloc)/(pf+pinfloc) --bad
              !num = num + this%material(im)%hydro%Cv*( min(max(this%material(im)%Ys(i,j,k),zero),one) *gm1*rhoE - abs(this%material(im)%Ys(i,j,k)) *pf - min(max(this%material(im)%Ys(i,j,k),zero),one) *this%material(im)%hydro%gam*pinfloc)/(pf+pinfloc) --bad
              !num = num + this%material(im)%hydro%Cv*( abs(this%material(im)%Ys(i,j,k)) *gm1*rhoE - min(max(this%material(im)%Ys(i,j,k),zero),one) *pf - abs(this%material(im)%Ys(i,j,k)) *this%material(im)%hydro%gam*pinfloc)/(pf+pinfloc) -- long but pressure oscillations
              !num = num + this%material(im)%hydro%Cv*( this%material(im)%Ys(i,j,k) *gm1*rhoE - abs(this%material(im)%Ys(i,j,k)) *pf - this%material(im)%Ys(i,j,k) *this%material(im)%hydro%gam*pinfloc)/(pf+pinfloc) -- long but bad
              !num = num + this%material(im)%hydro%Cv*( abs(this%material(im)%Ys(i,j,k)) *gm1*rhoE - this%material(im)%Ys(i,j,k) *pf - abs(this%material(im)%Ys(i,j,k)) *this%material(im)%hydro%gam*pinfloc)/(pf+pinfloc)
            enddo

        endif

    end subroutine


    subroutine rootfind_nr_1d(this,pf,fparams,iparams)
        use decomp_2d, only: nrank
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(inout)            :: pf      !equilibrium pressure
        real(rkind), intent(in), dimension(:) :: fparams !1: mixture hydrostatic (internal-elastic) energy
        integer, intent(in), dimension(:)     :: iparams !1: PTeqb flag, 2-4: i,j,k
    
        integer     :: ii, imat, itmax = 1000
        !real(rkind) :: tol = 1.0d-8
        real(rkind) :: dpf, num, den, den_conv, urlx,rhomin,tmp
    
        !pfinitguess = pf
        do ii = 1, itmax
          call this%fnumden(pf,fparams,iparams,num,den)
          ! if(dabs(den)>1.0d-12) then
          if(dabs(den)>1.0d-15) then
             dpf = num/den !original
          else
            write(*,*) 'den very small, please check.', num, num/den
            !write(*,*) 'failure at proc ', nrank, ' at index ', iparams(2:4)
            write(*,*) 'interation ',ii,'failure at proc ', nrank, ' at index ', iparams(2)+this%decomp%yst(1)-1,iparams(3)+this%decomp%yst(2)-1,iparams(4)+this%decomp%yst(3)-1
            stop
          endif

          pf = pf - dpf !original

          ! check for convergence
          if(dabs(pf)>1.0d-12) then
            den_conv = dabs(pf)
          else
            den_conv = one
          endif
          if(dabs(dpf)/den_conv<1.0d-8) exit
        enddo
        if(ii==itmax+1) then
           write(*,*) 'Newtons method for pf did not converge. Check details.', iparams(1)
        endif
  
    end subroutine rootfind_nr_1d

    subroutine thick_calculations(this,rho,u,p,pthick,uthick,rhothick,Ys_thick,VF_thick,Ys_wiggle,VF_wiggle,dx)
        use decomp_2d, only: nrank
        use reductions, only : P_MINVAL,P_MAXVAL
        use decomp_2d, only: transpose_y_to_x,transpose_x_to_y,transpose_y_to_z,transpose_z_to_y
        use operators, only: gradFV_N2Fx, gradFV_N2Fy,gradFV_N2Fz,gradient
        use constants, only: zero,epssmall,eps,one,two,third,half, pi
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(in), dimension(:,:,:) :: rho,u,p      !equilibrium pressure
        real(rkind), intent(inout) :: pthick,uthick,rhothick,dx,Ys_wiggle,VF_wiggle,Ys_thick,VF_thick  !1: mixture hydrostatic(internal-elastic) energy
        integer     :: i, imat
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,3) :: gradYs, gradVF,gradRho,gradU,gradP
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)   :: H1a, H2a,H1Ys,H2Ys
        real(rkind) :: dpf, num, den, num_p, den_p, den_conv, pf_p, pf_n,pst,ps_n,num_n,den_n,pf_no,pf_po,num_po,num_no,den_po,den_no
        real(rkind), dimension(2) :: Yswig,VFwig
        do i = 1, this%ns
          
          call gradient(this%decomp,this%der,this%material(i)%Ys,gradYs(:,:,:,1),gradYs(:,:,:,2),gradYs(:,:,:,3))
          call gradient(this%decomp,this%der,this%material(i)%VF,gradVF(:,:,:,1),gradVF(:,:,:,2),gradVF(:,:,:,3))
          Ys_thick = 1.0 / (dx * P_MAXVAL(gradYs))
          VF_thick = 1.0 / (dx * P_MAXVAL(gradVF))

          where( (this%material(i)%VF - 1) .LE. 0)

              H1a = 0

          elsewhere

              H1a = 1.0

          endwhere

          where( (-1*this%material(i)%VF ) .LE. 0)

              H2a = 0

          elsewhere

              H2a = 1.0

          endwhere

          VFwig(1) = P_MAXVAL(abs(this%material(i)%VF - 1)*H1a)
          VFwig(2) = P_MAXVAL(abs(this%material(i)%VF)*H2a )
          VF_wiggle = 100*maxval(VFwig )

          where( (this%material(i)%Ys - 1) .LE. 0)

              H1Ys = 0

          elsewhere

              H1Ys = 1.0

          endwhere

          where( (-1*this%material(i)%Ys ) .LE. 0)

              H2Ys = 0

          elsewhere

              H2Ys = 1.0

          endwhere

          Yswig(1) = P_MAXVAL(abs(this%material(i)%Ys - 1)*H1a)
          Yswig(2) = P_MAXVAL(abs(this%material(i)%Ys)*H2a)
          Ys_wiggle = 100*maxval(Yswig)
        enddo        
        
        call gradient(this%decomp,this%der,p,gradP(:,:,:,1),gradP(:,:,:,2),gradP(:,:,:,3))
        call gradient(this%decomp,this%der,u,gradU(:,:,:,1),gradU(:,:,:,2),gradU(:,:,:,3))
        call gradient(this%decomp,this%der,rho,gradRho(:,:,:,1),gradRho(:,:,:,2),gradRho(:,:,:,3))
        pthick = (P_MAXVAL(p) - P_MINVAL(p)) / (dx*P_MAXVAL(gradP))
        uthick = (P_MAXVAL(u) - P_MINVAL(u)) / (dx*P_MAXVAL(gradU)) 
        rhothick = (P_MAXVAL(rho) - P_MINVAL(rho)) / (dx*P_MAXVAL(gradRho))        

    end subroutine

    subroutine rootfind_nr_1d_new(this,pf,fparams,iparams,pmin,pmax,icount,icount2,isub,nsubs) !,rhomin,pmin,pmax)
        use decomp_2d, only: nrank
        use reductions, only : P_MINVAL
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(inout)            :: pf      !equilibrium pressure
        real(rkind), intent(in), dimension(:) :: fparams !1: mixture hydrostatic (internal-elastic) energy
        real(rkind), intent(in)               :: pmin,pmax ! rhomin
        integer, intent(in), dimension(:)     :: iparams !1: PTeqb flag, 2-4: i,j,k
        integer, intent(inout)                :: icount,icount2
        integer, intent(in)                   :: isub,nsubs
    
        integer     :: ii, imat, itmax = 10000,iter
        real(rkind) :: dpf, num, den, num_p, den_p, den_conv, pf_p, pf_n, pst,ps_n,num_n,den_n,pf_no,pf_po,num_po,num_no,den_po,den_no
        character(len=100) :: str

        ! !set bisection bounds
        ! ! !pst = 1.D3
        ! ! pf_p = pmax*pst          !positive bound for bisection
        ! ! !pf_n = max(pmin/pst,eps) !negative bound for bisection -- original -- could make min > eps for better stability
        ! ! pf_n = max(max(pmin/pst,this%intSharp_pfloor),eps) !negative bound for bisection
        ! ! !pf_n = -0.5*pf_p !test
        ! ! ps_n = pf_n

        ! pst = pf!1.1
        ! ! pf_p = min(pmax*pst,pf*pst)          !positive bound for bisection
        ! ! pf_n = max(max(pmin/pst,pf/pst),eps) !negative bound for bisection

        ! !pf_p = pmax*pst          !positive bound for bisection
        ! pf_p = pmax*1.D8          !positive bound for bisection
        ! pf_n = max(pmin/1.1,eps) !negative bound for bisection

        ! !make sure a root exists in between bounds
        ! num_n = one; num_p = one
        ! call this%fnumden_new(pf_n,fparams,iparams,num_no,den_no)!initialize num_no
        ! call this%fnumden_new(pf_p,fparams,iparams,num_po,den_po)!initialize num_po
        ! do while( (num_n*num_p.GT.zero) )
        !    call this%fnumden_new(pf_n,fparams,iparams,num_n,den_n)
        !    call this%fnumden_new(pf_p,fparams,iparams,num_p,den_p)
        !    if (num_n*num_p.GT.zero) then !likely there are two roots (or unlikely no roots) in between bounds
        !       pf_p = pf_p/1.01
        !       pf_n = pf_n*1.01
        !       ! call this%fnumden_new(pf_n,fparams,iparams,num_n,den_n)
        !       ! call this%fnumden_new(pf_p,fparams,iparams,num_p,den_p)
        !       if (num_p*num_po.LT.zero) then !root in between p and po
        !          pf_p = pf_po
        !          pf_n = pf_p
        !          exit
        !       elseif (num_n*num_no.LT.zero) then !root in between n and no
        !          pf_p = pf_n
        !          pf_n = pf_no
        !          exit
        !       endif
        !          pf_po = pf_p
        !          pf_no = pf_n
        !          num_po = num_p
        !          num_no = num_n
        !    else
        !       exit
        !    endif
        !    if(pf_p.le.pf_n) then
        !    endif
        !    if(pf_p.le.pf_n) then
        !       print*, 'root not found in bisection'
        !       ! pf_p = pf_po
        !       ! pf_n = pf_no
        !       pf_p = pst
        !       pf_n = pst
        !       exit
        !    endif
        ! enddo

        
        pst = 10.
        pf_p = pmax*pst          !positive bound for bisection
        pf_n = max(max(pmin/pst,this%intSharp_pfloor),eps) !negative bound for bisection -- original -- could make min > eps for better stability


        !pf_p = pf*pst          !positive bound for bisection
        !pf_n = max(pf/pst,eps) !negative bound for bisection -- original -- could make min > eps for better stability

        ps_n = pf_n

        !print warning when p negative before relaxation
        !if(pf.lt.zero) then
        if(pf.le.eps) then
           icount = icount + one
           ! if((icount .le. two).and.(isub.eq.one)) then
           if((icount .le. 2)) then
              write(str,'(1(ES8.1E2,1X),3(I5,1X),(I1))') pf,iparams(2)+this%decomp%yst(1)-1,iparams(3)+this%decomp%yst(2)-1,iparams(4)+this%decomp%yst(3)-1,isub
              print*,'p negative before relaxation '//trim(str) !set larger intSharp_pfloor to correct -- default zero
              if(icount .eq. 2) then
                 print*, 'not printing additional negative p before relaxation'
              endif
           endif
        endif

        ! !good new
        ! do ii = 1, itmax

        !    !bisection iteration
        !    pf = 0.5*(pf_p+pf_n)

        !    call this%fnumden_new(pf,fparams,iparams,num,den)
        !    call this%fnumden_new(pf_p,fparams,iparams,num_p,den_p)
        !    call this%fnumden_new(pf_n,fparams,iparams,num_n,den_n)
        !    !if (sign(one,num).EQ.sign(one,num_p)) then
        !    if (num*num_p.LT.zero) then !root exists between pf and pf_p
        !       pf_n = pf
        !    elseif (num*num_n.LT.zero) then !root exists between pf and pf_n
        !       pf_p = pf
        !    else !zero or two roots exist between pf and pf_p (or pf and pf_n) -- first assume zero then minimuze abs(num)

        !       if (den*den_p.LT.zero) then !min exists between pf and pf_p
        !          !pf_n = pf
        !          pf_n = 0.5*(pf+pf_n)
        !       elseif (den*den_n.LT.zero) then !min exists between pf and pf_n
        !          !pf_p = pf
        !          pf_p = 0.5*(pf+pf_p)
        !       else !two roots or minima exist between pf_n and pf_p -- tighten bounds
        !          pf_po = pf_p
        !          pf_no = pf_n
        !          num_po = num_p
        !          num_no = num_n
        !          pf_p = pf_p-0.01*(pf_p-pf_n)
        !          pf_n = pf_n+0.01*(pf_p-pf_n)
        !          call this%fnumden_new(pf_p,fparams,iparams,num_p,den_p)
        !          call this%fnumden_new(pf_n,fparams,iparams,num_n,den_n)
        !          if(abs(num_p).ge.abs(num_po)) then
        !             pf_p = pf_po
        !             num_p = num_po
        !          endif
        !          if(abs(num_n).ge.abs(num_no)) then
        !             pf_n = pf_no
        !             num_n = num_no
        !          endif
           
        !       endif
        !    endif

       

        !       !    !do while( (den*den_p.LT.zero) .and. (den*den_n.LT.zero))
        !       ! !do while( (den*den_p.LT.zero) .and. (den*den_n.LT.zero))
        !       ! !do while( (num*num_p.GT.zero) .and. (num*num_n.GT.zero))
        !       ! do while (.TRUE.)
        !       !    ! pf_p = pf_p-0.01*pf
        !       !    ! pf_n = pf_n+0.01*pf
        !       !    pf_po = pf_p
        !       !    pf_no = pf_n
        !       !    num_po = num_p
        !       !    num_no = num_n
        !       !    pf_p = pf_p-0.01*(pf_p-pf_n)
        !       !    pf_n = pf_n+0.01*(pf_p-pf_n)
        !       !    call this%fnumden_new(pf_p,fparams,iparams,num_p,den_p)
        !       !    call this%fnumden_new(pf_n,fparams,iparams,num_n,den_n)
        !       !    if(abs(num_p).ge.abs(num_po)) then
        !       !       pf_p = pf_po
        !       !       num_p = num_po
        !       !    endif
        !       !    if(abs(num_n).ge.abs(num_no)) then
        !       !       pf_n = pf_no
        !       !       num_n = num_no
        !       !    endif

        !       !    if(num*num_p.LT.zero) then
        !       !       pf_p = pf_po
        !       !       pf_n = pf_p
        !       !       exit
        !       !    endif
        !       !    if(num*num_n.LT.zero) then
        !       !       pf_p = pf_n
        !       !       pf_n = pf_no
        !       !       exit
        !       !    endif
                 
        !       !    if(pf_p/pf_n-one.lt.1.D-8) then
        !       !       pf = 0.5*(pf_p+pf_n)
        !       !       print*,'no root found within bisection bounds',pf
        !       !       exit
        !       !    endif
        !       ! enddo
        !       ! ! endif
        !   !endif
           
        !    dpf = zero
        !    if(pf_p/pf_n-one.lt.1.D-8) then
        !       pf = 0.5*(pf_p+pf_n)
        !       exit
        !    endif

        ! enddo


        !original
        do ii = 1, itmax

           !bisection iteration
           pf = 0.5*(pf_p+pf_n)

           call this%fnumden_new(pf,fparams,iparams,num,den)
           call this%fnumden_new(pf_p,fparams,iparams,num_p,den_p)
           !if (sign(one,num).EQ.sign(one,num_p)) then
           if (num*num_p.GT.zero) then
              pf_p = pf
           else
              pf_n = pf
           endif

           dpf = zero
           if( abs(pf_p/pf_n)-one.lt.1.0d-15 ) exit

        enddo

        if(pf/ps_n-one.lt.1.0d-15 ) then
           icount2 = icount2 + one
           !if((icount2 .le. two).and.(isub.eq.one)) then
           if((icount2 .le. two)) then
              write(str,'(1(ES8.1E2,1X),3(I5,1X),(I1))') pf,iparams(2)+this%decomp%yst(1)-1,iparams(3)+this%decomp%yst(2)-1,iparams(4)+this%decomp%yst(3)-1,isub
              print*,'pressure converging to floor '//trim(str) !set lower intSharp_pfloor to correct -- default zero
              if(icount2 .eq. two) then
                 print*, 'not printing additional pressure converging to floor'
              endif
           endif

           ! if(isub.eq.nsubs) then
           !    stop
           ! endif

        endif

        !check for convergence failure
        if(ii==itmax+1) then
           write(*,*) 'Bisection method for pf did not converge'
           write(str,'(3(ES8.1E2,1X),4(I5,1X))') dpf,pf,den,iparams(2)+this%decomp%yst(1)-1,iparams(3)+this%decomp%yst(2)-1,iparams(4)+this%decomp%yst(3)-1
           write(*,*) 'Peq: dpf,pf,den,i,j,k '//trim(str)
           !print*,'pf,pf_n,pf_p,ps_n,ps_p',pf,pf_n,pf_p,max(pmin/pst,eps),pmax*pst
           print*,'pf,pf_n,pf_p,ps_n,ps_p',pf,pf_n,pf_p,ps_n,pmax*pst
           
           if( abs(pf_p/pf_n)-one.gt.1.D-3 ) then
              stop !very bad convergence
           endif

        endif

    end subroutine rootfind_nr_1d_new





!     subroutine rootfind_nr_1d_PT(this,pf,fparams,iparams,rhomin,pmin,pmax)
!         use decomp_2d, only: nrank
!         use reductions, only : P_MINVAL
!         class(solid_mixture), intent(inout)   :: this
!         real(rkind), intent(inout)            :: pf      !equilibrium pressure
!         real(rkind), intent(in), dimension(:) :: fparams !1: mixture hydrostatic (internal-elastic) energy
!         real(rkind), intent(in)               :: rhomin,pmin,pmax
!         integer, intent(in), dimension(:)     :: iparams !1: PTeqb flag, 2-4: i,j,k
    
!         integer     :: ii, imat, itmax = 1000,iter
!         !real(rkind) :: tol = 1.0d-8
!         real(rkind) :: dpf, num, den, den_conv, urlx,dpfold,pfold,pdefault,urlx0 = 1.01
!         character(len=100) :: str

!         pdefault = pf !initial guess -- volume fraction weighted

!         dpfold = 1.0/eps
!         pfold = zero

!         !pfinitguess = pf
!         do ii = 1, itmax
!           call this%fnumden(pf,fparams,iparams,num,den)

!           !check for small den: not necessary?
!           ! if(dabs(den)>1.0d-12) then
!           if(dabs(den)>1.0d-15) then
!              dpf = num/den !original
!           else
!             write(*,*) 'den very small, please check.', num, num/den
!             !write(*,*) 'failure at proc ', nrank, ' at index ', iparams(2:4)
!             write(*,*) 'interation ',ii,'failure at proc ', nrank, ' at index ', iparams(2)+this%decomp%yst(1)-1,iparams(3)+this%decomp%yst(2)-1,iparams(4)+this%decomp%yst(3)-1
!             stop
!           endif

!           !dpf = num/den
!           !!pf = pf - dpf !original
          
!           urlx0 = 1.01
!           urlx = urlx0


!           iter = zero
!           !do while ( (abs(pf).lt.pmin*1.5) .OR. (abs(pf).gt.pmax/1.5) ) !avoid extreme pressure oscillations
!           !do while ( (abs(pf).lt.pmin) .OR. (abs(pf).gt.pmax) .OR. (abs(dpf).gt.abs(dpfold)) ) !avoid extreme pressure oscillations from too large step size
!           do while ( (abs(dpf).gt.abs(dpfold)) ) !avoid extreme pressure oscillations from too large step size
!              iter = iter + one
!              !urlx = one+0.5*(urlx-one)
!              urlx = one+0.99*(urlx-one)
!              pf = max(min(pf-dpf,pf*urlx),pf/urlx)
!              ! check for convergence
!              if(dabs(pf)>1.0d-12) then
!                 den_conv = dabs(pf)
!              else
!                 den_conv = one
!              endif
!              if(dabs(dpf)/den_conv<1.0d-8) goto 10
             

!              call this%fnumden(pf,fparams,iparams,num,den)
!              dpf = num/den

!              !if( (urlx.lt.one+1.0D-12) ) then
!              if( (urlx.lt.one+1.0D-12) .OR. (abs(den).lt.1.0D-12) ) then !take the large outer step if: 1) no root is found by succesive underrelaxation .OR. 2) den->0 without passing over a root, need to jump over den=0
!                 ! dpf = dpfold
!                 ! pf = pfold

!                 dpf = dpfold*0.999
!                 pf = max(min(pf-dpf,pf*urlx),pf/urlx)

!                 exit!goto 10
!              endif
!           enddo
!           pfold = pf
!           dpfold = dpf


!           ! iter = zero
!           ! do while ( (abs(den).lt.rhomin*1.0D-2) ) !run away from zero density
!           !    iter = iter + one
!           !    urlx = one+2.0*(urlx-one)
!           !    pf = max(min(pf-dpf,pf*urlx),pf/urlx)
!           !    ! check for convergence
!           !    if(dabs(pf)>1.0d-12) then
!           !       den_conv = dabs(pf)
!           !    else
!           !       den_conv = one
!           !    endif
!           !    if(dabs(dpf)/den_conv<1.0d-8) goto 10
             

!           !    call this%fnumden(pf,fparams,iparams,num,den)
!           !    dpf = num/den

!           !    if(urlx.ge.two) goto 10
!           ! urlx  = urlx0
          
!           ! iter = zero
!           ! !do while(abs(dpf).gt.abs(1.1*dpfold))
!           ! !do while( (abs(dpf).gt.abs(1.1*dpfold)) .OR. (abs(dpf).gt.0.1*pf) )
!           ! do while( (abs(dpf).gt.abs(1.1*dpfold)) .OR. (abs(den).lt.rhomin*1.0D-2) )
!           !    iter = iter + one
!           !    urlx = one+0.5*(urlx-one)

!           !    pf = max(min(pf-dpf,pf*urlx),pf/urlx)
!           !    ! check for convergence
!           !    if(dabs(pf)>1.0d-12) then
!           !       den_conv = dabs(pf)
!           !    else
!           !       den_conv = one
!           !    endif
!           !    if(dabs(dpf)/den_conv<1.0d-8) goto 10 !exit
             

!           !    call this%fnumden(pf,fparams,iparams,num,den)
!           !    dpf = num/den
!           !    !if(urlx.lt.one+1.0D-6) then
!           !    if(urlx.lt.one+1.0D-12) then
!           !       write(*,*) 'Exiting Peq tight iteration: pf,dpf,pfold,dpfold,iter,num,den',pf,dpf,pfold,dpfold,iter,num,den

!           !       !pf = max(min(pf,pmax),pmin)
!           !       !pf = pdefault
!           !       pf = pfold
                
!           !       iter = zero
!           !       do while(urlx.lt.urlx0)
!           !          iter = iter + one
!           !          urlx = one+2.0*(urlx-one)

!           !          pf = max(min(pf-dpf,pf*urlx),pf/urlx)
!           !          ! check for convergence
!           !          if(dabs(pf)>1.0d-12) then
!           !             den_conv = dabs(pf)
!           !          else
!           !             den_conv = one
!           !          endif
!           !          if(dabs(dpf)/den_conv<1.0d-8) goto 10 !exit

!           !          call this%fnumden(pf,fparams,iparams,num,den)
!           !          dpf = num/den
!           !       enddo
!           !       !write(*,*) 'Exiting Peq loose iteration: pf,dpf,pfold,dpfold,iter,num,den',pf,dpf,pfold,dpfold,iter,num,den

!           !       exit
!           !    endif
!           ! enddo
!           ! pfold = pf
!           ! dpfold = dpf

!           pf = max(min(pf-dpf,pf*urlx),pf/urlx)

! 10        continue

!           ! check for convergence
!           if(dabs(pf)>1.0d-12) then
!             den_conv = dabs(pf)
!           else
!             den_conv = one
!           endif
!           !if(dabs(dpf)/den_conv<1.0d-8) exit
!           if(dabs(dpf)/den_conv<1.0d-8) exit
!         enddo
!         if(ii==itmax+1) then
!            !write(*,*) 'Newtons method for pf did not converge. Check details.', iparams(1)
!            !write(*,*) 'Newtons method for pf did not converge: den, urlx, i,j,k', den,urlx,iparams(2)+this%decomp%yst(1)-1,iparams(3)+this%decomp%yst(2)-1,iparams(4)+this%decomp%yst(3)-1
!         write(str,'(4(ES8.1E2,1X),4(I5,1X))') dpf,pf,den,1.0-urlx,iparams(2)+this%decomp%yst(1)-1,iparams(3)+this%decomp%yst(2)-1,iparams(4)+this%decomp%yst(3)-1
!         write(*,*) 'Peq: dpf,pf,den,urlx,i,j,k '//trim(str)
!         endif

!         ! !overwrite with bounds values
!         ! pf = max(min(pf,pmax),pmin)

!     end subroutine rootfind_nr_1d_PT


end module
